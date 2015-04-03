(ns clj-gff3.core
  (:require [clj-biosequence.core :as bs]
            [clojure.string :refer [trim split]]
            [fs.core :refer [file?]]
            [clojure.set :refer [difference]]))

(declare gather-consecutive)

(defprotocol gffResolved
  (resolved? [this]
    "Returns true if entry denotes all forward references are
    resolved."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; entries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord gffEntry [seqid source type start
                     end score strand phase
                     attributes]
  gffResolved
  (resolved? [this] false))

(defn start
  "Returns the start field of a gffEntry as an integer."
  [e]
  (Integer/parseInt (:start e)))

(defn end
  "Returns the end field of a gffEntry as an integer."
  [e]
  (Integer/parseInt (:end e)))

(defn cds?
  "Returns true if 'type' field of gffEntry equals 'CDS'."
  [e]
  (= "CDS" (:type e)))

(defn gene?
  "Returns true if 'type' field of gffEntry equals 'gene'."
  [e]
  (= "gene" (:type e)))

(defn entry-range 
  "Returns the range of an entry as a set of integers."
  [e]
  (set (range (start e) (+ 1 (end e)))))

(defn entry-difference
  "Takes an entry and a collection of entries and provides a list of
  vectors of start and end values for ranges of the first entry that
  are not found in the list of entries."
  [e & es]
  (->> (apply difference (entry-range e) (map entry-range es))
       gather-consecutive
       (map #(vector (apply min %) (apply max %)))))

(defn- parse-gff-entry
  [line]
  (map->gffEntry (zipmap [:seqid :source :type :start
                          :end :score :strand :phase
                          :attributes]
                         (let [f (split line #"\t")
                               a (->> (split (last f) #";")
                                      (map #(let [f (split % #"=")]
                                              (vector (keyword (first f))
                                                      (split (last f) #","))))
                                      (into {}))]
                           (conj (vec (butlast f)) a)))))

(defn get-attribute
  "Returns an attribute value corresponding to the key
  argument. Always returns a vector."
  [entry key]
  (get-in entry [:attributes key]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; directives
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord gffDirective [name data]
  gffResolved
  (resolved? [this] (= :resolved (:name this))))

(defn seq-region?
  "Returns true if gffDirective is a sequence-region directive."
  [d]
  (= "sequence-region" (:name d)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; reader
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- gff-seq-helper
  [line]
  (let [l (trim line)]
    (cond (= "###" l)
          (->gffDirective :resolved nil)
          (= '(\# \#) (take 2 l))
          (let [[n d] (rest (re-find #"\#\#([^\s]+)\s+(.+)" l))]
            (->gffDirective n d))
          :else
          (parse-gff-entry l))))

(defrecord gffReader [strm fastas]
  bs/biosequenceReader
  (biosequence-seq [this]
    (if (:fastas this)
      (bs/biosequence-seq (:fastas this))
      (throw (Throwable. "No sequence stream. If file contains sequences initialise with 'init-gff-file' with :alphabet keyword."))))
  (get-biosequence [this accession]
    (if (:fastas this)
      (bs/get-biosequence (:fastas this) accession)
      (throw (Throwable. "No sequence stream. If file contains sequences initialise in 'init-gff-file' with :alphabet keyword."))))
  java.io.Closeable
  (close [this]
    (.close ^java.io.BufferedReader (:strm this))
    (if (:fastas this)
      (.close (:fastas this)))))

(defn gff-seq
  "Lazy sequence of all entries and directives in a gff3 file."
  [gffreader]
  (->> (line-seq (:strm gffreader))
       (take-while #(not (or (= % "##FASTA")
                             (= \> (first %)))))
       (filter #(not (or (= "" (trim %))
                         (re-find #"^\#[^#]" %))))
       (map gff-seq-helper)))

(defn- tokenise
  [m]
  (let [l (:remaining m)
        r (drop-while #(not (or (gene? %)
                                (resolved? %)))
                      (rest l))
        y (let [es (drop-while #(not (or (gene? %)
                                         (seq-region? %))) l)
                fi (take-while #(or (gene? %)
                                    (seq-region? %)) es)]
            (concat fi
                    (take-while #(not (or (gene? %)
                                          (resolved? %)))
                                (rest es))))]
    (if (seq y)
      {:yield y :remaining r}
      {:end true})))

(defn gene-seq
  "Returns a lazy list of vectors containing entries describing
  genes."
  [gffreader]
  (let [e (drop-while #(not (= "gene" (:type %)))
                      (gff-seq gffreader))]
    (->> {:remaining e}
         (iterate tokenise)
         rest
         (take-while #(not (contains? % :end)))
         (map #(vec (:yield %))))))

(defn get-fasta
  "Returns a fastaSequence (see clj-biosequence) from a gff3 file
  containing fasta sequences."
  [gffreader accession]
  (bs/get-biosequence gffreader accession))

(defn fasta-seq
  "Returns a lazy list of fasta sequence records (see clj-biosequence)
  in a gff file."
  [reader]
  (bs/biosequence-seq reader))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord gffFile [file alphabet opts])

(extend gffFile
  bs/biosequenceFile
  bs/default-biosequence-file)

(defn gff-reader
  [gff-file]
  (let [i (if (:alphabet gff-file)
            (bs/bs-reader (bs/index-biosequence-file
                           (apply bs/init-fasta-file
                                  (bs/bs-path gff-file)
                                  (:alphabet gff-file)
                                  (:opts gff-file)))))]
    (->gffReader (apply bs/bioreader (bs/bs-path gff-file)
                        (:opts gff-file))
                 i)))

(defn init-gff-file
  "Initialises a gffFile. If the keyword argument is specified with an
  argument compatible with the clj-biosequence package then the reaser
  can also provide a stream of fasta sequences that are in the
  file (see clj-biosequence)."
  [path & opts]
  (if (file? path)
    (->gffFile path (:alphabet (apply hash-map opts))
               (concat '(:junk true) opts))
    (throw (Throwable. (str "File not found: " path)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- gather-consecutive
  [lst]
  (let [f (atom true)]
    (loop [l (sort lst)
           curr ()
           acc ()]
      (if-not (second l)
        (if (= (- (first l) (first curr)) 1)
          (cons (reverse (cons (first l) curr)) acc)
          (cons (list (first l)) (cons (reverse curr) acc)))
        (cond @f
              (do (reset! f false)
                  (recur (rest l) (list (first l)) acc))
              (= (- (first l) (first curr)) 1)
              (recur (rest l) (cons (first l) curr) acc)
              :else
              (recur (rest l) (list (first l)) (cons (reverse curr) acc)))))))
