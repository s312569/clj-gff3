(ns clj-gff3.core
  (:require [clojure.string :refer [trim split]]
            [fs.core :refer [file?]]
            [clojure.set :refer [difference]]
            [clojure.java.io :as io]
            [clj-fasta.core :as fs]
            [bioreader.core :as br]
            [clojure.string :as st]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; data structures
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord gffEntry [seqid source type start end score strand phase attributes])
(defrecord gffDirective [name data])
(defrecord gffReader [resolved gff-strm fasta-strm]
  java.io.Closeable
  (close [this]
    (.close ^java.io.BufferedReader (:fasta-strm this))
    (.close ^java.io.BufferedReader (:gff-strm this))))

(defmulti gff-string (fn [e] (class e)))
(defmulti gff-directive (fn [n d] n))

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

(defn- parse-gff-entry
  [line]
  (map->gffEntry (zipmap [:seqid :source :type :start
                          :end :score :strand :phase
                          :attributes]
                         (let [f (split line #"\t")
                               a (->> (split (last f) #";")
                                      (map #(let [f (split % #"=")]
                                              (vector (keyword (st/lower-case (first f)))
                                                      (split (last f) #","))))
                                      (into {}))]
                           (conj (vec (butlast f)) a)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; entries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod gff-string gffEntry
  [e]
  (let [k [:seqid :source :type :start :end :score :strand
           :phase]
        a (->> (map (fn [[k v]]
                      (str (name k) "=" (->> v (interpose ",") (apply str))))
                    (:attributes e))
               (interpose ";")
               (apply str))
        l (->> (map #(get e %) k)
               (interpose \tab)
               (apply str))]
    (str l \tab a)))

(defn start
  "Returns the start field of a gffEntry as an integer."
  [e]
  (Integer/parseInt (:start e)))

(defn end
  "Returns the end field of a gffEntry as an integer."
  [e]
  (Integer/parseInt (:end e)))

(defn score
  "Returns the score field as a float or nil if none."
  [e]
  (if-not (= (:score e) ".")
    (Float/parseFloat (:score e))))

(defn phase
  "Returns the phase as an integer (0, 1 or 2) or nil if none."
  [e]
  (if-not (= (:phase e) ".")
    (Integer/parseInt (:phase e))))

(defn seqid
  "Returns the seqid of an entry."
  [e]
  (:seqid e))

(defn source
  "Returns the source of an entry."
  [e]
  (:source e))

(defn type
  "Returns the type of an entry."
  [e]
  (:type e))

(defn strand
  "Returns the strand of an entry."
  [e]
  (:strand e))

(defn cds?
  "Returns true if 'type' field of gffEntry equals 'CDS'."
  [e]
  (= "CDS" (:type e)))

(defn gene?
  "Returns true if 'type' field of gffEntry equals 'gene'."
  [e]
  (= "gene" (:type e)))

(defn exon?
  [e]
  (= (:type e) "exon"))

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

(defn get-attribute
  "Returns an attribute value corresponding to the key
  argument. Always returns a vector."
  [entry key]
  (get-in entry [:attributes key]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; directives
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod gff-string gffDirective
  [e]
  (if (= (:name e) :resolved)
    "###"
    (str "##" (:name e) \tab (:data e))))

(defn directive? [e]
  (instance? gffDirective e))

(defn resolved?
  "Returns true if directive is the resolved directive."
  [e]
  (= :resolved (:name e)))

(defn seq-region?
  "Returns a vector of seqid, start and end if directive is a
  seq-region directive. False if not."
  [d]
  (if (= "sequence-region" (:name d))
    (let [f (st/split (:data d) #"\t")]
      (vector (first f) (Integer/parseInt (second f)) (Integer/parseInt (last f))))
    false))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; reader
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod fs/fasta-line-seq gffReader
  [r]
  (line-seq (:fasta-strm r)))

(defn gff-reader
  "Returns a gffReader that can be used with fasta-seq and gff-seq."
  [gff-file]
  (->gffReader (atom nil) (br/bioreader gff-file) (br/bioreader gff-file)))

(defn gff-seq
  "Lazy sequence of all entries and directives in a gff3 file."
  [gffreader]
  (->> (line-seq (:gff-strm gffreader))
       (take-while #(not (or (= % "##FASTA")
                             (= \> (first %)))))
       (filter #(not (or (= "" (trim %))
                         (re-find #"^\#[^#]" %))))
       (map #(let [l (trim %)]
               (cond (= "###" l)
                     (do (reset! (:resolved gffreader) true)
                         (->gffDirective :resolved nil))
                     (= '(\# \#) (take 2 l))
                     (let [[n d] (rest (re-find #"\#\#([^\s]+)\s+(.+)" l))]
                       (->gffDirective n d))
                     :else
                     (parse-gff-entry l))))))

(defn reader-resolved?
  "Returns true if all forward references in the gff file have been resolved."
  [gff-reader]
  (:resolved gff-reader))

