# clj-gff3

A gff3 parsing library for Clojure

## Installation

To get from Clojars put the following in your project.clj file:

```clojure
[clj-gff "0.1.0"]
```

To use require in your app:

```clojure
(ns my-app.core
    (:require [clj-gff.core :as gf]))
```

API docs can be found [here](http://s312569.github.io/clj-gff3/)

## Usage

Firstly initialise a gff file (only gff3 is supported):

```clojure
user> (def tf (init-gff-file "/path/to/test-file.gff"))
```

Then open using the `gff-reader` function and access entries using
`gff-seq` which provides a lazy list of `gffEntry` and `gffDirective`
records:

```clojure
user> (with-open [r (gff-reader tf)]
                 (doall (take 3 (gff-seq r))))
(#clj_gff3.core.gffDirective{:name "gff-version", :data "3"}
#clj_gff3.core.gffDirective{:name "sequence-region", :data "KI657455 1
1890151"} #clj_gff3.core.gffEntry{:seqid "KI657455", :source
"WormBase_imported", :type "gene", :start "25", :end "387", :score
".", :strand "-", :phase ".", :attributes {:ID ["gene:NECAME_00001"],
:Name ["NECAME_00001"], :Note ["NOVEL protein_coding"], :Alias
["NECAME_00001"]}})
```

`gffDirective` records correspond to lines beginning with `##`. When
filtering entries the `resolved?` function returns true when the
resolved directive (###) is encountered and works on both `gffEntry`
and `gffDirective` records.

`gffEntry` records are the entries, comment lines (beginning with #)
are ignored. `gffEntry` records have keys corresponding to the defined
fields of the gff3 spec (e.g. :seqid :source :type :start
etc). Attributes can be accessed using the :attributes keyword or the
`get-attribute` function which takes a `gffEntry` and an attribute
name (e.g. 'ID', 'Name' etc).

## Fasta sequences

The gff3 format allows fasta sequences to follow annotation entries
and the functions `fasta-seq` and `get-fasta` provide access to these
sequences. The former returns a lazy list of all fasta sequences in
the file and the latter gets a specific fasta sequence. They both
return a fastaSequence record so you will need to include
`clj-biosequence.core` to do stuff with them (see
[clj-biosequence](https://github.com/s312569/clj-biosequence)).

```clojure
user> (require '[clj-biosequence.core :as bs])
nil
user> (with-open [r (gff-reader tf)]
        (first (fasta-seq r)))	
#clj_biosequence.core.fastaSequence{:acc "sp|C1IC49|3FN5_WALAE",
 :description "Three finger toxin Wa-V OS=Walterinnesia aegyptia PE=1
 SV=1", :alphabet :iupacAminoAcids, :sequence [\M \K \T \L \L \L \T \L
 \V \L \V \T \I \V \C \L \D \L \G \Y \T \L \T \C \L \I \C \P \K \K \Y
 \C \N \Q \V \H \T \C \R \N \G \E \N \L \C \I \K \T \F \Y \E \G \N \L
 \L \G \K \Q \F \K \R \G \C \A \A \T \C \P \E \A \R \P \R \E \I \V \E
 \C \C \S \R \D \K \C \N \H]}
user> (with-open [r (gff-reader tf)]
        (get-fasta r "sp|C1IC49|3FN5_WALAE"))
#clj_biosequence.core.fastaSequence{:acc "sp|C1IC49|3FN5_WALAE",
:description "Three finger toxin Wa-V OS=Walterinnesia aegyptia PE=1
SV=1", :alphabet :iupacAminoAcids, :sequence [\M \K \T \L \L \L \T \L
\V \L \V \T \I \V \C \L \D \L \G \Y \T \L \T \C \L \I \C \P \K \K \Y
\C \N \Q \V \H \T \C \R \N \G \E \N \L \C \I \K \T \F \Y \E \G \N \L
\L \G \K \Q \F \K \R \G \C \A \A \T \C \P \E \A \R \P \R \E \I \V \E
\C \C \S \R \D \K \C \N \H]}
```

## License

Copyright © 2015 Jason Mulvenna

Distributed under the Eclipse Public License either version 1.0.
