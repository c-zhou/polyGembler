# *polyGembler*: pseudomolecule construction combining *de novo* assembly and genetic mapping
## user manual and guide

--------------------

## Overview
*polyGembler* is program that constructs genetic linakge maps or chromosomal-scale pseudomolecules combining *de novo* assembly and genetic mapping. The method assumes availability of genome-wide genotyping data such as [GBS](https://en.wikipedia.org/wiki/Genotyping_by_sequencing) and array data, collected on a F1 outbred family, as well as high coverage (i.e. greater than 30X) whole genome sequence data on a reference sample, or alternatively the availability of a set of reference contigs or scaffolds. By mapping marker set to contigs *polyGembler* infers contig haplotypes for each sample. Contig haplotypes are then used to infer linkage groups corresponding to chromosomes as well as the optimal ordering of contigs within these chromosomes. *polyGembler* consists of three major modules, namely *variant detection*, *recombination frequency estimation* and *genetic mapping*.

<img src="https://github.com/c-zhou/polyGembler/raw/master/pipeline_flowchart.png" width=600/>

*polyGembler* also provides a tool for simulating outcrossed F1 mapping population GBS data. It uses the software [PedigreeSim V2.0](https://www.wur.nl/en/show/Software-PedigreeSim.htm) to simulate the full-sib family genomes and imitates the [GBS protocal](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0019379) to generate GBS reads.

*polyGembler* is written in Java. It has a user friendly design. All the tasks could be done with an one-liner command require no further user interactions. It has been multithreaded where possible.


## Installation

#### Download from the releases
The releases are executable jar files with all dependencies, which could be run directly.

#### Compile from source code
Java build tool either [Apache Ant](http://ant.apache.org/) or [Apache Maven](https://maven.apache.org/) is required.

Quick installation guide:

    $ git clone https://github.com/c-zhou/polyGembler.git
    $ cd polyGembler
    
    $ ant
    or
    $ mvn clean package

Both of them will generate a executable jar file *polyGembler-${version}-jar-with-dependencies.jar*. For [Apache Ant](http://ant.apache.org/), the jar file is under *polyGembler/dist* directory. For [Apache Maven](https://maven.apache.org/), the jar file is under *polyGembler/target* directory.
