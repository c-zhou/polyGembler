# *polyGembler*: pseudomolecule construction combining *de novo* assembly and genetic mapping
## user manual and guide

--------------------

## Overview
*polyGembler* is program that constructs genetic linakge maps or chromosomal-scale pseudomolecules combining *de novo* assembly and genetic mapping. The method assumes availability of genome-wide genotyping data such as [GBS](https://en.wikipedia.org/wiki/Genotyping_by_sequencing) and array data, collected on a F1 outbred family, as well as high coverage (i.e. greater than 30X) whole genome sequence data on a reference sample, or alternatively the availability of a set of reference contigs or scaffolds. By mapping marker set to contigs *polyGembler* infers contig haplotypes for each sample. Contig haplotypes are then used to infer linkage groups corresponding to chromosomes as well as the optimal ordering of contigs within these chromosomes. PolyGembler consists of three major modules, namely *variant detection*, *recombination frequency estimation* and *genetic mapping*.

*polyGembler* also provides a tool for simulating outcrossed F1 mapping population GBS data. It uses the software [PedigreeSim V2.0](https://www.wur.nl/en/show/Software-PedigreeSim.htm) to simulate the full-sib family genomes and imitates the [GBS protocal](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0019379) to generate GBS reads.

*polyGembler* is written in Java. It is design to minimize the user interaction. 
