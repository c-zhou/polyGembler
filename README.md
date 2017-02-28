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

Both of them will generate an executable jar file *polyGembler-${version}-jar-with-dependencies.jar*. For [Apache Ant](http://ant.apache.org/), the jar file is under *polyGembler/dist* directory. For [Apache Maven](https://maven.apache.org/), the jar file is under *polyGembler/target* directory.

The code has been tested with Java 7 and 8. The default version is Java 8. The [Apache Maven](https://maven.apache.org/) will fail if the system default Java version is not 8. This could be solved by either changing the system Java version or modifying the *pom.xml* file at line 22-23.

## Quick Start
A quick start to run the software. Modules listed below are independent from each other.

#### 1. Simulate GBS data for an outcrossed F1 mapping population
*Input*&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;a reference fasta file with all chromosomes

*Output*&nbsp;&nbsp;&nbsp;a zipped fastq file with all GBS reads

*Command*

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar popsimulation -r ${reference.fa} -t 32 -p 4 -c 200 -o ${pop_out_dir} -n 192
    $ java -jar polyGembler-${version}-jar-with-dependencies.jar gbssimulation -f ${pop_out_dir}/Sc1 -t 32 -m 5 -s 5 -o ${gbs_out_dir}

The first step simulates a tetraploid (-p option) F1 mapping population with 192 samples (-n option) from the reference genome ${reference.fa} (-r option). The total genetic length of the chromosomes is 200cm (-c option). It uses 32 CPUs (-t option). The output sample genomes (FASTA files) will be in ${pop_out_dir}/Sc1 (-o option). The second step takes the simulated mapping population genomes in the first step as input (-f option). The average sequencing depth of coverage for each copy of the chromosome is 5 (-m option) and the standard deviation is 5 (-s option). It uses 32 CPUs (-t option). The output GBS data (FASTQ file) will be in ${gbs_out_dir} (-o option). The program write GBS reads for each sample seperately as a [gzipped](https://en.wikipedia.org/wiki/Gzip) file. If you need to put all GBS reads together, simply use the [cat](https://en.wikipedia.org/wiki/Cat_(Unix)) command.

    $ cat ${gbs_out_dir}/*.gz > ${merged_gbs_file}.gz

#### 2. Run variant detection module for GBS data

*Input*&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;GBS FASTQ files, GBS key file, reference genome/assembly

*Output*&nbsp;&nbsp;&nbsp;VCF file

*Command*

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar gbspileup -i ${gbs_fastq_dir} -k ${gbs_key_file} -p 2 -t 32 -f ${assembly.fa} -o ${vd_out_dir}

To run this command, the following tools should be installed and added to system path.

* [samtools](http://samtools.sourceforge.net/)
* [bwa](http://bio-bwa.sourceforge.net/)
* [freebayes](https://github.com/ekg/freebayes)

This command runs the whole variant detection pipeline from the GBS reads demultiplexing to variant calling with [freebayes](https://github.com/ekg/freebayes). The GBS FASTQ file(s) should be found by the program in directory ${gbs_fastq_dir} (-i option). The GBS key file is provided with -k option. [Here](https://bytebucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline/Pipeline_Testing_key.txt?rev=19fcce52374296f80ed568fbdd752f8cffabb2ae) is a sample key file from [Buckler Lab](http://www.maizegenetics.net/tassel). The ploidy of the genome is specified with -p option. The reference genome assembly is provided with -f option. The output will be in ${vd_out_dir} directory. It should noted that the program will create ${vd_out_dir} if it is not existed. 

You may skip the freebayes variant calling step with -z option, and instead use other variant detection pipeline such as [samtools mpileup](http://samtools.sourceforge.net/mpileup.shtml) and [GATK](https://software.broadinstitute.org/gatk/). The BAM files for each sample can be found under the ${vd_out_dir}/bam directory. Running with -z option will not require [freebayes](https://github.com/ekg/freebayes) installed.
    
    $ java -jar polyGembler-${version}-jar-with-dependencies.jar gbspileup -i ${gbs_fastq_dir} -k ${gbs_key_file} -p 2 -t 32 -f ${assembly.fa} -o ${vd_out_dir} -z

#### 3. Run genetic linkage map or pseudomolecule construnction pipeline

*Input*&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;VCF file, reference genome assembly

*Output*&nbsp;&nbsp;&nbsp;genetic linkage maps or pseudomolecules

*Command*

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar gbspileup -i ${in_vcf_file} -a ${assembly_fasta_file} -f ${parent_sample_1}:${parent_sample_2} -p 2 -t 32 -o ${map_out_dir}
    
The command takes the VCF file (-i option) and the assembly fasta file (-a option) as input. Two parental samples of the full-sib family should be specified with -f option and use ":" as delimiter. This is a diploid genome so set -p option as default 2. It uses 32 CPUs. The output files will be found under ${map_out_dir} directory (-o option).

#### 4. Run haplotype phasing algorithm

*Input*&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;preprocessed file

*Output*&nbsp;&nbsp;&nbsp;inferred haplotypes for each sample

*Command*

Prepare a VCF file for haplotype phasing.

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar datapreparation -i ${in_vcf_file} -p 2 -q 10 -f 0.1 -m 0.5 -u 3000 -o ${out_zip_dir}
    
The VCF file is provided with -i option. The program will filter the variants according to the parameters specified by user. Here the program is told that it is dealing with a diploid genome (-p option). The minimum quality score of a base is set 10 (-q option). Bases with quility scores lower than 10 will be treated as missing. The lower bound of minor allele frequency is set to 0.1. The rare variants called from a F1 full-sib family is very likely caused by sequencing errors. The maximum missing data rate across a locus is set to 0.5. Loci with more than 50% missing genotypes will be filtered out. The upper bound of total allele depth (DP field in VCF files) for a position is set to 3000 (-u option). If the DP field exceeds 3000, the variant will be filtered out. This is used to remove ambiguous variants caused by copy number variantion. It should be noted here that the program only deals with bi-allelic SNPs currently. The output is a zipped file with the same prefix as the input VCF file. The location of the output file is specified with -o option.
    
Run the haplotype phasing algorithm.

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar haplotyper -i ${in_zip_file} -c ${contig_str_id} -p 2 -f ${parent_sample_1}:${parent_sample_2} -L -o ${out_file_dir}
    
The program takes the zipped file generated in the preprocess step and a contig/scaffold string id as input. User also need to specify the ploidy, the two parental sample names. -L option tells the program to use the phred-scaled likelihood scores, which is default. Otherwise user could run with -D option which utilises allele depth or -G option which utilises genotype information. The output file directory is specified with -o option.

You can run multiple contigs/scaffolds simultaneously. 

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar haplotyper -i ${in_zip_file} -c ${contig_str_id}:${contig_str_id2}:${contig_str_id3} -r true:false:false -s 0.05:0.1 -p 2 -f ${parent_sample_1}:${parent_sample_2} -L -o ${out_file_dir}

This is very simular to running with the single contig/scaffold. Mutiple contigs/scaffolds are provided with string ids seperated by ":". As there could be different concatenation directions, -r option specifies if each contig/scaffold is reversed or not. The default is not. The distances between the adjacent contigs/scaffolds are initialised with -s option, otherwise the program will generate them randomly. The distances could either be the recombination frequencies or the physical distances. The program will check these numbers. If all of them are smaller than 0.5, then they will be taken as recombination frequencies.

## Parameter options
#### 1. Simulate GBS data for an outcrossed F1 mapping population
#### 2. Run variant detection module for GBS data
#### 3. Run genetic linkage map or pseudomolecule construnction pipeline
#### 4. Run haplotype phasing algorithm

## Details about the output files
#### 1. Simulate GBS data for an outcrossed F1 mapping population
#### 2. Run variant detection module for GBS data
#### 3. Run genetic linkage map or pseudomolecule construnction pipeline
#### 4. Run haplotype phasing algorithm

## Citing polyGembler
