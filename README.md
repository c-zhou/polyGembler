# *polyGembler*: pseudomolecule construction combining *de novo* assembly and genetic mapping
## user manual and guide

--------------------

## Overview
*polyGembler* is program that constructs genetic linakge maps or chromosomal-scale pseudomolecules combining *de novo* assembly and genetic mapping. The method assumes availability of genome-wide genotyping data such as [GBS](https://en.wikipedia.org/wiki/Genotyping_by_sequencing) and array data, collected on a F1 outbred family, as well as high coverage (i.e. greater than 30X) whole genome sequence data on a reference sample, or alternatively the availability of a set of reference contigs or scaffolds. By mapping marker set to contigs *polyGembler* infers contig haplotypes for each sample. Contig haplotypes are then used to infer linkage groups corresponding to chromosomes as well as the optimal ordering of contigs within these chromosomes. *polyGembler* consists of three major modules, namely *variant detection*, *recombination frequency estimation* and *genetic mapping*.

<img src="https://github.com/c-zhou/polyGembler/raw/master/pipeline_flowchart.png" width=600/>

*polyGembler* also provides a tool for simulating outcrossed F1 mapping population GBS data. It uses the software [PedigreeSim V2.0](https://www.wur.nl/en/show/Software-PedigreeSim.htm) to simulate the full-sib family genomes and imitates the [GBS protocol](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0019379) to generate GBS reads.

*polyGembler* is written in Java. It has a user friendly design. All the tasks could be done with an one-liner command require no further manual intervention. It has been multithreaded where possible.


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
*Input*&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;a reference FASTA file with all chromosomes

*Output*&nbsp;&nbsp;&nbsp;a zipped FASTQ file with all GBS reads

*Command*

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar popsimulation -r ${reference.fa} -t 32 -p 4 -c 200 -o ${pop_out_dir} -n 192
    $ java -jar polyGembler-${version}-jar-with-dependencies.jar gbssimulation -f ${pop_out_dir}/Sc1 -t 32 -m 5 -s 5 -o ${gbs_out_dir}

The first step simulates a tetraploid (-p option) F1 mapping population with 192 samples (-n option) from the reference genome ${reference.fa} (-r option). The total genetic length of the chromosomes is 200cm (-c option). It uses 32 CPUs (-t option). The output sample genomes (FASTA files) will be in ${pop_out_dir}/Sc1 (-o option). The second step takes the simulated mapping population genomes in the first step as input (-f option). The average sequencing depth of coverage for each copy of the chromosome is 5 (-m option) and the standard deviation is 5 (-s option). It uses 32 CPUs (-t option). The output GBS data (FASTQ file) will be in ${gbs_out_dir} (-o option). The program write GBS reads for each sample separately as a [gzipped](https://en.wikipedia.org/wiki/Gzip) file. If you need to put all GBS reads together, simply use the [cat](https://en.wikipedia.org/wiki/Cat_(Unix)) command.

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

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar gembler -i ${in_vcf_file} -a ${assembly_fasta_file} -f ${parent_sample_1}:${parent_sample_2} -p 2 -t 32 -o ${map_out_dir}
    
The command takes the VCF file (-i option) and the assembly FASTA file (-a option) as input. Two parental samples of the full-sib family should be specified with -f option and use ":" as delimiter. This is a diploid genome so set -p option as default 2. It uses 32 CPUs. The output files will be found under ${map_out_dir} directory (-o option).

#### 4. Run haplotype phasing algorithm

*Input*&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;preprocessed file

*Output*&nbsp;&nbsp;&nbsp;inferred haplotypes for each sample

*Command*

Prepare a VCF file for haplotype phasing.

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar datapreparation -i ${in_vcf_file} -p 2 -q 10 -f 0.1 -m 0.5 -u 3000 -o ${out_zip_dir}
    
The VCF file is provided with -i option. The program will filter the variants according to the parameters specified by user. Here the program is told that it is dealing with a diploid genome (-p option). The minimum quality score of a base is set 10 (-q option). Bases with quality scores lower than 10 will be treated as missing. The lower bound of minor allele frequency is set to 0.1. The rare variants called from a F1 full-sib family is very likely caused by sequencing errors. The maximum missing data rate across a locus is set to 0.5. Loci with more than 50% missing genotypes will be filtered out. The upper bound of total allele depth (DP field in VCF files) for a position is set to 3000 (-u option). If the DP field exceeds 3000, the variant will be filtered out. This is used to remove ambiguous variants caused by copy number variation. It should be noted here that the program only deals with biallelic SNPs currently. The output is a zipped file with the same prefix as the input VCF file. The location of the output file is specified with -o option.
    
Run the haplotype phasing algorithm.

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar haplotyper -i ${in_zip_file} -c ${contig_str_id} -p 2 -f ${parent_sample_1}:${parent_sample_2} -L -o ${out_file_dir}
    
The program takes the zipped file generated in the preprocess step and a contig/scaffold string id as input. User also need to specify the ploidy, the two parental sample names. -L option tells the program to use the phred-scaled likelihood scores, which is default. Otherwise user could run with -D option which utilises allele depth or -G option which utilises genotype information. The output file directory is specified with -o option.

You can run multiple contigs/scaffolds simultaneously. 

    $ java -jar polyGembler-${version}-jar-with-dependencies.jar haplotyper -i ${in_zip_file} -c ${contig_str_id}:${contig_str_id2}:${contig_str_id3} -r true:false:false -s 0.05:0.1 -p 2 -f ${parent_sample_1}:${parent_sample_2} -L -o ${out_file_dir}

This is very similar to running with the single contig/scaffold. Multiple contigs/scaffolds are provided with string ids separated by ":". As there could be different concatenation directions, -r option specifies if each contig/scaffold is reversed or not. The default is not. The distances between the adjacent contigs/scaffolds are initialised with -s option, otherwise the program will generate them randomly. The distances could either be the recombination frequencies or the physical distances. The program will check these numbers. If all of them are smaller than 0.5, then they will be taken as recombination frequencies.

## More parameter options
#### 1. Six main pipelines, each triggered by the executable jar file
<pre>
popsimulation                    Simulate a full-sib mapping population.
gbssimulation                    Simulate GBS data.
gbspileup                        Variant calling from GBS data.
datapreparation                  Prepare data for haplotype phasing.
haplotyper                       Contig/scaffold haplotype construction from a mapping population.
gembler                          Run PolyGembler pipeline to construct genetic linkage maps/pseudomolecules.
</pre>

#### 2. Simulate an outcrossed F1 mapping population (popsimulation)
<pre>
-r/--reference                   Reference (fasta file format).
-n/--pop-size                    Population size including parents (default 96).
-p/--ploidy                      Copy number of chromosomes (default 2).
-c/--centimorgan                 Total genetic length to simulate. Assume the physical length
                                 and genetic length has linear correlation (default calculated
                                 from the reference, 1cM per 1Mbp).
-t/--threads                     Number of threads (default 1).
-s/--run-id                      Unique run id (default Sc1).
-o/--prefix                      Output directory (default current directory).
</pre>

#### 3. Simualte GBS data from the F1 mapping population (gbssimulation)
<pre>
-f/--fasta-file                  Directory contains genome fasta files to be sequenced. 
-e/--enzyme                      Enzyme(default PstI). 
-l/--library                     GBS protocol library preparation file (default null). 
-t/--threads                     Number of threads (default 1).
-b/--barcode-file                GBS protocol barcode file (default null).
-m/--avg-depth                   Depth of coverage (default 5).
-s/--sdev                        Standard deviation of depth of coverage (default 5).
-S/--random-seed                 Random seed (default system nano time).
-q/--quality-file                Markov chain parameter file for quality scores (default null). 
-o/--output-prefix               Output directory (default current directory).
</pre>

#### 4. Run variant detection module for GBS data (gbspileup)
<pre>
-i/--input-fastq                 Input directory containing FASTQ files in text or gzipped text.
                                 NOTE: Directory will be searched recursively and should
                                 be written WITHOUT a slash after its name.
-k/--key-file                    Key file listing barcodes distinguishing the samples.
-e/--enzyme                      Enzyme used to create the GBS library, if it differs from the one 
                                 listed in the key file.
-q/--min-qualS                   Minimum quality score (default is 10).
-p/--ploidy                      Ploidy for variant calling (default is 2).
                                 NOTE: You may call variant as diploid and the program will
                                 fit a binomial model to call genotypes and genotype
                                 qualities from allele depth.
-t/--threads                     Threads (default is 1).
-T/--trim-leading                The length of leading fragments to trim off.
-f/--reference                   The reference genome (in fasta format).
-z/--skip-freebayes              Skip the variant calling with freebayes (default not).
-o/--prefix                      Output directory to contain .cnt files (one per FASTQ file, defaults 
                                 to input directory).
</pre>

#### 5. Run data preparation for haplotype phasing algorithm (datapreparation)
<pre>
-i/--vcf		         Input VCF file.
-s/--id                          Unique id of this run (default: input VCF file name prefix).
-p/--ploidy			 Ploidy of genome (default 2).
                                 NOTE: If you called variant as diploid, then the program will
                                 fit a binomial model to call genotypes and genotype qualities
                                 from allele depth with the ploidy specified here.
-l/--min-depth		         Minimum depth to keep a SNP (DP).
-u/--max-depth		         Maximum depth to keep a SNP (DP).
-q/--min-qual  		         Minimum quality to keep a SNP (QUAL).
-f/--min-maf		         Minimum minor allele frequency to keep a SNP (default 0.1).
-m/--max-missing	         Maximum proportion of missing data to keep a SNP (default 0.5).
-o/--prefix			 Prefix for output files (default: input VCF file folder).
</pre>

#### 6. Run haplotype phasing algorithm (haplotyper)
<pre>
-i/--input                       Input zipped file.
-o/--prefix                      Output file location.
-ex/--experiment-id              Common prefix of haplotype files for this experiment.
-c/--scaffold                    The scaffold/contig/chromosome id will run.
-x/--max-iter                    Maxmium rounds for EM optimization (default 100).
-p/--ploidy                      Ploidy of genome (default 2).
-f/--parent                      Parent samples (separated by a \":\").
-s/--initial-separation          Initialisations of distances between the adjacent scaffolds 
                                 if multiple scaffolds will be jointly inferred. The separation
                                 could be either physical distances or recombination frequencies, 
                                 i.e., if all values provided is below 0.5, the
                                 program will take them as recombination frequencies.
                                 Distances should be separated by \":\".
-r/--reverse                     Take either 'true' or 'false', indicating whetherr the
                                 scaffold is reversed before inferring haplotypes. Multiple
                                 scaffolds are separated by \":\".
-G/--genotype                    Use genotypes to infer haplotypes. Mutually exclusive with
                                 option -D/--allele-depth and -L/--genetype likelihood.
-D/--allele-depth                Use allele depth to infer haplotypes. Mutually exclusive
                                 with option -G/--genotype and -L/--genetype likelihood.
-L/--genotype-likelihood         Use genotype likelihoods to infer haplotypes. Mutually
                                 exclusive with option -G/--genotype and -L/--allele-depth 
                                 (default).
-b/--segmental-kmeans            Use Viterbi training instead of Baum-Welch algorithm.
-e/--train-exp                   Re-estimate transition probabilities between founder/parental
                                 haplotypes at each step.
-S/--random-seed                 Random seed for this run.
-pp/--print-plot                 Plot the hidden Markov model.
-sp/--save-plot                  Save the plot as a pdf file. The file name should be provided here.
</pre>

#### 7. Run genetic linkage map or pseudomolecule construction pipeline (gembler)
<pre>
Common:
-i/--input-vcf                   Input VCF file.
-o/--prefix                      Output file location, create the directory if not exist.
-p/--ploidy                      Ploidy of genome (default 2).
-S/--random-seed                 Random seed for this run.
-t/--threads                     Threads (default 1).
-rlib/--R-external-libs          External library paths that you want R to search for packages.
                                 This could be useful if you are root users and install R
                                 packages in directories other than default.
                                 Multiple paths separated by ':' could be provided.

Data preparation:
-l/--min-depth                   Minimum depth to keep a SNP (DP).
-u/--max-depth                   Maximum depth to keep a SNP (DP).
-q/--min-qual                    Minimum quality to keep a SNP (QUAL).
-mf/--min-maf                    Minimum minor allele frequency to keep a SNP (default 0.1).
-mm/--max-missing                Maximum proportion of missing data to keep a SNP (default 0.5).

Haplotype inferring:
-x/--max-iter                    Maxmium rounds for EM optimization (default 100).
-f/--parent                      Parent samples (separated by a \":\").
-G/--genotype                    Use genotypes to infer haplotypes. Mutually exclusive with 
                                 option -D/--allele-depth and -L/--genetype likelihood.
-D/--allele-depth                Use allele depth to infer haplotypes. Mutually exclusive 
                                 with option -G/--genotype and -L/--genetype likelihood.
-L/--genotype-likelihood         Use genotype likelihoods to infer haplotypes. Mutually 
                                 exclusive with option -G/--genotype and -L/--allele-depth 
                                 (default).
-c/--min-snp-count               Minimum number of SNPs on a scaffold to run.
-r/--repeat                      Repeat haplotype inferring for multiple times as EM algorithm 
                                 could be trapped in local optima. The program takes three values 
                                 from here, i.e., for scaffold, for superscaffold and for genetic 
                                 linkage map refinement. Three values should be separated by ',' 
                                 (default 30,30,10).
-rr/--refinement-round           Number of rounds to refine pseudomelecules (default 10.)
   
Recombination frequency estimation:
-nb/--best                       The most likely nb haplotypes will be used (default 10).
-phi/--skew-phi                  For a haplotype inference, the frequencies of parental 
                                 haplotypes need to be in the interval [1/phi, phi], 
                                 otherwise will be discarded (default 2).
-nd/--drop                       At least nd haplotype inferences are required for calculation.		

Pseudomolecule construction:
-a/--input-assembly              Input assembly fasta file.
-frac/--frac-thresold            Lower threshold the genetic linkage map covers to construct
                                 pseudomolecules (default 0.8).
-gz/--genome-size                The estimated genome size (default size of the reference assembly).
</pre>

## Details about the output files
#### 1. Simulate an outcrossed F1 mapping population (popsimulation)
<pre>
a. ${out_dir}/Sc1/*.fasta.gz     Gzipped files with genome sequences. One for each sample.
b. ${out_dir}/*.*                Meta files used and generated by PedigreeSim V2.0.
</pre>

#### 2. Simualte GBS data from the F1 mapping population (gbssimulation)
<pre>
a. ${out_dir}/*.gz               Gzipped files with GBS reads. One for each sample.
b. ${out_dir}/*_key.txt          The key file for the GBS.
</pre>

#### 3. Run variant detection module for GBS data (gbspileup)
<pre>
a. ${out_dir}/tags               Encoded tag sequences called from the GBS FASTQ file(s).
b. ${out_dir}/mergedTags         Tag sequence files are merged into one.
c. ${out_dir}/tagFastq           The merged tag sequence file is decode to generate FASTQ file.
d. ${out_dir}/tagBam             BAM file generated from mapping tag FASTQ file to reference.
e. ${out_dir}/bam                Tag BAM file is distributed to generate BAM files for each sample.
f. ${out_dir}/bed                Reference is split for multithreading variant calling.
g. ${out_dir}/ref                Split reference files according to BED files.
h. ${out_dir}/splitBam           BAM file for each sample is split according to the BED files.
i. ${out_dir}/freebayes          Multithreading variant calling with freebayes.
    > out.vcf                    The resulted VCF file.
</pre>

#### 4. Run data preparation for haplotype phasing algorithm (datapreparation)
<pre>
a. ${out_dir}/*.recode.vcf       A preprocessed VCF files. Variants are filtered.
b. ${out_dir}/*.recode.zip       A zipped file can be recognised by the haplotype phasing algorithm.
</pre>

#### 5. Run haplotype phasing algorithm (haplotyper)
<pre>
a. ${out_dir}/*.zip              A zipped file contains all related results for haplotype phasing.
    > stderr_true                A log file.
    > snp_*.txt                  SNPs that used in this run.
    > results_hmm
       > emissionModel.txt       Emission probabilities of each allele for each parental haplotype
                                 at all loci.
       > transitionModel.txt     Transition probabilities or recombination frequencies between parental 
                                 haplotypes between adjacent loci.
    > phasedStates/*txt          Inferred inheritance pattern or haplotypes for all samples. 
</pre>

#### 6. Run genetic linkage map or pseudomolecule construnction pipeline (gembler)
<pre>
a. ${out_dir}/results            Result files including genetic linkage map (.mct), pseudomolecules (may 
                                 not generated) (.fa), and a log file for the genetic linkage map (.log).
b. ${out_dir}/data               Preprocessed data.
c. ${out_dir}/meta               Meta files generated by the program.
d. ${out_dir}/single_hap_infer   Haplotype phasing results for each single contig/scaffold.
e. ${out_dir}/2nn_hap_infer      Haplotype phasing results for superscaffolds.
f. ${out_dir}/refine_hap_infer   Haplotype phasing results for genetic linkage map refinement.
</pre>

## Citing polyGembler
