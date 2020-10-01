#!/bin/bash

#### This is an example for running Tassel 5 GBS v2 Pipeline.
#### See details in https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline.

## sbatch --nodes=1 --ntasks=6 --mem=64000 --time=1-00:00:00 --partition=snowy,physical --job-name=gbs2 --export=ploidy=2 runGBSv2.sh
## sbatch --nodes=1 --ntasks=6 --mem=96000 --time=1-00:00:00 --partition=snowy,physical --job-name=gbs4 --export=ploidy=4 runGBSv2.sh

export PATH=~/sw/tassel-5-standalone/:$PATH              ## Path to TASSEL-GBS. Edit this before run this script
export JAVA_TOOL_OPTIONS=-Dfile.encoding=UTF8
## export JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=$(pwd)
module load bwa/0.7.17                                   ## make sure BWA is available in the evironment. Edit this as needed before run this script

if [ ! ${ploidy} ]; then echo "Error: need to specify the ploidy level" && exit 1; fi

if [ ${ploidy} -eq 2 ]
then
	mem="64G"
fi
if [ ${ploidy} -eq 4 ]
then
        mem="96G"
fi

gbsOut="gbs_"${ploidy}
tasselOut="tassel_"${ploidy}

threads="6"

bwa="bwa"
inDir="../data/${gbsOut}"
keyFile=`ls ../data/${gbsOut}/*_key.txt`
ref="../data/ctg.fa"
outDir="../data/${tasselOut}"

if [ ! -d ${inDir} ]; then echo "${inDir} not existed"; exit 1; fi
if [ ! -f ${keyFile} ]; then echo "${keyFile} not existed"; exit 1; fi
if [ ! -f ${ref} ]; then echo "${ref} not existed"; exit 1; fi
if [ ! -d ${outDir} ]; then mkdir ${outDir}; fi

run_pipeline.pl -Xmx${mem} -fork1 -GBSSeqToTagDBPlugin -e ApeKI -i ${inDir} -db ${outDir}/tmp.db -k ${keyFile} -kmerLength 90 -minKmerL 20 -c 3 -mxKmerNum 1000000000 -batchSize 4 -endPlugin -runfork1 >${outDir}/GBSSeqToTagDBPlugin.log 2>&1

run_pipeline.pl -Xmx${mem} -fork1 -TagExportToFastqPlugin -db ${outDir}/tmp.db -o ${outDir}/tags.fq.gz -c 1 -endPlugin -runfork1 >${outDir}/TagExportToFastqPlugin.log 2>&1

${bwa} index ${ref}

${bwa} mem -t ${threads} ${ref} ${outDir}/tags.fq.gz > ${outDir}/tags.sam 2>${outDir}/aln.log

run_pipeline.pl -Xmx${mem} -fork1 -SAMToGBSdbPlugin -i ${outDir}/tags.sam -db ${outDir}/tmp.db -mapper bwaMem -aProp 0.0 -aLen 0 -endPlugin -runfork1 >${outDir}/SAMToGBSdbPlugin.log 2>&1

run_pipeline.pl -Xmx${mem} -fork1 -DiscoverySNPCallerPluginV2 -db ${outDir}/tmp.db -ref ${ref} -maxTagsCutSite 1024 -endPlugin -runfork1 >${outDir}/DiscoverySNPCallerPluginV2.log 2>&1

run_pipeline.pl -Xmx${mem} -fork1 -ProductionSNPCallerPluginV2 -db ${outDir}/tmp.db -e ApeKI -i ${inDir} -k ${keyFile} -batchSize 4 -kmerLength 90 -o ${outDir}/out.vcf -endPlugin -runfork1 >${outDir}/ProductionSNPCallerPluginV2.log 2>&1

