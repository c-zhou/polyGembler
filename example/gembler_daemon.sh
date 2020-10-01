#!/bin/bash

#### Need a HPC with slurm to run this script
#### You can either run this script locally or submit it as slurm batch script
#### For example:
#### sbatch --nodes=1 --ntasks=1 --mem=8000 --time=3-00:00:00 --job-name=gembler gembler_daemon.sh
#### Required R packages: argparse, TSP, MDSMap, igraph, doParallel, foreach
#### Optional R packages: unixtools
#### Required R packages if running SNP filtering: updog, XNomial
#### Make sure java and R/Rscript is available in the environment

#### An example for simulating GBS data.
## jobid1=$(sbatch --nodes=1 --ntasks=16 --mem=32000 --time=1-00:00:00 --partition=snowy,physical --job-name=f12 --export=ploidy=2 simulateF1GBS.sh | tail -1 | awk '{print $NF}')
#### An example for running Tassel 5 GBS v2 Pipeline.
#### See details in https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline.
## jobid2=$(sbatch --nodes=1 --ntasks=6 --mem=64000 --time=1-00:00:00 --partition=snowy,physical --job-name=gbs2 --export=ploidy=2 --dependency=afterok:${jobid1} runGBSv2.sh | tail -1 | awk '{print $NF}')
#### For SNP filtering: 1. call genotypes from allele depth with the R package updog using f1 model; 2. perform a multinomial exact test to test the segregation ratio.
#### This step is not mandatory. It might filter out a lot of SNPs especially when data coverage is low.
## jobid3=$(sbatch --nodes=1 --ntasks=16 --mem=32000 --time=1-00:00:00 --partition=snowy,physical --job-name=filter2 --export=ploidy=2 --dependency=afterok:${jobid2} SNPfiltering.sh | tail -1 | awk '{print $NF}')


######################################################
########## EDIT THESE PARAMETERS BEFORE RUN ##########
######################################################
module load java/1.8.0_241 r/3.6.2                                # load java and R module. Edit this as needed before run this script

inVCF="data/out2_filtered.vcf.gz"                                 # Input VCF file
inFASTA="data/ctg.fa.gz"                                          # Input FASTA file for contigs and scaffolds
ploidy=2                                                          # Ploidy
fhaps="P1:P2"                                                     # Two parents for the F1 mapping population
minSNP=5                                                          # Minimum number of SNPs a contig/scaffold harbors to be included in linkage analysis
field="-D"                                                        # Field used for haplotype phasing: '-D' for allele depth and '-G' for genotype
nr=10                                                             # Number of repeat runs for haplotype phasing
lodScore=3                                                        # LOD score threshold for linkage analysis
asmErrRF="0.1"                                                    # Will be detected as an assembly error if the RF estimation between two neighboring SNPs was larger than this value
experimentId="zzz"                                                # An arbitrary idenfier for running haplotype phasing algorithm
polyGembler="./polyGembler-1.0-jar-with-dependencies.jar"         # Absolute path to the polyGembler executable jar file
slurmOpts="--partition=snowy,physical --nodes=1 --time=12:00:00"  # Universal options for slurm jobs, NEED to change it to fit your system
nnSS=1                                                            # Refine RF esmitimation by construction of super scaffolds, set to 0 to skip
polishRound=3                                                     # Round of polishing for each linkage group, set to zero to skip polishing

outDir="results_$(date '+%Y-%m-%d-%H:%M:%S')"                     # Output dir for this run
mkdir ${outDir}
logDir="${outDir}/logs"                                           # Contains all log files of slurm jobs
mkdir ${logDir}

######################################################
########## NO NEED TO CHANGE THE CODE BELOW ##########
######################################################

log() {
    mesg=$1
    echo "[$(date)] ${mesg}"
}

#### prepare data for running polyGembler
log "Prepare data for running polyGembler"
datPref="out1"
datZip="${outDir}/${datPref}.zip"
ctgList="${outDir}/${datPref}_ctg.txt"
java -jar ${polyGembler} datapreparation -i ${inVCF} -s ${datPref} -o ${outDir}

#### get the list of contigs/scaffolds that will be included in linkage analysis
unzip -p ${datZip} contig >${ctgList}
nScf2Run=$(cat ${ctgList} | awk -v m=${minSNP} '{if($2>=m) n++}END{print n}')

#### initial run of the polyGembler haplotype phasing
log "Initial run of the polyGembler haplotype phasing"
outs="${outDir}/h1"
mkdir ${outs}
msg=$(sbatch --ntasks=1 --mem=3000 --wait $slurmOpts --job-name=haplo --export=polyGembler=${polyGembler},ctgList=${ctgList},datZip=${datZip},id=${experimentId},ploidy=${ploidy},fhaps=${fhaps},nr=${nr},field=${field},outs=${outs},method="haplotyper" --array=1-${nScf2Run} gembler_wrapper.sh)
exitCode1=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
if [ ${exitCode1} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode1=${exitCode1}). Check log."; exit 1; fi
jobId=$(echo ${msg} | awk '{print $NF}')
mv slurm-${jobId}*.out ${logDir}/
log "Slurm job 1 check-in: ${jobId}"

#### detect assembly errors
log "Detect assembly errors"
threads=8
datPref="out2"
outPref="${outDir}/${datPref}"
msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=asmerr --export=polyGembler=${polyGembler},inHaps=${outs},id=${experimentId},inVCF=${inVCF},asmErrRF=${asmErrRF},threads=${threads},outPref=${outPref},method="asmerr" gembler_wrapper.sh)
exitCode2=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
if [ ${exitCode2} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode2=${exitCode2}). Check log."; exit 1; fi
jobId=$(echo ${msg} | awk '{print $NF}')
mv slurm-${jobId}*.out ${logDir}/
log "Slurm job 2 check-in: ${jobId}"

#### process assembly errors if any
nAsmErr=$(cat ${outPref}".err" | wc -l)
if [ ${nAsmErr} -gt 0 ]
then
    ## merge VCF files for old and new contigs/scaffolds
    mv ${outPref}".vcf" ${outPref}"_tmp.vcf"
    awk '!/^#/{s[$1]++}END{for(i in s) print i"\t"s[i]}' ${outPref}"_tmp.vcf"  | sort -k2r,2n -k1,1 >${outPref}"_ctg.txt"
    if [[ ${inVCF} =~ .gz$ ]]
    then
        zcat ${inVCF} | grep -v "^#" >> ${outPref}"_tmp.vcf"
    else
        cat ${inVCF} | grep -v "^#" >> ${outPref}"_tmp.vcf"
    fi
    
    ## construct data file used by polyGembler and update related file names
    java -jar ${polyGembler} datapreparation -i ${outPref}"_tmp.vcf" -s ${datPref} -o ${outDir}
    rm ${outPref}"_tmp.vcf"
    datZip=${outPref}".zip"
    ctgList=${outPref}"_ctg.txt"
    nScf2Run=$(cat ${ctgList} | awk -v m=${minSNP} '{if($2>=m) n++}END{print n}')

    ## move haplotype phasing results for the misassembled contigs/scaffolds out of dir
    outsErr="${outDir}/herr"
    mkdir ${outsErr}
    for f in $(cat ${outPref}".err" | awk '{print $1}' | sort -u)
    do
        mv ${outs}/zzz.${f}.* ${outsErr}/
    done

    ## run haplotype phasing for new contigs/scaffolds
    log "Run haplotype phasing for new contigs/scaffolds"
    msg=$(sbatch --ntasks=1 --mem=3000 --wait $slurmOpts --job-name=haplo --export=polyGembler=${polyGembler},ctgList=${ctgList},datZip=${datZip},id=${experimentId},ploidy=${ploidy},fhaps=${fhaps},nr=${nr},field=${field},outs=${outs},method="haplotyper" --array=1-${nScf2Run} gembler_wrapper.sh)
    exitCode3=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
    if [ ${exitCode3} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode3=${exitCode3}). Check log."; exit 1; fi
    jobId=$(echo ${msg} | awk '{print $NF}')
    mv slurm-${jobId}*.out ${logDir}/
    log "Slurm job 3 check-in: ${jobId}"
fi

if [ ${nnSS} -gt 0 ]
then
    #### construct superscaffold using nearest-neighbor joining
    log "Construct superscaffold using nearest-neighbor joining"
    threads=8
    outPref="${outDir}/out3"
    
    ## two-point RF estimation
    msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=tprf --export=polyGembler=${polyGembler},inHaps=${outs},id=${experimentId},threads=${threads},outPref=${outPref},method="twopoint" gembler_wrapper.sh)
    exitCode4=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
    if [ ${exitCode4} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode4=${exitCode4}). Check log."; exit 1; fi
    jobId=$(echo ${msg} | awk '{print $NF}')
    mv slurm-${jobId}*.out ${logDir}/
    log "Slurm job 4 check-in: ${jobId}"

    ## superscaffold construciton
    msg=$(sbatch --ntasks=1 --mem=8000 --wait $slurmOpts --job-name=ss --export=polyGembler=${polyGembler},tpRfFile="${outPref}.txt",outPref=${outPref},method="superscaffold" gembler_wrapper.sh)
    exitCode5=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
    if [ ${exitCode5} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode5=${exitCode5}). Check log."; exit 1; fi
    jobId=$(echo ${msg} | awk '{print $NF}')
    mv slurm-${jobId}*.out ${logDir}/
    log "Slurm job 5 check-in: ${jobId}"

    #### run haplotyp phasing for superscaffolds
    log "Run haplotyp phasing for superscaffolds"
    outs="${outDir}/h2"
    mkdir ${outs}
    ssList="${outPref}.nns"
    nScf2Run=$(cat ${ssList} | wc -l)
    msg=$(sbatch --ntasks=1 --mem=3000 --wait $slurmOpts --job-name=haplo --export=polyGembler=${polyGembler},ssList=${ssList},datZip=${datZip},id=${experimentId},ploidy=${ploidy},fhaps=${fhaps},nr=${nr},field=${field},outs=${outs},method="haplotyper" --array=1-${nScf2Run} gembler_wrapper.sh)
    exitCode6=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
    if [ ${exitCode6} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode6=${exitCode6}). Check log."; exit 1; fi
    jobId=$(echo ${msg} | awk '{print $NF}')
    mv slurm-${jobId}*.out ${logDir}/
    log "Slurm job 6 check-in: ${jobId}"
fi

#### construct genetic linakge maps
log "Construct genetic linakge maps"
threads=8
outPref="${outDir}/out4"
## two-point RF estimation
msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=tprf --export=polyGembler=${polyGembler},inHaps=${outs},id=${experimentId},threads=${threads},outPref=${outPref},method="twopoint" gembler_wrapper.sh)
exitCode7=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
if [ ${exitCode7} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode7=${exitCode7}). Check log."; exit 1; fi
jobId=$(echo ${msg} | awk '{print $NF}')
mv slurm-${jobId}*.out ${logDir}/
log "Slurm job 7 check-in: ${jobId}"

## single-point RF estimation
msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=sprf --export=polyGembler=${polyGembler},inHaps=${outs},id=${experimentId},threads=${threads},outPref=${outPref},method="singlepoint" gembler_wrapper.sh)
exitCode8=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
if [ ${exitCode8} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode8=${exitCode8}). Check log."; exit 1; fi
jobId=$(echo ${msg} | awk '{print $NF}')
mv slurm-${jobId}*.out ${logDir}/
log "Slurm job 8 check-in: ${jobId}"

## linkage mapping
msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=lmap --export=polyGembler=${polyGembler},tpRfFile="${outPref}.txt",spRfFile="${outPref}.map",lodScore=${lodScore},threads=${threads},outPref=${outPref},method="map" gembler_wrapper.sh)
exitCode9=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
if [ ${exitCode9} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode9=${exitCode9}). Check log."; exit 1; fi
jobId=$(echo ${msg} | awk '{print $NF}')
mv slurm-${jobId}*.out ${logDir}/
log "Slurm job 9 check-in: ${jobId}"

#### prepare files for polishing
if [ ${polishRound} -gt 0 ]
then
    log "Prepare files for polishing"
    lgOutDir="${outDir}/hlg"
    mkdir ${lgOutDir}
    nLg=$(cat "${outPref}.par" | wc -l)
fi

for k in $(seq 1 ${polishRound})
do
    log "Polishing round ${k}"
    ## make config file
    conf="${outDir}/tmp.conf"
    cat "${outDir}/out$((3+k)).par" | awk -v k=${k} -v n=${nr} -v o=${lgOutDir} '{for(i=0;i<n;i++) {print $0"\t"1"\t"o"/"k"_"NR}}' >${conf}
    for f in $(cut -f3 ${conf} | sort -u); do mkdir $f; done
    
    ## run haplotype phasing with the config file
    nScf2Run=$(cat ${conf} | wc -l)
    msg=$(sbatch --ntasks=1 --mem=3000 --wait $slurmOpts --job-name=haplo --export=polyGembler=${polyGembler},config=${conf},datZip=${datZip},id=${experimentId},ploidy=${ploidy},fhaps=${fhaps},nr=${nr},field=${field},outs=${outs},method="haplotyper" --array=1-${nScf2Run} gembler_wrapper.sh)
    exitCode10=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
    if [ ${exitCode10} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode10=${exitCode10}). Check log."; exit 1; fi
    jobId=$(echo ${msg} | awk '{print $NF}')
    mv slurm-${jobId}*.out ${logDir}/
    log "Slurm job 10 check-in: ${jobId}"

    for j in $(seq 1 ${nLg})
    do
        threads=8
        outPref="${outDir}/out$((4+k))_${j}"
        inHaps="${lgOutDir}/${k}_${j}"
        ## two-point RF estimation
        msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=tprf --export=polyGembler=${polyGembler},inHaps=${inHaps},id=${experimentId},threads=${threads},outPref=${outPref},method="twopoint" gembler_wrapper.sh)
        exitCode11=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
        if [ ${exitCode11} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode11=${exitCode11}). Check log."; exit 1; fi
        jobId=$(echo ${msg} | awk '{print $NF}')
        mv slurm-${jobId}*.out ${logDir}/
        log "Slurm job 11 check-in: ${jobId}"

        ## single-point RF estimation
        msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=sprf --export=polyGembler=${polyGembler},inHaps=${inHaps},id=${experimentId},threads=${threads},outPref=${outPref},method="singlepoint" gembler_wrapper.sh)
        exitCode12=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
        if [ ${exitCode12} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode12=${exitCode12}). Check log."; exit 1; fi
        jobId=$(echo ${msg} | awk '{print $NF}')
        mv slurm-${jobId}*.out ${logDir}/
        log "Slurm job 12 check-in: ${jobId}"

        ## linkage mapping
        msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=lmap --export=polyGembler=${polyGembler},tpRfFile="${outPref}.txt",spRfFile="${outPref}.map",lodScore=${lodScore},threads=${threads},outPref=${outPref},opts="-1",method="map" gembler_wrapper.sh)
        exitCode13=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
        if [ ${exitCode13} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode13=${exitCode13}). Check log."; exit 1; fi
        jobId=$(echo ${msg} | awk '{print $NF}')
        mv slurm-${jobId}*.out ${logDir}/
        log "Slurm job 13 check-in: ${jobId}"
    done

    ## merge results for linkage groups
    outPref="${outDir}/out$((4+k))"
    outHaps="${lgOutDir}/round_${k}"
    mkdir ${outHaps}
    for j in $(seq 1 ${nLg})
    do
        cat "${outPref}_${j}.par" >>"${outPref}.par"
        head -1 "${outPref}_${j}.mct" | sed "s|LG1|LG${j}|g" >>"${outPref}.mct"
        tail -n +2 "${outPref}_${j}.mct" >>"${outPref}.mct"
        mv ${lgOutDir}/${k}_${j}/* ${outHaps}/
        ## clean data
        rm -r ${outPref}_${j}.* ${lgOutDir}/${k}_${j}
    done
    
    ## estimate RF all together
    threads=8
    ## two-point RF estimation
    msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=tprf --export=polyGembler=${polyGembler},inHaps=${outHaps},id=${experimentId},threads=${threads},outPref=${outPref},method="twopoint" gembler_wrapper.sh)
    exitCode14=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
    if [ ${exitCode14} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode14=${exitCode14}). Check log."; exit 1; fi
    jobId=$(echo ${msg} | awk '{print $NF}')
    mv slurm-${jobId}*.out ${logDir}/
    log "Slurm job 14 check-in: ${jobId}"

    ## single-point RF estimation
    msg=$(sbatch --ntasks=${threads} --mem=16000 --wait $slurmOpts --job-name=sprf --export=polyGembler=${polyGembler},inHaps=${outHaps},id=${experimentId},threads=${threads},outPref=${outPref},method="singlepoint" gembler_wrapper.sh)
    exitCode15=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
    if [ ${exitCode15} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode15=${exitCode15}). Check log."; exit 1; fi
    jobId=$(echo ${msg} | awk '{print $NF}')
    mv slurm-${jobId}*.out ${logDir}/
    log "Slurm job 15 check-in: ${jobId}"
done

#### prepare final results
log "Prepare final results"
wd=$(pwd)
cd ${outDir}
hi=$(ls | grep -P "^out.*\.mct$" | cut -d'.' -f1 | sed 's/out//g' | sort -r -n | head -1)
ln -s out2.err final.err
ln -s out${hi}.txt final.txt
ln -s out${hi}.map final.map
ln -s out${hi}.par final.par
ln -s out${hi}.mct final.mct
cd ${wd}

#### construct pseudomolecules
log "Construct pseudomolecules"
outPref="${outDir}/final"
msg=$(sbatch --ntasks=1 --mem=8000 --wait $slurmOpts --job-name=chroms --export=polyGembler=${polyGembler},mctFile="${outDir}/final.mct",inFasta=${inFASTA},errFile="${outDir}/final.err",outPref=${outPref},method="chromosomer" gembler_wrapper.sh)
exitCode16=$(sacct -j ${msg##* } -o ExitCode | tail -n +3 | awk -v FS=":" 'BEGIN{code=0}{if($1>code) code=$1}END{print code}')
if [ ${exitCode16} -gt 0 ]; then echo "One or more slurm jobs with non zero exit status (exitCode16=${exitCode16}). Check log."; exit 1; fi
jobId=$(echo ${msg} | awk '{print $NF}')
mv slurm-${jobId}*.out ${logDir}/
log "Slurm job 16 check-in: ${jobId}"

