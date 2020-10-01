#!/bin/bash

#### make sure java and R/Rscript is available
module load java/1.8.0_241 r/3.6.2 # load java and R module. Edit this as needed before run this script

method=$(echo ${method} | tr '[:upper:]' '[:lower:]')

case ${method} in
"haplotyper")
    #### run haplotype phasing for each scaffold listed in the $ctgList
    ## sbatch --nodes=1 --ntasks=1 --mem=3000 --time=12:00:00 --partition=physical,snowy --job-name=haplo --export=polyGembler=${polyGembler},ctgList=${},datZip=${},id=${},ploidy=${},fhaps=${},nr=${},field=${},outs=${},method="haplotyper" --array=1-${} gembler_wrapper.sh
    #### run haplotype phasing for each superscaffold listed in the $ssList
    ## sbatch --nodes=1 --ntasks=1 --mem=3000 --time=12:00:00 --partition=physical,snowy --job-name=haplo --export=polyGembler=${polyGembler},ssList=${},datZip=${},id=${},ploidy=${},fhaps=${},nr=${},field=${},outs=${},method="haplotyper" --array=1-${} gembler_wrapper.sh
    #### run haplotype phasing according to config file. Config file contains three columns: scf\tnr\touts
    ## sbatch --nodes=1 --ntasks=1 --mem=3000 --time=12:00:00 --partition=physical,snowy --job-name=haplo --export=polyGembler=${polyGembler},config=${},datZip=${},id=${},ploidy=${},fhaps=${},nr=${},field=${},outs=${},method="haplotyper" --array=1-${} gembler_wrapper.sh
    if [ ${ctgList} ]
    then
        scf=$(head -$SLURM_ARRAY_TASK_ID ${ctgList} | tail -1 | awk '{print $1}')
        scfOpts="-c ${scf}"
    elif [ ${ssList} ]
    then
        scfOpts=$(head -$SLURM_ARRAY_TASK_ID ${ssList} | tail -1)
    elif [ ${config} ]
    then
        scfOpts=$(head -$SLURM_ARRAY_TASK_ID ${config} | tail -1 | cut -f1)
        nr=$(head -$SLURM_ARRAY_TASK_ID ${config} | tail -1 | cut -f2)
        outs=$(head -$SLURM_ARRAY_TASK_ID ${config} | tail -1 | cut -f3)
    else
        exit 1
    fi
    for r in $(seq 1 ${nr})
    do
        java -Xmx3G -jar ${polyGembler} haplotyper -i ${datZip} --experiment-id ${id} --max-iter 1000 --ploidy ${ploidy} --parent ${fhaps} ${scfOpts} ${field} -o ${outs}
    done
    ;;
"asmerr")
    #### detect assembly errors
    ## sbatch --nodes=1 --ntasks=8 --mem=16000 --time=12:00:00 --partition=snowy,physical --job-name=asmerr --export=polyGembler=${polyGembler},inHaps=${},id=${},inVCF=${},asmErrRF=${},threads=${},outPref=${},method="asmerr" gembler_wrapper.sh
    java -jar ${polyGembler} asmerr -i ${inHaps} --experiment-id ${id} -vi ${inVCF} -r ${asmErrRF} -t ${threads} -o ${outPref}
    ;;
"singlepoint")
    #### estimate RFs between neighboring markers on the scaffold
    ## sbatch --nodes=1 --ntasks=8 --mem=16000 --time=12:00:00 --partition=snowy,physical --job-name=sprf --export=polyGembler=${polyGembler},inHaps=${},id=${},threads=${},outPref=${},method="singlepoint" gembler_wrapper.sh
    java -jar ${polyGembler} singlepoint -i ${inHaps} --experiment-id ${id} -t ${threads} -o ${outPref}
    ;;
"twopoint")
    #### estimate RFs between each pair of scaffolds
    ## sbatch --nodes=1 --ntasks=8 --mem=16000 --time=12:00:00 --partition=snowy,physical --job-name=tprf --export=polyGembler=${polyGembler},inHaps=${},id=${},threads=${},outPref=${},method="twopoint" gembler_wrapper.sh
    java -jar ${polyGembler} twopoint -i ${inHaps} --experiment-id ${id} -t ${threads} -o ${outPref}
    ;;
"map")
    #### linkage mapping: grouping and ordering
    ## sbatch --nodes=1 --ntasks=8 --mem=16000 --time=12:00:00 --partition=snowy,physical --job-name=lmap --export=polyGembler=${polyGembler},tpRfFile=${},spRfFile=${},lodScore=${},threads=${},outPref=${},method="map" gembler_wrapper.sh
    #### linkage mapping: no grouping, ordering only
    ## sbatch --nodes=1 --ntasks=8 --mem=16000 --time=12:00:00 --partition=snowy,physical --job-name=lmap --export=polyGembler=${polyGembler},tpRfFile=${},spRfFile=${},lodScore=${},threads=${},outPref=${},opts="-1",method="map" gembler_wrapper.sh
    java -jar ${polyGembler} map -i ${tpRfFile} -m ${spRfFile} -l ${lodScore} -t ${threads} -o ${outPref} ${opts}
    ;;
"superscaffold")
    #### build super scaffolds using nearest neighbour joining
    ## sbatch --nodes=1 --ntasks=1 --mem=8000 --time=12:00:00 --partition=snowy,physical --job-name=ss --export=polyGembler=${polyGembler},tpRfFile=${},outPref=${},method="superscaffold" gembler_wrapper.sh
    java -jar ${polyGembler} superscaffold -i ${tpRfFile} -o ${outPref}
    ;;
"chromosomer")
    #### construct pseudomolecules from the genetic linkage maps
    ## sbatch --nodes=1 --ntasks=1 --mem=8000 --time=12:00:00 --partition=snowy,physical --job-name=chroms --export=polyGembler=${polyGembler},mctFile=${},inFasta=${},errFile=${},outPref=${},method="chromosomer" gembler_wrapper.sh
    java -jar ${polyGembler} chromosomer -i ${mctFile} -a ${inFasta} -e ${errFile} -o ${outPref}
    ;;
*)
    echo "Undefined method: ${method}!!!"
    ;;
esac


