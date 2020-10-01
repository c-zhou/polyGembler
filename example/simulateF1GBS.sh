#!/bin/sh

#### This is an example for simulating GBS data for an outbred mapping population.

## sbatch --nodes=1 --ntasks=16 --mem=32000 --time=1-00:00:00 --partition=snowy,physical --job-name=f12 --export=ploidy=2 simulateF1GBS.sh
## sbatch --nodes=1 --ntasks=16 --mem=32000 --time=1-00:00:00 --partition=snowy,physical --job-name=f14 --export=ploidy=4 simulateF1GBS.sh

module load java/1.8.0_241   # load java module. Edit this as needed before run this script

dataDir="data"
ref="${dataDir}/ref.fa.gz"
polyGembler="polyGembler-1.0-jar-with-dependencies.jar" ## Path to polyGembler executable jar file, please edit this before run this script

faOut="f1_"${ploidy}
gbsOut="gbs_"${ploidy}

if [ ! -d ${dataDir}/${faOut} ]; then mkdir ${dataDir}/${faOut}; fi
if [ ! -d ${dataDir}/${gbsOut} ]; then mkdir ${dataDir}/${gbsOut}; fi

java -jar ${polyGembler} popsimulation -r ${ref} -n 192 -p $ploidy -t 16 -o ${dataDir}/${faOut}
java -jar ${polyGembler} gbssimulation -f ${dataDir}/${faOut}/Sc1 -e ApeKI -t 16 -m 10 -s 5 -o ${dataDir}/${gbsOut}

f=`ls ${dataDir}/${gbsOut} | grep "key" | cut -d'_' -f1-2`

cat ${dataDir}/${gbsOut}/*.gz >"$f"_fastq.txt.gz

rm ${dataDir}/${gbsOut}/*.gz && mv "$f"_fastq.txt.gz ${dataDir}/${gbsOut}/

