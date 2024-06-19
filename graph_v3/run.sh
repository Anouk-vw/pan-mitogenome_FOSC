#!/bin/bash

set -e  # To prevent a repeat of the big scare of 04/05/2023

function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}
eval $(parse_yaml config.yaml "CONF_")

if [ ${CONF_multiple_chromosomes} == 1 ]
then
  echo "Running sequence partitioning"
if [ ! -d ${CONF_sample}seqpart/ ]
then
  python3 scripts/combine.py ${CONF_sample}
  mkdir ${CONF_sample}seqpart/
  mv ${CONF_sample}combined.fa ${CONF_sample}/seqpart/combined.fa
  wd=${CONF_sample}seqpart/
  bgzip -@ 4 ${wd}combined.fa
  #mash
  #mash dist ${wd}combined.fa.gz ${wd}combined.fa.gz -s 3000 -i > ${wd}distances.tsv
  #python3 scripts/mash2net.py -m ${wd}distances.tsv
  #python3 scripts/net2communities.py \
  #  -e ${wd}distances.tsv.edges.list.txt \
  #  -w ${wd}distances.tsv.edges.weights.txt \
  #  -n ${wd}distances.tsv.vertices.id2name.txt
  ##add wfmash 3-08-2023
  #gunzip ${wd}combined.fa.gz
  samtools faidx ${wd}combined.fa.gz
  wfmash ${wd}combined.fa.gz -p 90 -n 3 -t 10 -m >${wd}mapping.paf
  echo ${wd}mapping.paf
  python3 scripts/paf2net.py -p ${wd}mapping.paf
  python3 scripts/net2communities.py \
    -e ${wd}mapping.paf.edges.list.txt \
    -w ${wd}mapping.paf.edges.weights.txt \
    -n ${wd}mapping.paf.vertices.id2name.txt\
    --plot
  python3 scripts/net2communities.py \
    -e ${wd}mapping.paf.edges.list.txt \
    -w ${wd}mapping.paf.edges.weights.txt \
    -n ${wd}mapping.paf.vertices.id2name.txt
    
else
  #18/04/2024: modify to allow for predifined cliques
  wd='/home/anouk/anouk2/Black_Sigatoka/graph/cliquing/set_parameters/95-2500-50000-0.01/'
  #${CONF_sample}seqpart/

fi

  echo "Indexing data"
  #18/04/2024 modify to allow for predifined cliques (mapped -> out)
  ncommunities=$(ls ${wd} | grep out.paf.edges.weights.txt.community | wc -l)

  seq 0 $ncommunities | while read i; do
    echo "community $i"
    #18/04/2024 modify to allow for predifined cliques (mapped -> out)
    samtools faidx ${wd}combined.fa.gz $(cat ${wd}out.paf.edges.weights.txt.community.$i.txt) | \
    bgzip -@ 4 -c > ${wd}community.$i.fa.gz
  done

  echo "Sequence partitioning finished."
  for (( i = 1; i < $ncommunities; i++))
  do
    echo "Analysing community: ${i}"
    sed "s#sample:.*#sample: ${wd}community.${i}.fa.gz#g" config.yaml > temp.yaml && mv temp.yaml config.yaml
    sed "s#runid:.*#runid: ${CONF_runid}/community${i}#g" config.yaml > temp.yaml && mv temp.yaml config.yaml

    echo "Initialising..."
    python3 scripts/init.py

    echo "Running snakemake..."
    snakemake #-p --forcerun --cores
  done

else # no sequence partitioning
  echo "Initialising..."
  python3 scripts/init.py

  echo "Running snakemake..."
  snakemake -p --forcerun --cores
fi






echo "Done! Results can be found in ${PWD}/output/${CONF_runid}"

