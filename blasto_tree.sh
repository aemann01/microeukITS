#############################
#Blastocystis Reference Tree
#############################
#first get outgroups from pr2 database
cd /Users/mann/github/microeukITS
mkdir blastTree
#get references from ncbi using this query: txid376790[Organism:exp] (internal transcribed spacer 1) NOT (mitochondrial OR mitochondria OR complete genome OR whole genome shotgun) 
#sort and cluster
cd blastTree
grep ">" ref.fa | sed 's/ /_/g' | sed 's/_/\t/' | sed 's/>//' > ref.tax
sed -i 's/ .*//' ref.fa
sed -i '/^$/d' ref.fa
usearch -sortbylength ref.fa -fastaout ref.sort.fa -minseqlength 500
usearch -cluster_smallmem ref.sort.fa -id 0.99 -centroids ref.clustered.fa -uc ref.clusters.uc
mafft --reorder --auto ref.clustered.fa > ref.align.fa
trimal -in ref.align.fa -out ref.trimal.fa -gt 0.3 -st 0.001
#upload to cluster to build tree
scp ref.trimal.fa mann@zoology.ubc.ca:/parfreylab/mann/temp/
cd /parfreylab/mann/temp
raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -s ref.trimal.fa -n blast.tre
