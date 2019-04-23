#############################
#Blastocystis Reference Tree
#############################
cd /Users/mann/github/microeukITS
mkdir blastTree
#get references from ncbi using this query: txid376790[Organism:exp] (internal transcribed spacer 1) NOT (mitochondrial OR mitochondria OR complete genome OR whole genome shotgun) 
#sort and cluster
cd blastRefTree
grep ">" ref.fa | sed 's/ /_/g' | sed 's/_/\t/' | sed 's/>//' > ref.tax
sed -i 's/ .*//' ref.fa
sed -i '/^$/d' ref.fa
usearch -sortbylength ref.fa -fastaout ref.sort.fa -minseqlength 200
usearch -cluster_smallmem ref.sort.fa -id 0.99 -centroids ref.clustered.fa -uc ref.clusters.uc
mafft --reorder --auto ref.clustered.fa > ref.align.fa
trimal -in ref.align.fa -out ref.trimal.fa -gt 0.3 -st 0.001
#upload to cluster to build tree
scp ref.trimal.fa mann@zoology.ubc.ca:/parfreylab/mann/temp/
cd /parfreylab/mann/temp
raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -s ref.trimal.fa -n blast.tre
#now do epa placement of your reads with this reference tree
#first filter out your blastocystis reads
cd github/microeukITS/blastRefTree/
grep "Blastocystis" microeuk.txt | awk '{print $1}' > blastocystis/blasto.ids
cd blastocystis
filter_fasta.py -f /Volumes/histolytica/its/vsearch/repset.nonchim.fa -o blasto.fa -s blasto.ids
align_seqs.py -i blasto.fa -t blastRefTree/ref.align.fa -o htes_aligned -p 80
cat htes_aligned/blasto_aligned.fasta blastRefTree/ref.align.fa > htes_aligned/htes_aligned_curated.fa
filter_alignment.py -i htes_aligned/htes_aligned_curated.fa -g 0.99 -s -e 0.0001 -o htes_aligned/
scp htes_aligned/htes_aligned_curated_pfiltered.fasta mann@zoology.ubc.ca:/parfreylab/mann/temp
rm *epa.tre
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n epa.tre -s htes_aligned_curated_pfiltered.fasta -t RAxML_bestTree.blast.tre
sed 's/QUERY__//g' RAxML_labelledTree.epa.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placement_tree_epa.tre
#now on local machine again
mkdir placementTree
scp mann@zoology.ubc.ca:/parfreylab/mann/temp/RA*epa.tre placementTree/
#ok now that we have a good placement blastocystis tree make a maximum likelihood tree to get structure
#raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -s htes_aligned_curated.trimal.fa -n full.tre

#phyloseq analysis of blastocystis reads
library("phyloseq");packageVersion("phyloseq")
library(ggplot2)
library(ape)
otu <- read.table("microeuk_nohash.txt", sep="\t", header=T, row.names=1)
map <- as.data.frame(read.table("mapPhyseq.txt", sep="\t", header=T, row.names=1))
map <- sample_data(map)
tre <- read.tree("blastocystis/RAxML_bipartitions.final.tre")
otutable <- otu_table(otu, taxa_are_rows=T)
physeq <- phyloseq(otutable)
physeq <- merge_phyloseq(physeq, tre, map)

#cleaned tree in archaeoptyrx, removed messy outliers, long branches
#pull node names from tree, get subtype information to add to taxonomy

#now want to build a tree with all ncbi blastocystis its regions, see if we can't get some more info on the unknown subtypes
usearch -sortbylength all_ncbi_blast_its.fasta -fastaout all_ncbi_blast_its.sort.fa -minseqlength 200
usearch -cluster_smallmem all_ncbi_blast_its.sort.fa -id 0.99 -centroids all_ncbi_blast_its.clustered.fa -uc all_ncbi_blast_its.clusters.uc
mafft --reorder --auto all_ncbi_blast_its.clustered.fa > all_ncbi_blast_its.align.fa
trimal -in all_ncbi_blast_its.align.fa -out all_ncbi_blast_its.trimal.fa -gt 0.90 -st 0.001
scp all_ncbi_blast_its.trimal.fa mann@zoology.ubc.ca:/parfreylab/mann/temp
raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -s all_ncbi_blast_its.trimal.fa -n all_ncbi.tre
#do placement tree with this as reference
align_seqs.py -i blasto.fa -t all_ncbi_blast_its.align.fa -o htes_aligned_allNCBI -p 80
#do we get more sequence aligned here? We do, still not sure if can type though
cat htes_aligned_allNCBI/blasto_aligned.fasta all_ncbi_blast_its.align.fa > htes_aligned_allNCBI/htes_algined_curated.fa
filter_alignment.py -i htes_aligned_allNCBI/htes_algined_curated.fa -g 0.99 -s -e 0.001 -o htes_aligned_allNCBI/
scp htes_aligned_allNCBI/htes_algined_curated_pfiltered.fasta mann@zoology.ubc.ca:/parfreylab/mann/temp/
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n allNCBI.epa.tre -s htes_algined_curated_pfiltered.fasta -t RAxML_bestTree.all_ncbi.tre
sed 's/QUERY__//g' RAxML_labelledTree.allNCBI.epa.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placement_tree_allNCBI.tre
scp mann@zoology.ubc.ca:/parfreylab/mann/temp/RAxML_placement_tree_allNCBI.tre .