#!/bin/bash
# Quality filter, merge, and process raw data
# Use: ./preprocessing.sh
# Downloaded demultiplexed data from deepthought, moved to histolytica in folder called raw

VSEARCH=/usr/local/bin/vsearch

#########################
#Trim primers & adapters
#########################
cd /Volumes/histolytica/its/raw
mkdir ../trim
echo "Trimming adapter and primer sequences"
ls *R1* | sed 's/_L001.*//' | parallel -j4 'cutadapt -a AATGATACGGCGACCACCGAGATCTAC -A CAAGCAGAAGACGGCATACGAGAT --trim-n -o ../trim/{}.R1.trim.fastq -p ../trim/{}.R2.trim.fastq {}_L001_R1_001.fastq.gz {}_L001_R2_001.fastq.gz 1>../trim/{}.adapt.trim.info'
cd ../trim
ls *R1* | sed 's/.R1.trim.fastq//' | parallel -j4 'cutadapt -a TCCGTAGGTGAACCTGCGG  -A CGGCTGCGTTCTTCATCGATGC  -o {}.R1.trimprimer.fastq -p {}.R2.trimprimer.fastq  {}.R1.trim.fastq {}.R2.trim.fastq 1> {}.primer.trim.info'
rm *trim.fastq

#########################
#Read merging
#########################
echo "Paired end read merging"
mkdir ../pear
#merge paired end reads
ls *R1* | sed 's/.R1.trimprimer.fastq//' | parallel -j4 'pear -f {}.R1.trimprimer.fastq -r {}.R2.trimprimer.fastq -o ../pear/{} -q 30 -n 50 -t 50 1>../pear/{}.out'
rm *trimprimer.fastq
cd ../pear
#get stats on merge rates
ls *out | sed 's/_.*//' > list
ls *out | while read line; do grep "Assembled reads ...................:" $line | grep -v "file" | awk -F":" '{print $2}' | awk -F" " '{print $4}' | sed 's/(//' | sed 's/)//' ; done > pmerge
paste list pmerge > mergerate.txt
rm list pmerge
#get rid of unmerged reverse and discarded read files
rm *discarded.fastq *unassembled.reverse.fastq
#remove empty files
find . -empty -delete

#########################
#Convert to fastq
#########################
echo "Fastq to fasta"
#convert to fasta and add sequence headers to both the merged and unmerged forward reads
ls *assembled.fastq | sed 's/.fastq//' | parallel -j4 'fastq_to_fasta -i {}.fastq -o {}.fasta -Q 33'
ls *unassembled.forward.fastq | sed 's/.fastq//' | parallel -j4 'fastq_to_fasta -i {}.fastq -o {}.fasta -Q 33'
#remove fastq
rm *fastq
#to keep fastq and zip up uncomment below
#ls *fastq | parallel 'gzip {}' &

echo "Adding headers and cleaning up"
#add headers
ls *assembled.fasta | sed 's/.assembled.fasta//' | while read line; do awk '/^>/{print ">'$line'."; next} {print}' < $line.assembled.fasta > $line.assembled.headers.fa; done
rm *assembled.fasta
ls *unassembled.forward.fasta | sed 's/.unassembled.forward.fasta//' | while read line; do awk '/^>/{print ">'$line'."; next} {print}' < $line.unassembled.forward.fasta > $line.unassembled.forward.headers.fa; done
rm *unassembled.forward.fasta
ls *assembled.headers.fa | sed 's/.assembled.headers.fa//' | while read line; do cat $line.assembled.headers.fa $line.unassembled.forward.headers.fa > $line.full.fa; done
ls *full.fa | sed 's/.full.fa//' | while read line; do awk '/^>/{$0=$0""(++i)}1' $line.full.fa > $line.full.fix.fa; done
cat *full.fix.fa > combinedseq.fa
#clean up files
rm *full.fa
rename 's/full.fix.fa/full.fa/' *full.fix.fa
rm *headers* *full*
cd ..
mkdir vsearch

#########################
#Denovo ref database 
#########################
#get representative sequences from dereplicated and sorted sequences
echo "Generating representative sequences file"
vsearch --derep_full pear/combinedseq.fa --output pear/combinedseq.derep.fa --sizeout --minuniquesize 5
vsearch --sortbylength pear/combinedseq.derep.fa --sizeout --output pear/combinedseq.sort.fa
vsearch --cluster_unoise pear/combinedseq.sort.fa --relabel OTU_ --centroids vsearch/repset.fa  
rm pear/combinedseq.sort.fa pear/combinedseq.derep.fa

#########################
#Chimera check
#########################
echo "Chimera check"
vsearch --uchime_denovo vsearch/repset.fa --chimeras vsearch/repset.chims.fa --nonchimeras vsearch/repset.nonchim.fa
rm vsearch/repset.fa

#########################
#OTU picking
#########################
echo "Picking OTUs"
vsearch --usearch_global pear/combinedseq.fa --db vsearch/repset.nonchim.fa --strand plus --uc vsearch/map.uc --maxaccepts 10 --id 0.99
rm pear/combinedseq.fa.gz
gzip pear/combinedseq.fa 
biom from-uc -i vsearch/map.uc -o otus.biom
biom summarize-table -i otus.biom -o otus.sum

#########################
#Taxonomy assignment
#########################
echo "Assigning taxonomy"
mkdir tax
echo "Get taxonomy assignments for repset using BLAST"
##NOTE: running these on the lab server, cruncher 1 
#lower pid threshold to get more results, should be resolved with megan
#below line for parallel, can't set outfmt columns
#cat repset.nonchim.fa | parallel --block 50k --recstart '>' --pipe blastn -qcov_hsp_perc 0.90 -evalue 1e-10 -max_target_seqs 500 -perc_identity 0.90 -db /parfreylab/shared/databases/NCBI_NT/NT_06.08.2018/nt -outfmt 6 -query - > tax/blast.out
blastn -evalue 1e-10 -max_target_seqs 500 -perc_identity 0.90 -db /parfreylab/shared/databases/NCBI_NT/NT_06.08.2018/nt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -query repset.nonchim.fa -out blast.out
#remove hits that are less than 95% similar to the target, query length must be more than 90% when aligned to subject sequence
awk '$4/$13 >= 0.95 && $3 >= 0.99 {print $0}' blast.out > goodblast.out
#open in megan to get taxonomy
#clean up taxonomy file
sed -i 's/"root;cellular organisms;//' tax/taxonomy.txt
sed -i 's/"//' tax/taxonomy.txt
sed -i 's/ /_/g' tax/taxonomy.txt
#add taxonomy to otu table
biom add-metadata -i otus.biom -o otus.wtax.biom --observation-metadata-fp tax/taxonomy.txt --observation-header OTUID,taxonomy --sc-separated taxonomy
#rarify
single_rarefaction.py -i otus.wtax.biom -o otus.d1k.biom -d 1000
#filter out those with no taxonomy and non microeuks
filter_otus_from_otu_table.py -i otus.d1k.biom -o microeuk.biom --negate_ids_to_exclude -e tax/microeuk_tax.txt -n 10
biom summarize-table -i microeuk.biom -o microeuk.sum
biom convert -i microeuk.biom -o microeuk.txt --table-type="OTU table" --to-tsv --header-key taxonomy

# Using ITSoneDB
# alias uclust=$VSEARCH
# assign_taxonomy.py -i vsearch/repset.nonchim.fa -t ~/refDB/ITSoneDB_rep_seq_1.131.tax -r ~/refDB/ITSoneDB_rep_seq_1.131.fasta -m uclust -o tax/






