# Gut microbiome of C57BL/6 mice exposed to forest and urban soil
Sequences, processing files, otu tables and metadata for the 16S  SoilGut project. Data was processed using QIIME 1.9.0
##Sequencing data Raw MiSeq sequencing data can be found in the 0_raw_reads folder.

##Library split Libraries were split using split_libraries_fastq.py script in QIIME 1.9.0
```
split_libraries_fastq.py -i 4Rds_R1_001.fastq -b combined_barcodes_fixedheaders.fastq -m vmf.txt --barcode_type 16 --store_demultiplexed_fastq -q 19 -o slout_R1_q20

```
##Chimera filtering Chimeras were identified with the identify_chimeric_seqs.py script (using USEARCH 6.1 and Greengenes 13_8), and removed with the filter_fasta.py script in QIIME 1.9.0
```
identify_chimeric_seqs.py -i seqs.fna -m usearch61 -o usearch_checked_chimeras/ -r 97_otus.fasta

filter_fasta.py -f seqs.fna -o seqs_chimeras61_filtered.fna -s usearch61_checked_chimeras/chimeras.txt -n 

```
##OTU picking Open-reference OTU picking was performed at 97% identity using Greengenes version 13_8, with usearch61 using the pick_open_reference_otus.py script in QIIME 1.9.0
```
pick_open_reference_otus.py -i seqs_chimeras61_filtered.fna -m usearch61 -o otus_chimerafiltered_usearch61_openref

```
##Filter out low depth sequences, filtering out low depth samples, with counts/sample lower than 943 using the filter_samples_from_otu_table.py script in QIIME 1.9.0
```
filter_samples_from_otu_table.py -i otu_table_mc2_w_tax_no_pynast_failures.biom -o otu_943cutoff.biom -n 943

```
