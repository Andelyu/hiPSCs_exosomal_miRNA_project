#!/bin/bash

#-------------prepare-----------------------#

ref_pre_miRNA='/data1/miRNA/hsa_hairpin.fa'
ref_mature_miRNA='/data1/miRNA/hsa_mature.fa'
ref_other_miRNA='/data1/miRNA/other_mature.fa'
refgenome='/data1/genomics/genome/hg38/hg38.fa'
index_refgenome='/data1/program/db/hsa_genome/hg38'
species='hsa'
adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA 

cmd_dir='/software/mirDeep2/mirdeep2.0.1.2/bin'

#-----------1. align to refgenome by bowtie---------#
outpath=/data1/miRNA
configfile='mapper_config.txt'
filename=ips_exosome_miRNA_all

mkdir $outpath/mapper_out; cd $outpath/mapper_out
$cmd_dir/mapper.pl $outpath/$configfile -d -e -i -j -h -k $adapter -l 18 -m -q -r 5 -p $index_refgenome \
            -s ${filename}_mappered_read.fa -t ${filename}_mapout.rtf -v -n -o 16 2>mapper_report.log
:<<EOF
paremeters:
-e: input fq
-i: convert rna to dna alphabet
-k: clip 3' adatper
-l: min seq len
-m: collapse reads
-p: map to genome
-q: allow 1 mismatch
-r: defualt 5,a read is allowed to map up to this number of positions in the genome
-t: mappping outfile name
-s: processed reads outfile name,  collapsed reads
-v: progress report
-o: number of threads to use for bowtie
EOF


#----2. microRNA detaction and quantification----#
if [ $ec == 0 ];then
	cd ../
	if [ -e miRDeep2 ];then rm -r miRDeep2;fi
	mkdir miRDeep2; cd miRDeep2
	$cmd_dir/miRDeep2.pl ../mapper/${filename}_mappered_read.fa $refgenome ../mapper/${filename}_mapout.rtf $ref_mature_miRNA $ref_other_miRNA  $ref_pre_miRNA -t hsa -v 2>miRDepp2_report.log

fi

:<<EOF
paremeters of miRDeep2.pl
input file：  A fasta file with deep sequencing reads,
                          a fasta file of the corresponding genome,
                          a file of mapped reads to the genome in miRDeep2 arf format,
                          an optional fasta file with known miRNAs of the analysing species，
              and an option fasta file of known miRNAs of related species.

-a int        minimum read stack height that triggers analysis.
              Using this option disables automatic estimation
              of the optimal value.

-t species    species being analyzed - this is used to link to
              the appropriate UCSC browser
-v            remove directory with temporary files
-s file       File with known miRBase star sequences
-P            use this switch if mature_ref_miRNAs contain miRBase
              v18 identifiers (5p and 3p) instead of previous
              ids from v17
EOF








