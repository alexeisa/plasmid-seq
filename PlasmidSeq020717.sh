#!/bin/sh

### PlasmidSeq script by A.Saveliev, 2017

### Assumes that bbmap is in PATH, NOVOPlasty in home directory, denovo seed fasta in /media/sf_vm_shared/seeds

### Script must find project folders in /media/sf_vm_shared/PlasmidSeq directory. Each folder must have 3 files: one plasmid .fa, a pair of MiSeq 2x250 bp fastq.gz's

cd /media/sf_vm_shared/PlasmidSeq/

### set up loop
for i in $(ls /media/sf_vm_shared/PlasmidSeq/)

do
 
  cd /media/sf_vm_shared/PlasmidSeq/$i/

  # unset variables
  FBLEN=""
  FBSEQ=""
  LOTNUMBER=""
  SAMPLENUMBER=""
  READ1=""
  READ2=""
  
  # remove BBMap fasta index if exists 
  rm -rf ref
  
  # assign 2 variables  based on fastq.gz names
  LOTNUMBER=$(echo *.gz | awk -F'[_]' '{print $1}')
  SAMPLENAME=$(echo *.gz | awk -F'[_]' '{print $1 "_" $2}')
  
  # create a new fasta reference file with a simple name
  cp *.fa $LOTNUMBER.fasta
  
  # extract fasta sequence, remove control characters, write fasta body to file
  sed  '1d' $LOTNUMBER.fasta | tr -d '\n \r' > $LOTNUMBER.body.fasta
  
  # assign variable to fasta body DNA length
  FBLEN=$(cat $LOTNUMBER.body.fasta | wc -m)
  
  # assign variable to fasta body DNA sequence
  FBSEQ=$(cat $LOTNUMBER.body.fasta)
  
  # shift fasta by 500 bp (move tail to head)
  rm -f $LOTNUMBER.movetail.fasta
  grep '>' $LOTNUMBER.fasta > $LOTNUMBER.movetail.fasta
  echo ${FBSEQ:$FBLEN-500:$FBLEN} | tr -d '\n \r' >> $LOTNUMBER.movetail.fasta
  echo ${FBSEQ:0:$FBLEN-500} >> $LOTNUMBER.movetail.fasta
  
  # perform read pre-processing in BBDuk
  bash /home/alex/bbmap/bbduk.sh overwrite=t ref=/home/alex/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 in1="$SAMPLENAME"_L001_R1_001.fastq.gz in2="$SAMPLENAME"_L001_R2_001.fastq.gz   out1="$SAMPLENAME"_L001_R1_001.clean.fastq.gz out2="$SAMPLENAME"_L001_R2_001.clean.fastq.gz minlen=250
  
  # gunzip pre-processed fastq files, use them for alignment and denovo assembly
  gunzip -kf *.clean.fastq.gz
  
  # perform alignment in BBMap, write ref, sam, bam, bai. Change maxindel value as needed. 
  bash /home/alex/bbmap/bbmap.sh overwrite=t in1="$SAMPLENAME"_L001_R1_001.clean.fastq.gz in2="$SAMPLENAME"_L001_R2_001.clean.fastq.gz ref=$LOTNUMBER.movetail.fasta outm=$LOTNUMBER.movetail.mapped.sam maxindel=2000 outu=$LOTNUMBER.movetail.unmapped.sam 
  
  echo ""
  echo -e "\e[1;34mStarting bamscript from BBTools as embedded...\e[0m"
  echo ""
  
  #BBTools bs.sh inserted below
  
  echo "Note: This script is designed to run with the amount of memory detected by BBMap."
  echo "      If Samtools crashes, please ensure you are running on the same platform as BBMap,"
  echo "      or reduce Samtools' memory setting (the -m flag)."
  echo "Note: Please ignore any warnings about 'EOF marker is absent'; this is a bug in samtools that occurs when using piped input."
  samtools view -bShu $LOTNUMBER.movetail.mapped.sam | samtools sort -m 2G -@ 3 - $LOTNUMBER.movetail.mapped_sorted
  samtools index $LOTNUMBER.movetail.mapped_sorted.bam
  echo "Note: Please ignore any warnings about 'EOF marker is absent'; this is a bug in samtools that occurs when using piped input."
  samtools view -bShu $LOTNUMBER.movetail.unmapped.sam | samtools sort -m 2G -@ 3 - $LOTNUMBER.movetail.unmapped_sorted
  samtools index $LOTNUMBER.movetail.unmapped_sorted.bam
 
 #BBTools bs.sh ends
  
  echo ""
  echo -e "\e[1;34mFinished bamscript from BBTools\e[0m"
  echo ""
	
  # perform variant calling in CallVariants, write as VCF
  bash /home/alex/bbmap/callvariants.sh overwrite=t in=$LOTNUMBER.movetail.mapped.sam ref=$LOTNUMBER.movetail.fasta out=$LOTNUMBER.movetail.vcf rarity=0.01 maf=0.01

  ### process VCF, build simplified variant table
  
  # extract VCF header, write to file
  head -n +52 $LOTNUMBER.movetail.vcf > $LOTNUMBER.movetail.header.txt
  
  # extract VCF body, write to file
  tail -n +53 $LOTNUMBER.movetail.vcf > $LOTNUMBER.movetail.body.txt
  
  # restore original fasta coordinates, write to new txt file as VCF body
  awk '$2>0 {$2=$2-500} 1'  OFS='\t' $LOTNUMBER.movetail.body.txt | awk '$2<1 {$2='$FBLEN'+$2} 1' OFS='\t' > $LOTNUMBER.unshift.vcf.txt
  
  # append VCF header and new rearranged VCF body, write to rearranged file. This is to archive a full VCF with the restored coordinates.
  cat $LOTNUMBER.movetail.header.txt $LOTNUMBER.unshift.vcf.txt > $LOTNUMBER.unshift.vcf
  
  ### process new rearranged VCF body
  
  # replace semicolons with tabs, write to tab-delimited txt file
  tr ';' '\t' <  $LOTNUMBER.unshift.vcf.txt >  $LOTNUMBER.unshift.vcf.tab.txt
  
  # extract specific columns, write to a new txt file
  cut -f1-2,4-6,11,30-34  $LOTNUMBER.unshift.vcf.tab.txt >  $LOTNUMBER.unshift.vcf.tab.simple.txt
  
  # sort by revised allele frequency, write to a new txt file
  sort -rk 9 $LOTNUMBER.unshift.vcf.tab.simple.txt > $LOTNUMBER.unshift.vcf.tab.simple.sortedbyRAF.txt
  
  ### build final report
  
  # append new column names to report 
  echo -e "PLASMID\tPOS\tREF\tALT\tQUAL\tTYPE\tDEPTH\tAF\tRAF\tSB\tDP4" >> $LOTNUMBER.FINAL.report.txt
  
  # append new variant body, sorted by revised allele frequency 
  cat $LOTNUMBER.unshift.vcf.tab.simple.sortedbyRAF.txt >> $LOTNUMBER.FINAL.report.txt
  
  # remove keys, keep numbers only
  sed -i 's/TYP=//g; s/DP=//g; s/RAF=//g; s/AF=//g; s/SB=//g'  $LOTNUMBER.FINAL.report.txt 
 

  # comment out the next line to debug
  #rm -f $LOTNUMBER.unshift.vcf.txt $LOTNUMBER.unshift.vcf.tab.txt $LOTNUMBER.unshift.vcf.tab.simple.txt $LOTNUMBER.unshift.vcf.tab.simple.sortedbyRAF.txt $LOTNUMBER.movetail.header.txt $LOTNUMBER.movetail.body.txt $LOTNUMBER.body.fasta
 
  # remove BBMap index 
  rm -rf ref
  
  ### start denovo assembly in NOVOPlasty
  
  # create variables from pre-processed fasta file names

  READ1=$(basename *_R1_001.clean.fastq)
  READ2=$(basename *_R2_001.clean.fastq)
  
  # clean up and create a new config file for NOVOPlasty
 
  rm -f config_plasmid.txt
  touch config_plasmid.txt
  
  # assign NOVOPLasty config values, seed location, preprocessed fasta files, append all to config_plasmid.txt file
  echo "Project name         = $LOTNUMBER" > config_plasmid.txt
  echo "Insert size          = 350" >> config_plasmid.txt
  echo "Insert size aut      = yes" >> config_plasmid.txt
  echo "Read Length          = 250" >> config_plasmid.txt
  echo "Type                 = mito" >> config_plasmid.txt
  echo "Genome Range         = "$((FBLEN-1000))" "-" "$((FBLEN+1000))"" >> config_plasmid.txt
  echo "K-mer                = 39"  >> config_plasmid.txt
  echo "Insert Range         = 1.6" >> config_plasmid.txt
  echo "Insert Range strict  = 1.2" >> config_plasmid.txt
  echo "Single/Paired        = PE" >> config_plasmid.txt
  echo "Max memory           = " >> config_plasmid.txt
  echo "Coverage Cut off     = 1000" >> config_plasmid.txt
  echo "Extended log         = 0" >> config_plasmid.txt
  echo "Combined reads       = " >> config_plasmid.txt
  echo "Forward reads        = " $READ1 >> config_plasmid.txt
  echo "Reverse reads        = " $READ2 >> config_plasmid.txt
  echo "Seed Input           = /media/sf_vm_shared/seeds/seed.fasta" >> config_plasmid.txt
  
  # run NOVOPlasty
  perl  /home/alex/NOVOPlasty/NOVOPlasty2.6.1.pl -c config_plasmid.txt


done

echo ""
echo "==============ANALYSIS COMPLETED==============="
echo ""