#!/bin/bash
#
#  Full pipeline to process fastq files   
#
# Usage:
# sh TSS.sh  /path/to/samples_config.csv name_of_job 
#
#
# CHECKLIST BEFORE USING:
#

## - Need to have the appropriate fa files in $HOME/genomes/ref.fa
# - HISAT2 setup in genomes folder

#
# - Installed version of cutadapt with path specified
# - Have RTEventsCounter.py in a folder in $HOME/scripts/RTEventsCounter.py
# - Need to have run R, and have locally installed dplyr and tidyr.


# - If applicable, need to have the appropriate primer files saved in $HOME/primers.dir
#Each line of the file looks like:
# Target,primer number,sequence
# 
# Tips:
#
# - make sure that the primer RNA = fa name = name within the fa file = RNA name in configuration file.
# - paste the csv into a text editor to make sure that there are no hidden characters
# - for a new RNA, paste the .fa file into an editor, make sure it has the correct RNA name
#and then in your genomes folder, use:  hisat2-build file.fa base-name 
#
###################################################################################
# 
#	Files needed:
#	1. Config file (.csv)
#	2. Reference file (.fa, with indexes built by Hisat2 or Bowtie2 as required)
#		Note: as written below, should be in /genomes/ folder
#	3. Primer list (if applicable)
# 		Note: as written below, should be in /primers.dir/ folder
#	4. RTEventsCounter.py
#	5. conf.py (for RTEventsCounter.py)
#	6. 1910_datawrangling_Terbium.R
# 	7. Terbium_data_analysis.R

# Set your home dir:

HOME="/home/NAME"

# Loading FastQC
module purge
module load FastQC/0.11.5-Java-1.8.0_121
module load HISAT2/2.1.0-foss-2016b
module load cutadapt/1.9.1-foss-2016b-Python-2.7.13
# Read in the configuration script.

lines=$(cat $1 | sed 's/ //g' | awk 'NR > 1 {print}') 
config=$1


# Use the job name 
 
mastername=$2


# Collect the files in one master folder:

mkdir $HOME/"$mastername".dir
cd $HOME/"$mastername".dir
mkdir "$mastername".tobackup.dir
cp $HOME/scripts/RTEventsCounter.py .
cp $HOME/scripts/conf.py .
cp $HOME/scripts/1910_datawrangling_Terbium.R .
cp $HOME/scripts/Terbium.sh .
cp $HOME/scripts/"$1" .


# First cycle through and make all the relavent alignments:
#  Note that this section can be commented out if you already have sam files.

for line in $lines
do
    # Set the sample name:

    sample=$(echo "$line" | sed 's/,/, /g' | awk 'BEGIN {FS = ","} {print $1}')
    echo "...processing " $sample
    

    # Set the link to the YCGA files:

    #samplelink=$(echo "$line" | sed 's/,/, /g' | awk 'BEGIN {FS = ","} {print $2}' | \
    samplelink=$(echo "$line" | sed 's/,/, /g' | awk '{print $2}' | sed 's/[;,]/ /g')


    # Set the reference name:

    ref=$(echo "$line" | sed 's/,/, /g' | awk 'BEGIN {FS = ","} {print $3}' | \
	sed 's/ //g')

     
    # Set the primer file:

    primerCoordinates=$(echo "$line" | sed 's/,/, /g' | awk 'BEGIN {FS = ","} {print $4}' | \
	sed 's/ //g')
    

    # Create a subfolder

    mkdir "$sample".dir
    cd "$sample".dir

        # copy data into dir and create individual r1 and r2 fastq files:

        for i in $samplelink; do cp "$i"/*.gz . ; done &&
	gunzip *.gz &&
	cat *_R1_*fastq > "$sample".R1.fastq &&
	cat *_R2_*fastq > "$sample".R2.fastq &&
	
	rm *_R[12]_*fastq # clean up as we go

	#fastqc "$sample"R1.fastq
	#fastqc "$sample"R2.fastq
        
	# trim reads #
	   # -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
           # -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \


        cutadapt \
	    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	    -A AGATCGGAAGAGCGTCGTGTAG \
	    -o "$sample".t.r1.fastq -p "$sample".t.r2.fastq "$sample".R1.fastq "$sample".R2.fastq &&
 
       
	mv "$sample".R[12].fastq ../"$mastername".tobackup.dir/. # save fastq for backup
	mv "$sample".R[12]_fastqc.html ../"$mastername".tobackup.dir/. # save the fastqc files for backup

        # align reads, trim 6 from forward, which contains (N)6
	# reverse (read 2) is reverse complement, starts at 5' end of cDNA (3' end of alignment)
	# forward (read 1) is sense, starts at 3' end of cDNA with 6nt indexed adaptor and N6 UMI
	# bowtie2 \
        #        -p 16 --very-sensitive-local -X 1000 --fr -5 3 \
        #       --no-mixed --no-unal \
        #         -S "$sample".sam -x $ref \
        #       -1 "$sample".t.r1.fastq -2 "$sample".t.r2.fastq &&

	hisat2 \
             -p 1 --fr -5 6 \
              --no-unal --no-mixed --reorder \
             -S "$sample".sam -x /home/sp956/genomes/$ref \
             -1 "$sample".t.r1.fastq -2 "$sample".t.r2.fastq &&

	# now run FastQC on the sam file
	
	#fastqc "$sample".sam

	# align reads to mouse 18s as a proxy for non-specific background
	# add "--no-discordant" after "--no-unal" to remove all reads < 75nt
	# this happens because of read overlap, causing read 2 to align 5' of read 1

	#hisat2 -p 1 --fr -5 6 --no-unal --no-mixed --reorder \            
        #      -S test_18s.sam -x /home/sp956/genomes/s \
        #      -1 "$sample".t.r1.fastq -2 "$sample".t.r2.fastq &&
 
	rm "$sample".t.r[12].fastq  # clean up
	
	# make bam file for backup/storage:

       samtools view -b -S -o "$sample".bam "$sample".sam &&
       samtools sort "$sample".bam "$sample".sorted &&
       samtools index "$sample".sorted.bam
	
	rm "$sample".bam
	mv "$sample".sorted* ../"$mastername".tobackup.dir/.
	mv "$sample"_fastqc.html ../"$mastername".tobackup.dir/.
	cd .. # go back out to the main directory
 

		# NOTE: check conf.py for Random_primers = True/False
	# If using random primers, remove primerCoordinates line below (leave the &&)
	
	python RTEventsCounter.py "$sample".dir/"$sample".sam \
	    ~/genomes/"$ref".fa && #\
	#    ~/primers.dir/"$primerCoordinates" &&
	
	#rm -r "$sample".dir/ #clean up

	echo "done with sample " $sample

done

# backup the fastq files and the bam alignments:

ls "$mastername".tobackup.dir > "$mastername".backupFileList.txt
tar -cvzf "$mastername".tar.gz "$mastername".tobackup.dir

# Load R and process RTEoutput 

module load R

#input in config file and output of RTEventsCounter file (sample SP082220.csv)
    Rscript 1910_datawrangling_Terbium"$mastername" "$mastername".Rout

# Now using file to final processing
Rscript Terbium_data_analysis.R


# And done!

echo "finished with" $mastername
