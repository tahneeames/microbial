
#installing software 

########################

#download miniforge to install with 
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh

########################

#checkM for bacterial tree-based single-copy ortholog analysis
#create checkM env 
conda create -n checkm python=3.9 #create new env named checkm installing python 3.9
conda activate checkm #activate environment 
conda install -c bioconda numpy matplotlib pysam #c bioconda goes to bioconda repository 
conda install -c bioconda hmmer prodigal pplacer
pip3 install checkm-genome
#install databases
wget https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz?download=1

tar -xvzf
#set root 
checkm data setRoot .checkm/

########################

##gtdb_tk for identity /tree building 
#install 
conda create -n gtdbtk-2.7.2 -c conda-forge -c bioconda gtdbtk=2.7.2
conda activate gtdbtk-2.7.2
#run database download 
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
#write path to database
conda env config vars set GTDBTK_DATA_PATH="$home/miniforge3/envs/gtdbtk-2.7.2/share/gtdbtk-2.7.2/db"

########################

#DRAM for functional annotation 
#install DRAM 
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/environment.yaml
conda env create -f environment.yaml -n DRAM

#setup databases
DRAM-setup.py prepare_databases --output_dir DRAM_data --kegg_loc kegg.pep
#seems like the script doesn't work so I need to fix it 
find / -type f -name "DRAM-setup.py"

#replace line with version off github 


########################################################################################################
########################################################################################################

#create new directory for new samples 
mkdir magFiles
cd magFiles
mkdir MAGs
mkdir checkM_magFiles
#grab samples 
sftp #cluster where the MAGs are at 
lcd $home/magFiles
cd #sample directory 
get *.fa #retrieve MAGs for the sample 

#make analyze directory 
mkdir $home/magFiles/checkM_magFiles/analyze

########################

#make job with other samples 
nano checkM.slurm


#!/bin/bash
#
##SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@address.edu 
#SBATCH --partition=normal 
#SBATCH --output=checkM_%j.txt
#SBATCH --error=checkM_%jerr.txt
#SBATCH --job-name="checkM"
##SBATCH --mem=187G
##SBATCH --array
##SBATCH --dependency=[dependency type:job id]
##SBACTH --requeue 

######################################################

conda activate checkm 

checkm lineage_wf -t 8 -x fa $home/magFiles/MAGs $home/magFiles/checkM_magFiles #metagenome file location - output directory

checkm analyze -t 8 -x fa $home/magFiles/checkM_magFiles/lineage.ms $home/magFiles/MAGs $home/magFiles/checkM_magFiles/analyze #lineage file from previous step - metagenome file location - output folder

checkm qa -f checkM_qa.tsv --tab_table -t 8 $home/magFiles/checkM_magFiles/lineage.ms $home/magFiles/checkM_magFiles/analyze #lineage file from lineage_wf step - analyze folder output from analyze step 


########################


#run 
sbatch checkM.slurm


########################################################################


#processing QA into high, med, and low quality bins for downstream steps 
cd $home/magFiles
#make new folders to sort checkM files 
mkdir high
mkdir med 
mkdir low  

#sort the high >90% complete <5% contaminated 
awk -F'\t' '$12 > 89.999 && $13 < 5.001' checkM_qa.tsv > checkMqa_high.tsv
awk -F'\t' '{print $1 ".fa"}' checkMqa_high.tsv > checkMqa_highList.txt 

#medium 50-89% complete <10% contaminated 
awk -F'\t' '$12 > 49.999 && $12 < 89.999 && $13 < 10.001' checkM_qa.tsv > checkMqa_med.tsv
awk -F'\t' '{print $1 ".fa"}' checkMqa_med.tsv > checkMqa_medList.txt 

#low <50% complete >10% contaminated 
awk -F'\t' '$12 < 49.999 || $13 > 9.999' checkM_qa.tsv > checkMqa_low.tsv
awk -F'\t' '{print $1 ".fa"}' checkMqa_low.tsv > checkMqa_lowList.txt 

#move the files from the MAGs directory into these folders based on above crireria 
cd MAGs 
while read file; do mv "$file" $home/magFiles/high; done < $home/magFiles/checkMqa_highList.txt 
while read file; do mv "$file" $home/magFiles/med; done < $home/magFiles/checkMqa_medList.txt  
while read file; do mv "$file" $home/magFiles/low; done < $home/magFiles/checkMqa_lowList.txt 


########################################################################


#GTDB-TK for bacterial classification 

#new directory for it in magFiles directory 
#make one for the high-quality samples 
cd $home/magFiles
mkdir gtdbtk_high  
cd gtdbtk_high
mkdir identify 
mkdir align 
mkdir classify 

#make another for medium quality samples 
cd $home/magFiles
mkdir gtdbtk_med
cd gtdbtk_med
mkdir identify 
mkdir align 
mkdir classify 

########################

#slurm script for high quality samples 
nano gtdbtk_high.slurm

#!/bin/bash
#
##SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@address.edu 
#SBATCH --partition=normal 
#SBATCH --output=gtdb_high_%j.txt
#SBATCH --error=gtdb_high_%jerr.txt
#SBATCH --job-name="gtdb_tk_high"
##SBATCH --mem=187G
##SBATCH --array
##SBATCH --dependency=[dependency type:job id]
##SBATCH --requeue 

######################################################

conda activate gtdbtk-2.7.2

gtdbtk identify --genome_dir $home/magFiles/high --out_dir $home/magFiles/gtdbtk_high/identify --extension fa --cpus 8

gtdbtk align --identify_dir $home/magFiles/gtdbtk_high/identify --out_dir $home/magFiles/gtdbtk_high/align --cpus 8

gtdbtk classify --genome_dir $home/magFiles/high --align_dir $home/magFiles/gtdbtk_high/align --out_dir $home/magFiles/gtdbtk_high/classify --extension fa --cpus 8


########################

#run 
sbatch gtdbtk_high.slurm


########################

#for medium quality
nano gtdbtk_med.slurm

#!/bin/bash
#
##SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@address.edu 
#SBATCH --partition=normal 
#SBATCH --output=gtdb_med_%j.txt
#SBATCH --error=gtdb_med_%jerr.txt
#SBATCH --job-name="gtdb_tk_MED"
##SBATCH --mem=187G
##SBATCH --array
##SBATCH --dependency=[dependency type:job id]
##SBATCH --requeue 

######################################################

conda activate gtdbtk-2.7.2

gtdbtk identify --genome_dir $home/magFiles/med --out_dir $home/magFiles/gtdbtk_med/identify --extension fa --cpus 8

gtdbtk align --identify_dir $home/magFiles/gtdbtk_med/identify --out_dir $home/magFiles/gtdbtk_med/align --cpus 8

gtdbtk classify --genome_dir $home/magFiles/med --align_dir $home/magFiles/gtdbtk_med/align --out_dir $home/magFiles/gtdbtk_med/classify --extension fa --cpus 8


















