#!/bin/bash

#Written RAWIII
#Dec 17, 2018

#SBATCH --partition=cahnrs_bigmem		    ### Partition 
#SBATCH --job-name=pipeline	    		    ### Job Name
#SBATCH --output=pipeline_181217_out.txt    	### File in which to store job output
#SBATCH --error=pipeline_181217_err.txt     	### File in which to store job error messages
#SBATCH --time=7-0:00:00       			    ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               		    ### Node count required for the job
#SBATCH --nodes=1               			### Node count required for the job
#SBATCH --cpus-per-task=60

#load modules
module load miniconda3
source activate /opt/apps/conda/envs/megahit

#remove adapters/trim
bbduk.sh -Xmx1g in1=F1_all_R1.fastq.gz in2=F1_all_R2.fastq.gz out1=F1_all_trim_R1.fastq out2=F1_all_trim_R2.fastq ref=~/bbmap/resources/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo 
bbduk.sh -Xmx1g in1=F2_all_R1.fastq.gz in2=F2_all_R2.fastq.gz out1=F2_all_trim_R1.fastq out2=F2_all_trim_R2.fastq ref=~/bbmap/resources/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo
bbduk.sh -Xmx1g in1=F3_all_R1.fastq.gz in2=F3_all_R2.fastq.gz out1=F3_all_trim_R1.fastq out2=F3_all_trim_R2.fastq ref=~/bbmap/resources/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo 
bbduk.sh -Xmx1g in1=F4_all_R1.fastq.gz in2=F4_all_R2.fastq.gz out1=F4_all_trim_R1.fastq out2=F4_all_trim_R2.fastq ref=~/bbmap/resources/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo
bbduk.sh -Xmx1g in1=U1_all_R1.fastq.gz in2=U1_all_R2.fastq.gz out1=U1_all_trim_R1.fastq out2=U1_all_trim_R2.fastq ref=~/bbmap/resources/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo 
bbduk.sh -Xmx1g in1=U2_all_R1.fastq.gz in2=U2_all_R2.fastq.gz out1=U2_all_trim_R1.fastq out2=U2_all_trim_R2.fastq ref=~/bbmap/resources/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo
bbduk.sh -Xmx1g in1=U3_all_R1.fastq.gz in2=U3_all_R2.fastq.gz out1=U3_all_trim_R1.fastq out2=U3_all_trim_R2.fastq ref=~/bbmap/resources/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo 
bbduk.sh -Xmx1g in1=U4_all_R1.fastq.gz in2=U4_all_R2.fastq.gz out1=U4_all_trim_R1.fastq out2=U4_all_trim_R2.fastq ref=~/bbmap/resources/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo

#decon remove phiX
bbduk.sh -Xmx1g in1=F1_all_trim_R1.fastq in2=F1_all_trim_R1.fastq out1=F1_all_decon_R1.fastq out2=F1_all_decon_R2.fastq qtrim=r trimq=25 maq=25 minlen=50 ref=~/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=F1_all_trim-decon-stats.txt
bbduk.sh -Xmx1g in1=F2_all_trim_R1.fastq in2=F2_all_trim_R1.fastq out1=F2_all_decon_R1.fastq out2=F2_all_decon_R2.fastq qtrim=r trimq=25 maq=25 minlen=50 ref=~/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=F2_all_trim-decon-stats.txt
bbduk.sh -Xmx1g in1=F3_all_trim_R1.fastq in2=F3_all_trim_R1.fastq out1=F3_all_decon_R1.fastq out2=F3_all_decon_R2.fastq qtrim=r trimq=25 maq=25 minlen=50 ref=~/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=F3_all_trim-decon-stats.txt
bbduk.sh -Xmx1g in1=F4_all_trim_R1.fastq in2=F4_all_trim_R1.fastq out1=F4_all_decon_R1.fastq out2=F4_all_decon_R2.fastq qtrim=r trimq=25 maq=25 minlen=50 ref=~/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=F4_all_trim-decon-stats.txt
bbduk.sh -Xmx1g in1=U1_all_trim_R1.fastq in2=U1_all_trim_R1.fastq out1=U1_all_decon_R1.fastq out2=U1_all_decon_R2.fastq qtrim=r trimq=25 maq=25 minlen=50 ref=~/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=U1_all_trim-decon-stats.txt
bbduk.sh -Xmx1g in1=U2_all_trim_R1.fastq in2=U2_all_trim_R1.fastq out1=U2_all_decon_R1.fastq out2=U2_all_decon_R2.fastq qtrim=r trimq=25 maq=25 minlen=50 ref=~/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=U2_all_trim-decon-stats.txt
bbduk.sh -Xmx1g in1=U3_all_trim_R1.fastq in2=U3_all_trim_R1.fastq out1=U3_all_decon_R1.fastq out2=U3_all_decon_R2.fastq qtrim=r trimq=25 maq=25 minlen=50 ref=~/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=U3_all_trim-decon-stats.txt
bbduk.sh -Xmx1g in1=U4_all_trim_R1.fastq in2=U4_all_trim_R1.fastq out1=U4_all_decon_R1.fastq out2=U4_all_decon_R2.fastq qtrim=r trimq=25 maq=25 minlen=50 ref=~/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=U4_all_trim-decon-stats.txt

#count reads in raw files
awk '{s++}END{print s/4}' F1_all_R1.fastq
awk '{s++}END{print s/4}' F2_all_R1.fastq
awk '{s++}END{print s/4}' F3_all_R1.fastq
awk '{s++}END{print s/4}' F4_all_R1.fastq
awk '{s++}END{print s/4}' U1_all_R1.fastq
awk '{s++}END{print s/4}' U2_all_R1.fastq
awk '{s++}END{print s/4}' U3_all_R1.fastq
awk '{s++}END{print s/4}' U4_all_R1.fastq
awk '{s++}END{print s/4}' F1_all_R2.fastq
awk '{s++}END{print s/4}' F2_all_R2.fastq
awk '{s++}END{print s/4}' F3_all_R2.fastq
awk '{s++}END{print s/4}' F4_all_R2.fastq
awk '{s++}END{print s/4}' U1_all_R2.fastq
awk '{s++}END{print s/4}' U2_all_R2.fastq
awk '{s++}END{print s/4}' U3_all_R2.fastq
awk '{s++}END{print s/4}' U4_all_R2.fastq

#count reads in trim files
awk '{s++}END{print s/4}' F1_all_trim_R1.fastq
awk '{s++}END{print s/4}' F2_all_trim_R1.fastq
awk '{s++}END{print s/4}' F3_all_trim_R1.fastq
awk '{s++}END{print s/4}' F4_all_trim_R1.fastq
awk '{s++}END{print s/4}' U1_all_trim_R1.fastq
awk '{s++}END{print s/4}' U2_all_trim_R1.fastq
awk '{s++}END{print s/4}' U3_all_trim_R1.fastq
awk '{s++}END{print s/4}' U4_all_trim_R1.fastq
awk '{s++}END{print s/4}' F1_all_trim_R2.fastq
awk '{s++}END{print s/4}' F2_all_trim_R2.fastq
awk '{s++}END{print s/4}' F3_all_trim_R2.fastq
awk '{s++}END{print s/4}' F4_all_trim_R2.fastq
awk '{s++}END{print s/4}' U1_all_trim_R2.fastq
awk '{s++}END{print s/4}' U2_all_trim_R2.fastq
awk '{s++}END{print s/4}' U3_all_trim_R2.fastq
awk '{s++}END{print s/4}' U4_all_trim_R2.fastq


#count reads in decon files
awk '{s++}END{print s/4}' F1_all_decon_R1.fastq
awk '{s++}END{print s/4}' F2_all_decon_R1.fastq
awk '{s++}END{print s/4}' F3_all_decon_R1.fastq
awk '{s++}END{print s/4}' F4_all_decon_R1.fastq
awk '{s++}END{print s/4}' U1_all_decon_R1.fastq
awk '{s++}END{print s/4}' U2_all_decon_R1.fastq
awk '{s++}END{print s/4}' U3_all_decon_R1.fastq
awk '{s++}END{print s/4}' U4_all_decon_R1.fastq
awk '{s++}END{print s/4}' F1_all_decon_R2.fastq
awk '{s++}END{print s/4}' F2_all_decon_R2.fastq
awk '{s++}END{print s/4}' F3_all_decon_R2.fastq
awk '{s++}END{print s/4}' F4_all_decon_R2.fastq
awk '{s++}END{print s/4}' U1_all_decon_R2.fastq
awk '{s++}END{print s/4}' U2_all_decon_R2.fastq
awk '{s++}END{print s/4}' U3_all_decon_R2.fastq
awk '{s++}END{print s/4}' U4_all_decon_R2.fastq

#run megahit assembly
megahit -m 0.99 -l 350 -1 F1_all_decon_R1.fastq -2 F1_all_decon_R2.fastq --k-min 21 --k-max 121 --out-dir F1_megahit_decon_dir
megahit -m 0.99 -l 350 -1 F2_all_decon_R1.fastq -2 F2_all_decon_R2.fastq --k-min 21 --k-max 121 --out-dir F2_megahit_decon_dir
megahit -m 0.99 -l 350 -1 F3_all_decon_R1.fastq -2 F3_all_decon_R2.fastq --k-min 21 --k-max 121 --out-dir F3_megahit_decon_dir
megahit -m 0.99 -l 350 -1 F4_all_decon_R1.fastq -2 F4_all_decon_R2.fastq --k-min 21 --k-max 121 --out-dir F4_megahit_decon_dir
megahit -m 0.99 -l 350 -1 U1_all_decon_R1.fastq -2 U1_all_decon_R2.fastq --k-min 21 --k-max 121 --out-dir U1_megahit_decon_dir
megahit -m 0.99 -l 350 -1 U2_all_decon_R1.fastq -2 U2_all_decon_R2.fastq --k-min 21 --k-max 121 --out-dir U2_megahit_decon_dir
megahit -m 0.99 -l 350 -1 U3_all_decon_R1.fastq -2 U3_all_decon_R2.fastq --k-min 21 --k-max 121 --out-dir U3_megahit_decon_dir
megahit -m 0.99 -l 350 -1 U4_all_decon_R1.fastq -2 U4_all_decon_R2.fastq --k-min 21 --k-max 121 --out-dir U4_megahit_decon_dir

#copy raw contigs and subsample >5k
mkdir contigs_5k
cd F1_megahit_decon_dir
cp *fasta* >F1_megahit_raw-assem.fasta
reformat.sh in=F1_megahit_raw-assem.fasta out=F1_megahit_5k.fasta minlength=5000 aqhist=F1_megahit_hist_5k.txt
mv *5k* ../contigs_5k
cd ../F2_megahit_decon_dir
cp *fasta* >F2_megahit_raw-assem.fasta
reformat.sh in=F2_megahit_raw-assem.fasta out=F2_megahit_5k.fasta minlength=5000 aqhist=F2_megahit_hist_5k.txt
mv *5k* ../contigs_5k
cd ../F3_megahit_decon_dir
cp *fasta* >F3_megahit_raw-assem.fasta
reformat.sh in=F3_megahit_raw-assem.fasta out=F3_megahit_5k.fasta minlength=5000 aqhist=F3_megahit_hist_5k.txt
mv *5k* ../contigs_5k
cd ../F4_megahit_decon_dir
cp *fasta* >F4_megahit_raw-assem.fasta
reformat.sh in=F4_megahit_raw-assem.fasta out=F4_megahit_5k.fasta minlength=5000 aqhist=F4_megahit_hist_5k.txt
mv *5k* ../contigs_5k
cd ../U1_megahit_decon_dir
cp *fasta* >U1_megahit_raw-assem.fasta
reformat.sh in=U1_megahit_raw-assem.fasta out=U1_megahit_5k.fasta minlength=5000 aqhist=U1_megahit_hist_5k.txt
mv *5k* ../contigs_5k
cd ../U2_megahit_decon_dir
cp *fasta* >U2_megahit_raw-assem.fasta
reformat.sh in=U2_megahit_raw-assem.fasta out=U2_megahit_5k.fasta minlength=5000 aqhist=U2_megahit_hist_5k.txt
mv *5k* ../contigs_5k
cd ../U3_megahit_decon_dir
cp *fasta* >U3_megahit_raw-assem.fasta
reformat.sh in=U3_megahit_raw-assem.fasta out=U3_megahit_5k.fasta minlength=5000 aqhist=U3_megahit_hist_5k.txt
mv *5k* ../contigs_5k
cd ../U4_megahit_decon_dir
cp *fasta* >U4_megahit_raw-assem.fasta
reformat.sh in=U4_megahit_raw-assem.fasta out=U4_megahit_5k.fasta minlength=5000 aqhist=U4_megahit_hist_5k.txt
mv *5k* ../contigs_5k
cd ../contigs_5k

mv *U* >U_all_5k.fasta 
mv *F* >F_all_5k.fasta 

#load metawrap
module load miniconda3
source activate /opt/apps/conda/envs/metawrap/

#Run metawrap binning and refinement
#binning
metawrap binning -a ../assemblies/subsampled/U_all_5k.fasta -o U_metawrap -t 28 --metabat2 --maxbin2 --concoct ../../data/decon/U1_all_decon_R1.fastq ../../data/decon/U1_all_decon_R2.fastq ../../data/decon/U1_all_decon_R1.fastq ../../data/decon/U2_all_decon_R1.fastq ../../data/decon/U2_all_decon_R2.fastq ../../data/decon/U3_all_decon_R1.fastq ../../data/decon/U3_all_decon_R2.fastq ../../data/decon/U4_all_decon_R1.fastq ../../data/decon/U4_all_decon_R2.fastq
metawrap binning -a ../assemblies/subsampled/F_all_5k.fasta -o F_metawrap -t 28 --metabat2 --maxbin2 --concoct ../../data/decon/F1_all_decon_R1.fastq ../../data/decon/F1_all_decon_R2.fastq ../../data/decon/F1_all_decon_R1.fastq ../../data/decon/F2_all_decon_R1.fastq ../../data/decon/F2_all_decon_R2.fastq ../../data/decon/F3_all_decon_R1.fastq ../../data/decon/F3_all_decon_R2.fastq ../../data/decon/F4_all_decon_R1.fastq ../../data/decon/F4_all_decon_R2.fastq

#Refine bins
metawrap bin_refinement -o bin_refinement -t 28 -A U_metawrap/metabat2_bins/ -B U_metawrap/maxbin2_bins/ -C U_metawrap/concoct_bins/ -c 5 -x 10
metawrap bin_refinement -o bin_refinement -t 28 -A F_metawrap/metabat2_bins/ -B F_metawrap/maxbin2_bins/ -C F_metawrap/concoct_bins/ -c 5 -x 10

#Load Prokka
module load miniconda3
source activate /opt/apps/conda/envs/prokka

#Run Prokka
prokka F1_megahit_5k.fasta --cpus 28 --outdir F1_megahit_5k --prefix F1_megahit_5k --rfam --metagenome
prokka F2_megahit_5k.fasta --cpus 28 --outdir F2_megahit_5k --prefix F2_megahit_5k --rfam --metagenome
prokka F3_megahit_5k.fasta --cpus 28 --outdir F3_megahit_5k --prefix F3_megahit_5k --rfam --metagenome
prokka F4_megahit_5k.fasta --cpus 28 --outdir F4_megahit_5k --prefix F4_megahit_5k --rfam --metagenome
prokka U1_megahit_5k.fasta --cpus 28 --outdir U1_megahit_5k --prefix U1_megahit_5k --rfam --metagenome
prokka U2_megahit_5k.fasta --cpus 28 --outdir U2_megahit_5k --prefix U2_megahit_5k --rfam --metagenome
prokka U3_megahit_5k.fasta --cpus 28 --outdir U3_megahit_5k --prefix U3_megahit_5k --rfam --metagenome
prokka U4_megahit_5k.fasta --cpus 28 --outdir U4_megahit_5k --prefix U4_megahit_5k --rfam --metagenome

#EggNOG mapper
python /data/friesen/bin/eggnog-mapper/emapper.py -i F1_megahit_5k.faa --output F1_megahit_5k_eggnog -m diamond --cpu 28 --annotate_hits_table --go_evidence
python /data/friesen/bin/eggnog-mapper/emapper.py -i F2_megahit_5k.faa --output F2_megahit_5k_eggnog -m diamond --cpu 28 --annotate_hits_table --go_evidence
python /data/friesen/bin/eggnog-mapper/emapper.py -i F3_megahit_5k.faa --output F3_megahit_5k_eggnog -m diamond --cpu 28 --annotate_hits_table --go_evidence
python /data/friesen/bin/eggnog-mapper/emapper.py -i F4_megahit_5k.faa --output F4_megahit_5k_eggnog -m diamond --cpu 28 --annotate_hits_table --go_evidence
python /data/friesen/bin/eggnog-mapper/emapper.py -i U1_megahit_5k.faa --output U1_megahit_5k_eggnog -m diamond --cpu 28 --annotate_hits_table --go_evidence
python /data/friesen/bin/eggnog-mapper/emapper.py -i U2_megahit_5k.faa --output U2_megahit_5k_eggnog -m diamond --cpu 28 --annotate_hits_table --go_evidence
python /data/friesen/bin/eggnog-mapper/emapper.py -i U3_megahit_5k.faa --output U3_megahit_5k_eggnog -m diamond --cpu 28 --annotate_hits_table --go_evidence
python /data/friesen/bin/eggnog-mapper/emapper.py -i U4_megahit_5k.faa --output U4_megahit_5k_eggnog -m diamond --cpu 28 --annotate_hits_table --go_evidence

#GTDHTK bin assignment
module load miniconda3
source activate /opt/apps/conda/envs/gtdbtk
export GTDBTK_DATA_PATH=/opt/apps/data/gtdbtk/release86/

gtdbtk classify_wf --cpus 28 -x fasta --genome_dir /data/friesen/clean/CQ/results/bins/refinement/refined --out_dir gtdbtk_output

#Run CAZy on contigs
diamond blastp --db /data/friesen/bin/databases/cazy/CAZY_2018_07_31_diamond.dmnd --query F1_megahit_5k.faa -o F1_megahit_5k_cazy18 -t /data/friesen/clean/CQ --more-sensitive -f 6 -k 5 -b 1 -p 24
diamond blastp --db /data/friesen/bin/databases/cazy/CAZY_2018_07_31_diamond.dmnd --query F2_megahit_5k.faa -o F2_megahit_5k_cazy18 -t /data/friesen/clean/CQ --more-sensitive -f 6 -k 5 -b 1 -p 24
diamond blastp --db /data/friesen/bin/databases/cazy/CAZY_2018_07_31_diamond.dmnd --query F3_megahit_5k.faa -o F3_megahit_5k_cazy18 -t /data/friesen/clean/CQ --more-sensitive -f 6 -k 5 -b 1 -p 24
diamond blastp --db /data/friesen/bin/databases/cazy/CAZY_2018_07_31_diamond.dmnd --query F4_megahit_5k.faa -o F4_megahit_5k_cazy18 -t /data/friesen/clean/CQ --more-sensitive -f 6 -k 5 -b 1 -p 24
diamond blastp --db /data/friesen/bin/databases/cazy/CAZY_2018_07_31_diamond.dmnd --query U1_megahit_5k.faa -o U1_megahit_5k_cazy18 -t /data/friesen/clean/CQ --more-sensitive -f 6 -k 5 -b 1 -p 24
diamond blastp --db /data/friesen/bin/databases/cazy/CAZY_2018_07_31_diamond.dmnd --query U2_megahit_5k.faa -o U2_megahit_5k_cazy18 -t /data/friesen/clean/CQ --more-sensitive -f 6 -k 5 -b 1 -p 24
diamond blastp --db /data/friesen/bin/databases/cazy/CAZY_2018_07_31_diamond.dmnd --query U3_megahit_5k.faa -o U3_megahit_5k_cazy18 -t /data/friesen/clean/CQ --more-sensitive -f 6 -k 5 -b 1 -p 24
diamond blastp --db /data/friesen/bin/databases/cazy/CAZY_2018_07_31_diamond.dmnd --query U4_megahit_5k.faa -o U4_megahit_5k_cazy18 -t /data/friesen/clean/CQ --more-sensitive -f 6 -k 5 -b 1 -p 24
