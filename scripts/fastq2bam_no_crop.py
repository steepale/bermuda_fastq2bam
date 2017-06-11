import sys
import os
import re
from shutil import copyfile

# reference files
config_file = './data/config.txt'
trim_adaptor_file = '/opt/software/Trimmomatic/0.33/adapters/TruSeq3-PE-2.fa'
galgal5_ref = './data/Galgal5/genome.fa'
galgal5_ref_index = galgal5_ref.replace('.fa', '.fa.fai') 
galgal5_ref_dict = galgal5_ref.replace('.fa', '.dict')
dbsnp = './data/dbsnp/snp/organisms/chicken_9031/VCF/all_sorted.vcf'

# Program files
trimmomatic = 'java -Xmx10g -jar $TRIM/trimmomatic'
picard_readgroups = 'java -Xmx10g -jar $PICARD/AddOrReplaceReadGroups.jar'
picard_sortsam = 'java -Xmx10g -jar $PICARD/SortSam.jar'
picard_markdups = 'java -Xmx10g -jar $PICARD/MarkDuplicates.jar'
picard_index = 'java -Xmx10g -jar $PICARD/BuildBamIndex.jar'
picard_dict = 'java -Xmx10g -jar $PICARD/CreateSequenceDictionary.jar'
gatk = 'java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar'

# Input sample
in_sample = sys.argv[1]

# Conditional variable
sample_type = 'whole_genome_sequence_blood'

# Iterate through the config file and grab appropriate sample identifiers and place in dictionary data structure
sample2ids = {}
sample2lanes = {}
sample2run = {}
sample2bcl = {}

# Iterate through config file and fill dictionaries
for line in open(config_file):
	if line[0] != '#':
		line = line.rstrip()
		col = line.split('\t')
		sample_id = col[0]
		date = col[1]
		instrument_id = col[2]
		run_number = col[3]
		flowcell = col[4]
		bcl_conversion = col[5]
		volume = col[8]
		barcode = col[9]
		lane = col[6]
		read = col[7]
		# Different dictionaries and id relationships
		if sample_id not in sample2ids.keys():
			sample2ids[sample_id] = barcode
		if sample_id not in sample2run.keys():
			sample2run[sample_id] = set()
			sample2run[sample_id].add(date+';'+instrument_id+';'+run_number+';'+flowcell)
		elif sample_id in sample2run.keys():
			sample2run[sample_id].add(date+';'+instrument_id+';'+run_number+';'+flowcell)
		if sample_id not in sample2lanes.keys():
			sample2lanes[sample_id] = set()
			sample2lanes[sample_id].add(lane)
		elif sample_id in sample2lanes.keys():
			sample2lanes[sample_id].add(lane)
		if sample_id not in sample2bcl.keys():
			sample2bcl[sample_id] = set()
			sample2bcl[sample_id].add(bcl_conversion)
		elif sample_id in sample2bcl.keys():
			sample2bcl[sample_id].add(bcl_conversion)

# Convert the sets to lists
sample2lanes[in_sample] = list(sample2lanes[in_sample])
sample2run[in_sample] = list(sample2run[in_sample])
sample2bcl[in_sample] = list(sample2bcl[in_sample])
# Create empty set of bams files for each lane to be filled in for loop
bams_per_lane_set = set()

### Run a fastqc analysis on individual sequencing runs before preprocessing to assess read quality and adaptor contamination
for run in sample2run[in_sample]:
	for lane in sample2lanes[in_sample]:
		for bcl_conversion in sample2bcl[in_sample]:
			date = run.split(';')[0]
			instrument_id = run.split(';')[1]
			run_number = run.split(';')[2]
			flowcell = run.split(';')[3]
			barcode = sample2ids[in_sample]
			# copy files to new folder and rename if they exist
			fastq_dir = './dev/fastq_all_samples/'
			r1_file_src = './dev/'+in_sample+'/02-FASTQ/'+date+'_'+instrument_id+'_'+run_number+'_'+flowcell+'/'+in_sample+'_'+bcl_conversion+'_'+lane+'_R1_001.fastq.gz'
			r1_file_dst = fastq_dir+in_sample+'_'+bcl_conversion+'_'+date+'_'+instrument_id+'_'+run_number+'_'+flowcell+'_'+lane+'_R1_001.fastq.gz'
			r2_file_src = './dev/'+in_sample+'/02-FASTQ/'+date+'_'+instrument_id+'_'+run_number+'_'+flowcell+'/'+in_sample+'_'+bcl_conversion+'_'+lane+'_R2_001.fastq.gz'
			r2_file_dst = fastq_dir+in_sample+'_'+bcl_conversion+'_'+date+'_'+instrument_id+'_'+run_number+'_'+flowcell+'_'+lane+'_R2_001.fastq.gz'
			if not os.path.exists(fastq_dir):
				os.makedirs(fastq_dir)
			if os.path.exists(r1_file_src) and os.path.exists(r2_file_src):
				copyfile(r1_file_src, r1_file_dst)
				copyfile(r2_file_src, r2_file_dst)
				out_dir = './analysis/fastqc_before_trim/'
				# Check if output directory exists, if not create it
				if not os.path.exists(out_dir):
					os.makedirs(out_dir)
				# Run FastQC command
				fastqc_cmd = 'fastqc -o'+' '+out_dir+' '+r1_file_dst+' '+r2_file_dst
				print('### Run a fastqc analysis before preprocessing to assess read quality and adaptor contamination')
				os.system(fastqc_cmd)

### Trim off adaptor sequences with trimmomatic and low quality reads
r1_file = "./dev/fastq_merged/"+in_sample+"_R1_merged.fastq.gz"
r2_file = "./dev/fastq_merged/"+in_sample+"_R2_merged.fastq.gz"
trimm_r1_paired_file = './dev/trimmomatic/'+in_sample+'_merged_R1_'+'paired_trimmomatic.fastq.gz'
trimm_r2_paired_file = './dev/trimmomatic/'+in_sample+'_merged_R2_'+'paired_trimmomatic.fastq.gz'
trimm_r1_unpaired_file = './dev/trimmomatic/'+in_sample+'_merged_R1_'+'unpaired_trimmomatic.fastq.gz'
trimm_r2_unpaired_file = './dev/trimmomatic/'+in_sample+'_merged_R2_'+'unpaired_trimmomatic.fastq.gz'
out_dir = './dev/trimmomatic/'
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
# Run command
trimmomatic_cmd = trimmomatic+' '+'PE -threads 2'+' '+r1_file+' '+r2_file+' '+trimm_r1_paired_file+' '+trimm_r1_unpaired_file+' '+trimm_r2_paired_file+' '+trimm_r2_unpaired_file+' '+'ILLUMINACLIP:'+trim_adaptor_file+':2:30:10'
print('### Trim off adaptor sequences with trimmomatic and low quality reads')
os.system(trimmomatic_cmd)

### Perform a FastQC analysis after adapter trimming
out_dir = './analysis/fastqc_post_trimmomatic/'
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
# Run command
fastqc_cmd = 'fastqc -o'+' '+out_dir+' '+trimm_r1_paired_file+' '+trimm_r2_paired_file+' '+trimm_r1_unpaired_file+' '+trimm_r2_unpaired_file
print('### Perform a FastQC analysis after adapter trimming')
os.system(fastqc_cmd)

### Trim off remaining low quality reads with sickle
sickle_r1_paired_file = './dev/sickle/'+in_sample+'_merged_R1_'+'paired_sickle.fastq.gz'
sickle_r2_paired_file = './dev/sickle/'+in_sample+'_merged_R2_'+'paired_sickle.fastq.gz'
sickle_singles_pe_file = './dev/sickle/'+in_sample+'_merged_'+'singles_pe_sickle.fastq.gz'
sickle_singles_se_file = './dev/sickle/'+in_sample+'_merged_'+'singles_se_sickle.fastq.gz'
out_dir = './dev/sickle/'
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
# Run command
sickle_pe_cmd = 'sickle pe -f '+trimm_r1_paired_file+' -r '+trimm_r2_paired_file+' -t sanger -o '+sickle_r1_paired_file+' -p '+sickle_r2_paired_file+' -s '+sickle_singles_pe_file+' -q 20 -l 50 -g'
sickle_se_cmd = 'sickle se -f '+trimm_r1_unpaired_file+\
' -t sanger -o '+sickle_singles_se_file+' -q 30 -l 50 -g'
print('### Trim off remaining low quality reads with sickle')
os.system(sickle_pe_cmd)
os.system(sickle_se_cmd)

### Perform FastQC post read quality trim via sickle
out_dir = './analysis/fastqc_post_sickle/'
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
# Run command
fastqc_cmd = 'fastqc -o'+' '+out_dir+' '+\
sickle_r1_paired_file+' '+sickle_r2_paired_file+' '+\
sickle_singles_pe_file+' '+sickle_singles_se_file
print('### Perform FastQC post read quality trim via sickle')
os.system(fastqc_cmd)
