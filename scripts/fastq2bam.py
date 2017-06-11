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
trimmomatic = 'java -Xmx40g -jar $TRIM/trimmomatic'
picard_readgroups = 'java -Xmx40g -jar $PICARD/AddOrReplaceReadGroups.jar'
picard_sortsam = 'java -Xmx40g -jar $PICARD/SortSam.jar'
picard_markdups = 'java -Xmx40g -jar $PICARD/MarkDuplicates.jar'
picard_index = 'java -Xmx40g -jar $PICARD/BuildBamIndex.jar'
picard_dict = 'java -Xmx40g -jar $PICARD/CreateSequenceDictionary.jar'
gatk = 'java -Xmx40g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar'

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
			fastq_dir = './data/fastq_all_samples/'
			r1_file_src = './data/'+in_sample+'/02-FASTQ/'+date+'_'+instrument_id+'_'+run_number+'_'+flowcell+'/'+in_sample+'_'+bcl_conversion+'_'+lane+'_R1_001.fastq.gz'
			r1_file_dst = fastq_dir+in_sample+'_'+bcl_conversion+'_'+date+'_'+instrument_id+'_'+run_number+'_'+flowcell+'_'+lane+'_R1_001.fastq.gz'
			r2_file_src = './data/'+in_sample+'/02-FASTQ/'+date+'_'+instrument_id+'_'+run_number+'_'+flowcell+'/'+in_sample+'_'+bcl_conversion+'_'+lane+'_R2_001.fastq.gz'
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
				print('\n'+'### Run a fastqc analysis before preprocessing to assess read quality and adaptor contamination')
				os.system(fastqc_cmd)

### Trim off adaptor sequences with trimmomatic and low quality reads
r1_file = "./data/fastq_merged/"+in_sample+"_R1_merged.fastq.gz"
r2_file = "./data/fastq_merged/"+in_sample+"_R2_merged.fastq.gz"
trimm_r1_paired_file = './data/trimmomatic/'+in_sample+'_merged_R1_'+'paired_trimmomatic.fastq.gz'
trimm_r2_paired_file = './data/trimmomatic/'+in_sample+'_merged_R2_'+'paired_trimmomatic.fastq.gz'
trimm_r1_unpaired_file = './data/trimmomatic/'+in_sample+'_merged_R1_'+'unpaired_trimmomatic.fastq.gz'
trimm_r2_unpaired_file = './data/trimmomatic/'+in_sample+'_merged_R2_'+'unpaired_trimmomatic.fastq.gz'
out_dir = './data/trimmomatic/'
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
# Run command
trimmomatic_cmd = trimmomatic+' '+'PE -threads 4'+' '+r1_file+' '+r2_file+' '+trimm_r1_paired_file+' '+trimm_r1_unpaired_file+' '+trimm_r2_paired_file+' '+trimm_r2_unpaired_file+' '+'ILLUMINACLIP:'+trim_adaptor_file+':2:30:10 CROP:144'
print('\n'+'### Trim off adaptor sequences with trimmomatic and low quality reads')
os.system(trimmomatic_cmd)

### Perform a FastQC analysis after adapter trimming
out_dir = './analysis/fastqc_post_trimmomatic/'
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
# Run command
fastqc_cmd = 'fastqc -o'+' '+out_dir+' '+trimm_r1_paired_file+' '+trimm_r2_paired_file+' '+trimm_r1_unpaired_file+' '+trimm_r2_unpaired_file
print('\n'+'### Perform a FastQC analysis after adapter trimming')
os.system(fastqc_cmd)

### Trim off remaining low quality reads with sickle
sickle_r1_paired_file = './data/sickle/'+in_sample+'_merged_R1_'+'paired_sickle.fastq.gz'
sickle_r2_paired_file = './data/sickle/'+in_sample+'_merged_R2_'+'paired_sickle.fastq.gz'
sickle_singles_pe_file = './data/sickle/'+in_sample+'_merged_'+'singles_pe_sickle.fastq.gz'
sickle_singles_se_file = './data/sickle/'+in_sample+'_merged_'+'singles_se_sickle.fastq.gz'
out_dir = './data/sickle/'
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
# Run command
sickle_pe_cmd = 'sickle pe -f '+trimm_r1_paired_file+' -r '+trimm_r2_paired_file+' -t sanger -o '+sickle_r1_paired_file+' -p '+sickle_r2_paired_file+' -s '+sickle_singles_pe_file+' -q 20 -l 50 -g'
sickle_se_cmd = 'sickle se -f '+trimm_r1_unpaired_file+' -t sanger -o '+sickle_singles_se_file+' -q 30 -l 50 -g'
print('\n'+'### Trim off remaining low quality reads with sickle')
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
print('\n'+'### Perform FastQC post read quality trim via sickle')
os.system(fastqc_cmd)

### Read mapping and alignment with bwa
out_dir = './data/bwa/'
bwa_sam = './data/bwa/'+in_sample+'_bwa_nrg_yet.sam'
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
# Check if reference is indexed
if not os.path.exists(galgal5_ref_index):
	idx_cmd = 'samtools faidx '+galgal5_ref
	os.system(idx_cmd)
# Run command
# Reminder: May want ot use -M command to make compatible with picard tools
bwa_cmd = 'bwa mem -t 4 -T 30 '+galgal5_ref+' '+sickle_r1_paired_file+' '+sickle_r2_paired_file+' > '+bwa_sam
print('\n'+'### Read mapping and alignment with bwa')
os.system(bwa_cmd)

### Remove trimmed reads
os.remove(trimm_r1_paired_file)
os.remove(trimm_r2_paired_file)
os.remove(trimm_r1_unpaired_file)
os.remove(trimm_r2_unpaired_file)

###Add ReadGroups using picard
barcode = sample2ids[in_sample].replace(' (', '_').replace(')', '')
out_dir = './data/post_alignment/'
pic_sam = './data/post_alignment/'+in_sample+'_bwa_readgroups.sam'
rgid = in_sample
rgpl = "ILLUMINA"
rgpu = barcode
rgsm = in_sample
rglb = sample_type
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
# Run the command
pic_cmd = picard_readgroups+' INPUT='+bwa_sam+' OUTPUT='+pic_sam+' RGID='+rgid+' RGPL='+rgpl+' RGPU='+rgpu+' RGSM='+rgsm+' RGLB='+rglb
print('\n'+'###Add readgroups using picard')
os.system(pic_cmd)
# Remove sam file
os.remove(bwa_sam)

### Sort the sam file
rg_bam = './data/post_alignment/'+in_sample+'_bwa_readgroups.bam'
sort_cmd = picard_sortsam+' INPUT='+pic_sam+' OUTPUT='+rg_bam+' SORT_ORDER=coordinate'
print('\n'+'### Sort the sam file')
os.system(sort_cmd)
# Remove read group sam
os.remove(pic_sam)

### Mark duplicates within the bam file
marked_bam = './data/post_alignment/'+in_sample+'_bwa_marked.bam'
metrics = './data/post_alignment/'+in_sample+'_dedup_metrics.list'
dupes_cmd = picard_markdups+' INPUT='+rg_bam+' OUTPUT='+marked_bam+' METRICS_FILE='+metrics+' MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000'
print('\n'+'### Mark duplicates within the bam file')
os.system(dupes_cmd)
# Remove sorted bam
os.remove(rg_bam)

### Index the marked by lane bam file
indexed_bam = './data/post_alignment/'+in_sample+'_bwa_marked.bai'
# The bam index is not automatically rewritten if the file exists. Check if file exists. If so, delete it and create a new one
if os.path.exists(indexed_bam):
	os.remove(indexed_bam)
index_bam_cmd = picard_index+' INPUT='+marked_bam+' OUTPUT='+indexed_bam
print('\n'+'### Index the marked bam file')
os.system(index_bam_cmd)

### Indel realignment
# A fasta dict file of reference is required. Check if exists, if not, then generate it
ref_dict_cmd = picard_dict+' R='+galgal5_ref+' O='+galgal5_ref_dict
if not os.path.exists(galgal5_ref_dict):
	print('\n'+'### Generate dictionary of contigs and length for reference')
	os.system(ref_dict_cmd)
intervals = './data/post_alignment/'+in_sample+'_intervals.list'
realigned_bam = './data/post_alignment/'+in_sample+'_realigned.bam'
realigned_bam_index = './data/post_alignment/'+in_sample+'_realigned.bai'
realign_target_cmd = gatk+' -T RealignerTargetCreator -R '+galgal5_ref+' -I '+marked_bam+' -o '+intervals
realign_cmd = gatk+' -T IndelRealigner -R '+galgal5_ref+' -I '+marked_bam+' -targetIntervals '+intervals+' -o '+realigned_bam
realign_index_cmd = picard_index+' INPUT='+realigned_bam+' OUTPUT='+realigned_bam_index
print('\n'+'### Create realignment targets')
os.system(realign_target_cmd)
print('\n'+'### Indel realignment and index output')
os.system(realign_cmd)
os.system(realign_index_cmd)
# Remove the marked bam and index
os.remove(marked_bam)
os.remove(indexed_bam)

### Perform base quality recalibration with gatk
## Analyze patterns of covariation in the sequence dataset
recal_table = './data/post_alignment/'+in_sample+'_recalibration_data.table'
bqsr_1_cmd = gatk+' -T BaseRecalibrator -R '+galgal5_ref+ ' -I '+realigned_bam+' -knownSites '+dbsnp+' -o '+recal_table
print('\n'+'## Analyze patterns of covariation in the sequence dataset')
os.system(bqsr_1_cmd)

## Perform a second pass to analyze covariation remaining after recalibration
post_recal_table = './data/post_alignment/'+in_sample+'_post_recalibration_data.table'
bqsr_2_cmd = gatk+' -T BaseRecalibrator -R '+galgal5_ref+' -I '+realigned_bam+' -knownSites '+dbsnp+' -BQSR '+recal_table+' -o '+post_recal_table
print('\n'+'## Perform a second pass to analyze covariation remaining after recalibration')
os.system(bqsr_2_cmd)

## Generate before/after plots
out_dir = './analysis/bsqr/'
# Check if output directory exists, if not create it
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
recal_plots = './analysis/bsqr/'+in_sample+'_recalibration_plots.pdf'
bqsr_3_cmd = gatk+' -T AnalyzeCovariates -R '+galgal5_ref+' -before '+recal_table+' -after '+post_recal_table+' -plots '+recal_plots
print('\n'+'## Generate before/after plots')
os.system(bqsr_3_cmd)

# If an error is expressed the gsalib may need to be installed 
# see http://gatkforums.broadinstitute.org/gatk/discussion/1244 to install package locally in R

# Check if output directory exists, if not create it
out_dir = './data/final_bam/'
if not os.path.exists(out_dir):
	os.makedirs(out_dir)

## Apply the recalibration to your sequence data
bqsr_bam = out_dir+in_sample+'_bwa_rg_dedupped_realigned_bqsr.bam'
bqsr_4_cmd = gatk+' -T PrintReads -R '+galgal5_ref+' -I '+realigned_bam+' -BQSR '+recal_table+' -o '+bqsr_bam
print('\n'+'## Apply the recalibration to your sequence data')
os.system(bqsr_4_cmd)

### Index the final bam file
bqsr_bam_index = out_dir+in_sample+'_bwa_rg_dedupped_realigned_bqsr.bai'
if os.path.exists(bqsr_bam_index):
	os.remove(bqsr_bam_index)
# Run command
final_bam_index_cmd = picard_index+' INPUT='+bqsr_bam+' OUTPUT='+bqsr_bam_index
print('\n'+'### Index the final bam file')
os.system(final_bam_index_cmd)

# Remove the dedupped realigned merged bam and index
#os.remove(realigned_bam)
#os.remove(realigned_bam_index)

### Finished Script
print('Fin'+'\n')
