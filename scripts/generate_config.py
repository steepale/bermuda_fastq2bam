import sys
import os

# reference files
lanes_file = './data/D.Wright_16_01_lanes_info.txt'
lib_file = './data/D.Wright_16_01_library_info.txt'
sample_file = './data/D.Wright_16_01_sample_info.txt'
int_config_file = './data/int_config.txt'

# outfile
outfile = open('./data/config.txt', 'w')

# Write header for output file
outfile.write("#ngi_name\tdate\tunknown-id-1\tunknown-id-2\tflowcell\tbcl-conversion-id\tlane\tread\tvolume\tbarcode\n")

# Create dictionary of sample id's to appropriate annotations
sample_anns2barcode = {}

for line in open(lib_file):
	if not line.startswith('NGI'):
		line = line.rstrip()
		col = line.split('\t')
		ngi_name = col[0]
		barcode = col[1]
		status = col[4]
		if status == 'PASSED':
			sample_anns2barcode[ngi_name] = barcode

# Iterate through reference files and sample strings to fill in config file
for line in open(int_config_file):
	if not line.startswith('#'):
		line = line.rstrip()
		col = line.split('\t')
		# Apply variables
		ngi_name = col[0]
		date = col[1]
		unknown_1 = col[2]
		unknown_2 = col[3]
		flowcell = col[4]
		bcl_conversion = col[5]
		lane = col[6]
		read = col[7]
		volume = col[8]
		if ngi_name in sample_anns2barcode.keys():
			barcode = sample_anns2barcode[ngi_name]
		else:
			print('WARNING')
			print(ngi_name+'\t'+date+'\t'+unknown_1+'\t'+unknown_2+'\t'+flowcell+'\t'+bcl_conversion+'\t'+lane+'\t'+read+'\t'+volume)
		# Print appropriate values to outfile
		outfile.write(ngi_name+'\t'+date+'\t'+unknown_1+'\t'+unknown_2+'\t'+flowcell+'\t'+bcl_conversion+'\t'+lane+'\t'+read+'\t'+volume+'\t'+barcode+'\n')
# Close outfile to protect EOF
outfile.close()
