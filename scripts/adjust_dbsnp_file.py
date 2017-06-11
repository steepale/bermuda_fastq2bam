import sys
import os
import re

infile = sys.argv[1]

outfile = open(sys.argv[2], 'w')

for line in open(infile):
	if line[0] == '#':
		outfile.write(line)
	if line[0] != '#':
		line = line.rstrip()
		col = line.split('\t')
		ref = col[3]
		alt = col[4]
		info = col[7]
		if re.search('-', ref):
			ref = 'N'
		if re.search('-', alt):
			alt = 'N'
		if re.search(' ', info):
			info = info.replace(' ', '_') 
		outfile.write(col[0]+'\t'+col[1]+'\t'+col[2]+'\t'+ref+'\t'+alt+'\t'+col[5]+'\t'+col[6]+'\t'+info+'\n')
outfile.close()
