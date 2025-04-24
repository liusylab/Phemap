#add eaf info;filter ADD;OR to BETA
import sys
import os
import math
import gzip
info =gzip.open(sys.argv[1],'r')
hl = info.readline().decode('utf8').strip().split('\t')
snpdict = {}
for line in info:
	line = line.decode('utf8').strip().split('\t')
	key = '%s_%s'%(line[hl.index('chr')],line[hl.index('bp')])
	snpdict[key] = [line[hl.index('EAF')],line[hl.index('INFO_SCORE')]]
info.close()

plinkfile = open(sys.argv[2],'r')
hl = plinkfile.readline().strip().split('\t')
out = open(sys.argv[3],'w')
#out.write('#CHROM\tPOSITION\tA1\tA2\tFRQ\tINFO\tBETA\tSE\tP\n')
for line in plinkfile:
	line = line.strip().split('\t')
	if line[hl.index('TEST')] == 'ADD':
		CHR = 'chr%s'%line[hl.index('#CHROM')]
		if line[hl.index('A1')] == line[hl.index('ALT')]:
			A2 = line[hl.index('REF')]
		else:
			A2 = line[hl.index('ALT')]
		k = '%s_%s'%(CHR,line[hl.index('POS')])
		m = sys.argv[2].split('.')[-1]
		if m == 'logistic' :
			if line[hl.index('OR')] == 'NA' or line[hl.index('OR')] == '0':
				beta = 'NA'
				SE='NA'
			else:
				beta = str(math.log(float(line[hl.index('OR')])))
				SE=line[hl.index('LOG(OR)_SE')]
				#line[hl.index('SE')] = line[hl.index('LOG(OR)_SE')]
		elif m == 'linear' :
			if line[hl.index('BETA')] == 'NA':
				beta = 'NA'
				SE='NA'
			else:
				beta = line[hl.index('BETA')]
				SE=line[hl.index('SE')]
		if k not in snpdict:
			frq = '0'
			info = '0'
		else:
			frq = snpdict[k][0]
			info = snpdict[k][1]
		o = [CHR,line[hl.index('POS')],line[hl.index('ID')],line[hl.index('A1')],A2,frq,info,beta,SE,line[hl.index('P')]]
		out.write('%s\n'%('\t'.join(o)))
plinkfile.close()
out.close()
