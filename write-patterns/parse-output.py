##
# parse-output.py - Created by Timothy Morey on 12/7/2012
#
# This file contains a script that analyzes the output of our batch runs, 
# gathering the relevant pieces of information and saving them to a single
# csv file.
#

import os
import sys

basename = 'wp'
if len(sys.argv) > 1:
	basename = sys.argv[1]
        
infiles = []
for filename in os.listdir('.'):
	if filename.startswith('o.' + basename):
		infiles.append(filename)

outfile = open('o.' + basename + '.csv', 'w')
outfile.write('jobid,N,pattern,stripesize,stripecount,time\n')
jobid = ''
N = ''
time = ''
stripesize = ''
stripecount = ''

for filename in infiles:
	print 'Reading', filename

	s = filename.replace('o.' + basename, '', 1)

	if s.startswith('_N'):
		s = s[2:]
		dotpos = s.find('.')
		scrpos = s.find('_')
		if scrpos < 0 or scrpos > dotpos:
			N = s[0:dotpos]
			s = s[dotpos:]
		else:
			N = s[0:scrpos]
			s = s[scrpos:]

	jobid = s[s.find('.')+1:]

	for line in open(filename):
		if line.startswith('Stripe size: '):
			stripesize = line[13:-2]
		if line.startswith('Stripe count: '):
			stripecount = line[14:-2]
		if line.startswith('Stripe-aligned write took '):
			time = line[26:-3]
			outfile.write(','.join([jobid, N, 'stripe-aligned', stripesize, stripecount, time]) + '\n')
		if line.startswith('Contiguous write took '):
			time = line[22:-3]
			outfile.write(','.join([jobid, N, 'contiguous', stripesize, stripecount, time]) + '\n')

outfile.close()
print 'Done'
