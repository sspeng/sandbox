##
# gather.py - Created by Timothy Morey on 10/17/2012
#
# This file contains a script that analyzes the output of many ranger batch
# runs.  It is currently tuned for doing the n,w parameter sweeps.
#

import os
import sys

basename = ''
if len(sys.argv) > 1:
	basename = sys.argv[1]

infiles = []
for filename in os.listdir('.'):
	if filename.startswith('o.' + basename):
		infiles.append(filename)

outfile = open('o.' + basename + '.csv', 'w')
outfile.write('n,w,i,format,init_time,init%,step_time,step%,out_time,out%,t_real,t_user,t_sys\n')
n = '?'
w = '?'
format = '?'
t_real = '?'
t_user = '?'
t_sys = '?'
i = '?'
init_time = '?'
init_pct = '?'
step_time = '?'
step_pct = '?'
out_time = '?'
out_pct = '?'
ignore_run = False

for filename in infiles:
	print 'Reading', filename

	s = filename.replace('o.' + basename, '', 1)
	if s.startswith('_n'):
		s = s[2:]
		n = s[0:s.find('.')]
		s = s[s.find('_'):]
		print 'Found n:', n

	for line in open(filename):
		if line.startswith('=== I = '):
			i = line[8:].strip()
			print 'Found i:', i

		if line.startswith('=== PISM_OFORMAT = '):
			format = line[18:].strip()
			print 'Found format:', format

		if line.startswith('#(spinup.sh)  running pclimate'):
			ignore_run = True

		if line.startswith(' 1:  Initialization:'):
			init_time = line[21:31].strip()
			init_pct = line[31:37].strip()
			print 'Found init_time:', init_time
			print 'Found init_pct:', init_pct

		if line.startswith(' 2:        Stepping:'):
			step_time = line[21:31].strip()
			step_pct = line[31:37].strip()
			print 'Found step_time:', step_time
			print 'Found step_pct:', step_pct

		if line.startswith(' 3:   PrimaryOutput:'):
			out_time = line[21:31].strip()
			out_pct = line[31:37].strip()
			print 'Found out_time:', out_time
			print 'Found out_pct:', out_pct

		if line.startswith('real '):
			t_real = line[5:].strip()
			print 'Found t_real:', t_real

		if line.startswith('user '):
			t_user = line[5:].strip()
			print 'Found t_user:', t_user

		if line.startswith('sys '):
			t_sys = line[4:].strip()
			print 'Found t_sys:', t_sys
			
			if ignore_run:
				ignore_run = False
			else:	
				outfile.write(','.join([n, w, i, format, 
				                        init_time, init_pct, 
				                        step_time, step_pct, 
				                        out_time, out_pct, 
				                        t_real, t_user, t_sys]))
				outfile.write('\n')

outfile.close()
print 'Done'
