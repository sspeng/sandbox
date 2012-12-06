##
# gather.py - Created by Timothy Morey on 10/17/2012
#
# This file contains a script that analyzes the output of our batch runs, 
# gathering the relevant pieces of information and saving them to a single
# csv file.
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
outfile.write('jobid,n,N,w,format,init_time,init%,step_time,step%,out_time,out%,t_real,t_user,t_sys\n')
jobid = ''
n = ''
N = ''
w = ''
format = ''
t_real = ''
t_user = ''
t_sys = ''
init_time = ''
init_pct = ''
step_time = ''
step_pct = ''
out_time = ''
out_pct = ''
ignore_run = False

for filename in infiles:
	print 'Reading', filename

	s = filename.replace('o.' + basename, '', 1)
	if s.startswith('_n'):
		s = s[2:]
                dotpos = s.find('.')
                scrpos = s.find('_')
                if scrpos < 0 or scrpos > dotpos:
                        n = s[0:dotpos]
                        s = s[dotpos:]
                else:
                        n = s[0:scrpos]
                        s = s[scrpos:]
		print 'Found n:', n

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
		print 'Found N:', N

	if s.startswith('_w'):
		s = s[2:]
		w = s[0:s.find('.')]
		print 'Found w:', w

	jobid = s[s.find('.')+1:]

	for line in open(filename):
		if line.lower().startswith('=== i = '):
			i = line[8:].strip()
			print 'Found i:', i

                elif line.lower().startswith('=== n = '):
                        n = line[8:].strip()
                        print 'Found n:', n

		elif line.startswith('=== PISM_OFORMAT = '):
			format = line[18:].strip()
			print 'Found format:', format

                elif line.startswith('-o_format '):
                        format = line[10:].strip()
                        print 'Found format:', format

		elif line.startswith('#(spinup.sh)  running pclimate'):
			ignore_run = True

		elif line.startswith(' 1:  Initialization:'):
			init_time = line[21:31].strip()
			init_pct = line[31:37].strip()
			print 'Found init_time:', init_time
			print 'Found init_pct:', init_pct

		elif line.startswith(' 2:        Stepping:'):
			step_time = line[21:31].strip()
			step_pct = line[31:37].strip()
			print 'Found step_time:', step_time
			print 'Found step_pct:', step_pct

		elif line.startswith(' 3:   PrimaryOutput:'):
			out_time = line[21:31].strip()
			out_pct = line[31:37].strip()
			print 'Found out_time:', out_time
			print 'Found out_pct:', out_pct

		elif line.startswith('real '):
			t_real = line[5:].strip()
			print 'Found t_real:', t_real

		elif line.startswith('user '):
			t_user = line[5:].strip()
			print 'Found t_user:', t_user

		elif line.startswith('sys '):
			t_sys = line[4:].strip()
			print 'Found t_sys:', t_sys
			
			if ignore_run:
				ignore_run = False
			else:	
				outfile.write(','.join([jobid, n, N, w, format, 
				                        init_time, init_pct, 
				                        step_time, step_pct, 
				                        out_time, out_pct, 
				                        t_real, t_user, t_sys]))
				outfile.write('\n')

outfile.close()
print 'Done'
