# -*- coding: utf-8 
"""
 script that runs an experiment, creates a new subdirectory with a copy of 
 the relevant code and parameters, stores the results (plots) and updates 
 the thesis homepage on githubwith this information
"""
DEBUG = True

import os, sys, time, re #,argparse??
t = time.gmtime()

datetime = '%02d%02d%d_%d%d'%(t.tm_mday,t.tm_mon,t.tm_year,t.tm_hour,t.tm_min)
parent_path = '~/uio/thesis/doc/results/experiment_%s'%datetime
code_path = parent_path+'/code'
parameter_path = parent_path+'/parameters'
result_path = parent_path+ '/results'

commandline_arguments = []

with open('uio/thesis/doc/web/index.html') as f:
	html = f.read()
pattern = r'</table>'
match = re.search(pattern,html)
print match

if not DEBUG:
	os.system('mkdir %s'%parent_path)
	os.system('mkdir %s'%parameter_path)
	os.system('mkdir %s'%result_path)
	os.system('mkdir %s'%code_path)

# os.makedir(parent_path)
# os.makedir(code_path)
# os.makedir(parameter_path)
# os.makedir(result_path)

# os.system('python run.py') #save results in result_path! they should be ignored by git
