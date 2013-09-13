# -*- coding: utf-8 
"""
 script that runs an experiment, creates a new subdirectory with a copy of 
 the relevant code and parameters, stores the results (plots) and updates 
 the thesis homepage on githubwith this information
"""
DEBUG = True
import os, sys, time, re #,argparse??
t = time.gmtime()
def right_split(s,delimiter):
	i = len(s)-1
	while i>0:
		if s[i] == delimiter:
			return s[:i]
		i -= 1

this_dir = right_split(os.getcwd(),'/')

datetime = '%02d%02d%d_%d%d'%(t.tm_mday,t.tm_mon,t.tm_year,t.tm_hour,t.tm_min)

html_code = """<h3><a name="%s" class="anchor" 
href="#%s"><span class="octicon octicon-link">
</span></a>New experiment %s.</h3>
"""%(datetime,datetime,time.ctime())

parent_path = this_dir +'/doc/results/experiment_%s'%datetime
code_path = parent_path+'/code'
parameter_path = parent_path+'/parameters'
result_path = parent_path+ '/results'


commandline_arguments = []

print this_dir
f = open(this_dir +'/doc/web/index.html','r')
html = f.read()
f.close()
pattern = r'</table>'
match = re.search(pattern,html)
part = match.span()[-1]+1

sum_html = html[:part]+html_code+html[part:]
f = open(this_dir +'/doc/web/index.html','w')
f.write(sum_html)
f.close()

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
