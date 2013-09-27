# -*- coding: utf-8 
"""
 script that runs an experiment, creates a new subdirectory with a copy of 
 the relevant code and parameters, stores the results (plots) and updates 
 the thesis homepage on githubwith this information
"""
import os, sys, time, re, numpy as np, matplotlib.pyplot as mpl #,argparse??
from combine import MultiscaleSolver

def right_split(s,delimiter):
	i = len(s)-1
	while i>0:
		if s[i] == delimiter:
			return s[:i]
		i -= 1

DEBUG = False
plot = True
save = True
gitpush = True

########################
# - Make directories - #
########################

this_dir = right_split(os.getcwd(),'/')
t = time.gmtime()
datetime = '%02d%02d%d_%d%d'%(t.tm_mday,t.tm_mon,t.tm_year,t.tm_hour,t.tm_min)

# This will be what is added to the webpage after a new run
html_code = """<h3><a name="%s" class="anchor" 
href="#%s"><span class="octicon octicon-link">
</span></a>New experiment %s.</h3>
"""%(datetime,datetime,time.ctime())

parent_path = this_dir +'/doc/results/experiment_%s'%datetime
code_path = parent_path+'/code'
parameter_path = parent_path+'/parameters'
result_path = parent_path+ '/results'

if not DEBUG:
	os.system('mkdir %s'%parent_path)
	os.system('mkdir %s'%parameter_path)
	os.system('mkdir %s'%result_path)
	os.system('mkdir %s'%code_path)

########################
# -- Run simulation -- #
########################

t = 0
T = 20
def setup_plot():
	mpl.ion()
	fig  = mpl.figure()
	ax = fig.add_subplot(111,projection='3d')
	ax.set_autoscaley_on(False)
	return fig,ax

"2D"
nx = 11; ny =11
Up = np.zeros((nx,ny))
Up[(nx)/2:,(nx)/2:] = 1
area = [[0.3,0.3],[0.5,0.5]]
if plot:
	fig,ax = setup_plot()
	X,Y = np.meshgrid(np.linspace(0,1,nx),np.linspace(0,1,ny))
mesh = [np.linspace(0,1,nx),np.linspace(0,1,ny)]
test = MultiscaleSolver(mesh)
test.AddWalkArea(area)
test.setInitialCondition(Up)
wframe = ax.plot_wireframe(X,Y,test.Up)
if plot: mpl.draw()
while t<T:
	if plot: ax.collections.remove(wframe)
	test.Solve()
	if save: test.SaveState(path=result_path)
	if plot:
		wframe = ax.plot_wireframe(X,Y,test.Up)
		mpl.draw()
	# time.sleep(1)
	t+=1

########################
# -- update website -- #
########################

commandline_arguments = []

# print this_dir
f = open(this_dir +'/doc/web/index.html','r')
html = f.read()
f.close()
pattern = r'</table>' 	# There is an empty table at the top of the html-file
match = re.search(pattern,html)
part = match.span()[-1]+1

sum_html = html[:part]+html_code+html[part:]

if not DEBUG:
	f = open(this_dir +'/doc/web/index.html','w')
	f.write(sum_html)
	f.close()


if gitpush:
	os.system('git add .')
	os.system('git commit -am "Ran new experiment"')
	os.system('git checkout gh-pages')
	os.system('git merge master')
	os.system('git push')
	os.system('git checkout master')
# os.makedir(parent_path)
# os.makedir(code_path)
# os.makedir(parameter_path)
# os.makedir(result_path)

# os.system('python run.py') #save results in result_path! they should be ignored by git
