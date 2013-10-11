# -*- coding: utf-8 
"""
 script that runs an experiment, creates a new subdirectory with a copy of 
 the relevant code and parameters, stores the results (plots) and updates 
 the thesis homepage on github with this information
"""
import os, sys, time, re, numpy as np, matplotlib.pyplot as mpl #,argparse??
import glob
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
gitpush = False
add_text_to_web = False
mode = 'test'


########################
# - Make directories - #
########################

this_dir = right_split(os.getcwd(),'/')
t = time.gmtime()
datetime = '%02d%02d%d_%02d%02d'%(t.tm_mday,t.tm_mon,t.tm_year,t.tm_hour,t.tm_min)
url = 'https://raw.github.com/fepettersen/thesis/master/doc/results'+'/experiment_%s/results'%datetime
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

tofile = 1 if save else 0;
x0 = 0.3
x1 = 0.5
y0 = 0.3
y1 = 0.5
n = 21
T = 50
filename = "Using_walk"

os.system('g++ *.cpp -O2 main_walk') 	#compile for good measure
os.system('./main_walk %d %f %f %f %f %d %d %s %s'%(tofile,x0,x1,y0,y1,n,T,result_path,filename))
filename = "Without_walk"
x0 = x1 = y0 = y1 = 0
# os.system('./main_walk %d %f %f %f %f %d %d %s %s'%(tofile,x0,x1,y0,y1,n,T,result_path,filename))
###########################
# -- Visualize results -- #
###########################

images = []
if plot:
	"2D"
	counter = 0
	def setup_plot():
		mpl.ion()
		fig  = mpl.figure()
		ax = fig.add_subplot(111,projection='3d')
		ax.set_autoscaley_on(False)
		return fig,ax

	fig,ax = setup_plot()
	X,Y = np.meshgrid(np.linspace(0,1,n),np.linspace(0,1,n))
	for step in sorted(glob.glob(result_path+'/*.txt')):
		img = np.loadtxt(step)
		wframe = ax.plot_wireframe(X,Y,img)
		mpl.draw()
		time.sleep(0.2)
		if counter==0 or counter==int(T/2) or counter==T-1:
			if not DEBUG: 
				mpl.savefig(result_path+'/from_simulation%s_%d.eps'%(datetime,counter))
				mpl.savefig(result_path+'/from_simulation%s_%d.png'%(datetime,counter))
			images.append('/from_simulation%s_%d.png'%(datetime,counter))
		# time.sleep(1)
		ax.collections.remove(wframe)
		counter+=1
# no_walk = []
# for step in sorted(glob.glob(result_path+'/W*.txt')):
# 	no_walk.append(np.loadtxt(step))
# walk = []
# for step in sorted(glob.glob(result_path+'/U*.txt')):
# 	walk.append(np.loadtxt(step))
# error = []
# for i in range(len(walk)):
# 	error.append(np.max(np.abs(no_walk[i]-walk[i])))
########################
# -- update website -- #
########################

# This will be what is added to the webpage after a new run
# html_code = """<h3><a name="%s" class="anchor" 
# href="#%s"><span class="octicon octicon-link">
# </span></a>New experiment %s.</h3>
# <img src="%s" height="42" width="42">
# """%(datetime,datetime,time.ctime(),url +images[0])

html_code = """<h3><a name="%s" class="anchor" 
href="#%s"><span class="octicon octicon-link">
</span></a>New experiment %s.</h3>

"""%(datetime,datetime,time.ctime())

if add_text_to_web:
	explanaiton = raw_input('Add description (optional):  ')
	html_code += '<br>'+explanaiton

if not DEBUG:
	"Write the values of all parameters to a .txt file"
	os.system('cp *.cpp *.h new_experiment.py %s'%code_path)
	f = open(parameter_path+'/parameters.txt','w')
	# f.write("This is the comparison between using and not unsing walk. Should be placed somewhere else! \n")
	# f.write("[")
	# for line in error:
	# 	to_write = str(line)+', '
	# 	f.write(to_write)
	# 	f.write("]")
	f.close()


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

if mode =='test':
	os.system('rm -rf %s'%parent_path)
	print 'Done testing, all files removed'
else:
	print 'Done, results in folder: \n'
	print result_path
if gitpush:
	"Add, commit and push the results to github -- Doesnt work because of directory..."
	os.system('cd %s'%this_dir)
	os.system('git checkout master') 	#Force branch master!
	os.system('git add .')
	os.system('git commit -am "Ran new experiment"')
	os.system('git checkout gh-pages')
	os.system('git merge master')
	os.system('git push')
	os.system('git checkout master')
	# os.system('cd code')
