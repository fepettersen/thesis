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

def setup_plot():
	mpl.ion()
	fig  = mpl.figure()
	ax = fig.add_subplot(111,projection='3d')
	return fig,ax

class Experiment:
	"""Usage:
	this_dir = right_split(os.getcwd(),'/')
	obj = Experiment(this_dir,True)
	Hc = [i for i in xrange(10)]
	obj.SetupRun(x0,x1,y0,y1,n,T)
	obj.RunDeterministic()
	for i in Hc:
		obj.RunSimulation(i)
	obj.CalculateError()
	obj.PlotError()
	obj.Finish()
	"""
	def __init__(self,path,visualize,DEBUG=False,save=True):
		self.debug = DEBUG
		self.visualize = visualize
		this_dir = path

		t = time.gmtime()
		self.datetime = '%02d%02d%d_%02d%02d'%(t.tm_mday,t.tm_mon,t.tm_year,t.tm_hour,t.tm_min)
		self.url = 'https://raw.github.com/fepettersen/thesis/master/doc/results'+'/experiment_%s/results'%self.datetime
		self.parent_path = this_dir +'/doc/results/experiment_%s'%self.datetime
		self.code_path = self.parent_path+'/code'
		self.parameter_path = self.parent_path+'/parameters'
		self.result_path = self.parent_path+ '/results'
		self.tofile = 1 if save else 0;
		if not DEBUG:
			os.system('mkdir %s'%self.parent_path)
			os.system('mkdir %s'%self.parameter_path)
			os.system('mkdir %s'%self.result_path)
			os.system('mkdir %s'%self.code_path)
		
		self.walk = []
		self.no_walk = []
		self.T = 11
		self.n = 11
		self.x0 = self.x1 = self.y0 = self.y1 = 0
		self.runcounter =0
		self.samples = 0

	def compile(self):
		os.system('g++ *.cpp -O2 -o main_walk')

	def SetupRun(self,x0,x1,y0,y1,n,T,filename="tmp_results"):
		self.x0 = x0; self.x1 = x1; self.y0 = y0; self.y1 = y1
		self.n = n; self.T = T
		self.filename = filename

	def RunDeterministic(self):
		# self.walk = [0]*self.T
		os.system('./main_walk %d %f %f %f %f %d %d %s %s %f'%(self.tofile,0,0,0,0,self.n,self.T,self.result_path,"Deterministic",0))
		for step in sorted(glob.glob(self.result_path+'/Deterministic*.bin')):
			self.walk.append(np.fromfile(step,sep=" "))
		# print 'RunDeterministic: shape(walk[0]): ',np.shape(self.walk[0])

	def RunSimulation(self,Hc,nsamples=30):
		tofile = self.tofile
		x0 = self.x0; x1 = self.x1; y0 = self.y0; y1 = self.y1
		n = self.n; T = self.T
		filename = self.filename
		result_path = self.result_path
		# Deterministic run first
		tmp = [0]*T
		for i in range(nsamples):
			os.system('./main_walk %d %f %f %f %f %d %d %s %s %f'%(tofile,x0,x1,y0,y1,n,T,result_path,filename,Hc))
			a = 0
			for step in sorted(glob.glob(result_path+'/tmp_results*.bin')):
				tmp[a] += np.fromfile(step,sep=" ")
				a += 1
		for i in xrange(T):
			tmp[i] /= nsamples

		self.no_walk.append(tmp)
		self.runcounter += 1
		self.samples = nsamples

	def CalculateError(self):
		self.error = [np.zeros(self.T)]*self.runcounter
		for i in xrange(self.runcounter):
			for j in xrange(self.T):
				self.error[i][j] = np.max(np.abs(self.walk[j]-self.no_walk[i][j]))
		# print self.walk

	def SaveError(self):
		fname = result_path+'/error.txt'
		f = np.asanyarray(self.error)
		np.savetxt(fname,f)

	def PlotError(self,save=True):
		# fig = mpl.figure()
		# ax = fig.add_subplot(111)
		# mpl.hold('on')
		# for i in range(self.runcounter):
		# 	ax.plot(self.error[i])
		# mpl.show()
		mpl.hold('on')
		for i in range(self.runcounter):
			mpl.plot(self.error[i])
		print np.where(self.error[0]!=self.error[1])
		if save:
			mpl.savefig(self.result_path+'/errorplot.png')
			mpl.savefig(self.result_path+'/errorplot.eps')
		mpl.show()

	def UpdateWebpageDefault(self):
		datetime = self.datetime
		f = open(this_dir +'/doc/web/index.html','r')
		html = f.read()
		f.close()
		pattern = r'</table>' 	# There is an empty table at the top of the html-file
		match = re.search(pattern,html)
		part = match.span()[-1]+1

		html_code = """<h3><a name="%s" class="anchor" 
		href="#%s"><span class="octicon octicon-link">
		</span></a>New experiment %s.</h3><br>
		Ran new experiment using %d different values of 
		Hc, averaging over %d samples. Results can be found 
		in doc/results/experiment_%s.
		"""%(datetime,datetime,time.ctime(),self.runcounter,self.samples,self.datetime)
		
		sum_html = html[:part]+html_code+html[part:]
		f = open(this_dir +'/doc/web/index.html','w')
		f.write(sum_html)
		f.close()

	def UpdateSpecial(self):
		datetime = self.datetime
		f = open(this_dir +'/doc/web/index.html','r')
		html = f.read()
		f.close()
		pattern = r'</table>' 	# There is an empty table at the top of the html-file
		match = re.search(pattern,html)
		part = match.span()[-1]+1
		
		title = raw_input('Title on webpage (html is taken care of)>> ')
		html_code = """<h3><a name="%s" class="anchor" 
		href="#%s"><span class="octicon octicon-link">
		</span></a>%s</h3>
		"""%(datetime,datetime,title)
		explanaiton = raw_input('Add description (optional):  ')
		html_code += '<br>'+explanaiton+'\n'

		sum_html = html[:part]+html_code+html[part:]
		f = open(this_dir +'/doc/web/index.html','w')
		f.write(sum_html)
		f.close()

	def Finish(self):
		if self.debug:
			os.system('rm -rf %s'%self.parent_path)

DEBUG = False
plot = True
save_files = True
add_text_to_web = False
mode = 'test'


this_dir = right_split(os.getcwd(),'/')

x0 = y0 = 0.6
x1 = y1 = 0.7
n = 21
T = 151
nsamples = 30
dx = 1.0/(n-1)
dt = dx**2/5.0
Hc = [5/dt,200/dt]

run = Experiment(this_dir,plot,DEBUG,save_files)
run.compile()
run.SetupRun(x0,x1,y0,y1,n,T)
run.RunDeterministic()
for i in Hc:
	print "Hc = %g"%i
	run.RunSimulation(i,nsamples)

run.CalculateError()
run.PlotError()
# run.UpdateSpecial()

run.Finish()

