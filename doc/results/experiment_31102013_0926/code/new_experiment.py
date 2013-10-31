# -*- coding: utf-8 
"""
 script that runs an experiment, creates a new subdirectory with a copy of 
 the relevant code and parameters, stores the results (plots) and updates 
 the thesis homepage on github with this information
"""
import os, sys, time, re, numpy as np, matplotlib.pyplot as mpl #,argparse??
import glob
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.axes3d import Axes3D

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
	def __init__(self,path, DEBUG=False, save=True):
		self.debug = DEBUG
		this_dir = path

		t = time.gmtime()
		self.datetime = '%02d%02d%d_%02d%02d'%(t.tm_mday,t.tm_mon,t.tm_year,t.tm_hour,t.tm_min)
		self.url = 'https://raw.github.com/fepettersen/thesis/master/doc/results'+'/experiment_%s/results'%self.datetime
		self.parent_path = this_dir +'/doc/results/experiment_%s'%self.datetime
		self.code_path = self.parent_path+'/code'
		self.parameter_path = self.parent_path+'/parameters'
		self.result_path = self.parent_path+ '/results'
		self.tofile = 1 if save else 0;

		os.system('mkdir %s'%self.parent_path)
		os.system('mkdir %s'%self.parameter_path)
		os.system('mkdir %s'%self.result_path)
		os.system('mkdir %s'%self.code_path)
		os.system('cp *.cpp *.h new_experiment.py %s'%self.code_path)
		
		self.walk = []
		self.no_walk = []
		self.T = 11
		self.n = 11
		self.x0 = self.x1 = self.y0 = self.y1 = 0
		self.runcounter =0
		self.samples = 0

	def compile(self):
		os.system('g++ *.cpp -O2 -o main_walk')

	def SetupRun(self,x0,x1,y0,y1,m,n,T,dt,filename="a"):
		self.x0 = x0; self.x1 = x1; self.y0 = y0; self.y1 = y1
		self.m = m; self.n = n; self.T = T; self.dt = dt
		self.filename = filename

	def RunDeterministic(self):
		# self.walk = [0]*self.T
		os.system('./main_walk %d %f %f %f %f %d %d %d %s %s %f %f'%(self.tofile,0,0,0,0,self.m,self.n,self.T,self.result_path,"Deterministic",0,self.dt))
		for step in sorted(glob.glob(self.result_path+'/Deterministic*.txt')):
			self.walk.append(np.fromfile(step,sep=" "))

	def RunSimulation(self,Hc):
		tofile = self.tofile
		x0 = self.x0; x1 = self.x1; y0 = self.y0; y1 = self.y1
		m = self.m; n = self.n; T = self.T
		filename = self.filename
		result_path = self.result_path
		tmp = [0]*T
		os.system('./main_walk %d %f %f %f %f %d %d %d %s %s %f %f'%(tofile,x0,x1,y0,y1,m,n,T,result_path,filename,Hc,self.dt))
		self.runcounter += 1

	def CalculateError(self,Hc,exact=None):
		self.legends = Hc
		self.error = []
		tmp = np.zeros(self.T)
		if exact is None:
			for i in xrange(self.runcounter):
				for j in xrange(self.T):
					self.error[i][j] = np.max(np.absolute(self.walk[j]-self.no_walk[i][j]))
		else:
			if self.n <= 1:
				X = np.linspace(0,1,self.m)
				Y = np.zeros(self.m)
			else:
				X,Y = np.meshgrid(np.linspace(0,1,self.m),np.linspace(0,1,self.n))
			for a in range(len(Hc)):
				i = Hc[a]
				for j in xrange(self.T): 
					step = np.loadtxt(self.result_path+'/results_FE_Hc%d_n%03d.txt'%(i+1,j))
					tmp[j] = np.max(np.linalg.norm(self.exact(X,Y,(j+1)*self.dt)-step))
				self.error.append(tmp.copy())



	def SaveError(self,header=None):
		fname = self.result_path+'/error.txt'
		f = np.asanyarray(self.error)
		if header is None:
			np.savetxt(fname,f)
		else:
			np.savetxt(fname,f)

	def PlotError(self,save=True):
		# mpl.hold('on')
		color = ['b-o','r-o','k-o','c-o','g-o','m-o','b-x','r-x','k-x','c-x','g-x','m-x']
		for i in range(len(self.error)):
			# mpl.plot(np.log(self.error[i]/self.dt))
			mpl.plot(self.error[i],color[i],label='Hc = %d'%(self.legends[i]+1))
			# print self.error[i]
		mpl.legend(loc=4)
		mpl.xlabel('timestep no.')
		mpl.ylabel('max(abs(simulation-exact))')
		mpl.title('Error plot')
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
		Hc. Results can be found 
		in doc/results/experiment_%s.
		"""%(datetime,datetime,time.ctime(),self.runcounter,self.datetime)
		
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


	def Visualize(self,path=None,filename=None):
		if path is None:
			path = self.result_path
		if filename is None:
			filename = '/results'
		im = []
		if self.n<=1:
			fig = mpl.figure()
			x = np.linspace(0,1,self.m)
			counter = 1
			for step in sorted(glob.glob(path+filename+'*.txt')):
				# f = open(step,'r')
				# tmp = np.asarray(f.readlines())
				# tmp = tmp.astype(np.float32)
				tmp = np.loadtxt(step)
				im.append(mpl.plot(x,(self.exact(x,np.zeros(self.m),counter*self.dt)-tmp),'b-'))
				# im.append(mpl.plot(x,tmp,'b-'))
				# print np.max(np.absolute(self.exact(x,np.zeros(self.m),counter*self.dt)-tmp))
				counter += 1
			ani = animation.ArtistAnimation(fig,im)
			mpl.show()
		else:
			X,Y = np.meshgrid(np.linspace(0,1,self.m),np.linspace(0,1,self.n))
			mpl.ion()
			fig = mpl.figure()
			ax = fig.add_subplot(111,projection='3d')
			counter = 1
			for step in sorted(glob.glob(path+filename+'*.txt')):
				tmp = np.loadtxt(step)
				wframe = ax.plot_wireframe(X,Y,(self.exact(X,Y,(counter*self.dt))-tmp))
				mpl.draw()
				if counter==1:
					pass
					# ax.set_autoscaley_on(False)
				ax.collections.remove(wframe)
				counter +=1
			# ax.plot_wireframe(X,Y,self.exact(X,Y,0))
			# mpl.draw()
			# time.sleep(1)


	def exact(self,x,y,t):
		return 0


	def Finish(self):
		if self.debug:
			os.system('rm -rf %s'%self.parent_path)



def f(x,y,t):
	return np.exp(-t*np.pi**2)*np.cos(np.pi*x) +1
		# return np.ones(np.shape(x))*1.5

if __name__ == '__main__':
	DEBUG = False
	save_files = True
	mode = 'test'


	this_dir = right_split(os.getcwd(),'/')

	x0 = 0.6
	y0 = 0.6
	x1 = 0.7
	y1 = 0.7
	m = 21
	n = 1
	T = 51	
	dx = 1.0/(m-1)
	dy = 1.0/(n-1) if n>1 else 0
	dt = dx*dy/5.0 if n>1 else dx**2/3.0
	Hc = [5/dt]
	name = '/home/fredriep/Dropbox/uio/thesis/doc/results/experiment_18102013_1337/results/'
	print dt

	run = Experiment(this_dir,DEBUG,save_files)
	run.exact = f
	run.compile()
	run.SetupRun(x0,x1,y0,y1,m,n,T,dt)
	# run.RunDeterministic()
	for i in Hc:
		print "Hc = %g"%i
		run.RunSimulation(i)
	time.sleep(1)
	run.CalculateError(Hc,exact=True)

	run.PlotError()
	run.SaveError(header="max(abs(error)) for manifactured solution u(x,t) = exp(-t*pi**2*cos(pi*x) in 1D. Hc = %g"%Hc[0])
	# run.UpdateSpecial()
	run.Visualize()
	run.Finish()

	# (self.tofile,0,0,0,0,self.n,self.T,self.result_path,"Deterministic",0))