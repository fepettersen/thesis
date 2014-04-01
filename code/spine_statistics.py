"""Script that reads the spine info file and gathers statistics in a separate file
structure of statistics file is:
position in um		seconds from av. to conc.increase		neck length		neck width		(total time 	no. of particles)#repeat
"""

import re


class Parser:
	"""docstring for Parser"""
	def __init__(self, infile):
		self.infile = open(infile,'r')
		self.timestamp = 0
		self.dx = self.dt = 0
		self.information = {}

	def read(self):
		for line in self.infile.readlines():
			tmp = re.search(r't\s?=\s?(\d+)\s+(\d+)',line)
			if tmp:
				self.timestamp = int(tmp.groups()[0])
				if self.timestamp == 1:
					# read next int(tmp.groups()[1]) + 1 lines
					pass
			else:
				tmp = re.search(r'^dx\s?=\s?(\d+\.?\d+)\s*dt\s?=\s?(\d+\.?\d+)',line)
				if tmp:
					dx,dt = tmp.groups()
					self.dx = float(dx); self.dt = float(dt)
				elif self.timestamp==1:
					pass
				else:
					self.parse(line)
		print 'done'
		self.infile.close()

	def parse(self,line):
		info = re.search(r'^spine\s\w+\s?(\d+)',line)
		if info:
			pos = int(info.groups()[0])
			self.information[pos] = [pos*self.dx,self.timestamp*self.dt]
		else:
			spline = line.split()
			# print spline
			pos = int(spline[-2])
			try:
				a = self.information[pos]
			except KeyError:
				self.information[pos] = [pos*self.dx,(self.timestamp-1)*self.dt]
			if len(self.information[pos])>3:
				# information has already been written, just update
				if self.information[pos][-1] != int(spline[0]):
					self.information[pos].append(self.timestamp*self.dt)
					self.information[pos].append(int(spline[0]))
			else:
				self.information[pos][1] = self.timestamp*self.dt - self.information[pos][1]
				self.information[pos].append(float(spline[1])) 				# neck length
				self.information[pos].append(float(spline[-1])*self.dx) 	# neck width
				self.information[pos].append(self.timestamp*self.dt)		# total time
				self.information[pos].append(int(spline[0]))				# number of particles

	def update_statisticsfile(self):
		debug = False
		statfile = open('../doc/results/statistics/spine_stats.txt','a')
		for key in self.information:
			string = ''
			for bla in self.information[key]:
				string += (str(bla)+'  ')
			if not debug:
				statfile.write(string+'\n')
		

balle = Parser('spine_info.txt')
balle.read()
balle.update_statisticsfile()
print balle.information
