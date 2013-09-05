# -*- coding: utf-8 
# Module to combine (diffusion) PDE solver and random walk module

class MultiscaleSolver:
	"""docstring for MultiscaleSolver
	Combine a normal diffusion PDE solver with random walk model for 
	diffusion in certain areas"""
	
	def __init__(self, arg):
		# super(MultiscaleSolver, self).__init__()
		self.arg = arg
		