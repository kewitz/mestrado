# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz

Created on Thu May 15 11:50:45 2014
"""
import numpy as np
from multiprocessing import Pool, Array

# Macros
indexed = lambda l, offset=0: zip(np.arange(len(l))+offset,l)

# Constantes
f = np.float64  # Default format.
pi = f(np.pi)  # Double precision pi.
e = np.exp(f(1))  # Exponential as unit.

class Yee:
	""" Classe de cálculo FDTD 1D. """
	def __init__(self, verbose = True):
		self.vp = f(299792458.0)
		self.eps = f(8.854187817E-12)
		self.mu = f(pi * 4E-7)
		self.sigma = f(5E-15)
		self.verbose = verbose

	def setFreq(self, fop):
		"""
		Ajusta os parâmetros da simulação de acordo com a frequência (Hz) de
		operação desejada.
		"""
		# Define constantes:
		self.fop = fop
		self.lamb = self.vp/fop
		self.t0 = 1.0/fop
		self.tal = self.t0/(2*np.sqrt(2*np.log(2)))
		self.dx = self.lamb/10.0
		self.dy = self.lamb/10.0
		self.dtal = self.dx/2.0
		self.dt = self.dtal/self.vp
		
		# Domínios
		self.xd = np.arange(0,self.lamb*10,self.dx,dtype=np.float64)
		self.yd = np.arange(0,self.lamb*10,self.dy,dtype=np.float64)
		
		self.bound = {
				'Ez': np.ones((len(self.xd),len(self.yd))),
				'Hx': np.ones((len(self.xd),len(self.yd))),
				'Hy': np.ones((len(self.xd),len(self.yd)))
			}	
		
	def setLambda(self, lamb):
		"""
		Ajusta os parâmetros da simulação de acordo com o comprimento de on-
		da (m) de operação desejado.
		"""
		self.setFreq(self.vp/lamb)
		
	def makeDomains(self, iteractions, skip):
		self.td = np.arange(0,self.dt*iteractions,self.dt*skip,dtype=np.float64)
		self.Ez = np.zeros((len(self.td),len(self.xd),len(self.yd)))
		self.Hx = np.zeros((len(self.td),len(self.xd),len(self.yd)))
		self.Hy = np.zeros((len(self.td),len(self.xd),len(self.yd)))

	def run(self, fx, t=500, skip=1):
		""" Executa o FDTD para os domínios configurados. """
		self.makeDomains(t, skip)
		# Constantes
		z = np.sqrt(self.mu/self.eps)
		CEy = z*self.dtal/self.dy
		CEx = z*self.dtal/self.dx
		CH = (1.0/z)*self.dtal/self.dx
		
		st, sx, sy = (len(self.td),len(self.xd),len(self.yd))
		
		pool = Pool(processes=4)
		
		# Iterações no domínio do tempo
		for k, t in indexed(self.td[:-1]):
			# Excitação
			fx(k,t,self)

			
			space = np.meshgrid(range(sx-1), range(sy-1))
			space = zip(space[0].flatten(), space[1].flatten())
			threads = len(space)
			
			# Propagação de H no espaço
			params = zip(space,[(k,self.Ez,self.Hx,CH)]*threads)
			Hy = np.array(pool.map(makeHy, params)).reshape((sx-1,sy-1))
			Hx = np.array(pool.map(makeHx, params)).reshape((sx-1,sy-1))
			self.Hx[k,0:-1,0:-1] = Hx*self.bound['Hx'][0:-1,0:-1]
			self.Hy[k,0:-1,0:-1] = Hy*self.bound['Hy'][0:-1,0:-1]
			
			for ix, x in indexed(self.xd[1:-1]):
				ys = range(len(self.yd[1:-1]))
				Ez = map(lambda iy: self.Ez[k,ix,iy] 
					+ (CEx * ( self.Hy[k,ix,iy] - self.Hy[k,ix,iy-1] ))
					- (CEy * ( self.Hx[k,ix,iy] - self.Hx[k,ix-1,iy] ))
					, ys)
				self.Ez[k+1,ix,1:-1] = Ez*self.bound['Ez'][ix,1:-1]

def makeHx(params):
	xy,params = params
	ix,iy = xy
	k, Ez, Hx, CH = params
	Hx = -CH * (Ez[k,ix+1,iy]-Ez[k,ix,iy])
	if k > 0:
		Hx += Hx[k-1,ix,iy]
	return Hx

def makeHy(params):
	xy,params = params
	ix,iy = xy
	k, Ez, Hx, CH = params
	Hy = -CH * (Ez[k,ix+1,iy]-Ez[k,ix,iy])
	if k > 0:
		Hy += Hy[k-1,ix,iy]
	return Hy
	