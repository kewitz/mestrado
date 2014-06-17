# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz

Created on Thu May 15 11:50:45 2014
fas"""
import numpy as np
from ctypes import *

# Macros
indexed = lambda l, offset=0: zip(np.arange(len(l))+offset,l)

fast = cdll.LoadLibrary('./src/Yee2D/yee.so')
#fast = cdll.LoadLibrary('./yee.so')

# Constantes
f = np.float64  # Default format.
pi = f(np.pi)  # Double precision pi.
e = np.exp(f(1))  # Exponential as unit.

class FYee:
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
		sx,sy = (len(self.xd),len(self.yd))

		# malloc.
		cHx = ( c_double * (sx * sy))()
		cHy = ( c_double * (sx * sy))()
		cEz = ( c_double * (sx * sy))()
		
		oldHx = ( c_double * (sx * sy))()
		oldHy = ( c_double * (sx * sy))()
		oldEz = ( c_double * (sx * sy))()

		# Iterações no domínio do tempo
		for k, t in indexed(self.td[:-1]):
			# Excitação
			fx(k,t,self)
			pack(cEz,self.Ez[k,:,:])

			# Propagação de H no espaço
				# void H (int k, unsigned int sx, unsigned int sy, double CH, double * Ez, double * Hx, double * Hy)
			pack(oldHx,self.Hx[k,:,:])
			pack(oldHy,self.Hy[k,:,:])
			fast.calcH(c_int(k), c_int(skip), c_int(sx), c_int(sy), c_double(CH), byref(cEz), byref(cHx), byref(cHy), byref(oldHx), byref(oldHy))
			self.Hx[k,:,:] = unpack(cHx,(sx,sy))*self.bound['Hx']
			self.Hy[k,:,:] = unpack(cHy,(sx,sy))*self.bound['Hy']
			
			# for i, x in indexed(self.xd[0:-1]):
			# 	for j, y in indexed(self.yd[0:-1]):
			# 		if self.bound['Hy'][i,j] == 1: self.Hy[k,i,j] = (CH/self.dx) * (self.Ez[k,i+1,j] - self.Ez[k,i,j])
			# 		if self.bound['Hx'][i,j] == 1: self.Hx[k,i,j] = (-CH/self.dx) * (self.Ez[k,i,j+1] - self.Ez[k,i,j])
			# 		if k > 0:
			# 			if self.bound['Hx'][i,j] == 1: self.Hx[k,i,j] += self.Hx[k-1,i,j]
			# 			if self.bound['Hy'][i,j] == 1: self.Hy[k,i,j] += self.Hy[k-1,i,j]

			# Propagação de E no espaço
				#int calcE (unsigned int sx, unsigned int sy, double CE, double * Ez, double * Hx, double * Hy)
			pack(oldEz,self.Ez[k,:,:])
			fast.calcE(c_int(skip), c_int(sx), c_int(sy), c_double(CEx), c_double(CEy), byref(cEz), byref(cHx), byref(cHy), byref(oldEz))
			self.Ez[k+1,:,:] = unpack(cEz,(sx,sy))*self.bound['Ez']
			# for i, x in indexed(self.xd[1:-1],1):
			# 	for j, y in indexed(self.yd[1:-1],1):
			# 		if self.bound['Ez'][i,j] == 1: self.Ez[k+1,i,j] = self.Ez[k,i,j] + (CE/self.dx) * (self.Hy[k,i,j] - self.Hy[k,i-1,j]) - (CE/self.dy) * (self.Hx[k,i,j] - self.Hx[k,i,j-1])

def pack(mem, array):
	flat = array.flatten().tolist()
	size = len(flat)
	for i in range(size):
		mem[i] = flat[i]

def unpack(mem,shape):
	flat = np.array([i for i in mem])
	return flat.reshape(shape)