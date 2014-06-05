# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz

Created on Thu May 15 11:50:45 2014
"""
import numpy as np
import multiprocessing

# Snippets and Constants
indexed = lambda l, offset=0: zip(np.arange(len(l))+offset,l)
f = np.float64  # Default format.
pi = f(np.pi)  # Double precision pi.
e = np.exp(f(1))  # Exponential as unit.
threads = multiprocessing.cpu_count() + 1  # Thread Count

class Yee:
    """Classe de cálculo FDTD 1D."""
    def __init__(self, verbose = True):
        self.vp = f(299792458)
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
        self.t0 = 1/fop
        self.dx = self.lamb/10.0
        self.dy = self.lamb/10.0
        self.dtal = np.sqrt(self.dx**2 + self.dy**2)
        self.dt = self.dtal/self.vp

        self.tal = self.t0/(2*np.sqrt(2*np.log(2)))
        self.setDomain()

    def setDomain(self):
        self.xd = np.arange(0,self.lamb*10,self.dx,dtype=np.float64)
        self.yd = np.arange(0,self.lamb*10,self.dy,dtype=np.float64)
        self.td = np.arange(0,self.dt*140,self.dt,dtype=np.float64)
        self.Ez = np.zeros((len(self.td),len(self.xd),len(self.yd)))
        self.Hx = np.zeros((len(self.td),len(self.xd),len(self.yd)))
        self.Hy = np.zeros((len(self.td),len(self.xd),len(self.yd)))
        self.bound = {
                'Ez': np.ones((len(self.xd),len(self.yd))),
                'Hx': np.ones((len(self.xd),len(self.yd))),
                'Hy': np.ones((len(self.xd),len(self.yd)))
            }
    
    def setLambda(self, lamb):
        """
        Ajusta os parâmetros da simulação de acordo com o comprimento de onda
        (m) de operação desejado.
        """
        self.setFreq(self.vp/lamb)

    def run(self, fx):
        """
        Executa o FDTD para os domínios configurados.
        """
        # Constantes
        z = np.sqrt(self.mu/self.eps)
        CE = z*self.dtal
        CH = (1/z)*self.dtal
        
        # Iterações no domínio do tempo
        for k, t in indexed(self.td[:-1]):
            # Excitação
            fx(k,t,self)

            # Propagação de H no espaço
            for i, x in indexed(self.xd[0:-1]):
                for j, y in indexed(self.yd[0:-1]):
                    if self.bound['Hy'][i,j] == 1: self.Hy[k,i,j] = (CH/self.dx) * (self.Ez[k,i+1,j] - self.Ez[k,i,j])
                    if self.bound['Hx'][i,j] == 1: self.Hx[k,i,j] = (-CH/self.dx) * (self.Ez[k,i,j+1] - self.Ez[k,i,j])
                    if k > 0:
                        if self.bound['Hx'][i,j] == 1: self.Hx[k,i,j] += self.Hx[k-1,i,j]
                        if self.bound['Hy'][i,j] == 1: self.Hy[k,i,j] += self.Hy[k-1,i,j]

            # Propagação de E no espaço
            for i, x in indexed(self.xd[1:-1],1):
                for j, y in indexed(self.yd[1:-1],1):
                    if self.bound['Ez'][i,j] == 1: self.Ez[k+1,i,j] = self.Ez[k,i,j] + (CE/self.dx) * (self.Hy[k,i,j] - self.Hy[k,i-1,j]) -\
            			(CE/self.dy) * (self.Hx[k,i,j] - self.Hx[k,i,j-1])