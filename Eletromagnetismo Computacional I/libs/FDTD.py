# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz

Created on Thu May 15 11:50:45 2014
"""
import numpy as np

# Snippets and Constants
indexed = lambda l, offset=0: zip(np.arange(len(l))+offset,l)
f = np.float64  # Default format.
pi = f(np.pi)  #Double precision pi.
e = np.exp(f(1))  # Exponential as unit.

class FDTD:
    """Classe de cálculo FDTD 1D."""
    def __init__(self, verbose = True):
        self.vp = f(299792458)
        self.eps = f(8.854187817E-12)
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
        self.tc = 1/fop
        self.dx = self.lamb/10.0
        self.dt = (self.dx)/(self.vp)
        self.tal = self.tc/(2*np.sqrt(2*np.log(2)))
    
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
        dtdx = (self.vp * self.dt) / self.dx
        
        # Define domínios:
        self.xd = np.arange(0,self.lamb*10,self.dx,dtype=np.float64)
        self.td = np.arange(0,self.dt*700,self.dt,dtype=np.float64)
        self.st = np.zeros((len(self.td),len(self.xd)))
        if self.verbose: print "DT = %.3E\nDX = %.3E\n" \
            "Space/Time Domain: %i/%ipoints" \
            % (self.dt,self.dx, len(self.xd), len(self.td))
        
        # Iterações no domínio do tempo
        for k, t in indexed(self.td[:-1]):
            # Excitação            
            self.st[k,0] = fx(t)
            # Propagação no espaço
            for i, x in indexed(self.xd[1:-1],1):
                self.st[k+1,i] = dtdx*(self.st[k,i-1]+self.st[k,i+1])\
                    + 2*(1-dtdx)*self.st[k,i] - self.st[k-1,i]
            # Condição de contorno absorvente
            self.st[k+1,-1] = self.st[k, -2]

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
        self.tc = 1/fop
        self.ds = self.lamb/10.0
        self.dt = (self.ds)/(self.vp)
        self.tal = self.tc/(2*np.sqrt(2*np.log(2)))
    
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
        dtds = -1 * self.dt / self.ds
        
        # Define domínios:
        self.sd = np.arange(0,self.lamb*10,self.ds,dtype=np.float64)
        self.td = np.arange(0,self.dt*700,self.dt,dtype=np.float64)
        self.et = np.zeros((len(self.td),len(self.sd)))
        self.ht = np.zeros((len(self.td),len(self.sd)))
        
        if self.verbose: print "DT = %.3E\nDS = %.3E\n" \
            "Space/Time Domain: %i/%ipoints" \
            % (self.dt,self.ds, len(self.sd), len(self.td))
        
        # Iterações no domínio do tempo
        for k, t in indexed(self.td[:-1]):
            # Excitação            
            self.et[k,0] = fx(t)
            # Propagação de H no espaço
            for i, x in indexed(self.sd[0:-1]):
                self.ht[k,i] = (dtds/self.mu) * (self.et[k,i+1] - self.et[k,i])
                if k > 0: self.ht[k,i] += self.ht[k-1,i]
            # Propagação de E no espaço
            for i, x in indexed(self.sd[1:-1],1):
                self.et[k+1,i] = (dtds/self.eps) * (self.ht[k,i] - self.ht[k,i-1])+\
                    self.et[k,i]
            # Condição de contorno absorvente no campo elétrico
            self.et[k+1,-1] = 0