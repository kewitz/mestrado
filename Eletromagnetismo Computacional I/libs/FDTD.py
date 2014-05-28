# -*- coding: utf-8 -*-
"""
The MIT License (MIT)

Copyright (c) 2014 Leonardo Kewitz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Created on Thu May 15 11:50:45 2014
@author: leo
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
        self.dx = self.lamb/40.0
        self.dt = (self.dx)/(self.vp)
        self.tal = self.tc/(2*np.sqrt(2*np.log(2)))
    
    def setLambda(self, lamb):
        self.setFreq(f(1)/lamb)

    def run(self):
        """
        Executa o FDTD para os domínios configurados.
        """
        # Constantes
        dtdx = (self.vp * self.dt) / self.dx
        width = (2*np.power(self.tal,2))
        omega = 6*pi*self.fop
        # Define domínios:
        self.xd = np.arange(0,self.lamb*10,self.dx,dtype=np.float64)
        self.td = np.arange(0,self.dt*700,self.dt,dtype=np.float64)
        self.st = np.zeros((len(self.td),len(self.xd)))
        if self.verbose: print "DT = %.3E\nDX = %.3E\n" \
            "Space/Time Domain: %i/%ipoints" \
            % (self.dt,self.dx, len(self.xd), len(self.td))
        
        for k, t in indexed(self.td[:-1]):
            # Excitação
            exc = np.sin(omega*t)\
                * np.exp(-np.power(t-2*self.tc,2) / width)
            self.st[k,0] = exc
            # Propagação
            for i, x in indexed(self.xd[1:-1],1):
                self.st[k+1,i] = dtdx*(self.st[k,i-1]+self.st[k,i+1])\
                    + 2*(1-dtdx)*self.st[k,i] - self.st[k-1,i]
            self.st[k+1,-1] = self.st[k, -2]

