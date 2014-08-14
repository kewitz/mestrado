# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz

Created on Thu May 15 11:50:45 2014
"""
import numpy as np
import sympy as sp

# Snippets and Constants
irange = lambda a: zip(range(len(a)),a)

class MOM:
    def __init__(self,lmn,gm,fx):
        """
        Inicia uma instância dos métodos dos momentos dada a expressão que
        define Lmn e Gm.

        Parâmetros
        ----------
        lmn, gm : string
            Função utilizada para definir `l[m,n]` e `g[m]` respectivamente.
        """
        lmn = sp.sympify(lmn.replace("N","ns")) if type(lmn) == str else lmn
        gm = sp.sympify(gm.replace("N","ns")) if type(gm) == str else gm
        fx = sp.sympify(fx.replace("N","ns")) if type(fx) == str else fx

        self.lmn = lmn
        self.gm = gm
        self.func = fx

        assert dir(self.lmn).count('evalf') == 1, "lmn must be string or Sympy expression."
        assert dir(self.gm).count('evalf') == 1, "gm must be string or Sympy expression."
        assert dir(self.func).count('evalf') == 1, "fx must be string or Sympy expression."

    def getAlphas(self,ns):
        """
        Retorna os valores de alpha dado a quantidade de `ns`.

        Parâmetros
        ----------
        ns : int
            Número de N.

        Exemplo
        -------
        >>> lmn = "( (m/n)**(n-1) )*n*(n+1)"
        >>> gm = "4 - 12*( (m / N)**2 )"
        >>> alphas(2,lmn,gm)
        matrix([[ 5.],
            [-3.]])
        """
        x, m, n, N = sp.symbols("x m n N")
        ms = ns

        L = np.matrix(np.zeros((ms,ns)))
        G = np.matrix(np.zeros((ns)))
        for i in range(ms):
            m = float(i+1.0)
            G[0,i] = self.gm.evalf(subs={'m':m,'ns':ns})
            for j in range(ns):
                n = float(j+1.0)
                L[i,j] = self.lmn.evalf(subs={'m':m,'ns':ns,'n':n})
        self.alphas = np.linalg.solve(L,G.T)
        return self

    def solveFor(self,domain):
        """
        Resolve a função fx para o domínio especificado.

        Example
        -------

        >>> lmn = "( (m/n)**(n-1) )*n*(n+1)"
        >>> gm = "4 - 12*( (m / N)**2 )"
        >>> x = np.arange(0,1,.1)
        >>> p1 = solve(alphas(1,lmn,gm),"x*(1-x**n)",x)
        array([ 0.  , -0.36, -0.64, -0.84, -0.96, -1.  , -0.96, -0.84, -0.64, -0.36])
        """
        ns = len(self.alphas)
        result = np.zeros(domain.shape)
        for ix,x in irange(domain):
            result[ix] = float()
            for ia, a in irange(self.alphas):
                result[ix] += a*(self.func.evalf(subs={'ns':ns,'n':ia+1,'x':x}))
        return result
