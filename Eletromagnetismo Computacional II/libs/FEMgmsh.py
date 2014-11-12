# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz
"""
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.patches import Polygon

nodesInElementType = {1: 2, 2: 3}


class Node(object):
    def __init__(self, *args):
        assert len(args[0]) == 4, "{0} need 4 values".format(args)
        i, x, y, z = args[0]
        self.i = int(i)
        self.x, self.y = float(x), float(y)


class Element(object):
    def __init__(self, *args, **kwargs):
        x = args[0]
        i, typ, tags = x[:3]
        self.i, self.dim = int(i), int(typ)
        self.tags = [int(a) for a in x[3:3+int(tags)]]
        if self.dim in nodesInElementType:
            non = nodesInElementType[self.dim]*-1
            if 'nodes' in kwargs:
                self.nodes = [kwargs['nodes'][int(a)-1] for a in x[non:]]
            else:
                self.nodes = [int(a) for a in x[non:]]


class Mesh(object):
    def __init__(self, verbose=False, **kwargs):
        """
        Read a Gmsh `file` and parse all nodes and elements.
        """
        assert 'file' in kwargs, "file not specified."
        # Read and parse the Gmsh file.
        with open(kwargs['file']) as f:
            x = f.read()
            ns, ne = x.find('$Nodes\n'), x.find('$EndNodes\n')
            nodes = map(lambda n: n.split(), x[ns+7:ne].split('\n')[1:-1])
            es, ee = x.find('$Elements\n'), x.find('$EndElements\n')
            elements = map(lambda n: n.split(), x[es+10:ee].split('\n')[1:-1])
        # Map nodes and Elements.
        self.nodes = map(lambda x: Node(x), nodes)
        self.elements = map(lambda x: Element(x, nodes=self.nodes), elements)
        # Verbosity...
        if verbose:
            print "Done parsing {0} nodes and {1} elements."\
                  .format(len(nodes), len(elements))

    def getStiff(self, calc, **kwargs):
        """
        Calculate the stiffness Matrix for this mesh.
        """
        elements = [e for e in self.elements if e.dim == 2]
        ns = len(self.nodes)
        self.stiff = zeros((ns, ns))
        for e in elements:
            calc(e)
            for li in range(3):
                gi = e.nodes[li].i - 1
                for lj in range(3):
                    gj = e.nodes[lj].i - 1
                    self.stiff[gi, gj] += e.M[li, lj]
        return self.stiff

    def plotMesh(self, **kwargs):
        fig, ax = plt.subplots()
        if 'result' in kwargs:
            self.plotResult(result=kwargs['result'])
        for e in filter(lambda el: el.dim == 2, self.elements):
            p = Polygon(map(lambda no: (no.x, no.y), e.nodes), closed=True,
                        fill=False, ec='k', linewidth=.5, alpha=.6)
            ax.add_patch(p)
        ax.axis('equal')
        plt.show()

    def plotResult(self, **kwargs):
        x = []
        y = []
        z = []
        for n in self.nodes:
            x.append(n.x)
            y.append(n.y)
            z.append(kwargs['result'][n.i-1])
        t = tri.Triangulation(x, y)
        #plt.tricontour(t, z, 15, linewidths=0.5)
        plt.tricontourf(t, z, 15, cmap=plt.cm.rainbow)
        plt.colorbar()
