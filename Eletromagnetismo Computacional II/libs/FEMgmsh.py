# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz
"""
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.patches import Polygon


class Node(object):
    def __init__(self, *args):
        assert len(args[0]) == 4, "{0} need 4 values".format(args)
        i, x, y, z = args[0]
        self.i = int(i)
        self.x, self.y = float(x), float(y)
        self.p = array([self.x, self.y])


class Element(object):
    def __init__(self, *args, **kwargs):
        x = args[0]
        i, typ, ntags = x[:3]
        self.i, self.dim = int(i), int(typ)
        self.tags = [int(a) for a in x[3:3+int(ntags)]]
        if 'nodes' in kwargs:
            self.nodes = [kwargs['nodes'][int(a)-1] for a in x[3+int(ntags):]]
        else:
            self.nodes = [int(a) for a in x[3+int(ntags):]]


class Mesh(object):
    verbose = False

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
            self.verbose = verbose
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
            indexes = [(n1.i - 1, n2.i-1) for n1 in e.nodes for n2 in e.nodes]
            for i, p in enumerate(indexes):
                x, y = p
                self.stiff[x, y] += e.M.flat[i]
        return self.stiff

    def nodesOnLine(self, tags, indexOnly=False):
        tags = [tags] if type(tags) is int else tags
        r = [ns for tag in tags for el in self.elements
             for ns in el.nodes if el.dim == 1
             and tag in el.tags]
        if indexOnly:
            r = [n.i-1 for n in r]
        return list(set(r))

    def elementsByTag(self, tags):
        tags = [tags] if type(tags) is int else tags
        e = [el for tag in tags for el in self.elements if el.dim == 2
             and tag in el.tags]
        return list(set(e))

    def triangulate(self):
        points = array([(n.x, n.y) for n in self.nodes])
        x = points[:, 0]
        y = points[:, 1]
        return tri.Triangulation(x, y)

    def plotMesh(self, **kwargs):
        fig, ax = plt.subplots()
        if 'result' in kwargs:
            self.plotResult(result=kwargs['result'])
        elements = kwargs['elements'] if 'elements' in kwargs\
            else filter(lambda el: el.dim == 2, self.elements)
        for e in elements:
            c = e.color if 'color' in e.__dict__ else 'k'
            p = Polygon(map(lambda no: (no.x, no.y), e.nodes), closed=True,
                        fill=False, ec=c, linewidth=.5, alpha=.6)
            ax.add_patch(p)
        ax.axis('equal')
        plt.show()

    def plotResult(self, **kwargs):
        t = self.triangulate()
        plt.tricontourf(t, kwargs['result'], 15, cmap=plt.cm.rainbow)
        plt.colorbar()
        plt.show()
