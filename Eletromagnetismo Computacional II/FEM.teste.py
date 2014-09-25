# -*- coding: utf-8 -*-
from libs import FEM

m = FEM.Mesh("""/home/leo/Documents/Master/Eletromagnetismo Computacional 2/vanti/teste.malha""")
m.solve()
m.plotMesh()
