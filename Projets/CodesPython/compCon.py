# décomposer en composantes connexes :
#   on passe en arguments :
#               -mesh : le fichier du mesh (fichier gmsh)
#               -lvlSet : le fichier level set qu'on va utiliser
#   on sort :
#               - un fichier Res+"nom du mesh de depart"
#               - un fichier texte "nbComp" avec le nombre de composantes interieures et le nombre de composantes extérieures

import inspect;
import numpy as np;
from mesh import *;
import argparse;


parser=argparse.ArgumentParser();
parser.add_argument('-mesh', help='Name of the mesh you want to deal with');
parser.add_argument('-lvlSet', help='Name of the levelset you want to deal with');
args=parser.parse_args();

nomMesh=args.mesh;
nomLvlSet=args.lvlSet;

print(nomMesh);
print(nomLvlSet);


print("load ",nomMesh);
print("load ",nomLvlSet);


M=Mesh(nomMesh);
P=P1Function(M,nomLvlSet);
#P.testIsSolid();

[compInt,compExt]=P.connectedComponent();
#print("nb de composantes intérieures : ",compInt)
#print("nb de composantes extérieures : ",compExt)

nomSortie="Res"+nomMesh;
M.saveFormat2(nomSortie);



f=open("nbComp","w");
#f.write("nb composantes interieures :"+"\n")
f.write(str(compInt)+"\n")
#f.write("nb composantes exterieures :"+"\n")
f.write(str(compExt))

