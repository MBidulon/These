import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt;

parser=argparse.ArgumentParser();
parser.add_argument('-nomVar',nargs='*');
parser.add_argument('-nomVarNonOkIt0',nargs='*');
parser.add_argument('-fichier',help='fichier de stockage des infos');

args=parser.parse_args();

nbVar=len(args.nomVar);
nbVar=nbVar-1
data=np.loadtxt(args.fichier)
abscisse=data[:,0]
abscisse2=data[1::1,0]

for i in range(nbVar-1) :
    variable=args.nomVar[i+1]
    if(variable in args.nomVarNonOkIt0):
        tab=data[1::1,i+1]
        figure=plt.figure(i+1)
        plt.plot(abscisse2,tab)
        plt.title("Evolution de "+variable)
        plt.xlabel("iteration")
        plt.ylabel(variable)
        s=variable+".png"
        figure.savefig(s)
    else :
        tab=data[:,i+1]
        figure=plt.figure(i+1)
        plt.plot(abscisse,tab)
        plt.title("Evolution de "+variable)
        plt.xlabel("iteration")
        plt.ylabel(variable)
        s=variable+".png"
        figure.savefig(s)



