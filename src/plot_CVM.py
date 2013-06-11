#!/usr/bin/python
from sys import argv,stdout
from math import *
from numpy import *
from numpy import linalg as lalg
from matplotlib import pyplot as plt
import random


def plotcvm(args):

    diff=None
    R=0
    for arg in args:
        if arg.startswith("-D"):
            diff=arg[2:]
        if arg.startswith("-R"):
            row=int(arg[2:])

    f=args[1]
    M=load(f)
    N=len(M[0])
    C=zeros((N,N))
    MM=zeros((N,N))
    ff=open("var","w")
    fff=open("cvm","w")
    i=0
    while i < N:
        ff.write("%g\n"%M[i][i])
        i+=1

    ff.close()
    print "var out"
    i=2000
    while i < N:
        j=0
        while j < N:
            fff.write("%d %d %g\n"%(i+1,j+1,M[i][j]))
            j+=1
        fff.write("\n")
        i+=1e999

    fff.close()
    print "cvm out"
    i=0
    while i < N:
        j=0
        while j < N:
            C[i][j]=M[i][j]/sqrt(M[i][i]*M[j][j])
            j+=1
        i+=1


    #save("clip",MM)
    #save("clip_i",lalg.inv(MM))


    plt.imshow(abs(C))
    plt.xlim(0,N)
    plt.ylim(0,N)
    plt.figure()
    plt.imshow(M)
    plt.xlim(0,N)
    plt.ylim(0,N)
    plt.show()




if __name__ == "__main__":

    plotcvm(argv)


