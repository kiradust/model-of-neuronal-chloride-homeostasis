# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 17:29:22 2016

@author: Kira
@title: Figure 4
"""
from plm_singlecomp_withkcc2 import plm, zp
import matplotlib.pyplot as plt
from figure_2or3 import minifig, minithreefig

sym=['-',':','--','-.']

def f4a(init_x=[30e-3,100e-3,200e-3,250e-3]):   
    plt.figure()
    for i in range(len(init_x)):
        endcl=plm(xinit=init_x[i],tt=50)
        plt.subplot(2,1,1)
        plt.plot(endcl[11],endcl[20],'b'+sym[i])
        plt.subplot(2,1,2)
        plt.plot(endcl[11],endcl[10],'k'+sym[i])
    plt.show()
    
    return
    
def f4b(init_x=range(1,251,19)):
    endv=[]
    endk=[]
    endcl=[]
    endw =[]
    
    for i in init_x:
        end=plm(xinit=i)
        endv.append(end[9])
        endk.append(end[6])
        endcl.append(end[7])
        endw.append(end[21])

    minifig([init_x,endcl,endk,[],[],endv,endw])
    plt.show()
    
    return
    
def f4c(gX=1e-9):
    dex=plm(gx=gX,xt=50,tt=200)
    minithreefig([dex[11][1:-1],dex[14][1:-1],dex[13][1:-1],dex[16][1:-1],dex[10][1:-1],dex[15][1:-1]],'b')
    return
    
def f4e(Z=range(-120,-50)):
    zee=zp(Z=Z)
    newx=[]
    for i in zee[9]:
        newx.append(i*6)
    minifig([Z,zee[3],zee[2],[],[],zee[5],zee[11]])
    plt.plot(Z,newx,'b')
    plt.show()
    
    return