# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 17:29:22 2016

@author: Kira
@title: Figure 4
"""
from plm_singlecomp_withkcc2 import plm, zp, xcolor
import matplotlib.pyplot as plt
from figure_2or3 import minifig, minithreefig
from pylab import rcParams
rcParams['figure.figsize'] = 8,8

sym=['-',':','--','-.']

def f4a(init_x=[30e-3,100e-3,200e-3,250e-3]):   
    plt.figure()
    for i in range(len(init_x)):
        endcl=plm(xinit=init_x[i],tt=50, k_init=0)
        plt.subplot(2,1,1)
        plt.plot(endcl[11][1:-1],endcl[20][1:-1],'m'+sym[i])
        plt.subplot(2,1,2)
        plt.plot(endcl[11][0:-1],endcl[10][0:-1],'k'+sym[i])
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
    
def f4c(gX=1e-9,tt=500,xt=50):
    dex=plm(gx=gX,xt=xt,tt=tt,k_init=0.118215377)
    minithreefig([dex[11][1:-1],dex[14][1:-1],dex[13][1:-1],dex[16][1:-1],dex[10][1:-1],dex[15][1:-1]],xcolor)
    return
    
def f4e(Z=range(-120,-50)):
    zee=zp(Z=Z)
    newx=[]
    for i in zee[9]:
        newx.append(i)
    minifig([Z,zee[3],zee[2],[],[],zee[5],zee[11]])
    plt.plot(Z,newx,xcolor)
    plt.show()
    
    return
    
def f4f():
    dez=plm(gx=5e-9,xt=50,two=1,tt=200)
    minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
    return
    
def f4d(f=1e-3):
    dxe=plm(graph=1,gx=0,xt=50,two=0,tt=200,f4d=f)
    minithreefig([dxe[11][1:-1],dxe[14][1:-1],dxe[13][1:-1],dxe[16][1:-1],dxe[10][1:-1],dxe[23][1:-1]],'k')
    return