# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:48:58 2016

@author: Kira
@title: Figure 1
"""
from plm_singlecomp_withkcc2 import plm, F, zplm
from plotting import clcolor
from matplotlib import gridspec
import matplotlib.pyplot as plt
from pylab import rcParams
import numpy as np
rcParams['figure.figsize'] = 8,8
from plotting import clcolor, kcolor, xcolor,nacolor,wcolor

sym=['-',':','--','-.']

def f1b(init_cl=[1e-3,15e-3,50e-3,90e-3]):
    leg=[]
    print "Figure 1B"
    plt.figure()
    for i in range(len(init_cl)):
        endcl=plm(clinit=init_cl[i],tt=1800,osmofix=False,k_init=0)
        plt.subplot(2,1,1)
        plt.plot(endcl[11][13:-1],endcl[17][13:-1],color=clcolor,linestyle=sym[i])
        plt.subplot(2,1,2)
        plt.plot(endcl[11][13:-1],endcl[10][13:-1],'k'+sym[i])
        leg.append(str(init_cl[i]*1000)+' mM')
    #plt.savefig('f1b.eps')
    plt.show()
    return
        
def f1c(new=0,l='-',g1=0,g2=0,g3=0):
    g=1
    if new!=0:
        g=2
    print "\nFigure 1C"
    offpump=plm(graph=g,ton=3000,toff=9000,tt=12000,title='f1c.eps',neww=new,ls=l,a0=g1,a1=g2,a2=g3)
    return offpump
    
def f1d(time=25000):
    T=[-7000,-6000,-5000,-4500,-4000,-3500,-3000,-2000,-1000,0,1000,2000]
    ti=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    
    # numerical solutions
    for k in T:
        q=10**(k/1000.0)/F
        # optimisations (alternative is to start sims at any values and allow each to run for a long time)
        if k==-7000:
            a=plm(p=q,tt=time,graph=0,k_init=0,xinit=30e-3,clinit=120e-3,na_init=140e-3,f1d=True)
        elif k==-6000: 
            a=plm(p=q,tt=time,graph=0,k_init=0,xinit=36e-3,clinit=113e-3,na_init=140e-3,f1d=True)
        else:
            a=plm(p=q,tt=time,graph=0,k_init=0,xinit=75e-3,clinit=75e-3,na_init=135e-3,f1d=True)
        for i in range(25):
            ti[i].append(a[i])
    
    # parametric solutions
    molinit=plm(gx=1e-8,xt=25,tt=100,two=1,paratwo=True,moldelt=0)
    para=zplm(molinit=molinit)
    
    # plotting
    print "\nFigure 1D"
    gs = gridspec.GridSpec(3, 1, height_ratios=[1.5, 1, 1]) 
    plt.subplot(gs[0])
    plt.plot(para[0],para[8],color=clcolor,linestyle='-')
    plt.plot(para[0],para[7],color=kcolor,linestyle='-')
    plt.plot(para[0],para[6],color=nacolor,linestyle='-')
    plt.plot(para[0],para[9],color=xcolor,linestyle='-')
    plt.scatter(T,ti[0],color=nacolor,marker='o')
    plt.scatter(T,ti[1],color=kcolor,marker='o')
    plt.scatter(T,ti[2],color=clcolor,marker='o')
    plt.scatter(T,ti[3],color=xcolor,marker='o')
    plt.subplot(gs[1])
    plt.plot(para[0],para[10],'k-')
    plt.plot(T,ti[4],'ko')
    plt.subplot(gs[2])
    plt.plot(para[0],para[11],color=wcolor,linestyle='-')
    plt.plot(T,ti[21],'ko')
    #plt.savefig('f1d.eps')
    plt.show()
    return ti