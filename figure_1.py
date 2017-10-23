# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:48:58 2016

@author: Kira
@title: Figure 1
"""
from plm_singlecomp_withkcc2 import plm, checkpara, F
from plotting import clcolor
import matplotlib.pyplot as plt
from pylab import rcParams
import numpy as np
rcParams['figure.figsize'] = 8,8

sym=['-',':','--','-.']

def f1b(init_cl=[1e-3,15e-3,50e-3,90e-3]):
    leg=[]    
    plt.figure()
    for i in range(len(init_cl)):
        endcl=plm(clinit=init_cl[i],tt=1800,osmofix=False,k_init=0)
        plt.subplot(2,1,1)
        plt.plot(endcl[11],endcl[17],color=clcolor,linestyle=sym[i])
        plt.subplot(2,1,2)
        plt.plot(endcl[11],endcl[10],'k'+sym[i])
        leg.append(str(init_cl[i]*1000)+' mM')
    #plt.legend(leg)
    plt.savefig('f1b.eps')
    plt.show()
    return
        
def f1c(new=0,l='-',g1=0,g2=0,g3=0):
    g=1
    if new!=0:
        g=2
    offpump=plm(graph=g,ton=3000,toff=9000,tt=12000,title='f1c.eps',neww=new,ls=l,a0=g1,a1=g2,a2=g3)
    return offpump
    
def f1d():
    c=checkpara()
    # "relative volume" is a lie --> assumes same initial [X]
    return