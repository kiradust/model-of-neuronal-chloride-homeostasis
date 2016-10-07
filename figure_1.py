# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:48:58 2016

@author: Kira
@title: Figure 1
"""
from plm_singlecomp_withkcc2 import plm, checkpara, clcolor
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 8,8

sym=['-',':','--','-.']

def f1b(init_cl=[1e-3,15e-3,50e-3,90e-3]):
    leg=[]    
    plt.figure()
    for i in range(len(init_cl)):
        endcl=plm(clinit=init_cl[i],tt=100)
        plt.subplot(2,1,1)
        plt.plot(endcl[11],endcl[17],color=clcolor,linestyle=sym[i])
        plt.subplot(2,1,2)
        plt.plot(endcl[11],endcl[10],'k'+sym[i])
        leg.append(str(init_cl[i]*1000)+' mM')
    #plt.legend(leg)
    plt.show()
    
    return
        
def f1c():
    offpump=plm(graph=1,ton=150,toff=350,tt=500)
    return
    
def f1d():
    c=checkpara()
    # "relative volume" is a lie --> assumes same initial [X]
    return
    