# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 12:11:36 2016

@author: Kira
@title: Figure 2/3
"""
from plm_singlecomp_withkcc2 import plm, zplm, z, ose, default_P, F, nae, ke, cle, xe, R
from plotting import minifig, minithreefig, twoaxes,wcolor, kcolor, clcolor, xcolor
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams
rcParams['figure.figsize'] = 8,8

def delta_gs(Gk=[50],Gna=[50],Gkcc=[10],Gcl=[10],molinit=0):
    vm=[]
    cli=[]
    nai=[]
    ki=[]
    xi=[]
    w=[]
    ek=[]
    ecl=[]
    ev=[]
    exi=[]
    ena=[]
    chosen=[]
    q=10**(default_P/1000.0)/(R*F)
    for i in range(max(len(Gcl),len(Gkcc),len(Gk),len(Gna))):
        for a in Gk, Gna, Gkcc, Gcl:
            if len(a)<=i:
                a.append(a[-1])
            if len(a)>i+1:
                chosen=a
        gkcc=Gkcc[i]*1.0e-9
        gna=Gna[i]*1.0e-10
        gk=Gk[i]*1.0e-9
        gcl=Gcl[i]*1.0e-9
        molinit=plm(gx=1e-8,xt=25,tt=100,two=1,paratwo=True,moldelt=0)
        if gk*gcl+gkcc*gcl+gk*gkcc !=0:
            beta=1.0/(gk*gcl+gkcc*gcl+gk*gkcc)
        else:
            print "yep"
        
        if z==-1:
            theta=0.5*ose/(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))
        else:
            theta=(-z*ose+np.sqrt(z**2*ose**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))    
        v=(-np.log(theta))*R
        vm.append(v)
        nai.append(nae*np.exp(-v/R-3*q/gna))
        ki.append(ke*np.exp(-v/R+2*q*(gcl+gkcc)*beta))
        cli.append(cle*np.exp(+v/R-2*q*gkcc*beta))
        xi.append(ose-nai[-1]-cli[-1]-ki[-1])
        w.append((molinit)/xi[-1])
        
        ek.append(1000*R*np.log(ke/ki[-1]))
        ena.append(1000*R*np.log(nae/nai[-1]))
        ecl.append(1000*R*np.log(cli[-1]/cle))
        exi.append(z*1000*R*np.log(xe/xi[-1]))
        ev.append(1000.0*v)
    
    return chosen,ecl,ek,ena,exi,ev,w

def f2a():
    minifig(delta_gs(Gk=range(10,400),Gna=[50],Gkcc=[10],Gcl=[10]))
    plt.show()
    minifig(delta_gs(Gna=range(1,400),Gk=[50],Gkcc=[10],Gcl=[10]))
    plt.show()
    minifig(delta_gs(Gcl=range(0,400),Gk=[50],Gna=[50],Gkcc=[10]))
    plt.show()
    minifig(delta_gs(Gcl=range(1,400),Gk=[50],Gna=[50],Gkcc=[0]))
    plt.show()
    #minifig(delta_gs(Gkcc=d_range,Gk=[50],Gna=[50],Gcl=[10]))
    return
    
def f3a():
    dg=plm(tk=25,tt=100)
    minithreefig([dg[11][1:-1],dg[14][1:-1],dg[13][1:-1],dg[16][1:-1],dg[10][1:-1],dg[24][1:-1]],'k')
    plt.show()
    return
    
def f3c():
    minifig(delta_gs(Gkcc=range(0,400),Gna=[50],Gk=[50],Gcl=[10]))
    plt.show()
    return
    
def f3b():
    gp=zplm()
    minifig([gp[0],gp[3],gp[2],[],[],gp[5],gp[11]])
    plt.show()
    return