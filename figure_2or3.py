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

def delta_gs(Gk=[50],Gna=[10],Gkcc=[10],Gcl=[10],molinit=0):
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
    df=[]
    chosen=[]
    q=10**(default_P/10000.0)/(R*F)
    for i in range(max(len(Gcl),len(Gkcc),len(Gk),len(Gna))):
        for a in Gk, Gna, Gkcc, Gcl:
            if len(a)<=i:
                a.append(a[-1])
            if len(a)>i+1:
                chosen=a
        gkcc=Gkcc[i]*1.0e-8
        gna=Gna[i]*1.0e-8
        gk=Gk[i]*1.0e-8
        gcl=Gcl[i]*1.0e-8
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
        df.append(ev[-1]-ecl[-1])
    
    return np.log10(chosen),ecl,ek,ena,df,ev,w

def f2a():
    minifig(delta_gs(Gk=range(1,1000),Gna=[20],Gkcc=[20],Gcl=[20]),x=np.log10(70))
    plt.savefig('f2a.eps')
    plt.show()
    minifig(delta_gs(Gna=range(1,1000),Gk=[70],Gkcc=[20],Gcl=[20]),x=np.log10(20))
    plt.savefig('f2b.eps')
    plt.show()
    minifig(delta_gs(Gcl=range(1,1000),Gk=[70],Gna=[20],Gkcc=[20]),x=np.log10(20),yl=[[-100,-60],[1.9e-12,2.05e-12]])
    plt.savefig('f2c.eps')
    plt.show()
    minifig(delta_gs(Gcl=range(1,1000),Gk=[70],Gna=[20],Gkcc=[0]),x=np.log10(20),yl=[[-100,-60],[1.9e-12,2.05e-12]])
    plt.savefig('f2d.eps')
    plt.show()
    #minifig(delta_gs(Gkcc=d_range,Gk=[50],Gna=[50],Gcl=[10]))
    return
    
def f3a():
    dg=plm(tk=60,tt=300)
    minithreefig([dg[11][1:-1],dg[14][1:-1],dg[13][1:-1],dg[16][1:-1],dg[10][1:-1],dg[24][1:-1]],'k',yl=[[-100,-70],[1.9e-12,2.05e-12],[0,2.4e-6]])
    plt.savefig('f3a.eps')
    plt.show()
    print dg[24][-1]
    print dg[14][0]
    print dg[14][-1]
    return
    
def f3b():
    dg=delta_gs(Gkcc=range(1,1000),Gna=[20],Gk=[70],Gcl=[20])
    minithreefig([dg[0],dg[1],dg[2],dg[5],dg[-1],dg[4]],'k',x=np.log10(20),yl=[[-100,-60],[1.9e-12,2.05e-12],[0,24]])
    plt.savefig('f3b.eps')
    plt.show()
    return
    
def f3d():
    molint=plm(gx=1e-8,xt=25,tt=100,two=1,paratwo=True,moldelt=0)
    gp=zplm(molinit=molint)
    minifig([gp[0],gp[3],gp[2],[],[],gp[5],gp[11]],x=0)
    plt.savefig('f3d.eps')
    plt.show()
    return