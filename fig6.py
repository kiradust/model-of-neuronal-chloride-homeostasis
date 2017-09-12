from plm_singlecomp_withkcc2 import gkcc, gcl, gk, gna, ose, cle, nae, ke, xe, gkcc, F, R, default_P, z, beta
from plotting import clcolor, kcolor, nacolor, xcolor
import numpy as np
import matplotlib.pyplot as plt

def osmoplm(os=range(-100,100),molinit=0):
    nai=[]
    ki=[]
    cli=[]
    xi=[]
    vm=[]
    oi=[]
    ena=[]
    ek=[]
    ecl=[]
    exi=[]
    ev=[]
    w=[]
    
    q=10**(default_P/10000)/(F*R)
    
    for o in os:
        theta=(-z*(ose+o/1000.0)+np.sqrt(z**2*(ose+o/1000.0)**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))    
        v=(-np.log(theta))*R
        vm.append(v)
        nai.append(nae*np.exp(-v/R-3*q/gna))
        ki.append(ke*np.exp(-v/R+2*q*(gcl+gkcc)*beta))
        cli.append(cle*np.exp(+v/R-2*q*gkcc*beta))
        xi.append(ose-nai[-1]-cli[-1]-ki[-1])
        
        ek.append(1000*R*np.log(ke/ki[-1]))
        ena.append(1000*R*np.log(nae/nai[-1]))
        ecl.append(1000*R*np.log(cli[-1]/cle))
        exi.append(z*1000*R*np.log(xe/xi[-1]))
        ev.append(1000.0*v)
    
    plt.figure()
    plt.plot(os,ecl,color=clcolor)
    plt.plot(os,ek,color=kcolor)
    plt.plot(os,ena,color=nacolor)
    plt.plot(os,xi,color=xcolor)
    plt.plot(os,ev,'k--')
    plt.ylabel('mV')
    plt.xlabel('osmolarity')
    plt.show()
    
    return oi, ena, ek, ecl, exi, ev, nai, ki, cli, xi, vm, w

a = osmoplm()