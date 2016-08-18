# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 22:21:56 2016

@author: Kira
"""
import numpy as np
import matplotlib.pyplot as plt

P=range(-8000,-4700)
gna=3e-9
gk=5e-8
gcl=1e-8
gkcc=2e-8 #Doyon uses between 0.2 and 2 mS/cm2
R=25.69*1e-3
F=96485.0
nae=138e-3
cle=119e-3
ke=2.8e-3
z=-0.85
gamma=gna/gk
alpha=1.0/(gk*gcl+gcl*gkcc-gk*gkcc)
beta=1.0/(gk*gcl-gkcc*gcl+gk*gkcc)
ose=285e-3
xe=ose-ke-nae-cle
nai=[]
ki=[]
cli=[]
xi=[]
vm=[]
zi=[]
pi=[]
ena=[]
ek=[]
ecl=[]
exi=[]
ev=[]

for p in P:
    q=10**(p/1000.0)/(F*R)
    if z==-1:
        theta=0.5*ose/(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))
    else:
        theta=(-z*ose+np.sqrt(z**2*ose**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))    
    v=(-np.log(theta))*R
    vm.append(v)
    zi.append(nae*np.exp(-v/R-3*q/gna))
    nai.append(nae*np.exp(-v/R-3*q/gna))
    ki.append(ke*np.exp(-v/R+2*q*(gcl+gkcc)*beta))
    cli.append(cle*np.exp(+v/R-2*q*gkcc*beta))
    xi.append(ose-nai[-1]-cli[-1]-ki[-1])
    pi.append(1000.0*np.log10(F*R*q/(((nae*np.exp(-v/R-3*q/gna))/nae)**3)))
    
    ek.append(1000*R*np.log(ke/ki[-1]))
    ena.append(1000*R*np.log(nae/nai[-1]))
    ecl.append(1000*R*np.log(cli[-1]/cle))
    exi.append(z*1000*R*np.log(xe/xi[-1]))
    ev.append(1000.0*v)
    
plt.plot(P,nai,'r',P,ki,'c',P,cli,'g',P,xi,'b',P,vm,'k')
plt.legend(['Nai','Ki','Cli','Xi','Vm'],loc=2)
plt.show()

plt.plot(vm,pi)
plt.show()

plt.plot(pi,nai,'r',pi,ki,'c',pi,cli,'g',pi,xi,'b',pi,vm,'k')
plt.show()

plt.plot(pi,ecl,'g',pi,ek,'c',pi,ev,'k')
plt.show()
