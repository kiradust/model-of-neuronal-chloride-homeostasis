# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 12:11:36 2016

@author: Kira
@title: Figure 2/3
"""
from plm_singlecomp_withkcc2 import plm, zplm, z, ose, default_P, F, nae, ke, cle, xe, R, default_p, P, gna, gk, gcl, beta, xe1, gkcc
from plotting import minifig, minithreefig, twoaxes,wcolor, kcolor, clcolor, xcolor
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams
rcParams['figure.figsize'] = 8,8
prange=range(-60000,-33000)

# delta_gs keeps the effective pump rate constant, while delta_gs3 uses a constant p and P_eff is dependent on the final sodium concentration (checks back to find corresponding values)

def delta_gs(Gk=[70],Gna=[20],Gkcc=[20],Gcl=[20],molinit=0):
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
    q=10**(default_P/10000.0)
    
    for i in range(max(len(Gcl),len(Gkcc),len(Gk),len(Gna))):
        for a in Gk, Gna, Gkcc, Gcl:
            if len(a)<=i:
                a.append(a[-1])
            if len(a)>i+1:
                chosen=a
        gkcc=Gkcc[i]*1.0e-4/F
        gna=Gna[i]*1.0e-4/F
        gk=Gk[i]*1.0e-4/F
        gcl=Gcl[i]*1.0e-4/F
        molinit=plm(gx=1e-8,xt=25,tt=100,two=1,paratwo=True,moldelt=0)
        if gk*gcl+gkcc*gcl+gk*gkcc !=0:
            beta=1.0/(gk*gcl+gkcc*gcl+gk*gkcc)
        else:
            print "invalid conductances: division by 0"
        
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
        
    print ecl[-1], ek[-1], ev[-1]
    
    return np.log10(chosen),ecl,ek,ena,df,ev,w,kpflux,kaflux

def delta_gs3(Gk=[70],Gna=[20],Gkcc=[20],Gcl=[20],molinit=0):
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
    kpflux = []
    kaflux = []
    kaflux1 = []
    kaflux2 = []
    molinit=plm(gx=1e-8,xt=25,tt=100,two=1,paratwo=True,moldelt=0)
    
    # only operate on conductance which is changing
    for i in range(max(len(Gcl),len(Gkcc),len(Gk),len(Gna))):
        for a in Gk, Gna, Gkcc, Gcl:
            if len(a)<=i:
                a.append(a[-1])
            if len(a)>i+1:
                chosen=a
        gkcc=Gkcc[i]*1.0e-4/F
        gna=Gna[i]*1.0e-4/F
        gk=Gk[i]*1.0e-4/F
        gcl=Gcl[i]*1.0e-4/F
        found=False
        
        # run through all pump rates (analytical solution)
        for l in prange:
            q=np.float64(10**(l/10000.0)/(F*R))
            
            if gk*gcl+gkcc*gcl+gk*gkcc !=0:
                beta=1.0/(gk*gcl+gkcc*gcl+gk*gkcc)
            else:
                print "invalid conductances: division by 0"

            if z==-1:
                theta=0.5*ose/(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))
            else:
                theta=(-z*ose+np.sqrt(z**2*ose**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))    
            v=(-np.log(theta))*R
            # efficient coding for expression pi=np.log10(F*R*q/(((np.exp(-v/R-3*q/gna)))**3))
            pi=np.exp(-v/R-3*q/gna)
            pi=-3*np.log10(pi)
            pi+=np.log10(F*R*q)
            
            # match constant pump rates to determine which pair (l, na) produce (very close to) the default pump rate
            # append to arrays
            if np.abs(pi-default_p)<0.01 and found==False:
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
                found=True
                kpflux.append(gk*(v-ek[-1]/1000.0))
                kaflux.append(-(2*q*R+gkcc*(ek[-1]-ecl[-1])/1000.0))
                kaflux1.append(-(2*q*R))
                kaflux2.append(-(gkcc*(ek[-1]-ecl[-1])/1000.0))
        
        # if when checking for correct pair (l, na) none is found, notify (this will mean that no figure is generated)
        if found==False:
            print 'no match',i,pi,default_p
            # may be fixed by increasing bounds of prange (i.e. larger search space)

    return np.log10(chosen),ecl,ek,ena,ki,ev,w,kpflux,kaflux,kaflux1,kaflux2
    
def f2():
    print "Figure 2A"
    dga = delta_gs3(Gk=range(1,1000),Gna=[20],Gkcc=[20],Gcl=[20])
    minifig(dga,x=np.log10(70),yl=[[-100,80],[1.9e-12,2.05e-12]])
    plt.show()
    print "\nFigure 2B"
    dga = delta_gs3(Gna=range(1,1000),Gk=[70],Gkcc=[20],Gcl=[20])
    minifig(dga,x=np.log10(20),yl=[[-100,100],[1.9e-12,2.05e-12]])
    plt.show()
    print "\nFigure 2C"
    #minifig(delta_gs3(Gcl=range(1,1000),Gk=[70],Gna=[20],Gkcc=[20]),x=np.log10(20),yl=[[-100,-60],[1.9e-12,2.05e-12]])
    dga = delta_gs3(Gcl=range(1,1000),Gk=[70],Gna=[20],Gkcc=[20])
    minifig(dga,x=np.log10(20),yl=[[-100,80],[1.9e-12,2.05e-12]])
    plt.show()
    print "\nFigure 2D"
    #minifig(delta_gs3(Gcl=range(1,1000),Gk=[70],Gna=[20],Gkcc=[0]),x=np.log10(20),yl=[[-100,-60],[1.9e-12,2.05e-12]])
    dga = delta_gs3(Gcl=range(1,1000),Gk=[70],Gna=[20],Gkcc=[0])
    minifig(dga,x=np.log10(20),yl=[[-100,80],[1.9e-12,2.05e-12]])
    plt.show()
    return
    
def f3a():
    dg=plm(tk=120,tt=600)
    print "Figure 3A"
    a1,a2,a3=minithreefig([dg[11][4:-1],dg[14][4:-1],dg[13][4:-1],dg[16][4:-1],dg[10][4:-1],dg[24][4:-1]],'k',yl=[[-100,-70],[1.92e-12,1.98e-12],[0,6e-7]])
    #plt.savefig('f3a.eps')
    plt.show()
    print dg[24][-1]
    print dg[14][0]
    print dg[17][0]
    print dg[17][-1]
    print dg[14][-1]
    return dg
    
def f3b():
    dg=delta_gs3(Gkcc=range(1,1000),Gna=[20],Gk=[70],Gcl=[20])
    print "\nFigure 3B"
    minithreefig([dg[0],dg[1],dg[2],dg[5],dg[-1],dg[4]],'k',x=np.log10(20),yl=[[-100,-60],[1.9e-12,2.05e-12],[0,24]])
    #minithreefig([10**dg[0],dg[1],dg[2],dg[5],dg[-1],dg[4]],'k',x=np.log10(20),yl=[[-100,-60],[1.9e-12,2.05e-12],[0,24]])
    #plt.savefig('f3b.eps')
    plt.show()
    return dg[0],dg[4], dg[1], dg[5]
    
def f3d(new=0,title='f3d.eps'): #doubles as f6f when new==1
    molint=plm(gx=1e-8,xt=25,tt=100,two=1,paratwo=True,moldelt=0)
    if new==0:
        gp=zplm(molinit=molint)
        print "\nFigure 3D"
        minifig([gp[0],gp[3],gp[2],[],[],gp[5],gp[11]],x=0)
    else:
        gp=zplm()
        gp_nokcc=zplm(gkcc=0)
        print "\nFigure 6F"
        plt.figure()
        plt.plot(gp[-1],np.subtract(gp[5],gp[3]),color='k')
        plt.plot(gp_nokcc[-1],np.subtract(gp_nokcc[5],gp_nokcc[3]),linestyle='--',color='k')
    #plt.savefig(title)
    plt.show()
    return

def k_vm(): #supplementary figure *note: pump rate not corrected
    vm=[]
    for i in range(1,500):
        i=i*1e-4
        xe1=-1*(cle-nae-ke)
        xe=xe1*0.2
        ose = xe1+cle+nae+i
        q=np.float64(10**(default_P/10000.0)/(F*R))
        theta=(-z*ose+np.sqrt(z**2*ose**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+i*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+i*np.exp(2*q*(gcl+gkcc)*beta))))    
        vm.append((-np.log(theta))*R*1000.0)
        
    plt.plot(np.log10(range(1,500))-1,vm,color=kcolor)
    plt.xlabel("log (extracellular potassium concentration) in mM")
    plt.ylabel("membrane potential (mV)")
    plt.savefig('k_vm.eps')
    plt.show()
    return