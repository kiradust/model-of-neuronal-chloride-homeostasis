# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 15:43:20 2016

@author: Kira

PLM model including KCC2

:::FUNCTIONS:::

plm(p,graph)
==> runs a time series plm run
p := desired pump rate
graph := {1 iff a graph is desired as output; and any other value otherwise}

zplm(z,gkcc,gcl)
==> runs the parametric solution over log pump rates in P for impermeant anion charge of z
z := charge of impermeant anions
gkcc := desired KCC2 conductance constant {default is 0}
gcl := desired Cl- conductance constant {default is 5e-8}

checkpara() 
==> checks the parametric solution over log pump rates in P against time series runs at points in P

zp(Z,p,gkcc) 
==> runs the parametric solution at the pump rate defined by p for different z values
Z := desired array of impermeant anion charges, multiplied by 100 (i.e. input as range(start*100, end*100))
p := initial log pump value satisfying the pump rate given by 10**(p)/(F*R) {default should be -5}
gkcc := kcc2 rate (NOT scaled i.e. add e-8 etc to it)

kcc2p(G,p,z)
==> runs the parametric solution at the pump rate defined by p for different gkcc2 values
G := desired array of gkcc2 (conductance through KCC2 co-transporter) values * 10**(10)
p := initial log pump value satisfying the pump rate given by 10**(p)/(F*R) {default should be -5}
z := charge of impermeant anions

kcc2para(G,z)
==> runs the parametric solution over log pump rates in P for different gkcc2s, for fixed impermeant anion charge of z
G := desired array of gkcc2 (conductance through KCC2 co-transporter) values [should be discrete here] * 10**(10)
z := charge of impermeant anions

gclp(G,p,z,gkcc2)
==> runs the parametric solution at the pump rate defined by p for different gcl values and a given gkcc2 conductance
G := desired array of gcl (chloride conductance) values * 10**(10)
p := initial log pump value satisfying the pump rate given by 10**(p)/(F*R) {default should be -5}
z := charge of impermeant anions
gkcc2 := desired KCC2 conductance constant

gclpara(G,z,gkcc2)
==> runs the parametric solution over log pump rates in P for different gcls, for fixed impermeant anion charge of z
G := desired array of gcl (chloride conductance) values * 10**(10)
z := charge of impermeant anions
gkcc2 := desired KCC2 conductance constant

deltax(X,z,p,gkcc)
==> runs the time series solution plm with a flux in impermeant anions allowed between 1000 and 3000 seconds
X := desired conductances of X (gx)
p := initial pump rate {default should be e-5/F}
z := charge of impermeant anions
gkcc := desired KCC2 conductance constant
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pylab import rcParams
rcParams['figure.figsize'] = 8,8
from plotting import clcolor, kcolor, xcolor,nacolor,wcolor

#constants
R=26.725*1e-3
F=96485.0 #R (RT/F) in Volts, where F is Faraday's constant in C/mol, and T is 37 deg C
n=200 #points to plot 
gna=1e-8
gk=5e-8
gcl=1e-8 #gna,gk,gcl: conductances in mS/cm^2 conv to S/dm^2 (10^-3/10^-2) - corrected for neuron
gkcc=1e-8 #1 is 'high' (Doyon) - use Chris's?
ck=2
cna=3 #cna,ck: pump (ATPase) stoichiometries
rad=0.5*1e-5 #radius in um convert to dm
length=130*1e-5 #length in um converted to dm
nao=145e-3
clo=119e-3
ko=3.5e-3 #nao,clo,ko: extracellular concentrations (mM converted to M)
z=-0.85 #intracellular (and extracellular) charge of impermeant anions
pkcc=1e-8
gamma=gna/gk
beta=1.0/(gk*gcl+gkcc*gk+gcl*gkcc)
nae=nao
ke=ko
cle=clo
xe1=-1*(cle-nae-ke)
xe=xe1*0.2
ose=xe1+cle+nae+ke
P=range(-44279,-44278) 
P=range(-70000,-42500)
default_p=-2.500
default_P=-44279.0
vw=0.018 #partial molar volume of water, dm3/mol
pw=0.0015 #osmotic permeability, biological membrane (muscle? unknown), dm s
km=6*10**(-6) #extensional rigidity of RBC at 23 deg, Mohandas and Evans (1994), N/m

def plm(p=(10**(default_p))/(F),graph=0,pkcc=gkcc,gx=0,xt=100000,os_init=ose,clinit=5.2e-3,toff=150000,ton=150000,tt=200,xinit=154.9e-3,two=0,xe=xe,f4d=0,ke=ke,n=1800,k_init=103.9e-3,tk=100000,ratio=0.98,xend=120,osmofix=False,paratwo=False,moldelt=1e-13,xflux=0,z=z,dz=0,Zx=-1,ztarget=-100,length=length,areascale=1,rad=rad,title='fig.eps',neww=0):
    #create plotting arrays
    Vm=[]
    K=[]
    Na=[]
    Cl=[]
    W=[]
    X=[]
    time=[]
    Cl2=[]
    X2=[]
    Na2=[]
    K2=[]
    z_delt=[]
    xe_delt=[]
    gkcc_delt=[]
    
    dt=1e-3 #zero time, dt time step
    ts=tt/n #plotting timestep 
    ctr=1 #counter for plotting points
    t=0 #real time
    sw=1 #switch for ATPase action 
    
    w=np.pi*rad**2*length #initial volume in liters
    sa=2*np.pi*rad*(length)
    w1=w #initial volume stored for graphing later
    Ar=2.0/rad #area constant (F and H method)
    if areascale==0 or areascale==1:
        Ar=sa/w
    C=2e-4 #capacitance (F/dm^2)
    FinvCAr=F/(C*Ar) #(F/C*area scaling constant)
    sarest=sa
    
    na=33e-3
    x=xinit
    #cl=((os_init-na-k)*z+na+k)/(1+z)
    cl=clinit
    #x=(cl-na-k)/z #na,k,cl,x: intracellular starting concentrations
    k=k_init
        
    if osmofix==True:
        if xinit==0:
            x=(os_init-2*cl)/(1-z)
        else:
            cl=(os_init+(z-1)*x)/2.0
            print cl
    
    if k_init==0:
        k=cl-z*x-na
    print "k_init: "+str(k)
    print "ose: "+str(k+cl+x+na)
    print "z_aim: "+str(ztarget) +" with zflux of "+str(Zx)
    xm=x*ratio
    xtemp=x*(1-ratio)
    zxm=z
    zx=z
    cle=clo
    pd=-20
    
    if two==1:
        zx=Zx
        zxm=(z*x-zx*xtemp)/xm
        if paratwo==True:
            return (w*xinit)
        
    while t < tt: #loop over time              
        V=FinvCAr*(na+k-cl+z*x) #voltage
        
        #update arrays for plotting
        if t>=(ctr-1)*ts:
            K.append(1000*R*np.log(ke/k))
            K2.append(1000*k)
            Na.append(1000*R*np.log(nao/na))
            Na2.append(1000*na)
            Cl.append(1000*R*np.log(cl/cle))
            Cl2.append(cl*1000.0)
            X.append(z*1000*R*np.log(xe/x))
            X2.append(1000*(x))
            W.append(w)
            Vm.append(1000*V)
            time.append(t)
            z_delt.append(z)
            xe_delt.append(xe)
            gkcc_delt.append(pkcc)
            ctr+=1
        
        if tk+180>t>tk:
            pkcc += 1e-12    #control switch for gkkc ramp

        if dz!=0 and xt<t<xt+120 and xtemp>0 and xm>0:
            xtemp+=dz
            xm-=dz

        if two==1:
            z=(zxm*xm+zx*xtemp)/(xm+xtemp)
        
        if f4d!=0:
            if xt+120>t>xt:
                xe+=f4d*4e-4
                cle-=f4d*4e-4*1
        
        if (toff>t) and (t>ton):
            sw=0
        elif t>toff:
            sw=1
            if pd<default_p:
                pd+=5e-6
                p=(10**(pd))/(F)
            jp=p*(na/nao)**3
        else:
            jp=p*(na/nao)**3 #cubic pump rate update (dependent on sodium gradient)
        
            # constant/linear made this worse :(
        #kcc2
        #jkcc2=50.0*pkcc*(ke*cle-k*cl) #Fraser and Huang
        jkcc2=pkcc*(K[ctr-2]-Cl[ctr-2])/1000.0 #Doyon

        #ionic flux equations
        dna=-dt*Ar*(gna*(V-R*np.log(nao/na))+cna*jp*sw) 
        dk=-dt*Ar*(gk*(V-R*np.log(ke/k))-ck*jp*sw-jkcc2)
        dcl=dt*Ar*(gcl*(V+R*np.log(cle/cl))+jkcc2) #dna,dk,dcl: increase in intracellular ion conc during time step dt
        dx=-dt*Ar*zx*(gx*(V-R/zx*np.log(xe/(xtemp))))
        na+=dna
        k+=dk
        cl+=dcl #increment concentrations
        
        if xend==0 and (t>xt):
            if (np.abs(x*w-xinit*w1)<moldelt) and (abs((np.abs(z)-np.abs(ztarget)))>0.001) and (min(z,zx)<=ztarget<=max(z,zx)):
                if xflux==0:
                    xtemp+=dx
                    tt=t+180
                else:
                    xtemp+=xflux
                    tt=t+180
            else:
                if (min(z,zx)<=ztarget<=max(z,zx)):
                    print 'anions stopped diffusing at '+str(t)
                    xend=1
                else:
                    if xflux!=0 and tt<1000 and xtemp>0 and (min(zxm,zx)<=ztarget<=max(zxm,zx)):
                        xtemp-=xflux
                        tt=t+50
                    else:
                        print 'anions stopped diffusing at '+str(t)
                        xend=1
                
        if xt+xend>t>xt:
            if xflux!=0:
                xtemp+=xflux
            else:
                xtemp+=dx 
            
        #update volume
        x=xm+xtemp
        osi=na+k+cl+x #intracellular osmolarity 
        ose=nae+ke+cle+xe+xe1*0.8
        w2=(w*osi)/ose #update volume
        
        if neww==1:
            dt=1e-3
            w2=w+dt*(vw*pw*sa*(osi-ose)+1e-7*pw*km*(sarest-sa)/sarest)
        #length=w2/(np.pi*rad**2)
        
        #correct ionic concentrations by volume change
        na=(na*w)/w2
        k=(k*w)/w2
        cl=(cl*w)/w2
        x=(x*w)/w2
        xm=(xm*w)/w2
        xtemp=(xtemp*w)/w2
        w=w2
        if areascale==1:
            rad=np.sqrt(w/(np.pi*length))
            sa=2*np.pi*rad*(length)
            Ar=sa/w
            FinvCAr=F/(C*Ar)
        elif areascale==0:
            Ar=sa/w
            FinvCAr=F/(C*Ar)
            
        t+=dt
        
    #plot if asked    
    if graph==1:
        gs = gridspec.GridSpec(3, 1, height_ratios=[1.5, 1, 1]) 
        plt.figure()
        plt.subplot(gs[0])
        plt.plot (time,Cl2,color=clcolor)
        plt.plot(time,K2,color=kcolor)
        plt.plot(time,X2,color=xcolor)
        plt.plot(time,Na2,color=nacolor)
        plt.subplot(gs[1])
        plt.plot(time,Vm,'k')
        plt.plot (time,Cl,color=clcolor)
        plt.plot(time,K,color=kcolor)
        plt.subplot(gs[2])
        plt.plot(time,W,color=wcolor,label='relative volume')
        plt.savefig(title)
        plt.show()
    
    print 'na', na, 'k', k, 'cl', cl, 'x', x, 'vm', V, 'cle', cle, 'ose', ose, 'osi', osi, 'deltx', x*w-xinit*w1
    print 'w', w, 'radius', rad
    return na, k, cl, x, V, Na[-1], K[-1], Cl[-1], X[-1], Vm[-1], W, time, Na, K, Cl, X, Vm, Cl2, Na2, K, X2, w, z_delt, xe_delt, gkcc_delt

def zplm(z=z,gkcc=gkcc,gcl=gcl,gna=gna,gk=gk,molinit=0):
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
    w=[]
    #beta=1.0/(gk*gcl-gkcc*gcl+gk*gkcc)
    for p in P:
        q=10**(p/10000.0)/(F*R)
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
        if molinit != 0:
            w.append((molinit)/xi[-1])
        else:
            w.append(0.1549/xi[-1])
    
    plt.figure()
    plt.plot(pi,ecl,color=clcolor)
    plt.plot(pi,ek,color=kcolor)
    #plt.plot(pi,ena,color=nacolor)
    #plt.plot(pi,xi,color=xcolor)
    plt.plot(pi,ev,'k--')
    plt.ylabel('mV')
    plt.xlabel('pump rate')
    plt.savefig('pump_mV.eps')
    plt.show()
    
    return pi, ena, ek, ecl, exi, ev, nai, ki, cli, xi, vm, w
    
def checkpara(time=80000):
    ti=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    T=[-7000,-6000,-5500,-5000,-4500,-4000,-3500,-2500,-1500,-500,500]
    
    for k in T:
        q=10**(k/1000.0)/F
        if k>-4000:
            time=1000
        elif k>-5500:
            time=5000
        a=plm(p=q,tt=time)
        print len(a)
        for i in range(25):
            ti[i].append(a[i])
    
    molinit=plm(gx=1e-8,xt=25,tt=100,two=1,paratwo=True,moldelt=0)
    para=zplm(molinit=molinit)
    gs = gridspec.GridSpec(3, 1, height_ratios=[1.5, 1, 1]) 
    plt.subplot(gs[0])
    plt.plot(para[0],para[8],color=clcolor,linestyle='-')
    plt.plot(para[0],para[7],color=kcolor,linestyle='-')
    plt.plot(para[0],para[6],color=nacolor,linestyle='-')
    plt.plot(para[0],para[9],color=xcolor,linestyle='-')
    plt.plot(T,ti[0],'ro')
    plt.plot(T,ti[1],'go')
    plt.plot(T,ti[2],'bo')
    plt.plot(T,ti[3],'mo')
    plt.subplot(gs[1])
    plt.plot(para[0],para[10],'k-')
    plt.plot(T,ti[4],'ko')
    plt.subplot(gs[2])
    plt.plot(para[0],para[11],color=wcolor,linestyle='-')
    plt.plot(T,ti[21],'ko')
    plt.savefig('checkpara.eps')
    plt.show()
    return ti

def zp(Z,p=default_P/10000.0,gkcc=gkcc,graph=0,molinit=0,moldelt=0):
    nai=[]
    ki=[]
    cli=[]
    xi=[]
    vm=[]
    pi=[]
    ena=[]
    ek=[]
    ecl=[]
    exi=[]
    ev=[]
    w=[]
    
    q=10**(p)/(F*R)
    
    for u in Z:
        z=u/100.0        
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
        pi.append(np.log10(F*R*q/(((np.exp(-v/R-3*q/gna)))**3)))
        
        ek.append(1000*R*np.log(ke/ki[-1]))
        ena.append(1000*R*np.log(nae/nai[-1]))
        ecl.append(1000*R*np.log(cli[-1]/cle))
        exi.append(z*1000*R*np.log(xe/xi[-1]))
        ev.append(1000.0*v)
        
        if molinit != 0:
            w.append((molinit+moldelt)/xi[-1])
        else:
            w.append(0.155157217542/xi[-1])
    
    if graph ==1:
        plt.figure()
        plt.plot(Z,nai,'r',Z,ki,'c',Z,cli,'g',Z,xi,'b',Z,vm,'k')
        plt.title('parametric plot: ion concentrations and membrane potential over z values at log pump rate of '+str(p))
        plt.xlabel('100.z')
        plt.ylabel('V                    M')
        plt.show()
        
        plt.figure()
        plt.plot(Z,ena,'r',Z,ek,'c',Z,ecl,'g',Z,exi,'b',Z,ev,'k')
        plt.title('parametric plot: ionic reversals and membrane potential over z values at log pump rate of '+str(p))
        plt.xlabel('100.z')
        plt.ylabel('mV')    
        plt.show()
    return pi, ena, ek, ecl, exi, ev, nai, ki, cli, xi, vm, w, z

def kcc2p(G,p=default_P,z=-0.85):
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
    q=10**(p)/(F*R)
    
    for u in G:
        gkcc=u*1e-10
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
    
    #plt.figure()
    #plt.plot(G,nai,'r',G,ki,'c',G,cli,'g',G,xi,'b',G,vm,'k')
    #plt.title('parametric plot: ion concentrations and membrane potential over gkcc2 conductances at log pump rate of '+str(p))
    #plt.xlabel('S (e10)')
    #plt.ylabel('V                    M')
    #plt.show()
    
    #plt.figure()
    #plt.plot(G,ena,'r',G,ek,'c',G,ecl,'g',G,exi,'b',G,ev,'k')
    #plt.title('parametric plot: ionic reversals and membrane potential over gkcc2 conductances at log pump rate of '+str(p))
    #plt.xlabel('S (e10)')
    #plt.ylabel('mV')    
    #plt.show()
    return pi, ena, ek, ecl, exi, ev, nai, ki, cli, xi, vm

def kcc2para(G,z=-0.85):
    N=[]
    K=[]
    C=[]
    X=[]
    V=[]
    sym=['-',':','--','-.']
    plt.figure()
    
    for g in range(len(G)):
        a=zplm(z,G[g]*1e-10,gcl)
        N.append(a[1])
        K.append(a[2])
        C.append(a[3])
        X.append(a[4])
        V.append(a[5])
        plt.plot(a[0],V[g],'k'+sym[g],a[0],K[g],'c'+sym[g],a[0],C[g],'g'+sym[g]) #,a[0],X[g],'b'+sym[g],a[0],N[g],'r'+sym[g]

    plt.title('parametric plot: ionic reversals and membrane potential over log pump rates for different gkcc2 conductances')
    plt.xlabel('1000.F.log(pump rate)')
    plt.ylabel('mV')    
    plt.show()
    
    return a[0], N, K, C, X, V

def gclp(G,p=-5,z=-0.85,gkcc=0e-8):
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
    q=10**(p)/(F*R)
    
    for u in G:
        gcl=u*1e-10    
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
    
    plt.figure()
    plt.plot(G,nai,'r',G,ki,'c',G,cli,'g',G,xi,'b',G,vm,'k')
    plt.title('parametric plot: ion concentrations and membrane potential over gcl conductances at log pump rate of '+str(p))
    plt.xlabel('S (e10)')
    plt.ylabel('V                    M')
    plt.show()
    
    plt.figure()
    plt.plot(G,ena,'r',G,ek,'c',G,ecl,'g',G,exi,'b',G,ev,'k')
    plt.title('parametric plot: ionic reversals and membrane potential over gcl conductances at log pump rate of '+str(p))
    plt.xlabel('S (e10)')
    plt.ylabel('mV')    
    plt.show()
    return pi, ena, ek, ecl, exi, ev, nai, ki, cli, xi, vm

def gclpara(G,z=-0.85,gkcc=0e-8):
    N=[]
    K=[]
    C=[]
    X=[]
    V=[]
    sym=['-',':','--','-.']
    plt.figure()
    
    for g in range(len(G)):
        a=zplm(z,gkcc,G[g]*1e-10)
        N.append(a[1])
        K.append(a[2])
        C.append(a[3])
        X.append(a[4])
        V.append(a[5])
        plt.plot(a[0],K[g],'c'+sym[g],a[0],C[g],'g'+sym[g],a[0],V[g],'k'+sym[g]) #a[0],N[g],'r'+sym[g],a[0],X[g],'b'+sym[g]

    plt.title('parametric plot: ionic reversals and membrane potential over log pump rates for different gcl conductances')
    plt.xlabel('S (e10)')
    plt.ylabel('mV')    
    plt.show()
    
    return a[0], N, K, C, X, V

def deltax(GX,z=-0.85,pip=1e-5/F,gkcc=[0e-8],two=0):
    N=[]
    K=[]
    C=[]
    X=[]
    V=[]
    W=[]
    sym=['-',':','--','-.']
    plt.figure()    
    for x in range(len(GX)):
        if len(gkcc)==1:
            a=plm(p=pip,graph=0,pkcc=gkcc[0],gx=GX[x]*1e-8)
        else:
            a=plm(p=pip,graph=0,pkcc=gkcc[x],gx=GX[x]*1e-8)
        N.append(a[12])
        K.append(a[13])
        C.append(a[14])
        X.append(a[15])
        V.append(a[16])
        W.append(a[10])
        plt.plot(a[11],K[x],'c'+sym[x],a[11],C[x],'g'+sym[x],a[11],V[x],'k'+sym[x]) #a[0],N[g],'r'+sym[g],a[0],X[g],'b'+sym[g]
        
    plt.title('changing absolute concentration of impermeant anions: time series') 
    plt.xlabel('time in seconds')
    plt.ylabel('mM / mV')  
    plt.show()
    
    plt.figure()
    for x in range(len(GX)):
        plt.plot(a[11],W[x],'b'+sym[x])
    plt.title('changing absolute concentration of impermeant anions: time series for relative volume')
    plt.xlabel('time in seconds')
    plt.ylabel('relative volume')
    plt.show()
    
    return a[11], C, K, V, W, X, N

def ecl(cli):
    return R*np.log(cli/clo)