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

zp(Z,p) 
==> runs the parametric solution at the pump rate defined by p for different z values
Z := desired array of impermeant anion charges, multiplied by 100 (i.e. input as range(start*100, end*100))
p := initial log pump value satisfying the pump rate given by 10**(p)/(F*R) {default should be -5}

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

#constants
R=25.69*1e-3
F=96485.0 #R (RT/F) in Volts, where F is Faraday's constant in C/mol
n=200 #points to plot 
gna=3e-9
gk=5e-8
gcl=5e-8 #gna,gk,gcl: conductances in mS/cm^2 conv to S/dm^2 (10^-3/10^-2) - corrected for neuron
gkcc=0e-8
ck=2
cna=3 #cna,ck: pump (ATPase) stoichiometries
rad=5*1e-5 #radius in um convert to dm
length=100*1e-5 #length in um converted to dm
nao=138e-3
clo=119e-3
ko=2.8e-3 #nao,clo,ko: extracellular concentrations (mM converted to M)
z=-0.85 #intracellular (and extracellular) charge of impermeant anions
pkcc=0.0
gamma=gna/gk
alpha=1.0/(gk*gcl+gcl*gkcc-gk*gkcc)
beta=1.0/(gk*gcl-gkcc*gcl+gk*gkcc)
nae=nao
ke=ko
cle=clo
xe1=-1*(cle-nae-ke)
xe=xe1/2.0
ose=xe1+cle+nae+ke
P=range(-8000,-4800)

def plm(p=1e-5/F,graph=0,pkcc=pkcc,gx=0,xt=1000):
    #create plotting arrays
    Vm=[]
    K=[]
    Na=[]
    Cl=[]
    W=[]
    X=[]
    time=[]
    
    dt=1e-3 #zero time, dt time step
    tt=4000.0 #total time
    if p<10**(-6)/F:
        tt=30000.0
    ts=tt/n #plotting timestep
    tit=int(round(tt/dt)) #tt total time, tit - total number of steps 
    ton=0
    toff=15000 #times when pump turned on & off 
    ctr=1 #counter for plotting points
    t=0 #real time
    sw=0 #switch for ATPase action 
    
    w=np.pi*rad**2*length #initial volume in liters
    w1=w #initial volume stored for graphing later
    Ar=4e6 #area constant (F and H method)
    C=7e-6 #capacitance (F/dm^2)
    FinvCAr=F/(C*Ar) #(F/C*area scaling constant)
    
    na=50e-3
    k=80e-3
    cl=((ose-na-k)*z+na+k)/(1+z)
    x=(cl-na-k)/z #na,k,cl,x: intracellular starting concentrations
    print cl
    
    for i in range(2,tit): #loop over time
        if (toff>t) and (t>ton):
            sw=1
        else:
            sw=0
        if t>toff:
            sw=0     #control switch
        
        V=FinvCAr*(na+k-cl+z*x) #voltage
        
        #update arrays for plotting
        if t>=(ctr-1)*ts:
            K.append(1000*R*np.log(ko/k))
            Na.append(1000*R*np.log(nao/na))
            Cl.append(1000*R*np.log(cl/clo))
            X.append(z*1000*R*np.log(xe/x))
            W.append(w/w1)
            Vm.append(1000*V)
            time.append(t)
            ctr+=1
        
        jp=p*(na/nao)**3 #cubic pump rate update (dependent on sodium gradient)
        
        #kcc2
        #jkcc2=sw*(gk*pkcc*(k*clo-k*cl)) #Fraser and Huang
        jkcc2=sw*gk*pkcc*(K[ctr-2]-Cl[ctr-2])/10000.0 #Doyon

        #ionic flux equations
        dna=-dt*Ar*(gna*(V-R*np.log(nao/na))+cna*jp) 
        dk=-dt*Ar*(gk*(V-R*np.log(ko/k))-ck*jp+jkcc2*sw)
        dcl=dt*Ar*(gcl*(V+R*np.log(clo/cl))-jkcc2*sw) #dna,dk,dcl: increase in intracellular ion conc during time step dt
        dx=-dt*Ar*(gx*(V-R/z*np.log(xe/x)))
        na+=dna
        k+=dk
        cl+=dcl #increment concentrations
        if xt+1000>t>xt:
            x+=dx 
        
        #update volume
        osi=na+k+cl+x #intracellular osmolarity 
        w2=(w*osi)/ose #update volume 
        
        #correct ionic concentrations by volume change
        na=(na*w)/w2
        k=(k*w)/w2
        cl=(cl*w)/w2
        x=(x*w)/w2
        w=w2
        t+=dt
    
    #plot if asked    
    if graph==1:
        plt.figure
        plt.subplot(2,1,1)
        plt.plot (time,K,'-c',time,Na,'-r', time,Cl,'-g',time,X,'-b',time,Vm,'--')
        plt.legend(['E_K','E_Na','E_Cl','E_X','V_m'])
        plt.subplot(2,1,2)
        plt.plot(time,W,label='relative volume')
        plt.legend()
        plt.show()
        
    print na, k, cl, x, V
        
    return Na[-1], K[-1], Cl[-1], X[-1], Vm[-1], na, k, cl, x, V, W[10:-1], time[10:-1], Na[10:-1], K[10:-1], Cl[10:-1], X[10:-1], Vm[10:-1]

def zplm(z,gkcc=0,gcl=5e-8):
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
    
    return pi, ena, ek, ecl, exi, ev, nai, ki, cli, xi, vm
    
def checkpara():
    ti=[[],[],[],[],[],[],[],[],[],[],[]]
    T=[-8000,-7000,-6000,-5500,-5000,-4500,-4000,-3000,-2000,-1000]
    
    for k in T:
        q=10**(k/1000.0)/F
        a=plm(p=q)
        for i in range(10):
            ti[i].append(a[i])
    
    para=zplm(z,0,gcl)
    
    plt.plot(para[0],para[1],'r',para[0],para[2],'c',para[0],para[3],'g',para[0],para[4],'b',para[0],para[5],'k',T,ti[0],'-or',T,ti[1],'-oc',T,ti[2],'-og',T,ti[3],'-ob',T,ti[4],'-ok')
    plt.title('parametric plot vs time series runs: ion concentrations and membrane potential over log(cubic pump rate)')
    plt.xlabel('1000.F.log(pump rate)')
    plt.ylabel('mV')
    plt.show()
    return

def zp(Z,p):
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
    
    for u in Z:
        z=u/100.0        
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
    return pi, ena, ek, ecl, exi, ev, nai, ki, cli, xi, vm

def kcc2p(G,p=-5,z=-0.85):
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
    
    plt.figure()
    plt.plot(G,nai,'r',G,ki,'c',G,cli,'g',G,xi,'b',G,vm,'k')
    plt.title('parametric plot: ion concentrations and membrane potential over gkcc2 conductances at log pump rate of '+str(p))
    plt.xlabel('S (e10)')
    plt.ylabel('V                    M')
    plt.show()
    
    plt.figure()
    plt.plot(G,ena,'r',G,ek,'c',G,ecl,'g',G,exi,'b',G,ev,'k')
    plt.title('parametric plot: ionic reversals and membrane potential over gkcc2 conductances at log pump rate of '+str(p))
    plt.xlabel('S (e10)')
    plt.ylabel('mV')    
    plt.show()
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

def deltax(GX,z=-0.85,pip=1e-5/F,gkcc=[0e-8]):
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
    
    return a[11], N, K, C, X, V, W