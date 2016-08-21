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

zplm(z)
==> runs the parametric solution overlog  pump rates in P for impermeant anion charge of z
z := charge of impermeant anions

checkpara() 
==> checks the parametric solution over log pump rates in P against time series runs at points in P

zpara(Z,p) 
==> runs the parametric solution at the pump rate defined by p for different z values
Z := desired array of impermeant anion charges, multiplied by 100 (i.e. input as range(start*100, end*100))
p := log pump value satisfying the pump rate given by 10**(p)/(F*R)
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
len=100*1e-5 #length in um converted to dm
nao=138e-3
clo=119e-3
ko=2.8e-3 #nao,clo,ko: extracellular concentrations (mM converted to M)
z=-0.8 #intracellular charge of impermeant anions
pkcc=0.0
gamma=gna/gk
alpha=1.0/(gk*gcl+gcl*gkcc-gk*gkcc)
beta=1.0/(gk*gcl-gkcc*gcl+gk*gkcc)
ose=285e-3
nae=nao
ke=ko
cle=clo
xe=ose-ke-nae-cle
P=range(-8000,-4800)

def plm(p,graph):
    #create plotting arrays
    Vm=[]
    K=[]
    Na=[]
    Cl=[]
    W=[]
    X=[]
    time=[]
    
    dt=1e-3 #zero time, dt time step
    tt=10000.0 #total time
    if p<10**(-6)/F:
        tt=30000.0
    ts=tt/n #plotting timestep
    tit=int(round(tt/dt)) #tt total time, tit - total number of steps 
    ton=0
    toff=15000 #times when pump turned on & off 
    ctr=1 #counter for plotting points
    t=0 #real time
    sw=0 #switch for ATPase action 
    
    w=np.pi*rad**2*len #initial volume in liters
    w1=w #initial volume stored for graphing later
    Ar=4e6 #area constant (F and H method)
    C=7e-6 #capacitance (F/dm^2)
    FinvCAr=F/(C*Ar) #(F/C*area scaling constant)
    
    na=50e-3
    k=80e-3
    cl=((ose-na-k)*z+na+k)/(1+z)
    x=(cl-na-k)/z #na,k,cl,x: intracellular starting concentrations
    
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
        jkcc2=sw*(gk*pkcc*(k*clo-k*cl)) #Fraser and Huang
        #jkcc2=sw*gk*pkcc*(K[ctr-2]-Cl[ctr-2])/10000.0 #Doyon

        #ionic flux equations
        dna=-dt*Ar*(gna*(V-R*np.log(nao/na))+cna*jp) 
        dk=-dt*Ar*(gk*(V-R*np.log(ko/k))-ck*jp+jkcc2*sw)
        dcl=dt*Ar*(gcl*(V+R*np.log(clo/cl))-jkcc2*sw) #dna,dk,dcl: increase in intracellular ion conc during time step dt
        na+=dna
        k+=dk
        cl+=dcl #increment concentrations
        
        #update volume
        osi=na+k+cl+x #intracellular osmolarity 
        w2=(w*osi)/ose #update volume 
        
        #correct ionic concentrations by volume change
        na=(na*w)/w2
        k=(k*w)/w2
        cl=(cl*w)/w2
        x=(x*w)/w2
        w=w2
        
        t=t+dt
    
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
        
    return Na[-1], K[-1], Cl[-1], X[-1], Vm[-1], na, k, cl, x, V

def zplm(z):
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
        a=plm(q,0)
        for i in range(10):
            ti[i].append(a[i])
    
    para=zplm(z)
    
    plt.plot(para[0],para[1],'r',para[0],para[2],'c',para[0],para[3],'g',para[0],para[4],'b',para[0],para[5],'k',T,ti[0],'-or',T,ti[1],'-oc',T,ti[2],'-og',T,ti[3],'-ob',T,ti[4],'-ok')
    plt.title('parametric plot vs time series runs: ion concentrations and membrane potential over log(cubic pump rate)')
    plt.xlabel('1000.F.log(pump rate)')
    plt.ylabel('mV')
    plt.show()
    return

def zpara(Z,p):
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
    
    plt.plot(Z,nai,'r',Z,ki,'c',Z,cli,'g',Z,xi,'b',Z,vm,'k')
    plt.title('parametric plot: ion concentrations and membrane potential over z values at log pump rate of '+str(p))
    plt.xlabel('100.z')
    plt.ylabel('V                    M')
    plt.show()
    
    plt.plot(Z,ena,'r',Z,ek,'c',Z,ecl,'g',Z,exi,'b',Z,ev,'k')
    plt.title('parametric plot: ionic reversals and membrane potential over z values at log pump rate of '+str(p))
    plt.xlabel('100.z')
    plt.ylabel('mV')    
    plt.show()
    return pi, ena, ek, ecl, exi, ev, nai, ki, cli, xi, vm