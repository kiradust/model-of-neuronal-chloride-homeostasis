from plm_singlecomp_withkcc2 import gkcc, gcl, gk, gna, ose, cle, nae, ke, xe, gkcc, F, R, default_P, z, beta, plm, default_p
from plotting import clcolor, kcolor, nacolor, xcolor, twoaxes
import numpy as np
import matplotlib.pyplot as plt
from fig4 import f4e

def f6b(moldelt=0):
    print "\nFigure 6B"
    XF=[-1.20,-1.15,-1.1,-1.05,-1.0,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.501]
    XFp=[]
    cl=[]
    vm=[]
    k=[]
    w=[]
    Z=[]
    x=[]
    df=[]
    newpr=[]
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
    df2=[]
    nai2=[]
    pi, zi, zee, newx = f4e(moldelt=moldelt)
    
    for i in range(len(zi)):
        for b in XF:
            if zi[i]/100.0==b:
                #XFp.append(pi[i])
                XFp.append(default_p)
    
    for a in range(len(XFp)):
        dez=plm(z=XF[a],tt=1000,p=(10**(XFp[a]))/F,graph=0,osmofix=True,xinit=0)
        newpr.append(dez[-3])
        cl.append(dez[7])
        k.append(dez[6])
        vm.append(dez[9])
        w.append(dez[10][-1])
        Z.append(dez[22][-1]*100)
        x.append(dez[3])
        df.append(dez[9]-dez[7])
    
    pi=[]
    pi1=[]
    pi2=[]
    df3=[]
    
    for m,o in enumerate(newpr,0):
        q=10**(default_P/10000.0)/(F*R)
        pi1.append(q*R)
        z=XF[m]
        theta=(-z*(ose)+np.sqrt(z**2*(ose)**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))
        v=(-np.log(theta))*R
        nai2.append(nae*np.exp(-v/R-3*q/gna))
        pi2.append(np.log10(F*R*q/(((np.exp(-v/R-3*q/gna)))**3)))
        df3.append(1000.0*v-1000*R*np.log(cle*np.exp(+v/R-2*q*gkcc*beta)/cle))

    for m,o in enumerate(newpr,0):
        q=10**(o)/(F*R)
        z=XF[m]
        theta=(-z*(ose)+np.sqrt(z**2*(ose)**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))
        v=(-np.log(theta))*R
        nai.append(nae*np.exp(-v/R-3*q/gna))
        ki.append(ke*np.exp(-v/R+2*q*(gcl+gkcc)*beta))
        cli.append(cle*np.exp(+v/R-2*q*gkcc*beta))
        xi.append(ose-nai[-1]-cli[-1]-ki[-1])
        pi.append(q*R)
        
        ek.append(1000*R*np.log(ke/ki[-1]))
        ena.append(1000*R*np.log(nae/nai[-1]))
        ecl.append(1000*R*np.log(cli[-1]/cle))
        exi.append(z*1000*R*np.log(xe/xi[-1]))
        ev.append(1000.0*v)
        df2.append(ev[-1]-ecl[-1])
        
    Z=XF[0:14]
    twoaxes(Z,pi,pi1,nai)
    plt.figure()
    plt.plot(Z,df2,'k')
    plt.plot(Z,df3,'k--')
    #plt.savefig('f6b.eps')
    plt.show()
    return

def f6d(os=range(-100,100)):
    print "\nFigure 6D"
    
    nai=[]
    ki=[]
    cli=[]
    ecl=[]
    ev=[]
    pi=[]
    os2=[]
    os3=[]
    
    q=10**(default_P/10000.0)/(F*R)
    
    for o in os:
        theta=(-z*(ose+o/10000.0)+np.sqrt(z**2*(ose+o/10000.0)**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))
        v=(-np.log(theta))*R
        nai.append(nae*np.exp(-v/R-3*q/gna))
        ki.append(ke*np.exp(-v/R+2*q*(gcl+gkcc)*beta))
        cli.append(cle*np.exp(+v/R-2*q*gkcc*beta))
        pi.append(q*R)
        ecl.append(1000.0*R*np.log(cli[-1]/cle))
        ev.append(1000.0*v)
    
    on=[-0.1,-0.075,-0.05,-0.025,-0.01,0,0.01,0.025,0.05,0.075,0.1]
    p_plm_osmo=[]
    
    for a,i in enumerate(os,0):
        for b in on:
            if a==b*1000:
                p_plm_osmo.append(pi[i])
    
    for m in (2,7):
        if m==1:
            shape='*'
        elif m==2:
            shape='s'
        elif m==3:
            shape='o'
        for l,o in enumerate(on,0):
            d=plm(neww=m,osmofix=True,os_choose=o/10.0,tt=1000,areascale=1,xflux=o*1e-5,xt=180,xend=180,graph=0)
            oh=(d[-2]-d[-1])*10000
            os2.append(10**(d[-3])/(F))
            os3.append(oh)
    
    nai2=[]
    ki2=[]
    cli2=[]
    xi2=[]
    ecl2=[]
    ev2=[]
    pi2=[]
    
    for x,o in enumerate(os3,0):
        q=os2[x]/R
        theta=(-z*(ose+o/10000.0)+np.sqrt(z**2*(ose+o/10000.0)**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))
        v=(-np.log(theta))*R
        nai2.append(nae*np.exp(-v/R-3*q/gna))
        ki2.append(ke*np.exp(-v/R+2*q*(gcl+gkcc)*beta))
        cli2.append(cle*np.exp(+v/R-2*q*gkcc*beta))
        pi2.append(np.log10(F*R*q/(((np.exp(-v/R-3*q/gna)))**3)))
        ecl2.append(1000.0*R*np.log(cli2[-1]/cle))
        ev2.append(1000.0*v)
    
    twoaxes(os,pi,os2,nai,os3)
    plt.figure()
    plt.plot(os,np.array(ev)-np.array(ecl),color='k')
    plt.plot(os3,np.array(ev2)-np.array(ecl2),'k--')
    #plt.savefig('f6d.eps')
    plt.show()
    return 