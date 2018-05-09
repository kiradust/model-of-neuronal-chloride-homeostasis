# -*- coding: utf-8 -*-
"""
@author: Kira
@title: Figure 5
"""

from plm_singlecomp_withkcc2 import plm, zp, F, default_p, default_P, R, ose, cle, nae, ke, gna, gk, gcl, gkcc, beta, xe
from plotting import minifig, minithreefig, minifigtwoaxes,twoaxes,xcolor,clcolor,wcolor,kcolor,nacolor
import matplotlib.pyplot as plt
from pylab import rcParams
from fig4 import f4e
import numpy as np
rcParams['figure.figsize'] = 8,8
from matplotlib import gridspec

sym=['-k',':k','--k','-.k']

def f5a(new=0,title='f5a.eps',dz=2.5e-7,tt=1800,ham=0): #doubles as f6a when new==1
    dez=plm(dz=dz,two=1,xt=360,tt=tt,ztarget=-1)
    delta=[]
    if new==0 and ham==0:
        print "Figure 5A"
        a0,a1,a2=minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],xcolor,yl=[[-100,-70],[1.8e-12,2.4e-12],[-1.1,-0.8]])
    else:
        print "Figure 6A"
        a0,a1,a2=minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[18][1:-1],dez[22][1:-1]],xcolor,yl=[[-100,-70],[13,19],[-1.1,-0.8]])
        print (dez[16][-1]-dez[14][-1])
        
        delta=plm(dz=dz,two=1,xt=360,tt=tt,ztarget=-1,neww=4*new,hamada=ham)
        a0.plot(delta[11][1:-1],delta[14][1:-1],color=clcolor,linestyle='--')
        a0.plot(delta[11][1:-1],delta[13][1:-1],color=kcolor,ls='--')
        a0.plot(delta[11][1:-1],delta[16][1:-1],'k',ls='--')
        a1.plot(delta[11][1:-1],delta[18][1:-1],color=nacolor,ls='--') #nai
        a2.plot(delta[11][1:-1],delta[22][1:-1],color=xcolor,ls='--')
        print (delta[16][-1]-delta[14][-1])
        
        if ham==1:
            delta=plm(dz=dz,two=1,xt=360,tt=tt,ztarget=-1,neww=4,hamada=ham)
            a0.plot(delta[11][1:-1],delta[14][1:-1],color=clcolor,linestyle='-.')
            a0.plot(delta[11][1:-1],delta[13][1:-1],color=kcolor,ls='-.')
            a0.plot(delta[11][1:-1],delta[16][1:-1],'k',ls='-.')
            a1.plot(delta[11][1:-1],delta[18][1:-1],color=nacolor,ls='-.') #nai
            a2.plot(delta[11][1:-1],delta[22][1:-1],color=xcolor,ls='-.')
            print (delta[16][-1]-delta[14][-1])
    plt.savefig(title)
    plt.show()
    return dez, delta

def f5b(moldelt=0):
    print "\nFigure 5B"
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
        #dez=plm(two=1,xt=10,dz=3e-6,Zx=a,tt=5000)
        dez=plm(z=XF[a],tt=2000,p=(10**(XFp[a]))/F,graph=0,osmofix=True,xinit=0)
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
        q=10**(o)/(F*R)
        z=XF[m]
        theta=(-z*(ose)+np.sqrt(z**2*(ose)**2+4*(1-z**2)*cle*np.exp(-2*q*gkcc*beta)*(nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))/(2*(1-z)*((nae*np.exp(-3*q/gna)+ke*np.exp(2*q*(gcl+gkcc)*beta))))
        v=(-np.log(theta))*R
        nai.append(nae*np.exp(-v/R-3*q/gna))
        ki.append(ke*np.exp(-v/R+2*q*(gcl+gkcc)*beta))
        cli.append(cle*np.exp(+v/R-2*q*gkcc*beta))
        xi.append(ose-nai[-1]-cli[-1]-ki[-1])
        #pi.append(np.log10(F*R*q/(((np.exp(-v/R-3*q/gna)))**3)))
        pi.append(q*R)
        
        ek.append(1000*R*np.log(ke/ki[-1]))
        ena.append(1000*R*np.log(nae/nai[-1]))
        ecl.append(1000*R*np.log(cli[-1]/cle))
        exi.append(z*1000*R*np.log(xe/xi[-1]))
        ev.append(1000.0*v)
        df2.append(ev[-1]-ecl[-1])

    XF.pop()
    plt.figure()
    gs = gridspec.GridSpec(3,1,height_ratios=[1.5,0.5,0.5])
    a0=plt.subplot(gs[0])
    a0.plot(Z,ecl,color=clcolor)
    a0.plot(Z,cl,'bo')
    a0.plot(Z,k,'go')
    a0.plot(Z,vm,'ko')
    a0.plot(Z,ek,color=kcolor)
    a0.plot(Z,ev,'k')
    a0.set_ylim([-110,-60])
    a2=plt.subplot(gs[1])
    a2.plot(Z,xi,color='m') #conc X
    a2.plot(Z,x,'mo')
    a2.set_ylim([0.11,0.2])
    a3=plt.subplot(gs[2])
    a3.plot(Z,df,'ko')
    a3.plot(Z,df2,'k-')
    a3.set_ylim([9,12])
    #plt.savefig('f5b.eps')
    plt.show()
    return

def f5c(ttt=3000,ratio=0.8,md=1e-12,new=0,title='f5c.eps'):
    print "\nFigure 5C"
    dez=plm(gx=1e-8,xt=360,two=1,xend=0,moldelt=md,xflux=4e-8,ztarget=-1,tt=ttt,Zx=-1.5)
    print(len(dez[16]))
    a0,a1,a2=minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k',yl=[[-100,-70],[1e-12,3.3e-12],[-1.1,-0.8]])
    print (dez[16][-1]-dez[14][-1])
    if new==1:
        delta=plm(gx=1e-8,xt=240,two=1,xend=0,moldelt=md,xflux=0.3*1e-6,ztarget=-1,tt=ttt,neww=1)
        a0.plot(delta[11][1:-1],delta[14][1:-1],color=clcolor,linestyle='--')
        a0.plot(delta[11][1:-1],delta[13][1:-1],color=kcolor,ls='--')
        a0.plot(delta[11][1:-1],delta[16][1:-1],'k',ls='--')
        a1.plot(delta[11][1:-1],delta[10][1:-1],color=wcolor,ls='--') #volume
        a2.plot(delta[11][1:-1],delta[22][1:-1],color=xcolor,ls='--')
        print (delta[16][-1]-delta[14][-1])
    #plt.savefig(title)
    plt.show()
    return dez

def f5d():
    print "\nFigure 5D"
    w=[[],[],[],[]]
    z=[[],[],[],[]]
    ZX=[-0.5,-1,-2,-3]
    ZT=[-0.5,-0.55,-0.6,-0.85,-0.9,-0.95,-1,-1.5,-1.75,-1.9,-1.94,-1.945,-2.0,-2.5,-2.75,-2.9,-2.94,-2.945,-3.0]
    plt.figure()
    for a in ZX:
        for b in ZT:
            dez=plm(xend=0,two=1,xt=10,xflux=1e-7,ztarget=b,Zx=a,moldelt=1,ratio=0.1)
            #minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
            w[abs(int(a))].append(np.log10(dez[10][-1]))
            z[abs(int(a))].append(dez[22][-1])
            print dez[14][100]-dez[14][-1]
            print (dez[10][-1])/dez[10][0]
        plt.plot(z[abs(int(a))],w[abs(int(a))],sym[abs(int(a))]+'o')
    #plt.savefig('f5d.eps')
    plt.show()
    return w,z