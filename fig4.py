# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 17:29:22 2016

@author: Kira
@title: Figure 4
"""
from plm_singlecomp_withkcc2 import plm, zp, F
from plotting import minifig, minithreefig, minifigtwoaxes,twoaxes,xcolor,clcolor,wcolor,kcolor,nacolor
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pylab import rcParams
import numpy as np
from figure_1 import f1c
rcParams['figure.figsize'] = 8,8

sym=['-',':','--','-.']

def f4a(init_x=[40e-3,80e-3,120e-3,160e-3]):
    print "Figure 4A"
    
    plt.figure()
    for i in range(len(init_x)):
        endcl=plm(xinit=init_x[i],tt=1800, k_init=0,osmofix=True)
        plt.subplot(2,1,1)
        plt.plot(endcl[11][1:-1],endcl[20][1:-1],color=xcolor,linestyle=sym[i])
        plt.subplot(2,1,2)
        plt.plot(endcl[11][0:-1],endcl[10][0:-1],'k'+sym[i])
    #plt.savefig('f4a.eps')
    plt.show()
    
    return endcl[20], endcl
    
def f4b(init_x=range(25,586,40),new=0,l='-',title='f4b.eps',a=0,b=0):
    endv=[]
    endk=[]
    endcl=[]
    endw =[]
    
    print "\nFigure 4B"
    for i in init_x:
        end=plm(xinit=i*1e-3,k_init=0,tt=2000,osmofix=False,neww=new,ls=l)
        endv.append(end[9])
        endk.append(end[6])
        endcl.append(end[7])
        endw.append(end[21])
        print end[9]
        print end[6]
        print end[7]
    
    plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 1]) 
    if a==0:
        a=plt.subplot(gs[0])
    a.plot(init_x,endcl,color=clcolor,linestyle=l)
    a.plot(init_x,endcl,'bo',linestyle=l)
    a.plot(init_x,endk,color=kcolor,linestyle=l)
    a.plot(init_x,endk,'ro',linestyle=l)
    a.plot(init_x,endv,'k',linestyle=l)
    a.plot(init_x,endv,'ko',linestyle=l)
    a.set_ylim([-100,-70])
    if b==0:
        b=plt.subplot(gs[1])
    b.plot(init_x,endw,color=wcolor,linestyle=l)
    b.plot(init_x,endw,'ko',linestyle=l)
    b.set_ylim([0,8e-12])
    #plt.savefig(title)
    plt.show()
    
    return a,b
    
def f4c(gX=1e-8,tt=3600,xt=360,xend=420,xflux=4e-7,new=0,title='f4c.eps',ham=0): #doubles as f6c when new!=0
    dex=plm(gx=gX,xt=xt,tt=tt,xflux=xflux,xend=xend,graph=0,hamada=0)
    delta=[]
    delta1=[]
    
    if new==0:
        print "\nFigure 4C"
        ax0,ax1,ax2=minithreefig([dex[11][1:-1],dex[14][1:-1],dex[13][1:-1],dex[16][1:-1],dex[10][1:-1],dex[20][1:-1]],xcolor,yl=[[-100,-70],[1.9e-12,2.5e-12],[154,157]])
    
    else:
        print "\nFigure 6C"
        ax0,ax1,ax2=minithreefig([dex[11][1:-1],dex[14][1:-1],dex[13][1:-1],dex[16][1:-1],dex[18][1:-1],dex[10][1:-1]],'k',yl=[[-100,-70],[13,19],[1.8e-12,2.2e-12]])

        delta=plm(gx=gX,xt=xt,tt=tt,xflux=xflux,xend=xend,neww=3,graph=0,hamada=ham)
        ax0.plot(delta[11][1:-1],delta[14][1:-1],color=clcolor,linestyle='--')
        ax0.plot(delta[11][1:-1],delta[13][1:-1],color=kcolor,ls='--')
        ax0.plot(delta[11][1:-1],delta[16][1:-1],'k',ls='--')
        print (delta[16][-1]-delta[14][-1])
        print (delta[16][350]-delta[14][350])
        print len(delta[16])
        ax1.plot(delta[11][1:-1],delta[18][1:-1],color=nacolor,ls='--') #nai
        ax2.plot(delta[11][1:-1],delta[10][1:-1],color='k',ls='--') #volume
        
        delta1=plm(gx=gX,xt=xt,tt=tt,xflux=xflux,xend=xend,neww=5,graph=0,hamada=ham)
        if ham==1:
            delta2=plm(gx=gX,xt=xt,tt=tt,xflux=xflux,xend=xend,graph=0,hamada=ham)
            ax0.plot(delta2[11][1:-1],delta2[14][1:-1],color=clcolor,linestyle=':')
            ax0.plot(delta2[11][1:-1],delta2[13][1:-1],color=kcolor,ls=':')
            ax0.plot(delta2[11][1:-1],delta2[16][1:-1],'k',ls=':')
            ax1.plot(delta2[11][1:-1],delta2[18][1:-1],color=nacolor,ls=':') #nai
            ax2.plot(delta2[11][1:-1],delta2[10][1:-1],color='k',ls=':') #volume
            print (delta2[16][-1]-delta2[14][-1])
            print len(delta2[16])
        ax0.plot(delta1[11][1:-1],delta1[14][1:-1],color=clcolor,linestyle='-.')
        ax0.plot(delta1[11][1:-1],delta1[13][1:-1],color=kcolor,ls='-.')
        ax0.plot(delta1[11][1:-1],delta1[16][1:-1],'k',ls='-.')
        ax1.plot(delta1[11][1:-1],delta1[18][1:-1],color=nacolor,ls='-.') #nai
        ax2.plot(delta1[11][1:-1],delta1[10][1:-1],color='k',ls='-.') #volume
        print (delta1[16][-1]-delta1[14][-1])
        print len(delta1[16])
    plt.savefig(title)
    plt.show()
    return dex,delta,delta1
    
def f4d(f=2e-3,new=0,title='f4d.eps'):
    print "\nFigure 4D"
    
    #dxe=plm(gx=0,xt=1600,two=0,tt=3000,f4d=f,neww=0,graph=1,p=(10**(-3))/(F),pkcc=1e-3/F)
    dxe=plm(gx=0,xt=360,two=0,tt=1800,f4d=f,neww=0)
    
    a0,a1,a2=minithreefig([dxe[11][1:-1],dxe[14][1:-1],dxe[13][1:-1],dxe[16][1:-1],dxe[10][1:-1],dxe[23][1:-1]],xcolor,yl=[[-100,-70],[1.85e-12,2.0e-12],[0,0.09]])
    print (dxe[16][-1]-dxe[14][-1])
    #print (dxe[16][350]-dxe[14][350])
    
    if new==1:
        delta=plm(gx=0,xt=360,two=0,tt=1800,f4d=f,neww=1,graph=0)
        a0.plot(delta[11][1:-1],delta[14][1:-1],color=clcolor,linestyle='--')
        a0.plot(delta[11][1:-1],delta[13][1:-1],color=kcolor,ls='--')
        a0.plot(delta[11][1:-1],delta[16][1:-1],'k',ls='--')
        a1.plot(delta[11][1:-1],delta[10][1:-1],color=wcolor,ls='--') #volume
        a2.plot(delta[11][1:-1],delta[23][1:-1],color=xcolor,ls='--') #conc X
    plt.savefig(title)
    plt.show()
    return dxe

def f4e(Z=range(-120,-50),moldelt=1e-12): # extra function needed in figure 5
    molinit=plm(gx=1e-8,xt=25,tt=100,two=1,paratwo=True,moldelt=moldelt)
    zee=zp(Z=Z,molinit=molinit,moldelt=moldelt)
    newx=[]
    for i in zee[9]:
        newx.append(i)
    return zee[0],Z,zee,newx

# earlier ideas for potential supplementary figures

def sf4fa():
    XF=[5e-14,1e-13,2e-13]
    GX=[-10,-9,-8]
    col=['mo','bo','go']
    plt.figure()
    for a in XF:
        for g in GX:
            dez=plm(tt=2500,moldelt=a,two=1,gx=10**g,xt=25,xend=0,ratio=0.8)
            #minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
            plt.plot(g,dez[22][-1],col[XF.index(a)])
    plt.show()
    return
    
def sf4fb():
    XF=[5e-14,1e-13,2e-13]
    rat=[0.2,0.4,0.6,0.8]
    GX=[-10,-9,-8]
    col=['mo','bo','go','ro']
    plt.figure()
    for s in XF:
        for r in rat:
            dez=plm(tt=500,moldelt=s,two=1,gx=10**(-8),xt=25,xend=0,ratio=r)
            #minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
            plt.plot(r,dez[22][-1],col[XF.index(s)])
            #plt.ylim((-0.83,-0.85))
    plt.show()
    return

def sf4c(GX=[5e-10,1e-9,5e-9,7e-9,1e-8,2e-8],tt=600,xt=25,ratio=0.98,xend=0):
    deltecl=[]
    maxdeltecl=[]
    deltw=[]
    deltx=[]
    for i in GX:
        dex=plm(gx=i,xt=xt,tt=tt,ratio=ratio,xend=xend)
        deltecl.append(dex[14][-1]-dex[14][1])
        maxdeltecl.append(max(np.absolute((dex[14]-dex[14][1]))))
        deltw.append((dex[10][-1])/dex[10][1])
        deltx.append(max(np.absolute((dex[20]-dex[20][1]))))
        minithreefig([dex[11][1:-1],dex[14][1:-1],dex[13][1:-1],dex[16][1:-1],dex[10][1:-1],dex[20][1:-1]],xcolor)
    print maxdeltecl
    twoaxes(GX,deltecl,maxdeltecl,deltx,deltw)
    return