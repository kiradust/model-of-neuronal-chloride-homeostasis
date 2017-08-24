from plm_singlecomp_withkcc2 import plm, zp, F
from plotting import minifig, minithreefig, minifigtwoaxes,twoaxes,xcolor,clcolor,wcolor,kcolor
import matplotlib.pyplot as plt
from pylab import rcParams
from fig4 import f4e
import numpy as np
rcParams['figure.figsize'] = 8,8
from matplotlib import gridspec

sym=['-k',':k','--k','-.k']

def f5a():
    dez=plm(dz=2.5e-7,two=1,xt=180,tt=540,ztarget=-1)
    minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k',yl=[[-100,-70],[1.8e-12,3.3e-12],[-0.95,-0.8]])
    plt.savefig('f5a.eps')
    print (dez[16][-1]-dez[14][-1])
    print (dez[16][8000]-dez[14][8000])
    return

def f5b(moldelt=0):
    XF=[-1.20,-1.15,-1.1,-1.05,-1.0,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.501]
    XFp=[]
    cl=[]
    vm=[]
    k=[]
    w=[]
    z=[]
    x=[]
    df=[]
    pi, zi, zee, newx = f4e(moldelt=moldelt)
    
    for i in range(len(zi)):
        for b in XF:
            if zi[i]/100.0==b:
                XFp.append(pi[i])
    for a in range(len(XFp)):
        #dez=plm(two=1,xt=10,dz=3e-6,Zx=a,tt=5000)
        print XFp[a]
        dez=plm(z=XF[a],tt=500,p=(10**(XFp[a]))/F,graph=1,osmofix=True,xinit=0)
        cl.append(dez[7])
        k.append(dez[6])
        vm.append(dez[9])
        w.append(dez[10][-1])
        z.append(dez[22][-1]*100)
        x.append(dez[3])
        df.append(dez[9]-dez[7])
    pi, zi, zee, newx = f4e(moldelt=moldelt)
    
    plt.figure()
    gs = gridspec.GridSpec(3,1,height_ratios=[1.5,0.5,0.5])
    a0=plt.subplot(gs[0])
    a0.plot(zi,zee[3],color=clcolor)
    a0.plot(z,cl,'bo')
    a0.plot(z,k,'go')
    a0.plot(z,vm,'ko')
    a0.plot(zi,zee[2],color=kcolor)
    a0.plot(zi,zee[5],'k')
    a0.set_ylim([-110,-60])
    #a1=plt.subplot(gs[1])
    #a1.plot(zi,zee[11],color='k') #vol
    #a1.plot(z,w,'ko')
    #a1.set_ylim([0,5e-12])
    a2=plt.subplot(gs[1])
    a2.plot(zi,newx,color='m') #conc X
    a2.plot(z,x,'mo')
    a2.set_ylim([0.11,0.2])
    a3=plt.subplot(gs[2])
    a3.plot(z,df)
    a3.set_ylim([0,21])
    plt.savefig('f5b_2.eps')
    plt.show()
    return

def f5c(ratio=0.8,md=1e-12):
    dez=plm(gx=1e-8,xt=180,two=1,xend=0,moldelt=md,xflux=0.3*1e-6,ztarget=-0.9,tt=540)
    minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k',yl=[[-100,-70],[1.8e-12,3.3e-12],[-0.95,-0.8]])
    plt.savefig('f5c.eps')
    print (dez[16][-1]-dez[14][-1])
    print (dez[16][8000]-dez[14][8000])
    return

def f5d():
    w=[[],[],[],[]]
    z=[[],[],[],[]]
    ZX=[-0.5,-1,-2,-3]
    ZT=[-0.5,-0.55,-0.6,-0.85,-0.9,-0.95,-1,-1.5,-1.75,-1.9,-1.94,-1.945,-2.0,-2.5,-2.75,-2.9,-2.94,-2.945,-3.0]
    plt.figure()
    for a in ZX:
        for b in ZT:
            dez=plm(xend=0,two=1,xt=10,xflux=1e-6,ztarget=b,Zx=a,moldelt=1,ratio=0.1)
            #minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
            w[abs(int(a))].append(np.log10(dez[10][-1]))
            z[abs(int(a))].append(dez[22][-1])
            print dez[14][-1]
            print (dez[10][-1])/dez[10][0]
        plt.plot(z[abs(int(a))],w[abs(int(a))],sym[abs(int(a))]+'o')
    plt.savefig('f5e.eps')
    plt.show()
    return w,z