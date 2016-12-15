from plm_singlecomp_withkcc2 import plm, zp
from plotting import minifig, minithreefig, minifigtwoaxes,twoaxes,xcolor,clcolor,wcolor,kcolor
import matplotlib.pyplot as plt
from pylab import rcParams
import numpy as np
from fig4 import f4e
rcParams['figure.figsize'] = 8,8

def f5a():
    dez=plm(dz=1e-5,two=1,xt=10)
    minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
    return

def f5b(moldelt=0):
    XF=[-1.2,-1.1,-1.0,-0.9,-0.85,-0.8,-0.7,-0.6,-0.55]
    cl=[]
    vm=[]
    k=[]
    w=[]
    z=[]
    x=[]
    for a in XF:
        dez=plm(tt=200,z=a)
        #minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
        cl.append(dez[7])
        k.append(dez[6])
        vm.append(dez[9])
        w.append(dez[10][-1])
        z.append(dez[22][-1]*100)
        x.append(dez[3])
    ax0, ax1, ax2 = f4e(moldelt=moldelt)
    ax0.plot(z,cl,'bo')
    ax0.plot(z,k,'go')
    ax0.plot(z,vm,'ko')
    ax2.plot(z,x,'mo')
    ax1.plot(z,w,'ko')
    plt.show()
    return

def f5c(ratio=0.8,md=1e-13):
    dez=plm(gx=1e-8,xt=25,tt=200,two=1,xend=0,moldelt=md,ratio=ratio)
    minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
    return

def f5d():
    return
