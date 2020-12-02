import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import decay_rate_cpt
from decay_rate_cpt import P_t_fs, Acp
import plot_utils
from plot_utils import plot_osc, plot_acp, fold_times, pp

# Theoretical Decay Rate
def plot_decrate_cpt(dm=17.757,dg=0.085,gs=0.664,
                     afs=0,rez=0,imz=0,
                     xmin=0,xmax=5,y_osc=2,y_mix=1,
                     name='plot_cpv.eps',save=False,
                     b_f=True,
                     bbar_f=True,
                     bbar_fbar=True,
                     b_fbar=True,
                     asymm='MIX'):

    # fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(30,10))
    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(20,10))
    # Decay Rate Equations
    t = np.linspace(xmin,xmax,pp)
    B_f_t = P_t_fs(t=t,qt=1,qf=1,dm=dm,dg=dg,gs=gs,rez=rez,imz=imz,afs=afs) if b_f else False
    Bbar_f_t = P_t_fs(t=t,qt=-1,qf=1,dm=dm,dg=dg,gs=gs,rez=rez,imz=imz,afs=afs) if bbar_f else False
    B_fbar_t = P_t_fs(t=t,qt=1,qf=-1,dm=dm,dg=dg,gs=gs,rez=rez,imz=imz,afs=afs) if b_fbar else False
    Bbar_fbar_t = P_t_fs(t=t,qt=-1,qf=-1,dm=dm,dg=dg,gs=gs,rez=rez,imz=imz,afs=afs) if bbar_fbar else False
    plot_osc(ax1,t,B_f_t,Bbar_f_t,B_fbar_t,Bbar_fbar_t,xmin,xmax,ymax=y_osc,title='')
    # CP Asymmetries
    if asymm == 'CPT':
        Amix = Acp(first=B_f_t,second=Bbar_fbar_t)
        acol = 'blue'
    elif asymm == 'CP':
        Amix = Acp(first=Bbar_f_t,second=B_fbar_t)
        acol = 'red'
    elif asymm == 'MIX':
        Amix = Acp(first=B_f_t,second=Bbar_f_t)
        acol = 'black'
    plot_acp(ax=ax2,t=t,Acp=Amix,xmin=xmin,xmax=xmax,ymin=-y_mix,ymax=y_mix,lc=acol,ytitle=r'$A_{'+asymm+r'}$')
    fig.tight_layout()
    if save: fig.savefig(name)
    plt.show()
