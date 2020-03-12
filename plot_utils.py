import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import mplhep as hep
plt.style.use(hep.style.ATLAS)
plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rc('text', usetex=True)

import decay_rate
from decay_rate import C_f, C_fbar, S_f, S_fbar, D_f, D_fbar
from decay_rate import B_f, B_fbar, Bbar_f, Bbar_fbar
from decay_rate import Amix_f, Amix_fbar
from detector_effects import eff_pow

lw = 2 # line width
fs = 40 # font size
pp = 1000 # points plot
fl = 30 # legend title font size

# legend titles
def bmix_leg(dm,dg,gs):
    leg = r'$\Delta m_{s}$ = '+'{:.3f}'.format(dm) + r' ps$^{-1}$' + '\n'
    leg += r'$\Delta \Gamma_{s}$ = '+'{:.3f}'.format(dg) + r' ps$^{-1}$' + '\n'
    leg += r'$\Gamma_{s}$ = '+'{:.3f}'.format(gs) + r' ps$^{-1}$'
    return leg

def cpv_leg(r,delta,gamma,beta):
    leg = r'$r$ = '+'{:.1f}'.format(r) + '\n'
    leg += r'$\delta$ = '+'{:.0f}'.format(delta) + r'$^{\circ}$' + '\n'
    leg += r'$\gamma$ = '+'{:.0f}'.format(gamma) + r'$^{\circ}$' + '\n'
    leg += r'$\beta_{s}$ = '+'{:.0f}'.format(beta) + r' mrad'
    return leg

def acc_leg(a,n,b,beta,cut):
    leg = r'$a$ = '+'{:.1f}'.format(a) + '\n'
    leg += r'$n$ = '+'{:.1f}'.format(n) + '\n'
    leg += r'$b$ = '+'{:.2f}'.format(b) + '\n'
    leg += r'$\beta$ = '+'{:.2f}'.format(beta) + '\n'
    leg += r'$t_{cut}$ = '+'{:.1f}'.format(cut)
    return leg

def res_leg(sigma_t):
    leg = r'$\sigma_{t}$ = ' + '{:.3f} ps'.format(sigma_t)
    return leg

def tag_leg(omega,d_omega,eff,d_eff):
    leg = r'$\omega_{tag}$ = ' + '{:.2f}'.format(omega) + '\n'
    leg += r'$\Delta\omega_{tag}$ = ' + '{:.2f}'.format(d_omega) + '\n'
    leg += r'$\varepsilon_{tag}$ = ' + '{:.2f}'.format(eff) + '\n'
    leg += r'$\Delta\varepsilon_{tag}$ = ' + '{:.2f}'.format(d_eff)
    return leg

def asymm_leg(a_prod,a_det):
    leg = r'$a_{prod}$ = ' + '{:.2f}'.format(a_prod) + '\n'
    leg += r'$a_{det}$ = ' + '{:.2f}'.format(a_det)
    return leg

# plot axes
def set_ax_labels(ax,xtitle,ytitle,title,xpos=[0.95, -0.070],ypos=[-0.08, 0.9]):
    ax.set_xlabel(xtitle,fontsize = fs)
    ax.set_ylabel(ytitle,fontsize = fs)
    ax.tick_params(labelsize = fs)
    ax.set_title(title,fontsize = fs)
    ax.xaxis.set_label_coords(xpos[0],xpos[1])
    ax.yaxis.set_label_coords(ypos[0],ypos[1])

def set_ax_lim(ax,xlim=False,ylim=False):
    if xlim: ax.set_xlim(xlim[0],xlim[1])
    if ylim: ax.set_ylim(ylim[0],ylim[1])
    ax.locator_params(nbins=5)

def angles(xmin,xmax,slope):
    xcoo = np.linspace(xmin,xmax,pp)
    ycoo = slope*xcoo
    return [xcoo,ycoo]

# Oscillation Plot
def plot_osc(ax,t,B_f_t,Bbar_f_t,B_fbar_t,Bbar_fbar_t,xmin=0,xmax=5,ymin=0,ymax=2,title='Decay Rate',leghead=''):
    ax.plot(t,B_f_t,color='blue',linewidth=lw,label=r'$B\to f$')
    ax.plot(t,Bbar_f_t,color='red',linewidth=lw,label=r'$\bar{B}\to f$')
    ax.plot(t,Bbar_fbar_t,color='blue',linewidth=lw,linestyle='--',label=r'$\bar{B}\to \bar{f}$')
    ax.plot(t,B_fbar_t,color='red',linewidth=lw,linestyle='--',label=r'$B\to \bar{f}$')
    set_ax_labels(ax,'t [ps]',r'$d\Gamma/dt$',title)
    set_ax_lim(ax,[xmin,xmax],[ymin,ymax])
    legend = ax.legend(loc='upper right',fontsize=fs,fancybox=True,title=leghead)
    plt.setp(legend.get_title(),fontsize=fl)

def plot_untag(ax,t,Untag_f_t,Untag_fbar_t,xmin=0,xmax=5,ymin=0,ymax=2,title='Decay Rate',leghead=''):
    ax.plot(t,Untag_f_t,color='black',linewidth=lw,label=r'$f$')
    ax.plot(t,Untag_fbar_t,color='red',linestyle='--',linewidth=lw,label=r'$\overline{f}$')
    set_ax_labels(ax,'t [ps]',r'$d\Gamma/dt$',title)
    set_ax_lim(ax,[xmin,xmax],[ymin,ymax])
    legend = ax.legend(loc='upper right',fontsize=fs,fancybox=True,title=leghead)
    plt.setp(legend.get_title(),fontsize=fl)

# Acceptance Plot
def plot_acc(ax,t,acc,xmin=0,xmax=0,ymin=0,ymax=2,
            title='Decay Rate',ytitle=r'$d\Gamma/dt$',
            leghead='',legcoo='upper right'):
    ax.plot(t,acc,color='blue',linewidth=lw,label=r'$B+\overline{B}$')
    set_ax_labels(ax,'t [ps]',ytitle,title)
    set_ax_lim(ax,[xmin,xmax],[ymin,ymax])
    legend = ax.legend(loc=legcoo,fontsize=fs,fancybox=True,title=leghead)
    plt.setp(legend.get_title(),fontsize=fl)

# Time Resolution Plot
def plot_res(ax,sigma_t,xmin=-0.35,xmax=0.35,title='Time Resolution',ytitle='',leghead=''):
    dt = np.linspace(xmin,xmax,pp)
    label = r'$G(\Delta t;\sigma_t)$'
    if sigma_t>0:
        # res_func = 1./(np.sqrt(2*np.pi)*sigma_t)*np.exp(-0.5*np.power(dt/sigma_t,2))
        res_func = np.exp(-0.5*np.power(dt/sigma_t,2))
        ax.plot(dt,res_func,color='blue',linewidth=lw,label=label)
    else:
        line = np.linspace(0,10,pp)
        zero = np.zeros(pp)
        ax.plot(zero,line,color='blue',linewidth=lw,label=label)
    set_ax_labels(ax,r'$\Delta t$ [ps]',ytitle,title,xpos=[0.9, -0.070])
    set_ax_lim(ax,[xmin,xmax],[0,1.1])
    # ax.set_ylim(bottom=0)
    legend = ax.legend(loc='upper right',fontsize=30,fancybox=True,title=leghead)
    plt.setp(legend.get_title(),fontsize=fl)

# Mixing Asymmetry
def plot_amix(ax,t,Amix_f_t,Amix_fbar_t,xmin,xmax,ymin=-1,ymax=1,
              title='Mixing Asymmetries',xtitle='t [ps]',xtitle_pos = [0.95, -0.070],
              leghead=''):
    ax.plot(t,Amix_f_t,linewidth=lw,color='black',label=r'$f$')
    ax.plot(t,Amix_fbar_t,linewidth=lw,color='red',linestyle='--',label=r'$\bar{f}$')
    set_ax_labels(ax,xtitle,r'$A_{\mathrm{mix}}$',title,xtitle_pos)
    set_ax_lim(ax,[xmin,xmax],[ymin,ymax])
    legend = ax.legend(loc='lower left',fontsize=fs,fancybox=True,title=leghead)
    plt.setp(legend.get_title(),fontsize=fl)

def fold_times(xmin,xmax,dm):
    first_osc = math.ceil(xmin/(2*np.pi/dm))
    last_osc = int(xmax/(2*np.pi/dm))
    nosc = last_osc-first_osc
    if nosc>0:
        xmin_fold = first_osc*2*np.pi/dm
        xmax_fold = last_osc*2*np.pi/dm
        t = np.linspace(xmin_fold,xmax_fold,pp*nosc)
        return t
    else:
        return [0]

# Constraints On Gamma
def plot_gamma(ax,r,delta,gamma,D_f_val,S_f_val,D_fbar_val,S_fbar_val):
    line = np.linspace(-1,1,pp)
    zero = np.zeros(pp)
    ax.plot(zero,line,color='black')
    ax.plot(line,zero,color='black')
    ## Interference
    if r!=0:
        circle = plt.Circle((0, 0), 2*r/(1+r**2), color='cyan',fill=False,label=r'$\sqrt{1-C_f^2}$',linewidth=lw)
        ax.add_patch(circle)
    ## Gamma
    ax.plot(angles(0 if gamma<=np.pi/2 else -1,1 if gamma<=np.pi/2 else 0,np.tan(gamma))[0],
            angles(0 if gamma<=np.pi/2 else -1,1 if gamma<=np.pi/2 else 0,np.tan(gamma))[1],
            color='orange',label=r'$\gamma-2\beta_s$',linewidth=lw)
    ## Delta
    ax.plot(angles(0 if (delta<=np.pi/2 or delta>3*np.pi/2) else -1,1 if (delta<=np.pi/2 or delta>3*np.pi/2) else 0,np.tan(delta))[0],
            angles(0 if (delta<=np.pi/2 or delta>3*np.pi/2) else -1,1 if (delta<=np.pi/2 or delta>3*np.pi/2) else 0,np.tan(delta))[1],
            color='blue',label=r'$\delta$',linewidth=lw)
    ## CP Parameters
    if D_f_val!=0:
        ax.plot(angles(0 if -D_f_val>0 else -D_f_val,-D_f_val if -D_f_val>0 else 0,S_f_val/D_f_val)[0],
                angles(0 if -D_f_val>0 else -D_f_val,-D_f_val if -D_f_val>0 else 0,S_f_val/D_f_val)[1],
                color='green',linewidth=lw)
    if D_fbar_val!=0:
        ax.plot(angles(0 if -D_fbar_val>0 else -D_fbar_val,-D_fbar_val if -D_fbar_val>0 else 0,S_fbar_val/D_fbar_val)[0],
                angles(0 if -D_fbar_val>0 else -D_fbar_val,-D_fbar_val if -D_fbar_val>0 else 0,S_fbar_val/D_fbar_val)[1],
                color='brown',linewidth=lw)
    ax.plot([-D_f_val],[-S_f_val],marker='o',color='green',markersize=15,linewidth=0,label=r'$(-D_f,-S_f)$')
    ax.plot([-D_fbar_val],[-S_fbar_val],marker='o',color='brown',markersize=15,linewidth=0,label=r'$(-D_{\bar{f}},-S_{\bar{f}})$')
    legend = ax.legend(loc='lower left',fontsize=fs,fancybox=True)
    set_ax_lim(ax,[-1,1],[-1,1])
    set_ax_labels(ax,r'$2r/(1+r^2)\cos((\gamma-2\beta_s) \pm \delta)$',
                     r'$2r/(1+r^2)\sin((\gamma-2\beta_s) \pm \delta)$',
                     r'Constraints on $\gamma$',
                     xpos=[0.7, -0.070],ypos=[-0.08, 0.7])
