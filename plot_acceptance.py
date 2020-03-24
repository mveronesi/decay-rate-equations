import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from decay_rate_compact import P_t, Amix_qf, Afold_qf
from detector_effects import eff_pow
import plot_utils
from plot_utils import plot_acc, plot_amix, fold_times, pp
from plot_utils import bmix_leg, cpv_leg, acc_leg

# Decay Time Acceptance
# a -
# n -
# b -
# beta -
# cutoff -
def plot_acceptance(a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                    dm,dg,gs,r=0,delta=0,gamma=0,beta=0,k=1,
                    xmin=0,xmax=5,y_osc=2,y_mix=1,
                    name='plot.eps',save=False,
                    fold_amix=True):

    fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(30,10))
    # phases
    delta_rad = delta*np.pi/180
    gamma_rad = gamma*np.pi/180
    beta_rad = beta/1000
    # Decay Rate Equations
    t = np.linspace(xmin,xmax,pp)
    eff_acc = eff_pow(t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc)
    B_f_obs_t = eff_acc*P_t(t=t,qt=1,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,k=k)
    Bbar_f_obs_t = eff_acc*P_t(t=t,qt=-1,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,k=k)
    B_fbar_obs_t = eff_acc*P_t(t=t,qt=1,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,k=k)
    Bbar_fbar_obs_t = eff_acc*P_t(t=t,qt=-1,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,k=k)
    tot_dec = B_f_obs_t + Bbar_f_obs_t + B_fbar_obs_t + Bbar_fbar_obs_t
    plot_acc(ax1,t,tot_dec,xmin,xmax,ymin=0,ymax=y_osc,leghead='')
    plot_acc(ax2,t,eff_acc,xmin,xmax,ymin=0,ymax=1,title='Acceptance',ytitle=r'$\varepsilon(t)$',
                leghead=acc_leg(a_acc,n_acc,b_acc,beta_acc,cutoff_acc),legcoo='lower right')
    # Mixing Asymmetry
    xmin_mix = max(np.power(b_acc,1./n_acc)/a_acc,cutoff_acc)
    t_fold = fold_times(xmin_mix,xmax,dm)
    if fold_amix and (len(t_fold)>1):
        Amix_f_t_fold = Afold_qf(t=t_fold,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,k=k)
        Amix_fbar_t_fold = Afold_qf(t=t_fold,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,k=k)
        t_osc = np.linspace(0,2*np.pi/dm,pp)
        plot_amix(ax3,t_osc,Amix_f_t_fold,Amix_fbar_t_fold,0,2*np.pi/dm,
                  title='Folded Asymmetries',xtitle=r't modulo $2\pi/\Delta m_{s}$ [ps]',
                  xtitle_pos=[0.7,-0.07],ymin=-y_mix,ymax=y_mix)
    else:
        Amix_f_t = Amix_qf(t=t,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,k=k)
        Amix_fbar_t = Amix_qf(t=t,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,k=k)
        plot_amix(ax3,t,Amix_f_t,Amix_fbar_t,xmin,xmax,ymin=-y_mix,ymax=y_mix)

    # Plot
    fig.tight_layout()
    if save: fig.savefig(name)
    plt.show()
