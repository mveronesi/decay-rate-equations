import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import decay_rate_compact
from decay_rate_compact import P_t, Amix_qf, Afold_qf
from decay_rate_compact import C_qf, S_qf, A_qf
import plot_utils
from plot_utils import plot_osc, plot_amix, plot_gamma, fold_times, pp


# Theoretical Decay Rate
# delta - degrees
# gamma - degrees
# beta - mrad
def plot_decrate(dm=17.757,dg=0.085,gs=0.664,
                 r=0,delta=0,gamma=0,beta=0,
                 xmin=0,xmax=5,y_osc=2,y_mix=1,
                 name='plot.eps',save=False,
                 qt_plus=True,qt_min=True,
                 qf_plus=True,qf_min=True,
                 fold_amix=True):
    # ymin=0
    # ymax=2
    fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(30,10))
    # phases
    delta_rad = delta*np.pi/180
    gamma_rad = gamma*np.pi/180
    beta_rad = beta/1000
    # CP coefficients
    C_f_val = C_qf(r,+1)
    C_fbar_val = C_qf(r,-1)
    A_f_val = A_qf(r,delta_rad,gamma_rad,beta_rad,+1)
    A_fbar_val = A_qf(r,delta_rad,gamma_rad,beta_rad,-1)
    S_f_val = S_qf(r,delta_rad,gamma_rad,beta_rad,+1)
    S_fbar_val = S_qf(r,delta_rad,gamma_rad,beta_rad,-1)
    # Decay Rate Equations
    t = np.linspace(xmin,xmax,pp)
    B_f_t = P_t(t=t,qt=1,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad) if (qt_plus and qf_plus) else False
    Bbar_f_t = P_t(t=t,qt=-1,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad) if (qt_min and qf_plus) else False
    B_fbar_t = P_t(t=t,qt=1,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad) if (qt_plus and qf_min) else False
    Bbar_fbar_t = P_t(t=t,qt=-1,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad) if (qt_min and qf_min) else False
    plot_osc(ax1,t,B_f_t,Bbar_f_t,B_fbar_t,Bbar_fbar_t,xmin,xmax,ymax=y_osc)
    # Mixing Asymmetries
    t_fold = fold_times(xmin,xmax,dm)
    if fold_amix and (len(t_fold)>1):
        Amix_f_t_fold = Afold_qf(t=t_fold,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad) if qf_plus else False
        Amix_fbar_t_fold = Afold_qf(t=t_fold,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad) if qf_min else False
        t_osc = np.linspace(0,2*np.pi/dm,pp)
        plot_amix(ax2,t_osc,Amix_f_t_fold,Amix_fbar_t_fold,0,2*np.pi/dm,
                  title='Folded Asymmetries',xtitle=r't modulo $2\pi/\Delta m_{s}$ [ps]',
                  xtitle_pos=[0.7,-0.07],ymin=-y_mix,ymax=y_mix)
    else:
        Amix_f_t = Amix_qf(t=t,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad) if qf_plus else False
        Amix_fbar_t = Amix_qf(t=t,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad) if qf_min else False
        plot_amix(ax2,t,Amix_f_t,Amix_fbar_t,xmin,xmax,ymin=-y_mix,ymax=y_mix)
    # Constraints on Gamma
    plot_gamma(ax3,r,delta_rad,gamma_rad,A_f_val,S_f_val,A_fbar_val,S_fbar_val)

    # Plot
    fig.tight_layout()
    if save: fig.savefig(name)
    plt.show()
