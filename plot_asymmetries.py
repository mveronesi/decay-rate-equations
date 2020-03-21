import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from decay_rate_compact import P_t, Amix_qf, Afold_qf
from detector_effects import eff_pow
import plot_utils
from plot_utils import plot_osc, plot_untag, plot_amix, fold_times, pp
from plot_utils import bmix_leg, cpv_leg, asymm_leg

def plot_asymmetries(a_prod,a_det,
                     omega,d_omega,eff,d_eff,
                     sigma_t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                     dm,dg,gs,r=0,delta=0,gamma=0,beta=0,
                     xmin=0,xmax=5,y_tag=0.4,y_untag=0.2,y_mix=1,
                     name='plot.eps',save=False,
                     k_acc=True,k_res=True,k_tag=True,
                     b_f=True,
                     bbar_f=True,
                     bbar_fbar=True,
                     b_fbar=True,
                     fold_amix=True):
    fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(30,10))
    # phases
    delta_rad = delta*np.pi/180
    gamma_rad = gamma*np.pi/180
    beta_rad = beta/1000
    # Decay Rate Equations
    t = np.linspace(xmin,xmax,pp)
    eff_acc = eff_pow(t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc) if k_acc else 1
    sigma_t = sigma_t if k_res else 0
    (omega, d_omega, eff, d_eff) = (omega, d_omega, eff, d_eff) if k_tag else (0,0,1,0)
    B_f_t = eff_acc*P_t(t=t,qt=1,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,
                                sigma_t=sigma_t,eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,
                                a_prod=a_prod,a_det=a_det) if b_f else False
    Bbar_f_t = eff_acc*P_t(t=t,qt=-1,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,
                                sigma_t=sigma_t,eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,
                                a_prod=a_prod,a_det=a_det) if bbar_f else False
    B_fbar_t = eff_acc*P_t(t=t,qt=1,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,
                                sigma_t=sigma_t,eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,
                                a_prod=a_prod,a_det=a_det) if b_fbar else False
    Bbar_fbar_t = eff_acc*P_t(t=t,qt=-1,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,
                                sigma_t=sigma_t,eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,
                                a_prod=a_prod,a_det=a_det) if bbar_fbar else False
    plot_osc(ax1,t,B_f_t,Bbar_f_t,B_fbar_t,Bbar_fbar_t,xmin,xmax,ymax=y_tag,title='Tagged Decay Rate',leghead='')
    Untag_f_t = eff_acc*P_t(t=t,qt=0,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,
                                sigma_t=sigma_t,eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,
                                a_prod=a_prod,a_det=a_det) if (b_f or bbar_f) else False
    Untag_fbar_t = eff_acc*P_t(t=t,qt=0,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,
                                sigma_t=sigma_t,eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,
                                a_prod=a_prod,a_det=a_det) if (b_fbar or bbar_fbar) else False
    plot_untag(ax2,t,Untag_f_t,Untag_fbar_t,xmin,xmax,ymax=y_untag,title='Untagged Decay Rate',leghead=asymm_leg(a_prod,a_det))
    # Mixing Asymmetry
    xmin_mix = max(np.power(b_acc,1./n_acc)/a_acc,cutoff_acc)
    t_fold = fold_times(xmin_mix,xmax,dm)
    if fold_amix and (len(t_fold)>1):
        Amix_f_t_fold = Afold_qf(t=t_fold,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,
                                        sigma_t=sigma_t,eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,
                                        a_prod=a_prod,a_det=a_det) if (b_f or bbar_f) else False
        Amix_fbar_t_fold = Afold_qf(t=t_fold,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,
                                        sigma_t=sigma_t,eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,
                                        a_prod=a_prod,a_det=a_det) if (b_fbar or bbar_fbar) else False
        t_osc = np.linspace(0,2*np.pi/dm,pp)
        plot_amix(ax3,t_osc,Amix_f_t_fold,Amix_fbar_t_fold,0,2*np.pi/dm,
                  title='Folded Asymmetries',xtitle=r't modulo $2\pi/\Delta m_{s}$ [ps]',
                  xtitle_pos=[0.7,-0.07],ymin=-y_mix,ymax=y_mix)
    else:
        Amix_f_t = Amix_qf(t=t,qf=1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,sigma_t=sigma_t,
                            eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,a_prod=a_prod,a_det=a_det) if (b_f or bbar_f) else False
        Amix_fbar_t = Amix_qf(t=t,qf=-1,dm=dm,dg=dg,gs=gs,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,sigma_t=sigma_t,
                            eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,a_prod=a_prod,a_det=a_det) if (b_fbar or bbar_fbar) else False
        plot_amix(ax3,t,Amix_f_t,Amix_fbar_t,xmin,xmax,ymin=-y_mix,ymax=y_mix)
    # Plot
    fig.tight_layout()
    if save: fig.savefig(name)
    plt.show()
