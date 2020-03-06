from detector_effects import B_f_obs, B_fbar_obs, Bbar_f_obs, Bbar_fbar_obs
from detector_effects import Amix_f_obs_fold, Amix_fbar_obs_fold
from detector_effects import eff_pow
import plot_utils
from plot_utils import plot_osc, plot_res, plot_amix, fold_times, pp
from plot_utils import bmix_leg, cpv_leg, res_leg
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

# Decay Time Resolution
# sigma_t - average resolution model
def plot_resolution(sigma_t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                    dm,dg,gs,r=0,delta=0,gamma=0,beta=0,
                    xmin=0,xmax=5,y_osc=2,y_mix=1,
                    name='plot.eps',save=False):

    fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(30,10))
    # phases
    delta_rad = delta*np.pi/180
    gamma_rad = gamma*np.pi/180
    beta_rad = beta/1000
    # Decay Rate Equations
    t = np.linspace(xmin,xmax,pp)
    B_f_obs_t = B_f_obs(t,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t)
    Bbar_f_obs_t = Bbar_f_obs(t,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t)
    B_fbar_obs_t = B_fbar_obs(t,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t)
    Bbar_fbar_obs_t = Bbar_fbar_obs(t,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t)
    plot_osc(ax1,t,B_f_obs_t,Bbar_f_obs_t,B_fbar_obs_t,Bbar_fbar_obs_t,xmin,xmax,ymax=y_osc,leghead=bmix_leg(dm,dg,gs))
    plot_res(ax2,sigma_t,leghead=res_leg(sigma_t))
    # Mixing Asymmetry
    xmin_mix = max(np.power(b_acc,1./n_acc)/a_acc,cutoff_acc)
    t_fold = fold_times(xmin_mix,xmax,dm)
    Amix_f_obs_t_fold = Amix_f_obs_fold(t_fold,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t)
    Amix_fbar_obs_t_fold = Amix_fbar_obs_fold(t_fold,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t)
    t_osc = np.linspace(0,2*np.pi/dm,pp)
    plot_amix(ax3,t_osc,Amix_f_obs_t_fold,Amix_fbar_obs_t_fold,0,2*np.pi/dm,
              title='Folded Asymmetries',xtitle=r't modulo $2\pi/\Delta m_{s}$ [ps]',
              xtitle_pos=[0.7,-0.07],ymin=-y_mix,ymax=y_mix,leghead=cpv_leg(r,delta,gamma,beta))
    # Plot
    fig.tight_layout()
    if save: fig.savefig(name)
    plt.show()
