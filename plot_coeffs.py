from decay_rate import exp_dec
from decay_rate import C_f, C_fbar, S_f, S_fbar, D_f, D_fbar
from detector_effects import eff_pow, dil_res
from eff_coeffs import k_cosh, k_sinh, k_cos, k_sin
from plot_utils import plot_func, pp
from plot_utils import cpv_leg, C_f_leg, S_f_leg, D_f_leg
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

def plot_coeffs(qt,qf,
                dm=17.757,dg=0.085,gs=0.664,
                r=0.4,delta=10,gamma=60,beta=0,
                a_acc=1.5,n_acc=1.5,b_acc=0.05,beta_acc=0.03,cutoff_acc=0.2,
                sigma_t=0.035,
                omega_tag=0.35,d_omega_tag=0,eff_tag=0.8,d_eff_tag=0,
                a_prod_asym=0,a_det_asym=0,
                xmin=0,xmax=5,
                y_cosh=0.5,y_sinh=0.01,y_cos=0.05,y_sin=0.05,
                name='plot.eps',save=False,
                k_dec=True,
                k_acc=True,
                k_res=True,
                k_tag=True,
                k_prod=True,
                k_det=True):
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(25,25))
    # phases
    delta_rad = delta*np.pi/180
    gamma_rad = gamma*np.pi/180
    beta_rad = beta/1000
    # detector effects
    t = np.linspace(xmin,xmax,pp)
    dec = exp_dec(t,gs) if k_dec else 1
    acc = eff_pow(t,a=a_acc,n=n_acc,b=b_acc,beta=beta_acc,cutoff=cutoff_acc) if k_acc else 1
    res = dil_res(sigma_t,dm) if k_res else 1
    (omega, d_omega, eff, d_eff) = (omega_tag, d_omega_tag, eff_tag, d_eff_tag) if k_tag else (0,0,1,0)
    a_prod = a_prod_asym if k_prod else 0
    a_det = a_det_asym if k_det else 0
    # trigonometric functions
    cosh_t = np.cosh(0.5*dg*t)
    sinh_t = np.sinh(0.5*dg*t)
    cos_t = np.cos(0.5*dm*t)
    sin_t = np.sin(0.5*dm*t)
    # labels
    title_sinh = r'$(D_f,D_{\overline{f}})$ = '+'({:.2f},{:.2f})'.format(D_f(r,delta,gamma,beta),D_fbar(r,delta,gamma,beta))
    title_cos = r'$(C_f,C_{\overline{f}})$ = '+'({:.2f},{:.2f})'.format(C_f(r),C_fbar(r))
    title_sin = r'$(S_f,S_{\overline{f}})$ = '+'({:.2f},{:.2f})'.format(S_f(r,delta,gamma,beta),S_fbar(r,delta,gamma,beta))
    title_cpv = r'$r$ = '+'{:.1f}'.format(r)
    title_cpv += r', $\delta$ = '+'{:.0f}'.format(delta) + r'$^{\circ}$'
    title_cpv += r', $\gamma$ = '+'{:.0f}'.format(gamma) + r'$^{\circ}$'
    title_cpv += r', $\beta_{s}$ = '+'{:.0f}'.format(beta) + r' mrad'
    q_list = [(1,1),(-1,1),(-1,-1),(1,-1),(0,1),(0,-1)]
    col_list = ['blue','red','blue','red','black','black']
    style_list = ['-','-','--','--','-','--']
    for i in range(6):
        (qt,qf) = q_list[i]
        col=col_list[i]
        style=style_list[i]
        # Effective coefficients
        coeff_cosh = k_cosh(qt,qf,omega,d_omega,eff,d_eff,a_prod,a_det)
        coeff_sinh = k_sinh(qt,qf,r,delta_rad,gamma_rad,beta_rad,omega,d_omega,eff,d_eff,a_prod,a_det)
        coeff_cos = k_cos(qt,qf,r,omega,d_omega,eff,d_eff,a_prod,a_det)
        coeff_sin = k_sin(qt,qf,r,delta_rad,gamma_rad,beta_rad,omega,d_omega,eff,d_eff,a_prod,a_det)
        # Effective CP coefficients
        A_cosh =  acc*dec*coeff_cosh*cosh_t
        B_sinh =  acc*dec*coeff_sinh*sinh_t
        C_cos =  acc*dec*res*coeff_cos*cos_t
        D_sin =  acc*dec*res*coeff_sin*sin_t
        plot_func(ax1,t,A_cosh,xmin=0,xmax=5,ymin=0,ymax=y_cosh,col=col,style=style,ytitle=r'$A\cosh(\Delta \Gamma_s t / 2)$',
                    title=title_cpv,leghead=r'$(q_t, q_f)$',label=r'({},{})'.format(qt,qf))
        plot_func(ax2,t,B_sinh,xmin=0,xmax=5,ymin=-y_sinh,ymax=y_sinh,col=col,style=style,ytitle=r'$B\sinh(\Delta \Gamma_s t / 2)$',
                    title=title_sinh,leghead=r'$(q_t, q_f)$',label=r'({},{})'.format(qt,qf))
        plot_func(ax3,t,C_cos,xmin=0,xmax=5,ymin=-y_cos,ymax=y_cos,col=col,style=style,ytitle=r'$C\cos(\Delta m_s t)$',
                    title=title_cos,leghead=r'$(q_t, q_f)$',label=r'({},{})'.format(qt,qf))
        plot_func(ax4,t,D_sin,xmin=0,xmax=5,ymin=-y_sin,ymax=y_sin,col=col,style=style,ytitle=r'$D\sin(\Delta m_s t)$',
                    title=title_sin,leghead=r'$(q_t, q_f)$',label=r'({},{})'.format(qt,qf))
    # Plot
    fig.tight_layout()
    if save: fig.savefig(name)
    plt.show()
