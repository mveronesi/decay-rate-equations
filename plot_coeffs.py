import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from itertools import product
from decay_rate import exp_dec
from detector_effects import eff_pow, dil_res
from decay_rate_compact import A_qf, C_qf, S_qf
from decay_rate_compact import k_cosh, k_sinh, k_cos, k_sin
from plot_utils import plot_func, pp
from plot_utils import cpv_leg, C_f_leg, S_f_leg, A_f_leg
from plot_utils import acc_leg, asymm_leg, tag_leg, res_leg

def plot_coeffs(dm=17.757,dg=0.085,gs=0.664,
                r=0.4,delta=10,gamma=60,beta=0,k=1,
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
                k_asymm=True,
                b_f=True,
                bbar_f=True,
                bbar_fbar=True,
                b_fbar=True,
                u_f=True,
                u_fbar=True):
    # fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(25,25))
    fig, (ax1,ax2,ax3,ax4) = plt.subplots(1,4, figsize=(40,10))
    # phases
    delta_rad = delta*np.pi/180
    gamma_rad = gamma*np.pi/180
    beta_rad = beta/1000
    # detector effects
    t = np.linspace(xmin,xmax,pp)
    dec = exp_dec(t,gs) if k_dec else 1
    acc = eff_pow(t,a=a_acc,n=n_acc,b=b_acc,beta=beta_acc,cutoff=cutoff_acc) if k_acc else 1
    sigma_t = sigma_t if k_res else 0
    (omega, d_omega, eff, d_eff) = (omega_tag, d_omega_tag, eff_tag, d_eff_tag) if k_tag else (0,0,1,0)
    a_prod = a_prod_asym if k_asymm else 0
    a_det = a_det_asym if k_asymm else 0
    # trigonometric functions
    cosh_t = np.cosh(0.5*dg*t)
    sinh_t = np.sinh(0.5*dg*t)
    cos_t = np.cos(0.5*dm*t)
    sin_t = np.sin(0.5*dm*t)
    # labels
    title_sinh = r'$(A^{\Delta\Gamma}_f,A^{\Delta\Gamma}_{\overline{f}})$ = '
    title_sinh +='({:.2f},{:.2f})'.format(A_qf(r,delta_rad,gamma_rad,beta_rad,k=k,qf=+1),A_qf(r,delta_rad,gamma_rad,beta_rad,k=k,qf=-1))
    title_cos = r'$(C_f,C_{\overline{f}})$ = '+'({:.2f},{:.2f})'.format(C_qf(r,+1),C_qf(r,-1))
    title_sin = r'$(S_f,S_{\overline{f}})$ = '+'({:.2f},{:.2f})'.format(S_qf(r,delta_rad,gamma_rad,beta_rad,k=k,qf=+1),S_qf(r,delta_rad,gamma_rad,beta_rad,k=k,qf=-1))
    title_cpv = r'$r$ = '+'{:.1f}'.format(r)
    title_cpv += r', $\delta$ = '+'{:.0f}'.format(delta) + r'$^{\circ}$'
    title_cpv += r', $\gamma$ = '+'{:.0f}'.format(gamma) + r'$^{\circ}$'
    title_cpv += r', $\beta_{s}$ = '+'{:.0f}'.format(beta) + r' mrad'
    q_list = [(1,1),(-1,1),(-1,-1),(1,-1),(0,1),(0,-1)]
    col_list = ['blue','red','blue','red','black','black']
    style_list = ['-','-','--','--','-','--']

    # qt_s = []
    # qf_s = []
    # if qt_plus: qt_s.append(+1)
    # if qt_nul: qt_s.append(0)
    # if qt_min: qt_s.append(-1)
    # if qf_plus: qf_s.append(+1)
    # if qf_min: qf_s.append(-1)

    q_prod = []
    if b_f: q_prod.append((1,1))
    if bbar_f: q_prod.append((-1,1))
    if bbar_fbar: q_prod.append((-1,-1))
    if b_fbar: q_prod.append((1,-1))
    if u_f: q_prod.append((0,1))
    if u_fbar: q_prod.append((0,-1))

    # leg heads
    leg_cosh = acc_leg(a_acc,n_acc,b_acc,beta_acc,cutoff_acc) if k_acc else ''
    leg_sinh = asymm_leg(a_prod,a_det) if k_asymm else ''
    leg_cos = tag_leg(omega,d_omega,eff,d_eff) if k_tag else ''
    leg_sin = res_leg(sigma_t) if k_res else ''
    for i in range(6):
        # if q_list[i] in product(qt_s,qf_s):
        if q_list[i] in q_prod:
            (qt,qf) = q_list[i]
            col=col_list[i]
            style=style_list[i]
            # Effective coefficients
            coeff_cosh = k_cosh(qt,qf,eff,d_eff,omega,d_omega,a_prod,a_det)
            coeff_sinh = k_sinh(qt,qf,r,delta_rad,gamma_rad,beta_rad,k,eff,d_eff,omega,d_omega,a_prod,a_det)
            coeff_cos = k_cos(qt,qf,r,dm,sigma_t,eff,d_eff,omega,d_omega,a_prod,a_det)
            coeff_sin = k_sin(qt,qf,r,delta_rad,gamma_rad,beta_rad,k,dm,sigma_t,eff,d_eff,omega,d_omega,a_prod,a_det)
            # Effective CP coefficients
            A_cosh =  acc*dec*coeff_cosh*cosh_t
            B_sinh =  acc*dec*coeff_sinh*sinh_t
            C_cos =  acc*dec*coeff_cos*cos_t
            D_sin =  acc*dec*coeff_sin*sin_t
            if qt>0: tag_lab = r'$B$'
            elif qt<0: tag_lab = r'$\overline{B}$'
            else: tag_lab = 'U'
            ch_lab = r'$f$' if qf>0 else r'$\overline{f}$'
            plot_func(ax1,t,A_cosh,xmin=xmin,xmax=xmax,ymin=0,ymax=y_cosh,col=col,style=style,ytitle=r'$A\cosh(\Delta \Gamma_s t / 2)$',
                        title=title_cpv,leghead='',label=r'{}$\to${}'.format(tag_lab,ch_lab),ypos=[-0.15, 0.8])
            plot_func(ax2,t,B_sinh,xmin=xmin,xmax=xmax,ymin=-y_sinh,ymax=y_sinh,col=col,style=style,ytitle=r'$B\sinh(\Delta \Gamma_s t / 2)$',
                        title=title_sinh,leghead='',label=r'{}$\to${}'.format(tag_lab,ch_lab),ypos=[-0.15, 0.8])
            plot_func(ax3,t,C_cos,xmin=xmin,xmax=xmax,ymin=-y_cos,ymax=y_cos,col=col,style=style,ytitle=r'$C\cos(\Delta m_s t)$',
                        title=title_cos,leghead='',label=r'{}$\to${}'.format(tag_lab,ch_lab),ypos=[-0.15, 0.85])
            plot_func(ax4,t,D_sin,xmin=xmin,xmax=xmax,ymin=-y_sin,ymax=y_sin,col=col,style=style,ytitle=r'$D\sin(\Delta m_s t)$',
                        title=title_sin,leghead='',label=r'{}$\to${}'.format(tag_lab,ch_lab),ypos=[-0.15, 0.85])
    # Plot
    fig.tight_layout()
    if save: fig.savefig(name)
    plt.show()
