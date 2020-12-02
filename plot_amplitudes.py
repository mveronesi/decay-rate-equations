import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from itertools import product
from decay_rate import exp_dec
from decay_rate_tagging import A_qf, C_qf, S_qf
from decay_rate_tagging import k_cosh, k_sinh, k_cos, k_sin
from detector_effects import eff_pow, dil_res
from plot_utils import plot_func, pp
from plot_utils import cpv_leg, C_f_leg, S_f_leg, A_f_leg
from plot_utils import acc_leg, asymm_leg, tag_leg, res_leg

def plot_amplitudes(dm=17.757,dg=0.085,gs=0.664,
                    r=0,delta=0,gamma=0,beta=0,
                    rez=0,imz=0,afs=0,
                    omega=0,d_omega=0,eff=1,d_eff=0,
                    a_prod=0,a_det=0,sigma_t=0,
                    a_acc=1.5,n_acc=1.5,b_acc=0.05,beta_acc=0.03,cutoff_acc=0.2,
                    xmin=0,xmax=5,
                    y_cosh=0.5,y_sinh=0.01,y_cos=0.05,y_sin=0.05,y_dg=1,
                    name='plot.eps',save=False,
                    b_f=True,
                    bbar_f=True,
                    bbar_fbar=True,
                    b_fbar=True,
                    u_f=True,
                    u_fbar=True,
                    tagging='MESON', # MESON/TAGGER
                    k_acc = False,
                    k_tag = False,
                    k_sum_f =False,
                    k_sum_fb =False
                    # asymm='DECRATE' # DECRATE/MIX/CP/CPT
                    ):

    if tagging == 'MESON':
        qp=+1
    elif tagging == 'TAGGER':
        qp=-1
    # turn off flavour tag
    if not k_tag:
        omega=0
        d_omega=0
        eff=1
        d_eff=0
        a_prod=0
        a_det=0

    # fig, (ax1,ax2,ax3,ax4) = plt.subplots(1,4, figsize=(40,10))
    fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1,5, figsize=(50,10))
    # phases
    delta_rad = delta*np.pi/180
    gamma_rad = gamma*np.pi/180
    beta_rad = beta/1000
    # detector effects
    t = np.linspace(xmin,xmax,pp)
    dec = exp_dec(t,gs)
    dil = dil_res(sigma_t=sigma_t,dm=dm)
    acc = eff_pow(t=t,a=a_acc,n=n_acc,b=b_acc,beta=beta_acc,cutoff=cutoff_acc) if k_acc else 1
    # trigonometric functions
    cosh_t = np.cosh(0.5*dg*t)
    sinh_t = np.sinh(0.5*dg*t)
    cos_t = np.cos(0.5*dm*t)
    sin_t = np.sin(0.5*dm*t)
    # labels
    # title_sinh = r'$(A^{\Delta\Gamma}_f,A^{\Delta\Gamma}_{\overline{f}})$ = '
    # title_sinh +='({:.2f},{:.2f})'.format(A_qf(r,delta_rad,gamma_rad,beta_rad,qf=+1),A_qf(r,delta_rad,gamma_rad,beta_rad,qf=-1))
    # title_cos = r'$(C_f,C_{\overline{f}})$ = '+'({:.2f},{:.2f})'.format(C_qf(r,+1),C_qf(r,-1))
    # title_sin = r'$(S_f,S_{\overline{f}})$ = '+'({:.2f},{:.2f})'.format(S_qf(r,delta_rad,gamma_rad,beta_rad,qf=+1),S_qf(r,delta_rad,gamma_rad,beta_rad,qf=-1))
    # title_cpv = r'$r$ = '+'{:.1f}'.format(r)
    # title_cpv += r', $\delta$ = '+'{:.0f}'.format(delta) + r'$^{\circ}$'
    # title_cpv += r', $\gamma$ = '+'{:.0f}'.format(gamma) + r'$^{\circ}$'
    # title_cpv += r', $\beta_{s}$ = '+'{:.0f}'.format(beta) + r' mrad'
    q_list = [(1,1),(-1,1),(-1,-1),(1,-1),(0,1),(0,-1)]
    col_list = ['blue','red','blue','red','black','black']
    style_list = ['-','-','--','--','-','--']

    q_prod = []
    if b_f: q_prod.append((1,1))
    if bbar_f: q_prod.append((-1,1))
    if bbar_fbar: q_prod.append((-1,-1))
    if b_fbar: q_prod.append((1,-1))
    if u_f: q_prod.append((0,1))
    if u_fbar: q_prod.append((0,-1))

    # leg heads
    # leg_cosh = acc_leg(a_acc,n_acc,b_acc,beta_acc,cutoff_acc) if k_acc else ''
    # leg_sinh = asymm_leg(a_prod,a_det) if k_asymm else ''
    # leg_cos = tag_leg(omega,d_omega,eff,d_eff) if k_tag else ''
    # leg_sin = res_leg(sigma_t) if k_res else ''
    sum_f = 0
    sum_fb = 0
    for i in range(6):
        if q_list[i] in q_prod:
            (qt,qf) = q_list[i]
            col=col_list[i]
            style=style_list[i]
            # Effective coefficients
            coeff_cosh = acc*k_cosh(qt=qt,qf=qf,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,rez=rez,imz=imz,afs=afs,
                                eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,a_prod=a_prod,a_det=a_det,qp=qp)
            coeff_sinh = acc*k_sinh(qt=qt,qf=qf,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,rez=rez,imz=imz,afs=afs,
                                eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,a_prod=a_prod,a_det=a_det,qp=qp)
            coeff_cos = acc*dil*k_cos(qt=qt,qf=qf,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,rez=rez,imz=imz,afs=afs,
                              eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,a_prod=a_prod,a_det=a_det,qp=qp)
            coeff_sin = acc*dil*k_sin(qt=qt,qf=qf,r=r,delta=delta_rad,gamma=gamma_rad,beta=beta_rad,rez=rez,imz=imz,afs=afs,
                              eff=eff,d_eff=d_eff,omega=omega,d_omega=d_omega,a_prod=a_prod,a_det=a_det,qp=qp)
            # Effective CP coefficients
            A_cosh = dec*coeff_cosh*cosh_t
            B_sinh = dec*coeff_sinh*sinh_t
            C_cos =  dec*coeff_cos*cos_t
            D_sin =  dec*coeff_sin*sin_t
            dG_dt = A_cosh + B_sinh + C_cos + D_sin

            if q_list[i][1]>0: sum_f += dG_dt
            else: sum_fb += dG_dt

            if qt>0: tag_lab = r'$B$'
            elif qt<0: tag_lab = r'$\overline{B}$'
            else: tag_lab = 'U'
            ch_lab = r'$f$' if qf>0 else r'$\overline{f}$'
            plot_func(ax1,t,A_cosh,xmin=xmin,xmax=xmax,ymin=0,ymax=y_cosh,col=col,style=style,ytitle=r'$k_{\cosh}^{q_t,q_f}\cosh(\Delta \Gamma_s t / 2)$ [a.u.]',
                        title='',leghead='',label=r'{}$\to${}'.format(tag_lab,ch_lab),ypos=[-0.15, 0.675])
            plot_func(ax2,t,B_sinh,xmin=xmin,xmax=xmax,ymin=-y_sinh,ymax=y_sinh,col=col,style=style,ytitle=r'$k_{\sinh}^{q_t,q_f}\sinh(\Delta \Gamma_s t / 2)$  [a.u.]',
                        title='',leghead='',label=r'{}$\to${}'.format(tag_lab,ch_lab),ypos=[-0.15, 0.675])
            plot_func(ax3,t,C_cos,xmin=xmin,xmax=xmax,ymin=-y_cos,ymax=y_cos,col=col,style=style,ytitle=r'$k_{\cos}^{q_t,q_f}\cos(\Delta m_s t)$ [a.u.]',
                        title='',leghead='',label=r'{}$\to${}'.format(tag_lab,ch_lab),ypos=[-0.15, 0.72])
            plot_func(ax4,t,D_sin,xmin=xmin,xmax=xmax,ymin=-y_sin,ymax=y_sin,col=col,style=style,ytitle=r'$k_{\sin}^{q_t,q_f}\sin(\Delta m_s t)$ [a.u.]',
                        title='',leghead='',label=r'{}$\to${}'.format(tag_lab,ch_lab),ypos=[-0.15, 0.72])
            plot_func(ax5,t,dG_dt,xmin=xmin,xmax=xmax,ymin=0,ymax=y_dg,col=col,style=style,ytitle=r'$d\Gamma^{q_t,q_f}/dt$ [a.u.]',
                        title='',leghead='',label=r'{}$\to${}'.format(tag_lab,ch_lab),ypos=[-0.15, 0.79])
    if k_sum_f:
        plot_func(ax5,t,sum_f,xmin=xmin,xmax=xmax,ymin=0,ymax=y_dg,col='gray',style='-',ytitle=r'$d\Gamma^{q_t,q_f}/dt$ [a.u.]',
                    title='',leghead='',label=r'Sum $f$',ypos=[-0.15, 0.79])
    if k_sum_fb:
        plot_func(ax5,t,sum_fb,xmin=xmin,xmax=xmax,ymin=0,ymax=y_dg,col='gray',style='--',ytitle=r'$d\Gamma^{q_t,q_f}/dt$ [a.u.]',
                    title='',leghead='',label=r'Sum $\bar{f}$',ypos=[-0.15, 0.79])
    # Plot
    fig.tight_layout()
    if save: fig.savefig(name)
    plt.show()
