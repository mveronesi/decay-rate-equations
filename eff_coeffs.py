from decay_rate import C_f, C_fbar
from decay_rate import S_f, S_fbar
from decay_rate import A_f, A_fbar
import numpy as np

# P(t)~A*cosh(dGamma*t/2)+B*sinh(dGamma*t/2)+C*cos(dm*t)+D*sin(dm*t)
def k_cosh(qt,qf,omega=0,d_omega=0,eff=1,d_eff=1,a_prod=0,a_det=0):
    k_det = 1+qf*a_det
    if qt!=0:
        mistag = omega+qt*d_omega*0.5
        eff = eff+qt*d_eff*0.5
        k_tag = ((1-mistag)*(1+qt*a_prod)+mistag*(1-qt*a_prod))*eff
    else:
        k_tag = (1-(eff+qt*d_eff*0.5))*(1+qt*a_prod) + (1-(eff-qt*d_eff*0.5))*(1-qt*a_prod)
    return k_tag*k_det

def k_sinh(qt,qf,r,delta,gamma,beta,omega=0,d_omega=0,eff=1,d_eff=1,a_prod=0,a_det=0):
    k_det = 1+qf*a_det
    k_cpv = A_f(r,delta,gamma,beta) if qf>0 else A_fbar(r,delta,gamma,beta)
    if qt!=0:
        mistag = omega+qt*d_omega*0.5
        eff = eff+qt*d_eff*0.5
        k_tag = ((1-mistag)*(1+qt*a_prod)+mistag*(1-qt*a_prod))*eff
    else:
        k_tag = (1-(eff+qt*d_eff*0.5))*(1+qt*a_prod) + (1-(eff-qt*d_eff*0.5))*(1-qt*a_prod)
    return k_tag*k_det*k_cpv

def k_cos(qt,qf,r,omega=0,d_omega=0,eff=1,d_eff=1,a_prod=0,a_det=0):
    k_det = 1+qf*a_det
    k_cpv = C_f(r) if qf>0 else C_fbar(r)
    if qt!=0:
        mistag = omega+qt*d_omega*0.5
        eff = eff+qt*d_eff*0.5
        k_tag = qt*((1-mistag)*(1+qt*a_prod)-mistag*(1-qt*a_prod))*eff
    else:
        k_tag = (1-(eff+qt*d_eff*0.5))*(1+qt*a_prod) - (1-(eff-qt*d_eff*0.5))*(1-qt*a_prod)
    return k_tag*k_det*k_cpv

def k_sin(qt,qf,r,delta,gamma,beta,omega=0,d_omega=0,eff=1,d_eff=1,a_prod=0,a_det=0):
    k_det = 1+qf*a_det
    k_cpv = S_f(r,delta,gamma,beta) if qf>0 else S_fbar(r,delta,gamma,beta)
    if qt!=0:
        mistag = omega+qt*d_omega*0.5
        eff = eff+qt*d_eff*0.5
        k_tag = -qt*((1-mistag)*(1+qt*a_prod)-mistag*(1-qt*a_prod))*eff
    else:
        k_tag = -(1-(eff+qt*d_eff*0.5))*(1+qt*a_prod) + (1-(eff-qt*d_eff*0.5))*(1-qt*a_prod)
    return k_tag*k_det*k_cpv
