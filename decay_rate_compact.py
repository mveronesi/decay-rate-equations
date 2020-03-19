import numpy as np

# Decay Rate Coefficients
# r - decay amplitudes ratio
# delta - strong phase
# gamma - weak phase
# beta - mixing phase
def C_qf(r,qf):
    return qf*(1-r**2)/(1+r**2)

def S_qf(r,delta,gamma,beta,qf):
    return qf*2*r*np.sin(delta-qf*(gamma-2*beta))/(1+r**2)

def A_qf(r,delta,gamma,beta,qf):
    return -2*r*np.cos(delta-qf*(gamma-2*beta))/(1+r**2)

# resolution
def dil(sigma_t,dm):
    return np.exp(-(sigma_t**2)*(dm**2)/2)

# tagging
def mistag(qt,omega,d_omega):
    return omega+0.5*qt*d_omega

def tageff(qt,eff,d_eff):
    return eff+0.5*qt*d_eff

# effective CP coeffs
def k_cosh(qt,qf,
           eff,d_eff,omega,d_omega,a_prod,a_det):
    if qt==0:
        k_tag = ( ( 1 - tageff(+1,eff,d_eff) )*(1-a_prod)
                + ( 1 - tageff(-1,eff,d_eff) )*(1+a_prod) )*(1-qf*a_det)
    else:
        k_tag = tageff(qt,eff,d_eff)*( ( 1 - mistag(qt,omega,d_omega) )*(1-qt*a_prod)
                                           + mistag(qt,omega,d_omega)*(1+qt*a_prod) )*(1-qf*a_det)
    return k_tag

def k_sinh(qt,qf,r,delta,gamma,beta,
           eff,d_eff,omega,d_omega,a_prod,a_det):
    A_f = A_qf(r,delta,gamma,beta,qf)
    if qt==0:
        k_tag = ( ( 1 - tageff(+1,eff,d_eff) )*(1-a_prod)
                + ( 1 - tageff(-1,eff,d_eff) )*(1+a_prod) )*(1-qf*a_det)
    else:
        k_tag = tageff(qt,eff,d_eff)*( ( 1 - mistag(qt,omega,d_omega) )*(1-qt*a_prod)
                                           + mistag(qt,omega,d_omega)*(1+qt*a_prod) )*(1-qf*a_det)
    return k_tag*A_f

def k_cos(qt,qf,r,dm,sigma_t,
          eff,d_eff,omega,d_omega,a_prod,a_det):
    C_f = C_qf(r,qf)
    k_dil = dil(sigma_t,dm)
    if qt==0:
        k_tag = ( ( 1 - tageff(+1,eff,d_eff) )*(1-a_prod)
                - ( 1 - tageff(-1,eff,d_eff) )*(1+a_prod) )*(1-qf*a_det)
    else:
        k_tag = qt*tageff(qt,eff,d_eff)*( ( 1 - mistag(qt,omega,d_omega) )*(1-qt*a_prod)
                                              - mistag(qt,omega,d_omega)*(1+qt*a_prod) )*(1-qf*a_det)
    return k_tag*k_dil*C_f

def k_sin(qt,qf,r,delta,gamma,beta,dm,sigma_t,
          eff,d_eff,omega,d_omega,a_prod,a_det):
    S_f = S_qf(r,delta,gamma,beta,qf)
    k_dil = dil(sigma_t,dm)
    if qt==0:
        k_tag = ( - ( 1 - tageff(+1,eff,d_eff) )*(1-a_prod)
                  + ( 1 - tageff(-1,eff,d_eff) )*(1+a_prod) )*(1-qf*a_det)
    else:
        k_tag = -qt*tageff(qt,eff,d_eff)*( ( 1 - mistag(qt,omega,d_omega) )*(1-qt*a_prod)
                                               - mistag(qt,omega,d_omega)*(1+qt*a_prod) )*(1-qf*a_det)
    return k_tag*k_dil*S_f


def P_t(t,qt,qf,
        dm,dg,gs,
        r,delta,gamma,beta,
        sigma_t=0,
        eff=1,d_eff=0,omega=0,d_omega=0,
        a_prod=0,a_det=0):
    A = k_cosh(qt,qf,eff,d_eff,omega,d_omega,a_prod,a_det)
    B = k_sinh(qt,qf,r,delta,gamma,beta,eff,d_eff,omega,d_omega,a_prod,a_det)
    C = k_cos(qt,qf,r,dm,sigma_t,eff,d_eff,omega,d_omega,a_prod,a_det)
    D = k_sin(qt,qf,r,delta,gamma,beta,dm,sigma_t,eff,d_eff,omega,d_omega,a_prod,a_det)
    return np.exp(-gs*t)*(A*np.cosh(0.5*dg*t)+B*np.sinh(0.5*dg*t)+C*np.cos(dm*t)+D*np.sin(dm*t))

def Amix_qf(t,qf,
            dm,dg,gs,
            r,delta,gamma,beta,
            sigma_t=0,
            eff=1,d_eff=0,omega=0,d_omega=0,
            a_prod=0,a_det=0):
    first = P_t(t,qf,qf,dm,dg,gs,r,delta,gamma,beta,sigma_t,eff,d_eff,omega,d_omega,a_prod,a_det)
    second = P_t(t,-qf,qf,dm,dg,gs,r,delta,gamma,beta,sigma_t,eff,d_eff,omega,d_omega,a_prod,a_det)
    num = first - second
    den = first + second
    return num/den

def Afold_qf(t,qf,
             dm,dg,gs,
             r,delta,gamma,beta,
             sigma_t=0,
             eff=1,d_eff=0,omega=0,d_omega=0,
             a_prod=0,a_det=0):
    nosc = int((t[-1] - t[0])/(2*np.pi/dm))
    return np.sum(np.split(Amix_qf(t,qf,dm,dg,gs,r,delta,gamma,beta,sigma_t,eff,d_eff,omega,d_omega,a_prod,a_det),nosc),axis=0)/nosc
