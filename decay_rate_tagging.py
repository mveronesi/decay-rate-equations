import numpy as np

# cp coefficients
def C_qf(r,qf):
    return qf*(1-r**2)/(1+r**2)

def S_qf(r,delta,gamma,beta,qf=1):
    return qf*2*r*np.sin(delta-qf*(gamma-2*beta))/(1+r**2)

def A_qf(r,delta,gamma,beta,qf=1):
    return -2*r*np.cos(delta-qf*(gamma-2*beta))/(1+r**2)

# tagging
def mistag(omega,d_omega,qt):
    return omega+0.5*qt*d_omega

def tageff(eff,d_eff,qt):
    return eff+0.5*qt*d_eff

# effective decay rate coefficients
def k_cosh(qt,qf,r,delta,gamma,beta,rez,imz,afs,
           eff,d_eff,omega,d_omega,a_prod,a_det,
           qp=+1 # property of the tagger -1 vs. meson +1
           ):
    A_f = A_qf(r=r,delta=delta,gamma=gamma,beta=beta,qf=qf)
    S_f = S_qf(r=r,delta=delta,gamma=gamma,beta=beta,qf=qf)
    if qt!=0:
        k_tag = 0.25*( tageff(eff=eff,d_eff=d_eff,qt=qt)*( 1 - mistag(omega=omega,d_omega=d_omega,qt=qt) )*(1-qt*a_prod)*(1+qt*rez*A_f-imz*S_f)*(1-int(float(qt-qf)/2)*afs)
                      +tageff(eff=eff,d_eff=d_eff,qt=-qt*qp)*mistag(omega=omega,d_omega=d_omega,qt=-qt*qp)*(1+qt*a_prod)*(1-qt*rez*A_f-imz*S_f)*(1+int(float(qt+qf)/2)*afs)
                     )*(1-qf*a_det)
    else:
        k_tag = 0.25*( (1-tageff(eff=eff,d_eff=d_eff,qt=-1))*(1-a_prod)*(1+rez*A_f-imz*S_f)*(1-int(float(1-qf)/2)*afs)
                      +(1-tageff(eff=eff,d_eff=d_eff,qt=+1))*(1+a_prod)*(1-rez*A_f-imz*S_f)*(1+int(float(1+qf)/2)*afs)
                     )*(1-qf*a_det)
    return k_tag

def k_sinh(qt,qf,r,delta,gamma,beta,rez,imz,afs,
           eff,d_eff,omega,d_omega,a_prod,a_det,
           qp=+1 # property of the tagger -1 vs. meson +1
           ):
    A_f = A_qf(r=r,delta=delta,gamma=gamma,beta=beta,qf=qf)
    C_f = C_qf(r=r,qf=qf)
    if qt!=0:
        k_tag = 0.25*( tageff(eff=eff,d_eff=d_eff,qt=qt)*( 1 - mistag(omega=omega,d_omega=d_omega,qt=qt) )*(1-qt*a_prod)*(A_f+rez*(C_f+qt))*(1-int(float(qt-qf)/2)*afs)
                      +tageff(eff=eff,d_eff=d_eff,qt=-qt*qp)*mistag(omega=omega,d_omega=d_omega,qt=-qt*qp)*(1+qt*a_prod)*(A_f+rez*(C_f-qt))*(1+int(float(qt+qf)/2)*afs)
                     )*(1-qf*a_det)
    else:
        k_tag = 0.25*( (1-tageff(eff=eff,d_eff=d_eff,qt=-1))*(1-a_prod)*(A_f+rez*(C_f+1))*(1-int(float(1-qf)/2)*afs)
                      +(1-tageff(eff=eff,d_eff=d_eff,qt=+1))*(1+a_prod)*(A_f+rez*(C_f-1))*(1+int(float(1+qf)/2)*afs)
                     )*(1-qf*a_det)
    return k_tag

def k_cos(qt,qf,r,delta,gamma,beta,rez,imz,afs,
          eff,d_eff,omega,d_omega,a_prod,a_det,
          qp=+1 # property of the tagger -1 vs. meson +1
          ):
    A_f = A_qf(r=r,delta=delta,gamma=gamma,beta=beta,qf=qf)
    S_f = S_qf(r=r,delta=delta,gamma=gamma,beta=beta,qf=qf)
    C_f = C_qf(r=r,qf=qf)
    if qt!=0:
        k_tag = 0.25*( tageff(eff=eff,d_eff=d_eff,qt=qt)*( 1 - mistag(omega=omega,d_omega=d_omega,qt=qt) )*(1-qt*a_prod)*(qt*C_f-qt*rez*A_f+imz*S_f)*(1-int(float(qt-qf)/2)*afs)
                      +tageff(eff=eff,d_eff=d_eff,qt=-qt*qp)*mistag(omega=omega,d_omega=d_omega,qt=-qt*qp)*(1+qt*a_prod)*(-qt*C_f+qt*rez*A_f+imz*S_f)*(1+int(float(qt+qf)/2)*afs)
                     )*(1-qf*a_det)
    else:
        k_tag = 0.25*( (1-tageff(eff=eff,d_eff=d_eff,qt=-1))*(1-a_prod)*(C_f-rez*A_f+imz*S_f)*(1-int(float(1-qf)/2)*afs)
                      +(1-tageff(eff=eff,d_eff=d_eff,qt=+1))*(1+a_prod)*(-C_f+rez*A_f+imz*S_f)*(1+int(float(1+qf)/2)*afs)
                     )*(1-qf*a_det)
    return k_tag

def k_sin(qt,qf,r,delta,gamma,beta,rez,imz,afs,
          eff,d_eff,omega,d_omega,a_prod,a_det,
          qp=+1 # property of the tagger -1 vs. meson +1
          ):
    S_f = S_qf(r=r,delta=delta,gamma=gamma,beta=beta,qf=qf)
    C_f = C_qf(r=r,qf=qf)
    if qt!=0:
        k_tag = 0.25*( tageff(eff=eff,d_eff=d_eff,qt=qt)*( 1 - mistag(omega=omega,d_omega=d_omega,qt=qt) )*(1-qt*a_prod)*(-qt*S_f+imz*(C_f+qt))*(1-int(float(qt-qf)/2)*afs)
                      +tageff(eff=eff,d_eff=d_eff,qt=-qt*qp)*mistag(omega=omega,d_omega=d_omega,qt=-qt*qp)*(1+qt*a_prod)*(qt*S_f+imz*(C_f-qt))*(1+int(float(qt+qf)/2)*afs)
                     )*(1-qf*a_det)
    else:
        k_tag = 0.25*( (1-tageff(eff=eff,d_eff=d_eff,qt=-1))*(1-a_prod)*(-S_f+imz*(C_f+1))*(1-int(float(1-qf)/2)*afs)
                      +(1-tageff(eff=eff,d_eff=d_eff,qt=+1))*(1+a_prod)*(S_f+imz*(C_f-1))*(1+int(float(1+qf)/2)*afs)
                     )*(1-qf*a_det)
    return k_tag

# Compute asymmetry
def Acp(first,second):
    np.seterr(divide='ignore', invalid='ignore')
    num = first - second
    den = first + second
    return num/den
