import numpy as np

# Decay Rate Coefficients
# r - decay amplitudes ratio
# delta - strong phase
# gamma - weak phase
# beta - mixing phase
def C_f(r):
    return (1-r**2)/(1+r**2)

def C_fbar(r):
    return -(1-r**2)/(1+r**2)

def S_f(r,delta,gamma,beta):
    return 2*r*np.sin(delta-(gamma-2*beta))/(1+r**2)

def S_fbar(r,delta,gamma,beta):
    return -2*r*np.sin(delta+(gamma-2*beta))/(1+r**2)

def D_f(r,delta,gamma,beta):
    return -2*r*np.cos(delta-(gamma-2*beta))/(1+r**2)

def D_fbar(r,delta,gamma,beta):
    return -2*r*np.cos(delta+(gamma-2*beta))/(1+r**2)

# Decay Rate Equations
# dm - mass difference
# dg - decay width difference
# gs - decay width average
# sigma_t - decay time resolution
def exp_dec(t,gs):
    return np.exp(-gs*t)

def dec_rate(t,dm,dg,gs,D,C,S,sigma_t=0):
    if sigma_t>0:
        dil = np.exp(-(sigma_t**2)*(dm**2)/2)
    else:
        dil = 1
    return np.exp(-gs*t)*(np.cosh(dg*t/2)+D*np.sinh(dg*t/2)+dil*C*np.cos(dm*t)+dil*S*np.sin(dm*t))

def B_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t=0):
    C_f_val = C_f(r)
    S_f_val = S_f(r,delta,gamma,beta)
    D_f_val = D_f(r,delta,gamma,beta)
    return dec_rate(t,dm,dg,gs,D_f_val,C_f_val,-S_f_val,sigma_t=sigma_t)

def Bbar_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t=0):
    C_f_val = C_f(r)
    S_f_val = S_f(r,delta,gamma,beta)
    D_f_val = D_f(r,delta,gamma,beta)
    return dec_rate(t,dm,dg,gs,D_f_val,-C_f_val,S_f_val,sigma_t=sigma_t)

def B_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t=0):
    C_fbar_val = C_fbar(r)
    S_fbar_val = S_fbar(r,delta,gamma,beta)
    D_fbar_val = D_fbar(r,delta,gamma,beta)
    return dec_rate(t,dm,dg,gs,D_fbar_val,C_fbar_val,-S_fbar_val,sigma_t=sigma_t)

def Bbar_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t=0):
    C_fbar_val = C_fbar(r)
    S_fbar_val = S_fbar(r,delta,gamma,beta)
    D_fbar_val = D_fbar(r,delta,gamma,beta)
    return dec_rate(t,dm,dg,gs,D_fbar_val,-C_fbar_val,S_fbar_val,sigma_t=sigma_t)

# Mixing Asymmetries
def Amix_f(t,dm,dg,gs,r,delta,gamma,beta):
    num = B_f(t,dm,dg,gs,r,delta,gamma,beta)-Bbar_f(t,dm,dg,gs,r,delta,gamma,beta)
    den = B_f(t,dm,dg,gs,r,delta,gamma,beta)+Bbar_f(t,dm,dg,gs,r,delta,gamma,beta)
    return num/den
def Amix_fbar(t,dm,dg,gs,r,delta,gamma,beta):
    num = Bbar_fbar(t,dm,dg,gs,r,delta,gamma,beta)-B_fbar(t,dm,dg,gs,r,delta,gamma,beta)
    den = Bbar_fbar(t,dm,dg,gs,r,delta,gamma,beta)+B_fbar(t,dm,dg,gs,r,delta,gamma,beta)
    return num/den

# Folded Asymmetries
def Amix_f_fold(t,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad):
    nosc = int((t[-1] - t[0])/(2*np.pi/dm))
    return np.sum(np.split(Amix_f(t,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad),nosc),axis=0)/nosc

def Amix_fbar_fold(t,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad):
    nosc = int((t[-1] - t[0])/(2*np.pi/dm))
    return np.sum(np.split(Amix_fbar(t,dm,dg,gs,r,delta_rad,gamma_rad,beta_rad),nosc),axis=0)/nosc
