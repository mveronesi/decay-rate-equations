import numpy as np
from decay_rate import B_f, Bbar_f, B_fbar, Bbar_fbar

# Acceptance
def eff_pow(t,a,n,b,beta,cutoff):
    return np.where((t<cutoff) | ((np.power((a*t),n)-b)<0),0,(1-1/(1+np.power((a*t),n)-b))*(1-beta*t))

# Time Resolution
def dil_res(sigma_t,dm):
    if sigma_t>0:
        return np.exp(-(sigma_t**2)*(dm**2)/2)
    else:
        return 1

# Observed Decay Rates
# dm,dg,gs - B mixing
# r,delta,gamma,beta - CPV
# t,a,b,beta,cutoff - acceptance
# sigma_t - decay time resolution
# omega_b(ar), eff_b(ar) - flavour tagging
# a_prod, a_det - production and detection asymmetries
def B_f_obs(t,dm,dg,gs,
            r,delta,gamma,beta,
            a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
            sigma_t=0,
            omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
            a_prod=0,a_det=0,
            k_acc=True,
            k_res=True,
            k_tag=True,
            k_asymm=True):
    acc = eff_pow(t,a=a_acc,n=n_acc,b=b_acc,beta=beta_acc,cutoff=cutoff_acc) if k_acc else 1
    sigma_t = sigma_t if k_res else 0
    (omega_b, omega_bar, eff_b, eff_bar) = (omega_b, omega_bar, eff_b, eff_bar) if k_tag else (0,0,1,1)
    (a_prod, a_det) = (a_prod,a_det) if k_asymm else (0,0)
    return acc*eff_b*((1-omega_b)*(1+a_prod)*(1+a_det)*B_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+omega_b*(1-a_prod)*(1+a_det)*Bbar_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))
    # return eff_pow(t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc)*eff_b*((1-omega_b)*(1+a_prod)*(1+a_det)*B_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+omega_b*(1-a_prod)*(1+a_det)*Bbar_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))

def Bbar_f_obs(t,dm,dg,gs,
               r,delta,gamma,beta,
               a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
               sigma_t=0,
               omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
               a_prod=0,a_det=0,
               k_acc=True,
               k_res=True,
               k_tag=True,
               k_asymm=True):
    acc = eff_pow(t,a=a_acc,n=n_acc,b=b_acc,beta=beta_acc,cutoff=cutoff_acc) if k_acc else 1
    sigma_t = sigma_t if k_res else 0
    (omega_b, omega_bar, eff_b, eff_bar) = (omega_b, omega_bar, eff_b, eff_bar) if k_tag else (0,0,1,1)
    (a_prod, a_det) = (a_prod,a_det) if k_asymm else (0,0)
    return acc*eff_bar*((1-omega_bar)*(1-a_prod)*(1+a_det)*Bbar_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+omega_bar*(1+a_prod)*(1+a_det)*B_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))
    # return eff_pow(t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc)*eff_bar*((1-omega_bar)*(1-a_prod)*(1+a_det)*Bbar_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+omega_bar*(1+a_prod)*(1+a_det)*B_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))

def B_fbar_obs(t,dm,dg,gs,
               r,delta,gamma,beta,
               a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
               sigma_t=0,
               omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
               a_prod=0,a_det=0,
               k_acc=True,
               k_res=True,
               k_tag=True,
               k_asymm=True):
    acc = eff_pow(t,a=a_acc,n=n_acc,b=b_acc,beta=beta_acc,cutoff=cutoff_acc) if k_acc else 1
    sigma_t = sigma_t if k_res else 0
    (omega_b, omega_bar, eff_b, eff_bar) = (omega_b, omega_bar, eff_b, eff_bar) if k_tag else (0,0,1,1)
    (a_prod, a_det) = (a_prod,a_det) if k_asymm else (0,0)
    return acc*eff_b*((1-omega_b)*(1+a_prod)*(1-a_det)*B_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+omega_b*(1-a_prod)*(1-a_det)*Bbar_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))
    # return eff_pow(t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc)*eff_b*((1-omega_b)*(1+a_prod)*(1-a_det)*B_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+omega_b*(1-a_prod)*(1-a_det)*Bbar_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))

def Bbar_fbar_obs(t,dm,dg,gs,
                  r,delta,gamma,beta,
                  a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                  sigma_t=0,
                  omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
                  a_prod=0,a_det=0,
                  k_acc=True,
                  k_res=True,
                  k_tag=True,
                  k_asymm=True):
    acc = eff_pow(t,a=a_acc,n=n_acc,b=b_acc,beta=beta_acc,cutoff=cutoff_acc) if k_acc else 1
    sigma_t = sigma_t if k_res else 0
    (omega_b, omega_bar, eff_b, eff_bar) = (omega_b, omega_bar, eff_b, eff_bar) if k_tag else (0,0,1,1)
    (a_prod, a_det) = (a_prod,a_det) if k_asymm else (0,0)
    return acc*eff_bar*((1-omega_bar)*(1-a_prod)*(1-a_det)*Bbar_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+omega_bar*(1+a_prod)*(1-a_det)*B_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))
    # return eff_pow(t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc)*eff_bar*((1-omega_bar)*(1-a_prod)*(1-a_det)*Bbar_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+omega_bar*(1+a_prod)*(1-a_det)*B_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))

def Untagged_f_obs(t,dm,dg,gs,
                   r,delta,gamma,beta,
                   a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                   sigma_t=0,
                   omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
                   a_prod=0,a_det=0,
                   k_acc=True,
                   k_res=True,
                   k_tag=True,
                   k_asymm=True):
    acc = eff_pow(t,a=a_acc,n=n_acc,b=b_acc,beta=beta_acc,cutoff=cutoff_acc) if k_acc else 1
    sigma_t = sigma_t if k_res else 0
    (omega_b, omega_bar, eff_b, eff_bar) = (omega_b, omega_bar, eff_b, eff_bar) if k_tag else (0,0,1,1)
    (a_prod, a_det) = (a_prod,a_det) if k_asymm else (0,0)
    return acc*((1-eff_b)*(1+a_prod)*(1+a_det)*B_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+(1-eff_bar)*(1-a_prod)*(1+a_det)*Bbar_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))
    # return eff_pow(t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc)*((1-eff_b)*(1+a_prod)*(1+a_det)*B_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+(1-eff_bar)*(1-a_prod)*(1+a_det)*Bbar_f(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))

def Untagged_fbar_obs(t,dm,dg,gs,
                      r,delta,gamma,beta,
                      a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                      sigma_t=0,
                      omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
                      a_prod=0,a_det=0,
                      k_acc=True,
                      k_res=True,
                      k_tag=True,
                      k_asymm=True):
    acc = eff_pow(t,a=a_acc,n=n_acc,b=b_acc,beta=beta_acc,cutoff=cutoff_acc) if k_acc else 1
    sigma_t = sigma_t if k_res else 0
    (omega_b, omega_bar, eff_b, eff_bar) = (omega_b, omega_bar, eff_b, eff_bar) if k_tag else (0,0,1,1)
    (a_prod, a_det) = (a_prod,a_det) if k_asymm else (0,0)
    return acc*((1-eff_b)*(1+a_prod)*(1-a_det)*B_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+(1-eff_bar)*(1-a_prod)*(1-a_det)*Bbar_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))
    # return eff_pow(t,a_acc,n_acc,b_acc,beta_acc,cutoff_acc)*((1-eff_b)*(1+a_prod)*(1-a_det)*B_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t)+(1-eff_bar)*(1-a_prod)*(1-a_det)*Bbar_fbar(t,dm,dg,gs,r,delta,gamma,beta,sigma_t))

# Observed Asymmetries
def Amix_f_obs(t,dm,dg,gs,
               r,delta,gamma,beta,
               a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
               sigma_t=0,
               omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
               a_prod=0,a_det=0,
               k_acc=True,
               k_res=True,
               k_tag=True,
               k_asymm=True):
    num = B_f_obs(t,dm,dg,gs,r,delta,gamma,beta,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t,omega_b,omega_bar,eff_b,eff_bar,a_prod,a_det,k_acc,k_res,k_tag,k_asymm)-Bbar_f_obs(t,dm,dg,gs,r,delta,gamma,beta,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t,omega_b,omega_bar,eff_b,eff_bar,a_prod,a_det,k_acc,k_res,k_tag,k_asymm)
    den = B_f_obs(t,dm,dg,gs,r,delta,gamma,beta,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t,omega_b,omega_bar,eff_b,eff_bar,a_prod,a_det,k_acc,k_res,k_tag,k_asymm)+Bbar_f_obs(t,dm,dg,gs,r,delta,gamma,beta,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t,omega_b,omega_bar,eff_b,eff_bar,a_prod,a_det,k_acc,k_res,k_tag,k_asymm)
    return num/den

def Amix_fbar_obs(t,dm,dg,gs,
                  r,delta,gamma,beta,
                  a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                  sigma_t=0,
                  omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
                  a_prod=0,a_det=0,
                  k_acc=True,
                  k_res=True,
                  k_tag=True,
                  k_asymm=True):
    num = Bbar_fbar_obs(t,dm,dg,gs,r,delta,gamma,beta,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t,omega_b,omega_bar,eff_b,eff_bar,a_prod,a_det,k_acc,k_res,k_tag,k_asymm)-B_fbar_obs(t,dm,dg,gs,r,delta,gamma,beta,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t,omega_b,omega_bar,eff_b,eff_bar,a_prod,a_det,k_acc,k_res,k_tag,k_asymm)
    den = Bbar_fbar_obs(t,dm,dg,gs,r,delta,gamma,beta,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t,omega_b,omega_bar,eff_b,eff_bar,a_prod,a_det,k_acc,k_res,k_tag,k_asymm)+B_fbar_obs(t,dm,dg,gs,r,delta,gamma,beta,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t,omega_b,omega_bar,eff_b,eff_bar,a_prod,a_det,k_acc,k_res,k_tag,k_asymm)
    return num/den

# Folded Asymmetries
def Amix_f_obs_fold(t,dm,dg,gs,
                    r,delta,gamma,beta,
                    a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                    sigma_t=0,
                    omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
                    a_prod=0,a_det=0,
                    k_acc=True,
                    k_res=True,
                    k_tag=True,
                    k_asymm=True):
    nosc = int((t[-1] - t[0])/(2*np.pi/dm))
    return np.sum(np.split(Amix_f_obs(t,dm,dg,gs,r,delta,gamma,beta,
                                      a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                                      sigma_t,omega_b,omega_bar,eff_b,eff_bar,
                                      a_prod,a_det,k_acc,k_res,k_tag,k_asymm),nosc),axis=0)/nosc

def Amix_fbar_obs_fold(t,dm,dg,gs,
                       r,delta,gamma,beta,
                       a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                       sigma_t=0,
                       omega_b=0,omega_bar=0,eff_b=1,eff_bar=1,
                       a_prod=0,a_det=0,
                       k_acc=True,
                       k_res=True,
                       k_tag=True,
                       k_asymm=True):
    nosc = int((t[-1] - t[0])/(2*np.pi/dm))
    return np.sum(np.split(Amix_fbar_obs(t,dm,dg,gs,r,delta,gamma,beta,
                                         a_acc,n_acc,b_acc,beta_acc,cutoff_acc,
                                         sigma_t,omega_b,omega_bar,eff_b,eff_bar,
                                         a_prod,a_det,k_acc,k_res,k_tag,k_asymm),nosc),axis=0)/nosc
