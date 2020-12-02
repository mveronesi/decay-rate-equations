import numpy as np

# CP Coefficients
# r - decay amplitudes ratio
# delta - strong phase
# gamma - weak phase
# beta - mixing phase

# Flavour specific
def k_cosh_fs(rez, imz, qt, qf):
    if qt == qf:
        return 1 + rez**2 + imz**2
    else:
        return 1

def k_cos_fs(rez, imz, qt, qf):
    if qt == qf:
        return 1 - rez**2 - imz**2
    else:
        return -1

def k_sinh_fs(rez, imz, qt, qf):
    if qt == qf:
        return qf*rez
    else:
        return 0

def k_sin_fs(rez, imz, qt, qf):
    if qt == qf:
        return qf*imz
    else:
        return 0

def P_t_fs(t,qt,qf,
           dm,dg,gs,
           rez,imz,
           afs):
    A = k_cosh_fs(rez, imz, qt, qf)
    B = k_sinh_fs(rez, imz, qt, qf)
    C = k_cos_fs(rez, imz, qt, qf)
    D = k_sin_fs(rez, imz, qt, qf)
    if qt==qf:
        pre_k = 1
    else:
        pre_k = (1-qt*afs)*np.sqrt((1-rez**2+imz**2)**2+4*(rez**2)*(imz**2))
    return pre_k*0.5*np.exp(-gs*t)*(A*np.cosh(0.5*dg*t)+B*np.sinh(0.5*dg*t)+C*np.cos(dm*t)+D*np.sin(dm*t))

def Acp(first,second):
    np.seterr(divide='ignore', invalid='ignore')
    num = first - second
    den = first + second
    return num/den
    # if (0 in den):
    #     return num/den
    # else:
    #     return den
