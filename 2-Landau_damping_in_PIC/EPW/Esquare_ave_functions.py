import sys
import numpy as np
import matplotlib.pyplot as plt

def total_field_energy(T,e1_D_xt):
    esquare = (e1_D_xt)**2
    n_x = np.linspace(0, len(x)-1,len(x))
    n_x = n_x.astype(int)
    n_x = n_x.tolist()
    n_t = np.linspace(0, len(T)-1,len(T))
    n_t = n_t.astype(int)
    n_t = n_t.tolist()
    TFE = []
    for i in n_t:
        integral = np.trapz(esquare[i,n_x], x)
        TFE.append(integral)
    return TFE

def total_field_energy_ave(T,e1_D_xt):
    T_input = T.tolist()
    T_ave = T_input[3:-3]
    TFEA = []
    TFE_ave = total_field_energy(T_input,e1_D_xt)
    for i in range(0,len(T_ave)):
        coarse_grain = (1/7)*(0.3*TFE_ave[i]+TFE_ave[i+1]+1.4*TFE_ave[i+2]+1.6*TFE_ave[i+3]+1.4*TFE_ave[i+4]+TFE_ave[i+5]+0.3*TFE_ave[i+6])
        TFEA.append(coarse_grain)
    return TFEA

def T_ave(t):
    T_input = t.tolist()
    return T_input[3:-3]

print(total_field_energy_ave(t,e1_D_xt))

def logderivativeEsquare(T,e1_D_xt):
    return np.gradient(np.log(total_field_energy_ave(T,e1_D_xt)))