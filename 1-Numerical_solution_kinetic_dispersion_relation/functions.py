import numpy as np
from plasmapy.dispersion import plasma_dispersion_func_deriv as zprime

def ion_dispersion(w,k,ti_over_te,mi_over_me = 1836.):
    return 1 - 1/(2*k**2)*zprime(1/np.sqrt(2)*(w/k))-1/(2*k**2)*1/ti_over_te*zprime(1/np.sqrt(2)*w/k*np.sqrt(mi_over_me/ti_over_te))