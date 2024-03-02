import numpy as np
from plasmapy.dispersion import plasma_dispersion_func_deriv as zprime

def ion_dispersion(w,k,ti_over_te,mi_over_me = 1836.):
    return 1 - 1/(2*k**2)*zprime(1/np.sqrt(2)*(w/k))-1/(2*k**2)*1/ti_over_te*zprime(1/np.sqrt(2)*w/k*np.sqrt(mi_over_me/ti_over_te))

def dfMdv(v):
    return 1./np.sqrt(2*np.pi) * (-v) * np.exp(-0.5*v**2)

def dfBOTdv(v, v0, nb):
    return (1-nb)*dfMdv(v) + 0.5*nb*(dfMdv(v-v0)) + 0.5*nb*(dfMdv(v+v0))

def D_analytic(wr, wi, k0, v0, nb):
    vphi = wr/k0
    eps_v = 1e-3
    v = np.linspace(-12, vphi-eps_v,2048)
    integral = np.trapz(dfBOTdv(v, v0, nb) / (v-(wr+1j*wi)/k0), v)
    v = np.linspace(vphi+eps_v,12, 2048)
    integral += np.trapz(dfBOTdv(v, v0, nb) / (v-(wr+1j*wi)/k0), v)
    d = 1.-1./(k0**2) *integral
    # deforming contour under pole
    if wi == 0:
        d -= 1j*np.pi/k0**2 * dfBOTdv((wr+1j*wi)/k0, v0, nb)
    elif wi<0:
        d -= 2j*np.pi/k0**2 * dfBOTdv((wr+1j*wi)/k0, v0, nb)
    return d

def find_w_roots(k0, v0, nb, method = D_analytic):
    Wr = np.linspace(0,1.7,128)
    Wi = np.linspace(-0.3,0.3,101)

    D_arr = np.zeros((len(Wi), len(Wr)), dtype = complex)
    if (method == D_analytic):
        for i in range(len(Wi)):
            for j in range(len(Wr)):
                D_arr[i,j] = method(Wr[j], Wi[i], k0=k0, v0 = v0, nb = nb)
    elif (method == z_func):
        for i in range(len(Wi)):
            for j in range(len(Wr)):
                D_arr[i,j] = method(Wr[j], Wi[i], k0=k0)

    
    D_arr = np.abs(D_arr)
    zeros = np.where(D_arr<2*D_arr.min())
    return Wr[zeros[1]][:], Wi[zeros[0]][:]
