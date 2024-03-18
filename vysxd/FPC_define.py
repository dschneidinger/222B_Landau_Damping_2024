# Run Imports
import numpy as np
from vysxd_define import *
from vysxd_analysis import *
from scipy.ndimage import gaussian_filter

def FPC(directory):
    p1x1_files = np.sort(os.listdir(f'{directory}/MS/PHA/p1x1/electrons/'))
    p1x1 = vysxd_get_data('EPW-alves/MS/PHA/p1x1/electrons/p1x1-electrons-000000.h5')
    A_result = np.empty((len(p1x1.X),len(p1x1.Y),len(p1x1_files)-1))
    B_result = np.empty((len(p1x1.X),len(p1x1.Y),len(p1x1_files)-1))
    for i in range(len(p1x1_files)-1):
        p1x1 = vysxd_get_data('EPW-alves/MS/PHA/p1x1/electrons/' + p1x1_files[i])

        e1_D_xt = get_osiris_quantity_1d(f'{directory}/MS/FLD/e1/')[0]

        e1_D_xt = e1_D_xt[i]
        e1_D_xt = np.array([e1_D_xt]*len(p1x1.Y))

        dfdr = np.gradient(p1x1.DATA,axis=1)
        dfdv = np.gradient(p1x1.DATA,axis=0)
        e_times_dfdv = np.multiply(e1_D_xt,dfdv)

        vsquared_arr = np.transpose(np.array([p1x1.Y**2]*len(p1x1.X)))
        vcubed_arr = np.transpose(np.array([p1x1.Y**3]*len(p1x1.X)))

        A = 1/2*np.multiply(vcubed_arr,dfdr)
        B = -1/2*np.multiply(vsquared_arr, e_times_dfdv)

        A_result[:,:,i] = np.transpose(A)
        B_result[:,:,i] = np.transpose(B)
    return A_result, B_result