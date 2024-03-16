import sys
sys.path.append('../../../vysxd')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

from vysxd_define import *
from vysxd_analysis import *

#if you are pointing to a particular folder add the / 
def Movie(directory, frame_rate = 5): 
    n = frame_rate
    # Create a sorted list of filenames you will be analyzing
    p1x1_files = np.sort(os.listdir(f'{directory}MS/PHA/p1x1/electrons/')) 
# Create a directory for the animation, if there isn't one already
    if not (os.path.isdir(f'{directory}figures-distribution-function')): 
        os.makedirs(f'{directory}figures-distribution-function')
 # Define the meta data for the movie
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='phase-space-animation', artist='Matplotlib',
                comment='visualizing the phase space evolution of the distribution function') # Describe what the animation is
    writer = FFMpegWriter(fps= n, metadata=metadata) # you can adjust the fps here.

# Initialize the movie
    fig = plt.figure()

    plt.xlabel(r'$x [c/\omega_{pe}]$') # Might need to rework this if the axes aren't static
    plt.ylabel(r'$v/c$')

    plt.tight_layout()
# Update the frames for the movie
    with writer.saving(fig, f"{directory}figures-distribution-function/distribution-function-movie.mp4", dpi=400):
        for i in range(len(p1x1_files)):
            p1x1 = vysxd_get_data(f'{directory}MS/PHA/p1x1/electrons/' + p1x1_files[i]) # Pull the phase space data
            plt.clf() # This clears the figure so you don't make a million colorbars lol
            plt.title(r'$t = ' + str(p1x1.TIME[0])+'\omega_{pe}^{-1}$')
            plt.imshow(np.log(-p1x1.DATA), origin='lower', extent=[p1x1.X[0], p1x1.X[-1],p1x1.Y[0],p1x1.Y[-1]], aspect='auto', vmin = -6)
            plt.colorbar(label=r'$f_e(x,v,t)$')
            writer.grab_frame()
