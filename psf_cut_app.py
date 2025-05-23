from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.utils.exceptions import AstropyUserWarning
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
import numpy as np
from lmfit.models import GaussianModel
from scipy.interpolate import griddata
import sys, getopt,os


arglist = sys.argv[1:]

data_dir = Path(arglist[0])
#write_loc = Path(arglist[1])

# take the names of all folders within the specified directory
folder_list = [folder for folder in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, folder))] 

# create an empty list of strings to hold the data parsed from the name of the directory
run = [None]*len(folder_list)
anode = [None]*len(folder_list)
detectorMode= [None]*len(folder_list)
temp = [None]*len(folder_list)
sweep = [None]*len(folder_list)
vTheta = [None]*len(folder_list)
hTheta = [None]*len(folder_list)


window = 100 #  number of pixels to either side of the peak
imPixelSizeX = 0.0149 # pixel size degrees
imPixelSizeY = 0.0149 # pixel size degrees
# take the folder name and parse it to keep track of run number, the anode used, the detector mode, what type of sweep
# and the commanded angular position of the gimbal




for a in range(len(folder_list)): 
    run[a], anode[a], detectorMode[a], temp[a], sweep[a], vTheta[a], hTheta[a] =  folder_list[a].split('_')
    vTheta[a] = vTheta[a].replace('N','-')
    hTheta[a] = hTheta[a].replace("N","-")
    
    data_path = Path(f"c:/Users/ajo5182/Documents/astro/y2024-12-09/{folder_list[a]}/Analysis/imaging_analysis_20241213/image_reconstruction.fits.gz")
    d = fits.open(data_path) # open fits file
    data = d[0].data # data contents of the fits file
    shape =  data.shape
    
    ylab = imPixelSizeY * np.linspace((shape[0]-1)/-2,(shape[0]-1)/2, shape[0]) # create array from -23 deg to 23 deg, centered on zero
    xlab = imPixelSizeX * np.linspace((shape[1]-1)/-2,(shape[1]-1)/2, shape[1]) # create array from -41 deg to 41 deg, centered on zero

    peak = np.max(data)
    peak_loc = np.unravel_index(np.argmax(data),data.shape)
    
    data_windowed = data[peak_loc[0]-window:peak_loc[0]+window,peak_loc[1]-window:peak_loc[1]+window]

    if peak_loc[1]-window < 0:
        xfit_xcut = 1* xlab[0:peak_loc[1]+window]
    
    else:
        xfit_xcut = 1* xlab[peak_loc[1]-window:peak_loc[1]+window]

    yfit_xcut = data[peak_loc[0]-window:peak_loc[0]+window,peak_loc[1]]

    if yfit_xcut.shape > xfit_xcut.shape :
        yfit_xcut = np.delete(yfit_xcut,0)

    #xfit_xcut = 1* xlab[peak_loc[1]-window:peak_loc[1]+window]
    #yfit_xcut = data[peak_loc[0]-window:peak_loc[0]+window,peak_loc[1]]
    
    Xmodel  = GaussianModel()
    params = Xmodel.guess(yfit_xcut, xfit_xcut)
    Xresult = Xmodel.fit(yfit_xcut, params,x = xfit_xcut)
    #print(Xresult.params['center'].value)

print('done')