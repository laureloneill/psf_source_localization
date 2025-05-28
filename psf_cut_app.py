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

Xresult = [None]*len(folder_list)
xfwhm = [None]*len(folder_list)
xCenter =[None]*len(folder_list)
xSigma = [None]*len(folder_list)
xAmplitude = [None]*len(folder_list)
xHeight = [None]*len(folder_list)

Yresult = [None]*len(folder_list)
yfwhm = [None]*len(folder_list)
yCenter =[None]*len(folder_list)
ySigma = [None]*len(folder_list)
yAmplitude = [None]*len(folder_list)
yHeight = [None]*len(folder_list)

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
    
    fwhm = Xresult.params["fwhm"]
    xfwhm[a] = fwhm.value
    center = Xresult.params["center"]
    xCenter[a] =center.value
    sigma = Xresult.params["sigma"]
    xSigma[a] = sigma.value
    amplitude = Xresult.params["amplitude"]
    xAmplitude[a] = amplitude.value
    height = Xresult.params["height"]  
    xHeight[a] = height.value

    
np.savez("x_fits",run, anode, detectorMode,temp,sweep, hTheta,vTheta,xfwhm,xCenter,xSigma,xAmplitude,xHeight  )
print('done x ')



for a in range(len(folder_list)): 
  
    data_path = Path(f"c:/Users/ajo5182/Documents/astro/y2024-12-09/{folder_list[a]}/Analysis/imaging_analysis_20241213/image_reconstruction.fits.gz")
    d = fits.open(data_path) # open fits file
    data = d[0].data # data contents of the fits file
    shape =  data.shape
    
    ylab = imPixelSizeY * np.linspace((shape[0]-1)/-2,(shape[0]-1)/2, shape[0]) # create array from -23 deg to 23 deg, centered on zero
    xlab = imPixelSizeX * np.linspace((shape[1]-1)/-2,(shape[1]-1)/2, shape[1]) # create array from -41 deg to 41 deg, centered on zero

    peak = np.max(data)
    peak_loc = np.unravel_index(np.argmax(data),data.shape)
    
    data_windowed = data[peak_loc[0]-window:peak_loc[0]+window,peak_loc[1]-window:peak_loc[1]+window]

    xfit_ycut =ylab[peak_loc[0]-window:peak_loc[0]+window]

    if peak_loc[1]-window < 0:
        yfit_ycut = data[peak_loc[0],0:peak_loc[1]+window]
    else:
        yfit_ycut = data[peak_loc[0],peak_loc[1]-window:peak_loc[1]+window]

    if yfit_ycut.shape < xfit_ycut.shape :
        xfit_ycut = np.delete(xfit_ycut,0)

    #xfit_xcut = 1* xlab[peak_loc[1]-window:peak_loc[1]+window]
    #yfit_xcut = data[peak_loc[0]-window:peak_loc[0]+window,peak_loc[1]]
    
    Ymodel  = GaussianModel()
    params = Ymodel.guess(yfit_ycut, xfit_ycut)
    Yresult = Ymodel.fit(yfit_ycut, params,x = xfit_ycut)
    
    fwhm = Yresult.params["fwhm"]
    yfwhm[a] = fwhm.value
    center = Yresult.params["center"]
    yCenter[a] =center.value
    sigma = Yresult.params["sigma"]
    ySigma[a] = sigma.value
    amplitude = Yresult.params["amplitude"]
    yAmplitude[a] = amplitude.value
    height = Yresult.params["height"]  
    yHeight[a] = height.value

np.savez("y_fits",run, anode, detectorMode,temp,sweep, hTheta,vTheta,yfwhm,yCenter,ySigma,yAmplitude,yHeight  )
