import numpy as np
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.io import fits
import sys, getopt,os



def fit_rotated_2d_gaussian(data, x=None, y=None, plot_result=False):
    """
    Fit a rotated 2D Gaussian to 2D data using Astropy's Gaussian2D model.

    Parameters:
        data : 2D numpy array
            The input data array to fit.
        x, y : 2D numpy arrays, optional
            Meshgrid coordinates corresponding to data. If None, will be auto-generated.
        plot_result : bool
            Whether to plot the original data and fit result.

    Returns:
        fitted_model : Gaussian2D
            The best-fit Gaussian2D model.
    """
    # Generate coordinate grids if not provided
    ny, nx = data.shape
    if x is None or y is None:
        y, x = np.mgrid[:ny, :nx]

    # Estimate initial parameters
    amplitude_init = np.max(data)
    x_mean_init = x[data == amplitude_init][0]
    y_mean_init = y[data == amplitude_init][0]

    # Initial guess for model
    gauss_init = models.Gaussian2D(amplitude=amplitude_init, x_mean=x_mean_init,
                                   y_mean=y_mean_init, x_stddev=5, y_stddev=5, theta=0)

    # Fitting with Levenberg-Marquardt algorithm
    fitter = fitting.LevMarLSQFitter()
    fitted_model = fitter(gauss_init, x, y, data)

    if plot_result:
        fit_data = fitted_model(x, y)
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        axes[0].imshow(data, origin='lower', cmap='viridis')
        axes[0].set_title("Original Data")
        axes[1].imshow(fit_data, origin='lower', cmap='viridis')
        axes[1].set_title("Fitted Gaussian")
        plt.tight_layout()
        plt.show()

    return fitted_model

############################################################################################

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

theta = [None]*len(folder_list)
xMean = [None]*len(folder_list)
yMean =[None]*len(folder_list)
xSigma = [None]*len(folder_list)
ySigma = [None]*len(folder_list)
Amplitude = [None]*len(folder_list)


window = 100 #  number of pixels to either side of the peak
imPixelSizeX = 0.0149 # pixel size degrees
imPixelSizeY = 0.0149 # pixel size degrees

for a in range(len(folder_list)): 
    run[a], anode[a], detectorMode[a], temp[a], sweep[a], vTheta[a], hTheta[a] =  folder_list[a].split('_')
    vTheta[a] = vTheta[a].replace('N','-')
    hTheta[a] = hTheta[a].replace("N","-")
    
    data_path = Path(f"c:/Users/ajo5182/Documents/astro/y2024-12-09/{folder_list[a]}/Analysis/imaging_analysis_20241213/image_reconstruction.fits.gz")
    d = fits.open(data_path) # open fits file
    data = d[0].data # data contents of the fits file
    shape =  data.shape
    
    y = imPixelSizeY * np.linspace((shape[0]-1)/-2,(shape[0]-1)/2, shape[0]) # create array from -23 deg to 23 deg, centered on zero
    x = imPixelSizeX * np.linspace((shape[1]-1)/-2,(shape[1]-1)/2, shape[1]) # create array from -41 deg to 41 deg, centered on zero

    peak = np.max(data)
    peak_loc = np.unravel_index(np.argmax(data),data.shape)
    
    if peak_loc[1]-window < 0:
        x_window = 1* x[0:peak_loc[1]+window]
    
    else:
        x_window = 1* x[peak_loc[1]-window:peak_loc[1]+window]


    y_window =y[peak_loc[0]-window:peak_loc[0]+window]
    
    if y_window.shape > x_window.shape :
        y_window = np.delete(y_window,0)

    X,Y = np.meshgrid(x_window,y_window)
    windowed_data = data[peak_loc[0]-window:peak_loc[0]+window, peak_loc[1]-window:peak_loc[1]+window]

    plt.contourf(X, Y, windowed_data, cmap=plt.cm.gist_earth_r)

    fitted = fit_rotated_2d_gaussian(windowed_data, X, Y, plot_result=False)

    xMean[a] = fitted.x_mean.value
    yMean[a] = fitted.y_mean.value
    xSigma[a] = fitted.x_stddev.value
    ySigma[a] = fitted.y_stddev.value
    Amplitude[a] = fitted.amplitude.value
    theta[a] =fitted.theta.value

np.savez("2dfit",run, anode, detectorMode, temp, sweep, hTheta, vTheta, xMean, xSigma, yMean, ySigma, Amplitude, theta)

    


   

print('done')