import numpy as np
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.io import fits
import sys, getopt,os
import scipy.stats as stats




def fit_rotated_2d_gaussian(data, x=None, y=None, plot_result=False,iter = None):
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
                                   y_mean=y_mean_init, x_stddev=0.089, y_stddev=0.069, theta=0)

    # Fitting with Levenberg-Marquardt algorithm
    fitter = fitting.LevMarLSQFitter()
    fitted_model = fitter(gauss_init, x, y, data)
    covariance_matrix = fitter.fit_info['param_cov']
    if fitted_model.x_fwhm > 1 or fitted_model.y_fwhm >1 :
        print(iter)

    if plot_result:
        fit_data = fitted_model(x, y)
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        axes[0].imshow(data, origin='lower', cmap='viridis')
        axes[0].set_title("Original Data")
        axes[1].imshow(fit_data, origin='lower', cmap='viridis')
        axes[1].set_title("Fitted Gaussian")
        plt.tight_layout()
        plt.show()

    return fitted_model,covariance_matrix

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
cov = [None]*len(folder_list)
validPeak = [None]*len(folder_list)
pval = [np.nan]*len(folder_list)

window =10 #  number of pixels to either side of the peak
imPixelSizeX = 0.0149 # pixel size degrees
imPixelSizeY = 0.0149 # pixel size degrees

for a in range(len(folder_list)): 
    numTerms = len(folder_list[a].split('_'))
    if numTerms == 7:
        run[a], anode[a], detectorMode[a], temp[a], sweep[a], vTheta[a], hTheta[a] =  folder_list[a].split('_')
    #    vTheta[a] = vTheta[a].replace('N','-')
    #    hTheta[a] = hTheta[a].replace("N","-")
    elif numTerms == 6:
        run[a], anode[a], detectorMode[a], temp[a], sweep[a], vTheta[a] = folder_list[a].split('_')
        vTheta[a] = 'NA'
        hTheta[a] = 'NA'
    elif numTerms == 4:
        run[a], anode[a], detectorMode[a], temp[a] = folder_list[a].split('_') 
    elif numTerms == 3:
        run[a],vTheta[a], hTheta[a] = folder_list[a].split('_')  

    if os.path.isfile(f"{data_dir}/{folder_list[a]}/Analysis/imaging_analysis_20241213/Figures/imaging_source.png"):
    
        data_path = Path(f"{data_dir}/{folder_list[a]}/Analysis/imaging_analysis_20241213/image_reconstruction.fits.gz")
        try:
            d = fits.open(data_path) # open fits file
        
            data = d[0].data # data contents of the fits file
            shape =  data.shape
    
            y = imPixelSizeY * np.linspace((shape[0]-1)/-2,(shape[0]-1)/2, shape[0]) # create array from -23 deg to 23 deg, centered on zero
            x = imPixelSizeX * np.linspace((shape[1]-1)/-2,(shape[1]-1)/2, shape[1]) # create array from -41 deg to 41 deg, centered on zero

            peak = np.max(data)
            peak_loc = np.unravel_index(np.argmax(data),data.shape)

            x_window =x[peak_loc[1]-window:peak_loc[1]+window]
            y_window =y [peak_loc[0]-window:peak_loc[0]+window]
            X,Y = np.meshgrid(x[peak_loc[1]-window:peak_loc[1]+window],y[peak_loc[0]-window:peak_loc[0]+window])
            windowed_data = data[peak_loc[0]-window:peak_loc[0]+window, peak_loc[1]-window:peak_loc[1]+window]
            
            shapiro_result = stats.shapiro(windowed_data, axis=0)
            pval[a] = np.mean(shapiro_result.pvalue)

            '''
            if peak_loc[0]-window < 0:
                x_window = 1* x[0:peak_loc[0]+window]
            else:
                x_window = 1* x[peak_loc[0]-window:peak_loc[0]+window]
            if peak_loc[1]-window < 0:
                y_window = 1* x[0:peak_loc[1]+window]
            else:
                y_window = 1* x[peak_loc[1]-window:peak_loc[1]+window]

            #y_window =y[peak_loc[0]-window:peak_loc[0]+window]
    
            if y_window.shape < x_window.shape :
                x_window = np.delete(x_window,0)

            X,Y = np.meshgrid(x_window,y_window)

            #windowed_data = data[x_window, y_window]
            if peak_loc[1]-window < 0:
                lefty = 0
            else:
                lefty = peak_loc[1]-window

            windowed_data = data[peak_loc[0]-window:peak_loc[0]+window, lefty:peak_loc[1]+window]
            '''
            #plt.contourf(X, Y, windowed_data, cmap=plt.cm.gist_earth_r)
    
            fitted,covariance_matrix = fit_rotated_2d_gaussian(windowed_data, -X, Y, plot_result=False,iter=folder_list[a])

            validPeak[a] = True
            xMean[a] = fitted.x_mean.value
            yMean[a] = fitted.y_mean.value
            xSigma[a] = fitted.x_stddev.value
            ySigma[a] = fitted.y_stddev.value
            Amplitude[a] = fitted.amplitude.value
            theta[a] =fitted.theta.value
            cov[a] = covariance_matrix
            runout = run[a]+anode[a]+detectorMode[a]+sweep[a]+hTheta[a]+vTheta[a]
        except:
            validPeak[a] = "corrupt Fits"
    else:
        validPeak[a] = False  
  

np.savez("data/test",run=run, anode=anode, detectorMode=detectorMode, temp=temp, 
         sweep=sweep, hTheta=hTheta, vTheta=vTheta, 
         xMean=xMean, xSigma=xSigma, 
         yMean=yMean, ySigma=ySigma,
         Amplitude=Amplitude, theta=theta,
         validPeak=validPeak)
np.savez("data/testcov", np.array(cov, dtype=object),allow_pickle = True)
np.savez("data/test_pval", pval = pval)

    


   
print('\a')
print('done')