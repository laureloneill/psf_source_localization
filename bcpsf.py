import numpy as np
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
import os
import scipy.stats as stats
import json
import pandas as pd
from xraysky.scripts import events_imaging


def get_peaksig_json(json_path):
    with open(json_path, 'r') as file:
        json_data = json.load(file)
    
    peaksig = json_data["peaksig"]

    return peaksig

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

def distFromCentroid(observation,plotHist=False):
    # computes the centroid of a finite set of points, and returns the distance that each point is from that centroid
    # observation -> pandas data frame consisting of the data from fitting a gaussian 
    
    xCoord = observation['xCenter'].to_numpy(dtype = np.float32)
    yCoord = observation['yCenter'].to_numpy(dtype = np.float32)

    centroidX = np.mean(xCoord)
    centroidY= np.mean(yCoord)
    dist = np.sqrt((centroidX-xCoord)**2+(centroidY-yCoord)**2)
    if plotHist == True:
        plt.hist(dist)
    return dist

def analyze_localization_directory(data_dir,window=10,plot = False ):
    data_dir = Path(data_dir)
    # take the names of all folders within the specified directory
    file_list = [file for file in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, file))]
    
    theta = [None]*len(file_list)
    xMean = [None]*len(file_list)
    yMean =[None]*len(file_list)
    xSigma = [None]*len(file_list)
    ySigma = [None]*len(file_list)
    Amplitude = [None]*len(file_list)
    peakSig = [None]*len(file_list)
    cov = [None]*len(file_list)
    validPeak = [None]*len(file_list)
    pval = [np.nan]*len(file_list) 
    

    for a in range(len(file_list)): 
        #
        #
        # numTerms = len(file_list[a].split('_'))
        #if numTerms == 7:
        #    run[a], anode[a], detectorMode[a], temp[a], sweep[a], vTheta[a], hTheta[a] =  file_list[a].split('_')

        #elif numTerms == 6:
        #    run[a], anode[a], detectorMode[a], temp[a], sweep[a], vTheta[a] = file_list[a].split('_')
        #    vTheta[a] = 'NA'
        #    hTheta[a] = 'NA'
        #elif numTerms == 4:
        #    run[a], anode[a], detectorMode[a], temp[a] = file_list[a].split('_') 
        #elif numTerms == 3:
        #    run[a],vTheta[a], hTheta[a] = file_list[a].split('_')
        #
        if os.path.isfile(f"{data_dir}/{file_list[a]}"):
            data_path = Path(f"{data_dir}/{file_list[a]}")
        try:
            d = fits.open(data_path) # open fits file
        
            data = d[0].data # data contents of the fits file
            #shape =  data.shape

            ny, nx = data.shape
            ym, xm = np.mgrid[:ny, :nx]

            w = WCS(d[0].header)

            x,y= w.array_index_to_world_values(ym,xm)
            x = (x+(180+360)) % 360 -180

            peak = np.max(data)
            peak_loc = np.unravel_index(np.argmax(data),data.shape)
            window = 10

            windowed_data = data[peak_loc[0]-window:peak_loc[0]+window, peak_loc[1]-window:peak_loc[1]+window]
            x_window = x[peak_loc[0]-window:peak_loc[0]+window, peak_loc[1]-window:peak_loc[1]+window]
            y_window = y[peak_loc[0]-window:peak_loc[0]+window, peak_loc[1]-window:peak_loc[1]+window]

            fitted,covariance_matrix = fit_rotated_2d_gaussian(windowed_data, x_window, y_window, plot_result=plot,iter=file_list[a])
        
            validPeak[a] = True
            xMean[a] = fitted.x_mean.value
            yMean[a] = fitted.y_mean.value
            xSigma[a] = fitted.x_stddev.value
            ySigma[a] = fitted.y_stddev.value
            Amplitude[a] = fitted.amplitude.value
            theta[a] =fitted.theta.value
            cov[a] = covariance_matrix
        except:
            validPeak[a] = "corrupt Fits"
        else:
            validPeak[a] = False  
        ######## Adjust whats in between the ########################################
        imager= events_imaging.BC_Imaging()
        peaks = imager.imager.findpeaks(data) # Index 5 has the significance
        if len(peaks) == 1: 
            peakSig = peaks[5]
        if len(peaks) == 2:
            peakSig = peaks[0]
            peakSig = peakSig[5]
        ###########################################################################
    result = {
        "xCenter" : xMean,
        "yCenter" : yMean,
        "xSigma" : xSigma,
        "ySigma" : ySigma,
        "peak" : Amplitude,
        "peakSig" : peakSig,
        "theta" : theta
        #'covMatrix' : cov
    }

    result = pd.DataFrame(result)
    return result, cov


def analyze_localization(eventlistin ,window=10,plot = False ):
    if type(eventlistin) == "str":
        eventlistin = fits.open(Path(eventlistin))
        data = d[0].data
        #w = WCS(d[0].header)

    if isinstance(eventlistin,(fits.HDUList)):
        data = eventlistin[0].data # data contents of the fits file
        #w = WCS(eventlistin[0].header)

    if isinstance(eventlistin,(fits.BinTableHDU)):
        data = eventlistin.data # data contents of the 
    

    ny, nx = data.shape
    ym, xm = np.mgrid[:ny, :nx]
    
    if isinstance(eventlistin,(fits.HDUList)):
        w = WCS(eventlistin[0].header)
    elif isinstance(eventlistin,(fits.BinTableHDU)):
        w = WCS(eventlistin.header)
    else:
        message = "Something about eventlistin is incorrect, couldnt assign w properly"
        print(message)
        return 
    
    x,y= w.array_index_to_world_values(ym,xm)
    x = (x+(180+360)) % 360 -180

    peak = np.max(data)
    peak_loc = np.unravel_index(np.argmax(data),data.shape)
    window = 10

    windowed_data = data[peak_loc[0]-window:peak_loc[0]+window, peak_loc[1]-window:peak_loc[1]+window]
    x_window = x[peak_loc[0]-window:peak_loc[0]+window, peak_loc[1]-window:peak_loc[1]+window]
    y_window = y[peak_loc[0]-window:peak_loc[0]+window, peak_loc[1]-window:peak_loc[1]+window]

    fitted,covariance_matrix = fit_rotated_2d_gaussian(windowed_data, x_window, y_window, plot_result=plot)

    validPeak = True
    xMean = fitted.x_mean.value
    yMean = fitted.y_mean.value
    xSigma = fitted.x_stddev.value
    ySigma = fitted.y_stddev.value
    Amplitude = fitted.amplitude.value
    theta =fitted.theta.value
    cov = covariance_matrix
    
    imager= events_imaging.BC_Imaging()
    peaks = imager.imager.findpeaks(data) # Index 5 has the significance
    if len(peaks) == 1: 
        peakSig = peaks[5]
    if len(peaks) == 2:
        peakSig = peaks[0]
        peakSig = peakSig[5]

    result = {
        "xCenter" : xMean,
        "yCenter" : yMean,
        "xSigma" : xSigma,
        "ySigma" : ySigma,
        "peak" : Amplitude,
        "peakSig" : peakSig,
        "theta" : theta,
        'covMatrix' : cov
    }

    #result = pd.DataFrame(result)
    return result