import shutil
from pathlib import Path 
from astropy.io import fits
from astropy.table import Table
import sys, getopt,os
import numpy as np


#def fitsmove(rawFitsDirectory):
#    directory = Path(rawFitsDirectory)

def stackFits(fitsDirectory):
    data_dir = Path(fitsDirectory)

    det_sns= ['s23056','s23196','s23197','s23200']

    for a in range(len(det_sns)):
        file_list = [file for file in os.listdir(f"{data_dir}/{det_sns[a]}") if os.path.isfile(os.path.join(f"{data_dir}/{det_sns[a]}", file))] 
        blank = np.zeros((550,550))
        
        #outfits = fits.PrimaryHDU(data=blank)
        
        for b in range(len(file_list)):
            path = Path(str(data_dir)+'/'+det_sns[a] +'/'+str(file_list[b]))
            d = fits.open(path)
            data = d[0].data # data contents of the fits file
            outFileName = data_dir/Path(Path(file_list[b]).stem).stem / "_combined"
            fits.append(outFileName,data)

            print('done')


#stackFits("Z:\Astro_BlackCAT\BlackCAT_Calibration_Data\LC_Calibration_Data\Raw Long Cell Data\y2024-12-09\BC001_Al_FF_243K_HORI_0_0")

def evntlist_in_browser(event_list):
    if type(event_list) == 'pathlib.PosixPath' or 'pathlib.WindowsPath':
        event_list = fits.open(event_list, memmap=True)
        event_data = Table(event_list[1].data)
        event_data.show_in_browser()
    elif type(event_list) != 'astropy.io.fits.fitsrec.FITS_rec':  
        event_data = Table(event_list[1].data)
        event_data.show_in_browser()

def random_sample_events_list(event_list, sampleNumber, seed, outname, save = True):
    # event_list: string or path to the event list you want to sample
    # sampleNumber: number of samples from the event list you want
    # outname: name to save the sub-sampled event list too
    # seed: seed for the random number generator
    # save: to save the file to the location in outname
    d = fits.open(event_list)
    data = d[1].data
    end = len(data)
    np.random.seed(seed)
    sample_index = np.random.randint(0,end,sampleNumber)
    d[1].data = d[1].data[sample_index]

    if save == True: 
        d.writeto(outname)
    
    return d