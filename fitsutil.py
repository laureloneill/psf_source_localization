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
    
    # given an event list this returns an event list with containing a given number
    # of subsamples from the event list
    
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

def subsample_eventlist(event_list, numberOfSamples, outname, outdir):
    # function to take a large event list and split it into smaller event lists of a defined number events, and write that to a directory
    # currently it saves it to a directory so that the current analysis pipeline for blackcat can be used, but you could just as easily output each list
    #
    # output is an eventlist with the name "outname{how many times it ran}.fits" in the directory specified in outdir

    # event_list -> path the the event list you want to split
    # numberOfSamples -> number of samples per new event list EG enter 50000 if you want 50000 events per file
    # outname -> basename for the files you'll get out
    # outdir -> directory that you will write the new fits files to, make sure it exists before you start the function

    # TO DO:
    #   1) have it create the directory if it doesn't exist already


    outdir = Path(outdir)
    if outdir.exists == False:
        os.mkdir(outdir)

    d = fits.open(event_list)
    data = d[1].data
    cols = d[1].columns
    size = len(data)

    subSampleHDU = fits.BinTableHDU.from_columns(cols)

    for n in range(int(size/numberOfSamples)):
        subSampleHDU.data = d[1].data[n*int(size/numberOfSamples):(n+1)*int(size/numberOfSamples)]
        fitsName = Path(outname+str(n)+'.fits')

        subSampleHDU.writeto(outdir/fitsName,overwrite=True)

        


