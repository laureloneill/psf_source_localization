import shutil
from pathlib import Path 
from astropy.io import fits
from astropy.table import Table
import sys, getopt,os
import numpy as np
from xhcd.core import Spectrum


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

def random_sample_events_list(input, sampleNumber, outdir=None, eventlist=None, save = False):
    
    # given an event list this returns an event list with containing a given number
    # of subsamples from the event list
    
    # event_list: string or path to the event list you want to sample
    # sampleNumber: number of samples from the event list you want
    # oudir: directory to save the subsample to
    # eventlist: name to save the sub-sampled event list too
    # 
    # save: to save the file to the location in eventlist
    
    if type(input) == 'str':
        directory = Path(input)
        d = fits.open(directory/Path(eventlist))
    else: 
        d = input
    if type(d) == "astropy.io.fits.hdu.hdulist.HDUList":
        data = d[1].data
        cols = d[1].columns
    else: 
        data = d.data
        cols = d.columns
    
    end = len(data)
    if sampleNumber > end:
        sampleNumber = end
        print("number of samples requested is greater than number of events")
    population = range(0,end)
    sample_index = np.random.choice(population, sampleNumber,replace = False)

    subSampleHDU = fits.BinTableHDU.from_columns(cols)
    subSampleHDU.data = data[sample_index]

    if save == True:
        if os.path.isdir(Path(outdir)) == False  :
            os.makedirs(Path(outdir))
        subSampleHDU.writeto(outdir/Path(eventlist),overwrite=True)
    
    return subSampleHDU

def subsample_eventlist(data_dir, numberOfLists, outdir):
    # function to take a large event list and split it into a defined number of snaller event lists, and write that to a directory
    # currently it saves it to a directory so that the current analysis pipeline for blackcat can be used, but you could just as easily output each list
    #
    # output is an eventlist with the name "outname{how many times it ran}.fits" in the directory specified in outdir

    # data_dir -> path the the event list collection you want to split
    # numberOfLists -> number of subsets of the event list you want to make 
    # outdir -> directory that you will write the new fits files to, make sure it exists before you start the function

    # TO DO:
    #   1) have it create the directory if it doesn't exist already



    file_list = [file for file in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, file))] 

    for a in range(len(file_list)):
        d = fits.open(data_dir/Path(file_list[a]))
        data = d[1].data 
        cols = d[1].columns
        size = len(data)
        numberOfSamples = int(size/numberOfLists)
        subSampleHDU = fits.BinTableHDU.from_columns(cols)
        for n in range(numberOfLists):
            subSampleHDU.data = d[1].data[n*numberOfSamples:(n+1)*numberOfSamples]
            out = outdir/Path('subsample'+str(n))
            if os.path.isdir((out)) == False:
                os.makedirs(out)
                #print('made the directory')
            #else:
            #    print("didn't make the directory")

            subSampleHDU.writeto(out/Path(file_list[a]),overwrite=True)

    '''
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
    '''

def numEvents(dataPath):
    d = fits.open(Path(dataPath))
    size = len(d[1].data)
    return size  

def splitEventList(eventlistin,splitThreshold,parameter, save = False):
    #if type(eventlistin) == 'pathlib.PosixPath' or 'pathlib.WindowsPath':
    #    d = fits.open(eventlistin)

    if type(eventlistin) == "str":
        d = fits.open(Path(eventlistin))
   
    else:
        d =  eventlistin
    #if type(d) == 'astropy.io.fits.hdu.hdulist.HDUList':
    cols = d[1].columns
    aboveThreshold = fits.BinTableHDU.from_columns(cols)
    belowThreshold = fits.BinTableHDU.from_columns(cols)
    aboveThreshold.data = d[1].data[d[1].data[str(parameter)]>splitThreshold]
    belowThreshold.data = d[1].data[d[1].data[str(parameter)]<splitThreshold]
    #else: 
     #   cols = d.columns
     #   aboveThreshold.data = d.data[d.data[str(parameter)]>splitThreshold]
     #   belowThreshold.data = d.data[d.data[str(parameter)]<splitThreshold]
   
    aboveThreshold = fits.BinTableHDU.from_columns(cols)
    belowThreshold = fits.BinTableHDU.from_columns(cols)
    #aboveThreshold.data = d[1].data[d[1].data[str(parameter)]>splitThreshold]
    #belowThreshold.data = d[1].data[d[1].data[str(parameter)]<splitThreshold]

    return aboveThreshold, belowThreshold

def halfEventList(eventListIn):
    # splits an eventlist in half
    # only accepts hdu list or path to one

    if type(eventListIn)=="str":
        d = fits.open(eventListIn)
    else:
        d = eventListIn
    if type(d) == "astropy.io.fits.hdu.hdulist.HDUList":
        middle =int(len(d[1]))
        cols = d[1].columns
        firstHalfDat = fits.BinTableHDU.from_columns(cols)
        secondHalfDat = fits.BinTableHDU.from_columns(cols)
        firstHalfDat.data = d[1].data[0:middle]
        secondHalfDat.data = d[1].data[middle:len(d[1])]
    else:
        middle =int(len(d.data))
        cols = d.columns
        firstHalfDat = fits.BinTableHDU.from_columns(cols)
        secondHalfDat = fits.BinTableHDU.from_columns(cols)
        firstHalfDat.data = d.data[0:middle]
        secondHalfDat.data = d.data[middle:len(d.data)]
    
    firstHalfDat = fits.BinTableHDU.from_columns(cols)
    secondHalfDat = fits.BinTableHDU.from_columns(cols)

    


    return firstHalfDat, secondHalfDat