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

folder_list = [folder for folder in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, folder))]
#dir_list = os.listdir(data_dir)
#print(folder_list)
run = np.zeros(len(folder_list))
for a in range(len(folder_list)):
    run[a], anode[a], runtype[a], temp[a], sweep[a], vTheta[a], hTheta[a] =  folder_list[a].split('_')

print(htheta)