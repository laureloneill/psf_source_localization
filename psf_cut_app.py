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
import sys, getopt


arglist = sys.argv[1:]

data_dir = Path(arglist[0])
write_loc = Path(arglist[1])

folder_list = [folder for folder in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, folder))]