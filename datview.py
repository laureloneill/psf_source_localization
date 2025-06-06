
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from pathlib import Path
from astropy.io import fits
import numpy as np

#datPath = Path("C:/Users/ajo5182/Documents/astro/psf/psf_source_localization/data/2dfit_test.npz")
datPath = Path("C:/Users/ajo5182/Documents/astro/psf/psf_source_localization/data/2dfit_2024_12_09.npz")
data = np.load(datPath, allow_pickle=True)

xSigma = np.array(data['xSigma'])
ySigma = np.array(data['ySigma'])
angle = np.array(data['theta'])
xMean = np.array(data['xMean'])
yMean = np.array(data['yMean'])
valid_peak = np.array(data['validPeak'])
run_name = (np.array(data['run'])+'_'+np.array(data['anode'])+'_'+np.array(data['detectorMode'])+'_'+
            np.array(data['sweep'])+'_'+np.array(data["vTheta"])+'_'+np.array(data["hTheta"]))

xSigma = [0 if item is None else item for item in xSigma]
ySigma = [0 if item is None else item for item in ySigma]
xMean = [0 if item is None else item for item in xMean]
yMean = [0 if item is None else item for item in yMean]
angle = [0 if item is None else item for item in angle]

xfwhm = 2.35*np.array(xSigma)
#ySigma = [0 if item is None else item for item in ySigma]
yfwhm = 2.35*np.array(ySigma)

imPixelSizeX = 0.0149 # pixel size degrees
imPixelSizeY = 0.0149 # pixel size degrees
shapey = 3137
shapex = 5585
y = imPixelSizeY * np.linspace((shapey-1)/-2,(shapey-1)/2, shapey) # create array from -23 deg to 23 deg, centered on zero
x = imPixelSizeX * np.linspace((shapex-1)/-2,(shapex-1)/2, shapex) # create array from -41 deg to 41 deg, centered on zero

ells =  [Ellipse(xy= (1*xMean[i],-1*yMean[i]),width=3*xfwhm[i],height=3*yfwhm[i], angle= angle[i]) for i in range(len(xfwhm))]
i = 0
fig, ax = plt.subplots()
ax.set(xlim=(-41, 41), ylim=(-23,23 ), aspect="equal")
for e in ells[10:11]:
    #print(i)
    #ax.legend(e,run_name)
    #ax.add_artist(e)
    ax.add_patch(e)
    e.set(clip_box=ax.bbox, alpha=1, facecolor=np.random.rand(3), label=run_name[i])
    i = i+1
ax.legend()
plt.show()
