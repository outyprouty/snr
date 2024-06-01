import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

from scipy.ndimage import gaussian_filter
from scipy.signal import resample

from glob import glob

from photutils.detection import DAOStarFinder

def extractRadialData(subFrame, xC, yC):
    #Get matrix of integer indices associated with subFrame
    y, x = np.indices((subFrame.shape))

    #Generate matrix of radius values
    r = np.sqrt((x - xC)**2 + (y - yC)**2)

    #Force integer values (np.sqrt gives floats)
    r = r.astype(int)

    #Generate a histogram of radius bin, each weighed by corresponding counts
    weightedRadiusHistogram = np.bincount(r.ravel(), weights=subFrame.ravel())
    unweightedRadiusHistogram = np.bincount(r.ravel())

    #Get average for each radius bin
    averageCountsPerRadiusBin = weightedRadiusHistogram / unweightedRadiusHistogram
    return averageCountsPerRadiusBin




fileLocation = "Z:\\ObservatoryData\\Sharpcap\\2024-05-31"

directories = glob(fileLocation+"\\*")

darks = fileLocation+"\\Darks_150s\\*\\*.fits"

flatdarks = fileLocation+"\\Darks_25ms\\*\\*.fits"

flats = fileLocation+"\\flats_v_25ms\\*\\*.fits"

lights = fileLocation+"\\izar_v\\*\\*.fits"

# Import all calibration frames and light frames into their own arrays

darkData = [fits.open(d)[0].data for d in glob(darks)]
flatDarkData = [fits.open(d)[0].data for d in glob(flatdarks)]
flatData = [fits.open(d)[0].data for d in glob(flats)]
lightData = [fits.open(d)[0].data for d in glob(lights)]

# Show/Save example of each type
#plt.imshow(lightData[0])
#plt.show()

# Perform Calibration

## Read in all flat darks
#DONE

## Average all flat darks into master_flatdark
master_flatDark = np.mean(flatDarkData, axis=0)


## Read in all flats
#DONE

## subtract master_faltdark from each flat
## Average result into master_flat * Ct
master_flat = np.mean(flatData-master_flatDark, axis=0)

## Identify the constant Ct
Ct = np.mean(master_flat)

## Divide by Ct to get master_flat
master_flat /= Ct



## Read in all darks
#DONE

## average all darks into master_dark
master_dark = np.mean(darkData, axis=0)

## read in all lights
#DONE

## subtract master dark from each light and divide result by master flat
## average results into Science Frame

scienceFrame = np.mean((lightData-master_dark)/master_flat, axis=0)

meanLight = np.mean(scienceFrame)
stdLight = np.std(scienceFrame)

## Find sources
finder = DAOStarFinder(threshold=meanLight, fwhm=30, exclude_border=True, peakmax=50000, brightest=5)
sources = finder(scienceFrame)
xc = int(sources[0][1])
yc = int(sources[0][2])

print(xc, yc)

sfSize = 75

radius = 50
subFrame = scienceFrame[yc-sfSize:yc+sfSize,xc-sfSize:xc+sfSize]
# plt.imshow(subFrame, vmin=meanLight-stdLight, vmax=meanLight+stdLight, cmap='gray')
# plt.plot(sfSize,sfSize,marker='x',ms=10, color='r')

# plt.figure()
# plt.plot(subFrame[50,:], '.')
# plt.show()


Y, X = np.ogrid[:sfSize*2, :sfSize*2]

dist = np.sqrt((X-sfSize)**2 + (Y-sfSize)**2)
ones = np.ones((sfSize*2, sfSize*2))
counts = subFrame[dist<radius]

radial_data_raw = extractRadialData(subFrame, xC=sfSize, yC=sfSize)[:sfSize]

nPix = np.sum(ones, where=dist<radius)

print(f"So={np.sum(counts):0.2f}")
print(f"nPix={nPix}")
print(f"t=0.15s")

darkcounts = (master_dark[yc-sfSize:yc+sfSize,xc-sfSize:xc+sfSize])[dist<radius]

print(f"Sd={np.sum(darkcounts):0.2f}")

skyCounts = (scienceFrame[900-sfSize:900+sfSize,1100-sfSize:1100+sfSize])[dist<radius]


print(f"Ss={np.sum(skyCounts):0.2f}")


## Use Science Frame to identify sky Signal

## Use Science Frame to identify Source Signal

## Use master Dark to identify Dark Signal

"""
So=7567105.99
nPix=7825.0  
t=0.15
Sd=135563.20 
Ss=3322.33

>>> from numpy import sqrt 
>>> Q = 0.79              
>>> SNR = So*Q*t/sqrt(So*Q*t + Ss*Q*t + Sd*Q*t) 
>>> SNR
938.3713267294833
>>> SNR = So*t/sqrt(So*t + Ss*t + Sd*t)         
>>> SNR
1055.7502262792425
>>> SNR = So/sqrt(So + Ss + Sd)                 
>>> SNR
2725.9353627561454
>>>

"""

