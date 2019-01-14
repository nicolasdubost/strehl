import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# from matplotlib import pyplot as pl
# from matplotlib.pyplot import rc
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text',usetex=True)
# rc('image',interpolation='nearest')
from numpy.fft import fft2,ifft2,fftshift,ifftshift
from scipy.optimize import curve_fit
# import FITS

'''
Gets Strehl from a 2D image. Design to be used in darcplot

import liveStrehl as ss
strehl = ss.getStrehl(data)
text = ("%.2f"%strehl,750,500,"white",2)
'''

## USER INPUT -----------
wvl = 655e-9        # wavelength in meters
wvl = 633e-9        # wavelength in meters
f = 225e-3          # EFL of nsfi camera
d0 = 10.54e-3       # pupil diamter in pupil
up = 6.45e-6        # camera pixel size in meters
nd = 60             # pupil size in array cells
xpxls = 1388        # pixels across camera
ypxls = 1029
y0 = 500            # top left corner of search window for rpsf.max()
x0 = 750
yw = 50             # width for search window
xw = 50
## ----------------------
def getStrehl(image,calibrated=False):
    '''
    full 2D image
    use calibrated==True when Background has been subtracted from image
    use calibrated==False when image are raw pixles
    '''
    # image = image.reshape((ypxls,xpxls)) # Real PSF
    ym,xm=np.unravel_index(image[y0:y0+yw,x0:x0+xw].argmax(),(yw,xw))
    print 'ym,xm = ',(ym,xm)
    ymax,xmax = (y0+ym,x0+xm)
    print 'ymax,xmax = ',(ymax,xmax)
    rpsf = image[ymax-nt/2:ymax+nt/2,xmax-nt/2:xmax+nt/2].astype(np.float32)
    if not calibrated:
        rpsf -= rpsf[:10,:10].mean()
    cpsf = centrePSF(rpsf)
    # import pylab as pl
    # pl.imshow(cpsf,interpolation='nearest')

    # Strehl calculation
    denominator = rpsf.sum()/tsum
    cptp = cpsf.max()
    print 'rpsf[:10,:10].mean() = ',rpsf[:10,:10].mean()
    cnumerator = cptp+0.
    return cnumerator/denominator


def centrePSF(psf):
    '''
    only works on a psf already cropped around it's center.
    in other words, this centering only applies a +/- 1/2 pixel shift
    '''
    yn,xn = psf.shape
    print psf.shape
    fpsf = fft2(ifftshift(psf))
    p0 = [0,0]
    sigma = 1/(np.abs(fpsf)+fpsf[0,0].real*1e-5)
    bounds = ([-np.pi/xn,-np.pi/yn],[np.pi/xn,np.pi/yn])
    if 'centrePSFxdata' not in globals():
        xlader = np.arange(xn)-xn/2
        ylader = np.arange(yn)-yn/2
        xg,yg = np.meshgrid(xlader,ylader)
        global centrePSFxdata
        centrePSFxdata = (ifftshift(xg),ifftshift(yg))
    ydata = np.angle(fpsf)
    popt,pcov = curve_fit(tiptilt,centrePSFxdata,ydata.ravel(),p0=p0,sigma=sigma.ravel())
    # popt,pcov = curve_fit(tiptilt,centrePSFxdata,ydata.ravel(),p0=p0,sigma=sigma.ravel(),bounds=bounds)
    nfpsf = fpsf*np.exp(-1j*tiptilt(centrePSFxdata,tiltx=popt[0],tilty=popt[1]).reshape((yn,xn)))
    # print 'xtilt = %.2fpi'%(popt[0]*xn/np.pi)
    # print 'ytilt = %.2fpi'%(popt[1]*yn/np.pi)
    return fftshift(ifft2(nfpsf).real)
def tiptilt((xg,yg),tiltx=0,tilty=0):
    return (xg*tiltx + yg*tilty).ravel()


## Theoretical PSF ------------------
p = wvl*f/(d0*up)       # padding to theor. PSF has same pixel scale as in cam images
nt = int(round(p*nd))   # total number of pixels
lader = np.arange(nt)-(nt-1)/2.
xg,yg = np.meshgrid(lader,lader)
xtilt_test = xg.astype(np.float64)/nt
fix = -1    # to make fitted and marechal's strehl match
pupil = np.sqrt(np.square(xg)+np.square(yg))<nd/2.+fix
tpsf = fftshift((np.abs(fft2(pupil))**2))
norm = tpsf.max()
tpsf = tpsf/norm
tsum = tpsf.sum()
## ----------------------------------
