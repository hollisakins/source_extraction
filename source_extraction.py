from astropy.io import fits
from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os 
import warnings 
import sep 

warnings.catch_warnings() 
warnings.simplefilter("ignore")

class Field(object):
    def __init__(self,filename):
        self.filename = filename
        print('*   Processing %s' % filename)

    def openFits(self):
        with fits.open(self.filename) as hdulist:
            hdr = hdulist[0].header
            img = hdulist[0].data
            img = np.array(img,dtype='<f4')
        return hdr,img

    def Source(self): # gathers source extraction data from .SRC file
        src = np.loadtxt(self.filename.replace('.fits','.SRC'))
        objects = src[:,0:2] # pixel X,Y coordinates of the objects in question
        X_pos = src[:,0]
        Y_pos = src[:,1]
        A = src[:,8]
        B = src[:,9]
        theta = src[:,10]
        print('    Gathered source data for %s' % self.filename)
        return {'obj':objects,'X':X_pos,'Y':Y_pos,'A':A,'B':B,'theta':theta}

    def Convert(self): # converts obj list in pixel coordinate to RA-dec coordinates
        hdr, img = self.openFits()
        w = wcs.WCS(hdr)
        objects = self.source['obj']
        world = w.wcs_pix2world(objects, 1)
        print('    Converted coordinates to RA/Dec for %s' % self.filename)
        return world

    def Photometry(self):
        hdr, img = self.openFits()
        bkg = sep.Background(img)
        img_sub = img - bkg
        flux, fluxerr, flag = sep.sum_circle(img_sub, self.source['X'], self.source['Y'], 10.0, err=bkg.globalrms, gain=1.0)
        magnitudes = -2.5*np.log(flux)

        # median_i = np.median(magnitudes)
        # median_c = np.median(catalog)
        # dif = median_c - median_i
        # magnitudes += dif
        return magnitudes

    def Plot(self):
        hdr, img = self.openFits()
        proj = wcs.WCS(hdr)
        fig = plt.figure(figsize=(13,10)) 
        ax = fig.add_subplot(111,projection=proj)
        m, s = np.mean(img), np.std(img)
        im = ax.imshow(img, interpolation='nearest', cmap='gray',
                    vmin=m-s, vmax=m+s, origin='lower')

        overlay = ax.get_coords_overlay('fk5')
        overlay.grid(color='white', ls='dotted')
        overlay[0].set_axislabel('Right Ascension (J2000)')
        overlay[1].set_axislabel('Declination (J2000)')

        # plot an ellipse for each object
        for i in range(len(self.source['X'])):
            e = Ellipse(xy=(self.source['X'][i], self.source['Y'][i]),
                        width=6*self.source['A'][i],
                        height=6*self.source['B'][i],
                        angle=self.source['theta'][i])
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        name = self.filename.replace('Calibrated Images/','TestImages/')
        name = name.replace('.fits','.jpg')
        plt.savefig(name)
        print('    Created plot of %s' % self.filename)

    def Extract(self):
        self.source = self.Source()
        world = self.Convert()
        magnitudes = self.Photometry()
        return {'obj':world,'mag':magnitudes}
        
path_to_files = 'Calibrated Images/'
list_of_files = [f for f in os.listdir(path_to_files) if os.path.isfile(os.path.join(path_to_files,f)) and not f.startswith('.') and f.endswith('.fits')]

for filename in list_of_files:
    name = path_to_files+filename
    f = Field(name)
    output = f.Extract()
    for i in range(len(output['obj'])):
        print('RA/Dec: %s, Instrumental Mag: %s' %(output['obj'][i],output['mag'][i]))
    f.Plot()
    # end with list of WCS coordinates and calibrated magnitudes