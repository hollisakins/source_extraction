from astropy.io import fits
from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os
import warnings

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
        return hdr,img

    def Source(self): # gathers source extraction data from .SRC file
        src = np.loadtxt(self.filename.replace('.fits','.SRC'))
        objects = src[:,0:2] # pixel X,Y coordinates of the objects in question
        magnitudes = src[:,2]
        X_pos = src[:,0]
        Y_pos = src[:,1]
        A = src[:,8]
        B = src[:,9]
        theta = src[:,10]
        self.source = {'obj':objects,'mag':magnitudes,'X':X_pos,'Y':Y_pos,'A':A,'B':B,'theta':theta}
        print('    Gathered source data for %s' % self.filename)

    def Convert(self): # converts obj list in pixel coordinate to RA-dec coordinates
        self.Source()
        hdr, img = self.openFits()
        w = wcs.WCS(hdr)
        # w.wcs.print_contents()
        objects = self.source['obj']
        instrumental = self.source['mag']
        world = w.wcs_pix2world(objects, 1)
        ### add into this function the ability to take the "world" and "magnitudes" arrays 
        ### and calibrate to a star catalog, based on the median value
        magnitudes = instrumental
        # median_i = np.median(instrumental)
        # median_c = np.median(catalog)
        # dif = median_c - median_i
        # magnitudes = instrumental + dif
        print('    Converted coordinates to RA/Dex for %s' % self.filename)
        return {'obj':world,'mag':magnitudes}


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

path_to_files = 'Calibrated Images/'
list_of_files = [f for f in os.listdir(path_to_files) if os.path.isfile(os.path.join(path_to_files,f)) and not f.startswith('.') and f.endswith('.fits')]

for filename in list_of_files:
    name = path_to_files+filename
    f = Field(name)
    o = f.Convert()
    # for i in range(len(o['obj'])):
    #     print('RA/Dec: %s, Instrumental Mag: %s' %(o['obj'][i],o['mag'][i]))
    f.Plot()
    # end with list of WCS coordinates and calibrated magnitudes