from astropy.io import fits
from astropy import wcs
from astroquery.vizier import Vizier 
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import math 
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os 
import warnings 
import sep 
import csv

warnings.catch_warnings() 
warnings.simplefilter('ignore')

def writeData(output):
    with open('sources.csv', 'a') as outfile:
        writer = csv.writer(outfile)
        global columnsWritten
        if columnsWritten==False:
            writer.writerow(output.keys())
            columnsWritten = True
        writer.writerows(zip(*output.values()))


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

        ### perform aperture photometry
        hdr, img = self.openFits()
        bkg = sep.Background(img)
        img_sub = img - bkg
        flux, fluxerr, flag = sep.sum_circle(img_sub, self.source['X'], self.source['Y'], 10.0, err=bkg.globalrms, gain=1.0)
        magnitudes = -2.5*np.log(flux)
        median_i = np.median(magnitudes)


        ### retrieve magnitudes from catalog
        time = hdr['DATE-OBS']
        time = datetime.strptime(time, '%Y-%m-%dT%H:%M:%S.%f')
        # filt = hdr['FILTER']
        filt = 'V'
        objects = self.world


        v = Vizier(columns=['UCAC4','+_r','RAJ2000','DEJ2000','Bmag','Vmag','rmag'])
        output = {'id':[],'RA_C':[],'DEC_C':[],'RA_M':[],'DEC_M':[],'DIF':[],'MAG_R':[],'MAG_V':[],'MAG_B':[],'CMAG_R':[],'CMAG_V':[],'CMAG_B':[],'DATETIME':[],'IMGNAME':[]}
        cmags = []
        misfires = 0
        
        for n in range(len(objects)):
            catalog = 'UCAC4'
            result = v.query_region(coord.SkyCoord(ra=objects[n,0], dec=objects[n,1],
            unit=(u.degree, u.degree), frame='fk5'),radius='2s',catalog=catalog)
            try:
                result = result[0]
            except:
                print('No star match within 2 arcseconds')
                misfires += 1 
                output['id'].append('nan')
                output['RA_C'].append('nan')
                output['DEC_C'].append('nan')
                output['RA_M'].append(objects[n,0])
                output['DEC_M'].append(objects[n,1])
                output['DIF'].append('nan')
                output['DATETIME'].append(time)
                output['IMGNAME'].append(self.filename)
                cmags.append('nan')
                continue

            ids = np.array(result['UCAC4'],str)
            ra = np.array(result['RAJ2000'],float)
            dec = np.array(result['DEJ2000'],float) # catalog RA and Dec
            dif = np.array(result['_r'],float)
            fluxtype = filt+'mag'
            if filt=='R':
                fluxtype = 'rmag'
            flux = np.array(result[fluxtype],float)

            for i in range(len(ids)):
                if dif[i] <= 2 and i==np.argmin(dif): # min residual value and less than 2 arcsec off
                    print('Star match in catalog %s, mag %s, residual %s arcsec' % (catalog,flux[i],dif[i]))
                    output['id'].append(ids[i])
                    output['RA_C'].append(ra[i])
                    output['DEC_C'].append(dec[i])
                    output['RA_M'].append(objects[n,0])
                    output['DEC_M'].append(objects[n,1])
                    output['DIF'].append(dif[i])
                    output['DATETIME'].append(time)
                    output['IMGNAME'].append(self.filename)
                    cmags.append(flux[i])

                    
                else:
                    print('No star match within 2 arcseconds')
                    misfires += 1
                    continue


        print('output %s stars' % len(output['id']))
        print('output %s unique stars' % len(set(output['id'])))
        print('Missed %s objects' % misfires)

        cmags_nonan = [k for k in cmags if not math.isnan(float(k))]
        median_c = np.median(np.array(cmags_nonan))

        d = median_c - median_i
        magnitudes += d
        for n in range(len(objects)):
            for j in ['R','V','B']:
                magtype = 'MAG_'+j
                if j==filt:
                    output[magtype].append(magnitudes[n])
                else:
                    output[magtype].append('---')
                magtype = 'CMAG_'+j
                if j==filt:
                    output[magtype].append(cmags[n])
                else:
                    output[magtype].append('---')

        return output

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
        self.world = self.Convert()
        out = self.Photometry()
        return out
        
path_to_files = 'Calibrated Images/'
list_of_files = [f for f in os.listdir(path_to_files) if os.path.isfile(os.path.join(path_to_files,f)) and not f.startswith('.') and f.endswith('.fits')]

columnsWritten = False

for filename in list_of_files:
    name = path_to_files+filename
    f = Field(name)
    output = f.Extract()
    writeData(output)
    # f.Plot()
    # end with list of WCS coordinates and calibrated magnitudes