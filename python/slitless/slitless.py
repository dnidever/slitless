#!/usr/bin/env python

import os
import time
import numpy as np
from dlnpyutils import utils as dln

def slitless():
     
    # Synthetic slitless spectral images 
     
     
    # -catalog of sources with 
    #  -X/Y position in the image 
    #  -magnitude 
    #  -stellar parameters (Teff, logg, [Fe/H]) 
    # 
    # -spatial fwhm 
    # -transmission curve 
    # -spectral resolution 
    # -orientation anges 
     
    # background, assume it is very low for the moment (this is mostly 
    #   true for space data) 
     
    # make two slitless spectra images at two different angles 
    # and one direct imaging image 
     
    # Parameters 
    fwhm = 2.5 
    dispersion = 1.0  # A per pixel 
    angle = 90.0  # 0.0 
    wr = [10000,11000] 
    w0 = np.mean(wr) 
    nx = 2048 
    ny = 2048 
     
    # Transmission/throughput curve
    trans = np.zeros(wr[1]-wr[0]+1,dtype=np.dtype([('wave',float),('flux',float)]))
    trans['wave'] = np.arange(wr[1]-wr[0]+1)+wr[0] 
    # -97.389234    0.0056913391   1.5786481e-06  -1.1699067e-10 
    tcoef = [-97.389234,    0.0056913391,   1.5786481e-06,  -1.1699067e-10] 
    trans['flux'] = np.polyval(np.flip(tcoef)trans['wave'])
     
    # Star catalog
    dt = [('id',int),('x',float),('y',float),('flux',float),('teff',float),('logg',float),('feh',float)]
    nstars = 500
    tab = np.zeros(nstars,dtype=np.dtype(dt))
    #tab[0] = {x:600.,y:400.,flux:1000.0,teff:3750.,logg:3.0,feh:-1.0} 
    #tab[1] = {x:700.,y:1500.,flux:10000.0,teff:3500.,logg:2.5,feh:-0.5} 
    #tab[2] = {x:1200.,y:700.,flux:5000.,teff:4000.,logg:2.5,feh:-0.5} 
    #tab[3] = {x:1300.,y:1400.,flux:2000.,teff:4500.,logg:3.0,feh:0.0} 
    tab['x'] = np.random.rand(nstars)*(nx-200)+100 
    tab['x'] = np.round(tab['x']).astype(int)
    tab['y'] = np.random.rand(nstars)*(ny-200)+100 
    tab['y'] = np.round(tab['y']).astype(int)
    tab['flux'] = np.random.rand(nstars)*1e4+200 
    teff = np.array([3500,3750,4000,4500]) 
    logg = np.array([2.5,3.0,2.5,3.0]))
    feh = np.array([-1.0,-0.5,-0.5,0.0])
    tab['teff'] = teff[np.round(np.random.rand(nstars)*3).astype(int)]
    tab['logg'] = logg[np.round(np.random.rand(nstars)*3).astype(int)]
    tab['feh'] = feh[np.round(np.random.rand(nstars)*3).astype(int)]
     
    # Initialize the image 
    im = np.zeros((nx,ny),float)
     
    # Synthetic spectra 
    specdir = '/Users/nidever/synspec/winter2017/grid2/' 
     
    # Loop over the stars 
    for i in range(nstars): 
        tab1 = tab[i] 
         
        # NO dispersion 
        #-------------- 
        If dispersion == 0.0: 
             
            # Generate the spatial Gaussian (FWHM) kernel 
            xr = [ np.maximum(np.floor(tab1['x']-5*fwhm),0), np.minimum(np.ceil(tab1['x']+5*fwhm),(nx-1)) ] 
            yr = [ np.maximum(np.floor(tab1['y']-5*fwhm),0), np.minimum(np.ceil(tab1['y']+5*fwhm),(ny-1)) ] 
            xx = (np.arange(xr[1]-xr[0]+1).astype(float)+xr[0]).reshape(-1,1) + np.zeros(yr[1]-yr[0]+1).reshape(-1,1)
            yy = np.zeros(xr[1]-xr[0]+1).reshape(-1,1) + (np.arange(yr[1]-yr[0]+1).astype(float)+yr[0]).reshape(-1,1)
            subim = exp(-0.5*((xx-tab1['x'])**2+(yy-tab1['y'])**2)/(fwhm/2.35)**2) 
            # Scale by flux 
            subim /= np.sum(subim) 
            subim *= tab1.flux 
            # Add to image 
            im[xr[0]:xr[1],yr[0]:yr[1]] += subim 
             
        # Dispersion 
        #----------- 
        else: 
             
            # Load the synthetic spectrum for this star 
            if tab1['feh'] >= 0.0: 
                msign='+' 
            else: 
                msign='-' 
            synfile = specdir+'spec.t'+string(tab1['teff'],format='(I4)')+'g'+string(tab1['logg'],format='(F4.2)')+
                  'm'+msign+string(abs(tab1['feh']),format='(F4.2)')+'a+0.00.fits' 
            if os.path.exists(synfile)==False:
                print(synfile+' NOT FOUND')
                import pdb; pdb.set_trace()
            synstr = Table.read(synfile)
            print(i+1,' ',synfile)
             
            # Scale the spectrum by the magnitude/flux 
            dum = dln.closest(w0,synstr.wave,ind=wind) 
            synstr['flux'] *= tab1['flux']/synstr['flux'][wind]
             
            # Generate the spatial Gaussian (FWHM) kernel 
            nkernel = np.ceil(fwhm*2.5) 
            if nkernel % 2 == 0 : 
                nkernel += 1 
            xkernel = np.arange(nkernel)-nkernel/2 
            kernel = np.exp(-0.5*xkernel**2/(fwhm/2.35)**2) 
            kernel /= np.sum(kernel)  # normalize 
             
            if angle == 0.0: 
                # Figure out the X/Y center and how far the spectrum stretches 
                # Traces are STRAIGHT (for now) 
                xr = (wr-w0)/dispersion + tab1['x'] 
                nxspec = xr[1]-xr[0]+1 
                wave = np.arange(nxspec)*dispersion+wr[0] 
                 
                # Bin the spectrum for the pixels 
                osamp = 4 
                xx = np.arange(nxspec*osamp).astype(float)/osamp+xr[0] 
                wavefine = np.arange(nxspec*osamp).astype(float)/osamp*dispersion+wr[0] 
                fluxfine = interp1d(wavefine,synstr['wave'],synstr.['flux']) 
                # multiply times the transmission 
                transfine = interp(wavefine,trans['wave'],trans['flux']) 
                binflux = dln.rebin(fluxfine*transfine,nxspec) 
                binxx = dln.rebin(xx,nxspec) 
                 
                # Spread the spectrum by the appropriate spatial Gaussian (FWHM) 
                binflux2d = binflux.reshape(1,-1) * kernel.reshape(-1,1)
                 
                # Trim overlap 
                if xr[0]<0: 
                    gd, = np.where(binxx >= 0) 
                    binflux2d = binflux2d[gd,:] 
                    xr[0] = 0 
                if xr[1]> nx-1: 
                    gd, = np.where(binxx >= 0) 
                    binflux2d = binflux2d[gd,:] 
                    xr[1] = nx-1 
                 
                # Add it to the image 
                yr = [-nkernel/2,nkernel/2]+tab1['y'] 
                im[xr[0]:xr[1],yr[0]:yr[1]] += binflux2d 
             
            elif angle == 90.0: 
                # Figure out the X/Y center and how far the spectrum stretches 
                # Traces are STRAIGHT (for now) 
                yr = (wr-w0)/dispersion + tab1['y'] 
                nyspec = yr[1]-yr[0]+1 
                wave = np.arange(nyspec)*dispersion+wr[0] 
                 
                # Bin the spectrum for the pixels 
                osamp = 4 
                yy = np.arange(nyspec*osamp).astype(float)/osamp+yr[0] 
                wavefine = np.arange(nyspec*osamp).astype(float)/osamp*dispersion+wr[0] 
                fluxfine = interp(wavefine,synstr['wave'],synstr['flux']) 
                # multiply times the transmission 
                transfine = interp(wavefine,trans['wave'],trans['flux']) 
                binflux = dln.rebin(fluxfine*transfine,nyspec) 
                binyy = dln.rebin(yy,nyspec) 
                 
                # Spread the spectrum by the appropriate spatial Gaussian (FWHM)
                binflux2d = binflux.reshape(1,-1) * kernel.reshape(-1,1)                  
                 
                # Trim overlap 
                if yr[0]<0: 
                    gd, = np.where(binyy >= 0) 
                    binflux2d = binflux2d[gd,:] 
                    yr[0] = 0 
                if yr[1]>ny-1: 
                    gd, = np.where(binyy >= 0) 
                    binflux2d = binflux2d[gd,:] 
                    yr[1] = ny-1 
                 
                # Add it to the image 
                xr = [-nkernel/2,nkernel/2]+tab1['x']
                im[xr[0]:xr[1],yr[0]:yr[1]] += binflux2d.T
             
    import pdb; pdb.set_trace() 
     
    return im
 
