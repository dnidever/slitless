#!/usr/bin/env python

import os
import time
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table
from dlnpyutils import utils as dln

def image(pars,startab):
    """ Make a regular imaging image, no dispersion"""
     
    # Initialize the image 
    im = np.zeros((pars['nx'],pars['ny']),float)
     
    # Loop over the stars 
    for i in range(len(startab)):
        tab1 = startab[i] 
        # Generate the spatial Gaussian (FWHM) kernel 
        xr = [ int(np.maximum(np.floor(tab1['x']-5*pars['fwhm']),0)),
               int(np.minimum(np.ceil(tab1['x']+5*pars['fwhm']),(pars['nx']-1))) ] 
        yr = [ int(np.maximum(np.floor(tab1['y']-5*pars['fwhm']),0)),
               int(np.minimum(np.ceil(tab1['y']+5*pars['fwhm']),(pars['ny']-1))) ]
        xx = (np.arange(xr[1]-xr[0]+1).astype(float)+xr[0]).reshape(1,-1) + np.zeros(yr[1]-yr[0]+1).reshape(-1,1)
        yy = np.zeros(xr[1]-xr[0]+1).reshape(1,-1) + (np.arange(yr[1]-yr[0]+1).astype(float)+yr[0]).reshape(-1,1)
        subim = np.exp(-0.5*((xx-tab1['x'])**2+(yy-tab1['y'])**2)/(pars['fwhm']/2.35)**2) 
        # Scale by flux 
        subim /= np.sum(subim) 
        subim *= tab1['flux'] 
        # Add to image
        im[xr[0]:xr[1]+1,yr[0]:yr[1]+1] += subim
            
    return im

            
def spectrum(pars,startab):
     
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
    #pars = {}
    #pars['fwhm'] = 2.5 
    #pars['dispersion'] = 1.0  # A per pixel 
    #pars['angle'] = 90.0  # 0.0 
    #pars['wr'] = [10000,11000] 
    #pars['w0'] = np.mean(pars['wr']) 
    #pars['nx'] = 2048 
    #pars['ny'] = 2048 
     
    # Transmission/throughput curve
    #trans = np.zeros(pars['wr'][1]-pars['wr'][0]+1,dtype=np.dtype([('wave',float),('flux',float)]))
    #trans['wave'] = np.arange(pars['wr'][1]-pars['wr'][0]+1)+wars['wr'][0] 
    ## -97.389234    0.0056913391   1.5786481e-06  -1.1699067e-10 
    #tcoef = [-97.389234,    0.0056913391,   1.5786481e-06,  -1.1699067e-10] 
    #trans['flux'] = np.polyval(np.flip(tcoef),trans['wave'])
    #pars['trans'] = trans
    
    # Star catalog
    #dt = [('id',int),('x',float),('y',float),('flux',float),('teff',float),('logg',float),('feh',float)]
    #nstars = 500
    #tab = np.zeros(nstars,dtype=np.dtype(dt))
    ##tab[0] = {x:600.,y:400.,flux:1000.0,teff:3750.,logg:3.0,feh:-1.0} 
    ##tab[1] = {x:700.,y:1500.,flux:10000.0,teff:3500.,logg:2.5,feh:-0.5} 
    ##tab[2] = {x:1200.,y:700.,flux:5000.,teff:4000.,logg:2.5,feh:-0.5} 
    ##tab[3] = {x:1300.,y:1400.,flux:2000.,teff:4500.,logg:3.0,feh:0.0} 
    #tab['x'] = np.random.rand(nstars)*(pars['nx']-200)+100 
    #tab['x'] = np.round(tab['x']).astype(int)
    #tab['y'] = np.random.rand(nstars)*(pars['ny']-200)+100 
    #tab['y'] = np.round(tab['y']).astype(int)
    #tab['flux'] = np.random.rand(nstars)*1e4+200 
    #teff = np.array([3500,3750,4000,4500]) 
    #logg = np.array([2.5,3.0,2.5,3.0])
    #feh = np.array([-1.0,-0.5,-0.5,0.0])
    #tab['teff'] = teff[np.round(np.random.rand(nstars)*3).astype(int)]
    #tab['logg'] = logg[np.round(np.random.rand(nstars)*3).astype(int)]
    #tab['feh'] = feh[np.round(np.random.rand(nstars)*3).astype(int)]
     
    # Initialize the image 
    im = np.zeros((pars['nx'],pars['ny']),float)
     
    # Synthetic spectra 
    specdir = '/Users/nidever/synspec/winter2017/grid2/' 
     
    # Loop over the stars 
    for i in range(len(startab)):
        tab1 = startab[i] 
             
        # Load the synthetic spectrum for this star 
        if tab1['feh'] >= 0.0: 
            msign='+' 
        else: 
            msign='-' 
        synbase = 'spec.t%4dg%4.2fm%1s%4.2fa+0.00.fits' % (tab1['teff'],tab1['logg'],msign,np.abs(tab1['feh']))
        synfile = specdir+synbase
        if os.path.exists(synfile)==False:
            print(synfile+' NOT FOUND')
            import pdb; pdb.set_trace()
        synstr = Table.read(synfile)
        for n in synstr.colnames:   # lowercase names
            synstr[n].name = n.lower()
        print(i+1,' ',synfile)
        
        # Scale the spectrum by the magnitude/flux
        
        dum,wind = dln.closest(synstr['wave'],pars['w0'])
        synstr['flux'] *= tab1['flux']/synstr['flux'][wind]
             
        # Generate the spatial Gaussian (FWHM) kernel 
        nkernel = int(np.ceil(pars['fwhm']*2.5))
        if nkernel % 2 == 0: 
            nkernel += 1 
        xkernel = np.arange(nkernel)-nkernel/2 
        kernel = np.exp(-0.5*xkernel**2/(pars['fwhm']/2.35)**2) 
        kernel /= np.sum(kernel)  # normalize 
             
        if pars['angle'] == 0.0: 
            # Figure out the X/Y center and how far the spectrum stretches 
            # Traces are STRAIGHT (for now) 
            xr = (pars['wr']-pars['w0'])/pars['dispersion'] + tab1['x']
            xr = np.array(xr).astype(int)
            nxspec = int(xr[1]-xr[0]+1)
            wave = np.arange(nxspec)*pars['dispersion']+pars['wr'][0] 
                 
            # Bin the spectrum for the pixels 
            osamp = 4 
            xx = np.arange(nxspec*osamp).astype(float)/osamp+xr[0] 
            wavefine = np.arange(nxspec*osamp).astype(float)/osamp*pars['dispersion']+pars['wr'][0] 
            fluxfine = interp1d(synstr['wave'],synstr['flux'])(wavefine)
            # multiply times the transmission 
            transfine = interp1d(pars['trans']['wave'],pars['trans']['flux'])(wavefine)
            binflux = dln.rebin(fluxfine*transfine,nxspec) 
            #binxx = dln.rebin(xx,nxspec)
            binxx = np.arange(nxspec)+xr[0]
                 
            # Spread the spectrum by the appropriate spatial Gaussian (FWHM) 
            binflux2d = binflux.reshape(1,-1) * kernel.reshape(-1,1)
            
            # Trim overlap 
            if xr[0]<0:
                gd, = np.where(binxx >= 0) 
                binflux2d = binflux2d[:,gd] 
                xr[0] = 0 
            if xr[1]>(pars['nx']-1):
                gd, = np.where(binxx <= (pars['nx']-1)) 
                binflux2d = binflux2d[:,gd] 
                xr[1] = pars['nx']-1 
                 
            # Add it to the image
            yr = np.array([-(nkernel//2),nkernel//2])+int(tab1['y'])
            im[yr[0]:yr[1]+1,xr[0]:xr[1]+1] += binflux2d 
             
        elif pars['angle'] == 90.0: 
            # Figure out the X/Y center and how far the spectrum stretches 
            # Traces are STRAIGHT (for now) 
            yr = (pars['wr']-pars['w0'])/pars['dispersion'] + tab1['y']
            yr = np.array(yr).astype(int)
            nyspec = int(yr[1]-yr[0]+1 )
            wave = np.arange(nyspec)*pars['dispersion']+pars['wr'][0] 
                 
            # Bin the spectrum for the pixels 
            osamp = 4 
            yy = np.arange(nyspec*osamp).astype(float)/osamp+yr[0] 
            wavefine = np.arange(nyspec*osamp).astype(float)/osamp*pars['dispersion']+pars['wr'][0]
            fluxfine = interp1d(synstr['wave'],synstr['flux'])(wavefine)
            # multiply times the transmission 
            transfine = interp1d(pars['trans']['wave'],pars['trans']['flux'])(wavefine)
            binflux = dln.rebin(fluxfine*transfine,nyspec) 
            #binyy = dln.rebin(yy,nyspec)
            binyy = np.arange(nyspec)+yr[0]
                 
            # Spread the spectrum by the appropriate spatial Gaussian (FWHM)
            binflux2d = binflux.reshape(-1,1) * kernel.reshape(1,-1)                  
            
            # Trim overlap 
            if yr[0]<0:
                gd, = np.where(binyy >= 0) 
                binflux2d = binflux2d[gd,:] 
                yr[0] = 0 
            if yr[1]>(pars['ny']-1):
                gd, = np.where(binyy <= (pars['ny']-1)) 
                binflux2d = binflux2d[gd,:] 
                yr[1] = pars['ny']-1 
                 
            # Add it to the image 
            xr = np.array([-(nkernel//2),nkernel//2])+int(tab1['x'])
            im[yr[0]:yr[1]+1,xr[0]:xr[1]+1] += binflux2d
     
    return im
 

def example():
    """ Run example of synth code."""


    pars = {}
    pars['fwhm'] = 2.5 
    pars['dispersion'] = 1.0  # A per pixel 
    pars['angle'] = 0
    pars['wr'] = [10000,11000] 
    pars['w0'] = np.mean(pars['wr'])
    pars['nx'] = 2048 
    pars['ny'] = 2048 
     
    # Transmission/throughput curve
    trans = np.zeros(pars['wr'][1]-pars['wr'][0]+2,dtype=np.dtype([('wave',float),('flux',float)]))
    trans['wave'] = np.arange(pars['wr'][1]-pars['wr'][0]+2)+pars['wr'][0] 
    # -97.389234    0.0056913391   1.5786481e-06  -1.1699067e-10 
    tcoef = [-97.389234,    0.0056913391,   1.5786481e-06,  -1.1699067e-10] 
    trans['flux'] = np.polyval(np.flip(tcoef),trans['wave'])
    pars['trans'] = trans
    
    # Star catalog
    dt = [('id',int),('x',float),('y',float),('flux',float),('teff',float),('logg',float),('feh',float)]
    nstars = 500
    tab = np.zeros(nstars,dtype=np.dtype(dt))
    tab['x'] = np.random.rand(nstars)*(pars['nx']-200)+100 
    tab['x'] = np.round(tab['x']).astype(int)
    tab['y'] = np.random.rand(nstars)*(pars['ny']-200)+100 
    tab['y'] = np.round(tab['y']).astype(int)
    tab['flux'] = np.random.rand(nstars)*1e4+200 
    teff = np.array([3500,3750,4000,4500]) 
    logg = np.array([2.5,3.0,2.5,3.0])
    feh = np.array([-1.0,-0.5,-0.5,0.0])
    tab['teff'] = teff[np.round(np.random.rand(nstars)*3).astype(int)]
    tab['logg'] = logg[np.round(np.random.rand(nstars)*3).astype(int)]
    tab['feh'] = feh[np.round(np.random.rand(nstars)*3).astype(int)]

    im = image(pars,tab)
    pars['angle'] = 0    
    specim0 = spectrum(pars,tab)
    pars['angle'] = 90
    specim90 = spectrum(pars,tab)

    return im,specim0,specim90
