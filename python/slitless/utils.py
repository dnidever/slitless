#!/usr/bin/env python

import os
import time
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table
from dlnpyutils import utils as dln

class Projection():

    def __init__(self,pars):
        for n in ['wr','w0','nx','ny','angle','dispersion']:
            if hasattr(pars,n)==False:
                raise ValueError(n+' not found in input parameters')
            setattr(self,pars[n])

        # Trace is linear

        # Length of the trace (in pixels)
        tracelenpix = (wr[1]-wr[0])/dispersion
        self.tracelenpix = tracelenpix
        
        # Vertical
        if angle % 180 == 0:
            self.vertical = True
            self.pslp = np.inf
            self.poffset = 0
            self.x2wslp = np.inf
            self.w2xslp = np.nan
            if angle % 360 == 0:
                self.w2yslp = dispersion
            else:
                self.w2yslp = -dispersion                
                
        # At an angle
        else:            
            self.vertical = False

            # Trace position
            # y = m*x + b
            pslp = np.tan(np.deg2rad(angle))
            poffset = y0-pslp*x0
            self.pslp = pslp
            self.poffset = poffset
            
            # Wavelength along the trace
            # w = m2*x + b2
            x2wslp = np.sqrt(1+pslp**2)*dispersion
            #x2woffset = w0-x2wslp*x0
            self.x2wslp = x2wslp
            
            # Solve for beginning/ending X values
            #x1 = (wr[0]-woffset)/wslp
            #x2 = (wr[1]-woffset)/wslp
            
            # Wave to X, reverse of above
            # x = m*w + b
            w2xslp = 1/x2wslp
            #w2xoffset = -x2woffset/x2wslp
            self.w2xslp = w2xslp
            
            # Wave to Y
            # y = m*w + b
            # use w->x and then x->y
            w2yslp = w2xslp * pslp
            #w2yoffset = w2xoffset * pslp + poffset
            self.w2yslp = w2yslp
            
            #wave = np.linspace(wr[0],wr[1],dispersion)
            #x = w2xslp * wave + w2xoffset
            #y = w2yslp * wave + w2yoffset

    def __call__(self,coo):
        x0 = coo[0]
        y0 = coo[1]
        wave = np.linspace(self.wr[0],self.wr[1],self.dispersion)
        
        # Vertical
        if self.angle % 180 == 0:
            # Trace position
            x = np.zeros(len(wave),float)+x0
            # Wave to Y
            # Y = w2yslp*W + w2yoffset
            w2yoffset = y0-1/self.dispersion*self.w0
            y = wave/self.dispersion+w2yoffset
                
        # At an angle
        else:            
            # Trace position
            # Y = pslp*X + poffset
            poffset = y0-self.pslp*x0
            
            # Wavelength along the trace
            # W = x2wslp*X + x2woffset
            x2woffset = self.w0-self.x2wslp*x0
            
            # Solve for beginning/ending X values
            x1 = (wr[0]-x2woffset)/self.x2wslp
            x2 = (wr[1]-x2woffset)/self.x2wslp
            
            # Wave to X, reverse of above
            # X = w2xslp*W + w2xoffset
            w2xoffset = -x2woffset/self.x2wslp
            
            # Wave to Y
            # Y = w2yslp*W + w2yoffset
            # use w->x and then x->y
            w2yoffset = w2xoffset * self.pslp + poffset
            
            x = self.w2xslp * wave + w2xoffset
            y = self.w2yslp * wave + w2yoffset
        
        return Trace(w,x,y)
        
    def wave2xy(self,wave):
        pass

class Trace():

    def __init__(self,wave,x,y):
        self.npix = len(wave)
        self.pathlength = np.sum((x[0:-1]-x[1:])**2+(y[0:-1]-y[1:])**2)
        self.wave = wave
        self.x = x
        self.y = y

    def wave2xy(self,w):
        x = interp1d(self.wave,self.x)(w)
        y = interp1d(self.wave,self.y)(w)        
        return x,y

    
def trace(pars,coo):
    """ Make the trace for a star."""

    x0 = coo[0]
    y0 = coo[1]
    wr = pars['wr']
    w0 = pars['w0']
    nx = pars['nx']
    ny = pars['ny']
    angle = pars['angle']
    dispersion = pars['dispersion']
    
    dt = [('x',float),('y',float),('wave',float)]
    info = np.zeros(1000,dtype=np.dtype(dt))

    # Trace is linear

    # Length of the trace (in pixels)
    tracelenpix = (wr[1]-wr[0])/dispersion
    
    # Vertical
    if angle % 180 == 0:
        pass

    # At an angle
    else:
        # Trace position
        # y = m*x + b
        pslp = np.tan(np.deg2rad(angle))
        poffset = y0-pslp*x0

        # Wavelength along the trace
        # w = m2*x + b2
        x2wslp = np.sqrt(1+pslp**2)*dispersion
        x2woffset = w0-x2wslp*x0

        # Solve for beginning/ending X values
        x1 = (wr[0]-woffset)/wslp
        x2 = (wr[1]-woffset)/wslp

        # Wave to X, reverse of above
        # x = m*w + b
        w2xslp = 1/x2wslp
        w2xoffset = -x2woffset/x2wslp

        # Wave to Y
        # y = m*w + b
        # use w->x and then x->y
        w2yslp = w2xslp * pslp
        w2yoffset = w2xoffset * pslp + poffset
        
        wave = np.linspace(wr[0],wr[1],dispersion)
        x = w2xslp * wave + w2xoffset
        y = w2yslp * wave + w2yoffset

            
    return im

