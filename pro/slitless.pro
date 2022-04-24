pro slitless

  ;; Synthetic slitless spectral images


  ;; -catalog of sources with
  ;;  -X/Y position in the image
  ;;  -magnitude
  ;;  -stellar parameters (Teff, logg, [Fe/H])
  ;;
  ;; -spatial fwhm
  ;; -transmission curve
  ;; -spectral resolution
  ;; -orientation anges

  ;; background, assume it is very low for the moment (this is mostly
  ;;   true for space data)
  
  ;; make two slitless spectra images at two different angles
  ;; and one direct imaging image

  ;; Parameters
  fwhm = 2.5
  dispersion = 1.0  ;; A per pixel
  angle = 90.0  ;; 0.0
  wr = [10000,11000]
  w0 = mean(wr)  
  nx = 2048
  ny = 2048

  ;; Transmission/throughput curve
  trans = replicate({wave:0.0,flux:0.0},wr[1]-wr[0]+1)
  trans.wave = findgen(wr[1]-wr[0]+1)+wr[0]
  ;; -97.389234    0.0056913391   1.5786481e-06  -1.1699067e-10
  tcoef = [-97.389234,    0.0056913391,   1.5786481e-06,  -1.1699067e-10]
  trans.flux = poly(trans.wave,tcoef)

  ;; Star catalog
  schema = {x:0.0,y:0.0,flux:0.0,teff:0.0,logg:0.0,feh:0.0}
  nstars = 500
  str = replicate(schema,nstars)
  ;str[0] = {x:600.,y:400.,flux:1000.0,teff:3750.,logg:3.0,feh:-1.0}
  ;str[1] = {x:700.,y:1500.,flux:10000.0,teff:3500.,logg:2.5,feh:-0.5}
  ;str[2] = {x:1200.,y:700.,flux:5000.,teff:4000.,logg:2.5,feh:-0.5}
  ;str[3] = {x:1300.,y:1400.,flux:2000.,teff:4500.,logg:3.0,feh:0.0}
  str.x = randomu(seed,nstars)*(nx-200)+100
  str.x = round(str.x)
  str.y = randomu(seed,nstars)*(ny-200)+100
  str.y = round(str.y)
  str.flux = randomu(0,nstars)*1e4+200
  teff = [3500,3750,4000,4500]
  logg = [2.5,3.0,2.5,3.0]
  feh = [-1.0,-0.5,-0.5,0.0]
  str.teff = teff[round(randomu(0,nstars)*3)]
  str.logg = logg[round(randomu(0,nstars)*3)]
  str.feh = feh[round(randomu(0,nstars)*3)]
  
  ;; Initialize the image
  im = fltarr(nx,ny)

  ;; Synthetic spectra
  specdir = '/Users/nidever/synspec/winter2017/grid2/'
  
  ;; Loop over the stars
  FOR i=0,nstars-1 do begin
     str1 = str[i]
     
     ;; NO dispersion
     ;;--------------
     If dispersion eq 0.0 then begin

       ;; Generate the spatial Gaussian (FWHM) kernel
       xr = [floor(str1.x-5*fwhm)>0,ceil(str1.x+5*fwhm)<(nx-1)]
       yr = [floor(str1.y-5*fwhm)>0,ceil(str1.y+5*fwhm)<(ny-1)]
       xx = (findgen(xr[1]-xr[0]+1)+xr[0])#replicate(1,yr[1]-yr[0]+1)
       yy = replicate(1,xr[1]-xr[0]+1)#(findgen(yr[1]-yr[0]+1)+yr[0])
       subim = exp(-0.5*((xx-str1.x)^2+(yy-str1.y)^2)/(fwhm/2.35)^2)
       ;; Scale by flux
       subim /= total(subim)
       subim *= str1.flux
       ;; Add to image
       im[xr[0]:xr[1],yr[0]:yr[1]] += subim
       
     ;; Dispersion
     ;;-----------
     Endif else begin
     
       ;; Load the synthetic spectrum for this star
       if str1.feh ge 0.0 then msign='+' else msign='-'
       synfile = specdir+'spec.t'+string(str1.teff,format='(I4)')+'g'+string(str1.logg,format='(F4.2)')+$
                         'm'+msign+string(abs(str1.feh),format='(F4.2)')+'a+0.00.fits'
       if file_test(synfile) eq 0 then stop,synfile+' NOT FOUND'
       synstr = mrdfits(synfile,1,/silent)
       print,i+1,' ',synfile

       ;; Scale the spectrum by the magnitude/flux
       dum = closest(w0,synstr.wave,ind=wind)
       synstr.flux *= str1.flux/synstr[wind].flux

       ;; Generate the spatial Gaussian (FWHM) kernel
       nkernel = ceil(fwhm*2.5)
       if nkernel mod 2 eq 0 then nkernel += 1
       xkernel = lindgen(nkernel)-nkernel/2
       kernel = exp(-0.5*xkernel^2/(fwhm/2.35)^2)
       kernel /= total(kernel)  ;; normalize
     
       if angle eq 0.0 then begin
         ;; Figure out the X/Y center and how far the spectrum stretches
         ;; Traces are STRAIGHT (for now)
         xr = (wr-w0)/dispersion + str1.x
         nxspec = xr[1]-xr[0]+1
         wave = lindgen(nxspec)*dispersion+wr[0]
     
         ;; Bin the spectrum for the pixels
         osamp = 4
         xx = findgen(nxspec*osamp)/osamp+xr[0]
         wavefine = findgen(nxspec*osamp)/osamp*dispersion+wr[0]
         fluxfine = spline(synstr.wave,synstr.flux,wavefine)
         ;; multiply times the transmission
         transfine = spline(trans.wave,trans.flux,wavefine)
         binflux = rebin(fluxfine*transfine,nxspec)
         binxx = rebin(xx,nxspec)
       
         ;; Spread the spectrum by the appropriate spatial Gaussian (FWHM)
         binflux2d = binflux # kernel

         ;; Trim overlap
         if xr[0] lt 0 then begin
           gd = where(binxx ge 0,ngd)
           binflux2d = binflux2d[gd,*]
           xr[0] = 0 
         endif
         if xr[1] gt nx-1 then begin
           gd = where(binxx ge 0,ngd)          
           binflux2d = binflux2d[gd,*]
           xr[1] = nx-1
         endif
       
         ;; Add it to the image
         yr = [-nkernel/2,nkernel/2]+str1.y
         im[xr[0]:xr[1],yr[0]:yr[1]] += binflux2d       
       endif

       if angle eq 90.0 then begin
         ;; Figure out the X/Y center and how far the spectrum stretches
         ;; Traces are STRAIGHT (for now)
         yr = (wr-w0)/dispersion + str1.y
         nyspec = yr[1]-yr[0]+1
         wave = lindgen(nyspec)*dispersion+wr[0]
     
         ;; Bin the spectrum for the pixels
         osamp = 4
         yy = findgen(nyspec*osamp)/osamp+yr[0]
         wavefine = findgen(nyspec*osamp)/osamp*dispersion+wr[0]
         fluxfine = spline(synstr.wave,synstr.flux,wavefine)
         ;; multiply times the transmission       
         transfine = spline(trans.wave,trans.flux,wavefine)
         binflux = rebin(fluxfine*transfine,nyspec)       
         binyy = rebin(yy,nyspec)
       
         ;; Spread the spectrum by the appropriate spatial Gaussian (FWHM)
         binflux2d = binflux # kernel

         ;; Trim overlap
         if yr[0] lt 0 then begin
           gd = where(binyy ge 0,ngd)
           binflux2d = binflux2d[gd,*]
           yr[0] = 0 
         endif
         if yr[1] gt ny-1 then begin
           gd = where(binyy ge 0,ngd)          
           binflux2d = binflux2d[gd,*]
           yr[1] = ny-1
         endif
       
         ;; Add it to the image
         xr = [-nkernel/2,nkernel/2]+str1.x
         im[xr[0]:xr[1],yr[0]:yr[1]] += transpose(binflux2d)
       endif

    Endelse  ; dispersion
       
  ENDFOR  ;; star loop

  
  stop


end
