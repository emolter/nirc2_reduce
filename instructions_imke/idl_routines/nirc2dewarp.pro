;+
; NAME:
;	NIRC2DEWARP
;
; PURPOSE:
;	Remove distortion from a NIRC2 image.
;
; EXPLANATION:
;	All three NIRC2 cameras have known distortion.  This function returns
;	a rectified image, with each pixel in im_out bilinearly interpolated
;	at the correct location in im_in.  This function is NOT INTENDED FOR
;	PRECISE PHOTOMETRIC ANALYSIS, since the variations in the effective
;	area of pixels has not been addressed.
;
;	A discussion of the NIRC2 optical quality is provided in the pre-ship
;	testing document, available at:
;	http://www.keck.hawaii.edu/realpublic/inst/nirc2/preship_testing.pdf
;       and updated investigation of the distortion can be found at:
;       http:// URL
;
; CALLING SEQUENCE:
;	im_out = NIRC2DEWARP(im_in [,hd=hd, camera=camera, version=version])
;
; INPUTS:
;	IM_IN - Input image
;
; OPTIONAL INPUTS:
;	None
;
; OPTIONAL KEYWORD PARAMETERS:
;	HD - Input image FITS header.  Used only to determine camera.
;
;	CAMERA - Used only if HD is not specified, should be set to
;	  'narrow', 'medium', or 'wide'. Note medium camera does not
;         have a current solution, pre-ship will be used.
;
;       VERSION - Version of the solution to be used. Default is 'new,
;              but 'preship' is also available.
;
; OUTPUTS:
;	IM_OUT - Distortion-corrected image
;
; SIDE EFFECTS:
;	Giddyness
;
; RESTRICTIONS:
;	Sub-images (less than 1024x1024 pixels) are assumed to be centered
;	on the NIRC2 array.  De-centered subimages will not be correctly
;	rectified!
;
; PROCEDURE:
;	Note that there is some mathematical sloppyness here.  Instead of
;	inverting the transformation matrix, we have simply applied the
;	negative of the forward transformation to find the locations at
;	which to interpolate the input image.
;
; EXAMPLE:
;	IDL> im_in = READFITS('n0001.fits',hd)
;	IDL> im_out = NIRC2DEWARP(im_in,hd=hd)
;
; PROCEDURES CALLED:
;	BILINEAR()
;
; MODIFICATION HISTORY:
;	Writen by A. Bouchez, Keck Observatory, July 2004.
;       Revised by P. B. Cameron, Caltech, March 2007.
;
;-
FUNCTION	NIRC2DEWARP, im_in, hd=hd, camera=camera, version=version

;;; 1. Cubic polynomial coefficients from pre-ship testing document

  sz = SIZE(im_in)
  x0 = sz[1]/2.		; X origin of polynomial
  y0 = sz[2]/2.		; Y origin of polynomial

                                ;The offset terms are irrelevant for
                                ;correcting distortion (which is
                                ;simply relative). We set them here to
                                ;zero for the new solution, but
                                ;include them for the preship solution
                                ;since they were used previously.

    if ( KEYWORD_SET(version) ) then begin
      vers=version
  endif else begin
      vers = 'new'
  endelse

  if( vers eq 'preship') then begin
      a = DBLARR(3,15)
      a[*,0] = [-2.7200d-01, -3.12530d-01,  4.0600d-01]
      a[*,1] = [ 1.0009d+00,  9.99382d-01,  1.0008d+00]
      a[*,2] = [ 5.0800d-03,  1.23970d-03, -3.2400d-03]
      a[*,3] = [-5.5200d-06, -2.80460d-06, -1.7500d-06]
      a[*,4] = [ 7.7000d-07,  1.20020d-06, -5.6000d-07]
      a[*,5] = [ 5.3700d-06,  1.08800d-05,  1.0000d-05]
      a[*,6] = [-8.6000d-09, -9.61700d-09, -6.8000d-09]
      a[*,7] = [-8.0000d-10, -1.51400d-09,  8.0000d-10]
      a[*,8] = [-6.1000d-09, -5.05300d-09, -3.6000d-09]
      a[*,9] = [-4.4000d-09,  3.90000d-10,  9.0000d-10]
      
      b = DBLARR(3,15)
      b[*,0] = [ 1.6800d-01,  3.87000d-01,  3.4000d-02]
      b[*,1] = [ 1.6000d-04, -1.40000d-04,  2.1000d-04]
      b[*,2] = [ 1.0008d+00,  9.95730d-01,  9.9477d-01]
      b[*,3] = [ 2.1000d-07,  7.60000d-07, -4.4000d-07]
      b[*,4] = [-9.9300d-06, -8.27000d-06, -9.2000d-07]
      b[*,5] = [ 1.7800d-06,  1.86000d-06, -1.5700d-06]
      b[*,6] = [-1.6000d-09,  3.00000d-10,  2.0000d-10]
      b[*,7] = [-2.8000d-09, -1.60000d-09,  2.6000d-09]
      b[*,8] = [-4.0000d-10,  8.00000d-10,  5.0000d-10]
      b[*,9] = [-1.2000d-08, -6.40000d-09,  5.9000d-09]
  endif else begin
      a = DBLARR(3,15)
      a[*,0] = [         0d0, -3.12530d-01,         0d0]
      a[*,1] = [ 1.00116d+00,  9.99382d-01, 1.00258d+00]
      a[*,2] = [ 1.96010d-03,  1.23970d-03,-1.07553d-03]
      a[*,3] = [-3.14182d-06, -2.80460d-06,-8.44603d-07]
      a[*,4] = [-3.08788d-06,  1.20020d-06,-8.66153d-07]
      a[*,5] = [ 5.59549d-06,  1.08800d-05, 1.03161d-05]
      a[*,6] = [-6.04724d-09, -9.61700d-09,-2.80832d-09]
      a[*,7] = [ 6.53706d-10, -1.51400d-09, 2.11712d-12]
      a[*,8] = [-3.99943d-09, -5.05300d-09,-8.90944d-10]
      a[*,9] = [-1.63796d-09,  3.90000d-10,-3.14757d-11]
      a[*,10]= [-1.12204d-11,          0d0,         0d0]
      a[*,11]= [ 1.35291d-11,          0d0,         0d0]
      a[*,12]= [ 5.73354d-12,          0d0,         0d0]
      a[*,13]= [ 3.37186d-12,          0d0,         0d0]
      a[*,14]= [-8.50332d-13,          0d0,         0d0]
      
      b = DBLARR(3,15)
      b[*,0] = [         0d0,  3.87000d-01,         0d0]
      b[*,1] = [ 1.75403d-03, -1.40000d-04,-1.39196d-03]
      b[*,2] = [ 1.00129d+00,  9.95730d-01, 9.96806d-01]
      b[*,3] = [-1.91215d-06,  7.60000d-07,-2.81300d-07]
      b[*,4] = [-1.16438d-05, -8.27000d-06,-5.77556d-07]
      b[*,5] = [-2.51404d-06,  1.86000d-06,-1.44116d-06]
      b[*,6] = [ 6.03978d-11,  3.00000d-10, 5.30766d-10]
      b[*,7] = [-3.24360d-09, -1.60000d-09, 4.39543d-09]
      b[*,8] = [-3.55861d-09,  8.00000d-10,-3.07567d-10]
      b[*,9] = [-8.66718d-09, -6.40000d-09, 7.84907d-09]
      b[*,10]= [ 5.17175d-12,          0d0,         0d0]
      b[*,11]= [ 2.66960d-12,          0d0,         0d0]
      b[*,12]= [ 5.29880d-12,          0d0,         0d0]
      b[*,13]= [ 7.83758d-12,          0d0,         0d0]
      b[*,14]= [ 4.80703d-12,          0d0,         0d0]
  endelse
 

;;; 2. Camera definition

  if KEYWORD_SET(hd) then cam = STRTRIM(STRUPCASE(SXPAR(hd,'CAMNAME')),2)
  if KEYWORD_SET(camera) then cam = STRTRIM(STRUPCASE(camera),2)
  if not KEYWORD_SET(cam) then begin
    cam_input = ''
    READ,"Define camera ('narrow','medium',or 'wide'): ",cam_input
    cam = STRTRIM(STRUPCASE(cam_input),2)
  endif

  case cam of
    'NARROW':i=0
    'MEDIUM':i=1
    'WIDE':i=2
  endcase

;;; 3. Set up x and y coordinate arrays
  xobs = REBIN(FINDGEN(sz[1],1),sz[1],sz[2])
  x1 = xobs - x0
  yobs = REBIN(FINDGEN(1,sz[2]),sz[1],sz[2])
  y1 = yobs - y0

;;; 4. Perform forward transformation (x -> x')
  xp = 512.+POLYSOL(x1,y1,a[i,*])
  yp = 512.+POLYSOL(x1,y1,b[i,*])

;;; 5. Bilinearly interpolate output image at negative offset locations
  ix = 2*xobs - xp
  jy = 2*yobs - yp
  im_out = BILINEAR(im_in,ix,jy)

RETURN, im_out
END

;#########################
FUNCTION POLYSOL, x, y, p

n = p[0] + p[1]*x + p[2]*y + p[3]*x^2. + p[4]*x*y + p[5]*y^2. + $
     p[6]*x^3. + p[7]*x^2*y + p[8]*x*y^2 + p[9]*y^3 + $
     p[10]*x^4. + p[11]*x^3.*y + p[12]*x^2.*y^2. + p[13]*x*y^3. + p[14]*y^4.

RETURN,n

END
