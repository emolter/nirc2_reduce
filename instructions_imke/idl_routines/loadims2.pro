pro Loadims, istart, iend, out, head0, imlist, ncoaddlist, $
             step=step, numlist=numlist, $
             fname=fname, fdir=fdir, ext=ext, $
             coadd=coadd, $
             gemini=gemini,lick=lick, ctio = ctio, ircal = ircal, scam = scam, $
             silent = silent

;+
; loads FITS images with sequential or irregularly spaced numbers and
; divides each by number of coadds unless /coadd set
;
; loads the header of the first image into 'head0'
; if header is loaded and images are divided by # of coadds, then
;  the coadd header keyword is changed to 1 and the old # of coadds
;  is saved in a new field with '_OLD' appended to it
;
; default loading convention is set to working with NIRC images
;  (eg, 's12345.fits')
;
; INPUTS
;   istart   starting file number
;   iend     ending file number
;
; OUTPUTS
;   out      output images
;   head0    header of the first image 
;            (with coadd keyword possibly changed - see above)
;   imlist   string array with list of images loaded
;   ncoaddlist  integer array with # of coadds for each image
;
; KEYWORD PARAMETERS
;   fname    file prefix name
;   fdir     file directory
;   ext      file suffix name
;   /coadd   preserve the number of coadds, don't divide
;   step     increment for file numbers
;   numlist  list of file numbers (can be irregularly spaced)
;   /lick    use name convention for Lick LIRC2 (e.g. 'd1234.ccd') 
;   /gemini  use name convention for UCLA Gemini camera ('01noa123.fts') 
;   /ctio    use name convention for CTIO's CIRIM ('ctio1234.fits')
;              *and* set /coadd since CIRIM coaverages
;   /ircal   use name conventions for Lick IRCAL (e.g., 'ircal1234.fits')
;
; HISTORY: Written by M. Liu (UCB) 07/03/94
; 06/28/95 (MCL): added /lick
; 12/01/95 (MCL): added 'step' keyword
; 04/30/96 (MCL): prints object name as files are loaded
; 10/04/96 (MCL): added 'numlist' keyword
; 02/25/97 (MCL): added header loading ('head0' variable)
;                   changes header if divide by # of coadds
; 05/13/97 (MCL): if load an image w/only 1 coadd, won't alter header 
;                   (good for already reduced images)
; 05/29/97 (MCL): corrected bug if no coadd header keyword present
; 12/11/97 (MCL): checks if files exist
; 04/02/98 (MCL): added 'imlist' output
; 01/24/99 (MCL): added 'ncoaddlist' output array (used by BIAS.PRO)
; 07/07/99 (MCL): added /ircal, doesn't divide by # of reads as it should
; 09/26/99 (MCL): added /scam for NIRSPEC SCAM
;-

;on_error, 2

; default directory and suffix for file names - CHANGE TO WHAT YOU WANT
if keyword_set(scam) then begin
   if keyword_set(fname) eq 0 then fname = '23sei'
    if keyword_set(ext) eq 0 then ext = '.fits'
    cchead = 'COADDS'
endif else begin
    if keyword_set(fname) eq 0 then fname = 's'
    if keyword_set(ext) eq 0 then ext = '.fits'
    cchead = 'M56COAD0'
endelse
if not(keyword_set(step)) then step = 1
if keyword_set(fdir) eq 0 then fdir=''


; user info
if n_params() lt 3 then begin
    print, 'loadims, istart, iend, out, (head0), (imlist), $
    print, '         [step='+strc(step)+'], [numlist=], $
    print, '         [fname='+strc(fname)+'], [fdir='+strc(fdir)+'], [ext='+strc(ext)+'], $
    print, '         [coadd],  $
    print, '         [gemini], [lick], [ctio], [scam] $
    print, '         [silent]'
    return
endif

ncoadd=indgen(1)
ncoadd=1

; check if dividing by coadds
if keyword_set(silent) eq 0 then $
  if keyword_set(coadd) then $
  message, 'not dividing each frame by number of coadds', /info $
else begin
    message, '* dividing each frame by number of coadds *', /info
    message, '  (header variable will be accordingly altered)', /info
endelse


; if desired, use a specific list of image numbers
; (good for non-sequential files)
if keyword_set(numlist) then begin
    sz =  size(numlist)
    if sz(0) ne 1 then begin
        message, 'numlist must be a 1-d array!', /info
        retall
    endif else $
      n = n_elements(numlist) 
endif else n = (iend-istart)/step + 1 
if (n le 0) then begin
    message, 'file order is backwards!'
    return
endif


; begin loading loop

for i = 0, n-1 do begin

    hh = ''
    if keyword_set(numlist) then nn = numlist(i) else nn = i*step+istart
    
;   construct the file name
    if keyword_set(scam) then begin
        file = strc(fname+'0')
        index = strc(string(nn))
        strput, file, index, strlen(file)-strlen(index)
    endif else begin
        file = strc(fname+'0')
        index = strc(string(nn))
        strput, file, index, strlen(file)-strlen(index)
    endelse

    file = fdir + file + ext
    if keyword_set(silent) eq 0 then $
      print, format = '($,"  opening ",A)', file
   
;   get first image, initialize output arrays
    if (i eq 0) then begin
        if ((findfile(file))(0) eq '') then begin
            print
            message, 'file does not exist = "'+file+'"', /info
            retall
        endif

        out = readfits(file, hh, /silent)
        sz = size(out)
        out = fltarr(sz(1), sz(2), n)
        ncoaddlist = intarr(n)
        imlist = strarr(n)
            head0 = hh
            if not(keyword_set(coadd)) then extract_head,hh,'COADDS',ncoadd
	    print,ncoadd
        
    ; load files 
    if ((findfile(file))(0) eq '') then begin
        print
        message, 'file does not exist = "'+file+'"', /info
        retall
    endif
	endif

    out(*, *, i) = readfits(file, hh, /silent)
    imlist(i) = file
            if not(keyword_set(coadd)) then extract_head,hh,'COADDS',ncoadd
	    print,ncoadd

    ncoaddlist(i) = ncoadd
   
; divide by number of coadds
          out(*, *, i) = (out(*, *, i)) / ncoadd
        if keyword_set(silent) eq 0 then $
          if (ncoadd gt 1) then print,format = '($," dividing by")' $
        else print, form = '($," not dividing by")'
    
    
endfor


end
