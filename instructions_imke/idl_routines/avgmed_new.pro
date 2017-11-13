pro avgmed_new,out,images,filename=filename,log=log, $
		noweight=noweight,nocoadd=nocoadd,average=average, $
		evenmed=evenmed,maxrej=maxrej,minrej=minrej, $
		silent=silent

; Median averages FITS images whose names are given in a file.
; All frames are weighted by their medians, unless /noweight is set.
; If /log, then writes weights to a file as well and the options chosen.
;
; 5/20/94 MCL (UCB)
;
; All frames read from file are normalized by their # of coadds 
; unless /nocoadd is set
;
; 6/27/94
;
; Added ability to read images from IDL array (in a very inelegant way)
;   if want to work with IDL array, just call avgmed like:
;	IDL> avgmed, out, array
;   to read stuff from a file, just define the 'filename' parameter:
;	IDL> avgmed, out, filename='file.list'
; 7/01/94
;
; significant modifications:
;
;   1. 	modified median averaging so images are scaled to the average
;    	of all the medians, instead of the median of the first image
;   	-> this bombs if an image has a median = 0
;
;   2. 	when median filtering if have <= 6 images, will take the mean
;	of the middle 2 images when there is an even (2,4,6) number of images
;	(only if /evenmed set since this can be quite S-L-O-W); this probably
;	make little difference in practice
;
;   3.  added 'minrej' and 'maxrej': will toss these number of values
;	out before median averaging (no effect for arithmatic avg)
;	note: this is done AFTER any weighting.
;
;   4.  major editing of program to make it nicer
;
; 12/13/94
;
; added minrej/maxrej option for arithmatic averaging as well
;
; 6/2/95 MCL
;
; modified maxrej/minrej to toss more than one pixel
;
; I think this all works ... hard to test the entire parameter space
;   esp. have not tested the new version of the file reader 
;   (though could just replace it with 'loadimf')
;
; 'avgmed_new' = testing if will run faster if don't make
; 	a copy of the original images, instead doing scaling
;	when do median filtering - runs about 10% faster for
;	a stack of 30 images with simple median averaging
;
; 7/2/95 MCL
;
; corrected error in calculating the proper image weights
; (was incorrectly scaling the weights to 1.0 by multiplying by the
;   average value instead of dividing by it!)
; 12/17/96 MCL
; 
; the above bug fix was wrong - no idea why I thought it was right!
; 01/29/97


if (n_params() eq 0) then begin
   print,'avgmed_new,out,images,[filename=],[log=],'
   print,'       [noweight],[nocoadd],[average],'
   print,'       [evenmed],[maxrej=],[minrej=],[silent]'
   return
endif


; beginning messages
if keyword_set(silent) eq 0 then begin

  print,'--> did you remember to subtract the bias? <--'

  if (keyword_set(noweight) eq 1) then print,'* no weighting for averaging *' $
     else print,'* median weighting for averaging *'

  if (keyword_set(filename) eq 1) then begin
     if (keyword_set(nocoadd) eq 1) then $
	print,'* not dividing each frame by number of coadds *' $
     else print,'* dividing each frame by number of coadds *'
  endif

endif


;
; --------------------------------------------------
; (1) figure out how many images:
; --------------------------------------------------
;     (a) either from an input file, not counting blank lines 
;
if keyword_set(filename) eq 1 then begin

   print,format='($,"* reading images from file ",A,": ")',filename 
   openr,unit0,filename,/get_lun
   n = 0
   file = ' ' 
   while eof(unit0) ne 1 do begin
 	readf,unit0,file
	if (strlen(file) ne 0) then n = n+1
   endwhile
   print,strc(n),' files *'
   free_lun,unit0

;
; (b) or from IDL array
;
endif else begin 
  sz = size(images)
  xs = sz(1)
  ys = sz(2)
  n = sz(3)
  imgs = fltarr(xs,ys,n)
  if keyword_set(silent) eq 0 then $
  	print,'* reading ',strc(n),' images from IDL array *'
endelse


;
; set up output array and make sure things are ok
;
out = fltarr(xs,ys)


;
; once set up maxrej and minrej to handle more than one pixel,
;    these statements will be useful:
;
;if keyword_set(maxrej) and keyword_set(minrej) and $
;  (maxrej+minrej) gt n then message,'(minrej+maxrej) too large!'
;if keyword_set(minrej) eq 0 then minrej = 0 $
;  else if minrej gt n then message,'minrej too large!' 
;if keyword_set(maxrej) eq 0 then maxrej = 0 $
; else if maxrej gt n then message,'maxrej too large!'
	

;
; --------------------------------------------------
; (2) read in images 
; --------------------------------------------------;

; (2a) from a file
;
if keyword_set(filename) eq 1 then begin

   openr,unit0,filename,/get_lun

;  get the images
   for i =0,(n-1) do begin

	readf,unit0,file

;	first image defines array sizes
	if (i eq 0) then begin

	  out = readfits(file,head,/silent)
	  sz = size(out)
	  xs = sz(1)
	  ys = sz(2)
	  imgs = fltarr(xs,ys,n)
	  imgs(*,*,0) = out
	  print,'image size = ',strc(sz(1)),' x ',strc(sz(2))

	endif else begin

	  imgs(*,*,i) = readfits(file,head,/silent)

	endelse

;	normalize by # of coadds if desired
	if (keyword_set(nocoadd) eq 0) then begin 
	  extract_head,head,'M56COAD0',coadd
	  imgs(*,*,i) = imgs(*,*,i) / float(coadd)
	endif

   endfor
   free_lun,unit0


;
; (2b) or images already in a 3d array, just get the medians
;
endif else begin

;;;  imgs = images			; don't change the originals!

endelse


;
; --------------------------------------------------
; (3) median weighting so medians for all images
;     are equal to the average of all the image medians
; --------------------------------------------------
;
wtlist = fltarr(n)+1.0
if keyword_set(noweight) eq 0 then begin

; 	get medians
	medlist = fltarr(n,/nozero)		
	for k=0,n-1 do begin
		medlist(k) = median(images(*,*,k))
	endfor

;	scale the images
	mm = total(medlist) / n
	for k = 0,n-1 do begin
		wtlist(k) = mm / medlist(k)   
;;;		imgs(*,*,k) = temporary(imgs(*,*,k)) * wtlist(k)
	endfor

;	fancy two column printing of weights
	if keyword_set(silent) eq 0 then begin
	  for k = 0,n/2-1 do begin
		print,'image ',strc(k),' weight = ',strc(wtlist(k)),$
			'      image ',strc(k+n/2),' weight = ', $
			strc(wtlist(k+n/2))
	  endfor
	  if (n-1)/2 eq n/2 then $
		print,'image ',strc(n-1),' weight = ',strc(wtlist(n-1))
	endif

endif


;
; --------------------------------------------------
; (4) median averaging
; --------------------------------------------------
;
if keyword_set(average) eq 0 then begin

  print,'* median averaging frames *' 

;
; (4a) evenmed: if number of images is even,
;	instead of using IDL's median (the higher value), take average
;	of the two values straddling the center
;
  if (n/2 eq n/2.0) and keyword_set(evenmed) then begin
    print,'* even # of frames: using modified median *'
    for j = 0,(ys-1) do $
	for i = 0,(xs-1) do begin
	  cut = images(i,j,*)*wtlist
	  cut = cut(sort(cut))
	  out(i,j) = (cut((n-1)/2) + cut(n/2)) / 2.0
	  if i eq 40 and j eq 40 then $
		print,'compare ',out(i,j),' with ',median(images(i,j,*))
	endfor

;
; (4b) minrej & maxrej: toss the highest and/or lowest pixels before
;	median averaging
;
  endif else if keyword_set(maxrej) or keyword_set(minrej) then begin
    aa = 0
    zz = n-1
    if keyword_set(maxrej) then begin
	print,'* maxrej: tossing highest pixel *'
	zz = n-2
    endif
    if keyword_set(minrej) then begin
	print,'* minrej: tossing lowest pixel *'
	aa = 1
    endif
    for j = 0,(ys-1) do begin
	for i = 0,(xs-1) do begin
	  cut = images(i,j,*)*wtlist
	  cut = cut(sort(cut))
;	  cut = images(i,j,sort(images(i,j,*)))
	  out(i,j) = median(cut(aa:zz))
	endfor
    endfor

;
; (4c) just median average each pixel
;
  endif else begin
    for j = 0,(ys-1) do $
	for i = 0,(xs-1) do begin
	  out(i,j) = median(images(i,j,*)*wtlist)
	endfor
  endelse


;
; --------------------------------------------------
; (5) arithmatic averaging
; --------------------------------------------------
;
endif else begin

  print,'* arithmatic averaging frames *'

;
; (5a) with minrej and/or maxrej before averaging
;
  if keyword_set(maxrej) or keyword_set(minrej) then begin
    aa = 0
    zz = n-1
    nn = n
    if keyword_set(maxrej) then begin
	print,'* maxrej: tossing highest ',strc(maxrej),' pixels *'
	zz = zz - maxrej
	nn = n - maxrej
    endif
    if keyword_set(minrej) then begin
	print,'* minrej: tossing lowest ',strc(minrej),' pixels *'
	aa = minrej
	nn = n - minrej
    endif
    for j = 0,(ys-1) do begin
	for i = 0,(xs-1) do begin
	  cut = images(i,j,*)*wtlist
	  cut = cut(sort(cut))
	  out(i,j) = total(cut(aa:zz))/nn
	endfor
    endfor

;
; (5b) plain arithmatic averaging
;
  endif else begin
    for j = 0,(ys-1) do $
	for i = 0,(xs-1) do $
	  out(i,j) = total(images(i,j,*)*wtlist)/n
  endelse

endelse


;
; --------------------------------------------------
; (6) write to log file
; --------------------------------------------------
;
if keyword_set(log) then begin

   print,'writing to log file: ',log
   openw,unitlog,log,/get_lun
   printf,'** log file from avgmed.pro **'

   if keyword_set(filename) eq 1 then begin
	openr,unit0,filename,/get_lun
   	file = ' ' 
   	while eof(unit0) ne 1 do begin
 		readf,unit0,file
		if (strlen(file) ne 0) then n = n+1
	endwhile
	printf,unitlog,file,'	',strc(wtlist(i))  
   	free_lun,unit0
   endif else for i=0,(n-1) do printf,unitlog,'image ',strc(i),$
	'	',strc(wtlist(i))  

   if keyword_set(nocoadd) then $
	printf,unitlog,'* not normalized to counts/coadd *' $
   else printf,unitlog,'* normalized to counts/coadd *'

   if keyword_set(noweight) then $
	printf,unitlog,'* no weighting for frames *' $
   else printf,unitlog,'* median weighting for frames *'

   if keyword_set(average) eq 0 then $
	printf,unitlog,'* median averaged frames *' $
   else printf,unitlog,'* arithmatic averaged frames *'

   close,unitlog
endif

end

