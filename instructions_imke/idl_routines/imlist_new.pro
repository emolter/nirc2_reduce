pro imlist_new,istart,iend,$
               infile = infile,outfile = outfile, $           
               fname = fname,fdir = fdir,ext = ext,$
               gemini = gemini,lick = lick, ctio  = ctio, ircal = ircal, $
               scam = scam, lws=lws,nirspec=nirspec, kcam=kcam, $
               silent=silent

;hr's notes:
;ext=filenameextentions ...e.g. '.fits' or '.fits.gz'
; a simple example:
;imlist_new,60,63,outfile=testlog,fdir='/data/hroe/lws/23sep99/',fname='lws',ext='.fits.gz',/lws

;or
; imlist_new,1,479,outfile='/data/hroe/nirspec/23sep99/23sep99.speclog',fdir='/data/hroe/nirspec/23sep99/spec/',fname='23ses',ext='.fits.gz',/nirspec

;can also do a list of files (useful when there are breaks in numbering
;sequence:
;$ls /home/hroe/obsdata/lws/23nov99UT/*.fits.gz > temp.in
;imlist_new,infile='temp.in',outfile='log.991123UT',/lws


;+
; Program to read and print out header info from images.  User
; specifies either a range of file numbers or a text file which
; contains the list of images (useful for non-consequtively numbered
; images).
;
; This is a re-write of the original IMLIST.PRO with
; considerably more flexible.  To change which keywords are printed,
; simply adjust the list in the 'kwlist' variable.  Note that it is 
; unable to make any special accomodations for particular keywords,
; e.g., if you want to divide the int. time in the header by 10
; (like we have to do for the LIRC2 data).
;
; Note that any floats/doubles are rounded to two decimal places
;
; Default uses the filename conventions for NIRC.
;
; USES
;    readfits, strdup, strc, strn (Goddard)
;
; CAVEATS
; Routine uses the values in the first file in the list to determine
; the size of each column, if later images have longer values this
; will lead to info being chopped off.  To fix this, set a minimum
; size for the corresponding element of the 'len' array, as I have
; done for OBJECT and some of the other keywords.
;
; NOTES
; - needs to be tuned slightly for diff telescope's header words
; - by default the 'OBJECT' keyword is allow to keep its max length
; - to print use:
;       enscript -r -fCourier8   (for Lick)
;       enscript -r -fCourier8   (for NIRC and SCAM)
;       enscript -r -fCourier9   (for Gemini)
;       enscript -r -fCourier9   (for CTIO)
;
; HISTORY
; Written by M. Liu (UCB): 02/13/95
; 07/25/95 (MCL): added /lick keyword option
; 09/16/95 (MCL): prints the new LIRC2 filter keyword info
; 09/20/95 (MCL): now handles a list of files given in 'filename'
;                 (i.e. absorbed 'imlistf.pro')
; 04/02/96 (MCL): for /lick, automatically divides the int. time by
;                 the appropriate power of 10 to give actual int. time
; 11/05/96 (MCL): new flexible scheme for header keywords (see above)
;                 removed the divide by 10 for Lick data
; 03/29/97 (MCL): established min length for some keywords (OBJECT, etc.)
; 10/24/97 (MCL): added /ctio
; 11/28/97 (MCL): fields of CTIO 4-m headers can be 70+ chars so use
;                   SXPAR_CTIO.PRO for these data (kludge)
; 07/11/98 (MCL): now checks if outfile already exists before opening
; 10/14/98 (MCL): small bug fix - if didn't find an image, was using
;                   previous image's header info. now it will crash.
; 07/07/99 (MCL): added /ircal option
; 09/26/99 (MCL): added /scam (for NIRSPEC SCAM) -> check Ra/Dec is ok!!
; 09/27/99 (HR): added /lws option (Keck LWS)
; 10/06/99 (HR): added /nirspec option (NIRSPEC spectra)
; 10/30/99 (HR): added /kcam option (AO/KCAM imaging)
;-


; Establish parameters
;   fdir    directory name (appended to all files in input file)
;   fname   base file name
;   ext     file extension name
;   kwlist  list of keywords to get from header
;   kwtitle list of titles for the keywords when printing
;
if not(keyword_set(fdir)) then fdir = ''
if keyword_set(lick) then begin
    kwlist = ['DATE-OBS','OBJECT','RA','DEC','TIME','HA','COADDS',$
              'EXPOSURE','FRONT','BACK','LENS']
    kwtitle = ['Date(UT)','Object','RA','Dec','Time(UT)','HA','Coadds',$
              'T(int)','Front','Back','Lens']
;    kwlist = ['DATE-OBS','OBJECT','RA','DEC','TIME','HA','EXPOSURE',$
;              'COADDS','FRONT','BACK','LENS']
;    kwtitle = ['Date(UT)','Object','RA','Dec','Time(UT)','HA','T(int)',$
;              'Coadds','Front','Back','Lens']
    if keyword_set(fname) eq 0 then fname = 'd'
    if keyword_set(ext) eq 0 then ext = '.ccd'
endif else if keyword_set(ircal) then begin
    kwlist = ['DATE-OBS','OBJECT','COADDTIM','NCOADDS','NREADS', 'AIRMASS',  $
              'RA','DEC','TIME','FILTER1','FILTER2','APERTURE']
    kwtitle = ['Date(UT)','Object','Tint(ms)','Coadds','Reads','AM', $
               'RA','DEC','UT Time', 'Filter1','Filter2','Aperture']
    if keyword_set(fname) eq 0 then fname = 'ircal'
    if keyword_set(ext) eq 0 then ext = '.fits'
endif else if keyword_set(gemini) then begin
    kwlist = ['DATE_OBS','OBJECT','RA','DEC','TIME_OBS','HA','ITIME',$
              'COADDS','FILTER','SAMPMODE','N_MREADS','CHANNEL']
    kwtitle = ['Date','Object','RA','Dec','Time','HA','T_int', $
              'Coadd','Filt','SmpMd','#Rd','Ch']
    if keyword_set(fname) eq 0 then fname = '01noa'
    if keyword_set(ext) eq 0 then ext = '.fts'
endif else if keyword_set(ctio) then begin
    kwlist = ['DATE','OBJECT','RA','DEC','OFFSET','UT','HA','INT_S',$
              'COADDS','FILTER', 'LNRS']
    kwtitle = ['Date','Object','RA','Dec','Offset','UT', 'HA','T_int', $
              'Coadd','Filt', '#Rds']
    if keyword_set(fname) eq 0 then fname = 'ctio'
    if keyword_set(ext) eq 0 then ext = '.fits.'
endif else if keyword_set(scam) then begin   ; NIRSPEC SCAM
    kwlist = ['DATE-OBS', 'OBJECT','FILTER','ITIME', 'COADDS',  $
              'AIRMASS', 'RA','DEC','UTC', $
              'SCAMPA', 'SAMPMODE', 'MULTISCA']
    kwtitle = ['Date(UT)', 'Object','Filter', 'T(int)', 'Coadd',  $
               'AM', 'RA','Dec','Time(UT)',$
              'PA_scam','SmpMd','#Rds']
    if keyword_set(fname) eq 0 then fname = '23sei'
    if keyword_set(ext) eq 0 then ext = '.fits'
endif else if keyword_set(lws) then begin   ; LWS
    kwlist = ['DATE-OBS','UTC', 'FRAMENO','OBJNAME','FILNAME',  $
              'AIRMASS', 'NAXIS4','NAXIS6','CHPCOADD','FRMCOADD', $
              'FRMTIME','RA','DEC']
    kwtitle = ['Date(UT)','Time(UT)', '#','Object','Filter',   $
               'AM','#CHOPsets','#NODsets','CHPCOADD','FRMCOADD', $
               'FRMTIME','RA','Dec']
    if keyword_set(fname) eq 0 then fname = 'lws'
    if keyword_set(ext) eq 0 then ext = '.fits'
endif else if keyword_set(kcam) then begin   ; AO-kcam
    kwlist = ['DATE-OBS','UTC', 'FRAMENO','OBJECT','FILTER','OBKCNAME',  $
              'OBKSNAME','AIRMASS','ITIME','COADDS','WCDMSTAT', $
              'RA','DEC','PARANG']
    kwtitle = ['Date(UT)','Time(UT)', '#','Object','Filter','OBKCNAME',  $
               'OBKSNAME','AM','ITIME','COADD','WCDMSTAT', $
               'RA','Dec','Parang']
    if keyword_set(fname) eq 0 then fname = 'AOKCAM'
    if keyword_set(ext) eq 0 then ext = '.fits'
endif else if keyword_set(nirspec) then begin   ; NIRSPEC
    kwlist = ['FILENAME','UTC','OBJECT','FILNAME','SLITNAME','SLITPA', $
              'ITIME','COADDS','ECHLPOS','DISPPOS','AIRMASS', $
              'NEON','ARGON','KRYPTON','XENON','ETALON','FLAT', $
              'RA','DEC']
    kwtitle = ['Filename','Time(UT)','Object','FilterName',   $
               'Slit','SlitPA','ITIME','COADDS','ECHLPOS','DISPPOS', $
               'AIRMASS','N','A','K','X','E','F','RA','Dec']
    if keyword_set(fname) eq 0 then fname = '23ses'
    if keyword_set(ext) eq 0 then ext = '.fits'
endif else begin
    keck = 1
    kwlist = ['DATE-OBS', 'OBJECT','RA','DEC','UTC','AIRMASS','TINT',$
              'M56COAD0','FILTER','ROTPOSN','ROTPPOSN']
    kwtitle = ['Date(UT)', 'Object','RA','Dec','Time(UT)','AM','T(int)',$
              'Coadds','Filter','User PA','Phys PA']
    if keyword_set(fname) eq 0 then fname = 's'
    if keyword_set(ext) eq 0 then ext = '.fits'
endelse


; provide info
if n_params() lt 2 and not(keyword_set(infile)) then begin
    print,'imlist_new,istart,iend,[infile=],[outfile=]'
    print, '      [fname="', fname, '"],[ext="', ext, $
      '"],[fdir="', fdir, '"],' 
    print, '      [gemini],[lick],' 
    print, '      [silent]'
    return
endif


; if given input file, figure out how many images, 
; not counting blank lines in file or comment lines ('#')
if not(keyword_set(silent)) then print
if keyword_set(infile) then begin
    openr,inunit,infile,/get_lun
    file = ' ' 
    n = 0
    while eof(inunit) ne 1 do begin
        readf,inunit,file
        if (strlen(file) ne 0 and strmid(file,0,1) ne '#') then n = n+1
    endwhile
    if keyword_set(silent) eq 0 then print, '** listing ', $
      strc(n), ' images **'
    free_lun,inunit
    openr,inunit,infile,/get_lun
endif else begin
    n = iend - istart + 1
endelse


; initialize
kwlist = ['',kwlist]
kwtitle = ['Image',kwtitle]
nkw = n_elements(kwtitle)
len = intarr(nkw)           ; # of printed characters in each field
val = strarr(nkw)           ; value of keywords for a given image


; loop through image list
if not(keyword_set(silent)) then print
for i=0,n-1 do begin

    ; reset header variable so if program doesn't find an image, 
    ; it will crash instead of using the last image's header
    h = ''

    if keyword_set(infile) then begin
        got = 0
        repeat begin
            readf,inunit,file
            if (strlen(file) ne 0 and strmid(file,0,1) ne '#') then begin
                im = readfits(fdir+file,h,/silent)
                got = 1
                file = strmid(file,0,strpos(file,ext))
            endif
        endrep until (got)        

    endif else begin

;       construct the file name: NIRC format is default
        if keyword_set(lick) or keyword_set(ctio) then begin
            file = fname+strc(i+istart)
        endif else if keyword_set(gemini) then begin
            file = strc(fname+'0000')
            index = strc(string(i+istart))
            strput,file,index,strlen(file)-strlen(index)
        endif else if keyword_set(ircal) then begin
            file = strc(fname+'0000')
            index = strc(string(i+istart))
            strput,file,index,strlen(file)-strlen(index)
        endif else if keyword_set(scam) then begin
            file = strc(fname+'0000')
            index = strc(string(i+istart))
            strput,file,index,strlen(file)-strlen(index)
        endif else if keyword_set(lws) then begin
            file = strc(fname+'0000')
            index = strc(string(i+istart))
            strput,file,index,strlen(file)-strlen(index)
        endif else if keyword_set(kcam) then begin
            file = strc(fname+'0000')
            index = strc(string(i+istart))
            strput,file,index,strlen(file)-strlen(index)
        endif else if keyword_set(nirspec) then begin
            file = strc(fname+'0000')
            index = strc(string(i+istart))
            strput,file,index,strlen(file)-strlen(index)
        endif else begin
            file = strc(fname+'0000')
            index = strc(string(i+istart))
            strput,file,index,strlen(file)-strlen(index)
        endelse

        filename = fdir+file+ext
        im=readfits(filename,h,/silent)
        if n_elements(im) eq 1 then begin
            message, 'unable to open file "'+filename+'"', /info
            retall
        endif

    endelse


; listing starts with file name/number
  val(0) = file
  if (i eq 0) then begin 
      len(0) = strlen(file) > strlen(kwtitle(0))
      ss = len(0)
      hs = strlen(kwtitle(0))
      hrule = '  ' + strdup('=',ss)
      form = '("| ", A'+strc(ss)+', '
      print, $
        format='($,"  ",A'+strc(hs)+',"'+$
        strdup(' ', strc(ss-hs))+'")',kwtitle(0)
      if keyword_set(outfile) then begin
          if filecheck(outfile) then openw,outunit,outfile,/get_lun
          printf, outunit, $
            format='($,"  ",A'+strc(hs)+',"'+$
            strdup(' ', strc(ss-hs))+'")',kwtitle(0)
      endif
  endif


; now get the desired keywords
  for j=1,nkw-1 do begin

      vv = ''
      if keyword_set(ctio) then begin
          vv = sxpar_ctio(h,kwlist(j))
      endif else $
        vv = sxpar(h,kwlist(j))

;     for Keck, change RA and DEC from degs to standard form 
;     need to take care to handle values which are in sci notation -
;        use ROUND and STRC to truncate and get rid of extra zeros
      if keyword_set(keck) or keyword_set(scam) then begin
          if (kwlist(j) eq 'RA') then begin
              radec,vv,0.0,ihr,imin,rsec
              rsec = strc(round(rsec*10.)/10.)
              vv = strn(ihr,len=2,padc='0') + ':' + $
                strn(imin,len=2,padc='0') + ':' + $
                strn(rsec,len=4,padc='0')
          endif else if (kwlist(j) eq 'DEC') then begin
              radec,0.0,vv,hr,rmin,rsec,ideg,imn,asec
              asec = strc(round(asec*10.)/10.)
              if (asec lt 10) then asec='0'+strc(asec)
              vv = strn(ideg,len=2,padc='0') + ':' + $
                strn(imn,len=2,padc='0') + ':' + $
                strn(asec,len=4,padc='0') 
          endif
      endif

      val(j) = strc(vv)

;     round any floats of their excess zeros and to two decimal places,
;     putting a zero at the end if there's an exposed decimal pt
      sz = size(vv)
      typ = sz(0)+1
      if (sz(typ) eq 5) or (sz(typ) eq 4) then begin
          ss = strc(round(vv*100.)/100.)
;          ss = val(j)
          while (strmid(ss,strlen(ss)-1,1) eq '0') do $
            ss = strmid(ss,0,strlen(ss)-1)
          if (strmid(ss,strlen(ss)-1,1) eq '.') then ss = ss+'0'
          val(j) = ss
      endif 

;     initialize printing format and title heads
      if (i eq 0) then begin

          len(j) = strlen(val(j)) > strlen(kwtitle(j)) 
;          if (kwlist(j) eq 'OBJECT') then len(j) = strlen(vv)
          if keyword_set(keck) then begin
              if (kwlist(j) eq 'AIRMASS') then len(j) = len(j) > 5
              if (kwlist(j) eq 'OBJECT') then len(j) = 12
              if (kwlist(j) eq 'FILTER') then len(j) = 15
          endif else if keyword_set(ctio) then begin
              if (kwlist(j) eq 'OFFSET') then len(j) = len(j) > 13
          endif
          if (kwlist(j) eq 'AIRMASS') then len(j) = len(j) > 5
          if (kwlist(j) eq 'OBJECT') then len(j) = len(j) > 12
          if (kwlist(j) eq 'RA') then len(j) = len(j) > 10
          if (kwlist(j) eq 'DEC') then len(j) = len(j) > 9
          if (kwlist(j) eq 'HA') then len(j) = len(j) > 9
          ss = len(j)
          hs = strlen(kwtitle(j))
          hrule = hrule + '   ' + strdup('=',ss)
          form = form + '" | ", A'+strc(ss)+', '          
          print, $
            format='($,"   ",A'+strc(hs)+',"'+$
            strdup(' ', strc(ss-hs))+'")',kwtitle(j)
          if keyword_set(outfile) then $
            printf, outunit, $
            format='($,"   ",A'+strc(hs)+',"'+$
            strdup(' ', strc(ss-hs))+'")',kwtitle(j)
      endif 

  endfor
  
  if (i eq 0) then begin
      form = strmid(form,0,strlen(form)-2)+', " |")'
      print
      print,hrule
      if keyword_set(outfile) then begin
          printf,outunit
          printf,outunit,hrule
      endif
  endif

  print,format=form,val
  if keyword_set(outfile) then $
    printf, outunit,format=form,val

endfor

if keyword_set(outfile) then free_lun,outunit
end



