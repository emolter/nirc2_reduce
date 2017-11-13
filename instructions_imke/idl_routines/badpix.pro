pro badpix,imin,imout,hpmask
; applies the badpixel map to the image

if n_params() eq 0 then begin
        print,'badpix,imin,imout,hpmask'
        return
endif

s=size(imin)
; s(1) and s(2) is size array

if(s(0) eq 2) then s(3)=1
print,s(0),s(1),s(2),s(3)

 temp=fltarr(s(1),s(2),s(3))
 imout=fltarr(s(1),s(2),s(3))
; first set all bad points to 0.0

 temp=imin

  for k=0,s(3)-1 do begin
; print,k
    for i=0,s(1)-1 do begin
      for j=0,s(2)-1 do begin
         if(hpmask(i,j) eq 0) then temp(i,j,k)=0.0
          if(hpmask(i,j) eq 0) then begin
            if(i ge 2 and i le s(1)-3) then begin
              if(j ge 2 and j le s(2)-3) then begin
                 temp(i,j,k)=median(temp(i-2:i+2,j-2:j+2,k))
         endif
        endif
		if(i eq 0 or i eq s(1)-1) then begin
		  if(j ge 1 and j le s(2)-2) then temp(i,j,k)=median(temp(I,j-1:j+1,k))
		endif
		if(j eq 0 or j eq s(2)-1) then begin
		  if(i ge 1 and i le s(1)-2) then temp(i,j,k)=median(temp(I-1:i+1,j,k))
		endif
       endif
      endfor
    endfor
  endfor

   imout=temp

 end
