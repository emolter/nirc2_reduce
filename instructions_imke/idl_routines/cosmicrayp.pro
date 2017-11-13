pro cosmicray,imin,imout
; determines if pixel is way positive and replaces with median neighborhood
s=size(imin)
 imout=imin
 if(s(0) eq 2) then s(3)=1
print,s(1),s(2),s(3)
  for k=0,s(3)-1 do begin
   for i=2,s(2)-3 do begin
    for j=2,s(1)-3 do begin
      x1=median(imin(j-2:j+2,i-2:i+2,k))
	if(imin(j,i,k) gt 20.0 and imin(j,i,k) gt 5.0*x1)then imout(j,i,k)=x1
	if(imin(j,i,k) gt 20.0 and imin(j,i,k) gt 5.0*x1) then print,j,i,k,imin(j,i,k),x1
;	if(imin(j,i,k) gt 100.0 and imin(j,i,k) gt 3.0*x1)then imout(j,i,k)=x1
;	if(imin(j,i,k) gt 100.0 and imin(j,i,k) gt 3.0*x1) then print,j,i,k,imin(j,i,k),x1
   endfor
   endfor
 endfor
 return
 end
