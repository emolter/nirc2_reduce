pro cosmicneg,imin,imout
; determines if pixel is way negative and replaces with median neighborhood
s=size(imin)
 imout=imin
 if(s(0) eq 2) then s(3)=1
print,s(1),s(2),s(3)
  for k=0,s(3)-1 do begin
   for i=2,s(2)-3 do begin
    for j=4,s(1)-3 do begin
      x1=median(imin(j-2:j+2,i-2:i+2,k))
;        if(imin(j,i,k) lt -.20 and imin(j,i,k) lt 3.*x1 ) then  print,j,i,k,imin(j,i,k),x1
;	if(imin(j,i,k) lt -.2 and imin(j,i,k) lt 3.*x1 ) then imout(j,i,k)=x1
	if(x1 gt 20.0 and imin(j,i,k) lt 0.5*x1) then imout(j,i,k)=x1
	if(x1 gt 20. and imin(j,i,k) lt 0.5*x1) then  print,j,i,k,imin(j,i,k),x1
	if(imin(j,i,k) gt 100.0 and imin(j,i,k) lt 0.8*x1) then imout(j,i,k)=x1
	if(imin(j,i,k) gt 100. and imin(j,i,k) lt 0.8*x1) then  print,j,i,k,imin(j,i,k),x1
;	if(imin(j,i,k) gt 100.0 and imin(j,i,k) gt 2.0*x1) then imout(j,i,k)=x1
;	if(imin(j,i,k) gt 100. and imin(j,i,k) gt 2.0*x1) then  print,j,i,k,imin(j,i,k),x1
	if(imin(j,i,k) lt -50.0 and x1 gt 0) then imout(j,i,k)=x1
	if(imin(j,i,k) lt -50. and x1 gt 0) then  print,j,i,k,imin(j,i,k),x1
;	if(imin(j,i,k) lt -100.0) then imout(j,i,k)=x1
;	if(imin(j,i,k) lt -100.) then  print,j,i,k,imin(j,i,k),x1

     endfor
   endfor
 endfor
 return
 end
