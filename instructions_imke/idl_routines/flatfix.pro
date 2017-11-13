pro flatfix,imin,imout
; make zero's in flat equal to 1; will be corrected with badpixel map later
s = size(imin)
im=imin
imout=im
print,s(1),s(2)
 for i=0,s(1)-1 do begin
   for j=0,s(2)-1 do begin
    if(imin(i,j) lt 0.01) then imout(i,j)=1.0
   endfor
 endfor

return
end

