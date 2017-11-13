pro nirc2dewarp_ext,imin,imout
; program for narrow camera only
; doesn't seem to work on 512x512 images

s=size(imin)
if (s(0) eq 2) then s(3)=1
print,s(1),s(2),s(3)

imout=imin

for i=0,s(3)-1 do begin
im=imin(*,*,i)
out=nirc2dewarp(im,hd=hd,camera='narrow')

imout(*,*,i)=out(*,*)
endfor

return
end

