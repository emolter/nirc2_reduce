pro shif_fft,imin,imout,shif
;  input:  cube of images to shift
;  output: shifted image, last image combined
; shif is array with shifts from Henry's program

s=size(imin)
imout=fltarr(s(1),s(2),s(3))

; a positive xdiff,ydiff shifts Io to smaller numbers in pixels

for i=0,s(3)-1 do begin
xdiff=shif(i,0)
ydiff=shif(i,1)
im=imin(*,*,i)
imout(*,*,i)=poly_2d(im,[[xdiff,0],[1,0]],[[ydiff,1],[0,0]],1) 

endfor

return

end


