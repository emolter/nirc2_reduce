function find_offsets_imstack_fft, ims, verbose = verbose

;written 020503 by hroe from find_offset_fft

;return is [numim,2] where 2=x/y offsets to shift each image by

;note that program is a memory hog.  This is to save recalculating
;fft's more than necessary.

;note that currently all math is done as DOUBLE...later should change
;this to optional

s = size(ims)
numim = s[3]
xs = s[1]
ys = s[2]

ftims = dcomplexarr(xs, ys, numim)
for i = 0, numim-1 do ftims[*, *, i] = fft(ims[*, *, i], -1, /double)

Xoffsets = dblarr(numim, numim)
Yoffsets = dblarr(numim, numim)

for a = 0, numim-1 do for b = a, numim-1 do begin 
 cc = double(fft(ftims[*, *, a]*conj(ftims[*, *, b]), 1, /double))
 where_xy, cc, where(cc eq max(cc)), xoff, yoff
 if xoff[0] ge xs/2 then xoff[0] = xoff[0] - xs
 if yoff[0] ge ys/2 then yoff[0] = yoff[0] - ys
 Xoffsets[a, b] = xoff[0]
 Yoffsets[a, b] = yoff[0]
 Xoffsets[b, a] = -xoff[0]
 Yoffsets[b, a] = -yoff[0]
 if keyword_set(verbose) then print, "  done with ", a, b
endfor

xoffs = solve_for_offsets(xoffsets)
yoffs = solve_for_offsets(yoffsets)

return, [[xoffs - xoffs[0]],[yoffs - yoffs[0]]]

end
