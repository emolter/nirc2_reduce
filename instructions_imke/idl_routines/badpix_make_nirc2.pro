pro hpmake,imin,imout,hpmask,max,min
; mask out high point in frames; look with plothist at levels; make mask
; works only on normalized images (1 average)

if n_params() eq 0 then begin
        print,'hpmake,imin,imout,hpmask,max,min'
        return
endif
s=size(imin)
; s(1) and s(2) is size array

print,s(1),s(2),s(3)

k1=220  ;to trace bad stripe
k2=370

 temp=fltarr(s(1),s(2))

imout=fltarr(s(1),s(2))
mask1=intarr(s(1),s(2))
hpmask=intarr(s(1),s(2))
 n=s(1)
 zy=findgen(4)
im=imin
 x=median(imin)
 print,x
; y=1.25
; z=0.75
 y=max
 z=min
 print,y,z
imout=im
 for i=2,n-3 do begin
   for j=2,(n-3) do begin
    if(imin(i,j) gt y*median(imin(i-2:i+2,j-2:j+2))) then imout(i,j)=0.0
    if(imin(i,j) lt z*median(imin(i-2:i+2,j-2:j+2))) then imout(i,j)=0.0
   endfor
 endfor
; add other edges
 for i=0,1 do begin
   for j=2,(n-3) do begin
	if(imin(i,j) gt y*median(imin(0:i+2,j-2:j+2))) then imout(i,j)=0.0
	endfor
   for j=2,(n-3) do begin
	if(imin(i,j) lt z*median(imin(0:i+2,j-2:j+2))) then imout(i,j)=0.0
	endfor
 endfor
 for i=(n-3),(n-1) do begin
 for j=2,(n-3) do  begin
	if(imin(i,j) gt y*median(imin(n-3:n-1,j-2:j+2))) then imout(i,j)=0.0
	endfor
 for j=2,(n-3) do begin 
	if(imin(i,j) lt z*median(imin(n-3:n-1,j-2:j+2))) then imout(i,j)=0.0
	endfor
 endfor
 for j=0,1 do begin
 for i=2,(n-3) do begin
	if(imin(i,j) gt y*median(imin(i-2:i+2,0:1))) then imout(i,j)=0.0
	endfor
 for i=2,(n-3) do begin
	if(imin(i,j) lt z*median(imin(i-2:i+2,0:1))) then imout(i,j)=0.0
	endfor
 endfor
 for j=(n-3),(n-1) do begin
 for i=2,(n-3) do begin
	if(imin(i,j) gt y*median(imin(i-2:i+2,n-3:n-1))) then imout(i,j)=0.0
	endfor
 for i=2,(n-3) do begin
	if(imin(i,j) lt z*median(imin(i-2:i+2,n-3:n-1))) then imout(i,j)=0.0
	endfor
 endfor

 for j=0,k2 do begin
 for i=0,k1 do begin
	if(imin(i,j) lt 0.97) then imout(i,j)=0.0
	endfor
endfor
imout(213,360)=0.
imout(212,359)=0.
imout(212,358)=0.
imout(211,357)=0.
imout(211,356)=0.
imout(210,355)=0.
imout(210,354)=0.
imout(209,354)=0.
imout(209,353)=0.
imout(208,352)=0.
imout(208,351)=0.
imout(207,350)=0.
imout(206,349)=0.
imout(206,348)=0.
imout(205,347)=0.
imout(205,346)=0.
imout(204,345)=0.
imout(203,344)=0.
imout(203,343)=0.
imout(202,342)=0.
imout(202,341)=0.
imout(201,340)=0.
imout(200,339)=0.
imout(200,338)=0.
imout(199,337)=0.
imout(199,336)=0.
imout(198,335)=0.
imout(197,334)=0.
imout(197,333)=0.
imout(196,332)=0.
imout(196,331)=0.
imout(195,330)=0.
imout(194,329)=0.
imout(194,328)=0.
imout(193,327)=0.
imout(193,326)=0.
imout(192,325)=0.
imout(191,324)=0.
imout(191,323)=0.
imout(190,322)=0.
imout(190,321)=0.
imout(189,320)=0.
imout(188,319)=0.
imout(188,318)=0.
imout(187,317)=0.
imout(186,316)=0.
imout(186,315)=0.
imout(185,314)=0.
imout(185,313)=0.
imout(184,312)=0.
imout(184,311)=0.
imout(183,311)=0.
imout(183,310)=0.
imout(182,309)=0.
imout(182,308)=0.
imout(181,307)=0.
imout(181,306)=0.
imout(180,305)=0.
imout(179,304)=0.
imout(179,303)=0.
imout(178,302)=0.
imout(177,301)=0.
imout(177,300)=0.
imout(176,299)=0.
imout(176,298)=0.
imout(175,297)=0.
imout(174,296)=0.
imout(174,295)=0.
imout(173,294)=0.
imout(173,293)=0.
imout(172,292)=0.
imout(171,291)=0.
imout(171,290)=0.
imout(170,289)=0.
imout(170,288)=0.
imout(169,287)=0.
imout(168,286)=0.
imout(168,285)=0.
imout(167,284)=0.
imout(167,283)=0.
imout(166,282)=0.
imout(165,281)=0.
imout(165,280)=0.
imout(164,279)=0.
imout(164,278)=0.
imout(163,277)=0.
imout(162,277)=0.
imout(162,276)=0.
imout(162,275)=0.
imout(161,274)=0.
imout(161,273)=0.
imout(160,272)=0.
imout(159,271)=0.
imout(159,270)=0.
imout(158,269)=0.
imout(158,268)=0.
imout(157,267)=0.
imout(157,268)=0.
imout(156,266)=0.
imout(156,265)=0.
imout(155,264)=0.
imout(155,263)=0.
imout(154,263)=0.
imout(154,262)=0.
imout(153,261)=0.
imout(153,260)=0.
imout(152,259)=0.
imout(151,258)=0.
imout(151,257)=0.
imout(150,256)=0.
imout(150,255)=0.
imout(149,254)=0.
imout(149,253)=0.
imout(148,253)=0.
imout(148,252)=0.
imout(147,251)=0.
imout(147,250)=0.
imout(146,249)=0.
imout(145,248)=0.
imout(145,247)=0.
imout(144,246)=0.
imout(144,245)=0.
imout(144,35)=0.
imout(143,244)=0.
imout(142,243)=0.
imout(142,242)=0.
imout(141,241)=0.
imout(140,239)=0.
imout(139,238)=0.
imout(139,237)=0.
imout(138,236)=0.
imout(138,235)=0.
imout(137,234)=0.
imout(136,233)=0.
imout(135,231)=0.
imout(135,230)=0.
imout(134,229)=0.
imout(133,228)=0.
imout(133,227)=0.
imout(132,226)=0.
imout(132,225)=0.
imout(131,224)=0.
imout(130,223)=0.
imout(130,222)=0.
imout(129,221)=0.
imout(129,220)=0.
imout(127,218)=0.
imout(127,217)=0.
imout(126,215)=0.
imout(126,216)=0.
imout(125,215)=0.
imout(125,214)=0.
imout(124,213)=0.
imout(124,212)=0.
imout(123,211)=0.
imout(123,210)=0.
imout(122,210)=0.
imout(122,209)=0.
imout(121,208)=0.
imout(121,207)=0.
imout(120,205)=0.
imout(119,205)=0.
imout(119,204)=0.
imout(118,203)=0.
imout(118,201)=0.
imout(117,201)=0.
imout(116,200)=0.
imout(116,199)=0.
imout(115,198)=0.
imout(115,197)=0.
imout(114,196)=0.
imout(113,195)=0.
imout(113,194)=0.
imout(112,193)=0.
imout(111,191)=0.
imout(110,190)=0.
imout(110,189)=0.
imout(109,188)=0.
imout(109,187)=0.
imout(108,186)=0.
imout(107,185)=0.
imout(107,184)=0.
imout(106,183)=0.
imout(106,182)=0.
imout(105,181)=0.
imout(104,180)=0.
imout(104,179)=0.
imout(103,178)=0.
imout(103,177)=0.
imout(102,176)=0.
imout(101,175)=0.
imout(101,174)=0.
imout(100,173)=0.
imout(100,172)=0.
imout(99,172)=0.
imout(99,171)=0.
imout(98,170)=0.
imout(98,169)=0.
imout(97,168)=0.
imout(97,167)=0.
imout(96,167)=0.
imout(96,166)=0.
imout(95,165)=0.
imout(95,164)=0.
imout(94,163)=0.
imout(94,162)=0.
imout(93,161)=0.
imout(92,160)=0.
imout(92,159)=0.
imout(91,158)=0.
imout(91,157)=0.
imout(90,157)=0.
imout(90,156)=0.
imout(89,155)=0.
imout(89,154)=0.
imout(87,152)=0.
imout(87,151)=0.
imout(86,150)=0.
imout(86,149)=0.
imout(85,148)=0.
imout(84,147)=0.
imout(84,146)=0.
imout(83,145)=0.
imout(83,144)=0.
imout(82,143)=0.
imout(81,142)=0.
imout(81,141)=0.
imout(80,140)=0.
imout(79,138)=0.
imout(78,137)=0.
imout(78,136)=0.
imout(77,135)=0.
imout(77,134)=0.
imout(76,133)=0.
imout(75,132)=0.
imout(75,131)=0.
imout(74,130)=0.
imout(74,129)=0.
imout(73,128)=0.
imout(71,125)=0.
imout(71,124)=0.
imout(70,124)=0.
imout(70,123)=0.
imout(69,122)=0.
imout(69,121)=0.
imout(68,120)=0.
imout(68,119)=0.
imout(67,118)=0.
imout(67,118)=0.
imout(66,117)=0.
imout(66,116)=0.
imout(65,115)=0.
imout(65,114)=0.
imout(64,113)=0.
imout(63,112)=0.
imout(63,111)=0.
imout(62,110)=0.
imout(61,109)=0.
imout(60,107)=0.
imout(60,106)=0.
imout(59,105)=0.
imout(58,104)=0.
imout(58,103)=0.
imout(57,102)=0.
imout(57,101)=0.
imout(56,100)=0.
imout(55,99)=0.
imout(55,98)=0.
imout(54,97)=0.
imout(54,96)=0.
imout(52,94)=0.
imout(52,93)=0.
imout(51,92)=0.
imout(51,91)=0.
imout(50,90)=0.
imout(49,89)=0.
imout(49,88)=0.
imout(48,87)=0.
imout(48,86)=0.
imout(47,86)=0.
imout(47,85)=0.
imout(46,84)=0.
imout(46,83)=0.
imout(45,82)=0.
imout(45,81)=0.
imout(44,81)=0.
imout(44,80)=0.
imout(43,79)=0.
imout(43,78)=0.
imout(42,77)=0.
imout(42,76)=0.
imout(41,76)=0.
imout(41,75)=0.
imout(40,74)=0.
imout(39,72)=0.
imout(38,71)=0.
imout(38,70)=0.
imout(37,69)=0.
imout(37,68)=0.
imout(36,67)=0.
imout(35,66)=0.
imout(35,65)=0.
imout(34,64)=0.
imout(34,63)=0.
imout(33,62)=0.
imout(32,61)=0.
imout(32,60)=0.
imout(31,59)=0.
imout(31,58)=0.
imout(30,57)=0.
imout(29,56)=0.
imout(29,55)=0.
imout(28,54)=0.
imout(28,53)=0.
imout(27,52)=0.
imout(26,51)=0.
imout(26,50)=0.
imout(25,49)=0.
imout(25,48)=0.
imout(24,48)=0.
imout(24,47)=0.
imout(23,46)=0.
imout(23,45)=0.
imout(22,44)=0.
imout(22,43)=0.
imout(21,42)=0.
imout(20,41)=0.
imout(20,40)=0.
imout(19,39)=0.
imout(19,38)=0.
imout(18,38)=0.
imout(18,37)=0.
imout(17,36)=0.
imout(17,35)=0.
imout(16,34)=0.
imout(16,33)=0.
imout(15,33)=0.
imout(15,32)=0.
imout(14,31)=0.
imout(14,30)=0.
imout(13,29)=0.
imout(12,28)=0.
imout(12,27)=0.
imout(11,26)=0.
imout(11,25)=0.
imout(10,24)=0.
imout(9,23)=0.
imout(9,22)=0.
imout(7,19)=0.
imout(7,11)=0.
imout(7,10)=0.


; makemask for proper weighting

mask1=abs(imout)/(abs(imout)>.0001)   ; =1 for nonzero =0 for zero
hpmask=mask1


return
end

