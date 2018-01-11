load'wavelet.ijs'
temp =. (64#(0,1,2,_1)) + (1 o. i.256) + (1 o. (i.256 ) % o.1 ) NB. interesting TS
temp =. 256$(1,2,3,4)
cwav=.(4;'db4')dwtLc temp
wav=. 4 'db4'  dwtL temp
(> wav) - >cwav NB. c version and J version tie out
'db4' idwtLc cwav
'db4' idwtL wav

wv=. 'db8' 
jj=. ts 'b=: 2 wv dwtL i.2^20'
c=. ts 'b=: (2;wv) dwtLc i.2^20'
jj%c

wt =. 'db4'
'hpf lpf'=.wdict wt
temp =. 16$(1,2,3,4)
wav=.(1;wt)dwtL temp
yl=.ydec=.>1{wav
wv=.wavl=.>0{wav
wt idwtLc  wav


oddx =: ] {~ ([: (] #~ 0 1 $~ #) [: i. [: # ]) -/ [: i. [: # [

jj=. ts 'b=: 2 wv dwtL i.2^20'
jj=.ts 'a=.wv idwtL 2 wv dwtL i.2^20'
c=.ts 'a=.wv idwtLc (2;wv) dwtLc i.2^20'

NB. test modwtLc
temp =. (64#(0,1,2,_1)) + (1 o. i.256) + (1 o. (i.256 ) % o.1 ) NB. interesting TS
cwav =. (4;'db4') modwtLc temp
'db4' imodwtLc cwav