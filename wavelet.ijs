NB.coclass 'jwavelet'
NB. dwt_z_ =: dwtL_jwavelet_
NB. idwt_z_ =: idwtL_jwavelet_


NB. following the waveslim package in R
NB.* kind dwt y
NB. -returns nextY;Wavelet
oddx=: ] {~ ([: (] #~ 0 1 $~ #) [: i. [: # ]) -/ [: i. [: # [

dwt=: 4 : 0
'hpf lpf'=.wdict x
 yvals=. hpf oddx y
 (yvals +/ . * hpf);yvals +/ . * lpf
)

NB. * lev kind dwtL y
NB. -returns level decomposed wavelet of y
NB. -with structure wav1;wav2;...;wavL;yL
dwtL =: 1 : 0
:
  'wn yn'=. m dwt y
  wn; (<:x) (m dwtL)^:(x>1) yn
)
   
NB. drot reorients and rotates wavelet reduced data for mult to mother wvs
NB. 13 : '((i.(#y)) +/ (0,((_1 * #y) + }.i.(#{.x)))) { y'  
drot =:   ] {~ ([: i. [: # ]) +/ 0 , (_1 * [: # ]) + [: }. [: i. [: # [: {. [
NB. filtrot reorients mother wavelet for multiplication to reduced data
NB. 13 : '|: |. ((2 %~ (#x)),2) $ |. y'  
filtrot =: [: |: [: |. (2 ,~ 2 %~ [: # [) $ [: |. ]
NB. reducer adds it all together as a simple dyad
NB. 13 : ' (2 # (x drot y)) * ((2*#y)$x)'
reducer=: (2 # drot) * [ $~ 2 * [: # ]

NB.* kind wavelet y
NB. -returns nextY;Wavelet
idwt=: 4 : 0
 'wv yl' =. y 
 'hpf lpf'=. filtrot each wdict x
 +/"1 (lpf reducer yl) + hpf reducer wv
)


NB. * kind&idwtl wavelet
NB. -returns fully reconstructed series from stacked wavelet
idwtL=: 4 : 0
 idw =. [: < [: x&idwt ,
 > idw/y
)

NB. get the high pass from the low pass, return a box list of both
HpLp =: ] ;~ |. * _1 ^ [: i. #
NB. *wdict 'db4'
NB. -dictionary of wavelet coefficients
NB. -defined wavelets: 
NB. -haar,w4
NB. -minimum bandwidth/mb4,8,16,24
NB. -Daubechies extremal phase db4,6,8,12,16,20
NB. -Daubechies least asymmetric la8,16,20
NB. -fk4,6,8,14
wdict=: 3 : 0
 select. y
  case. 'haar' do.
   HpLp 2 # %: % ] 2
  case. 'w4' do.
   ((_1, 3, 3, _1) % 8);   (_1, 3, _3, 1)%8 
  case. 'mb4' do.
   HpLp 0.4801755, 0.8372545, 0.2269312, _0.1301477
  case. 'mb8' do.
   HpLp _0.1673619, 0.01847751, 0.5725771, 0.7351331, 0.2947855, _0.1108673, 0.007106015, 0.06436345
  case. 'mb16' do.
   HpLp _0.0130277, 0.02173677, 0.1136116, _0.0577657, _0.2278359, 0.1188725, 0.6349228, 0.6701646, 0.2345342, _0.05656657, _0.01987986, 0.05474628, _0.02483876, _0.04984698, 0.009620427, 0.005765899
  case. 'mb24' do.
   HpLp _2.132706e_05, 0.0004745736, 0.0007456041, _0.004879053, _0.001482995, 0.04199576, _0.002658282, _0.006559513, 0.1019512, 0.1689456, 0.1243531, 0.1949147, 0.4581101, 0.6176385, 0.2556731, _0.3091111, _0.3622424, _0.004575448, 0.1479342, 0.01027154, _0.01644859, _0.002062335, 0.001193006, 5.361301e_05
  case. 'db4' do.
   HpLp 0.482962913144534, 0.836516303737808, 0.224143868042013,  _0.12940952255126
  case. 'db6' do.
   HpLp 0.332670552950083, 0.806891509311093, 0.459877502118491, _0.135011020010255, _0.0854412738820267, 0.0352262918857096
  case. 'db8' do.
   HpLp 0.230377813307443, 0.714846570548406, 0.630880767935879, _0.0279837694166834, _0.187034811717913, 0.0308413818353661, 0.0328830116666778, _0.0105974017850021
  case. 'db16' do.
   HpLp 0.0544158422431049, 0.312871590914303, 0.67563073629729, 0.585354683654191, _0.0158291052563816, _0.28401554296157, 0.0004724845739124, 0.128747426620484, _0.0173693010018083, _0.0440882539307952, 0.0139810279173995, 0.0087460940474061, _0.0048703529934518, _0.000391740373377, 0.0006754494064506, _0.0001174767841248
  case. 'fk4' do.
   HpLp 0.653927555569765, 0.753272492839487, 0.0531792287790598, _0.0461657148152177
  case. 'fk6' do.
   HpLp 0.42791503242231, 0.812919643136907, 0.356369511070187, _0.146438681272577, _0.0771777574069701, 0.0406258144232379
  case. 'fk8' do.
   HpLp 0.3492381118638, 0.782683620384065, 0.475265135079471, _0.0996833284505732, _0.15997809743403, 0.0431066681065162, 0.0425816316775818, _0.0190001788537359
  case. 'fk14' do.
   HpLp 0.260371769291396, 0.686891477239599, 0.611554653959511, 0.0514216541421191, _0.245613928162192, _0.0485753390858553, 0.124282560921513, 0.0222267396224631, _0.0639973730391417, _0.00507437254997285, 0.029779711590379, _0.00329747915270872, _0.00927061337444824, 0.00351410097043596
  case. 'la8' do.
   HpLp _0.0757657147893567, _0.0296355276459604, 0.497618667632563, 0.803738751805386, 0.297857795605605, _0.0992195435769564, _0.0126039672622638, 0.0322231006040782
  case. 'la16' do.
   HpLp 0.0544158422431049, 0.312871590914303, 0.67563073629729, 0.585354683654191, _0.0158291052563816, _0.28401554296157, 0.0004724845739124, 0.128747426620484, _0.0173693010018083, _0.0440882539307952, 0.0139810279173995, 0.0087460940474061, _0.0048703529934518, _0.000391740373377, 0.0006754494064506,  _0.0001174767841248
  case. 'la20' do. 
   HpLp 0.000770159809103, 9.56326707837e_05, _0.0086412992759401, _0.0014653825833465, 0.0459272392237649, 0.0116098939129724, _0.159494278857531, _0.0708805358108615, 0.471690666842659, 0.769510037014339, 0.383826761225382, _0.0355367403054689, _0.0319900568281631, 0.049994972079156, 0.0057649120455518, _0.020354939803946, _0.000804358934537, 0.0045931735836703, 5.7036084339e_05, _0.0004593294205481
 end.
)
 





NB. NB. note, while this is reasonably efficient, it isn't very general. There also isn't a defined
NB. NB. idft, and though I can guess at a solution, it's better to start from a general solution.


NB. debugging code to figure out how wavelets work ....
ts=: 6!:2, 7!:2@]
cd=: 15!:0
typeof=: 3!:0
LDWT=: '/home/scott/src/lugos/jstuff/wavelet/src/libdwt.so'

dwtc=: 4 : 0
 'hpf lpf'=.wdict x
 Newy=.((#y)%2)#(2.2-2.2)
 Wout=.((#y)%2)#(2.2-2.2)
 cmd=. LDWT, ' dwt n *d i i *d *d *d *d'
 cmd cd (y+2.2-2.2);(<.(#y));(#hpf);hpf;lpf;Wout;Newy
 Wout;Newy
)

dwtLc=: 4 : 0
 'lev k' =. x
 'wn yn'=. k dwtc y
 wn; (((<:lev);k)&dwtLc^:(lev>1) yn)
)

JVAL=: 1
modwtc =: 4 : 0
 'hpf lpf'=.wdict x
NB. Newy=.((#y)%2)#(2.2-2.2)
NB.  Wout=.((#y)%2)#(2.2-2.2)
 Newy=.(#y)#(2.2-2.2)
 Wout=.(#y)#(2.2-2.2)
 cmd=. LDWT, ' modwt n *d i i i *d *d *d *d'
 cmd cd (y + 2.2 - 2.2);(<.(#y));(JVAL);(#hpf);hpf;lpf;Wout;Newy
 Wout;Newy
)

modwtLc =: 4 : 0
 'lev k'=. x
 'wn yn'=. k modwtc y
 wn; (((<:lev);k)&modwtLc^:(lev>1) yn)
)

imodwtc =: 4 : 0
 'wv yl' =. y
 'hpf lpf'=.wdict x
 Newy=.(#yl)#(2.2-2.2) 
 cmd=. LDWT, ' imodwt n *d i i i *d *d *d *d'
 cmd cd wv;yl;(#yl);(JVAL);(#hpf);hpf;lpf;Newy
 Newy;Wout
)

imodwtLc=: 4 : 0
idw =. [: < [: x&imodwtc ,
 > idw/y
)

idwtc=: 4 : 0
 'wv yl'=. y
 'hpf lpf'=.wdict x
 Newy=.((#yl)*2)#(2.2-2.2) 
 cmd=. LDWT, ' idwt n *d *d i i *d *d *d'
 cmd cd wv;yl;(#yl);(#hpf);hpf;lpf;Newy
 Newy
)


idwtLc =: 4 : 0
 idw =. [: < [: x&idwtc ,
 > idw/y
)









