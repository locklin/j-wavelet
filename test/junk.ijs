load'wavelet.ijs'

dwtL=: 4 : 0
 'lev k'=.x
 'wn yn'=. k dwt y
 wn; (((<:lev);k)&dwtL^:(lev>1) yn)
)


dwtdumb=: 4 : 0
 'hpf lpf'=.wdict x
 Wout=.Newy=. '' NB. ((#y)%2)#0
 for_tt. i.((#y)%2) do.
  u =. 1+ 2 * tt 
NB.  smoutput 'outer=', (":u) ,' t=', ":tt
  W=. ({.hpf) * (u{y)
  N=. ({.lpf) * (u{y)
   for_n. (1+i.((#hpf)-1)) do.
    u =. <: u 
    if. u < 0 do.
     u=. (#y) -1
    end.
NB.    smoutput 'inner=',(":u),' t=',(":tt),' n=',":n
    W=.W+(n{hpf)* (u{y)
    N=.N+(n{lpf)* (u{y)
   end.
   Wout=.Wout,W
   Newy=.Newy,N
   W=.N=.0
  end.
Newy;Wout
)

dwtdL =: 4 : 0
 'lev k'=.x
 'yn wn'=. k dwtdumb y
 wn; (((<:lev);k)&dwtdL^:(lev>1) yn)
)

NB. Daubechies wavelets, thx to Roger Hui on the J-list
NB. http://jsoftware.2058.n7.nabble.com/J-vs-Matlab-vs-NumPy-for-speed-by-wavelet-computing-td8390.html
D4a=: 3 : 0
C1=.  1.7320508075688772
C2=.  0.4330127018922193
C3=. -0.066987298107780702
C4=.  0.51763809020504137
C5=.  1.9318516525781364
 e=.  y#~(#y)$1 0
 d=. (y#~(#y)$0 1) - C1*e
 s=. e+(C2*d)+C3*1|.d
 , (D4a^:(2<#y) C5*s) ,. C4*d+_1|.s
) 



NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.NB.
NB. Filters
NB. takes the last power of 2 rows of y
lastPow2 =: ] {.~ _1 * 2 ^ [: <. 2 ^.  $
firstPow2=: ] {.~ 2 ^ [: <. 2 ^.  $

NB. *(level;'kind') wavSamp yin
NB. -wavelet samples of last 2^n of yin
NB. -default
wavSamp =: 3 : 0
 (6;'db4') wavSamp y
:
'lev knd'=.x
yy=. lastPow2 y
wav=.((lev;knd)&dwtL yy)
(_1*2^((2^.(#yy)) - lev)){. wav
)

NB. same as above, but with the first 2^n of yin
wavSamp2 =: 3 : 0
 (6;'db4') wavSamp y
:
'lev knd'=.x
yy=. firstPow2 y
wav=.((lev;knd)&dwtL yy)
(_1*2^((2^.(#yy)) - lev)){. wav
)

NB. subsample, returns 1st 2n wavelet decomp, then last 2^2 wdecomp
waveletSS =: 4 : 0
NB.(x&wavSamp2"1 |: y),(x&wavSamp"1 |: y)
(x&wavSamp"1 |: y)
)


temp =. (64#(0,1,2,_1)) + (1 o. i.256) + (1 o. (i.256 ) % o.1 ) NB. interesting TS
temp =. 256$(1,2,3,4)
cwav=.(4;'db4')dwtLc temp
wav=.(4;'db4')dwtL temp
(> wav) - >cwav


wt =. 'db4'
'hpf lpf'=.wdict wt
temp =. 16$(1,2,3,4)
wav=.(1;wt)dwtL temp
yl=.ydec=.>1{wav
wv=.wavl=.>0{wav
wt idwtLc  wav

NB. db4 
NB. outer m=0 n=1 i=1 j=0 u=0
NB. outer m=2 n=3 i=1 j=0 u=1
NB. outer m=4 n=5 i=1 j=0 u=2
NB. outer m=6 n=7 i=1 j=0 u=3
NB. outer m=8 n=9 i=1 j=0 u=4
NB. outer m=10 n=11 i=1 j=0 u=5
NB. outer m=12 n=13 i=1 j=0 u=6
NB. outer m=14 n=15 i=1 j=0 u=7

NB. inner m=0 n=1 i=3 j=2 u=1
NB. inner m=2 n=3 i=3 j=2 u=2
NB. inner m=4 n=5 i=3 j=2 u=3
NB. inner m=6 n=7 i=3 j=2 u=4
NB. inner m=8 n=9 i=3 j=2 u=5
NB. inner m=10 n=11 i=3 j=2 u=6
NB. inner m=12 n=13 i=3 j=2 u=7
NB. inner m=14 n=15 i=3 j=2 u=0

dro =: 2 # ]
outlh =. 16$(1,0){hpf    NB.  outer loop
outll =. 16$(1,0){lpf    NB.  outer loop
inlh=. 16$(3,2){hpf   NB. inner loop
inll=. 16$(3,2){lpf   NB. inner loop
outer =. ((dro ydec) * outll) + ((dro wavl) * outlh) NB. outer loop
13 : ' (2|. 2 # i.(#y)){y '
drot =: ] {~ 2 |. 2 # ]
inner =. ((drot ydec) * inll) + ((drot wavl) * inlh) NB. outer loop
outer + inner
NB. ja, ja, jawol! -but this only works for size 4 wavelets...


drotg =: ' (2 # i.(#y)){y '
((i.16) +/ (0,_2)) { (2#i.8)
((i.8) +/ (0,_1)) { i.8
((i.8) +/ (0 _7)) {i.8
((i.8) +/ (0 _7 _6)) {i.8
13 : '((i.(#y)) +/ (0,((_1 * #y) + }.i.(#x)))) { y'
drotx =: ] {~ ([: i. [: # ]) +/ 0 , (_1 * [: # ]) + [: }. [: i. [: # [
13 : '((i.(#y)) +/ (0,((_1 * #y) + }.i.(#{.x)))) { y'
drotx =:   ] {~ ([: i. [: # ]) +/ 0 , (_1 * [: # ]) + [: }. [: i. [: # [: {. [

|: |.  3 2 $ |.i.6
|: |.  2 2 $ |.i.4
13 : '|: |. ((2 %~ (#x)),2) $ |. y'
filtrot =: [: |: [: |. (2 ,~ 2 %~ [: # [) $ [: |. ]
13 : ' (2 # (x drotx y)) * ((2*#y)$x)'
reducer=: (2 # drotx) * [ $~ 2 * [: # ]

idwt=: 4 : 0
 'wv yl' =. y 
 'hpf lpf'=. filtrot each wdict x
 +/"1 (lpf reducer yl) + (hpf reducer wv)
)
 

NB. db6
wt =. 'db6'
'hpf lpf'=.wdict wt
temp =. 16$(1,2,3,4)
wav=.(1;wt)dwtL temp
wt idwtLc  wav
ydec=.>1{wav
wavl=.>0{wav

NB. outer m=0 n=1 i=1 j=0 u=0
NB. outer m=2 n=3 i=1 j=0 u=1
NB. outer m=4 n=5 i=1 j=0 u=2
NB. outer m=6 n=7 i=1 j=0 u=3
NB. outer m=8 n=9 i=1 j=0 u=4
NB. outer m=10 n=11 i=1 j=0 u=5
NB. outer m=12 n=13 i=1 j=0 u=6
NB. outer m=14 n=15 i=1 j=0 u=7
outlh =. 16$(1,0){hpf    NB.  outer loop
outll =. 16$(1,0){lpf    NB.  outer loop
outer =. ((dro ydec) * outll) + ((dro wavl) * outlh) NB. outer loop

NB. inner m=0 n=1 i=3 j=2 u=1
NB. inner m=0 n=1 i=5 j=4 u=2
NB. inner m=2 n=3 i=3 j=2 u=2
NB. inner m=2 n=3 i=5 j=4 u=3
NB. inner m=4 n=5 i=3 j=2 u=3
NB. inner m=4 n=5 i=5 j=4 u=4
NB. inner m=6 n=7 i=3 j=2 u=4
NB. inner m=6 n=7 i=5 j=4 u=5
NB. inner m=8 n=9 i=3 j=2 u=5
NB. inner m=8 n=9 i=5 j=4 u=6
NB. inner m=10 n=11 i=3 j=2 u=6
NB. inner m=10 n=11 i=5 j=4 u=7
NB. inner m=12 n=13 i=3 j=2 u=7
NB. inner m=12 n=13 i=5 j=4 u=0
NB. inner m=14 n=15 i=3 j=2 u=0
NB. inner m=14 n=15 i=5 j=4 u=1
inlh=. 16$(3,2){hpf   NB. inner loop
inll=. 16$(3,2){lpf   NB. inner loop


NB. db6
wt =. 'db8'
'hpf lpf'=.wdict wt
temp =. 16$(1,2,3,4)
wav=.(1;wt)dwtL temp
wt idwtLc  wav
ydec=.>1{wav
wavl=.>0{wav

NB. outer m=0 n=1 i=1 j=0 u=0
NB. outer m=2 n=3 i=1 j=0 u=1
NB. outer m=4 n=5 i=1 j=0 u=2
NB. outer m=6 n=7 i=1 j=0 u=3
NB. outer m=8 n=9 i=1 j=0 u=4
NB. outer m=10 n=11 i=1 j=0 u=5
NB. outer m=12 n=13 i=1 j=0 u=6
NB. outer m=14 n=15 i=1 j=0 u=7

NB. inner m=0 n=1 i=3 j=2 u=1
NB. inner m=0 n=1 i=5 j=4 u=2
NB. inner m=0 n=1 i=7 j=6 u=3

NB. inner m=2 n=3 i=3 j=2 u=2
NB. inner m=2 n=3 i=5 j=4 u=3
NB. inner m=2 n=3 i=7 j=6 u=4

NB. inner m=4 n=5 i=3 j=2 u=3
NB. inner m=4 n=5 i=5 j=4 u=4
NB. inner m=4 n=5 i=7 j=6 u=5

NB. inner m=6 n=7 i=3 j=2 u=4
NB. inner m=6 n=7 i=5 j=4 u=5
NB. inner m=6 n=7 i=7 j=6 u=6

NB. inner m=8 n=9 i=3 j=2 u=5
NB. inner m=8 n=9 i=5 j=4 u=6
NB. inner m=8 n=9 i=7 j=6 u=7

NB. inner m=10 n=11 i=3 j=2 u=6
NB. inner m=10 n=11 i=5 j=4 u=7
NB. inner m=10 n=11 i=7 j=6 u=0

NB. inner m=12 n=13 i=3 j=2 u=7
NB. inner m=12 n=13 i=5 j=4 u=0
NB. inner m=12 n=13 i=7 j=6 u=1

NB. inner m=14 n=15 i=3 j=2 u=0
NB. inner m=14 n=15 i=5 j=4 u=1
NB. inner m=14 n=15 i=7 j=6 u=2



wt =. 'db8'
'hpf lpf'=.wdict wt
temp =. 32$(1,2,3,4)
wav=.(1;wt)dwtL temp
wt idwtLc  wav
ydec=.>1{wav
wavl=.>0{wav


NB. outer m=0 n=1 i=1 j=0 u=0
NB. outer m=2 n=3 i=1 j=0 u=1
NB. outer m=4 n=5 i=1 j=0 u=2
NB. outer m=6 n=7 i=1 j=0 u=3
NB. outer m=8 n=9 i=1 j=0 u=4
NB. outer m=10 n=11 i=1 j=0 u=5
NB. outer m=12 n=13 i=1 j=0 u=6
NB. outer m=14 n=15 i=1 j=0 u=7
NB. outer m=16 n=17 i=1 j=0 u=8
NB. outer m=18 n=19 i=1 j=0 u=9
NB. outer m=20 n=21 i=1 j=0 u=10
NB. outer m=22 n=23 i=1 j=0 u=11
NB. outer m=24 n=25 i=1 j=0 u=12
NB. outer m=26 n=27 i=1 j=0 u=13
NB. outer m=28 n=29 i=1 j=0 u=14
NB. outer m=30 n=31 i=1 j=0 u=15


NB. inner m=0 n=1 i=3 j=2 u=1
NB. inner m=0 n=1 i=5 j=4 u=2
NB. inner m=0 n=1 i=7 j=6 u=3
NB. inner m=2 n=3 i=3 j=2 u=2
NB. inner m=2 n=3 i=5 j=4 u=3
NB. inner m=2 n=3 i=7 j=6 u=4
NB. inner m=4 n=5 i=3 j=2 u=3
NB. inner m=4 n=5 i=5 j=4 u=4
NB. inner m=4 n=5 i=7 j=6 u=5
NB. inner m=6 n=7 i=3 j=2 u=4
NB. inner m=6 n=7 i=5 j=4 u=5
NB. inner m=6 n=7 i=7 j=6 u=6
NB. inner m=8 n=9 i=3 j=2 u=5
NB. inner m=8 n=9 i=5 j=4 u=6
NB. inner m=8 n=9 i=7 j=6 u=7
NB. inner m=10 n=11 i=3 j=2 u=6
NB. inner m=10 n=11 i=5 j=4 u=7
NB. inner m=10 n=11 i=7 j=6 u=8
NB. inner m=12 n=13 i=3 j=2 u=7
NB. inner m=12 n=13 i=5 j=4 u=8
NB. inner m=12 n=13 i=7 j=6 u=9
NB. inner m=14 n=15 i=3 j=2 u=8
NB. inner m=14 n=15 i=5 j=4 u=9
NB. inner m=14 n=15 i=7 j=6 u=10
NB. inner m=16 n=17 i=3 j=2 u=9
NB. inner m=16 n=17 i=5 j=4 u=10
NB. inner m=16 n=17 i=7 j=6 u=11
NB. inner m=18 n=19 i=3 j=2 u=10
NB. inner m=18 n=19 i=5 j=4 u=11
NB. inner m=18 n=19 i=7 j=6 u=12
NB. inner m=20 n=21 i=3 j=2 u=11
NB. inner m=20 n=21 i=5 j=4 u=12
NB. inner m=20 n=21 i=7 j=6 u=13
NB. inner m=22 n=23 i=3 j=2 u=12
NB. inner m=22 n=23 i=5 j=4 u=13
NB. inner m=22 n=23 i=7 j=6 u=14
NB. inner m=24 n=25 i=3 j=2 u=13
NB. inner m=24 n=25 i=5 j=4 u=14
NB. inner m=24 n=25 i=7 j=6 u=15
NB. inner m=26 n=27 i=3 j=2 u=14
NB. inner m=26 n=27 i=5 j=4 u=15
NB. inner m=26 n=27 i=7 j=6 u=0
NB. inner m=28 n=29 i=3 j=2 u=15
NB. inner m=28 n=29 i=5 j=4 u=0
NB. inner m=28 n=29 i=7 j=6 u=1
NB. inner m=30 n=31 i=3 j=2 u=0
NB. inner m=30 n=31 i=5 j=4 u=1
NB. inner m=30 n=31 i=7 j=6 u=2





NB. outer m=4 n=5 i=1 j=0 u=2
NB. outer m=6 n=7 i=1 j=0 u=3
NB. outer m=8 n=9 i=1 j=0 u=4
NB. outer m=10 n=11 i=1 j=0 u=5
NB. outer m=12 n=13 i=1 j=0 u=6
NB. outer m=14 n=15 i=1 j=0 u=7
NB. outer m=16 n=17 i=1 j=0 u=8
NB. outer m=18 n=19 i=1 j=0 u=9
NB. outer m=20 n=21 i=1 j=0 u=10
NB. outer m=22 n=23 i=1 j=0 u=11
NB. outer m=24 n=25 i=1 j=0 u=12
NB. outer m=26 n=27 i=1 j=0 u=13
NB. outer m=28 n=29 i=1 j=0 u=14
NB. outer m=30 n=31 i=1 j=0 u=15

NB. outer m=0 n=1 i=1 j=0 u=0
NB. inner m=0 n=1 i=3 j=2 u=1
NB. inner m=0 n=1 i=5 j=4 u=2
NB. inner m=0 n=1 i=7 j=6 u=3

NB. outer m=2 n=3 i=1 j=0 u=1
NB. inner m=2 n=3 i=3 j=2 u=2
NB. inner m=2 n=3 i=5 j=4 u=3
NB. inner m=2 n=3 i=7 j=6 u=4
NB. inner m=4 n=5 i=3 j=2 u=3
NB. inner m=4 n=5 i=5 j=4 u=4
NB. inner m=4 n=5 i=7 j=6 u=5
NB. inner m=6 n=7 i=3 j=2 u=4
NB. inner m=6 n=7 i=5 j=4 u=5
NB. inner m=6 n=7 i=7 j=6 u=6
NB. inner m=8 n=9 i=3 j=2 u=5
NB. inner m=8 n=9 i=5 j=4 u=6
NB. inner m=8 n=9 i=7 j=6 u=7
NB. inner m=10 n=11 i=3 j=2 u=6
NB. inner m=10 n=11 i=5 j=4 u=7
NB. inner m=10 n=11 i=7 j=6 u=8
NB. inner m=12 n=13 i=3 j=2 u=7
NB. inner m=12 n=13 i=5 j=4 u=8
NB. inner m=12 n=13 i=7 j=6 u=9
NB. inner m=14 n=15 i=3 j=2 u=8
NB. inner m=14 n=15 i=5 j=4 u=9
NB. inner m=14 n=15 i=7 j=6 u=10
NB. inner m=16 n=17 i=3 j=2 u=9
NB. inner m=16 n=17 i=5 j=4 u=10
NB. inner m=16 n=17 i=7 j=6 u=11
NB. inner m=18 n=19 i=3 j=2 u=10
NB. inner m=18 n=19 i=5 j=4 u=11
NB. inner m=18 n=19 i=7 j=6 u=12
NB. inner m=20 n=21 i=3 j=2 u=11
NB. inner m=20 n=21 i=5 j=4 u=12
NB. inner m=20 n=21 i=7 j=6 u=13
NB. inner m=22 n=23 i=3 j=2 u=12
NB. inner m=22 n=23 i=5 j=4 u=13
NB. inner m=22 n=23 i=7 j=6 u=14
NB. inner m=24 n=25 i=3 j=2 u=13
NB. inner m=24 n=25 i=5 j=4 u=14
NB. inner m=24 n=25 i=7 j=6 u=15
NB. inner m=26 n=27 i=3 j=2 u=14
NB. inner m=26 n=27 i=5 j=4 u=15
NB. inner m=26 n=27 i=7 j=6 u=0
NB. inner m=28 n=29 i=3 j=2 u=15
NB. inner m=28 n=29 i=5 j=4 u=0
NB. inner m=28 n=29 i=7 j=6 u=1
NB. inner m=30 n=31 i=3 j=2 u=0
NB. inner m=30 n=31 i=5 j=4 u=1
NB. inner m=30 n=31 i=7 j=6 u=2





NB. outer m=0 n=1 i=1 j=0 u=0
NB. inner m=26 n=27 i=7 j=6 u=0
NB. inner m=28 n=29 i=5 j=4 u=0
NB. inner m=30 n=31 i=3 j=2 u=0

NB. outer m=2 n=3 i=1 j=0 u=1
NB. inner m=0 n=1 i=3 j=2 u=1
NB. inner m=28 n=29 i=7 j=6 u=1
NB. inner m=30 n=31 i=5 j=4 u=1


NB. outer m=4 n=5 i=1 j=0 u=2
NB. inner m=0 n=1 i=5 j=4 u=2
NB. inner m=2 n=3 i=3 j=2 u=2
NB. inner m=30 n=31 i=7 j=6 u=2

NB. outer m=6 n=7 i=1 j=0 u=3
NB. inner m=0 n=1 i=7 j=6 u=3
NB. inner m=2 n=3 i=5 j=4 u=3
NB. inner m=4 n=5 i=3 j=2 u=3

NB. outer m=8 n=9 i=1 j=0 u=4
NB. inner m=2 n=3 i=7 j=6 u=4
NB. inner m=4 n=5 i=5 j=4 u=4
NB. inner m=6 n=7 i=3 j=2 u=4

NB. outer m=10 n=11 i=1 j=0 u=5
NB. inner m=4 n=5 i=7 j=6 u=5
NB. inner m=6 n=7 i=5 j=4 u=5
NB. inner m=8 n=9 i=3 j=2 u=5

NB. outer m=12 n=13 i=1 j=0 u=6
NB. inner m=6 n=7 i=7 j=6 u=6
NB. inner m=8 n=9 i=5 j=4 u=6
NB. inner m=10 n=11 i=3 j=2 u=6

NB. outer m=14 n=15 i=1 j=0 u=7
NB. inner m=8 n=9 i=7 j=6 u=7
NB. inner m=10 n=11 i=5 j=4 u=7
NB. inner m=12 n=13 i=3 j=2 u=7

NB. outer m=16 n=17 i=1 j=0 u=8
NB. inner m=10 n=11 i=7 j=6 u=8
NB. inner m=12 n=13 i=5 j=4 u=8
NB. inner m=14 n=15 i=3 j=2 u=8

NB. outer m=18 n=19 i=1 j=0 u=9
NB. inner m=12 n=13 i=7 j=6 u=9
NB. inner m=14 n=15 i=5 j=4 u=9
NB. inner m=16 n=17 i=3 j=2 u=9

NB. outer m=20 n=21 i=1 j=0 u=10
NB. inner m=14 n=15 i=7 j=6 u=10
NB. inner m=16 n=17 i=5 j=4 u=10
NB. inner m=18 n=19 i=3 j=2 u=10

NB. outer m=22 n=23 i=1 j=0 u=11
NB. inner m=16 n=17 i=7 j=6 u=11
NB. inner m=18 n=19 i=5 j=4 u=11
NB. inner m=20 n=21 i=3 j=2 u=11

NB. outer m=24 n=25 i=1 j=0 u=12
NB. inner m=18 n=19 i=7 j=6 u=12
NB. inner m=20 n=21 i=5 j=4 u=12
NB. inner m=22 n=23 i=3 j=2 u=12

NB. outer m=26 n=27 i=1 j=0 u=13
NB. inner m=20 n=21 i=7 j=6 u=13
NB. inner m=24 n=25 i=3 j=2 u=13
NB. inner m=22 n=23 i=5 j=4 u=13

NB. outer m=28 n=29 i=1 j=0 u=14
NB. inner m=22 n=23 i=7 j=6 u=14
NB. inner m=24 n=25 i=5 j=4 u=14
NB. inner m=26 n=27 i=3 j=2 u=14

NB. outer m=30 n=31 i=1 j=0 u=15
NB. inner m=24 n=25 i=7 j=6 u=15
NB. inner m=26 n=27 i=5 j=4 u=15
NB. inner m=28 n=29 i=3 j=2 u=15







































13 : ' (1|. 2 # i.(#y)){y '
drot =: ] {~ 1 |. 2 # [: i. #
tmpy=.drot ydec
tmpw=.drot wavl
13 : '((#y)$ (i.(#x))){x'
[ {~ ([: # ]) $ [: i. [: # [
takx=: 2 # [: (] #~ 0 1 $~ #) [
NB. index acting on ydec;
 1|. 2 # i.32 NB. index for ydec 128
NB. hpf/lpf index for evens 
inlh=. 16$(3,2){hpf   NB. inner loop
inll=. 16$(3,2){lpf   NB. inner loop
NB. hpf/lpf index for odds
outlh =. 16$(1,0){hpf    NB.  outer loop
outll =. 16$(1,0){lpf    NB.  outer loop


outer=. ( (tmpy * outll) + tmpw * outlh) 
inner=.  (tmpy * inll) + tmpw * inlh


temp =. 128$ (1,2,3,4)
wav=.(1;'db4')dwtL temp
'db4'idwtLc  wav


