c    def.inc
c
c    Copyright (c) 1991-2000, Ilan Bar-On
c
c    You may distribute under the terms of either the GNU General Public
c    License or the Artistic License, as specified in the README file.
c
c     	Epsilons
      	real*8 mach_eps
      	parameter(mach_eps =2.0d-16) 

c     	Constants
      	real*8 ZERO,ONE,TWO
      	parameter ( ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0 ) 

      	complex*16 CZERO,CONE,CIMG
      	parameter ( 
     $ CZERO=(0.0D0,0.0D0),
     $ CONE =(1.0D0,0.0D0),
     $ CIMG =(0.0D0,1.0D0) 
     $       )


 

