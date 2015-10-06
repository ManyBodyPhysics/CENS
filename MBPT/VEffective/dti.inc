c    dti.inc
c
c    Copyright (c) 1991-2000, Ilan Bar-On
c
c    You may distribute under the terms of either the GNU General Public
c    License or the Artistic License, as specified in the README file.
c
      real*4 dtime
      real*4 tt_time(2), sec  
      common /dtimers/ tt_time, sec

      real*4 nflops
      common /flops/ nflops

      character*70 header, nheader, kheader
      common /titles/ header, nheader, kheader

      

