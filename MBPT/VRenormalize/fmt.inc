c    fmt.inc
c
c    Copyright (c) 1991-2000, Ilan Bar-On
c
c    You may distribute under the terms of either the GNU General Public
c    License or the Artistic License, as specified in the README file.
c
c	Formats
	
100   format (4g19.8)
200   format (1x,'(',I4,')',a20,g12.2)
300   format (A,3g10.2)
399   format (a20, 2x, I4)
400   format (2x,'Max |x|=',g10.1,4x,'|lap-new|=',g10.1,4x,'err=',g10.1)
