c    drv.inc
c
c    Copyright (c) 1991-2000, Ilan Bar-On
c
c    You may distribute under the terms of either the GNU General Public
c    License or the Artistic License, as specified in the README file.
c
c	Compare Lapack
	logical cmp_lap_trd, cmp_lap_dns
	common /lap/ cmp_lap_trd, cmp_lap_dns



c	Check Condition of Eigenvectors 
	logical chk_trd_cnd, chk_dns_cnd
	common /cnd/ chk_trd_cnd, chk_dns_cnd


c	Normalize Eigenvectors
	logical set_normal_on
	common /normal/ set_normal_on

	logical whole_set
	common /whole/ whole_set 

c	Norms
      	real*8 a_norm,t_norm 
       	common /mat_norm/ a_norm,t_norm


c	Growth
      	REAL*8 max_cpx
      	integer cpx_ind
      	common /GROWTH_CPX/ max_cpx, cpx_ind

      	REAL*8 min_rho
      	integer rho_ind
      	common /GROWTH_RHO/ min_rho,rho_ind


c	Fail
      	integer no_fail
      	common /FAIL/ no_fail

      






