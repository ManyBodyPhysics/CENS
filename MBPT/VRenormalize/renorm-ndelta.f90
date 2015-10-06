!             Program block renorm-ndelta.f90
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!             LAST UPGRADE : April 2009, 
!             Sets up a nucleon-delta interaction with pions and rho mesons
!
!
SUBROUTINE  setup_ndelta
  USE partial_waves
  USE constants
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  IMPLICIT NONE

END SUBROUTINE setup_ndelta





!	 **************************************************
!	 *    Here we set up the NNND, or the  NNDD       *
!        *    transition potentials and g-matrices        *
!        *    The isospin-spin algebra is set up          *
!        *    by the function SIFAC, whereas the masses,  *
!        *    form factors and cutoffs are set up in      *
!        *    COEFFS. Similarly, the coupling constants   *
!        *    are set up in COUPLFA. Observe that in the  *
!        *    the coupling constants differ in NNND and   *
!        *    the NNDD cases. For the ND-mes              *
!        *    coupling constants we have used the quark   *
!        *    model (non-rel) whereas in the NN case we   *
!        *    follow Machleidt in Adv. Nucl. Phys. vol 19 *
!        *    The function FACTOR sets up factorials to   *
!        *    be used in the calculation of Legendre      *
!        *    functions of the second kind                *
!	 **************************************************

      SUBROUTINE vnnnd(nmesh,ifree,ndelta)
      IMPLICIT NONE
      COMMON/gnn/grm(25,6,25,6,2,2,5),q(25),u(25)

      corr = 1.0d0/(hbarc**3)              !  potential in units of mev**-2
      CALL factor


!                 ****************************
!                 *   set up pot. for nnnd   *
!                 *   and search for nn-chan *
!                 ****************************


        it = 1

        CALL couplfa(1)
        CALL sifac(1, 1, 1, 2)
        DO j = 0, 4
           DO is = 0, 1
	      CALL jtest(j, l1, l2, tester, is, it)
              IF(MOD(l1+is+it,2).EQ.0) CYCLE
	      IF(tester) THEN
	        CALL g_nnnd(j,is,l1,l2,nmesh,ifree)
	      ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE vnnd



      SUBROUTINE vnnnd(j,is,l1,l2,nmesh,ifree)
      IMPLICIT REAL*8(a-h,o-z)
      COMMON/gnn/grm(25,6,25,6,2,2,5),q(25),u(25)
      COMMON/nnnd/gnnnd(25,6,25,6,2,2,5)

      DIMENSION vnnnd(25, 25, 8)
      REAL*8 hbarom/197.327d0/             ! mev*fermi
      INTEGER itabnn(3,8)
     :   /0,0,1,0,1,1,1,1,3,2,2,3,2,1,8,3,3,3,4,4,2,4,3,2/
      INTEGER itabnd(2, 8, 8)
      itabnd(1,1,1) = 2
      itabnd(2,1,1) = 2
      itabnd(1,1,2) = 1
      itabnd(2,1,2) = 1
      itabnd(1,1,3) = 1
      itabnd(2,1,3) = 1
      itabnd(1,2,3) = 2
      itabnd(2,2,3) = 1
      itabnd(1,3,3) = 2
      itabnd(2,3,3) = 3
      itabnd(1,1,4) = 2
      itabnd(2,1,4) = 0
      itabnd(1,2,4) = 2
      itabnd(2,2,4) = 2
      itabnd(1,3,4) = 2
      itabnd(2,3,4) = 4
      itabnd(1,1,5) = 1
      itabnd(2,1,5) = 1
      itabnd(1,2,5) = 2
      itabnd(2,2,5) = 1
      itabnd(1,3,5) = 1
      itabnd(2,3,5) = 3
      itabnd(1,4,5) = 2
      itabnd(2,4,5) = 3
      itabnd(1,5,5) = 1
      itabnd(2,5,5) = 1
      itabnd(1,6,5) = 2
      itabnd(2,6,5) = 1
      itabnd(1,7,5) = 1
      itabnd(2,7,5) = 3
      itabnd(1,8,5) = 2
      itabnd(2,8,5) = 3
      itabnd(1,1,6) = 1
      itabnd(2,1,6) = 3
      itabnd(1,2,6) = 2
      itabnd(2,2,6) = 1
      itabnd(1,3,6) = 2
      itabnd(2,3,6) = 3
      itabnd(1,1,7) = 2
      itabnd(2,1,7) = 4
      itabnd(1,2,7) = 2
      itabnd(2,2,7) = 2
      itabnd(1,1,8) = 1
      itabnd(2,1,8) = 3
      itabnd(1,2,8) = 2
      itabnd(2,2,8) = 3

      corr  = 1.0d0/(hbarom**3)              ! units in (mev-(2))*fm

!            *********   format declarations     *********

 200  FORMAT(//'no nn - ndelta transition potential for jtot, spin,l',
     :         3i5///)
 201  FORMAT(//'nnnd free transition potential for jtot =',i3,2x,
     : 'lnn =',i3,2x,'snn =',i3,2x,'lnd =',i3,2x,'snd =',i3,2x,d12.6)

!            *********************************************
!            *    search for possible nnnd channels      *
!            *********************************************

      llnn = l1                 !  l1 <= l2  (always!!!!)
      DO i =1, 8
         IF(j /= itabnn(1,i)) CYCLE
         IF(llnn /= itabnn(2,i)) CYCLE
         icas = i
         GOTO 20
      ENDDO
      WRITE(6, 200) j, is, llnn
      RETURN

 20   CONTINUE
      nchan = itabnn(3, icas)

!            *****************************************
!            *  setup of free transition potential   *
!            *****************************************      

      DO i = 1, nchan
         IF(i > 4) llnn = l2                ! coupled channel 3p2-3f2
         isnd = itabnd(1, i, icas)
         llnd = itabnd(2, i, icas)
         DO k1 = 1, nmesh
            ak1 = q(k1)/hbarom
            DO k2 = 1, nmesh
               ak2 = q(k2)/hbarom
               sumcoe=0.d0
               DO icoeff = 1, 2
                  IF (icoeff == 1) CALL coeffs(1)
                  IF(icoeff == 2)  CALL coeffs(2)
                  sumcoe=sumcoe+0.5d0*corr*(   &
                      vrho(ak1,llnn,is,ak2,llnd,isnd,j,1) &
                     + vpion(ak1,llnn,is,ak2,llnd,isnd,j,1))
               ENDDO
               vnnnd(k2,k1,i)=sumcoe
             ENDDO
         ENDDO
      ENDDO

      DO i = 1, nchan
         isnd = itabnd(1, i, icas)
         llnd = itabnd(2, i, icas)
	 llnn = l1
	 indx = 4
	 ll1 = l2
	 IF(i > 4) THEN
	   llnn = l2
	   indx =-4
	   ll1 = l1
	 ENDIF
	 DO k1 = 1, nmesh
	    DO k2 = 1, nmesh
		 gnnnd(k2,llnd+1,k1,llnn+1,isnd,is+1,j+1)=
     :                       vnnnd(k2,k1,i)
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE vnnnd




!             
!                 This function calculates the non-relativistic    
!                 transition potential < bb  v bb > for pions      
!                 coupled to  nn or nd . contributions from tensor 
!                 and central force only. the necessary spin and   
!                 isospin has been set up in routine sifac.        
!                 vt is the tensor part, vc the central part       
!                 arising from the spin-spin term.                 
!             

      REAL(DP) FUNCTION vpion(qa, la, isa, qb, lb, isb, j ,it)
      USE constants
      USE IsobarNucleonVariables
      IMPLICIT NONE
!      COMMON/coupls/gpion,grho
!      COMMON/spinrm/ smo(4), srm2(4,4)
!      COMMON/pion/fmpion,pcut1,pcut2,pform1,pform2
      vpion=0.0d0
      it1 = it + 1
      ia = isa + 1
      ib = isb + 1
      al = la
      bl = lb
      sa = isa
      sb = isb
      tj = j
      vt = 0.0d0
      vc = 0.0d0
      ff = smo(it1)*srm2(ia, ib)
      IF(ABS(ff)-1.0d-8 > 0.d0) THEN
	ip = isa + j + lb
	p = 1 + 2*(ip/2*2 - ip)
	f = ff*p*sixj(al, sa, tj, sb, bl, 2.0d0)
	vt = -f*3.0d0*(potten(la,qa,lb,qb,fmpion)  &
            -pform2*potten(la,qa,lb,qb,pcut1)      &
            +pform1*potten(la,qa,lb,qb,pcut2))
      ENDIF
      IF(((isa-isb) == 0).AND.((la-lb) == 0)) THEN
	f = smo(ia)*smo(it1)
	IF(ABS(f)-1.0d-8 > 0d0) THEN
	  vc = f*(potcen(la,qa,qb,fmpion)*(fmpion**2)
              -pform2*potcen(la,qa,qb,pcut1)*(pcut1**2) &
              +pform1*potcen(la,qa,qb,pcut2)*(pcut2**2))  &
	ENDIF
      ENDIF
      vpion = (vt + vc)*gpion

      END FUNCTION vpion


!             *******************************************************
!             *    this function calculates the non-relativistic    *
!             *    transition potential < bb  v bb > for rho mesons *
!             *    coupled to  nn or nd .                           *
!             *******************************************************

      REAL(DP) FUNCTION vrho(qa, la, isa, qb, lb, isb, j ,it)
      IMPLICIT NONE
      COMMON/coupls/gpion,grho
      COMMON/spinrm/ smo(4), srm2(4,4)
      COMMON/rho/fmrho,rcut1,rcut2,rform1,rform2
      vrho=0.0d0
      it1 = it + 1
      ia = isa + 1
      ib = isb + 1
      al = la
      bl = lb
      sa = isa
      sb = isb
      tj = j
      vt = 0.0d0
      vc1 = 0.0d0
      ff = smo(it1)*srm2(ia, ib)
      IF(ABS(ff)-1.0d-8 > 0.d0) THEN
	ip = isa + j + lb
	p = 1 + 2*(ip/2*2 - ip)
	f = ff*p*sixj(al, sa, tj, sb, bl, 2.0d0)/4.d0
	vt= f*(potten(la,qa,lb,qb,fmrho) &
            -rform2*potten(la,qa,lb,qb,rcut1) &
            +rform1*potten(la,qa,lb,qb,rcut2))
      ENDIF
      IF(((isa-isb) == 0).AND.((la-lb) == 0)) THEN
	f = smo(ia)*smo(it1)/6.d0
	IF(ABS(f)-1.0d-8 > 0d0) THEN
	  vc1=f*(potcen(la,qa,qb,fmrho)*(fmrho**2) &
             -rform2*potcen(la,qa,qb,rcut1)*(rcut1**2) &
             +rform1*potcen(la,qa,qb,rcut2)*(rcut2**2))
	ENDIF
      ENDIF
      vrho = (vt + vc1 )*grho

      END FUNCTION vrho



!       *************************************************************
!       *   this function sets up the spin-isospin algebra in a non *
!       *   relativistic calculation. ni = 1 for nucleons and 2 for *
!       *   deltas.                                                 *
!       *************************************************************

      SUBROUTINE sifac(n1,n2,n3,n4)
      IMPLICIT NONE
      COMMON/spinrm/smo(4), srm2(4,4)
      DIMENSION rs(2,2)
      DATA rs/2.449489743, 2.0d0, -2.0d0, 7.745966692/
      s1 = n1 - 0.50d0
      s2 = n2 - 0.50d0
      s3 = n3 - 0.50d0
      s4 = n4 - 0.50d0
      DO is = 1, 4
         s = is - 1
         ip = n2 + n3 + is
         p = 1 + 2*(ip/2*2 - ip)
         smo(is) = p*sixj(s1,s2,s,s4,s3,1.0d0)*rs(n1,n3)*rs(n2,n4)
      ENDDO
      DO isa = 1, 4
         sa = isa - 1
         DO isb = 1, 4
            srm2(isa,isb) = 0.0d0
            IF(isa+isb-4 < 0) CYCLE
            IF(ABS(isa-isb)-2 > 0) CYCLE
              sb = isb - 1
              srm2(isa, isb)=SQRT((sa+sa+1.0d0)*5.0d0*(sb+sb+1.0d0))&
                 *anjs(s1,s2,sa,s3,s4,sb,1.0d0,1.0d0,2.0d0)*rs(n1,n3)* &
                    rs(n2,n4)
         ENDDO
      ENDDO

      END SUBROUTINE sifac
!              *****    end of function  sifac    *****              

!       *************************************************************
!       *   this function sets up the necessary form factors        *
!       *   cutoffs and masses for the various mesons.              *
!       *   the definition of the form factors are given in eq.     *
!       *   a.28 of r. machleidt, adv. nucl. phys. 19               *
!       *   the masses are given in table a.3 of the same work.     *
!       *   legend to the variables:                                *
!       *   fm**  = mass for the various mesons divided by hbar*c   *
!       *   *cut1 = cutoff + epsilon(10 mev)   -----  "  -------    *
!       *   *cut2 = cutoff - epsilon                                *
!       *   *form1 = (cut1**2 - mesmass**2)/(cut2**2 - cut1**2)     *
!       *   if ntype = 1, then the meson mass is given by           *
!       *                                                           *      
!       *                  m_(meson)                                *
!       *   elseif ntype = 2, we have an effective meson mass       *
!       *                                                           *
!       *      sqrt( m_(meson)*(m_(delta)-m_(nucleon)+m_(meson)))   *     
!       *                                                           *
!       *************************************************************

      SUBROUTINE coeffs(ntype)
      IMPLICIT NONE
      COMMON/pion/fmpion,pcut1,pcut2,pform1,pform2
      COMMON/rho/fmrho,rcut1,rcut2,rform1,rform2

!              ********************
!              *      masses      *
!              ********************

      IF(ntype == 1) THEN
         fmpion = 0.6995
         fmrho = 3.897
      ELSEIF(ntype == 2) THEN
         fmpion = 1.236
         fmrho =  4.5799
      ENDIF

!              ********************
!              *     cutoffs      *
!              ********************

      pcut1 = 6.03d0                        ! 1.2 gev
      pcut2 = 6.13d0 
      rcut1 = 6.54d0                         ! 1.3 gev
      rcut2 = 6.64d0

!              ********************
!              *  form factors    *
!              ********************

      pform1 = (pcut1**2 - fmpion**2)/(pcut2**2 - pcut1**2)
      pform2 = (pcut2**2 - fmpion**2)/(pcut2**2 - pcut1**2)
      rform1 = (rcut1**2 - fmrho**2)/(rcut2**2 - rcut1**2)
      rform2 = (rcut2**2 - fmrho**2)/(rcut2**2 - rcut1**2)

      END SUBROUTINE coeffs



!      *************************************************************
!       *   this function sets up the necessary coupling constants  *
!       *   for the various mesons.                                 *
!       *   legend to coupling constants  g**2(mes)/4*pi            *
!       *                  -----------------------------            *
!       *                       nn             nd                   *
!       *                  -----------------------------            *
!       *           g_pi       14.6        72/25*14.6               *
!       *           g_rho      0.95        72/25*g_rho*7.1**2       *
!       *************************************************************

      SUBROUTINE couplfa(itype)
      IMPLICIT NONE
      COMMON/coupls/gpion,grho
      REAL(DP) hbarc/197.327d0/           ! particle-data group 1988
      REAL(DP) quarkm/2.88d0/
      REAL(DP) gvgt/50.41d0/
      tmass = (938.9/hbarc)**2          ! particle-data group 1988
      IF(itype == 1) THEN
	gpion = 14.6*sqrt(quarkm)*hbarc/12.d0/tmass
	grho =  0.95*sqrt(quarkm)*gvgt*hbarc/tmass
      ELSEIF(itype == 2) THEN
        gpion = 14.6*quarkm*hbarc/12.0d0/tmass
        grho  = 0.95*quarkm*gvgt*hbarc/tmass
      ENDIF

      END SUBROUTINE couplfa



!       *************************************************************
!       *   this function calculates the central contribution       *
!       *************************************************************

      REAL(DP) FUNCTION potcen(la,qa,qb,ac)
      IMPLICIT NONE
      DATA f/0.318309886d0/          ! 1/pi
      x = (qa*qa+qb*qb+ac*ac)/(2.0d0*qa*qb)
      potcen = f*qlx(la,x)/(qa*qb)

      END FUNCTION potcen



!       *************************************************************
!       *   this function calculates the tensor contribution       *
!       *************************************************************

      REAL(DP) FUNCTION potten(la,qa,lb,qb,ac)
      IMPLICIT NONE
      DATA f/0.8164965809/          ! square root of 2/3
      DATA sr5/2.236067977/         ! square root of 5
      DATA pi/0.318309886d0/        ! 1/pi
      potten = 0.0d0
      a = 0.0d0
      b = 0.0d0
      x = (qa*qa+qb*qb+ac*ac)/(2.0d0*qa*qb)
      al = la
      bl = lb
      IF(la-lb) 20, 10, 20
 20   IF(la-lb-2) 21,11,21
 21   IF(la-lb+2) 500, 12, 500
 10   IF(la) 500, 500, 13
 13   	ma = la -1
	mb = la + 1
	am = ma
	bm = mb
	b=2.0d0*sr5*((am+am+1.0d0)*
     :       dreij(1.0d0,am,al,0.0d0,0.0d0,0.0d0)**2
     :      *sixj(al,bl,2.0d0,1.0d0,1.0d0,am)*qlx(ma,x)+
     :       (bm+bm+1.d0)*dreij(1.0d0,bm,al,0.0d0,0.0d0,0.0d0)**2
     :       *sixj(al,bl,2.0d0,1.0d0,1.0d0,bm)*qlx(mb,x))*qa*qb
        GOTO 100
 11     ma = la -1
        GOTO 15
 12     ma = la +1
 15     am = ma
	b = 2.d0*sr5*qa*qb*(am+am+1.d0)
     :      *dreij(1.0d0,am,al,0.0d0,0.0d0,0.0d0)
     :      *dreij(1.0d0,am,bl,0.0d0,0.0d0,0.0d0)*qlx(ma,x)
     :      *sixj(al,bl,2.0d0,1.0d0,1.0d0,am)
 100  pa = 1 + 2*(la/2*2-la)
      pb = 1 + 2*(lb/2*2-lb)
      a = f*(qa*qa*pb*qlx(lb,x)+qb*qb*pa*qlx(la,x))
     :      *dreij(2.0d0,al,bl,0.0d0,0.0d0,0.0d0)
      potten = pi*(a+b)*dsqrt((al+al+1.d0)*(bl+bl+1.d0))/(qa*qb)
 500  RETURN
      END
!              ******  end of function peten   ******



!             *******************************************************
!             *    legendre function of the second kind. see        *
!             *    haftel and tabakin nucl. phys. a145 (1970) p. 1  *
!             *    for further references                           *
!             ********************************************************

      REAL(DP) FUNCTION qlx(l,x)
      IMPLICIT REAL*8(a-h,o-z)

      qlx = 0.0d0
      IF(x-100.0d0 >= 0.d0) THEN           !  for large x values.
	y = 0.0d0
	DO n = 1, 5
	   yd = fact(n)*dfact(l+n)*((2.0d0*x*x)**(n-1))
	   y = y + fact(l+n+n-1)/yd
        ENDDO
	qlx = y/(x**(l+1))
      ELSEIF(x-100.d0 < 1.0d0) THEN
	a = 0.50d0*LOG((x+1.0d0)/(x-1.0d0))
	IF(l < 0) THEN
	  qlx = 0.0d0
	ELSEIF(l == 0) THEN
	  qlx = a
	ELSEIF(l > 0) THEN
	  b = a
	  a = b*x-1.0d0
	  IF(l-1 < 0) THEN
	    qlx = 0.0d0
	  ELSEIF(l-1 == 0) THEN
	    qlx = a
	  ELSEIF(l-1 > 0) THEN
	    DO k =2, l
	       c = b
	       b = a
	       gk = 1.0d0/dfloat(k)
	       a = (2.0d0-gk)*x*b - (1.0d0 - gk)*c
            ENDDO
	    qlx = a
	  ENDIF
	ENDIF
      ENDIF

      END function qlx



!         **********************************************************
!         *   this function finds the possible nn channels which   *
!         *   relate to the nd coupling if isospin is equal to 1   *
!         *   or equal to 0.                                       *
!         *   it returns the values l1  and l2.                    *
!         *   if tester =.false. no nnnd potential will be calc.   *
!         *   if l1 = l2 no coupled nn channels.                   *
!         *   note that   l1 <= l2                                 *
!         **********************************************************

      subroutine jtest(j, l1, l2, tester, is, it)
      implicit none
      logical tester
      tester =.false.
      if(it.eq.1) then
        if((j.eq.0).and.(is.eq.0))then
           l1 = 0
           l2 = 0
           tester = .true.
        elseif((j.eq.0).and.(is.eq.1))then
           l1 = 1
           l2 = l1
           tester = .true.
        elseif((j.eq.1).and.(is.eq.1))then
           l1 = 1
           l2 = l1
           tester = .true.
        elseif((j.eq.2).and.(is.eq.0))then
           l1 = 2
           l2 = l1
           tester =.true.
        elseif((j.eq.2).and.(is.eq.1))then
           l1 = 1
           l2 = 3
           tester = .true.
        elseif((j.eq.3).and.(is.eq.1))then
           l1 = 3
           l2 = l1
           tester = .true.
        elseif((j.eq.4).and.(is.eq.0))then
           l1 = 4
           l2 = l1
           tester = .true.
        elseif((j.eq.4).and.(is.eq.1))then
           l1 = 3
           l2 = l1 
           tester = .true.
        else
           tester = .false.
        endif
      elseif(it.eq.0) then
        if((j.eq.1).and.(is.eq.0)) then
           l1 = 1
           l2 = l1
           tester = .true.
        elseif((j.eq.1).and.(is.eq.1)) then
           l1 = 0
           l2 = 2
           tester = .true.
        elseif((j.eq.2).and.(is.eq.1)) then
           l1 = 2
           l2 = l1
           tester = .true.
        elseif((j.eq.3).and.(is.eq.0)) then
           l1 = 3
           l2 = l1
           tester = .true.
        elseif((j.eq.3).and.(is.eq.1)) then
           l1 = 2
           l2 = 4
           tester = .true.
        elseif((j.eq.4).and.(is.eq.1)) then
           l1 = 4
           l2 = l1
           tester = .true.
        else
           tester = .false.
        endif
      endif
      return
      end subroutine jtest

