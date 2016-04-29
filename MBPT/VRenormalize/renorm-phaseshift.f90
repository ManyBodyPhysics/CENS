!             Program block renorm-phaseshift.f90
!
!             Authors:  Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!             LAST UPGRADE : Added December 2011
!
!                        
!
!                    Begin setup of phaseshift code
!
SUBROUTINE  setup_vlowk
  USE partial_waves
  USE constants
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  IMPLICIT NONE
  INTEGER :: loop, number_orbits, max_conf

  !     reserve space in memory for mesh point and h.o. wave function  arrays
  ALLOCATE ( ra (n_rel), wra (n_rel));   ALLOCATE ( krel (n_rel), wkrel (n_rel))
  ALLOCATE ( rnlr (n_rel, 0:lmax, 0:nmax) )
  ALLOCATE(coulomb_relcom(0:nmax, 0:nmax, 0:lmax))
  !     set up Coulomb interaction in relative coordinates using Brody-Jacob-Moshinsky
  coulomb_relcom = 0.0_dp  
  !  The argonne v18 interaction contains the coulomb + electromagnetic corrections 
  !  explicitely
  IF ( type_of_pot /= 'argonnev18') THEN
     IF ( coulomb_included =='coulomb') CALL coulomb_integral
  ENDIF
  CALL make_configurations_relcm
  !     set up mesh points for relative coordinates
  CALL vlowk_mesh                   
  !     find all partial waves with given jmin and jmax
  CALL setup_channels          
  !
  !     The configurations defined by the relative quantum numbers nl are then
  !     used to diagonalize the deuteron in momentum space 
  !     and obtain a model space effective interaction using a similarity
  !     transformation. This is done by the function setup_com_vlowk.
  !     Note that the interaction does not depend on the CoM coordinates
  !     and thus the problem can be separated out. The final effective 
  !     interaction is stored for transformation to the lab system.
  ! 
  CALL setup_com_vlowk
  !     With the final interaction in a harmonic oscillator basis defined by
  !      < nlNLJT_ZS | Veff | n'l'N'LJT_ZS > we perform the final transformation
  !     to the lab system in the function final_vlow_labsystem
  DEALLOCATE(coulomb_relcom)
  DEALLOCATE ( ra, wra, krel, wkrel) ; DEALLOCATE ( rnlr) 

END SUBROUTINE setup_vlowk
!
!           Set up the interaction in the relative system using
!           a similarity transformation in momentum space. Quantum numbers are 
!           < nl JT_ZS | Veff | n'l'JT_ZS >
!
SUBROUTINE setup_com_vlowk
  USE single_particle_orbits
  USE configurations
  USE constants
  USE relcm_gmatrix
  USE partial_waves
  IMPLICIT NONE
  INTEGER ::  pq_confs, loop, mspace, bra, ket
  REAL(KIND = 8), ALLOCATABLE :: v2body(:,:)

  DO loop=1,no_channels   
     pq_confs = rel_conf%nconfs_rel(loop)
     WRITE(6,*) ' Total states and model states for channel: ', loop, pq_confs
     IF ( pq_confs <= 0 ) CYCLE
     ALLOCATE(v2body(pq_confs, pq_confs))
     v2body = 0.
     !     Setup the interaction
     CALL vlowk_channel(loop,v2body)
     DO ket = 1, pq_confs
        DO bra = 1, pq_confs
           v_com_rel(loop,bra,ket) = v2body(bra,ket)
        ENDDO
     ENDDO
     !  then we set up the final CoM and relative coordinate effective interaction
     !      free space
     DEALLOCATE(v2body)
  ENDDO

END SUBROUTINE setup_com_vlowk
!
!     Obtain the bare interaction in a harmonic oscillator 
!     basis plus the kinetic energy and the Coulomb part. It contains
!     also the CoM correction to the relative coordinates. The latter depends
!     on the mass number of the nucleus
! 
SUBROUTINE vlowk_channel(i,vint)
  USE wave_functions
  USE relcm_gmatrix
  USE partial_waves
  USE configurations
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i  ! loop variable over NN channels
  INTEGER :: ncoup, n,np, bra, ket, k1, k2, a, c
  INTEGER :: la, lb, jang, lim1, lim2, lim3, lim4
  REAL(KIND = 8) :: e_coulomb, vsum
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: heff(:,:), vzz(:,:)
  REAL(KIND = 8), ALLOCATABLE :: vkk(:,:)
  REAL(KIND = 8),  INTENT(INOUT):: &
       vint(rel_conf%nconfs_rel(i),rel_conf%nconfs_rel(i))

  ncoup=1
  IF ( orb_lrel_max(i) == jang_rel(i) ) ncoup = 1
  IF ( orb_lrel_max(i) /= jang_rel(i) ) ncoup = 2
  ALLOCATE(vkk (ncoup*n_rel,ncoup*n_rel))
  ALLOCATE(vzz (ncoup*n_rel,ncoup*n_rel))
  vzz = 0.0D0
  jang=jang_rel(i)
  !     setup the NN interaction
  vkk = 0.0D0
  CALL nocorepotential_interface(ncoup,vkk,i)
  vzz = vkk
  WRITE(6,'(12H Channel Nr:,I3,7H l_min:,I3,7H l_max:,I3,3H J:,I3,3H S:,I3,4H Tz:,I3)') &
       i, orb_lrel_min(i), orb_lrel_max(i), jang_rel(i), spin_rel(i), iso(i)
  ALLOCATE(heff (ncoup*n_k1,ncoup*n_k1))
  heff =0.0D0
  IF ((type_of_pot == 'OPEP').OR.(type_of_pot == 'Tensorinteraction')   &
       .OR.(type_of_pot == 'LSinteraction') ) THEN
     heff = vzz
  ELSEIF (type_of_renormv =='vbare') THEN
      heff = vzz
  ELSE
     ! get similarity transformed interaction
     CALL vlowk_mtx(ncoup,vzz,heff,i)
  ENDIF
  !     make now transformation to h.o. basis in the rel and cm system
  !     loop over all cm and rel coordinate configurations
  DO bra =1, rel_conf%nconfs_rel(i)
     k1=bra
     a=rel_conf%rel_ab(i,k1)
     n=relcm_sp_data%nrel(a)
     la=relcm_sp_data%lrel(a)
     IF ((n+n+la) > nlmax) CYCLE
     DO ket =1,  rel_conf%nconfs_rel(i) 
        k2=ket
        c=rel_conf%rel_ab(i,k2)
        np=relcm_sp_data%nrel(c)        
        lb=relcm_sp_data%lrel(c)
        !  No dependence of the bare interaction upon the CoM momenta
        !  The Hamiltonian is also diagonal in L and N 
        vsum = 0.
        IF ((np+np+lb) > nlmax) CYCLE
        IF ( ncoup == 1) THEN
           lim1=1; lim2=n_k1 ; lim3=1 ; lim4=n_k1
        ELSEIF ( ncoup == 2 ) THEN
           IF ( (la == lb).AND. ( jang > la) ) THEN
              lim1=1; lim2=n_k1 ; lim3=1 ; lim4=n_k1
           ELSEIF ( (la == lb).AND. ( jang < la) ) THEN
              lim1=1+n_k1; lim2=n_k1+n_k1 ; lim3=1+n_k1 ; lim4=n_k1+n_k1
           ELSEIF ( la >  lb ) THEN
              lim1=1+n_k1; lim2=n_k1+n_k1 ; lim3=1 ; lim4=n_k1
           ELSEIF ( la <  lb ) THEN
              lim1=1; lim2=n_k1 ; lim3=1+n_k1 ; lim4=n_k1+n_k1
           ENDIF
        ENDIF
        CALL vlowk_hosc(rnlr(:,lb,np), rnlr(:,la,n),REAL(heff(lim1:lim2,lim3:lim4)),bra,ket,vsum)
        !       Here we set up the Coulomb part in an oscillator basis. 
        e_coulomb = 0.
        IF ( ( la == lb ).AND.( iso(i) == -1 ))  THEN
           e_coulomb = coulomb_relcom(n,np,la) 
        ENDIF
        !  only potential energy : coulomb + V_NN
        vint(ket,bra)=vsum+e_coulomb
     ENDDO
  ENDDO
  DEALLOCATE(heff); DEALLOCATE(vkk); DEALLOCATE(vzz)

END SUBROUTINE vlowk_channel
!
!       Compute the T(G)-mtx for the actual channel 
!
SUBROUTINE vlowk_mtx(ncoup,vzz,heff,ichan)
  USE wave_functions
  USE partial_waves
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncoup, ichan
  COMPLEX*16, DIMENSION(ncoup*n_k,ncoup*n_k), INTENT(IN) :: vzz
  COMPLEX*16, DIMENSION(ncoup*n_k1,ncoup*n_k1), INTENT(INOUT) :: heff
  INTEGER :: i, j, ntot, i1, np, nq, j1,j2, tisoz 
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: ham, cvec, temp
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) ::  ceig
  COMPLEX*16 :: sum, input_energy

  ! dimension of vectors and matrices
  ntot = ncoup*n_k
  ALLOCATE( ham(ntot, ntot), cvec(ntot,ntot), ceig(ntot), temp(ntot,ntot)) 
  ! setup hamiltonian to be diagonalized
  ham = DCMPLX(0.0D0,0.0D0) 
  cvec = dcmplx(0.0d0, 0.0d0)
  ceig = dcmplx(0.0d0, 0.0d0)

  CALL complex_hamiltonian(ham,vzz, ntot, ichan)
  temp = ham
  CALL vlowkdiag_exact( temp, cvec, ceig, ntot )
  !
  ! construct renormalized nn-interaction via the Lee-Suzuki sim. transform
  !
  ! size of model space     : np = n_k1*ncoup
  ! size of complement space: nq = n_k2*ncoup
  np = n_k1*ncoup; nq = n_k2*ncoup
  CALL effective_int( np, nq, ntot, ncoup, cvec, ceig, heff )   
  !
  ! subtract diagonal elements to obtain Veff, 
  ! and compare exact with effective interactions: 
  !
  tisoz = iso(ichan) 
  DO i = 1, np
     i1 = i
     IF ( i > n_k1 ) i1 = i-n_k1
     heff(i,i) = heff(i,i) - ( krel(i1) * krel(i1) )/p_mass(tisoz)
     DO j = 1, np
        j1 = j
        IF ( j > n_k1 ) j1 = j-n_k1
        heff(i,j) = heff(i,j)/SQRT( wkrel(i1)*wkrel(j1) )/krel(i1)/krel(j1)
     ENDDO
  ENDDO
  DEALLOCATE(ham, temp, cvec );   DEALLOCATE(ceig)

END SUBROUTINE vlowk_mtx
!
! complex scaled hamiltonian in momentum representation 
!
SUBROUTINE complex_hamiltonian(h,vzz,ntot,ichan)
  USE wave_functions
  USE partial_waves
  USE constants
  IMPLICIT NONE
  REAL(KIND = 8) :: delta
  INTEGER, INTENT(IN) :: ntot, ichan
  COMPLEX*16, DIMENSION(ntot,ntot), INTENT(INOUT) :: h
  COMPLEX*16, DIMENSION(ntot,ntot), INTENT(IN) :: vzz
  INTEGER :: i, j, i1, i2, tisoz
  COMPLEX*16 :: wi

  tisoz = iso(ichan) 
  h = 0.
  DO i = 1, ntot 
     i1 = i
     IF ( i > n_k ) i1 = i-n_k
     h(i,i) = ( krel(i1) * krel(i1) )/p_mass(tisoz)  + & 
          (krel(i1)) * (krel(i1)) * wkrel(i1) * Vzz(i,i) 
     DO j = 1, ntot 
        i2 = j
        IF ( j > n_k ) i2 = j-n_k
        IF (i /= j ) THEN
           h(i,j) = SQRT( wkrel(i2) * wkrel(i1) ) *  krel(i1) * krel(i2) * Vzz(i,j)   
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE complex_hamiltonian
!
! calculate effective interaction for gamow shell model calcs
!
SUBROUTINE effective_int( np, nq, ntot, ncoup, cvec, ceig, heff )   
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER, INTENT(in) ::  np, nq, ntot, ncoup
  INTEGER :: npp, k1, i1,i,k2,nqq,j
  INTEGER ::  model(np), orb, orbit_no(ntot)
  COMPLEX*16, DIMENSION(ntot,ntot), INTENT(in) :: cvec
  COMPLEX*16, DIMENSION(ntot,ntot) :: cvec_temp
  COMPLEX*16, DIMENSION(ntot), INTENT(in) :: ceig
  COMPLEX*16 :: cvec_pp(np,np), cvec_qp(nq,np), cvec_model(nq,np) 
  COMPLEX*16 :: ceig_p(np), ceig_model(np)
  REAL(kind = 8),  DIMENSION(ntot) ::  temp
  REAL(kind = 8) :: cvec_max, overlap(ntot)
  COMPLEX*16 :: e_a, e_b
  COMPLEX*16, DIMENSION(np,np), INTENT(inout) :: heff
  REAL*8 :: a1,a2,b1,b2

  ! calculate P-space overlap of all eigenvectors 
  ! loop over all eigen vectors
  DO orb = 1, ntot
     overlap(orb) = 0.
     DO i = 1, n_k1
        IF ( ncoup == 1 )THEN
           overlap(orb) = overlap(orb) + ABS( cvec(i,orb) )**2
        ELSEIF( ncoup == 2 ) THEN
           overlap(orb) = overlap(orb) + ABS( cvec(i,orb) )**2 + ABS( cvec(n_k + i,orb) )**2
        END IF
     END DO
     orbit_no(orb) = orb
  END DO
  ! sort all overlaps and corresponding orbit numbers 
  CALL eigenvalue_sort(overlap,orbit_no, ntot)
  cvec_pp = 0.; cvec_qp = 0.
  ! Set up eigenvalues and eigenvectors for model space and excluded space
  IF ( ncoup == 1 ) THEN
     DO i=1, np
        ! loop over all model space coefficients of exact eigenvectors |k>
        ! 
        k1 = orbit_no(i)
        DO j = 1, n_k
           IF ( j <= n_k1 ) THEN
              cvec_pp(j,i) = cvec(j,k1)
              ceig_p(i) = ceig(k1)
           ELSEIF ( j > n_k1 ) THEN
              cvec_qp(j-n_k1,i) = cvec(j,k1)
           ENDIF
        ENDDO
     ENDDO
  END IF
  IF ( ncoup == 2 ) THEN
     DO i=1, np
        ! loop over all model space coefficients of exact eigenvectors |k>
        ! 
        k1 = orbit_no(i)
        ceig_p(i) = ceig(k1)
        DO j = 1, n_k1
           cvec_pp(j,i) = cvec(j,k1)
           cvec_pp(j+n_k1,i) = cvec(j+n_k,k1)
        ENDDO
     ENDDO
     DO i=1, np
        k1 = orbit_no(i)
        DO j = 1, n_k2
           cvec_qp(j,i) = cvec(n_k1+j,k1)
           cvec_qp(j+n_k2,i) = cvec(n_k+n_k1+j,k1)
        END DO
     END DO
  END IF
  heff = 0.
  CALL lee_suzuki2( cvec_pp, cvec_qp, ceig_p, np, nq, heff )
  !  subtract diagonal terms

END SUBROUTINE effective_int
!
! eigenvalue sort
! sort cmplx vector real(a(1))<real(a(2)) < ... < real(a(n))
!
SUBROUTINE eigenvalue_sort(A,B, n)
  IMPLICIT NONE
  INTEGER :: i, j, n
  REAL(kind = 8), DIMENSION(n), INTENT(INOUT) :: A
  INTEGER, DIMENSION(n), INTENT(inout) :: B
  REAL(kind = 8) :: temp, temp1
  REAL(kind = 8), DIMENSION(n) :: temp2
  INTEGER :: orb

  DO i = 1, n
     DO j = 1, n
        IF ( ABS( A(i) )  > ABS( A(j) ) ) THEN
           temp = A(i)
           A(i) = A(j) 
           A(j) = temp
           orb = B(i) 
           B(i) = B(j)
           B(j) = orb
        END IF
     END DO
  END DO

END SUBROUTINE eigenvalue_sort
!
! Lee Suzuki similarity transformation 
!
SUBROUTINE lee_suzuki2( cvec_pp, cvec_qp, ceig, np, nq, heff )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: np, nq
  COMPLEX*16, DIMENSION(np,np), INTENT(IN) :: cvec_pp
  COMPLEX*16, DIMENSION(nq,np), INTENT(IN)  :: cvec_qp
  COMPLEX*16, DIMENSION(np), INTENT(IN) :: ceig
  COMPLEX*16, DIMENSION(np,np) :: temp, omega2, sim1, sim2, omega_2inv, u, u_inv, & 
       cvec_pp_inv
  COMPLEX*16, DIMENSION(nq,np) :: omega
  COMPLEX*16, DIMENSION(np) :: ceig_p, omega2_eig, eigval2
  COMPLEX*16, DIMENSION(np,np) :: eigen_vec, vl, heff_rhs
  COMPLEX*16, DIMENSION(np,np), INTENT(INOUT) ::  heff
  REAL(KIND = 8), DIMENSION(2*np) :: rwork
  COMPLEX*16, DIMENSION(10000) :: work1
  COMPLEX*16 :: d, sum1, temp1(np,np), temp2(np,nq), norm, determinant
  INTEGER :: i_p,j_p, j_q, ii, jj, k, i_q , a_p, a_pp 
  INTEGER :: i, lda, ldb, ldvl, ldvr, info, lwork, ilo , ihi
  CHARACTER*1 :: jobvl, jobvr, balanc, sense
  REAL(KIND = 8), DIMENSION(np) :: scale, rconde, rcondv
  REAL(KIND = 8) :: abnrm
  INTEGER :: ipiv(np) ,j
  REAL(KIND = 8) :: a1, a2, b1, b2

  balanc = 'n';  jobvl = 'n' ;  jobvr = 'v';  sense = 'n';  lda = np
  ldvl = 1;  ldvr = np;  lwork = 10000
  eigen_vec = 0. 
  temp1 = TRANSPOSE(cvec_pp) 
  CALL zgeev( jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  determinant = PRODUCT(omega2_eig(:))
  !  write(6,*) 'check determinant', determinant
  ! the P->Q space transformation matrix, omega 
  cvec_pp_inv = cvec_pp
  CALL cmplxmatinv(cvec_pp_inv, np, d)
  DO i_p = 1, np
     DO i_q = 1, nq
        omega(i_q,i_p) = SUM( cvec_pp_inv(:,i_p)*cvec_qp(i_q,:) )
     ENDDO
  ENDDO
  ! non-hermitian effective interaction
  ! setup 2-p effective interaction in P-space
  heff = 0.
  DO i_p = 1, np
     DO j_p = 1, np 
        heff(i_p,j_p) = SUM( cvec_pp(i_p,:)*ceig(:)*cvec_pp(j_p,:) ) 
        sum1 = 0.
        DO k = 1, np
           DO i_q = 1, nq
              sum1 = sum1 + cvec_pp(i_p,k) * ceig(k) * cvec_qp(i_q,k) * omega(i_q,j_p) 
           ENDDO
        ENDDO
        heff(i_p,j_p) = heff(i_p,j_p) + sum1
     ENDDO
  ENDDO
  ! organizing the matrix (P(1 + omega * omega)P)
  omega2 = MATMUL(TRANSPOSE(cvec_pp_inv),cvec_pp_inv)
  ! calculate sqrt and inverse sqrt of matrix (P(1 + omega * omega)P)
  CALL sqrtmat( omega2, U, U_inv ,np ) 
  heff_rhs =  MATMUL( U, MATMUL( heff,u_inv ) ) 
  ! check if heff is symmetrized:
  heff = 0.
  heff = heff_rhs
  DO i_p = 1, np
     DO j_p = 1, np
        ! make heff manifestly symmetric
        IF ( i_p /= j_p ) THEN
           IF ( ABS(heff(i_p,j_p)- heff(j_p,i_p)) < 1.E-6 ) CYCLE 
           !           WRITE(6,*) 'sym test', heff(i_p,j_p), heff(j_p,i_p)
        ENDIF
     ENDDO
  ENDDO
  ! diagonalize 2p-effective shell model hamiltonian
  CALL zgeev( jobvl, jobvr, np, heff_rhs, lda, ceig_p, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  !   compare spectrum from exact and P-space diagonalization
  WRITE(6,*) 'Compare model space two-body spectrum with exact spectrum:' 
  DO i_p = 1, np
     a1 = REAL( ceig_p(i_p))
     a2 = AIMAG( ceig_p(i_p))
     b1 = REAL( ceig(i_p) )
     b2 = AIMAG( ceig(i_p) ) 
     WRITE(6,*) A1, B1
  ENDDO

END SUBROUTINE lee_suzuki2
!
! eigenvalue sort
! sort cmplx vector real(a(1))<real(a(2)) < ... < real(a(n))
!
SUBROUTINE eigenvalue_sort_cmplx(A, n)
  IMPLICIT NONE
  INTEGER :: i, j, n
  COMPLEX*16, DIMENSION(n), INTENT(INOUT) :: A
  COMPLEX*16 :: temp, temp1
  COMPLEX*16, DIMENSION(n) :: temp2

  DO i = 1, n
     DO j = 1, n
        IF ( ABS( A(i) )  > ABS( A(j) ) ) THEN
           temp = A(i)
           A(i) = A(j) 
           A(j) = temp
        END IF
     END DO
  END DO

END SUBROUTINE eigenvalue_sort_cmplx
!
!fast diagonalizing of complex symmetric matrices
!
SUBROUTINE vlowkdiag_exact( h, cvec, ceig, n )
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(in) :: h
  INTEGER :: k, np
  INTEGER :: p, i,j, i1, kvec, lwork, option, il ,iu, info
  REAL(kind = 8) :: A(n*(n+1)), work(300*n)
  REAL(kind = 8) :: thresh
  COMPLEX*16 :: cvl, cvu
  LOGICAL :: flag(n)  
  COMPLEX*16 :: ceig(n), sum1
  COMPLEX*16 :: cvec(n,n)
  REAL(kind = 8) ::  vl,vu

  lwork = 300*n
  thresh = 30.0
  kvec = n
  option = 4
  info = 1
  i1 = 0
  DO i =  1, n
     DO j =  1, i
        i1 = i1 + 1
        a(i1) = DBLE(h(j,i))
        a(i1+n*(n+1)/2) = AIMAG(h(j,i))
     END DO
  END DO
  CALL cs(n,a,ceig,kvec,cvec,lwork,work,thresh, option,il,iu,cvl,cvu,flag,info)

END SUBROUTINE vlowkdiag_exact
!
!          Obtain the bare potential in oscillator basis
!          Final V is in units of MeV or eV
!          Time consuming part. Parallel version available.
!
SUBROUTINE vlowk_hosc(wave_bra, wave_ket, a,bra,ket,vsum)
  USE wave_functions
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  REAL(KIND = 8), INTENT(INOUT) :: vsum
  INTEGER, INTENT(IN) :: bra, ket
  INTEGER :: i, j
  REAL(KIND = 8), DIMENSION(n_k1,n_k1), INTENT(IN) :: a
  REAL(KIND = 8) :: sum1, sum2, hbarc3
  REAL(KIND = 8), DIMENSION(n_rel), INTENT(IN) :: wave_bra, wave_ket 

  hbarc3=hbarc**3
  sum1=0.
  DO i=1, n_k1
     sum2=0.
     DO j=1, n_k1
        sum2=sum2+wave_ket(j)*a(j,i)
     ENDDO
     sum1=sum1+sum2*wave_bra(i)
  ENDDO
  vsum = sum1*hbarc**3

END SUBROUTINE vlowk_hosc
