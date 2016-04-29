!             Program block renorm-vkrg.f90
!
!             Authors:  Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!             LAST UPGRADE : May 2009, p and n masses different
!
!                        
!
!                    Begin setup of renormalization group based interaction
!                    Here we have coded Eq. (6) of Bogner et al,
!                    arXiv:nucl-th0611045.   It is a simple set of
!                    coupled, first-order differential equations.
!                    This is the momentum space version
!
SUBROUTINE  setup_vkrg
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
  !     set up Coulomb interaction in relative coordinates 
  coulomb_relcom = 0.0_dp  
  !  The argonne v18 interaction contains the coulomb + electromagnetic corrections 
  !  explicitely
  IF ( type_of_pot /= 'argonnev18') THEN
     IF ( coulomb_included =='coulomb') CALL coulomb_integral
  ENDIF
  !     set up quantum numbers nl and NL relative and cm system
  CALL find_spdata_relcm(number_orbits)
  relcm_sp_data%max_value=number_orbits
  CALL allocate_relcm_array(relcm_sp_data, number_orbits) 
  CALL make_configurations_relcm
  !     set up mesh points for relative coordinates
  CALL vkrgk_mesh                   
  !     setup h.o. wave functions, only relative coordinates needed
  CALL nocoreho_wfunction                
  !     find all partial waves with given jmin and jmax
  CALL setup_channels          
  !
  !     Note: first we find all possible configurations in the relative system
  !     only based on h.o. quantum number nl for all partial waves
  !     and allocate array for the given partial wave
  !     all configs are set up in the call to setup_configurations_rel
  !     First we search for the number of configurations in each channel
  !
  ALLOCATE ( rel_conf%nconfs_rel(no_channels) ) 
  ALLOCATE ( rel_conf%nconfs_relmodel(no_channels) ) 
  CALL number_configurations_rel(rel_conf)
  max_conf=0
  DO loop = 1, no_channels
     IF ( max_conf <  rel_conf%nconfs_rel(loop) ) THEN
        max_conf= rel_conf%nconfs_rel(loop)
     ENDIF
  ENDDO
  ALLOCATE( rel_conf%rel_ab(no_channels, max_conf))
  ALLOCATE ( v_com_rel(no_channels,max_conf,max_conf))
  v_com_rel = 0. 
  CALL setup_configurations_rel(rel_conf)      
  ! 
  !     Then we find all possible configurations in the relative and cm system
  !     based on h.o. quantum number nlNL for all partial waves
  !     and allocate array for the given partial wave
  !     all configs are set up in the call to setup_configurations_cm
  !     First we search for the number of configurations in each channel
  !
  ALLOCATE ( relcm_conf%nconfs_relcm(no_channels) ) 
  CALL number_configurations_relcm(relcm_conf)
  max_conf=0
  DO loop = 1, no_channels
     IF ( max_conf <  relcm_conf%nconfs_relcm(loop) ) THEN
        max_conf= relcm_conf%nconfs_relcm(loop)
     ENDIF
  ENDDO
  ALLOCATE( relcm_conf%relcm_ab(no_channels, max_conf+max_conf))
  CALL setup_configurations_relcm(relcm_conf)      
  !
  !     The configurations defined by the relative quantum numbers nl are then
  !     used to diagonalize the deuteron in momentum space 
  !     and obtain a model space effective interaction using a similarity
  !     transformation. This is done by the function setup_com_vkrgk.
  !     Note that the interaction does not depend on the CoM coordinates
  !     and thus the problem can be separated out. The final effective 
  !     interaction is stored for transformation to the lab system.
  !
  CALL setup_com_vkrgk
  !     With the final interaction in a harmonic oscillator basis defined by
  !      < nlNLJT_ZS | Veff | n'l'N'LJT_ZS > we perform the final transformation
  !     to the lab system in the function final_vkrg_labsystem
  CALL final_vkrg_labsystem
  DEALLOCATE(coulomb_relcom)
  DEALLOCATE ( ra, wra, krel, wkrel) ; DEALLOCATE ( rnlr) 
  DEALLOCATE (rel_conf%rel_ab) 
  DEALLOCATE (rel_conf%nconfs_rel ) ;   DEALLOCATE (rel_conf%nconfs_relmodel ) 
  CALL deallocate_relcm_array(relcm_sp_data)
  DEALLOCATE (relcm_conf%relcm_ab) 
  DEALLOCATE (relcm_conf%nconfs_relcm ) 

END SUBROUTINE setup_vkrg
!
!           Set up the effective interaction in the lab-frame using the
!           mtx elements in the rel-CoM frame
!
SUBROUTINE final_vkrg_labsystem
  USE single_particle_orbits  
  USE configurations
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs 
  REAL(KIND = 8), ALLOCATABLE :: lab_to_relcoeff(:,:)
  INTEGER, ALLOCATABLE :: lab_to_relconf(:,:)
  INTEGER, ALLOCATABLE :: lab_to_relnumber(:)
  INTEGER ::  p_parity, ang_mom, isospin_z, max_coeff, pq_confs
  REAL(KIND=8), ALLOCATABLE::  gna(:,:)
  REAL(KIND=8), ALLOCATABLE:: twobody_com(:,:), twobody_r2(:,:), twobody_p2(:,:)

  !     loop over isospin projection
  DO isospin_z=itzmin,itzmax 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=j_lab_min,j_lab_max
           !     find all possible configurations, large space and model space 
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,gmatrix_configs)
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
           pq_confs=gmatrix_configs%number_confs
           ALLOCATE(gmatrix_configs%config_ab(pq_confs+pq_confs) )
           CALL setup_gmatrix_configurations &
                (ang_mom,p_parity,isospin_z,gmatrix_configs)
           !     find max possible number of transformation coeffs
           CALL find_ncoeffs(isospin_z,ang_mom,gmatrix_configs,max_coeff)           
           !     allocate space for various arrays needed in transform rel-cm -> lab 
           ALLOCATE(lab_to_relconf(pq_confs,max_coeff), &
                lab_to_relcoeff(pq_confs,max_coeff))
           ALLOCATE(lab_to_relnumber(pq_confs))
           ALLOCATE(twobody_p2(pq_confs, pq_confs))
           ALLOCATE(twobody_r2(pq_confs, pq_confs))
           ALLOCATE(twobody_com(pq_confs, pq_confs))
           !     setup transformation coefficients for oscillator basis
           !     transformations from the c.m. frame to the lab frame
           CALL mosh_transf(isospin_z,ang_mom,gmatrix_configs, &
                lab_to_relcoeff, lab_to_relconf, lab_to_relnumber, &
                max_coeff)
           ALLOCATE(gna(pq_confs, pq_confs))
           gna=0.0D0
           twobody_com = 0.D0; twobody_p2 = 0.0D0; twobody_r2 = 0.0D0
           !     Performs HO transformation from rel and cm coordinates to
           !     lab system.
           CALL vkrgk_free(gna,twobody_com,twobody_r2, twobody_p2,lab_to_relcoeff, lab_to_relconf,&
                lab_to_relnumber,max_coeff,gmatrix_configs,isospin_z)
           CALL vkrgk_print(isospin_z,p_parity,ang_mom,gna,twobody_com,twobody_r2,twobody_p2,gmatrix_configs)
           !      free space
           DEALLOCATE(gmatrix_configs%config_ab)
           DEALLOCATE(twobody_com)
           DEALLOCATE(twobody_r2,twobody_p2)
           DEALLOCATE(lab_to_relcoeff, lab_to_relconf )
           DEALLOCATE (lab_to_relnumber); DEALLOCATE(gna)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE final_vkrg_labsystem
!
!                    Print out veff-mtx
!
SUBROUTINE vkrgk_print(it,ip,ij,gna,twobody_com,twobody_r2, twobody_p2,gmatrix_configs)
  USE constants
  USE configurations
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, ijd, it, ip, ij, ia, ib, ic,id
  REAL(KIND=8), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(IN)  :: gna
  REAL(KIND=8), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(IN) :: twobody_com,twobody_r2, twobody_p2
  ijd=ij+ij
  DO i=1,gmatrix_configs%number_confs
     ia= gmatrix_configs%config_ab(i*2-1)
     ib= gmatrix_configs%config_ab(i*2)
     DO j=i,gmatrix_configs%number_confs
        ic= gmatrix_configs%config_ab(j+j-1)
        id= gmatrix_configs%config_ab(j+j)
        WRITE(7,1001) it, ip, ijd, ia, ib, ic, id, REAL(gna(j,i)), &
             twobody_com(j,i), twobody_r2(j,i), twobody_p2(j,i)       
     ENDDO
  ENDDO
1001 FORMAT(7I4,5X,4(5X,E12.6))
END SUBROUTINE vkrgk_print
!
!          Find veff-mtx and make transformation from CoM and
!          Rel system to lab system
!
SUBROUTINE vkrgk_free(gfree,twobody_com,twobody_r2, twobody_p2,lab_to_relcoeff, lab_to_relconf,&
     lab_to_relnumber,max_coeff,gmatrix_configs,isospin_z)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE partial_waves
  USE relcm_gmatrix
  IMPLICIT NONE
  TYPE (configuration_descriptor) , INTENT(IN)  :: gmatrix_configs
  INTEGER :: k2,k1,i,j,nlc1,nlc2, conf1, conf2, lcm1, lcm2, chan1, chan2, as, ls, ns,&
       max_coeff, isospin_z, ncm1, ncm2, lr1, lr2, nr1, nr2, bra, ket, state 
  REAL(KIND=8), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(INOUT) :: gfree
  REAL(KIND=8) :: dcoem, e_com, e_r2, e_rcm, e_p2, e_p2com
  REAL(KIND=8), DIMENSION(gmatrix_configs%number_confs,max_coeff), &
       INTENT(IN)  :: lab_to_relcoeff
  INTEGER, DIMENSION(gmatrix_configs%number_confs,max_coeff), &
       INTENT(IN)  :: lab_to_relconf
  INTEGER, DIMENSION(gmatrix_configs%number_confs), &
       INTENT(IN)  :: lab_to_relnumber
  REAL(KIND=8), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(INOUT) :: twobody_com, twobody_r2, twobody_p2

  gfree = 0.
  DO k2=1,gmatrix_configs%number_confs
     DO k1=1,gmatrix_configs%number_confs
        DO j=1,lab_to_relnumber(k2)
           nlc2=lab_to_relconf(k2,j)
           conf2=nlc2/1000;  chan2=nlc2-conf2*1000
           lcm2=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan2,conf2+conf2))
           ncm2=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan2,conf2+conf2))
           DO i=1,lab_to_relnumber(k1)
              nlc1=lab_to_relconf(k1,i)
              conf1=nlc1/1000;  chan1=nlc1-conf1*1000
              IF ( chan1 /= chan2 ) CYCLE
              lcm1=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan1,conf1+conf1))
              IF ( lcm1 /= lcm2 ) CYCLE   ! cm orbital mom must be equal
              ncm1=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan1,conf1+conf1))
              lr1=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan1,conf1+conf1-1))
              nr1=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan1,conf1+conf1-1))
              lr2=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan2,conf2+conf2-1))
              nr2=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan2,conf2+conf2-1))
              ! find bra and ket  states for nl and n'l'
              e_com = 0.0D0; e_r2 = 0.0D0; e_rcm = 0.0D0; e_p2=0.0D0; e_p2com= 0.0D0
              !             explicit computation of CoM and rel part
              IF ( (ncm1== ncm2) .AND.( lr1 == lr2 )) THEN
                 IF ( nr1 == nr2-1) e_r2 = -SQRT(nr2*(nr2+lr2+0.5))
                 IF ( nr1 == nr2-1) e_p2 = SQRT(nr2*(nr2+lr2+0.5))
                 IF ( nr1 == nr2) e_r2 =  2*nr2+lr2+1.5
                 IF ( nr1 == nr2) e_p2 =  2*nr2+lr2+1.5
                 IF ( nr1 == nr2+1) e_r2 = -SQRT((nr2+1.)*(nr2+lr2+1.5))
                 IF ( nr1 == nr2+1) e_p2 = SQRT((nr2+1.)*(nr2+lr2+1.5))
              ENDIF
              IF ( (nr1 == nr2) .AND. ( lr1 == lr2 ) ) THEN
                 IF ( ncm1 == ncm2-1) e_rcm = -SQRT(ncm2*(ncm2+lcm2+0.5))
                 IF ( ncm1 == ncm2-1) e_p2com = SQRT(ncm2*(ncm2+lcm2+0.5))
                 IF ( ncm1 == ncm2) e_rcm =  2*ncm2+lcm2+1.5
                 IF ( ncm1 == ncm2) e_p2com =  2*ncm2+lcm2+1.5
                 IF ( ncm1 == ncm2+1) e_rcm = -SQRT((ncm2+1.)*(ncm2+lcm2+1.5))
                 IF ( ncm1 == ncm2+1) e_p2com = SQRT((ncm2+1.)*(ncm2+lcm2+1.5))
              ENDIF
              IF ( ( lr1 == lr2 ).AND.(ncm1== ncm2).AND.( nr1 == nr2)) THEN
                 e_com =  2*ncm1+lcm1-2*nr2-lr2
              ENDIF
              bra = 0; ket = 0
              DO state =1, rel_conf%nconfs_rel(chan1)
                 as=rel_conf%rel_ab(chan1,state)
                 ns=relcm_sp_data%nrel(as)
                 ls=relcm_sp_data%lrel(as)
                 IF ( (ns == nr1) .AND. (ls == lr1)) ket = state
                 IF ( (ns == nr2) .AND. (ls == lr2)) bra = state
              ENDDO
              dcoem=lab_to_relcoeff(k1,i)*lab_to_relcoeff(k2,j)
              !  The interaction is diagonal in N
              !  this is not the case for  the G-matrix. Note also that the
              !  CoM part plus the twobody term -m\omega^2(r_i-r_j)^2/A is 
              !  independent of N, N' and L.
              IF ( ncm1 == ncm2 ) THEN  
                 twobody_com(k1,k2)=twobody_com(k1,k2)+dcoem*e_com
                 gfree(k1,k2)= gfree(k1,k2)+dcoem*v_com_rel(chan1,bra,ket) 
                 twobody_r2(k1,k2) = twobody_r2(k1,k2)+dcoem*e_r2
              ENDIF
              ! The p_ip_j correction depends however on N and N' and n and n'
              twobody_p2(k1,k2)=twobody_p2(k1,k2)+dcoem*(e_p2com-e_p2)*0.5D0  
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE vkrgk_free
!
!           Set up the interaction in the relative system using
!           a similarity transformation in momentum space. Quantum numbers are 
!           < nl JT_ZS | Veff | n'l'JT_ZS >
!
SUBROUTINE setup_com_vkrgk
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
     CALL vkrgk_channel(loop,v2body)
     DO ket = 1, pq_confs
        DO bra = 1, pq_confs
           v_com_rel(loop,bra,ket) = v2body(bra,ket)
        ENDDO
     ENDDO
     DEALLOCATE(v2body)
  ENDDO

END SUBROUTINE setup_com_vkrgk
!
!     Obtain the bare interaction in a harmonic oscillator 
!     basis plus the kinetic energy and the Coulomb part. It contains
!     also the CoM correction to the relative coordinates. The latter depends
!     on the mass number of the nucleus
! 
SUBROUTINE vkrgk_channel(i,vint)
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
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: heff(:,:), vzz(:,:)
  REAL(KIND = 8), ALLOCATABLE :: vkk(:,:)
  REAL(KIND = 8),  INTENT(INOUT):: &
       vint(rel_conf%nconfs_rel(i),rel_conf%nconfs_rel(i))

  integer :: jxx, sxx
  jxx=jang_rel(i)
  sxx = spin_rel(i)

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
  ALLOCATE(heff (ncoup*n_rel,ncoup*n_rel))
  heff =0.0D0
  ! get similarity transformed interaction  with soft cutoff
  CALL vkrgk_mtx(ncoup,vzz,heff,i)
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
           lim1=1; lim2=n_rel ; lim3=1 ; lim4=n_rel
        ELSEIF ( ncoup == 2 ) THEN
           IF ( (la == lb).AND. ( jang > la) ) THEN
              lim1=1; lim2=n_rel ; lim3=1 ; lim4=n_rel
           ELSEIF ( (la == lb).AND. ( jang < la) ) THEN
              lim1=1+n_rel; lim2=n_rel+n_rel ; lim3=1+n_rel ; lim4=n_rel+n_rel
           ELSEIF ( la >  lb ) THEN
              lim1=1+n_rel; lim2=n_rel+n_rel ; lim3=1 ; lim4=n_rel
           ELSEIF ( la <  lb ) THEN
              lim1=1; lim2=n_rel ; lim3=1+n_rel ; lim4=n_rel+n_rel
           ENDIF
        ENDIF
        CALL vkrgk_hosc(rnlr(:,lb,np), rnlr(:,la,n),heff(lim1:lim2,lim3:lim4),bra,ket,vsum)
        !       Here we set up the Coulomb part in an oscillator basis. 
        e_coulomb = 0.
        IF ( ( la == lb ).AND.( iso(i) == -1 ))  THEN
           e_coulomb = coulomb_relcom(n,np,la) 
        ENDIF
        !  only potential energy : coulomb + V_NN
        vint(ket,bra)=vsum+e_coulomb
!        write(6,'(6I4,2x,e20.8)') n,np,la,lb,sxx,jxx,vsum
     ENDDO
  ENDDO
  DEALLOCATE(heff); DEALLOCATE(vkk); DEALLOCATE(vzz)

END SUBROUTINE vkrgk_channel
!
!       Compute the T(G)-mtx for the actual channel 
!
SUBROUTINE vkrgk_mtx(ncoup,vzz,heff,ichan)
  USE wave_functions
  USE partial_waves
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncoup, ichan
  REAL(KIND=8), DIMENSION(ncoup*n_rel,ncoup*n_rel), INTENT(IN) :: vzz
  REAL(KIND=8), DIMENSION(ncoup*n_rel,ncoup*n_rel), INTENT(INOUT) :: heff
  INTEGER(KIND=4) :: i, j, n_ode, ij, iflag, i1, j1, tisoz
  INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: iwork
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: rg_vec, work
  REAL(KIND=8) :: relerr, abserr, k_lambda, start_point
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: cvec, cvecv, temp, tempv
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) ::  ceig, ceigv

  INTERFACE
     SUBROUTINE derivatives(k_lambda, v, dv, tiso)
       USE constants
       INTEGER ::  ij, i, j, i1, j2, kk, k, tiso
       REAL(DP) :: sum, k_lambda, v(:), dv(:), k1, k2, p
       REAL(DP), ALLOCATABLE :: vij(:,:)
     END SUBROUTINE derivatives
  END INTERFACE

  relerr = 1.0e-8_dp; abserr =  1.0e-8_dp 
  !  cutoff in MeV
  k_lambda  = k_cutoff*hbarc
  ! dimension of vectors and matrices
  n_total = ncoup*n_rel; n_ode = n_total*(n_total+1)/2
  ! total dim used by derivatives function
  ALLOCATE(rg_vec(n_ode)); ALLOCATE( work(100+21*n_ode)); ALLOCATE (iwork(5))
  rg_vec = 0.0_dp; work = 0.0_dp;  iwork = 0
  !  Transform the matrix -V- to a one-dim vector of dim ntot*(ntot+1)/2
  !  which is also the number of ordinary differential equations

  ALLOCATE(cvecv(n_total,n_total), ceigv(n_total), tempv(n_total,n_total)) 
  ! setup hamiltonian to be diagonalized
  tempv = dcmplx(0.0d0, 0.0d0); cvecv = dcmplx(0.0d0, 0.0d0); ceigv = dcmplx(0.0d0, 0.0d0)
  tisoz = iso(ichan) 
  DO i = 1, n_total
     i1 = i
     IF ( i > n_rel ) i1 = i-n_rel
     DO j = 1, n_total
        j1 = j
        IF ( j > n_rel ) j1 = j-n_rel
        tempv(j,i) = vzz(j,i)*SQRT( wkrel(i1)*wkrel(j1)) *krel(i1)*krel(j1)
     ENDDO
     tempv(i,i) = tempv(i,i) + ( krel(i1) * krel(i1) )/p_mass(tisoz)
  ENDDO
  !  Diagonalize  to test  if one gets the deuteron
  CALL vlowkdiag_exact( tempv, cvecv, ceigv, n_total )
  ij = 0  
  DO i = 1, n_total
     DO j = i, n_total
        ij = ij + 1
        rg_vec(ij) = vzz(j,i)
     ENDDO
  ENDDO
  iflag = 1
  start_point = krel(n_rel)
  CALL ode(derivatives,n_ode,rg_vec,start_point,k_lambda,relerr,abserr,iflag,work,iwork, tisoz)
  WRITE(6,*) 'iflag=', iflag
  IF ( iflag /= 2) WRITE(6,*) 'Error in ode, iflag not equal 2'
  !  Now transform back and get final effective interaction
  ij = 0  
  DO i = 1, n_total
     DO j = i, n_total
        ij = ij + 1
        heff(j,i) = rg_vec(ij)
        heff(i,j) = rg_vec(ij)
     ENDDO
  ENDDO
  DEALLOCATE(rg_vec); DEALLOCATE(work); DEALLOCATE(iwork)

  ALLOCATE(cvec(n_total,n_total), ceig(n_total), temp(n_total,n_total)) 
  ! setup hamiltonian to be diagonalized
  temp = dcmplx(0.0d0, 0.0d0); cvec = dcmplx(0.0d0, 0.0d0); ceig = dcmplx(0.0d0, 0.0d0)
  DO i = 1, n_total
     i1 = i
     IF ( i > n_rel ) i1 = i-n_rel
     DO j = 1, n_total
        j1 = j
        IF ( j > n_rel ) j1 = j-n_rel
        temp(j,i) = heff(j,i)*SQRT( wkrel(i1)*wkrel(j1)) *krel(i1)*krel(j1)
     ENDDO
     temp(i,i) = temp(i,i) + ( krel(i1) * krel(i1) )/p_mass(tisoz)
  ENDDO
  !  Diagonalize  to test  if one gets the deuteron
  CALL vlowkdiag_exact( temp, cvec, ceig, n_total )
  DO i = 1, n_total
     WRITE(6,*) i, REAL(ceig(i)), REAL(ceigv(i))
  ENDDO
  DEALLOCATE(ceig); DEALLOCATE(cvec, temp)
  DEALLOCATE(ceigv); DEALLOCATE(cvecv, tempv)
  DO i = 1, n_total
    write(6,*) krel(i), heff(i,i), vzz(i,i)
  ENDDO

END SUBROUTINE vkrgk_mtx
!
!fast diagonalization of complex symmetric matrices
!
SUBROUTINE vkrgkdiag_exact( h, cvec, ceig, n )
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

END SUBROUTINE vkrgkdiag_exact
!
!          Obtain the bare potential in oscillator basis
!          Final V is in units of MeV or eV
!          Time consuming part. Parallel version available.
!
SUBROUTINE vkrgk_hosc(wave_bra, wave_ket, a,bra,ket,vsum)
  USE wave_functions
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  REAL(KIND = 8), INTENT(INOUT) :: vsum
  INTEGER, INTENT(IN) :: bra, ket
  INTEGER :: i, j
  REAL(KIND = 8), DIMENSION(n_rel,n_rel), INTENT(IN) :: a
  REAL(KIND = 8) :: sum1, sum2, hbarc3
  REAL(KIND = 8), DIMENSION(n_rel), INTENT(IN) :: wave_bra, wave_ket 

  hbarc3=hbarc**3
  sum1=0.
  DO i=1, n_rel
     sum2=0.
     DO j=1, n_rel
        sum2=sum2+wave_ket(j)*a(j,i)
     ENDDO
     sum1=sum1+sum2*wave_bra(i)
  ENDDO
  vsum = sum1*hbarc3

END SUBROUTINE vkrgk_hosc


SUBROUTINE derivatives(k_lambda, v, dv, tisoz)
  USE wave_functions
  USE constants
  IMPLICIT NONE
  INTEGER ::  ij, i, j, i1, j2, kk, k, i2, tisoz
  REAL(DP) :: sum, k_lambda, v(:), dv(:), k1, k2, p
  REAL(DP), ALLOCATABLE :: vij(:,:)

  ALLOCATE(vij(n_total,n_total)) 
  ij = 0
  DO j = 1, n_total
     DO i = j, n_total
        ij = ij + 1
        vij(i,j) = v(ij)
        vij(j,i) = v(ij)
     ENDDO
  ENDDO
  ij = 0
  DO i = 1, n_total
     i1 = i
     IF ( i > n_rel ) i1 = i-n_rel
     k1 = (krel(i1)**2)/p_mass(tisoz)
     DO j = i, n_total
        i2 = j
        IF ( j > n_rel ) i2 = j-n_rel
        k2 = (krel(i2)**2)/p_mass(tisoz)
        ij = ij + 1
        sum = 0_dp
        DO k = 1, n_total
           kk = k
           IF(k > n_rel) kk = k - n_rel 
           p  = (krel(kk)**2)
           sum = sum +  (k1+k2-2.0_dp*p/p_mass(tisoz))*p*wkrel(kk)*vij(j,k)*vij(k,i) 
        ENDDO
        dv(ij) =  sum - (k2-k1)*(k2-k1)*vij(j,i)
     ENDDO
  ENDDO
  dv = -4.0_dp*dv*(p_mass(tisoz)**2)/(k_lambda**5)
  DEALLOCATE(vij)

END SUBROUTINE derivatives

!
!  This function is taken from www.netlib.org
!
SUBROUTINE ode(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork,isotz)
  IMPLICIT REAL*8(a-h,o-z)
  !
  !   double precision subroutine ode integrates a system of neqn
  !   first order ordinary differential equations of the form:
  !             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
  !             y(i) given at  t .
  !   the subroutine integrates from  t  to  tout .  on return the
  !   parameters in the call list are set for continuing the integration.
  !   the user has only to define a new value  tout  and call  ode  again.
  !
  !   the differential equations are actually solved by a suite of codes
  !   de ,  step , and  intrp .  ode  allocates virtual storage in the
  !   arrays  work  and  iwork  and calls  de .  de  is a supervisor which
  !   directs the solution.  it calls on the routines  step  and  intrp
  !   to advance the integration and to interpolate at output points.
  !   step  uses a modified divided difference form of the adams pece
  !   formulas and local extrapolation.  it adjusts the order and step
  !   size to control the local error per unit step in a generalized
  !   sense.  normally each call to  step  advances the solution one step
  !   in the direction of  tout .  for reasons of efficiency  de
  !   integrates beyond  tout  internally, though never beyond
  !   t+10*(tout-t), and calls  intrp  to interpolate the solution at
  !   tout .  an option is provided to stop the integration at  tout  but
  !   it should be used only if it is impossible to continue the
  !   integration beyond  tout .
  !
  !   this code is completely explained and documented in the text,
  !   computer solution of ordinary differential equations:  the initial
  !   value problem  by l. f. shampine and m. k. gordon.
  !
  !   the parameters represent:
  !      f -- double precision subroutine f(t,y,yp) to evaluate
  !                derivatives yp(i)=dy(i)/dt
  !      neqn -- number of equations to be integrated (integer*4)
  !      y(*) -- solution vector at t                 (real*8)
  !      t -- independent variable                    (real*8)
  !      tout -- point at which solution is desired   (real*8)
  !      relerr,abserr -- relative and absolute error tolerances for local
  !           error test (real*8).  at each step the code requires
  !             dabs(local error) .le. dabs(y)*relerr + abserr
  !           for each component of the local error and solution vectors
  !      iflag -- indicates status of integration     (integer*4)
  !      work(*)  (real*8)  -- arrays to hold information internal to
  !      iwork(*) (integer*4)    which is necessary for subsequent calls
  !
  !   first call to ode --
  !
  !   the user must provide storage in his calling program for the arrays
  !   in the call list,
  !      y(neqn), work(100+21*neqn), iwork(5),
  !   declare  f  in an external statement, supply the double precision
  !   subroutine f(t,y,yp)  to evaluate
  !      dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn))
  !   and initialize the parameters:
  !      neqn -- number of equations to be integrated
  !      y(*) -- vector of initial conditions
  !      t -- starting point of integration
  !      tout -- point at which solution is desired
  !      relerr,abserr -- relative and absolute local error tolerances
  !      iflag -- +1,-1.  indicator to initialize the code.  normal input
  !           is +1.  the user should set iflag=-1 only if it is
  !           impossible to continue the integration beyond  tout .
  !   all parameters except  f ,  neqn  and  tout  may be altered by the
  !   code on output so must be variables in the calling program.
  !
  !   output from  ode  --
  !
  !      neqn -- unchanged
  !      y(*) -- solution at  t
  !      t -- last point reached in integration.  normal return has
  !           t = tout .
  !      tout -- unchanged
  !      relerr,abserr -- normal return has tolerances unchanged.  iflag=3
  !           signals tolerances increased
  !      iflag = 2 -- normal return.  integration reached  tout
  !            = 3 -- integration did not reach  tout  because error
  !                   tolerances too small.  relerr ,  abserr  increased
  !                   appropriately for continuing
  !            = 4 -- integration did not reach  tout  because more than
  !                   500 steps needed
  !            = 5 -- integration did not reach  tout  because equations
  !                   appear to be stiff
  !            = 6 -- invalid input parameters (fatal error)
  !           the value of  iflag  is returned negative when the input
  !           value is negative and the integration does not reach  tout ,
  !           i.e., -3, -4, -5.
  !      work(*),iwork(*) -- information generally of no interest to the
  !           user but necessary for subsequent calls.
  !
  !   subsequent calls to  ode --
  !
  !   subroutine  ode  returns with all information needed to continue
  !   the integration.  if the integration reached  tout , the user need
  !   only define a new  tout  and call again.  if the integration did not
  !   reach  tout  and the user wants to continue, he just calls again.
  !   the output value of  iflag  is the appropriate input value for
  !   subsequent calls.  the only situation in which it should be altered
  !   is to stop the integration internally at the new  tout , i.e.,
  !   change output  iflag=2  to input  iflag=-2 .  error tolerances may
  !   be changed by the user before continuing.  all other parameters must
  !   remain unchanged.
  !
  !***********************************************************************
  !*  subroutines  de  and  step  contain machine dependent constants. *
  !*  be sure they are set before using  ode .                          *
  !***********************************************************************
  !
  LOGICAL start,phase1,nornd
  DIMENSION y(neqn),work(1),iwork(5)
  DATA ialpha,ibeta,isig,iv,iw,ig,iphase,ipsi,ix,ih,ihold,istart,  &
       itold,idelsn/1,13,25,38,50,62,75,76,88,89,90,91,92,93/

  INTERFACE
     SUBROUTINE f(lam, v, dv, tisoz)
       USE constants
       INTEGER ::  ij, i, j, ii, jj, kk, k, tisoz
       REAL(DP) :: ans, lam, v(:), dv(:), k1, k2, p
       REAL(DP), ALLOCATABLE :: vij(:,:)
     END SUBROUTINE f
  END INTERFACE

  iyy = 100
  iwt = iyy + neqn
  ip = iwt + neqn
  iyp = ip + neqn
  iypout = iyp + neqn
  iphi = iypout + neqn
  IF(iabs(iflag) /= 1) THEN
     start = work(istart) > 0.0d0
     phase1 = work(iphase) > 0.0d0
     nornd = iwork(2) /= -1
  ELSE
     CALL de(f,neqn,y,t,tout,relerr,abserr,iflag,work(iyy), & 
          work(iwt),work(ip),work(iyp),work(iypout),work(iphi),  &
          work(ialpha),work(ibeta),work(isig),work(iv),work(iw),work(ig), &
          phase1,work(ipsi),work(ix),work(ih),work(ihold),start,  &
          work(itold),work(idelsn),iwork(1),nornd,iwork(3),iwork(4),  &
          iwork(5),isotz)
  ENDIF
  work(istart) = -1.0d0
  IF(start) work(istart) = 1.0d0
  work(iphase) = -1.0d0
  IF(phase1) work(iphase) = 1.0d0
  iwork(2) = -1
  IF(nornd) iwork(2) = 1

END SUBROUTINE  ode


SUBROUTINE de(f,neqn,y,t,tout,relerr,abserr,iflag,  &
     yy,wt,p,yp,ypout,phi,alpha,beta,sig,v,w,g,phase1,psi,x,h,hold, &
     start,told,delsgn,ns,nornd,k,kold,isnold,isotz)
  IMPLICIT REAL*8(a-h,o-z)
  !
  !   ode  merely allocates storage for  de  to relieve the user of the
  !   inconvenience of a long call list.  consequently  de  is used as
  !   described in the comments for  ode .
  !
  !   this code is completely explained and documented in the text,
  !   computer solution of ordinary differential equations:  the initial
  !   value problem  by l. f. shampine and m. k. gordon.
  !
  LOGICAL stiff,crash,start,phase1,nornd
  DIMENSION y(neqn),yy(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),  &
       ypout(neqn),psi(12),alpha(12),beta(12),sig(13),v(12),w(12),g(13)

  INTERFACE
     SUBROUTINE f(lam, v, dv, tisoz)
       USE constants
       INTEGER ::  ij, i, j, ii, jj, kk, k, tisoz
       REAL(DP) :: ans, lam, v(:), dv(:), k1, k2, p
       REAL(DP), ALLOCATABLE :: vij(:,:)
     END SUBROUTINE f
  END INTERFACE
  !
  !***********************************************************************
  !*  the only machine dependent constant is based on the machine unit   *
  !*  roundoff error  u  which is the smallest positive number such that *
  !*  1.0+u .gt. 1.0 .  u  must be calculated and  fouru=4.0*u  inserted *
  !*  in the following data statement before using  de .  the routine    *
  !*  machin  calculates  u .  fouru  and  twou=2.0*u  must also be      *
  !*  inserted in subroutine  step  before calling  de .                 *
  !     data fouru/.888d-15/                                              ***
  !***********************************************************************
  !
  !   the constant  maxnum  is the maximum number of steps allowed in one
  !   call to  de .  the user may change this limit by altering the
  !   following statement
  DATA maxnum/5000000/


  fouru = 4.0 * d1mach(4)                                          ! ***
  IF(neqn <  1) go to 10
  IF(t ==  tout) go to 10
  IF(relerr  <  0.0d0  .OR.  abserr .LT. 0.0d0) go to 10
  eps = dmax1(relerr,abserr)
  IF(eps  <=  0.0d0) go to 10
  IF(iflag ==  0) go to 10
  isn = isign(1,iflag)
  iflag = iabs(iflag)
  IF(iflag ==  1) go to 20
  IF(t /=  told) go to 10
  IF(iflag .GE. 2  .AND.  iflag .LE. 5) go to 20
10 iflag = 6
  RETURN
  !
  !   on each call set interval of integration and counter for number of
  !   steps.  adjust input error tolerances to define weight vector for
  !   subroutine  step
  !
20 del = tout - t
  absdel = dabs(del)
  tend = t + 10.0d0*del
  IF(isn .LT. 0) tend = tout
  nostep = 0
  kle4 = 0
  stiff = .FALSE.
  releps = relerr/eps
  abseps = abserr/eps
  IF(iflag .EQ. 1) go to 30
  IF(isnold .LT. 0) go to 30
  IF(delsgn*del .GT. 0.0d0) go to 50
  !
  !   on start and restart also set work variables x and yy(*), store the
  !   direction of integration and initialize the step size
  !
30 start = .TRUE.
  x = t
  DO l = 1,neqn
     yy(l) = y(l)
  ENDDO
  delsgn = dsign(1.0d0,del)
  h = dsign(dmax1(dabs(tout-x),fouru*dabs(x)),tout-x)
  !
  !   if already past output point, interpolate and return
  !
50 IF(dabs(x-t) .LT. absdel) go to 60
  CALL intrp(x,yy,tout,y,ypout,neqn,kold,phi,psi)
  iflag = 2
  t = tout
  told = t
  isnold = isn
  RETURN
  !
  !   if cannot go past output point and sufficiently close,
  !   extrapolate and return
  !
60 IF(isn .GT. 0  .OR.  dabs(tout-x) .GE. fouru*dabs(x)) go to 80
  h = tout - x
  CALL f(x,yy,yp,isotz)
  DO l = 1,neqn
     y(l) = yy(l) + h*yp(l)
  ENDDO
  iflag = 2
  t = tout
  told = t
  isnold = isn
  RETURN
  !
  !   test for too many steps
  !
80 IF(nostep .LT. maxnum) go to 100
  iflag = isn*4
  IF(stiff) iflag = isn*5
  DO l = 1,neqn
     y(l) = yy(l)
  ENDDO
  t = x
  told = t
  isnold = 1
  RETURN
  !
  !   limit step size, set weight vector and take a step
  !
100 h = dsign(dmin1(dabs(h),dabs(tend-x)),h)
  DO l = 1,neqn
     wt(l) = releps*dabs(yy(l)) + abseps
  ENDDO
  CALL step(x,yy,f,neqn,h,eps,wt,start,  &
       hold,k,kold,crash,phi,p,yp,psi,         &
       alpha,beta,sig,v,w,g,phase1,ns,nornd,isotz)
  !
  !   test for tolerances too small
  !
  IF(.NOT.crash) go to 130
  iflag = isn*3
  relerr = eps*releps
  abserr = eps*abseps
  DO l = 1,neqn
     y(l) = yy(l)
  ENDDO
  t = x
  told = t
  isnold = 1
  RETURN
  !
  !   augment counter on number of steps and test for stiffness
  !
130 nostep = nostep + 1
  kle4 = kle4 + 1
  IF(kold .GT. 4) kle4 = 0
  IF(kle4 .GE. 50) stiff = .TRUE.
  go to 50

END SUBROUTINE  de



  SUBROUTINE step(x,y,f,neqn,h,eps,wt,start,  &
       hold,k,kold,crash,phi,p,yp,psi,  &
       alpha,beta,sig,v,w,g,phase1,ns,nornd,isotz)
    IMPLICIT REAL*8(a-h,o-z)
    !
    !   double precision subroutine  step
    !   integrates a system of first order ordinary
    !   differential equations one step, normally from x to x+h, using a
    !   modified divided difference form of the adams pece formulas.  local
    !   extrapolation is used to improve absolute stability and accuracy.
    !   the code adjusts its order and step size to control the local error
    !   per unit step in a generalized sense.  special devices are included
    !   to control roundoff error and to detect when the user is requesting
    !   too much accuracy.
    !
    !   this code is completely explained and documented in the text,
    !   computer solution of ordinary differential equations:  the initial
    !   value problem  by l. f. shampine and m. k. gordon.
    !
    !
    !   the parameters represent:
    !      x -- independent variable             (real*8)
    !      y(*) -- solution vector at x          (real*8)
    !      yp(*) -- derivative of solution vector at  x  after successful
    !           step                             (real*8)
    !      neqn -- number of equations to be integrated (integer*4)
    !      h -- appropriate step size for next step.  normally determined by
    !           code                             (real*8)
    !      eps -- local error tolerance.  must be variable  (real*8)
    !      wt(*) -- vector of weights for error criterion   (real*8)
    !      start -- logical variable set .true. for first step,  .false.
    !           otherwise                        (logical*4)
    !      hold -- step size used for last successful step  (real*8)
    !      k -- appropriate order for next step (determined by code)
    !      kold -- order used for last successful step
    !      crash -- logical variable set .true. when no step can be taken,
    !           .false. otherwise.
    !   the arrays  phi, psi  are required for the interpolation subroutine
    !   intrp.  the array p is internal to the code.  all are real*8
    !
    !   input to  step
    !
    !      first call --
    !
    !   the user must provide storage in his driver program for all arrays
    !   in the call list, namely
    !
    !     dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
    !
    !   the user must also declare  start  and  crash  logical variables
    !   and  f  an external subroutine, supply the subroutine  f(x,y,yp)
    !   to evaluate
    !      dy(i)/dx = yp(i) = f(x,y(1),y(2),...,y(neqn))
    !   and initialize only the following parameters:
    !      x -- initial value of the independent variable
    !      y(*) -- vector of initial values of dependent variables
    !      neqn -- number of equations to be integrated
    !      h -- nominal step size indicating direction of integration
    !           and maximum size of step.  must be variable
    !      eps -- local error tolerance per step.  must be variable
    !      wt(*) -- vector of non-zero weights for error criterion
    !      start -- .true.
    !
    !   step  requires the l2 norm of the vector with components
    !   local error(l)/wt(l)  be less than  eps  for a successful step.  the
    !   array  wt  allows the user to specify an error test appropriate
    !   for his problem.  for example,
    !      wt(l) = 1.0  specifies absolute error,
    !            = dabs(y(l))  error relative to the most recent value of
    !                 the l-th component of the solution,
    !            = dabs(yp(l))  error relative to the most recent value of
    !                 the l-th component of the derivative,
    !            = dmax1(wt(l),dabs(y(l)))  error relative to the largest
    !                 magnitude of l-th component obtained so far,
    !            = dabs(y(l))*relerr/eps + abserr/eps  specifies a mixed
    !                 relative-absolute test where  relerr  is relative
    !                 error,  abserr  is absolute error and  eps =
    !                 dmax1(relerr,abserr) .
    !
    !      subsequent calls --
    !
    !   subroutine  step  is designed so that all information needed to
    !   continue the integration, including the step size  h  and the order
    !   k , is returned with each step.  with the exception of the step
    !   size, the error tolerance, and the weights, none of the parameters
    !   should be altered.  the array  wt  must be updated after each step
    !   to maintain relative error tests like those above.  normally the
    !   integration is continued just beyond the desired endpoint and the
    !   solution interpolated there with subroutine  intrp .  if it is
    !   impossible to integrate beyond the endpoint, the step size may be
    !   reduced to hit the endpoint since the code will not take a step
    !   larger than the  h  input.  changing the direction of integration,
    !   i.e., the sign of  h , requires the user set  start = .true. before
    !   calling  step  again.  this is the only situation in which  start
    !   should be altered.
    !
    !   output from  step
    !
    !      successful step --
    !
    !   the subroutine returns after each successful step with  start  and
    !   crash  set .false. .  x  represents the independent variable
    !   advanced one step of length  hold  from its value on input and  y
    !   the solution vector at the new value of  x .  all other parameters
    !   represent information corresponding to the new  x  needed to
    !   continue the integration.
    !
    !      unsuccessful step --
    !
    !   when the error tolerance is too small for the machine precision,
    !   the subroutine returns without taking a step and  crash = .true. .
    !   an appropriate step size and error tolerance for continuing are
    !   estimated and all other information is restored as upon input
    !   before returning.  to continue with the larger tolerance, the user
    !   just calls the code again.  a restart is neither required nor
    !   desirable.
    !
    LOGICAL start,crash,phase1,nornd
    DIMENSION y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
    DIMENSION alpha(12),beta(12),sig(13),w(12),v(12),g(13),  &
         gstr(13),two(13)

    INTERFACE
       SUBROUTINE f(lam, v, dv, tisoz)
         USE constants
         INTEGER ::  ij, i, j, ii, jj, kk, k, tisoz
         REAL(DP) :: ans, lam, v(:), dv(:), k1, k2, p
         REAL(DP), ALLOCATABLE :: vij(:,:)
       END SUBROUTINE f
    END INTERFACE
    !***********************************************************************
    !*  the only machine dependent constants are based on the machine unit *
    !*  roundoff error  u  which is the smallest positive number such that *
    !*  1.0+u .gt. 1.0  .  the user must calculate  u  and insert          *
    !*  twou=2.0*u  and  fouru=4.0*u  in the data statement before calling *
    !*  the code.  the routine  machin  calculates  u .                    *
    !     data twou,fouru/.444d-15,.888d-15/                                ***
    !***********************************************************************
    DATA two/2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,256.0d0,  &
         512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/
    DATA gstr/0.500d0,0.0833d0,0.0417d0,0.0264d0,0.0188d0,0.0143d0, &
         0.0114d0,0.00936d0,0.00789d0,0.00679d0,0.00592d0,0.00524d0,  &
         0.00468d0/

    twou = 2.0 * d1mach(4)                                           ! ***
    fouru = 2.0 * twou                                               ! ***
    !
    !   if step size is too small, determine an acceptable one
    !
    crash = .TRUE.
      IF(dabs(h) .GE. fouru*dabs(x)) go to 5
      h = dsign(fouru*dabs(x),h)
      RETURN
    5 p5eps = 0.5d0*eps
!
!   if error tolerance is too small, increase it to an acceptable value
!
      round = 0.0d0
      DO l = 1,neqn
      round = round + (y(l)/wt(l))**2
      ENDDO
      round = twou*dsqrt(round)
      IF(p5eps .GE. round) go to 15
      eps = 2.0*round*(1.0d0 + fouru)
      RETURN
   15 crash = .FALSE.
      g(1)=1.0d0
      g(2)=0.5d0
      sig(1)=1.0d0
      IF(.NOT.start) go to 99
!
!   initialize.  compute appropriate step size for first step
!
      CALL f(x,y,yp,isotz)
      sum = 0.0d0
      DO l = 1,neqn
        phi(l,1) = yp(l)
        phi(l,2) = 0.0d0
        sum = sum + (yp(l)/wt(l))**2
      ENDDO
      sum = dsqrt(sum)
      absh = dabs(h)
      IF(eps <  16.0d0*sum*h*h) absh = 0.25d0*dsqrt(eps/sum)
      h = dsign(dmax1(absh,fouru*dabs(x)),h)
      hold = 0.0d0
      k = 1
      kold = 0
      start = .FALSE.
      phase1 = .TRUE.
      nornd = .TRUE.
      IF(p5eps .GT. 100.0d0*round) go to 99
      nornd = .FALSE.
      DO l = 1,neqn
       phi(l,15) = 0.0d0
      ENDDO
   99 ifail = 0
!       ***     end block 0     ***
!
!       ***     begin block 1     ***
!   compute coefficients of formulas for this step.  avoid computing
!   those quantities not changed when step size is not changed.
!                   ***
!
  100 kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
!
!   ns is the number of steps taken with size h, including the current
!   one.  when k.lt.ns, no coefficients change
!
      IF(h .NE. hold) ns = 0
      IF(ns.LE.kold)   ns=ns+1
      nsp1 = ns+1
      IF (k .LT. ns) go to 199
!
!   compute those components of alpha(*),beta(*),psi(*),sig(*) which
!   are changed
!
      beta(ns) = 1.0d0
      realns = ns
      alpha(ns) = 1.0d0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0d0
      IF(k <  nsp1) go to 110
      DO i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
        sig(i+1) = reali*alpha(i)*sig(i)
      ENDDO
  110 psi(k) = temp1
!
!   compute coefficients g(*)
!
!   initialize v(*) and set w(*).  g(2) is set in data statement
!
      IF(ns >  1) go to 120
      DO iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0d0/temp3
        w(iq) = v(iq)
      ENDDO
      go to 140
!
!   if order was raised, update diagonal part of v(*)
!
  120 IF(k <=  kold) go to 130
      temp4 = k*kp1
      v(k) = 1.0d0/temp4
      nsm2 = ns-2
      IF(nsm2 .LT. 1) go to 130
      DO j = 1,nsm2
        i = k-j
       v(i) = v(i) - alpha(j+1)*v(i+1)
      ENDDO
!
!   update v(*) and set w(*)
!
  130 limit1 = kp1 - ns
      temp5 = alpha(ns)
      DO iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
        w(iq) = v(iq)
      ENDDO
      g(nsp1) = w(1)
!
!   compute the g(*) in the work vector w(*)
!
  140 nsp2 = ns + 2
      IF(kp1 .LT. nsp2) go to 199
      DO 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        DO 145 iq = 1,limit2
  145     w(iq) = w(iq) - temp6*w(iq+1)
  150   g(i) = w(1)
  199   CONTINUE
!       ***     end block 1     ***
!
!       ***     begin block 2     ***
!   predict a solution p(*), evaluate derivatives using predicted
!   solution, estimate local error at order k and errors at orders k,
!   k-1, k-2 as if constant step size were used.
!                   ***
!
!   change phi to phi star
!
      IF(k .LT. nsp1) go to 215
      DO 210 i = nsp1,k
        temp1 = beta(i)
        DO 205 l = 1,neqn
  205     phi(l,i) = temp1*phi(l,i)
  210   CONTINUE
!
!   predict solution and differences
!
  215 DO 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0d0
  220   p(l) = 0.0d0
      DO 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        DO 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
  225     phi(l,i) = phi(l,i) + phi(l,ip1)
  230   CONTINUE
      IF(nornd) go to 240
      DO 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
  235   phi(l,16) = (p(l) - y(l)) - tau
      go to 250
  240 DO 245 l = 1,neqn
  245   p(l) = y(l) + h*p(l)
  250 xold = x
      x = x + h
      absh = dabs(h)
      CALL f(x,p,yp,isotz)
!
!   estimate errors at orders k,k-1,k-2
!
      erkm2 = 0.0d0
      erkm1 = 0.0d0
      erk = 0.0d0
      DO 265 l = 1,neqn
        temp3 = 1.0d0/wt(l)
        temp4 = yp(l) - phi(l,1)
        IF(km2)265,260,255
  255   erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
  260   erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
  265   erk = erk + (temp4*temp3)**2
      IF(km2)280,275,270
  270 erkm2 = absh*sig(km1)*gstr(km2)*dsqrt(erkm2)
  275 erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
  280 temp5 = absh*dsqrt(erk)
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      knew = k
!
!   test if order should be lowered
!
      IF(km2)299,290,285
  285 IF(dmax1(erkm1,erkm2) .LE. erk) knew = km1
      go to 299
  290 IF(erkm1 .LE. 0.5d0*erk) knew = km1
!
!   test if step successful
!
  299 IF(err .LE. eps) go to 400
!       ***     end block 2     ***
!
!       ***     begin block 3     ***
!   the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
!   if third consecutive failure, set order to one.  if step fails more
!   than three times, consider an optimal step size.  double error
!   tolerance and return if estimated step size is too small for machine
!   precision.
!                   ***

!   restore x, phi(*,*) and psi(*)

      phase1 = .FALSE.
      x = xold
      DO 310 i = 1,k
        temp1 = 1.0d0/beta(i)
        ip1 = i+1
        DO 305 l = 1,neqn
  305     phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
  310   CONTINUE
      IF(k .LT. 2) go to 320
      DO 315 i = 2,k
  315   psi(i-1) = psi(i) - h
!
!   on third failure, set order to one.  thereafter, use optimal step
!   size
!
  320 ifail = ifail + 1
      temp2 = 0.5d0
      IF(ifail - 3) 335,330,325
  325 IF(p5eps .LT. 0.25d0*erk) temp2 = dsqrt(p5eps/erk)
  330 knew = 1
  335 h = temp2*h
      k = knew
      IF(dabs(h) .GE. fouru*dabs(x)) go to 340
      crash = .TRUE.
      h = dsign(fouru*dabs(x),h)
      eps = eps + eps
      RETURN
  340 go to 100
!       ***     end block 3     ***
!
!       ***     begin block 4     ***
!   the step is successful.  correct the predicted solution, evaluate
!   the derivatives using the corrected solution and update the
!   differences.  determine best order and step size for next step.
!                   ***
  400 kold = k
      hold = h
!
!   correct and evaluate
!
      temp1 = h*g(kp1)
      IF(nornd) go to 410
      DO 405 l = 1,neqn
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
  405   phi(l,15) = (y(l) - p(l)) - rho
      go to 420
  410 DO 415 l = 1,neqn
  415   y(l) = p(l) + temp1*(yp(l) - phi(l,1))
  420 CALL f(x,y,yp,isotz)
!
!   update differences for next step
!
      DO 425 l = 1,neqn
        phi(l,kp1) = yp(l) - phi(l,1)
  425   phi(l,kp2) = phi(l,kp1) - phi(l,kp2)
      DO 435 i = 1,k
        DO 430 l = 1,neqn
  430     phi(l,i) = phi(l,i) + phi(l,kp1)
  435   CONTINUE
!
!   estimate error at order k+1 unless:
!     in first phase when always raise order,
!     already decided to lower order,
!     step size not constant so estimate unreliable
!
      erkp1 = 0.0d0
      IF(knew .EQ. km1  .OR.  k .EQ. 12) phase1 = .FALSE.
      IF(phase1) go to 450
      IF(knew .EQ. km1) go to 455
      IF(kp1 .GT. ns) go to 460
      DO 440 l = 1,neqn
  440   erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
      erkp1 = absh*gstr(kp1)*dsqrt(erkp1)
!
!   using estimated error at order k+1, determine appropriate order
!   for next step
!
      IF(k .GT. 1) go to 445
      IF(erkp1 .GE. 0.5d0*erk) go to 460
      go to 450
  445 IF(erkm1 .LE. dmin1(erk,erkp1)) go to 455
      IF(erkp1 .GE. erk  .OR.  k .EQ. 12) go to 460
!
!   here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would have
!   been lowered in block 2.  thus order is to be raised
!
!   raise order
!
  450 k = kp1
      erk = erkp1
      go to 460
!
!   lower order
!
  455 k = km1
      erk = erkm1
!
!   with new order determine appropriate step size for next step
!
  460 hnew = h + h
      IF(phase1) go to 465
      IF(p5eps .GE. erk*two(k+1)) go to 465
      hnew = h
      IF(p5eps .GE. erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0d0/temp2)
      hnew = absh*dmax1(0.5d0,dmin1(0.9d0,r))
      hnew = dsign(dmax1(hnew,fouru*dabs(x)),h)
  465 h = hnew
      RETURN
!       ***     end block 4     ***

  END SUBROUTINE  step


  SUBROUTINE intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)
    IMPLICIT REAL*8(a-h,o-z)
    !
    !   the methods in subroutine  step  approximate the solution near  x
    !   by a polynomial.  subroutine  intrp  approximates the solution at
    !   xout  by evaluating the polynomial there.  information defining this
    !   polynomial is passed from  step  so  intrp  cannot be used alone.
    !
    !   this code is completely explained and documented in the text,
    !   computer solution of ordinary differential equations:  the initial
    !   value problem  by l. f. shampine and m. k. gordon.
    !
    !   input to intrp --
    !
    !   all floating point variables are double precision
    !   the user provides storage in the calling program for the arrays in
    !   the call list
    DIMENSION y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)
    !   and defines
    !      xout -- point at which solution is desired.
    !   the remaining parameters are defined in  step  and passed to  intrp
    !   from that subroutine
    !
    !   output from  intrp --
    !
    !      yout(*) -- solution at  xout
    !      ypout(*) -- derivative of solution at  xout
    !   the remaining parameters are returned unaltered from their input
    !   values.  integration with  step  may be continued.
    !
    DIMENSION g(13),w(13),rho(13)
    DATA g(1)/1.0d0/,rho(1)/1.0d0/
    !
    hi = xout - x
    ki = kold + 1
    kip1 = ki + 1
    !
    !   initialize w(*) for computing g(*)
    !
    DO  i = 1,ki
       temp1 = i
       w(i) = 1.0d0/temp1
    ENDDO
    term = 0.0d0
    !
    !   compute g(*)
    !
    DO j = 2,ki
       jm1 = j - 1
       psijm1 = psi(jm1)
       gamma = (hi + term)/psijm1
       eta = hi/psijm1
       limit1 = kip1 - j
       DO  i = 1,limit1
          w(i) = gamma*w(i) - eta*w(i+1)
       ENDDO
       g(j) = w(1)
       rho(j) = gamma*rho(jm1)
       term = psijm1
    ENDDO
    !
    !   interpolate
    !
    DO  l = 1,neqn
       ypout(l) = 0.0d0
       yout(l) = 0.0d0
    ENDDO
    DO  j = 1,ki
       i = kip1 - j
       temp2 = g(i)
       temp3 = rho(i)
       DO l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
          ypout(l) = ypout(l) + temp3*phi(l,i)
       ENDDO
    ENDDO
    DO l = 1,neqn
       yout(l) = y(l) + hi*yout(l)
    ENDDO

  END SUBROUTINE  intrp


  DOUBLE PRECISION FUNCTION d1mach (i)
    IMPLICIT NONE
    INTEGER :: i
    DOUBLE PRECISION :: b, x
    !***begin prologue  d1mach
    !***purpose  return floating point machine dependent constants.
    !***library   slatec
    !***category  r1
    !***type      single precision (d1mach-s, d1mach-d)
    !***keywords  machine constants
    !***author  fox, p. a., (bell labs)
    !           hall, a. d., (bell labs)
    !           schryer, n. l., (bell labs)
    !***description
    !
    !   d1mach can be used to obtain machine-dependent parameters for the
    !   local machine environment.  it is a function subprogram with one
    !   (input) argument, and can be referenced as follows:
    !
    !        a = d1mach(i)
    !
    !   where i=1,...,5.  the (output) value of a above is determined by
    !   the (input) value of i.  the results for various values of i are
    !   discussed below.
    !
    !   d1mach(1) = b**(emin-1), the smallest positive magnitude.
    !   d1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
    !   d1mach(3) = b**(-t), the smallest relative spacing.
    !   d1mach(4) = b**(1-t), the largest relative spacing.
    !   d1mach(5) = log10(b)
    !
    !   assume single precision numbers are represented in the t-digit,
    !   base-b form
    !
    !              sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
    !
    !   where 0 .le. x(i) .lt. b for i=1,...,t, 0 .lt. x(1), and
    !   emin .le. e .le. emax.
    !
    !   the values of b, t, emin and emax are provided in i1mach as
    !   follows:
    !   i1mach(10) = b, the base.
    !   i1mach(11) = t, the number of base-b digits.
    !   i1mach(12) = emin, the smallest exponent e.
    !   i1mach(13) = emax, the largest exponent e.
    !
    !
    !***references  p. a. fox, a. d. hall and n. l. schryer, framework for
    !                 a portable library, acm transactions on mathematical
    !                 software 4, 2 (june 1978), pp. 177-188.
    !***routines called  xermsg
    !***revision history  (yymmdd)
    !   790101  date written
    !   960329  modified for fortran 90 (be after suggestions by ehg)      
    !***end prologue  d1mach
    !      
    x = 1.0d0
    b = RADIX(x)
    SELECT CASE (i)
    CASE (1)
       d1mach = b**(MINEXPONENT(x)-1) ! the smallest positive magnitude.
    CASE (2)
       d1mach = HUGE(x)               ! the largest magnitude.
    CASE (3)
       d1mach = b**(-DIGITS(x))       ! the smallest relative spacing.
    CASE (4)
       d1mach = b**(1-DIGITS(x))      ! the largest relative spacing.
    CASE (5)
       d1mach = LOG10(b)
    CASE default
       WRITE (*, fmt = 9000)
9000   FORMAT ('1error    1 in d1mach - i out of bounds')
       STOP
    END SELECT

  END FUNCTION  d1mach
