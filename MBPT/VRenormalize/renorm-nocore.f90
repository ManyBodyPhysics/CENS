!             Program block renorm-nocore.f90
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!             LAST UPGRADE : September 2007, Similarity transformed interaction
!             no energy dependence. Print out for many-body perturbation theory
!             no-core shell-model calculations and coupled-cluster codes in
!             both m-scheme and j-scheme. The interaction includes both Coulomb and (r_i-r_j)^2/2m
!             if A > 2, see Eq. (93) of Prog. Part. Nucl. Phys. vol 53 pages 419-500. 
!             The final interaction is tailored to a triangular Pauli operator only.
!             It contains also the RG approach in oscillator space of Bogner et al, Phys. Rev. C75
!             (2007) 061001, see Eq. (6), but modified to an oscillator basis.
!                        
!
!                    Begin setup of similarity transformation
!
SUBROUTINE  setup_nocore
  USE partial_waves
  USE constants
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  IMPLICIT NONE
  INTEGER :: loop, number_orbits, max_conf

  !     reserve space in memory for mesh point and h.o. wave function  arrays
  ALLOCATE ( ra (n_rel), wra (n_rel))
  ALLOCATE ( rnlr (n_rel, 0:lmax, 0:nmax) )
  !     set up quantum numbers nl and NL relative and cm system
  ALLOCATE(coulomb_relcom(0:nmax, 0:nmax, 0:lmax))
  !     set up Coulomb interaction in relative coordinates
  coulomb_relcom = 0.0_dp  
  !  The argonne v18 interaction contains the coulomb + electromagnetic corrections 
  !  explicitely
  SELECT CASE (physical_system) 
  CASE('nuclear_physics')
     IF ( type_of_pot /= 'argonnev18') THEN
        IF ( coulomb_included =='coulomb') CALL coulomb_integral
     ENDIF
  CASE('atomic_physics')
     CALL electroncoulomb_integral
  END SELECT

  CALL find_spdata_relcm(number_orbits)
  relcm_sp_data%max_value=number_orbits
  CALL allocate_relcm_array(relcm_sp_data, number_orbits) 
  CALL make_configurations_relcm
  SELECT CASE (physical_system) 
  CASE('nuclear_physics')
     !     set up mesh points for relative coordinates
     CALL nocore_mesh                   
     !     setup h.o. wave functions, only relative coordinates needed
     CALL nocoreho_wfunction                
  END SELECT
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
  !     used to diagonalize the deuteron with an A-dependent CoM correction
  !     and obtain a model space effective interaction using a similarity
  !     transformation. This is done by the function setup_com_interaction.
  !     Note that the interaction does not depend on the CoM coordinates
  !     and thus the problem can be separated out. The final effective 
  !     interaction is stored for transformation to the lab system.
  !
  CALL setup_com_interaction
  !     With the final interaction in a harmonic oscillator basis defined by
  !      < nlNLJT_ZS | Veff | n'l'N'LJT_ZS > we perform the final transformation
  !     to the lab system in the function final_veff_labsystem
  CALL final_veff_labsystem
  DEALLOCATE(coulomb_relcom)
  DEALLOCATE ( ra, wra) ; DEALLOCATE ( rnlr) 
  DEALLOCATE (rel_conf%rel_ab) 
  DEALLOCATE (rel_conf%nconfs_rel ) ;   DEALLOCATE (rel_conf%nconfs_relmodel ) 
  CALL deallocate_relcm_array(relcm_sp_data)
  DEALLOCATE (relcm_conf%relcm_ab) 
  DEALLOCATE (relcm_conf%nconfs_relcm ) 

END SUBROUTINE setup_nocore
!
!           Set up the effective interaction in the lab-frame using the
!           mtx elements in the rel-CoM frame
!
SUBROUTINE final_veff_labsystem  
  USE single_particle_orbits
  USE configurations
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs 
  REAL(KIND=8), ALLOCATABLE :: lab_to_relcoeff(:,:)
  INTEGER, ALLOCATABLE :: lab_to_relconf(:,:)  
  INTEGER, ALLOCATABLE :: lab_to_relnumber(:)
  INTEGER ::  p_parity, ang_mom, isospin_z, max_coeff, pq_confs
  COMPLEX*16, ALLOCATABLE::  gna(:,:)
  REAL(KIND=8) :: dij
  REAL(KIND=8), ALLOCATABLE:: twobody_com(:,:), twobody_r2(:,:), twobody_p2(:,:)

  !     loop over isospin projection
  DO isospin_z=itzmin,itzmax 
     !     loop over parity values, here positive parity is 0, negative         
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=j_lab_min,j_lab_max
           !     find all possible configurations, large space and model space 
           CALL number_nocore_confs &
                (ang_mom,p_parity,isospin_z,gmatrix_configs)
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
           pq_confs=gmatrix_configs%number_confs
           ALLOCATE(gmatrix_configs%config_ab(pq_confs+pq_confs) )
           CALL setup_nocore_configurations &
                (ang_mom,p_parity,isospin_z,gmatrix_configs)
           !     find max possible number of transformation coeffs
           CALL find_ncoeffs(isospin_z,ang_mom,gmatrix_configs,max_coeff)           
           !     allocate space for various arrays needed in transform rel-cm -> lab 
           ALLOCATE(lab_to_relconf(pq_confs,max_coeff), &
                lab_to_relcoeff(pq_confs,max_coeff))          
           ALLOCATE(twobody_p2(pq_confs, pq_confs))
           ALLOCATE(twobody_r2(pq_confs, pq_confs))
           ALLOCATE(twobody_com(pq_confs, pq_confs))
           ALLOCATE(lab_to_relnumber(pq_confs))
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
           CALL nocore_free(gna,twobody_com,twobody_r2, twobody_p2,lab_to_relcoeff, lab_to_relconf,&
                lab_to_relnumber,max_coeff,gmatrix_configs,isospin_z)
           CALL veff_print(isospin_z,p_parity,ang_mom,gna,twobody_com,twobody_r2, twobody_p2,gmatrix_configs)
           !      free space
           DEALLOCATE(gmatrix_configs%config_ab)
           DEALLOCATE(twobody_com)
           DEALLOCATE(twobody_r2,twobody_p2)
           DEALLOCATE(lab_to_relcoeff, lab_to_relconf )
           DEALLOCATE (lab_to_relnumber); DEALLOCATE(gna)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE final_veff_labsystem
!
!                    Print out veff-mtx
!
SUBROUTINE veff_print(it,ip,ij,gna,twobody_com,twobody_r2, twobody_p2,gmatrix_configs)
  USE constants
  USE configurations
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, ijd, it, ip, ij, ia, ib, ic,id
  COMPLEX*16, DIMENSION(gmatrix_configs%number_confs, &
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
        SELECT CASE (type_of_renormv) 
        CASE('no-core')
           WRITE(7,1001) it, ip, ijd, ia, ib, ic, id, REAL(gna(j,i)), &
                twobody_com(j,i), twobody_r2(j,i), twobody_p2(j,i)       
        CASE('v-nrg')
           WRITE(7,1001) it, ip, ijd, ia, ib, ic, id, & 
                REAL(gna(j,i)), &
                twobody_com(j,i), twobody_r2(j,i), twobody_p2(j,i)       
        END SELECT
     ENDDO
  ENDDO
1001 FORMAT(7I4,5X,4(5X,E12.6))
END SUBROUTINE veff_print
!
!          Find veff-mtx and make transformation from CoM and
!          Rel system to lab system
!
SUBROUTINE nocore_free(gfree,twobody_com,twobody_r2, twobody_p2,lab_to_relcoeff, lab_to_relconf,&
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
  COMPLEX*16, DIMENSION(gmatrix_configs%number_confs, &
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

END SUBROUTINE nocore_free
!
!           Set up the interaction in the relative system using
!           a similarity transformation or flow renormalization group approach. 
!           Quantum numbers are 
!           < nl JT_ZS | Veff | n'l'JT_ZS >
!
SUBROUTINE setup_com_interaction
  USE single_particle_orbits
  USE configurations
  USE constants
  USE relcm_gmatrix
  USE partial_waves
  IMPLICIT NONE
  INTEGER ::  p_confs, pq_confs, loop, mspace, bra, ket
  COMPLEX*16, ALLOCATABLE :: cvec_pp(:,:)
  COMPLEX*16, ALLOCATABLE :: cvec_qp(:,:), heff(:,:)
  COMPLEX*16, ALLOCATABLE::  ceig_p(:) 
  REAL(KIND=8), ALLOCATABLE :: v2body(:,:), kinetic_energy(:,:)
  COMPLEX*16, ALLOCATABLE :: gna(:,:)
  DO loop=1,no_channels   
     pq_confs = rel_conf%nconfs_rel(loop)
     p_confs = rel_conf%nconfs_relmodel(loop)
     WRITE(6,*) ' Total states and model states for channel: ', loop, pq_confs, p_confs
     IF ( pq_confs <= 0 ) CYCLE
     IF ( ( p_confs == 0) .OR. ( p_confs > pq_confs)) CYCLE
     ALLOCATE(v2body(pq_confs, pq_confs)); ALLOCATE(kinetic_energy(pq_confs, pq_confs))
     v2body = 0.D0; kinetic_energy = 0.D0
     !     Setup the interaction to diagonalize in the big space
     !     It returns also the kinetic energy
     CALL nocore_channel(loop,v2body,kinetic_energy)
     IF ((type_of_pot == 'OPEP').OR.(type_of_pot == 'Tensorinteraction')   &
          .OR.(type_of_pot == 'LSinteraction') ) THEN
        DO ket = 1, p_confs
           DO bra = 1, p_confs
              v_com_rel(loop,bra,ket) = v2body(bra,ket)
           ENDDO
        ENDDO
     ELSE
        SELECT CASE (type_of_renormv) 
        CASE('no-core')
           ALLOCATE(cvec_pp(p_confs, p_confs)); ALLOCATE(cvec_qp(pq_confs-p_confs,p_confs))
           ALLOCATE(ceig_p(p_confs)); ALLOCATE(heff(p_confs, p_confs))
           ALLOCATE(gna(pq_confs, pq_confs))
           heff = 0.D0; gna = 0.0D0
           gna = kinetic_energy+v2body
           !      Get model space eigenvalues, eigenvectors of 
           !      model space and excluded space and
           !      model space configurations which match the corresponding ones
           !      of the large space
           CALL eigenvalues_large(cvec_pp,cvec_qp,ceig_p,gna,pq_confs,p_confs)
           !      setup 2p-effective interaction in P-space using
           !      the Lee-Suzuki similarity transformation
           CALL lee_suzuki( cvec_pp, cvec_qp, ceig_p, p_confs, pq_confs-p_confs, heff)
           DO ket = 1, p_confs
              DO bra = 1, p_confs
                 v_com_rel(loop,bra,ket) = heff(bra,ket)-kinetic_energy(bra,ket)
              ENDDO
           ENDDO
           DEALLOCATE(heff); DEALLOCATE(cvec_pp); DEALLOCATE(cvec_qp); DEALLOCATE(ceig_p); DEALLOCATE(gna)
        CASE('v-nrg')
           ALLOCATE(heff(pq_confs, pq_confs))
           heff = 0.D0
           CALL vnrgk_mtx(pq_confs,v2body,heff,kinetic_energy,loop)
           DO ket = 1, p_confs
              DO bra = 1, p_confs
                 v_com_rel(loop,bra,ket) = heff(bra,ket)
              ENDDO
           ENDDO
           DEALLOCATE(heff)
        END SELECT
     ENDIF
     !      free space
     DEALLOCATE(v2body, kinetic_energy)
  ENDDO

END SUBROUTINE setup_com_interaction
!
!     Obtain the bare interaction in a harmonic oscillator 
!     basis plus the kinetic energy and the Coulomb part. It contains
!     also the CoM correction to the relative coordinates. The latter depends
!     on the mass number of the nucleus
! 
SUBROUTINE nocore_channel(i,vint,kinetic_energy)
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
  REAL(KIND=8) ::  e_kin, e_r2, e_coulomb, vsum
  REAL(KIND=8), ALLOCATABLE :: vkk(:,:)
  REAL(KIND=8),  INTENT(INOUT):: &
       vint(rel_conf%nconfs_rel(i),rel_conf%nconfs_rel(i)), &
       kinetic_energy(rel_conf%nconfs_rel(i),rel_conf%nconfs_rel(i))

  ncoup=1
  IF ( orb_lrel_max(i) == jang_rel(i) ) ncoup = 1
  IF ( orb_lrel_max(i) /= jang_rel(i) ) ncoup = 2
  ALLOCATE(vkk (ncoup*n_rel,ncoup*n_rel))
  jang=jang_rel(i)
  !     setup the NN interaction
  vkk = 0.; vint = 0; kinetic_energy = 0.

  SELECT CASE (physical_system) 
  CASE('nuclear_physics')
     CALL nocorepotential_interface(ncoup,vkk,i)
  END SELECT
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
        vsum = 0.; e_kin =0.; e_r2 = 0.
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
        CALL nocore_hosc(rnlr(:,lb,np), rnlr(:,la,n),vkk(lim1:lim2,lim3:lim4),bra,ket,vsum)
        !       Here we set up the Coulomb part and the h.o. energy in an
        !       oscillator basis. We have only diagonal terms for the h.o. energy
        e_coulomb = 0.
        SELECT CASE (physical_system) 
        CASE('nuclear_physics')
           IF ( ( la == lb ).AND.( iso(i) == -1 ))  THEN
              e_coulomb = coulomb_relcom(n,np,la) 
           ENDIF
        CASE('atomic_physics')
           IF (la == lb )  THEN
              e_coulomb = coulomb_relcom(n,np,la) 
           ENDIF
        END SELECT
        ! h.o.  energy from relative motion with CoM correction, note not diagonal
        ! in n and np
        IF ( ( la == lb )) THEN
           IF ( n == np-1) THEN
              e_r2 = -SQRT(np*(np+lb+0.5)) ! harmonic oscillator potential
              e_kin = SQRT(np*(np+lb+0.5)) ! harmonic oscillator kinetic
           ENDIF
           IF ( n == np) THEN 
              e_kin =  2*np+lb+1.5 ! harmonic oscillator kinetic
              e_r2 =  2*np+lb+1.5  ! harmonic oscillator potential
           ENDIF
           IF ( n == np+1) THEN
              e_r2 = -SQRT((np+1.)*(np+lb+1.5)) ! harmonic oscillator potential
              e_kin = SQRT((np+1.)*(np+lb+1.5)) ! harmonic oscillator kinetic
           ENDIF
        ENDIF
        ! for the ho case listed here it contains kinetic energy +
        ! ho potential energy in relative coordinates
        ! this means that you must use HO sp energies in SM studies.
        SELECT CASE (physical_system) 
        CASE('nuclear_physics')
           IF ( mass_nucleus == 0)  THEN 
              kinetic_energy(ket,bra) = 0.5*hbar_omega*(e_kin+ e_r2)
              vint(ket,bra)=vsum+e_coulomb
           ELSE 
              kinetic_energy(ket,bra) = 0.5*hbar_omega*(e_kin+ e_r2)
              vint(ket,bra)=vsum-hbar_omega*e_r2/mass_nucleus+e_coulomb
           ENDIF
           !       We use Kvaal's notation for the Hamiltonian, see his thesis, UiO 2008
        CASE('atomic_physics')
           kinetic_energy(ket,bra) = 0.5_dp*e_kin+0.5_dp*e_r2
           vint(ket,bra)=SQRT(0.5_dp)*e_coulomb
        END SELECT
     ENDDO
  ENDDO

END SUBROUTINE nocore_channel
!
!          Obtain the bare potential in oscillator basis
!          Final V is in units of MeV or eV
!          Time consuming part. Parallel version available.
!
SUBROUTINE nocore_hosc(wave_bra, wave_ket, a,bra,ket,vsum)
  USE wave_functions
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  REAL(KIND=8), INTENT(INOUT) :: vsum
  INTEGER, INTENT(IN) :: bra, ket
  INTEGER :: i, j
  REAL(KIND=8), DIMENSION(n_rel,n_rel), INTENT(IN) :: a
  REAL(KIND=8) :: sum1, sum2, hbarc3
  REAL(KIND=8), DIMENSION(n_rel), INTENT(IN) :: wave_bra, wave_ket 

  hbarc3=hbarc**3
  sum1=0.
  DO i=1, n_rel
     sum2=0.
     DO j=1, n_rel
        sum2=sum2+wave_ket(j)*a(j,i)
     ENDDO
     sum1=sum1+sum2*wave_bra(i)
  ENDDO
  vsum = sum1*hbarc**3

END SUBROUTINE nocore_hosc
!
!                 Set up h.o. wf for rel system
!
SUBROUTINE nocoreho_wfunction
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER :: n, l, i, j
  REAL(KIND = 8)  :: cx(0:200), factor, z_rel, xp, ph, oscl_r, sum_rel, contrib

  oscl_r=oscl*SQRT(2.)        ! Oscillator parameter for relative
  DO n=0,nmax
     ph=(-1.D0)**n
     DO l=0,lmax
        sum_rel=0.
        factor = 0.5D0*((n+l+2)*LOG(2.D0)+fac(n)-dfac(2*n+2*l+1)-0.5D0*LOG(pi))
        factor = EXP(factor)
        DO i=1,n_rel
           z_rel= ra(i)*ra(i)*oscl_r*oscl_r
           CALL laguerre_general( n, l+0.5D0, z_rel, cx )
           xp = EXP(-z_rel*0.5)*((ra(i)*oscl_r)**l)*cx(n)
           rnlr(i,l,n) = xp*(wra(i)*(ra(i)**2))*ph*factor*(oscl_r**(1.5D0))     ! rel wf
           contrib = xp*ph*factor*(oscl_r**(1.5D0)) 
           sum_rel=sum_rel+ wra(i)*(contrib*ra(i))**2
        ENDDO
        WRITE(6,*) 'Norm rel ho wf n,l : ', n, l, sum_rel
     ENDDO
  ENDDO

END SUBROUTINE nocoreho_wfunction
!
!Fast diagonalization of real(kind=8) symmetric matrices
!
SUBROUTINE diag_exact(cvec, ceig, n, h)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: h
  INTEGER :: k, np
  INTEGER :: p, i,j, i1, kvec, lwork, option, il ,iu, info
  REAL(KIND=8) :: a(n*(n+1)), work(30*n), thresh, vl, vu
  LOGICAL :: flag(n)  
  COMPLEX*16 :: cvu, cvl
  COMPLEX*16, INTENT(INOUT) :: ceig(n), cvec(n,n)

  cvec = 0.; ceig = 0.
  lwork = 30*n
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
     ENDDO
  ENDDO
  CALL cs(n,a,ceig,kvec,cvec,lwork,work,thresh, option,il,iu,cvl,cvu,flag,info)

END SUBROUTINE diag_exact
!
! Get eigenvalues for large space, eigenvectors for model space and excluded space
!
SUBROUTINE eigenvalues_large(cvec_pp,cvec_qp,ceig_p,hamilton,n,np)
  USE single_particle_orbits
  USE configurations
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, np
  INTEGER :: k1, i, nqq
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: cvec
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: ceig
  COMPLEX*16, INTENT(INOUT)  :: ceig_p(np)
  COMPLEX*16, INTENT(INOUT) :: cvec_pp(np,np)
  COMPLEX*16, INTENT(INOUT) :: cvec_qp(n-np, np) 
  COMPLEX*16, INTENT(IN) :: hamilton(n,n)
  REAL(KIND=8) :: cvec_max
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: temp

  ! diagonalize 2p-effective shell model hamiltonian, this takes a symmetric 
  ! matrix as input.
  ALLOCATE(cvec(n,n)); ALLOCATE(ceig(n))
  cvec = 0.; ceig = 0.
  CALL diag_exact(cvec, ceig, n, hamilton)
  cvec_pp = 0.; cvec_qp = 0.
  ! Set up eigenvalues and eigenvectors for model space and excluded space
  DO k1=1, np
     ! loop over all model space coefficients of exact eigenvectors |k>
     nqq = 0
     DO i = 1, n
        IF ( i <= np ) THEN
           cvec_pp(i,k1) = cvec(i,k1) 
           ceig_p(k1) = ceig(k1) 
        ELSEIF ( i > np ) THEN
           nqq = nqq + 1
           cvec_qp(i-np,k1) = cvec(i,k1) 
        ENDIF
     ENDDO
  ENDDO
  DEALLOCATE(ceig); DEALLOCATE(cvec)

END SUBROUTINE eigenvalues_large
!
! Lee Suzuki similarity transformation 
!
SUBROUTINE lee_suzuki( cvec_pp, cvec_qp, ceig, np, nq, heff)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: np, nq
  COMPLEX*16, DIMENSION(np,np), INTENT(IN) :: cvec_pp
  COMPLEX*16, DIMENSION(nq,np), INTENT(IN)  :: cvec_qp
  COMPLEX*16, DIMENSION(np), INTENT(IN) :: ceig
  COMPLEX*16, DIMENSION(np,np), INTENT(INOUT) ::  heff
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: omega2, omega_2inv, u, u_inv, cvec_pp_inv, &
       eigen_vec, vl, heff_rhs, temp1, omega, temp2
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: ceig_p, omega2_eig
  REAL(KIND=8), DIMENSION(2*np) :: rwork
  COMPLEX*16, DIMENSION(10000) :: work1
  COMPLEX*16 :: d, sum1, norm
  INTEGER :: i_p,j_p, j_q, ii, jj, k, i_q , a_p, a_pp 
  INTEGER :: i, lda, ldb, ldvl, ldvr, info, lwork, ilo , ihi
  CHARACTER*1 :: jobvl, jobvr, balanc, sense
  REAL(KIND=8), DIMENSION(np) :: scale, rconde, rcondv
  REAL(KIND=8) :: abnrm, a1, a2, b1, b2
  INTEGER :: ipiv(np) ,j

  balanc = 'n';  jobvl = 'n' ;  jobvr = 'v';  sense = 'n';  lda = np
  ldvl = 1;  ldvr = np;  lwork = 10000
  ALLOCATE(omega2(np,np), omega_2inv(np,np),cvec_pp_inv(np,np)) 
  ALLOCATE(eigen_vec(np,np), vl(np,np), temp1(np,np)) 
  ALLOCATE(omega(nq,np),temp2(np,nq)); ALLOCATE(ceig_p(np), omega2_eig(np))
  eigen_vec = (0.D0,0.D0) 
  temp1 = TRANSPOSE(cvec_pp) 
  CALL zgeev( jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
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
  heff = (0.0D0,0.0D0)
  DO i_p = 1, np
     DO j_p = 1, np 
        heff(i_p,j_p) = SUM( cvec_pp(i_p,:)*ceig(:)*cvec_pp(j_p,:) ) 
        sum1 = (0.D0,0.D0)
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
  ALLOCATE(u(np,np), u_inv(np,np),heff_rhs(np,np)) 
  u = (0.D0,0.D0); u_inv =  (0.D0,0.D0)
  CALL sqrtmat( omega2, u, u_inv ,np ) 
  heff_rhs =  MATMUL( u, MATMUL( heff,u_inv ) ) 
  DEALLOCATE(u, u_inv) 
  heff = (0.D0,0.D0)
  heff = heff_rhs
  ! diagonalize 2p-effective shell model hamiltonian
  CALL zgeev( jobvl, jobvr, np, heff_rhs, lda, ceig_p, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  ! compare spectrum from exact and P-space diagonalization
  WRITE(6,*) 'Compare model space two-body spectrum with exact spectrum:' 
  DO i_p = 1, np
     a1 = REAL( ceig_p(i_p))
     a2 = AIMAG( ceig_p(i_p))
     b1 = REAL( ceig(i_p) )
     b2 = AIMAG( ceig(i_p) ) 
     WRITE(6,*) A1, A2, B1, B2
  ENDDO
  DEALLOCATE(omega2, omega_2inv, cvec_pp_inv) 
  DEALLOCATE(eigen_vec, vl, heff_rhs, temp1); DEALLOCATE(omega) 
  DEALLOCATE(temp2); DEALLOCATE(ceig_p, omega2_eig)

END SUBROUTINE lee_suzuki
!
!
!          This routine contains old-fashioned common blocks
!          acting as interface between the effective interaction program
!          and the potential routines of the Bonn type, Nijmegen or 
!          Argonne groups. 
!          vkk is in units of MeV^-2,  
!          ra are mesh points in rel coordinates, 
!          units of fm^-1
!          in partial wave basis, v(6) 
!          means
!          v(1)= singlet uncoupled
!          V(2)= triplet uncoupled
!          V(3)= triplet coupled, < J+1 | V | J+1>
!          V(4)= triplet coupled, < J-1 | V | J-1>
!          V(5)= triplet coupled, < J+1 | V | J-1>
!          V(6)= triplet coupled, < J-1 | V | J+1>
!
SUBROUTINE nocorepotential_interface(ncoup,vkk,ichan)
  USE wave_functions
  USE relcm_gmatrix
  USE partial_waves
  USE constants
  USE idaho_chiral_potential
  IMPLICIT NONE
  INTEGER :: i,j, ncoup, k, l, jxx, ichan, inn, isospin_tz, spin, &
       n1, ix, iy
  REAL(KIND=8)  :: v,xmev,ymev, c, q
  CHARACTER (LEN=4) :: label
  COMMON /cnn/ inn
  COMMON /cpot/ v(6),xmev,ymev
  common /cpts/   q(97),c,n1,ix,iy
  COMMON /cstate/ jxx, heform, sing, trip, coup, endep, label
  ! n2lo500-pounders read/write stuff
  common /crdwrt/ kread,kwrite
  INTEGER :: kread, kwrite
  !
  ! chp stuff
  !
  INTEGER :: tz, T
  !
  LOGICAL :: sing, trip, coup, heform, endep
  REAL(KIND=8), DIMENSION(ncoup*n_rel,ncoup*n_rel), INTENT(OUT) :: vkk
  REAL(KIND=8) :: v00, v11, v12, v21, v22, fmult, besl, bs, b11, b12, b21, b22
  REAL(KIND=8) :: xb(500), wb(500), v8(500,2,2), temp(500,2,2)
  INTEGER :: j1,l1,l2,s1,isot,m,nt,tz1,tz2,lpot

  kread=5
  kwrite=6

  heform=.FALSE.
  jxx=jang_rel(ichan)
  sing=spin_rel(ichan) == 0
  trip=spin_rel(ichan) == 1
  coup=(ncoup == 2) 
  isospin_tz=iso(ichan)
  spin = spin_rel(ichan)
  SELECT CASE (type_of_pot)
  CASE('Idaho-A')
     inn = 1
  CASE('Idaho-B')
     inn = 2
  CASE('n3lo')
     SELECT CASE ( isospin_tz)
     CASE (-1)
        IF ( (csb == 'csb').AND.(cib == 'cib'))       inn = 1    !  pp case = 1
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    !  pp case = nn case
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    !  pp case = pn case
        tz1 = 1; tz2 = 1
     CASE (0)
        inn = 2    !  pn case = 2  if all inn=2, no CSB or ISB
        tz1 = 1; tz2 = -1
     CASE ( 1)
        IF ( (csb == 'csb').AND.(cib == 'cib'))       inn = 3    !  nn case = 3
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    !  nn case = pp case 
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    !  nn case = pn case
        tz1 = -1; tz2 = -1
     END SELECT
  CASE('chiral-potential')
     SELECT CASE ( isospin_tz)
     CASE (-1)
        T=1
        IF ( (csb == 'csb').AND.(cib == 'cib'))       tz = -1    !  pp case = 1
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    tz = +1    !  pp case = nn case
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) tz =  0    !  pp case = pn case
     CASE (0)
        IF ( ncoup == 1) THEN
           l = jxx
        ELSEIF ( ncoup == 2) THEN
           l = jxx-1
        END IF
        IF (MOD(l+spin+1,2) == 0) T=0
        IF (MOD(l+spin+0,2) == 0) T=1
        tz = 0    !  pn case = 2  if all inn=2, no CSB or ISB
     CASE ( 1)
        T=1
        IF ( (csb == 'csb').AND.(cib == 'cib'))       tz = +1    !  nn case = 3
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    tz = -1    !  nn case = pp case 
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) tz =  0    !  nn case = pn case
     END SELECT
  CASE('n3lo3b')
     SELECT CASE ( isospin_tz)
     CASE (-1)
        IF ( (csb == 'csb').AND.(cib == 'cib'))       inn = 1    !  pp case = 1
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    !  pp case = nn case
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    !  pp case = pn case
        tz1 = 1; tz2 = 1
     CASE (0)
        inn = 2    !  pn case = 2  if all inn=2, no CSB or ISB
        tz1 = 1; tz2 = -1
     CASE ( 1)
        IF ( (csb == 'csb').AND.(cib == 'cib'))       inn = 3    !  nn case = 3
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    !  nn case = pp case 
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    !  nn case = pn case
        tz1 = -1; tz2 = -1
     END SELECT
  CASE('CD-bonn')
     SELECT CASE ( isospin_tz)
     CASE (-1)
        IF ( (csb == 'csb').AND.(cib == 'cib'))       inn = 1    
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    
        tz1 = 1; tz2 = 1
     CASE (0)
        inn = 2    !  pn case = 2  if all inn=2, no CSB or ISB
        tz1 = 1; tz2 = -1
     CASE ( 1)
        IF ( (csb == 'csb').AND.(cib == 'cib'))       inn = 3    
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    
        tz1 = -1; tz2 = -1
     END SELECT
  END SELECT
  IF ( (type_of_pot == 'argonnev8').OR.(type_of_pot == 'argonnev18').OR.   &
       (type_of_pot == 'OPEP').OR.(type_of_pot == 'Tensorinteraction')   &
       .OR.(type_of_pot == 'LSinteraction').OR.(type_of_pot == 'MalflietTjon') ) THEN
     SELECT CASE ( isospin_tz)
     CASE (-1)
        fmult = 2.D0/ACOS(-1.D0)/p_mass(-1)/hbarc
        tz1 = 1; tz2 = 1
     CASE (0)
        tz1 = 1; tz2 = -1
        fmult = 2.D0/ACOS(-1.D0)/p_mass(0)/hbarc
     CASE ( 1)
        tz1 = -1; tz2 = -1
        fmult = 2.D0/ACOS(-1.D0)/p_mass(1)/hbarc
     END SELECT
     lpot = 1  ! switch for argonne potentials, lpot 1 means v18
     IF (type_of_pot == 'argonnev8') lpot = 2   ! lpot = 2 means v8
     IF (type_of_pot == 'OPEP') lpot = 3   ! lpot = 3 means OPEP from V18
     IF (type_of_pot == 'Tensorinteraction') lpot = 4   ! lpot = 4 means full tensor from V18
     IF (type_of_pot == 'LSinteraction') lpot = 5   ! lpot = 5 means full LS from V18
     IF (type_of_pot == 'MalflietTjon') lpot = 6   ! lpot = 6 Malfliet-Tjon
     !  sett mesh points in r-space
     nt = 300
     CALL gauss_legendre(0.D0,10.D0,xb,wb,nt)
     v8 = 0.
     j1 = jang_rel(ichan); s1 = spin_rel(ichan)
     IF ( ncoup == 1) THEN
        l1 = jang_rel(ichan)
        l2 = jang_rel(ichan)
     ELSEIF ( ncoup == 2) THEN
        l1 = jang_rel(ichan)-1
        l2 = jang_rel(ichan)+1
     ENDIF
     IF ( ABS(isospin_tz) == 1) THEN
        isot = 1
        CALL argonne_pots(j1,l1,s1,isot,tz1,tz2,v8,xb,nt,lpot)
     ELSEIF ( isospin_tz == 0) THEN
        v8 = 0.
        !  now check partial waves in good isospin formalism
        DO isot = 0, 1
           IF(MOD(l1+s1+isot,2) == 0) CYCLE
           temp = 0.
           CALL argonne_pots(j1,l1,s1,isot,tz1,tz2,temp,xb,nt,lpot)
           v8 = v8+temp
        ENDDO
     ENDIF
     DO i=1,n_rel
        IF (ncoup == 2)   k=i+n_rel 
        DO j=1,n_rel
           IF (ncoup == 2) l=j+n_rel 
           IF (ncoup == 1 ) THEN
              v00 = 0.
              DO m=1,nt
                 bs=besl(xb(m)*ra(i),l1)*besl(xb(m)*ra(j),l1)
                 v00=v00+bs*v8(m,1,1)*wb(m)*xb(m)**2
              ENDDO
              vkk(j,i) = v00*fmult
           ELSEIF (ncoup ==2  ) THEN 
              v11=0.d0
              v12=0.d0
              v21=0.d0
              v22=0.d0
              DO m=1,nt
                 c=wb(m)*xb(m)**2
                 b11=besl(xb(m)*ra(i),l1)
                 b12=besl(xb(m)*ra(i),l2)
                 b21=besl(xb(m)*ra(j),l1)
                 b22=besl(xb(m)*ra(j),l2)
                 v11=v11+c*b11*b21*v8(m,1,1)
                 v12=v12+c*b11*b22*v8(m,1,2)
                 v21=v21+c*b12*b21*v8(m,2,1)
                 v22=v22+c*b12*b22*v8(m,2,2)
              ENDDO
              vkk(j,i)=v11*fmult
              vkk(j,k)=v21*fmult
              vkk(l,i)=v12*fmult
              vkk(l,k)= v22*fmult
           ENDIF
        ENDDO
     ENDDO
  ELSE
     DO i=1,n_rel
        xmev=ra(i)*hbarc
        IF (ncoup == 2)   k=i+n_rel 
        DO j=1,n_rel
           ymev=ra(j)*hbarc
           IF (ncoup == 2) l=j+n_rel 
           SELECT CASE (type_of_pot)
           CASE('CD-bonn')
              CALL cdbonn
           CASE('Idaho-A')
              CALL idaho
           CASE('Idaho-B')
              CALL idaho
           CASE('chiral-potential')
              CALL chp(xmev,ymev,coup,spin,jxx,T,tz,v)
           CASE('n3lo')
              CALL idaho_n3lo
           CASE('n3lo3b')
              CALL n3lo3b
           END SELECT
           IF (sing ) THEN
              vkk(j,i)=v(1)
           ELSEIF ((trip).AND.(.NOT.coup )  ) THEN 
              vkk(j,i) = v(2)
           ELSEIF (coup ) THEN
              vkk(j,i)=v(4)
              vkk(j,k)=v(5)
              vkk(l,i)=v(6)
              vkk(l,k)= v(3)
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  
END SUBROUTINE nocorepotential_interface
!
!       Compute the renormalized interaction in oscillator space
!
SUBROUTINE vnrgk_mtx(nconfs,vzz,heff, kinetic_energy,ichan)
  USE wave_functions
  USE partial_waves
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nconfs, ichan
  REAL(KIND=8), DIMENSION(nconfs,nconfs), INTENT(IN) :: vzz, kinetic_energy
  REAL(KIND=8), DIMENSION(nconfs,nconfs), INTENT(INOUT) :: heff
  INTEGER(KIND=4) :: i, j, n_ode, ij, iflag, i1, j1, isospin_tz
  INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: iwork
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: rg_vec, work
  REAL(KIND=8) :: relerr, abserr, k_lambda, start_point
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: cvec, temp
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) ::  ceig

  INTERFACE
     SUBROUTINE oscillator_derivative(k_lambda, v, dv)
       USE wave_functions
       USE constants
       INTEGER ::  ij, i, j, i1, j2, kk, k
       REAL(DP) :: sum, k_lambda, v(:), dv(:), k1, k2, p
       REAL(DP), ALLOCATABLE :: vij(:,:)
     END SUBROUTINE oscillator_derivative
  END INTERFACE

  isospin_tz=iso(ichan)
  relerr = 1.0e-8_dp; abserr =  1.0e-8_dp 
  !  cutoff in MeV
  !   oscillator energy of model space limit
  k_lambda  = (nlmax_model+1.5_dp)*hbar_omega
  ! dimension of vectors and matrices
  n_total = nconfs; n_ode = n_total*(n_total+1)/2
  ! total dim used by derivatives function
  ALLOCATE(rg_vec(n_ode)); ALLOCATE( work(100+21*n_ode)); ALLOCATE (iwork(5))
  ALLOCATE ( ho_onebodyenergy(n_total))
  rg_vec = 0.0_dp; work = 0.0_dp;  iwork = 0
  !  Transform the matrix -V- to a one-dim vector of dim ntot*(ntot+1)/2
  !  which is also the number of ordinary differential equations
  ij = 0  
  DO i = 1, n_total
     ho_onebodyenergy(i) = kinetic_energy(i,i)
     DO j = i, n_total
        ij = ij + 1
        rg_vec(ij) = vzz(j,i)
     ENDDO
  ENDDO
  iflag = 1
  !   oscillator energy of max space, nlmas typically around 200-300 (2n+l)

  start_point = (nlmax+1.5_dp)*hbar_omega
  CALL ode(oscillator_derivative,n_ode,rg_vec,start_point,k_lambda,relerr,abserr,iflag,work,iwork,isospin_tz)
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
  temp = kinetic_energy + heff
  !  Diagonalize  to test  if one gets the deuteron
  CALL vlowkdiag_exact( temp, cvec, ceig, n_total )
  DO i = 1, n_total
     WRITE(6,*) i, REAL(ceig(i))
  ENDDO
  DEALLOCATE(ceig); DEALLOCATE(cvec, temp); DEALLOCATE(ho_onebodyenergy)
  DO i = 1, n_total
     write(6,*) i, heff(i,i), vzz(i,i)
  ENDDO

END SUBROUTINE vnrgk_mtx

SUBROUTINE oscillator_derivative(k_lambda, v, dv, isotz)
  USE wave_functions
  USE constants
  IMPLICIT NONE
  INTEGER ::  ij, i, j, i1, j2, kk, k, i2, isotz
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
     k1 = ho_onebodyenergy(i) 
     DO j = i, n_total
        k2 = ho_onebodyenergy(j) 
        ij = ij + 1
        sum = 0_dp
        DO k = 1, n_total
           p = ho_onebodyenergy(k) 
           sum = sum +  (k1+k2-2.0_dp*p)*vij(j,k)*vij(k,i) 
        ENDDO
        dv(ij) =  sum - (k2-k1)*(k2-k1)*vij(j,i)
     ENDDO
  ENDDO
  dv = -2.0_dp*dv/(k_lambda**3)
  DEALLOCATE(vij)

END SUBROUTINE oscillator_derivative
!
!                 Coulomb interaction in coordinate space for electrons
!                 Strictly in atomic units, with the potential by a strength factor
SUBROUTINE electroncoulomb_integral
  USE constants
  USE wave_functions
  USE relcm_gmatrix
  IMPLICIT NONE
  INTEGER :: n1, n2, lr, i, nrel_max
  REAL(KIND=8) :: int_sum, xr, xp, factor1, factor2
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: rr, wrr
  REAL(KIND=8) :: cx(0:200)

  ! Oscillator parameter for relative
  nrel_max = 1000
  ALLOCATE ( rr(nrel_max ), wrr(nrel_max ))
  CALL gauss_legendre(0.d0, 50.d0, rr, wrr, nrel_max )
  DO lr = 0, lmax, 1
     DO n1 = 0, nmax, 1
        factor1 = 0.5D0*((n1+lr+1)*LOG(2.D0)+fac(n1)-dfac(2*n1+2*lr+1)-0.5D0*LOG(pi))
        factor1 = EXP(factor1)
        DO n2 = 0, nmax,1
           IF ( n2 == n1) THEN
              factor2 = factor1
           ELSE
              factor2 = 0.5D0*((n2+lr+1)*LOG(2.D0)+fac(n2)-dfac(2*n2+2*lr+1)-0.5D0*LOG(pi))
              factor2 = EXP(factor2)
           ENDIF
           int_sum = 0.D0
           DO i=1,nrel_max
              CALL laguerre_general( n1, lr+0.5D0, rr(i), cx )
              xp = cx(n1)*exp(-rr(i)*0.5)*(rr(i)**(lr*0.5))
              IF ( n1 == n2) THEN
                 xr = xp
              ELSE
                 CALL laguerre_general( n2, lr+0.5D0, rr(i), cx )
                 xr = cx(n2)*exp(-rr(i)*0.5)*(rr(i)**(lr*0.5))
              ENDIF
              int_sum=int_sum+wrr(i)*xp*xr
           ENDDO
           ! Coulomb energy in units of hbar_omega
           coulomb_relcom(n2, n1, lr) = atomic_strength*int_sum*factor1*factor2
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE ( rr, wrr)

END SUBROUTINE electroncoulomb_integral
