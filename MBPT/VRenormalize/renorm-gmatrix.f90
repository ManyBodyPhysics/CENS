!             Program block renorm-gmatrix.f
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!
!                    Begin g-matrix code
!
SUBROUTINE  setup_g_matrix
  USE partial_waves
  USE constants
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  USE jobs_mpi
  IMPLICIT NONE
  INTEGER :: loop, number_orbits, max_conf

  my_rank_mpi=0 
  !     reserve space in memory for various arrays
  ALLOCATE ( ra (n_rel), wra (n_rel));  ALLOCATE ( rgkk (n_rel), wgkk (n_rel))
  ALLOCATE ( rca (n_cm_mom), wrca (n_cm_mom))
  ALLOCATE ( hol (n_rel, 0:lmax, 0:nmax) )
  ALLOCATE ( chol (n_rel, 0:lmax, 0:nmax) )
  ALLOCATE ( rnlr (n_rel, 0:lmax, 0:nmax) )
  ALLOCATE ( rnlc (n_cm_mom, 0:lmax, 0:nmax))
  ALLOCATE(coulomb_relcom(0:nmax, 0:nmax, 0:lmax))
  !     set up Coulomb interaction in relative coordinates
  coulomb_relcom = 0.0_dp  
  !  The argonne v18 interaction contains the coulomb + electromagnetic corrections 
  !  explicitely
  IF ( type_of_pot /= 'argonnev18') THEN
     IF ( coulomb_included =='coulomb') CALL coulomb_integral
  ENDIF
  !     set up quantum numbers for nl and NL relative and cm system
  CALL find_spdata_relcm(number_orbits)
  relcm_sp_data%max_value=number_orbits
  CALL allocate_relcm_array(relcm_sp_data, number_orbits) 
  CALL make_configurations_relcm
  !     set up mesh points in lab frame
  CALL rel_mesh                   
  !     set up mesh points in cm frame 
  CALL cm_mesh                    
  !     setup ho wave functions
  CALL ho_wfunction                
  !     find all partial waves with given jmin jmax
  CALL setup_channels          
  !
  !     Find all possible configurations in the relative and cm system
  !     based on h.o. quantum number nlNL for all partial waves
  !     and allocate array for the given partial wave
  !     all configs are set up in the call to setup_configurations_cm
  !     First we search for the number of configurations in each channel
  !
  ALLOCATE ( relcm_conf%nconfs_relcm(no_channels) ) 
  CALL number_configurations_relcm(relcm_conf)
  max_conf=0
  length_config_all = 0
  DO loop = 1, no_channels
     IF ( max_conf <  relcm_conf%nconfs_relcm(loop) ) THEN
        max_conf= relcm_conf%nconfs_relcm(loop)
     ENDIF
     length_config_all = length_config_all +relcm_conf%nconfs_relcm(loop) * & 
          relcm_conf%nconfs_relcm(loop)
  ENDDO
  WRITE(6,*) 'Total number of relative and com frame matrix elements', length_config_all
  ALLOCATE( relcm_conf%relcm_ab(no_channels, max_conf+max_conf))
  CALL setup_configurations_relcm(relcm_conf)      
  !     propagator in k-space for g-matrix without pauli operator
  ALLOCATE ( propagator(n_rel,n_cm_mom,n_startenergy_g) )
  ALLOCATE ( osc_propagator(n_rel,n_cm_mom,n_startenergy_g) )
  ALLOCATE ( wave_cm(n_cm_mom)) 
  CALL k_space_propagator
  !
  !     allocate space for the free g-matrix in the rel-cm system
  !
  ALLOCATE ( gtf_tot_all ( 4, length_config_all, n_startenergy_g ) )
  gtf_tot_all = 0.0_dp
  DO loop=1,no_channels    ! loop over channels in rel&cm system  
     CALL g_channel(loop)
  ENDDO
  !     find total number of matrix elements in lab frame to be stored 
  CALL setupg
  !     deallocation of arrays
  DEALLOCATE ( gtf_tot_all)
  DEALLOCATE (relcm_conf%relcm_ab) 
  DEALLOCATE (relcm_conf%nconfs_relcm ) ; DEALLOCATE (osc_propagator)
  DEALLOCATE ( rnlc ) 
  DEALLOCATE ( wave_cm ); ;   DEALLOCATE ( rgkk, wgkk) ; DEALLOCATE ( rca, wrca )
  DEALLOCATE(coulomb_relcom)
  CALL deallocate_relcm_array(relcm_sp_data)
  DEALLOCATE ( ra, wra) 
  DEALLOCATE ( hol, rnlr); DEALLOCATE ( chol )

END SUBROUTINE setup_g_matrix
!
!           Set up the G-mtx in the lab-frame using the
!           mtx elements in rel-cm frame, NO bhf calculation
!
SUBROUTINE setupg
  USE single_particle_orbits
  USE configurations
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs
  REAL(DP), ALLOCATABLE :: lab_to_relcoeff(:,:)
  INTEGER, ALLOCATABLE :: lab_to_relconf(:,:)
  INTEGER, ALLOCATABLE :: lab_to_relnumber(:)
  INTEGER ::  p_parity, ang_mom, isospin_z, max_coeff, n_confs, ie, bra, ket
  REAL(DP), ALLOCATABLE :: gna(:,:,:), coulomb(:,:)
  REAL(DP), ALLOCATABLE:: twobody_com(:,:), twobody_r2(:,:), twobody_p2(:,:)
  REAL(DP), ALLOCATABLE :: gfree_r(:,:,:), gfree(:,:,:), gfree_e(:,:,:), &
       gfree_ee(:,:,:), pep(:,:,:), a(:,:)
  REAL(DP) :: d
  INTEGER :: ix, iy
  !     loop over isospin projection
  DO isospin_z=itzmin,itzmax 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=j_lab_min,j_lab_max
           !     find all possible configurations 
           CALL number_gmatrix_confs &
                (ang_mom,p_parity,isospin_z,gmatrix_configs)
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
           n_confs=gmatrix_configs%number_confs
           ALLOCATE(gmatrix_configs%config_ab(n_confs+n_confs) )
           CALL setup_gmatrix_configurations &
                (ang_mom,p_parity,isospin_z,gmatrix_configs)
           !     find max possible number of transformation coeffs
           CALL find_ncoeffs(isospin_z,ang_mom,gmatrix_configs,max_coeff)
           !     allocate space in heap for various arrays needed in transform rel-cm -> lab 
           ALLOCATE(lab_to_relconf(n_confs,max_coeff), &
                lab_to_relcoeff(n_confs,max_coeff))
           ALLOCATE(lab_to_relnumber(n_confs))
           !     setup transformation coefficients for oscillator basis
           !     transformations from the c.m. frame to the lab frame
           CALL mosh_transf(isospin_z,ang_mom,gmatrix_configs, &
                lab_to_relcoeff, lab_to_relconf, lab_to_relnumber, &
                max_coeff)
           ALLOCATE(gfree_e(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(gfree_ee(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(gfree_r(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(gna(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(gfree(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(pep(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(coulomb(n_confs, n_confs))
           ALLOCATE(twobody_p2(n_confs, n_confs))
           ALLOCATE(twobody_r2(n_confs, n_confs))
           ALLOCATE(twobody_com(n_confs, n_confs))
           gfree_e=0.0_dp ; gfree=0.0_dp ; gna=0.0_dp; pep=0.0_dp
           gfree_r=0.0_dp; gfree_ee=0.0_dp; coulomb = 0.0_dp
           twobody_com = 0.D0; twobody_p2 = 0.0D0; twobody_r2 = 0.0D0
           !     find the g-matrix with no pauli operator
           CALL gmtx_free(gfree,gfree_e,gfree_r,gfree_ee,coulomb, &
                twobody_com,twobody_r2, twobody_p2, lab_to_relcoeff, lab_to_relconf,&
                lab_to_relnumber,max_coeff,gmatrix_configs,isospin_z)
           DEALLOCATE(lab_to_relcoeff, lab_to_relconf )
           DEALLOCATE (lab_to_relnumber)
           !     setup energy deno for part with pauli operator  
           CALL pspace_energy(isospin_z,ang_mom,pep,gmatrix_configs)
           ! build first 1/e +1/eG_F1/e and invert it           
           !     final g-matrix, setup matrix to invert
           ALLOCATE(a(n_confs, n_confs))
           DO ie=1, n_startenergy_g  
              a=0.0_dp
              a(:,:) = REAL(gfree_ee(:,:,ie))+REAL(pep(:,:,ie))
              !     Invert matrix
              CALL matinv(a,n_confs,d)
              gfree_ee(:,:,ie) = a(:,:) 
           ENDDO
           DEALLOCATE(pep,a)
           ! Build final G-matrix
           DO ie=1, n_startenergy_g
              gfree_ee(:,:,ie)=MATMUL(gfree_ee(:,:,ie),gfree_e(:,:,ie))
              gna(:,:,ie)=gfree(:,:,ie)-MATMUL(gfree_r(:,:,ie),gfree_ee(:,:,ie))
              IF (isospin_z == -1)  THEN
                 gna(:,:,ie)=gna(:,:,ie)+coulomb(:,:)
              ENDIF
           ENDDO
           DEALLOCATE(gfree, gfree_e, gfree_r, gfree_ee)
           CALL gprint(isospin_z,p_parity,ang_mom,gna,twobody_com,twobody_r2, twobody_p2,gmatrix_configs)
           !     free space
           DEALLOCATE(gmatrix_configs%config_ab)
           DEALLOCATE(coulomb); DEALLOCATE(twobody_com)
           DEALLOCATE(twobody_r2,twobody_p2)
           DEALLOCATE(gna)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE setupg
!
!                    Print out G-mtx
!
SUBROUTINE gprint(it,ip,ij,gna,twobody_com,twobody_r2, twobody_p2,gmatrix_configs)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, ijd, it, ip, ij, ia, ib, ic,id,ie
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(IN) :: twobody_com,twobody_r2, twobody_p2
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_g), INTENT(IN)  :: gna

  ijd=ij+ij
  DO i=1,gmatrix_configs%number_confs
     ia=gmatrix_configs%config_ab(i*2-1)
     ib=gmatrix_configs%config_ab(i*2)
     DO j=i,gmatrix_configs%number_confs
        ic=gmatrix_configs%config_ab(j+j-1)
        id=gmatrix_configs%config_ab(j+j)
        !  write final G-matrix with no BHF calculation
        WRITE(7,1001) it, ip, ijd,ia, ib, ic, id,((gna(j,i,ie)),ie=1,n_startenergy_g), &
                                                  twobody_com(j,i), twobody_r2(j,i), &
                                                  twobody_p2(j,i)       
     ENDDO
  ENDDO
1001 FORMAT(7I4,5X,15(E12.6,1X))

END SUBROUTINE gprint
!
!          Find  free G-mtx, and contributions to
!          Pauli corrections
!
SUBROUTINE gmtx_free(gfree,gfree_e,gfree_r,gfree_ee,coulomb,&
     twobody_com,twobody_r2, twobody_p2,lab_to_relcoeff, lab_to_relconf,&
     lab_to_relnumber,max_coeff,gmatrix_configs,isospin_z)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE partial_waves
  USE relcm_gmatrix
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor) , INTENT(IN)  :: gmatrix_configs
  INTEGER :: k2,k1,i,j,nlc1,nlc2, conf1, conf2, lcm1, lcm2, chan1, chan2, &
       max_coeff, n_optimal, find_address_all, ncm1, ncm2, nr1, lr1, nr2, lr2, isospin_z 
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_g), INTENT(INOUT) :: &
       gfree, gfree_e, gfree_r, gfree_ee
  REAL(DP) :: dcoem, e_com, e_r2, e_rcm, e_p2, e_p2com
  REAL(DP), DIMENSION(gmatrix_configs%number_confs,max_coeff), &
       INTENT(IN)  :: lab_to_relcoeff
  INTEGER, DIMENSION(gmatrix_configs%number_confs,max_coeff), &
       INTENT(IN)  :: lab_to_relconf
  INTEGER, DIMENSION(gmatrix_configs%number_confs), &
       INTENT(IN)  :: lab_to_relnumber
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(INOUT) :: coulomb
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(INOUT) :: twobody_com, twobody_r2, twobody_p2


  gfree = 0.; gfree_e = 0. ; gfree_r = 0.; gfree_ee = 0.
  twobody_com = 0.0D0; twobody_r2=0.0D0;  twobody_p2=0.0D0
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
              ncm1=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan1,conf1+conf1))
              IF ( lcm1 /= lcm2 ) CYCLE   ! cm orbital mom must be equal
              lr1=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan1,conf1+conf1-1))
              nr1=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan1,conf1+conf1-1))
              lr2=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan2,conf2+conf2-1))
              nr2=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan2,conf2+conf2-1))
              dcoem=lab_to_relcoeff(k1,i)*lab_to_relcoeff(k2,j)
              n_optimal = find_address_all (chan1, conf1, conf2, 0, 'tot')
              gfree(k1,k2,:)= gfree(k1,k2,:)+ gtf_tot_all( 1, n_optimal, :)*dcoem       
              gfree_e(k1,k2,:)= gfree_e(k1,k2,:)+ gtf_tot_all(2, n_optimal, :)*dcoem
              gfree_r(k1,k2,:)= gfree_r(k1,k2,:)+gtf_tot_all(3,n_optimal,:)*dcoem
              gfree_ee(k1,k2,:)= gfree_ee(k1,k2,:)+gtf_tot_all(4,n_optimal, :)*dcoem
              !  Coulomb correction to be added perturbatively  at the end       
              IF ( (ncm1 == ncm2).AND.( lr1 == lr2 ).AND.( isospin_z == -1 ))  THEN
                 coulomb(k1,k2) = coulomb(k1,k2)+coulomb_relcom(nr1,nr2,lr1)*dcoem
              ENDIF
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
              IF ( ncm1 == ncm2 ) THEN  
                 twobody_com(k1,k2)=twobody_com(k1,k2)+dcoem*e_com
                 twobody_r2(k1,k2) = twobody_r2(k1,k2)+dcoem*e_r2
              ENDIF
              ! The p_ip_j correction depends however on N and N' and n and n'
              twobody_p2(k1,k2)=twobody_p2(k1,k2)+dcoem*(e_p2com-e_p2)*0.5D0
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE gmtx_free
!
!                    obtain the free G-mtx for negative and positive energies
!
SUBROUTINE g_channel(i)
  USE relcm_gmatrix
  USE partial_waves
  USE configurations
  USE single_particle_orbits
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i  ! loop variable over NN channels
  INTEGER :: ncoup,nc,ncp,lc, n,np, bra, ket, k1, k2, a, b, k3, k4, c, d
  INTEGER :: cm, ev, la, lb, lcp, jang, lim1, lim2, lim3, lim4, ii, i1
  REAL(DP) :: k0, e0
  REAL(DP), ALLOCATABLE ::  vkk(:,:)
  REAL(DP), ALLOCATABLE ::  krel_temp(:)
  REAL(DP), ALLOCATABLE ::  gfree(:,:,:,:)
  REAL(DP), ALLOCATABLE ::  gfree_temp(:,:), g_bhf(:,:)

  ncoup=1
  IF ( orb_lrel_max(i) == jang_rel(i) ) ncoup = 1
  IF ( orb_lrel_max(i) /= jang_rel(i) ) ncoup = 2
  ALLOCATE(vkk (ncoup*n_rel,ncoup*n_rel))
  ALLOCATE(g_bhf (ncoup*n_rel,ncoup*n_rel) )
  jang=jang_rel(i)
  !     setup the NN interaction
  CALL potential_interface(ncoup,vkk,i)
  ALLOCATE(gfree ( ncoup*n_rel,ncoup*n_rel,n_startenergy_g,n_cm_mom))
  !     construct free G-matrix with no Pauli operator and in k-space
  WRITE(6,*) 'Setting up g-matrix for partial wave quantum numbers:'
  WRITE(6,*) 'J,S,Tz:', jang_rel(i), spin_rel(i), iso(i)
  gfree = 0.0_dp
  IF ((type_of_pot == 'OPEP').OR.(type_of_pot == 'Tensorinteraction')   &
       .OR.(type_of_pot == 'LSinteraction') ) THEN
     DO cm=1, n_cm_mom
        DO ev=1, n_startenergy_g
           gfree(:,:,ev,cm) = vkk(:,:) 
        ENDDO
     ENDDO
  ELSE
     DO cm=1, n_cm_mom
        DO ev=1, n_startenergy_g
           g_bhf = 0.0_dp
           CALL g_mtx(ncoup,vkk,g_bhf,ev,cm)
           CALL g_interpolate(ncoup,g_bhf)
           gfree(:,:,ev,cm) = g_bhf(:,:) 
        ENDDO
     ENDDO
  ENDIF
  DEALLOCATE(vkk);   DEALLOCATE(g_bhf)
  !     make now transformation to h.o. basis in the rel and cm system
  !     loop over all cm and rel coordinate configurations
  DO bra =1, relcm_conf%nconfs_relcm(i)
     k2=bra+bra ; k1=k2-1
     b=relcm_conf%relcm_ab(i,k2)
     a=relcm_conf%relcm_ab(i,k1)
     n=relcm_sp_data%nrel(a)
     la=relcm_sp_data%lrel(a)
     nc=relcm_sp_data%nrel(b)
     lc=relcm_sp_data%lrel(b)
     IF ((nc+nc+lc+n+n+la) > nlmax) CYCLE
     DO ket =1,  relcm_conf%nconfs_relcm(i) 
        k4=ket+ket ;  k3=k4-1
        d=relcm_conf%relcm_ab(i,k4)
        c=relcm_conf%relcm_ab(i,k3)
        np=relcm_sp_data%nrel(c)        
        lb=relcm_sp_data%lrel(c)
        ncp=relcm_sp_data%nrel(d)
        lcp=relcm_sp_data%lrel(d)
        IF ( lcp /= lc ) CYCLE
        IF ((ncp+ncp+lc+np+np+lb) > nlmax) CYCLE
        wave_cm(:)=rnlc(:,lc,nc)*rnlc(:,lc,ncp)
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
        DO ev=1, n_startenergy_g
           CALL g_ho_osc(lb, np, la, n, gfree(lim1:lim2,lim3:lim4,ev,:),bra,ket,i,ev)
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE(gfree)

END SUBROUTINE g_channel
!
!          Obtain the free g-mat in oscillator basis
!          Final G is in units of MeV or eV
!          Time consuming part. Less time-consuming algo is welcome!
!
SUBROUTINE g_ho_osc(lb, np, la,n, a, bra, ket, chan, ie)
  USE wave_functions
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: bra, ket, chan, ie, lb, np, la,n 
  INTEGER :: cm, i, j, n_optimal, find_address_all 
  REAL(DP) :: c1
  REAL(DP), DIMENSION(n_rel,n_rel,n_cm_mom), INTENT(IN) :: a
  REAL(DP) :: sum1, sum2, sum3, sum4, sum5, sum6
  REAL(DP), DIMENSION(n_rel) :: wave_bra, wave_ket

  sum1=0.0_dp; sum3=0.0_dp; sum5=0.0_dp; sum6=0.0_dp         
  wave_bra(:) = rnlr(:,lb,np); wave_ket = rnlr(:,la,n)
  DO cm=1, n_cm_mom 
     c1=wave_cm(cm)
     DO i=1, n_rel
        sum2=0.0_dp; sum4=0.0_dp
        DO j=1, n_rel
           sum2=sum2+wave_ket(j)*a(j,i,cm)
           sum4=sum4+wave_ket(j)*a(j,i,cm)*osc_propagator(j,cm,ie)
        ENDDO
        sum1=sum1+sum2*wave_bra(i)*c1
        sum3=sum3+sum4*c1*wave_bra(i)
        sum5=sum5+sum2*c1*wave_bra(i)*osc_propagator(i,cm,ie)
        sum6=sum6+sum4*c1*wave_bra(i)*osc_propagator(i,cm,ie)
     ENDDO
  ENDDO
  n_optimal = find_address_all (chan, bra, ket, 0, 'tot')
  gtf_tot_all(1,n_optimal,ie) = sum1
  gtf_tot_all(2,n_optimal,ie) = sum3
  gtf_tot_all(3,n_optimal,ie) = sum5
  gtf_tot_all(4,n_optimal,ie) = sum6

END SUBROUTINE g_ho_osc
!
!          This routine contains old-fashioned common blocks
!          acting as interface between the effective interaction program
!          and the potential routines of the Bonn type, Nijmegen or
!          Argonne groups.
!          vkk is in units of MeV^-2,
!          rgkk are mesh points in rel coordinates,
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
SUBROUTINE potential_interface(ncoup,vkk,ichan)
  USE wave_functions
  USE relcm_gmatrix
  USE partial_waves
  USE constants
  IMPLICIT NONE
  INTEGER :: i,j, ncoup, k, l, jxx, ichan, inn, isospin_tz, spin, n1, ix, iy
  REAL(DP)  :: v,xmev,ymev, c, q
  CHARACTER (LEN=4) :: label
  COMMON /cnn/ inn
  COMMON /cpts/   q(97),c,n1,ix,iy
  COMMON /cpot/ v(6),xmev,ymev
  COMMON /cstate/ jxx, heform, sing, trip, coup, endep, label
  LOGICAL :: sing, trip, coup, heform, endep
  REAL(DP), DIMENSION(ncoup*n_rel,ncoup*n_rel), INTENT(OUT) :: vkk
  REAL(DP):: v00, v11, v12, v21, v22, fmult, besl, bs, b11, b12, b21, b22
  REAL(DP) :: xb(500), wb(500), v8(500,2,2), temp(500,2,2)
  INTEGER :: j1,l1,l2,s1,isot,m,nt,tz1,tz2,lpot

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
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    !  pp case = 1
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    !  pp case = 1
        tz1 = 1; tz2 = 1
     CASE (0)
        inn = 2    !  pn case = 2  if all inn=2, no CSB or ISB
        tz1 = 1; tz2 = -1
     CASE ( 1)
        IF ( (csb == 'csb').AND.(cib == 'cib'))       inn = 3    !  nn case = 3
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    !  nn case = 3
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    !  nn case = 3
        tz1 = -1; tz2 = -1
     END SELECT
  CASE('CD-bonn')
     SELECT CASE ( isospin_tz)
     CASE (-1)
        IF ( (csb == 'csb').AND.(cib == 'cib'))       inn = 1    !  pp case = 1
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    !  pp case = 1
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    !  pp case = 1
        tz1 = 1; tz2 = 1
     CASE (0)
        inn = 2    !  pn case = 2  if all inn=2, no CSB or ISB
        tz1 = 1; tz2 = -1
     CASE ( 1)
        IF ( (csb == 'csb').AND.(cib == 'cib'))       inn = 3    !  nn case = 3
        IF ( (csb == 'no-csb').AND.(cib == 'cib'))    inn = 3    !  nn case = 3
        IF ( (cib == 'no-cib').AND.(csb == 'no-csb')) inn = 2    !  nn case = 3
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
        fmult = 2.D0/ACOS(-1.D0)/p_mass(0)/hbarc
        tz1 = 1; tz2 = -1
     CASE ( 1)
        fmult = 2.D0/ACOS(-1.D0)/p_mass(1)/hbarc
        tz1 = -1; tz2 = -1
     END SELECT
     lpot = 1  ! switch for argonne potentials, lpot 1 means v18
     IF (type_of_pot == 'argonnev8') lpot = 2   ! lpot = 2 means v8
     IF (type_of_pot == 'OPEP') lpot = 3   ! lpot = 3 means OPEP from V18
     IF (type_of_pot == 'Tensorinteraction') lpot = 4   ! lpot = 4 means full tensor from V18
     IF (type_of_pot == 'LSinteraction') lpot = 5   ! lpot = 5 means full LS from V18
     IF (type_of_pot == 'MalflietTjon') lpot = 6   ! lpot = 5 means full LS from V18
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
                 bs=besl(xb(m)*rgkk(i),l1)*besl(xb(m)*rgkk(j),l1)
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
                 b11=besl(xb(m)*rgkk(i),l1)
                 b12=besl(xb(m)*rgkk(i),l2)
                 b21=besl(xb(m)*rgkk(j),l1)
                 b22=besl(xb(m)*rgkk(j),l2)
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
        xmev=rgkk(i)*hbarc
        ix = i
        IF (ncoup == 2)   k=i+n_rel
        DO j=1,n_rel
           ymev=rgkk(j)*hbarc
           iy = j
           IF (ncoup == 2) l=j+n_rel
           SELECT CASE (type_of_pot)
           CASE('CD-bonn')
              CALL cdbonn
           CASE('n3lo')
              CALL n3lo
           CASE('Idaho-A')
              CALL idaho
           CASE('Idaho-B')
              CALL idaho
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

END SUBROUTINE potential_interface
!
!                 Set up h.o. wf for rel cm system and lab frame
!                 It computes also a complex
!
SUBROUTINE ho_wfunction
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER :: n, l, i, j
  REAL(DP) :: ph, oscl_r, oscl_cm, sum_rel, sum_cm
  REAL(DP) :: hbarc3
  REAL(DP)  :: cx(0:200), factor, z_lab, z_rel, z_cm, xp, zc_lab, xpc, sum_relc

  hbarc3=hbarc**3
  oscl_r=oscl*SQRT(2.)        ! Oscillator parameter for relative
  oscl_cm=oscl*SQRT(0.5)      ! Oscillator parameter for cm points
  DO n=0,nmax
     ph=(-1.D0)**n
     DO l=0,lmax
        factor = 0.5D0*((n+l+2)*LOG(2.D0)+fac(n)-dfac(2*n+2*l+1)-0.5D0*LOG(pi))
        factor = EXP(factor)
        sum_rel=0.0_dp; sum_relc = 0.0_dp
        DO i=1,n_rel
           !  real ho wave function
           z_lab= ra(i)*ra(i)*oscl*oscl; cx = 0.0_dp
           CALL laguerre_general( n, l+0.5D0, z_lab, cx )
           xp = EXP(-z_lab*0.5D0)*((ra(i)*oscl)**l)*cx(n)
           hol(i,l,n) = xp*ph*factor*(oscl**(1.5D0)) ! lab wf
           z_rel= ra(i)*ra(i)*oscl_r*oscl_r; cx = 0.0_dp
           CALL laguerre_general( n, l+0.5D0, z_rel, cx )
           xp = EXP(-z_rel*0.5)*((ra(i)*oscl_r)**l)*cx(n)
           rnlr(i,l,n) = xp*(wra(i)*(ra(i)**2))*ph*factor*(oscl_r**(1.5D0))     ! rel w
           sum_rel=sum_rel+ wra(i)*(hol(i,l,n)*ra(i))**2
           !  complex HO wave functions
           zc_lab= ra(i)*ra(i)*oscl*oscl; cx = 0.0_dp
           CALL laguerre_general( n, l+0.5D0, zc_lab, cx )
           xpc = EXP(-zc_lab*0.5D0)*((ra(i)*oscl)**l)*cx(n)
           chol(i,l,n) = xpc*ph*factor*(oscl**(1.5D0)) ! lab wf
           sum_relc=sum_relc+ wra(i)*((chol(i,l,n)*ra(i))**2)
        ENDDO
        WRITE(6,'(21H Norm cm ho wf n,l : ,2I3,2X,3F12.7)') n, l, sum_rel, sum_relc
        sum_cm=0.
        DO j=1,n_cm_mom
           z_cm=rca(j)*rca(j)*oscl_cm*oscl_cm; cx = 0.0_dp
           CALL laguerre_general( n, l+0.5D0, z_cm, cx )
           xp = EXP(-z_cm*0.5D0)*((rca(j)*oscl_cm)**l)*cx(n)
           rnlc(j,l,n)=xp*ph*factor*(oscl_cm**(1.5D0))   ! cm wf
           sum_cm=sum_cm+wrca(j)*(rnlc(j,l,n)*rca(j))**2
           rnlc(j,l,n)= rnlc(j,l,n)*rca(j)*SQRT(wrca(j)*hbarc3)
        ENDDO
        WRITE(6,'(21H Norm cm ho wf n,l : ,2I3,2X,F12.7)') n, l, sum_cm
     ENDDO
  ENDDO

END SUBROUTINE ho_wfunction
!
!                 Set Coulomb interaction
!                 in coordinate space, nuclear case, returns in MeV
!
SUBROUTINE coulomb_integral
  USE constants
  USE wave_functions
  USE relcm_gmatrix
  IMPLICIT NONE
  INTEGER :: n1, n2, lr, i, nrel_max
  REAL(DP) :: oscl_r, int_sum, xr, xp, z, factor1, factor2
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rr, wrr
  REAL(DP) :: cx(0:200)

  oscl_r=oscl*SQRT(2.0D0)        ! Oscillator parameter for relative
  nrel_max = 400
  ALLOCATE ( rr(nrel_max ), wrr(nrel_max ))
  CALL gauss_legendre(0.d0, 20.d0, rr, wrr, nrel_max )
  DO lr = 0, lmax, 1
     DO n1 = 0, nmax, 1
        factor1 = 0.5D0*((n1+lr+2)*LOG(2.D0)+fac(n1)-dfac(2*n1+2*lr+1)-0.5D0*LOG(pi))
        factor1 = EXP(factor1)
        DO n2 = 0, nmax,1
           IF ( n2 == n1) THEN
              factor2 = factor1
           ELSE
              factor2 = 0.5D0*((n2+lr+2)*LOG(2.D0)+fac(n2)-dfac(2*n2+2*lr+1)-0.5D0*LOG(pi))
              factor2 = EXP(factor2)
           ENDIF
           int_sum = 0.D0
           DO i=1,nrel_max
              z= rr(i)/oscl_r
              CALL laguerre_general( n1, lr+0.5D0, z*z, cx )
              xp = cx(n1)*EXP(-z*z*0.5)*(z**lr)
              IF ( n1 == n2) THEN 
                 xr = xp
              ELSE
                 CALL laguerre_general( n2, lr+0.5D0, z*z, cx )
                 xr = cx(n2)*EXP(-z*z*0.5)*(z**lr)
              ENDIF
              int_sum=int_sum+wrr(i)*rr(i)*xp*xr
           ENDDO
           !  Coulomb energy in MeV or eV
           coulomb_relcom(n2, n1, lr) = int_sum*factor1*factor2*1.439965183D0/(oscl_r**3)
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE ( rr, wrr)

END SUBROUTINE coulomb_integral
!
!         Propagator used in g-matrix code
!         in units of MeV^2 
!
SUBROUTINE k_space_propagator
  USE wave_functions
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  INTEGER :: ev, cm
  REAL(DP) :: xsum

  propagator = 0.0_dp; osc_propagator = 0.0_dp
  DO ev=1, n_startenergy_g     
     DO cm=1, n_cm_mom
        xsum=-e_start_g(ev)+0.25*rca(cm)*rca(cm)*hb2ip
        propagator(:,cm,ev)=(hbarc**3)*wgkk(:)*rgkk(:)*rgkk(:)/ &   
             (xsum+(rgkk(:)*rgkk(:)*hb2ip))
        osc_propagator(:,cm,ev)=1.0_dp/ &
             (xsum+(ra(:)*ra(:)*hb2ip))
     ENDDO
  ENDDO

END SUBROUTINE k_space_propagator
!
!        Function used to get the free g-matrix in k-space for
!        starting energies less than E_d
!
SUBROUTINE g_mtx(ncoup,vkk,g_bhf,ev,cm)
  USE wave_functions
  USE constants
  IMPLICIT NONE
  INTEGER :: i, ncoup, i1,ev,cm, j, naux, iopt, ndim
  REAL(DP), DIMENSION(ncoup*n_rel,ncoup*n_rel) ::  g_bhf
  REAL(DP), DIMENSION(ncoup*n_rel,ncoup*n_rel) ::  vkk
  REAL(DP), ALLOCATABLE :: a(:,:),aux(:)
  REAL(DP) :: det(2), rcond

  iopt=0
  ndim=ncoup*n_rel
  naux=200*ndim
  ALLOCATE ( aux(naux)) 
  ALLOCATE(a (ncoup*n_rel,ncoup*n_rel))
  a = 0.
  DO i=1,ndim
     i1=i
     IF(i > n_rel) i1=i-n_rel
     DO j = 1, ndim
        a(j,i)=vkk(j,i)*propagator(i1,cm,ev)
     ENDDO
     a(i,i)=a(i,i)+1.D0
  ENDDO
  !     Provided library function from Numerical Recipes
  CALL matinv(a,ndim)
  !     Essl library call
  !      CALL dgeicd(a, ndim, ndim, iopt, rcond, det, aux, naux)
  g_bhf=MATMUL(a,vkk)
  DEALLOCATE(a); DEALLOCATE (aux)

END SUBROUTINE g_mtx

!
!            Setup pspace_energy  which enters in P(e)^-1P. 
!            It is in turn
!            multiplied with the free G-matrix and then
!            used to invert the final matrix in the function
!            gfinal
!

SUBROUTINE pspace_energy(itz,ij,pep,gmatrix_configs)
  USE single_particle_orbits
  USE constants
  USE relcm_gmatrix
  USE configurations
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: itz,ij,k1,k2,ie,i1,i2,i3,i4,n1,n2,n3,n4,l1,l2,l3,l4, &
       j1,j2,j3,j4,j,k,it,it_min,it_max
  REAL(DP) :: ev, fnorm, dij
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_g), INTENT(INOUT):: pep
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: deno
  REAL(DP) :: e, e_dir,e_exc,c3, c4, c1, c2 
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: int_factor, kin_energy

  ALLOCATE (deno(n_rel, n_rel) )
  ALLOCATE (int_factor(n_rel), kin_energy(n_rel) )
  int_factor(:) = wra(:)*ra(:)*ra(:)
  SELECT CASE ( itz)
  CASE (-1)
     kin_energy(:) = 0.5*ra(:)*ra(:)*hbarc*hbarc/p_mass(-1)
  CASE (0)
     kin_energy(:) = 0.5*ra(:)*ra(:)*hbarc*hbarc/p_mass(0)
  CASE ( 1)
     kin_energy(:) = 0.5*ra(:)*ra(:)*hbarc*hbarc/p_mass(1)
  END SELECT
  DO ie=1,n_startenergy_g
     ev=e_start_g(ie)
     deno=0.
     DO j=1,n_rel
        c1=int_factor(j) ; c2=kin_energy(j)
        DO k=1,n_rel
           deno(k,j)=int_factor(k)*c1/( ev-c2-kin_energy(k) )
        ENDDO
     ENDDO
     DO k2=1,gmatrix_configs%number_confs
        DO k1=1,k2
           i1=gmatrix_configs%config_ab(k1+k1-1)
           i2=gmatrix_configs%config_ab(k1+k1)
           i3=gmatrix_configs%config_ab(k2+k2-1)
           i4=gmatrix_configs%config_ab(k2+k2)
           n1=all_orbit%nn(i1)
           n2=all_orbit%nn(i2)
           n3=all_orbit%nn(i3)
           n4=all_orbit%nn(i4)
           l1=all_orbit%ll(i1)
           l2=all_orbit%ll(i2)
           l3=all_orbit%ll(i3)
           l4=all_orbit%ll(i4)
           j1=all_orbit%jj(i1)
           j2=all_orbit%jj(i2)
           j3=all_orbit%jj(i3)
           j4=all_orbit%jj(i4)
           e=0.D0
           e_dir=0.d0
           e_exc=0.D0
           IF ( itz /= 0) THEN
              fnorm = 1.d0/ dij(i1,i2) /  dij(i3,i4)
           ELSEIF( itz == 0 ) THEN
              fnorm = 1.d0
           ENDIF
           !          direct part
           IF((l1==l3).AND.(l2==l4).AND.(j1==j3).AND.(j2==j4))THEN
              DO j=1,n_rel
                 c3=chol(j,l1,n1)*chol(j,l3,n3)
                 DO k=1,n_rel
                    e_dir=e_dir+chol(k,l2,n2)*chol(k,l4,n4)*deno(k,j)*c3
                 ENDDO
              ENDDO
           ENDIF
           !          exchange part
           IF((l1==l4).AND.(l2==l3).AND.(j1==j4).AND.(j2==j3) .AND. (itz /= 0) )THEN
              DO j=1,n_rel
                 c4=chol(j,l1,n1)*chol(j,l4,n4)
                 DO k=1,n_rel
                    e_exc=e_exc+chol(k,l2,n2)*chol(k,l3,n3)*c4*deno(k,j)
                 ENDDO
              ENDDO
           ENDIF
           e= (e_dir-e_exc*((-1.0_dp)**((2*ij-j3-j4)/2)))*fnorm
           pep(k1,k2,ie)= e
           pep(k2,k1,ie)= e
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE ( deno)
  DEALLOCATE (int_factor, kin_energy )

END SUBROUTINE pspace_energy
!
!           Setup moshinsky brackets for transf
!           rel-cm frame ----> lab frame
!
!
SUBROUTINE mosh_transf(it,ij,gmatrix_configs,lab_to_relcoeff, &
     lab_to_relconf,lab_to_relnumber,max_coeff)
  USE single_particle_orbits
  USE configurations
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN) :: gmatrix_configs
  INTEGER, INTENT(IN) :: it, ij, max_coeff 
  INTEGER :: k1,k2,i,j,n1,l1,j1d,n2,l2,j2d,k, number
  REAL(DP), DIMENSION(gmatrix_configs%number_confs,max_coeff), & 
       INTENT(INOUT) :: lab_to_relcoeff
  INTEGER, DIMENSION(gmatrix_configs%number_confs,max_coeff), &
       INTENT(INOUT) :: lab_to_relconf 
  INTEGER, DIMENSION(gmatrix_configs%number_confs), &
       INTENT(INOUT) :: lab_to_relnumber 
  REAL(DP), DIMENSION(max_coeff) :: coeb
  INTEGER, DIMENSION(max_coeff) :: mqnb

  DO k=1,gmatrix_configs%number_confs              
     k2=k*2
     k1=k2-1
     i=gmatrix_configs%config_ab(k1)
     j=gmatrix_configs%config_ab(k2)
     n1=all_orbit%nn(i) 
     l1=all_orbit%ll(i)
     j1d=all_orbit%jj(i)
     n2=all_orbit%nn(j)
     l2=all_orbit%ll(j)
     j2d=all_orbit%jj(j)
     CALL transf_mbs(n1,l1,j1d,n2,l2,j2d,it,ij,mqnb,coeb,max_coeff,number)
     lab_to_relcoeff(k,:)=coeb(:)
     lab_to_relconf(k,:)=mqnb(:)
     lab_to_relnumber(k)=number
  ENDDO

END SUBROUTINE  mosh_transf
!
!           Explicit evaluation of transf coeff
!
SUBROUTINE transf_mbs(n1,l1,j1d,n2,l2,j2d,itz,ij,mqnb,coeb,max_coeff,number)
  USE single_particle_orbits
  USE partial_waves
  USE ang_mom_functions
  USE configurations
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n1,l1,j1d,n2,l2,j2d,itz,ij,max_coeff
  INTEGER, DIMENSION(max_coeff), INTENT(OUT) :: mqnb
  REAL(DP), DIMENSION(max_coeff), INTENT(OUT) :: coeb
  INTEGER :: ld, l, n, lcd, lambd, isign, lambdd,  a, b, k1, k2, &
       l2d, l1d, nemx, ijd, l_min, l_max, &
       jlsd, nch, isd, jls, ispin,nc, lc, nconfs
  INTEGER, INTENT(OUT) :: number
  REAL(DP) :: bk, w6, w9, s , co, cst, delta
  LOGICAL triag

  mqnb=0 ; coeb=0.
  l1d=l1*2 ;   l2d=l2*2 ;   ijd=ij*2 ;   nemx=2*n1+l1+2*n2+l2
  number=0
  !     loop over possible partial wave channels which specify J, l, S and Tz
  DO nch=1,no_channels
     IF ( itz /= iso(nch) ) CYCLE
     l_min=orb_lrel_min(nch)
     l_max=orb_lrel_max(nch)
     ispin=spin_rel(nch) ; isd=2*ispin
     jls=jang_rel(nch) ; jlsd=2*jls
     !     loop over rel and cm configs specifying n, l, N and L
     DO nconfs =1, relcm_conf%nconfs_relcm(nch)
        k2=nconfs*2
        k1=k2-1
        b=relcm_conf%relcm_ab(nch,k2)
        a=relcm_conf%relcm_ab(nch,k1)
        n=relcm_sp_data%nrel(a)
        l=relcm_sp_data%lrel(a)
        IF ( (l /= l_max) .AND. ( l /= l_min ) ) CYCLE
        nc=relcm_sp_data%nrel(b)
        lc=relcm_sp_data%lrel(b)
        !     pauli test for identical particles in partial waves
        IF((ABS(itz) == 1).AND. ((-1)**(l+ispin+1) > 0)) CYCLE
        !     get right sign from isospin clebsch-gordan coeff for the partial waves
        !     since we have for Tz=0 always a coupling order | (pn) J >
        !     this means that the T=0 contrib goes like  -(1-(-)^(l+spin))/sqrt(2)
        !     and l+s have to be odd       
        cst = 1.0_dp
        IF ( itz /= 0 ) THEN
           cst = SQRT(2.D0)/SQRT( 1.d0+ delta(n1,n2)*delta(l1,l2)*delta(j1d,j2d) )
        ELSEIF ( (itz == 0) ) THEN
           IF ( (MOD(l+ispin,2) /= 0 ) )  cst= -1.0_dp
        ENDIF
        IF(triag(jls,ij,lc)) CYCLE
        IF((-1)**(l1+l2+l+lc) < 0) CYCLE
        IF((2*n+l+2*nc+lc) /= nemx) CYCLE
        IF((2*n+l+2*nc+lc) > nlmax) CYCLE
        IF(triag(jls,l,ispin)) CYCLE
        ld=l*2 ; lcd=lc*2
        co=0.D0
        DO lambd = ABS(l1-l2), l1+l2
           lambdd=lambd*2
           IF(triag(lambd,l,lc)) CYCLE
           IF(triag(ij,lambd,ispin)) CYCLE
           bk= gmosh(n,l,nc,lc,n1,l1,n2,l2,lambd,1.D0)
!           IF( ABS (bk) == 0.D0) CYCLE
           w9=snj(l1d,1,j1d,l2d,1,j2d,lambdd,isd,ijd)
!           IF ( ABS( w9) == 0.D0) CYCLE
           w9=w9*SQRT((j1d+1.)*(j2d+1.)*(lambdd+1.)*(isd+1.))
           isign=(lcd+ld+ijd+isd)/2
           s=1.-2.*MOD(isign,2)
           w6=s*sjs(lcd,ld,lambdd,isd,ijd,jlsd)
           w6=w6*SQRT(FLOAT((lambdd+1)*(jlsd+1)))
           co=co+cst*bk*w6*w9*(-1.)**(l+lambd+jls+ij)
        ENDDO
        IF(DABS(co) <= 1.D-8) co = 0.0_dp!THEN
           number=number+1
           mqnb(number)=nconfs*1000+nch
           coeb(number)=co
!        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE transf_mbs
!
!           Finds total number of transformation  coeffs 
!           for given J, parity and isospin projection
!
SUBROUTINE find_ncoeffs(it,ij,this,max_coeff)
  USE single_particle_orbits
  USE configurations
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN) :: this
  INTEGER, INTENT(IN) :: it, ij 
  INTEGER :: k1,k2,i,j,n1,l1,n2,l2,number, k
  INTEGER, INTENT(OUT) :: max_coeff
  max_coeff=0
  DO k=1, this%number_confs              
     k2=k*2
     k1=k2-1
     i=this%config_ab(k1)
     j=this%config_ab(k2)
     n1=all_orbit%nn(i) 
     l1=all_orbit%ll(i)
     n2=all_orbit%nn(j)
     l2=all_orbit%ll(j)
     CALL count_trans_mbs(n1,l1,n2,l2,it,ij,number)
     IF ( max_coeff < number ) max_coeff = number
  ENDDO

END SUBROUTINE find_ncoeffs
!
!           counts max possible number  of transf coeff
!
SUBROUTINE count_trans_mbs(n1,l1,n2,l2,itz,ij,number)
  USE single_particle_orbits
  USE partial_waves
  USE configurations
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: number
  INTEGER, INTENT(IN) :: n1,l1,n2,l2,itz,ij
  INTEGER :: l, n, nemx, l_min, l_max, nconfs , nch, jls, & 
       ispin, nc, lc, a, b,  k1, k2 
  LOGICAL triag

  nemx=2*n1+l1+2*n2+l2
  number=0
  DO nch=1,no_channels
     IF ( itz /= iso(nch) ) CYCLE
     l_min=orb_lrel_min(nch)
     l_max=orb_lrel_max(nch)
     ispin=spin_rel(nch)
     jls=jang_rel(nch)
     DO nconfs =1, relcm_conf%nconfs_relcm(nch)
        k2=nconfs*2
        k1=k2-1
        b=relcm_conf%relcm_ab(nch,k2)
        a=relcm_conf%relcm_ab(nch,k1)
        n=relcm_sp_data%nrel(a)
        l=relcm_sp_data%lrel(a)
        IF ( (l /= l_max) .AND. ( l /= l_min ) ) CYCLE
        nc=relcm_sp_data%nrel(b)
        lc=relcm_sp_data%lrel(b)
        IF((ABS(itz)==1).AND.((-1)**(l+ispin+1) > 0)) CYCLE
        IF(triag(jls,ij,lc)) CYCLE
        IF((-1)**(l1+l2+l+lc) < 0) CYCLE
        IF((2*n+l+2*nc+lc) /= nemx) CYCLE
        IF((2*n+l+2*nc+lc) > nlmax) CYCLE
        IF(triag(jls,l,ispin)) CYCLE
        number=number+1
     ENDDO
  ENDDO

END SUBROUTINE count_trans_mbs
!
!      Interpolates free-gmatrix for negative starting energies from a mesh
!      with k \in [ 0, infty] to k \in [ 0, limit], where limit is set 
!      by the max extension of the ho wave function with most nodes
!
SUBROUTINE g_interpolate(ncoup,g_bhf)
  USE wave_functions
  USE constants
  IMPLICIT NONE
  INTEGER :: i, ncoup, i1,ev,cm, j, ndim, j1, lim1, lim2, lim3, lim4
  REAL(DP), DIMENSION(ncoup*n_rel,ncoup*n_rel), INTENT(INOUT) ::  g_bhf
  REAL(DP), ALLOCATABLE :: a(:,:), b(:,:)
  REAL(DP) :: mtx_element

  ndim=ncoup*n_rel

  ALLOCATE(a (ncoup*n_rel,ncoup*n_rel))
  ALLOCATE(b (ncoup*n_rel,ncoup*n_rel))
  a = 0.
  b = DBLE(g_bhf)
  DO i=1,ndim
     lim1 = 1; lim2 = n_rel 
     i1=i
     IF(i > n_rel) THEN 
        i1=i-n_rel
        lim1 = n_rel+1
        lim2 = n_rel+n_rel
     ENDIF
     DO j = 1, ndim
        lim3 = 1; lim4 = n_rel
        j1 = j
        IF(j > n_rel ) THEN 
           j1 = j -n_rel
           lim3 = n_rel+1
           lim4 = n_rel+n_rel
        ENDIF
        CALL lagrange_2dim(ra(j1),ra(i1),b(lim3:lim4,lim1:lim2),mtx_element)
        a(j, i) = mtx_element
     ENDDO
  ENDDO
  g_bhf = a

  DEALLOCATE(a, b )

END SUBROUTINE g_interpolate



     
