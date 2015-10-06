!             Program block bhf-coupledcluster.f    
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F90 
!             LAST UPGRADE : September 2006
!                            
!             Program to calculate effective interactions 
!             for finite nuclei and atoms within the spherical shell
!             model. The program runs in the proton-neutron
!             formalism. 
!             The output is taylored for coupled cluster
!             calculations. 
!
!     The two-body effective interaction in J-scheme is set up here
!     This interaction is projected onto a smaller space, defined by the
!     specific model space
!
SUBROUTINE coupled_cluster_twobody
  USE constants
  USE configurations
  USE single_particle_orbits
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs
  INTEGER :: p_parity, isospin_z, ang_mom, bra_side, ket_side
  INTEGER :: sp_state_a, sp_state_b, sp_state_c, sp_state_d, pq_confs, p_confs
  REAL(DP) :: norm, dij, particle_mass
  REAL(DP), DIMENSION( n_startenergy_veff) :: ans
  REAL(DP), DIMENSION( n_startenergy_g+3) :: additional
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: sp_energy 
  COMPLEX*16, ALLOCATABLE::  heff(:,:), cvec_pp(:,:), cvec_qp(:,:), q_box(:,:), cspenergy(:,:)
  REAL(DP), ALLOCATABLE::  s_box(:,:,:)
  COMPLEX*16, ALLOCATABLE::  ceig_p(:)

  CALL onebody_contribution
  DO isospin_z=itzmin,itzmax
     DO p_parity=0,1
        DO ang_mom=j_lab_min,j_lab_max
           !     find all possible configurations, large space and model space
           CALL large_number_confs &
                (ang_mom,p_parity,isospin_z,gmatrix_configs)
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
           pq_confs=gmatrix_configs%number_confs
           ALLOCATE(gmatrix_configs%config_ab(pq_confs+pq_confs) )
           CALL number_configurations_model(ang_mom,p_parity,isospin_z,p_confs)
           gmatrix_configs%model_space_confs = p_confs
           IF ( p_confs <= 0 ) CYCLE
           CALL setup_configurations(ang_mom,p_parity,isospin_z,gmatrix_configs)
           ALLOCATE ( q_box(pq_confs,pq_confs) )
           ALLOCATE ( sp_energy(pq_confs,pq_confs) );  ALLOCATE ( cspenergy(pq_confs,pq_confs) )
           ALLOCATE ( s_box(pq_confs,pq_confs,n_startenergy_veff) )
           q_box=0.0_dp; s_box = 0.0_dp; sp_energy = 0.0_dp
           CALL setup_onebody_2_twobody(isospin_z,p_parity,ang_mom, sp_energy,s_box, gmatrix_configs)
           DO bra_side=1, gmatrix_configs%number_confs
              sp_state_a=gmatrix_configs%config_ab(bra_side*2-1)
              sp_state_b=gmatrix_configs%config_ab(bra_side*2)
              DO ket_side=1, gmatrix_configs%number_confs  
                 sp_state_c=gmatrix_configs%config_ab(ket_side*2-1)
                 sp_state_d=gmatrix_configs%config_ab(ket_side*2)
                 ans = 0.0_dp
                 norm=1./dij(sp_state_a,sp_state_b)/dij(sp_state_c,sp_state_d)
                 CALL diag1(sp_state_a, sp_state_b, sp_state_c, sp_state_d, ang_mom,ans)
                 q_box(bra_side,ket_side)= sp_energy(bra_side,ket_side) + ans(n_startenergy_veff/2+1)*norm
              ENDDO
           ENDDO
           ALLOCATE(heff(p_confs, p_confs))
           heff = (0.0_dp,0.0_dp)
           cspenergy  = sp_energy
           ALLOCATE(cvec_pp(p_confs, p_confs))
           ALLOCATE(cvec_qp(pq_confs-p_confs,p_confs))
           ALLOCATE(ceig_p(p_confs))
           !      Get model space eigenvalues, eigenvectors of 
           !      model space and excluded space and
           !      model space configurations which match the corresponding ones 
           !      of the large space
           CALL eigenvalues_large_maxvector(cvec_pp,cvec_qp,ceig_p,q_box,pq_confs,p_confs,gmatrix_configs)
           !           CALL eigenvalues_large(cvec_pp,cvec_qp,ceig_p,q_box,pq_confs,p_confs)
           !      setup 2p-effective interaction in square P-space using           
           !      the Okubo similarity transformation
           CALL lee_suzuki( cvec_pp, cvec_qp, ceig_p, p_confs, pq_confs-p_confs, heff, cspenergy)
           DEALLOCATE(cvec_pp); DEALLOCATE(cvec_qp); DEALLOCATE(ceig_p)
           DO bra_side=1, p_confs
              sp_state_a=gmatrix_configs%config_ab(bra_side*2-1)
              sp_state_b=gmatrix_configs%config_ab(bra_side*2)
              DO ket_side=bra_side, p_confs
                 sp_state_c=gmatrix_configs%config_ab(ket_side*2-1)
                 sp_state_d=gmatrix_configs%config_ab(ket_side*2)
                 CALL pphhmtx(sp_state_a, sp_state_b,sp_state_c, sp_state_d, ang_mom, additional)
                 additional = additional/dij(sp_state_a,sp_state_b)/dij(sp_state_c,sp_state_d)
                 WRITE(14,'(7I4,5X,4(E12.6,2X))') &
                      isospin_z,p_parity,2*ang_mom,sp_state_a  ,sp_state_b, sp_state_c, sp_state_d,  &
                      REAL(heff(ket_side,bra_side)), additional(n_startenergy_g+1), additional(n_startenergy_g+2),  &
                      additional(n_startenergy_g+3)
              ENDDO
           ENDDO
           DEALLOCATE (q_box,s_box); DEALLOCATE(heff); DEALLOCATE(sp_energy); DEALLOCATE(cspenergy)
           DEALLOCATE ( gmatrix_configs%config_ab) 
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE coupled_cluster_twobody
!
!   Prints out single-particle data after self-consistency.
!
SUBROUTINE ccsd_data_model
  USE constants
  USE single_particle_orbits
  USE wave_functions 
  USE relcm_gmatrix
  IMPLICIT NONE
  REAL(DP) :: ekin, hcom, particle_mass
  INTEGER :: a,na,np,npp,la,ja,tza,c,nc,lc,jc,tzc,k

  WRITE(11,*) hbar_omega
  np = 0 
  DO a = 1, all_orbit%total_orbits
     na = all_orbit%nn(a)
     la = all_orbit%ll(a)
     ja = all_orbit%jj(a)
     tza = all_orbit%itzp(a)
     IF ( all_orbit%model_space(a) == 'outside' ) CYCLE
     np = np + 1
     ! write single-particle data to disk for ccsd-calcs
     WRITE(11,'(1X,I3,1X,I3,1X,I3,1X,I3,1X,I3,5X,2(E16.10,1X))') np, na, la, ja, tza, DBLE( all_orbit%e(a) ), 0.d0
  ENDDO
  !
  ! calculate kinetic energies for ccsd-calcs 
  ! <a|k^2/m|c> diagonal in l,j,tz
  !
  np = 0
  DO a = 1, all_orbit%total_orbits
     na = all_orbit%nn(a) !- 1
     la = all_orbit%ll(a)
     ja = all_orbit%jj(a)
     tza = all_orbit%itzp(a)
     IF ( all_orbit%model_space(a) == 'outside' ) CYCLE
     np = np + 1
     npp = 0
     DO c = 1, all_orbit%total_orbits
        nc = all_orbit%nn(c) ! - 1
        lc = all_orbit%ll(c)
        jc = all_orbit%jj(c)
        tzc = all_orbit%itzp(c)
        IF ( all_orbit%model_space(c) == 'outside' ) CYCLE
        npp = npp + 1
        ekin = 0.d0
        IF ( la == lc .AND. ja == jc .AND. tza == tzc ) THEN
           SELECT CASE (physical_system) 
           CASE('nuclear_physics')
              particle_mass = p_mass(tza)
           CASE('atomic_physics')
              particle_mass = e_mass
           END SELECT
           IF ( hf_iterations > 0 ) THEN
              DO k=1,n_rel
                 ekin=ekin+(wave_function(k,c)*wave_function(k,a))*wra(k)*(ra(k)**4)*0.5*hbarc*hbarc/particle_mass
              ENDDO
           ELSE
              DO k=1,n_rel
                 ekin=ekin+(hol(k,la,na)*hol(k,lc,nc))*wra(k)*(ra(k)**4)*0.5*hbarc*hbarc/particle_mass
              ENDDO
           ENDIF
        ENDIF
        IF ( a == c ) THEN
           hcom = DBLE( all_orbit%e(a) )
        ELSE
           hcom = 0.d0
        ENDIF
        WRITE(10,'(1X,I3,1X,I3,5X,4(E16.10,1X))') np,npp,a_factor*DBLE(ekin), 0.d0, hcom, 0.d0
     ENDDO
  ENDDO

END SUBROUTINE ccsd_data_model
