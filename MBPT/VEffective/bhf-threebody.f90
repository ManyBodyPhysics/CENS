!             Program block bhf-threebody.f    
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO, Norway 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F90/F95 
!             LAST UPGRADE : June 2008
!                            
!             Program to calculate effective 
!             three-body interactions for nuclei within the spherical shell
!             model. The program runs in the proton-neutron
!             formalism. 
!             Syntax throughout is that of free source format,
!             additional compiler options may be necessary
!
!     The three-body effective interaction in J-scheme is set up here
!
SUBROUTINE threebody_contribution
  USE constants
  USE configurations
  USE single_particle_orbits
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: mspace_configs
  REAL(DP),  DIMENSION(n_startenergy_veff ) :: two_d, d1a, d2a, d3a, d4a, d5a, d6a,d7a, d8a, d9a, &
       d1b, d2b, d3b, d4b, d5b, d6b,d7b, d8b, d9b
  INTEGER :: p_parity, isospin_z, ang_mom, bra_side, ket_side, jt_max
  INTEGER :: a, b, c, d, e, f, iph
  INTEGER :: nconfs, j_ab, j_de, i, nrot, np, ia
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) ::   vmtx, vmtx_vect
  REAL(DP), ALLOCATABLE, DIMENSION(:) ::   vmtx_eigen
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nmap

  ! Set up of output numbers for model space, starting with 1.
  ALLOCATE( nmap(all_orbit%total_orbits) )
  np = 0
  DO ia = 1, all_orbit%total_orbits
     IF(all_orbit%model_space(ia)=='outside')CYCLE
     np = np + 1
     nmap(ia) = np 
  ENDDO
  ! Find maximum value of three-particle J  
  jt_max = 0
  DO i = 1, all_orbit%total_orbits
     IF (all_orbit%jj(i) > jt_max) jt_max = all_orbit%jj(i) 
  ENDDO
  jt_max = 3*jt_max
  ! Loop over configs to find total interaction  for the given model space
  ! Note all possible values of isospin projection included by default
  DO isospin_z=-3,3,2
     DO ang_mom=1, jt_max, 2
        DO p_parity=0,1
           CALL number_threebody_configurations_model(ang_mom,p_parity,isospin_z,nconfs)
           mspace_configs%number_confs = nconfs
           IF ( mspace_configs%number_confs <= 0 ) CYCLE
           ALLOCATE ( mspace_configs%config_ab(4*nconfs) ) 
           CALL setup_threebody_configurations_model(ang_mom,p_parity,isospin_z,mspace_configs)
           ALLOCATE(vmtx(nconfs,nconfs)); ALLOCATE(vmtx_vect(nconfs,nconfs)); ALLOCATE(vmtx_eigen(nconfs))
           WRITE(6,'(42HNumber of configurations, J, Tz and parity,2x,4I4)') nconfs, ang_mom, isospin_z, p_parity
           DO bra_side=1, nconfs
              j_ab=mspace_configs%config_ab(bra_side*4)
              a=mspace_configs%config_ab(bra_side*4-1)
              b=mspace_configs%config_ab(bra_side*4-2)
              c=mspace_configs%config_ab(bra_side*4-3)
              DO ket_side=1, nconfs
                 j_de=mspace_configs%config_ab(ket_side*4)
                 d=mspace_configs%config_ab(ket_side*4-1)
                 e=mspace_configs%config_ab(ket_side*4-2)
                 f=mspace_configs%config_ab(ket_side*4-3)
                 CALL  two_2_three(a,b,c,d,e,f,j_ab,j_de,ang_mom,two_d)
                 CALL  three_body_1a(a,b,c,d,e,f,j_ab,j_de,ang_mom,d1a)
                 CALL  three_body_2a(a,b,c,d,e,f,j_ab,j_de,ang_mom,d2a)
                 CALL  three_body_3a(a,b,c,d,e,f,j_ab,j_de,ang_mom,d3a)
                 CALL  three_body_4a(a,b,c,d,e,f,j_ab,j_de,ang_mom,d4a)
                 CALL  three_body_5a(a,b,c,d,e,f,j_ab,j_de,ang_mom,d5a)
                 CALL  three_body_6a(a,b,c,d,e,f,j_ab,j_de,ang_mom,d6a)
                 CALL  three_body_7a(a,b,c,d,e,f,j_ab,j_de,ang_mom,d7a)
                 CALL  three_body_8a(a,b,c,d,e,f,j_ab,j_de,ang_mom,d8a)
                 CALL  three_body_9a(a,b,c,d,e,f,j_ab,j_de,ang_mom,d9a)
                 CALL  three_body_1b(a,b,c,d,e,f,j_ab,j_de,ang_mom,d1b)
                 CALL  three_body_2b(a,b,c,d,e,f,j_ab,j_de,ang_mom,d2b)
                 CALL  three_body_3b(a,b,c,d,e,f,j_ab,j_de,ang_mom,d3b)
                 CALL  three_body_4b(a,b,c,d,e,f,j_ab,j_de,ang_mom,d4b)
                 CALL  three_body_5b(a,b,c,d,e,f,j_ab,j_de,ang_mom,d5b)
                 CALL  three_body_6b(a,b,c,d,e,f,j_ab,j_de,ang_mom,d6b)
                 CALL  three_body_7b(a,b,c,d,e,f,j_ab,j_de,ang_mom,d7b)
                 CALL  three_body_8b(a,b,c,d,e,f,j_ab,j_de,ang_mom,d8b)
                 CALL  three_body_9b(a,b,c,d,e,f,j_ab,j_de,ang_mom,d9b)
                 vmtx(bra_side,ket_side) = two_d(n_startenergy_veff/2+1)&
                 + d1a(n_startenergy_veff/2+1)+d2a(n_startenergy_veff/2+1) &
                      +d3a(n_startenergy_veff/2+1)+d4a(n_startenergy_veff/2+1)  &
                      +d5a(n_startenergy_veff/2+1)+d6a(n_startenergy_veff/2+1)  &
                      +d7a(n_startenergy_veff/2+1)+d8a(n_startenergy_veff/2+1)  &
                      +d9a(n_startenergy_veff/2+1)+d1b(n_startenergy_veff/2+1)  &
                      +d2b(n_startenergy_veff/2+1)+d3b(n_startenergy_veff/2+1)  &
                      +d4b(n_startenergy_veff/2+1)+d5b(n_startenergy_veff/2+1)  & 
                      +d6b(n_startenergy_veff/2+1)+d7b(n_startenergy_veff/2+1)  &
                      +d8b(n_startenergy_veff/2+1)+d9b(n_startenergy_veff/2+1)
                 WRITE(6,'(10I4,5X,3F12.6)') &
                      a, b, c, d, e, f, j_ab, j_de, ang_mom, isospin_z, vmtx(bra_side,ket_side)
              ENDDO
           ENDDO
           DEALLOCATE ( mspace_configs%config_ab); DEALLOCATE(vmtx); DEALLOCATE(vmtx_vect)
           DEALLOCATE(vmtx_eigen)  
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE threebody_contribution




                    
                        
