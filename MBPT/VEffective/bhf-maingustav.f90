!             Program block bhf-main.f    
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F90/F95 
!             LAST UPGRADE : June 2008
!                            
!             Program to calculate effective interactions 
!             for finite nuclei within the spherical shell
!             model. The program runs in the proton-neutron
!             formalism. It allows also to perform calculations
!             for electrons confined by a spherical 3dim harmonic
!             oscillator interaction interacting via the Coulomb
!             interaction. This runs strictly in atomic units only
!             Syntax throughout is that of free source format,
!             additional compiler options may be necessary
!             The code computes effective interactions for
!             the nuclear shell model, 3Dim quantum dots, coupled-cluster input and 
!             and the binding energy of closed-shell nuclei.
!             It allows one to 
!             perform Hartree-Fock calculations as well.
!
!     Main program starts here
!
PROGRAM  bhf_effective
  USE ang_mom_functions
  USE constants
  USE wave_functions
  USE inifile
  IMPLICIT NONE
  CHARACTER (LEN=120) :: infilename, outputfile, renorminteraction_file
  CHARACTER (LEN=120) :: HFrenorminteraction_file, spdata_file, HFspdata_file, vpp_file, vpn_file, vnn_file
  CHARACTER (LEN=120) :: CCinteraction_file, CCspdata_file, CCkinetic_file
  LOGICAL :: fail

  !     Read in and organize all relevant data, the file bhf.ini contains all data
  infilename   = 'bhf.ini'
  call ini_open(infilename, 5, fail, .false.)
  IF (fail) STOP 'Error opening parameter file, probably wrong file, use bhf.ini as name'
  !     open the output file
  outputfile = ini_read_string('output_run')
  OPEN(unit=6,file=outputfile)
  renorminteraction_file = Ini_Read_String('renorminteraction_file')
  OPEN(UNIT=7,FILE=renorminteraction_file)
  spdata_file = Ini_Read_String('spdata_file')
  OPEN(UNIT=12,FILE=spdata_file)
  order_of_interaction = Ini_Read_String('order_of_interaction')
  type_of_interaction = Ini_Read_String('type_of_interaction')
  type_of_operator = Ini_Read_String('type_of_operator')
  lambda = Ini_Read_Int('lambda')
  hf_iterations = Ini_Read_Int('hf_iterations')
  physical_system = ini_read_string('physical_system')
  n_body_int = ini_read_string('manybody_interaction')

  IF ( hf_iterations > 0 ) THEN
     HFrenorminteraction_file = Ini_Read_String('HFrenorminteraction_file')
     HFspdata_file = Ini_Read_String('HFspdata_file')
     OPEN(UNIT=9,FILE=HFrenorminteraction_file)
     OPEN(UNIT=8,FILE=HFspdata_file)
  ENDIF

  SELECT CASE (physical_system) 
  CASE('nuclear_physics')
     IF ( type_of_interaction == 'open-diagrams') THEN
        vpp_file = Ini_Read_String('veff_pp')
        vpn_file = Ini_Read_String('veff_pn')
        vnn_file = Ini_Read_String('veff_nn')
        OPEN(UNIT=13,FILE=vpp_file)
        OPEN(UNIT=14,FILE=vpn_file)
        OPEN(UNIT=15,FILE=vnn_file)
     ENDIF
  CASE('atomic_physics')
     IF ( type_of_interaction == 'open-diagrams') THEN
        vnn_file = Ini_Read_String('veff_nn')
        OPEN(UNIT=15,FILE=vnn_file)
     ENDIF
  END SELECT

  IF ( type_of_interaction == 'coupled-cluster') THEN
     CCkinetic_file = Ini_Read_String('CCkinetic_file')
     CCinteraction_file = Ini_Read_String('CCinteraction_file')
     CCspdata_file = Ini_Read_String('CCspdata_file')
     OPEN (UNIT=10,FILE=CCkinetic_file)
     OPEN (UNIT=11,FILE=CCspdata_file)
     OPEN (UNIT=14,FILE=CCinteraction_file)
  ENDIF
  CALL read_data
  !    Setup max number of mesh points needed
  CALL setup_ho_cutoff
  !     Factorials for 3j, 6j and 9j symbols and for the moshinsky transf coeffs
  CALL commons_to_angmom  
  !     Setup of effective interaction depending on which
  !     type of effective interaction and physical system is desired
  SELECT CASE (type_of_interaction)
  CASE ('core-diagrams')
     IF ( hf_iterations > 0 )   THEN
        CALL  setup_hfmatrix
     ELSEIF ( hf_iterations == 0) THEN
        !     reserve space in memory for various arrays
        ALLOCATE ( ra (n_rel), wra (n_rel));  ALLOCATE ( rgkk (n_rel), wgkk (n_rel))
        ALLOCATE ( hol (n_rel, 0:lmax, 0:nmax) )
        !     set up mesh points in lab frame
        CALL rel_mesh                   
        !     setup ho wave functions
        CALL ho_wfunction                
        CALL core_diagrams
        DEALLOCATE ( rgkk, wgkk, ra, wra) ; DEALLOCATE ( hol)
     ENDIF
  CASE ('open-diagrams')
     IF ( hf_iterations > 0 )   THEN
        CALL  setup_hfmatrix
     ELSEIF ( hf_iterations == 0) THEN
        !     reserve space in memory for various arrays
        ALLOCATE ( ra (n_rel), wra (n_rel));  ALLOCATE ( rgkk (n_rel), wgkk (n_rel))
        ALLOCATE ( hol (n_rel, 0:lmax, 0:nmax) )
        !     set up mesh points in lab frame
        CALL rel_mesh                   
        !     setup ho wave functions
        CALL ho_wfunction                
        IF (n_body_int == 'onebody') THEN 
           CALL onebody_contribution
           CALL effective_operator
        ELSEIF (n_body_int == 'twobody') THEN 
           CALL onebody_contribution
           CALL twobody_contribution
        ELSEIF (n_body_int == 'threebody') THEN 
           CALL onebody_contribution
           CALL twobody_contribution
           CALL threebody_contribution
        ENDIF
        DEALLOCATE ( rgkk, wgkk, ra, wra) ; DEALLOCATE ( hol)
     ENDIF
  CASE ('coupled-cluster')
     IF ( hf_iterations > 0 ) THEN
        CALL  setup_hfmatrix
     ELSEIF ( hf_iterations == 0) THEN
        !     reserve space in memory for various arrays
        ALLOCATE ( ra (n_rel), wra (n_rel));  ALLOCATE ( rgkk (n_rel), wgkk (n_rel))
        ALLOCATE ( hol (n_rel, 0:lmax, 0:nmax) )
        !     set up mesh points in lab frame
        CALL rel_mesh                   
        !     setup ho wave functions
        CALL ho_wfunction  
        CALL ccsd_data_model
        SELECT CASE ( order_of_interaction )
        CASE ('zero')
            CALL coupled_cluster_passthrough
        CASE ('first')
           CALL coupled_cluster_twobody
        CASE ('second')
           CALL onebody_contribution
           CALL twobody_contribution
        CASE ('third')
           CALL onebody_contribution
           CALL twobody_contribution
        END SELECT
        DEALLOCATE ( rgkk, wgkk, ra, wra) ; DEALLOCATE ( hol)
     ENDIF
  END SELECT

END PROGRAM bhf_effective
!
!    Read in single-particle data
!
SUBROUTINE read_data
  USE constants
  USE single_particle_orbits  
  USE stored_bare_interaction
  USE inifile
  IMPLICIT NONE
  CHARACTER :: text
  INTEGER :: i, j, k, nlarge_max, nmodel_max, neutrons, protons, modeloscillator
  CHARACTER(LEN= 10) space, model
  INTEGER :: number_twobody_elements, number_pp_elements, number_pn_elements, number_nn_elements 

  n_startenergy_veff = Ini_Read_Int('n_startenergy_veff')
  IF ( n_startenergy_veff < 0) THEN
     WRITE(6,*) 'Number of starting energies less than zero'
     STOP 'startingenergyveffwrong'
  ENDIF
  !  always odd number of starting energies
  IF ( MOD(n_startenergy_veff,2) == 0)n_startenergy_veff =  n_startenergy_veff+1
  !     starting energy at which the effective interactions is evaluated
  starting_energy = Ini_Read_Double('starting_energy')
  ALLOCATE ( wcn (n_startenergy_veff) ) 
  IF ( n_startenergy_veff == 1) THEN
     wcn(1) = 0.0_dp
  ELSE
     !  center the starting energies around the chosen starting energy
     DO i=1, n_startenergy_veff
        wcn(i)=-1.*i+(n_startenergy_veff/2+1)
     ENDDO
  ENDIF
  !  get the number of hbaromega excitations
  number_homega_exct = Ini_Read_Int('hbaromega_excitations')
  !  For the atomic physics case, mass_nucleus = number_electrons
  READ(12,*)
  READ(12,'(A65,I30)') text, mass_nucleus 
  SELECT CASE (physical_system) 
  CASE('atomic_physics')
     number_electrons = mass_nucleus
  END SELECT
  adependency = Ini_Read_String('A-dependence')
  a_factor = 1.0_dp
  keep_originalenergies = Ini_Read_String('keep_original_spenergies')
  IF ( adependency =='yes') a_factor = (1.0_dp-1.0_dp/mass_nucleus)
  READ(12,'(A30,E12.6,2X,E12.6)')  text, oscl, hbar_omega
  READ(12,*); READ(12,*); READ(12,*)
  READ(12,'(A56,I10)') text, nlarge_max
  READ(12,'(A56,I10)') text, nmodel_max
  READ(12,'(A39,I15)') text, all_orbit%total_orbits 
  READ(12,*) 
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits) 
  DO j = 1, all_orbit%total_orbits
     !     READ(12,'(A7,5(I4,2X),2X,E12.6,2X,A10,2X,A10)') text, i, all_orbit%nn(j), all_orbit%ll(j), &
     READ(12,*) text, i, all_orbit%nn(j), all_orbit%ll(j), &
          all_orbit%jj(j), all_orbit%itzp(j), all_orbit%nshell(j), all_orbit%e(j),&
          all_orbit%evalence(j), all_orbit%orbit_status(j), all_orbit%model_space(j)
     all_orbit%e_original(j) = all_orbit%e(j)
  ENDDO
  j_lab_max = 0; i = 0; nmax = 0; lmax = 0
  DO j = 1, all_orbit%total_orbits
     IF ( all_orbit%jj(j) > j_lab_max) j_lab_max = all_orbit%jj(j)
     IF ( all_orbit%ll(j) > lmax) lmax = all_orbit%ll(j)
     IF ( all_orbit%nn(j) > nmax) nmax = all_orbit%nn(j)
  ENDDO
  ! Here we fix itzmin, itzmax and j_lab_min
  SELECT CASE (physical_system) 
  CASE('nuclear_physics')
     itzmin = -1; itzmax=1; j_lab_min=0
  CASE('atomic_physics')
     itzmin = 1; itzmax=1; j_lab_min=0
  END SELECT

  READ(7,*); READ(7,*)
  READ(7,'(A20,A20)')text, type_of_renormv
  IF ( type_of_renormv =='no-core') THEN
     nlmax = 2*nmax+lmax
  ELSEIF ( type_of_renormv =='v-nrg') THEN
     nlmax = 2*nmax+lmax
  ELSE
     IF ( nlarge_max == nmodel_max) THEN
        nlmax  = nmodel_max
     ELSEIF ( nlarge_max > nmodel_max) THEN
        nlmax = 2*nmax
     ENDIF
  ENDIF
  READ(7,'(A38,I4)') text, n_startenergy_g
  IF ( n_startenergy_g == 0) n_startenergy_g = 1
  ALLOCATE(e_start_g(n_startenergy_g))
  READ(7,*)  (e_start_g(i), i=1,n_startenergy_g )
  READ(7,'(A38,4I15 )') text, number_twobody_elements, number_pp_elements, number_pn_elements, number_nn_elements
  SELECT CASE (physical_system) 
  CASE('atomic_physics')
     number_nn_elements = number_twobody_elements
     number_pp_elements = 0
     number_pn_elements = 0
  END SELECT
  READ(7,*); READ(7,*)
  ! Then read the matrix elements
  CALL read_gmatrix(number_twobody_elements, number_pp_elements, number_pn_elements, number_nn_elements)
  ! Now print info to new files for the Hartree-Fock output, single-particle data and interaction data
  IF ( hf_iterations > 0 ) THEN
     WRITE(8,'(68H   ----> Oscillator parameters, Model space and single-particle data)')
     SELECT CASE (physical_system) 
     CASE('nuclear_physics')
        WRITE(8,'(65HMass number A of chosen nucleus (important for CoM corrections): ,I4)') mass_nucleus 
     CASE('atomic_physics')
        WRITE(8,'(65HNumber of electrons (important for CoM corrections)            : ,I4)') number_electrons 
     END SELECT
     WRITE(8,'(30HOscillator length and energy: ,E12.6,2X,E12.6)')  oscl, hbar_omega
     WRITE(8,*) 'Max value of l in lab fram= ', lmax
     WRITE(8,*) 'Max value of n in lab frame :', nmax
     WRITE(8,*) 'Max value of 2*n + l+ cm 2*N +L for large space:', nlarge_max
     WRITE(8,*) 'Max value of 2*n + l+ cm 2*N +L for model space:', nmodel_max
     WRITE(8,*) 'Min and max value of partial wave ang. mom', jmin, jmax
     WRITE(8,*) 'Total number of single-particle orbits', all_orbit%total_orbits 
     WRITE(8,'(55HLegend:         n   l  2j  tz  energy  particle or hole)') 
     WRITE(9,'(34H   ----> Interaction part         )')
     WRITE(9,'(34HNucleon-Nucleon interaction model:,A20)') type_of_pot
     WRITE(9,'(20HType of calculation:,A20)') type_of_renormv
     WRITE(9,'(38HNumber and value of starting energies:,I4)')  n_startenergy_g  
     WRITE(9,'(10(1X,E12.6) )') (e_start_g(i), i=1,n_startenergy_g )
     WRITE(9,'(38HTotal number of twobody matx elements:,4I15 )') number_twobody_elements, number_pp_elements, &
          number_pn_elements, number_nn_elements 
     WRITE(9,'(91HMatrix elements with the following legend, NOTE no hbar_omega/A for Hcom, p_ip_j and r_ir_j)')
     WRITE(9,'(103H  Tz Par  2J   a   b   c   d          <ab|V|cd>        <ab|Hcom|cd>     <ab|r_ir_j|cd>   <ab|p_ip_j|cd>)') 
  ENDIF

END SUBROUTINE read_data
!
!                  SUBROUTINE to fix the cutoff in mesh points
!                  for the harmonic oscillator wave functions
!                  which depend on the oscillator parameter.    
!
SUBROUTINE setup_ho_cutoff
  USE single_particle_orbits
  USE constants
  USE partial_waves
  USE wave_functions  
  IMPLICIT NONE
  INTEGER :: h,nh, lh, iq, number_of_iterations, int_points
  REAL(DP) :: sigma, norm, oscl_r
  REAL(DP) :: qmin, qmax, sum_hf, sum_norm(0:lmax)
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q
  REAL(DP) :: z_rel, factor, xp, ph, cx(0:200), contrib

  sigma = 1.0_dp; norm = 1.0_dp; int_points = 10
  number_of_iterations = 0
  nh = nmax; oscl_r=oscl*SQRT(2.)
  !  for momentum space if we set z= 0.5*(oscl_r*k)**2 > 60-70, then
  !  EXP(-60) < 10E-27. Making it twice as large ensures that we account
  !  for extensions due to the Laguerre polynoms which depend on z**(2n).
  cutoff = 2*SQRT(60.0_dp)/oscl
  qmin = 0.0D0; qmax = cutoff
  DO WHILE( (number_of_iterations < 20) .AND. (ABS(sigma) > 1E-4) )         
     ALLOCATE ( q_points(int_points), weight_q(int_points) )
     CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
     sum_norm  = 0.0D0
     DO h=0, lmax
        lh = h
        sum_hf = 0.
        factor = 0.5D0*((nh+lh+2)*LOG(2.D0)+fac(nh)-dfac(2*nh+2*lh+1)-0.5D0*LOG(pi))
        factor = EXP(factor)
        DO iq=1, int_points
           z_rel= q_points(iq)*q_points(iq)*oscl_r*oscl_r
           CALL laguerre_general( nh, lh+0.5D0, z_rel, cx )
           xp = EXP(-z_rel*0.5)*((q_points(iq)*oscl_r)**lh)*cx(nh)
           contrib = xp*factor*(oscl_r**(1.5D0)) 
           sum_hf=sum_hf+ weight_q(iq)*(contrib*q_points(iq))**2
        ENDDO
        sum_norm(h) = sum_hf
     ENDDO
     sigma = 0.D0
     DO h = 0, lmax
        sigma = sigma +ABS( sum_norm(h)-norm)
     ENDDO
     sigma = sigma/(lmax+1)
     number_of_iterations = number_of_iterations +1 
     !     WRITE(6,*) 'Sigma for this iteration', sigma, number_of_iterations
     DEALLOCATE ( weight_q, q_points)
     IF (ABS(sigma) > 1E-4) THEN 
        int_points = int_points+10
     ENDIF
  ENDDO
  n_rel = int_points
  n_cm_mom = int_points
  WRITE(6,'(41H New HO cutoff and # integration points: ,E12.6,2X,I5)') &
       cutoff, n_rel

END SUBROUTINE setup_ho_cutoff
