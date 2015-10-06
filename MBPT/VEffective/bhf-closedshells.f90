!
!     setup of energy from closed shell core
!
SUBROUTINE core_diagrams
  USE constants
  USE single_particle_orbits
  USE wave_functions
  IMPLICIT NONE
  REAL(DP) :: ekin, first_order, sum_kin, closed_core_E, second_order, second_order_HF, &
       third_order1, third_order2, third_order3, correlation_energy, second_orderrecoupl, &
       third_orderhf, hf31, hf32, tdahf, rpa1hf, rpa2hf, particle_mass
  INTEGER :: h, n_hole, l_hole, j_hole, k
  REAL(DP) ::  d10(n_startenergy_veff), d11(n_startenergy_veff), d12(n_startenergy_veff), d13(n_startenergy_veff), &
       d14(n_startenergy_veff), d15(n_startenergy_veff)
  ekin=0.0_dp; first_order = 0.0_dp; second_order = 0.0_dp; closed_core_E = 0.0_dp
  third_order1 = 0.0_dp;  third_order2 =  0.0_dp;  third_order3 = 0.0_dp; second_order_HF = 0.0_dp  
  third_orderhf= 0.0_dp;  hf31= 0.0_dp; hf32= 0.0_dp; tdahf = 0.0_dp;  rpa1hf = 0.0_dp;  rpa2hf  = 0.0_dp 
  correlation_energy = 0.0_dp
  DO h=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
     n_hole = all_orbit%nn(h)
     l_hole = all_orbit%ll(h)
     j_hole = all_orbit%jj(h)
     sum_kin=0.0_dp
     SELECT CASE (physical_system) 
     CASE('nuclear_physics')
        particle_mass = p_mass(all_orbit%itzp(h))
        IF ( hf_iterations == 0 ) THEN
           DO k=1,n_rel
              sum_kin=sum_kin+(hol(k,l_hole,n_hole)**2)*wra(k)*(ra(k)**4)*0.5*hbarc*hbarc/particle_mass
           ENDDO
        ELSE
           DO k=1,n_rel
              sum_kin=sum_kin+(wave_function(k,h)**2)*wra(k)*(ra(k)**4)*0.5*hbarc*hbarc/particle_mass
           ENDDO
        ENDIF
     CASE('atomic_physics')
        IF ( hf_iterations == 0 ) THEN
           DO k=1,n_rel
              sum_kin=sum_kin+(hol(k,l_hole,n_hole)**2)*wra(k)*(ra(k)**4)*0.5
           ENDDO
        ELSE
           DO k=1,n_rel
              sum_kin=sum_kin+(wave_function(k,h)**2)*wra(k)*(ra(k)**4)*0.5
           ENDDO
        ENDIF
     END SELECT

     SELECT CASE (type_of_renormv)
     CASE ('no-core')
        ekin=ekin+2*sum_kin*(j_hole+1)
     CASE ('v-nrg')
        ekin=ekin+2*sum_kin*(j_hole+1)
     CASE ('vlowk')
        ekin=ekin+a_factor*sum_kin*(j_hole+1)
     CASE ('v-krg')
        ekin=ekin+a_factor*sum_kin*(j_hole+1)
     CASE ('g-matrix')
        ekin=ekin+a_factor*sum_kin*(j_hole+1)
     END SELECT
  ENDDO
  CALL core_E_first_order(first_order)
  WRITE(6,*) 'FIRST-ORDER=    ', first_order
  closed_core_E=ekin+first_order
  IF (order_of_interaction  == 'second') THEN
     CALL core_E_second_order(second_order)
     WRITE(6,*) 'SECOND-ORDER 2p-2h=    ', second_order
     IF (hf_iterations < 1) THEN
        CALL core_E_second_order_HF(second_order_HF)  ! two HF insertions
        WRITE(6,*) 'SECOND-ORDER HF=    ', second_order_HF
     ENDIF
     WRITE(6,*) 'TOTAL SECOND-ORDER=    ', second_order+second_order_HF
     closed_core_E=closed_core_E+second_order+second_order_HF
     correlation_energy = second_order+second_order_HF
  ELSEIF ((order_of_interaction  == 'third')) THEN
     CALL core_E_second_order(second_order)
     WRITE(6,*) 'SECOND-ORDER 2p-2h=    ', second_order
     IF (hf_iterations < 1) THEN
        CALL core_E_second_order_HF(second_order_HF)  ! two HF insertions
        WRITE(6,*) 'SECOND-ORDER HF=    ', second_order_HF
     ENDIF
     WRITE(6,*) 'TOTAL SECOND-ORDER=    ', second_order+second_order_HF
     CALL core_E_third_order1(third_order1)
     CALL core_E_third_order2(third_order2)
     CALL core_E_third_order3(third_order3)
     WRITE(6,*) 'THIRD-ORDER 4p-2h=    ', third_order1
     WRITE(6,*) 'THIRD-ORDER 4h-2p=    ', third_order2
     WRITE(6,*) 'THIRD-ORDER 3p-3h=    ', third_order3
     IF (hf_iterations < 1) THEN
        ! get third-order diagrams with HF insertions, using one-body diagrams
        DO h=1, all_orbit%total_orbits 
           IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
           d10 =0.0_dp; d11 = 0.0_dp; d12 = 0.0_dp
           d13 =0.0_dp; d14 = 0.0_dp; d15 = 0.0_dp
           CALL diagram_10(h,h,d10)
           CALL diagram_11(h,h,d11)
           CALL diagram_12(h,h,d12)
           CALL diagram_13(h,h,d13)
           CALL diagram_14(h,h,d14)
           CALL diagram_15(h,h,d15)
           third_orderhf = third_orderhf + (all_orbit%jj(h)+1)*(d10((n_startenergy_veff/2+1))+&
                d11((n_startenergy_veff/2+1))+d12((n_startenergy_veff/2+1))+d13((n_startenergy_veff/2+1)) &
                +d14((n_startenergy_veff/2+1))+d15((n_startenergy_veff/2+1)))
        ENDDO
        ! adding tda and rpa diagrams + 3 HF insertions (not from one-body code)
        CALL tda_E_third_order_HF(tdahf);  CALL rpa1_E_third_order_HF(rpa1hf)
        CALL rpa2_E_third_order_HF(rpa2hf); CALL hf31_E_third_order_HF(hf31); 
        CALL hf32_E_third_order_HF(hf32) 
        third_orderhf = third_orderhf+tdahf+rpa1hf+rpa2hf+hf31+hf32
        WRITE(6,*) 'THIRD-ORDER HF=    ', third_orderhf
     ENDIF
     closed_core_E=closed_core_E+second_order+third_order1+third_order2+third_order3+&
          third_orderhf+second_order_HF
     correlation_energy = second_order+third_order1+third_order2+third_order3+second_order_HF+third_orderhf
  ENDIF
  WRITE(6,*) 'Total energy per particle= ', closed_core_E/mass_nucleus
  WRITE(6,*) 'Total energy= ', closed_core_E
  WRITE(6,*) 'Kinetic energy= ', ekin
  WRITE(6,*) 'Correlation energy= ', correlation_energy
  WRITE(6,*) 'Correlation energy per particle= ', correlation_energy/mass_nucleus

END SUBROUTINE core_diagrams
