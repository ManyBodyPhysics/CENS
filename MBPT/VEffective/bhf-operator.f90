!             Program block bhf-operator.f
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F90 
!             LAST UPGRADE : oct 1997
!
!             This program block sets up the one-body 
!             effective operator with folded diagrams as well, 
!             Diagrams up to second  order in the interaction G in
!             standard Rayleigh-Schroedinger theory can be
!             included.
!                        
!
!     Sets up all diagrams through second order in the interaction G
!     The type of effective operator is defined by the input values of 
!     the variable type_of_operator
!
SUBROUTINE effective_operator
  USE bare_operator_value
  USE wave_functions
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  REAL(DP) :: q
  INTEGER :: a, c, iph, iz, ia, r_loop, number_operators, j_initial, j_final
  REAL(DP) :: corepol, phase
  REAL(DP), DIMENSION(n_startenergy_veff) :: effoper_diagram_1, &
       effoper_diagram_2,effoper_diagram_3,effoper_diagram_4, &
       effoper_diagram_5,effoper_diagram_6,effoper_diagram_7, &
       effoper_diagram_8,effoper_diagram_9, effoper_diagram_10, &
       effoper_diagram_11,effoper_diagram_12,effoper_diagram_13, &
       effoper_diagram_14,effoper_diagram_15,effoper_diagram_16, &
       effoper_diagram_17,effoper_diagram_18,effoper_diagram_19, &
       effoper_diagram_20,effoper_diagram_21,effoper_diagram_22, &
       effoper_diagram_23,effoper_diagram_24,effoper_diagram_25, &
       effoper_diagram_26,effoper_diagram_27,effoper_diagram_28, &
       effoper_diagram_29, effoper_diagram_30,effoper_diagram_31, &
       effoper_diagram_32,eff_box, effoper_diagram_23f, &
       effoper_diagram_24f, effoper_diagram_27f,effoper_diagram_28f
  LOGICAL triag
  !     set up the bare operator depending on the type wanted 
  ALLOCATE( bare_operator(all_orbit%total_orbits, all_orbit%total_orbits))
  DO a = 1, all_orbit%total_orbits
     DO c= a, all_orbit%total_orbits
        bare_operator(a,c) =0.0_dp
     ENDDO
  ENDDO
  SELECT CASE ( type_of_operator ) 
  CASE ('none') 
     RETURN
  CASE ('gamow-teller')
     CALL gamtel
  CASE ('electro')
     CALL em_matrix_element
  CASE ('magnetic')
     CALL magnetic_matrix_element
  END SELECT
  WRITE(6,*)' Effective operator for model space states, bare and higher order'

  DO a = 1, all_orbit%total_orbits
     IF (all_orbit%model_space(a) == 'outside') CYCLE         
     DO c= a, all_orbit%total_orbits
        phase=iph(lambda-(all_orbit%jj(a)-all_orbit%jj(c))/2)*SQRT(2.*lambda+1.)
        !            CALL one_body_diagram_phase(a,c,phase)
        IF (all_orbit%model_space(c) == 'outside') CYCLE
        eff_box = 0.0_dp   
        SELECT CASE ( order_of_interaction )
        CASE ('first')
           CALL effective_operator_1(a,c,effoper_diagram_1)
           CALL effective_operator_2(a,c,effoper_diagram_2)
           eff_box=effoper_diagram_1+effoper_diagram_2
        CASE ('second')
           CALL effective_operator_1(a,c,effoper_diagram_1)
           CALL effective_operator_2(a,c,effoper_diagram_2)
           CALL effective_operator_3(a,c,effoper_diagram_3)
           CALL effective_operator_4(a,c,effoper_diagram_4)
           CALL effective_operator_5(a,c,effoper_diagram_5)
           CALL effective_operator_6(a,c,effoper_diagram_6)
           CALL effective_operator_7(a,c,effoper_diagram_7)
           CALL effective_operator_8(a,c,effoper_diagram_8)
           CALL effective_operator_9(a,c,effoper_diagram_9)
           CALL effective_operator_10(a,c,effoper_diagram_10)
           CALL effective_operator_11(a,c,effoper_diagram_11)
           CALL effective_operator_12(a,c,effoper_diagram_12)
           CALL effective_operator_13(a,c,effoper_diagram_13)
           CALL effective_operator_14(a,c,effoper_diagram_14)
           CALL effective_operator_15(a,c,effoper_diagram_15)
           CALL effective_operator_16(a,c,effoper_diagram_16)
           CALL effective_operator_17(a,c,effoper_diagram_17)
           CALL effective_operator_18(a,c,effoper_diagram_18)
           CALL effective_operator_19(a,c,effoper_diagram_19)
           CALL effective_operator_20(a,c,effoper_diagram_20)
           CALL effective_operator_21(a,c,effoper_diagram_21)
           CALL effective_operator_22(a,c,effoper_diagram_22)
           CALL effective_operator_23(a,c,effoper_diagram_23)
           CALL effective_operator_24(a,c,effoper_diagram_24)
           CALL effective_operator_23folded(a,c,effoper_diagram_23f)
           CALL effective_operator_24folded(a,c,effoper_diagram_24f)
           CALL effective_operator_25(a,c,effoper_diagram_25)
           CALL effective_operator_26(a,c,effoper_diagram_26)
           CALL effective_operator_27(a,c,effoper_diagram_27)
           CALL effective_operator_28(a,c,effoper_diagram_28)
           CALL effective_operator_27folded(a,c,effoper_diagram_27f)
           CALL effective_operator_28folded(a,c,effoper_diagram_28f)
           CALL effective_operator_29(a,c,effoper_diagram_29)
           CALL effective_operator_30(a,c,effoper_diagram_30)
           CALL effective_operator_31(a,c,effoper_diagram_31)
           CALL effective_operator_32(a,c,effoper_diagram_32)
           eff_box=effoper_diagram_1+effoper_diagram_2+ &
                effoper_diagram_3+effoper_diagram_4+ &
                effoper_diagram_5+effoper_diagram_6+ &
                effoper_diagram_7+effoper_diagram_8+ &
                effoper_diagram_9+effoper_diagram_10 +&
                effoper_diagram_11+effoper_diagram_12 + &
                effoper_diagram_13+effoper_diagram_14+ &
                effoper_diagram_15+effoper_diagram_16+ &
                effoper_diagram_17+effoper_diagram_18+ &
                effoper_diagram_19+effoper_diagram_20 + &
                effoper_diagram_21+ effoper_diagram_22+ &
                effoper_diagram_23+effoper_diagram_24+ &
                effoper_diagram_25+effoper_diagram_26+ &
                effoper_diagram_27+effoper_diagram_28+ & 
                effoper_diagram_29+effoper_diagram_30+ &
                effoper_diagram_31+effoper_diagram_32 + &
                effoper_diagram_23f+effoper_diagram_24f+ &
                effoper_diagram_27f+effoper_diagram_28f
        CASE ('third')
           WRITE(6,*) 'Warning, effective operator only to second order'  
           CALL effective_operator_1(a,c,effoper_diagram_1)
           CALL effective_operator_2(a,c,effoper_diagram_2)
           CALL effective_operator_3(a,c,effoper_diagram_3)
           CALL effective_operator_4(a,c,effoper_diagram_4)
           CALL effective_operator_5(a,c,effoper_diagram_5)
           CALL effective_operator_6(a,c,effoper_diagram_6)
           CALL effective_operator_7(a,c,effoper_diagram_7)
           CALL effective_operator_8(a,c,effoper_diagram_8)
           CALL effective_operator_9(a,c,effoper_diagram_9)
           CALL effective_operator_10(a,c,effoper_diagram_10)
           CALL effective_operator_11(a,c,effoper_diagram_11)
           CALL effective_operator_12(a,c,effoper_diagram_12)
           CALL effective_operator_13(a,c,effoper_diagram_13)
           CALL effective_operator_14(a,c,effoper_diagram_14)
           CALL effective_operator_15(a,c,effoper_diagram_15)
           CALL effective_operator_16(a,c,effoper_diagram_16)
           CALL effective_operator_17(a,c,effoper_diagram_17)
           CALL effective_operator_18(a,c,effoper_diagram_18)
           CALL effective_operator_19(a,c,effoper_diagram_19)
           CALL effective_operator_20(a,c,effoper_diagram_20)
           CALL effective_operator_21(a,c,effoper_diagram_21)
           CALL effective_operator_22(a,c,effoper_diagram_22)
           CALL effective_operator_23(a,c,effoper_diagram_23)
           CALL effective_operator_24(a,c,effoper_diagram_24)
           CALL effective_operator_23folded(a,c,effoper_diagram_23f)
           CALL effective_operator_24folded(a,c,effoper_diagram_24f)
           CALL effective_operator_25(a,c,effoper_diagram_25)
           CALL effective_operator_26(a,c,effoper_diagram_26)
           CALL effective_operator_27(a,c,effoper_diagram_27)
           CALL effective_operator_28(a,c,effoper_diagram_28)
           CALL effective_operator_27folded(a,c,effoper_diagram_27f)
           CALL effective_operator_28folded(a,c,effoper_diagram_28f)
           CALL effective_operator_29(a,c,effoper_diagram_29)
           CALL effective_operator_30(a,c,effoper_diagram_30)
           CALL effective_operator_31(a,c,effoper_diagram_31)
           CALL effective_operator_32(a,c,effoper_diagram_32)
           eff_box=effoper_diagram_1+effoper_diagram_2+ &
                effoper_diagram_3+effoper_diagram_4+ &
                effoper_diagram_5+effoper_diagram_6+ &
                effoper_diagram_7+effoper_diagram_8+ &
                effoper_diagram_9+effoper_diagram_10 +&
                effoper_diagram_11+effoper_diagram_12 + &
                effoper_diagram_13+effoper_diagram_14+ &
                effoper_diagram_15+effoper_diagram_16+ &
                effoper_diagram_17+effoper_diagram_18+ &
                effoper_diagram_19+effoper_diagram_20 + &
                effoper_diagram_21+ effoper_diagram_22+ &
                effoper_diagram_23+effoper_diagram_24+ &
                effoper_diagram_25+effoper_diagram_26+ &
                effoper_diagram_27+effoper_diagram_28+ & 
                effoper_diagram_29+effoper_diagram_30+ &
                effoper_diagram_31+effoper_diagram_32 + &
                effoper_diagram_23f+effoper_diagram_24f+ &
                effoper_diagram_27f+effoper_diagram_28f      
        END SELECT
        corepol = eff_box(n_startenergy_veff/2+1)*phase
        WRITE(6,'(2I3,2X,2F12.6,2X,I3)')a, c, bare_operator(a,c)*phase,corepol
     ENDDO
  ENDDO
  DEALLOCATE (bare_operator)

END SUBROUTINE effective_operator
!
!     calculates the reduced gamow-teller operat. 
!     from glaudemans  &  brussard, eq. 12.17
!     this follows the convention of Edmonds. That means
!     that the reduced matrix element used in the evaluation
!     of the various diagrams has a factor
!     (-)^(ji+lambda-jf)/sqrt(2lambda+1), see foonote 5
!     of Kuo et al. Ann. Phys. 132 (1981) 237
!
SUBROUTINE gamtel
  USE bare_operator_value
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: iph, nphase, ii, if, nf, lf, jf, itzf, ni, li, ji, itzi
  REAL(DP) :: factr

  DO if=1, all_orbit%total_orbits
     nf=all_orbit%nn(if)
     lf=all_orbit%ll(if)
     jf=all_orbit%jj(if)
     itzf=all_orbit%itzp(if)
     DO ii=1, all_orbit%total_orbits
        ni=all_orbit%nn(ii)
        li=all_orbit%ll(ii)
        ji=all_orbit%jj(ii)
        itzi=all_orbit%itzp(ii)
        IF(lf /= li) CYCLE ; IF (ni /= nf) CYCLE ; IF (itzi == itzf ) CYCLE
        !     expression of glaudemans and Brussard, eq. 12.17 
        nphase = (jf +3)/2+lf
        factr = 6.*SQRT((jf+1.)*(ji+1.))
        bare_operator(if,ii)=factr*sjs(1,1,2*lambda,ji,jf,2*lf)*iph(nphase)
        !     multiply with new factors
        bare_operator(if,ii)=bare_operator(if,ii)/SQRT(2.*lambda+1.)*iph(lambda+(ji-jf)/2)
     ENDDO
  ENDDO

END SUBROUTINE gamtel
!
!
!     calculates the reduced electromagnetic operator for given multipolarity
!     lambda provided in input file 
!     from glaudemans  &  brussard, eqs. 10.55 & 10.56 but with no isoscalar 
!     or isovector part. Note also typo in 10.56. Use eq. A.3.e5 
!     convention used is that of Edmonds or the expressions in the appendix
!     of Ring & Schuck, eq. B.81 and B.82. That means
!     that the reduced matrix element used in the evaluation
!     of the various diagrams has a factor
!     (-)^(ji+lambda-jf)/sqrt(2lambda+1), see foonote 5
!     of Kuo et al. Ann. Phys. 132 (1981) 237
!     The bare operator is only given for protons, since e=0 for neutrons
!
SUBROUTINE em_matrix_element
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE wave_functions
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: iph, nphase, ii, if, nf, lf, jf, itzf, ni, li, ji, itzi, i_mesh
  REAL(DP) :: factr, sum_rel
  LOGICAL triag
  !     loop over final sp state

  DO if=1, all_orbit%total_orbits
     nf=all_orbit%nn(if)
     lf=all_orbit%ll(if)
     jf=all_orbit%jj(if)
     itzf=all_orbit%itzp(if)
     !     only protons
     IF (itzf == 1) CYCLE
     !     loop over initial & final state
     DO ii=1, all_orbit%total_orbits
        ni=all_orbit%nn(ii)
        li=all_orbit%ll(ii)
        ji=all_orbit%jj(ii)
        itzi=all_orbit%itzp(ii)
        !     only protons
        IF (itzi == 1) CYCLE
        IF ( itzi /= itzf) CYCLE
        IF ( triag ( jf, ji, 2*lambda) ) CYCLE
        IF ( iph(li+lf) /= iph(lambda) ) CYCLE
        !     get h.o. integral
        IF ( hf_iterations == 0 ) THEN
           sum_rel=0.
           DO i_mesh=1,n_rel
              sum_rel=sum_rel+hol(i_mesh,lf,nf)*hol(i_mesh,li,ni)*wra(i_mesh)*(ra(i_mesh)**(2+lambda))
           ENDDO
        ELSE
           sum_rel=0.
           DO i_mesh=1,n_rel
              sum_rel=sum_rel+wave_function(i_mesh,if)*wave_function(i_mesh,ii)*wra(i_mesh)*(ra(i_mesh)**(2+lambda))
           ENDDO
        ENDIF
        !     expression of Ring & Schuck., Eq. B.81
        nphase= iph((jf-1)/2)*(1+iph(lf+li+lambda))
        factr=sum_rel*0.25*SQRT((jf+1.)*(ji+1.)*(2.*lambda+1.))&
             /SQRT(4.*ACOS(-1.))
        bare_operator(if,ii)=factr*tjs(jf,2*lambda,ji,-1,0,+1)*nphase
        !     multiply with iph(lambda+(ji-jf)/2)/SQRT(2.*lambda+1.)
        bare_operator(if,ii)=bare_operator(if,ii)/SQRT(2.*lambda+1.)*iph(lambda+(ji-jf)/2)
     ENDDO
  ENDDO

END SUBROUTINE em_matrix_element
!
!     calculates the reduced magnetic operator for given multipolarity
!     lambda provided in input file 
!     from glaudemans  &  brussard, eqs. 10.55 & 10.56 but with no isoscalar 
!     or isovector part. Note also typo in 10.56. Use eq. A.3.e5 
!     convention used is that of Edmonds or the expressions in the appendix
!     of Ring & Schuck, eq. B.81 and B.82. That means
!     that the reduced matrix element used in the evaluation
!     of the various diagrams has a factor
!     (-)^(ji+lambda-jf)/sqrt(2lambda+1), see foonote 5
!     of Kuo et al. Ann. Phys. 132 (1981) 237
!     no bohr-magneton in the expression
!     
!
SUBROUTINE magnetic_matrix_element
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE wave_functions
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: iph, nphase, ii, if, nf, lf, jf, itzf, ni, li, ji, itzi, i_mesh
  REAL(DP) :: factr, sum_rel, kappa
  REAL(DP), DIMENSION(-1:1) :: gs 
  REAL(DP), DIMENSION(-1:1) :: gl
  DATA gs(-1:1) /5.586D0,0.D0,-3.826D0/
  DATA gl(-1:1) /1.D0,0.D0,0.D0/
  LOGICAL triag

  !     loop over final sp state
  DO if=1, all_orbit%total_orbits
     nf=all_orbit%nn(if)
     lf=all_orbit%ll(if)
     jf=all_orbit%jj(if)
     itzf=all_orbit%itzp(if)
     !     loop over initial & final state
     DO ii=1, all_orbit%total_orbits
        ni=all_orbit%nn(ii)
        li=all_orbit%ll(ii)
        ji=all_orbit%jj(ii)
        itzi=all_orbit%itzp(ii)
        IF ( itzi /= itzf) CYCLE
        IF ( triag ( jf, ji, 2*lambda) ) CYCLE
        IF ( iph(li+lf) /= iph(lambda) ) CYCLE
        !     get h.o. integral
        IF ( hf_iterations == 0 ) THEN
           sum_rel=0.
           DO i_mesh=1,n_rel
              sum_rel=sum_rel+hol(i_mesh,lf,nf)*hol(i_mesh,li,ni)*wra(i_mesh)*(ra(i_mesh)**(1+lambda))
           ENDDO
        ELSE
           sum_rel=0.
           DO i_mesh=1,n_rel
              sum_rel=sum_rel+wave_function(i_mesh,if)*wave_function(i_mesh,ii)*wra(i_mesh)*(ra(i_mesh)**(1+lambda))
           ENDDO
        ENDIF
        !     expression of Ring & Schuck eq. B.82 
        nphase= iph((jf-1)/2)*(1+iph(lf+li+lambda))
        kappa=(ji+1.)/2.*iph((2*li+ji+1)/2) + &
             (jf+1.)/2.*iph((2*lf+jf+1)/2)
        factr=sum_rel*0.5*SQRT((jf+1.)*(ji+1.)*(2.*lambda+1.))&
             /SQRT(4.*ACOS(-1.))*(lambda-kappa)* &
             (0.5*gs(itzi)-gl(itzi)*(1.+kappa/(lambda+1.)) )
        bare_operator(if,ii)=factr*tjs(jf,2*lambda,ji,-1,0,+1)*nphase
        !     multiply with iph(lambda+(ji-jf)/2)/SQRT(2.*lambda+1.)
        bare_operator(if,ii)=bare_operator(if,ii)/SQRT(2.*lambda+1.)*iph(lambda+(ji-jf)/2)
        write(6,*) if, ii,  bare_operator(if,ii)
     ENDDO
  ENDDO

END SUBROUTINE magnetic_matrix_element

