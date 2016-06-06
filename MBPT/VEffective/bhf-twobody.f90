!             Program block bhf-twobody.f    
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO, Norway 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F90/F95 
!             LAST UPGRADE : Oct 2007
!                            
!             Program to calculate effective interactions 
!             for nuclei within the spherical shell
!             model. The program runs in the proton-neutron
!             formalism. 
!             Syntax throughout is that of free source format,
!             additional compiler options may be necessary
!             The code computes effective interactions for
!             the nuclear shell model and allows one to 
!             perform Hartree-Fock calculations as well.
!
!     The two-body effective interaction in J-scheme is set up here
!
SUBROUTINE twobody_contribution
  USE constants
  USE configurations
  USE single_particle_orbits
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: mspace_configs
  INTEGER :: p_parity, isospin_z, ang_mom, bra_side, ket_side
  INTEGER :: sp_state_a, sp_state_b, sp_state_c, sp_state_d, np, ia
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nmap
  INTEGER :: nconfs
  REAL(DP) :: u_pot, norm, dij
  REAL(DP), DIMENSION( n_startenergy_veff) :: two_body, one_body, dg1, dg2, dg3, dg4, d3rd
  REAL(DP), DIMENSION( n_startenergy_g) :: additional
  REAL(DP),  ALLOCATABLE, DIMENSION(:,:,:):: q_box, s_box
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: final_v, folded_twobody, folded_onebody, nonhermitian_v, sp_energy 
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: first, corepol, holeladder, particleladder, third, totalunfolded

  ALLOCATE( nmap(all_orbit%total_orbits) )
  np = 0
  DO ia = 1, all_orbit%total_orbits
     IF(all_orbit%model_space(ia)=='outside')CYCLE
     np = np + 1
     nmap(ia) = np 
  ENDDO
  WRITE(6,*) 'Two-body diagrams for < (ab) J | Veff | (cd) J>'   
  WRITE(6,'(127H     a        b        c        d    2Tz 2J    first         corepol     4p2h        &
       2p-ladder   third       total  Q    Folded)')   
  DO isospin_z=itzmin,itzmax
     DO ang_mom=j_lab_min,j_lab_max
        DO p_parity=0,1
           CALL number_configurations_model &
                (ang_mom,p_parity,isospin_z,nconfs)
           mspace_configs%number_confs = nconfs
           IF ( mspace_configs%number_confs <= 0 ) CYCLE
           ALLOCATE ( mspace_configs%config_ab(2*nconfs) ) 
           CALL setup_configurations_model &
                (ang_mom,p_parity,isospin_z,mspace_configs)
           ALLOCATE ( q_box(nconfs,nconfs,n_startenergy_veff) )
           ALLOCATE ( s_box(nconfs,nconfs,n_startenergy_veff) )  
           ALLOCATE ( folded_twobody(nconfs,nconfs) )
           ALLOCATE ( folded_onebody(nconfs,nconfs) )      
           ALLOCATE ( first(nconfs,nconfs) );   ALLOCATE ( third(nconfs,nconfs) )      
           ALLOCATE ( corepol(nconfs,nconfs) );   ALLOCATE ( holeladder(nconfs,nconfs) )      
           ALLOCATE ( particleladder(nconfs,nconfs) );   ALLOCATE ( totalunfolded(nconfs,nconfs) )      
           ALLOCATE ( final_v(nconfs,nconfs) );  ALLOCATE (nonhermitian_v(nconfs,nconfs) )
           ALLOCATE ( sp_energy(nconfs,nconfs) ); sp_energy = 0.0_dp
           final_v=0.; q_box=0.0_dp; s_box = 0.0_dp
           folded_onebody = 0.0_dp;   folded_twobody = 0.0_dp
           first = 0.0_dp; corepol = 0.0_dp; holeladder = 0.0_dp; particleladder = 0.0_dp
           third = 0.0_dp; totalunfolded = 0.0_dp;
           !  Note that veff is in principle non-hermitian. It can be hermitized via a similarity
           !  transformation
           !  We compute first the onebody part transformed to a two-body part
           CALL setup_onebody_2_twobody(isospin_z,p_parity,ang_mom, sp_energy,s_box, mspace_configs)
           DO bra_side=1, mspace_configs%number_confs
              sp_state_a=mspace_configs%config_ab(bra_side*2-1)
              sp_state_b=mspace_configs%config_ab(bra_side*2)
              DO ket_side=1, mspace_configs%number_confs  
                 sp_state_c=mspace_configs%config_ab(ket_side*2-1)
                 sp_state_d=mspace_configs%config_ab(ket_side*2)
                 two_body=0.0_dp
                 CALL get_twobody_diagrams(sp_state_a,  &
                      sp_state_b, sp_state_c, sp_state_d, &
                      ang_mom, two_body,dg1, dg2, dg3, dg4, d3rd)
                 q_box(bra_side,ket_side,:) = two_body(:)+s_box(bra_side,ket_side,:)
                 totalunfolded(bra_side,ket_side)=two_body(n_startenergy_veff/2+1)
                 first(bra_side,ket_side)=dg1(n_startenergy_veff/2+1)
                 corepol(bra_side,ket_side)=dg2(n_startenergy_veff/2+1)
                 holeladder(bra_side,ket_side)=dg3(n_startenergy_veff/2+1)
                 particleladder(bra_side,ket_side)=dg4(n_startenergy_veff/2+1)
                 third(bra_side,ket_side)=d3rd(n_startenergy_veff/2+1)
              ENDDO
           ENDDO
           IF ( n_startenergy_veff <= 1 ) THEN
              nonhermitian_v = q_box(:,:,n_startenergy_veff/2+1)
              sp_energy(:,:) = sp_energy(:,:) + s_box(:,:,n_startenergy_veff/2+1) 
              !  make the interaction hermitian
              final_v = nonhermitian_v
              CALL hermitize(nonhermitian_v,sp_energy,nconfs,final_v)
           ELSE    
              CALL vfolded_diagrams(q_box,nconfs,folded_twobody)
              CALL vfolded_diagrams(s_box,nconfs,folded_onebody)
              nonhermitian_v = folded_twobody-folded_onebody
              !  make the interaction hermitian
              CALL hermitize(nonhermitian_v,sp_energy,nconfs,final_v)
           ENDIF
           DO bra_side=1, mspace_configs%number_confs
              sp_state_a=mspace_configs%config_ab(bra_side*2-1)
              sp_state_b=mspace_configs%config_ab(bra_side*2)
              DO ket_side=bra_side, mspace_configs%number_confs
                 sp_state_c=mspace_configs%config_ab(ket_side*2-1)
                 sp_state_d=mspace_configs%config_ab(ket_side*2)
                 CALL pphhmtx(sp_state_a, sp_state_b,sp_state_c, sp_state_d, ang_mom, additional)
                 additional = additional/dij(sp_state_a,sp_state_b)/dij(sp_state_c,sp_state_d)
                 SELECT CASE (type_of_interaction)
                 CASE ('coupled-cluster')
                    WRITE(14,'(7I4,5X,4(F12.6,2X))') &
                         isospin_z,p_parity,2*ang_mom,nmap(sp_state_a)  ,nmap(sp_state_b), nmap(sp_state_c), nmap(sp_state_d),  &
                         final_v(ket_side,bra_side)!, additional(n_startenergy_g+1), additional(n_startenergy_g+2),  &
!                         additional(n_startenergy_g)
                 CASE ('open-diagrams')
!                    WRITE(14+isospin_z,'(14I3,2x,F12.6,2X,4F12.6)') all_orbit%nn(sp_state_a), &
                    WRITE(14+isospin_z,'(14I3,2x,F12.6,2X,2F14.10)') all_orbit%nn(sp_state_a), &
                         all_orbit%ll(sp_state_a), &
                         all_orbit%jj(sp_state_a), &
                         all_orbit%nn(sp_state_b), &
                         all_orbit%ll(sp_state_b), &
                         all_orbit%jj(sp_state_b), &
                         all_orbit%nn(sp_state_c), &
                         all_orbit%ll(sp_state_c), &
                         all_orbit%jj(sp_state_c), &
                         all_orbit%nn(sp_state_d), &
                         all_orbit%ll(sp_state_d), &
                         all_orbit%jj(sp_state_d), &
                         isospin_z*2, ang_mom*2, &                       
                         final_v(ket_side,bra_side)!, additional(n_startenergy_g+1)!, additional(n_startenergy_g+2),  &
!                         additional(n_startenergy_g)
                    WRITE(6,'(14I3,2x,F12.6,2X,7F12.6)') all_orbit%nn(sp_state_a), &
                         all_orbit%ll(sp_state_a), &
                         all_orbit%jj(sp_state_a), &
                         all_orbit%nn(sp_state_b), &
                         all_orbit%ll(sp_state_b), &
                         all_orbit%jj(sp_state_b), &
                         all_orbit%nn(sp_state_c), &
                         all_orbit%ll(sp_state_c), &
                         all_orbit%jj(sp_state_c), &
                         all_orbit%nn(sp_state_d), &
                         all_orbit%ll(sp_state_d), &
                         all_orbit%jj(sp_state_d), &
                         isospin_z*2, ang_mom*2, first(ket_side,bra_side),  corepol(ket_side,bra_side), &
                         holeladder(ket_side,bra_side), particleladder(ket_side,bra_side), &
                         third(ket_side,bra_side), totalunfolded(ket_side,bra_side), final_v(ket_side,bra_side)
                 END SELECT
              ENDDO
           ENDDO
           DEALLOCATE (q_box, s_box)
           DEALLOCATE ( mspace_configs%config_ab) 
           DEALLOCATE( final_v, folded_onebody,folded_twobody,sp_energy,nonhermitian_v)
           DEALLOCATE ( first, third, corepol)
           DEALLOCATE ( holeladder, particleladder, totalunfolded )      
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE( nmap )

END SUBROUTINE twobody_contribution
!
!      Sets up onebody contributions
!
SUBROUTINE onebody_contribution
  USE single_particle_orbits
  USE onebody_diagrams
  USE constants
  USE wave_functions
  IMPLICIT NONE
  REAL(DP), DIMENSION(n_startenergy_veff) :: onebody_diagram_1, &
       onebody_diagram_2,onebody_diagram_3, onebody_diagram_5, &
       onebody_diagram_6,onebody_diagram_7,onebody_diagram_8, &
       onebody_diagram_9,onebody_diagram_4, onebody_diagram_10, &
       onebody_diagram_11,onebody_diagram_12,onebody_diagram_13, &
       onebody_diagram_14,onebody_diagram_15,onebody_diagram_16, &
       onebody_diagram_17,onebody_diagram_18,onebody_diagram_19, &
       onebody_diagram_20,onebody_diagram_21
  REAL(DP), DIMENSION(1,1,n_startenergy_veff) :: sbox
  REAL(DP) :: sum_kin, particle_mass
  REAL(DP), DIMENSION(all_orbit%total_orbits,all_orbit%total_orbits) :: kinenergy
  INTEGER :: a, c, iph, k

  ALLOCATE (  one_body_terms &
       (all_orbit%total_orbits, &
       all_orbit%total_orbits,n_startenergy_veff) )
  ALLOCATE (  one_body_folded &
       (all_orbit%total_orbits, &
       all_orbit%total_orbits) )
  kinenergy = 0.0_dp
  WRITE(6,*) 'a, c,  onebody part, kinetic energy and  total sp energy'
  DO a = 1, all_orbit%total_orbits
     IF (all_orbit%model_space(a) == 'outside') CYCLE         
     DO c= a, all_orbit%total_orbits
        !        CALL one_body_diagram_phase(a,c,phase)
        sbox=0.0_dp
        IF (all_orbit%model_space(c) == 'outside') CYCLE         
        IF(iph(all_orbit%ll(a)) /= iph(all_orbit%ll(c))) CYCLE
        IF ( all_orbit%jj(a) /= all_orbit%jj(c)) CYCLE
        IF(all_orbit%itzp(a) /= all_orbit%itzp(c) ) CYCLE
        IF((itzmin == 1).AND.(all_orbit%itzp(a) == -1)) CYCLE
        IF((itzmax == -1).AND.(all_orbit%itzp(a) == 1) ) CYCLE
        SELECT CASE (physical_system) 
        CASE('nuclear_physics')
           particle_mass = p_mass(all_orbit%itzp(a))
           sum_kin=0.0_dp
           IF ( hf_iterations == 0 ) THEN
              DO k=1,n_rel
                 sum_kin=sum_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))*hol(k,all_orbit%ll(c),all_orbit%nn(c))  &
                      *wra(k)*(ra(k)**4)*0.5*hbarc*hbarc/particle_mass
              ENDDO
           ELSE
              DO k=1,n_rel
                 sum_kin=sum_kin+wave_function(k,a)*wave_function(k,c)*wra(k)*(ra(k)**4)*0.5*hbarc*hbarc/particle_mass
              ENDDO
           ENDIF
        CASE('atomic_physics')
           sum_kin=0.0_dp
           IF ( hf_iterations == 0 ) THEN
              DO k=1,n_rel
                 sum_kin=sum_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))*hol(k,all_orbit%ll(c),all_orbit%nn(c))  &
                      *wra(k)*(ra(k)**4)*0.5
              ENDDO
           ELSE
              DO k=1,n_rel
                 sum_kin=sum_kin+wave_function(k,a)*wave_function(k,c)*wra(k)*(ra(k)**4)*0.5
              ENDDO
           ENDIF
        END SELECT
        SELECT CASE (type_of_renormv)
        CASE ('no-core')
           IF ( a == c) THEN
              sum_kin=2*sum_kin
           ELSE
              sum_kin = 0.0_dp
           ENDIF
        CASE ('v-nrg')
           IF ( a == c) THEN
              sum_kin=2*sum_kin
           ELSE
              sum_kin = 0.0_dp
           ENDIF
        CASE ('vlowk')
           sum_kin = a_factor*sum_kin
        CASE ('v-krg')
           sum_kin = a_factor*sum_kin
        CASE ('g-matrix')
           sum_kin = sum_kin*a_factor
        END SELECT
        kinenergy(c,a) = sum_kin
        onebody_diagram_1 = 0.0_dp
        SELECT CASE ( order_of_interaction )
        CASE ('first')
           CALL diagram_1(a,c,onebody_diagram_1)
           one_body_terms(a,c,:)=onebody_diagram_1(:)
        CASE ('second')
           CALL diagram_1(a,c,onebody_diagram_1)
           CALL diagram_2(a,c,onebody_diagram_2)
           CALL diagram_3(a,c,onebody_diagram_3)
           one_body_terms(a,c,:)= &
                (onebody_diagram_1(:)+ &
                onebody_diagram_2(:)+ &
                onebody_diagram_3(:))
        CASE ('third')
           CALL diagram_1(a,c,onebody_diagram_1)
           CALL diagram_2(a,c,onebody_diagram_2)
           CALL diagram_3(a,c,onebody_diagram_3)
           CALL diagram_4(a,c,onebody_diagram_4)
           CALL diagram_5(a,c,onebody_diagram_5)
           CALL diagram_6(a,c,onebody_diagram_6)
           CALL diagram_7(a,c,onebody_diagram_7)
           CALL diagram_8(a,c,onebody_diagram_8)
           CALL diagram_9(a,c,onebody_diagram_9)
           CALL diagram_10(a,c,onebody_diagram_10)
           CALL diagram_11(a,c,onebody_diagram_11)
           CALL diagram_12(a,c,onebody_diagram_12)
           CALL diagram_13(a,c,onebody_diagram_13)
           CALL diagram_14(a,c,onebody_diagram_14)
           CALL diagram_15(a,c,onebody_diagram_15)
           CALL diagram_16(a,c,onebody_diagram_16)
           CALL diagram_17(a,c,onebody_diagram_17)
           CALL diagram_18(a,c,onebody_diagram_18)
           CALL diagram_19(a,c,onebody_diagram_19)
           CALL diagram_20(a,c,onebody_diagram_20)        
           CALL diagram_21(a,c,onebody_diagram_21)
           one_body_terms(a,c,:)= &
                (onebody_diagram_1(:)+ &
                onebody_diagram_2(:)+ &
                onebody_diagram_3(:)+onebody_diagram_4(:)+onebody_diagram_5(:)+ &
                onebody_diagram_6(:)+onebody_diagram_7(:)+onebody_diagram_8(:)+ &
                onebody_diagram_9(:)+onebody_diagram_10(:)+onebody_diagram_11(:)+ &
                onebody_diagram_12(:)+onebody_diagram_13(:)+onebody_diagram_14(:)+ &
                onebody_diagram_15(:)+onebody_diagram_16(:)+onebody_diagram_17(:)+ &
                onebody_diagram_18(:)+onebody_diagram_19(:)+onebody_diagram_20(:)+ &
                onebody_diagram_21(:))
        END SELECT
        one_body_terms(c,a,:)=one_body_terms(a,c,:)
        sbox(1,1,:)=one_body_terms(a,c,:)

        WRITE(6,'(2I4,2X,3F12.6)') a , c, one_body_terms(a,c,n_startenergy_veff/2+1), &
             kinenergy(c,a),one_body_terms(a,c,n_startenergy_veff/2+1)+kinenergy(c,a)
     ENDDO
  ENDDO

END SUBROUTINE onebody_contribution
!
!     All two body diagrams through third order in perturbation
!     theory can be computed here. The order of the interaction is
!     specified by the variable -order_of_interaction-
!     The terms are normalized
!
SUBROUTINE get_twobody_diagrams(a,b,c,d,jtot,twobd,d1, d2, d3, d4, third)
  USE constants
  IMPLICIT NONE
  REAL(DP) :: phase_ab, phase_cd, phase_fact, norm , dij
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(INOUT) :: twobd, third, d1, d2, d3, d4
  REAL(DP), DIMENSION(n_startenergy_veff) :: dg1,dg2,dg3,& 
       dg4, dsum1,dsum2,dsum4, dg10,dg22,dg23a, dg23b,&
       dsum5,dg11a,dg11b,drpa, second
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: i

  norm=1./dij(a,b)/dij(c,d)
  twobd=0.; dg1 = 0.; dg2 = 0.; dg3 = 0.; dg4 = 0.;
  dsum1 = 0.; dsum2 = 0.; dsum4= 0.; drpa=0.; dsum5 =0.;
  dg11a = 0.; dg11b = 0.; dg22 = 0.; dg23a = 0.; dg23b= 0.
  second = 0.0_dp; third = 0.0_dp
  CALL diag1(a,b,c,d,jtot,dg1)
  SELECT CASE ( order_of_interaction)
  CASE ('second')
     CALL core_polarization_sum(a,b,c,d,jtot,dg2)
     CALL hole_hole_ladder(a,b,c,d,jtot,dg3)
     CALL particle_particle_ladder(a,b,c,d,jtot,dg4)
     second = dg2+dg3+dg4
  CASE ( 'third') 
     CALL core_polarization_sum(a,b,c,d,jtot,dg2)
     CALL hole_hole_ladder(a,b,c,d,jtot,dg3)
     CALL particle_particle_ladder(a,b,c,d,jtot,dg4)
     second = dg2+dg3+dg4
     CALL number_conserving_set(a,b,c,d,jtot,dsum1)
     CALL diag10(a,b,c,d,jtot,dg10)
     CALL diag11a(a,b,c,d,jtot,dg11a)
     CALL diag11b(a,b,c,d,jtot,dg11b)
     CALL ladder_corepol(a,b,c,d,jtot,dsum2)
     CALL tda_rpa_diagrams(a,b,c,d,jtot,drpa)
     CALL screening_thirdorder(a,b,c,d,jtot,dsum4)
     CALL diag22(a,b,c,d,jtot,dg22)
     CALL diag23(a,b,c,d,jtot,dg23a,dg23b)
     CALL onebody_insertions(a,b,c,d,jtot,dsum5)
     third=dg10+dsum1+dsum2+dsum4+ &
          dg22+dg23a+dg23b+dsum5+dg11a+dg11b+drpa
  END SELECT
  twobd = (dg1+second+third)*norm
  third = third*norm
  d1 = dg1*norm; d2 = dg2*norm; d3 = dg3*norm; d4 = dg4*norm

END SUBROUTINE get_twobody_diagrams
!
!     Using an unfolded one-body, or two-body or three-body set
!     of diagrams, this routine sets up the 
!     contribution from folded diagrams at a given starting energy 
!     This function is taylored to a degenerate model space and uses
!     the standard approach of iterating with Q-box derivatives
!     For references see M. Hjorth-Jensen et al, Physics Reports 261 (1995) 125
!     It yields a non-hermitian effective interaction.
!
SUBROUTINE vfolded_diagrams(unfolded_diagrams,nconfs,final_v)
  USE constants
  USE configurations
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nconfs
  INTEGER :: number_of_derivatives, number_of_iterations
  INTEGER :: bra, ket,i 
  REAL(DP), DIMENSION(nconfs,nconfs), INTENT(INOUT) :: final_v
  REAL(DP), DIMENSION(nconfs,nconfs,n_startenergy_veff), INTENT(IN) :: unfolded_diagrams
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: qsum, qq, aa, veff, tempq, xqsum
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:):: dq

  ALLOCATE(dq(nconfs,nconfs,0:n_startenergy_veff-1))
  ALLOCATE(qsum(nconfs,nconfs),qq(nconfs,nconfs),aa(nconfs,nconfs),veff(nconfs,nconfs),&
       tempq(nconfs,nconfs),xqsum(nconfs,nconfs))
  !     calculate first the derivatives of the unfolded diagrams
  CALL derivate_nonfolded(unfolded_diagrams,dq,nconfs)
  qsum=0.0_dp
  veff(:,:)=unfolded_diagrams(:,:,n_startenergy_veff/2+1)
  final_v=0.0_dp
  DO number_of_iterations=1,n_startenergy_veff-1
     xqsum=0.0_dp
     DO number_of_derivatives=1,n_startenergy_veff-2                    
        aa(:,:)=dq(:,:,number_of_derivatives)          
        IF(number_of_derivatives == 2) THEN
           qq=veff
        ELSEIF(number_of_derivatives >= 3) THEN
           qq=tempq
        ENDIF
        IF (number_of_derivatives == 1) THEN
           !     Q-box * first der. of Q 
           qsum=MATMUL(aa,veff)
        ELSEIF (number_of_derivatives >= 2) THEN
           !     Prods of Q-boxes
           tempq=MATMUL(qq,veff)
           !     sec. der or higher der * products of q-boxes
           qsum=MATMUL(aa,tempq)
        ENDIF
        xqsum=xqsum+qsum
     ENDDO
     IF(number_of_iterations == 1) THEN
        veff=veff+xqsum
     ELSE
        veff(:,:)=unfolded_diagrams(:,:,n_startenergy_veff/2+1)+xqsum(:,:)
     ENDIF
     final_v=veff
  ENDDO
  DEALLOCATE(qsum,qq,aa,veff,tempq,xqsum); DEALLOCATE(dq)

END SUBROUTINE vfolded_diagrams
!
!     Using an unfolded or folded one-body, or two-body or three-body set
!     of diagrams, this routine sets up a
!     hermitian effective interaction.
!     This function is taylored to a degenerate model space and uses
!     the standard approach of iterating with Q-box derivatives
!     For references see M. Hjorth-Jensen et al, Physics Reports 261 (1995) 125
!     
!
SUBROUTINE hermitize(nonhermitian_v,sp_energy,nconfs,final_v)
  USE constants
  USE configurations
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nconfs
  INTEGER :: bra, ket,i 
  REAL(DP), DIMENSION(nconfs,nconfs), INTENT(INOUT) :: final_v
  REAL(DP), DIMENSION(nconfs,nconfs), INTENT(IN) :: nonhermitian_v, sp_energy
  COMPLEX*16, ALLOCATABLE :: cvec_pp(:,:)
  COMPLEX*16, ALLOCATABLE :: cvec_qp(:,:), heff(:,:)
  COMPLEX*16, ALLOCATABLE::  ceig_p(:) 
  COMPLEX*16, ALLOCATABLE :: gna(:,:)
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) ::  cspenergy

  ! make the interaction hermitian through the Lee-Suzuki similarity
  ! transformation
  ALLOCATE(cvec_pp(nconfs, nconfs))
  ALLOCATE(cvec_qp(0,nconfs)); ALLOCATE(gna(nconfs, nconfs))
  ALLOCATE(ceig_p(nconfs)); ALLOCATE(heff(nconfs, nconfs))
  ALLOCATE ( cspenergy(nconfs,nconfs) ); cspenergy = (0.0_dp,0.0_dp)
  cspenergy = sp_energy
  gna = nonhermitian_v+cspenergy
  heff = (0.0_dp,0.0_dp)
  !      Get model space eigenvalues, eigenvectors of 
  !      model space and excluded space and
  !      model space configurations which match the corresponding ones
  !      of the large space
  CALL eigenvalues_large(cvec_pp,cvec_qp,ceig_p,gna,nconfs,nconfs)
  !      setup 2p-effective interaction in P-space using
  !      the Lee-Suzuki similarity transformation
  CALL lee_suzuki( cvec_pp, cvec_qp, ceig_p, nconfs, nconfs-nconfs, heff, cspenergy )
  final_v   = heff
  !      free space
  DEALLOCATE(cvec_pp); DEALLOCATE(cvec_qp); DEALLOCATE(ceig_p)
  DEALLOCATE(gna); DEALLOCATE(heff)

END SUBROUTINE hermitize
!
!     This routine calculates the derivative of the x-box
!
SUBROUTINE derivate_nonfolded(q,dq,nconfs)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nconfs
  INTEGER :: bra, ket, ien
  REAL(DP), DIMENSION(nconfs,nconfs,n_startenergy_veff), INTENT(IN) :: q
  REAL(DP), DIMENSION(nconfs,nconfs,0:n_startenergy_veff-1),INTENT(OUT):: dq
  REAL(DP), DIMENSION(n_startenergy_veff) :: en,x,dx1,&
       dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,fct
  REAL(DP) :: der, w 

  IF ( n_startenergy_veff <= 1) RETURN
  DO ien=1, n_startenergy_veff
     en(ien)=starting_energy-(2.*ien-1.)+n_startenergy_veff
  ENDDO
  fct(1)=1.
  DO ien=2, n_startenergy_veff
     fct(ien)=fct(ien-1)/FLOAT(ien)
  ENDDO
  dq=0.0_dp
  dq(:,:,0)=q(:,:,n_startenergy_veff/2+1)
  DO bra=1, nconfs
     DO ket=bra, nconfs
        x(:)=q(bra,ket,:)
        IF ( x(n_startenergy_veff/2+1) == 0. ) CYCLE 
        dx1=0.0_dp; dx2=0.0_dp;dx3=0.0_dp;dx4=0.0_dp;dx5=0.0_dp;dx6=0.0_dp;dx7=0.0_dp
        dx8=0.0_dp;dx9=0.0_dp
        !     First derivative
        DO ien=1, n_startenergy_veff
           w=starting_energy-ien*1.d0+n_startenergy_veff/2.0_dp+1.0_dp
           CALL deriv(w,en,x,der)
           dx1(ien)=der
        ENDDO
        dq(bra,ket,1)=dx1(n_startenergy_veff/2+1)
        dq(ket,bra,1)=dq(bra,ket,1)
        IF ( n_startenergy_veff <= 2) CYCLE
        !     second derivative
        w=starting_energy+4.999d0
        DO ien=1, n_startenergy_veff
           CALL deriv(w,en,dx1,der)
           dx2(ien)=der
           w=w-0.9998d0
        ENDDO
        dq(bra,ket,2)=dx2(n_startenergy_veff/2+1)*fct(2)
        dq(ket,bra,2)=dq(bra,ket,2)
        IF ( n_startenergy_veff <= 3) CYCLE
        !     third derivative
        w=starting_energy+4.998d0
        DO ien=1, n_startenergy_veff
           CALL deriv(w,en,dx2,der)
           dx3(ien)=der
           w=w-0.9996d0
        ENDDO
        dq(bra,ket,3)=dx3(n_startenergy_veff/2+1)*fct(3)
        dq(ket,bra,3)=dq(bra,ket,3)
        IF ( n_startenergy_veff <= 4) CYCLE
        !     fourth derivative
        w=starting_energy+4.995d0
        DO ien=1, n_startenergy_veff
           CALL deriv(w,en,dx3,der)
           dx4(ien)=der
           w=w-0.999d0
        ENDDO
        dq(bra,ket,4)=dx4(n_startenergy_veff/2+1)*fct(4)
        dq(ket,bra,4)=dq(bra,ket,4)
        IF ( n_startenergy_veff <= 5) CYCLE
        !     fifth derivative
        w=starting_energy+4.990d0
        DO ien=1, n_startenergy_veff
           CALL deriv(w,en,dx4,der)
           dx5(ien)=der
           w=w-0.998d0
        ENDDO
        dq(bra,ket,5)=dx5(n_startenergy_veff/2+1)*fct(5)
        dq(ket,bra,5)=dq(bra,ket,5)
        IF ( n_startenergy_veff <= 6) CYCLE
        !     sixth derivative
        w=starting_energy+4.988d0
        DO ien=1, n_startenergy_veff
           CALL deriv(w,en,dx5,der)
           dx6(ien)=der
           w=w-0.9976d0
        ENDDO
        dq(bra,ket,6)=dx6(n_startenergy_veff/2+1)*fct(6)
        dq(ket,bra,6)=dq(bra,ket,6)
        IF ( n_startenergy_veff <= 7) CYCLE
        !     seventh derivative
        w=starting_energy+4.986d0
        DO ien=1, n_startenergy_veff
           CALL deriv(w,en,dx6,der)
           dx7(ien)=der
           w=w-0.9972d0
        ENDDO
        dq(bra,ket,7)=dx7(n_startenergy_veff/2+1)*fct(7)
        dq(ket,bra,7)=dq(bra,ket,7)
        IF ( n_startenergy_veff <= 8) CYCLE
        !     eight derivative
        w=starting_energy+4.984d0
        DO ien=1, n_startenergy_veff
           CALL deriv(w,en,dx7,der)
           dx8(ien)=der
           w=w-0.9968d0
        ENDDO
        dq(bra,ket,8)=dx8(n_startenergy_veff/2+1)*fct(8)
        dq(ket,bra,8)=dq(bra,ket,8)
        IF ( n_startenergy_veff <= 9) CYCLE
        !     ninth derivative
        w=starting_energy+4.982d0
        DO ien=1, n_startenergy_veff
           CALL deriv(w,en,dx8,der)
           dx9(ien)=der
           w=w-0.9964d0
        ENDDO
        dq(bra,ket,9)=dx9(n_startenergy_veff/2+1)*fct(9)
        dq(ket,bra,9)=dq(bra,ket,9)
     ENDDO
  ENDDO

END SUBROUTINE derivate_nonfolded
!
!     Compute derivative
!
SUBROUTINE deriv(x,w,f,der)
  USE constants
  IMPLICIT NONE   
  REAL(DP) :: x,w,f,p_ip,p_i,elp, el,der,a,b
  INTEGER :: j, k, i
  DIMENSION w(n_startenergy_veff), f(n_startenergy_veff), & 
       p_ip(n_startenergy_veff), elp(n_startenergy_veff), &
       el(n_startenergy_veff)

  el=0.0_dp; elp=0.0_dp ; p_ip=1.0_dp; p_i=1.0_dp; der=0.0_dp
  DO k=1,n_startenergy_veff
     IF (x /= w(k)) p_i=p_i*(x-w(k))
     p_ip(k)=PRODUCT(w(k)-w, MASK=(w(k) /= w) )
  ENDDO
  outer_loop :  DO k=1,n_startenergy_veff
     b=0.
     DO j=1,n_startenergy_veff
        IF (j == k) CYCLE 
        IF (x == w(k)) THEN 
           a=(x-w(j))*p_ip(k)
           b=b+1./a
        ELSEIF (x == w(j)) THEN
           elp(k)=p_i/((x-w(k))*p_ip(k)) ; CYCLE outer_loop
        ELSE
           a=(x-w(k))*(x-w(j))*p_ip(k)
           b=b+1.0_dp/a
        ENDIF
     ENDDO
     elp(k)=p_i*b
  ENDDO outer_loop
  DO i=1, n_startenergy_veff
     der=der+elp(i)*f(i)
  ENDDO

END SUBROUTINE deriv
!
!  Setup single-particle energies and single-particle s_box and transform it to 
!  two-body form
!
SUBROUTINE setup_onebody_2_twobody(it,ip,ij, sp_energy,s_box, gmatrix_configs)
  USE single_particle_orbits
  USE onebody_diagrams
  USE constants
  USE wave_functions
  USE configurations
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(in)  :: gmatrix_configs
  INTEGER :: i,j, ijd, it, ip, ij, ia, ib, ic,id,iabcd,ie,iab, icd, itemp
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(INOUT)  :: sp_energy
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_veff), INTENT(INOUT)  :: s_box
  REAL(DP), DIMENSION(n_startenergy_veff) :: one_1
  REAL(DP) :: upot1, upot2, upot3, upot4, dij, particle_mass
  REAL(DP), DIMENSION(n_startenergy_veff) :: onebd1, onebd2, onebd3, onebd4, s_dir, s_ex
  REAL(DP) :: ekin1, ekin2, ekin3, ekin4, e_dir, e_ex, sbox1, sbox2, sbox3, sbox4
  INTEGER :: iph, n_confs, a,b,c,d,na,la,ja,tza,nb,lb,jb,tzb,nc,lc,jc,tzc,nd,ld,jd,tzd, i1, j1, k

  n_confs = gmatrix_configs%number_confs
  ! setup single-particle energy and two-body S-box
  sp_energy = 0.0_dp; s_box = 0.0_dp  
  DO i1=1,gmatrix_configs%number_confs
     a=gmatrix_configs%config_ab(i1*2-1)
     b=gmatrix_configs%config_ab(i1*2)
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)
     nb=all_orbit%nn(b)
     lb=all_orbit%ll(b)
     jb=all_orbit%jj(b)
     tzb = all_orbit%itzp(b)
     DO j1=i1,gmatrix_configs%number_confs
        c=gmatrix_configs%config_ab(j1+j1-1)
        d=gmatrix_configs%config_ab(j1+j1)
        nc=all_orbit%nn(c) 
        lc=all_orbit%ll(c)
        jc=all_orbit%jj(c)
        tzc = all_orbit%itzp(c)
        nd=all_orbit%nn(d)
        ld=all_orbit%ll(d)
        jd=all_orbit%jj(d)
        tzd = all_orbit%itzp(d)
        ! term 1
        ekin1 = 0.0_dp; ekin2 = 0.0_dp; ekin3 = 0.0_dp ; ekin4 = 0.0_dp
        upot1 = 0.0_dp; upot2 = 0.0_dp; upot3 = 0.0_dp ; upot4 = 0.0_dp
        onebd1 = 0.0_dp; onebd2 = 0.0_dp; onebd3 = 0.0_dp ; onebd4 = 0.0_dp
        IF ( b == d .AND. la == lc .AND. ja == jc .AND. tza == tzc ) THEN
           IF (hf_iterations == 0) THEN
              IF ( na == nc-1) ekin1 = SQRT(nc*(nc+lc+0.5))
              IF ( na == nc)   ekin1 =  2.*nc+lc+1.5
              IF ( na == nc+1) ekin1 = SQRT((nc+1.)*(nc+lc+1.5))
              ekin1 = ekin1*hbar_omega*0.50_dp
           ELSE
              DO k=1,n_rel
                 ekin1=ekin1+wave_function(k,a)*wave_function(k,c)*wra(k)*(ra(k)**4)
              ENDDO
              SELECT CASE (physical_system) 
              CASE('nuclear_physics')
                 particle_mass = p_mass(all_orbit%itzp(a))
              CASE('atomic_physics')
                 particle_mass = e_mass
              END SELECT
              ekin1 = ekin1*0.5*hbarc*hbarc/particle_mass
           ENDIF
           onebd1(:)=one_body_terms(a,c,:)
           CALL diagram_1(a,c,one_1)
           upot1=one_1(n_startenergy_veff/2+1)
        ENDIF
        ! term 2
        IF ( a == c .AND. lb == ld .AND. jb == jd .AND. tzb == tzd ) THEN
           IF (hf_iterations == 0) THEN
              IF ( nb == nd-1) ekin2 = SQRT(nd*(nd+ld+0.5))
              IF ( nb == nd)   ekin2 =  2.*nd+ld+1.5
              IF ( nb == nd+1) ekin2 = SQRT((nd+1.)*(nd+ld+1.5))
              ekin2 = ekin2*hbar_omega*0.50_dp
           ELSE
              DO k=1,n_rel
                 ekin2=ekin2+wave_function(k,b)*wave_function(k,d)*wra(k)*(ra(k)**4)
              ENDDO
              SELECT CASE (physical_system) 
              CASE('nuclear_physics')
                 particle_mass = p_mass(all_orbit%itzp(b))
              CASE('atomic_physics')
                 particle_mass = e_mass
              END SELECT
              ekin2 = ekin2*0.5*hbarc*hbarc/particle_mass
           ENDIF
           onebd2(:)=one_body_terms(b,d,:)
           CALL diagram_1(b,d,one_1)
           upot2=one_1(n_startenergy_veff/2+1)
        ENDIF
        ! term 3
        IF ( b == c .AND. la == ld .AND. ja == jd .AND. tza == tzd ) THEN
           IF (hf_iterations == 0) THEN
              IF ( na == nd-1) ekin3 = SQRT(nd*(nd+ld+0.5))
              IF ( na == nd)   ekin3 =  2.*nd+ld+1.5
              IF ( na == nd+1) ekin3 = SQRT((nd+1.)*(nd+ld+1.5))
              ekin3 = ekin3*hbar_omega*0.50_dp
           ELSE
              DO k=1,n_rel
                 ekin3=ekin3+wave_function(k,a)*wave_function(k,d)*wra(k)*(ra(k)**4)
              ENDDO
              SELECT CASE (physical_system) 
              CASE('nuclear_physics')
                 particle_mass = p_mass(all_orbit%itzp(a))
              CASE('atomic_physics')
                 particle_mass = e_mass
              END SELECT
              ekin3 = ekin3*0.5*hbarc*hbarc/particle_mass
           ENDIF
           onebd3(:)=one_body_terms(a,d,:)
           CALL diagram_1(a,d,one_1)
           upot3=one_1(n_startenergy_veff/2+1)
        ENDIF
        ! term 4
        IF ( a == d .AND. lb == lc .AND. jb == jc .AND. tzb == tzc ) THEN
           IF (hf_iterations == 0) THEN
              IF ( nb == nc-1) ekin4 = SQRT(nc*(nc+lc+0.5))
              IF ( nb == nc)   ekin4 =  2.*nc+lc+1.5
              IF ( nb == nc+1) ekin4 = SQRT((nc+1.)*(nc+lc+1.5))
              ekin4 = ekin4*hbar_omega*0.50_dp
           ELSE
              DO k=1,n_rel
                 ekin4=ekin4+wave_function(k,b)*wave_function(k,c)*wra(k)*(ra(k)**4)
              ENDDO
              SELECT CASE (physical_system) 
              CASE('nuclear_physics')
                 particle_mass = p_mass(all_orbit%itzp(b))
              CASE('atomic_physics')
                 particle_mass = e_mass
              END SELECT
              ekin4 = ekin4*0.5*hbarc*hbarc/particle_mass
           ENDIF
           onebd4(:)=one_body_terms(b,c,:)
           CALL diagram_1(b,c,one_1)
           upot4=one_1(n_startenergy_veff/2+1)
        ENDIF
        !  first sp energies

        upot1= 0.0_dp; upot2= 0.0_dp; upot3= 0.0_dp; upot4= 0.0_dp;
        e_dir = 0.0_dp; e_ex = 0.0_dp
        e_dir = (a_factor*(ekin1 + ekin2) + upot1 + upot2)/dij(a,b)/dij(c,d) 
        e_ex = (a_factor*(ekin3 + ekin4) + upot3 + upot4 )/dij(a,b)/dij(c,d) 
        IF( it /= 0 ) sp_energy(i1,j1) = ( e_dir -iph( (2*ij-jc-jd)/2)*e_ex )
        IF( it == 0 ) sp_energy(i1,j1) = ( a_factor*(ekin1 + ekin2)+upot1+upot2 ) 
        sp_energy(j1,i1) = sp_energy(i1,j1) 
        ! then s_box with subtracted HF term
        s_dir = 0.0_dp; s_ex = 0.0_dp
        s_dir(:) = (onebd1(:) + onebd2(:) - upot1 - upot2)/dij(a,b)/dij(c,d) 
        s_ex = (onebd3 + onebd4(:) - upot3 - upot4 )/dij(a,b)/dij(c,d) 
        ! p_1^2 + p_2^2        
        IF( it /= 0 ) s_box(i1,j1,:) = ( s_dir(:) -iph( (2*ij-jc-jd)/2)*s_ex(:) )
        IF( it == 0 ) s_box(i1,j1,:) = ( onebd1(:)+onebd2(:) -upot1 -upot2 ) 
     ENDDO
  ENDDO

END SUBROUTINE setup_onebody_2_twobody
