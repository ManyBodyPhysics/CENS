!             Program block bhf-matrix.f
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!
!                    Begin bhf-matrix code
!
SUBROUTINE  setup_hfmatrix
  USE constants
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  IMPLICIT NONE
  INTEGER :: number_orbits
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: coeffs  

  !     reserve space in memory for various arrays
  ALLOCATE ( ra (n_rel), wra (n_rel));  ALLOCATE ( rgkk (n_rel), wgkk (n_rel))
  ALLOCATE ( hol (n_rel, 0:lmax, 0:nmax) )
  !     set up mesh points in lab frame
  CALL rel_mesh                   
  !     setup ho wave functions
  CALL ho_wfunction                
  !    Expansion coefficients for sp wave functions. 
  ALLOCATE(coeffs(all_orbit%total_orbits,all_orbit%total_orbits))
  ALLOCATE ( bhf_hol (n_rel,all_orbit%total_orbits))
  ALLOCATE ( wave_function (n_rel,all_orbit%total_orbits))
  coeffs = 0.0_dp; wave_function = 0.0_dp
  !  perform the HF calculation
  CALL brueckner_hartree_fock(coeffs,all_orbit%total_orbits)
  !  update the g-matrix
  CALL setupg_bhf(coeffs,all_orbit%total_orbits)
  SELECT CASE (type_of_interaction)
  CASE ('core-diagrams')
     CALL core_diagrams
  CASE ('coupled-cluster')
     CALL ccsd_data_model
     SELECT CASE ( order_of_interaction )
     CASE ('first')
        CALL coupled_cluster_twobody
     CASE ('second')
        CALL onebody_contribution
        CALL twobody_contribution
     CASE ('third')
        CALL onebody_contribution
        CALL twobody_contribution
     END SELECT
  CASE ('open-diagrams')
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
  END SELECT
  DEALLOCATE(coeffs) ;    DEALLOCATE ( bhf_hol,wave_function )
  DEALLOCATE ( rgkk, wgkk, ra, wra) ; DEALLOCATE ( hol)

END SUBROUTINE setup_hfmatrix
!
!           Set up the BHF G-mtx in the lab-frame 
!
SUBROUTINE setupg_bhf(coeffs,ncoeffs)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE hoosc_gmatrix
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs
  REAL(DP), ALLOCATABLE :: bhf_coeff(:,:)
  INTEGER, INTENT(IN) :: ncoeffs
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  INTEGER ::  p_parity, ang_mom, isospin_z, n_confs, ie, bra, ket
  REAL(DP), ALLOCATABLE :: gna(:,:,:),  temp(:,:,:)

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
           ALLOCATE(temp(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(gna(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(bhf_coeff(n_confs, n_confs))
           gna=0.0_dp; bhf_coeff = 0.0_dp; temp = 0.0_dp
           CALL fetch_matrix(isospin_z,p_parity,ang_mom,gna,gmatrix_configs)
           CALL bhf_coefficients(ang_mom,isospin_z,n_confs,ncoeffs,gmatrix_configs,bhf_coeff,coeffs)
           DO ie=1,n_startenergy_g
              temp(:,:,ie) = MATMUL(gna(:,:,ie),TRANSPOSE(bhf_coeff(:,:)))
              gna(:,:,ie) = MATMUL(bhf_coeff(:,:),temp(:,:,ie))
           ENDDO
           !     update table of matrix elements 
           CALL update(isospin_z,p_parity,ang_mom,gna,gmatrix_configs)
           !     free space in heap
           DEALLOCATE(gmatrix_configs%config_ab)
           DEALLOCATE(bhf_coeff)
           DEALLOCATE(gna, temp)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE setupg_bhf
!
!                    Update out G-mtx
!
SUBROUTINE update(it,ip,ij,gna,gmatrix_configs)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  USE stored_bare_interaction
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, ijd, it, ip, ij, ia, ib, ic,id,ie
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_g), INTENT(IN)  :: gna

  ijd=ij+ij
  DO i=1,gmatrix_configs%number_confs
     ia=gmatrix_configs%config_ab(i*2-1)
     ib=gmatrix_configs%config_ab(i*2)
     DO j=i,gmatrix_configs%number_confs
        ic=gmatrix_configs%config_ab(j+j-1)
        id=gmatrix_configs%config_ab(j+j)
        CALL replace_g(ia,ib, ic, id, it, ip, ij ,gna(j,i,:))
        WRITE(9,'(7I4,5X,10(5X,E12.6))') it, ip, ijd, ia, ib, ic, id, (gna(j,i,ie),ie=1,n_startenergy_g)
     ENDDO
  ENDDO

END SUBROUTINE update
!
!                    Get just the G-mtx to be used in various HF iterations
!
SUBROUTINE fetch_matrix(it,ip,ij,gna,gmatrix_configs)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  USE stored_bare_interaction
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, it, ip, ij, ia, ib, ic,id
  REAL(DP) :: norm, dij
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_g), INTENT(INOUT)  :: gna
  REAL(DP), DIMENSION(n_startenergy_g) :: ans

  DO i=1,gmatrix_configs%number_confs
     ia=gmatrix_configs%config_ab(i*2-1)
     ib=gmatrix_configs%config_ab(i*2)
     DO j=i,gmatrix_configs%number_confs
        ic=gmatrix_configs%config_ab(j+j-1)
        id=gmatrix_configs%config_ab(j+j)
        CALL pphhmtx(ia,ib,ic,id,ij,ans)
        norm=1.0_dp/dij(ia,ib)/dij(ic,id)
        gna(j,i,:)=ans(:)*norm
        gna(i,j,:) = ans(:)*norm
     ENDDO
  ENDDO

END SUBROUTINE fetch_matrix
!
!                 Set up h.o. wf for rel cm system and lab frame
!                 It computes also a complex
!
SUBROUTINE ho_wfunction
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER :: n, l, i, j
  REAL(DP) :: ph, sum_rel
  REAL(DP)  :: cx(0:200), factor, z_lab,xp

  DO n=0,nmax
     ph=(-1.D0)**n
     DO l=0,lmax
        factor = 0.5D0*((n+l+2)*LOG(2.D0)+fac(n)-dfac(2*n+2*l+1)-0.5D0*LOG(pi))
        factor = EXP(factor)
        sum_rel=0.0_dp
        DO i=1,n_rel
           !  real ho wave function
           z_lab= ra(i)*ra(i)*oscl*oscl; cx = 0.0_dp
           CALL laguerre_general( n, l+0.5D0, z_lab, cx )
           xp = EXP(-z_lab*0.5D0)*((ra(i)*oscl)**l)*cx(n)
           hol(i,l,n) = xp*ph*factor*(oscl**(1.5D0)) ! lab wf
           sum_rel=sum_rel+ wra(i)*(hol(i,l,n)*ra(i))**2
        ENDDO
        WRITE(6,'(21H Norm cm ho wf n,l : ,2I3,2X,3F12.7)') n, l, sum_rel
     ENDDO
  ENDDO

END SUBROUTINE ho_wfunction
!
!     Brueckner-Hartree-Hock self-consistency
!
SUBROUTINE brueckner_hartree_fock(coeffs, ncoeffs)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: sp_configs
  INTEGER, INTENT(IN) :: ncoeffs
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(INOUT)  :: coeffs  
  INTEGER :: a, c, ket , bra, nrot, hole_states, hf_iter, max_hfiter
  INTEGER :: la, lc, na, nc, k, i
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: hf_vect
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: hf_mtx
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: hf_eigen
  REAL(DP), DIMENSION(n_rel) ::  sum_wf
  REAL(DP) :: e_kin, hf, sigma
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: int_factor, kin_energy
  ALLOCATE (int_factor(n_rel), kin_energy(n_rel) )
  int_factor(:)=wra(:)*ra(:)*ra(:)

  !     initialize bhf harmonic oscillator sp wave function
  bhf_hol=0.0_dp; max_hfiter = 100
  ALLOCATE (hf_mtx(all_orbit%total_orbits,all_orbit%total_orbits))
  ALLOCATE (hf_vect(all_orbit%total_orbits,all_orbit%total_orbits))
  ALLOCATE (hf_eigen(all_orbit%total_orbits))
  hf_vect = 0.0_dp
  !  for first iteration coeffs has only diagonal values equal 1
  DO bra = 1, all_orbit%total_orbits
     coeffs(bra,bra) = 1.0_dp
  ENDDO
  hf_iter = 0; sigma = 1.0_dp
  DO WHILE ((hf_iter <= max_hfiter).AND.(ABS(sigma) > 1E-8))
     !     set up bhf matrix to be diagonalized
     hf_mtx = 0.0_dp; hf_vect = 0.0_dp
     DO bra = 1, all_orbit%total_orbits
        a = bra
        DO ket= bra, all_orbit%total_orbits
           c = ket
           IF(all_orbit%ll(a) /= all_orbit%ll(c)) CYCLE
           IF ( all_orbit%jj(a) /= all_orbit%jj(c)) CYCLE
           IF(all_orbit%itzp(a) /= all_orbit%itzp(c) ) CYCLE
           hf = 0.0_dp               
           SELECT CASE (physical_system) 
           CASE('nuclear_physics')
              kin_energy(:)=0.5*ra(:)*ra(:)*hbarc*hbarc/p_mass(all_orbit%itzp(c))
           CASE('atomic_physics')
              kin_energy(:)=0.5*ra(:)*ra(:)
           END SELECT
           CALL  diagram_HF(a, c, coeffs, hf_iter, hf, ncoeffs)
           !  compute the kinetic energy or unperturbed one-body H
           e_kin = 0.0_dp
           SELECT CASE (type_of_renormv)
           CASE ('no-core')
              IF ( a == c) THEN
                 DO k=1,n_rel
                    e_kin=e_kin+2*hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                         hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                         int_factor(k)*kin_energy(k)
                 ENDDO
              ENDIF
           CASE ('v-nrg')
              IF ( a == c) THEN
                 DO k=1,n_rel
                    e_kin=e_kin+2*hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                         hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                         int_factor(k)*kin_energy(k)
                 ENDDO
              ENDIF
           CASE ('vlowk')
              DO k=1,n_rel
                 e_kin=e_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                      hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                      int_factor(k)*a_factor*kin_energy(k)
              ENDDO
           CASE ('v-krg')
              DO k=1,n_rel
                 e_kin=e_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                      hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                      int_factor(k)*a_factor*kin_energy(k)
              ENDDO
           CASE ('g-matrix')
              DO k=1,n_rel
                 e_kin=e_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                      hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                      int_factor(k)*a_factor*kin_energy(k)
              ENDDO
           END SELECT
           hf_mtx(bra,ket) = e_kin+hf
           hf_mtx(ket,bra) = hf_mtx(bra,ket)
        ENDDO
     ENDDO
     !     obtain the BHF coefficients 
     CALL matrix_diag(hf_mtx,all_orbit%total_orbits, all_orbit%total_orbits,&
          hf_eigen, hf_vect,nrot)
     !     set up the new bhf harmonic  oscillator wave function
     !     Note memory stride
     sigma = 0.0_dp
     DO bra = 1, all_orbit%total_orbits
        a = bra
        la = all_orbit%ll(a)
        na = all_orbit%nn(a)
        sum_wf = 0.0_dp
        DO ket= 1, all_orbit%total_orbits
           c = ket
           lc = all_orbit%ll(c)
           nc = all_orbit%nn(c)
           sum_wf(:)=sum_wf(:)+hf_vect(ket,bra)*hol(:,lc,nc)
           coeffs(a,c) = hf_vect(bra,ket)
        ENDDO
        bhf_hol(:,a)=sum_wf(:)
        ! set up of new single-particle energy 
        sigma = sigma +ABS(all_orbit%e(a) - hf_eigen(bra))
        all_orbit%e(a) = hf_eigen(bra)                    
     ENDDO
     sigma = sigma/all_orbit%total_orbits
     WRITE(6,*) 'iteration nr and sigma', hf_iter, sigma
     DO i=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(i) /= 'hole') CYCLE         
        WRITE(6,'(7HNumber:,6(I3,2X),2X,E20.10)') i, all_orbit%nn(i), all_orbit%ll(i), &
             all_orbit%jj(i), &
             all_orbit%nshell(i), all_orbit%itzp(i), all_orbit%e(i)
     ENDDO
     hf_iter = hf_iter+1
  ENDDO
  DO i=1, all_orbit%total_orbits
     IF ( keep_originalenergies == 'no') THEN 
        IF (all_orbit%model_space(i) == 'inside')  THEN
           all_orbit%evalence(i) = all_orbit%e(i)
        ELSE
           all_orbit%evalence(i) = 0.0_dp
        ENDIF
     ELSEIF ( keep_originalenergies == 'yes') THEN
        all_orbit%e(i) = all_orbit%e_original(i)
     ENDIF
     WRITE(8,'(7HNumber:,5(I4,2X),2X,E12.6,2X,E12.6,2X,A10,2X,A10)') i, all_orbit%nn(i), all_orbit%ll(i), &
          all_orbit%jj(i), all_orbit%itzp(i), all_orbit%e(i), all_orbit%evalence(i), &
          all_orbit%orbit_status(i), all_orbit%model_space(i)
  ENDDO
  DEALLOCATE ( hf_mtx)
  DEALLOCATE ( hf_eigen)
  DEALLOCATE ( hf_vect)
  !     new harmonic oscillator wave function
  wave_function=bhf_hol
  DEALLOCATE (int_factor, kin_energy )

END SUBROUTINE brueckner_hartree_fock
!
!    The Hartree-Fock diagram
!
SUBROUTINE diagram_HF(a,c,coeffs, iteration,onebody_diagram_HF, ncoeffs)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c, iteration, ncoeffs
  INTEGER :: j_min, j_max, jph, h, h1, h2, i
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  REAL(DP) :: val, ang_mom_factor, w
  REAL(DP), INTENT(INOUT) :: onebody_diagram_HF
  REAL(DP), DIMENSION(n_startenergy_g) :: ans

  onebody_diagram_HF=0.0_dp
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO h1=1, all_orbit%total_orbits
        IF  (all_orbit%jj(h) /= all_orbit%jj(h1)) CYCLE
        IF  (all_orbit%ll(h) /= all_orbit%ll(h1)) CYCLE
        IF  (all_orbit%itzp(h) /= all_orbit%itzp(h1)) CYCLE
        DO h2=1, all_orbit%total_orbits
           IF  (all_orbit%jj(h) /= all_orbit%jj(h2)) CYCLE
           IF  (all_orbit%ll(h) /= all_orbit%ll(h2)) CYCLE
           IF  (all_orbit%itzp(h) /= all_orbit%itzp(h2)) CYCLE
           SELECT CASE (type_of_renormv)
           CASE ('no-core')
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 onebody_diagram_HF=onebody_diagram_HF+ans(1)*ang_mom_factor*coeffs(h2,h)*coeffs(h1,h)
              ENDDO
           CASE ('v-nrg')
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 onebody_diagram_HF=onebody_diagram_HF+ans(1)*ang_mom_factor*coeffs(h2,h)*coeffs(h1,h)
              ENDDO
           CASE ('vlowk')
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 onebody_diagram_HF=onebody_diagram_HF+ans(1)*ang_mom_factor*coeffs(h1,h)*coeffs(h2,h)
              ENDDO
           CASE ('v-krg')
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 onebody_diagram_HF=onebody_diagram_HF+ans(1)*ang_mom_factor*coeffs(h1,h)*coeffs(h2,h)
              ENDDO
           CASE ('g-matrix')
              w = all_orbit%e(h2)+(all_orbit%e(c)+all_orbit%e(a))/2
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 CALL interpolate(w,e_start_g,ans,val)
                 onebody_diagram_HF=onebody_diagram_HF+val*ang_mom_factor*coeffs(h1,h)*coeffs(h2,h)
              ENDDO
           END SELECT
        ENDDO
     ENDDO
  ENDDO

END  SUBROUTINE diagram_HF

SUBROUTINE bhf_coefficients(ang_mom,itz,n_confs,ncoeffs,gmatrix_configs,bhf_coeff,coeffs)
  USE single_particle_orbits
  USE configurations
  USE constants
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER, INTENT(IN) :: ncoeffs, n_confs, ang_mom, itz
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  REAL(DP), DIMENSION(n_confs,n_confs), INTENT(INOUT)  :: bhf_coeff
  REAL(DP) :: dij, fnorm
  REAL(DP) :: e_dir, e_exc, e
  INTEGER :: ket , bra, a,b, la, ja, lb, jb, p,q, lq, lp, jp, jq

  fnorm = 0.d0
  DO bra = 1, n_confs
     a = gmatrix_configs%config_ab(bra+bra-1)
     b = gmatrix_configs%config_ab(bra+bra)
     la=all_orbit%ll(a); jb=all_orbit%jj(b)
     ja=all_orbit%jj(a); lb=all_orbit%ll(b)
     DO ket = 1, n_confs
        p = gmatrix_configs%config_ab(ket+ket-1)
        q = gmatrix_configs%config_ab(ket+ket)
        lp = all_orbit%ll(p); jp=all_orbit%jj(p)
        lq = all_orbit%ll(q); jq=all_orbit%jj(q)
        IF ( itz /= 0) THEN
           fnorm = 1.d0/ dij(a,b) /  dij(p,q)
        ELSEIF( itz == 0 ) THEN
           fnorm = 1.d0
        ENDIF
        e_dir=0.0_dp; e_exc=0.D0; e = 0.0_dp
        ! direct term
        IF ( (la == lp).AND.( lb == lq ).AND.( ja == jp ).AND.( jb == jq )) THEN
           e_dir = coeffs(p,a)*coeffs(q,b)
        ENDIF
        ! exchange term
        IF ( (la == lq).AND.( lb == lp ).AND.( ja == jq ).AND.( jb == jp ) .AND. itz /= 0 ) THEN
           e_exc = coeffs(p,b)*coeffs(q,a)
        ENDIF
        e= (e_dir-e_exc*((-1.0_dp)**((2*ang_mom-jp-jq)/2)))*fnorm
        bhf_coeff(bra,ket) = e
     ENDDO
  ENDDO

END SUBROUTINE bhf_coefficients

!
!     This function returns the matrix for V
!
SUBROUTINE pphhmtx(ja,jb,jc,jd,jt,gmtpn)
  USE stored_bare_interaction
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  LOGICAL TRIAG
  INTEGER, INTENT(IN) :: ja,jb,jc,jd,jt
  INTEGER :: iph
  REAL(DP), INTENT(INOUT), DIMENSION(n_startenergy_g) :: gmtpn
  REAL(DP) :: dij, delta, xxx

  gmtpn=0.
  IF(2*all_orbit%nn(ja)+all_orbit%ll(ja) +2*all_orbit%nn(jb)+all_orbit%ll(jb) > nlmax ) RETURN
  IF(2*all_orbit%nn(jc)+all_orbit%ll(jc) +2*all_orbit%nn(jd)+all_orbit%ll(jd) > nlmax ) RETURN
  IF((all_orbit%itzp(ja)+all_orbit%itzp(jb)) /= &
       (all_orbit%itzp(jc)+all_orbit%itzp(jd))) RETURN
  IF((-1)**(all_orbit%ll(ja)+all_orbit%ll(jb)) /=  &
       (-1)**(all_orbit%ll(jc)+all_orbit%ll(jd))) RETURN
  IF((ja == jb).AND.(MOD(jt,2)/=0)) RETURN       
  IF((jc == jd).AND.(MOD(jt,2)/=0)) RETURN       
  IF(triag(all_orbit%jj(ja),all_orbit%jj(jb),2*jt)) RETURN
  IF(triag(all_orbit%jj(jc),all_orbit%jj(jd),2*jt)) RETURN
  CALL mtx_elements(ja,jb,jc,jd,jt,gmtpn)
  !  Test purpose, leave as is
  !  IF( (all_orbit%itzp(ja)+all_orbit%itzp(jb) == 0 ) ) THEN
  !     gmtpn = delta(ja,jc)*delta(jb,jd)
  !  ELSE
  !     gmtpn = (delta(ja,jc)*delta(jb,jd)-delta(ja,jd)*delta(jb,jc)*iph((all_orbit%jj(ja)+all_orbit%jj(jb))/2-jt))/dij(ja,jb)/dij(jc,jd)
  !  ENDIF

  !  gmtpn = gmtpn/dij(ja,jb)/dij(jc,jd)

END SUBROUTINE pphhmtx
!
!     Calculates the crosscoupled matrix element type 1
!
SUBROUTINE cross_coupled_mtxel1(ja,jb,jc,jd,jtot,cross)
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE 
  REAL(DP),  DIMENSION (n_startenergy_g) :: ans
  INTEGER :: ja,jb,jc,jd,jtot,jbra_min, jket_min, jbra_max ,jket_max, &
       jt, j_min, j_max, iph
  REAL(DP) :: fnorm, ang_mom_factor
  REAL(DP), DIMENSION(n_startenergy_g), INTENT(OUT) :: cross

  cross=0.
  fnorm=SQRT(2.*jtot+1.)*iph((all_orbit%jj(ja)+all_orbit%jj(jd))/2+jtot)
  IF(iph(all_orbit%ll(ja)+all_orbit%ll(jb)) & 
       /=iph(all_orbit%ll(jc)+all_orbit%ll(jd))) RETURN
  IF((all_orbit%itzp(ja)+all_orbit%itzp(jb)) /= &
       (all_orbit%itzp(jc)+all_orbit%itzp(jd))) RETURN
  jbra_min=ABS(all_orbit%jj(ja)-all_orbit%jj(jb))
  jbra_max=all_orbit%jj(ja)+all_orbit%jj(jb)
  jket_min=ABS(all_orbit%jj(jc)-all_orbit%jj(jd))
  jket_max=all_orbit%jj(jc)+all_orbit%jj(jd)
  j_max=MIN(jbra_max,jket_max)/2
  j_min=MAX(jbra_min,jket_min)/2
  IF(j_min > j_max) RETURN
  DO jt=j_min,j_max
     ang_mom_factor=sjs(all_orbit%jj(jc),all_orbit%jj(ja),2*jtot, &
          all_orbit%jj(jb),all_orbit%jj(jd),2*jt)* &
          (2.*jt+1.)*fnorm*iph(jt)
     CALL pphhmtx(ja,jb,jc,jd,jt,ans)
     cross=cross+ans*ang_mom_factor
  ENDDO

END SUBROUTINE cross_coupled_mtxel1
!
!     Calculates the crosscoupled matrix element type 2
!
SUBROUTINE cross_coupled_mtxel2(ja,jb,jc,jd,jtot,cross)
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE 
  REAL(DP), DIMENSION (n_startenergy_g) :: ans
  INTEGER :: ja,jb,jc,jd,jtot, jbra_min, jket_min, jbra_max ,jket_max, &
       jt, j_min, j_max, iph
  REAL(DP) :: fnorm, ang_mom_factor
  REAL(DP), DIMENSION(n_startenergy_g), INTENT(OUT) :: cross

  cross=0.
  fnorm=SQRT(2.*jtot+1.)*iph((all_orbit%jj(ja)+all_orbit%jj(jd))/2+ &
       jtot+all_orbit%jj(jb))
  IF(iph(all_orbit%ll(ja)+all_orbit%ll(jb)) & 
       /=iph(all_orbit%ll(jc)+all_orbit%ll(jd))) RETURN
  IF((all_orbit%itzp(ja)+all_orbit%itzp(jb)) /= &
       (all_orbit%itzp(jc)+all_orbit%itzp(jd))) RETURN
  jbra_min=ABS(all_orbit%jj(ja)-all_orbit%jj(jb))
  jbra_max=all_orbit%jj(ja)+all_orbit%jj(jb)
  jket_min=ABS(all_orbit%jj(jc)-all_orbit%jj(jd))
  jket_max=all_orbit%jj(jc)+all_orbit%jj(jd)
  j_max=MIN(jbra_max,jket_max)/2
  j_min=MAX(jbra_min,jket_min)/2
  IF(j_min > j_max) RETURN
  DO jt=j_min,j_max
     ang_mom_factor=sjs(all_orbit%jj(jc),all_orbit%jj(jb),2*jtot, &
          all_orbit%jj(ja),all_orbit%jj(jd),2*jt)* &
          (2.*jt+1.)*fnorm
     CALL pphhmtx(ja,jb,jc,jd,jt,ans)
     cross=cross+ans*ang_mom_factor
  ENDDO

END SUBROUTINE cross_coupled_mtxel2
!
!Fast diagonalization of real(kind=8) symmetric matrices
!
SUBROUTINE diag_exact(cvec, ceig, n, h)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: h
  INTEGER :: k, np, p, i,j, i1, kvec, lwork, option, il ,iu, info
  REAL(DP) :: a(n*(n+1)), work(30*n), thresh, vl, vu
  LOGICAL :: flag(n)  
  COMPLEX*16 :: cvu, cvl
  COMPLEX*16, INTENT(INOUT) :: ceig(n), cvec(n,n)

  cvec = 0.0_dp; ceig = 0.0_dp
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
!Fast diagonalization of real(kind=8) symmetric matrices
!
SUBROUTINE diag_exactzgeev(cvec, ceig, n, h)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: h
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: cvec
  COMPLEX*16, DIMENSION(n), INTENT(INOUT) :: ceig
  REAL(DP), DIMENSION(2*n) :: rwork
  COMPLEX*16, DIMENSION(n,n) :: vl
  COMPLEX*16, DIMENSION(10000) :: work1
  INTEGER :: i, lda, ldb, ldvl, ldvr, lwork, info
  CHARACTER*1 :: jobvl, jobvr, balanc, sense

  balanc = 'n';  jobvl = 'n' ;  jobvr = 'v';  sense = 'n';  lda = n
  ldvl = 1;  ldvr = n;  lwork = 10000
  cvec = 0.0_dp; ceig = 0.0_dp
  CALL zgeev( jobvl, jobvr, n, h, lda, ceig, vl, ldvl, cvec, ldvr, &
       work1, lwork, rwork, info )

END SUBROUTINE diag_exactzgeev
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
  REAL(DP) :: cvec_max
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: temp

  ! diagonalize 2p-effective shell model hamiltonian, this takes a symmetric 
  ! matrix as input.
  ALLOCATE(cvec(n,n)); ALLOCATE(ceig(n))
  cvec = 0.0_dp; ceig = 0.0_dp
  CALL diag_exact(cvec, ceig, n, hamilton)
  cvec_pp = 0.0_dp; cvec_qp = 0.0_dp
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

SUBROUTINE eigenvalues_large_maxvector(cvec_pp,cvec_qp,ceig_p,hamilton,n,np,gmatrix_configs)
  USE single_particle_orbits
  USE configurations
  USE constants 
  IMPLICIT NONE
  TYPE (configuration_descriptor) , INTENT(IN)  :: gmatrix_configs
  INTEGER, INTENT(IN) :: n, np
  INTEGER :: k1, k2, i, j,nqq, npp, i1, i2
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: cvec, cvec_temp
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: ceig
  COMPLEX*16, INTENT(INOUT)  :: ceig_p(np)
  COMPLEX*16, INTENT(INOUT) :: cvec_pp(np,np)
  COMPLEX*16, INTENT(INOUT) :: cvec_qp(n-np, np) 
  COMPLEX*16, INTENT(IN) :: hamilton(n,n)
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: temp
  REAL(KIND=8) :: cvec_max
  INTEGER, ALLOCATABLE, DIMENSION(:) :: model
  ! diagonalize 2p-effective shell model hamiltonian, this takes a symmetric 
  ! matrix as input.
  ALLOCATE(cvec(n,n)); ALLOCATE(ceig(n)); ALLOCATE(model(np)); ALLOCATE(temp(n))
  ALLOCATE(cvec_temp(n,n))
  cvec = 0.; ceig = 0.
  ! Diagonalize in big space first
  CALL diag_exactzgeev(cvec, ceig, n, hamilton)
  cvec_pp = 0.; cvec_qp = 0.; model = 0; ceig_p = 0.
  ! Set up eigenvalues and eigenvectors for model space and excluded space
  ! using the largest overlap vector
  cvec_temp = cvec
  npp = 0
  DO k1=1, n
     i1=gmatrix_configs%config_ab(k1+k1-1)
     i2=gmatrix_configs%config_ab(k1+k1)
     ! loop over all model space coefficients of exact eigenvectors |k>
     IF(all_orbit%model_space(i1) == 'inside' .AND.   &
          all_orbit%model_space(i2) == 'inside' ) THEN
        temp = ABS( cvec_temp(k1,:) )
        cvec_max =  MAXVAL(temp(:))
        DO i = 1, n
           ! find the np eigenvectors |k>, cvec, with largest overlap with
           ! two particle model space vectors.
           IF (ABS(cvec_temp(k1,i)) == cvec_max) THEN
              ! store the i'th element of the np eigenvectors |k> and eigenvalues
              npp = npp + 1; model(npp) = i
              cvec_temp(:,i) = 0.
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  ! setup P and Q space overlap matrices: cvec_pp, cvec_qp
  DO i = 1, np
     k2 = i
     npp = 0; nqq = 0
     DO k1=1, n
        i1=gmatrix_configs%config_ab(k1+k1-1)
        i2=gmatrix_configs%config_ab(k1+k1)
        IF ( all_orbit%model_space(i1) == 'inside' .AND. &
             all_orbit%model_space(i2) == 'inside' ) THEN
           npp = npp + 1
           cvec_pp(npp,i) = cvec(k1,k2)
           ceig_p(i) = ceig(k2)
        ELSE
           nqq = nqq + 1
           cvec_qp(nqq,i) = cvec(k1,k2)
        ENDIF
     ENDDO
  ENDDO
  DEALLOCATE(ceig); DEALLOCATE(cvec); DEALLOCATE(model); DEALLOCATE(temp)
  DEALLOCATE(cvec_temp)

END SUBROUTINE eigenvalues_large_maxvector
!
! Lee Suzuki similarity transformation 
!
SUBROUTINE lee_suzuki( cvec_pp, cvec_qp, ceig, np, nq, heff, sp_energy )
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: np, nq
  COMPLEX*16, DIMENSION(np,np), INTENT(IN) :: cvec_pp
  COMPLEX*16, DIMENSION(nq,np), INTENT(IN)  :: cvec_qp
  COMPLEX*16, DIMENSION(np), INTENT(IN) :: ceig
  COMPLEX*16, DIMENSION(np,np), INTENT(INOUT) ::  heff
  COMPLEX*16, DIMENSION(np+nq,np+nq), INTENT(IN) :: sp_energy
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: omega2, omega_2inv, u, u_inv, cvec_pp_inv, &
       eigen_vec, vl, heff_rhs, temp1, omega, temp2
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: ceig_p, omega2_eig
  REAL(DP), DIMENSION(2*np) :: rwork
  COMPLEX*16, DIMENSION(10000) :: work1
  COMPLEX*16 :: d, sum1, norm
  INTEGER :: i_p,j_p, j_q, ii, jj, k, i_q , a_p, a_pp 
  INTEGER :: i, lda, ldb, ldvl, ldvr, info, lwork, ilo , ihi
  CHARACTER*1 :: jobvl, jobvr, balanc, sense
  REAL(DP), DIMENSION(np) :: scale, rconde, rcondv
  REAL(DP) :: abnrm, a1, a2, b1, b2
  INTEGER :: ipiv(np) ,j

  balanc = 'n';  jobvl = 'n' ;  jobvr = 'v';  sense = 'n';  lda = np
  ldvl = 1;  ldvr = np;  lwork = 10000
  ALLOCATE(omega2(np,np), omega_2inv(np,np),cvec_pp_inv(np,np)) 
  ALLOCATE(eigen_vec(np,np), vl(np,np), temp1(np,np)) 
  ALLOCATE(omega(nq,np)); ALLOCATE(ceig_p(np), omega2_eig(np))
  ALLOCATE(temp2(np,np))
  eigen_vec = (0.D0,0.D0) 
  temp1 = TRANSPOSE(cvec_pp) 
  CALL zgeev( jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  ! the P->Q space transformation matrix, omega 
  cvec_pp_inv = cvec_pp
  CALL cmplxmatinv(cvec_pp_inv, np, d)
  !  set up of omega
  DO i_p = 1, np
     DO i_q = 1, nq
        omega(i_q,i_p) = SUM( cvec_pp_inv(:,i_p)*cvec_qp(i_q,:) )
     ENDDO
  ENDDO
  !  set up  modified single-particle energy part O(P+QwP)
  temp2 = (0.D0,0.D0)
  DO i_p = 1, np
     DO j_p = 1, np 
        sum1 = (0.D0,0.D0) 
        DO i_q = 1, nq
           sum1 = sum1 +sp_energy(j_p,np+i_q)*omega(i_q,i_p) 
        ENDDO
        temp2(j_p,i_p) = sum1
     ENDDO
     temp2(i_p,i_p) =  temp2(i_p,i_p) + sp_energy(i_p,i_p)
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
  heff_rhs =  MATMUL( U, MATMUL( heff,u_inv ) ) 
  temp2 =  MATMUL( U, MATMUL( temp2,u_inv ) )
  DEALLOCATE(u, u_inv) 
  ! check if heff is symmetrized:
  heff = (0.D0,0.D0)
  heff = heff_rhs
  ! diagonalize 2p-effective shell model hamiltonian
  CALL zgeev( jobvl, jobvr, np, heff_rhs, lda, ceig_p, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  ! compare spectrum from exact and P-space diagonalization
  !  WRITE(6,*) 'Compare model space two-body spectrum with exact spectrum:' 
  !  DO i_p = 1, np
  !     a1 = 0.0; a2 = 0.0_dp
  !     a1 = REAL( ceig_p(i_p))
  !     a2 = IMAG( ceig_p(i_p))
  !     b1 = REAL( ceig(i_p) )
  !     b2 = IMAG( ceig(i_p) ) 
  !     WRITE(6,'(4E20.6)') A1, A2, B1, B2
  !  ENDDO
  heff = heff - temp2
  DEALLOCATE(omega2, omega_2inv, cvec_pp_inv) 
  DEALLOCATE(eigen_vec, vl, heff_rhs, temp1); DEALLOCATE(omega) 
  DEALLOCATE(temp2); DEALLOCATE(ceig_p, omega2_eig)

END SUBROUTINE lee_suzuki
