!             Program block renorm-modules.f90   
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95
!             LAST UPGRADE : April 2007
!             This program block contains the definiton of
!             all modules used by the various program blocks.
!
!    This module contains all constants and declarations 
!    of variables read in by the function read_data. These
!    variables are used by many functions.

MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  ! min and max isospin projection
  INTEGER, PUBLIC :: itzmin, itzmax
  ! min and max total two-body angular momentum in lab frame
  INTEGER, PUBLIC :: j_lab_min, j_lab_max
  ! min and max total two-body angular momentum in Rel-CoM frame
  INTEGER, PUBLIC :: jmin, jmax
  ! max value of 2n+l = nlmax in rel and CoM frame. lmax is max relative l
  INTEGER, PUBLIC :: nlmax, lmax, nmax, nlmax_model, lab_lmax, lab_nmax
  ! number of integration points in relative and CoM momenta
  INTEGER , PUBLIC :: n_rel, n_cm_mom, n_total
  ! number of starting energies used to compute the g-matrix
  INTEGER, PUBLIC :: n_startenergy_g
  ! number of starting energies used to compute the Q-box, typically 11
  ! can then compute up to the tenth derivative of the Q-box
  INTEGER, PUBLIC :: n_startenergy_veff
  ! number of excitations (using harmonic oscillator picture) for diagrams
  INTEGER, PUBLIC :: number_homega_exct
  ! A of closed shell core
  INTEGER, PUBLIC :: mass_nucleus
  ! N as number of electrons
  INTEGER, PUBLIC :: number_electrons
  ! define the physical system, 3dim electrons or nucleons, labelled by atomic_physics or nuclear_physics
  CHARACTER (LEN=100), PUBLIC :: physical_system
  ! define whether onebody, twobody or threebody interaction
  CHARACTER (LEN=100), PUBLIC :: n_body_int
  !    type of  renormalization and coulomb
  CHARACTER (LEN=100), PUBLIC :: type_of_renormv, coulomb_included
  !    determine the type of Pauli operator
  CHARACTER (LEN=100), PUBLIC :: pauli_operator
  ! the nn potential, cd-bonn, idaho-a, idaho-b etc, with csb and cib options
  CHARACTER (LEN=100), PUBLIC  :: type_of_pot, csb, cib
  ! test for square Pauli operator
  LOGICAL, PUBLIC :: square_calculation
  ! arrays of starting energies
  REAL(DP), PUBLIC, ALLOCATABLE :: e_start_g(:)
  ! which starting energy is used to compute the effective interaction
  REAL(DP), PUBLIC :: starting_energy
  ! oscillator energy and oscillator length
  REAL(DP), PUBLIC :: hbar_omega, oscl, cutoff, atomic_strength
  ! nucleon masses for different isospin channels
  REAL(DP), PUBLIC, DIMENSION(-1:1), PARAMETER:: p_mass = &
       (/938.27231_dp, 938.918725_dp, 939.56563_dp/)
! If you want to have an average mass, use the variables here, given by the 
! average neutron and proton masses.
!       (/938.918725_dp, 938.918725_dp, 938.918725_dp, 1232.0_dp, 1232.0_dp/)
  REAL(DP) , PARAMETER, PUBLIC :: p_massave =   938.918725_dp  !   2006 value of average (m_p+m_n)/2
  REAL(DP) , PARAMETER, PUBLIC :: e_mass =   510998.918_dp  !   2006 value for electron mass  in eV
  REAL(DP), PARAMETER, PUBLIC ::  bohr_r=0.05291772108_dp  ! Bohr radius in nm
  REAL(DP), PARAMETER, PUBLIC :: theta_rot = 0.0_dp! 0.125_dp
  REAL(DP), PARAMETER, PUBLIC :: hbarc = 197.326968_dp    !  2006 value
  REAL(DP), PARAMETER, PUBLIC :: hb2ip = hbarc*hbarc/p_massave ! only for nuclei
  REAL(DP), PUBLIC, PARAMETER :: pi = 3.141592741012573_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_2 = 1.570796370506287_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_4 = 0.7853981852531433_dp
END MODULE constants

!     Definition of single particle data

MODULE single_particle_orbits
  USE constants
  TYPE, PUBLIC :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, DIMENSION(:), POINTER :: nn, ll, jj, itzp, nshell,mvalue, baryon_spin
     CHARACTER (LEN=10), DIMENSION(:), POINTER :: orbit_status, model_space
     REAL(DP), DIMENSION(:), POINTER :: e, evalence
  END TYPE single_particle_descript
  TYPE, PUBLIC :: rel_cm_data
     INTEGER, DIMENSION(:), POINTER :: nrel, lrel
     INTEGER :: max_value
  END TYPE rel_cm_data
  TYPE (rel_cm_data), PUBLIC :: relcm_sp_data
  TYPE (single_particle_descript), PUBLIC :: all_orbit, neutron_data, &
       proton_data, mscheme_basis, electron_data
CONTAINS
  SUBROUTINE allocate_sp_array(this_array,n)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    IF (ASSOCIATED (this_array%nn) ) DEALLOCATE(this_array%nn)
    ALLOCATE(this_array%nn(n))
    IF (ASSOCIATED (this_array%ll) ) DEALLOCATE(this_array%ll)
    ALLOCATE(this_array%ll(n))
    IF (ASSOCIATED (this_array%jj) ) DEALLOCATE(this_array%jj)
    ALLOCATE(this_array%jj(n))
    IF (ASSOCIATED (this_array%itzp) ) DEALLOCATE(this_array%itzp)
    ALLOCATE(this_array%itzp(n))
    IF (ASSOCIATED (this_array%mvalue) ) DEALLOCATE(this_array%mvalue)
    ALLOCATE(this_array%mvalue(n))
    IF (ASSOCIATED (this_array%e) ) DEALLOCATE(this_array%e)
    ALLOCATE(this_array%e(n))
    IF (ASSOCIATED (this_array%evalence) ) DEALLOCATE(this_array%evalence)
    ALLOCATE(this_array%evalence(n))
    IF (ASSOCIATED (this_array%nshell) ) DEALLOCATE(this_array%nshell)
    ALLOCATE(this_array%nshell(n))
    IF (ASSOCIATED (this_array%baryon_spin) ) DEALLOCATE(this_array%baryon_spin)
    ALLOCATE(this_array%baryon_spin(n))
    IF (ASSOCIATED (this_array%orbit_status) ) DEALLOCATE(this_array%orbit_status)
    ALLOCATE(this_array%orbit_status(n))
    IF (ASSOCIATED (this_array%model_space) ) DEALLOCATE(this_array%model_space)
    ALLOCATE(this_array%model_space(n))
    !           blank all characters and zero all other values
    DO i= 1, n
       this_array%model_space(i)= ' '
       this_array%orbit_status(i)= ' '
       this_array%e(i)=0.0_dp
       this_array%baryon_spin(i)=0
       this_array%evalence(i)=0.0_dp
       this_array%nn(i)=0
       this_array%ll(i)=0
       this_array%jj(i)=0
       this_array%nshell(i)=0
       this_array%itzp(i)=0
       this_array%mvalue(i)=0
    ENDDO

  END SUBROUTINE allocate_sp_array

  SUBROUTINE deallocate_sp_array(this_array)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%nn) ; DEALLOCATE(this_array%ll)
    DEALLOCATE(this_array%jj) ;DEALLOCATE(this_array%itzp)
    DEALLOCATE(this_array%evalence) ;DEALLOCATE(this_array%mvalue)
    DEALLOCATE(this_array%e) ;DEALLOCATE(this_array%nshell)
    DEALLOCATE(this_array%orbit_status); DEALLOCATE(this_array%model_space)
    DEALLOCATE(this_array%baryon_spin)
  END SUBROUTINE deallocate_sp_array

  SUBROUTINE allocate_relcm_array(this_array,n)
    TYPE (rel_cm_data), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    IF (ASSOCIATED (this_array%nrel) ) DEALLOCATE(this_array%nrel)
    ALLOCATE(this_array%nrel(n))
    IF (ASSOCIATED (this_array%lrel) ) DEALLOCATE(this_array%lrel)
    ALLOCATE(this_array%lrel(n))
    !           blank all characters and zero all other values
    DO i= 1, n
       this_array%nrel(i)=0
       this_array%lrel(i)=0
    ENDDO

  END SUBROUTINE allocate_relcm_array

  SUBROUTINE deallocate_relcm_array(this_array)
    TYPE (rel_cm_data), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%nrel) ; DEALLOCATE(this_array%lrel)

  END SUBROUTINE deallocate_relcm_array

END MODULE single_particle_orbits

!     This module sets up the mpi start and end points
!     for where to search for Relative and cm configs
!     In addition we include the variables size_mpi
!     and my_rank_mpi given in the main function through
!     the calls to
!      CALL MPI_INIT ( err_mpi )
!      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, size_mpi, err_mpi )
!      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank_mpi, err_mpi )
!      which initialize MPI

MODULE jobs_mpi     
  INTEGER, PUBLIC :: length_config_all, length_config_mpi
  INTEGER, PUBLIC :: size_mpi, my_rank_mpi
  TYPE, PUBLIC :: jobs_mpi_descriptor    
     !        MPI modification
     INTEGER, DIMENSION(:), POINTER :: loop_start, loop_end
     INTEGER, DIMENSION(:,:), POINTER :: loop_start_all, loop_end_all
  END TYPE jobs_mpi_descriptor
  TYPE (jobs_mpi_descriptor), PUBLIC :: jobs

END MODULE jobs_mpi

!     the g-matrix in the rel and c.m. frame
MODULE relcm_gmatrix
  USE constants
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: v_com_rel(:,:,:)
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: gtf_tot_all(:,:,:), &
       gtf_tot_rank(:,:,:), &
       gtf_buffer(:,:,:)
END MODULE relcm_gmatrix

!     setup all partial wave combinations from jmin to jmax
!     defined on input, as default use jmin = 0 and jmax = 10
MODULE partial_waves
  USE constants
  INTEGER, PUBLIC :: no_channels, nchans
  PARAMETER ( nchans=200)
  INTEGER,PUBLIC :: orb_lrel_min(nchans),spin_rel(nchans), &
       jang_rel(nchans),iso(nchans), orb_lrel_max(nchans)
CONTAINS

  !
  !          Setup partial waves.
  !          The max value of the partial wave J is def
  !          in the input file and transferred as
  !          jmin and jmax
  !
  SUBROUTINE setup_channels
    IMPLICIT NONE
    INTEGER :: j_ang, lorb_min, lorb_max, it, i_spin,l_orb
    LOGICAL :: triag
    no_channels=0

    DO j_ang=jmin,jmax
       DO i_spin=0,1
          lorb_min=j_ang ; lorb_max=j_ang+i_spin
          DO l_orb=lorb_min,lorb_max
             DO it=itzmin, itzmax                 
                IF(.NOT.triag(j_ang,l_orb,i_spin)) THEN
                   IF ( ABS( it ) == 1) THEN
                      !     Pauli principle test for identical particles, pp or nn
                      IF(MOD(l_orb+1+i_spin,2) /= 0) THEN
                         IF ( l_orb == j_ang ) THEN
                            no_channels=no_channels+1
                            orb_lrel_max(no_channels)=l_orb
                            orb_lrel_min(no_channels)=l_orb
                         ELSEIF( l_orb /= j_ang ) THEN
                            IF ( j_ang > 0 ) THEN
                               no_channels=no_channels+1
                               orb_lrel_max(no_channels)=j_ang+i_spin
                               orb_lrel_min(no_channels)=j_ang-i_spin
                            ELSEIF ( j_ang == 0) THEN
                               no_channels=no_channels+1 
                               orb_lrel_max(no_channels)=j_ang+i_spin
                               orb_lrel_min(no_channels)=j_ang+i_spin
                            ENDIF
                         ENDIF
                         spin_rel(no_channels)=i_spin
                         jang_rel(no_channels)=j_ang
                         iso(no_channels)=it
                      ENDIF
                   ELSEIF (it == 0) THEN
                      !     For T_z=0 all possible waves are included, no restrictions
                      IF ( l_orb == j_ang ) THEN
                         no_channels=no_channels+1
                         orb_lrel_max(no_channels)=l_orb
                         orb_lrel_min(no_channels)=l_orb
                      ELSEIF( l_orb /= j_ang ) THEN
                         IF ( j_ang > 0 ) THEN
                            no_channels=no_channels+1
                            orb_lrel_max(no_channels)=j_ang+i_spin
                            orb_lrel_min(no_channels)=j_ang-i_spin
                         ELSEIF ( j_ang == 0) THEN
                            no_channels=no_channels+1
                            orb_lrel_max(no_channels)=j_ang+i_spin
                            orb_lrel_min(no_channels)=j_ang+i_spin
                         ENDIF
                      ENDIF
                      spin_rel(no_channels)=i_spin
                      jang_rel(no_channels)=j_ang
                      iso(no_channels)=it
                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE setup_channels

  !
  !          Setup partial waves for the Nucleon-Delta part
  !          The isospin value is only one, since this channel will be matched
  !          with  an NN state. No NNDeltaDelta interaction is studied. 
  !
  SUBROUTINE setup_channelsNDelta
    IMPLICIT NONE
    INTEGER :: j_ang, lorb_min, lorb_max, it, i_spin,l_orb, tndelta_min, tndelta_max
    LOGICAL :: triag
    no_channels=0

    tndelta_min = 1; tndelta_max = 1
    DO j_ang=jmin,jmax
       DO i_spin=1,2   !  We can couple channels with total spin = 1 and 2. The NN pot is no longer
          ! diagonal in S
          lorb_min=j_ang ; lorb_max=j_ang+i_spin
          DO l_orb=lorb_min,lorb_max
             DO it=tndelta_min, tndelta_max                 
                IF(.NOT.triag(j_ang,l_orb,i_spin)) THEN
                   IF ( l_orb == j_ang ) THEN
                      no_channels=no_channels+1
                      orb_lrel_max(no_channels)=l_orb
                      orb_lrel_min(no_channels)=l_orb
                   ELSEIF( l_orb /= j_ang ) THEN
                      IF ( j_ang > 0 ) THEN
                         no_channels=no_channels+1
                         orb_lrel_max(no_channels)=j_ang+i_spin
                         orb_lrel_min(no_channels)=j_ang-i_spin
                      ELSEIF ( j_ang == 0) THEN
                         no_channels=no_channels+1
                         orb_lrel_max(no_channels)=j_ang+i_spin
                         orb_lrel_min(no_channels)=j_ang+i_spin
                      ENDIF
                   ENDIF
                   spin_rel(no_channels)=i_spin
                   jang_rel(no_channels)=j_ang
                   iso(no_channels)=it
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE setup_channelsNDelta

END MODULE partial_waves

!           
!     This module contains the angular momentun functions
!     and transformation coefficients when going from
!     lab system  <--> cm system
!
MODULE ang_mom_functions
  USE constants
  integer, private, parameter:: maxjj=200
  REAL(DP), PRIVATE :: f_mb(maxjj),g_mb(maxjj),w_mb(maxjj),sfact(100),dfact(100)
  INTEGER, PRIVATE :: kh(4*maxjj)
  REAL(DP), PRIVATE :: q(maxjj,maxjj), cn(0:maxjj+1,0:maxjj+1)

CONTAINS
  !
  ! factorials for 3j,6j and 9j symbols
  ! for moshinsky trans brackets and for
  ! vector brackets
  !
  SUBROUTINE commons_to_angmom
    IMPLICIT NONE
    INTEGER :: l, k, i, j
    REAL(DP) :: a , sq_pi, fj, tfj, fk, s
    ! 3j, 6j and 9j symbols
    kh=1
    kh(2*maxjj) =0
    DO l=1,maxjj
       q(l,1)=1.0d0
       q(l,l)=1.0d0
       kh(l+l+2*maxjj)=0
    ENDDO
    DO l=2,maxjj-1
       DO k=2,l
          q(l+1,k)=q(l,k-1)+q(l,k)
       ENDDO
    ENDDO
    ! Moshinsky brackets
    f_mb(1)=0.
    g_mb(1)=LOG(0.5D0)
    w_mb(1)=0.
    DO i=2,maxjj
       a=i-1
       f_mb(i)=f_mb(i-1)+LOG(a)
       g_mb(i)=g_mb(i-1)+LOG(a+0.5D0)
       w_mb(i)=LOG(a+a+1.)
    ENDDO
    ! spherical harmonics
    cn=0.
    sq_pi=1./SQRT(2.*pi)
    DO j=0,maxjj+1
       cn(0,j)=SQRT(0.5*(2.*j+1.))
    ENDDO
    DO j=1,maxjj+1
       tfj=2.*j
       cn(j,j)=cn(j-1,j-1)*SQRT((tfj+1.)/tfj)
    ENDDO
    DO j=0,maxjj+1
       fj=FLOAT(j)
       DO k=1,j-1
          fk=FLOAT(k)
          cn(k,j)=cn(k-1,j)*SQRT((fj+fk)*(fj-fk+1.))*0.5/fk
       ENDDO
    ENDDO
    cn=cn*sq_pi
    ! Legendre functions of the second kind
    sfact(1) = 1.0_dp
    dfact(1) = 1.0_dp
    DO i = 1,30
       s = i
       sfact(i+1) = sfact(i)*s
       dfact(i+1) = (s+s+1.0d0)*dfact(i)
    ENDDO

  END SUBROUTINE commons_to_angmom
  !
  ! calculates 3j-symbols
  !
  REAL(DP) FUNCTION tjs(j_a,j_b,j_c,m_a,m_b,m_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,m_a,m_b,m_c
    INTEGER :: ja, jb, jc, mb, ma, mc, la, lb, lc, lt, ld, ja2, jb2, &
         jc2, i, k0, k1, k, ip
    REAL(DP) :: x, fn, p

    tjs=0.
    ja=(j_a+m_a)/2+1
    ma=(j_a-m_a)/2+1
    jb=(j_b+m_b)/2+1
    mb=(j_b-m_b)/2+1
    jc=(j_c+m_c)/2+1
    mc=(j_c-m_c)/2+1
    la=(j_b+j_c-j_a)/2+1
    lb=(j_c+j_a-j_b)/2+1
    lc=(j_a+j_b-j_c)/2+1
    lt=(j_a+j_b+j_c)/2+1
    ld=MIN(ja,jb,jc,ma,mb,mc,la,lb,lc)
    IF(((m_a+m_b+m_c) <= 0).AND.(ld > 0)) THEN
       ja2=j_a+m_a
       jb2=j_b+m_b
       jc2=j_c+m_c
       i=ja2+jb2+jc2-ja2/2*2-jb2/2*2-jc2/2*2
       IF(i == 0) then
          fn=q(ja+ma-1,lc)*q(jb+mb-1,lc)/(q(lt,jc+mc-1)*q(lt+1,2) &
               *q(ja+ma-1,ja)*q(jb+mb-1,jb)*q(jc+mc-1,jc))
          k0=MAX(0,lc-ja,lc-mb)+1
          k1=MIN(lc,ma,jb)
          x=0.
          DO k=k0,k1
             x=-x-q(lc,k)*q(lb,ma-k+1)*q(la,jb-k+1)
          ENDDO
          ip=k1+lb+jc
          p=1-2*(ip-ip/2*2)
          tjs=p*x*SQRT(fn)
       ENDIF
    ENDIF

  END FUNCTION tjs
  !
  ! calculates 6j-symbols
  !
  REAL(DP) FUNCTION sjs(j_a,j_b,j_c,l_a,l_b,l_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,l_a,l_b,l_c
    INTEGER :: ja,jb,jc,la,lb,lc,i,mt,ma,mb,mc,na,nb,nc,ka,&
         kb,kc,l,l0,l1, ihlp, idum
    REAL(DP) :: x, fs, fss

    ihlp=2*maxjj-1

    sjs=0.0d0
    ja=j_a + 1
    jb=j_b + 1
    jc=j_c + 1
    la=l_a + 1
    lb=l_b + 1
    lc=l_c + 1
    i=kh(ja+jb-jc+ihlp)+kh(jb+jc-ja+ihlp)+kh(jc+ja-jb+ihlp)+kh(ja+lb-lc+ihlp) &
         +kh(lb+lc-ja+ihlp)+kh(lc+ja-lb+ihlp)+kh(la+jb-lc+ihlp)+kh(jb+lc-la+ihlp)&
         +kh(lc+la-jb+ihlp)+kh(la+lb-jc+ihlp)+kh(lb+jc-la+ihlp)+kh(jc+la-lb+ihlp)
    IF(i <= 0) THEN
       mt=(j_a+j_b+j_c)/2 + 2
       ma=(j_a+l_b+l_c)/2+ 2
       mb=(l_a+j_b+l_c)/2+ 2
       mc=(l_a+l_b+j_c)/2+ 2
       na=mt-ja
       nb=mt-jb
       nc=mt-jc
       ka=ma-lc
       kb=mb-lc
       kc=mc-jc

       idum=max(mt,ja+1,nc,ma,ka,mb,la+1,kb,mc,kc)

       if(idum.gt.maxjj) then
          write(6,*) 'increase maxjj in MODULE ang_mom_functions from', maxjj, 'to', idum
          stop
       end if

       fss=q(mt,ja+1)*q(ja,nc)/(q(ma,ja+1)*q(ja,ka)*q(mb,la+1)* &
            q(la,kb)*q(mc,la+1)*q(la,kc))
       fs=SQRT(fss)/(l_a + 1.)
       l0=MAX(mt,ma,mb,mc)+1
       l1=MIN(ma+na,mb+nb,mc+nc)
       x=0.
       DO l=l0,l1
          x=-x+q(l-1,mt)*q(na,l-ma)*q(nb,l-mb)*q(nc,l-mc)
       ENDDO
       sjs=-(1+2*(l1/2*2-l1))*fs*x
    ENDIF

  END FUNCTION sjs
  !
  ! calculates ninej-symbols
  !
  REAL(DP) FUNCTION snj (ia,ib,ie,ic,id,if,ig,ih,it)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ia,ib,ie,ic,id,if,ig,ih,it
    INTEGER :: ja,jb,je,jc,jd,jf,jg,jh,jt,i,la,ld,ma,mc,na,nb,le,lf,&
         lg,me,mf,mg,ne,nf,ng,lx,mx,nx,jsi,jsf, js,is,lb, lc, &
         mb, ly, my,ny,l,l0,m0,n0,l1,m1,n1,m,n,ihx, ihlp
    REAL(DP) :: x, fn, fd, ps, fs, u, y, z, ud, p

    ihlp=2*maxjj-1

    snj=0.
    ja=ia+1
    jb=ib+1
    jc=ic+1
    jd=id+1
    je=ie+1
    jf=IF+1
    jg=ig+1
    jh=ih+1
    jt=it+1
    i=kh(ja+jb-je+ihlp)+kh(jb+je-ja+ihlp)+kh(je+ja-jb+ihlp)+kh(jc+jd-jf+ihlp) &
         +kh(jd+jf-jc+ihlp)+kh(jf+jc-jd+ihlp)+kh(jg+jh-jt+ihlp)+kh(jh+jt-jg+ihlp)&
         +kh(jt+jg-jh+ihlp)+kh(ja+jc-jg+ihlp)+kh(jc+jg-ja+ihlp)+kh(jg+ja-jc+ihlp)&
         +kh(jb+jd-jh+ihlp)+kh(jd+jh-jb+ihlp)+kh(jh+jb-jd+ihlp)+kh(je+jf-jt+ihlp)&
         +kh(jf+jt-je+ihlp)+kh(jt+je-jf+ihlp)
    IF(i <= 0) THEN
       la=(ie+IF+it)/2+2
       ld=(ig+ih+it)/2+2
       ma=(ia+ic+ig)/2+2
       mc=(IF+ic+id)/2+2
       na=(ib+id+ih)/2+2
       nb=(ib+ie+ia)/2+2
       le=(ie+IF-it)/2+1
       lf=(IF+it-ie)/2+1
       lg=(it+ie-IF)/2+1
       me=(ia+ic-ig)/2+1
       mf=(ic+ig-ia)/2+1
       mg=(ig+ia-ic)/2+1
       ne=(ib+id-ih)/2+1
       nf=(id+ih-ib)/2+1
       ng=(ih+ib-id)/2+1
       lx=(it+ig-ih)/2+1
       mx=(ic+id-IF)/2+1
       nx=(ib+ie-ia)/2+1
       fn=q(la,jt+1)*q(jt,lg)*q(ma,jc+1)*q(jc,mf)*q(na,jb+1)*q(jb,ne)
       fd=q(ld,jt+1)*q(jt,lx)*q(mc,jc+1)*q(jc,mx)*q(nb,jb+1)*q(jb,nx)
       jsi=MAX(ABS(je-jh),ABS(jg-jf),ABS(ja-jd))+1
       jsf=MIN(je+jh,jg+jf,ja+jd)-1
       ps=-1-2*(jsi/2*2-jsi)
       fs=ps*SQRT(fn/fd)/FLOAT((ig+1)*(ie+1))
       u=0.
       DO js=jsi,jsf,2
          is=js-1
          lb=(ie+ih+is)/2+2
          lc=(ig+IF+is)/2+2
          mb=(ia+id+is)/2+2
          ly=(ie+ih-is)/2+1
          my=(ig+IF-is)/2+1
          ny=(ia-id+is)/2+1
          ud=q(lb,je+1)*q(je,ly)*q(lc,jg+1)*q(jg,my)*q(mb,js+1)*q(js,ny)
          l0=MAX(la,lb,lc,ld)+1
          m0=MAX(ma,mb,mc,lc)+1
          n0=MAX(na,nb,mb,lb)+1
          l1=MIN(le+ld,lf+lb,lg+lc)
          m1=MIN(me+lc,mf+mb,mg+mc)
          n1=MIN(ne+lb,nf+nb,ng+mb)
          x=0.
          DO l=l0,l1
             x=-x-q(l-1,la)*q(le,l-ld)*q(lf,l-lb)*q(lg,l-lc)
          ENDDO
          y=0.
          DO m=m0,m1
             y=-y-q(m-1,ma)*q(me,m-lc)*q(mf,m-mb)*q(mg,m-mc)
          ENDDO
          z=0.
          DO n=n0,n1
             z=-z-q(n-1,na)*q(ne,n-lb)*q(nf,n-nb)*q(ng,n-mb)
          ENDDO
          ihx=l1+m1+n1
          p=1+2*(ihx/2*2-ihx)
          u=u+p*x*y*z/ud
       ENDDO
       snj=u*fs
    ENDIF

  END FUNCTION snj

  !
  ! This routine calculates the moshinsky vector bracket
  ! Note that D=mass1/mass2
  ! Ref m.sotona and m.gmitro comp.phys.comm 3(1972)53
  !
  REAL(DP) FUNCTION gmosh &
       (n,l,nc,lc,n1,l1,n2,l2,lr,d)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n,l,nc,lc,n1,l1,n2,l2,lr
    REAL(DP), INTENT(IN) :: d
    INTEGER :: ip,ixf,ix, iyi, iyf, j1f,j2,k1i,k1f,m1f,iy,m2f,k2, &
         m2,m2i,m1,j1,k2f,k2i,k1
    REAL(DP) :: dl, d1l, bb, ba, anorm, y, p, bc, cfac, bm , &
         sm, s, sxy, bxy

    gmosh=0.
    IF(n+n+nc+nc+l+lc-n1-n1-n2-n2-l1-l2 /= 0 ) RETURN
    IF(l+lc-lr < 0 ) RETURN
    IF(l1+l2-lr < 0 ) RETURN
    IF(ABS(l-lc)-lr > 0 ) RETURN
    IF(ABS(l1-l2)-lr > 0 ) RETURN
    DL=LOG(D)
    D1L=LOG(D+1.)
    bb=f_mb(n1+1)+f_mb(n2+1)+f_mb(n+1)-f_mb(nc+1)+ &
         g_mb(n1+l1+1)+g_mb(n2+l2+1) &
         -g_mb(n+l+1)-g_mb(nc+lc+1)
    ba=w_mb(l1+1)+w_mb(l2+1)+w_mb(lc+1)+w_mb(l+1)+ &
         f_mb(l1+l2-lr+1)+f_mb(l+lc+lr+2) &
         +f_mb(l+lc-lr+1)+f_mb(lc+lr-l+1)+ &
         f_mb(lr+l-lc+1)-f_mb(l1+l2+lr+2) &
         -f_mb(l1+lr-l2+1)-f_mb(l2+lr-l1+1)-DBLE(l)*d1l
    ip=lr+n+n1+n2
    p=1+2*(ip/2*2-ip)
    anorm=p*EXP(0.5D0*(bb+ba))
    y=0.
    j1f=l+1
    DO j1=1,j1f
       j2=l+2-j1
       k1i=ABS(l1-j1+1)+1
       k1f=l1+j1
       DO k1=k1i,k1f,2
          m1f=n1-(j1+k1-l1)/2+2
          IF(m1f-1 < 0 ) CYCLE
          k2i=MAX(ABS(l2-j2+1),ABS(lc-k1+1))+1
          k2f=MIN(l2+j2,lc+k1)
          IF(k2i-k2f > 0 ) CYCLE
          DO k2=k2i,k2f,2
             m2f=n2-(j2+k2-l2)/2+2
             IF(m2f-1 < 0 ) CYCLE
             ip=j2-1+(l1+k1+j1+l2+k2+j2)/2
             p=1+2*(ip/2*2-ip)
             bc=0.5D0*(DBLE(k1+j2-2)*dl-DBLE(k1+k2-2)*d1l) &
                  +f_mb(k1+l1-j1+1)+f_mb(k1+k2-lc-1)+ &
                  f_mb(k2+l2-j2+1)-f_mb(k1+l1+j1)-f_mb(k1+k2+lc)- &
                  f_mb(k2+l2+j2)+w_mb(k1)+w_mb(k2)+f_mb((k1+l1+j1)/2)+ &
                  f_mb((k1+k2+lc)/2)+f_mb((k2+l2+j2)/2)- &
                  f_mb((k1+l1-j1)/2+1)-f_mb((l1+j1-k1)/2+1)- &
                  f_mb((j1+k1-l1)/2)-f_mb((k1+k2-lc)/2)- &
                  f_mb((k2+lc-k1)/2+1)-f_mb((lc+k1-k2)/2+1) &
                  -f_mb((k2+l2-j2)/2+1)-f_mb((l2+j2-k2)/2+1)- &
                  f_mb((j2+k2-l2)/2)
             cfac=p*EXP(bc)
             sxy=0.
             ixf=MIN(k1+k1,k1+k2-lc)-1
             DO ix=1,ixf
                iyi=MAX(1,ix+j1+l2-k1-lr)
                iyf=MIN(l2+l2+1,l1+l2-lr+1,l2+lc+ix-k1-j2+2)
                IF(iyi-iyf > 0 ) CYCLE
                DO iy=iyi,iyf
                   ip=ix+iy
                   p=1+2*(ip/2*2-ip)
                   bxy=f_mb(k1+k1-ix)+f_mb(l2+l2-iy+2)+ &
                        f_mb(k2+lc-k1+ix)+f_mb(l1+lr-l2+iy) &
                        -f_mb(ix)-f_mb(iy)-f_mb(k1+k2-lc-ix)- &
                        f_mb(l1+l2-lr-iy+2)-f_mb(k1-l2+lr-j1+iy-ix+1)- &
                        f_mb(l2-k1+lc-j2+ix-iy+3)
                   sxy=sxy+p*EXP(bxy)
                ENDDO
             ENDDO
             s=cfac*sxy
             sm=0.
             DO m1=1,m1f
                m2i=MAX(1,nc-m1-(k1+k2-lc)/2+3)
                IF(m2i-m2f > 0 ) CYCLE
                DO m2=m2i,m2f
                   ip=m1+m2
                   p=1+2*(ip/2*2-ip)
                   bm=DBLE(m1-1)*DL-DBLE(m1+m2-2)*d1l+g_mb(1) &
                        +g_mb(m1+m2+(k1+k2+lc)/2-2)-g_mb(k1+m1-1)- &
                        g_mb(k2+m2-1)+f_mb(m1+m2+(k1+k2-lc)/2-2)- &
                        f_mb(m1)-f_mb(m2)-f_mb(n1-m1-(j1+k1-l1)/2+3)- &
                        f_mb(n2-m2-(j2+k2-l2)/2+3) &
                        -f_mb(m1+m2-nc+(k1+k2-lc)/2-2)
                   sm=sm+p*EXP(bm)
                ENDDO
             ENDDO
             y=y+s*sm
          ENDDO
       ENDDO
    ENDDO
    gmosh=anorm*y

  END FUNCTION gmosh
  ! This routine calculates the vector bracket
  ! allowing for a transformation from rel and com coordinates
  ! to the lab frame for plane waves. Present version assumes
  ! identical masses only. Dimensionless.
  !

  REAL(DP) FUNCTION vector_trcoefficients(ak,akk,ak1,ak2,l,ll,lam,l1,l2,tisoz)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ak,akk,ak1,ak2
    INTEGER , INTENT(IN) :: l,ll,lam,l1,l2,tisoz
    INTEGER :: ic1, ic2, m, mm, mu, mmu, mi, mp, ix, ixx, md
    REAL(DP) :: sak, dak, tmass, xm, xm1, xm2, xm3, x, xx, &
         xa, xb, xc, aiii, xd, psa, sl2, sll, sggn,&
         sumin, sign, sl, sl1, dr, dr1, dr2

    vector_trcoefficients=0.0_dp
    sak=ak1+ak2
    dak=DABS(ak1-ak2)
    IF (akk > sak) RETURN
    IF (akk < dak) RETURN
    ! tmass=p_mass(it1)+p_mass(it2)
    tmass=p_mass(tisoz)+p_mass(tisoz)
    ! xm=p_mass(it1)/tmass
    xm=p_mass(tisoz)/tmass
    ! xm1=(p_mass(it2)/tmass)**2
    xm1=(p_mass(tisoz)/tmass)**2
    xm2=xm**2
    ! xm3=p_mass(it1)*p_mass(it2)/(tmass**2)
    xm3=p_mass(tisoz)*p_mass(tisoz)/(tmass**2)
    x=(ak1*ak1-ak*ak-akk*akk*xm2)/(2.D0*xm*ak*akk)
    xx=x*x
    IF (xx > 1.0_dp) RETURN
    ic1=l1+l2+l+ll
    ic2=ic1-2*(ic1/2)
    IF (ic2 /= 0) RETURN
    aiii=0.0_dp
    xa=(ak1*ak1+ak*ak-akk*akk*xm2)/(2.d0*ak*ak1)
    xb=(ak1*ak1+akk*akk*xm2-ak*ak)/(2.D0*xm*ak1*akk)
    xc=(ak1*ak1*xm1+ak2*ak2*xm2-ak*ak)/(2.d0*ak1*ak2*xm3)
    md=0
    xd=1.0_dp
    sl1=spherical_harmonics(md,l1,xd)
    IF (sl1 == 0.0_dp) RETURN
    DO m=-l,l
       sl=spherical_harmonics(m,l,xa)
       IF (sl == 0.0_dp) CYCLE
       DO mm=-ll,ll
          mu=m+mm
          IF (IABS(mu) > l2) CYCLE
          mmu=-mu
          dr1=tjs(l+l,ll+ll,lam+lam,m+m,mm+mm,mmu+mmu)
          dr2=tjs(l1+l1,l2+l2,lam+lam,md+md,mu+mu,mmu+mmu)
          dr=dr1*dr2
          sign=1.
          mi=m-l-l1+l2+ll
          mp=mi-2*(mi/2)
          IF(mp /= 0) sign=-1.0_dp
          sll= spherical_harmonics(mm,ll,xb)
          IF (sll == 0.0_dp) CYCLE
          sl2= spherical_harmonics(mu,l2,xc)
          IF (sl2 == 0.0_dp) CYCLE
          psa=sl*sll*sl1*sl2
          sumin=0.
          sumin=dr*psa*sign
          aiii=aiii+dr*psa*sign
       ENDDO
    ENDDO
    vector_trcoefficients=16.0d0*pi*pi*aiii
    sggn=1.
    ix=(l1+l2-l-ll)/2
    ixx=ix-2*(ix/2)
    IF (ixx /= 0) sggn=-1.
    vector_trcoefficients=vector_trcoefficients*sggn

  END FUNCTION vector_trcoefficients

  ! Spherical harmonics from Num. Recipes

  REAL(DP) FUNCTION spherical_harmonics(m1,l,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m1, l
    REAL(DP), INTENT(IN) :: x
    REAL(DP), DIMENSION(0:51) :: y
    INTEGER :: iphase, m, j
    REAL(DP) :: fj, z, fac, div, sum, a, b, c
    spherical_harmonics=0.
    m=IABS(m1)
    IF(m.LT.0) m=-m1
    y(0)=1.
    IF(l.EQ.0) THEN
       sum=y(0)
    ELSE
       a=m-l
       b=l+m+1
       c=m+1
       z=0.5_dp-x*0.5_dp
       DO j=1,l-m+1
          fj=j-1
          y(j)=y(j-1)*(a+fj)*(b+fj)*z
          div=(c+fj)*(fj+1.)
          y(j)=y(j)/div
       ENDDO
       IF(m > 0) then
          fac=(1.-x*x)**m
          fac=SQRT(fac)
       ELSE
          fac=1.
       ENDIF
       sum=0.
       DO j=0,l-m
          sum=sum+y(j)
       ENDDO
       iphase=m
       IF(m1.LT.0) then
          iphase=0
       ENDIF
       sum=sum*fac*((-1)**iphase)
    ENDIF
    spherical_harmonics=cn(m,l)*sum

  END FUNCTION spherical_harmonics

END MODULE ang_mom_functions
!
!     Modules specific to the g-matrix calculation and effective operators
!     only
!     arrays containing mesh points and harmonic oscillator wave functions
!     In addition, routines for various wave functions are also included
!
MODULE wave_functions  
  USE constants
  INTEGER , PUBLIC :: n_k1, n_k, n_k2
  REAL(DP), ALLOCATABLE, PUBLIC :: ra(:),wra(:), krel(:), wkrel(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rgkk(:),wgkk(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rca(:),wrca(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: hol(:,:,:)
  ! ho onebody energies used for RG flow equations
  REAL(DP), PUBLIC, ALLOCATABLE :: ho_onebodyenergy(:)
  REAL(DP), PUBLIC :: k_cutoff, k_max
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: chol(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: bhf_hol(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rnlc(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rnlr(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: propagator(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: osc_propagator(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: coulomb_relcom(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC:: wave_cm(:)
CONTAINS

  !
  !     H.O. functions using Kummers function   
  !
  REAL(DP) FUNCTION rnl(n,l,z)
    IMPLICIT NONE
    INTEGER :: lll, nn
    INTEGER, INTENT(IN) :: l, n
    REAL(DP) :: y, dl, gamfaa, dfll, gamfab, dfnn
    REAL(DP), INTENT(IN) :: z

    rnl=0. ; y=0.5_dp*z*z
    IF(y > 60.0_dp) RETURN
    dl = l
    IF((ABS(z) < 1.0d-6) .AND. (l == 0)) rnl = 1.0_dp
    IF( ABS(z) > 1.0d-6) rnl = (z**l) * EXP(-y) * hypkum(n,dl+1.5_dp,z*z)
    gamfaa = 0.5_dp * SQRT(pi)
    IF(l /= 0) THEN
       DO lll = 1, l
          dfll = lll - 1
          gamfaa = gamfaa * (dfll + 1.5_dp)
       ENDDO
    ENDIF
    gamfab = gamfaa
    IF(n /= 0) THEN
       dfll = dl + 0.5_dp
       DO nn = 1, n
          dfnn = nn
          gamfab = gamfab * ((dfnn + dfll) / dfnn)
       ENDDO
    ENDIF
    rnl = rnl * (SQRT(2.0_dp * gamfab) / gamfaa)

  END FUNCTION rnl
  !
  !     Kummers function, Abramowitz & Stegun   
  !     exp. 13.1.2. a(there) equals (-n)       
  !  
  REAL(DP) FUNCTION hypkum(n,b,z)
    IMPLICIT NONE
    INTEGER :: nmax, nf
    INTEGER, INTENT(IN)  :: n
    REAL(DP) :: af, bf, zf, term, dfnf, xadd, sum
    REAL(DP), INTENT(IN) :: b, z

    IF(n < 0) WRITE (6,*)' error exit in hypkum ',  n,b,z
    hypkum = 1.0_dp
    IF(n == 0) RETURN
    nmax = n ; af = - n ; bf = b ; zf = z ; sum = 1.0 ; term = 1.0_dp
    DO nf = 1, nmax
       dfnf = nf
       xadd = dfnf - 1.0_dp
       term = term * ((af + xadd) / (bf + xadd)) * (zf / dfnf)
       IF(ABS(term) <  1.0d-12) EXIT
       sum = sum + term
    ENDDO
    hypkum = sum

  END FUNCTION hypkum
  !  This function sets up the recursive relation
  !  for the associated Legendre polynomials
  REAL(DP) FUNCTION legendre_polynomials(l, m, x)
    IMPLICIT NONE
    REAL(DP) ::  fact,pll,pmm,pmmp1,somx2
    REAL(DP), INTENT(IN)  :: x
    INTEGER ::  i,ll
    INTEGER, INTENT(IN) :: l, m

    !  check whether m, l and x are ok
    IF((M < 0).OR.(M > L).OR.(ABS(X) > 1.)) THEN
       WRITE(6,*) 'bad arguments', m, l, x; RETURN
    ENDIF
    !  calculate now pmm as starting point for iterations
    pmm=1.0
    IF (m > 0) THEN
       somx2=SQRT((1.0-x)*(1.0+x))
       fact=1.0_dp;
       DO i=1, m
          pmm = -fact*somx2*pmm
          fact = fact+2.0_dp
       ENDDO
    ENDIF
    !  if l == m we do not need to use recursion relation
    IF (l == m) THEN
       legendre_polynomials=pmm
       !  recursive relation for associated Legendre polynomials
    ELSE
       pmmp1=x*(2*m+1)*pmm
       !  analytical formula for the case l == m+1
       IF (l == (m+1)) THEN
          legendre_polynomials=pmmp1
       ELSE
          DO ll=m+2, l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          legendre_polynomials= pll
       ENDIF
    ENDIF

  END FUNCTION legendre_polynomials
  !
  !
  !      This routine calculates gauss-legendre mesh points and weights      
  !      input:                                                              
  !      x1   : lower limit of the integration interval                      
  !      x2   : upper limit ---------- "" -------------                      
  !      n    : the desired number of mesh points                            
  !      output :                                                            
  !      x     : gauss-legendre mesh points on the interval (x1,x2)          
  !      w     : the corresponding weights                                   
  !      From  : Numerical recipes
  !      F90 version : M. Hjorth-Jensen
  !
  SUBROUTINE gauss_legendre(x1,x2,x,w,n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j, m
    REAL(DP), INTENT(IN) :: x1, x2
    REAL(DP), INTENT(INOUT) :: x, w
    REAL(DP) :: eps
    DIMENSION :: x(n), w(n)
    PARAMETER (eps=3.D-14)
    REAL(DP) :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5_dp*(x2+x1)
    xl=0.5_dp*(x2-x1)
    DO i=1,m
       z1=0.
       z=COS(pi*(i-.25_dp)/(n+.5_dp))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.0_dp
          p2=0.0_dp
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.0_dp*j-1.)*z*p2-(j-1.0_dp)*p3)/j
          ENDDO
          pp=n*(z*p1-p2)/(z*z-1.)
          z1=z
          z=z-p1/pp
       ENDDO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO

  END SUBROUTINE gauss_legendre
  !
  !             Set up of c.m. mesh and weights
  !
  SUBROUTINE cm_mesh
    IMPLICIT NONE
    INTEGER :: i
    REAL(DP), DIMENSION(n_cm_mom) :: u,s

    CALL gauss_legendre(0.0_dp,cutoff,u,s,n_cm_mom)
    rca=u
    wrca=s

  END SUBROUTINE cm_mesh
  !
  !             Set up of relative mesh and weights
  !
  SUBROUTINE nocore_mesh
    IMPLICIT NONE
    INTEGER :: i
    REAL(KIND = 8) :: pih,u,s,xx,c
    PARAMETER (c=0.75)
    DIMENSION u(n_rel), s(n_rel)
    !   set up points for interpolation with HO basis, limited extent
    !   of the wave function
    CALL gauss_legendre(0.D0,cutoff,u,s,n_rel)
    ra=u
    wra=s

  END SUBROUTINE nocore_mesh
  !
  !             Set up of k mesh and weights for vfree obtained
  !             through a Lee-Suzuki similarity transformation
  !
  SUBROUTINE vlowk_mesh
    IMPLICIT NONE
    INTEGER :: i
    REAL(KIND=8) :: u1(n_k1), s1(n_k1),u2(n_k2), s2(n_k2)

    !   mesh points in units of fm-1
    CALL gauss_legendre(0.D0,k_cutoff,u1,s1,n_k1)
    CALL gauss_legendre(k_cutoff, k_max ,u2,s2,n_k2)
    DO i= 1, n_k
       IF ( i <= n_k1 ) THEN
          ra(i) = u1(i)
          wra(i)= s1(i)
       ELSEIF ( i > n_k1 ) THEN
          ra(i) = u2(i-n_k1)
          wra(i)= s2(i-n_k1)
       ENDIF
    ENDDO
    krel = ra*hbarc
    wkrel = wra*hbarc
    WRITE(6,*) 'points for ho oscillator basis:',  n_rel
    WRITE(6,'(5D16.8)') ra
    WRITE(6,'(5D16.8)') wra

  END SUBROUTINE vlowk_mesh
  !
  !             Set up of relative mesh and weights
  !
  SUBROUTINE rel_mesh
    IMPLICIT NONE
    INTEGER :: i
    REAL(DP) :: u,s,xx,c
    PARAMETER (c=0.75_dp)
    DIMENSION u(n_rel), s(n_rel)
    !   set up mesh points for G-mat calc with x \in [0,\infty]
    !   mapping of Gauss-Legendre mesh points u \in [-1,1]
    CALL gauss_legendre(-1.0_dp,1.0_dp,u,s,n_rel)
    DO i=1,n_rel
       xx=pi_4*(u(i)+1.0_dp); rgkk(i)=DTAN(xx)*c
       wgkk(i)=pi_4*c/DCOS(xx)**2*s(i)
    ENDDO
    WRITE(6,*) 'points for gkk matrix calculation:',  n_rel
    WRITE(6,'(5D16.8)') rgkk
    WRITE(6,'(5D16.8)') wgkk
    !   set up points for interpolation with HO basis, limited extent
    CALL gauss_legendre(0.0_dp,cutoff,u,s,n_rel)
    ra=u
    wra=s
    !    ra=rgkk
    !    wra=wgkk
    WRITE(6,*) 'points for ho oscillator basis:',  n_rel
    WRITE(6,'(5D16.8)') ra
    WRITE(6,'(5D16.8)') wra

  END SUBROUTINE rel_mesh

  !
  !             Set up of k mesh and weights for vfree obtained
  !             through the solution of RG flow equations
  !
  SUBROUTINE vkrgk_mesh               
    IMPLICIT NONE
    INTEGER :: i
    REAL(KIND=8) :: u(n_k), s(n_k)

    !   mesh points in units of fm-1
    CALL gauss_legendre(0.D0, k_max ,u,s,n_k)
    DO i= 1, n_k
       ra(i) = u(i)
       wra(i)= s(i)
    ENDDO
    krel = ra*hbarc
    wkrel = wra*hbarc
    WRITE(6,*) 'points for ho oscillator basis:',  n_rel
    WRITE(6,'(5D16.8)') ra
    WRITE(6,'(5D16.8)') wra

  END SUBROUTINE vkrgk_mesh


  SUBROUTINE laguerre_general( n, alpha, x, cx )
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL (dp ) ::  alpha
    REAL ( dp ) :: cx(0:n)
    INTEGER :: i
    REAL ( dp ), INTENT(IN) ::  x

    IF ( alpha <= -1.0D+00 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'LAGUERRE_GENERAL - Fatal error!'
       WRITE ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
       WRITE ( *, '(a)' ) '  but ALPHA must be greater than -1.'
       STOP
    END IF
    IF ( n < 0 ) THEN
       RETURN
    END IF
    cx(0) = 1.0D+00
    IF ( n == 0 ) THEN
       RETURN
    END IF
    cx(1) = 1.0D+00 + alpha - x
    DO i = 2, n
       cx(i) = ( ( REAL ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( REAL (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
            / REAL (     i,     kind = 8 )
    END DO

  END SUBROUTINE laguerre_general

  SUBROUTINE laguerre_complex( n, alpha, x, cx )
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL (dp ) ::  alpha
    COMPLEX ( DPC ) :: cx(0:n)
    INTEGER :: i
    COMPLEX ( DPC ), INTENT(IN) ::  x

    IF ( alpha <= -1.0D+00 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'LAGUERRE_GENERAL - Fatal error!'
       WRITE ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
       WRITE ( *, '(a)' ) '  but ALPHA must be greater than -1.'
       STOP
    END IF
    IF ( n < 0 ) THEN
       RETURN
    END IF
    cx(0) = 1.0D+00
    IF ( n == 0 ) THEN
       RETURN
    END IF
    cx(1) = 1.0D+00 + alpha - x
    DO i = 2, n
       cx(i) = ( ( REAL ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( REAL (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
            / REAL (     i,     kind = 8 )
    END DO

  END SUBROUTINE laguerre_complex

  DOUBLE PRECISION FUNCTION  fac(m)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: m
    INTEGER :: i

    fac = 0.0D0
    IF(m == 0) RETURN
    DO i=1,m
       fac=fac+LOG(FLOAT(i))
    ENDDO

  END FUNCTION  fac

  DOUBLE PRECISION FUNCTION  dfac(m)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: m
    INTEGER :: i

    IF (MOD(m,2).NE.1) STOP 'wrong argument to dfac'
    dfac = 0.0D0
    IF (m == 1)RETURN
    DO i=3,m,2
       dfac=dfac+LOG(FLOAT(i))
    ENDDO
  END FUNCTION  dfac

END MODULE wave_functions

!     module which defines configuration type, general structure
!     either for lab frame case or relcm system
MODULE configurations
  USE constants
  USE single_particle_orbits
  USE partial_waves
  USE ang_mom_functions
  TYPE configuration_descriptor        
     INTEGER :: number_confs
     INTEGER, DIMENSION(:), POINTER :: config_ab
  END TYPE configuration_descriptor
  TYPE, PUBLIC ::  configuration_relcm        
     INTEGER, DIMENSION(:), POINTER :: nconfs_relcm
     INTEGER, DIMENSION(:,:), POINTER :: relcm_ab, rel_ab 
     INTEGER, DIMENSION(:), POINTER :: nconfs_rel, nconfs_relmodel
  END TYPE configuration_relcm
  TYPE (configuration_relcm), PUBLIC :: relcm_conf, rel_conf

CONTAINS

  !                       
  !     setting up all configurations for given J, Tz and parity for
  !     the g-matrix model space

  SUBROUTINE number_gmatrix_confs(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         itza, itzb
    LOGICAL triag
    nconfs=0
    DO a=1, all_orbit%total_orbits
       na = all_orbit%nn(a)
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b=1, b_end
          nb = all_orbit%nn(b)
          lb=all_orbit%ll(b)
          itzb=all_orbit%itzp(b)
          jb=all_orbit%jj(b)
          IF ((all_orbit%model_space(a) == 'outside').AND.  &
               (all_orbit%model_space(b) == 'outside') ) CYCLE
          IF ((square_calculation).AND.(all_orbit%model_space(a) == 'inside').AND.  &
               (all_orbit%model_space(b) == 'outside') ) CYCLE
          IF ((square_calculation).AND.(all_orbit%model_space(a) == 'outside').AND.  &
               (all_orbit%model_space(b) == 'inside') ) CYCLE
          IF ( na+na+la+nb+nb+lb > nlmax) CYCLE
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
       ENDDO
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE number_gmatrix_confs
  !
  !
  !
  SUBROUTINE setup_gmatrix_configurations(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         k1, k2, itza, itzb
    LOGICAL triag
    nconfs=0
    DO a = 1, all_orbit%total_orbits
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       na = all_orbit%nn(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b = 1, b_end
          lb=all_orbit%ll(b)
          jb=all_orbit%jj(b)
          nb = all_orbit%nn(b)
          itzb=all_orbit%itzp(b)
          IF ( na+na+la+nb+nb+lb > nlmax) CYCLE
          IF ((all_orbit%model_space(a) == 'outside').AND.  &
               (all_orbit%model_space(b) == 'outside') ) CYCLE
          IF ((square_calculation).AND.(all_orbit%model_space(a) == 'inside').AND.  &
               (all_orbit%model_space(b) == 'outside') ) CYCLE
          IF ((square_calculation).AND.(all_orbit%model_space(a) == 'outside').AND.  &
               (all_orbit%model_space(b) == 'inside') ) CYCLE
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
          k2=nconfs*2
          k1=k2-1
          this%config_ab(k1)=b
          this%config_ab(k2)=a
       ENDDO
    ENDDO
    IF ( nconfs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF
    CALL sort_configs(this,nconfs)
    SELECT CASE (itz )
    CASE ( -1)
       WRITE(6,*) nconfs,' proton-proton configurations for J',ij
    CASE (1)
       IF ( physical_system == 'nuclear_physics') THEN
          WRITE(6,*) nconfs,' neutron-neutron configurations for J ',ij
       ELSEIF ( physical_system == 'atomic_physics') THEN
          WRITE(6,*) nconfs,' electron-electron configurations for J ',ij
       ENDIF
    CASE (0)
       WRITE(6,*) nconfs,' proton-neutron configurations for J ',ij
    END SELECT

  END SUBROUTINE setup_gmatrix_configurations
  !                       
  !     setting up all configurations for given J, Tz and parity for
  !     the nocore triangular  model space
  !
  SUBROUTINE number_nocore_confs(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         itza, itzb
    LOGICAL triag
    nconfs=0
    DO a=1, all_orbit%total_orbits
       na = all_orbit%nn(a)
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b=1, b_end
          nb = all_orbit%nn(b)
          lb=all_orbit%ll(b)
          itzb=all_orbit%itzp(b)
          jb=all_orbit%jj(b)
          IF ( na+na+la+nb+nb+lb > nlmax_model) CYCLE
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
       ENDDO
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE number_nocore_confs
  !
  !
  !
  SUBROUTINE setup_nocore_configurations(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         k1, k2, itza, itzb
    LOGICAL triag
    nconfs=0
    DO a = 1, all_orbit%total_orbits
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       na = all_orbit%nn(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b = 1, b_end
          lb=all_orbit%ll(b)
          jb=all_orbit%jj(b)
          nb = all_orbit%nn(b)
          itzb=all_orbit%itzp(b)
          IF ( na+na+la+nb+nb+lb > nlmax_model) CYCLE
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
          k2=nconfs*2
          k1=k2-1
          this%config_ab(k1)=b
          this%config_ab(k2)=a
       ENDDO
    ENDDO
    IF ( nconfs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF
    CALL sort_configs(this,nconfs)
    SELECT CASE (itz )
    CASE ( -1)
       WRITE(6,*) nconfs,' proton-proton configurations for J',ij
    CASE (1)
       IF ( physical_system == 'nuclear_physics') THEN
          WRITE(6,*) nconfs,' neutron-neutron configurations for J ',ij
       ELSEIF ( physical_system == 'atomic_physics') THEN
          WRITE(6,*) nconfs,' electron-electron configurations for J ',ij
       ENDIF
    CASE (0)
       WRITE(6,*) nconfs,' proton-neutron configurations for J ',ij
    END SELECT

  END SUBROUTINE setup_nocore_configurations


  !
  !     setting up all configurations for given J, Tz and parity for
  !     the single particle model space
  !
  SUBROUTINE number_sp_confs(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, a, nconfs

    nconfs=0
    DO a=1, all_orbit%total_orbits
       IF ( MOD(all_orbit%ll(a),2) /= ipar ) CYCLE
       IF ( all_orbit%jj(a) /= ij) CYCLE
       IF ( all_orbit%itzp(a) /= itz ) CYCLE
       nconfs=nconfs+1
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE number_sp_confs
  !
  !     Allocates space for the sp configurations and sets them up
  !
  SUBROUTINE setup_sp_configurations(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: a, ij, ipar, itz, nconfs
    nconfs=0
    DO a = 1, all_orbit%total_orbits
       IF ( MOD(all_orbit%ll(a),2) /= ipar ) CYCLE
       IF ( all_orbit%jj(a) /= ij) CYCLE
       IF ( all_orbit%itzp(a) /= itz ) CYCLE
       nconfs=nconfs+1
       this%config_ab(nconfs)=a
    ENDDO
    IF ( nconfs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF

  END SUBROUTINE setup_sp_configurations
  !
  !           Sort the two-particle configurations as
  !           ascending series
  !
  SUBROUTINE sort_configs(this,n_confs)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: i, j, tempa, tempb, n_confs

    DO i = 1, n_confs
       DO j =i, n_confs
          IF( (this%config_ab(2*i-1)*1000 +this%config_ab(2*i)) >  &
               (this%config_ab(2*j-1)*1000 +this%config_ab(2*j)) ) THEN 
             tempa=this%config_ab(2*i-1)
             tempb=this%config_ab(2*i)
             this%config_ab(2*i-1)=this%config_ab(2*j-1)
             this%config_ab(2*i)=this%config_ab(2*j)
             this%config_ab(2*j-1)=tempa
             this%config_ab(2*j)=tempb                             
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE sort_configs
  !
  !     setting up all configurations for given J, Tz, spin
  !     and parity for the rel coordinates only by the
  !     number channel. First find the number of configurations.
  !
  SUBROUTINE number_configurations_rel(this)
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, n_rel, a, nconfs, channel, l_min, l_max, nconfs_model

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0; nconfs_model = 0
       DO a=1, relcm_sp_data%max_value
          n_rel=relcm_sp_data%nrel(a)
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min) ) CYCLE
          IF (2*n_rel+l_rel > nlmax) CYCLE
          IF ( 2*n_rel+l_rel <= nlmax_model) THEN
             nconfs_model = nconfs_model +1
          ENDIF
          nconfs=nconfs+1
       ENDDO
       this%nconfs_rel(channel)=nconfs              
       this%nconfs_relmodel(channel)=nconfs_model              
    ENDDO

  END SUBROUTINE number_configurations_rel
  !
  !     Now we find the possible configurations
  !
  SUBROUTINE setup_configurations_rel(this)
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, n_rel, a, channel, l_min, l_max, nconfs

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0
       DO a=1, relcm_sp_data%max_value
          n_rel=relcm_sp_data%nrel(a)
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min ) ) CYCLE
          IF ( 2*n_rel+l_rel > nlmax) CYCLE
          nconfs=nconfs+1
          this%rel_ab(channel,nconfs)=a
       ENDDO
       IF ( nconfs /= this%nconfs_rel(channel) ) THEN
          WRITE(6,*) ' Error in configuration allocation ' ; STOP
       ENDIF
    ENDDO

  END SUBROUTINE setup_configurations_rel
  !
  !     setting up all configurations for given J, Tz, spin
  !     and parity for the rel and cm system specified by the
  !     number channel. First find the number of configurations.
  !
  SUBROUTINE number_configurations_relcm(this)
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, l_cm, n_rel, n_cm, a, b, &
         nconfs, channel, l_min, l_max

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0
       DO a=1, relcm_sp_data%max_value
          n_rel=relcm_sp_data%nrel(a)
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min) ) CYCLE
          DO b=1, relcm_sp_data%max_value
             n_cm=relcm_sp_data%nrel(b)
             l_cm=relcm_sp_data%lrel(b)
             IF ( 2*n_cm+l_cm+2*n_rel+l_rel > nlmax) CYCLE
             nconfs=nconfs+1
          ENDDO
       ENDDO
       this%nconfs_relcm(channel)=nconfs              
    ENDDO

  END SUBROUTINE number_configurations_relcm
  !
  !     Now we find the possible configurations
  !
  SUBROUTINE setup_configurations_relcm(this)
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, l_cm, n_rel, n_cm, a, b, &
         channel, l_min, l_max, k1, k2, nconfs

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0
       DO a=1, relcm_sp_data%max_value
          n_rel=relcm_sp_data%nrel(a)
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min ) ) CYCLE
          DO b=1, relcm_sp_data%max_value
             n_cm=relcm_sp_data%nrel(b)
             l_cm=relcm_sp_data%lrel(b)
             IF ( 2*n_cm+l_cm+2*n_rel+l_rel > nlmax) CYCLE
             nconfs=nconfs+1
             k2=nconfs*2
             k1=k2-1
             this%relcm_ab(channel,k1)=a
             this%relcm_ab(channel,k2)=b
          ENDDO
       ENDDO
       IF ( nconfs /= this%nconfs_relcm(channel) ) THEN
          WRITE(6,*) ' Error in configuration allocation ' ; STOP
       ENDIF
    ENDDO

  END SUBROUTINE setup_configurations_relcm

  !
  !     Finding all possible values for rel and cm variables n, N, l, L
  !
  SUBROUTINE find_spdata_relcm(a)
    IMPLICIT NONE
    INTEGER , INTENT(OUT) :: a
    INTEGER :: n, l

    a=0
    DO n=0, nmax
       DO l=0, lmax
          a=a+1
       ENDDO
    ENDDO

  END SUBROUTINE find_spdata_relcm
  !
  !     setting up all quantum numbers for relative and
  !     c.m coordinates in the block relcm_sp_data
  !
  SUBROUTINE make_configurations_relcm
    IMPLICIT NONE
    INTEGER :: a, n, l

    a=0
    DO n=0, nmax
       DO l=0, lmax
          a=a+1
          relcm_sp_data%nrel(a)=n
          relcm_sp_data%lrel(a)=l
       ENDDO
    ENDDO

  END SUBROUTINE make_configurations_relcm

END MODULE configurations
!
!     Function to calculate phase factors    
!
INTEGER FUNCTION iph(n)
  IMPLICIT NONE
  INTEGER :: n
  iph=(-1)**n
END FUNCTION iph
!
!     Function to check # osc. excitations   
! 
LOGICAL FUNCTION dencheck(i)
  USE constants
  IMPLICIT NONE
  INTEGER :: i 
  DENCHECK = ((ABS(i) == 0).OR.(ABS(i) > number_homega_exct))

END FUNCTION dencheck
!
!     Function to check triangular relations     
!
LOGICAL FUNCTION triag(i,j,k)
  IMPLICIT NONE
  INTEGER :: i, j, k
  triag = ((i-j-k)*(i-ABS(j-k)) > 0)

END FUNCTION triag
!
!      Function to calculate norm of g-mat    
!
REAL(KIND=8) FUNCTION dij(ja,jb)
  USE constants
  IMPLICIT NONE
  INTEGER :: ja, jb
  IF(ja == jb ) THEN
     dij=SQRT(2.0_dp)
  ELSE
     dij=1.0_dp
  ENDIF

END FUNCTION dij
!
!      Function to calculate norm of g-mat
!
REAL(DP) FUNCTION delta(ja,jb)
  USE constants
  IMPLICIT NONE
  INTEGER :: ja, jb
  IF(ja == jb ) THEN
     delta=1.d0
  ELSE
     delta=0.d0
  ENDIF

END FUNCTION delta
!
!     swaps values of 2 integers
!
SUBROUTINE SWAP(a, b)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: a, b
  INTEGER :: c

  c = a; a = b; b = c

END SUBROUTINE SWAP
!
!            Routines to do mtx inversion, from Numerical
!            Recepies, Teukolsky et al. Routines included
!            below are MATINV, LUDCMP and LUBKSB. See chap 2
!            of Numerical Recipes for further details
!            Recoded in FORTRAN 90 by M. Hjorth-Jensen
!
SUBROUTINE matinv(a,n)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  REAL(DP), DIMENSION(n,n), INTENT(INOUT)  :: a
  REAL(DP), ALLOCATABLE :: y(:,:)
  REAL(DP) :: d
  INTEGER, ALLOCATABLE :: indx(:)

  ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
  y=0.
  !     setup identity matrix
  DO i=1,n
     y(i,i)=1.
  ENDDO
  !     LU decompose the matrix just once
  CALL  lu_decompose(a,n,indx,d)

  !     Find inverse by columns
  DO j=1,n
     CALL lu_linear_equation(a,n,indx,y(:,j))
  ENDDO
  !     The original matrix a was destroyed, now we equate it with the inverse y 
  a=y

  DEALLOCATE ( y ); DEALLOCATE ( indx )

END SUBROUTINE matinv

!     Given an NxN matrix A(N,N), this routine replaces it by the LU 
!     decomposed one, where the matrix elements are stored in the same 
!     matrix A. The array indx is  an output vector which records the row
!     permutation effected by the partial pivoting. d is the determinant
!
SUBROUTINE lu_decompose(a,n,indx,d)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, i, j, k, imax
  REAL(DP) :: sum , tiny, aamax, dum, d
  REAL(DP), DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  REAL(DP), ALLOCATABLE :: vv(:)

  tiny=1.0e-20
  ALLOCATE ( vv(n) )
  D=1.
  DO i=1,n
     aamax=0.
     DO j=1,n
        IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
     ENDDO
     !     Zero is the largest element
     IF (aamax == 0.) STOP 'Singular matrix.'
     !     No nonzero largest element
     vv(i)=1./aamax
  ENDDO
  !     loop over columns
  DO j=1,n
     !     solves equation 2.3.12 except for i=j of Numerical Recipes
     IF (j > 1) THEN
        DO i=1,j-1
           sum=a(i,j)
           IF (i > 1)THEN
              DO k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
              ENDDO
              a(i,j)=sum
           ENDIF
        ENDDO
     ENDIF
     !    start searching for largest pivot element
     aamax=0.
     DO i=j,n
        sum=a(i,j)
        IF (j > 1)THEN
           DO k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           ENDDO
           a(i,j)=sum
        ENDIF
        dum=vv(i)*ABS(sum)
        IF (dum >= aamax) THEN
           imax=i
           aamax=dum
        ENDIF
     ENDDO
     !    interchange of rows
     IF (j /= imax)THEN
        DO k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        ENDDO
        !    change of parity for determinant
        d=-d
        vv(imax)=vv(j)
     ENDIF
     indx(j)=imax
     IF(j /= n) THEN
        IF(a(j,j) == 0.) a(j,j)=tiny
        dum=1./a(j,j)
        DO i=j+1,n
           a(i,j)=a(i,j)*dum
        ENDDO
     ENDIF
     !    set up determinant
     d=d*a(j,j)
  ENDDO
  IF(a(n,n) == 0.)  a(n,n)=tiny
  DEALLOCATE ( vv)

END SUBROUTINE lu_decompose

!     Solves set of linear equations Ax=b, A is input as an LU decompomsed
!     matrix and indx keeps track of the permutations of the rows. b is input
!     as the right-hand side vector b and returns the solution x. A, n and indx
!     are not modified by this routine. This function takes into that b can contain
!     many zeros and is therefore suitable for matrix inversion


SUBROUTINE lu_linear_equation(a,n,indx,b)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, ii, ll, i, j
  REAL(DP) :: sum 
  REAL(DP), DIMENSION(n,n) :: a
  REAL(DP), DIMENSION(n) :: b
  INTEGER, DIMENSION(n) :: indx

  ii=0
  !     First we solve equation 2.3.6 of numerical recipes 
  DO i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     IF (ii /= 0)THEN
        DO j=ii,i-1
           sum=sum-a(i,j)*b(j)
        ENDDO
     ELSEIF (sum /= 0.) THEN
        ii=i
     ENDIF
     b(i)=sum
  ENDDO
  !     then we solve equation 2.3.7
  DO i=n,1,-1
     sum=b(i)
     IF (i < n) THEN
        DO j=i+1,n
           sum=sum-a(i,j)*b(j)
        ENDDO
     ENDIF
     !     store a component of the solution x in the same place as b
     b(i)=sum/a(i,i)
  ENDDO

END SUBROUTINE lu_linear_equation
! 
! routine for calculating the principal square root ( and inverse ) 
! of a general matrix AA, i.e. solving the matrix equation: AA - A = 0. 
! The routine is built on the property : 
!
!            | 0  AA |      | 0    A | 
! sign(B) = (|       |) =  (|        |)
!            | 1   0 |      | A^-1 0 |
! 
! the sign of the matrix AA is calculated using Newtons iteration method,
! see Higham et. al. ref. Numerical Algorithms 15 (1997) 227-242 
! 
SUBROUTINE sqrtmat( aa, a, a_inv ,n ) 
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: aa
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: a, a_inv
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: x, x0, x1, x2
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: i_mat, temp
  INTEGER :: i,j,k
  COMPLEX*16 :: d

  ALLOCATE( x(n+n,n+n), x0(n+n,n+n), x1(n+n,n+n),x2(n+n,n+n))
  ALLOCATE( i_mat(n,n),temp(n,n))
  ! setup real identity matrix only
  i_mat = (0.D0,0.D0)
  DO i = 1, n
     i_mat(i,i) = (1.d0, 0.d0)
  ENDDO
  DO i = 1, 2*n
     DO j = 1, 2*n
        x0(j,i) = (0.d0, 0.d0)
        x1(j,i) = (0.d0, 0.d0)
        x2(j,i) = (0.d0, 0.d0)
        x(j,i) = (0.d0, 0.d0)
     ENDDO
  ENDDO
  DO i = 1, n
     DO j = 1, n
        temp(j,i) = (0.d0, 0.d0)
     ENDDO
  ENDDO
  DO i = n+1, 2*n
     DO j = 1, n
        x0(j,i) = aa(j,i-n)
     ENDDO
  ENDDO
  DO i = 1, n
     DO j = n+1, 2*n
        x0(j,i) = i_mat(j-n,i)
     ENDDO
  ENDDO
  k = 0
  DO WHILE( MAXVAL(ABS(temp-aa)) > 1.D-14 .AND.  k < 1000 )
     x1 = x0 
     x2 = x0 
     CALL cmplxmatinv(x2,n+n,d)
     x = 0.5d0 * ( x1 + x2 ) 
     x0 = x 
     k = k + 1
     DO i = 1, n
        DO j =1, n
           a(i,j) = x(i,j+n)
           a_inv(i,j) = x(i+n,j) 
        ENDDO
     ENDDO
     temp = MATMUL( a,a )
  ENDDO
  DEALLOCATE(i_mat,temp); DEALLOCATE(x,x0,x1,x2)

END SUBROUTINE sqrtmat
!
!    F90 program library, adapted from Numerical Recipes
!    All functions have been translated to F90 from F77
!    This is the complex*16 version of the 
!    routines to do matrix inversion, from Numerical
!    Recipes, Teukolsky et al. Routines included
!    below are MATINV, LUDCMP and LUBKSB. See chap 2
!    of Numerical Recipes for further details
!
SUBROUTINE cmplxmatinv(a,n,d)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT)  :: a
  COMPLEX*16, ALLOCATABLE :: y(:,:)
  COMPLEX*16 :: d
  INTEGER, ALLOCATABLE :: indx(:)

  ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
  y=0.
  !     setup identity matrix
  DO i=1,n
     y(i,i)=(1.d0, 0.d0) 
  ENDDO
  !     LU decompose the matrix just once
  CALL  cmplxlu_decompose(a,n,indx,d)

  !     Find inverse by columns
  DO j=1,n
     CALL cmplxlu_linear_equation(a,n,indx,y(:,j))
  ENDDO
  !     The original matrix a was destroyed, now we equate it with the inverse y 
  a=y

  DEALLOCATE ( y ); DEALLOCATE ( indx )

END SUBROUTINE cmplxmatinv
!
!     Given an NxN matrix A(N,N), this routine replaces it by the LU 
!     decomposed one, where the matrix elements are stored in the same 
!     matrix A. The array indx is  an output vector which records the row
!     permutation effected by the partial pivoting. d is the determinant
!
SUBROUTINE cmplxlu_decompose(a,n,indx,d)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, i, j, k, imax
  COMPLEX*16 :: sum, dum, tiny, d, aamax
  COMPLEX*16, DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  COMPLEX*16, ALLOCATABLE :: vv(:)

  tiny= ( 1.0D-20, 0.d0 ) 
  ALLOCATE ( vv(n) )
  D=1.
  DO i=1,n
     aamax=(0.D0,0.D0)
     DO j=1,n
        IF (ABS(a(i,j)) > ABS(aamax) ) aamax=ABS(a(i,j))
     ENDDO
     !     Zero is the largest element
     IF (aamax == (0.D0,0.D0)) STOP 'Singular matrix.'
     !     No nonzero largest element
     vv(i)=1./aamax
  ENDDO
  !     loop over columns
  DO j=1,n
     !     solves equation 2.3.12 except for i=j of Numerical Recipes
     IF (j > 1) THEN
        DO i=1,j-1
           sum=a(i,j)
           IF (i > 1)THEN
              DO k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
              ENDDO
              a(i,j)=sum
           ENDIF
        ENDDO
     ENDIF
     !    start searching for largest pivot element
     aamax=(0.D0,0.D0)
     DO i=j,n
        sum=a(i,j)
        IF (j > 1)THEN
           DO k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           ENDDO
           a(i,j)=sum
        ENDIF
        dum=vv(i)*ABS(sum)
        IF (ABS( dum ) >= ABS( aamax) ) THEN
           imax=i
           aamax=dum
        ENDIF
     ENDDO
     !    interchange of rows
     IF (j /= imax)THEN
        DO k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        ENDDO
        !    change of parity for determinant
        d=-d
        vv(imax)=vv(j)
     ENDIF
     indx(j)=imax
     IF(j /= n) THEN
        IF(a(j,j) == 0.) a(j,j)=tiny
        dum=1./a(j,j)
        DO i=j+1,n
           a(i,j)=a(i,j)*dum
        ENDDO
     ENDIF
     !    set up determinant
     d=d*a(j,j)
  ENDDO
  IF(a(n,n) == (0.d0,0.d0) )  a(n,n)=tiny
  DEALLOCATE ( vv)

END SUBROUTINE cmplxlu_decompose
!
!     Solves set of linear equations Ax=b, A is input as an LU decompomsed
!     matrix and indx keeps track of the permutations of the rows. b is input
!     as the right-hand side vector b and returns the solution x. A, n and indx
!     are not modified by this routine. This function takes into that b can contain
!     many zeros and is therefore suitable for matrix inversion
!
SUBROUTINE cmplxlu_linear_equation(a,n,indx,b)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, ii, ll, i, j
  COMPLEX*16 :: sum 
  COMPLEX*16, DIMENSION(n,n) :: a
  COMPLEX*16, DIMENSION(n) :: b
  INTEGER, DIMENSION(n) :: indx

  ii=0
  !     First we solve equation 2.3.6 of numerical recipes 
  DO i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     IF (ii /= 0)THEN
        DO j=ii,i-1
           sum=sum-a(i,j)*b(j)
        ENDDO
     ELSEIF (sum /= 0.) THEN
        ii=i
     ENDIF
     b(i)=sum
  ENDDO
  !     then we solve equation 2.3.7
  DO i=n,1,-1
     sum=b(i)
     IF (i < n) THEN
        DO j=i+1,n
           sum=sum-a(i,j)*b(j)
        ENDDO
     ENDIF
     !     store a component of the solution x in the same place as b
     b(i)=sum/a(i,i)
  ENDDO

END SUBROUTINE cmplxlu_linear_equation
!
!     Does the lagrangian interpolation                    
!     see abramowitz & stegun, p.878 & p.882            
!     x=pt. at which you evaluate                       
!     w= interpolation points                           
!     f= values of f at w's                             
!
SUBROUTINE interpolate(x,w,f,val)
  USE constants
  IMPLICIT NONE   
  REAL(DP) :: x, w, f, val, pip, p_i
  INTEGER :: i, k
  DIMENSION w(n_startenergy_g), f(n_startenergy_g), pip(n_startenergy_g)

  pip=1.
  p_i=1.0d0
  val=0.0d0
  DO i=1,n_startenergy_g
     IF (x == w(i)) val=f(i)
  ENDDO
  IF(val == 0.) THEN
     DO i=1,n_startenergy_g
        p_i=p_i*(x-w(i))
     ENDDO
     DO k=1,n_startenergy_g
        pip(k)=PRODUCT(w(k)-w, mask=(w(k) /= w) )
     ENDDO
     DO i=1,n_startenergy_g
        val = val + p_i*f(i)/((x-w(i))*pip(i))
     ENDDO
  ENDIF

END SUBROUTINE interpolate
!    
!                Two-dimensional lagrange interpol           
!     
SUBROUTINE lagrange_2dim(x,y,a,fxy)
  USE wave_functions
  USE constants
  IMPLICIT NONE
  INTEGER :: nx0, ny0, i0, j0, k, l, m
  REAL(DP), INTENT(INOUT) :: fxy
  REAL(DP), INTENT(IN) :: a(n_rel,n_rel)
  REAL(DP), INTENT(IN) :: x, y
  REAL(DP), DIMENSION(n_rel) :: q
  REAL(DP) :: fxk, p
  fxy=0.d0
  nx0=n_rel
  ny0=n_rel
  q = rgkk
  IF(X > 0.5*(q(nx0-2)+q(nx0-1))) THEN
     i0=nx0-3
  ELSE 
     i0=-1
400  i0=i0+1
     IF (x > 0.5*(q(i0+2)+q(i0+3))) GOTO 400
  ENDIF
  IF(y > 0.5*(q(ny0-2)+q(ny0-1))) THEN
     j0=ny0-3
  ELSE 
     j0=-1
410  j0=j0+1
     IF (y > 0.5*(q(j0+2)+q(j0+3))) GOTO 410
  ENDIF
  DO k=1,3
     fxk=0.
     DO l=1,3
        p=1.d0
        DO m=1,3
           IF (m /= l) p=p*(x-q(i0+m))/(q(i0+l)-q(i0+m))
        ENDDO
        fxk=fxk+p*a(i0+l,j0+k)
     ENDDO
     p=1.d0
     DO l=1,3
        IF (k /= l) p=p*(y-q(j0+l))/(q(j0+k)-q(j0+l))
     ENDDO
     fxy=fxy+fxk*p
  ENDDO

END SUBROUTINE lagrange_2dim
!
!     Function to find which number does the element from old arrays gtf_* have in gtf_tot_rank
!     To be called from slave processes
!
INTEGER FUNCTION find_address_rank ( num_channel, num_bra, num_ket )
  USE configurations
  USE jobs_mpi
  IMPLICIT NONE
  INTEGER num_channel, num_bra, num_ket, j

  IF (size_mpi == 1) THEN
     IF ( length_config_mpi /= length_config_all) THEN
        WRITE(6,*) 'Error counting length for size_mpi = 1'
        STOP
     ENDIF
  ENDIF
  !     The general formula for 1 process is:
  !     (i,j,k,l) = (n2, l)
  !      n2 = ( E_{j=1}^{i-1} K(j)^2 ) + ( j(i) - 1 ) * K(i) + k (i)
  !
  !     i = num_channel
  !     j = num_bra 
  !     k = num_ket, K(i) = MAX ( k(i) ) = conf_number
  find_address_rank = 0
  DO j = 1, num_channel - 1
     find_address_rank = find_address_rank + &
          ( jobs%loop_end(j) - jobs%loop_start(j) + 1) * relcm_conf%nconfs_relcm(j)
  ENDDO
  find_address_rank = find_address_rank + &
       ( num_bra - jobs%loop_start(num_channel)    ) * &
       relcm_conf%nconfs_relcm(j) + num_ket                     
  IF (find_address_rank > length_config_mpi) THEN
     WRITE (6,*) 'Error in find_address_rank ',find_address_rank ,' ', length_config_mpi
     WRITE (6,*) 'channel = ',num_channel, 'bra = ', num_bra  ,'ket = ', num_ket
  ENDIF
  IF (find_address_rank <= 0 ) THEN
     !          WRITE(6,*) 'Error in find_address_rank - less then or zero', find_address_rank
  ENDIF
  !	  WRITE (6,*) find_address_rank ,' ', length_config_mpi

END FUNCTION find_address_rank

!     Function to find which number does the element from old arrays gtf_* have in gtf_tot_all
!     
!     To be called from master process! rank does not matter when call with which = 'tot'
!           
INTEGER FUNCTION find_address_all ( num_channel, num_bra, num_ket, rank, which )
  USE configurations
  USE partial_waves
  USE jobs_mpi
  IMPLICIT NONE
  CHARACTER ( LEN = 3 ) which
  INTEGER num_channel, num_bra, num_ket, j, rank

  !     The general formula for 1 process is:
  !     (i,j,k,l) = (n2, l)
  !      n2 = ( E_{j=1}^{i-1} K(j)^2 ) + ( j(i) - 1 ) * K(i) + k (i)
  !
  !     i = num_channel
  !     j = num_bra 
  !     k = num_ket, K(i) = MAX ( k(i) ) = conf_number
  find_address_all = 0
  SELECT CASE (which)
  CASE ('rnk')
     DO j = 1, num_channel - 1
        find_address_all = find_address_all + &
             ( jobs%loop_end_all(rank + 1, j) - jobs%loop_start_all(rank + 1, j) + 1) * &
             relcm_conf%nconfs_relcm(j)
     ENDDO
     find_address_all = find_address_all + &
          ( num_bra - jobs%loop_start_all(rank + 1, num_channel)  ) * &
          relcm_conf%nconfs_relcm(j) + num_ket                     
     length_config_mpi = 0
     DO j = 1, no_channels
        length_config_mpi = length_config_mpi + &
             ( jobs%loop_end_all(rank + 1, j) - jobs%loop_start_all(rank + 1, j) + 1 ) * & 
             relcm_conf%nconfs_relcm(j)          
     ENDDO
     IF (find_address_all > length_config_mpi) THEN
        WRITE (6,*) 'Error in find_address_all ',find_address_all ,' ', length_config_mpi, &
             'which = ', which
        WRITE (6,*) 'channel = ',num_channel, 'bra = ', num_bra  ,'ket = ', num_ket
     ENDIF
  CASE ('tot')
     DO j = 1, num_channel - 1
        find_address_all = find_address_all + &
             relcm_conf%nconfs_relcm(j) * relcm_conf%nconfs_relcm(j)
     ENDDO
     find_address_all = find_address_all + &
          ( num_bra - 1) * relcm_conf%nconfs_relcm(j) + num_ket                     
     IF (find_address_all > length_config_all) THEN
        WRITE (6,*) 'Error in find_address_all ',find_address_all ,' ', length_config_all, &
             'which = ', which
        WRITE (6,*) 'channel = ',num_channel, 'bra = ', num_bra  ,'ket = ', num_ket
     ENDIF
  CASE DEFAULT
     WRITE (6,*) 'Could not recognize value of which in find_address_all'
  END SELECT

END FUNCTION find_address_all

!Module to read in name/value pairs from a file, with each line of the form line 'name = value'

MODULE inifile
  IMPLICIT NONE
  PUBLIC
  INTEGER, PARAMETER :: ini_max_name_len = 128
  INTEGER, PARAMETER :: ini_max_string_len = 1024
  LOGICAL :: ini_fail_on_not_found = .FALSE.
  LOGICAL :: ini_echo_read = .FALSE.
  TYPE tnamevalue
     !no known way to make character string pointers..
     CHARACTER(ini_max_name_len)  :: name
     CHARACTER(ini_max_string_len):: value
  END TYPE tnamevalue

  TYPE tnamevalue_pointer
     TYPE(tnamevalue), POINTER :: p
  END TYPE tnamevalue_pointer

  TYPE tnamevaluelist
     INTEGER count
     INTEGER delta
     INTEGER capacity
     TYPE(tnamevalue_pointer), DIMENSION(:), POINTER :: items
  END TYPE tnamevaluelist

  TYPE tinifile
     LOGICAL slashcomments
     TYPE (tnamevaluelist) :: l, readvalues
  END TYPE tinifile

  TYPE(tinifile) :: defini

CONTAINS

  SUBROUTINE tnamevaluelist_init(l)
    TYPE (tnamevaluelist) :: l

    l%count = 0
    l%capacity = 0
    l%delta = 128
    NULLIFY(l%items)

  END SUBROUTINE tnamevaluelist_init

  SUBROUTINE tnamevaluelist_clear(l)
    TYPE (tnamevaluelist) :: l
    INTEGER i, status

    DO i=l%count,1,-1
       DEALLOCATE (l%items(i)%p, stat = status)
    END DO
    DEALLOCATE (l%items, stat = status)
    CALL tnamevaluelist_init(l)

  END SUBROUTINE tnamevaluelist_clear

  SUBROUTINE tnamevaluelist_valueof(l, aname, avalue)
    TYPE (tnamevaluelist) :: l
    CHARACTER(len=*), INTENT(in) :: aname
    CHARACTER(len=*) :: avalue
    INTEGER i

    DO i=1, l%count
       IF (l%items(i)%p%name == aname) THEN
          avalue = l%items(i)%p%value 
          RETURN
       END IF
    END DO
    avalue = ''

  END SUBROUTINE tnamevaluelist_valueof

  SUBROUTINE tnamevaluelist_add(l, aname, avalue)
    TYPE (tnamevaluelist) :: l
    CHARACTER(len=*), INTENT(in) :: aname, avalue

    IF (l%count == l%capacity) CALL tnamevaluelist_setcapacity(l, l%capacity + l%delta)
    l%count = l%count + 1
    ALLOCATE(l%items(l%count)%p)
    l%items(l%count)%p%name = aname
    l%items(l%count)%p%value = avalue

  END SUBROUTINE tnamevaluelist_add

  SUBROUTINE tnamevaluelist_setcapacity(l, c)
    TYPE (tnamevaluelist) :: l
    INTEGER c
    TYPE(tnamevalue_pointer), DIMENSION(:), POINTER :: tmpitems

    IF (l%count > 0) THEN
       IF (c < l%count) STOP 'tnamevaluelist_setcapacity: smaller than count'
       ALLOCATE(tmpitems(l%count))
       tmpitems = l%items(1:l%count)
       DEALLOCATE(l%items)
       ALLOCATE(l%items(c))
       l%items(1:l%count) = tmpitems
       DEALLOCATE(tmpitems)
    ELSE
       ALLOCATE(l%items(c))
    END IF
    l%capacity = c

  END SUBROUTINE tnamevaluelist_setcapacity

  SUBROUTINE tnamevaluelist_delete(l, i)
    TYPE (tnamevaluelist) :: l
    INTEGER, INTENT(in) :: i

    DEALLOCATE(l%items(i)%p)
    IF (l%count > 1) l%items(i:l%count-1) = l%items(i+1:l%count)
    l%count = l%count -1

  END SUBROUTINE tnamevaluelist_delete

  SUBROUTINE ini_namevalue_add(ini,ainline)
    TYPE(tinifile) :: ini
    CHARACTER (len=*), INTENT(in) :: ainline
    INTEGER eqpos, slashpos, lastpos
    CHARACTER (len=LEN(ainline)) :: aname, s, inline

    inline=TRIM(ADJUSTL(ainline))
    eqpos = SCAN(inline,'=')
    IF (eqpos/=0 .AND. inline(1:1)/='#' .AND. inline(1:7) /= 'comment' ) THEN
       aname = TRIM(inline(1:eqpos-1))
       s = ADJUSTL(inline(eqpos+1:)) 
       IF (ini%slashcomments) THEN
          slashpos=SCAN(s,'/')
          IF (slashpos /= 0) THEN
             s  = s(1:slashpos-1)
          END IF
       END IF
       lastpos=LEN_TRIM(s)
       IF (lastpos>1) THEN
          IF (s(1:1)=='''' .AND. s(lastpos:lastpos)=='''') THEN
             s = s(2:lastpos-1)
          END IF
       END IF
       CALL tnamevaluelist_add(ini%l, aname, s)
    END IF

  END SUBROUTINE ini_namevalue_add


  SUBROUTINE ini_open(filename, unit_id,  error, slash_comments)
    CHARACTER (len=*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: unit_id
    LOGICAL, OPTIONAL, INTENT(out) :: error
    LOGICAL, OPTIONAL, INTENT(in) :: slash_comments
    LOGICAL aerror

    CALL tnamevaluelist_init(defini%l)
    CALL tnamevaluelist_init(defini%readvalues)
    IF (PRESENT(slash_comments)) THEN
       CALL ini_open_file(defini,filename,unit_id,aerror,slash_comments)
    ELSE
       CALL ini_open_file(defini,filename,unit_id,aerror)
    END IF
    IF (PRESENT(error)) THEN
       error = aerror
    ELSE
       IF (aerror) THEN
          WRITE (*,*) 'ini_open: error opening file ' // TRIM(filename)
          STOP
       END IF
    END IF

  END SUBROUTINE ini_open


  SUBROUTINE ini_open_file(ini, filename, unit_id,  error, slash_comments)
    TYPE(tinifile) :: ini
    CHARACTER (len=*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: unit_id
    LOGICAL, INTENT(out) :: error
    LOGICAL, OPTIONAL, INTENT(in) :: slash_comments
    CHARACTER (len=120) :: inline

    CALL tnamevaluelist_init(ini%l)
    CALL tnamevaluelist_init(ini%readvalues)

    IF (PRESENT(slash_comments)) THEN
       ini%slashcomments = slash_comments
    ELSE
       ini%slashcomments = .FALSE.
    END IF
    OPEN(unit=unit_id,file=filename,form='formatted',status='old', err=500)
    DO 
       READ (unit_id,'(a)',END=400) inline
       IF (inline == 'end') EXIT;
       IF (inline /= '') CALL ini_namevalue_add(ini,inline) 
    END DO
400 CLOSE(unit_id)
    error=.FALSE.
    RETURN
500 error=.TRUE.

  END SUBROUTINE ini_open_file

  SUBROUTINE ini_open_fromlines(ini, lines, numlines, slash_comments)
    TYPE(tinifile) :: ini
    INTEGER, INTENT(in) :: numlines
    CHARACTER (len=*), DIMENSION(numlines), INTENT(in) :: lines
    LOGICAL, INTENT(in) :: slash_comments
    INTEGER i

    CALL tnamevaluelist_init(ini%l)
    ini%slashcomments = slash_comments
    DO i=1,numlines
       CALL ini_namevalue_add(ini,lines(i))
    END DO

  END  SUBROUTINE ini_open_fromlines

  SUBROUTINE ini_close

    CALL ini_close_file(defini)

  END SUBROUTINE ini_close


  SUBROUTINE ini_close_file(ini)
    TYPE(tinifile) :: ini

    CALL tnamevaluelist_clear(ini%l)
    CALL tnamevaluelist_clear(ini%readvalues)

  END  SUBROUTINE ini_close_file


  FUNCTION ini_read_string(key, notfoundfail) RESULT(avalue)
    CHARACTER (len=*), INTENT(in) :: key
    LOGICAL, OPTIONAL, INTENT(in) :: notfoundfail
    CHARACTER(len=ini_max_string_len) :: avalue

    IF (PRESENT(notfoundfail)) THEN
       avalue = ini_read_string_file(defini, key, notfoundfail)
    ELSE
       avalue = ini_read_string_file(defini, key)
    END IF

  END FUNCTION ini_read_string


  FUNCTION ini_read_string_file(ini, key, notfoundfail) RESULT(avalue)
    TYPE(tinifile) :: ini
    CHARACTER (len=*), INTENT(in) :: key
    LOGICAL, OPTIONAL, INTENT(in) :: notfoundfail
    CHARACTER(len=ini_max_string_len) :: avalue

    CALL tnamevaluelist_valueof(ini%l, key, avalue)
    IF (avalue/='') THEN
       CALL  tnamevaluelist_add(ini%readvalues, key, avalue)
       IF (ini_echo_read) WRITE (*,*) TRIM(key)//' = ',TRIM(avalue)
       RETURN
    END IF
    IF (ini_fail_on_not_found) THEN
       WRITE(*,*) 'key not found : '//key
       STOP
    END IF
    IF (PRESENT(notfoundfail)) THEN
       IF (notfoundfail) THEN
          WRITE(*,*) 'key not found : '//key
          STOP
       END IF
    END IF

  END FUNCTION ini_read_string_file


  FUNCTION ini_read_int(key, default)
    INTEGER, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    INTEGER ini_read_int

    IF (PRESENT(default)) THEN
       ini_read_int = ini_read_int_file(defini, key, default)
    ELSE
       ini_read_int = ini_read_int_file(defini, key)
    END IF

  END FUNCTION ini_read_int


  FUNCTION ini_read_int_file(ini, key, default)
    TYPE(tinifile) :: ini
    INTEGER ini_read_int_file
    INTEGER, OPTIONAL, INTENT(in) :: default
    CHARACTER  (len=*), INTENT(in) :: key
    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini, key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_int_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       IF (VERIFY(TRIM(s),'-+0123456789') /= 0) GOTO 10
       READ (s,*, err = 10) ini_read_int_file
    END IF
    RETURN
10  WRITE (*,*) 'error reading integer for key: '//key
    STOP

  END FUNCTION ini_read_int_file

  FUNCTION ini_read_double(key, default)
    DOUBLE PRECISION, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    DOUBLE PRECISION ini_read_double

    IF (PRESENT(default)) THEN
       ini_read_double = ini_read_double_file(defini, key, default)
    ELSE
       ini_read_double = ini_read_double_file(defini, key)
    END IF

  END FUNCTION ini_read_double


  FUNCTION ini_read_double_file(ini,key, default)
    TYPE(tinifile) :: ini
    DOUBLE PRECISION ini_read_double_file 
    DOUBLE PRECISION, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini,key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_double_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       READ (s,*, err=10) ini_read_double_file
    END IF
    RETURN
10  WRITE (*,*) 'error reading double for key: '//key
    STOP

  END FUNCTION ini_read_double_file


  FUNCTION ini_read_real(key, default)
    REAL, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    REAL ini_read_real

    IF (PRESENT(default)) THEN
       ini_read_real = ini_read_real_file(defini, key, default)
    ELSE
       ini_read_real = ini_read_real_file(defini, key)
    END IF

  END FUNCTION ini_read_real

  FUNCTION ini_read_real_file(ini,key, default)
    TYPE(tinifile) :: ini
    REAL ini_read_real_file 
    REAL, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini,key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_real_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       READ (s,*, err=10) ini_read_real_file
    END IF
    RETURN
10  WRITE (*,*) 'error reading double for key: '//key
    STOP

  END FUNCTION ini_read_real_file


  FUNCTION ini_read_logical(key, default)
    LOGICAL, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    LOGICAL ini_read_logical

    IF (PRESENT(default)) THEN
       ini_read_logical = ini_read_logical_file(defini, key, default)
    ELSE
       ini_read_logical = ini_read_logical_file(defini, key)
    END IF

  END FUNCTION ini_read_logical

  FUNCTION ini_read_logical_file(ini, key, default)
    TYPE(tinifile) :: ini
    LOGICAL ini_read_logical_file
    LOGICAL, OPTIONAL, INTENT(in) :: default
    CHARACTER  (len=*), INTENT(in) :: key

    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini,key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_logical_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       IF (VERIFY(TRIM(s),'10tf') /= 0) GOTO 10  
       READ (s,*, err = 10) ini_read_logical_file
    END IF

    RETURN

10  WRITE (*,*) 'error reading logical for key: '//key
    STOP
  END FUNCTION ini_read_logical_file


  SUBROUTINE ini_savereadvalues(afile,unit_id)
    CHARACTER(len=*)  :: afile
    INTEGER, INTENT(in) :: unit_id

    CALL ini_savereadvalues_file(defini, afile, unit_id)

  END SUBROUTINE ini_savereadvalues


  SUBROUTINE ini_savereadvalues_file(ini, afile, unit_id)
    TYPE(tinifile) :: ini
    CHARACTER(len=*), INTENT(in) :: afile
    INTEGER, INTENT(in) :: unit_id
    INTEGER i

    OPEN(unit=unit_id,file=afile,form='formatted',status='replace', err=500)

    DO i=1, ini%readvalues%count
       WRITE (unit_id,'(a)') TRIM(ini%readvalues%items(i)%p%name) // ' = ' &
            //TRIM(ini%readvalues%items(i)%p%value)

    END DO

    CLOSE(unit_id)
    RETURN

500 WRITE(*,*) 'ini_savereadvalues_file: error creating '//TRIM(afile)

  END SUBROUTINE ini_savereadvalues_file

END MODULE inifile





