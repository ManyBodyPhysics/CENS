!             Program block heff-diagrams.f
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F90  
!             LAST UPGRADE : October 2007
!
!             This program block sets up all types of diagrams for
!             one-body, two-body, three-body, effective operators and closed core diagrams.
!             Three-body diagrams are only to second order in the interaction
!
!     Begin one-body diagrams
!
!
!
!
!     HF diagram
!
SUBROUTINE diagram_1(a,c,onebody_diagram_1)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: j_min, j_max, jph, h, i
  REAL(DP), DIMENSION(n_startenergy_veff) :: w 
  REAL(DP) :: val, ang_mom_factor
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_1
  REAL(DP), DIMENSION(n_startenergy_g) :: ans

  onebody_diagram_1=0.
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     w=all_orbit%evalence(c)+all_orbit%e(h)+wcn
     DO jph=j_min,j_max
        ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
        CALL pphhmtx(a,h,c,h,jph,ans) ; IF ( ans(1) == 0.0_dp) CYCLE
        DO i=1, n_startenergy_veff
           call interpolate(w(i),e_start_g,ans,val)
           onebody_diagram_1(i)=onebody_diagram_1(i)+val*ang_mom_factor
        ENDDO
     ENDDO
  ENDDO
END  SUBROUTINE diagram_1
!
!     2p1h diagram
!
SUBROUTINE diagram_2(a,c,onebody_diagram_2)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h, p1, p2, i, nshell1, nshell2, iosc
  REAL(DP) :: val1, val2, factr
  REAL(DP), DIMENSION(n_startenergy_veff) :: w , de 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL DENCHECK

  onebody_diagram_2=0.
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     nshell1=all_orbit%nshell(h)+all_orbit%nshell(c)
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/(all_orbit%jj(a)+1.)
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) == 'hole') CYCLE         
           DO p2=1, all_orbit%total_orbits
              IF (all_orbit%orbit_status(p2) == 'hole') CYCLE         
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              iosc=nshell2-nshell1
              IF(dencheck(iosc))CYCLE
              de=all_orbit%evalence(c)+all_orbit%e(h)- &
                   all_orbit%e(p1)-all_orbit%e(p2)+wcn
              CALL pphhmtx(a,h,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              CALL pphhmtx(p1,p2,c,h,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
              w=all_orbit%evalence(c)+all_orbit%e(h)+wcn
              DO  i=1, n_startenergy_veff
                 CALL interpolate(w(i),e_start_g,ans1,val1)
                 CALL interpolate(w(i),e_start_g,ans2,val2)
                 onebody_diagram_2(i)=onebody_diagram_2(i)+ &
                      factr*0.5*val1*val2/de(i)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_2
!
!     2h1p diagram
!
SUBROUTINE diagram_3(a,c,onebody_diagram_3)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, p, i, h1, h2, nshell1, nshell2, iosc
  REAL(DP) :: val1, val2, factr
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2 , de 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_3
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL DENCHECK

  onebody_diagram_3=0.
  DO p=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p) == 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/(all_orbit%jj(a)+1.)
        DO h1=1,all_orbit%total_orbits
           IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE         
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE         
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
              iosc=nshell2-nshell1
              IF(dencheck(iosc)) CYCLE
              de=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p)-all_orbit%evalence(a)+wcn
              CALL pphhmtx(h1,h2,c,p,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              CALL pphhmtx(a,p,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                   all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              DO i=1, n_startenergy_veff
                 CALL  interpolate(w1(i),e_start_g,ans1,val1)
                 CALL  interpolate(w2(i),e_start_g,ans2,val2)
                 onebody_diagram_3(i)=onebody_diagram_3(i)- & 
                      factr*0.5d0*val1*val2/de(i)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_3



SUBROUTINE diagram_4(a,c,onebody_diagram_4)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h, p1, p2, p3, p4, &
       nshell1, nshell2, nshell3, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w, de 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_4
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_4=0.
  DO h=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
     nshell1=all_orbit%nshell(h)+all_orbit%nshell(c)
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              iosc1=nshell2-nshell1
              IF (dencheck(iosc1)) CYCLE
              CALL pphhmtx(a,h,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              den1=all_orbit%evalence(c)+all_orbit%e(h)-all_orbit%e(p1)-all_orbit%e(p2)
              DO p3=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                 DO p4=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p4) == 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(p3)+all_orbit%nshell(p4)
                    iosc2=nshell1-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%evalence(c)+all_orbit%e(h)-  &
                         all_orbit%e(p3)-all_orbit%e(p4)
                    CALL pphhmtx(p1,p2,p3,p4,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p3,p4,c,h,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w=all_orbit%evalence(c)+all_orbit%e(h)+wcn
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w(i),e_start_g,ans1,val1)
                       CALL interpolate(w(i),e_start_g,ans2,val2)
                       CALL interpolate(w(i),e_start_g,ans3,val3)
                       onebody_diagram_4(i)=onebody_diagram_4(i)+ &
                            0.25*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_4



SUBROUTINE diagram_5(a,c,onebody_diagram_5)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p1, p2, p3, h2, &
       nshell1, nshell2, nshell3, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_5
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_5=0.
  DO p3=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p3))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p3))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              nshell2=all_orbit%nshell(a)+all_orbit%nshell(p3)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(a)-all_orbit%e(p3)
              CALL pphhmtx(h1,h2,c,p3,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
              DO p1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
                 DO p2=1, all_orbit%total_orbits                    
                    IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(p1)+& 
                         all_orbit%nshell(p2)
                    iosc2=nshell1-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                         all_orbit%e(p1)-all_orbit%e(p2)
                    CALL pphhmtx(a,p3,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p1,p2,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w3=all_orbit%e(h1)+all_orbit%e(h2)+ &
                         all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_5(i)=onebody_diagram_5(i)- &
                            0.25*val1*factr*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_5



SUBROUTINE diagram_6(a,c,onebody_diagram_6)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p1, p2, p3, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de      
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_6
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_6=0.
  DO p2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE 
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p2))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
     nshell2=all_orbit%nshell(a)+all_orbit%nshell(p2)
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(a)-all_orbit%e(p2)
              DO p1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
                 CALL pphhmtx(a,p2,h1,h2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                 DO p3=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(c)+nshell1
                    nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p3)+ &
                         all_orbit%nshell(a)
                    iosc2=nshell4-nshell3
                    den2=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p3)+all_orbit%evalence(c)-all_orbit%evalence(a)
                    IF(dencheck(iosc2)) CYCLE
                    CALL pphhmtx(h1,h2,p1,p3,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p1,p3,c,p2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w1=all_orbit%e(h2)+all_orbit%e(h1)+wcn
                    w2=all_orbit%e(h2)+all_orbit%e(h1)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w3=w2
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_6(i)=onebody_diagram_6(i)- &
                            0.25d0*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_6


SUBROUTINE diagram_7(a,c,onebody_diagram_7)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p, h2, h3, h4, &
       nshell1, nshell2, nshell3, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_7
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_7=0.
  DO p=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE   
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE   
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p)- &
                   all_orbit%evalence(a)
              CALL pphhmtx(a,p,h1,h2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              DO h3=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE 
                 DO h4=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h3)+all_orbit%nshell(h4)
                    iosc2=nshell2-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%e(h3)+all_orbit%e(h4)-all_orbit%e(p)- &
                         all_orbit%evalence(a)
                    CALL pphhmtx(h1,h2,h3,h4,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(h3,h4,c,p,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w2=all_orbit%e(h3)+all_orbit%e(h4)-all_orbit%evalence(a)- &
                         all_orbit%e(p)+w1
                    w3=all_orbit%e(h3)+all_orbit%e(h4)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_7(i)=onebody_diagram_7(i)- &
                            0.25*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_7



SUBROUTINE diagram_8(a,c,onebody_diagram_8)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p1, p2, h3, h2, &
       nshell1, nshell2, nshell3, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_8
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_8=0.
  DO h3=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h3))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h3))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h3)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h3)+all_orbit%evalence(c)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              DO h1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                    iosc2=nshell3-nshell2
                    den2=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p2)
                    IF(dencheck(iosc2)) CYCLE
                    CALL pphhmtx(a,h3,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p1,p2,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(h1,h2,c,h3,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w1=all_orbit%evalence(c)+all_orbit%e(h3)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w3=w2+all_orbit%e(h3)+all_orbit%evalence(c)- &
                         all_orbit%e(p1)-all_orbit%e(p2)
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_8(i)=onebody_diagram_8(i)+ &
                            0.25*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_8



SUBROUTINE diagram_9(a,c,onebody_diagram_9)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p1, p2, h3, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP) :: val1, val2, val3, factr,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_9
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_9=0.
  DO h2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h2))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              iosc1=nshell2-nshell1
              den1=all_orbit%e(h2)+all_orbit%evalence(c)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              IF(dencheck(iosc1)) CYCLE
              DO h1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                 DO h3=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h1)+all_orbit%nshell(c)+ &
                         all_orbit%nshell(h3)
                    nshell4=all_orbit%nshell(a)+nshell2
                    iosc2=nshell4-nshell3
                    den2=all_orbit%e(h1)+all_orbit%e(h3)-all_orbit%e(p1)- &
                         all_orbit%e(p2)+all_orbit%evalence(c)-all_orbit%evalence(a)
                    IF(dencheck(iosc2)) CYCLE
                    CALL pphhmtx(h1,h3,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(a,h2,h1,h3,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p1,p2,c,h2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w1=all_orbit%e(h1)+all_orbit%e(h3)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=all_orbit%e(h2)+all_orbit%evalence(a)-all_orbit%e(p1)- &
                         all_orbit%e(p2)+w1
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_9(i)=onebody_diagram_9(i)+ &
                            0.25*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_9



SUBROUTINE diagram_10(a,c,onebody_diagram_10)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, h3, h2, par1, par2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, iph
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_10
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_10=0.
  DO h1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        IF(all_orbit%jj(h1).ne.all_orbit%jj(h2)) CYCLE
        IF(iph(all_orbit%ll(c)+all_orbit%ll(h2)) /= &
             iph(all_orbit%ll(a)+all_orbit%ll(h1))) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(h1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(h1))/2
        DO h3=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
           par1=iph(all_orbit%ll(h1)+all_orbit%ll(h3))
           par2=iph(all_orbit%ll(h2)+all_orbit%ll(h3))
           IF(par1 /= par2) CYCLE 
           j2min=ABS((all_orbit%jj(h3)-all_orbit%jj(h1))/2)
           j2max=(all_orbit%jj(h3)+all_orbit%jj(h1))/2
           DO jtot1=j1min,j1max
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(h1))+1.)
                 DO p1=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
                    DO p2=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                       nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h3)
                       nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                       iosc1=nshell2-nshell1
                       nshell3=all_orbit%nshell(h2)+all_orbit%nshell(h3)+ &
                            all_orbit%nshell(c)
                       nshell4=nshell2+all_orbit%nshell(a)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc1).or.dencheck(iosc2)) CYCLE
                       den1=all_orbit%e(h1)+all_orbit%e(h3)- &
                            all_orbit%e(p1)-all_orbit%e(p2)
                       den2=all_orbit%e(h2)+all_orbit%e(h3)- &
                            all_orbit%e(p1)-all_orbit%e(p2)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)
                       CALL pphhmtx(a,h1,c,h2,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(h2,h3,p1,p2,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(p1,p2,h1,h3,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%evalence(c)+all_orbit%e(h3)- &
                            all_orbit%e(p1)-all_orbit%e(p2)+wcn
                       w2=all_orbit%e(h2)+all_orbit%e(h3)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
                       w3=all_orbit%e(h1)+all_orbit%e(h3)+wcn
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_10(i)=onebody_diagram_10(i)- &
                               0.5*factr1*factr2*  &
                               val1*val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_10

SUBROUTINE diagram_11(a,c,onebody_diagram_11)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, p3, h2, par1, par2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, iph
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_11
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_11=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE    
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
        IF(all_orbit%jj(p1) /= all_orbit%jj(p2)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        DO p3=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
           par1=iph(all_orbit%ll(p1)+all_orbit%ll(p3))
           par2=iph(all_orbit%ll(p2)+all_orbit%ll(p3))
           IF(par1 /= par2) CYCLE 
           j2min=ABS((all_orbit%jj(p3)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(p3)+all_orbit%jj(p1))/2
           DO jtot1=j1min,j1max
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
                 DO h1=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                    DO h2=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                       nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                       nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p3)
                       iosc1=nshell2-nshell1
                       nshell3=all_orbit%nshell(c)+nshell1
                       nshell4=all_orbit%nshell(p3)+all_orbit%nshell(p2)+ &
                            all_orbit%nshell(a)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc1).or.dencheck(iosc2)) CYCLE
                       den1=all_orbit%e(h1)+all_orbit%e(h2)-  &
                            all_orbit%e(p1)-all_orbit%e(p3)
                       den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                            all_orbit%e(p2)-all_orbit%e(p3)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)
                       CALL pphhmtx(a,p2,c,p1,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(h1,h2,p2,p3,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(p1,p3,h1,h2,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%evalence(c)-all_orbit%e(p3)+wcn
                       w2=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
                       w3=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_11(i)=onebody_diagram_11(i) &
                               +0.5*factr1*factr2*val1* &
                               val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_11


SUBROUTINE diagram_12(a,c,onebody_diagram_12)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, p3, h2, par1, par2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, iph
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_12
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_12=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE    
        IF(all_orbit%jj(p1) /= all_orbit%jj(h1)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        nshell1=all_orbit%nshell(h1)
        nshell2=all_orbit%nshell(p1)
        iosc1=nshell2-nshell1
        IF(dencheck(iosc1)) CYCLE
        den1=all_orbit%e(h1)-all_orbit%e(p1)

        DO h2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE    
           par1=iph(all_orbit%ll(p1)+all_orbit%ll(h2))
           par2=iph(all_orbit%ll(h1)+all_orbit%ll(h2))
           IF(par1 /= par2) CYCLE
           j2min=ABS((all_orbit%jj(h2)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(h2)+all_orbit%jj(p1))/2
           DO jtot1=j1min,j1max
              CALL pphhmtx(a,h1,c,p1,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
                 DO p2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                    DO p3=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                       nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                       nshell4=all_orbit%nshell(p3)+all_orbit%nshell(p2)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc2)) CYCLE
                       den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                            all_orbit%e(p2)-all_orbit%e(p3)
                       CALL pphhmtx(p1,h2,p2,p3,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(p2,p3,h1,h2,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                       w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                       w3=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_12(i)=onebody_diagram_12(i)+0.5 &
                               *factr1*factr2*val1*val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_12


SUBROUTINE diagram_13(a,c,onebody_diagram_13)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, p3, h2, par1, par2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, iph
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_13
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_13=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE    
        IF(all_orbit%jj(p1) /= all_orbit%jj(h1)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
        nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
        iosc1=nshell2-nshell1
        den1=all_orbit%e(h1)+all_orbit%evalence(c)-all_orbit%evalence(a)-all_orbit%e(p1)
        IF(dencheck(iosc1)) CYCLE
        DO h2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE    
           par1=iph(all_orbit%ll(p1)+all_orbit%ll(h2))
           par2=iph(all_orbit%ll(h1)+all_orbit%ll(h2))
           IF(par1 /= par2) CYCLE
           j2min=ABS((all_orbit%jj(h2)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(h2)+all_orbit%jj(p1))/2
           DO jtot1=j1min,j1max
              CALL pphhmtx(a,p1,c,h1,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
                 DO p2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE 
                    DO p3=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE 
                       nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                       nshell4=all_orbit%nshell(p2)+all_orbit%nshell(p3)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc2)) CYCLE
                       den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                            all_orbit%e(p2)-all_orbit%e(p3)
                       CALL pphhmtx(h1,h2,p2,p3,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(p2,p3,p1,h2,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                       w2=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
                       w3=w2
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_13(i)=onebody_diagram_13(i)+ &
                               0.5*factr1*factr2*val1*val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_13




SUBROUTINE diagram_14(a,c,onebody_diagram_14)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, h3, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_14
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_14=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE    
        IF(all_orbit%jj(p1) /= all_orbit%jj(h1)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        nshell1=all_orbit%nshell(h1)
        nshell2=all_orbit%nshell(p1)
        iosc1=nshell2-nshell1
        den1=all_orbit%e(h1)-all_orbit%e(p1)
        IF(dencheck(iosc1)) CYCLE
        DO p2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE 
           j2min=ABS((all_orbit%jj(p2)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
           DO jtot1=j1min,j1max
              CALL pphhmtx(a,h1,c,p1,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE 
                    DO h3=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE    
                       nshell3=all_orbit%nshell(h2)+all_orbit%nshell(h3)
                       nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc2)) CYCLE
                       den2=all_orbit%e(h3)+all_orbit%e(h2)- &
                            all_orbit%e(p1)-all_orbit%e(p2)
                       CALL pphhmtx(p1,p2,h2,h3,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(h2,h3,h1,p2,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                       w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
                       w3=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%e(h3)-all_orbit%e(p1)+wcn
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_14(i)=onebody_diagram_14(i)- &
                               0.5*factr1*factr2*val1*val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_14
!
!
!
SUBROUTINE diagram_15(a,c,onebody_diagram_15)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, h3, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, ie1, ie2
  REAL(DP) :: val1, val2, val3, factr2, factr1
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1,w2,w3,den1, den2, &
       de,xmtx1,xmtx2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_15
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_15=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE 
        IF(all_orbit%jj(p1) /= all_orbit%jj(h1)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
        nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
        iosc1=nshell2-nshell1
        den1=wcn+all_orbit%e(h1)+all_orbit%evalence(c)-all_orbit%evalence(a)- &
             all_orbit%e(p1)
        IF(dencheck(iosc1)) CYCLE
        xmtx1=0.
        DO jtot1=j1min,j1max
           factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
           CALL pphhmtx(a,p1,c,h1,jtot1,ans1); IF(ans1(1).eq.0.) CYCLE 
           w1=all_orbit%evalence(c)+all_orbit%e(h1)+wcn
           DO ie1=1, n_startenergy_veff
              CALL interpolate(w1(ie1),e_start_g,ans1,val1)
              xmtx1(ie1)=xmtx1(ie1)+factr1*val1
           ENDDO
        ENDDO
        xmtx2=0.
        DO p2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE 
           j2min=ABS((all_orbit%jj(p2)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
           DO jtot2=j2min,j2max
              factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
              DO h2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE       
                 DO h3=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE 
                    nshell3=all_orbit%nshell(h3)+all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=wcn+all_orbit%e(h3)+all_orbit%e(h2)- &
                         all_orbit%e(p1)-all_orbit%e(p2)
                    CALL pphhmtx(h1,p2,h2,h3,jtot2,ans2); IF(ans2(1).eq.0.) CYCLE 
                    CALL pphhmtx(h2,h3,p1,p2,jtot2,ans3); IF(ans3(1).eq.0.) CYCLE 
                    w2=all_orbit%evalence(c)+all_orbit%e(h1)+all_orbit%e(h2)+ &
                         all_orbit%e(h3)-all_orbit%evalence(a)-all_orbit%e(p1)+wcn
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+all_orbit%e(h3)- &
                         all_orbit%evalence(a)+wcn
                    de=den1*den2
                    DO ie2=1, n_startenergy_veff
                       CALL interpolate(w2(ie2),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie2),e_start_g,ans3,val3)
                       xmtx2(ie2)=xmtx2(ie2)+factr2*val2*val3/de(ie2)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        onebody_diagram_15=onebody_diagram_15-0.5*xmtx1*xmtx2
     ENDDO
  ENDDO

END SUBROUTINE diagram_15
!
!
!
SUBROUTINE diagram_16(a,c,onebody_diagram_16)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, h, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_16
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_16=0.
  DO h=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              CALL cross_coupled_mtxel1(a,h1,h,p1,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                    nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                    iosc1=nshell2-nshell1
                    nshell3=all_orbit%nshell(h)+all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(a)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2).or.dencheck(iosc1))CYCLE
                    den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                         all_orbit%e(p1)-all_orbit%e(p2)
                    den2=all_orbit%e(h2)+all_orbit%e(h)- &
                         all_orbit%evalence(a)-all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(h,h2,c,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(p1,p2,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+  &
                         all_orbit%jj(h)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                         all_orbit%jj(h)-jtot)
                    w1=all_orbit%e(h)+all_orbit%e(h2)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w3=w2+all_orbit%e(h)-all_orbit%e(p2)
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_16(i)=onebody_diagram_16(i)- &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_16


SUBROUTINE diagram_17(a,c,onebody_diagram_17)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, h, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_17
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_17=0.
  DO h=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              CALL cross_coupled_mtxel1(h,p1,c,h1,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)+ &
                         all_orbit%nshell(h2)
                    nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)+ &
                         all_orbit%nshell(a)
                    iosc1=nshell2-nshell1
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(h)
                    nshell4=all_orbit%nshell(p2)+all_orbit%nshell(a)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2).or.dencheck(iosc1))CYCLE
                    den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p2)+all_orbit%evalence(c)-all_orbit%evalence(a)
                    den2=all_orbit%e(h2)+all_orbit%e(h)-all_orbit%evalence(a)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(h1,h2,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(a,p2,h,h2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+  &
                         all_orbit%jj(h)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                         all_orbit%jj(h)-jtot)
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=w1+all_orbit%e(h)-all_orbit%e(p2)
                    w3=all_orbit%e(h)+all_orbit%e(h2)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_17(i)=onebody_diagram_17(i)- &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_17


SUBROUTINE diagram_18(a,c,onebody_diagram_18)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, p, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_18
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_18=0.
  DO p=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              CALL cross_coupled_mtxel1(p,h1,c,p1,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                    nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                    iosc1=nshell2-nshell1
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(c)
                    nshell4=all_orbit%nshell(p)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2).or.dencheck(iosc1))CYCLE
                    den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p2)
                    den2=all_orbit%evalence(c)+all_orbit%e(h2)-all_orbit%e(p)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(p1,p2,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(a,h2,p,p2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+  &
                         all_orbit%jj(p)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(p)+all_orbit%jj(h1)+ &
                         all_orbit%jj(p2)-jtot)
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)- &
                         all_orbit%e(p2)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_18(i)=onebody_diagram_18(i)+ &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_18



SUBROUTINE diagram_19(a,c,onebody_diagram_19)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, p, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_19
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_19=0.
  DO p=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              CALL cross_coupled_mtxel1(a,p1,p,h1,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                    nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                    iosc1=nshell2-nshell1
                    nshell3=all_orbit%nshell(a)+all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(p)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2).or.dencheck(iosc1))CYCLE
                    den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p2)
                    den2=all_orbit%evalence(a)+all_orbit%e(h2)-all_orbit%e(p)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(h1,h2,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(p,p2,c,h2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+ &
                         all_orbit%jj(p)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(h1)+all_orbit%jj(p2)+ &
                         all_orbit%jj(p)-jtot)
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=w1+all_orbit%evalence(a)-all_orbit%e(p2)
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_19(i)=onebody_diagram_19(i)+ &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_19



SUBROUTINE diagram_20(a,c,onebody_diagram_20)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, p, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_20
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_20=0.
  DO p=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(c)+all_orbit%nshell(h1)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%evalence(c)+all_orbit%e(h1)-all_orbit%e(p)- &
                   all_orbit%e(p1)
              CALL cross_coupled_mtxel1(p,p1,c,h1,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(c)
                    nshell4=all_orbit%nshell(p)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%evalence(c)+all_orbit%e(h2)-all_orbit%e(p)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(a,h2,p,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(h1,p2,p1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+ &
                         all_orbit%jj(p)+all_orbit%jj(p2))/2+  &
                         all_orbit%jj(h1)+all_orbit%jj(p2)+ &
                         all_orbit%jj(p)-jtot)
                    w1=all_orbit%e(h2)+all_orbit%evalence(c)+wcn
                    w2=w1+all_orbit%e(h1)-all_orbit%e(p)
                    w3=all_orbit%evalence(c)+all_orbit%e(h1)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_20(i)=onebody_diagram_20(i)+ &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_20



SUBROUTINE diagram_21(a,c,onebody_diagram_21)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, h, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_21
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_21=0.
  DO h=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(h1)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h)+all_orbit%e(h1)-all_orbit%evalence(a)- &
                   all_orbit%e(p1)
              CALL cross_coupled_mtxel1(a,p1,h,h1,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h)+all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(a)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%e(h)+all_orbit%e(h2)-all_orbit%evalence(a)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(h1,p2,p1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(h,h2,c,p2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+ &
                         all_orbit%jj(h)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                         all_orbit%jj(h)-jtot)
                    w1=all_orbit%e(h)+all_orbit%e(h1)+wcn
                    w2=w1+all_orbit%e(h2)-all_orbit%evalence(a)
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+all_orbit%e(h)- &
                         all_orbit%evalence(a)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_21(i)=onebody_diagram_21(i)- &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_21
!
!
!     Begin two-body diagrams
!
!
!
!     G-matrix contribution
!
SUBROUTINE diag1(a,b,c,d,jtot,dg1)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: i
  REAL(DP) :: val
  REAL(DP), DIMENSION(n_startenergy_g) :: ans
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: dg1
  REAL(DP), DIMENSION(n_startenergy_veff ) :: w

  dg1=0.
  CALL pphhmtx(a,b,c,d,jtot,ans)
  IF ( ans (1) == 0.0_dp) RETURN
  IF ( type_of_interaction == 'coupled-cluster') THEN
     w=all_orbit%e(c)+all_orbit%e(d)+wcn
  ELSEIF ( type_of_interaction == 'open-diagrams') THEN
     w=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
  ENDIF
  DO i=1, n_startenergy_veff
     CALL interpolate(w(i),e_start_g,ans,val)
     dg1(i)=val
  ENDDO

END SUBROUTINE diag1
!
!     Two-body core-polarization diagram, setting up all possible contributions
!
SUBROUTINE core_polarization_sum(a,b,c,d,jtot,sum)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: iph
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: sum
  REAL(DP), DIMENSION(n_startenergy_veff ) :: d2a, d2b, d2c, d2d, dia2
  REAL(DP) ::       sgabcd, sgcdjt, sgabjt

  sgabjt=iph((all_orbit%jj(a)+all_orbit%jj(b))/2-jtot+1)      
  sgcdjt=iph((all_orbit%jj(c)+all_orbit%jj(d))/2-jtot+1)      
  sgabcd=sgcdjt*sgabjt
  sum=0.0_dp; d2a=0.0_dp;   d2b=0.0_dp;    d2c=0.0_dp;   d2d=0.
  CALL diag2(a,b,c,d,jtot,dia2)
  d2a=dia2
  IF(c == d) THEN
     d2b=sgcdjt*d2a
  ELSE
     CALL diag2(a,b,d,c,jtot,dia2)
     d2b=sgcdjt*dia2
  ENDIF
  IF(a == b) THEN
     d2c=sgabjt*d2b
  ELSE
     CALL diag2(b,a,d,c,jtot,dia2)
     d2c=sgabcd*dia2
  ENDIF
  IF(a == b)THEN
     d2d=sgabjt*d2a
  ELSE
     CALL diag2(b,a,c,d,jtot,dia2)
     d2d=sgabjt*dia2
  ENDIF
  sum=d2a+d2b+d2c+d2d

END SUBROUTINE core_polarization_sum
!
!     The core-polarization diagram itself
!
SUBROUTINE diag2(a,b,c,d,jtot,dia2)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: iad, jp, jh, jph, jphmin, jphmax, ie, &
       nshell1, nshell2, idiff, iphx,j1max, j2max, j1min, j2min, iph
  REAL(DP) :: val1, val2, sixj, sg
  REAL(DP), DIMENSION(n_startenergy_veff )  :: dia2, de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  iad=(all_orbit%jj(a)+all_orbit%jj(d))/2+jtot
  dia2=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  j2max=(all_orbit%jj(b)+all_orbit%jj(d))/2
  jphmin=MAX(j1min,j2min)
  jphmax=MIN(j1max,j2max)
  DO jph=jphmin,jphmax
     sixj=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     DO jh=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(jh) /= 'hole' ) CYCLE
        DO jp=1,  all_orbit%total_orbits 
           IF(all_orbit%orbit_status(jp) == 'hole') CYCLE
           nshell1=all_orbit%nshell(d)+all_orbit%nshell(jh)
           nshell2=all_orbit%nshell(jp)+all_orbit%nshell(b)
           idiff=nshell1-nshell2
           IF(dencheck(idiff)) CYCLE
           CALL cross_coupled_mtxel1(a,jh,c,jp,jph,ans1)
           CALL cross_coupled_mtxel1(jp,b,jh,d,jph,ans2)
           IF ((ans1(1) == 0.)) CYCLE
           IF ((ans2(1) == 0.)) CYCLE
           iphx=(all_orbit%jj(jh)+all_orbit%jj(jp))/2
           sg=iph(iad+iphx+all_orbit%jj(jh))*sixj
           de=all_orbit%evalence(d)+all_orbit%e(jh)- &
                all_orbit%e(jp)-all_orbit%evalence(b)+wcn
           w1=all_orbit%e(jh)+all_orbit%evalence(c)+all_orbit%evalence(d)- &
                all_orbit%evalence(b)+wcn
           w2=all_orbit%e(jh)+all_orbit%evalence(d)+wcn
           DO ie=1, n_startenergy_veff
              CALL interpolate(w1(ie),e_start_g,ans1,val1)
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              dia2(ie)=dia2(ie)+sg*val1*val2/de(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag2
!
!     Two-body 4p-2h diagram
!
SUBROUTINE hole_hole_ladder(a,b,c,d,jtot,dg3)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: jh1, jh2, ie, nshell1, nshell2, idiff
  REAL(DP) :: val1, val2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: dg3
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  dg3=0.
  DO jh1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(jh1) /= 'hole' ) CYCLE
     DO jh2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(jh2) /= 'hole' ) CYCLE
        nshell1=all_orbit%nshell(jh1)+all_orbit%nshell(jh2)
        nshell2=all_orbit%nshell(a)+all_orbit%nshell(b)
        idiff=nshell2-nshell1
        IF(dencheck(idiff)) CYCLE
        CALL pphhmtx(a,b,jh1,jh2,jtot,ans1)
        IF ((ans1(1) == 0.)) CYCLE
        CALL pphhmtx(jh1,jh2,c,d,jtot,ans2)
        IF ((ans2(1) == 0.)) CYCLE
        de=all_orbit%e(jh1)+all_orbit%e(jh2)- &
             all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
        w1=all_orbit%e(jh1)+all_orbit%e(jh2)+wcn
        w2=w1+all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%evalence(a)- &
             all_orbit%evalence(b)
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           dg3(ie)=dg3(ie)+val1*val2*0.5/de(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE hole_hole_ladder
!
!     Two-body particle-particle ladder diagram to second order in G
!
SUBROUTINE particle_particle_ladder(a,b,c,d,jtot,dg4)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: jh1, jh2, ie, nshell1, nshell2, idiff
  REAL(DP) :: val1, val2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: dg4
  REAL(DP), DIMENSION(n_startenergy_veff )  :: de, w
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  dg4=0.
  !     Loop over particle-particle intermediate configurations
  DO jh1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(jh1) == 'hole' ) CYCLE
     DO jh2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(jh2) == 'hole' ) CYCLE
        nshell1=all_orbit%nshell(jh1)+all_orbit%nshell(jh2)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(d)
        idiff=nshell2-nshell1
        IF(dencheck(idiff)) CYCLE
        CALL pphhmtx(a,b,jh1,jh2,jtot,ans1)
        IF ((ans1(1) == 0.)) CYCLE
        CALL pphhmtx(jh1,jh2,c,d,jtot,ans2)
        IF ((ans2(1) == 0.)) CYCLE
        de=all_orbit%evalence(c)+all_orbit%evalence(d)- &
             all_orbit%e(jh1)-all_orbit%e(jh2)+wcn
        w=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
        DO ie=1,n_startenergy_veff
           CALL interpolate(w(ie),e_start_g,ans1,val1)
           CALL interpolate(w(ie),e_start_g,ans2,val2)
           dg4(ie)=dg4(ie)+val1*val2*0.5/de(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE particle_particle_ladder
!
!     summuning all non-folded third-order diagrams
!
SUBROUTINE number_conserving_set(a,b,c,d,jtot,sum)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: iph
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: sum

  REAL(DP), DIMENSION(n_startenergy_veff ) :: dg4,d4a, &
       d4b,d4c,d4d,d4,dg5,d5a,d5b,d5c,d5d,d5,dg7, &
       d7a,d7b,d7c,d7d,d7,dg8,d8a, d8b,d8c,d8d,d8, &
       dg9a,dg9b,d91a,d91d,d92a,d92b,d9a,d9b, &
       dg6a,d61a,d61b,d6a,dg6b,d62a,d62d,d6b
  REAL(DP) ::       sgabcd, sgcdjt, sgabjt

  sgabjt=iph((all_orbit%jj(a)+all_orbit%jj(b))/2-jtot+1)      
  sgcdjt=iph((all_orbit%jj(c)+all_orbit%jj(d))/2-jtot+1)      
  sgabcd=sgcdjt*sgabjt
  sum=0.0_dp
  dg4=0.; d4a=0.;d4b=0.;d4c=0.;d4d=0.;d4=0.;
  dg5=0.;d5a=0.;d5b=0.;d5c=0.;d5d=0.;d5=0.;
  dg7=0.; d7a=0.;d7b=0.;d7c=0.;d7d=0.;d7=0.;
  dg8=0.;d8a=0.; d8b=0.;d8c=0.;d8d=0.;d8=0.;
  dg9a=0.;dg9b=0.;d91a=0.;d91d=0.;
  d92a=0.;d92b=0.;d9a=0.;d9b=0.;
  dg6a=0.;d61a=0.;d61b=0.;d6a=0.;
  dg6b=0.;d62a=0.;d62d=0.;d6b=0.;
  CALL diag4(a,b,c,d,jtot,d4)
  CALL diag5(a,b,c,d,jtot,d5)
  CALL diag6a(a,b,c,d,jtot,d6a)
  CALL diag6b(a,b,c,d,jtot,d6b)
  CALL diag7(a,b,c,d,jtot,d7)
  CALL diag8(a,b,c,d,jtot,d8)
  CALL diag9a(a,b,c,d,jtot,d9a)
  CALL diag9b(a,b,c,d,jtot,d9b)
  d4a=d4
  d5a=d5
  d61a=d6a
  d62a=d6b
  d7a=d7
  d8a=d8
  d91a=d9a
  d92a=d9b
  IF(c == d) then
     d4b=sgcdjt*d4a
     d5b=sgcdjt*d5a
     d61b=sgcdjt*d61a
     d7b=sgcdjt*d7a
     d8b=sgcdjt*d8a
     d92b=sgcdjt*d92a
  ELSE
     CALL diag4(a,b,d,c,jtot,d4)
     CALL diag5(a,b,d,c,jtot,d5)
     CALL diag6a(a,b,d,c,jtot,d6a)
     CALL diag7(a,b,d,c,jtot,d7)
     CALL diag8(a,b,d,c,jtot,d8)
     CALL diag9b(a,b,d,c,jtot,d9b)
     d4b=sgcdjt*d4
     d5b=sgcdjt*d5
     d61b=sgcdjt*d6a
     d7b=sgcdjt*d7
     d8b=sgcdjt*d8
     d92b=sgcdjt*d9b
  ENDIF
  IF(a == b) then
     d4c=sgabjt*d4b
     d5c=sgabjt*d5b
     d7c=sgabjt*d7b
     d8c=sgabjt*d8b
  ELSE
     CALL diag4(b,a,d,c,jtot,d4)
     CALL diag5(b,a,d,c,jtot,d5)
     CALL diag7(b,a,d,c,jtot,d7)
     CALL diag8(b,a,d,c,jtot,d8)
     d4c=sgabcd*d4
     d5c=sgabcd*d5
     d7c=sgabcd*d7
     d8c=sgabcd*d8
  ENDIF
  IF(a == b) then
     d4d=sgabjt*d4a
     d5d=sgabjt*d5a
     d62d=sgabjt*d62a
     d7d=sgabjt*d7a
     d8d=sgabjt*d8a
     d91d=sgabjt*d91a
  ELSE
     CALL diag4(b,a,c,d,jtot,d4)
     CALL diag5(b,a,c,d,jtot,d5)
     CALL diag6b(b,a,c,d,jtot,d6b)
     CALL diag7(b,a,c,d,jtot,d7)
     CALL diag8(b,a,c,d,jtot,d8)
     CALL diag9a(b,a,c,d,jtot,d9a)
     d4d=sgabjt*d4
     d5d=sgabjt*d5
     d62d=sgabjt*d6b
     d7d=sgabjt*d7
     d8d=sgabjt*d8
     d91d=sgabjt*d9a
  ENDIF
  dg4=d4a+d4b+d4c+d4d
  dg5=d5a+d5b+d5c+d5d
  dg6a=d61a+d61b
  dg6b=d62a+d62d
  dg7=d7a+d7b+d7c+d7d
  dg8=d8a+d8b+d8c+d8d
  dg9a=d91a+d91d
  dg9b=d92a+d92b
  sum=dg4+dg5+dg7+dg8+dg9a+dg9b+dg6a+dg6b
1001 FORMAT(2X,20HDIAGRAM 4            ,11F8.4)
1002 FORMAT(2X,20HDIAGRAM 5            ,11F8.4)
1003 FORMAT(2X,20HDIAGRAM 6(1)         ,11F8.4)
1004 FORMAT(2X,20HDIAGRAM 6(2)         ,11F8.4)
1005 FORMAT(2X,20HDIAGRAM 7            ,11F8.4)
1006 FORMAT(2X,20HDIAGRAM 8            ,11F8.4)
1007 FORMAT(2X,20HDIAGRAM 9(1)         ,11F8.4)
1008 FORMAT(2X,20HDIAGRAM 9(2)         ,11F8.4)
  WRITE(6,1001) DG4
  WRITE(6,1002) DG5
  WRITE(6,1003) DG6A
  WRITE(6,1004) DG6B
  WRITE(6,1005) DG7
  WRITE(6,1006) DG8
  WRITE(6,1007) DG9A
  WRITE(6,1008) DG9B

END SUBROUTINE number_conserving_set
!
!
!
SUBROUTINE ladder_corepol(a,b,c,d,jtot,sum)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: iph
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: sum

  REAL(DP), DIMENSION(n_startenergy_veff ) :: d121a, &
       d121b,d12a,dg12a,d122a,d122d,dg12b,d12b,d131a, &
       d131d,d132a,d132b,dg13a,dg13b,d13a,d13b,d141a, &
       d141d,d142a,d142b,dg14a,dg14b,d14a,d14b
  REAL(DP) ::       sgabcd, sgcdjt, sgabjt

  sgabjt=iph((all_orbit%jj(a)+all_orbit%jj(b))/2-jtot+1)      
  sgcdjt=iph((all_orbit%jj(c)+all_orbit%jj(d))/2-jtot+1)      
  sgabcd=sgcdjt*sgabjt
  sum=0.
  dg12a=0.d0
  d121a=0.d0
  d121b=0.d0
  dg12b=0.d0
  d122a=0.d0
  d122d=0.d0
  dg13a=0.d0
  d131a=0.d0
  d131d=0.d0
  dg13b=0.d0
  d132a=0.d0
  d132b=0.d0
  dg14a=0.d0
  d141a=0.d0
  d141d=0.d0
  dg14b=0.d0
  d142a=0.d0
  d142b=0.d0
  CALL diag12a(a,b,c,d,jtot,d12a)
  CALL diag12b(a,b,c,d,jtot,d12b)
  CALL diag13a(a,b,c,d,jtot,d13a)
  CALL diag13b(a,b,c,d,jtot,d13b)
  CALL diag14a(a,b,c,d,jtot,d14a)
  CALL diag14b(a,b,c,d,jtot,d14b)
  d121a=d12a
  d122a=d12b
  d131a=d13a
  d132a=d13b
  d141a=d14a
  d142a=d14b
  IF(c == d) then
     d121b=sgcdjt*d121a
     d132b=sgcdjt*d132a
     d142b=sgcdjt*d142a
  ELSE
     CALL diag12a(a,b,d,c,jtot,d12a)
     CALL diag13b(a,b,d,c,jtot,d13b)
     CALL diag14b(a,b,d,c,jtot,d14b)
     d121b=sgcdjt*d12a
     d132b=sgcdjt*d13b
     d142b=sgcdjt*d14b
  ENDIF
  IF(a == b) then
     d122d=sgabjt*d122a
     d131d=sgabjt*d131a
     d141d=sgabjt*d141a
  ELSE
     CALL diag12b(b,a,c,d,jtot,d12b)
     CALL diag13a(b,a,c,d,jtot,d13a)
     CALL diag14a(b,a,c,d,jtot,d14a)
     d122d=sgabjt*d12b
     d131d=sgabjt*d13a
     d141d=sgabjt*d14a
  ENDIF
  dg12a=d121a+d121b
  dg12b=d122a+d122d
  dg13a=d131a+d131d
  dg13b=d132a+d132b
  dg14a=d141a+d141d
  dg14b=d142a+d142b
  sum=dg12a+dg12b+dg13a+dg13b +dg14a+dg14b
1001 FORMAT(2X,20HDIAGRAM 12(1)        ,11F8.4)
1002 FORMAT(2X,20HDIAGRAM 12(2)        ,11F8.4)
1003 FORMAT(2X,20HDIAGRAM 13(1)        ,11F8.4)
1004 FORMAT(2X,20HDIAGRAM 13(2)        ,11F8.4)
1005 FORMAT(2X,20HDIAGRAM 14(1)        ,11F8.4)
1006 FORMAT(2X,20HDIAGRAM 14(2)        ,11F8.4)
  WRITE(6,1001) DG12A
  WRITE(6,1002) DG12B
  WRITE(6,1003) DG13A
  WRITE(6,1004) DG13B
  WRITE(6,1005) DG14A
  WRITE(6,1006) DG14B

END SUBROUTINE ladder_corepol
!
!
!
SUBROUTINE tda_rpa_diagrams(a,b,c,d,jtot,sum)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: iph
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: sum
  REAL(DP), DIMENSION(n_startenergy_veff ) ::  dg15, &
       d15a,d15b,d15c,d15d,d15,dg16a,d161a,d161b,d161c,d161d, &
       d16a,dg16b,d162a,d162b,d162c,d162d,d16b
  REAL(DP) ::       sgabcd, sgcdjt, sgabjt

  sgabjt=iph((all_orbit%jj(a)+all_orbit%jj(b))/2-jtot+1)      
  sgcdjt=iph((all_orbit%jj(c)+all_orbit%jj(d))/2-jtot+1)      
  sgabcd=sgcdjt*sgabjt
  sum=0.
  dg15=0.d0
  d15a=0.d0
  d15b=0.d0
  d15c=0.d0
  d15d=0.d0
  dg16a=0.d0
  d161a=0.d0
  d161b=0.d0
  d161c=0.d0
  d161d=0.d0
  dg16b=0.d0
  d162a=0.d0
  d162b=0.d0
  d162c=0.d0
  d162d=0.d0
  CALL diag15(a,b,c,d,jtot,d15)
  CALL diag16a(a,b,c,d,jtot,d16a)
  CALL diag16b(a,b,c,d,jtot,d16b)
  d15a=d15
  d161a=d16a
  d162a=d16b
  IF(c == d) then
     d15b=sgcdjt*d15a
     d161b=sgcdjt*d161a
     d162b=sgcdjt*d162a
  ELSE
     CALL diag15(a,b,d,c,jtot,d15)
     CALL diag16a(a,b,d,c,jtot,d16a)
     CALL diag16b(a,b,d,c,jtot,d16b)
     d15b=sgcdjt*d15
     d161b=sgcdjt*d16a
     d162b=sgcdjt*d16b
  ENDIF
  IF(a == b) then
     d15c=sgabjt*d15b
     d161c=sgabjt*d161b
     d162c=sgabjt*d162b
  ELSE
     CALL diag15(b,a,d,c,jtot,d15)
     CALL diag16a(b,a,d,c,jtot,d16a)
     CALL diag16b(b,a,d,c,jtot,d16b)
     d15c=sgabcd*d15
     d161c=sgabcd*d16a
     d162c=sgabcd*d16b
  ENDIF
  IF(a == b) then
     d15d=sgabjt*d15a
     d161d=sgabjt*d161a
     d162d=sgabjt*d162a
  ELSE
     CALL diag15(b,a,c,d,jtot,d15)
     CALL diag16a(b,a,c,d,jtot,d16a)
     CALL diag16b(b,a,c,d,jtot,d16b)
     d15d=sgabjt*d15
     d161d=sgabjt*d16a
     d162d=sgabjt*d16b
  ENDIF
  dg15=d15a+d15b+d15c+d15d
  dg16a=d161a+d161b+d161c+d161d
  dg16b=d162a+d162b+d162c+d162d
  sum=dg15+dg16a+dg16b
1001 FORMAT(2X,20H3RD-ORDER TDA       ,11F8.4)
1002 FORMAT(2X,20HDIAGRAM 16(1)       ,11F8.4)
1003 FORMAT(2X,20HDIAGRAM 16(2)       ,11F8.4)
  WRITE(6,1001) DG15
  WRITE(6,1002) DG16A
  WRITE(6,1003) DG16B

END SUBROUTINE tda_rpa_diagrams
!
!
!
SUBROUTINE screening_thirdorder(a,b,c,d,jtot,sum)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: iph
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: sum
  REAL(DP), DIMENSION(n_startenergy_veff ) :: dg17a, &
       d17a,d171a,d171b,d171c,d171d, dg17b,d17b,d172a,d172b, &
       d172c,d172d,dg18a,d18a,d181a,d181b,d181c,d181d, &
       dg18b,d18b,d182a,d182b,d182c,d182d,dg19a,d19a, &
       d191a,d191b,d191c,d191d,dg19b,d19b,d192a,d192b, &
       d192c,d192d,dg20a,d20a,d201a,d201b,d201c,d201d, &
       dg20b,d20b,d202a,d202b,d202c,d202d,dg21a,d211a, &
       d211d,d21a,dg21b,d212a,d212b,d21b
  REAL(DP) ::       sgabcd, sgcdjt, sgabjt

  sgabjt=iph((all_orbit%jj(a)+all_orbit%jj(b))/2-jtot+1)      
  sgcdjt=iph((all_orbit%jj(c)+all_orbit%jj(d))/2-jtot+1)      
  sgabcd=sgcdjt*sgabjt
  sum=0.
  CALL diag17a(a,b,c,d,jtot,d17a)
  CALL diag17b(a,b,c,d,jtot,d17b)
  CALL diag18a(a,b,c,d,jtot,d18a)
  CALL diag18b(a,b,c,d,jtot,d18b)
  CALL diag19a(a,b,c,d,jtot,d19a)
  CALL diag19b(a,b,c,d,jtot,d19b)
  CALL diag20a(a,b,c,d,jtot,d20a)
  CALL diag20b(a,b,c,d,jtot,d20b)
  CALL diag21a(a,b,c,d,jtot,d21a)
  CALL diag21b(a,b,c,d,jtot,d21b)
  d171a=d17a
  d172a=d17b
  d181a=d18a
  d182a=d18b
  d191a=d19a
  d192a=d19b
  d201a=d20a
  d202a=d20b
  d211a=d21a
  d212a=d21b
  IF(c == d) then
     d171b=sgcdjt*d171a
     d172b=sgcdjt*d172a
     d181b=sgcdjt*d181a
     d182b=sgcdjt*d182a
     d191b=sgcdjt*d191a
     d192b=sgcdjt*d192a
     d201b=sgcdjt*d201a
     d202b=sgcdjt*d202a
     d212b=sgcdjt*d212a
  ELSE
     CALL diag17a(a,b,d,c,jtot,d17a)
     CALL diag17b(a,b,d,c,jtot,d17b)
     CALL diag18a(a,b,d,c,jtot,d18a)
     CALL diag18b(a,b,d,c,jtot,d18b)
     CALL diag19a(a,b,d,c,jtot,d19a)
     CALL diag19b(a,b,d,c,jtot,d19b)
     CALL diag20a(a,b,d,c,jtot,d20a)
     CALL diag20b(a,b,d,c,jtot,d20b)
     CALL diag21b(a,b,d,c,jtot,d21b)
     d171b=sgcdjt*d17a
     d172b=sgcdjt*d17b
     d181b=sgcdjt*d18a
     d182b=sgcdjt*d18b
     d191b=sgcdjt*d19a
     d192b=sgcdjt*d19b
     d201b=sgcdjt*d20a
     d202b=sgcdjt*d20b
     d212b=sgcdjt*d21b
  ENDIF
  IF(a == b) then
     d171c=sgabjt*d171b
     d172c=sgabjt*d172b
     d181c=sgabjt*d181b
     d182c=sgabjt*d182b
     d191c=sgabjt*d191b
     d192c=sgabjt*d192b
     d201c=sgabjt*d201b
     d202c=sgabjt*d202b
  ELSE
     CALL diag17a(b,a,d,c,jtot,d17a)
     CALL diag17b(b,a,d,c,jtot,d17b)
     CALL diag18a(b,a,d,c,jtot,d18a)
     CALL diag18b(b,a,d,c,jtot,d18b)
     CALL diag19a(b,a,d,c,jtot,d19a)
     CALL diag19b(b,a,d,c,jtot,d19b)
     CALL diag20a(b,a,d,c,jtot,d20a)
     CALL diag20b(b,a,d,c,jtot,d20b)
     d171c=sgabcd*d17a
     d172c=sgabcd*d17b
     d181c=sgabcd*d18a
     d182c=sgabcd*d18b
     d191c=sgabcd*d19a
     d192c=sgabcd*d19b
     d201c=sgabcd*d20a
     d202c=sgabcd*d20b
  ENDIF
  IF(a == b) then
     d171d=sgabjt*d171a
     d172d=sgabjt*d172a
     d181d=sgabjt*d181a
     d182d=sgabjt*d182a
     d191d=sgabjt*d191a
     d192d=sgabjt*d192a
     d201d=sgabjt*d201a
     d202d=sgabjt*d202a
     d211d=sgabjt*d211a
  ELSE
     CALL diag17a(b,a,c,d,jtot,d17a)
     CALL diag17b(b,a,c,d,jtot,d17b)
     CALL diag18a(b,a,c,d,jtot,d18a)
     CALL diag18b(b,a,c,d,jtot,d18b)   
     CALL diag19a(b,a,c,d,jtot,d19a)
     CALL diag19b(b,a,c,d,jtot,d19b)
     CALL diag20a(b,a,c,d,jtot,d20a)
     CALL diag20b(b,a,c,d,jtot,d20b)
     CALL diag21a(b,a,c,d,jtot,d21a)
     d171d=sgabjt*d17a
     d172d=sgabjt*d17b
     d181d=sgabjt*d18a
     d182d=sgabjt*d18b
     d191d=sgabjt*d19a
     d192d=sgabjt*d19b
     d201d=sgabjt*d20a
     d202d=sgabjt*d20b
     d211d=sgabjt*d21a
  ENDIF
  dg17a=d171a+d171b+d171c+d171d
  dg17b=d172a+d172b+d172c+d172d
  dg18a=d181a+d181b+d181c+d181d
  dg18b=d182a+d182b+d182c+d182d
  dg19a=d191a+d191b+d191c+d191d
  dg19b=d192a+d192b+d192c+d192d
  dg20a=d201a+d201b+d201c+d201d
  dg20b=d202a+d202b+d202c+d202d
  dg21a=d211a+d211d
  dg21b=d212a+d212b
  sum=dg17a+dg17b+dg18a+dg18b+dg19a+dg19b+dg20a+dg20b+dg21a+dg21b
1001 FORMAT(2X,20HDIAGRAM 17(1)        ,11F8.4)
1002 FORMAT(2X,20HDIAGRAM 17(2)        ,11F8.4)
1003 FORMAT(2X,20HDIAGRAM 18(1)        ,11F8.4)
1004 FORMAT(2X,20HDIAGRAM 18(2)        ,11F8.4)
1005 FORMAT(2X,20HDIAGRAM 19(1)        ,11F8.4)
1006 FORMAT(2X,20HDIAGRAM 19(2)        ,11F8.4)
1007 FORMAT(2X,20HDIAGRAM 20(1)        ,11F8.4)
1008 FORMAT(2X,20HDIAGRAM 20(2)        ,11F8.4)
1009 FORMAT(2X,20HDIAGRAM 21(1)        ,11F8.4)
1010 FORMAT(2X,20HDIAGRAM 21(2)        ,11F8.4)
  WRITE(6,1001) DG17A
  WRITE(6,1002) DG17B
  WRITE(6,1003) DG18A
  WRITE(6,1004) DG18B
  WRITE(6,1005) DG19A
  WRITE(6,1006) DG19B
  WRITE(6,1007) DG20A
  WRITE(6,1008) DG20B
  WRITE(6,1009) DG21A
  WRITE(6,1010) DG21B

END SUBROUTINE screening_thirdorder
!
!   These are additional diagrams with onebody insertions
!
SUBROUTINE onebody_insertions(a,b,c,d,jtot,sum)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: iph
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: sum
  REAL(DP), DIMENSION(n_startenergy_veff ) :: dg27a, &
       d271a,d271b,d27a,dg27b,d272a,d272d,d27b, &
       dg26a, d261a,d261b,d26a,dg26b,d262a,d262d,d26b, &
       dg24a,d241a,d241b,d24a,dg24b,d242a,d242d, &
       d24b,dg25a,d251a,d251b,d25a,dg25b,d252a, &
       d252d,d25b

  REAL(DP) ::       sgabcd, sgcdjt, sgabjt

  sgabjt=iph((all_orbit%jj(a)+all_orbit%jj(b))/2-jtot+1)      
  sgcdjt=iph((all_orbit%jj(c)+all_orbit%jj(d))/2-jtot+1)      
  sgabcd=sgcdjt*sgabjt
  sum=0.
  CALL diag24a(a,b,c,d,jtot,d24a)
  CALL diag24b(a,b,c,d,jtot,d24b)
  CALL diag25a(a,b,c,d,jtot,d25a)
  CALL diag25b(a,b,c,d,jtot,d25b)
  CALL diag26a(a,b,c,d,jtot,d26a)
  CALL diag26b(a,b,c,d,jtot,d26b)
  CALL diag27a(a,b,c,d,jtot,d27a)
  CALL diag27b(a,b,c,d,jtot,d27b)
  d241a=d24a
  d242a=d24b
  d251a=d25a
  d252a=d25b
  d261a=d26a
  d262a=d26b
  d271a=d27a
  d272a=d27b
  IF(c == d) then
     d241b=sgcdjt*d241a
     d251b=sgcdjt*d251a
     d261b=sgcdjt*d261a
     d271b=sgcdjt*d271a
  ELSE
     CALL diag24a(a,b,d,c,jtot,d24a)
     CALL diag25a(a,b,d,c,jtot,d25a)
     CALL diag26a(a,b,d,c,jtot,d26a)
     CALL diag27a(a,b,d,c,jtot,d27a)
     d241b=sgcdjt*d24a
     d251b=sgcdjt*d25a
     d261b=sgcdjt*d26a
     d271b=sgcdjt*d27a
  ENDIF
  IF(a == b) then
     d242d=sgabjt*d242a
     d252d=sgabjt*d252a
     d262d=sgabjt*d262a
     d272d=sgabjt*d272a
  ELSE
     CALL diag24b(b,a,c,d,jtot,d24b)
     CALL diag25b(b,a,c,d,jtot,d25b)
     CALL diag26b(b,a,c,d,jtot,d26b)
     CALL diag27b(b,a,c,d,jtot,d27b)
     d242d=sgabjt*d24b
     d252d=sgabjt*d25b
     d262d=sgabjt*d26b
     d272d=sgabjt*d27b
  ENDIF
  dg24a=d241a+d241b
  dg24b=d242a+d242d
  dg25a=d251a+d251b
  dg25b=d252a+d252d
  dg26a=d261a+d261b
  dg26b=d262a+d262d 
  dg27a=d271a+d271b
  dg27b=d272a+d272d                                     
  sum=dg24a+dg24b+dg25a+dg25b+dg26a+dg26b+dg27a+dg27b
1001 FORMAT(2X,20HDIAGRAM 24(1)        ,11F8.4)
1002 FORMAT(2X,20HDIAGRAM 24(2)        ,11F8.4)
1003 FORMAT(2X,20HDIAGRAM 25(1)        ,11F8.4)
1004 FORMAT(2X,20HDIAGRAM 25(2)        ,11F8.4)
1005 FORMAT(2X,20HDIAGRAM 26(1)        ,11F8.4)
1006 FORMAT(2X,20HDIAGRAM 26(2)        ,11F8.4)
1007 FORMAT(2X,20HDIAGRAM 27(1)        ,11F8.4)
1008 FORMAT(2X,20HDIAGRAM 27(2)        ,11F8.4)
  WRITE(6,1001) DG24A
  WRITE(6,1002) DG24B
  WRITE(6,1003) DG25A
  WRITE(6,1004) DG25B
  WRITE(6,1005) DG26A
  WRITE(6,1006) DG26B
  WRITE(6,1007) DG27A
  WRITE(6,1008) DG27B

END SUBROUTINE onebody_insertions
!
!
!
SUBROUTINE diag4(c,d,a,b,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2, prdjt
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: g1, g2, g3
  INTEGER :: p2, h, p3, p1, nshell1, nshell2, idiff1, idiff2, nshell3, &
       nshell4,j1mx,j1mn,j2mx,j2mn,phmx2, phmn2, iph, &
       ne1, ne2, ne3, jj1, jj2, jj3, ic1, ph, j1max,j1min,&
       phmin, phmax, j3max, j3mn, j3mx, j3min,j2min, j2max, &
       jph, phmx3, phmn3, j1, j2, j3
  REAL(DP) :: val1, val2, val3, sg, sprd
  REAL(DP) :: gg1, gg2, gg3, ssj1, ssj2, ssj3
  DIMENSION :: gg1(0:30,n_startenergy_veff), &
       gg2(0:30,n_startenergy_veff), gg3(0:30,n_startenergy_veff), &
       ssj1(0:30,0:30),ssj2(0:30,0:30),ssj3(0:30,0:30)
  LOGICAL dencheck, triag

  diagram=0.
  sg=iph((all_orbit%jj(a)+all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(d))/2)
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE      
        j1mx=(all_orbit%jj(p1)+all_orbit%jj(p2))/2
        j1mn=ABS((all_orbit%jj(p1)-all_orbit%jj(p2))/2)
        j2mx=(all_orbit%jj(p2)+all_orbit%jj(d))/2
        j2mn=ABS((all_orbit%jj(p2)-all_orbit%jj(d))/2)
        phmx2=(all_orbit%jj(c)+all_orbit%jj(p2))/2
        phmn2=ABS((all_orbit%jj(c)-all_orbit%jj(p2))/2)
        DO p3=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE      
           j2max=MIN(j2mx,(all_orbit%jj(p3)+all_orbit%jj(b))/2)
           j2min=MAX(j2mn,ABS((all_orbit%jj(p3)-all_orbit%jj(b))/2))
           IF (j2max < j2min) CYCLE
           j3mx=(all_orbit%jj(p1)+all_orbit%jj(p3))/2
           j3mn=ABS((all_orbit%jj(p1)-all_orbit%jj(p3))/2)
           phmx3=(all_orbit%jj(a)+all_orbit%jj(p3))/2
           phmn3=ABS((all_orbit%jj(a)-all_orbit%jj(p3))/2)
           DO h=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p3)
              nshell3=all_orbit%nshell(h)+all_orbit%nshell(a)+ &
                   all_orbit%nshell(b)
              nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)+ &
                   all_orbit%nshell(d)
              idiff1=nshell1-nshell2
              idiff2=nshell3-nshell4
              IF(dencheck(idiff1).or.dencheck(idiff2))CYCLE
              den1=wcn+all_orbit%e(h)+all_orbit%evalence(a)-all_orbit%e(p1)- &
                   all_orbit%e(p3)
              den2=wcn+all_orbit%e(h)+all_orbit%evalence(a)+all_orbit%evalence(b)- &
                   all_orbit%evalence(d)-all_orbit%e(p1)-all_orbit%e(p2)
              de=den1*den2
              phmax=MIN((all_orbit%jj(h)+all_orbit%jj(p1))/2,phmx2,phmx3)
              phmin=MAX(ABS((all_orbit%jj(h)-all_orbit%jj(p1))/2),phmn2,phmn3)
              IF (phmax < phmin) CYCLE
              j1max=MIN(j1mx,(all_orbit%jj(c)+all_orbit%jj(h))/2)
              j1min=MAX(j1mn,ABS((all_orbit%jj(c)-all_orbit%jj(h))/2))
              IF (j1max < j1min) CYCLE
              j3max=MIN(j3mx,(all_orbit%jj(a)+all_orbit%jj(h))/2)
              j3min=MAX(j3mn,ABS((all_orbit%jj(a)-all_orbit%jj(h))/2))
              IF (j3max < j3min) CYCLE
              gg1=0.D0; ssj1=0.D0
              DO jj1=j1min,j1max
                 DO jph=phmin,phmax
                    ssj1(jj1,jph)=sjs(all_orbit%jj(h),all_orbit%jj(c), &
                         2*jj1,all_orbit%jj(p2),all_orbit%jj(p1),2*jph) &
                         *(2.*jph+1.)
                 ENDDO
                 CALL  pphhmtx(c,h,p1,p2,jj1,g1)
                 IF (g1(1) == 0.) CYCLE
                 w1=all_orbit%evalence(a)+all_orbit%e(h)+wcn+all_orbit%evalence(b)- &
                      all_orbit%evalence(d)
                 DO ne1=1,n_startenergy_veff
                    CALL  interpolate(w1(ne1),e_start_g,g1,val1)
                    gg1(jj1,ne1)=val1*(2*jj1+1.)
                 ENDDO
              ENDDO
              gg2=0.D0; ssj2=0.D0
              DO jj2=j2min,j2max
                 DO jph=phmin,phmax
                    IF(triag(jj2,jph,jtot)) CYCLE
                    ssj2(jj2,jph)=sjs(all_orbit%jj(c),2*jph, &
                         all_orbit%jj(p2),2*jj2,all_orbit%jj(d),2*jtot) &
                         *sjs(all_orbit%jj(a),2*jph,all_orbit%jj(p3), &
                         2*jj2,all_orbit%jj(b),2*jtot)
                 ENDDO
                 CALL pphhmtx(p2,d,p3,b,jj2,g2)
                 IF(g2(1) == 0.) CYCLE
                 w2=all_orbit%evalence(a)+all_orbit%evalence(b)+all_orbit%e(h)- &
                      all_orbit%e(p1)+wcn
                 DO ne2=1,n_startenergy_veff
                    CALL interpolate(w2(ne2),e_start_g,g2,val2)
                    gg2(jj2,ne2)=val2*(2*jj2+1.)
                 ENDDO
              ENDDO
              gg3=0.D0; ssj3=0.D0
              DO jj3=j3min,j3max
                 DO jph=phmin,phmax
                    ssj3(jj3,jph)=sjs(all_orbit%jj(p1),all_orbit%jj(p3), &
                         2*jj3,all_orbit%jj(a),all_orbit%jj(h),2*jph)
                 ENDDO
                 CALL pphhmtx(p1,p3,a,h,jj3,g3); IF(g3(1) == 0.) CYCLE
                 w3=all_orbit%e(h)+all_orbit%evalence(a)+wcn
                 DO ne3=1,n_startenergy_veff
                    CALL interpolate(w3(ne3),e_start_g,g3,val3)
                    gg3(jj3,ne3)=val3*(2*jj3+1.)
                 ENDDO
              ENDDO
              prdjt=0.
              DO j1=j1min,j1max
                 DO j2=j2min,j2max
                    DO j3=j3min,j3max
                       sprd=0.d0
                       DO ph=phmin,phmax
                          sprd=sprd+ssj1(j1,ph)*ssj2(j2,ph)*ssj3(j3,ph)
                       ENDDO
                       DO ic1=1,n_startenergy_veff
                          prdjt(ic1)=prdjt(ic1)+sprd*gg1(j1,ic1)* &
                               gg2(j2,ic1)*gg3(j3,ic1)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
              diagram=diagram+prdjt*sg/de
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag4
!
!
!
SUBROUTINE diag5(c,d,a,b,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2, prdjt
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: g1, g11, g2
  INTEGER :: p2, h1, h2, p1, nshell1, nshell2, idiff1, idiff2, nshell3, &
       nshell4,j1mx,j1mn,j2mx,j2mn,phmx2, phmn2, iph, ipph, ij1, &
       ne1, ne2, ne3, jj1, jj2, jj3, ij2, ij3, isgn,j1mn2, j1mx1, &
       j1mx3,j1mx2,j1mn1, j1mn3,j1min,j2min,j1max,j2max,j2mnn,j2mxx
  REAL(DP) :: val1, val2, val3, ffnorm, sprd, ssnj
  REAL(DP) :: gg1, gg2
  DIMENSION :: gg1(0:30,n_startenergy_veff), gg2(0:30,n_startenergy_veff)
  LOGICAL dencheck

  diagram=0.
  isgn=(all_orbit%jj(a)+all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(d))/2
  ffnorm=-iph(isgn)*0.5
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE      
        j1mx1=(all_orbit%jj(p1)+all_orbit%jj(p2))/2
        j1mn1=ABS((all_orbit%jj(p1)-all_orbit%jj(p2))/2)
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell1=all_orbit%nshell(h1)+all_orbit%nshell(a)+ &
                all_orbit%nshell(b)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)+ &
                all_orbit%nshell(d)
           idiff1=nshell1-nshell2
           IF(dencheck(idiff1))CYCLE
           den2=all_orbit%e(h1)+all_orbit%evalence(a)+all_orbit%evalence(b)- &
                all_orbit%e(p1)-all_orbit%e(p2)-all_orbit%evalence(d)+wcn
           j1mx2=(all_orbit%jj(c)+all_orbit%jj(h1))/2
           j1mn2=ABS((all_orbit%jj(c)-all_orbit%jj(h1))/2)
           j2mx=(all_orbit%jj(h1)+all_orbit%jj(b))/2
           j2mn=ABS((all_orbit%jj(h1)-all_orbit%jj(b))/2)
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              nshell3=all_orbit%nshell(h2)+all_orbit%nshell(a)
              nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff2=nshell3-nshell4
              IF(dencheck(idiff2))CYCLE
              den1=all_orbit%evalence(a)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)+wcn
              j1mx3=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              j1mn3=ABS((all_orbit%jj(a)-all_orbit%jj(h2))/2)
              j1max=MIN(j1mx1,j1mx2,j1mx3)
              j1min=MAX(j1mn1,j1mn2,j1mn3)
              IF (j1max < j1min) CYCLE
              j2mxx=(all_orbit%jj(h2)+all_orbit%jj(d))/2
              j2mnn=ABS((all_orbit%jj(h2)-all_orbit%jj(d))/2)
              j2max=MIN(j2mx,j2mxx)
              j2min=MAX(j2mn,j2mnn)
              IF (j2max < j2min) CYCLE
              de=den1*den2
              gg1=0.
              DO jj1=j1min,j1max
                 CALL pphhmtx(c,h1,p1,p2,jj1,g1)
                 CALL pphhmtx(p1,p2,a,h2,jj1,g11)
                 IF((g1(1) == 0.d0).or.(g11(1) == 0.d0)) CYCLE
                 w1=all_orbit%e(h1)+all_orbit%evalence(a)+all_orbit%evalence(b)- &
                      all_orbit%evalence(d)+wcn
                 w2=all_orbit%evalence(a)+all_orbit%e(h2)+wcn
                 DO ne1=1,n_startenergy_veff
                    CALL interpolate(w1(ne1),e_start_g,g1,val1)
                    CALL interpolate(w2(ne1),e_start_g,g11,val2)
                    gg1(jj1,ne1)=val1*val2*(2.*jj1+1.)
                 ENDDO
              ENDDO
              gg2=0.
              DO jj2=j2min,j2max
                 CALL pphhmtx(h2,d,h1,b,jj2,g2)
                 IF(g2(1) == 0.d0) CYCLE
                 w3=all_orbit%e(h2)+all_orbit%e(h1)+all_orbit%evalence(a)+ &
                      all_orbit%evalence(b)-all_orbit%e(p1)-all_orbit%e(p2)+wcn
                 DO ne2=1,n_startenergy_veff
                    CALL interpolate(w3(ne2),e_start_g,g2,val3)
                    gg2(jj2,ne2)=val3*(2.*jj2+1.)
                 ENDDO
              ENDDO
              prdjt=0.
              DO jj1=j1min,j1max
                 DO jj2=j2min,j2max
                    ssnj=  snj(all_orbit%jj(h2),2*jj1,all_orbit%jj(a), &
                         all_orbit%jj(d),all_orbit%jj(c),2*jtot, &
                         2*jj2,all_orbit%jj(h1),all_orbit%jj(b))
                    prdjt(:)=prdjt(:)+gg1(jj1,:)*ssnj*gg2(jj2,:)
                 ENDDO
              ENDDO
              diagram=diagram+ffnorm*prdjt/de
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag5
!
!
!
SUBROUTINE diag6a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, h2, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, pp, j_min, j_max, h1
  REAL(DP) :: val1, val2, val3, factr1, factr2
  LOGICAL dencheck

  diagram=0.
  factr1=-0.5/(all_orbit%jj(c)+1.)
  DO h2=1, all_orbit%total_orbits
     IF((all_orbit%jj(h2) /= all_orbit%jj(c)).or.  &
          (all_orbit%ll(h2) /= all_orbit%ll(c))) CYCLE
     IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(c))/2)
        j2min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j1max=(all_orbit%jj(h1)+all_orbit%jj(c))/2
        j2max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        j_min=max0(j1min,j2min)
        j_max=MIN(j1max,j2max)
        DO pp=j_min,j_max
           factr2=(2.*pp+1.)*factr1
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
                 nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff1=nshell1-nshell2
                 nshell3=all_orbit%nshell(h2)+all_orbit%nshell(d)+nshell1
                 nshell4=all_orbit%nshell(a)+all_orbit%nshell(b)+nshell2
                 idiff2=nshell3-nshell4
                 IF(dencheck(idiff1).or.dencheck(idiff2))CYCLE
                 den1=wcn+all_orbit%e(h1)+all_orbit%evalence(c)- &
                      all_orbit%e(p1)-all_orbit%e(p2)
                 den2=wcn+all_orbit%evalence(c)+all_orbit%evalence(d)+ &
                      all_orbit%e(h1)+all_orbit%e(h2)- &
                      all_orbit%e(p1)-all_orbit%e(p2)- &
                      all_orbit%evalence(a)-all_orbit%evalence(b)
                 CALL pphhmtx(a,b,h2,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
                 CALL pphhmtx(p1,p2,h1,c,pp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(h1,h2,p1,p2,pp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%e(h2)+all_orbit%evalence(d)+ &
                      all_orbit%e(h1)+all_orbit%evalence(c)- &
                      all_orbit%e(p1)-all_orbit%e(p2)+wcn
                 w2=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                 w3=w2+all_orbit%e(h2)+all_orbit%evalence(d)- &
                      all_orbit%evalence(a)-all_orbit%evalence(b)
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+factr2*val1*val2*val3/de(ie)   
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag6a
!
!
!
SUBROUTINE diag6b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, h2, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, pp, j_min, j_max, h1
  REAL(DP) :: val1, val2, val3, factr1, factr2
  LOGICAL dencheck

  diagram=0.
  factr1=-0.5/(all_orbit%jj(a)+1.)
  DO h2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
     IF((all_orbit%jj(h2) /= all_orbit%jj(a)).or. &
          (all_orbit%ll(h2) /= all_orbit%ll(a))) CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j1max=(all_orbit%jj(h1)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        DO pp=j_min,j_max
           factr2=(2*pp+1.)*factr1   
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                 nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff1=nshell1-nshell2
                 nshell3=all_orbit%nshell(c)+all_orbit%nshell(d)+nshell1
                 nshell4=all_orbit%nshell(h2)+all_orbit%nshell(b)+nshell2
                 idiff2=nshell3-nshell4
                 IF(dencheck(idiff1).or.dencheck(idiff2)) CYCLE
                 den1=wcn+all_orbit%e(h1)+all_orbit%e(h2)- &
                      all_orbit%e(p1)-all_orbit%e(p2)
                 den2=den1+all_orbit%evalence(c)+all_orbit%evalence(d)- &
                      all_orbit%e(h2)-all_orbit%evalence(b)
                 CALL pphhmtx(h2,b,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
                 CALL pphhmtx(h1,a,p1,p2,pp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(p1,p2,h1,h2,pp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                      all_orbit%evalence(c)+all_orbit%evalence(d)- &
                      all_orbit%e(p1)-all_orbit%e(p2)+wcn
                 w2=all_orbit%evalence(c)+all_orbit%evalence(d)+all_orbit%e(h1)-&
                      all_orbit%evalence(b)+wcn
                 w3=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+factr2*val1*val2*val3/de(ie)   
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag6b
!
!
!
SUBROUTINE diag6f(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, jv1, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, pp, j_min, j_max, h1
  REAL(DP) :: val1, val2, val3, factr1, factr2
  LOGICAL dencheck

  diagram=0.
  factr1=-0.5/(all_orbit%jj(c)+1.)
  DO jv1=1, all_orbit%total_orbits
     IF(all_orbit%model_space(jv1) /= 'inside' ) CYCLE
     IF(all_orbit%jj(jv1).ne.all_orbit%jj(a)) CYCLE
     CALL pphhmtx(jv1,b,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(jv1))/2)
        j2min=ABS((all_orbit%jj(h1)-all_orbit%jj(a))/2)
        j1max=(all_orbit%jj(h1)+all_orbit%jj(jv1))/2
        j2max=(all_orbit%jj(h1)+all_orbit%jj(a))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        DO pp=j_min,j_max
           factr2=(2*pp+1)*factr1
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 nshell1=all_orbit%nshell(h1)+all_orbit%nshell(jv1)
                 nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff1=nshell1-nshell2
                 nshell3=all_orbit%nshell(c)+all_orbit%nshell(d)+all_orbit%nshell(h1)
                 nshell4=all_orbit%nshell(b)+nshell2
                 idiff2=nshell3-nshell4
                 IF(dencheck(idiff1).or.dencheck(idiff2)) CYCLE
                 den1=wcn+all_orbit%e(h1)+all_orbit%evalence(jv1)- &
                      all_orbit%e(p1)-all_orbit%e(p2)
                 den2=wcn+all_orbit%evalence(c)+all_orbit%evalence(d)+ &
                      all_orbit%e(h1)-all_orbit%e(p1)- &
                      all_orbit%e(p2)-all_orbit%evalence(b)
                 CALL pphhmtx(p1,p2,h1,jv1,pp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(h1,a,p1,p2,pp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%evalence(c)+all_orbit%evalence(d)+ &
                      all_orbit%e(h1)+all_orbit%evalence(jv1)- &
                      all_orbit%e(p1)-all_orbit%e(p2)+wcn
                 w2=all_orbit%e(h1)+all_orbit%evalence(jv1)+wcn
                 w3=all_orbit%evalence(c)+all_orbit%evalence(d)+ &
                      all_orbit%e(h1)-all_orbit%evalence(b)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)-0.5*factr2*val1*val2*val3/de(ie)   
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag6f
!
!
!
SUBROUTINE diag7(c,d,a,b,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2, prdjt
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: g1, g11, g2
  INTEGER :: p2, h1, h2, p1, nshell1, nshell2, idiff1, idiff2, nshell3, &
       nshell4,j1mx,j1mn,j2mx,j2mn,phmx2, phmn2, iph, ipph, ij1, &
       ne1, ne2, ne3, jj1, jj2, jj3, ij2, ij3, isgn,j1mn2, j1mx1, &
       j1mx3,j1mx2,j1mn1, j1mn3,j1min,j2min,j1max,j2max,j2mnn,j2mxx
  REAL(DP) :: val1, val2, val3, ffnorm, sprd, ssnj
  REAL(DP) :: gg1, gg2
  DIMENSION :: gg1(0:30,n_startenergy_veff), gg2(0:30,n_startenergy_veff)
  LOGICAL dencheck

  diagram=0.
  isgn=(all_orbit%jj(a)+all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(d))/2
  ffnorm=-iph(isgn)*0.5
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) /= 'hole' ) CYCLE      
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) /= 'hole' ) CYCLE      
        j1mx1=(all_orbit%jj(p1)+all_orbit%jj(p2))/2
        j1mn1=ABS((all_orbit%jj(p1)-all_orbit%jj(p2))/2)
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) == 'hole' ) CYCLE
           nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1))CYCLE
           den1=wcn+all_orbit%e(p1)+all_orbit%e(p2)-all_orbit%e(h1)-&
                all_orbit%evalence(c)
           j1mx2=(all_orbit%jj(c)+all_orbit%jj(h1))/2
           j1mn2=ABS((all_orbit%jj(c)-all_orbit%jj(h1))/2)
           j2mx=(all_orbit%jj(h1)+all_orbit%jj(b))/2
           j2mn=ABS((all_orbit%jj(h1)-all_orbit%jj(b))/2)
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) == 'hole' ) CYCLE
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(b)+ &
                   all_orbit%nshell(p2)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(c)+ &
                   all_orbit%nshell(d)
              idiff2=nshell3-nshell4
              IF(dencheck(idiff2))CYCLE
              den2=wcn+all_orbit%e(p1)+all_orbit%e(p2)+all_orbit%evalence(b)-&
                   all_orbit%e(h2)-all_orbit%evalence(c)-all_orbit%evalence(d)
              j1mx3=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              j1mn3=ABS((all_orbit%jj(a)-all_orbit%jj(h2))/2)
              j1max=MIN(j1mx1,j1mx2,j1mx3)
              j1min=MAX(j1mn1,j1mn2,j1mn3)
              IF (j1max < j1min) CYCLE
              j2mxx=(all_orbit%jj(h2)+all_orbit%jj(d))/2
              j2mnn=ABS((all_orbit%jj(h2)-all_orbit%jj(d))/2)
              j2max=MIN(j2mx,j2mxx)
              j2min=MAX(j2mn,j2mnn)
              IF (j2max < j2min) CYCLE
              de=den1*den2
              gg1=0.
              DO jj1=j1min,j1max
                 CALL pphhmtx(c,h1,p1,p2,jj1,g1)
                 CALL pphhmtx(p1,p2,a,h2,jj1,g11)
                 IF((g1(1) == 0.d0).or.(g11(1) == 0.d0)) CYCLE
                 w1=all_orbit%e(p1)+all_orbit%e(p2)+wcn
                 w2=w1+all_orbit%evalence(a)+all_orbit%evalence(b)-all_orbit%evalence(c)- &
                      all_orbit%evalence(d)
                 DO ne1=1,n_startenergy_veff
                    CALL interpolate(w1(ne1),e_start_g,g1,val1)
                    CALL interpolate(w2(ne1),e_start_g,g11,val2)
                    gg1(jj1,ne1)=val1*val2*(2.*jj1+1.)
                 ENDDO
              ENDDO
              gg2=0.
              DO jj2=j2min,j2max
                 CALL pphhmtx(h2,d,h1,b,jj2,g2)
                 IF(g2(1) == 0.d0) CYCLE
                 w3=all_orbit%e(p1)+all_orbit%e(p2)+all_orbit%evalence(b)- &
                      all_orbit%evalence(c)+wcn
                 DO ne2=1,n_startenergy_veff
                    CALL interpolate(w3(ne2),e_start_g,g2,val3)
                    gg2(jj2,ne2)=val3*(2.*jj2+1.)
                 ENDDO
              ENDDO
              prdjt=0.
              DO jj1=j1min,j1max
                 DO jj2=j2min,j2max
                    ssnj=  snj(all_orbit%jj(h2),2*jj1,all_orbit%jj(a), &
                         all_orbit%jj(d),all_orbit%jj(c),2*jtot, &
                         2*jj2,all_orbit%jj(h1),all_orbit%jj(b))
                    prdjt(:)=prdjt(:)+gg1(jj1,:)*ssnj*gg2(jj2,:)
                 ENDDO
              ENDDO
              diagram=diagram+ffnorm*prdjt/de
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag7
!
!
!
SUBROUTINE diag8(c,d,a,b,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2, prdjt
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: g1, g2, g3
  INTEGER :: p2, h, p3, p1, nshell1, nshell2, idiff1, idiff2, nshell3, &
       nshell4,j1mx,j1mn,j2mx,j2mn,phmx2, phmn2, iph, ipph, ij1, &
       ne1, ne2, ne3, jj1, jj2, jj3, ij2, ij3, ic1, ph, j1max,j1min,&
       phmin, phmax, phmxx, j3max, j3mn, j3mx, j3min,j2min, j2max, &
       jph, phmx3, phmn3, j1, j2, j3
  REAL(DP) :: val1, val2, val3, sg, sprd
  REAL(DP) :: gg1, gg2, gg3, ssj1, ssj2, ssj3
  DIMENSION :: gg1(0:30,n_startenergy_veff), &
       gg2(0:30,n_startenergy_veff), gg3(0:30,n_startenergy_veff), &
       ssj1(0:30,0:30),ssj2(0:30,0:30),ssj3(0:30,0:30)
  LOGICAL dencheck, triag

  diagram=0.
  sg=iph((all_orbit%jj(a)+all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(d))/2)
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) /= 'hole' ) CYCLE      
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) /= 'hole' ) CYCLE      
        j1mx=(all_orbit%jj(p1)+all_orbit%jj(p2))/2
        j1mn=ABS((all_orbit%jj(p1)-all_orbit%jj(p2))/2)
        j2mx=(all_orbit%jj(p2)+all_orbit%jj(d))/2
        j2mn=ABS((all_orbit%jj(p2)-all_orbit%jj(d))/2)
        phmx2=(all_orbit%jj(c)+all_orbit%jj(p2))/2
        phmn2=ABS((all_orbit%jj(c)-all_orbit%jj(p2))/2)
        DO p3=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p3) /= 'hole' ) CYCLE      
           j2max=MIN(j2mx,(all_orbit%jj(p3)+all_orbit%jj(b))/2)
           j2min=MAX(j2mn,ABS((all_orbit%jj(p3)-all_orbit%jj(b))/2))
           IF (j2max < j2min) CYCLE
           j3mx=(all_orbit%jj(p1)+all_orbit%jj(p3))/2
           j3mn=ABS((all_orbit%jj(p1)-all_orbit%jj(p3))/2)
           phmx3=(all_orbit%jj(a)+all_orbit%jj(p3))/2
           phmn3=ABS((all_orbit%jj(a)-all_orbit%jj(p3))/2)
           DO h=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h) == 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p3)+ &
                   all_orbit%nshell(b)
              nshell4=all_orbit%nshell(h)+all_orbit%nshell(c)+ &
                   all_orbit%nshell(d)
              idiff1=nshell1-nshell2
              idiff2=nshell3-nshell4
              IF(dencheck(idiff1).OR.dencheck(idiff2))CYCLE
              den1=wcn+all_orbit%e(p1)+all_orbit%e(p2)-all_orbit%e(h)-&
                   all_orbit%evalence(c)
              den2=wcn+all_orbit%e(p1)+all_orbit%e(p3)+all_orbit%evalence(b)-&
                   all_orbit%e(h)-all_orbit%evalence(c)-all_orbit%evalence(d)
              de=den1*den2
              phmax=MIN((all_orbit%jj(h)+all_orbit%jj(p1))/2,phmx2,phmx3)
              phmin=MAX(ABS((all_orbit%jj(h)-all_orbit%jj(p1))/2),phmn2,phmn3)
              IF (phmax < phmin) CYCLE
              j1max=MIN(j1mx,(all_orbit%jj(c)+all_orbit%jj(h))/2)
              j1min=MAX(j1mn,ABS((all_orbit%jj(c)-all_orbit%jj(h))/2))
              IF (j1max < j1min) CYCLE
              j3max=MIN(j3mx,(all_orbit%jj(a)+all_orbit%jj(h))/2)
              j3min=MAX(j3mn,ABS((all_orbit%jj(a)-all_orbit%jj(h))/2))
              IF (j3max < j3min) CYCLE
              gg1=0.D0; ssj1= 0.D0
              DO jj1=j1min,j1max
                 DO jph=phmin,phmax
                    ssj1(jj1,jph)=sjs(all_orbit%jj(h),all_orbit%jj(c), &
                         2*jj1,all_orbit%jj(p2),all_orbit%jj(p1),2*jph) &
                         *(2.*jph+1.)
                 ENDDO
                 CALL  pphhmtx(c,h,p1,p2,jj1,g1)
                 IF (g1(1) == 0.) CYCLE
                 w1=all_orbit%e(p1)+all_orbit%e(p2)+wcn
                 DO ne1=1,n_startenergy_veff
                    CALL  interpolate(w1(ne1),e_start_g,g1,val1)
                    gg1(jj1,ne1)=val1*(2*jj1+1.)
                 ENDDO
              ENDDO
              gg2=0.D0; ssj2 = 0.D0
              DO jj2=j2min,j2max
                 DO jph=phmin,phmax
                    IF(triag(jj2,jph,jtot)) CYCLE
                    ssj2(jj2,jph)=sjs(all_orbit%jj(c),2*jph, &
                         all_orbit%jj(p2),2*jj2,all_orbit%jj(d),2*jtot) &
                         *sjs(all_orbit%jj(a),2*jph,all_orbit%jj(p3), &
                         2*jj2,all_orbit%jj(b),2*jtot)
                 ENDDO
                 CALL pphhmtx(p2,d,p3,b,jj2,g2)
                 IF(g2(1) == 0.) CYCLE
                 w2=all_orbit%e(p1)+all_orbit%e(p2)+all_orbit%e(p3)- &
                      all_orbit%e(h)-all_orbit%evalence(c)+all_orbit%evalence(b)+wcn
                 DO ne2=1,n_startenergy_veff
                    CALL interpolate(w2(ne2),e_start_g,g2,val2)
                    gg2(jj2,ne2)=val2*(2*jj2+1.)
                 ENDDO
              ENDDO
              gg3=0.0D0; ssj3 = 0.D0
              DO jj3=j3min,j3max
                 DO jph=phmin,phmax
                    ssj3(jj3,jph)=sjs(all_orbit%jj(p1),all_orbit%jj(p3), &
                         2*jj3,all_orbit%jj(a),all_orbit%jj(h),2*jph)
                 ENDDO
                 CALL pphhmtx(p1,p3,a,h,jj3,g3); IF(g3(1) == 0.) CYCLE
                 w3=all_orbit%e(p1)+all_orbit%e(p3)+all_orbit%evalence(b)+ &
                      all_orbit%evalence(a)-all_orbit%evalence(c)-all_orbit%evalence(d)+wcn
                 DO ne3=1,n_startenergy_veff
                    CALL interpolate(w3(ne3),e_start_g,g3,val3)
                    gg3(jj3,ne3)=val3*(2*jj3+1.)
                 ENDDO
              ENDDO
              prdjt=0.
              DO j1=j1min,j1max
                 DO j2=j2min,j2max
                    DO j3=j3min,j3max
                       sprd=0.d0
                       DO ph=phmin,phmax
                          sprd=sprd+ssj1(j1,ph)*ssj2(j2,ph)*ssj3(j3,ph)
                       ENDDO
                       prdjt(:)=prdjt(:)+sprd*gg1(j1,:)*gg2(j2,:)*gg3(j3,:)
                    ENDDO
                 ENDDO
              ENDDO
              diagram=diagram+prdjt*sg/de
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag8
!
!
!
SUBROUTINE diag9a(a,b,c,d,jt,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jt  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, h2, p1, nshell1, nshell2, diff1, &
       diff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, pp, j_min, j_max, h1, iph
  REAL(DP) :: val1, val2, val3, factr1, factr2
  LOGICAL dencheck

  diagram=0.
  factr1=-0.5/(all_orbit%jj(a)+1.)
  DO p2=1, all_orbit%total_orbits
     IF((all_orbit%jj(p2) /= all_orbit%jj(a)).or. &
          (iph(all_orbit%ll(p2)) /= iph(all_orbit%ll(a))))CYCLE
     IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE      
        DO h2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
           j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
           j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              j2min=abs((all_orbit%jj(p1)-all_orbit%jj(a))/2)
              j2max=(all_orbit%jj(p1)+all_orbit%jj(a))/2
              j_min=MAX(j1min,j2min)
              j_max=MIN(j1max,j2max)
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
              diff1=nshell2-nshell1
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(d)+nshell1
              nshell4=all_orbit%nshell(b)+all_orbit%nshell(p2)+nshell2
              diff2=nshell3-nshell4
              IF(dencheck(diff1).or.dencheck(diff2)) CYCLE
              den1=wcn+all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%evalence(a)
              den2=all_orbit%evalence(c)+all_orbit%evalence(d)- &
                   all_orbit%evalence(b)-all_orbit%e(p2)+den1
              DO pp=j_min,j_max
                 factr2=(2*pp+1)*factr1
                 CALL pphhmtx(p2,b,c,d,jt,ans1);IF ((ans1(1) == 0.)) CYCLE
                 CALL pphhmtx(p1,a,h1,h2,pp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(h1,h2,p1,p2,pp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                      all_orbit%evalence(c)+all_orbit%evalence(d)- &
                      all_orbit%evalence(a)-all_orbit%e(p1)+wcn
                 w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                 w3=w1+all_orbit%e(p1)-all_orbit%evalence(b)
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+factr2*val1*val2*val3/de(ie)   
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag9a
!
!       
!
SUBROUTINE diag9b(a,b,c,d,jt,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jt  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, h2, p1, nshell1, nshell2, diff1, &
       diff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, pp, j_min, j_max, h1, iph
  REAL(DP) :: val1, val2, val3, factr1, factr2
  LOGICAL dencheck

  diagram=0.
  factr1=-0.5/(all_orbit%jj(c)+1.)
  DO p2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p2)== 'hole' ) CYCLE
     IF((all_orbit%jj(p2) /= all_orbit%jj(c)).or. &
          (iph(all_orbit%ll(p2)) /= iph(all_orbit%ll(c))))CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
        DO h2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
           j1min=abs((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
           j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              j2min=ABS((all_orbit%jj(p1)-all_orbit%jj(c))/2)
              j2max=(all_orbit%jj(p1)+all_orbit%jj(c))/2
              j_min=MAX(j1min,j2min)
              j_max=MIN(j1max,j2max)
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              diff1=nshell2-nshell1
              nshell3=all_orbit%nshell(p2)+all_orbit%nshell(d)+nshell1
              nshell4=all_orbit%nshell(b)+all_orbit%nshell(a)+nshell2
              diff2=nshell3-nshell4
              IF(dencheck(diff1).or.dencheck(diff2)) CYCLE
              den1=wcn+all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%e(p2)+all_orbit%evalence(d)- &
                   all_orbit%evalence(a)-all_orbit%evalence(b)+den1
              DO pp=j_min,j_max
                 factr2=(2*pp+1)*factr1
                 CALL pphhmtx(a,b,p2,d,jt,ans1);IF ((ans1(1) == 0.)) CYCLE
                 CALL pphhmtx(h1,h2,p1,c,pp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(p1,p2,h1,h2,pp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                      all_orbit%evalence(d)-all_orbit%e(p1)+wcn
                 w2=w1+all_orbit%evalence(c)-all_orbit%evalence(a)- &
                      all_orbit%evalence(b)+all_orbit%e(p1)
                 w3=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+factr2*val1*val2*val3/de(ie)   
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag9b
!
!       
!
SUBROUTINE diag9f(a,b,c,d,jt,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jt  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, h2, p1, nshell1, nshell2, diff1, &
       diff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, pp, j_min, j_max, h1, iph
  REAL(DP) :: val1, val2, val3, factr1, factr2
  LOGICAL dencheck

  diagram=0.
  factr1=-0.5/(all_orbit%jj(a)+1.)
  DO p2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
     IF((all_orbit%jj(p2) /= all_orbit%jj(a)).or. &
          (iph(all_orbit%ll(p2)) /= iph(all_orbit%ll(a))))CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE   
        DO h2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
           j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
           j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              j2min=ABS((all_orbit%jj(p1)-all_orbit%jj(a))/2)
              j2max=(all_orbit%jj(p1)+all_orbit%jj(a))/2
              j_min=MAX(j1min,j2min)
              j_max=MIN(j1max,j2max)
              j2min=ABS((all_orbit%jj(p1)-all_orbit%jj(a))/2)
              j2max=(all_orbit%jj(p1)+all_orbit%jj(a))/2
              j_min=MAX(j1min,j2min)
              j_max=MIN(j1max,j2max)
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
              diff1=nshell2-nshell1
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(d)+nshell1
              nshell4=all_orbit%nshell(b)+all_orbit%nshell(p2)+nshell2
              diff2=nshell3-nshell4
              IF(dencheck(diff1).or.dencheck(diff2)) CYCLE
              den1=wcn+all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%evalence(a)
              den2=all_orbit%evalence(c)+all_orbit%evalence(d)- &
                   all_orbit%evalence(b)-all_orbit%e(p2)+den1
              DO pp=j_min,j_max
                 factr2=(2*pp+1)*factr1
                 CALL pphhmtx(p2,b,c,d,jt,ans1);IF ((ans1(1) == 0.)) CYCLE
                 CALL pphhmtx(p1,a,h1,h2,pp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(h1,h2,p1,p2,pp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                      all_orbit%evalence(c)+all_orbit%evalence(d)- &
                      all_orbit%evalence(a)-all_orbit%e(p1)+wcn
                 w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                 w3=w1+all_orbit%e(p1)-all_orbit%evalence(b)
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)-factr2*val1*val2*val3/de(ie)   
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag9f
!
!
!
SUBROUTINE diag10(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: h3, h4, h1, h2, nshell1, nshell2, nshell3, nshell4, idiff1, idiff2, ie
  REAL(DP) :: val1, val2, val3
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
        nshell2=all_orbit%nshell(a)+all_orbit%nshell(b)
        idiff1=nshell2-nshell1
        IF(dencheck(idiff1)) CYCLE
        den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%evalence(a)- &
             all_orbit%evalence(b)+wcn
        CALL pphhmtx(a,b,h1,h2,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
        DO h3=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
           DO h4=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
              nshell3=all_orbit%nshell(h3)+all_orbit%nshell(h4)
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(b)
              idiff2=nshell3-nshell4
              IF(dencheck(idiff2)) CYCLE
              den2=all_orbit%e(h3)+all_orbit%e(h4)-all_orbit%evalence(a)- &
                   all_orbit%evalence(b)+wcn
              CALL pphhmtx(h1,h2,h3,h4,jtot,ans2); IF(ans2(1)==0.) CYCLE
              CALL pphhmtx(h3,h4,c,d,jtot,ans3);IF(ans3(1)==0.) CYCLE
              w2=w1+all_orbit%e(h3)+all_orbit%e(h4)-all_orbit%evalence(a)- &
                   all_orbit%evalence(b)
              w3=all_orbit%e(h3)+all_orbit%e(h4)+all_orbit%evalence(c)+&
                   all_orbit%evalence(d)-all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
              de=den1*den2
              DO ie=1,n_startenergy_veff
                 CALL interpolate(w1(ie),e_start_g,ans1,val1)
                 CALL interpolate(w2(ie),e_start_g,ans2,val2)
                 CALL interpolate(w3(ie),e_start_g,ans3,val3)
                 diagram(ie)=diagram(ie)+0.25*val1*val2*val3/de(ie)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag10
!
!
!
SUBROUTINE diag11a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p1, p2, p3, p4, nshell1, nshell2, nshell3, idiff1, idiff2, ie
  REAL(DP) :: val1, val2, val3
  LOGICAL dencheck

  diagram=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
        nshell1=all_orbit%nshell(p1)+all_orbit%nshell(p2)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(d)
        idiff1=nshell2-nshell1
        IF(dencheck(idiff1))CYCLE
        den1=wcn+all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%e(p1)- &
             all_orbit%e(p2)
        CALL pphhmtx(a,b,p1,p2,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO p3=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p3) /= 'hole' ) CYCLE
           DO p4=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p4) /= 'hole' ) CYCLE
              nshell3=all_orbit%nshell(p3)+all_orbit%nshell(p4)
              idiff2=nshell3-nshell1
              IF(dencheck(idiff2))CYCLE
              den2=wcn+all_orbit%e(p3)+all_orbit%e(p4)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              CALL pphhmtx(p1,p2,p3,p4,jtot,ans2); IF(ans2(1)==0.) CYCLE
              CALL pphhmtx(p3,p4,c,d,jtot,ans3);IF(ans3(1)==0.) CYCLE
              w1=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
              w2=all_orbit%e(p3)+all_orbit%e(p4)+wcn
              w3=den2 
              de=den1*den2
              DO ie=1,n_startenergy_veff
                 CALL interpolate(w1(ie),e_start_g,ans1,val1)
                 CALL interpolate(w2(ie),e_start_g,ans2,val2)
                 CALL interpolate(w3(ie),e_start_g,ans3,val3)
                 diagram(ie)=diagram(ie)+0.25*val1*val2*val3/de(ie)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag11a
!
!
!
SUBROUTINE diag11b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p1, p2, h1, h2, nshell1, nshell2, nshell3, nshell4, idiff1, idiff2, ie
  REAL(DP) :: val1, val2, val3
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        CALL pphhmtx(a,b,h1,h2,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
              nshell1=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              nshell2=all_orbit%nshell(c)+all_orbit%nshell(d)
              idiff1=nshell2-nshell1
              nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)+nshell2
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(b)+nshell1
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1).or.dencheck(idiff2))CYCLE
              den1=wcn+all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%evalence(a)- &
                   all_orbit%evalence(b)+all_orbit%evalence(c)+all_orbit%evalence(d)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=wcn+all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              CALL pphhmtx(h1,h2,p1,p2,jtot,ans2); IF(ans2(1)==0.) CYCLE
              CALL pphhmtx(p1,p2,c,d,jtot,ans3);IF(ans3(1)==0.) CYCLE
              w1=all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%e(p2)- &
                   all_orbit%e(p1)+all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%evalence(a)- &
                   all_orbit%evalence(b)+all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w3=all_orbit%evalence(d)+all_orbit%evalence(c)+wcn
              de=den1*den2
              DO ie=1,n_startenergy_veff
                 CALL interpolate(w1(ie),e_start_g,ans1,val1)
                 CALL interpolate(w2(ie),e_start_g,ans2,val2)
                 CALL interpolate(w3(ie),e_start_g,ans3,val3)
                 diagram(ie)=diagram(ie)+0.25*val1*val2*val3/de(ie)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag11b
!         
!
!
SUBROUTINE diag11fa(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p1, p2, p3, p4, nshell1, nshell2, nshell3, idiff1, idiff2, ie
  REAL(DP) :: val1, val2, val3
  LOGICAL dencheck

  diagram=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
        nshell1=all_orbit%nshell(p1)+all_orbit%nshell(p2)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(d)
        idiff1=nshell2-nshell1
        IF(dencheck(idiff1))CYCLE
        den1=wcn+all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%e(p1)- &
             all_orbit%e(p2)
        CALL pphhmtx(a,b,p1,p2,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO p3=1, all_orbit%total_orbits
           IF(all_orbit%model_space(p3) /= 'inside' ) CYCLE      
           DO p4=1, all_orbit%total_orbits
              IF(all_orbit%model_space(p4) /= 'inside' ) CYCLE      
              nshell3=all_orbit%nshell(p3)+all_orbit%nshell(p4)
              idiff2=nshell3-nshell1
              IF(dencheck(idiff2))CYCLE
              den2=wcn+all_orbit%evalence(p3)+all_orbit%evalence(p4)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              CALL pphhmtx(p1,p2,p3,p4,jtot,ans2); IF(ans2(1)==0.) CYCLE
              CALL pphhmtx(p3,p4,c,d,jtot,ans3);IF(ans3(1)==0.) CYCLE
              w1=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
              w2=all_orbit%evalence(p3)+all_orbit%evalence(p4)+wcn
              w3=den2
              de=den1*den2
              DO ie=1,n_startenergy_veff
                 CALL interpolate(w1(ie),e_start_g,ans1,val1)
                 CALL interpolate(w2(ie),e_start_g,ans2,val2)
                 CALL interpolate(w3(ie),e_start_g,ans3,val3)
                 diagram(ie)=diagram(ie)-0.25*val1*val2*val3/de(ie)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag11fa
!
!
!
SUBROUTINE diag11fb(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p1, p2, h1, h2, nshell1, nshell2, nshell3, nshell4, idiff1, idiff2, ie
  REAL(DP) :: val1, val2, val3
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
        nshell2=all_orbit%nshell(a)+all_orbit%nshell(b)
        idiff1=nshell2-nshell1
        IF(dencheck(idiff1)) CYCLE
        den1=wcn+all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%evalence(a)- &
             all_orbit%evalence(b)
        CALL pphhmtx(a,b,h1,h2,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%model_space(p1) /= 'inside' ) CYCLE      
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%model_space(p2) /= 'inside' ) CYCLE
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(d)+nshell1
              nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)+nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff2))CYCLE
              den2=wcn+den1+all_orbit%evalence(c)+all_orbit%evalence(d)- &
                   all_orbit%evalence(p1)-all_orbit%evalence(p2)
              CALL pphhmtx(h1,h2,p1,p2,jtot,ans2); IF(ans2(1)==0.) CYCLE
              CALL pphhmtx(p1,p2,c,d,jtot,ans3);IF(ans3(1)==0.) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1+all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%evalence(a)- &
                   all_orbit%evalence(b)
              w3=w2
              de=den1*den2
              DO ie=1,n_startenergy_veff
                 CALL interpolate(w1(ie),e_start_g,ans1,val1)
                 CALL interpolate(w2(ie),e_start_g,ans2,val2)
                 CALL interpolate(w3(ie),e_start_g,ans3,val3)
                 diagram(ie)=diagram(ie)-0.25*val1*val2*val3/de(ie)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag11fb
!
!
!
SUBROUTINE diag12a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: phmin, phmx, jph, p1, &
       p2, p3, h1, iph, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, j2max
  REAL(DP) :: val1, val2, val3, sixj, sg
  LOGICAL dencheck

  diagram=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(p1)-all_orbit%jj(c))/2)
        j2min=ABS((all_orbit%jj(p2)-all_orbit%jj(d))/2)
        j1max=(all_orbit%jj(p1)+all_orbit%jj(c))/2
        j2max=(all_orbit%jj(p2)+all_orbit%jj(d))/2
        phmin=MAX(j1min,j2min)
        phmx=MIN(j1max,j2max)
        nshell1=all_orbit%nshell(p1)+all_orbit%nshell(p2)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(d)
        idiff1=nshell2-nshell1
        IF(dencheck(idiff1)) CYCLE
        den2=all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%e(p1)-all_orbit%e(p2)+wcn
        CALL pphhmtx(a,b,p1,p2,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO jph=phmin,phmx
           sixj=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
                all_orbit%jj(p2),all_orbit%jj(p1),2*jph)
           DO p3=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE      
              DO h1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p3)
                 nshell4=all_orbit%nshell(c)+all_orbit%nshell(h1)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2))CYCLE
                 sg=iph((all_orbit%jj(p1)+all_orbit%jj(d)+ &
                      all_orbit%jj(p3)+all_orbit%jj(h1))/2+jtot &
                      +all_orbit%jj(h1))
                 den1=all_orbit%evalence(c)+all_orbit%e(h1)-all_orbit%e(p1) &
                      -all_orbit%e(p3)+wcn
                 CALL cross_coupled_mtxel1(p1,p3,c,h1,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(h1,p2,p3,d,jph,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
                 w2=all_orbit%evalence(c)+all_orbit%e(h1)+wcn
                 w3=w1+all_orbit%e(h1)-all_orbit%e(p1)
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+sixj*sg*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag12a
!
!
!
SUBROUTINE diag12b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: phmin, phmx, jph, p1, &
       p2, p3, h1, iph, nshell1, nshell2,  idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, j2max
  REAL(DP) :: val1, val2, val3, sixj, sg
  LOGICAL dencheck

  diagram=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(p1)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(p2)-all_orbit%jj(b))/2)
        j1max=(all_orbit%jj(p1)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(p2)+all_orbit%jj(b))/2
        phmin=MIN(j1min,j2min)
        phmx=MAX(j1max,j2max)
        nshell1=all_orbit%nshell(p1)+all_orbit%nshell(p2)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(d)
        idiff1=nshell2-nshell1
        IF(dencheck(idiff1)) CYCLE
        den1=all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%e(p1)-all_orbit%e(p2)+wcn
        CALL pphhmtx(p1,p2,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO jph=phmin,phmx
           sixj=sjs(all_orbit%jj(p1),all_orbit%jj(p2),2*jtot, &
                all_orbit%jj(b),all_orbit%jj(a),2*jph)
           DO p3=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE      
              DO h1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p3)+ &
                      all_orbit%nshell(b)
                 nshell4=all_orbit%nshell(h1)+nshell2
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=all_orbit%evalence(c)+all_orbit%evalence(d)+all_orbit%e(h1)- &
                      all_orbit%e(p3)-all_orbit%e(p1)-all_orbit%evalence(b)+wcn
                 sg=iph((all_orbit%jj(p3)+all_orbit%jj(p2)+ &
                      all_orbit%jj(a)+all_orbit%jj(h1))/2+jtot+ &
                      all_orbit%jj(h1))
                 CALL cross_coupled_mtxel1(a,h1,p1,p3,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(p3,b,h1,p2,jph,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
                 w2=w1+all_orbit%e(h1)-all_orbit%evalence(b)
                 w3=w1+all_orbit%e(h1)-all_orbit%e(p1)
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+sixj*sg*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag12b
!
!
!
SUBROUTINE diag13a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: phmin, phmx, jph, p1, &
       p2, h2, h1, iph, nshell1, nshell2,  idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, j2max
  REAL(DP) :: val1, val2, val3, sixj, sg
  LOGICAL dencheck

  diagram=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(p1)-all_orbit%jj(b))/2)
        j1max=(all_orbit%jj(h1)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(p1)+all_orbit%jj(b))/2
        phmin=MAX(j1min,j2min)
        phmx=MIN(j1max,j2max)
        CALL pphhmtx(h1,p1,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO jph=phmin,phmx
           sixj=sjs(all_orbit%jj(h1),all_orbit%jj(p1), &
                2*jtot,all_orbit%jj(b),all_orbit%jj(a),2*jph)
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
              DO h2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                 nshell2=all_orbit%nshell(a)+all_orbit%nshell(p2)
                 idiff1=nshell2-nshell1
                 nshell3=all_orbit%nshell(h2)+all_orbit%nshell(d)+ &
                      all_orbit%nshell(c)
                 nshell4=all_orbit%nshell(p1)+nshell2
                 idiff2=nshell3-nshell4
                 IF(dencheck(idiff1).or.dencheck(idiff2))CYCLE
                 den1=all_orbit%e(h2)+all_orbit%e(h1)-all_orbit%evalence(a)- &
                      all_orbit%e(p2)+wcn
                 den2=all_orbit%evalence(c)+all_orbit%evalence(d)+all_orbit%e(h2)- &
                      all_orbit%evalence(a)-all_orbit%e(p1)-all_orbit%e(p2)+wcn
                 CALL cross_coupled_mtxel1(p2,a,h2,h1,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(h2,b,p2,p1,jph,ans3);IF(ans3(1)==0.) CYCLE
                 sg=iph((all_orbit%jj(p1)+all_orbit%jj(p2)+all_orbit%jj(a) &
                      +all_orbit%jj(h2))/2+jtot+all_orbit%jj(h2))
                 w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)+ &
                      all_orbit%evalence(d)-all_orbit%evalence(a)-all_orbit%e(p2)+wcn
                 w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                 w3=all_orbit%e(h2)+all_orbit%evalence(c)+all_orbit%evalence(d)- &
                      all_orbit%evalence(a)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)-sixj*sg*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag13a
!
!
!
SUBROUTINE diag13b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: phmin, phmx, jph, p1, &
       p2, h2, h1, iph, nshell1, nshell2,  idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, j2max
  REAL(DP) :: val1, val2, val3,  sixj, sg
  LOGICAL dencheck

  diagram=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(c))/2)
        j2min=ABS((all_orbit%jj(p1)-all_orbit%jj(d))/2)
        j1max=(all_orbit%jj(h1)+all_orbit%jj(c))/2
        j2max=(all_orbit%jj(p1)+all_orbit%jj(d))/2
        phmin=MAX(j1min,j2min)
        phmx=MIN(j1max,j2max)
        CALL pphhmtx(a,b,h1,p1,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO jph=phmin,phmx      
           sixj=sjs(all_orbit%jj(h1),all_orbit%jj(p1),2*jtot, &
                all_orbit%jj(d),all_orbit%jj(c),2*jph)
           IF(sixj == 0.) CYCLE
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
              DO h2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 nshell1=all_orbit%nshell(d)+all_orbit%nshell(h2)
                 nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff1=nshell2-nshell1
                 nshell3=all_orbit%nshell(h2)+all_orbit%nshell(d)+ &
                      all_orbit%nshell(h1)
                 nshell4=all_orbit%nshell(p2)+all_orbit%nshell(a)+ &
                      all_orbit%nshell(b)
                 idiff2=nshell3-nshell4
                 IF(dencheck(idiff1).or.dencheck(idiff2))CYCLE
                 den1=all_orbit%e(h2)+all_orbit%evalence(d)-all_orbit%e(p1)- &
                      all_orbit%e(p2)+wcn
                 den2=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(d)- &
                      all_orbit%evalence(a)-all_orbit%evalence(b)-all_orbit%e(p2)+wcn
                 CALL cross_coupled_mtxel1(h2,h1,p2,c,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(p2,p1,h2,d,jph,ans3);IF(ans3(1)==0.) CYCLE
                 sg=iph((all_orbit%jj(h1)+all_orbit%jj(p2)+ &
                      all_orbit%jj(d)+all_orbit%jj(h2))/2+jtot+ &
                      all_orbit%jj(h2))
                 w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(d)- &
                      all_orbit%e(p2)+wcn
                 w2=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(d)+ &
                      all_orbit%evalence(c)-all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
                 w3=all_orbit%e(h2)+all_orbit%evalence(d)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)-sixj*sg*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag13b
!         
!
!
SUBROUTINE diag13f(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: phmin, phmx, jph, p1, &
       p2, h2, h1, nshell1, nshell2,  idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, j2max
  REAL(DP) :: val1, val2, val3, sixj
  LOGICAL dencheck

  diagram=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%model_space(h1) /= 'inside' ) CYCLE
        j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(p1)-all_orbit%jj(b))/2)
        j1max=(all_orbit%jj(h1)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(p1)+all_orbit%jj(b))/2
        phmin=MAX(j1min,j2min)
        phmx=MIN(j1max,j2max)
        CALL pphhmtx(h1,p1,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO jph=phmin,phmx
           sixj=sjs(all_orbit%jj(a),all_orbit%jj(b),2*jtot, &
                all_orbit%jj(p1),all_orbit%jj(h1),2*jph)
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE  
              DO h2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                 nshell2=all_orbit%nshell(a)+all_orbit%nshell(p2)
                 idiff1=nshell2-nshell1
                 nshell3=all_orbit%nshell(h2)+all_orbit%nshell(d)+ &
                      all_orbit%nshell(c)
                 nshell4=all_orbit%nshell(p1)+nshell2
                 idiff2=nshell3-nshell4
                 IF(dencheck(idiff1).or.dencheck(idiff2))CYCLE
                 den1=all_orbit%e(h2)+all_orbit%evalence(h1)-all_orbit%evalence(a)- &
                      all_orbit%e(p2)+wcn
                 den2=all_orbit%evalence(c)+all_orbit%e(h2)+all_orbit%evalence(d)- &
                      all_orbit%evalence(a)-all_orbit%e(p1)-all_orbit%e(p2)+wcn
                 CALL cross_coupled_mtxel1(a,p2,h1,h2,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(p1,p2,b,h2,jph,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%evalence(h1)+all_orbit%e(h2)+all_orbit%evalence(d)+ &
                      all_orbit%evalence(c)-all_orbit%evalence(a)-all_orbit%e(p2)+wcn
                 w2=all_orbit%evalence(h1)+all_orbit%e(h2)+wcn
                 w3=all_orbit%e(h2)+all_orbit%evalence(d)+all_orbit%evalence(c)- &
                      all_orbit%evalence(a)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)-sixj*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag13f
!
!
!
SUBROUTINE diag14a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: phmin, phmx, jph, p1, &
       h2, h3, h1, iph, nshell1, nshell2,  idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, j2max
  REAL(DP) :: val1, val2, val3, sixj, sg
  LOGICAL dencheck

  diagram=0.
  DO h2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE         
     DO h3=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
        nshell1=all_orbit%nshell(h3)+all_orbit%nshell(h2)
        nshell2=all_orbit%nshell(a)+all_orbit%nshell(b)
        idiff1=nshell2-nshell1
        IF(dencheck(idiff1)) CYCLE
        den1=all_orbit%e(h3)+all_orbit%e(h2)-all_orbit%evalence(b)-all_orbit%evalence(a)+wcn
        j1min=ABS((all_orbit%jj(h3)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(h2)-all_orbit%jj(b))/2)
        j1max=(all_orbit%jj(h3)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(h2)+all_orbit%jj(b))/2
        phmin=MAX(j1min,j2min)
        phmx=MIN(j1max,j2max)
        CALL pphhmtx(h2,h3,c,d,jtot,ans1); IF ((ans1(1) == 0.)) CYCLE
        DO jph=phmin,phmx
           sixj=sjs(all_orbit%jj(h3),all_orbit%jj(h2),2*jtot,all_orbit%jj(b),all_orbit%jj(a),2*jph)
           IF(sixj == 0.) CYCLE
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
              DO h1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h3)+all_orbit%nshell(h1)
                 nshell4=all_orbit%nshell(a)+all_orbit%nshell(p1)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=all_orbit%e(h3)+all_orbit%e(h1)-all_orbit%e(p1)- &
                      all_orbit%evalence(a)+wcn
                 CALL cross_coupled_mtxel1(a,p1,h3,h1,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(h1,b,p1,h2,jph,ans3); IF(ans3(1)==0.) CYCLE
                 sg=iph((all_orbit%jj(h1)+all_orbit%jj(p1)+all_orbit%jj(a)+ &
                      all_orbit%jj(h3))/2+all_orbit%jj(h2)+all_orbit%jj(h1))
                 w1=all_orbit%e(h3)+all_orbit%e(h2)+all_orbit%evalence(d)+ &
                      all_orbit%evalence(c)-all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
                 w2=all_orbit%e(h1)+all_orbit%e(h3)+wcn
                 w3=all_orbit%e(h2)+all_orbit%e(h1)+all_orbit%e(h3)- &
                      all_orbit%evalence(a)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+sixj*sg*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag14a
!
!
!
SUBROUTINE diag14b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: phmin, phmx, jph, p1, &
       h2, h3, h1, iph, nshell1, nshell2,  idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, j2max
  REAL(DP) :: val1, val2, val3, sixj, sg
  LOGICAL dencheck

  diagram=0.
  DO h2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE   
     DO h3=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
        nshell1=all_orbit%nshell(h3)+all_orbit%nshell(h2)
        nshell2=all_orbit%nshell(a)+all_orbit%nshell(b)
        idiff1=nshell2-nshell1
        IF(dencheck(idiff1)) CYCLE
        den1=all_orbit%e(h3)+all_orbit%e(h2)-all_orbit%evalence(b)-all_orbit%evalence(a)+wcn
        j1min=ABS((all_orbit%jj(h3)-all_orbit%jj(c))/2)
        j2min=ABS((all_orbit%jj(h2)-all_orbit%jj(d))/2)
        j1max=(all_orbit%jj(h3)+all_orbit%jj(c))/2
        j2max=(all_orbit%jj(h2)+all_orbit%jj(d))/2
        phmin=MAX(j1min,j2min)
        phmx=MIN(j1max,j2max)
        CALL pphhmtx(a,b,h2,h3,jtot,ans1); IF ((ans1(1) == 0.)) CYCLE
        DO jph=phmin,phmx
           sixj=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot,all_orbit%jj(h2),all_orbit%jj(h3),2*jph)
           IF(sixj == 0.) CYCLE
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
              DO h1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h3)+all_orbit%nshell(h1)+ &
                      all_orbit%nshell(d)
                 nshell4=all_orbit%nshell(a)+all_orbit%nshell(b)+ &
                      all_orbit%nshell(p1)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=all_orbit%e(h3)+all_orbit%e(h1)+all_orbit%evalence(d)- &
                      all_orbit%evalence(b)-all_orbit%evalence(a)-all_orbit%e(p1)+wcn
                 CALL cross_coupled_mtxel1(h2,h1,d,p1,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(h3,p1,c,h1,jph,ans3); IF(ans3(1)==0.) CYCLE
                 sg=iph((all_orbit%jj(h1)+all_orbit%jj(p1)+ &
                      all_orbit%jj(d)+all_orbit%jj(h2))/2)
                 w1=all_orbit%e(h2)+all_orbit%e(h3)+wcn
                 w2=all_orbit%e(h1)+all_orbit%e(h3)+all_orbit%evalence(d)+ &
                      all_orbit%evalence(c)-all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
                 w3=all_orbit%e(h2)+all_orbit%e(h1)+all_orbit%e(h3)+ &
                      all_orbit%evalence(d)-all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+sixj*sg*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag14b
!
!
!
SUBROUTINE diag15(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p1, p2, h1, h2, j1min, j2min, j2max, j1max, phmin, &
       phmx, nshell1, nshell2, nshell3, nshell4, idiff1, &
       idiff2, isgcd, jph, iph, ie
  REAL(DP) :: val1, val2, val3, sixj, sg
  LOGICAL dencheck

  diagram=0.
  J1MIN=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  J2MIN=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  J1MAX=(all_orbit%jj(a)+all_orbit%jj(c))/2
  J2MAX=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(j1min,j2min)
  phmx=MIN(j1max,j2max)
  DO jph=phmin,phmx
     sixj=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     IF(sixj == 0.) CYCLE
     sixj=sixj/SQRT(2.*jph+1.)
     DO p1=1,all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h1=1,all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell3=all_orbit%nshell(h1)+all_orbit%nshell(d)
           nshell4=all_orbit%nshell(p1)+all_orbit%nshell(b)
           idiff2=nshell4-nshell3
           IF(dencheck(idiff2)) CYCLE
           CALL cross_coupled_mtxel1(a,h1,c,p1,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           den1=wcn+all_orbit%e(h1)+all_orbit%evalence(d)- &
                all_orbit%evalence(b)-all_orbit%e(p1)
           DO p2=1,all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' )CYCLE
              DO h2=1,all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 nshell1=all_orbit%nshell(d)+all_orbit%nshell(h2)
                 nshell2=all_orbit%nshell(p2)+all_orbit%nshell(b)
                 idiff1=nshell2-nshell1
                 IF(dencheck(idiff1 ) ) CYCLE
                 isgcd=(all_orbit%jj(d)+all_orbit%jj(a)+ &
                      all_orbit%jj(p1)+all_orbit%jj(p2)+ &
                      all_orbit%jj(h1)+all_orbit%jj(h2))/2
                 sg=iph(isgcd+all_orbit%jj(h1)+all_orbit%jj(h2)&
                      +jtot+jph)
                 den2=wcn+all_orbit%evalence(d)+all_orbit%e(h2)- &
                      all_orbit%evalence(b)-all_orbit%e(p2)
                 CALL cross_coupled_mtxel1(p1,h2,h1,p2,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(p2,b,h2,d,jph,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%e(h1)+all_orbit%evalence(d)+ &
                      all_orbit%evalence(c)-all_orbit%evalence(b)+wcn
                 w2=all_orbit%e(h1)+all_orbit%e(h2)+ &
                      all_orbit%evalence(d)-all_orbit%evalence(b)+wcn
                 w3=all_orbit%e(h2)+all_orbit%evalence(d)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+sg*sixj*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag15
!
!
!
SUBROUTINE diag16a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p1, p2, h1, h2, j1min, j2min, j2max, j1max, phmin, &
       phmx, nshell1, nshell2, nshell3, nshell4, idiff1, &
       idiff2, isgcd, jph, iph, ie
  REAL(DP) :: val1, val2, val3, sixj, sg  
  LOGICAL dencheck

  diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j2min=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2max=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(j1min,j2min)
  phmx=MIN(j1max,j2max)
  if(phmx.lt.phmin) return
  DO jph=phmin,phmx
     sixj=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     sixj=sixj/sqrt(2.*jph+1.)
     DO p1=1,all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
        DO h1=1,all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell1=all_orbit%nshell(d)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(b)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1)) CYCLE
           CALL cross_coupled_mtxel1(a,h1,c,p1,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           den1=wcn+all_orbit%e(h1)+all_orbit%evalence(d)- &
                all_orbit%e(p1)-all_orbit%evalence(b)
           DO p2=1,all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' )CYCLE 
              DO h2=1,all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                 nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(h1)+all_orbit%e(h2)- &
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL cross_coupled_mtxel1(p1,p2,h1,h2,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(h2,b,p2,d,jph,ans3);IF(ans3(1)==0.) CYCLE
                 isgcd=(all_orbit%jj(d)+all_orbit%jj(a)+ &
                      all_orbit%jj(p1)+all_orbit%jj(p2)+ &
                      all_orbit%jj(h1)+all_orbit%jj(h2))/2
                 sg=iph(isgcd+all_orbit%jj(h1)+ &
                      all_orbit%jj(h2)+jtot+jph)
                 w1=all_orbit%e(h1)+all_orbit%evalence(d)+ &
                      all_orbit%evalence(c)-all_orbit%evalence(b)+wcn
                 w2=all_orbit%e(h2)+all_orbit%e(h1)+wcn
                 w3=all_orbit%e(h1)+all_orbit%e(h2)+ &
                      all_orbit%evalence(d)-all_orbit%e(p1)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+sg*sixj*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag16a
!
!
!
SUBROUTINE diag16b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p1, p2, h1, h2, j1min, j2min, j2max, j1max, phmin, &
       phmx, nshell1, nshell2, nshell3, nshell4, idiff1, &
       idiff2, isgcd, jph, iph, ie
  REAL(DP) :: val1, val2, val3, sixj, sg
  LOGICAL dencheck

  diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j2min=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2max=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(j1min,j2min)
  phmx=MIN(j1max,j2max)
  DO jph=phmin,phmx
     sixj=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     IF (sixj == 0.) CYCLE
     sixj=sixj/sqrt(2.*jph+1.)
     DO p1=1,all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE     
        DO h1=1,all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell1=all_orbit%nshell(c)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
           idiff1=nshell2-nshell1
           IF (dencheck(idiff1)) CYCLE
           den1=wcn+all_orbit%e(h1)+all_orbit%evalence(c)- &
                all_orbit%e(p1)-all_orbit%evalence(a)
           CALL cross_coupled_mtxel1(a,p1,c,h1,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           DO p2=1,all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' )CYCLE    
              DO h2=1,all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' )CYCLE
                 nshell3=all_orbit%nshell(h2)+all_orbit%nshell(d)+nshell1
                 nshell4=all_orbit%nshell(p2)+all_orbit%nshell(b)+nshell2
                 idiff2=nshell4-nshell3
                 IF (dencheck(idiff2)) CYCLE
                 den2=den1+all_orbit%e(h2)+all_orbit%evalence(d) &
                      -all_orbit%e(p2)-all_orbit%evalence(b)
                 CALL cross_coupled_mtxel1(h1,h2,p1,p2,jph,ans2); IF(ans2(1)==0.) CYCLE
                 CALL cross_coupled_mtxel1(p2,b,h2,d,jph,ans3);IF(ans3(1)==0.) CYCLE
                 isgcd=(all_orbit%jj(d)+all_orbit%jj(a)+ &
                      all_orbit%jj(p1)+all_orbit%jj(p2)+ &
                      all_orbit%jj(h1)+all_orbit%jj(h2))/2
                 sg=iph(isgcd+all_orbit%jj(h1)+all_orbit%jj(h2)&
                      +jtot+jph)
                 w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                 w2=all_orbit%e(h1)+all_orbit%e(h2)+ &
                      all_orbit%evalence(d)+all_orbit%evalence(c)- &
                      all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
                 w3=all_orbit%e(h2)+all_orbit%e(h1)+ &
                      all_orbit%evalence(c)+all_orbit%evalence(d)-all_orbit%evalence(a)&
                      -all_orbit%e(p1)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+sg*sixj*val1*val2*val3/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag16b
!
!
!      
SUBROUTINE diag17a(a,b,c,d,jtot,diagram)
  USE constants
  USE ang_mom_functions
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: acmin, bdmin, acmax, bdmax, ppmin, ppmax, jph, p1, &
       p2, h2, h1, iph, nshell1, nshell2, isgcd, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, phmin, phmx
  REAL(DP) :: val1, val2, val3, factr, sixj1, sixj2, sg
  LOGICAL dencheck

  diagram=0.
  acmin=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  bdmin=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  acmax=(all_orbit%jj(a)+all_orbit%jj(c))/2
  bdmax=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(acmin,bdmin)
  phmx=MIN(acMAX,bdMAX)
  DO jph=phmin,phmx
     sixj1=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     IF(sixj1 == 0.) CYCLE
     sixj1=sixj1*sqrt(2.*jph+1.)
     DO p1=1,all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE            
        DO h1=1,all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell1=all_orbit%nshell(c)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1)) CYCLE
           den1=wcn+all_orbit%e(h1)+all_orbit%evalence(c)- &
                all_orbit%e(p1)-all_orbit%evalence(a)
           j1min=ABS((all_orbit%jj(a)-all_orbit%jj(h1))/2)
           j2min=ABS((all_orbit%jj(c)-all_orbit%jj(p1))/2)
           j1max=(all_orbit%jj(a)+all_orbit%jj(h1))/2
           j2max=(all_orbit%jj(c)+all_orbit%jj(p1))/2
           ppmin=MAX(j1min,j2min)
           ppmax=MIN(j1MAX,j2MAX)
           IF(ppmax.lt.ppmin) CYCLE
           DO jpp=ppmin,ppmax
              sixj2=sjs(all_orbit%jj(c),all_orbit%jj(a),&
                   2*jph,all_orbit%jj(h1),all_orbit%jj(p1),2*jpp)
              IF(sixj2 == 0.) CYCLE
              DO p2=1,all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1,all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(c)
                    nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                    idiff2=nshell4-nshell3
                    IF(dencheck(idiff2)) CYCLE
                    den2=wcn+all_orbit%e(h2)+all_orbit%evalence(c)-&
                         all_orbit%e(p2)-all_orbit%e(p1)
                    CALL cross_coupled_mtxel1(h1,b,p1,d,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
                    CALL cross_coupled_mtxel2(a,h2,p2,h1,jpp,ans2); IF(ans2(1)==0.) CYCLE
                    CALL cross_coupled_mtxel2(p2,p1,c,h2,jpp,ans3);IF(ans3(1)==0.) CYCLE
                    isgcd=(all_orbit%jj(a)+all_orbit%jj(d)+ &
                         all_orbit%jj(p1)+all_orbit%jj(p2)+ &
                         all_orbit%jj(h1)+all_orbit%jj(h2))/2
                    sg=iph(isgcd+all_orbit%jj(h1)+all_orbit%jj(h2)+ &
                         all_orbit%jj(a)+jtot-jpp)
                    factr=sixj1*sixj2*sg
                    w1=all_orbit%e(h1)+all_orbit%evalence(d)+all_orbit%evalence(c)-&
                         all_orbit%evalence(a)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)-&
                         all_orbit%e(p1)+wcn
                    w3=all_orbit%e(h2)+all_orbit%evalence(c)+wcn
                    de=den1*den2
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)-factr*val1*val2*val3/de(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag17a
!
!
!
SUBROUTINE diag17b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: acmin, bdmin, acmax, bdmax, ppmin, ppmax, jph, p1, &
       p2, h2, h1, iph, nshell1, nshell2, isgcd, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, phmin, phmx
  REAL(DP) :: val1, val2, val3, factr, sixj1, sixj2, sg
  LOGICAL dencheck

  diagram=0.
  acmin=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  bdmin=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  acmax=(all_orbit%jj(a)+all_orbit%jj(c))/2
  bdmax=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(acmin,bdmin)
  phmx=MIN(acMAX,bdMAX)
  DO jph=phmin,phmx
     sixj1=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot,&
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     IF(sixj1 == 0.) CYCLE
     sixj1=sixj1*sqrt(2.*jph+1.)
     DO p1=1,all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h1=1,all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell1=all_orbit%nshell(d)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(b)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1)) CYCLE
           den1=wcn+all_orbit%e(h1)+all_orbit%evalence(d)- &
                all_orbit%e(p1)-all_orbit%evalence(b)
           CALL cross_coupled_mtxel1(p1,b,h1,d,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
           j2min=ABS((all_orbit%jj(c)-all_orbit%jj(h1))/2)
           j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
           j2max=(all_orbit%jj(c)+all_orbit%jj(h1))/2
           ppmin=MAX(j1min,j2min)
           ppmax=MIN(j1max,j2max)
           DO jpp=ppmin,ppMAX
              sixj2=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                   2*jph,all_orbit%jj(p1),all_orbit%jj(h1),2*jpp)
              IF(sixj2 == 0.) CYCLE
              DO p2=1,all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE      
                 DO h2=1,all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(d)+ &
                         all_orbit%nshell(c)
                    nshell4=all_orbit%nshell(p1)+all_orbit%nshell(b)+ &
                         all_orbit%nshell(p2)
                    idiff2=nshell4-nshell3
                    IF(dencheck(idiff2)) CYCLE
                    den2=wcn+all_orbit%e(h2)+all_orbit%evalence(c)+all_orbit%evalence(d)- &
                         all_orbit%e(p2)-all_orbit%e(p1)-all_orbit%evalence(b)
                    CALL cross_coupled_mtxel2(a,h2,p2,p1,jpp,ans2); IF(ans2(1)==0.) CYCLE
                    CALL cross_coupled_mtxel2(p2,h1,c,h2,jpp,ans3);IF(ans3(1)==0.) CYCLE
                    isgcd=(all_orbit%jj(a)+all_orbit%jj(h1)+all_orbit%jj(d)+ &
                         all_orbit%jj(p1)+all_orbit%jj(p2)+all_orbit%jj(h2))/2
                    sg=iph(isgcd+all_orbit%jj(h2)+2*all_orbit%jj(h1)+ &
                         all_orbit%jj(a)+all_orbit%jj(p1)+jtot-jpp)
                    factr=sixj1*sixj2*sg
                    w1=all_orbit%e(h1)+all_orbit%evalence(d)+wcn
                    w2=all_orbit%evalence(c)+all_orbit%e(h2)+all_orbit%evalence(d)- &
                         all_orbit%evalence(b)+wcn
                    w3=w2+all_orbit%e(h1)-all_orbit%e(p1)
                    de=den1*den2
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)-factr*val1*val2*val3/de(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag17b
!
!
!
SUBROUTINE diag18a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: acmin, bdmin, acmax, bdmax, ppmin, ppmax, jph, p1, &
       p2, h2, h1, iph, nshell1, nshell2, isgcd, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, phmin, phmx
  REAL(DP) :: val1, val2, val3, factr, sixj1, sixj2, sg
  LOGICAL dencheck

  diagram=0.
  acmin=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  bdmin=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  acmax=(all_orbit%jj(a)+all_orbit%jj(c))/2
  bdmax=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(acmin,bdmin)
  phmx=MIN(acmax,bdmax)
  DO jph=phmin,phmx
     sixj1=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     sixj1=sixj1*sqrt(2.*jph+1.)
     DO p1=1,all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
        DO h1=1,all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell1=all_orbit%nshell(c)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1)) CYCLE
           den1=wcn+all_orbit%e(h1)+all_orbit%evalence(c)-all_orbit%e(p1)-all_orbit%evalence(a)
           j1min=ABS((all_orbit%jj(a)-all_orbit%jj(h1))/2)
           j2min=ABS((all_orbit%jj(c)-all_orbit%jj(p1))/2)
           j1max=(all_orbit%jj(a)+all_orbit%jj(h1))/2
           j2max=(all_orbit%jj(c)+all_orbit%jj(p1))/2
           ppmin=MAX(j1min,j2min)
           ppmax=MIN(j1max,j2max)
           CALL cross_coupled_mtxel1(h1,b,p1,d,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           DO jpp=ppmin,ppmax
              sixj2=sjs(all_orbit%jj(c),all_orbit%jj(a),2*jph, &
                   all_orbit%jj(h1),all_orbit%jj(p1),2*jpp)
              IF(sixj2 == 0.) CYCLE
              DO p2=1,all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1,all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(p2)+all_orbit%nshell(a)
                    idiff2=nshell4-nshell3
                    IF(dencheck(idiff2)) CYCLE
                    den2=wcn+all_orbit%e(h2)+all_orbit%e(h1)-all_orbit%e(p2)-&
                         all_orbit%evalence(a)
                    CALL cross_coupled_mtxel2(a,p2,h2,h1,jpp,ans2); IF(ans2(1)==0.) CYCLE
                    CALL cross_coupled_mtxel2(h2,p1,c,p2,jpp,ans3);IF(ans3(1)==0.) CYCLE
                    isgcd=(all_orbit%jj(d)+all_orbit%jj(a)+all_orbit%jj(p2)+&
                         all_orbit%jj(h1)+all_orbit%jj(p1)+all_orbit%jj(h2))/2
                    sg=iph(isgcd+all_orbit%jj(h1)+all_orbit%jj(a)+ &
                         all_orbit%jj(p2)+jtot-jpp)
                    factr=sixj1*sixj2*sg
                    w1=all_orbit%e(h1)+all_orbit%evalence(d)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w3=w2+all_orbit%evalence(c)-all_orbit%evalence(a)
                    de=den1*den2
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)-factr*val1*val2*val3/de(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag18a
!
!
!
SUBROUTINE diag18b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: acmin, bdmin, acmax, bdmax, ppmin, ppmax, jph, p1, &
       p2, h2, h1, iph, nshell1, nshell2, isgcd, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, phmin, phmx
  REAL(DP) :: val1, val2, val3, factr, sixj1, sixj2, sg
  LOGICAL dencheck

  diagram=0.
  acmin=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  bdmin=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  acmax=(all_orbit%jj(a)+all_orbit%jj(c))/2
  bdmax=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(acmin,bdmin)
  phmx=MIN(acmax,bdmax)
  DO jph=phmin,phmx
     sixj1=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     IF(sixj1 == 0.) CYCLE
     sixj1=sixj1*sqrt(2.*jph+1.)
     DO p1=1,all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE         
        DO h1=1,all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell1=all_orbit%nshell(d)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(b)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1)) CYCLE
           den1=wcn+all_orbit%e(h1)+all_orbit%evalence(d)- &
                all_orbit%e(p1)-all_orbit%evalence(b)
           j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
           j2min=ABS((all_orbit%jj(c)-all_orbit%jj(h1))/2)
           j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
           j2max=(all_orbit%jj(c)+all_orbit%jj(h1))/2
           ppmin=MAX(j1min,j2min)
           ppmax=MIN(j1max,j2max)
           CALL cross_coupled_mtxel1(p1,b,h1,d,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           DO jpp=ppmin,ppmax
              sixj2=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                   2*jph,all_orbit%jj(p1),all_orbit%jj(h1),2*jpp)
              IF(sixj2 == 0.) CYCLE
              DO p2=1,all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE      
                 DO h2=1,all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h1)+all_orbit%nshell(d)+&
                         all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(p2)+all_orbit%nshell(b)+ &
                         all_orbit%nshell(a)
                    idiff2=nshell4-nshell3
                    IF(dencheck(idiff2)) CYCLE
                    den2=wcn+all_orbit%e(h2)+all_orbit%e(h1)+all_orbit%evalence(d)- &
                         all_orbit%e(p2)-all_orbit%evalence(a)-all_orbit%evalence(b)

                    CALL cross_coupled_mtxel2(h2,h1,c,p2,jpp,ans2); IF(ans2(1)==0.) CYCLE
                    CALL cross_coupled_mtxel2(a,p2,h2,p1,jpp,ans3);IF(ans3(1)==0.) CYCLE
                    isgcd=(all_orbit%jj(d)+all_orbit%jj(a)+all_orbit%jj(p1)+ &
                         all_orbit%jj(p2)+all_orbit%jj(h1)+all_orbit%jj(h2))/2
                    sg=iph(isgcd+2*all_orbit%jj(h1)+all_orbit%jj(p1)+ &
                         all_orbit%jj(p2)+all_orbit%jj(a)+jtot-jpp)
                    factr=sixj1*sixj2*sg
                    w1=all_orbit%e(h1)+all_orbit%evalence(d)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(d)+ &
                         all_orbit%evalence(c)-all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
                    w3=all_orbit%evalence(d)-all_orbit%evalence(b)+all_orbit%e(h2)+ &
                         all_orbit%e(h1)+wcn
                    de=den1*den2
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)-factr*val1*val2*val3/de(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag18b
!
!
!
SUBROUTINE diag19a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: acmin, bdmin, acmax, bdmax, ppmin, ppmax, jph, p1, &
       h3, h2, h1, iph, nshell1, nshell2, isgcd, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, phmin, phmx
  REAL(DP) :: val1, val2, val3, factr, sixj1, sixj2
  LOGICAL dencheck

  diagram=0.
  acmin=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  bdmin=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  acMAX=(all_orbit%jj(a)+all_orbit%jj(c))/2
  bdMAX=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(acmin,bdmin)
  phmx=MIN(acmax,bdmax)
  DO jph=phmin,phmx
     sixj1=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     IF(sixj1 == 0.) CYCLE
     sixj1=sixj1*sqrt(2.*jph+1.)
     DO p1=1,all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
        DO h1=1,all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
           j2min=ABS((all_orbit%jj(c)-all_orbit%jj(h1))/2)
           j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
           j2max=(all_orbit%jj(c)+all_orbit%jj(h1))/2
           ppmin=MAX(j1min,j2min)
           ppmax=MIN(j1max,j2max)
           nshell1=all_orbit%nshell(c)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1)) CYCLE
           den1=wcn+all_orbit%e(h1)+all_orbit%evalence(c)-all_orbit%e(p1)-all_orbit%evalence(a)
           CALL cross_coupled_mtxel1(h1,b,p1,d,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           DO jpp=ppmin,ppMAX
              sixj2=sjs(all_orbit%jj(c),all_orbit%jj(a),2*jph, &
                   all_orbit%jj(p1),all_orbit%jj(h1),2*jpp)
              IF(sixj2 == 0.) CYCLE
              sixj2=sixj2*(2.*jpp+1.)
              DO h2=1,all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE      
                 DO h3=1,all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(h3)
                    nshell4=all_orbit%nshell(p1)+all_orbit%nshell(a)
                    idiff2=nshell4-nshell3
                    IF(dencheck(idiff2)) CYCLE
                    den2=wcn+all_orbit%e(h3)+all_orbit%e(h2)- &
                         all_orbit%evalence(a)-all_orbit%e(p1)
                    CALL pphhmtx(a,p1,h2,h3,jpp,ans2); IF(ans2(1)==0.) CYCLE
                    CALL pphhmtx(h2,h3,c,h1,jpp,ans3);IF(ans3(1)==0.) CYCLE
                    isgcd=(all_orbit%jj(d)+all_orbit%jj(p1))/2+ &
                         all_orbit%jj(a)+2*all_orbit%jj(h1)
                    factr=sixj2*sixj1*0.5*iph(isgcd+jtot+jph+jpp)
                    w1=all_orbit%e(h1)+all_orbit%evalence(d)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
                    w3=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%e(h3)+ &
                         all_orbit%evalence(c)-all_orbit%evalence(a)-all_orbit%e(p1)+wcn
                    de=den1*den2
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+factr*val1*val2*val3/de(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag19a
!
!
!
SUBROUTINE diag19b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: acmin, bdmin, acmax, bdmax, ppmin, ppmax, jph, p1, &
       h3, h2, h1, iph, nshell1, nshell2, isgcd, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, phmin, phmx
  REAL(DP) :: val1, val2, val3, factr, sixj1, sixj2
  LOGICAL dencheck

  diagram=0.
  acmin=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  bdmin=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  acMAX=(all_orbit%jj(a)+all_orbit%jj(c))/2
  bdMAX=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(acmin,bdmin)
  phmx=MIN(acmax,bdmax)
  IF(phmx.lt.phmin) return
  DO jph=phmin,phmx
     sixj1=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     IF(sixj1 == 0.) CYCLE
     sixj1=sixj1*sqrt(2.*jph+1.)
     DO p1=1,all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h1=1,all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           j1min=ABS((all_orbit%jj(a)-all_orbit%jj(h1))/2)
           j2min=ABS((all_orbit%jj(c)-all_orbit%jj(p1))/2)
           j1MAX=(all_orbit%jj(a)+all_orbit%jj(h1))/2
           j2MAX=(all_orbit%jj(c)+all_orbit%jj(p1))/2
           ppmin=MAX(j1min,j2min)
           ppmax=MIN(j1max,j2max)
           nshell1=all_orbit%nshell(d)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(b)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1)) CYCLE
           den1=wcn+all_orbit%e(h1)+all_orbit%evalence(d)-all_orbit%evalence(b)-all_orbit%e(p1)
           CALL cross_coupled_mtxel1(p1,b,h1,d,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           DO jpp=ppmin,ppmax
              sixj2=sjs(all_orbit%jj(c),all_orbit%jj(a),2*jph, &
                   all_orbit%jj(h1),all_orbit%jj(p1),2*jpp)
              IF(sixj2 == 0.) CYCLE
              sixj2=sixj2*(2*jpp+1)
              DO h2=1,all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 DO h3=1,all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(d)+ &
                         all_orbit%nshell(h3)
                    nshell4=all_orbit%nshell(p1)+all_orbit%nshell(b)+ &
                         all_orbit%nshell(a)
                    idiff2=nshell4-nshell3
                    IF(dencheck(idiff2)) CYCLE
                    den2=wcn+all_orbit%e(h3)+all_orbit%e(h2)+ &
                         all_orbit%evalence(d)-all_orbit%e(p1)-all_orbit%evalence(a)- &
                         all_orbit%evalence(b)
                    CALL pphhmtx(a,h1,h2,h3,jpp,ans2); IF(ans2(1)==0.) CYCLE
                    CALL pphhmtx(h2,h3,c,p1,jpp,ans3);IF(ans3(1)==0.) CYCLE
                    isgcd=(all_orbit%jj(d)+all_orbit%jj(h1))/2+ &
                         all_orbit%jj(a)+all_orbit%jj(h1)+all_orbit%jj(p1)
                    factr=sixj2*sixj1*0.5*iph(isgcd+jtot+jph+jpp)
                    w1=all_orbit%e(h1)+all_orbit%evalence(d)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%e(h3)+&
                         all_orbit%evalence(d)-all_orbit%evalence(b)-all_orbit%e(p1)+wcn
                    w3=all_orbit%e(h2)+all_orbit%e(h3)+all_orbit%evalence(d)+ &
                         all_orbit%evalence(c)-all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
                    de=den1*den2
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+factr*val1*val2*val3/de(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag19b
!
!     diagram 20 A
!
SUBROUTINE diag20a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot 
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: acmin, bdmin, acmax, bdmax, phmin, phmx, jph, p1, &
       p2, p3, h1, iph, nshell1, nshell2, isgcd, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jjpp, jppmin, jppmax
  REAL(DP) :: val1, val2, val3, factr, sixj1, sixj2, sg
  LOGICAL dencheck

  diagram=0.
  acmin=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  bdmin=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  acmax=(all_orbit%jj(a)+all_orbit%jj(c))/2
  bdmax=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(acmin,bdmin)
  phmx=MIN(acmax,bdmax)
  DO jph=phmin,phmx
     sixj1=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     IF(sixj1 == 0.) CYCLE
     sixj1=sixj1*SQRT(2.*jph+1.)
     DO p1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell1=all_orbit%nshell(c)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1)) CYCLE
           den1=wcn+all_orbit%evalence(c)-all_orbit%evalence(a)+ & 
                all_orbit%e(h1)-all_orbit%e(p1)
           j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
           j2min=ABS((all_orbit%jj(c)-all_orbit%jj(h1))/2)
           j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
           j2max=(all_orbit%jj(c)+all_orbit%jj(h1))/2
           jppmin=MAX(j1min,j2min)
           jppmax=MIN(j1max,j2max)
           IF(jppmax < jppmin) CYCLE
           CALL cross_coupled_mtxel1(h1,b,p1,d,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           DO jjpp=jppmin,jppmax
              sixj2=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                   2*jph,all_orbit%jj(p1),all_orbit%jj(h1),2*jjpp)
              IF(sixj2 == 0.) CYCLE
              sixj2=sixj2*(2*jjpp+1)
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO p3=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(c)+all_orbit%nshell(h1)
                    nshell4=all_orbit%nshell(p2)+all_orbit%nshell(p3)
                    idiff2=nshell4-nshell3
                    IF(dencheck(idiff2)) CYCLE
                    den2=wcn+all_orbit%evalence(c)+all_orbit%e(h1)-all_orbit%e(p3)-&
                         all_orbit%e(p2)
                    w1=all_orbit%e(h1)+all_orbit%evalence(d)+all_orbit%evalence(c)-&
                         all_orbit%evalence(a)+wcn
                    w2=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                    w3=w2
                    de=den1*den2
                    CALL pphhmtx(a,p1,p2,p3,jjpp,ans2); IF(ans2(1)==0.) CYCLE
                    CALL pphhmtx(p2,p3,c,h1,jjpp,ans3);IF(ans3(1)==0.) CYCLE
                    isgcd=(all_orbit%jj(d)+all_orbit%jj(p1))/2+&
                         all_orbit%jj(a)+2*all_orbit%jj(h1)
                    sg=iph(isgcd+jtot+jph+jjpp)
                    factr=sixj2*sg*sixj1
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+0.5*factr*val1*val2*val3/de(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag20a
!
!
!
SUBROUTINE diag20b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: acmin, bdmin, acmax, bdmax, phmin, phmx, jph, p1, &
       p2, p3, h1, iph, nshell1, nshell2, isgcd, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jjpp, jppmin, jppmax
  REAL(DP) :: val1, val2, val3, factr, sixj1, &
       sixj2, sg                 
  LOGICAL dencheck

  diagram=0.
  acmin=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  bdmin=ABS((all_orbit%jj(b)-all_orbit%jj(d))/2)
  acmax=(all_orbit%jj(a)+all_orbit%jj(c))/2
  bdmax=(all_orbit%jj(b)+all_orbit%jj(d))/2
  phmin=MAX(acmin,bdmin)
  phmx=MIN(acmax,bdmax)
  DO jph=phmin,phmx
     sixj1=sjs(all_orbit%jj(c),all_orbit%jj(d),2*jtot, &
          all_orbit%jj(b),all_orbit%jj(a),2*jph)
     IF(sixj1 == 0.) CYCLE
     sixj1=sixj1*sqrt(2.*jph+1.)
     DO p1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           nshell1=all_orbit%nshell(d)+all_orbit%nshell(h1)
           nshell2=all_orbit%nshell(p1)+all_orbit%nshell(b)
           idiff1=nshell2-nshell1
           IF(dencheck(idiff1)) CYCLE
           den1=wcn+all_orbit%evalence(d)-all_orbit%evalence(b)+ &
                all_orbit%e(h1)-all_orbit%e(p1)
           j1min=ABS((all_orbit%jj(a)-all_orbit%jj(h1))/2)
           j2min=ABS((all_orbit%jj(c)-all_orbit%jj(p1))/2)
           j1max=(all_orbit%jj(a)+all_orbit%jj(h1))/2
           j2max=(all_orbit%jj(c)+all_orbit%jj(p1))/2
           jppmin=MAX(j1min,j2min)
           jppmax=MIN(j1max,j2max)
           IF(jppmax < jppmin) CYCLE
           CALL cross_coupled_mtxel1(p1,b,h1,d,jph,ans1);IF ((ans1(1) == 0.)) CYCLE
           DO jjpp=jppmin,jppmax
              sixj2=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                   2*jph,all_orbit%jj(h1),all_orbit%jj(p1),2*jjpp)
              IF(sixj2 == 0.) CYCLE
              sixj2=sixj2*float(2*jjpp+1)
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE      
                 DO p3=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h1)+all_orbit%nshell(d)+&
                         all_orbit%nshell(c)
                    nshell4=all_orbit%nshell(p2)+all_orbit%nshell(b)+&
                         all_orbit%nshell(p3)
                    idiff2=nshell4-nshell3
                    IF(dencheck(idiff2)) CYCLE
                    den2=wcn+all_orbit%evalence(d)+all_orbit%e(h1)+ &
                         all_orbit%evalence(c)-all_orbit%evalence(b)- &
                         all_orbit%e(p3)-all_orbit%e(p2)
                    CALL pphhmtx(a,h1,p2,p3,jjpp,ans2); IF(ans2(1)==0.) CYCLE
                    CALL pphhmtx(p2,p3,c,p1,jjpp,ans3);IF(ans3(1)==0.) CYCLE
                    isgcd=(all_orbit%jj(d)+all_orbit%jj(h1))/2+&
                         all_orbit%jj(a)+all_orbit%jj(h1)+ &
                         all_orbit%jj(p1)
                    sg=iph(isgcd+jtot+jph+jjpp)
                    factr=sixj2*sg*sixj1
                    w1=all_orbit%e(h1)+all_orbit%evalence(d)+wcn
                    w2=all_orbit%evalence(c)-all_orbit%evalence(b)+w1
                    w3=w2
                    de=den1*den2
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+0.5*factr*val1*val2*val3/de(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag20b
!
!
!
SUBROUTINE diag21a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: h2, h3, h1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, p
  REAL(DP) :: val1, val2, val3, that, den1a
  LOGICAL dencheck

  diagram=0.
  den1a=0.5/(all_orbit%jj(a)+1.)
  DO h1=1, all_orbit%total_orbits
     IF(all_orbit%jj(h1) /= all_orbit%jj(a)) CYCLE
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     nshell1=all_orbit%nshell(h1)
     nshell2=all_orbit%nshell(a)
     idiff1=nshell2-nshell1
     IF(dencheck(idiff1)) CYCLE
     den1=wcn+all_orbit%e(h1)-all_orbit%evalence(a)
     DO p=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(p)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(p)-all_orbit%jj(h1))/2)
        j1max=(all_orbit%jj(p)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(h1)+all_orbit%jj(p))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        CALL pphhmtx(h1,b,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO jpp=j_min,j_max
           that=float(2*jpp+1)
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2)/= 'hole') CYCLE      
              DO h3=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h2)+all_orbit%nshell(h3)
                 nshell4=all_orbit%nshell(p)+all_orbit%nshell(a)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(h2)+all_orbit%e(h3)- &
                      all_orbit%evalence(a)-all_orbit%e(p)
                 CALL pphhmtx(p,a,h2,h3,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(h2,h3,p,h1,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%e(h1)+all_orbit%evalence(c)+ &
                      all_orbit%evalence(d)-all_orbit%evalence(a)+wcn
                 w2=all_orbit%e(h3)+all_orbit%e(h2)+wcn
                 w3=all_orbit%e(h3)+all_orbit%e(h2)+ &
                      all_orbit%e(h1)-all_orbit%evalence(a)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag21a
!
!
!
SUBROUTINE diag21b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: h2, h3, h1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, p
  REAL(DP) :: val1, val2, val3, den1a, that
  LOGICAL dencheck

  diagram=0.
  den1a=0.5/(all_orbit%jj(d)+1.)
  DO h1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
     nshell2=all_orbit%nshell(a)+all_orbit%nshell(b)
     idiff1=nshell1-nshell2
     IF(dencheck(idiff1)) CYCLE
     den1=wcn+all_orbit%e(h1)+all_orbit%evalence(c)- &
          all_orbit%evalence(a)-all_orbit%evalence(b)
     IF(all_orbit%jj(h1) /= all_orbit%jj(d)) CYCLE
     DO p=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(p)-all_orbit%jj(d))/2)
        j2min=ABS((all_orbit%jj(p)-all_orbit%jj(h1))/2)
        j1max=(all_orbit%jj(p)+all_orbit%jj(d))/2
        j2max=(all_orbit%jj(h1)+all_orbit%jj(p))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        CALL pphhmtx(a,b,c,h1,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO jpp=j_min,j_max
           that=float(2*jpp+1)
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              DO h3=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h2)+all_orbit%nshell(h3)+&
                      all_orbit%nshell(c)
                 nshell4=all_orbit%nshell(p)+all_orbit%nshell(a)+&
                      all_orbit%nshell(b)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(h2)+all_orbit%e(h3)+ &
                      all_orbit%evalence(c)-all_orbit%evalence(b)- &
                      all_orbit%evalence(a)-all_orbit%e(p)
                 CALL pphhmtx(h1,p,h2,h3,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(h2,h3,d,p,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                 w2=all_orbit%e(h3)+all_orbit%e(h2)+ &
                      all_orbit%e(h1)+all_orbit%evalence(c)- &
                      all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
                 w3=all_orbit%e(h3)+all_orbit%e(h2)+ &
                      all_orbit%evalence(c)+all_orbit%evalence(d)- &
                      all_orbit%evalence(b)-all_orbit%evalence(a)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag21b
!
!
!
SUBROUTINE diag22(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, p3, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, p4
  REAL(DP) :: val1, val2, val3
  LOGICAL dencheck

  diagram=0.
  w1=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
  w2=w1
  w3=w1
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
        nshell1=all_orbit%nshell(p1)+all_orbit%nshell(p2)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(d)
        idiff1=nshell1-nshell2
        IF(dencheck(idiff1))CYCLE
        den1=wcn+all_orbit%evalence(c)+all_orbit%evalence(d)-& 
             all_orbit%e(p1)-all_orbit%e(p2)
        CALL pphhmtx(p1,p2,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
        DO p3=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
           DO p4=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p4) == 'hole' ) CYCLE
              nshell3=all_orbit%nshell(p3)+all_orbit%nshell(p4)
              idiff2=nshell3-nshell2
              IF(dencheck(idiff2))CYCLE
              den2=wcn+all_orbit%evalence(c)+all_orbit%evalence(d)- &
                   all_orbit%e(p3)-all_orbit%e(p4)
              CALL pphhmtx(a,b,p3,p4,jtot,ans2); IF(ans2(1)==0.) CYCLE
              CALL pphhmtx(p3,p4,p1,p2,jtot,ans3);IF(ans3(1)==0.) CYCLE
              de=den1*den2
              DO ie=1,n_startenergy_veff
                 CALL interpolate(w1(ie),e_start_g,ans1,val1)
                 CALL interpolate(w2(ie),e_start_g,ans2,val2)
                 CALL interpolate(w3(ie),e_start_g,ans3,val3)
                 diagram(ie)=diagram(ie)+0.25*val1*val2*val3/de(ie)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag22
!
!
!
SUBROUTINE diag23(a,b,c,d,jtot,dg23a,dg23b)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: w1, w2, w3, &
       den1, den2, den3, w4, w5, w6, de2, de1
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: dg23a, dg23b
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3, &
       ans4, ans5, ans6
  INTEGER :: p2, p1, h1, h2, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, idiff3 
  REAL(DP) :: val1, val2, val3, val4, val5, val6
  LOGICAL dencheck

  dg23a=0.
  dg23b=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE 
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              nshell2=all_orbit%nshell(a)+all_orbit%nshell(b)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              nshell4=all_orbit%nshell(c)+all_orbit%nshell(d)
              idiff1=nshell1-nshell2
              idiff2=nshell1-nshell3
              idiff3=nshell1+nshell4-nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2).or.dencheck(idiff3)) CYCLE
              den1=wcn+all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(a)-all_orbit%evalence(b)
              den2=wcn+all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den3=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%e(p2)+ &
                   all_orbit%evalence(d)+all_orbit%evalence(c)- &
                   all_orbit%evalence(a)-all_orbit%evalence(b)+wcn
              CALL pphhmtx(a,b,p1,p2,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
              CALL pphhmtx(a,b,h1,h2,jtot,ans2); IF(ans2(1)==0.) CYCLE
              CALL pphhmtx(h1,h2,p1,p2,jtot,ans3);IF(ans3(1)==0.) CYCLE
              CALL pphhmtx(p1,p2,h1,h2,jtot,ans4); IF(ans4(1)==0.) CYCLE
              CALL pphhmtx(h1,h2,c,d,jtot,ans5); IF(ans5(1)==0.) CYCLE
              CALL pphhmtx(p1,p2,c,d,jtot,ans6); IF(ans6(1)==0.) CYCLE
              de1=den1*den2
              de2=den1*den3
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w4=w1
              w5=w1+all_orbit%evalence(c)+all_orbit%evalence(d)- &
                   all_orbit%evalence(a)-all_orbit%evalence(b)
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w3=w1+all_orbit%evalence(c)+all_orbit%evalence(d)- &
                   all_orbit%evalence(a)-all_orbit%evalence(b)
              w6=w3
              DO ie=1,n_startenergy_veff
                 CALL interpolate(w1(ie),e_start_g,ans1,val1)
                 CALL interpolate(w2(ie),e_start_g,ans2,val2)
                 CALL interpolate(w3(ie),e_start_g,ans3,val3)
                 CALL interpolate(w4(ie),e_start_g,ans4,val4)
                 CALL interpolate(w5(ie),e_start_g,ans5,val5)
                 CALL interpolate(w6(ie),e_start_g,ans6,val6)
                 dg23a(ie)=dg23a(ie)+0.25*val1*val4*val5/de1(ie)
                 dg23b(ie)=dg23b(ie)+0.25*val2*val3*val6/de2(ie)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag23
!
!
!
SUBROUTINE diag24a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, p3, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, h
  REAL(DP) :: val1, val2, val3, that, den1a
  LOGICAL dencheck

  diagram=0.
  den1a=0.5/(all_orbit%jj(c)+1.)
  w1=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
  DO p3=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
     nshell1=all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p3)
     idiff1=nshell1-nshell2
     IF(dencheck(idiff1)) CYCLE
     den1=wcn+all_orbit%evalence(c)-all_orbit%e(p3)
     IF(all_orbit%jj(p3) /= all_orbit%jj(c).or. &
          all_orbit%ll(p3) /= all_orbit%ll(c)) CYCLE
     CALL pphhmtx(a,b,p3,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
     DO h=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h)-all_orbit%jj(c))/2)
        j2min=ABS((all_orbit%jj(h)-all_orbit%jj(p3))/2)
        j1max=(all_orbit%jj(h)+all_orbit%jj(c))/2
        j2max=(all_orbit%jj(p3)+all_orbit%jj(h))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        DO jpp=j_min,j_max
           that=2.*jpp+1.
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h)+all_orbit%nshell(c)
                 nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(h)+all_orbit%evalence(c)-&
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL pphhmtx(h,p3,p1,p2,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(p1,p2,h,c,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w2=all_orbit%evalence(c)+all_orbit%e(h)+wcn
                 w3=w2
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag24a
!
!
!
SUBROUTINE diag24b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, p3, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, h
  REAL(DP) :: val1, val2, val3, that, den1a
  LOGICAL dencheck

  diagram=0.
  den1a=0.5d0/(all_orbit%jj(a)+1.)
  DO p3=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
     nshell1=all_orbit%nshell(c)+all_orbit%nshell(d)
     nshell2=all_orbit%nshell(p3)+all_orbit%nshell(b)
     idiff1=nshell1-nshell2
     IF(dencheck(idiff1)) CYCLE
     den1=wcn+all_orbit%evalence(c)+all_orbit%evalence(d)- &
          all_orbit%evalence(b)-all_orbit%e(p3)
     IF(all_orbit%jj(p3) /= all_orbit%jj(a).or. &
          all_orbit%ll(p3) /= all_orbit%ll(a)) CYCLE
     CALL pphhmtx(p3,b,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
     DO h=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(h)-all_orbit%jj(p3))/2)
        j1max=(all_orbit%jj(h)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(p3)+all_orbit%jj(h))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        DO jpp=j_min,j_max
           that=2.*jpp+1.
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE      
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h)+all_orbit%nshell(c)+ &
                      all_orbit%nshell(d)
                 nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)+ &
                      all_orbit%nshell(b)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(h)+all_orbit%evalence(c)+ &
                      all_orbit%evalence(d)-all_orbit%evalence(b)- &
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL pphhmtx(h,p3,p1,p2,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(p1,p2,h,a,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
                 w2=w1+all_orbit%e(h)-all_orbit%evalence(b)
                 w3=w2
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag24b
!
!
!
SUBROUTINE diag25a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, h2, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, h1
  REAL(DP) :: val1, val2, val3, that, den1a
  LOGICAL dencheck

  diagram=0.
  den1a=-0.5/(all_orbit%jj(c)+1.)
  DO p2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
     IF(all_orbit%jj(p2) /= all_orbit%jj(c).or.all_orbit%ll(p2) /= &
          all_orbit%ll(c)) CYCLE
     nshell1=all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p2)
     idiff1=nshell1-nshell2
     IF(dencheck(idiff1)) CYCLE
     den1=wcn+all_orbit%evalence(c)-all_orbit%e(p2)
     CALL pphhmtx(a,b,p2,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE 
     DO p1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(p1)-all_orbit%jj(c))/2)
        j2min=ABS((all_orbit%jj(p1)-all_orbit%jj(p2))/2)
        j1max=(all_orbit%jj(p1)+all_orbit%jj(c))/2
        j2max=(all_orbit%jj(p1)+all_orbit%jj(p2))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        DO jpp=j_min,j_max
           that=float((2*jpp+1))   
           DO h1=1, all_orbit%total_orbits
              DO h2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                 nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(h1)+all_orbit%e(h2)-&
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL pphhmtx(h1,h2,p1,c,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(p1,p2,h1,h2,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
                 w2=all_orbit%evalence(c)+all_orbit%e(h1)+ &
                      all_orbit%e(h2)-all_orbit%e(p2)+wcn
                 w3=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
 
END SUBROUTINE diag25a
!
!
!
SUBROUTINE diag25b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, h2, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, h1
  REAL(DP) :: val1, val2, val3, that, den1a
  LOGICAL dencheck

  diagram=0.
  den1a=-0.5/(all_orbit%jj(a)+1.)
  DO p2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
     nshell1=all_orbit%nshell(c)+all_orbit%nshell(d)
     nshell2=all_orbit%nshell(p2)+all_orbit%nshell(b)
     idiff1=nshell1-nshell2
     IF(dencheck(idiff1)) CYCLE
     den1=wcn+all_orbit%evalence(c)+all_orbit%evalence(d)- &
          all_orbit%evalence(b)-all_orbit%e(p2)
     IF(all_orbit%jj(p2) /= all_orbit%jj(a).or. &
          all_orbit%ll(p2) /= all_orbit%ll(a)) CYCLE
     CALL pphhmtx(p2,b,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
     DO p1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(p1)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(p1)-all_orbit%jj(p2))/2)
        j1max=(all_orbit%jj(p1)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        DO jpp=j_min,j_max
           that=float((2*jpp+1))
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              DO h2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                 nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(h1)+all_orbit%e(h2)+&
                      all_orbit%evalence(c)+all_orbit%evalence(d)-&
                      all_orbit%evalence(b)-all_orbit%e(p1)-&
                      all_orbit%e(p2)-all_orbit%evalence(a)
                 CALL pphhmtx(p1,a,h1,h2,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(h1,h2,p1,p2,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
                 w2=w1+all_orbit%e(h1)+all_orbit%e(h2)-&
                      all_orbit%e(p2)-all_orbit%evalence(b)
                 w3=w1+all_orbit%e(h1)+all_orbit%e(h2)-&
                      all_orbit%evalence(a)-all_orbit%evalence(b)
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag25b
!
!
!
SUBROUTINE diag26a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, p3, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, h
  REAL(DP) :: val1, val2, val3, that, den1a
  LOGICAL dencheck

  diagram=0.
  den1a=0.5/(all_orbit%jj(c)+1.)
  DO p3=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p3) /= 'hole' ) CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(p3)
     nshell2=all_orbit%nshell(a)+all_orbit%nshell(b)
     idiff1=nshell1-nshell2
     IF(dencheck(idiff1)) CYCLE
     w1=all_orbit%e(p3)+all_orbit%evalence(d)+wcn
     den1=wcn+all_orbit%evalence(d)+all_orbit%e(p3)-all_orbit%evalence(a)-all_orbit%evalence(b)
     IF(all_orbit%jj(p3) /= all_orbit%jj(c).or. &
          all_orbit%ll(p3) /= all_orbit%ll(c)) CYCLE
     CALL pphhmtx(a,b,p3,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
     DO h=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h) == 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h)-all_orbit%jj(c))/2)
        j2min=ABS((all_orbit%jj(h)-all_orbit%jj(p3))/2)
        j1max=(all_orbit%jj(h)+all_orbit%jj(c))/2
        j2max=(all_orbit%jj(p3)+all_orbit%jj(h))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        DO jpp=j_min,j_max
           that=2.*jpp+1.
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) /= 'hole' ) CYCLE 
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h)+all_orbit%nshell(b)+all_orbit%nshell(a)
                 nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)+all_orbit%nshell(d)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(p1)+all_orbit%evalence(d)+&
                      all_orbit%e(p2)-all_orbit%e(h)-all_orbit%evalence(a)-all_orbit%evalence(b)
                 CALL pphhmtx(h,p3,p1,p2,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(p1,p2,h,c,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w2=wcn+all_orbit%e(p3)+all_orbit%e(p2)+all_orbit%e(p1) +&
                    all_orbit%evalence(d)-all_orbit%evalence(a)-all_orbit%evalence(b)
                 w3=wcn+all_orbit%evalence(c)+all_orbit%e(p2)+all_orbit%e(p1) +&
                    all_orbit%evalence(d)-all_orbit%evalence(a)-all_orbit%evalence(b)
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag26a
!
!
!
SUBROUTINE diag26b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, p3, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, h
  REAL(DP) :: val1, val2, val3, that, den1a
  LOGICAL dencheck

  diagram=0.
  den1a=0.5d0/(all_orbit%jj(a)+1.)
  DO p3=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p3) /= 'hole' ) CYCLE
     nshell1=all_orbit%nshell(a)
     nshell2=all_orbit%nshell(p3)
     idiff1=nshell1-nshell2
     IF(dencheck(idiff1)) CYCLE
     den1=wcn+all_orbit%e(p3)-all_orbit%evalence(a)
     IF(all_orbit%jj(p3) /= all_orbit%jj(a).or. &
          all_orbit%ll(p3) /= all_orbit%ll(a)) CYCLE
     CALL pphhmtx(p3,b,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
     w1=all_orbit%evalence(c)+all_orbit%evalence(d)+all_orbit%e(p3)-all_orbit%evalence(a)+wcn
     DO h=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h) == 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(h)-all_orbit%jj(p3))/2)
        j1max=(all_orbit%jj(h)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(p3)+all_orbit%jj(h))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        DO jpp=j_min,j_max
           that=2.*jpp+1.
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) /= 'hole' ) CYCLE      
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) /= 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h)+all_orbit%nshell(a)
                 nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(p1)+all_orbit%e(p2)-all_orbit%evalence(a)- &
                      all_orbit%e(h)
                 CALL pphhmtx(h,p3,p1,p2,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(p1,p2,h,a,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w2=wcn+all_orbit%e(p1)+all_orbit%e(p2)
                 w3=w2+all_orbit%e(p3)-all_orbit%evalence(a)
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag26b
!
!
!
SUBROUTINE diag27a(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, h2, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, h1
  REAL(DP) :: val1, val2, val3, that, den1a
  LOGICAL dencheck

  diagram=0.
  den1a=-0.5/(all_orbit%jj(c)+1.)
  DO h2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(h2)
     nshell2=all_orbit%nshell(a)+all_orbit%nshell(b)
     idiff1=nshell1-nshell2
     IF(dencheck(idiff1)) CYCLE
     den1=wcn+all_orbit%evalence(d)+all_orbit%e(h2)-all_orbit%evalence(a)-all_orbit%evalence(b)
     IF(all_orbit%jj(h2) /= all_orbit%jj(c).or.all_orbit%ll(h2) /= &
          all_orbit%ll(c)) CYCLE
     CALL pphhmtx(a,b,h2,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(c))/2)
        j2min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j1max=(all_orbit%jj(h1)+all_orbit%jj(c))/2
        j2max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        DO jpp=j_min,j_max
           that=(2*jpp+1.)
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h1)+all_orbit%nshell(c)
                 nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff2=nshell4+nshell2-nshell3-nshell1
                 IF(dencheck(idiff2)) CYCLE
                 den2=den1+all_orbit%e(h1)+ &
                      all_orbit%evalence(c)-all_orbit%e(p1)- &
                      all_orbit%e(p2)
                 CALL pphhmtx(p1,p2,h1,c,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(h1,h2,p1,p2,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%e(h2)+all_orbit%evalence(d)+wcn
                 w2=w1+all_orbit%e(h1)+all_orbit%evalence(c)- &
                      all_orbit%evalence(a)-all_orbit%evalence(b)
                 w3=w2
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag27a
!
!
!
SUBROUTINE diag27b(a,b,c,d,jtot,diagram)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot  
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2, w3, &
       den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: diagram
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  INTEGER :: p2, h2, p1, nshell1, nshell2, idiff1, &
       idiff2, ie, nshell3, nshell4, j1min, j1max, j2min, &
       j2max, jpp, j_min, j_max, h1
  REAL(DP) :: val1, val2, val3, that, den1a
  LOGICAL dencheck

  diagram=0.
  den1a=-0.5/(all_orbit%jj(a)+1.)
  DO h2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
     nshell1=all_orbit%nshell(h2)
     nshell2=all_orbit%nshell(a)
     idiff1=nshell1-nshell2
     IF(dencheck(idiff1)) CYCLE
     den1=wcn+all_orbit%e(h2)-all_orbit%evalence(a)
     IF(all_orbit%jj(h2) /= all_orbit%jj(a).or.all_orbit%ll(h2) /= &
          all_orbit%ll(a)) CYCLE
     CALL pphhmtx(h2,b,c,d,jtot,ans1);IF ((ans1(1) == 0.)) CYCLE
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
        j1min=ABS((all_orbit%jj(h1)-all_orbit%jj(a))/2)
        j2min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j1max=(all_orbit%jj(h1)+all_orbit%jj(a))/2
        j2max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        j_min=MAX(j1min,j2min)
        j_max=MIN(j1max,j2max)
        IF(j_max < j_min) CYCLE
        DO jpp=j_min,j_max
           that=float((2*jpp+1))   
           DO p1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE     
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                 nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff2=nshell4-nshell3
                 IF(dencheck(idiff2)) CYCLE
                 den2=wcn+all_orbit%e(h1)+all_orbit%e(h2)-&
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL pphhmtx(h1,a,p1,p2,jpp,ans2); IF(ans2(1)==0.) CYCLE
                 CALL pphhmtx(p1,p2,h1,h2,jpp,ans3);IF(ans3(1)==0.) CYCLE
                 w1=all_orbit%evalence(c)+all_orbit%evalence(d)+ &
                      all_orbit%e(h2)-all_orbit%evalence(a)+wcn
                 w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                 w3=w2
                 de=den1*den2
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w1(ie),e_start_g,ans1,val1)
                    CALL interpolate(w2(ie),e_start_g,ans2,val2)
                    CALL interpolate(w3(ie),e_start_g,ans3,val3)
                    diagram(ie)=diagram(ie)+that*val1*val2*val3*den1a/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diag27b

!
!     Begin closed core diagrams
!
!
!     first-order contribution to the core-energy
!
SUBROUTINE core_E_first_order(first_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min
  INTEGER :: h1, h2, ie
  REAL(DP) :: val, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w
  REAL(DP), INTENT(OUT) :: first_order
  REAL(DP), DIMENSION(n_startenergy_g) :: ans

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        DO jtot=j_min,j_max
           fact=2*jtot+1.
           CALL pphhmtx(h1,h2,h1,h2,jtot,ans)
           IF ((ans(1) == 0.)) CYCLE
           w=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           DO ie=1,n_startenergy_veff
              CALL interpolate(w(ie),e_start_g,ans,val)
              diagram(ie)=diagram(ie)+fact*val
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  first_order=diagram(n_startenergy_veff/2+1)*0.5

END SUBROUTINE core_E_first_order
!
!  The second order diagram, 2p-2h
!
SUBROUTINE core_E_second_order(second_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min
  INTEGER :: h1, h2, p1, p2, ie, nshell1, nshell2, idiff
  REAL(DP), INTENT(OUT) :: second_order
  REAL(DP) :: val1, val2, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w, de
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        w=all_orbit%e(h1)+all_orbit%e(h2)+wcn
        DO jtot=j_min,j_max
           fact=2*jtot+1.
           DO p1=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 CALL pphhmtx(p1,p2,h1,h2,jtot,ans1)
                 CALL pphhmtx(h1,h2,p1,p2,jtot,ans2)
                 IF ( ans1(1) == 0.) CYCLE; IF ( ans2(1) == 0.) CYCLE
                 de=all_orbit%e(h1)+all_orbit%e(h2)+wcn -&
                      all_orbit%e(p1)-all_orbit%e(p2)
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w(ie),e_start_g,ans1,val1)
                    CALL interpolate(w(ie),e_start_g,ans2,val2)
                    diagram(ie)=diagram(ie)+val1*val2*fact/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  second_order=diagram(n_startenergy_veff/2+1)*0.25

END SUBROUTINE core_E_second_order


!
!  The second order diagram in  a recoupled way, 2p-2h
!
SUBROUTINE core_E_second_orderrecoupl(second_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min, iph
  INTEGER :: h1, h2, p1, p2, ie, nshell1, nshell2, idiff
  REAL(DP), INTENT(OUT) :: second_order
  REAL(DP) :: val1, val2, sg
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w, de
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(p1))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
        DO jtot=j_min,j_max
           DO h2=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              w=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                 nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff=nshell2-nshell1
!                 IF(dencheck(idiff)) CYCLE
                 CALL cross_coupled_mtxel1(p1,p2,h1,h2,jtot,ans1)
                 CALL cross_coupled_mtxel1(h1,h2,p1,p2,jtot,ans2)
                 sg=iph( (all_orbit%jj(p1)+all_orbit%jj(p2) +  &
                      3*all_orbit%jj(h1)+3*all_orbit%jj(h2))/2)
                 IF ( ans1(1) == 0.) CYCLE; IF ( ans2(1) == 0.) CYCLE
                 de=all_orbit%e(h1)+all_orbit%e(h2)+wcn -&
                      all_orbit%e(p1)-all_orbit%e(p2)
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w(ie),e_start_g,ans1,val1)
                    CALL interpolate(w(ie),e_start_g,ans2,val2)
                    diagram(ie)=diagram(ie)+sg*val1*val2/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  second_order=diagram(n_startenergy_veff/2+1)*0.25

END SUBROUTINE core_E_second_orderrecoupl
!
!  The second order diagram with two HF insertions
!
SUBROUTINE core_E_second_order_HF(second_order)
  USE constants
  USE single_particle_orbits              
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, p1, ie, jtot
  REAL(DP), INTENT(OUT) :: second_order
  REAL(DP) :: val1, val2, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, de
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO p1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits
        IF (all_orbit%jj(h2) /= all_orbit%jj(p1)) CYCLE
        IF (all_orbit%ll(h2) /= all_orbit%ll(p1)) CYCLE
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        de=all_orbit%e(h2)+wcn-all_orbit%e(p1)
        DO h1=1, all_orbit%total_orbits 
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h3)+1.))
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              CALL cross_coupled_mtxel1(h2,h3,p1,h3,jtot,ans2)
              IF ( ans2(1) == 0.) CYCLE
              DO ie=1,n_startenergy_veff
                 CALL interpolate(w1(ie),e_start_g,ans1,val1)
                 CALL interpolate(w2(ie),e_start_g,ans2,val2)
                 diagram(ie)=diagram(ie)+val1*val2*fact/de(ie)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  second_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE core_E_second_order_HF
!
!  The third order diagram, 4p-2h
!
SUBROUTINE core_E_third_order1(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min
  INTEGER :: h1, h2, p1, p2, p3, p4, ie
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        w=all_orbit%e(h1)+all_orbit%e(h2)+wcn
        DO jtot=j_min,j_max
           fact=2*jtot+1.
           DO p1=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 de1=all_orbit%e(h1)+all_orbit%e(h2)+wcn-& 
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL pphhmtx(h1,h2,p1,p2,jtot,ans1)
                 IF ((ans1(1) == 0.)) CYCLE
                 DO p3=1, all_orbit%total_orbits 
                    IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                    DO p4=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(p4) == 'hole' ) CYCLE
                       de2=all_orbit%e(h1)+all_orbit%e(h2)+wcn-&  
                            all_orbit%e(p3)-all_orbit%e(p4)
                       CALL pphhmtx(p1,p2,p3,p4,jtot,ans2)
                       IF(ans2(1) == 0.) CYCLE
                       CALL pphhmtx(p3,p4,h1,h2,jtot,ans3)
                       IF(ans3(1) == 0.) CYCLE
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w(ie),e_start_g,ans1,val1)
                          CALL interpolate(w(ie),e_start_g,ans2,val2)
                          CALL interpolate(w(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)+fact*val1*val2*val3/de1(ie)/de2(ie)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)*0.125

END SUBROUTINE core_E_third_order1
!
!  The third order diagram, 2p-4h
!
SUBROUTINE core_E_third_order2(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min
  INTEGER :: h1, h2, p1, p2, h3, h4, ie
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
        DO jtot=j_min,j_max
           fact=2*jtot+1.
           DO p1=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 de1=all_orbit%e(h1)+all_orbit%e(h2)+wcn -&  
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL pphhmtx(p1,p2,h1,h2,jtot,ans1)
                 IF ((ans1(1) == 0.)) CYCLE
                 DO h3=1, all_orbit%total_orbits 
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                    DO h4=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                       de2=all_orbit%e(h3)+all_orbit%e(h4)+wcn -&
                            all_orbit%e(p1)-all_orbit%e(p2)
                       w2 = w1 + all_orbit%e(h3)+all_orbit%e(h4) -&
                            all_orbit%e(p1)-all_orbit%e(p2)
                       w3 = de2
                       CALL pphhmtx(h1,h2,h3,h4,jtot,ans2)
                       IF(ans2(1) == 0.) CYCLE
                       CALL pphhmtx(h3,h4,p1,p2,jtot,ans3)
                       IF(ans3(1) == 0.) CYCLE
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w1(ie),e_start_g,ans1,val1)
                          CALL interpolate(w2(ie),e_start_g,ans2,val2)
                          CALL interpolate(w3(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)+fact*val1*val2*val3/de1(ie)/de2(ie) 
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)*0.125

END SUBROUTINE core_E_third_order2
!
!  The third order diagram, 3p-3h
!
SUBROUTINE core_E_third_order3(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min, nshell4, iph
  INTEGER :: h1, h2, h3, p1, p2, p3, ie
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact, sg, deno
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(p1))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
        DO jtot=j_min,j_max
           fact=SQRT(2*jtot+1.)
           DO h2=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 de1=all_orbit%e(h1)+all_orbit%e(h2)+wcn -&
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL cross_coupled_mtxel1(p1,p2,h1,h2,jtot,ans1) 
                 IF ( ans1(1) == 0.0_dp) CYCLE
                 DO h3=1, all_orbit%total_orbits 
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                    DO p3=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                       de2=all_orbit%e(h1)+all_orbit%e(h3)+wcn -&
                            all_orbit%e(p1)-all_orbit%e(p3)
                       w2 = all_orbit%e(h3)+all_orbit%e(h1) +wcn
                       w3 = w2 + all_orbit%e(h2)-all_orbit%e(p1)
                       sg=iph((all_orbit%jj(p2)+all_orbit%jj(h2)+ &
                            all_orbit%jj(h3)+all_orbit%jj(h1)+  &
                            all_orbit%jj(p1)+all_orbit%jj(p3))/2+ &
                            all_orbit%jj(h2)+all_orbit%jj(h1)+ &
                            all_orbit%jj(h3)-jtot)
                       CALL cross_coupled_mtxel1(h1,h3,p1,p3,jtot,ans2) 
                       IF ( ans2(1) == 0.0_dp) CYCLE
                       CALL cross_coupled_mtxel1(p3,h2,h3,p2,jtot,ans3) 
                       IF ( ans3(1) == 0.0_dp) CYCLE
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w1(ie),e_start_g,ans1,val1)
                          CALL interpolate(w2(ie),e_start_g,ans2,val2)
                          CALL interpolate(w3(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)+sg*fact*val1*val2*val3/de1(ie)/de2(ie)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE core_E_third_order3
!
!  The third order TDA diagram with two HF insertions
!
SUBROUTINE tda_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, p1, p2, ie, idiff1, idiff2, jtot
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO h1=1, all_orbit%total_orbits   ! outer hole line
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h2=1, all_orbit%total_orbits 
           IF (all_orbit%jj(h2) /= all_orbit%jj(p1)) CYCLE
           IF (all_orbit%ll(h2) /= all_orbit%ll(p1)) CYCLE
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
           idiff1=all_orbit%nshell(h2)-all_orbit%nshell(p1)
           IF(dencheck(idiff1)) CYCLE
           de1=all_orbit%e(h2)+wcn-all_orbit%e(p1)
           w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 IF (all_orbit%jj(h3) /= all_orbit%jj(p2)) CYCLE
                 IF (all_orbit%ll(h3) /= all_orbit%ll(p2)) CYCLE
                 idiff2=all_orbit%nshell(h3)-all_orbit%nshell(p2)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h3)+wcn-all_orbit%e(p2)
                 CALL cross_coupled_mtxel1(h2,p2,p1,h3,jtot,ans2)
                 IF ( ans2(1) == 0.) CYCLE
                 DO h4=1, all_orbit%total_orbits   ! outer hole line
                    IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                    fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h4)+1.))
                    CALL cross_coupled_mtxel1(h3,h4,p2,h4,jtot,ans3)
                    IF ( ans3(1) == 0.) CYCLE
                    w3=all_orbit%e(h4)+all_orbit%e(h3)+wcn
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+val1*val2*val3*fact/de1(ie)/de2(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE tda_E_third_order_HF
!
!  The third order RPA-1 diagram with two HF insertions
!
SUBROUTINE rpa1_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, p1, p2, ie, idiff1, idiff2, jtot
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO h1=1, all_orbit%total_orbits   ! outer hole line
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h2=1, all_orbit%total_orbits 
           IF (all_orbit%jj(h2) /= all_orbit%jj(p1)) CYCLE
           IF (all_orbit%ll(h2) /= all_orbit%ll(p1)) CYCLE
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
           idiff1=all_orbit%nshell(h2)-all_orbit%nshell(p1)
           IF(dencheck(idiff1)) CYCLE
           de1=all_orbit%e(h2)+wcn-all_orbit%e(p1)
           w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 IF (all_orbit%jj(h3) /= all_orbit%jj(p2)) CYCLE
                 IF (all_orbit%ll(h3) /= all_orbit%ll(p2)) CYCLE
                 idiff2=all_orbit%nshell(h3)+all_orbit%nshell(h2)-&
                      all_orbit%nshell(p2)-all_orbit%nshell(p1)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h3)+all_orbit%e(h2)+wcn-all_orbit%e(p2)-all_orbit%e(p1)
                 CALL cross_coupled_mtxel1(h2,h3,p1,p2,jtot,ans2)
                 IF ( ans2(1) == 0.) CYCLE
                 DO h4=1, all_orbit%total_orbits   ! outer hole line
                    IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                    fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h4)+1.))
                    CALL cross_coupled_mtxel1(p2,h4,h3,h4,jtot,ans3)
                    IF ( ans3(1) == 0.) CYCLE
                    w3=all_orbit%e(h4)+w2-all_orbit%e(p1)
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+val1*val2*val3*fact/de1(ie)/de2(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE rpa1_E_third_order_HF
!
!  The third order RPA-2 diagram with two HF insertions 
!
SUBROUTINE rpa2_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, p1, p2, ie, idiff1, idiff2, jtot
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO h1=1, all_orbit%total_orbits   ! outer hole line
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h2=1, all_orbit%total_orbits 
           IF (all_orbit%jj(h2) /= all_orbit%jj(p1)) CYCLE
           IF (all_orbit%ll(h2) /= all_orbit%ll(p1)) CYCLE
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
           idiff1=all_orbit%nshell(h2)-all_orbit%nshell(p1)
           IF(dencheck(idiff1)) CYCLE
           de1=all_orbit%e(h2)+wcn-all_orbit%e(p1)
           w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h1,h2,h1,p1,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 IF (all_orbit%jj(h3) /= all_orbit%jj(p2)) CYCLE
                 IF (all_orbit%ll(h3) /= all_orbit%ll(p2)) CYCLE
                 idiff2=all_orbit%nshell(h3)+all_orbit%nshell(h2)-&
                      all_orbit%nshell(p2)-all_orbit%nshell(p1)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h3)+all_orbit%e(h2)+wcn-all_orbit%e(p2)-all_orbit%e(p1)
                 CALL cross_coupled_mtxel1(p1,p2,h2,h3,jtot,ans2)
                 IF ( ans2(1) == 0.) CYCLE
                 DO h4=1, all_orbit%total_orbits   ! outer hole line
                    IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                    fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h4)+1.))
                    CALL cross_coupled_mtxel1(h3,h4,p2,h4,jtot,ans3)
                    IF ( ans3(1) == 0.) CYCLE
                    w3=all_orbit%e(h4)+w2-all_orbit%e(p1)
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+val1*val2*val3*fact/de1(ie)/de2(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE rpa2_E_third_order_HF
!
!  The third order diagram with 3  HF insertions, 4 holes, 2 particles
!
SUBROUTINE hf31_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, p1, p2, ie, idiff1, idiff2, jtot, jtot2, j_min, j_max
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO p2=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        IF (all_orbit%jj(h2) /= all_orbit%jj(p2)) CYCLE
        IF (all_orbit%ll(h2) /= all_orbit%ll(p2)) CYCLE
        idiff1=all_orbit%nshell(h2)-all_orbit%nshell(p2)
        IF(dencheck(idiff1)) CYCLE
        de1=all_orbit%e(h2)+wcn-all_orbit%e(p2)
        DO h4=1, all_orbit%total_orbits 
           IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
           w1=all_orbit%e(h4)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h2,h4,p2,h4,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              DO p1=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
                 IF (all_orbit%jj(p1) /= all_orbit%jj(p2)) CYCLE
                 IF (all_orbit%ll(p1) /= all_orbit%ll(p2)) CYCLE
                 w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
                 idiff2=all_orbit%nshell(h2)-all_orbit%nshell(p1)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h2)+wcn-all_orbit%e(p1)
                 j_min=ABS((all_orbit%jj(p2)-all_orbit%jj(p1))/2)
                 j_max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
                 DO jtot2 = j_min, j_max
                    CALL cross_coupled_mtxel1(h3,p2,h3,p1,jtot2,ans2)
                    IF ( ans2(1) == 0.) CYCLE
                    DO h1=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                       fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h3)+1.) &
                            *(all_orbit%jj(h4)+1.)/(all_orbit%jj(p2)+1.))
                       CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot2,ans3)
                       IF ( ans3(1) == 0.) CYCLE
                       w3=all_orbit%e(h2)+all_orbit%e(h1)+wcn
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w1(ie),e_start_g,ans1,val1)
                          CALL interpolate(w2(ie),e_start_g,ans2,val2)
                          CALL interpolate(w3(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)+val1*val2*val3*fact/de1(ie)/de2(ie)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE hf31_E_third_order_HF
!
!  The third order diagram with 3  HF insertions, 5 holes, 1 particle
!
SUBROUTINE hf32_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, h5, p1, ie, idiff1, idiff2, jtot, jtot2, j_min, j_max
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO p1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO h4=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
        IF (all_orbit%jj(h4) /= all_orbit%jj(p1)) CYCLE
        IF (all_orbit%ll(h4) /= all_orbit%ll(p1)) CYCLE
        idiff1=all_orbit%nshell(h4)-all_orbit%nshell(p1)
        IF(dencheck(idiff1)) CYCLE
        de1=all_orbit%e(h4)+wcn-all_orbit%e(p1)
        DO h5=1, all_orbit%total_orbits 
           IF(all_orbit%orbit_status(h5) /= 'hole' ) CYCLE
           w1=all_orbit%e(h4)+all_orbit%e(h5)+wcn
           CALL cross_coupled_mtxel1(h4,h5,p1,h5,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              DO h2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 IF (all_orbit%jj(p1) /= all_orbit%jj(h2)) CYCLE
                 IF (all_orbit%ll(p1) /= all_orbit%ll(h2)) CYCLE
                 w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
                 idiff2=all_orbit%nshell(h2)-all_orbit%nshell(p1)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h2)+wcn-all_orbit%e(p1)
                 j_min=ABS((all_orbit%jj(h2)-all_orbit%jj(h4))/2)
                 j_max=(all_orbit%jj(h2)+all_orbit%jj(h4))/2
                 DO jtot2 = j_min, j_max
                    CALL cross_coupled_mtxel1(h3,h2,h3,h4,jtot2,ans2)
                    IF ( ans2(1) == 0.) CYCLE
                    DO h1=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                       fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h3)+1.) &
                            *(all_orbit%jj(h5)+1.)/(all_orbit%jj(p1)+1.))
                       CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot2,ans3)
                       IF ( ans3(1) == 0.) CYCLE
                       w3=all_orbit%e(h2)+all_orbit%e(h1)+wcn
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w1(ie),e_start_g,ans1,val1)
                          CALL interpolate(w2(ie),e_start_g,ans2,val2)
                          CALL interpolate(w3(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)-val1*val2*val3*fact/de1(ie)/de2(ie)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE hf32_E_third_order_HF

!
!
!     Begin three-body diagrams
!
!
!     Sets up three-body contrib with one-body and two-body diagrams only

SUBROUTINE two_2_three(a,b,c,d,e,f,j_ab,j_de,jtot,two_d)
  USE constants
  USE single_particle_orbits
  USE stored_qbox
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, j_ab, j_de, j_ac, jtot, iph, &
       j_ac_min, j_ac_max, j_bc_min, j_bc_max, i
  REAL(DP) :: term, val, dij   
  REAL(DP), DIMENSION(n_startenergy_g) :: q1
  REAL(DP), DIMENSION(n_startenergy_veff ) :: w
  REAL(DP),  DIMENSION(n_startenergy_veff ), INTENT(OUT) :: two_d
  REAL(DP), DIMENSION(n_startenergy_veff) :: dg1,dg2,dg3,dg4, &
       dg5,dg6,dg7,dg8,dg9
  LOGICAL triag

  two_d=0.0_dp; dg1=0.0_dp; dg2=0.0_dp; dg3=0.0_dp; dg4=0.0_dp; dg5=0.0_dp; dg6=0.0_dp; dg7=0.0_dp;
  dg8=0.0_dp; dg9=0.0_dp
  j_bc_min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j_bc_max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j_ac_min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j_ac_max=(all_orbit%jj(a)+all_orbit%jj(c))/2

  ! The three-particle starting energy
  w=all_orbit%evalence(d)+all_orbit%evalence(e)+all_orbit%evalence(f)+wcn
  IF(c == f) THEN
     IF(j_ab == j_de) THEN
        CALL pphhmtx(a,b,d,e,j_ab,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg1(i)=dg1(i)+val
        ENDDO
     ENDIF
  ENDIF

  IF(c == d) THEN
     term=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
          all_orbit%jj(f),jtot,2*j_ab)*SQRT((2.*j_ab+1.)*(2.*j_de+1.))&
          *iph((all_orbit%jj(e)+all_orbit%jj(f)+all_orbit%jj(d)+jtot)/2) &
          *iph((ABS(all_orbit%jj(d)-jtot)/2)+j_ab)
     CALL pphhmtx(a,b,e,f,j_ab,q1)
     DO i=1, n_startenergy_veff
        CALL interpolate(w(i),e_start_g,q1,val)
        dg2(i)=dg2(i)+val*term
     ENDDO
  ENDIF

  IF(c == e) THEN
     term=-sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,jtot, &
          all_orbit%jj(f),2*j_ab)*SQRT((2.*j_ab+1.)*(2.*j_de+1.)) &
          *iph((all_orbit%jj(e)+all_orbit%jj(f))/2+j_ab+j_de)
     CALL pphhmtx(a,b,d,f,j_ab,q1)
     DO i=1, n_startenergy_veff
        CALL interpolate(w(i),e_start_g,q1,val)
        dg3(i)=dg3(i)+val*term
     ENDDO
  ENDIF

  IF(b == f) THEN
     term=-sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
          all_orbit%jj(c),2*j_de)*SQRT((2.*j_ab+1.)*(2.*j_de+1.))&
          *iph((all_orbit%jj(b)+all_orbit%jj(c))/2+j_de+j_ab)
     CALL pphhmtx(a,c,d,e,j_de,q1)
     DO i=1, n_startenergy_veff
        CALL interpolate(w(i),e_start_g,q1,val)
        dg4(i)=dg4(i)+val*term
     ENDDO
  ENDIF

  IF(b == d) THEN
     DO j_ac=j_ac_min,j_ac_max
        IF(triag(2*j_ac,all_orbit%jj(b),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(d),jtot)) CYCLE
        term=sjs(all_orbit%jj(a),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(b),2*j_ab) &
             *sjs(all_orbit%jj(e),all_orbit%jj(f),2*j_ac,jtot, &
             all_orbit%jj(d),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1) &
             *iph((all_orbit%jj(d)+all_orbit%jj(e))/2-j_de+1)
        CALL pphhmtx(a,c,e,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg5(i)=dg5(i)+val*term
        ENDDO
     ENDDO
  ENDIF

  IF(b == e) THEN
     DO j_ac=j_ac_min,j_ac_max
        IF(triag(2*j_ac,all_orbit%jj(b),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(e),jtot)) CYCLE
        term=sjs(all_orbit%jj(a),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(b),2*j_ab) &
             *sjs(all_orbit%jj(d),all_orbit%jj(f),2*j_ac, &
             jtot,all_orbit%jj(e),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1)
        CALL pphhmtx(a,c,d,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg6(i)=dg6(i)+val*term
        ENDDO
     ENDDO
  ENDIF

  IF(a == f) THEN
     term=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,all_orbit%jj(c), &
          jtot,2*j_de)*SQRT((2.*j_ab+1.)*(2.*j_de+1.)) &
          *iph((all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(a)+jtot)/2)&
          *iph(ABS(all_orbit%jj(a)-jtot)/2+j_de)
     CALL pphhmtx(b,c,d,e,j_de,q1)
     DO i=1, n_startenergy_veff
        CALL interpolate(w(i),e_start_g,q1,val)
        dg7(i)=dg7(i)+val*term
     ENDDO
  ENDIF

  IF(a == d) THEN
     DO j_ac=j_bc_min,j_bc_max
        IF(triag(2*j_ac,all_orbit%jj(a),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(d),jtot)) CYCLE
        term=sjs(all_orbit%jj(b),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(a),2*j_ab) &
             *sjs(all_orbit%jj(e),all_orbit%jj(f),2*j_ac,jtot, &
             all_orbit%jj(d),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1) &
             *iph((all_orbit%jj(d)+all_orbit%jj(e)+all_orbit%jj(b)+ &
             all_orbit%jj(a))/2-j_ab-j_de)
        CALL pphhmtx(b,c,e,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg8(i)=dg8(i)+val*term
        ENDDO
     ENDDO
  ENDIF

  IF(a == e) THEN
     DO j_ac=j_bc_min,j_bc_max
        IF(triag(2*j_ac,all_orbit%jj(a),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(e),jtot)) CYCLE
        term=sjs(all_orbit%jj(b),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(a),2*j_ab) &
             *sjs(all_orbit%jj(d),all_orbit%jj(f),2*j_ac,jtot, &
             all_orbit%jj(e),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1) &
             *iph((all_orbit%jj(b)+all_orbit%jj(a))/2-j_ab+1)
        CALL pphhmtx(b,c,d,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg9(i)=dg9(i)+val*term
        ENDDO
     ENDDO
  ENDIF
  two_d=dg1+dg2+dg3+dg4+dg5+dg6+dg7+dg8+dg9

END SUBROUTINE two_2_three


!
!
!     Begin three-body diagrams
!
!
!     Sets up three-body contrib with one-body and two-body diagrams only

SUBROUTINE two_2_threex(a,b,c,d,e,f,j_ab,j_de,jtot,two_d)
  USE constants
  USE single_particle_orbits
  USE stored_qbox
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, j_ab, j_de, j_ac, jtot, iph, j_ac_min, j_ac_max, i
  REAL(DP) :: term, val, dij
  REAL(DP), DIMENSION(n_startenergy_g) :: q1
  REAL(DP), DIMENSION(n_startenergy_veff ) :: w
  REAL(DP),  DIMENSION(n_startenergy_veff ), INTENT(OUT) :: two_d
  REAL(DP), DIMENSION(n_startenergy_veff) :: dg1
  LOGICAL triag

  two_d = 0.0_dp; dg1 = 0.0_dp
  j_ac_min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j_ac_max=(all_orbit%jj(a)+all_orbit%jj(c))/2

  ! The three-particle starting energy
  w=all_orbit%evalence(d)+all_orbit%evalence(e)+all_orbit%evalence(f)+wcn
  
  IF(b == e) THEN
     DO j_ac=j_ac_min,j_ac_max
        IF(triag(2*j_ac,all_orbit%jj(b),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(e),jtot)) CYCLE
        term=sjs(all_orbit%jj(a),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(b),2*j_ab) &
             *sjs(all_orbit%jj(d),all_orbit%jj(f),2*j_ac, &
             jtot,all_orbit%jj(e),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1)
        CALL pphhmtx(a,c,d,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg1(i)=dg1(i)+val*term
        ENDDO
     ENDDO
  ENDIF
  two_d = dg1

END SUBROUTINE two_2_threex



!
!              Three-body diagram 1 with particle intermediate state
!

SUBROUTINE three_body_1a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(e)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(p)+all_orbit%jj(b)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2)
     CALL pphhmtx(a,p,d,e,j_de,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(d)+all_orbit%e(e)-all_orbit%e(p)-all_orbit%e(a)
     w1=all_orbit%e(d)+all_orbit%e(e)+wcn
     w2=all_orbit%e(f)-all_orbit%e(a)+all_orbit%e(d)+ &
          all_orbit%e(e)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)
        sixj2=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_de, &
             all_orbit%jj(f),jtot,2*j_bc)
        fact2=fact1*sixj1*sixj2*(2.*j_bc+1.)
        CALL pphhmtx(b,c,p,f,j_bc,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_1a
!
!              Three-body diagram 1 with hole intermediate state
!
SUBROUTINE three_body_1b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(b)+all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(f)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))*iph((all_orbit%jj(p)+ &
          all_orbit%jj(b)+all_orbit%jj(f)+all_orbit%jj(c))/2+all_orbit%jj(p))
     CALL pphhmtx(a,p,d,e,j_de,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(p)+all_orbit%e(f)-all_orbit%e(b)-all_orbit%e(c)
     w1=all_orbit%e(p)+all_orbit%e(f)+all_orbit%e(d)+all_orbit%e(e)- &
          all_orbit%e(b)-all_orbit%e(c)+wcn
     w2=all_orbit%e(p)+all_orbit%e(f)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c), jtot,2*j_bc)
        sixj2=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_de, &
             all_orbit%jj(f),jtot,2*j_bc)
        fact2=fact1*sixj1*sixj2*(2.*j_bc+1.)
        CALL pphhmtx(b,c,p,f,j_bc,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_1b



SUBROUTINE three_body_2a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(e)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(b)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph(ABS((all_orbit%jj(c)-all_orbit%jj(f))/2+j_ab-j_de))
     CALL pphhmtx(p,b,d,e,j_de,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     IF(ANS1(2) == 0.)CYCLE
     den=wcn+all_orbit%e(d)+all_orbit%e(e)-all_orbit%e(p)-all_orbit%e(b)
     w1=all_orbit%e(d)+all_orbit%e(e)+wcn
     w2=all_orbit%e(f)-all_orbit%e(b)+all_orbit%e(d)+ &
          all_orbit%e(e)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)
        sixj2=sjs(all_orbit%jj(p),all_orbit%jj(b),2*j_de,jtot, &
             all_orbit%jj(f),2*j_ac)
        fact2=fact1*sixj1*sixj2*(2.*j_ac+1.)
        CALL pphhmtx(a,c,p,f,j_ac,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_2a




SUBROUTINE three_body_2b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, &
       sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(a)+all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(f)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph(ABS((all_orbit%jj(c)-all_orbit%jj(f))/2+j_ab- &
          j_de+all_orbit%jj(p)))
     CALL pphhmtx(p,b,d,e,j_de,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(f)+all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)+ &
          all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(c)+wcn
     w2=all_orbit%e(f)+all_orbit%e(p)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)
        sixj2=sjs(all_orbit%jj(p),all_orbit%jj(b),2*j_de,jtot, &
             all_orbit%jj(f),2*j_ac)
        fact2=fact1*sixj1*sixj2*(2.*j_ac+1.)
        CALL pphhmtx(a,c,p,f,j_ac,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        IF(ans2(2) == 0.) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                den(ie)
        ENDDO
     ENDDO
  ENDDO


END SUBROUTINE three_body_2b


SUBROUTINE three_body_3a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, iph,  &
       nshell1, nshell2, idiff, jtot
  REAL(DP) :: val1, val2, fact                         
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2 
  REAL(DP), DIMENSION(n_startenergy_veff) :: den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(e)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(c)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de)* &
          sjs(all_orbit%jj(p),all_orbit%jj(f),2*j_ab,jtot, &
          all_orbit%jj(c),2*j_de)
     CALL pphhmtx(a,b,p,f,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     CALL pphhmtx(p,c,d,e,j_de,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
     den=wcn+all_orbit%e(d)+all_orbit%e(e)-all_orbit%e(p)-all_orbit%e(c)
     w2=all_orbit%e(d)+all_orbit%e(e)+wcn
     w1=all_orbit%e(f)-all_orbit%e(a)+all_orbit%e(d)+ &
          all_orbit%e(e)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        CALL interpolate(w2(ie),e_start_g,ans2,val2)
        three_diagram(ie)=three_diagram(ie)+fact*val1*val2/den(ie)
     ENDDO
  ENDDO

END SUBROUTINE three_body_3a




SUBROUTINE three_body_3b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, iph,  &
       nshell1, nshell2, idiff, jtot
  REAL(DP) :: val1, val2, fact                          
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(b)+all_orbit%nshell(a)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(f)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+ &
          all_orbit%jj(p))* &
          sjs(all_orbit%jj(p),all_orbit%jj(f),2*j_ab,jtot, &
          all_orbit%jj(c),2*j_de)
     CALL pphhmtx(a,b,p,f,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     CALL pphhmtx(p,c,d,e,j_de,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
     den=wcn+all_orbit%e(p)+all_orbit%e(f)-all_orbit%e(a)-all_orbit%e(b)
     w1=all_orbit%e(p)+all_orbit%e(f)+wcn
     w2=all_orbit%e(f)-all_orbit%e(a)+all_orbit%e(d)+all_orbit%e(e)- &
          all_orbit%e(b)+all_orbit%e(p)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        CALL interpolate(w2(ie),e_start_g,ans2,val2)
        three_diagram(ie)=three_diagram(ie)+fact*val1*val2/den(ie)
     ENDDO
  ENDDO

END SUBROUTINE three_body_3b


SUBROUTINE three_body_4a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ef, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j1max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(f)+all_orbit%nshell(e)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(c)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(p)+all_orbit%jj(c)+all_orbit%jj(f)+ &
          all_orbit%jj(e))/2)
     CALL pphhmtx(a,b,d,p,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(f)+all_orbit%e(e)-all_orbit%e(p)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(c)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_ef=j1min,j1max
        sixj1=sjs(all_orbit%jj(d),all_orbit%jj(p),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_ef)
        sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
             all_orbit%jj(f),jtot,2*j_ef)
        fact2=fact1*sixj1*sixj2*(2.*j_ef+1.)
        CALL pphhmtx(p,c,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_4a


SUBROUTINE three_body_4b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ef, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j1max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(a)+all_orbit%nshell(b)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(d)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(p)+all_orbit%jj(c)+all_orbit%jj(f)+ &
          all_orbit%jj(e))/2+all_orbit%jj(p))
     CALL pphhmtx(a,b,d,p,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(d)+all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(b)
     w1=all_orbit%e(d)+all_orbit%e(p)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+all_orbit%e(d)+ &
          all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(b)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_ef=j1min,j1max
        sixj1=sjs(all_orbit%jj(d),all_orbit%jj(p),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_ef)
        sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
             all_orbit%jj(f),jtot,2*j_ef)
        fact2=fact1*sixj1*sixj2*(2.*j_ef+1.)
        CALL pphhmtx(p,c,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_4b


SUBROUTINE three_body_5a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j1max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(c)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph(ABS((all_orbit%jj(c)-all_orbit%jj(f))/2+j_ab-j_de))
     CALL pphhmtx(a,b,p,e,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(d)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(c)+wcn
     w2=all_orbit%e(f)+all_orbit%e(d)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_df=j1min,j1max
        sixj1=sjs(all_orbit%jj(p),all_orbit%jj(e),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_df)
        sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,jtot, &
             all_orbit%jj(f),2*j_df)
        fact2=fact1*sixj1*sixj2*(2.*j_df+1.)
        CALL pphhmtx(p,c,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_5a




SUBROUTINE three_body_5b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j1max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(a)+all_orbit%nshell(b)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(e)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph(ABS((all_orbit%jj(c)-all_orbit%jj(f))/2+ &
          j_ab-j_de+all_orbit%jj(p)))
     CALL pphhmtx(a,b,p,e,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(p)+all_orbit%e(e)-all_orbit%e(a)-all_orbit%e(b)
     w1=all_orbit%e(e)+all_orbit%e(p)+wcn
     w2=all_orbit%e(e)+all_orbit%e(p)+all_orbit%e(f)+ &
          all_orbit%e(d)-all_orbit%e(a)-all_orbit%e(b)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_df=j1min,j1max
        sixj1=sjs(all_orbit%jj(p),all_orbit%jj(e),2*j_ab, &
             jtot,all_orbit%jj(c),2*j_df)
        sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
             jtot,all_orbit%jj(f),2*j_df)
        fact2=fact1*sixj1*sixj2*(2.*j_df+1.)
        CALL pphhmtx(p,c,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_5b



SUBROUTINE three_body_6a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot, j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(e)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2+j_de)
     den=wcn+all_orbit%e(d)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(a)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(a)+wcn
     w2=all_orbit%e(f)+all_orbit%e(d)+wcn
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)*(2.*j_bc+1.)
        CALL pphhmtx(b,c,e,p,j_bc,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_df=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                jtot,all_orbit%jj(f),2*j_df)
           sixj3=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_df, &
                all_orbit%jj(e),jtot,2*j_bc)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_df+1.)*iph(j_bc+j_df)
           CALL pphhmtx(a,p,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_6a



SUBROUTINE three_body_6b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph, j_df,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(b)+all_orbit%nshell(b)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(e)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(e)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2+j_de+all_orbit%jj(p))
     den=wcn+all_orbit%e(p)+all_orbit%e(e)-all_orbit%e(b)-all_orbit%e(c)
     w1=all_orbit%e(p)+all_orbit%e(e)+wcn
     w2=all_orbit%e(f)+all_orbit%e(d)+all_orbit%e(e)+ &
          all_orbit%e(p)-all_orbit%e(b)-all_orbit%e(c)+wcn
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)*(2.*j_bc+1.)
        CALL pphhmtx(b,c,e,p,j_bc,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_df=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                jtot,all_orbit%jj(f),2*j_df)
           sixj3=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_df, &
                all_orbit%jj(e),jtot,2*j_bc)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_df+1.)*iph(j_bc+j_df)
           CALL pphhmtx(a,p,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_6b


SUBROUTINE three_body_7a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(b)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_de+j_ab)
     den=wcn+all_orbit%e(d)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(b)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(b)+wcn
     w2=all_orbit%e(f)+all_orbit%e(d)+wcn
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)*(2.*j_ac+1.)
        CALL pphhmtx(a,c,p,e,j_ac,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_df=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,jtot, &
                all_orbit%jj(f),2*j_df)
           sixj3=sjs(all_orbit%jj(p),all_orbit%jj(e),2*j_ac,jtot, &
                all_orbit%jj(b),2*j_df)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_df+1.)
           CALL pphhmtx(p,b,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_7a


SUBROUTINE three_body_7b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(c)+all_orbit%nshell(a)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(e)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_de+j_ab+ &
          all_orbit%jj(p))
     den=wcn+all_orbit%e(p)+all_orbit%e(e)-all_orbit%e(c)-all_orbit%e(a)
     w1=all_orbit%e(p)+all_orbit%e(e)+wcn
     w2=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)+ &
          all_orbit%e(p)-all_orbit%e(c)-all_orbit%e(a)+wcn
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)*(2.*j_ac+1.)
        CALL pphhmtx(a,c,p,e,j_ac,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_df=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,jtot, &
                all_orbit%jj(f),2*j_df)
           sixj3=sjs(all_orbit%jj(p),all_orbit%jj(e),2*j_ac,jtot, &
                all_orbit%jj(b),2*j_df)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_df+1.)
           CALL pphhmtx(p,b,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_7b


SUBROUTINE three_body_8a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph, j_ef,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(e)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(e)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2+jtot)
     den=wcn+all_orbit%e(e)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(a)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(a)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+wcn
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)*(2.*j_bc+1.)
        CALL pphhmtx(b,c,d,p,j_bc,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_ef=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                all_orbit%jj(f),jtot,2*j_ef)
           sixj3=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_ef, &
                all_orbit%jj(d),jtot,2*j_bc)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_ef+1.)*iph(j_bc+j_ef)
           CALL pphhmtx(a,p,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_8a



SUBROUTINE three_body_8b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph, j_ef,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(b)+all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(d)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(e)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2+jtot+all_orbit%jj(p))
     den=wcn+all_orbit%e(p)+all_orbit%e(d)-all_orbit%e(b)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(p)+wcn
     w2=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)+ &
          all_orbit%e(p)-all_orbit%e(b)-all_orbit%e(c)+wcn
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)*(2.*j_bc+1.)
        CALL pphhmtx(b,c,d,p,j_bc,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_ef=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                all_orbit%jj(f),jtot,2*j_ef)
           sixj3=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_ef, &
                all_orbit%jj(d),jtot,2*j_bc)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_ef+1.)*iph(j_bc+j_ef)
           CALL pphhmtx(a,p,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                   den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_8b



SUBROUTINE three_body_9a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, j_ef, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(e)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(b)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(f)+ &
          all_orbit%jj(e))/2+j_ab)
     den=wcn+all_orbit%e(e)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(b)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(b)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+wcn
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)*(2.*j_ac+1.)
        CALL pphhmtx(a,c,d,p,j_ac,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_ef=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                all_orbit%jj(f),jtot,2*j_ef)
           sixj3=sjs(all_orbit%jj(d),all_orbit%jj(p),2*j_ac, &
                all_orbit%jj(b),jtot,2*j_ef)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_ef+1.)*iph(j_ac+j_ef)
           CALL pphhmtx(b,p,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                   den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO


END SUBROUTINE three_body_9a


SUBROUTINE three_body_9b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, j_ef, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(a)+all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(d)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(f)+ &
          all_orbit%jj(e))/2+j_ab+all_orbit%jj(p))
     den=wcn+all_orbit%e(d)+all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(p)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+all_orbit%e(d)+ &
          all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(c)+wcn         
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)*(2.*j_ac+1.)
        CALL pphhmtx(a,c,d,p,j_ac,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_ef=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,&
                all_orbit%jj(f),jtot,2*j_ef)
           sixj3=sjs(all_orbit%jj(d),all_orbit%jj(p),2*j_ac, &
                all_orbit%jj(b),jtot,2*j_ef)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_ef+1.)*iph(j_ac+j_ef)
           CALL pphhmtx(b,p,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                   den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_9b
!
!     Begin one-body effective operator diagrams
!
!
!     effective operator diagram with ph intermdediate state, first order in G
!
SUBROUTINE effective_operator_1(a,c,effoper_diagram_1)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, l, p,  nshell1, nshell2, idiff, iphase, iph
  REAL(DP) :: val, factr
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_1
  REAL(DP), DIMENSION(n_startenergy_g) :: ans
  LOGICAL dencheck
  effoper_diagram_1=0.
  DO h=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
     DO p=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
        nshell1=all_orbit%nshell(a)+all_orbit%nshell(p)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(h)
        idiff=nshell1-nshell2
        IF(dencheck(idiff)) CYCLE
        de=all_orbit%evalence(c)+all_orbit%e(h)-all_orbit%evalence(a)-all_orbit%e(p)+wcn
        IF(bare_operator(h,p) == 0.) CYCLE
        iphase=iph((3*all_orbit%jj(h)+all_orbit%jj(p))/2-lambda)
        factr=iphase*bare_operator(h,p)/SQRT(2.*lambda+1.)
        CALL cross_coupled_mtxel1(a,p,c,h,lambda,ans)
        IF(ans(2) == 0.) CYCLE
        w=all_orbit%evalence(c)+all_orbit%e(h)+wcn
        DO l=1,n_startenergy_veff
           CALL interpolate(w(l),e_start_g,ans,val)
           effoper_diagram_1(l)=effoper_diagram_1(l)+val*factr/de(l)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_1
!
!     effective operator diagram with ph intermdediate state, first order in G
!     hermitian conjugate of previous diagram
!
SUBROUTINE effective_operator_2(a,c,effoper_diagram_2)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, l, p,  nshell1, nshell2, idiff, iphase, iph
  REAL(DP) :: val, factr
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w
  REAL(DP), DIMENSION(n_startenergy_g) :: ans
  LOGICAL dencheck

  effoper_diagram_2=0.
  DO h=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
     DO p=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
        nshell1=all_orbit%nshell(a)+all_orbit%nshell(p)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(h)
        idiff=nshell1-nshell2
        IF(dencheck(idiff)) CYCLE
        IF(bare_operator(p,h) == 0.) CYCLE
        iphase=iph((3*all_orbit%jj(h)+all_orbit%jj(p))/2-lambda)
        factr=iphase*bare_operator(p,h)/SQRT(2.*lambda+1.)
        CALL cross_coupled_mtxel1(a,h,c,p,lambda,ans)
        IF(ans(2) == 0.) CYCLE
        de=all_orbit%evalence(a)+all_orbit%e(h)-all_orbit%evalence(c)-all_orbit%e(p)+wcn
        w=all_orbit%evalence(a)+all_orbit%e(h)+wcn
        DO l=1,n_startenergy_veff
           CALL interpolate(w(l),e_start_g,ans,val)
           effoper_diagram_2(l)=effoper_diagram_2(l)+val*factr/de(l)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_2
!
!
!
SUBROUTINE effective_operator_3(a,c,effoper_diagram_3)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: w1, w2, de
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_3
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_3=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(h2,p2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(a)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(a)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(c)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)- &
                   all_orbit%evalence(a)-all_orbit%e(p1)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)- &
                   all_orbit%evalence(a)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,p1,c,h1,lambda,ans1)
              CALL cross_coupled_mtxel1(h1,p2,p1,h2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3*all_orbit%jj(h2)+&
                   all_orbit%jj(p1)+all_orbit%jj(p2))/2) &
                   *bare_operator(h2,p2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              w2=w1+all_orbit%e(h2)-all_orbit%evalence(a)
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_3(l)=effoper_diagram_3(l)+& 
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_3
!
!
!
SUBROUTINE effective_operator_4(a,c,effoper_diagram_4)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_4
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_4=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE         
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE 
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(p2,h2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(c)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(c)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(a)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(a)- &
                   all_orbit%evalence(c)-all_orbit%e(p1)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)- &
                   all_orbit%evalence(c)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,h1,c,p1,lambda,ans1)
              CALL cross_coupled_mtxel1(p1,h2,h1,p2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3*all_orbit%jj(h2)+&
                   all_orbit%jj(p1)+all_orbit%jj(p2))/2) &
                   *bare_operator(p2,h2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(a)+wcn
              w2=w1+all_orbit%e(h2)-all_orbit%evalence(c)
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_4(l)=effoper_diagram_4(l)+&
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_4
!
!
!
subroutine effective_operator_5(a,c,effoper_diagram_5)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff )  :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_5
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_5=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE    
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE 
              IF(bare_operator(p2,h2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(a)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)-& 
                   all_orbit%evalence(a)-all_orbit%e(p1)
              den2=all_orbit%e(h1)+all_orbit%e(h2)-& 
                   all_orbit%e(p1)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,p1,c,h1,lambda,ans1)
              CALL cross_coupled_mtxel1(h1,h2,p1,p2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3* &
                   all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                   all_orbit%jj(p2))/2) &
                   *bare_operator(p2,h2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_5(l)=effoper_diagram_5(l)+ &
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_5
!
!
!
subroutine effective_operator_6(a,c,effoper_diagram_6)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_6
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_6=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE         
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE 
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(h2,p2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(c)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(a)- &
                   all_orbit%evalence(c)-all_orbit%e(p1)
              den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,h1,c,p1,lambda,ans1)
              CALL cross_coupled_mtxel1(p1,p2,h1,h2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3* &
                   all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                   all_orbit%jj(p2))/2) &
                   *bare_operator(h2,p2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(a)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_6(l)=effoper_diagram_6(l)+& 
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_6
!
!
!
subroutine effective_operator_7(a,c,effoper_diagram_7)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_7
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_7=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE         
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE            
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE 
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(p2,h2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p2)+all_orbit%nshell(c)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h2)+all_orbit%evalence(a)- &
                   all_orbit%evalence(c)-all_orbit%e(p2)
              den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,p1,c,h1,lambda,ans1)
              CALL cross_coupled_mtxel1(h1,h2,p1,p2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3*all_orbit%jj(h2)+&
                   all_orbit%jj(p1)+all_orbit%jj(p2))/2) &
                   *bare_operator(p2,h2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(a)+ &
                   all_orbit%e(h2)-all_orbit%evalence(c)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_7(l)=effoper_diagram_7(l)+ &
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_7
!
!
!
SUBROUTINE effective_operator_8(a,c,effoper_diagram_8)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_8
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_8=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE 
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE   
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(h2,p2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p2)+all_orbit%nshell(a)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h2)+all_orbit%evalence(c)- &
                   all_orbit%evalence(a)-all_orbit%e(p2)
              den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,h1,c,p1,lambda,ans1)
              CALL cross_coupled_mtxel1(p1,p2,h1,h2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3*all_orbit%jj(h2)+&
                   all_orbit%jj(p1)+all_orbit%jj(p2))/2) &
                   *bare_operator(h2,p2)/(2.*lambda+1.)
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w1=all_orbit%evalence(c)-all_orbit%evalence(a)+w2
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_8(l)=effoper_diagram_8(l)+ &
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_8
!
!
!
SUBROUTINE effective_operator_9(a,c,effoper_diagram_9)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, p1, l, p2, p3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_9
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_9=0.
  DO p3=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE 
     DO p2=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE        
        IF(bare_operator(p3,p2) == 0.) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p3)
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(h)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h)+all_orbit%evalence(c)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%e(h)+all_orbit%evalence(a)- &
                   all_orbit%e(p1)-all_orbit%e(p3)
              w1=all_orbit%e(h)+all_orbit%evalence(a)+wcn
              w2=all_orbit%e(h)+all_orbit%evalence(c)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p3))/2
              j2min=ABS(all_orbit%jj(h)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(p3))/2
              j2max=(all_orbit%jj(h)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p3), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL cross_coupled_mtxel2(a,h,p1,p3,jt,ans1)
                 CALL cross_coupled_mtxel2(p1,p2,c,h,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=iph((3*all_orbit%jj(p3)+3*all_orbit%jj(h)&
                      +all_orbit%jj(p2)+2*all_orbit%jj(a) &
                      +all_orbit%jj(p1))/2)*sixjj*iph(lambda-jt)&
                      *bare_operator(p3,p2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_9(l)=effoper_diagram_9(l)+&
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_9
!
!
!
SUBROUTINE effective_operator_10(a,c,effoper_diagram_10)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_10
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_10=0.
  DO h1=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE         
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        IF(bare_operator(h1,h2) == 0.) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO p2=1, all_orbit%total_orbits                  
              IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(h2)
              idiff1=nshell1-nshell2
              idiff2=nshell3-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%e(h2)+all_orbit%evalence(a)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(h1))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h2)+all_orbit%evalence(a)+wcn
              w2=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(h1),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL pphhmtx(a,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-0.5D0*iph((all_orbit%jj(h2)+&
                      all_orbit%jj(a))/2+jt)*sixjj*(2.*jt+1.) &
                      *bare_operator(h1,h2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_10(l)=effoper_diagram_10(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_10
!
!
!
SUBROUTINE effective_operator_11(a,c,effoper_diagram_11)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p, h1, l, h2, h3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_11
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_11=0.
  DO h3=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        IF(bare_operator(h2,h3) == 0.) CYCLE
        DO p=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h3)
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(p)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p)-all_orbit%evalence(a)
              den2=all_orbit%e(h1)+all_orbit%e(h3)- &
                   all_orbit%e(p)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(h3))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h3))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h3)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local      
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(h3),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL cross_coupled_mtxel2(a,p,h1,h2,jt,ans1)
                 CALL cross_coupled_mtxel2(h1,h3,c,p,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=iph((3*all_orbit%jj(h1)+3* &
                      all_orbit%jj(h2)+all_orbit%jj(p)&
                      +2*all_orbit%jj(a)+&
                      all_orbit%jj(h3))/2)*sixjj*iph(lambda-jt) &
                      *bare_operator(h2,h3)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_11(l)=effoper_diagram_11(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_11
!
!
!
SUBROUTINE effective_operator_12(a,c,effoper_diagram_12)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_12
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_12=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     DO p2=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
        IF(bare_operator(p2,p1) == 0.0) CYCLE
        DO h1=1,all_orbit%total_orbits
           IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell3-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(c)-all_orbit%e(p2)
              den2=all_orbit%e(h2)+all_orbit%e(h1)- &
                   all_orbit%e(p1)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p1))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h2)+all_orbit%e(h1)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p1), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL pphhmtx(a,p1,h1,h2,jt,ans1)
                 CALL pphhmtx(h1,h2,c,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=iph((all_orbit%jj(p1)+all_orbit%jj(a))/2+jt)
                 factr=-0.5D0*factr*(2.*jt+1.)*sixjj*bare_operator(p2,p1)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_12(l)=effoper_diagram_12(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_12
!
!
!
SUBROUTINE effective_operator_13(a,c,effoper_diagram_13)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, p1, l, p2, p3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_13
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_13=0.
  DO h=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE	 
     DO p3=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE	   
        IF(bare_operator(h,p3) == 0.) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) /= 'particle') CYCLE
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(p3)
              idiff1=nshell1-nshell2
              idiff2=nshell1-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h)+all_orbit%evalence(c)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%evalence(c)+all_orbit%e(h)- &
                   all_orbit%e(p3)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p3))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(h))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(p3))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%evalence(c)+all_orbit%e(h)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p3), &
                      all_orbit%jj(h),2*jt)
                 IF(sixjj == 0.) CYCLE
                 factr=0.5D0*iph(jt+(all_orbit%jj(a)+ &
                      all_orbit%jj(p3))/2)*sixjj &
                      *(2.*jt+1.)*bare_operator(h,p3)
                 CALL pphhmtx(a,p3,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_13(l)=effoper_diagram_13(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_13
!
!
!
SUBROUTINE effective_operator_14(a,c,effoper_diagram_14)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, p1, l, p2, p3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_14
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_14=0.
  DO h=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
     DO p3=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE
        IF(bare_operator(p3,h) == 0.) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE 
           DO p2=1, all_orbit%total_orbits
              IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE 
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(p3)
              idiff1=nshell1-nshell2
              idiff2=nshell1-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h)+all_orbit%evalence(a)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%evalence(a)+all_orbit%e(h)- &
                   all_orbit%e(p3)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p3))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(a)-all_orbit%jj(h))/2
              j1max=(all_orbit%jj(c)+all_orbit%jj(p3))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%evalence(a)+all_orbit%e(h)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h), &
                      all_orbit%jj(p3),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL pphhmtx(a,h,p1,p2,jt,ans1)
                 IF(ans1(2) == 0.) CYCLE
                 CALL pphhmtx(p1,p2,c,p3,jt,ans2)
                 IF(ans2(2) == 0.) CYCLE
                 factr=iph(jt+(all_orbit%jj(a)+3*all_orbit%jj(h)&
                      +2*all_orbit%jj(p3))/2)*sixjj &
                      *0.5*(2.*jt+1.)*bare_operator(p3,h)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_14(l)=effoper_diagram_14(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_14
!
!
!
SUBROUTINE effective_operator_15(a,c,effoper_diagram_15)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p, h1, l, h2, h3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_15
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_15=0.
  DO p=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
     DO h3=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
        IF(bare_operator(h3,p) == 0.0D0) CYCLE
        DO h1=1,all_orbit%total_orbits
           IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(h3)
              idiff1=nshell2-nshell1
              idiff2=nshell3-nshell1
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(a)-all_orbit%e(p)
              den2=all_orbit%evalence(c)-all_orbit%e(p)+ &
                   all_orbit%e(h3)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h3))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h1))/2
              j3min=ABS(all_orbit%jj(a)-all_orbit%jj(p))/2
              j1max=(all_orbit%jj(c)+all_orbit%jj(h3))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h2)+all_orbit%e(h1)+wcn
              w2=w1+all_orbit%evalence(c)+all_orbit%e(h3)- &
                   all_orbit%e(p)-all_orbit%evalence(a)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p), &
                      all_orbit%jj(h3),2*jt)
                 if(sixjj == 0.0) cycle
                 CALL pphhmtx(a,p,h1,h2,jt,ans1)
                 CALL pphhmtx(h1,h2,c,h3,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=0.5*iph(jt+(all_orbit%jj(a)+ &
                      all_orbit%jj(p))/2)*sixjj*(2.*jt+1.)*bare_operator(h3,p)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_15(l)=effoper_diagram_15(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_15
!
!
!
SUBROUTINE effective_operator_16(a,c,effoper_diagram_16)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, p, l, h2, h3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_16
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_16=0.
  DO p=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
     DO h3=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
        IF(bare_operator(p,h3) == 0.) CYCLE
        DO h1=1,all_orbit%total_orbits
           IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(h3)
              idiff1=nshell2-nshell1
              idiff2=nshell3-nshell1
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(c)-all_orbit%e(p)
              den2=all_orbit%evalence(a)-all_orbit%e(p)+ &
                   all_orbit%e(h3)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h3))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(p))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(h3))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%evalence(a)+all_orbit%e(h2)+ &
                   all_orbit%e(h3)+all_orbit%e(h1)- &
                   all_orbit%e(p)-all_orbit%evalence(c)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h3), &
                      all_orbit%jj(p),2*jt)
                 IF(sixjj == 0.) CYCLE
                 factr=iph(jt+(2*all_orbit%jj(p)+ &
                      3*all_orbit%jj(h3)+all_orbit%jj(a))/2)*sixjj &
                      *0.5*(2.*jt+1.)*bare_operator(p,h3)
                 CALL pphhmtx(a,h3,h1,h2,jt,ans1)
                 CALL pphhmtx(h1,h2,c,p,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_16(l)=effoper_diagram_16(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_16
!
!
!
SUBROUTINE effective_operator_17(a,c,effoper_diagram_17)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, p2, l, h2, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_17
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_17=0.
  DO p2=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        IF(bare_operator(h2,p2) == 0.0D0) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(p2)
              nshell4=all_orbit%nshell(c)+all_orbit%nshell(h2)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%evalence(a)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)- &
                   all_orbit%e(p2)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1+all_orbit%evalence(c)-all_orbit%evalence(a)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.0) CYCLE
                 CALL cross_coupled_mtxel2(a,p1,h1,h2,jt,ans1)
                 CALL cross_coupled_mtxel2(h1,p2,c,p1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-iph(lambda-jt+(all_orbit%jj(p1)+ &
                      3*all_orbit%jj(h1)+2*all_orbit%jj(a) &
                      +all_orbit%jj(p2)+3*all_orbit%jj(h2))/2)*sixjj &
                      *bare_operator(h2,p2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_17(l)=effoper_diagram_17(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_17
!
!
!
SUBROUTINE effective_operator_18(a,c,effoper_diagram_18)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_18
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_18=0.
  DO p2=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE      
        IF(bare_operator(p2,h2) == 0.0D0) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(p2)
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(h2)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%evalence(c)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)- &
                   all_orbit%e(p2)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(a)-all_orbit%jj(p2))/2
              j1max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w1=w2+all_orbit%evalence(a)-all_orbit%evalence(c)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(a),all_orbit%jj(c), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL cross_coupled_mtxel2(a,p1,h1,p2,jt,ans1)
                 CALL cross_coupled_mtxel2(h1,h2,c,p1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-iph(-lambda-jt+(all_orbit%jj(p1)+ &
                      3*all_orbit%jj(h1)+2*all_orbit%jj(a) &
                      +2*all_orbit%jj(c)+3*all_orbit%jj(p2)+&
                      3*all_orbit%jj(h2))/2)*sixjj*bare_operator(p2,h2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_18(l)=effoper_diagram_18(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_18
!
!
!
SUBROUTINE effective_operator_19(a,c,effoper_diagram_19)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_19
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_19=0.
  DO p2=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE            
        IF(bare_operator(h2,p2) == 0.0D0) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(h2)
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)- &
                   all_orbit%e(p2)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                   all_orbit%evalence(c)-all_orbit%e(p2)+wcn
              w2=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL cross_coupled_mtxel2(a,h1,p1,h2,jt,ans1)
                 CALL cross_coupled_mtxel2(p1,p2,c,h1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-iph(lambda-jt+(all_orbit%jj(p1)+ &
                      3*all_orbit%jj(h1)+2*all_orbit%jj(a) &
                      +all_orbit%jj(p2)+3*all_orbit%jj(h2))/2)*sixjj &
                      *bare_operator(h2,p2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_19(l)=effoper_diagram_19(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_19
!
!
!
SUBROUTINE effective_operator_20(a,c,effoper_diagram_20)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_20
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_20=0.
  DO p2=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        IF(bare_operator(p2,h2) == 0.0D0) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(h2)
              nshell4=all_orbit%nshell(c)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(a)-& 
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)- &
                   all_orbit%e(p2)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(JMIN_LOCAL > JMAX_LOCAL) CYCLE
              w2=all_orbit%e(h1)+all_orbit%e(h2)+ &
                   all_orbit%evalence(a)-all_orbit%e(p2)+wcn
              w1=all_orbit%e(h1)+all_orbit%evalence(a)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p2), &
                      all_orbit%jj(h2),2*jt)
                 IF(sixjj == 0.0D0) CYCLE
                 CALL cross_coupled_mtxel2(a,h1,p1,p2,jt,ans1)
                 CALL cross_coupled_mtxel2(p1,h2,c,h1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-iph(-lambda-jt+(all_orbit%jj(p1)+ &
                      3*all_orbit%jj(h1)+2*all_orbit%jj(a) &
                      +3*all_orbit%jj(p2)+5*all_orbit%jj(h2))/2)*sixjj &
                      *bare_operator(p2,h2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_20(l)=effoper_diagram_20(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_20
!
!
!
SUBROUTINE effective_operator_21(a,c,effoper_diagram_21)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_21
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_21=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(a,p1) == 0.) CYCLE
     fact=-0.5/(all_orbit%jj(c)+1.)*bare_operator(a,p1)* &
          iph((all_orbit%jj(c)-all_orbit%jj(p1))/2)
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE	    
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE  
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,c,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_21(l)=effoper_diagram_21(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_21
!
!
!
SUBROUTINE effective_operator_22(a,c,effoper_diagram_22)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_22
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_22=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(p1,c) == 0.) CYCLE
     fact=-bare_operator(p1,c)*0.5/(all_orbit%jj(a)+1.) * &
          iph((all_orbit%jj(p1)-all_orbit%jj(a)))
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p2)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(a,p2,h1,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_22(l)=effoper_diagram_22(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_22
!
!
!
SUBROUTINE effective_operator_23(a,c,effoper_diagram_23)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_23
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_23=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=-bare_operator(h1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((-all_orbit%jj(a)+all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell3
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%evalence(a)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans1)
                 CALL pphhmtx(a,h2,p1,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_23(l)=effoper_diagram_23(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_23
!
!
!
SUBROUTINE effective_operator_24(a,c,effoper_diagram_24)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_24
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_24=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=-bare_operator(a,h1)*0.5/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell3
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_24(l)=effoper_diagram_24(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_24
!
!
!
SUBROUTINE effective_operator_23folded(a,c,effoper_diagram_23f)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_23f
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_23f=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%model_space(h1) == 'outside') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=-bare_operator(h1,c)*0.25/(all_orbit%jj(a)+1.)* &
          iph((-all_orbit%jj(a)+all_orbit%jj(h1)))
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell3
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%evalence(a)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans1)
                 CALL pphhmtx(a,h2,p1,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_23f(l)=effoper_diagram_23f(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_23folded
!
!
!
SUBROUTINE effective_operator_24folded(a,c,effoper_diagram_24f)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_24f
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_24f=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%model_space(h1) == 'outside') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=-bare_operator(a,h1)*0.25/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell3
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_24f(l)=effoper_diagram_24f(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_24folded

!
!
!
SUBROUTINE effective_operator_25(a,c,effoper_diagram_25)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, p2, l, p3, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_25
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_25=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(p1,c) == 0.) CYCLE
     fact=bare_operator(p1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(p1)-all_orbit%jj(a)))
     DO p2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
        DO p3=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p3)
              nshell3=all_orbit%nshell(h1)+all_orbit%nshell(a)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(a)-all_orbit%e(p3)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(a)-all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h1))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p3))/2
              j3min=ABS(all_orbit%jj(p1)-all_orbit%jj(h1))/2
              j1max=(all_orbit%jj(p1)+all_orbit%jj(h1))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p3))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p2,p3,p1,h1,jt,ans1)
                 CALL pphhmtx(a,h1,p2,p3,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_25(l)=effoper_diagram_25(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_25
!
!
!
SUBROUTINE effective_operator_26(a,c,effoper_diagram_26)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, p2, l, p3, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_26
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_26=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(a,p1) == 0.) CYCLE
     fact=bare_operator(a,p1)*0.5/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(p1))/2)
     DO p2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
        DO p3=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p3)
              nshell3=all_orbit%nshell(h1)+all_orbit%nshell(c)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)-all_orbit%e(p3)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(c)-all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h1))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p3))/2
              j3min=ABS(all_orbit%jj(p1)-all_orbit%jj(h1))/2
              j1max=(all_orbit%jj(p1)+all_orbit%jj(h1))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p3))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)  		  
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,h1,p2,p3,jt,ans1)
                 CALL pphhmtx(p2,p3,c,h1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_26(l)=effoper_diagram_26(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_26
!
!
!
SUBROUTINE effective_operator_27(a,c,effoper_diagram_27)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_27
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_27=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=-bare_operator(h1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(h1)-all_orbit%jj(a))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans1)
                 CALL pphhmtx(a,h2,p1,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_27(l)=effoper_diagram_27(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_27
!
!
!
SUBROUTINE effective_operator_28(a,c,effoper_diagram_28)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_28
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_28=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=-bare_operator(a,h1)*0.5/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_28(l)=effoper_diagram_28(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_28
!
!
!
SUBROUTINE effective_operator_27folded(a,c,effoper_diagram_27f)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_27f
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_27f=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%model_space(h1) == 'outside') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=-bare_operator(h1,c)*0.25/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(h1)-all_orbit%jj(a))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans1)
                 CALL pphhmtx(a,h2,p1,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_27f(l)=effoper_diagram_27f(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_27folded
!
!
!
SUBROUTINE effective_operator_28folded(a,c,effoper_diagram_28f)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_28f
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_28f=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%model_space(h1) == 'outside') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=-bare_operator(a,h1)*0.25/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_28f(l)=effoper_diagram_28f(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_28folded
!
!
!
SUBROUTINE effective_operator_29(a,c,effoper_diagram_29)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_29
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_29=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(a,p1) == 0.0D0) CYCLE
     fact=-0.5/(all_orbit%jj(c)+1.)*bare_operator(a,p1)* &
          iph((all_orbit%jj(c)-all_orbit%jj(p1))/2)
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE	    
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE  
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(c)-all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)&
                   -all_orbit%e(p1)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,c,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_29(l)=effoper_diagram_29(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_29
!
!
!
SUBROUTINE effective_operator_30(a,c,effoper_diagram_30)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_30
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_30=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(p1,c) == 0.0D0) CYCLE
     fact=-bare_operator(p1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(a)-all_orbit%jj(p1)))
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(a)-all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1+all_orbit%evalence(a)-all_orbit%e(p1)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(a,p2,h1,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_30(l)=effoper_diagram_30(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_30
!
!
!
SUBROUTINE effective_operator_31(a,c,effoper_diagram_31)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, h2, l, h3, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_31
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_31=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=bare_operator(h1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(a)-all_orbit%jj(h1))/2)
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        DO h3=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
           DO p1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h3)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(a)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=-(all_orbit%e(p1)+all_orbit%evalence(a)-all_orbit%e(h3)- &
                   all_orbit%e(h2))
              den2=all_orbit%e(h1)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p1))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h3))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h3))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              w1=w2+all_orbit%e(h1)-all_orbit%evalence(a)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h2,h3,h1,p1,jt,ans1)
                 CALL pphhmtx(a,p1,h2,h3,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_31(l)=effoper_diagram_31(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_31
!
!
!
SUBROUTINE effective_operator_32(a,c,effoper_diagram_32)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, h2, l, h3, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_32
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_32=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=bare_operator(a,h1)*0.5/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        DO h3=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
           DO p1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h3)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(c)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=-(all_orbit%e(p1)+all_orbit%evalence(c)-all_orbit%e(h3)- &
                   all_orbit%e(h2))
              den2=all_orbit%e(h1)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p1))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h3))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h3))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h2)+all_orbit%e(h3)+wcn+ &
                   all_orbit%e(h1)-all_orbit%evalence(c)
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,p1,h2,h3,jt,ans1)
                 CALL pphhmtx(h2,h3,c,p1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_32(l)=effoper_diagram_32(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_32

