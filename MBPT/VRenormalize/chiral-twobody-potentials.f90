!
!     C H I R A L    T W O B O D Y    P O T E N T I A L S
!
!**********************************************************************
!
! PROGRAM PACKAGE TO COMPUTE ALL POSSIBLE CHIRAL TWOBODY POTENTIALS
! AT LEADING-ORDER (LO), NEXT-TO LEADING-ORDER (NLO)
! AND NEXT-TO-NEXT-TO LEADING-ORDER (NNLO)
!
! THIS CODE IS BASED ON THE WORK AND HELP OF RUPRECHT MACHLEIDT et al.
! AT MOSCOW UNIVERSITY, IDAHO.
!
! AUTHOR: Andreas Ekström
! ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
! EMAIL: jaeks@fys.uio.no
! LAST UPDATE: 2013-08-16
!
! N3LO-addition:
!   AUTHOR: Boris D. Carlsson
!   EMAIL: borisc@chalmers.se
!
! IMPORTANT: COMPILE WITH INTEL FORTRAN
! OTHER COMPILERS, SUCH AS PORTLAND GIVE
! WRONG RESULTS. WILL CORRECT THIS IN THE
! FUTURE.
!
!**********************************************************************
MODULE chp_aux
  
  IMPLICIT NONE
  
  ! basic mesh info
  TYPE, PUBLIC :: chp_mesh_info
     INTEGER  :: amount
     REAL(8) :: xmin, xmax
  END TYPE chp_mesh_info
  
  ! GAUSS-LEGENDRE MESH POINTS
  ! Gauss-Legendre mesh point x, corresponding integration weight w and corresponding x*x*w-value
  TYPE, PUBLIC :: chp_gauleg_mesh_point
     SEQUENCE
     REAL(8) :: x, w, xxw
  END TYPE chp_gauleg_mesh_point
  
  ! mesh points and weights in momentum space
  TYPE, PUBLIC :: chp_gauleg_mesh
     TYPE(chp_mesh_info)                                    :: info
     TYPE(chp_gauleg_mesh_point), DIMENSION(:), ALLOCATABLE :: pnt
  END TYPE chp_gauleg_mesh
  
  TYPE, PUBLIC :: chp_real_type
     REAL(8)           :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_real_type
  
  TYPE, PUBLIC :: chp_char2_type
     CHARACTER(LEN=2) :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_char2_type
  
  TYPE, PUBLIC :: chp_int_type
     INTEGER          :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_int_type

  TYPE, PUBLIC :: chp_chn_type
     INTEGER          :: S,j,T,tz
     LOGICAL          :: coup
  END TYPE chp_chn_type
  
  ! mathematical constants
  TYPE(chp_real_type), PARAMETER :: chp_pi=chp_real_type(ACOS(-1.0D0),'pi',.TRUE.)
  REAL(8)             , PARAMETER :: fourpi=4.0D0*chp_pi%val
  REAL(8)             , PARAMETER :: twopi =2.0D0*chp_pi%val
  REAL(8)             , PARAMETER :: pi2   = chp_pi%val*chp_pi%val
  REAL(8)             , PARAMETER :: pi_inv= 1D0 / chp_pi%val
  
  !--- standard formatted contacts to pw contacts conversion matrix
  REAL(8), PUBLIC :: LOST2PW(2,2)
  REAL(8), PUBLIC :: LOPW2ST(2,2)
  REAL(8), PUBLIC :: NLOST2PW(7,7)
  REAL(8), PUBLIC :: NLOPW2ST(7,7)
  REAL(8), PUBLIC :: N3LOST2PW(15,15)
  REAL(8), PUBLIC :: N3LOPW2ST(15,15)

CONTAINS
  !
  !  This function sets up the recursive relation
  !  for the associated Legendre polynomials
  !
  REAL(8) FUNCTION chp_legendre_polynomials(l, m, x)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: l, m
    REAL(8), INTENT(IN) :: x
    REAL(8)             :: fact,pll,pmm,pmmp1,somx2
    INTEGER              :: i,ll
    

    !  check whether m, l and x are ok
    IF((M < 0).OR.(M > L).OR.(ABS(X) > 1.0D0)) THEN
       WRITE(*,*) 'legendre_polynomials: bad arguments', m, l, x
       chp_legendre_polynomials = 0.0D0
       RETURN
    ENDIF

    !  calculate now pmm as starting point for iterations
    pmm=1.0
    IF (m > 0) THEN
       somx2=SQRT((1.0-x)*(1.0+x))
       fact=1.0D0;
       DO i=1, m
          pmm = -fact*somx2*pmm
          fact = fact+2.0D0
       ENDDO
    ENDIF

    !  if l == m we do not need to use recursion relation
    IF (l == m) THEN
       chp_legendre_polynomials=pmm

       !  recursive relation for associated Legendre polynomials
    ELSE
       pmmp1=x*(2*m+1)*pmm

       !  analytical formula for the case l == m+1
       IF (l == (m+1)) THEN
          chp_legendre_polynomials=pmmp1
       ELSE
          DO ll=m+2, l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          chp_legendre_polynomials= pll
       ENDIF
    ENDIF
    
  END FUNCTION chp_legendre_polynomials
  
  !
  !      This routine calculates gauss-legendre mesh points and weights      
  !      INPUT:                                                              
  !      mesh%info%xmin     : lower limit of the integration interval                      
  !      mesh%info%xmax     : upper limit ---------- "" -------------                      
  !      mesh%info%amount   : the desired number of mesh points                            
  !      OUTPUT:                                                            
  !      mesh%pnt(:)%x      : gauss-legendre mesh points on the interval (x1,x2)          
  !      mesh%pnt(:)%w      : the corresponding weights                                   
  !      FROM               : Numerical recipes
  !      F90 version        : M. Hjorth-Jensen
  !      Object interface   : M. Kartamyshev
  !
  SUBROUTINE chp_setup_gauleg_mesh (mesh)


    TYPE(chp_gauleg_mesh), INTENT(INOUT)        :: mesh
    INTEGER                                 :: i, j, m, n
    REAL(8)                                :: x1, x2
    REAL(8), DIMENSION(:), ALLOCATABLE     :: x, w
    REAL(8)                                :: p1,p2,p3,pp,xl,xm,z,z1
    REAL(8), PARAMETER                     :: EPS = 3.D-14
    
    ALLOCATE(x(mesh%info%amount))
    ALLOCATE(w(mesh%info%amount))
        
    ! allocate points and weights storages
    CALL chp_destroy_gauleg_mesh (mesh)
    
    IF (mesh%info%amount <= 0 .OR. mesh%info%xmin >= mesh%info%xmax) THEN
       WRITE(*,*) ': incorrect mesh info', mesh%info ; STOP
    ENDIF
    
    ALLOCATE( mesh%pnt( 1:mesh%info%amount ) )
    mesh%pnt(:) = chp_gauleg_mesh_point(0.0D0, 0.0D0, 0.0D0)
    
    ! set values of local variables
    x1 = mesh%info%xmin ; x2 = mesh%info%xmax; n = mesh%info%amount
    
    m=(n+1)/2
    xm=0.5D0*(x2+x1)
    xl=0.5D0*(x2-x1)
    DO i=1,m
       z1=0.0D0
       z=COS(chp_pi%val*(i - 0.25D0)/(n + 0.5D0))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.0D0
          p2=0.0D0
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.0D0*j-1.0D0)*z*p2-(j-1.0D0)*p3)/j
          END DO
          pp=n*(z*p1-p2)/(z*z-1.0D0)
          z1=z
          z=z-p1/pp
       END DO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0D0*xl/((1.0D0-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO
    
    ! set return values
    mesh%pnt(:)%x = x(:) ; mesh%pnt(:)%w = w(:) ; mesh%pnt(:)%xxw = x(:) * x(:) * w(:)

    DEALLOCATE(w)
    DEALLOCATE(x)

  END SUBROUTINE chp_setup_gauleg_mesh
  
  SUBROUTINE chp_destroy_gauleg_mesh (mesh)
    TYPE(chp_gauleg_mesh), INTENT(INOUT) :: mesh
    IF (ALLOCATED(mesh%pnt) ) DEALLOCATE (mesh%pnt)
  END SUBROUTINE chp_destroy_gauleg_mesh

  ! phase factor of type (-1)**arg
  FUNCTION chp_minus_power (arg) RESULT (res)
    INTEGER              :: res
    INTEGER, INTENT (IN) :: arg
    INTEGER :: exponent
    
    exponent = ABS(arg) ! (-1)**N = (-1)**(-N)
    
    SELECT CASE (MOD(exponent,2))
    CASE (0)
       res =  1
    CASE (1)
       res = -1
    END SELECT

  END FUNCTION chp_minus_power

  !
  !     Given an NxN matrix A(N,N), this routine replaces it by the LU
  !     decomposed one, where the matrix elements are stored in the same
  !     matrix A. The array indx is  an output vector which records the row
  !     permutation effected by the partial pivoting. d is the determinant
  !
  SUBROUTINE chp_lu_decompose(a,n,indx,d)
    
    IMPLICIT NONE
    INTEGER :: n, i, j, k, imax
    REAL(8) :: sum , tiny, aamax, dum, d
    REAL(8), DIMENSION(n,n) :: a
    INTEGER, DIMENSION(n) :: indx
    REAL(8), ALLOCATABLE :: vv(:)
    
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
    
  END SUBROUTINE chp_lu_decompose
  
  !     Solves set of linear equations Ax=b, A is input as an LU decompomsed
  !     matrix and indx keeps track of the permutations of the rows. b is input
  !     as the right-hand side vector b and returns the solution x. A, n and indx
  !     are not modified by this routine. This function takes into that b can contain
  !     many zeros and is therefore suitable for matrix inversion
  
  
  SUBROUTINE chp_lu_linear_equation(a,n,indx,b)
    
    IMPLICIT NONE
    INTEGER :: n, ii, ll, i, j
    REAL(8) :: sum
    REAL(8), DIMENSION(n,n) :: a
    REAL(8), DIMENSION(n) :: b
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
    
  END SUBROUTINE chp_lu_linear_equation
  
   !            Routines to do mtx inversion, from Numerical
  !            Recepies, Teukolsky et al. Routines included
  !            below are MATINV, LUDCMP and LUBKSB. See chap 2
  !            of Numerical Recipes for further details
  !            Recoded in FORTRAN 90 by M. Hjorth-Jensen
  !
  SUBROUTINE chp_matinv(a,n)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j
    REAL(8), DIMENSION(n,n), INTENT(INOUT)  :: a
    REAL(8), ALLOCATABLE :: y(:,:)
    REAL(8) :: d
    INTEGER, ALLOCATABLE :: indx(:)
    
    ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
    y=0.
    !     setup identity matrix
    DO i=1,n
       y(i,i)=1.
    ENDDO
    !     LU decompose the matrix just once
    CALL  chp_lu_decompose(a,n,indx,d)
    
    !     Find inverse by columns
    DO j=1,n
       CALL chp_lu_linear_equation(a,n,indx,y(:,j))
    ENDDO
    !     The original matrix a was destroyed, now we equate it with the inverse y
    a=y
    
    DEALLOCATE ( y ); DEALLOCATE ( indx )
    
  END SUBROUTINE chp_matinv
  
  SUBROUTINE set_contact_conversion_matrices
    
    LOST2PW = transpose(reshape((/ 1.0D0, -3.0D0,  &
         1.0D0,  1.0D0 /), shape(LOST2PW)))
    
    LOST2PW = fourpi * LOST2PW

    LOPW2ST = LOST2PW
    
    CALL chp_matinv(LOPW2ST,2)

    NLOST2PW = transpose(reshape((/ 1.0D0, 0.25D0, -3.0D0, -0.75D0, 0.0D0, -1.0D0, -0.25D0, &
                                    -2.0D0/3.0D0, 1.0D0/6.0D0, -2.0D0/3.0D0, 1.0D0/6.0D0, -2.0D0/3.0D0, 2.0D0, - 0.5D0, &
                                    -2.0D0/3.0D0, 1.0D0/6.0D0, 2.0D0, -0.5D0, 0.0D0, 2.0D0/3.0D0, -1.0D0/6.0D0, &
                                    -2.0D0/3.0D0, 1.0D0/6.0D0, -2.0D0/3.0D0, 1.0D0/6.0D0, -1.0D0/3.0D0, -4.0D0/3.0D0, 1.0D0/3.0D0, &
                                    1.0D0, 0.25D0, 1.0D0, 0.25D0, 0.0D0, 1.0D0/3.0D0, 1.0D0/12.0D0, &
                                    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -2.0D0*DSQRT(2.0D0)/3.0D0, -DSQRT(2.0D0)/6.0D0, &
                                    -2.0D0/3.0D0, 1.0D0/6.0D0, -2.0D0/3.0D0, 1.0D0/6.0D0, +1.0D0/3.0D0, 0.0D0, 0.0D0 /), shape(NLOST2PW)))
    NLOST2PW = fourpi * NLOST2PW

    NLOPW2ST = NLOST2PW
    
    CALL chp_matinv(NLOPW2ST,7)
    
    N3LOST2PW = transpose(reshape((/ &
                    1D0, 1D0/16D0, 0.25D0, 0D0, -3D0, -3D0/16D0, -0.75D0, 0D0, 0D0, 0D0, -1D0, -0.25D0, -0.25D0, -1D0/16D0, 0D0, &
                    10D0/3D0, 5D0/24D0, 1D0/6D0, 2D0/3D0, -10D0, -0.625D0, -0.5D0, -2D0, 0D0, 0D0, -10D0/3D0, -1D0/6D0, -1D0/6D0, -5D0/24D0, -2D0/3D0, &
                    -4D0/3D0, 1D0/12D0, 0D0, 0D0, -4D0/3D0, 1D0/12D0, 0D0, 0D0, -2D0/3D0, -1D0/6D0, 8D0/3D0, 1D0/3D0, -1D0/3D0, -1D0/6D0, 0D0, &
                    -4D0/3D0, 1D0/12D0, 0D0, 0D0, 4D0, -0.25D0, 0D0, 0D0, 0D0, 0D0, 4D0/3D0, 0D0, 0D0, -1D0/12D0, 0D0, &
                    -4D0/3D0, 1D0/12D0, 0D0, 0D0, -4D0/3D0, 1D0/12D0, 0D0, 0D0, -1D0/3D0, -1D0/12D0, -2D0, -1D0/6D0, 1D0/6D0, 1D0/8D0, 0D0, &
                    1D0, 1D0/16D0, 0.25D0, 0D0, 1D0, 1D0/16D0, 0.25D0, 0D0, 0D0, 0D0, 1D0/3D0, 1D0/12D0, 1D0/12D0, 1D0/48D0, 0D0, &
                    10D0/3D0, 5D0/24D0, 1D0/6D0, 2D0/3D0, 10D0/3D0, 5D0/24D0, 1D0/6D0, 2D0/3D0, 0D0, 0D0, 10D0/9D0, 1D0/18D0, 1D0/18D0, 5D0/72D0, 2D0/9D0, &
                    8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 2D0/5D0, -0.1D0, -4D0/9D0, 1D0/9D0, 1D0/9D0, -1D0/36D0, -16D0/45D0, &
                    0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, -DSQRT(2D0)*2D0/3D0, -DSQRT(2D0)/6D0, -DSQRT(2D0)/6D0, -DSQRT(2D0)/24D0, 0D0, &
                    0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, -DSQRT(2D0)*14D0/9D0, DSQRT(2D0)/18D0, DSQRT(2D0)/18D0, -DSQRT(2D0)*7D0/72D0, DSQRT(2D0)*2D0/9D0, &
                    8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, -8D0/5D0, -0.1D0, 0.4D0, 0.4D0, 0D0, 0D0, -8D0/15D0, 2D0/15D0, 2D0/15D0, -1D0/30D0, 2D0/15D0, &
                    8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 2D0/15D0, -1D0/30D0, 0.8D0, -0.2D0, -0.2D0, 0.05D0, 4D0/15D0, &
                    -4D0/3D0, 1D0/12D0, 0D0, 0D0, -4D0/3D0, 1D0/12D0, 0D0, 0D0, 1D0/3D0, 1D0/12D0, -2D0/15D0, 1D0/30D0, -1D0/30D0, 1D0/120D0, 0D0, &
                    0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, DSQRT(6D0)*4D0/15D0, -DSQRT(6D0)/15D0, DSQRT(6D0)/15D0, -DSQRT(6D0)/60D0, 0D0, &
                    8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, -4D0/15D0, 1D0/15D0, 0D0, 0D0, 0D0, 0D0, -2D0/15D0 &
                /), shape(N3LOST2PW)))
    N3LOST2PW = fourpi * N3LOST2PW
    N3LOPW2ST = N3LOST2PW
    CALL chp_matinv(N3LOPW2ST, 15)

  END SUBROUTINE set_contact_conversion_matrices
     
END MODULE chp_aux

! Ref: K. Erkelenz et al., NPA176, 413-432 (1971)
MODULE twobody_pwd
  
  USE chp_aux

  IMPLICIT NONE
  
  INTEGER, PRIVATE :: NZ, JMAX

  TYPE(chp_gauleg_mesh) :: zmesh
  !                 z,J
  REAL(8), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: legP
  !            z**l z,l
  REAL(8), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: zl
CONTAINS
  
  SUBROUTINE setup_twobody_pwd(set_NZ, set_JMAX)
    
    INTEGER :: set_NZ, set_JMAX
        
    NZ   = set_NZ
    JMAX = set_JMAX

    CALL setup_zmesh
    CALL setup_legendre_polynomials
    
  END SUBROUTINE setup_twobody_pwd
  
  SUBROUTINE print_pwd_numerics(unit)
    
    INTEGER, INTENT(IN) :: unit

    WRITE(unit,"(A8,I6)") 'JMAX:', JMAX
    WRITE(unit,"(A8,I6)") 'NZ  :', NZ

  END SUBROUTINE print_pwd_numerics

  SUBROUTINE setup_zmesh
    
    zmesh%info = chp_mesh_info(NZ, -1.0D0, +1.0D0)
    CALL chp_setup_gauleg_mesh(zmesh)
    
  END SUBROUTINE setup_zmesh

  SUBROUTINE setup_legendre_polynomials
    
    INTEGER :: J, iz
    REAL(8)  :: z
    REAL(8)  :: val, err

    IF (ALLOCATED(legP)) DEALLOCATE(legP)
    IF (ALLOCATED(zl)) DEALLOCATE(zl)
    ALLOCATE(legP(1:zmesh%info%amount,0:JMAX))
    ALLOCATE(zl(1:zmesh%info%amount,0:JMAX))
    DO J=0, JMAX
       DO iz=1, zmesh%info%amount
          
          z = zmesh%pnt(iz)%x
          legP(iz,J) = chp_legendre_polynomials(J,0,z)
          zl(iz,J) = z**J
       END DO

       val = sum(legP(:,J) * zmesh%pnt(:)%w)
       if(J == 0) then
           err = abs(val - 2)
       else
           err = abs(val)
       end if
       if(err > 1.e-13) then
           write(*,*) 'Legendre polynomial ', J, ' bad, error is ', err
           stop
       end if
    END DO
    
  END SUBROUTINE setup_legendre_polynomials
  
  !   c     pwd of central force
  ! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_c(W, CHN, pfinal, pinit, POT)
    
    REAL(8)            ,INTENT(INOUT) :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)            ,INTENT(INOUT) :: POT(1:6)
    INTEGER :: j
    
    j = CHN%j
    ! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 2.0D0*pwd_integral(W,0,j)
       
    END IF
    
    IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*pwd_integral(W,0,j)
       
       ! once here, return
       RETURN
    END IF
    
    ! coupled channels
    IF (CHN%coup) THEN
       
       !++Vj
       POT(3) = POT(3) + 2.0D0*pwd_integral(W,0,j+1)

       !--Vj
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*pwd_integral(W,0,j-1)
       
       POT(5) = POT(5) + 0.0D0
       
       POT(6) = POT(6) + 0.0D0

    END IF
    
  END SUBROUTINE pwd_c

  !   c     pwd of central force
  ! ref: E. Epelbaum et al., NPA 747 (2005), 362
  !SUBROUTINE pwd_c_E(W, CHN, pfinal, pinit, POT)
  !  
  !  REAL(8)            ,INTENT(INOUT) :: W(1:NZ)
  !  TYPE(chp_chn_type),INTENT(IN)    :: CHN
  !  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
  !  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
  !  INTEGER :: j
  !  
  !  j = CHN%j
  !  ! uncoupled singlet
  !  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
  !     
  !     POT(1) = POT(1) + 2.0D0*pwd_integral(W,0,j)
  !     
  !  END IF
  !  
  !  IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
  !     
  !     POT(2) = POT(2) + 2.0D0*pwd_integral(W,0,j)
  !     
  !     ! once here, return
  !     RETURN
  !  END IF
  !  
  !  ! coupled channels
  !  IF (CHN%coup) THEN
  !     
  !     POT(3) = POT(3) + 2.0D0*pwd_integral(W,0,j+1)
  !     IF (j==0) RETURN
  !     POT(4) = POT(4) + 2.0D0*pwd_integral(W,0,j-1)
  !
  !     POT(5) = POT(5) + 0.0D0
  !     
  !     POT(6) = POT(6) + 0.0D0
  !     
  !  END IF
  !  
  !END SUBROUTINE pwd_c_E

  !   s     spin-spin
  ! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_s(W, CHN, pfinal, pinit, POT)
    
    REAL(8)            ,INTENT(INOUT) :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)            ,INTENT(INOUT) :: POT(1:6)
    INTEGER :: j
    
    j = CHN%j
    ! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) - 6.0D0*pwd_integral(W,0,j)
       
    ELSE IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*pwd_integral(W,0,j)
       
    ! coupled channels
    ELSE IF (CHN%coup) THEN
       
       POT(3) = POT(3) + 2.0D0*pwd_integral(W,0,j+1)
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*pwd_integral(W,0,j-1)
       
       POT(5) = POT(5) + 0.0D0
       
       POT(6) = POT(6) + 0.0D0

    END IF
    
  END SUBROUTINE pwd_s

  !   s     spin-spin
  ! ref: E. Epelbaum et al., NPA 747 (2005), 362
  !SUBROUTINE pwd_s_E(W, CHN, pfinal, pinit, POT)
  !  
  !  REAL(8)            ,INTENT(INOUT) :: W(1:NZ)
  !  TYPE(chp_chn_type),INTENT(IN)    :: CHN
  !  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
  !  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
  !  INTEGER :: j
  !  
  !  j = CHN%j
  !  ! uncoupled singlet
  !  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
  !     
  !     POT(1) = POT(1) - 6.0D0*pwd_integral(W,0,j)
  !     
  !  ELSE IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
  !     
  !     POT(2) = POT(2) + 2.0D0*pwd_integral(W,0,j)
  !     
  !  ! coupled channels
  !  ELSE IF (CHN%coup) THEN
  !     
  !     POT(3) = POT(3) + 2.0D0*pwd_integral(W,0,j+1)
  !     IF (j==0) RETURN
  !     POT(4) = POT(4) + 2.0D0*pwd_integral(W,0,j-1)
  !     
  !     POT(5) = POT(5) + 0.0D0
  !     
  !     POT(6) = POT(5) + 0.0D0
  !
  !  END IF
  !  
  !END SUBROUTINE pwd_s_E

  !   LS    spin-orbit
  ! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_LS(W, CHN, pfinal, pinit, POT)
    
    REAL(8)            ,INTENT(IN)    :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)            ,INTENT(INOUT) :: POT(1:6)
    REAL(8)  :: jj, p, pp, p2, pp2, rj
    INTEGER :: j
    
    rj = REAL(CHN%j,kind=8)
    jj = 2.0D0*rj+1.0D0
    j = CHN%j
    
    pp  = pfinal; p  = pinit
    pp2 = pp*pp ; p2 = p*p
    ! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 0.0D0
       
    END IF
    
    IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*pp*p*( pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj
       
       ! once here, return
       RETURN
    END IF
    
    ! coupled channels
    IF (CHN%coup) THEN
       
       POT(3) = POT(3) + 2.0D0*pp*p*(rj+2.0D0)*( pwd_integral(W,0,j+2) - pwd_integral(W,0,j))/(2.0D0*rj+3.0D0)
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*pp*p*(rj-1.0D0)*( pwd_integral(W,0,j-2) - pwd_integral(W,0,j))/(2.0D0*rj-1.0D0)
       
       POT(5) = POT(5) + 0.0D0
       
       POT(6) = POT(6) + 0.0D0

    END IF
    
  END SUBROUTINE pwd_LS

  !  LS     spin-orbit
  ! ref: E. Epelbaum et al., NPA 747 (2005), 362
  !SUBROUTINE pwd_LS_E(W, CHN, pfinal, pinit, POT)
  !  
  !  REAL(8)            ,INTENT(IN)    :: W(1:NZ)
  !  TYPE(chp_chn_type),INTENT(IN)    :: CHN
  !  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
  !  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
  !  REAL(8)  :: jj, p, pp, p2, pp2, rj
  !  INTEGER :: j
  !  
  !  rj = REAL(CHN%j,kind=8)
  !  jj = 2.0D0*rj+1.0D0
  !  j = CHN%j
  !  
  !  pp  = pfinal; p  = pinit
  !  pp2 = pp*pp ; p2 = p*p
  !
  !  ! uncoupled singlet
  !  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
  !     
  !     POT(1) = POT(1) + 0.0D0
  !     
  !  END IF
  !  
  !  IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
  !     
  !     POT(2) = POT(2) + 4.0D0*pp*p*pwd_integral(W,1,j) - 2.0D0*pp*p*(pwd_integral(W,0,j-1) + &
  !          pwd_integral(W,0,j+1))
  !     
  !     ! once here, return
  !     RETURN
  !  END IF
  !  
  !  ! coupled channels
  !  IF (CHN%coup) THEN
  !     
  !     POT(3) = POT(3) - 2.0D0*pp*p*pwd_integral(W,0,j)+2.0D0*pp*p*pwd_integral(W,1,j+1)
  !     IF (j==0) RETURN
  !     POT(4) = POT(4) - 2.0D0*pp*p*pwd_integral(W,0,j)+2.0D0*pp*p*pwd_integral(W,1,j-1)
  !     
  !     POT(5) = POT(5) + 0.0D0
  !     
  !     POT(6) = POT(6) + 0.0D0
  !
  !  END IF
  !  
  !END SUBROUTINE pwd_LS_E
  
  !   sigL  sigma-L
  ! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_sigL(W, CHN, pfinal, pinit, POT)
    
    REAL(8)            ,INTENT(IN)    :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)            ,INTENT(INOUT) :: POT(1:6)
    REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
    INTEGER :: j
    
    rj = REAL(CHN%j,kind=8)
    jj = 2.0D0*rj+1.0D0
    jj1 = rj*(rj+1.0D0)
    j = CHN%j
    
    pp  = pfinal; p  = pinit
    pp2 = pp*pp ; p2 = p*p
    ! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 2.0D0*pp2*p2*(pwd_integral(W,2,j) - pwd_integral(W,0,j))
       
    END IF
    
    IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*pp2*p2*( -1.0D0*pwd_integral(W,0,j)  + &
            ((rj-1.0D0)*pwd_integral(W,1,j+1) + (rj+2.0D0)*pwd_integral(W,1,j-1))/jj)
       
       ! once here, return
       RETURN
    END IF
    
    ! coupled channels
    IF (CHN%coup) THEN
       
       POT(3) = POT(3) + 2.0D0*pp2*p2*( -1.0D0*pwd_integral(W,2,j+1)  + &
            ((2.0D0*rj+3.0D0)*pwd_integral(W,0,j+1) - (2.0D0)*pwd_integral(W,1,j))/jj)
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*pp2*p2*( -1.0D0*pwd_integral(W,2,j-1)  + &
            ((2.0D0*rj-1.0D0)*pwd_integral(W,0,j-1) + (2.0D0)*pwd_integral(W,1,j))/jj)

       POT(5) = POT(5) - DSQRT(jj1)*4.0D0*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj**2

       POT(6) = POT(6) - DSQRT(jj1)*4.0D0*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj**2

    END IF
    
  END SUBROUTINE pwd_sigL

  !   sigL  sigma-L
  ! ref: E. Epelbaum et al., NPA 747 (2005), 362
  !SUBROUTINE pwd_sigL_E(W, CHN, pfinal, pinit, POT)
  !  
  !  REAL(8)            ,INTENT(IN)    :: W(1:NZ)
  !  TYPE(chp_chn_type),INTENT(IN)    :: CHN
  !  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
  !  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
  !  REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
  !  INTEGER :: j
  !  
  !  rj = REAL(CHN%j,kind=8)
  !  jj = 2.0D0*rj+1.0D0
  !  jj1 = rj*(rj+1.0D0)
  !  j = CHN%j
  !  
  !  pp  = pfinal; p  = pinit
  !  pp2 = pp*pp ; p2 = p*p
  !  ! uncoupled singlet
  !  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
  !     
  !     POT(1) = POT(1) + 2.0D0*pp2*p2*(pwd_integral(W,2,j) -pwd_integral(W,0,j))
  !     
  !  END IF
  !  
  !  IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
  !     
  !     POT(2) = POT(2) - 2.0D0*pp2*p2*(pwd_integral(W,0,j) +3.0D0*pwd_integral(W,2,j) &
  !          -2.0D0*(pwd_integral(W,1,j-1) + pwd_integral(W,1,j+1)))
  !     
  !     ! once here, return
  !     RETURN
  !  END IF
  !  
  !  ! coupled channels
  !  IF (CHN%coup) THEN
  !     
  !     POT(3) = POT(3) + 4.0D0*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,1,j))/jj + &
  !          2.0D0*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,2,j+1))
  !     IF (j==0) RETURN
  !     POT(4) = POT(4) + 4.0D0*pp2*p2*(pwd_integral(W,1,j) - pwd_integral(W,0,j-1))/jj + &
  !          2.0D0*pp2*p2*(pwd_integral(W,0,j-1) - pwd_integral(W,2,j-1))
  !     
  !     POT(5) = POT(5) - 4.0D0*DSQRT(jj1)*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj**2
  !
  !     POT(6) = POT(6) - 4.0D0*DSQRT(jj1)*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj**2
  !
  !  END IF
  !  
  !END SUBROUTINE pwd_sigL_E

  !   T     tensor
  ! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_T(W, CHN, pfinal, pinit, POT)
    
    REAL(8)            ,INTENT(IN)    :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)            ,INTENT(INOUT) :: POT(1:6)
    REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
    INTEGER :: j
  
    rj  = REAL(CHN%j,kind=8)
    jj  = 2.0D0*rj+1.0D0
    jj1 = rj*(rj+1.0D0)
    j   = CHN%j
    
    pp  = pfinal; p  = pinit
    pp2 = pp*pp ; p2 = p*p
    ! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 2.0D0*(-(pp2 + p2) * pwd_integral(W,0,j) + &
            2.0D0*pp*p*pwd_integral(W,1,j))
       
    ELSE IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*( (pp2 + p2)*pwd_integral(W,0,j) - &
            2.0D0*pp*p*(rj*pwd_integral(W,0,j+1) + &
            (rj+1.0D0)*pwd_integral(W,0,j-1))/jj)
       
       ! coupled channels
    ELSE IF(CHN%coup) THEN
       
       POT(3) = POT(3) + 2.0D0*(-(pp2 + p2) * pwd_integral(W,0,j+1) + &
            2.0D0*pp*p*pwd_integral(W,0,j))/jj
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*(+(pp2 + p2) * pwd_integral(W,0,j-1) - &
            2.0D0*pp*p*pwd_integral(W,0,j))/jj
       
       POT(5) = POT(5) - 4.0D0*DSQRT(jj1)*(p2*pwd_integral(W,0,j+1) + &
            pp2*pwd_integral(W,0,j-1) - 2.0D0*pp*p* &
            pwd_integral(W,0,j) )/jj
       
       POT(6) = POT(6) - 4.0D0*DSQRT(jj1)*(p2*pwd_integral(W,0,j-1) + &
            pp2*pwd_integral(W,0,j+1) - 2.0D0*pp*p* &
            pwd_integral(W,0,j) )/jj
  
    END IF
    
  END SUBROUTINE pwd_T

  !   T     tensor
  ! ref: E. Epelbaum et al., NPA 747 (2005), 362
  !SUBROUTINE pwd_T_E(W, CHN, pfinal, pinit, POT)
  !
  !  REAL(8)            ,INTENT(IN)    :: W(1:NZ)
  !  TYPE(chp_chn_type),INTENT(IN)    :: CHN
  !  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
  !  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
  !  REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
  !  INTEGER :: j
  !
  !  rj  = REAL(CHN%j,kind=8)
  !  jj  = 2.0D0*rj+1.0D0
  !  jj1 = rj*(rj+1.0D0)
  !  j   = CHN%j
  !
  !  pp  = pfinal; p  = pinit
  !  pp2 = pp*pp ; p2 = p*p
  !
  !  ! uncoupled singlet
  !  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
  !
  !     POT(1) = POT(1)-2.0D0*(pp2+p2)*pwd_integral(W,0,j) + 4.0D0*pp*p*pwd_integral(W,1,j)
  !
  !  ELSE IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
  !
  !     POT(2) = POT(2) + 2.0D0*((pp2+p2)*pwd_integral(W,0,j) + 2.0D0*pp*p*pwd_integral(W,1,j) - &
  !          2.0D0*pp*p*(pwd_integral(W,0,j-1)+pwd_integral(W,0,j+1)))
  !
  !     ! coupled channels
  !  ELSE IF (CHN%coup) THEN
  !
  !     POT(3) = POT(3) + (4.0D0*pp*p*pwd_integral(W,0,j)-2.0D0*(pp2+p2)*pwd_integral(W,0,j+1))/jj
  !     IF (j==0) RETURN
  !     POT(4) = POT(4) + (-4.0D0*pp*p*pwd_integral(W,0,j)-2.0D0*(pp2+p2)*pwd_integral(W,0,j-1))/jj
  !
  !     POT(5) = POT(5) - DSQRT(jj1)*(-8.0D0*pp*p*pwd_integral(W,0,j) + 4.0D0*pp2*pwd_integral(W,0,j-1) + &
  !          4.0D0*p2*pwd_integral(W,0,j+1) )/jj
  !
  !     POT(6) = POT(6) - DSQRT(jj1)*(-8.0D0*pp*p*pwd_integral(W,0,j) + 4.0D0*pp2*pwd_integral(W,0,j+1) + &
  !          4.0D0*p2*pwd_integral(W,0,j-1) )/jj
  !
  !  END IF
  !
  !END SUBROUTINE pwd_T_E
  
  !   sigk momentum-tensor
  ! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_sigk(W, CHN, pfinal, pinit, POT)
    
    REAL(8)            ,INTENT(IN)    :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)            ,INTENT(INOUT) :: POT(1:6)
    REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
    INTEGER :: j
    
    rj = REAL(CHN%j,kind=8)
    jj = 2.0D0*rj+1.0D0
    jj1 = rj*(rj+1.0D0)
    j = CHN%j
    
    pp  = pfinal; p  = pinit
    pp2 = pp*pp ; p2 = p*p
    ! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 0.5D0*(-(pp2 + p2) * pwd_integral(W,0,j) - &
            2.0D0*pp*p*pwd_integral(W,1,j))
       
    END IF
    
    IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 0.5D0*( (pp2 + p2)*pwd_integral(W,0,j) + &
            2.0D0*pp*p*(rj*pwd_integral(W,0,j+1) + &
            (rj+1.0D0)*pwd_integral(W,0,j-1))/jj)
       
       ! once here, return
       RETURN
    END IF
    
    ! coupled channels
    IF (CHN%coup) THEN
       
       POT(3) = POT(3) + 0.5D0*(-(pp2 + p2) * pwd_integral(W,0,j+1) - &
            2.0D0*pp*p*pwd_integral(W,0,j))/jj
       IF (j==0) RETURN
       POT(4) = POT(4) + 0.5D0*(+(pp2 + p2) * pwd_integral(W,0,j-1) + &
            2.0D0*pp*p*pwd_integral(W,0,j))/jj
       
       POT(5) = POT(5) - DSQRT(jj1)*(p2*pwd_integral(W,0,j+1) + &
            pp2*pwd_integral(W,0,j-1) + 2.0D0*pp*p* &
            pwd_integral(W,0,j))/jj
       
       POT(6) = POT(6) - DSQRT(jj1)*(p2*pwd_integral(W,0,j-1) + &
            pp2*pwd_integral(W,0,j+1) + 2.0D0*pp*p* &
            pwd_integral(W,0,j))/jj

    END IF
    
  END SUBROUTINE pwd_sigk

  !   sigk     momentum-tensor
  ! ref: E. Epelbaum et al., NPA 747 (2005), 362
  !SUBROUTINE pwd_sigk_E(W, CHN, pfinal, pinit, POT)
  !  
  !  REAL(8)            ,INTENT(IN)    :: W(1:NZ)
  !  TYPE(chp_chn_type),INTENT(IN)    :: CHN
  !  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
  !  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
  !  REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
  !  INTEGER :: j
  !  
  !  rj = REAL(CHN%j,kind=8)
  !  jj = 2.0D0*rj+1.0D0
  !  jj1 = rj*(rj+1.0D0)
  !  j = CHN%j
  !  
  !  pp  = pfinal; p  = pinit
  !  pp2 = pp*pp ; p2 = p*p
  !
  !  ! uncoupled singlet
  !  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
  !     
  !     POT(1) = POT(1)-0.5D0*(pp2+p2)*pwd_integral(W,0,j) - pp*p*pwd_integral(W,1,j)
  !     
  !  END IF
  !  
  !  IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
  !     
  !     POT(2) = POT(2) + 0.5D0*(pp2+p2)*pwd_integral(W,0,j) - pp*p*pwd_integral(W,1,j) + &
  !          pp*p*(pwd_integral(W,0,j-1)+pwd_integral(W,0,j+1))
  !     
  !     ! once here, return
  !     RETURN
  !  END IF
  !  
  !  ! coupled channels
  !  IF (CHN%coup) THEN
  !     
  !     POT(3) = POT(3) + (-1.0D0*pp*p*pwd_integral(W,0,j)-0.5D0*(pp2+p2)*pwd_integral(W,0,j+1))/jj
  !     IF (j==0) RETURN
  !     POT(4) = POT(4) + (+1.0D0*pp*p*pwd_integral(W,0,j)+0.5D0*(pp2+p2)*pwd_integral(W,0,j-1))/jj
  !     
  !     POT(5) = POT(5) - (DSQRT(jj1)*(2.0D0*pp*p*pwd_integral(W,0,j) + pp2*pwd_integral(W,0,j-1) + &
  !          p2*pwd_integral(W,0,j+1)))/jj
  !
  !     POT(6) = POT(6) - (DSQRT(jj1)*(2.0D0*pp*p*pwd_integral(W,0,j) + pp2*pwd_integral(W,0,j+1) + &
  !          p2*pwd_integral(W,0,j-1)))/jj
  !
  !  END IF
  !  
  !END SUBROUTINE pwd_sigk_E
  
  FUNCTION pwd_integral(W,l,j) RESULT(res)
    
    REAL(8), INTENT(IN)  :: W(1:NZ)
    INTEGER, INTENT(IN) :: l, j
    REAL(8) :: res
    
    res = 0.0D0

    IF (j<0) RETURN
    
    ! not using precalculated z**l
    !res = chp_pi%val*SUM(W(:)*zmesh%pnt(:)%x**l*legP(:,j)*zmesh%pnt(:)%w)
    ! using precalculated z**l (atleast a factor of 2 in speedup)
    res = chp_pi%val*SUM(W(:)*zl(:,l)*legP(:,j)*zmesh%pnt(:)%w)
    
  END FUNCTION pwd_integral
  
END MODULE twobody_pwd

![1] = Machleidts Phys Rep
![2] = Epelbaum NPA
!
MODULE idaho_chiral_potential
  
  USE chp_aux
  USE twobody_pwd

  IMPLICIT NONE
  
  !MAXIMUM J
  INTEGER, PARAMETER, PRIVATE :: maximum_angular_momentum = 40 ! maximum angular momentum. 
  ! Z=cos(theta), theta angle between p' and p, i.e., final and initial momenta
  ! q = (p'-p)
  ! q^2 = p'^2 + p^2 - 2p'pZ
  INTEGER, PARAMETER, PRIVATE :: nof_theta_int_points   = 96
  INTEGER, PARAMETER, PRIVATE :: nof_2PE_2loop_int_x_points = 50
  INTEGER, PARAMETER, PRIVATE :: nof_2PE_2loop_int_z_points = 30
  
  ! FORTRAN UNITS
  INTEGER, PARAMETER, PRIVATE :: CHP_ERR = 6
  INTEGER, PARAMETER, PRIVATE :: CHP_SCR = 6
  
  ! chiral order definition
  INTEGER, PRIVATE, PARAMETER :: LO   = 0
  ! all contributions vanish at order 1
  ! due to parity and time-reversal  invariance
  INTEGER, PARAMETER :: NLO  = 2
  INTEGER, PARAMETER :: NNLO = 3
  INTEGER, PARAMETER :: N3LO = 4
  
  ! set the type of iterated two-pion exchange
  ! formulation to keep relativistic elastic
  ! unitarity. See [1] Eq. 4.25 and 4.25
  ! EM : Entem Machleidt, equivalent to using Eq. 4.26
  ! KW : Kaiser Brockmann Weise, equivalent to using Eq. 4.25
  
  TYPE(chp_char2_type), PRIVATE :: chp_itope
  
  REAL(8), PRIVATE  :: rirrel = -99.99D0
  INTEGER, PRIVATE :: iirrel = -99
  CHARACTER(LEN=2), PRIVATE :: cirrel = 'XX'
  ! masses, MeV
  TYPE(chp_real_type), PRIVATE :: chp_mnuc(-1:1) ! -1: pp , 0: pn , +1: nn
  TYPE(chp_real_type), PRIVATE :: chp_mpi(-1:2)  ! -1: pi-, 0: pi0, +1: pi+, 2: average
  ! chiral order
  TYPE(chp_int_type), PRIVATE :: chp_chiral_order ! LO, NLO, NNLO, N3LO
  
  ! renormalization stuff
  TYPE(chp_real_type), PRIVATE :: chp_ren ! SF or DR, with SF-cut in MeV
  ! regularization stuff
  TYPE(chp_real_type), PRIVATE :: chp_lambda ! regularization cutoff in MeV, typically ~ 500
  TYPE(chp_real_type), PRIVATE :: chp_regcut_1PE, chp_regcut_2PE
  ! potential coupling constants
  TYPE(chp_real_type), PRIVATE :: chp_gA  ! dimensionless
  TYPE(chp_real_type), PRIVATE :: chp_fpi ! MeV
  TYPE(chp_real_type), PRIVATE :: chp_fine_structure ! dimensionless

  ! optional features
  TYPE(chp_int_type) , PRIVATE :: chp_2PE_CSB_correct_mass ! 1 or 0
  TYPE(chp_int_type) , PRIVATE :: chp_use_2PE_2loop_int ! 1 or 0
  
  ! variables for calculating optional 2PE 2loop diagrams
  LOGICAL, PRIVATE :: chp_2PE_2loop_int_VTS_data_set
  LOGICAL, PRIVATE :: chp_2PE_2loop_int_WC_data_set
  TYPE(chp_gauleg_mesh), PRIVATE :: chp_2PE_2loop_int_mesh_x
  TYPE(chp_gauleg_mesh), PRIVATE :: chp_2PE_2loop_int_mesh_z
  REAL(8) :: chp_2PE_2loop_int_VTS_x_fact(nof_2PE_2loop_int_z_points)
  REAL(8) :: chp_2PE_2loop_int_WC_x_fact(nof_2PE_2loop_int_z_points)

  ! LECs ci
  TYPE(chp_real_type) chp_c1
  TYPE(chp_real_type) chp_c3
  TYPE(chp_real_type) chp_c4

  ! New at N3LO
  TYPE(chp_real_type) chp_c2
  TYPE(chp_real_type) chp_d1_plus_d2
  TYPE(chp_real_type) chp_d3
  TYPE(chp_real_type) chp_d5
  TYPE(chp_real_type) chp_d14_minus_d15

  ! CONTACT TERMS
  ! THEY CAN BE INPUT IN PWD FORMAT
  ! OR STANDARD FORMAT. THE CODE USES
  ! THE PWD FORMAT. WHEN PRINTED, THEY
  ! ARE CONVERTED AND PRINTED IN BOTH
  TYPE(chp_int_type), PRIVATE :: chp_contact_format! PW (partial-wave) or ST (standard)
  !
  ! @LO      PW         
  ! CLO(1) : Ctilde 1S0  CS
  ! CLO(2) : Ctilde 3S1  CT
  TYPE(chp_real_type), PRIVATE :: chp_CLO(1:2), chp_CLO_n(1:2)
  ! @NLO     PW
  ! CLO(1) : C1S0     C1
  ! CLO(2) : C3P0     C2
  ! CLO(3) : C1P1     C3
  ! CLO(4) : C3P1     C4
  ! CLO(5) : C3S1     C5
  ! CLO(6) : C3S1-3D1 C6
  ! CLO(7) : C3P1     C7
  TYPE(chp_real_type), PRIVATE :: chp_CNLO(1:7), chp_CNLO_n(1:7)
  !  at the NLO order we also have CIB LO contacts
  TYPE(chp_real_type), PRIVATE :: chp_CIB_CLO(-1:1,1:2)
  ! NOTE: Ctilde pp 3S1 = Ctilde pn 3S1 = Ctilde nn 3S1
  ! i.e., all CIB in 1S0 contact
  ! CIB_CLO(-1,1) : Ctilde pp 1S0  CS pp
  ! CIB_CLO(-1,2) : Ctilde pp 3S1  CT pp
  ! CIB_CLO( 0,1) : Ctilde pn 1S0  CS pn
  ! CIB_CLO( 0,2) : Ctilde pn 3S1  CT pn
  ! CIB_CLO(+1,1) : Ctilde nn 1S0  CS nn
  ! CIB_CLO(+1,2) : Ctilde nn 3S1  CT nn
  
  ! @N3LO       PW
  ! DN3LO( 1) : D1S0t     D1
  ! DN3LO( 2) : D1S0      D2
  ! DN3LO( 3) : D3P0      D3
  ! DN3LO( 4) : D1P1      D4
  ! DN3LO( 5) : D3P1      D5
  ! DN3LO( 6) : D3S1t     D6
  ! DN3LO( 7) : D3S1      D7
  ! DN3LO( 8) : D3D1      D8
  ! DN3LO( 9) : D3S1-3D1t D9
  ! DN3LO(10) : D3S1-3D1  D10
  ! DN3LO(11) : D1D2      D11
  ! DN3LO(12) : D3D2      D12
  ! DN3LO(13) : D3P2      D13
  ! DN3LO(14) : D3P2-3F2  D14
  ! DN3LO(15) : D3D3      D15
  TYPE(chp_real_type), PRIVATE :: chp_DN3LO(1:15), chp_DN3LO_n(1:15)

  ! LIST OF OPERATOR STRUCTURES
  
  REAL(8) :: Vc(nof_theta_int_points)
  REAL(8) :: Wc(nof_theta_int_points)
  REAL(8) :: Vs(nof_theta_int_points)
  REAL(8) :: Ws(nof_theta_int_points)
  REAL(8) :: VLS(nof_theta_int_points)
  REAL(8) :: WLS(nof_theta_int_points)
  REAL(8) :: VsigL(nof_theta_int_points)
  REAL(8) :: WsigL(nof_theta_int_points)
  REAL(8) :: VT(nof_theta_int_points)
  REAL(8) :: WT(nof_theta_int_points)
  REAL(8) :: Vsigk(nof_theta_int_points)
  REAL(8) :: Wsigk(nof_theta_int_points)
  
  ! DERIVED CONSTANTS THAT ARE REPEATEDLY USED, ONCE THE MASSES AND COUPLINGS ARE SET
  ! THESE ARE COMPUTED IN THE ROUTINE set_units_and_derive_constants.
  
  REAL(8) :: c1,c3,c4
  REAL(8) :: c2, d1_plus_d2, d3, d5, d14_minus_d15
  REAL(8) :: mnuc2(-1:1)   ! nucleon mass squared
  REAL(8) :: mnuc_inv(-1:1)! nucleon mass inversed
  REAL(8) :: mpi2(-1:2)    ! pion mass squared
  REAL(8) :: mpi3(-1:2)    ! pion mass cubed
  REAL(8) :: mpi4(-1:2)    ! mpi2 squared
  REAL(8) :: mpi5(-1:2)    ! mpi^5
  REAL(8) :: twompi(-1:2)  ! two times pion mass
  REAL(8) :: fourmpi2(-1:2)! four times pion mass squared
  
  REAL(8) :: sfr                 ! sfr cutoff
  REAL(8) :: sfr2                ! sfr cutoff squared
  REAL(8) :: sfr_heavyside(-1:2) ! THETA(sfr-twompi)

  REAL(8),PRIVATE :: iso(0:1)      ! expectation value of tau_1 \cdot tau_2
  REAL(8) :: fpi2          ! fpi squared
  REAL(8) :: fpi4          ! fpi2 squared
  REAL(8) :: fpi_inv       ! fpi inverse
  REAL(8) :: gA2           ! gA squared
  REAL(8) :: gA4           ! gA2 squared
  
  REAL(8) :: const(50)      

 CONTAINS
  ! uncoupled singlet
  ! POT(1) = <L=J  ,S=0,J|V|L=J  ,S=0,J>
  ! uncoupled triplet
  ! POT(2) = <L=J  ,S=1,J|V|L=J  ,S=1,J>
  ! coupled triplets
  ! POT(3) = <L=J+1,S=1,J|V|L=J+1,S=0,J>
  ! POT(4) = <L=J-1,S=1,J|V|L=J-1,S=1,J>
  ! POT(5) = <L=J+1,S=1,J|V|L=J-1,S=1,J>
  ! POT(6) = <L=J-1,S=1,J|V|L=J+1,S=1,J>
  !
  SUBROUTINE chp(pout, pin, coup, S, j, T, tz, POT)
    
    REAL(8) ,INTENT(IN)    :: pout   ! initial relative momentum
    REAL(8) ,INTENT(IN)    :: pin    ! final relative momentum  
    LOGICAL,INTENT(IN)    :: coup   ! IF true, coupled channel
    INTEGER,INTENT(IN)    :: S      ! total spin: 0 or 1
    INTEGER,INTENT(IN)    :: j      ! relative angular momentum
    INTEGER,INTENT(IN)    :: T      ! total isospin: 0 or 1
    INTEGER,INTENT(IN)    :: tz     ! isospin projection: -1, 0, 1 ! pp pn nn
    REAL(8) ,INTENT(INOUT) :: POT(6) ! potential in LSJ formalism ! MeV
    TYPE(chp_chn_type)    :: CHN    ! coup S, j, T, tz
    REAL(8)                :: vlo(6), vnlo(6), vnnlo(6), vn3lo(6), v_contact(6)
    !
    ! Z=cos(theta), theta angle between p' and p, i.e., final and initial momenta
    ! q = (p'-p)
    ! q^2 = p'^2 + p^2 - 2p'pZ
    REAL(8) :: pfinal, pinit
    REAL(8) :: pinit2 ! pinit^2
    REAL(8) :: pfinal2! pfinal^2
    REAL(8) :: q2(1:nof_theta_int_points) ! q^2
    REAL(8) :: q(1:nof_theta_int_points)  ! q
    REAL(8) :: k2(1:nof_theta_int_points) ! k^2
    ! loop functions
    REAL(8) :: loop_w(nof_theta_int_points)
    REAL(8) :: loop_s
    REAL(8) :: loop_L(nof_theta_int_points)
    REAL(8) :: loop_w_tz(nof_theta_int_points)
    REAL(8) :: loop_s_tz
    REAL(8) :: loop_L_tz(nof_theta_int_points)
    REAL(8) :: loop_wtilde2(nof_theta_int_points)
    REAL(8) :: loop_A(nof_theta_int_points)
    INTEGER :: imnuc_2PE
    !

    POT     = 0.0D0
    Vc      = 0.0D0
    Wc      = 0.0D0
    Vs      = 0.0D0
    Ws      = 0.0D0
    VLS     = 0.0D0
    WLS     = 0.0D0
    VsigL   = 0.0D0
    WsigL   = 0.0D0
    VT      = 0.0D0
    WT      = 0.0D0
    Vsigk   = 0.0D0
    Wsigk   = 0.0D0
    
    vlo = 0.0D0 ; vnlo = 0.0D0 ; vnnlo = 0.0D0 ; vn3lo = 0.0D0
    v_contact = 0.0D0
    
    loop_w       = 0.0D0
    loop_s       = 0.0D0
    loop_L       = 0.0D0
    loop_wtilde2 = 0.0D0
    loop_A       = 0.0D0

    IF (S==0 .AND. coup) THEN
       WRITE(CHP_ERR,"(A)")
       WRITE(CHP_ERR,"(A)") '   ****************************'
       WRITE(CHP_ERR,"(A)") '   error(chp): illegal channel'
       WRITE(CHP_ERR,"(A)") '   ****************************'
       WRITE(CHP_ERR,"(A)")
       RETURN
    END IF
    !
    ! set the channel quantum numbers
    CHN%j    = j
    CHN%S    = S
    CHN%T    = T
    CHN%tz   = tz
    CHN%coup = coup
    !
    !DERIVED QUANTITIES
    ! initial and final momenta squared
    pinit = pin
    pfinal = pout
    pinit2  = pinit*pinit
    pfinal2 = pfinal*pfinal
    
    ! squared momentum transfer
    q2 = pfinal2 + pinit2 - 2.0D0*pinit*pfinal*zmesh%pnt(:)%x
    ! momentum transfer
    q  = DSQRT(q2)

    !@LO
    !static one-pion exchange
    IF (chp_chiral_order%val == LO) THEN
       
       ! static OPE contribution
       ! using average pion mass only
       !
       WT = iso(T)*chp_one_pion_exchange(q2, 2)
       
       CALL pwd_T(WT, CHN, pfinal, pinit, vlo)
       
       ! add LO contact terms
       CALL chp_LO_contact_terms(CHN, pinit, pfinal, v_contact)

    ELSE
       
       !@NLO
       !NLO loop functions
       loop_w = chp_NLO_two_pion_exchange_loop_w(q2,2)
       IF (chp_ren%name == 'SF') THEN
          loop_s = chp_NLO_sfr_two_pion_exchange_loop_s(2)
          loop_L = chp_NLO_sfr_two_pion_exchange_loop_L(q, q2, loop_w, loop_s, 2)
       ELSE IF(chp_ren%name == 'DR') THEN
          loop_L = chp_NLO_dr_two_pion_exchange_loop_L(q,loop_w,2)
       END IF
       
       ! LO CIB one-pion exchanges
       ! [1] Eq. 4.77-4.79
       IF (tz == 0) THEN
           WT = -1.0D0*chp_one_pion_exchange(q2, 0) + &
                chp_minus_power(CHN%T+1)*2.0D0*chp_one_pion_exchange(q2,1)
           if(chp_chiral_order%val >= N3LO) then
              WT = WT + chp_minus_power(CHN%T+1)*2.0D0*chp_one_pion_exchange_gamma(q2, 1)
           end if
       ELSE
          WT = chp_one_pion_exchange(q2, 0)
       END IF
       
       ! add contacts to vlo
       CALL chp_CIB_LO_contact_terms(CHN, pinit, pfinal, v_contact)
       
       ! irreducible, or non-polynomial, NLO two-pion exchanges
       if(chp_chiral_order%val > N3LO) then
           if(chn%T == 0) then
               Wc = iso(T) * chp_NLO_two_pion_exchange_Wc(q2, loop_L, loop_w, 2)
           else
               loop_w_tz = chp_NLO_two_pion_exchange_loop_w(q2,chn%tz)
               IF (chp_ren%name == 'SF') THEN
                  loop_s_tz = chp_NLO_sfr_two_pion_exchange_loop_s(chn%tz)
                  loop_L_tz = chp_NLO_sfr_two_pion_exchange_loop_L(q, q2, loop_w_tz, loop_s_tz, chn%tz)
               ELSE IF(chp_ren%name == 'DR') THEN
                  loop_L_tz = chp_NLO_dr_two_pion_exchange_loop_L(q,loop_w_tz,chn%tz)
               END IF
               Wc = iso(T) * chp_NLO_two_pion_exchange_Wc(q2, loop_L_tz, loop_w_tz, chn%tz)
           end if
       else
           Wc = iso(T) * chp_NLO_two_pion_exchange_Wc(q2, loop_L, loop_w, 2)
       end if
       Vs =          chp_NLO_two_pion_exchange_Vs(q2, loop_L, 2)
       VT =          chp_NLO_two_pion_exchange_VT(loop_L, 2)
       
       ! project to LSJ basis and add to vnlo array
       CALL pwd_c(Wc, CHN, pfinal, pinit, vnlo)
       CALL pwd_s(Vs, CHN, pfinal, pinit, vnlo)
       CALL pwd_T(VT, CHN, pfinal, pinit, vnlo)
       ! project to LSJ basis and add to vlo
       CALL pwd_T(WT, CHN, pfinal, pinit, vlo)

       ! add contacts to vnlo
       CALL chp_NLO_contact_terms(CHN, pinit, pfinal, v_contact)
       
       !@NNLO
       !
       IF (chp_chiral_order%val > NLO) THEN
          
          if(chp_2PE_CSB_correct_mass%val == 0) then
              imnuc_2PE = 0
          else
              imnuc_2PE = tz
          end if
          ! NNLO loop functions
          IF (chp_ren%name == 'SF') THEN
             loop_wtilde2 = chp_NNLO_two_pion_exchange_loop_wtilde2(q2,2)
             loop_A = chp_NNLO_sfr_two_pion_exchange_loop_A(q, q2, 2)
          ELSE IF (chp_ren%name == 'DR') THEN
             loop_wtilde2 = chp_NNLO_two_pion_exchange_loop_wtilde2(q2,2)
             loop_A       = chp_NNLO_dr_two_pion_exchange_loop_A(q,2)
          END IF
          
          ! chp_itope type is set in the header of this module.
          Vc  =           chp_NNLO_two_pion_exchange_Vc (q2, loop_w, loop_A, loop_wtilde2, 2, imnuc_2PE)
          Wc  = iso(T) *  chp_NNLO_two_pion_exchange_Wc (q2, loop_w, loop_A, loop_wtilde2, 2, imnuc_2PE)
          VLS =           chp_NNLO_two_pion_exchange_VLS(            loop_A, loop_wtilde2, 2, imnuc_2PE)
          WLS = iso(T) *  chp_NNLO_two_pion_exchange_WLS(    loop_w, loop_A              , 2, imnuc_2PE)
          VT  =           chp_NNLO_two_pion_exchange_VT (    loop_w, loop_A, loop_wtilde2, 2, imnuc_2PE)
          Vs  =           chp_NNLO_two_pion_exchange_Vs (q2, loop_w, loop_A, loop_wtilde2, 2, imnuc_2PE)
          WT  = iso(T) *  chp_NNLO_two_pion_exchange_WT (q2, loop_w, loop_A              , 2, imnuc_2PE)
          Ws  = iso(T) *  chp_NNLO_two_pion_exchange_Ws (q2, loop_w, loop_A              , 2, imnuc_2PE)
          
          ! project to LSJ basis and add to vnnlo array
          
          CALL pwd_c(Vc, CHN, pfinal, pinit, vnnlo)

          CALL pwd_c(Wc, CHN, pfinal, pinit, vnnlo)

          CALL pwd_LS(VLS,  CHN, pfinal, pinit, vnnlo)

          CALL pwd_LS(WLS,  CHN, pfinal, pinit, vnnlo)

          CALL pwd_T(VT, CHN, pfinal, pinit, vnnlo)  

          CALL pwd_s(Vs, CHN, pfinal, pinit, vnnlo)

          CALL pwd_T(WT, CHN, pfinal, pinit, vnnlo)  

          CALL pwd_s(Ws, CHN, pfinal, pinit, vnnlo)            
          
       END IF
       
       !@N3LO
       IF (chp_chiral_order%val > NNLO)THEN
          k2 = 0.25D0 * (pfinal2 + pinit2) + 0.5D0*pinit*pfinal*zmesh%pnt(:)%x
          ! Add contacts to v3lo
          CALL chp_N3LO_contact_terms(CHN, pinit, pfinal, v_contact)

          Vc  =           chp_N3LO_2PE_Vc_1loop_ci2 (        loop_w, loop_L, loop_wtilde2,         2) + &
                          chp_N3LO_2PE_Vc_1loop_ci1 (    q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_Vc_1loop_ci0 (    q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_Vc_2loop     (    q2,                 loop_wtilde2, loop_A, 2)
          Wc  = iso(T) * (chp_N3LO_2PE_Wc_1loop_ci1 (    q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_Wc_1loop_ci0 (k2, q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_Wc_2loop_a   (    q2, loop_w, loop_L, loop_wtilde2,         2) + &
                          chp_N3LO_2PE_Wc_2loop_b   (    q2,                                       2))
          VLS =           chp_N3LO_2PE_VLS_1loop_ci1(        loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_VLS_1loop_ci0(    q2, loop_w, loop_L,                       2, imnuc_2PE)
          WLS = iso(T) * (chp_N3LO_2PE_WLS_1loop_ci1(    q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_WLS_1loop_ci0(    q2, loop_w, loop_L,                       2, imnuc_2PE))
          VT  =           chp_N3LO_2PE_VT_1loop_ci0 (k2, q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_VT_2loop_a   (        loop_w, loop_L,                       2) + &
                          chp_N3LO_2PE_VT_2loop_b   (    q2,                                       2)
          WT  = iso(T) * (chp_N3LO_2PE_WT_1loop_ci2 (        loop_w, loop_L,                       2) + &
                          chp_N3LO_2PE_WT_1loop_ci1 (    q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_WT_1loop_ci0 (    q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_WT_2loop     (        loop_w,                       loop_A, 2))
          Vs  =           chp_N3LO_2PE_Vs_1loop_ci0 (k2, q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_Vs_2loop_a   (    q2, loop_w, loop_L,                       2) + &
                          chp_N3LO_2PE_Vs_2loop_b   (    q2,                                       2)
          Ws  = iso(T) * (chp_N3LO_2PE_Ws_1loop_ci2 (    q2, loop_w, loop_L,                       2) + &
                          chp_N3LO_2PE_Ws_1loop_ci1 (    q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_Ws_1loop_ci0 (    q2, loop_w, loop_L,                       2, imnuc_2PE) + &
                          chp_N3LO_2PE_Ws_2loop     (    q2, loop_w,                       loop_A, 2))
          VsigL =         chp_N3LO_2PE_VsigL_1loop_ci0(              loop_L,                       2, imnuc_2PE)

          CALL pwd_c   (Vc   , CHN, pfinal, pinit, vn3lo)
          CALL pwd_c   (Wc   , CHN, pfinal, pinit, vn3lo)
          CALL pwd_LS  (VLS  , CHN, pfinal, pinit, vn3lo)
          CALL pwd_LS  (WLS  , CHN, pfinal, pinit, vn3lo)
          CALL pwd_T   (VT   , CHN, pfinal, pinit, vn3lo)  
          CALL pwd_T   (WT   , CHN, pfinal, pinit, vn3lo)  
          CALL pwd_s   (Vs   , CHN, pfinal, pinit, vn3lo)
          CALL pwd_s   (Ws   , CHN, pfinal, pinit, vn3lo)
          CALL pwd_sigL(VsigL, CHN, pfinal, pinit, vn3lo)            
       END IF
       

    END IF
    
    ! adopt minimal relativity
    vlo   = vlo   * const(1)*chp_mnuc(CHN%tz)%val/DSQRT(cm_energy(pinit,CHN%tz)*cm_energy(pfinal,CHN%tz))
    vnlo  = vnlo  * const(1)*chp_mnuc(CHN%tz)%val/DSQRT(cm_energy(pinit,CHN%tz)*cm_energy(pfinal,CHN%tz))
    vnnlo = vnnlo * const(1)*chp_mnuc(CHN%tz)%val/DSQRT(cm_energy(pinit,CHN%tz)*cm_energy(pfinal,CHN%tz))
    vn3lo = vn3lo * const(1)*chp_mnuc(CHN%tz)%val/DSQRT(cm_energy(pinit,CHN%tz)*cm_energy(pfinal,CHN%tz))
    v_contact = v_contact * const(1)*chp_mnuc(CHN%tz)%val/DSQRT(cm_energy(pinit,CHN%tz)*cm_energy(pfinal,CHN%tz))
    
    ! add the contributions computed above
    ! and regulate them according to their chiral order
    POT = vlo  *freg(pfinal,pinit,chp_regcut_1PE%val) + &
          vnlo *freg(pfinal,pinit,chp_regcut_2PE%val) + &
          vnnlo*freg(pfinal,pinit,chp_regcut_2PE%val) + &
          vn3lo*freg(pfinal,pinit,chp_regcut_2PE%val) + &
          v_contact
    
  END SUBROUTINE chp
  
  SUBROUTINE chp_setup_2PE_2loop_int_data(nof_x_points, nof_z_points)
    INTEGER, INTENT(IN) :: nof_x_points
    INTEGER, INTENT(IN) :: nof_z_points

    chp_2PE_2loop_int_mesh_x%info = chp_mesh_info(nof_x_points, 0.0D0, 1.0D0)
    CALL chp_setup_gauleg_mesh(chp_2PE_2loop_int_mesh_x)

    chp_2PE_2loop_int_mesh_z%info = chp_mesh_info(nof_z_points, 0.0D0, 1.0D0)
    CALL chp_setup_gauleg_mesh(chp_2PE_2loop_int_mesh_z)
  END SUBROUTINE

  ! static one pion exchange, [1] eq 4.5
  ! without isospin structure
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use,
  FUNCTION chp_one_pion_exchange(q2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = -1.0D0 * const(2)/ (q2 + mpi2(impi))
    
  END FUNCTION chp_one_pion_exchange
  
  FUNCTION chp_one_pion_exchange_gamma(q2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) :: res(1:nof_theta_int_points)
    REAL(8) :: b2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    
    b2 = q2 / mpi2(impi)
    res = -(0.25d0 * chp_fine_structure%val * gA2 / (chp_pi%val * fpi2 * mpi2(impi))) * (-((1-b2)**2 / (2.0d0 * b2**2 * (1 + b2))) * dlog(b2+1) + 1.0d0 / (2*b2))

  END FUNCTION chp_one_pion_exchange_gamma
  
  ! NLO loop function w [1] Eq 4.12 (DR and SFR)
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_two_pion_exchange_loop_w(q2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = DSQRT(fourmpi2(impi) + q2)
    
  END FUNCTION chp_NLO_two_pion_exchange_loop_w

  ! NLO SFR loop function s [2] Eq 2.16
  ! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s(impi) RESULT(res)
    
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = 0.0D0
    res = DSQRT(sfr2 - fourmpi2(impi))
    
  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s

  ! NLO dr loop function L [1] Eq 4.11
  ! with isospin structure
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_dr_two_pion_exchange_loop_L(q, w, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = w * dlog((w+q)/(twompi(impi) ) )/q
    
  END FUNCTION chp_NLO_dr_two_pion_exchange_loop_L

  ! NLO SFR loop function L [2] Eq 2.16
  ! q   : momentum transfer
  ! q2  : momentum transfer squared
  ! w,s : SFR loop functions
  ! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L(q, q2, w, s, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: s
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    IF (sfr_heavyside(impi) == 0.0D0) return
    
    res = w * dlog( (sfr2*w*w+q2*s*s+2.0D0*sfr*q*w*s)/( fourmpi2(impi)*(sfr2+q2) ) )/(2.0D0*q)
    
  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L
  
  ! NLO function Eq 4.9
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_two_pion_exchange_Wc(q2, L, w, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = -L *(fourmpi2(impi)*const(4) + q2*const(5) + const(6)*mpi4(impi)/(w*w))/const(3)
    
  END FUNCTION chp_NLO_two_pion_exchange_Wc
  
  ! NLO function Eq 4.10
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_two_pion_exchange_Vs(q2,L,impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = +const(7)*L*q2/const(8)
    
  END FUNCTION chp_NLO_two_pion_exchange_Vs

  ! NLO Vt function Eq 4.10
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_two_pion_exchange_VT(L,impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = -const(7)*L/const(8)
    
  END FUNCTION chp_NLO_two_pion_exchange_VT

  ! NNLO loop function wtilde SQUARED [1] Eq 4.20 (DR)
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use,
  FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2(q2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = 2.0D0*mpi2(impi) + q2
    
  END FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2

  ! NNLO loop function wtilde [1] Eq 4.19 (DR)
  ! q   : momentum transfer
  ! impi: determines which mpi to use,
  FUNCTION chp_NNLO_dr_two_pion_exchange_loop_A(q, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = datan(q/(twompi(impi)))/(2.0D0*q)
    
  END FUNCTION chp_NNLO_dr_two_pion_exchange_loop_A

  ! NNLO loop function wtilde [2] Eq 2.17 (SFR)
  ! q   : momentum transfer
  ! q2  : momentum transfer squared
  ! impi: determines which mpi to use,
  FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A(q, q2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    IF (sfr_heavyside(impi) == 0.0D0) return

    res = datan( q*(sfr-twompi(impi) )/(q2 + sfr*twompi(impi) ) )/(2.0D0*q)
    
  END FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A

  ! NNLO function Eq [2] 3.38
  ! q2    : momentum transfer squared
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  FUNCTION chp_NNLO_sfr_two_pion_exchange_Vc(q2, A, wt2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    
    res = -const(9)*(mpi2(impi)*const(11) - q2*c3) *wt2*A
    
  END FUNCTION chp_NNLO_sfr_two_pion_exchange_Vc

  ! NNLO function Eq [1] 4.13 and 4.21
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Vc(q2, w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = const(9)*(const(10)*mpi5(impi)/(chp_mnuc(imnuc)%val*w*w) - &
         (mpi2(impi)*const(11) - q2*c3 - q2*3.0D0*const(10)/chp_mnuc(imnuc)%val ) *wt2*A)
    
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope%val == 'EM') THEN
       res = res - const(12)*(chp_mpi(impi)%val*w*w+wt2*wt2*A)/chp_mnuc(imnuc)%val
    END IF
 
  END FUNCTION chp_NNLO_two_pion_exchange_Vc

  ! NNLO function Eq [1] 4.14 and 4.22
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Wc(q2, w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = const(13) * (3.0D0*gA2*mpi5(impi)/(w*w) - &
         (fourmpi2(impi) + 2.0D0*q2 - gA2*(fourmpi2(impi)+3.0D0*q2))*wt2*A)/chp_mnuc(imnuc)%val
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope%val == 'EM') THEN
       res = res + const(14)*(chp_mpi(impi)%val*w*w + wt2*wt2*A)/chp_mnuc(imnuc)%val
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_Wc
  
  ! NNLO function Eq [1] 4.15 and 4.23
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_VT(w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = 3.0D0*const(15)*wt2*A/chp_mnuc(imnuc)%val
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope%val == 'EM') THEN
       res = res + const(15)*(chp_mpi(impi)%val + w*w*A )/chp_mnuc(imnuc)%val
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_VT

  ! NNLO function Eq [1] 4.15 and 4.23
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Vs(q2, w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = -3.0D0*q2*const(15)*wt2*A/chp_mnuc(imnuc)%val
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope%val == 'EM') THEN
       res = res - 1.0D0*const(15)*q2*(chp_mpi(impi)%val + w*w*A)/chp_mnuc(imnuc)%val
    END IF

  END FUNCTION chp_NNLO_two_pion_exchange_Vs
  
  ! NNLO function Eq [1] 4.16 and 4.24
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_WT(q2, w, A, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = -1.0D0*const(16)*A*( (c4 + 1.0D0/(4.0D0*chp_mnuc(imnuc)%val))*w*w - &
         const(17)*(10.0D0*mpi2(impi) + 3.0D0*q2)/chp_mnuc(imnuc)%val)
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope%val == 'EM') THEN
       res = res - const(18)*(chp_mpi(impi)%val + w*w*A)/chp_mnuc(imnuc)%val
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_WT

  ! NNLO function Eq [1] 4.16 and 4.24
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Ws(q2, w, A, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = q2*const(16)*A*( (c4+1.0D0/(4.0D0*chp_mnuc(imnuc)%val))*w*w - &
         const(17)*(10.0D0*mpi2(impi) +3.0D0*q2)/chp_mnuc(imnuc)%val)
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope%val == 'EM') THEN
       res = res + const(18)*q2*(chp_mpi(impi)%val + w*w*A)/chp_mnuc(imnuc)%val
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_Ws

  ! NNLO function Eq [1] 4.17
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  FUNCTION chp_NNLO_two_pion_exchange_VLS(A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = const(19) * wt2*A/chp_mnuc(imnuc)%val
             
  END FUNCTION chp_NNLO_two_pion_exchange_VLS

  ! NNLO function Eq [1] 4.18
  ! w   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  FUNCTION chp_NNLO_two_pion_exchange_WLS(w, A, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = const(20)*w*w*A/chp_mnuc(imnuc)%val
             
  END FUNCTION chp_NNLO_two_pion_exchange_WLS

  ! N3LO two pion exchange functions
  ! k2    : mean momentum squared
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! L     : L loop function (depends on what regularization is used)
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use

  ! [1]Eq. D.1
  FUNCTION chp_N3LO_2PE_Vc_1loop_ci2 (        w, L, wt2,    impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = c2 * (1D0/6D0) * w * w + c3 * wt2 - c1 * fourmpi2(impi)
    res = ((3D0 / 16D0) * const(21)) * L * (res*res + c2*c2*w*w*w*W*(1D0/45D0))
  END FUNCTION chp_N3LO_2PE_Vc_1loop_ci2

  ! [1]Eq. D.4
  FUNCTION chp_N3LO_2PE_Vc_1loop_ci1 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (-gA2 * const(21) * mnuc_inv(imnuc) * (1D0 / 32D0)) * L * ( &
              (c2 - 6D0*c3) * q2 * q2 + 4D0 * (6D0*c1 + c2 - 3D0*c3)*q2*mpi2(impi) + &
              6D0 * (c2 - 2D0*c3)*mpi4(impi) + &
              24D0*(2D0*c1 + c3)*mpi4(impi)*mpi2(impi)/(w*w) )
  END FUNCTION chp_N3LO_2PE_Vc_1loop_ci1

  ! [1]Eq. D.9
  FUNCTION chp_N3LO_2PE_Vc_1loop_ci0 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -gA2*gA2*const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/32D0) * ( &
              L*(2D0*mpi4(impi)*mpi4(impi)/(w*w*w*w) + 8D0*mpi3(impi)*mpi3(impi)/(w*w) - &
                  q2*q2 - 2D0*mpi4(impi) ) + mpi3(impi)*mpi3(impi)*0.5D0/(w*w) )
  END FUNCTION chp_N3LO_2PE_Vc_1loop_ci0

  ! [1]Eq. D.18
  FUNCTION chp_N3LO_2PE_Vc_2loop     (    q2,       wt2, A, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 3D0*gA2*wt2*A*const(21)*const(2)*(1D0/256D0)*( &
              (mpi2(impi) + 2D0*q2)*(2D0*chp_mpi(impi)%val + wt2*A) + &
              4D0*gA2*chp_mpi(impi)%val*wt2)
  END FUNCTION chp_N3LO_2PE_Vc_2loop     

  ! [1]Eq. D.5
  FUNCTION chp_N3LO_2PE_Wc_1loop_ci1 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -c4*q2*L*const(21)*mnuc_inv(imnuc)*(1D0/192D0) * ( &
              gA2*(8D0*mpi2(impi) + 5D0*q2) + w*w)
  END FUNCTION chp_N3LO_2PE_Wc_1loop_ci1

  ! [1]Eq. D.10
  FUNCTION chp_N3LO_2PE_Wc_1loop_ci0 (k2, q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: k2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/768D0)*( &
              L*(8D0*gA2*(1.5D0*q2*q2 + 3D0*mpi2(impi)*q2 + 3D0*mpi4(impi) - &
                      6D0*mpi3(impi)*mpi3(impi)/(w*w) - k2*(8D0*mpi2(impi) + 5D0*q2)) + &
                  4D0*gA2*gA2*(k2*(20D0*mpi2(impi) + 7D0*q2 - 16D0*mpi4(impi)/(w*w)) + &
                      16D0*mpi4(impi)*mpi4(impi)/(w*w*w*w) + &
                      12D0*mpi3(impi)*mpi3(impi)/(w*w) - 4D0*mpi4(impi)*q2/(w*w) - &
                      5D0*q2*q2 - 6D0*mpi2(impi)*q2 - 6D0*mpi4(impi)) - 4D0*k2*w*w) + &
              16D0*gA2*gA2*mpi3(impi)*mpi3(impi)/(w*w))
  END FUNCTION chp_N3LO_2PE_Wc_1loop_ci0

  ! [1]Eq. D.20
  FUNCTION chp_N3LO_2PE_Wc_2loop_a   (    q2, w, L, wt2,    impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    ! const(21) = 1 / (pi^2 * f_pi^4)
    res = L*const(21)*fpi_inv*fpi_inv*pi_inv*pi_inv*(1D0/18432D0)*( &
              192D0*pi2*fpi2*w*w*d3*(2D0*gA2*wt2 - 0.6D0*(gA2 - 1D0)*w*w) + &
              (6D0*gA2*wt2 - (gA2 - 1D0)*w*w)*(384D0*pi2*fpi2*(wt2*d1_plus_d2 + &
                      4D0*mpi2(impi)*d5) + L*(4D0*mpi2(impi)*(1D0+2D0*gA2) + &
                          q2*(1D0 + 5D0*gA2)) - (q2*(1D0/3D0)*(5D0 + 13D0*gA2) + &
                              8D0*mpi2(impi)*(1D0 + 2D0*gA2))))
  END FUNCTION chp_N3LO_2PE_Wc_2loop_a   

  ! Helper function needed by chp_calculate_2PE_2loop_int_*_data functions
  FUNCTION chp_N3LO_2PE_2loop_int_helper(x) RESULT(res)
    REAL(8), INTENT(IN) :: x
    REAL(8) :: x2, x4, x2_inv, A
    REAL(8) :: res

    IF(x < 0.05) THEN
      ! Use Taylor expansion around 0 with terms up to x^6 for small x
      x2 = x*x
      x4 = x2*x2

      res = (((-8.0D0/315.0D0)*x4*x2 + (2.0D0/35.0D0)*x4) - 0.2*x2) - 4.0D0/3.0D0
    ELSE
      x2_inv = 1.0D0 / (x*x)
      A = sqrt(1 + x2_inv)

      res = x2_inv - (1+x2_inv)*A * log(x * (1 + A))
    END IF
  END FUNCTION chp_N3LO_2PE_2loop_int_helper

  ! [1]Eq. D.21 and D.22
  SUBROUTINE chp_calculate_2PE_2loop_int_WC_data
    INTEGER :: z_nr, x_nr, i, j
    REAL(8) :: fact, C1
    REAL(8) :: res, z, z2, zt, D1
    REAL(8) :: w, x, y, y2, y_z

    IF(chp_2PE_2loop_int_WC_data_set) return

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    x_nr = chp_2PE_2loop_int_mesh_x%info%amount
    fact = 1.0D0 / (2048.0D0 * pi2 * pi2)
    C1 = 2.0D0 * gA2 * gA2 / 3.0D0

!$OMP parallel do schedule(static) default(none) shared(z_nr, x_nr, fact, C1, chp_2PE_2loop_int_mesh_z, gA2, chp_2PE_2loop_int_mesh_x, chp_2PE_2loop_int_WC_x_fact) private(i, j, res, z, z2, zt, D1, w, x, y, y2, y_z)
    DO i = 1, z_nr
      res = 0
      z = chp_2PE_2loop_int_mesh_z%pnt(i)%x
      z2 = z * z
      zt = sqrt(1 - z2)
      D1 = gA2 * (2 - z2)

      DO j = 1, x_nr
        w = chp_2PE_2loop_int_mesh_x%pnt(j)%w
        x = chp_2PE_2loop_int_mesh_x%pnt(j)%x
        y = x * zt
        y2 = y * y
        y_z = y / z

        res = res + fact * w * zt * z * ((gA2 - 1) * y2 - D1) * (-y2 + 2*y * sqrt(z2 + y2) * log(y_z + sqrt(1 + y_z*y_z)) + C1 * (2 - y2 - z2) * (chp_N3LO_2PE_2loop_int_helper(y_z) + 5.0D0/6.0D0))
      END DO

      chp_2PE_2loop_int_WC_x_fact(i) = res
    END DO
!$OMP end parallel do

    chp_2PE_2loop_int_WC_data_set = .TRUE.
  END SUBROUTINE chp_calculate_2PE_2loop_int_WC_data

  FUNCTION chp_N3LO_2PE_Wc_2loop_b   (    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    INTEGER :: i
    REAL(8) :: fact, a1
    
    if(chp_use_2PE_2loop_int%val == 0) then
      res = 0
      return
    end if

    call chp_calculate_2PE_2loop_int_WC_data

    fact = fpi_inv**6
    a1 = 4.0D0 * mpi2(impi)

!$OMP parallel do schedule(static) default(none) shared(q2, impi, res, fact, a1, chp_2PE_2loop_int_mesh_z, chp_2PE_2loop_int_WC_x_fact) private(i)
    DO i = 1, nof_theta_int_points
      res(i) = fact * q2(i)**3 * sum(chp_2PE_2loop_int_mesh_z%pnt(:)%w * chp_2PE_2loop_int_WC_x_fact(:) / (a1 + chp_2PE_2loop_int_mesh_z%pnt(:)%x**2 * q2(i)))
    END DO
!$OMP end parallel do
  END FUNCTION chp_N3LO_2PE_Wc_2loop_b   

  ! [1]Eq. D.7
  FUNCTION chp_N3LO_2PE_VLS_1loop_ci1(        w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (c2*gA2*const(21)*mnuc_inv(imnuc)*(1D0/8D0))*w*w*L
  END FUNCTION chp_N3LO_2PE_VLS_1loop_ci1

  ! [1]Eq. D.13
  FUNCTION chp_N3LO_2PE_VLS_1loop_ci0(    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (gA2*gA2*const(21)*0.25D0*mnuc_inv(imnuc)*mnuc_inv(imnuc))*L*( &
              (11D0/32D0)*q2 + mpi4(impi)/(w*w))
  END FUNCTION chp_N3LO_2PE_VLS_1loop_ci0

  ! [1]Eq. D.8
  FUNCTION chp_N3LO_2PE_WLS_1loop_ci1(    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (-c4*const(21)*mnuc_inv(imnuc)*(1D0/48D0))*L*( &
              gA2*(8D0*mpi2(impi) + 5D0*q2) + w*w)
  END FUNCTION chp_N3LO_2PE_WLS_1loop_ci1

  ! [1]Eq. D.14
  FUNCTION chp_N3LO_2PE_WLS_1loop_ci0(    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/256D0))*L*( &
              16D0*gA2*(mpi2(impi) + 0.375D0*q2) + (4D0/3D0)*gA2*gA2*( &
                  4D0*mpi4(impi)/(w*w) - (11D0/4D0)*q2 - 9D0*mpi2(impi)) - w*w)
  END FUNCTION chp_N3LO_2PE_WLS_1loop_ci0

  ! [1]Eq. D.11
  FUNCTION chp_N3LO_2PE_VT_1loop_ci0 (k2, q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: k2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (gA2*gA2*const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/32D0))*L*( &
              k2 + 0.625D0*q2 + mpi4(impi)/(w*w))
  END FUNCTION chp_N3LO_2PE_VT_1loop_ci0

  ! [1]Eq. D.24
  FUNCTION chp_N3LO_2PE_VT_2loop_a   (        w, L,             impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (-gA2*const(21)*(1D0/32D0)*d14_minus_d15)*w*w*L
  END FUNCTION chp_N3LO_2PE_VT_2loop_a   

  ! [1]Eq. D.25
  SUBROUTINE chp_calculate_2PE_2loop_int_VTS_data
    INTEGER :: z_nr, x_nr, i, j
    REAL(8) :: fact
    REAL(8) :: res, z, zt
    REAL(8) :: w, x

    IF(chp_2PE_2loop_int_VTS_data_set) return

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    x_nr = chp_2PE_2loop_int_mesh_x%info%amount
    fact = -1.0D0 / (1024.0D0 * pi2 * pi2)

!$OMP parallel do schedule(static) default(none) shared(z_nr, x_nr, fact, chp_2PE_2loop_int_mesh_z, chp_2PE_2loop_int_mesh_x, chp_2PE_2loop_int_VTS_x_fact) private(i, j, res, z, zt, w, x)
    DO i = 1, z_nr
      res = 0
      z = chp_2PE_2loop_int_mesh_z%pnt(i)%x
      zt = sqrt(1 - z*z)

      DO j = 1, x_nr
        w = chp_2PE_2loop_int_mesh_x%pnt(j)%w
        x = chp_2PE_2loop_int_mesh_x%pnt(j)%x

        res = res + fact * w * (1 - x*x) * z * zt * (1 - z*z) * (chp_N3LO_2PE_2loop_int_helper(x*zt/z) - 1.0D0/6.0D0)
      END DO

      chp_2PE_2loop_int_VTS_x_fact(i) = res
    END DO
!$OMP end parallel do

    chp_2PE_2loop_int_VTS_data_set = .TRUE.
  END SUBROUTINE chp_calculate_2PE_2loop_int_VTS_data

  FUNCTION chp_N3LO_2PE_VT_2loop_b   (    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    INTEGER :: i
    REAL(8) :: fact, a1
    
    if(chp_use_2PE_2loop_int%val == 0) then
      res = 0
      return
    end if

    call chp_calculate_2PE_2loop_int_VTS_data

    fact = (chp_gA%val * fpi_inv)**6
    a1 = 4.0D0 * mpi2(impi)

!$OMP parallel do schedule(static) default(none) shared(q2, impi, res, fact, a1, chp_2PE_2loop_int_mesh_z, chp_2PE_2loop_int_VTS_x_fact) private(i)
    DO i = 1, nof_theta_int_points
      res(i) = fact * q2(i)**2 * sum(chp_2PE_2loop_int_mesh_z%pnt(:)%w * chp_2PE_2loop_int_VTS_x_fact(:) / (a1 + chp_2PE_2loop_int_mesh_z%pnt(:)%x**2 * q2(i)))
    END DO
!$OMP end parallel do
  END FUNCTION chp_N3LO_2PE_VT_2loop_b   

  ! [1]Eq. D.2
  FUNCTION chp_N3LO_2PE_WT_1loop_ci2 (        w, L,             impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (c4*c4*const(21)*(1D0/96D0))*w*w*L
  END FUNCTION chp_N3LO_2PE_WT_1loop_ci2

  ! [1]Eq. D.6
  FUNCTION chp_N3LO_2PE_WT_1loop_ci1 (q2,     w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (-c4*const(21)*mnuc_inv(imnuc)*(1D0/192D0))*L*( &
              gA2*(16D0*mpi2(impi) + 7D0*q2) - w*w)
  END FUNCTION chp_N3LO_2PE_WT_1loop_ci1

  ! [1]Eq. D.12
  FUNCTION chp_N3LO_2PE_WT_1loop_ci0 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (mnuc_inv(imnuc)*mnuc_inv(imnuc)*const(21)*(1D0/1536D0))*L*( &
              4D0*gA2*gA2*(7D0*mpi2(impi) + (17D0/4D0)*q2 + 4D0*mpi4(impi)/(w*w)) - &
              32D0*gA2*(mpi2(impi) + (7D0/16D0)*q2) + w*w)
  END FUNCTION chp_N3LO_2PE_WT_1loop_ci0

  ! [1]Eq. D.27
  FUNCTION chp_N3LO_2PE_WT_2loop     (        w,             A, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (gA2*gA2*const(21)*fpi_inv*fpi_inv*(1D0/2048D0))*w*w*A*( &
              w*w*A + 2D0*chp_mpi(impi)%val*(1D0 + 2D0*gA2))
  END FUNCTION chp_N3LO_2PE_WT_2loop     

  ! [1]Eq. D.11
  FUNCTION chp_N3LO_2PE_Vs_1loop_ci0 (k2, q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: k2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_VT_1loop_ci0(k2, q2, w, L, impi, imnuc)
  END FUNCTION chp_N3LO_2PE_Vs_1loop_ci0

  ! [1]Eq. D.24
  FUNCTION chp_N3LO_2PE_Vs_2loop_a   (    q2, w, L,             impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_VT_2loop_a(w, L, impi)
  END FUNCTION chp_N3LO_2PE_Vs_2loop_a   

  ! [1]Eq. D.25
  FUNCTION chp_N3LO_2PE_Vs_2loop_b   (    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_VT_2loop_b(q2, impi)
  END FUNCTION chp_N3LO_2PE_Vs_2loop_b   

  ! [1]Eq. D.2
  FUNCTION chp_N3LO_2PE_Ws_1loop_ci2 (    q2, w, L,             impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_WT_1loop_ci2(w, L, impi)
  END FUNCTION chp_N3LO_2PE_Ws_1loop_ci2

  ! [1]Eq. D.6
  FUNCTION chp_N3LO_2PE_Ws_1loop_ci1 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_WT_1loop_ci1(q2, w, L, impi, imnuc)
  END FUNCTION chp_N3LO_2PE_Ws_1loop_ci1

  ! [1]Eq. D.12
  FUNCTION chp_N3LO_2PE_Ws_1loop_ci0 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_WT_1loop_ci0(q2, w, L, impi, imnuc)
  END FUNCTION chp_N3LO_2PE_Ws_1loop_ci0

  ! [1]Eq. D.27
  FUNCTION chp_N3LO_2PE_Ws_2loop     (    q2, w,             A, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_WT_2loop(w, A, impi)
  END FUNCTION chp_N3LO_2PE_Ws_2loop     

  ! [1]Eq. D.15
  FUNCTION chp_N3LO_2PE_VsigL_1loop_ci0(         L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (gA2*gA2*const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/32D0))*L
  END FUNCTION chp_N3LO_2PE_VsigL_1loop_ci0


  SUBROUTINE chp_LO_contact_terms(chn, pfinal, pinit, POT)
    
    IMPLICIT NONE

    TYPE(chp_chn_type), INTENT(IN)    :: chn
    REAL(8)            , INTENT(INOUT) :: POT(6)
    REAL(8) , intent(in) :: pfinal, pinit
    REAL(8) :: p, pp

    p = pinit
    pp = pfinal

    !1S0
    IF (CHN%J == 0 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_CLO(1)%val * freg(pp,p,chp_CLO_n(1)%val)
    ELSE IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       !3S1
       POT(4)  = POT(4) + chp_CLO(2)%val * freg(pp,p,chp_CLO_n(2)%val)
    END IF

  END SUBROUTINE chp_LO_contact_terms
  
  SUBROUTINE chp_CIB_LO_contact_terms(chn, pfinal, pinit, POT)
    
    IMPLICIT NONE

    TYPE(chp_chn_type), INTENT(IN)    :: chn
    REAL(8)            , INTENT(INOUT) :: POT(6)
    REAL(8) , intent(in) :: pfinal, pinit
    REAL(8) :: p, pp
    
    p = pinit
    pp = pfinal

    ! CIB contacts
    !1S0
    IF (CHN%J == 0 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_CIB_CLO(CHN%tz,1)%val * freg(pp,p,chp_CLO_n(1)%val)
       !3S1
    ELSE IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(4)  = POT(4) + chp_CIB_CLO(CHN%tz,2)%val * freg(pp,p,chp_CLO_n(2)%val)
    END IF
    
  END SUBROUTINE chp_CIB_LO_contact_terms

  SUBROUTINE chp_NLO_contact_terms(chn, pfinal, pinit, POT)
    
    IMPLICIT NONE

    TYPE(chp_chn_type), INTENT(IN)    :: chn
    REAL(8)            , INTENT(IN)    :: pfinal, pinit
    REAL(8)            , INTENT(INOUT) :: POT(6)
    REAL(8) :: p, pp, ppp
    
    p  = pinit
    pp = pfinal
    ppp=p*pp

    ! NLO contacts
    !1S0
    IF (CHN%J == 0 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_CNLO(1)%val*(p*p+pp*pp) * freg(pp,p,chp_CNLO_n(1)%val)
    ELSE IF (CHN%J == 0 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       !3P0
       POT(3)  = POT(3) + chp_CNLO(2)%val*ppp * freg(pp,p,chp_CNLO_n(2)%val)
    ELSE IF (CHN%J == 1 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       !1P1
       POT(1) = POT(1) + chp_CNLO(3)%val*ppp * freg(pp,p,chp_CNLO_n(3)%val)
    ELSE IF (CHN%J == 1 .AND. CHN%S == 1 .AND. .NOT. CHN%coup) THEN
       !3P1
       POT(2) = POT(2) + chp_CNLO(4)%val*ppp * freg(pp,p,chp_CNLO_n(4)%val)
    ELSE IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       !3S1
       POT(4)  = POT(4) + chp_CNLO(5)%val*(p*p+pp*pp) * freg(pp,p,chp_CNLO_n(5)%val)
       !3S1-3D1
       POT(6)  = POT(6) + chp_CNLO(6)%val*pp*pp * freg(pp,p,chp_CNLO_n(6)%val)
       !3D1-3S1
       POT(5)  = POT(5) + chp_CNLO(6)%val*p*p * freg(pp,p,chp_CNLO_n(6)%val)
    ELSE IF (CHN%J == 2 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       !3P2
       POT(4)  = POT(4) + chp_CNLO(7)%val*p*pp * freg(pp,p,chp_CNLO_n(7)%val)
    END IF

  END SUBROUTINE chp_NLO_contact_terms

  SUBROUTINE chp_N3LO_contact_terms(chn, pfinal, pinit, POT)
    
    IMPLICIT NONE

    TYPE(chp_chn_type), INTENT(IN)    :: chn
    REAL(8)            , INTENT(IN)    :: pfinal, pinit
    REAL(8)            , INTENT(INOUT) :: POT(6)
    REAL(8) :: p, pp, ppp, p2, pp2
    
    p  = pinit
    p2 = p * p
    pp = pfinal
    pp2 = pp * pp
    ppp=p*pp

    ! N3LO contacts
    ! Eq [1] E.1, the contacts are numbered in the order they appear
    ! in the equation [1] E.1.
    ! FIXME: borisc: p and pp seems to be switched compared with what [1] uses. I do it this
    !                way to be consistent with chp_NLO_contact_terms
    !                UPDATE: On closer inspection, it seems that the function is called with
    !                        reversed arguments, so that pinit <-> pfinal, so everything is OK
    
    !1S0
    IF (CHN%J == 0 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_DN3LO(1)%val*(p2 * p2 + pp2 * pp2) * freg(pp,p,chp_DN3LO_n(1)%val) + chp_DN3LO(2)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(2)%val)
    END IF
    
    !3P0
    IF (CHN%J == 0 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(3) = POT(3) + chp_DN3LO(3)%val*ppp*(p2 + pp2) * freg(pp,p,chp_DN3LO_n(3)%val)
    END IF
    
    !1P1
    IF (CHN%J == 1 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_DN3LO(4)%val*ppp*(p2 + pp2) * freg(pp,p,chp_DN3LO_n(4)%val)
    END IF
    
    !3P1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. .NOT. CHN%coup) THEN
       POT(2) = POT(2) + chp_DN3LO(5)%val*ppp*(p2 + pp2) * freg(pp,p,chp_DN3LO_n(5)%val)
    END IF
    
    !3S1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(4) = POT(4) + chp_DN3LO(6)%val*(p2 * p2 + pp2 * pp2) * freg(pp,p,chp_DN3LO_n(6)%val) + chp_DN3LO(7)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(7)%val)
    END IF
    
    !3D1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(3) = POT(3) + chp_DN3LO(8)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(8)%val)
    END IF
    
    !3S1-3D1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(6) = POT(6) + chp_DN3LO(9)%val*(pp2 * pp2) * freg(pp,p,chp_DN3LO_n(9)%val) + chp_DN3LO(10)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(10)%val)
    END IF
    
    !3D1-3S1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(5) = POT(5) + chp_DN3LO(9)%val*(p2 * p2) * freg(pp,p,chp_DN3LO_n(9)%val) + chp_DN3LO(10)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(10)%val)
    END IF
    
    !1D2
    IF (CHN%J == 2 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_DN3LO(11)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(11)%val)
    END IF
    
    !3D2
    IF (CHN%J == 2 .AND. CHN%S == 1 .AND. .NOT. CHN%coup) THEN
       POT(2) = POT(2) + chp_DN3LO(12)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(12)%val)
    END IF
    
    !3P2
    IF (CHN%J == 2 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(4) = POT(4) + chp_DN3LO(13)%val*ppp*(p2 + pp2) * freg(pp,p,chp_DN3LO_n(13)%val)
    END IF
    
    !3P2-3F2
    IF (CHN%J == 2 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(6) = POT(6) + chp_DN3LO(14)%val*ppp*(pp2) * freg(pp,p,chp_DN3LO_n(14)%val)
    END IF
    
    !3F2-3P2
    IF (CHN%J == 2 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(5) = POT(5) + chp_DN3LO(14)%val*ppp*(p2) * freg(pp,p,chp_DN3LO_n(14)%val)
    END IF
    
    !3D3
    IF (CHN%J == 3 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       IF (chp_DN3LO_n(15)%val < 0.0D0) THEN
          POT(4) = POT(4) + chp_DN3LO(15)%val*(p2 * pp2) * 0.5D0 * (freg(pp,p,2.0D0) + freg(pp,p,3.0D0))
       ELSE
          POT(4) = POT(4) + chp_DN3LO(15)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(15)%val)
       END IF
    END IF

  END SUBROUTINE chp_N3LO_contact_terms


  ! cm energy defined on p.28
  FUNCTION cm_energy(p,tz) RESULT(res)

    REAL(8), INTENT(IN)  :: p
    INTEGER, INTENT(IN) :: tz
    REAL(8) :: res

    res = DSQRT(mnuc2(tz) + p*p)

  END FUNCTION cm_energy

  ! regulator, eq 4.63
  ! pfinal : final momentum
  ! pinit  : initial momentum
  ! n      : cutoff order
  ! LAMBDA is accessed from the chp constant chp_lambda
  FUNCTION freg(pfinal, pinit, n) RESULT(res)
    
    REAL(8) , INTENT(IN) :: pfinal, pinit
    REAL(8) , INTENT(IN) :: n
    REAL(8) :: res,exponent,lambda

    res = 0.0D0

    lambda = chp_lambda%val
    
    
    exponent = (pfinal/lambda)**(2.0D0*n) + &
               (pinit/lambda)**(2.0D0*n)
    
    res = dexp(-exponent)
    
  END FUNCTION freg

  SUBROUTINE initialize_chiral_potential
    
    CALL set_contact_conversion_matrices
    CALL setup_twobody_pwd(nof_theta_int_points, maximum_angular_momentum)
    CALL chp_setup_2PE_2loop_int_data(nof_2PE_2loop_int_x_points, nof_2PE_2loop_int_z_points)
    chp_mnuc(-1)     = chp_real_type(rirrel,'m_prot' ,.FALSE.)
    chp_mnuc(0)      = chp_real_type(rirrel,'m_nucl' ,.FALSE.)
    chp_mnuc(+1)     = chp_real_type(rirrel,'m_neut' ,.FALSE.)
    chp_mpi(-1)      = chp_real_type(rirrel,'m_pi-'  ,.FALSE.)
    chp_mpi(0)       = chp_real_type(rirrel,'m_pi0'  ,.FALSE.)
    chp_mpi(+1)      = chp_real_type(rirrel,'m_pi+'  ,.FALSE.)
    chp_mpi(+2)      = chp_real_type(rirrel,'m_pi '  ,.FALSE.)
    chp_chiral_order = chp_int_type(iirrel, 'nu'     , .FALSE.)
    chp_ren          = chp_real_type(rirrel, 'REN'   , .FALSE.)
    chp_lambda       = chp_real_type(rirrel, 'LAMBDA', .FALSE.)
    chp_regcut_1PE   = chp_real_type(rirrel, 'rcut 1', .FALSE.)
    chp_regcut_2PE   = chp_real_type(rirrel, 'rcut 2', .FALSE.)
    chp_itope        = chp_char2_type(cirrel, 'itope' , .FALSE.)
    chp_gA           = chp_real_type(rirrel, 'gA'    , .FALSE.)
    chp_fpi          = chp_real_type(rirrel, 'fpi'   , .FALSE.)
    chp_fine_structure = chp_real_type(rirrel, 'alpha'   , .FALSE.)
    chp_2PE_CSB_correct_mass = chp_int_type(0, '2PE CSB mass', .TRUE.)
    chp_use_2PE_2loop_int = chp_int_type(0, '2PE 2loop int', .TRUE.)
    chp_2PE_2loop_int_VTS_data_set = .FALSE.
    chp_2PE_2loop_int_WC_data_set  = .FALSE.
    
    ! N2LO LECs c1 c3 c4
    chp_c1           = chp_real_type(rirrel, 'c1'   , .FALSE.)
    chp_c3           = chp_real_type(rirrel, 'c3'   , .FALSE.)
    chp_c4           = chp_real_type(rirrel, 'c4'   , .FALSE.)

    ! N3LO LECs c2, d1+d2, d3, d5, d14-d15
    chp_c2            = chp_real_type(rirrel, 'c2'     , .FALSE.)
    chp_d1_plus_d2    = chp_real_type(rirrel, 'd1+d2'  , .FALSE.)
    chp_d3            = chp_real_type(rirrel, 'd3'     , .FALSE.)
    chp_d5            = chp_real_type(rirrel, 'd5'     , .FALSE.)
    chp_d14_minus_d15 = chp_real_type(rirrel, 'd14-d15', .FALSE.)

    ! LO contacts
    chp_CLO(1)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CLO(2)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CLO_n(1)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CLO_n(2)     = chp_real_type(rirrel, 'x', .FALSE.)
    
    ! NLO contacts
    chp_CIB_CLO(-1,1)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CIB_CLO(-1,2)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CIB_CLO( 0,1)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CIB_CLO( 0,2)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CIB_CLO(+1,1)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CIB_CLO(+1,2)       = chp_real_type(rirrel, 'x', .FALSE.)
    
    chp_CNLO(1)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO(2)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO(3)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO(4)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO(5)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO(6)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO(7)       = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(1)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(2)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(3)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(4)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(5)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(6)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(7)     = chp_real_type(rirrel, 'x', .FALSE.)

    ! N3LO contacts
    chp_DN3LO( 1) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO( 2) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO( 3) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO( 4) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO( 5) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO( 6) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO( 7) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO( 8) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO( 9) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO(10) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO(11) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO(12) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO(13) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO(14) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO(15) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 1) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 2) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 3) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 4) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 5) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 6) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 7) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 8) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 9) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(10) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(11) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(12) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(13) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(14) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(15) = chp_real_type(rirrel, 'x', .FALSE.)

  END SUBROUTINE initialize_chiral_potential
  
  SUBROUTINE chp_set_units_and_derive_constants
    
    USE chp_aux
    
    implicit none
    INTEGER :: i, row
    REAL(8), ALLOCATABLE, DIMENSION(:) :: tmp_contacts
    REAL(8) :: reg_val

    ! set the regulator cutoffs depending on the chiral order
    IF (.NOT. chp_regcut_1PE%set) THEN
       IF (chp_chiral_order%val == N3LO) THEN
          CALL chp_set_1PE_reg_par(4.0D0)
       ELSE
          CALL chp_set_1PE_reg_par(3.0D0)
       END IF
    END IF
    IF (.NOT. chp_regcut_2PE%set) THEN
       IF (chp_chiral_order%val == N3LO) THEN
          CALL chp_set_2PE_reg_par(2.0D0)
       ELSE
          CALL chp_set_2PE_reg_par(3.0D0)
       END IF
    END IF

    ! Set regulator parameter for LO contacts if not set
    DO i = 1, 2
       IF(chp_CLO_n(i)%set) CYCLE
       CALL chp_set_LO_contact_reg_par(i, 3.0D0)
    END DO
    ! Set regulator parameter for NLO contacts if not set
    IF (chp_chiral_order%val > LO) THEN
       DO i = 1, 7
          IF(chp_CNLO_n(i)%set) CYCLE
          IF(chp_chiral_order%val == N3LO) THEN
             reg_val = 2.0D0
          ELSE
             reg_val = 3.0D0
          END IF
          CALL chp_set_NLO_contact_reg_par(i, reg_val)
       END DO
    END IF
    IF (chp_chiral_order%val > NNLO) THEN
       DO i = 1, 15
          IF(chp_DN3LO_n(i)%set) CYCLE
          SELECT CASE(i)
          CASE(1)
             reg_val = 2.0D0
          CASE(2)
             reg_val = 2.0D0
          CASE(3)
             reg_val = 3.0D0
          CASE(4)
             reg_val = 2.0D0
          CASE(5)
             reg_val = 4.0D0
          CASE(6)
             reg_val = 2.0D0
          CASE(7)
             reg_val = 2.0D0
          CASE(8)
             reg_val = 2.0D0
          CASE(9)
             reg_val = 2.0D0
          CASE(10)
             reg_val = 2.0D0
          CASE(11)
             reg_val = 4.0D0
          CASE(12)
             reg_val = 2.0D0
          CASE(13)
             reg_val = 2.0D0
          CASE(14)
             reg_val = 4.0D0
          CASE(15)
             reg_val = -1.0D0
          END SELECT
          CALL chp_set_N3LO_contact_reg_par(i, reg_val)
       END DO
    END IF

    ! NNLO c1 c3 c4
    IF (chp_chiral_order%val > NLO) THEN

       c1 = chp_c1%val
       c3 = chp_c3%val
       c4 = chp_c4%val

       IF (.NOT. chp_c1%set) THEN
          WRITE(CHP_ERR,"(A)") 'for chiral_order>NLO you must set c1'
          STOP
       END IF

       IF (.NOT. chp_c3%set) THEN
          WRITE(CHP_ERR,"(A)") 'for chiral_order>NLO you must set c3'
          STOP
       END IF

       IF (.NOT. chp_c4%set) THEN
          WRITE(CHP_ERR,"(A)") 'for chiral_order>NLO you must set c4'
          STOP
       END IF
       
       IF (chp_chiral_order%val > NNLO) THEN
          c2 = chp_c2%val
          d1_plus_d2 = chp_d1_plus_d2%val
          d3 = chp_d3%val
          d5 = chp_d5%val
          d14_minus_d15 = chp_d14_minus_d15%val

          IF (.NOT. chp_c2%set) THEN
             WRITE(CHP_ERR,"(A)") 'for chiral_order>NNLO you must set c2'
             STOP
          END IF

          IF (.NOT. chp_d1_plus_d2%set) THEN
             WRITE(CHP_ERR,"(A)") 'for chiral_order>NNLO you must set d1+d2'
             STOP
          END IF

          IF (.NOT. chp_d3%set) THEN
             WRITE(CHP_ERR,"(A)") 'for chiral_order>NNLO you must set d3'
             STOP
          END IF

          IF (.NOT. chp_d5%set) THEN
             WRITE(CHP_ERR,"(A)") 'for chiral_order>NNLO you must set d5'
             STOP
          END IF

          IF (.NOT. chp_d14_minus_d15%set) THEN
             WRITE(CHP_ERR,"(A)") 'for chiral_order>NNLO you must set d14-d15'
             STOP
          END IF
        END IF
    END IF

    mnuc2(:)   = chp_mnuc(:)%val*chp_mnuc(:)%val
    mnuc_inv(:)= 1D0 / chp_mnuc(:)%val
    mpi2(:)    = chp_mpi(:)%val*chp_mpi(:)%val
    mpi3(:)    = chp_mpi(:)%val*mpi2(:)
    mpi4(:)    = mpi2*mpi2
    mpi5(:)    = mpi4*chp_mpi(:)%val
    twompi(:)  = 2.0D0*chp_mpi(:)%val
    fourmpi2(:)= 4.0D0*mpi2
    
    IF (chp_ren%name == 'SF') THEN
       
       sfr  = chp_ren%val
       sfr2 = sfr*sfr

       DO i=-1,2
          IF ( (sfr-twompi(i)) <  0.0D0 ) sfr_heavyside(i) = 0.0D0
          IF ( (sfr-twompi(i)) >= 0.0D0 ) sfr_heavyside(i) = 1.0D0
       END DO

    END IF
    
    gA2        = chp_gA%val*chp_gA%val
    gA4        = gA2*gA2
    
    fpi2       = chp_fpi%val*chp_fpi%val
    fpi4       = fpi2*fpi2
    fpi_inv    = 1D0 / chp_fpi%val

    iso(0) = -3.0D0
    iso(1) = +1.0D0
    
    const = 0.0D0
    !1: 1/(2pi)^3
    const(1) = 1.0D0/(twopi**3)
    !2: gA^2/(4*fpi^2)
    const(2) = gA2/(4.0D0*fpi2)
    !3: 384pi^2*fpi^4
    const(3) = 384.0D0*pi2*fpi4
    !4: 5gA^4-4gA^2-1
    const(4) = 5.0D0*gA4-4.0D0*gA2-1.0D0
    !5: 23gA^4-10gA^2-1
    const(5) = 23.0D0*gA4-10.0D0*gA2 - 1.0D0
    !6: 48gA^4
    const(6) = 48.0D0*gA4
    !7: 3gA^4
    const(7) = 3.0D0*gA4
    !8: 64pi^2fpi^4
    const(8) = 64.0D0*pi2*fpi4
    !9: 3gA^2/16pifpi^4
    const(9) = 3.0D0*gA2/(16.0D0*chp_pi%val*fpi4)
    !10: ga^2/16
    const(10) = ga2/16.0D0
    !11: 2.0D0*(2c1-c3)
    const(11) = 2.0D0*(2.0D0*c1-c3)
    !12 : const(7)/256pifpi^4
    const(12) = const(7)/(256.0D0*chp_pi%val*fpi4)
    !13: gA^2/(128pifpi^4)
    const(13) = gA2/(128.0D0*chp_pi%val*fpi4)
    !14: gA^4/(128pifpi^4)
    const(14) = gA4/(128.0D0*chp_pi%val*fpi4)
    !15: 3gA^4/(512pifpi^4)
    const(15) = 3.0D0*gA4/(512.0D0*chp_pi%val*fpi4)
    !16: gA2/(32pifpi^4)
    const(16) = gA2/(32.0D0*chp_pi%val*fpi4)
    !17: gA2/8
    const(17) = gA2/8.0D0
    !18: gA4/(256pifpi^4)
    const(18) = gA4/(256.0D0*chp_pi%val*fpi4)
    !19: 3gA4/(32pifpi^4)
    const(19) = 3.0D0*gA4/(32.0D0*chp_pi%val*fpi4)
    !20: const(16)*(1-gA2)
    const(20) = const(16)*(1.0D0-gA2)
    !21: 1 / (pi^2fpi^4)
    const(21) = 1D0 / (pi2 * fpi4)
    ! -----------------------------------------------------
    
    ! if the contacts were input in PW, transform to ST and print
    IF (chp_contact_format%name == 'PW') THEN
       
       !WRITE(CHP_SCR,*)
       !WRITE(CHP_SCR,"(A)") '   *************************'
       !WRITE(CHP_SCR,"(A)") '   CONTACTS INPUT FORMAT: PW'
       !WRITE(CHP_SCR,"(A)") '   IN THE ST FORMAT THEY ARE:'
       !WRITE(CHP_SCR,"(A)") '   *************************'
       !WRITE(CHP_SCR,*)
       
       ! LO CONTACTS
       IF (chp_chiral_order%val == LO) THEN
          
          IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
          ALLOCATE(tmp_contacts(1:2))
          tmp_contacts = 0.0D0

          !simple matrix vector multiplication
          
          DO row=1, 2
             tmp_contacts(row) = SUM(LOPW2ST(row,:)*chp_CLO(:)%val)
          END DO
          
          !WRITE(CHP_SCR,"(A12,F30.16)") 'CS', tmp_contacts(1)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'CT', tmp_contacts(2)

       END IF

       ! NLO CONTACTS
       IF (chp_chiral_order%val > LO) THEN
          
          ! do the CIB LO contacts
          DO i=-1,1
             
             IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
             ALLOCATE(tmp_contacts(1:2))
             tmp_contacts = 0.0D0
             
             !simple matrix vector multiplication
             
             DO row=1, 2
                tmp_contacts(row) = SUM(LOPW2ST(row,:)*chp_CIB_CLO(i,:)%val)
             END DO

             SELECT CASE (i)
             CASE(-1)
                !WRITE(CHP_SCR,"(A12,F30.16)") 'CSpp', tmp_contacts(1)
                !WRITE(CHP_SCR,"(A12,F30.16)") 'CTpp', tmp_contacts(2)
             CASE( 0)
                !WRITE(CHP_SCR,"(A12,F30.16)") 'CSpn', tmp_contacts(1)
                !WRITE(CHP_SCR,"(A12,F30.16)") 'CTpn', tmp_contacts(2)
             CASE(+1)
                !WRITE(CHP_SCR,"(A12,F30.16)") 'CSnn', tmp_contacts(1)
                !WRITE(CHP_SCR,"(A12,F30.16)") 'CTnn', tmp_contacts(2)
             END SELECT
          END DO

          ! do the NLO contacts
          
          IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
          ALLOCATE(tmp_contacts(1:7))
          tmp_contacts = 0.0D0
          
          !simple matrix vector multiplication
          
          DO row=1,7
             tmp_contacts(row) = SUM(NLOPW2ST(row,:)*chp_CNLO(:)%val)
          END DO
          
          !WRITE(CHP_SCR,"(A12,F30.16)") 'C1', tmp_contacts(1)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'C2', tmp_contacts(2)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'C3', tmp_contacts(3)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'C4', tmp_contacts(4)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'C5', tmp_contacts(5)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'C6', tmp_contacts(6)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'C7', tmp_contacts(7)
          
       END IF
       
       IF (chp_chiral_order%val > NNLO) THEN
          ! do the N3LO contacts
          
          IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
          ALLOCATE(tmp_contacts(1:15))
          tmp_contacts = 0.0D0
          
          !simple matrix vector multiplication
          
          DO row=1,15
             tmp_contacts(row) = SUM(N3LOPW2ST(row,:)*chp_DN3LO(:)%val)
          END DO
          
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D1' , tmp_contacts( 1)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D2' , tmp_contacts( 2)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D3' , tmp_contacts( 3)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D4' , tmp_contacts( 4)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D5' , tmp_contacts( 5)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D6' , tmp_contacts( 6)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D7' , tmp_contacts( 7)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D8' , tmp_contacts( 8)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D9' , tmp_contacts( 9)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D10', tmp_contacts(10)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D11', tmp_contacts(11)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D12', tmp_contacts(12)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D13', tmp_contacts(13)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D14', tmp_contacts(14)
          !WRITE(CHP_SCR,"(A12,F30.16)") 'D15', tmp_contacts(15)
        END IF
       
    END IF
       
    ! if the contacts were input in ST, transform to PW replace, and print
    IF (chp_contact_format%name == 'ST') THEN
       
       WRITE(CHP_SCR,*)
       WRITE(CHP_SCR,"(A)") '   *************************'
       WRITE(CHP_SCR,"(A)") '   CONTACTS INPUT FORMAT: ST'
       WRITE(CHP_SCR,"(A)") '   IN THE PW FORMAT THEY ARE:'
       WRITE(CHP_SCR,"(A)") '   *************************'
       WRITE(CHP_SCR,*)

       ! LO
       ! FIXME: borisc: What about higher orders? They don't seem to
       !                be converted to PW format.
       IF (chp_chiral_order%val == LO) THEN
          
          ! matrix vector multiplication
          
          IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
          ALLOCATE(tmp_contacts(1:2))
          tmp_contacts = 0.0D0

          !simple matrix vector multiplication
          
          DO row=1, 2
             tmp_contacts(row) = SUM(LOST2PW(row,:)*chp_CLO(:)%val)
          END DO
          
          WRITE(CHP_SCR,"(A12,F30.16)") '1S0', tmp_contacts(1)
          WRITE(CHP_SCR,"(A12,F30.16)") '3S1', tmp_contacts(2)
          
          ! replace the ST matrix elements with st ones
          
          DO i=1,2
             chp_CLO(i)%val = tmp_contacts(i)
          END DO
          
          ! change the format identifier
          chp_contact_format%name = 'PW'
       END IF

    END IF
    
    ! convert the contact parameters to units of [MeV]
    
    ! LO
    ! the LO contacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2
    !
    !IF (chp_chiral_order%val == LO) THEN
    !   
    !   DO i=1, 2
    !      chp_CLO(i)%val = chp_CLO(i)%val * 0.01D0
    !   END DO
    !   
    !END IF

    ! NLO
    ! the LO CIBcontacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2
    ! the NLO contacts are input in units of 10^4/GeV^4
    ! this = 10^-8/MeV^4
    ! the N3LO contacts are input in units of 10^4/GeV^6
    ! this = 10^(4-3*6)/MeV^6 = 10^-14/MeV^6
    
    !IF (chp_chiral_order%val > LO) THEN
    !   
    !   DO i=1, 2
    !      DO j=-1,1
    !         chp_CIB_CLO(j,i)%val = chp_CIB_CLO(j,i)%val * 0.01D0
    !      END DO
    !   END DO
    !   
    !   DO i=1,7
    !      chp_CNLO(i)%val = chp_CNLO(i)%val * 1.E-08
    !   END DO
    !
    !END IF
    
    !IF (chp_chiral_order%val > NNLO) THEN
    !   DO i=1,15
    !      chp_DN3LO(i)%val = chp_DN3LO(i)%val * 1.D-14
    !   END DO
    !END IF


    ! -----------------------------------------------------
    
    ! DO SOME CONSISTENCY CHECKS

    !IF (chp_gA%val == 1.29D0 .AND. chp_chiral_order%val == LO) THEN
    !   WRITE(CHP_ERR,"(A)")
    !   WRITE(CHP_ERR,"(A)") '   ****************************************'
    !   WRITE(CHP_ERR,"(A)") '   AT LO, YOU DONT HAVE TO CORRECT FOR THE '
    !   WRITE(CHP_ERR,"(A)") '   GOLDBERGER-TREIMAN DISCREPANCY. CHANGE'
    !   WRITE(CHP_ERR,"(A)") '   gA FROM 1.29 TO 1.276'
    !   WRITE(CHP_ERR,"(A)") '   ****************************************'
    !   WRITE(CHP_ERR,"(A)")
    !END IF
    
    IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)

  END SUBROUTINE chp_set_units_and_derive_constants

  SUBROUTINE chp_print_constants(title,unit)
    
    INTEGER :: unit,i,j
    CHARACTER(LEN=42) :: title

    WRITE(unit,"(A42)") '------------------------------------------'
    WRITE(unit,"(A42)") title
    WRITE(unit,"(A42)") '------------------------------------------'

    WRITE(unit,"(A)") 'CHIRAL POTENTIAL: NUMERICS'
    CALL print_pwd_numerics(unit)
    WRITE(unit,"(A)") 'CHIRAL POTENTIAL: CONSTANTS'
    CALL chp_print_rconst(unit, chp_mnuc(-1))
    CALL chp_print_rconst(unit, chp_mnuc(0))
    CALL chp_print_rconst(unit, chp_mnuc(+1))
    CALL chp_print_rconst(unit, chp_mpi(-1))
    CALL chp_print_rconst(unit, chp_mpi(0))
    CALL chp_print_rconst(unit, chp_mpi(+1))
    CALL chp_print_rconst(unit, chp_mpi(+2))
    WRITE(unit,"(A42)") '------------------------------------------'
    CALL chp_print_iconst(unit, chp_chiral_order)
    CALL chp_print_rconst(unit, chp_ren)
    CALL chp_print_rconst(unit, chp_lambda)
    CALL chp_print_rconst(unit, chp_regcut_1PE)
    CALL chp_print_rconst(unit, chp_regcut_2PE)
    CALL chp_print_cconst(unit, chp_itope)
    CALL chp_print_rconst(unit, chp_gA)
    CALL chp_print_rconst(unit, chp_fpi)
    CALL chp_print_rconst(unit, chp_fine_structure)
    WRITE(unit,"(A42)") '------------------------------------------'
    IF (chp_chiral_order%val == LO) THEN
       CALL chp_print_rconst(unit, chp_CLO(1))
       CALL chp_print_rconst(unit, chp_CLO(2))
    END IF
    IF (chp_chiral_order%val == NLO) THEN
       DO i=1,2
          DO j=-1,1
             CALL chp_print_rconst(unit, chp_CIB_CLO(j,i))
          END DO
       END DO
       DO i=1,7
          CALL chp_print_rconst(unit, chp_CNLO(i))
       END DO
    END IF
    IF (chp_chiral_order%val >= NNLO) THEN
       CALL chp_print_rconst(unit, chp_c1)
       CALL chp_print_rconst(unit, chp_c3)
       CALL chp_print_rconst(unit, chp_c4)
       DO i=1,2
          DO j=-1,1
             CALL chp_print_rconst(unit, chp_CIB_CLO(j,i))
          END DO
       END DO
       DO i=1,7
          CALL chp_print_rconst(unit, chp_CNLO(i))
       END DO
    END IF
    IF (chp_chiral_order%val >= N3LO) THEN
       DO i=1,15
          CALL chp_print_rconst(unit, chp_DN3LO(i))
       END DO
    END IF
    
  END SUBROUTINE chp_print_constants
  
  SUBROUTINE chp_print_rconst(unit, chp_const)
    
    INTEGER, INTENT(IN) :: unit
    TYPE(chp_real_type) :: chp_const
    
    IF (chp_const%set)  THEN
       WRITE(unit,"(A12,F30.16)") chp_const%name, chp_const%val
    ELSE
       WRITE(unit,"(A12,A30)") chp_const%name, '----------------'
    END IF
    
  END SUBROUTINE chp_print_rconst

  SUBROUTINE chp_print_cconst(unit, chp_const)
    
    INTEGER, INTENT(IN) :: unit
    TYPE(chp_char2_type) :: chp_const
    
    IF (chp_const%set)  THEN
       WRITE(unit,"(A12,A30)") chp_const%name, chp_const%val
    ELSE
       WRITE(unit,"(A12,A30)") chp_const%name, '----------------'
    END IF
    
  
END SUBROUTINE chp_print_cconst

  SUBROUTINE chp_print_iconst(unit, chp_const)
    
    INTEGER, INTENT(IN) :: unit
    TYPE(chp_int_type) :: chp_const
    
    IF (chp_const%set)  THEN
       WRITE(unit,"(A12,I30)") chp_const%name, chp_const%val
    ELSE
       WRITE(unit,"(A12,A30)") chp_const%name, '----------------'
    END IF
    
  END SUBROUTINE chp_print_iconst

  SUBROUTINE chp_set_mass_nucleon(set_mnuc)
    
    REAL(8), INTENT(IN) :: set_mnuc(-1:1)
    
    chp_mnuc(-1)%val  = set_mnuc(-1) ! proton mass
    chp_mnuc(-1)%set  = .TRUE.
    chp_mnuc(0)%val   = set_mnuc(0) ! nucleon mass
    chp_mnuc(0)%set   = .TRUE.
    chp_mnuc(+1)%val  = set_mnuc(+1) ! neutron mass
    chp_mnuc(+1)%set  = .TRUE.
    
  END SUBROUTINE chp_set_mass_nucleon
  
  SUBROUTINE chp_get_mass_nucleon(get_mnuc)
    
    REAL(8), INTENT(INOUT) :: get_mnuc(-1:1)
    
    IF (.NOT. chp_mnuc(-1)%set) THEN
       WRITE(CHP_ERR,"(A)") 'mnuc(-1) NOT SET'
       STOP
    END IF
    
    get_mnuc(-1) = chp_mnuc(-1)%val ! proton mass

    IF (.NOT. chp_mnuc(0)%set) THEN
       WRITE(CHP_ERR,"(A)") 'mnuc(0) NOT SET'
       STOP
    END IF
    
    get_mnuc(0) = chp_mnuc(0)%val ! nucleon mass

    IF (.NOT. chp_mnuc(+1)%set) THEN
       WRITE(CHP_ERR,"(A)") 'mnuc(+1) NOT SET'
       STOP
    END IF
    
    get_mnuc(+1) = chp_mnuc(+1)%val ! neutron mass
        
  END SUBROUTINE chp_get_mass_nucleon

  SUBROUTINE chp_set_mass_pion(set_mpi)
    
    REAL(8), INTENT(IN) :: set_mpi(-1:1)
    
    chp_mpi(-1)%val  = set_mpi(-1) ! pi-
    chp_mpi(-1)%set  = .TRUE.
    chp_mpi(0)%val   = set_mpi(0)  ! pi0
    chp_mpi(0)%set   = .TRUE.
    chp_mpi(+1)%val  = set_mpi(+1) ! pi+
    chp_mpi(+1)%set  = .TRUE.
    
    chp_mpi(+2)%val = (chp_mpi(+1)%val + chp_mpi(-1)%val + chp_mpi(0)%val)/3.0D0
    chp_mpi(+2)%set = .TRUE.
    
  END SUBROUTINE chp_set_mass_pion

  SUBROUTINE chp_get_mass_pion(get_mpi)
    
    REAL(8), INTENT(INOUT) :: get_mpi(-1:2)
    
    IF (.NOT. chp_mpi(-1)%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_mpi(-1) NOT SET'
       STOP
    END IF
    
    get_mpi(-1) = chp_mpi(-1)%val! pi-

    IF (.NOT. chp_mpi(01)%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_mpi(01) NOT SET'
       STOP
    END IF
    
    get_mpi(0) = chp_mpi(0)%val! pi0

    IF (.NOT. chp_mpi(+1)%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_mpi(+1) NOT SET'
       STOP
    END IF
    
    get_mpi(+1) = chp_mpi(+1)%val! pi+

    IF (.NOT. chp_mpi(+2)%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_mpi(+2) NOT SET'
       STOP
    END IF
    
    get_mpi(+2) = chp_mpi(+2)%val! average pion mass

  END SUBROUTINE chp_get_mass_pion
  
  SUBROUTINE chp_set_chiral_order(set_chiral_order)
    
    INTEGER, INTENT(IN) :: set_chiral_order
    INTEGER :: i
    
    IF (set_chiral_order == 1 .OR. set_chiral_order>4 .OR. set_chiral_order <0) THEN
       WRITE(CHP_ERR,"(A,I5)") 'error(set_chiral_order): illegal chiral order', set_chiral_order
       STOP
    END IF
    chp_chiral_order%val = set_chiral_order
    chp_chiral_order%set = .TRUE.
    chp_regcut_1PE%set = .FALSE.
    chp_regcut_2PE%set = .FALSE.
    DO i = 1, 2
       chp_CLO_n(i)%set = .FALSE.
    END DO
    DO i = 1, 7
       chp_CNLO_n(i)%set = .FALSE.
    END DO
    DO i = 1, 15
       chp_DN3LO_n(i)%set = .FALSE.
    END DO

  END SUBROUTINE chp_set_chiral_order

  SUBROUTINE chp_get_chiral_order(get_chiral_order)
    
    INTEGER, INTENT(INOUT) :: get_chiral_order
    
    IF (.NOT. chp_chiral_order%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_chiral_order NOT SET'
       STOP
    END IF
    
    get_chiral_order = chp_chiral_order%val
    
  END SUBROUTINE chp_get_chiral_order
  
  SUBROUTINE chp_set_itope(set_itope)
    
    CHARACTER(LEN=2), INTENT(IN) :: set_itope
    
    IF (set_itope /= 'KW' .AND. set_itope/= 'EM') THEN
       WRITE(CHP_ERR, "(A,A)") 'error(chp_set_reg): illegal iterated ope treatment', set_itope
       STOP
    END IF
    
    chp_itope%val = set_itope
    chp_itope%set  = .TRUE.
    
  END SUBROUTINE chp_set_itope

  SUBROUTINE chp_get_itope(get_itope)
    
    CHARACTER(LEN=2), INTENT(INOUT) :: get_itope
    
    IF (.NOT. chp_itope%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_itope NOT SET'
       STOP
    END IF
    
    get_itope = chp_itope%val
    
  END SUBROUTINE chp_get_itope
  
  SUBROUTINE chp_set_1PE_reg_par(val)
    REAL(8), INTENT(IN) :: val

    IF (.NOT. chp_chiral_order%set) THEN
       STOP 'Must set chiral order before setting regulator parameters'
    END IF
    chp_regcut_1PE%set = .TRUE.
    chp_regcut_1PE%val = val
  END SUBROUTINE chp_set_1PE_reg_par

  SUBROUTINE chp_set_2PE_reg_par(val)
    REAL(8), INTENT(IN) :: val

    IF (.NOT. chp_chiral_order%set) THEN
       STOP 'Must set chiral order before setting regulator parameters'
    END IF
    chp_regcut_2PE%set = .TRUE.
    chp_regcut_2PE%val = val
  END SUBROUTINE chp_set_2PE_reg_par

  SUBROUTINE chp_set_reg(set_reg, set_cutoff)
    
    CHARACTER(LEN=2), INTENT(IN) :: set_reg
    REAL(8), INTENT(IN) :: set_cutoff
    
    IF (set_reg /= 'SF' .AND. set_reg/= 'DR') THEN
       WRITE(CHP_ERR, "(A,A)") 'error(chp_set_reg): illegal regularization choice', set_reg
       STOP
    END IF
    
    chp_ren%val  = set_cutoff
    chp_ren%name = set_reg
    chp_ren%set  = .TRUE.
    
  END SUBROUTINE chp_set_reg
  
  SUBROUTINE chp_get_reg(get_reg_type,get_cutoff)
    
    CHARACTER(LEN=2), INTENT(INOUT) :: get_reg_type
    REAL(8), INTENT(INOUT) :: get_cutoff
    
    IF (.NOT. chp_ren%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_ren NOT SET'
       STOP
    END IF
    
    get_cutoff = chp_ren%val
    get_reg_type = chp_ren%name(1:2)
    
  END SUBROUTINE chp_get_reg
  
  SUBROUTINE chp_set_contact_format(set_cnt_fmt)
    
    CHARACTER(LEN=2), INTENT(IN) :: set_cnt_fmt
    
    IF (set_cnt_fmt /= 'PW' .AND. set_cnt_fmt/= 'ST') THEN
       WRITE(CHP_ERR, "(A,A)") 'error(chp_contact_format): illegal contact format', set_cnt_fmt
       STOP
    END IF
    
    chp_contact_format%name = set_cnt_fmt
    chp_contact_format%set  = .TRUE.
    
  END SUBROUTINE chp_set_contact_format

  SUBROUTINE chp_get_contact_format(get_cnt_fmt)
    
    CHARACTER(LEN=2), INTENT(INOUT) :: get_cnt_fmt
    
    IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_contact_format NOT SET'
       STOP
    END IF
    
    get_cnt_fmt = chp_contact_format%name(1:2)
    
  END SUBROUTINE chp_get_contact_format
  
  SUBROUTINE chp_set_LO_contact_reg_par(contact_no, val)
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val

    IF (.NOT. chp_chiral_order%set) THEN
       STOP 'Must set chiral order before setting regulator parameters'
    END IF
    chp_CLO_n(contact_no)%set = .TRUE.
    chp_CLO_n(contact_no)%val = val
  END SUBROUTINE chp_set_LO_contact_reg_par

  SUBROUTINE chp_set_LO_contact(contact_no, val)    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val
    
    IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'error(chp_set_LO_contact): you must set the format before the contact value'
       STOP
    END IF
    ! LO
    ! the LO contacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2
    !
    chp_CLO(contact_no)%val  = val* 0.01D0
    chp_CLO(contact_no)%set  = .TRUE.
    
    SELECT CASE(contact_no)
    CASE(1)
       IF (chp_contact_format%name == 'PW') chp_CLO(1)%name = '1S0'
       IF (chp_contact_format%name == 'ST') chp_CLO(1)%name = 'CS->1S0'
    CASE(2)
       IF (chp_contact_format%name == 'PW') chp_CLO(2)%name = '3S1'
       IF (chp_contact_format%name == 'ST') chp_CLO(2)%name = 'CT->3S1'
    CASE DEFAULT
       WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
       STOP
    END SELECT

  END SUBROUTINE chp_set_LO_contact

  SUBROUTINE chp_get_LO_contact(contact_no, val, name)    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(INOUT) :: val
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    

    ! LO
    ! the LO contacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2
    !
    
    IF (.NOT. chp_CLO(contact_no)%set) THEN
       WRITE(CHP_ERR,"(A,I5,A)") 'chp_CLO',contact_no, ' NOT SET'
       STOP
    END IF
    
    val = chp_CLO(contact_no)%val * 100.0D0
    name = chp_CLO(contact_no)%name

  END SUBROUTINE chp_get_LO_contact

  SUBROUTINE chp_set_CIB_LO_contact(contact_no, tz, val)    
    INTEGER, INTENT(IN) :: contact_no, tz
    REAL(8), INTENT(IN) :: val
    
    IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'error(chp_set_LO_contact): you must set the format before the contact value'
       STOP
    END IF

    chp_CIB_CLO(tz,contact_no)%val  = val*0.01D0
    chp_CIB_CLO(tz,contact_no)%set  = .TRUE.
    
    IF (tz == -1) THEN
       SELECT CASE(contact_no)
       CASE(1)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,1)%name = '1S0pp'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,1)%name = 'CS->1S0pp'
       CASE(2)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,2)%name = '3S1pp'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,2)%name = 'CT->3S1pp'
       CASE DEFAULT
          WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
          STOP
       END SELECT
    END IF

    IF (tz == 0) THEN
       SELECT CASE(contact_no)
       CASE(1)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,1)%name = '1S0pn'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,1)%name = 'CS->1S0pn'
       CASE(2)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,2)%name = '3S1pn'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,2)%name = 'CT->3S1pn'
       CASE DEFAULT
          WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
          STOP
       END SELECT
    END IF

    IF (tz == +1) THEN
       SELECT CASE(contact_no)
       CASE(1)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,1)%name = '1S0nn'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,1)%name = 'CS->1S0nn'
       CASE(2)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,2)%name = '3S1nn'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,2)%name = 'CT->3S1nn'
       CASE DEFAULT
          WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
          STOP
       END SELECT
    END IF
  END SUBROUTINE chp_set_CIB_LO_contact

  SUBROUTINE chp_get_CIB_LO_contact(contact_no, tz, val,name)    
    INTEGER, INTENT(IN) :: contact_no, tz
    REAL(8), INTENT(INOUT) :: val
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_CIB_CLO(tz,contact_no)%set) THEN
       WRITE(CHP_ERR,"(A,2I5,A)") 'chp_CIB_CLO ', tz, contact_no, ' NOT SET'
       STOP
    END IF

    ! the LO CIBcontacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2
    
    val = chp_CIB_CLO(tz,contact_no)%val*100.0D0
    name = chp_CIB_CLO(tz,contact_no)%name
    
  END SUBROUTINE chp_get_CIB_LO_contact
  
  SUBROUTINE chp_set_NLO_contact_reg_par(contact_no, val)
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val

    IF (.NOT. chp_chiral_order%set) THEN
       STOP 'Must set chiral order before setting regulator parameters'
    END IF
    chp_CNLO_n(contact_no)%set = .TRUE.
    chp_CNLO_n(contact_no)%val = val
  END SUBROUTINE chp_set_NLO_contact_reg_par

  SUBROUTINE chp_set_NLO_contact(contact_no, val)    
    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val
    
    IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'error(chp_set_NLO_contact): you must set the format before the contact value'
       STOP
    END IF
    
    ! the NLO contacts are input in units of 10^4/GeV^4
    ! this = 10^-8/MeV^4
    chp_CNLO(contact_no)%val  = val*1.E-08
    chp_CNLO(contact_no)%set  = .TRUE.
    
    SELECT CASE(contact_no)
    CASE(1)
       IF (chp_contact_format%name == 'PW') chp_CNLO(1)%name = '1S0'
       IF (chp_contact_format%name == 'ST') chp_CNLO(1)%name = 'C1->1S0'
    CASE(2)
       IF (chp_contact_format%name == 'PW') chp_CNLO(2)%name = '3P0'
       IF (chp_contact_format%name == 'ST') chp_CNLO(2)%name = 'C2->3P0'
    CASE(3)
       IF (chp_contact_format%name == 'PW') chp_CNLO(3)%name = '1P1'
       IF (chp_contact_format%name == 'ST') chp_CNLO(3)%name = 'C3->1P1'
    CASE(4)
       IF (chp_contact_format%name == 'PW') chp_CNLO(4)%name = '3P1'
       IF (chp_contact_format%name == 'ST') chp_CNLO(4)%name = 'C4->3P1'
    CASE(5)
       IF (chp_contact_format%name == 'PW') chp_CNLO(5)%name = '3S1'
       IF (chp_contact_format%name == 'ST') chp_CNLO(5)%name = 'C5->3S1'
    CASE(6)
       IF (chp_contact_format%name == 'PW') chp_CNLO(6)%name = '3S1-3D1'
       IF (chp_contact_format%name == 'ST') chp_CNLO(6)%name = 'C6->3S1-3D1'
    CASE(7)
       IF (chp_contact_format%name == 'PW') chp_CNLO(7)%name = '3P2'
       IF (chp_contact_format%name == 'ST') chp_CNLO(7)%name = 'C7->3P2'
    CASE DEFAULT
       WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
       STOP
    END SELECT
    
  END SUBROUTINE chp_set_NLO_contact

  SUBROUTINE chp_get_NLO_contact(contact_no, val, name)    
    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(INOUT) :: val
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    
    IF (.NOT. chp_CNLO(contact_no)%set) THEN
       WRITE(CHP_ERR,"(A,I5,A)") 'chp_CNLO ', contact_no, ' NOT SET'
       STOP
    END IF
    ! the NLO contacts are input in units of 10^4/GeV^4
    ! this = 10^-8/MeV^4
    val = chp_CNLO(contact_no)%val*1.0E+08
    name = chp_CNLO(contact_no)%name
    
  END SUBROUTINE chp_get_NLO_contact

  SUBROUTINE chp_set_N3LO_contact_reg_par(contact_no, val)
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val

    IF (.NOT. chp_chiral_order%set) THEN
       STOP 'Must set chiral order before setting regulator parameters'
    END IF
    chp_DN3LO_n(contact_no)%set = .TRUE.
    chp_DN3LO_n(contact_no)%val = val
  END SUBROUTINE chp_set_N3LO_contact_reg_par

  SUBROUTINE chp_set_N3LO_contact(contact_no, val)    
    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val
    
    IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'error(chp_set_N3LO_contact): you must set the format before the contact value'
       STOP
    END IF

    ! the N3LO contacts are input in units of 10^4/GeV^6
    ! this = 10^(4-3*6)/MeV^6 = 10^-14/MeV^6
    chp_DN3LO(contact_no)%val  = val * 1.D-14
    chp_DN3LO(contact_no)%set  = .TRUE.
    
    SELECT CASE(contact_no)
    CASE(1)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(1)%name = 't1S0'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(1)%name = 'D1->t1S0'
    CASE(2)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(2)%name = '1S0'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(2)%name = 'D2->1S0'
    CASE(3)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(3)%name = '3P0'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(3)%name = 'D3->3P0'
    CASE(4)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(4)%name = '1P1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(4)%name = 'D4->1P1'
    CASE(5)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(5)%name = '3P1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(5)%name = 'D5->3P1'
    CASE(6)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(6)%name = 't3S1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(6)%name = 'D6->t3S1'
    CASE(7)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(7)%name = '3S1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(7)%name = 'D7->3S1'
    CASE(8)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(8)%name = '3D1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(8)%name = 'D8->3D1'
    CASE(9)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(9)%name = 't3S1-3D1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(9)%name = 'D9->t3S1-3D1'
    CASE(10)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(10)%name = '3S1-3D1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(10)%name = 'D10->3S1-3D1'
    CASE(11)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(11)%name = '1D2'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(11)%name = 'D11->1D2'
    CASE(12)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(12)%name = '3D2'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(12)%name = 'D12->3D2'
    CASE(13)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(13)%name = '3P2'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(13)%name = 'D13->3P2'
    CASE(14)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(14)%name = '3P2-3F2'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(14)%name = 'D14->3P2-3F2'
    CASE(15)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(15)%name = '3D3'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(15)%name = 'D15->3D3'
    CASE DEFAULT
       WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_N3LO_contact): illegal contact nummber', contact_no
       STOP
    END SELECT
    
  END SUBROUTINE chp_set_N3LO_contact

  SUBROUTINE chp_get_N3LO_contact(contact_no, val, name)    
    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(INOUT) :: val
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    
    IF (.NOT. chp_DN3LO(contact_no)%set) THEN
       WRITE(CHP_ERR,"(A,I5,A)") 'chp_DN3LO ', contact_no, ' NOT SET'
       STOP
    END IF
    ! the N3LO contacts are input in units of 10^4/GeV^6
    ! this = 10^(4-3*6)/MeV^6 = 10^-14/MeV^6
    val = chp_DN3LO(contact_no)%val*1.0D+14
    name = chp_DN3LO(contact_no)%name
    
  END SUBROUTINE chp_get_N3LO_contact

  SUBROUTINE chp_set_c1(set_c1)
    
    REAL(8), INTENT(IN) :: set_c1

    chp_c1%val  = set_c1 * 1.0E-03 ! transform to MeV^-1
    chp_c1%set  = .TRUE.
    
  END SUBROUTINE chp_set_c1

  SUBROUTINE chp_get_c1(get_c1, name)
    
    REAL(8), INTENT(INOUT) :: get_c1
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_c1%set) THEN
       WRITE(CHP_ERR,"(A)") 'c1 NOT SET'
       STOP
    END IF
    
    get_c1 = chp_c1%val * 1000.0D0 ! transform to  GeV^-1
    name = chp_c1%name

  END SUBROUTINE chp_get_c1

  SUBROUTINE chp_set_c3(set_c3)
    
    REAL(8), INTENT(IN) :: set_c3
    
    chp_c3%val  = set_c3*1.0E-3
    chp_c3%set  = .TRUE.
    
  END SUBROUTINE chp_set_c3

  SUBROUTINE chp_get_c3(get_c3,name)
    
    REAL(8), INTENT(INOUT) :: get_c3
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_c3%set) THEN
       WRITE(CHP_ERR,"(A)") 'c3 NOT SET'
       STOP
    END IF
    get_c3 = chp_c3%val * 1000.0D0 ! transform to GeV^-1
    name = chp_c3%name
    
  END SUBROUTINE chp_get_c3
  
  SUBROUTINE chp_set_c4(set_c4)
    
    REAL(8), INTENT(IN) :: set_c4
    
    chp_c4%val  = set_c4*1.0E-3
    chp_c4%set  = .TRUE.
    
  END SUBROUTINE chp_set_c4

  SUBROUTINE chp_get_c4(get_c4,name)
    
    REAL(8), INTENT(INOUT) :: get_c4
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_c4%set) THEN
       WRITE(CHP_ERR,"(A)") 'c4 NOT SET'
       STOP
    END IF
    
    get_c4 = chp_c4%val * 1000.0D0 ! transform to GeV^-1
    name = chp_c4%name

  END SUBROUTINE chp_get_c4
  
  SUBROUTINE chp_set_c2(set_c2)
    
    REAL(8), INTENT(IN) :: set_c2
    
    chp_c2%val  = set_c2 * 1.0D-3
    chp_c2%set  = .TRUE.
    
  END SUBROUTINE chp_set_c2
  
  SUBROUTINE chp_get_c2(get_c2,name)
    
    REAL(8), INTENT(INOUT) :: get_c2
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_c2%set) THEN
       WRITE(CHP_ERR,"(A)") 'c2 NOT SET'
       STOP
    END IF
    
    get_c2 = chp_c2%val * 1000.0D0 ! transform to GeV^-1
    name = chp_c2%name

  END SUBROUTINE chp_get_c2

  SUBROUTINE chp_set_d1_plus_d2(set_d1_plus_d2)
    
    REAL(8), INTENT(IN) :: set_d1_plus_d2
    
    chp_d1_plus_d2%val  = set_d1_plus_d2 * 1.0D-6
    chp_d1_plus_d2%set  = .TRUE.
    
  END SUBROUTINE chp_set_d1_plus_d2
  
  SUBROUTINE chp_get_d1_plus_d2(get_d1_plus_d2,name)
    
    REAL(8), INTENT(INOUT) :: get_d1_plus_d2
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_d1_plus_d2%set) THEN
       WRITE(CHP_ERR,"(A)") 'd1+d2 NOT SET'
       STOP
    END IF
    
    get_d1_plus_d2 = chp_d1_plus_d2%val * 1.0D6 ! transform to GeV^-2
    name = chp_d1_plus_d2%name

  END SUBROUTINE chp_get_d1_plus_d2

  SUBROUTINE chp_set_d3(set_d3)
    
    REAL(8), INTENT(IN) :: set_d3
    
    chp_d3%val  = set_d3 * 1.0D-6
    chp_d3%set  = .TRUE.
    
  END SUBROUTINE chp_set_d3
  
  SUBROUTINE chp_get_d3(get_d3,name)
    
    REAL(8), INTENT(INOUT) :: get_d3
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_d3%set) THEN
       WRITE(CHP_ERR,"(A)") 'd3 NOT SET'
       STOP
    END IF
    
    get_d3 = chp_d3%val * 1.0D6 ! transform to GeV^-2
    name = chp_d3%name

  END SUBROUTINE chp_get_d3

  SUBROUTINE chp_set_d5(set_d5)
    
    REAL(8), INTENT(IN) :: set_d5
    
    chp_d5%val  = set_d5 * 1.0D-6
    chp_d5%set  = .TRUE.
    
  END SUBROUTINE chp_set_d5
  
  SUBROUTINE chp_get_d5(get_d5,name)
    
    REAL(8), INTENT(INOUT) :: get_d5
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_d5%set) THEN
       WRITE(CHP_ERR,"(A)") 'd5 NOT SET'
       STOP
    END IF
    
    get_d5 = chp_d5%val * 1.0D6 ! transform to GeV^-2
    name = chp_d5%name

  END SUBROUTINE chp_get_d5

  SUBROUTINE chp_set_d14_minus_d15(set_d14_minus_d15)
    
    REAL(8), INTENT(IN) :: set_d14_minus_d15
    
    chp_d14_minus_d15%val  = set_d14_minus_d15 * 1.0D-6
    chp_d14_minus_d15%set  = .TRUE.
    
  END SUBROUTINE chp_set_d14_minus_d15

  SUBROUTINE chp_get_d14_minus_d15(get_d14_minus_d15,name)
    
    REAL(8), INTENT(INOUT) :: get_d14_minus_d15
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_d14_minus_d15%set) THEN
       WRITE(CHP_ERR,"(A)") 'd14_minus_d15 NOT SET'
       STOP
    END IF
    
    get_d14_minus_d15 = chp_d14_minus_d15%val * 1.0D6 ! transform to GeV^-2
    name = chp_d14_minus_d15%name

  END SUBROUTINE chp_get_d14_minus_d15

  SUBROUTINE chp_set_lambda(set_lambda)
    
    REAL(8), INTENT(IN) :: set_lambda

    chp_lambda%val  = set_lambda
    chp_lambda%set  = .TRUE.

  END SUBROUTINE chp_set_lambda

  SUBROUTINE chp_get_lambda(get_lambda)
    
    REAL(8), INTENT(INOUT) :: get_lambda

    IF (.NOT. chp_lambda%set) THEN
       WRITE(CHP_ERR,"(A)") 'lambda NOT SET'
       STOP
    END IF
    
    get_lambda = chp_lambda%val
  
  END SUBROUTINE chp_get_lambda

  SUBROUTINE chp_set_gA(set_gA)
    
    REAL(8), INTENT(IN) :: set_gA

    ! This depends on gA, so make it unset when gA is changed
    IF(.NOT. chp_gA%set .OR. chp_gA%val /= set_gA) THEN
       chp_2PE_2loop_int_WC_data_set  = .FALSE.
    END IF
    chp_gA%val  = set_gA
    chp_gA%set  = .TRUE.

  END SUBROUTINE chp_set_gA

  SUBROUTINE chp_get_gA(get_gA)
    
    REAL(8), INTENT(INOUT) :: get_gA
    
    IF (.NOT. chp_gA%set) THEN
       WRITE(CHP_ERR,"(A)") 'gA NOT SET'
       STOP
    END IF
    
    get_gA = chp_gA%val
  
  END SUBROUTINE chp_get_gA

  SUBROUTINE chp_set_fpi(set_fpi)
    
    REAL(8), INTENT(IN) :: set_fpi

    chp_fpi%val  = set_fpi
    chp_fpi%set  = .TRUE.
    
  END SUBROUTINE chp_set_fpi

  SUBROUTINE chp_get_fpi(get_fpi)
    
    REAL(8), INTENT(INOUT) :: get_fpi

    IF (.NOT. chp_fpi%set) THEN
       WRITE(CHP_ERR,"(A)") 'fpi NOT SET'
       STOP
    END IF
    
    get_fpi = chp_fpi%val
    
  END SUBROUTINE chp_get_fpi

  SUBROUTINE chp_set_fine_structure(set_fine_structure)
    
    REAL(8), INTENT(IN) :: set_fine_structure

    chp_fine_structure%val  = set_fine_structure
    chp_fine_structure%set  = .TRUE.
    
  END SUBROUTINE chp_set_fine_structure

  SUBROUTINE chp_get_fine_structure(get_fine_structure)
    
    REAL(8), INTENT(INOUT) :: get_fine_structure

    IF (.NOT. chp_fine_structure%set) THEN
       WRITE(CHP_ERR,"(A)") 'fine_structure NOT SET'
       STOP
    END IF
    
    get_fine_structure = chp_fine_structure%val
    
  END SUBROUTINE chp_get_fine_structure

  SUBROUTINE chp_set_2PE_CSB_correct_mass(set_correct)
    
    INTEGER, INTENT(IN) :: set_correct

    chp_2PE_CSB_correct_mass%val  = set_correct
    chp_2PE_CSB_correct_mass%set  = .TRUE.
    
  END SUBROUTINE chp_set_2PE_CSB_correct_mass

  SUBROUTINE chp_get_2PE_CSB_correct_mass(get_correct)
    
    INTEGER, INTENT(INOUT) :: get_correct

    IF (.NOT. chp_2PE_CSB_correct_mass%set) THEN
       WRITE(CHP_ERR,"(A)") '2PE_CSB_correct_mass NOT SET'
       STOP
    END IF
    
    get_correct = chp_2PE_CSB_correct_mass%val
    
  END SUBROUTINE chp_get_2PE_CSB_correct_mass

  SUBROUTINE chp_set_2PE_2loop_int(set)
    
    INTEGER, INTENT(IN) :: set

    chp_use_2PE_2loop_int%val  = set
    chp_use_2PE_2loop_int%set  = .TRUE.
    
  END SUBROUTINE chp_set_2PE_2loop_int

  SUBROUTINE chp_get_2PE_2loop_int(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chp_use_2PE_2loop_int%set) THEN
       WRITE(CHP_ERR,"(A)") '2PE_2loop_int NOT SET'
       STOP
    END IF
    
    get = chp_use_2PE_2loop_int%val
    
  END SUBROUTINE chp_get_2PE_2loop_int

END MODULE idaho_chiral_potential

