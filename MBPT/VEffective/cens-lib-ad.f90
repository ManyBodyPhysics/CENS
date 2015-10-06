!!$****************************************************************
!!$ TaylUR:
!!$ An arbitrary-order auto-differentiation module
!!$ (revised version)
!!$
!!$ (C) copyright Georg M. von Hippel, 2005-2007
!!$
!!$****************************************************************

MODULE cens_lib_ad
  USE cens_kinds
  
!!$ Define useintrinsic if your compiler supports the IEEE modules
  !!!DEC$ DEFINE useintrinsic
  
  !DEC$ IF DEFINED (useintrinsic)
  USE, INTRINSIC :: ieee_arithmetic
  !DEC$ ENDIF
  
  IMPLICIT NONE
  
  PRIVATE
  
!!$****************************************************************
!!$ CONSTANTS AND VARIABLES
!!$****************************************************************
  
!!$ Kinds of floating point numbers to use
  
  !INTEGER, PARAMETER :: dr_kind = KIND(0.D0)
!  INTEGER, PARAMETER :: dc_kind = KIND((0.D0,0.D0))
  INTEGER, PARAMETER :: dr_kind = dp
  INTEGER, PARAMETER :: dc_kind = dp!KIND((rs_kind, rs_kind)) 
  INTEGER, PARAMETER :: sr_kind = KIND(0.)
  INTEGER, PARAMETER :: sc_kind = KIND((0.,0.))
  
!!$ Number of variables
  
  INTEGER, PARAMETER, PUBLIC :: N_taylor_vars = 1
  
!!$ Definition of Taylor expansion order
  
  INTEGER, PARAMETER, PUBLIC :: Max_taylor_order = 15
  
  
!!$ Order at which to cut off the Taylor expansion of quantities
  
  INTEGER, PUBLIC :: Taylor_order = Max_taylor_order
  
  
!!$  taylor object with an imaginary part greater than
!!$  Real_args_tol is passed to the real-argument only F95
!!$  functions, and Real_args_warn is set to .TRUE., a NaN
!!$  value is returned as the result.
  
  LOGICAL, PUBLIC :: Real_args_warn = .TRUE.
  REAL(kind=dr_kind), PUBLIC :: Real_args_tol = 1.D-15
  
!!$****************************************************************
!!$***************** DO NOT EDIT BELOW THIS POINT *****************
!!$****************************************************************

!!$****************************************************************
!!$ TYPES AND INTERFACES
!!$****************************************************************
  
!!$ Definition of taylor type
  
  PUBLIC taylor
  
  TYPE taylor
     COMPLEX(kind=dc_kind) :: drv(1:N_taylor_vars,0:Max_taylor_order)
  END TYPE taylor
  
!!$ Interface for taylor assignment
  
  PUBLIC ASSIGNMENT(=)
  
  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE assign_taylors, &
          assign_taylor_real, assign_taylor_sreal, &
          assign_taylor_complex, assign_taylor_scomplex, &
          assign_taylor_int, &
          assign_real_taylor, assign_sreal_taylor, &
          assign_complex_taylor, assign_scomplex_taylor
  END INTERFACE
  
!!$ Interfaces for arithmetic operations on taylors
  
   PUBLIC OPERATOR(+), OPERATOR(-), OPERATOR(*), OPERATOR(/), &
        OPERATOR(**)
   
   INTERFACE OPERATOR(+)
      MODULE PROCEDURE add_taylors, &
           add_taylor_real, add_real_taylor, &
           add_taylor_sreal, add_sreal_taylor, &
           add_taylor_complex, add_complex_taylor, &
           add_taylor_scomplex, add_scomplex_taylor, &
           add_int_taylor,add_taylor_int
   END INTERFACE
   
   INTERFACE OPERATOR(-)
      MODULE PROCEDURE neg_taylor, sub_taylors, &
           sub_real_taylor, sub_taylor_real, &
           sub_sreal_taylor, sub_taylor_sreal, &
           sub_taylor_complex, sub_complex_taylor, &
           sub_taylor_scomplex, sub_scomplex_taylor, &
           sub_int_taylor, sub_taylor_int
   END INTERFACE
   
   INTERFACE OPERATOR(*)
      MODULE PROCEDURE mult_taylors, &
           mult_taylor_real,mult_real_taylor, &
           mult_taylor_sreal,mult_sreal_taylor, &
           mult_taylor_complex,mult_complex_taylor, &
           mult_taylor_scomplex,mult_scomplex_taylor, &
           mult_taylor_int,mult_int_taylor
   END INTERFACE
   
   INTERFACE OPERATOR(/)
      MODULE PROCEDURE div_taylors, &
           div_taylor_real, div_real_taylor, &
           div_taylor_sreal, div_sreal_taylor, &
           div_taylor_complex, div_complex_taylor, &
           div_taylor_scomplex, div_scomplex_taylor, &
           div_taylor_int, div_int_taylor
   END INTERFACE
   
   INTERFACE OPERATOR(**)
      MODULE PROCEDURE power_taylor_int
   END INTERFACE

!!$ Interfaces for comparison operators on (the values of) taylors
   
   PUBLIC OPERATOR(==), OPERATOR(.IDENT.), &
        OPERATOR(/=), OPERATOR(.NIDENT.), &
        OPERATOR(<), OPERATOR(<=), &
        OPERATOR(>), OPERATOR(>=)
   
   INTERFACE OPERATOR(==)
      MODULE PROCEDURE eq_taylors, &
           eq_taylor_real, eq_real_taylor, &
           eq_taylor_sreal, eq_sreal_taylor, &
           eq_taylor_complex, eq_complex_taylor, &
           eq_taylor_scomplex, eq_scomplex_taylor, &
           eq_taylor_int, eq_int_taylor
   END INTERFACE
   
   INTERFACE OPERATOR(.IDENT.)
      MODULE PROCEDURE ident_taylors
   END INTERFACE
   
   INTERFACE OPERATOR(/=)
      MODULE PROCEDURE ne_taylors, &
           ne_taylor_real, ne_real_taylor, &
           ne_taylor_sreal, ne_sreal_taylor, &
           ne_taylor_complex, ne_complex_taylor, &
           ne_taylor_scomplex, ne_scomplex_taylor, &
           ne_taylor_int, ne_int_taylor
   END INTERFACE

   INTERFACE OPERATOR(.NIDENT.)
      MODULE PROCEDURE nident_taylors
   END INTERFACE
   
   INTERFACE OPERATOR(<) ! Note: compares real parts only
      MODULE PROCEDURE lt_taylors, &
           lt_taylor_real, lt_real_taylor, &
           lt_taylor_sreal, lt_sreal_taylor, &
           lt_taylor_int, lt_int_taylor
   END INTERFACE
   
   INTERFACE OPERATOR(<=) ! Note: compares real parts only
      MODULE PROCEDURE le_taylors, &
           le_taylor_real, le_real_taylor, &
           le_taylor_sreal, le_sreal_taylor, &
           le_taylor_int, le_int_taylor
   END INTERFACE

   INTERFACE OPERATOR(>) ! Note: compares real parts only
      MODULE PROCEDURE gt_taylors, &
           gt_taylor_real, gt_real_taylor, &
           gt_taylor_sreal, gt_sreal_taylor, &
           gt_taylor_int, gt_int_taylor
   END INTERFACE
   
   INTERFACE OPERATOR(>=) ! Note: compares real parts only
      MODULE PROCEDURE ge_taylors, &
           ge_taylor_real, ge_real_taylor, &
           ge_taylor_sreal, ge_sreal_taylor, &
           ge_taylor_int, ge_int_taylor
   END INTERFACE
   
!!$ Interfaces for mathematical functions on taylors
   
   PUBLIC SQRT, AIMAG, REAL, EXP, LOG, LOG10, COS, SIN, TAN, &
        COSH, SINH, TANH, ACOS, ASIN, ATAN, ATAN2
   
   INTERFACE SQRT
      MODULE PROCEDURE sqrt_taylor
   END INTERFACE
   
   INTERFACE AIMAG ! Note: Does NOT convert to real type
      MODULE PROCEDURE aimag_taylor
   END INTERFACE
   
   INTERFACE REAL ! Note: Does NOT convert to real type
      MODULE PROCEDURE real_taylor
   END INTERFACE
   
   INTERFACE EXP
      MODULE PROCEDURE exp_taylor
   END INTERFACE

   INTERFACE LOG
      MODULE PROCEDURE log_taylor
   END INTERFACE
   
   INTERFACE LOG10 ! Note: takes real part
      MODULE PROCEDURE log10_taylor
   END INTERFACE
   
   INTERFACE COS
      MODULE PROCEDURE cos_taylor
   END INTERFACE
   
   INTERFACE SIN
      MODULE PROCEDURE sin_taylor
   END INTERFACE
   
   INTERFACE TAN
      MODULE PROCEDURE tan_taylor
   END INTERFACE

   INTERFACE COSH
      MODULE PROCEDURE cosh_taylor
   END INTERFACE

   INTERFACE SINH
      MODULE PROCEDURE sinh_taylor
   END INTERFACE

   INTERFACE TANH
      MODULE PROCEDURE tanh_taylor
   END INTERFACE

   INTERFACE ACOS ! Note: takes real part
      MODULE PROCEDURE acos_taylor
   END INTERFACE
      
   INTERFACE ASIN ! Note: takes real part
      MODULE PROCEDURE asin_taylor
   END INTERFACE
   
   INTERFACE ATAN ! Note: takes real part
      MODULE PROCEDURE atan_taylor
   END INTERFACE

   INTERFACE ATAN2 ! Note: takes real part
      MODULE PROCEDURE atan2_taylors, &
          atan2_taylor_real,atan2_real_taylor, &
          atan2_taylor_sreal,atan2_sreal_taylor
   END INTERFACE

!!$  Interfaces for less differentiable functions on taylors

   PUBLIC ABS, SIGN, MAX, MIN, DIM, CONJG, MOD, MODULO, &
          AINT, ANINT, CEILING, FLOOR, INT, NINT

   INTERFACE ABS
      MODULE PROCEDURE abs_taylor
   END INTERFACE

   INTERFACE SIGN ! Note: takes real part
      MODULE PROCEDURE sign_taylors, &
           sign_taylor_real, sign_real_taylor, &
           sign_taylor_sreal, sign_sreal_taylor
   END INTERFACE

   INTERFACE MAX ! Note: takes real part
      MODULE PROCEDURE max_taylors, &
           max_taylor_real, max_real_taylor, &
           max_taylor_sreal, max_sreal_taylor
   END INTERFACE

   INTERFACE MIN ! Note: takes real part
      MODULE PROCEDURE min_taylors, &
           min_taylor_real, min_real_taylor, &
           min_taylor_sreal, min_sreal_taylor
   END INTERFACE

   INTERFACE DIM ! Note: takes real part
      MODULE PROCEDURE dim_taylors, &
           dim_taylor_real, dim_real_taylor, &
           dim_taylor_sreal, dim_sreal_taylor
   END INTERFACE

   INTERFACE CONJG ! Note: assumes independent variables are real
      MODULE PROCEDURE conjg_taylor
   END INTERFACE

   INTERFACE MOD ! Note: takes real part
      MODULE PROCEDURE mod_taylor_real, mod_taylor_sreal
   END INTERFACE

   INTERFACE MODULO ! Note: takes real part
      MODULE PROCEDURE modulo_taylor_real, modulo_taylor_sreal
   END INTERFACE
   
   INTERFACE AINT ! Note: takes real part
      MODULE PROCEDURE aint_taylor
   END INTERFACE

   INTERFACE ANINT ! Note: takes real part
      MODULE PROCEDURE anint_taylor
   END INTERFACE

   INTERFACE CEILING ! Note: takes real part
      MODULE PROCEDURE ceiling_taylor
   END INTERFACE

   INTERFACE FLOOR ! Note: takes real part
      MODULE PROCEDURE floor_taylor
   END INTERFACE

   INTERFACE INT ! Note: takes real part
      MODULE PROCEDURE int_taylor
   END INTERFACE

   INTERFACE NINT ! Note: takes real part
      MODULE PROCEDURE nint_taylor
   END INTERFACE

!!$  Interfaces for array reduction and vectorial functions on taylors

   PUBLIC MAXVAL, MINVAL, MAXLOC, MINLOC, SUM, PRODUCT, DOT_PRODUCT, &
          MATMUL

   INTERFACE MAXVAL ! Note: takes real part
      MODULE PROCEDURE maxval_taylors
   END INTERFACE

   INTERFACE MINVAL ! Note: takes real part
      MODULE PROCEDURE minval_taylors
   END INTERFACE

   INTERFACE MAXLOC ! Note: takes real part
      MODULE PROCEDURE maxloc_taylors
   END INTERFACE

   INTERFACE MINLOC ! Note: takes real part
      MODULE PROCEDURE minloc_taylors
   END INTERFACE

   INTERFACE SUM
      MODULE PROCEDURE sum_taylors
   END INTERFACE
    
   INTERFACE PRODUCT
      MODULE PROCEDURE product_taylors
   END INTERFACE

   INTERFACE DOT_PRODUCT
      MODULE PROCEDURE dot_taylors, &
              dot_taylors_reals, dot_reals_taylors, &
              dot_taylors_sreals, dot_sreals_taylors, &
              dot_taylors_complexes, dot_complexes_taylors, &
              dot_taylors_scomplexes, dot_scomplexes_taylors
   END INTERFACE

   INTERFACE MATMUL
      MODULE PROCEDURE matmul_taylors, matvec_taylors, vecmat_taylors, &
              matmul_taylors_reals, matmul_reals_taylors, &
              matvec_taylors_reals, matvec_reals_taylors, &
              vecmat_taylors_reals, vecmat_reals_taylors, &
              matmul_taylors_sreals, matmul_sreals_taylors, &
              matvec_taylors_sreals, matvec_sreals_taylors, &
              vecmat_taylors_sreals, vecmat_sreals_taylors, &
              matmul_taylors_complexes, matmul_complexes_taylors, &
              matvec_taylors_complexes, matvec_complexes_taylors, &
              vecmat_taylors_complexes, vecmat_complexes_taylors, &
              matmul_taylors_scomplexes, matmul_scomplexes_taylors, &
              matvec_taylors_scomplexes, matvec_scomplexes_taylors, &
              vecmat_taylors_scomplexes, vecmat_scomplexes_taylors
   END INTERFACE

!!$ Interfaces for user-defined functions

   PUBLIC DERIVATIVE, EXPANSION, GRADIENT, INDEPENDENT, LAPLACIAN, &
          VALUE, REALVALUE, IMAGVALUE

   INTERFACE INDEPENDENT
      MODULE PROCEDURE independent_real, independent_sreal, &
              independent_complex, independent_scomplex, &
              independent_int
   END INTERFACE
   
!!$****************************************************************
!!$ MISCELLANEOUS
!!$****************************************************************
   
!!$ Value of Pi for use in trigonometric functions
   
   REAL(kind=dr_kind), PARAMETER :: Pi = 3.14159265358979323844_dp

!!$ Persistent workspace for FDB

   !INTEGER,SAVE :: scratchi(700,1+Max_taylor_order,Max_taylor_order)
   !INTEGER,SAVE :: scratchr(700,Max_taylor_order)
 
   INTEGER,SAVE :: scratchi(1000,1+Max_taylor_order,Max_taylor_order)
   INTEGER,SAVE :: scratchr(1000,Max_taylor_order)
   
   INTEGER,SAVE :: have(Max_taylor_order)
   DATA have/Max_taylor_order*0/

   !DEC$ IF .NOT. DEFINED (useintrinsic)
!!$ Emulate the needed parts of the IEEE intrinsic module
   INTEGER, PARAMETER, PRIVATE :: ieee_quiet_nan = 1
   
   INTERFACE IEEE_VALUE
       MODULE PROCEDURE IEEE_VALUE_DOUBLE, IEEE_VALUE_SINGLE
   END INTERFACE

   PRIVATE IEEE_VALUE

   !DEC$ ENDIF

   


!!$ Power operator that fixes the 0**0 NAG compiler problem
   INTERFACE OPERATOR(.POW.)
       MODULE PROCEDURE DCPOWI
   END INTERFACE

   PRIVATE OPERATOR(.POW.)

CONTAINS

  !DEC$ IF .NOT. DEFINED (useintrinsic)
!!$ Dummy version of IEEE support function to return a NaN value
   PURE FUNCTION IEEE_VALUE_DOUBLE(dummy,type)

       IMPLICIT NONE

       REAL(kind=dr_kind), INTENT(IN) :: dummy
       INTEGER, INTENT(IN) :: type
       REAL(kind=dr_kind) :: ieee_value_double

       ieee_value_double = 0.D0/0.D0

   END FUNCTION IEEE_VALUE_DOUBLE
   PURE FUNCTION IEEE_VALUE_SINGLE(dummy,type)

       IMPLICIT NONE

       REAL(kind=sr_kind), INTENT(IN) :: dummy
       INTEGER, INTENT(IN) :: type
       REAL(kind=sr_kind) :: ieee_value_single

       ieee_value_single = 0.0/0.0
       
   END FUNCTION IEEE_VALUE_SINGLE
   !DEC$ ENDIF
   
!!$ Power operator that fixes the 0**0 NAG compiler problem
   ELEMENTAL FUNCTION DCPOWI(z,n)

       IMPLICIT NONE

       COMPLEX(kind=dc_kind), INTENT(IN) :: z
       INTEGER, INTENT(IN) :: n
       COMPLEX(kind=dc_kind) :: dcpowi

       IF(n==0) THEN
          dcpowi = 1._dp
          RETURN
       ELSE
          dcpowi = z**n
          RETURN
       ENDIF

   END FUNCTION DCPOWI


!!$****************************************************************
!!$ OVERLOADED OPERATORS
!!$****************************************************************

!!$****************************************************************
!!$ Assign two taylors
!!$****************************************************************

   ELEMENTAL SUBROUTINE assign_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor), INTENT(OUT) :: self

!!$ Code starts here

      self%drv(:,:) = other%drv(:,:)

   END SUBROUTINE assign_taylors

!!$****************************************************************
!!$ Assign an INTEGER to a taylor (derivatives all zero)
!!$****************************************************************

   ELEMENTAL SUBROUTINE assign_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(OUT) :: self
      INTEGER, INTENT(IN) :: other

!!$ Code starts here

      self%drv(:,:) = 0
      self%drv(:,0) = other

   END SUBROUTINE assign_taylor_int

!!$****************************************************************
!!$ Assign a REAL to a taylor (derivatives all zero)
!!$****************************************************************

   ELEMENTAL SUBROUTINE assign_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(OUT) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other

!!$ Code starts here

      self%drv(:,:) = 0
      self%drv(:,0) = other

   END SUBROUTINE assign_taylor_real

!!$****************************************************************
!!$ Assign a default kind REAL to a taylor (derivatives all zero)
!!$****************************************************************

   ELEMENTAL SUBROUTINE assign_taylor_sreal(self,other)

      IMPLICIT NONE
      
!!$ Passed variables

      TYPE(taylor), INTENT(OUT) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other

!!$ Code starts here

      self%drv(:,:) = 0
      self%drv(:,0) = other

   END SUBROUTINE assign_taylor_sreal

!!$****************************************************************
!!$ Assign a COMPLEX to a taylor (derivatives all zero)
!!$****************************************************************

   ELEMENTAL SUBROUTINE assign_taylor_complex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN) :: other
      TYPE(taylor), INTENT(OUT) :: self

!!$ Code starts here

      self%drv(:,:) = 0
      self%drv(:,0) = other

   END SUBROUTINE assign_taylor_complex

!!$****************************************************************
!!$ Assign a default COMPLEX to a taylor (derivatives all zero)
!!$****************************************************************

   ELEMENTAL SUBROUTINE assign_taylor_scomplex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=sc_kind), INTENT(IN) :: other
      TYPE(taylor), INTENT(OUT) :: self

!!$ Code starts here

      self%drv(:,:) = 0
      self%drv(:,0) = other

      END SUBROUTINE assign_taylor_scomplex

!!$****************************************************************
!!$ Assign the value of a taylor to a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part

   ELEMENTAL SUBROUTINE assign_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=dr_kind), INTENT(OUT) :: self

!!$ Code starts here

      self = other%drv(1,0)

      END SUBROUTINE assign_real_taylor

!!$****************************************************************
!!$ Assign the value of a taylor to a default kind REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part

   ELEMENTAL SUBROUTINE assign_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=sr_kind), INTENT(OUT) :: self

!!$ Code starts here

      self = other%drv(1,0)

      END SUBROUTINE assign_sreal_taylor

!!$****************************************************************
!!$ Assign the value of a taylor to a COMPLEX
!!$****************************************************************

   ELEMENTAL SUBROUTINE assign_complex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      COMPLEX(kind=dc_kind), INTENT(OUT) :: self

!!$ Code starts here

      self = other%drv(1,0)

      END SUBROUTINE assign_complex_taylor

!!$****************************************************************
!!$ Assign the value of a taylor to a default kind COMPLEX
!!$****************************************************************

   ELEMENTAL SUBROUTINE assign_scomplex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      COMPLEX(kind=sc_kind), INTENT(OUT) :: self

!!$ Code starts here

      self = other%drv(1,0)

      END SUBROUTINE assign_scomplex_taylor

!!$****************************************************************
!!$ The negative of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION neg_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: neg_taylor

!!$ Code starts here

      neg_taylor%drv(:,:) = -self%drv(:,:)

   END FUNCTION neg_taylor

!!$****************************************************************
!!$ Add two taylors
!!$****************************************************************

   ELEMENTAL FUNCTION add_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      TYPE(taylor) :: add_taylors

!!$ Code starts here

      add_taylors%drv(:,:) = self%drv(:,:) + other%drv(:,:)

   END FUNCTION add_taylors

!!$****************************************************************
!!$ Add a taylor and a REAL
!!$****************************************************************

   ELEMENTAL FUNCTION add_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor) :: add_taylor_real

!!$ Code starts here

      add_taylor_real%drv(:,:) = self%drv(:,:)
      add_taylor_real%drv(:,0) = self%drv(:,0) + other

   END FUNCTION add_taylor_real

!!$****************************************************************
!!$ Add a REAL and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION add_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor) :: add_real_taylor

!!$ Code starts here

      add_real_taylor%drv(:,:) = other%drv(:,:)
      add_real_taylor%drv(:,0) = other%drv(:,0) + self

   END FUNCTION add_real_taylor

!!$****************************************************************
!!$ Add a taylor and a default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION add_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor) :: add_taylor_sreal

!!$ Code starts here

      add_taylor_sreal%drv(:,:) = self%drv(:,:)
      add_taylor_sreal%drv(:,0) = self%drv(:,0) + other

   END FUNCTION add_taylor_sreal

!!$****************************************************************
!!$ Add a default kind REAL and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION add_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor) :: add_sreal_taylor

!!$ Code starts here

      add_sreal_taylor%drv(:,:) = other%drv(:,:)
      add_sreal_taylor%drv(:,0) = other%drv(:,0) + self

   END FUNCTION add_sreal_taylor

!!$****************************************************************
!!$ Add a taylor and a COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION add_taylor_complex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=dc_kind), INTENT(IN) :: other
      TYPE(taylor) :: add_taylor_complex

!!$ Code starts here

      add_taylor_complex%drv(:,:) = self%drv(:,:)
      add_taylor_complex%drv(:,0) = self%drv(:,0) + other

   END FUNCTION add_taylor_complex

!!$****************************************************************
!!$ Add a COMPLEX and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION add_complex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      COMPLEX(kind=dc_kind), INTENT(IN) :: self
      TYPE(taylor) :: add_complex_taylor

!!$ Code starts here

      add_complex_taylor%drv(:,:) = other%drv(:,:)
      add_complex_taylor%drv(:,0) = other%drv(:,0) + self

   END FUNCTION add_complex_taylor

!!$****************************************************************
!!$ Add a taylor and a default kind COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION add_taylor_scomplex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=sc_kind), INTENT(IN) :: other
      TYPE(taylor) :: add_taylor_scomplex

!!$ Code starts here

      add_taylor_scomplex%drv(:,:) = self%drv(:,:)
      add_taylor_scomplex%drv(:,0) = self%drv(:,0) + other

   END FUNCTION add_taylor_scomplex

!!$****************************************************************
!!$ Add a default kind COMPLEX and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION add_scomplex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      COMPLEX(kind=sc_kind), INTENT(IN) :: self
      TYPE(taylor) :: add_scomplex_taylor

!!$ Code starts here

      add_scomplex_taylor%drv(:,:) = other%drv(:,:)
      add_scomplex_taylor%drv(:,0) = other%drv(:,0) + self

   END FUNCTION add_scomplex_taylor

!!$****************************************************************
!!$ Add a taylor and an INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION add_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: other
      TYPE(taylor) :: add_taylor_int

!!$ Code starts here

      add_taylor_int%drv(:,:) = self%drv(:,:)
      add_taylor_int%drv(:,0) = self%drv(:,0) + other

   END FUNCTION add_taylor_int

!!$****************************************************************
!!$ Add an INTEGER and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION add_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      INTEGER, INTENT(IN) :: self
      TYPE(taylor) :: add_int_taylor

!!$ Code starts here

      add_int_taylor%drv(:,:) = other%drv(:,:)
      add_int_taylor%drv(:,0) = other%drv(:,0) + self

   END FUNCTION add_int_taylor

!!$****************************************************************
!!$ Subtract two taylors
!!$****************************************************************

   ELEMENTAL FUNCTION sub_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      TYPE(taylor) :: sub_taylors

!!$ Code starts here

      sub_taylors%drv(:,:) = self%drv(:,:) - other%drv(:,:)

   END FUNCTION sub_taylors

!!$****************************************************************
!!$ Subtract a taylor and a REAL
!!$****************************************************************

   ELEMENTAL FUNCTION sub_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor) :: sub_taylor_real

!!$ Code starts here

      sub_taylor_real%drv(:,:) = self%drv(:,:)
      sub_taylor_real%drv(:,0) = self%drv(:,0) - other

   END FUNCTION sub_taylor_real

!!$****************************************************************
!!$ Subtract a REAL and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION sub_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor) :: sub_real_taylor

!!$ Code starts here

      sub_real_taylor%drv(:,:) = -other%drv(:,:)
      sub_real_taylor%drv(:,0) = self - other%drv(:,0)

   END FUNCTION sub_real_taylor

!!$****************************************************************
!!$ Subtract a taylor and a default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION sub_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor) :: sub_taylor_sreal

!!$ Code starts here

      sub_taylor_sreal%drv(:,:) = self%drv(:,:)
      sub_taylor_sreal%drv(:,0) = self%drv(:,0) - other

   END FUNCTION sub_taylor_sreal

!!$****************************************************************
!!$ Subtract a default kind REAL and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION sub_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor) :: sub_sreal_taylor

!!$ Code starts here

      sub_sreal_taylor%drv(:,:) = -other%drv(:,:)
      sub_sreal_taylor%drv(:,0) = self - other%drv(:,0)

   END FUNCTION sub_sreal_taylor

!!$****************************************************************
!!$ Subtract a taylor and a COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION sub_taylor_complex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=dc_kind), INTENT(IN) :: other
      TYPE(taylor) :: sub_taylor_complex

!!$ Code starts here

      sub_taylor_complex%drv(:,:) = self%drv(:,:)
      sub_taylor_complex%drv(:,0) = self%drv(:,0) - other

   END FUNCTION sub_taylor_complex

!!$****************************************************************
!!$ Subtract a COMPLEX and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION sub_complex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      COMPLEX(kind=dc_kind), INTENT(IN) :: self
      TYPE(taylor) :: sub_complex_taylor

!!$ Code starts here

      sub_complex_taylor%drv(:,:) = -other%drv(:,:)
      sub_complex_taylor%drv(:,0) = self - other%drv(:,0)

   END FUNCTION sub_complex_taylor

!!$****************************************************************
!!$ Subtract a taylor and a default kind COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION sub_taylor_scomplex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=sc_kind), INTENT(IN) :: other
      TYPE(taylor) :: sub_taylor_scomplex

!!$ Code starts here

      sub_taylor_scomplex%drv(:,:) = self%drv(:,:)
      sub_taylor_scomplex%drv(:,0) = self%drv(:,0) - other

   END FUNCTION sub_taylor_scomplex

!!$****************************************************************
!!$ Subtract a default kind COMPLEX and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION sub_scomplex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      COMPLEX(kind=sc_kind), INTENT(IN) :: self
      TYPE(taylor) :: sub_scomplex_taylor

!!$ Code starts here

      sub_scomplex_taylor%drv(:,:) = -other%drv(:,:)
      sub_scomplex_taylor%drv(:,0) = self - other%drv(:,0)

   END FUNCTION sub_scomplex_taylor

!!$****************************************************************
!!$ Subtract a taylor and an INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION sub_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: other
      TYPE(taylor) :: sub_taylor_int

!!$ Code starts here

      sub_taylor_int%drv(:,:) = self%drv(:,:)
      sub_taylor_int%drv(:,0) = self%drv(:,0) - other

   END FUNCTION sub_taylor_int

!!$****************************************************************
!!$ Subtract an INTEGER and a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION sub_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      INTEGER, INTENT(IN) :: self
      TYPE(taylor) :: sub_int_taylor

!!$ Code starts here

      sub_int_taylor%drv(:,:) = -other%drv(:,:)
      sub_int_taylor%drv(:,0) = self - other%drv(:,0)

   END FUNCTION sub_int_taylor

!!$****************************************************************
!!$ Multiply two taylors using Leibniz's rule
!!$****************************************************************

   ELEMENTAL FUNCTION mult_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      TYPE(taylor) :: mult_taylors

!!$ Internal variables

      INTEGER :: n,k

!!$ Code starts here

      mult_taylors%drv(:,:) = 0

      DO n=0,Taylor_order
       mult_taylors%drv(:,n) = 0
       DO k=0,n
        mult_taylors%drv(:,n) = mult_taylors%drv(:,n) + &
           choose(n,k)*self%drv(:,k)*other%drv(:,n-k)
       ENDDO
      ENDDO
           
   END FUNCTION mult_taylors

!!$****************************************************************
!!$ Multiply a taylor by a REAL
!!$****************************************************************

   ELEMENTAL FUNCTION mult_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN)  :: other
      TYPE(taylor) :: mult_taylor_real

!!$ Code starts here

      mult_taylor_real%drv(:,:) = self%drv(:,:) * other

   END FUNCTION mult_taylor_real

!!$****************************************************************
!!$ Multiply a taylor by a REAL
!!$****************************************************************

   ELEMENTAL FUNCTION mult_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN)  :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: mult_real_taylor

!!$ Code starts here

      mult_real_taylor%drv(:,:) = self * other%drv(:,:)

   END FUNCTION mult_real_taylor

!!$****************************************************************
!!$ Multiply a taylor by a default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION mult_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN)  :: other
      TYPE(taylor) :: mult_taylor_sreal

!!$ Code starts here

      mult_taylor_sreal%drv(:,:) = self%drv(:,:) * other

   END FUNCTION mult_taylor_sreal

!!$****************************************************************
!!$ Multiply a taylor by a default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION mult_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN)  :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: mult_sreal_taylor

!!$ Code starts here

      mult_sreal_taylor%drv(:,:) = self * other%drv(:,:)

   END FUNCTION mult_sreal_taylor

!!$****************************************************************
!!$ Multiply a taylor by a COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION mult_taylor_complex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=dc_kind), INTENT(IN)  :: other
      TYPE(taylor) :: mult_taylor_complex

!!$ Code starts here

      mult_taylor_complex%drv(:,:) = self%drv(:,:) * other

   END FUNCTION mult_taylor_complex

!!$****************************************************************
!!$ Multiply a taylor by a COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION mult_complex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN)  :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: mult_complex_taylor

!!$ Code starts here

      mult_complex_taylor%drv(:,:) = self * other%drv(:,:)

   END FUNCTION mult_complex_taylor


!!$****************************************************************
!!$ Multiply a taylor by a default kind COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION mult_taylor_scomplex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=sc_kind), INTENT(IN)  :: other
      TYPE(taylor) :: mult_taylor_scomplex

!!$ Code starts here

      mult_taylor_scomplex%drv(:,:) = self%drv(:,:) * other

   END FUNCTION mult_taylor_scomplex

!!$****************************************************************
!!$ Multiply a taylor by a default kind COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION mult_scomplex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=sc_kind), INTENT(IN)  :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: mult_scomplex_taylor

!!$ Code starts here

      mult_scomplex_taylor%drv(:,:) = self * other%drv(:,:)

   END FUNCTION mult_scomplex_taylor

!!$****************************************************************
!!$ Multiply a taylor by a INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION mult_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN)  :: other
      TYPE(taylor) :: mult_taylor_int

!!$ Code starts here

      mult_taylor_int%drv(:,:) = self%drv(:,:) * other

   END FUNCTION mult_taylor_int

!!$****************************************************************
!!$ Multiply a taylor by a INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION mult_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN)  :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: mult_int_taylor

!!$ Code starts here

      mult_int_taylor%drv(:,:) = self * other%drv(:,:)

   END FUNCTION mult_int_taylor

!!$****************************************************************
!!$ Divide two taylors by each other
!!$****************************************************************

   ELEMENTAL FUNCTION div_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      TYPE(taylor) :: div_taylors

!!$ Code starts here

      div_taylors = self*(1._dp/other)

   END FUNCTION div_taylors

!!$****************************************************************
!!$ Divide a REAL by a taylor (using Leibniz's and chain rule)
!!$****************************************************************

   ELEMENTAL FUNCTION div_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: div_real_taylor

!!$  Internal variables

      INTEGER :: n,k

!!$ Code starts here

      div_real_taylor%drv(:,:) = 0

      div_real_taylor%drv(:,0) = self/(other%drv(:,0))
      DO n=1,Taylor_order
       div_real_taylor%drv(:,n) = 0._dp
       DO k=0,n-1
        div_real_taylor%drv(:,n) = div_real_taylor%drv(:,n) - &
           choose(n,k)*other%drv(:,n-k)*div_real_taylor%drv(:,k)
       ENDDO
       div_real_taylor%drv(:,n) = div_real_taylor%drv(:,n)/(other%drv(:,0))
      ENDDO
      
   END FUNCTION div_real_taylor

!!$****************************************************************
!!$ Divide a taylor by a REAL
!!$****************************************************************

   ELEMENTAL FUNCTION div_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: div_taylor_real

!!$ Code starts here
      
      div_taylor_real%drv(:,:) = self%drv(:,:)/other

   END FUNCTION div_taylor_real

!!$****************************************************************
!!$ Divide a default kind REAL by a taylor (using rules as before)
!!$****************************************************************

   ELEMENTAL FUNCTION div_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: div_sreal_taylor

!!$  Internal variables

      INTEGER :: n,k

!!$ Code starts here

      div_sreal_taylor%drv(:,:) = 0

      div_sreal_taylor%drv(:,0) = self/(other%drv(:,0))
      DO n=1,Taylor_order
       div_sreal_taylor%drv(:,n) = 0._dp
       DO k=0,n-1
        div_sreal_taylor%drv(:,n) = div_sreal_taylor%drv(:,n) - &
           choose(n,k)*other%drv(:,n-k)*div_sreal_taylor%drv(:,k)
       ENDDO
       div_sreal_taylor%drv(:,n) = div_sreal_taylor%drv(:,n)/(other%drv(:,0))
      ENDDO
      
   END FUNCTION div_sreal_taylor

!!$****************************************************************
!!$ Divide a taylor by a default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION div_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: div_taylor_sreal

!!$ Code starts here
      
      div_taylor_sreal%drv(:,:) = self%drv(:,:)/other

   END FUNCTION div_taylor_sreal

!!$****************************************************************
!!$ Divide a COMPLEX by a taylor (using Leibniz's and chain rule)
!!$****************************************************************

   ELEMENTAL FUNCTION div_complex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: div_complex_taylor

 !!$  Internal variables

      INTEGER :: n,k

!!$ Code starts here

      div_complex_taylor%drv(:,:) = 0

      div_complex_taylor%drv(:,0) = self/(other%drv(:,0))
      DO n=1,Taylor_order
       div_complex_taylor%drv(:,n) = 0._dp
       DO k=0,n-1
        div_complex_taylor%drv(:,n) = div_complex_taylor%drv(:,n) - &
           choose(n,k)*other%drv(:,n-k)*div_complex_taylor%drv(:,k)
       ENDDO
       div_complex_taylor%drv(:,n) = div_complex_taylor%drv(:,n)/(other%drv(:,0))
      ENDDO

   END FUNCTION div_complex_taylor

!!$****************************************************************
!!$ Divide a taylor by a COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION div_taylor_complex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN) :: other
      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: div_taylor_complex

!!$ Code starts here
      
      div_taylor_complex%drv(:,:) = self%drv(:,:)/other

   END FUNCTION div_taylor_complex

!!$****************************************************************
!!$ Divide a default kind COMPLEX by a taylor (using chain rule)
!!$****************************************************************

   ELEMENTAL FUNCTION div_scomplex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=sc_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: div_scomplex_taylor

 !!$  Internal variables

      INTEGER :: n,k

!!$ Code starts here

      div_scomplex_taylor%drv(:,:) = 0

      div_scomplex_taylor%drv(:,0) = self/(other%drv(:,0))
      DO n=1,Taylor_order
       div_scomplex_taylor%drv(:,n) = 0._dp
       DO k=0,n-1
        div_scomplex_taylor%drv(:,n) = div_scomplex_taylor%drv(:,n) - &
           choose(n,k)*other%drv(:,n-k)*div_scomplex_taylor%drv(:,k)
       ENDDO
       div_scomplex_taylor%drv(:,n) = div_scomplex_taylor%drv(:,n)/(other%drv(:,0))
      ENDDO

   END FUNCTION div_scomplex_taylor

!!$****************************************************************
!!$ Divide a taylor by a default kind COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION div_taylor_scomplex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=sc_kind), INTENT(IN) :: other
      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: div_taylor_scomplex

!!$ Code starts here
      
      div_taylor_scomplex%drv(:,:) = self%drv(:,:)/other

   END FUNCTION div_taylor_scomplex

!!$****************************************************************
!!$ divide a taylor by an INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION div_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: other
      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: div_taylor_int

!!$ Code starts here
      
      div_taylor_int%drv(:,:) = self%drv(:,:)/other

   END FUNCTION div_taylor_int

!!$****************************************************************
!!$ Divide an INTEGER by a taylor (using Leibniz's and chain rule)
!!$****************************************************************

   ELEMENTAL FUNCTION div_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      TYPE(taylor) :: div_int_taylor

!!$  Internal variables

      INTEGER :: n,k

!!$ Code starts here

      div_int_taylor%drv(:,:) = 0

      div_int_taylor%drv(:,0) = self/(other%drv(:,0))
      DO n=1,Taylor_order
       div_int_taylor%drv(:,n) = 0._dp
       DO k=0,n-1
        div_int_taylor%drv(:,n) = div_int_taylor%drv(:,n) - &
           choose(n,k)*other%drv(:,n-k)*div_int_taylor%drv(:,k)
       ENDDO
       div_int_taylor%drv(:,n) = div_int_taylor%drv(:,n)/(other%drv(:,0))
      ENDDO
      
   END FUNCTION div_int_taylor

!!$****************************************************************
!!$ Take an integer power of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION power_taylor_int(self,n)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: n
      TYPE(taylor) :: power_taylor_int

!!$  Internal variables

      INTEGER :: k

!!$ Code starts here

      power_taylor_int = REAL(1._dp,dr_kind)
      DO k=1,n
         power_taylor_int = power_taylor_int * self
      ENDDO

   END FUNCTION power_taylor_int

!!$****************************************************************
!!$ Check equality of taylors
!!$****************************************************************

   ELEMENTAL FUNCTION eq_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      LOGICAL :: eq_taylors

!!$ Code starts here

      eq_taylors = (self%drv(1,0) == other%drv(1,0))

   END FUNCTION eq_taylors

!!$****************************************************************
!!$ Check equality of taylor and REAL
!!$****************************************************************

   ELEMENTAL FUNCTION eq_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      LOGICAL :: eq_taylor_real

!!$ Code starts here

      eq_taylor_real = (self%drv(1,0) == other)

   END FUNCTION eq_taylor_real

!!$****************************************************************
!!$ Check equality of REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION eq_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: eq_real_taylor

!!$ Code starts here

      eq_real_taylor = (self == other%drv(1,0))

   END FUNCTION eq_real_taylor

!!$****************************************************************
!!$ Check equality of taylor and default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION eq_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      LOGICAL :: eq_taylor_sreal

!!$ Code starts here

      eq_taylor_sreal = (self%drv(1,0) == other)

   END FUNCTION eq_taylor_sreal

!!$****************************************************************
!!$ Check equality of default kind REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION eq_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: eq_sreal_taylor

!!$ Code starts here

      eq_sreal_taylor = (self == other%drv(1,0))

   END FUNCTION eq_sreal_taylor

!!$****************************************************************
!!$ Check equality of taylor and COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION eq_taylor_complex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=dc_kind), INTENT(IN) :: other
      LOGICAL :: eq_taylor_complex

!!$ Code starts here

      eq_taylor_complex = (self%drv(1,0) == other)

   END FUNCTION eq_taylor_complex

!!$****************************************************************
!!$ Check equality of COMPLEX and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION eq_complex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: eq_complex_taylor

!!$ Code starts here

      eq_complex_taylor = (self == other%drv(1,0))

   END FUNCTION eq_complex_taylor

!!$****************************************************************
!!$ Check equality of taylor and default kind COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION eq_taylor_scomplex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=sc_kind), INTENT(IN) :: other
      LOGICAL :: eq_taylor_scomplex

!!$ Code starts here

      eq_taylor_scomplex = (self%drv(1,0) == other)

   END FUNCTION eq_taylor_scomplex

!!$****************************************************************
!!$ Check equality of default kind COMPLEX and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION eq_scomplex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=sc_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: eq_scomplex_taylor

!!$ Code starts here

      eq_scomplex_taylor = (self == other%drv(1,0))

   END FUNCTION eq_scomplex_taylor

!!$****************************************************************
!!$ Check equality of taylor and INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION eq_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: other
      LOGICAL :: eq_taylor_int

!!$ Code starts here

      eq_taylor_int = (self%drv(1,0) == other)

   END FUNCTION eq_taylor_int

!!$****************************************************************
!!$ Check equality of INTEGER and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION eq_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: eq_int_taylor

!!$ Code starts here

      eq_int_taylor = (self == other%drv(1,0))

   END FUNCTION eq_int_taylor

!!$****************************************************************
!!$ Check inequality of taylors
!!$****************************************************************

   ELEMENTAL FUNCTION ne_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      LOGICAL :: ne_taylors

!!$ Code starts here

      ne_taylors = (self%drv(1,0) /= other%drv(1,0))

   END FUNCTION ne_taylors

!!$****************************************************************
!!$ Check inequality of taylor and REAL
!!$****************************************************************

   ELEMENTAL FUNCTION ne_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      LOGICAL :: ne_taylor_real

!!$ Code starts here

      ne_taylor_real = (self%drv(1,0) /= other)

   END FUNCTION ne_taylor_real

!!$****************************************************************
!!$ Check inequality of REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION ne_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: ne_real_taylor

!!$ Code starts here

      ne_real_taylor = (self /= other%drv(1,0))

   END FUNCTION ne_real_taylor

!!$****************************************************************
!!$ Check inequality of taylor and default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION ne_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      LOGICAL :: ne_taylor_sreal

!!$ Code starts here

      ne_taylor_sreal = (self%drv(1,0) /= other)

   END FUNCTION ne_taylor_sreal

!!$****************************************************************
!!$ Check inequality of default kind REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION ne_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: ne_sreal_taylor

!!$ Code starts here

      ne_sreal_taylor = (self /= other%drv(1,0))

   END FUNCTION ne_sreal_taylor

!!$****************************************************************
!!$ Check inequality of taylor and COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION ne_taylor_complex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=dc_kind), INTENT(IN) :: other
      LOGICAL :: ne_taylor_complex

!!$ Code starts here

      ne_taylor_complex = (self%drv(1,0) /= other)

   END FUNCTION ne_taylor_complex

!!$****************************************************************
!!$ Check inequality of COMPLEX and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION ne_complex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: ne_complex_taylor

!!$ Code starts here

      ne_complex_taylor = (self /= other%drv(1,0))

   END FUNCTION ne_complex_taylor

!!$****************************************************************
!!$ Check inequality of taylor and default kind COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION ne_taylor_scomplex(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=sc_kind), INTENT(IN) :: other
      LOGICAL :: ne_taylor_scomplex

!!$ Code starts here

      ne_taylor_scomplex = (self%drv(1,0) /= other)

   END FUNCTION ne_taylor_scomplex

!!$****************************************************************
!!$ Check inequality of default kind COMPLEX and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION ne_scomplex_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=sc_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: ne_scomplex_taylor

!!$ Code starts here

      ne_scomplex_taylor = (self /= other%drv(1,0))

   END FUNCTION ne_scomplex_taylor

!!$****************************************************************
!!$ Check inequality of taylor and INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION ne_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: other
      LOGICAL :: ne_taylor_int

!!$ Code starts here

      ne_taylor_int = (self%drv(1,0) /= other)

   END FUNCTION ne_taylor_int

!!$****************************************************************
!!$ Check inequality of INTEGER and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION ne_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: ne_int_taylor

!!$ Code starts here

      ne_int_taylor = (self /= other%drv(1,0))

   END FUNCTION ne_int_taylor

!!$****************************************************************
!!$ Check ordering of taylors
!!$****************************************************************

   ELEMENTAL FUNCTION lt_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      LOGICAL :: lt_taylors

!!$ Code starts here

      lt_taylors = (REAL(self%drv(1,0)) < REAL(other%drv(1,0)))

   END FUNCTION lt_taylors

!!$****************************************************************
!!$ Check ordering of taylor and REAL
!!$****************************************************************

   ELEMENTAL FUNCTION lt_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      LOGICAL :: lt_taylor_real

!!$ Code starts here

      lt_taylor_real = (REAL(self%drv(1,0)) < other)

   END FUNCTION lt_taylor_real

!!$****************************************************************
!!$ Check ordering of REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION lt_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: lt_real_taylor

!!$ Code starts here

      lt_real_taylor = (self < REAL(other%drv(1,0)))

   END FUNCTION lt_real_taylor

!!$****************************************************************
!!$ Check ordering of taylor and default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION lt_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      LOGICAL :: lt_taylor_sreal

!!$ Code starts here

      lt_taylor_sreal = (REAL(self%drv(1,0)) < other)

   END FUNCTION lt_taylor_sreal

!!$****************************************************************
!!$ Check ordering of default kind REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION lt_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: lt_sreal_taylor

!!$ Code starts here

      lt_sreal_taylor = (self < REAL(other%drv(1,0)))

   END FUNCTION lt_sreal_taylor

!!$****************************************************************
!!$ Check ordering of taylor and INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION lt_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: other
      LOGICAL :: lt_taylor_int

!!$ Code starts here

      lt_taylor_int = (REAL(self%drv(1,0)) < other)

   END FUNCTION lt_taylor_int

!!$****************************************************************
!!$ Check ordering of INTEGER and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION lt_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: lt_int_taylor

!!$ Code starts here

      lt_int_taylor = (self < REAL(other%drv(1,0)))

   END FUNCTION lt_int_taylor

!!$****************************************************************
!!$ Check ordering of taylors
!!$****************************************************************

   ELEMENTAL FUNCTION le_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      LOGICAL :: le_taylors

!!$ Code starts here

      le_taylors = (REAL(self%drv(1,0)) <= REAL(other%drv(1,0)))

   END FUNCTION le_taylors

!!$****************************************************************
!!$ Check ordering of taylor and REAL
!!$****************************************************************

   ELEMENTAL FUNCTION le_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      LOGICAL :: le_taylor_real

!!$ Code starts here

      le_taylor_real = (REAL(self%drv(1,0)) <= other)

   END FUNCTION le_taylor_real

!!$****************************************************************
!!$ Check ordering of REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION le_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: le_real_taylor

!!$ Code starts here

      le_real_taylor = (self <= REAL(other%drv(1,0)))

   END FUNCTION le_real_taylor

!!$****************************************************************
!!$ Check ordering of taylor and default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION le_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      LOGICAL :: le_taylor_sreal

!!$ Code starts here

      le_taylor_sreal = (REAL(self%drv(1,0)) <= other)

   END FUNCTION le_taylor_sreal

!!$****************************************************************
!!$ Check ordering of default kind REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION le_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: le_sreal_taylor

!!$ Code starts here

      le_sreal_taylor = (self <= REAL(other%drv(1,0)))

   END FUNCTION le_sreal_taylor

!!$****************************************************************
!!$ Check ordering of taylor and INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION le_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: other
      LOGICAL :: le_taylor_int

!!$ Code starts here

      le_taylor_int = (REAL(self%drv(1,0)) <= other)

   END FUNCTION le_taylor_int

!!$****************************************************************
!!$ Check ordering of INTEGER and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION le_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: le_int_taylor

!!$ Code starts here

      le_int_taylor = (self <= REAL(other%drv(1,0)))

   END FUNCTION le_int_taylor

!!$****************************************************************
!!$ Check ordering of taylors
!!$****************************************************************

   ELEMENTAL FUNCTION gt_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      LOGICAL :: gt_taylors

!!$ Code starts here

      gt_taylors = (REAL(self%drv(1,0)) > REAL(other%drv(1,0)))

   END FUNCTION gt_taylors

!!$****************************************************************
!!$ Check ordering of taylor and REAL
!!$****************************************************************

   ELEMENTAL FUNCTION gt_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      LOGICAL :: gt_taylor_real

!!$ Code starts here

      gt_taylor_real = (REAL(self%drv(1,0)) > other)

   END FUNCTION gt_taylor_real

!!$****************************************************************
!!$ Check ordering of REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION gt_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: gt_real_taylor

!!$ Code starts here

      gt_real_taylor = (self > REAL(other%drv(1,0)))

   END FUNCTION gt_real_taylor

!!$****************************************************************
!!$ Check ordering of taylor and default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION gt_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      LOGICAL :: gt_taylor_sreal

!!$ Code starts here

      gt_taylor_sreal = (REAL(self%drv(1,0)) > other)

   END FUNCTION gt_taylor_sreal

!!$****************************************************************
!!$ Check ordering of default kind REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION gt_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: gt_sreal_taylor

!!$ Code starts here

      gt_sreal_taylor = (self > REAL(other%drv(1,0)))

   END FUNCTION gt_sreal_taylor

!!$****************************************************************
!!$ Check ordering of taylor and INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION gt_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: other
      LOGICAL :: gt_taylor_int

!!$ Code starts here

      gt_taylor_int = (REAL(self%drv(1,0)) > other)

   END FUNCTION gt_taylor_int

!!$****************************************************************
!!$ Check ordering of INTEGER and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION gt_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: gt_int_taylor

!!$ Code starts here

      gt_int_taylor = (self > REAL(other%drv(1,0)))

   END FUNCTION gt_int_taylor

!!$****************************************************************
!!$ Check ordering of taylors
!!$****************************************************************

   ELEMENTAL FUNCTION ge_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      LOGICAL :: ge_taylors

!!$ Code starts here

      ge_taylors = (REAL(self%drv(1,0)) >= REAL(other%drv(1,0)))

   END FUNCTION ge_taylors

!!$****************************************************************
!!$ Check ordering of taylor and REAL
!!$****************************************************************

   ELEMENTAL FUNCTION ge_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      LOGICAL :: ge_taylor_real

!!$ Code starts here

      ge_taylor_real = (REAL(self%drv(1,0)) >= other)

   END FUNCTION ge_taylor_real

!!$****************************************************************
!!$ Check ordering of REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION ge_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: ge_real_taylor

!!$ Code starts here

      ge_real_taylor = (self >= REAL(other%drv(1,0)))

   END FUNCTION ge_real_taylor

!!$****************************************************************
!!$ Check ordering of taylor and default kind REAL
!!$****************************************************************

   ELEMENTAL FUNCTION ge_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      LOGICAL :: ge_taylor_sreal

!!$ Code starts here

      ge_taylor_sreal = (REAL(self%drv(1,0)) >= other)

   END FUNCTION ge_taylor_sreal

!!$****************************************************************
!!$ Check ordering of default kind REAL and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION ge_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: ge_sreal_taylor

!!$ Code starts here

      ge_sreal_taylor = (self >= REAL(other%drv(1,0)))

   END FUNCTION ge_sreal_taylor

!!$****************************************************************
!!$ Check ordering of taylor and INTEGER
!!$****************************************************************

   ELEMENTAL FUNCTION ge_taylor_int(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: other
      LOGICAL :: ge_taylor_int

!!$ Code starts here

      ge_taylor_int = (REAL(self%drv(1,0)) >= other)

   END FUNCTION ge_taylor_int

!!$****************************************************************
!!$ Check ordering of INTEGER and taylor
!!$****************************************************************

   ELEMENTAL FUNCTION ge_int_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: self
      TYPE(taylor), INTENT(IN) :: other
      LOGICAL :: ge_int_taylor

!!$ Code starts here

      ge_int_taylor = (self >= REAL(other%drv(1,0)))

   END FUNCTION ge_int_taylor


!!$****************************************************************
!!$ USER-DEFINED OPERATORS
!!$****************************************************************

!!$****************************************************************
!!$ Check identity of taylors
!!$****************************************************************

   ELEMENTAL FUNCTION ident_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      LOGICAL :: ident_taylors

!!$ Code starts here

      ident_taylors = ALL(self%drv == other%drv)

   END FUNCTION ident_taylors

!!$****************************************************************
!!$ Check non-identity of taylors
!!$****************************************************************

   ELEMENTAL FUNCTION nident_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      LOGICAL :: nident_taylors

!!$ Code starts here

      nident_taylors = ANY(self%drv /= other%drv)

   END FUNCTION nident_taylors


!!$****************************************************************
!!$ FORTRAN INTRINSICS
!!$****************************************************************

!!$****************************************************************
!!$ Take the square root of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION sqrt_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: sqrt_taylor

!!$  Internal variables

      INTEGER :: n,k
      COMPLEX(kind=dc_kind) :: sqrtterm

!!$ Code starts here

      sqrt_taylor%drv(:,:) = 0

      sqrtterm = SQRT(self%drv(1,0))
      sqrt_taylor%drv(:,0) = sqrtterm

      DO n=1,Taylor_order
       sqrt_taylor%drv(:,n) = self%drv(:,n)
       DO k=1,n/2
        sqrt_taylor%drv(:,n) = sqrt_taylor%drv(:,n) - &
           2.*choose(n,k)*sqrt_taylor%drv(:,n-k)*sqrt_taylor%drv(:,k)
       ENDDO
       IF(MODULO(n,2)==0) sqrt_taylor%drv(:,n) = sqrt_taylor%drv(:,n) + &
           choose(n,n/2)*sqrt_taylor%drv(:,n/2)**2
       sqrt_taylor%drv(:,n) = sqrt_taylor%drv(:,n)/(2.*sqrtterm)
      ENDDO

   END FUNCTION sqrt_taylor

!!$****************************************************************
!!$ Take the imaginary part of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION aimag_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: aimag_taylor

!!$ Code starts here

      aimag_taylor%drv(:,:) = AIMAG(self%drv(:,:))

   END FUNCTION aimag_taylor

!!$****************************************************************
!!$ Take the real part of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION real_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: real_taylor

!!$ Code starts here

      real_taylor%drv(:,:) = REAL(self%drv(:,:),dr_kind)

   END FUNCTION real_taylor

!!$****************************************************************
!!$ Take the exponential of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION exp_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: exp_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      COMPLEX(kind=dc_kind) :: external_drv(Max_taylor_order)

!!$ Code starts here

      exp_taylor%drv(:,:) = 0
      k(:) = 0

      exp_taylor%drv(:,0) = EXP(self%drv(:,0))
      external_drv(:) = exp_taylor%drv(1,0)

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         exp_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) exp_taylor%drv(i,p) = exp_taylor%drv(i,p) + &
           wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION exp_taylor

!!$****************************************************************
!!$ Take the natural logarithm of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION log_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: log_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      COMPLEX(kind=dc_kind) :: external_drv(Max_taylor_order)

!!$ Code starts here

      log_taylor%drv(:,:) = 0
      k(:) = 0

      log_taylor%drv(:,0) = LOG(self%drv(:,0))
      external_drv(:) = (/ ((-1)**(i-1)*fac(i-1)/self%drv(1,0)**i,&
                             i=1,Max_taylor_order) /)

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         log_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) log_taylor%drv(i,p) = log_taylor%drv(i,p) + &
          wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION log_taylor

!!$****************************************************************
!!$ Take the dekadic logarithm of a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION log10_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: log10_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      REAL(kind=dr_kind) :: external_drv(Max_taylor_order)

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        log10_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      log10_taylor%drv(:,:) = 0
      k(:) = 0

      log10_taylor%drv(:,0) = LOG10(REAL(self%drv(:,0),dr_kind))
      external_drv(:) = (/ ((-1)**(i-1)*fac(i-1)/self%drv(1,0)**i,&
                             i=1,Max_taylor_order) /)/LOG(10._dp)

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         log10_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) log10_taylor%drv(i,p) = log10_taylor%drv(i,p) + &
          wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION log10_taylor

!!$****************************************************************
!!$ Take the cosine of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION cos_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: cos_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      COMPLEX(kind=dc_kind) :: external_drv(Max_taylor_order)

!!$ Code starts here

      cos_taylor%drv(:,:) = 0
      k(:) = 0

      cos_taylor%drv(:,0) = COS(self%drv(:,0))
      external_drv(:) = (/ (COS(self%drv(1,0)+i*Pi/2._dp),&
                             i=1,Max_taylor_order) /)

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         cos_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) cos_taylor%drv(i,p) = cos_taylor%drv(i,p) + &
          wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION cos_taylor

!!$****************************************************************
!!$ Take the sine of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION sin_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: sin_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      COMPLEX(kind=dc_kind) :: external_drv(Max_taylor_order)

!!$ Code starts here

      sin_taylor%drv(:,:) = 0
      k(:) = 0

      sin_taylor%drv(:,0) = SIN(self%drv(:,0))
      external_drv(:) = (/ (SIN(self%drv(1,0)+i*Pi/2._dp),&
                             i=1,Max_taylor_order) /)

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         sin_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) sin_taylor%drv(i,p) = sin_taylor%drv(i,p) + &
          wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION sin_taylor

!!$****************************************************************
!!$ Take the tangent of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION tan_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: tan_taylor

!!$ Code starts here

      tan_taylor = SIN(self)/COS(self)

   END FUNCTION tan_taylor

!!$****************************************************************
!!$ Take the hyperbolic cosine of a COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION cosh_complex(self)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN) :: self
      COMPLEX(kind=dc_kind) :: cosh_complex

!!$ Code starts here

      cosh_complex = CMPLX(COSH(REAL(self,dr_kind))*COS(AIMAG(self)),&
                           SINH(REAL(self,dr_kind))*SIN(AIMAG(self)))

   END FUNCTION cosh_complex

!!$****************************************************************
!!$ Take the hyperbolic sine of a COMPLEX
!!$****************************************************************

   ELEMENTAL FUNCTION sinh_complex(self)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN) :: self
      COMPLEX(kind=dc_kind) :: sinh_complex

!!$ Code starts here

      sinh_complex = CMPLX(SINH(REAL(self,dr_kind))*COS(AIMAG(self)),&
                           COSH(REAL(self,dr_kind))*SIN(AIMAG(self)))

   END FUNCTION sinh_complex

!!$****************************************************************
!!$ Take the hyperbolic cosine of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION cosh_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: cosh_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      COMPLEX(kind=dc_kind) :: external_drv(Max_taylor_order)

!!$ Code starts here

      cosh_taylor%drv(:,:) = 0
      k(:) = 0

      cosh_taylor%drv(:,0) = cosh_complex(self%drv(:,0))
      external_drv(:) = (/ (sinh_complex(self%drv(1,0))*(1+(-1)**(i-1))+ &
                            cosh_complex(self%drv(1,0))*(1+(-1)**i), &
                             i=1,Max_taylor_order) /)/2._dp

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         cosh_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) cosh_taylor%drv(i,p) = cosh_taylor%drv(i,p) + &
          wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION cosh_taylor

!!$****************************************************************
!!$ Take the hyperbolic sine of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION sinh_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: sinh_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      COMPLEX(kind=dc_kind) :: external_drv(Max_taylor_order)

!!$ Code starts here

      sinh_taylor%drv(:,:) = 0
      k(:) = 0

      sinh_taylor%drv(:,0) = sinh_complex(self%drv(:,0))
      external_drv(:) = (/ (sinh_complex(self%drv(1,0))*(1+(-1)**i)+ &
                            cosh_complex(self%drv(1,0))*(1+(-1)**(i-1)), &
                             i=1,Max_taylor_order) /)/2._dp

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         sinh_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) sinh_taylor%drv(i,p) = sinh_taylor%drv(i,p) + &
          wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION sinh_taylor

!!$****************************************************************
!!$ Take the hyperbolic tangent of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION tanh_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: tanh_taylor

!!$ Code starts here

      tanh_taylor = SINH(self)/COSH(self)

   END FUNCTION tanh_taylor

!!$****************************************************************
!!$ Take the arccosine of a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION acos_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: acos_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      REAL(kind=dr_kind) :: external_drv(Max_taylor_order)
      TYPE(taylor) :: auxiliary

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        acos_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      acos_taylor%drv(:,:) = 0
      k(:) = 0

      acos_taylor%drv(:,0) = ACOS(REAL(self%drv(:,0),dr_kind))
      auxiliary%drv(:,:) = 0._dp
      auxiliary%drv(:,0) = self%drv(:,0)
      auxiliary%drv(1,1) = 1._dp
      auxiliary = -1/SQRT(1-auxiliary**2) 
      external_drv(:) = (/ (auxiliary%drv(1,i-1),i=1,Max_taylor_order) /)

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         acos_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) acos_taylor%drv(i,p) = acos_taylor%drv(i,p) + &
          wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION acos_taylor

!!$****************************************************************
!!$ Take the arcsine of a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION asin_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: asin_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      REAL(kind=dr_kind) :: external_drv(Max_taylor_order)
      TYPE(taylor) :: auxiliary

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        asin_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      asin_taylor%drv(:,:) = 0
      k(:) = 0

      asin_taylor%drv(:,0) = ASIN(REAL(self%drv(:,0),dr_kind))
      auxiliary%drv(:,:) = 0._dp
      auxiliary%drv(:,0) = self%drv(:,0)
      auxiliary%drv(1,1) = 1._dp
      auxiliary = 1/SQRT(1-auxiliary**2)
      external_drv(:) = (/ (auxiliary%drv(1,i-1),i=1,Max_taylor_order) /)

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         asin_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) asin_taylor%drv(i,p) = asin_taylor%drv(i,p) + &
          wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION asin_taylor

!!$****************************************************************
!!$ Take the arctangent of a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION atan_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: atan_taylor

!!$ Internal variables

      INTEGER :: p,m,i
      INTEGER, POINTER :: idx(:,:),wgt(:)
      INTEGER :: k(Max_taylor_order)
      REAL(kind=dr_kind) :: external_drv(Max_taylor_order)
      TYPE(taylor) :: auxiliary

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        atan_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      atan_taylor%drv(:,:) = 0
      k(:) = 0

      atan_taylor%drv(:,0) = ATAN(REAL(self%drv(:,0),dr_kind))
      auxiliary%drv(:,:) = 0._dp
      auxiliary%drv(:,0) = self%drv(:,0)
      auxiliary%drv(1,1) = 1._dp
      auxiliary = 1/(1+auxiliary**2)
      external_drv(:) = (/ (auxiliary%drv(1,i-1),i=1,Max_taylor_order) /)

      DO p=1,Taylor_order

         CALL fdb(p,idx,wgt)

         atan_taylor%drv(:,p) = 0._dp

         DO m=1,SIZE(wgt)
          k(1:p) = idx(m,2:p+1)
          FORALL(i=1:N_taylor_vars) atan_taylor%drv(i,p) = atan_taylor%drv(i,p) + &
          wgt(m)*external_drv(idx(m,1))*PRODUCT(self%drv(i,1:p).POW.k(1:p),mask=k(1:p)>0)
         ENDDO
         
         DEALLOCATE(idx,wgt)

      ENDDO

   END FUNCTION atan_taylor

!!$****************************************************************
!!$ Take the arctangent of a pair of taylors
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION atan2_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      TYPE(taylor) :: atan2_taylors

!!$ Code starts here

      IF(Real_args_warn.AND.(ABS(AIMAG(self)).GT.Real_args_tol&
     &   .OR.ABS(AIMAG(other)).GT.Real_args_tol)) THEN
        atan2_taylors%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(other>0.) THEN
         atan2_taylors = ATAN(self/other)
      ELSE
         atan2_taylors = SIGN(Pi,REAL(self%drv(1,0),dr_kind)) + &
              ATAN(self/other)
      ENDIF

   END FUNCTION atan2_taylors

!!$****************************************************************
!!$ Take the pairwise arctangent of a taylor and a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION atan2_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor) :: atan2_taylor_real

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        atan2_taylor_real%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(other>0.) THEN
         atan2_taylor_real = ATAN(self/other)
      ELSE
         atan2_taylor_real = SIGN(Pi,REAL(self%drv(1,0),dr_kind)) + &
              ATAN(self/other)
      ENDIF

   END FUNCTION atan2_taylor_real

!!$****************************************************************
!!$ Take the pairwise arctangent of a taylor and a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION atan2_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor) :: atan2_real_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        atan2_real_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(other>0.) THEN
         atan2_real_taylor = ATAN(self/other)
      ELSE
         atan2_real_taylor = SIGN(Pi,self) + ATAN(self/other)
      ENDIF

   END FUNCTION atan2_real_taylor

!!$****************************************************************
!!$ Take the pairwise arctangent of a taylor and a default REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION atan2_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor) :: atan2_taylor_sreal

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        atan2_taylor_sreal%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(other>0.) THEN
         atan2_taylor_sreal = ATAN(self/other)
      ELSE
         atan2_taylor_sreal = SIGN(Pi,REAL(self%drv(1,0),dr_kind)) + &
              ATAN(self/other)
      ENDIF

   END FUNCTION atan2_taylor_sreal

!!$****************************************************************
!!$ Take the pairwise arctangent of a taylor and a default REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION atan2_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor) :: atan2_sreal_taylor

!!$ Code starts here
 
      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        atan2_sreal_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(other>0.) THEN
         atan2_sreal_taylor = ATAN(self/other)
      ELSE
         atan2_sreal_taylor = SIGN(REAL(Pi,sr_kind),self) + ATAN(self/other)
      ENDIF

   END FUNCTION atan2_sreal_taylor

!!$****************************************************************
!!$ Take the absolute value of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION abs_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: abs_taylor

!!$ Code starts here

      abs_taylor = SQRT(CONJG(self)*self)

   END FUNCTION abs_taylor

!!$****************************************************************
!!$ Take the SIGN function of two taylors
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION sign_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self, other
      TYPE(taylor) :: sign_taylors

!!$ Code starts here

      IF(Real_args_warn.AND.(ABS(AIMAG(self)).GT.Real_args_tol&
     &   .OR.ABS(AIMAG(other)).GT.Real_args_tol)) THEN
        sign_taylors%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      sign_taylors = ABS(self)*sign(1._dp,REAL(other%drv(1,0),dr_kind))

      IF(other==0.) sign_taylors%drv(:,1:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)

   END FUNCTION sign_taylors

!!$****************************************************************
!!$ Take the SIGN function of a taylor and a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION sign_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor) :: sign_taylor_real

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        sign_taylor_real%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      sign_taylor_real = ABS(self)*sign(1._dp,other)

   END FUNCTION sign_taylor_real

!!$****************************************************************
!!$ Take the SIGN function of a REAL and a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION sign_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=dr_kind), INTENT(IN) :: self
      REAL(kind=dr_kind) :: sign_real_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        sign_real_taylor = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      sign_real_taylor = ABS(self)*sign(1._dp,REAL(other%drv(1,0),dr_kind))

   END FUNCTION sign_real_taylor

!!$****************************************************************
!!$ Take the SIGN function of a taylors and a default kind REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION sign_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor) :: sign_taylor_sreal

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        sign_taylor_sreal%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      sign_taylor_sreal = ABS(self)*sign(1.,other)

   END FUNCTION sign_taylor_sreal

!!$****************************************************************
!!$ Take the SIGN function of a default kind REAL and a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION sign_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=sr_kind), INTENT(IN) :: self
      REAL(kind=sr_kind) :: sign_sreal_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        sign_sreal_taylor = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      sign_sreal_taylor = ABS(self)*sign(1._dp,REAL(other%drv(1,0),dr_kind))

   END FUNCTION sign_sreal_taylor

!!$****************************************************************
!!$ Take the greater of two taylors
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION max_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self, other
      TYPE(taylor) :: max_taylors

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.(ABS(AIMAG(self)).GT.Real_args_tol&
     &   .OR.ABS(AIMAG(other)).GT.Real_args_tol)) THEN
        max_taylors%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self>other) THEN
         max_taylors = self
      ELSE
         max_taylors = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(self%drv(j,i)/=other%drv(j,i)) max_taylors%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION max_taylors

!!$****************************************************************
!!$ Take the greater of a taylor and a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION max_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor) :: max_taylor_real

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        max_taylor_real%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self>other) THEN
         max_taylor_real = self
      ELSE
         max_taylor_real = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(self%drv(j,i)/=0.) max_taylor_real%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION max_taylor_real

!!$****************************************************************
!!$ Take the greater of a taylor and a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION max_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor) :: max_real_taylor

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        max_real_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self>other) THEN
         max_real_taylor = self
      ELSE
         max_real_taylor = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(other%drv(j,i)/=0.) max_real_taylor%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION max_real_taylor

!!$****************************************************************
!!$ Take the greater of a taylor and a default kind REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION max_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor) :: max_taylor_sreal

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        max_taylor_sreal%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self>other) THEN
         max_taylor_sreal = self
      ELSE
         max_taylor_sreal = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(self%drv(j,i)/=0.) max_taylor_sreal%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION max_taylor_sreal

!!$****************************************************************
!!$ Take the lesser of a taylor and a default kind REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION max_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor) :: max_sreal_taylor

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        max_sreal_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self>other) THEN
         max_sreal_taylor = self
      ELSE
         max_sreal_taylor = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(other%drv(j,i)/=0.) max_sreal_taylor%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION max_sreal_taylor

!!$****************************************************************
!!$ Take the lesser of two taylors
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION min_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self, other
      TYPE(taylor) :: min_taylors

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.(ABS(AIMAG(self)).GT.Real_args_tol&
     &   .OR.ABS(AIMAG(other)).GT.Real_args_tol)) THEN
        min_taylors%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self<other) THEN
         min_taylors = self
      ELSE
         min_taylors = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(self%drv(j,i)/=other%drv(j,i)) min_taylors%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION min_taylors

!!$****************************************************************
!!$ Take the lesser of a taylor and a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION min_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor) :: min_taylor_real

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        min_taylor_real%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self<other) THEN
         min_taylor_real = self
      ELSE
         min_taylor_real = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(self%drv(j,i)/=0.) min_taylor_real%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION min_taylor_real

!!$****************************************************************
!!$ Take the lesser of a taylor and a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION min_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor) :: min_real_taylor

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        min_real_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self<other) THEN
         min_real_taylor = self
      ELSE
         min_real_taylor = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(other%drv(j,i)/=0.) min_real_taylor%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION min_real_taylor

!!$****************************************************************
!!$ Take the lesser of a taylor and a default kind REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION min_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor) :: min_taylor_sreal

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        min_taylor_sreal%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self<other) THEN
         min_taylor_sreal = self
      ELSE
         min_taylor_sreal = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(self%drv(j,i)/=0.) min_taylor_sreal%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION min_taylor_sreal

!!$****************************************************************
!!$ Take the lesser of a taylor and a default kind REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION min_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor) :: min_sreal_taylor

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        min_sreal_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      IF(self<other) THEN
         min_sreal_taylor = self
      ELSE
         min_sreal_taylor = other
      ENDIF

      IF(self==other) THEN
        DO j=1,N_taylor_vars ; DO i=1,Taylor_order ;
           IF(other%drv(j,i)/=0.) min_sreal_taylor%drv(j,i:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)
        ENDDO ; ENDDO
      ENDIF

   END FUNCTION min_sreal_taylor

!!$****************************************************************
!!$ Take the DIM function of two taylors
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION dim_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self,other
      TYPE(taylor) :: dim_taylors

!!$ Code starts here

      IF(Real_args_warn.AND.(ABS(AIMAG(self)).GT.Real_args_tol&
         .OR.ABS(AIMAG(other)).GT.Real_args_tol)) THEN
        dim_taylors%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      dim_taylors = MAX(self-other,0.)

   END FUNCTION dim_taylors

!!$****************************************************************
!!$ Take the DIM function of a taylor and a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION dim_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor) :: dim_taylor_real

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        dim_taylor_real%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      dim_taylor_real = MAX(self-other,0.)

   END FUNCTION dim_taylor_real

!!$****************************************************************
!!$ Take the DIM function of a taylor and a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION dim_real_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=dr_kind), INTENT(IN) :: self
      TYPE(taylor) :: dim_real_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        dim_real_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      dim_real_taylor = MAX(self-other,0.)

   END FUNCTION dim_real_taylor

!!$****************************************************************
!!$ Take the DIM function of a taylor and a default kind REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION dim_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor) :: dim_taylor_sreal

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        dim_taylor_sreal%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      dim_taylor_sreal = MAX(self-other,0.)

   END FUNCTION dim_taylor_sreal

!!$****************************************************************
!!$ Take the DIM function of a taylor and a default kind REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION dim_sreal_taylor(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other
      REAL(kind=sr_kind), INTENT(IN) :: self
      TYPE(taylor) :: dim_sreal_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(other)).GT.Real_args_tol) THEN
        dim_sreal_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      dim_sreal_taylor = MAX(self-other,0.)

   END FUNCTION dim_sreal_taylor

!!$****************************************************************
!!$ Take the complex conjugate of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION conjg_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: conjg_taylor

!!$ Code starts here

      conjg_taylor%drv(:,:) = CONJG(self%drv(:,:))

   END FUNCTION conjg_taylor

!!$****************************************************************
!!$ Take the remainder of a taylor modulo a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION mod_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor) :: mod_taylor_real

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        mod_taylor_real%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      mod_taylor_real%drv(:,0) = MOD(REAL(self%drv(:,0),dr_kind),other)
      mod_taylor_real%drv(:,1:) = self%drv(:,1:)
      IF(self==other) mod_taylor_real%drv(:,1:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)

   END FUNCTION mod_taylor_real

!!$****************************************************************
!!$ Take the remainder of a taylor modulo a REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION modulo_taylor_real(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dr_kind), INTENT(IN) :: other
      TYPE(taylor) :: modulo_taylor_real

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        modulo_taylor_real%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      modulo_taylor_real%drv(:,0) = MODULO(REAL(self%drv(:,0),dr_kind),other)
      modulo_taylor_real%drv(:,1:) = self%drv(:,1:)
      IF(self==other) modulo_taylor_real%drv(:,1:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)

   END FUNCTION modulo_taylor_real

!!$****************************************************************
!!$ Take the remainder of a taylor modulo a default REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION mod_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor) :: mod_taylor_sreal

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        mod_taylor_sreal%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      mod_taylor_sreal%drv(:,0) = MOD(REAL(self%drv(:,0),sr_kind),other)
      mod_taylor_sreal%drv(:,1:) = self%drv(:,1:)
      IF(self==other) mod_taylor_sreal%drv(:,1:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)

   END FUNCTION mod_taylor_sreal

!!$****************************************************************
!!$ Take the remainder of a taylor modulo a default REAL
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION modulo_taylor_sreal(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=sr_kind), INTENT(IN) :: other
      TYPE(taylor) :: modulo_taylor_sreal

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        modulo_taylor_sreal%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      modulo_taylor_sreal%drv(:,0) = MODULO(REAL(self%drv(:,0),sr_kind),other)
      modulo_taylor_sreal%drv(:,1:) = self%drv(:,1:)
      IF(self==other) modulo_taylor_sreal%drv(:,1:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)

   END FUNCTION modulo_taylor_sreal

!!$****************************************************************
!!$ Take the integer part of a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION aint_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: aint_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        aint_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      aint_taylor%drv(:,0) = AINT(REAL(self%drv(:,0)))
      aint_taylor%drv(:,1:) = 0
      IF(self==aint_taylor) aint_taylor%drv(:,1:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)

   END FUNCTION aint_taylor

!!$****************************************************************
!!$ Take the rounded value of a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION anint_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      TYPE(taylor) :: anint_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        anint_taylor%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      anint_taylor%drv(:,0) = ANINT(REAL(self%drv(:,0)))
      anint_taylor%drv(:,1:) = 0
      IF(self==anint_taylor-0.5) anint_taylor%drv(:,1:) = IEEE_VALUE(0.,IEEE_QUIET_NAN)

   END FUNCTION anint_taylor

!!$****************************************************************
!!$ Take the ceiling of a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION ceiling_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER :: ceiling_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        ceiling_taylor = -HUGE(0)
        RETURN
      ENDIF

      ceiling_taylor = CEILING(REAL(self%drv(1,0)))

   END FUNCTION ceiling_taylor

!!$****************************************************************
!!$ Take the floor of a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION floor_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER :: floor_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        floor_taylor = HUGE(0)
        RETURN
      ENDIF

      floor_taylor = FLOOR(REAL(self%drv(1,0)))

   END FUNCTION floor_taylor

!!$****************************************************************
!!$ Take the integer truncation of a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION int_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER :: int_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        int_taylor = HUGE(0)
        RETURN
      ENDIF

      int_taylor = INT(REAL(self%drv(1,0)))

   END FUNCTION int_taylor

!!$****************************************************************
!!$ Take the nearest integer to a taylor
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   ELEMENTAL FUNCTION nint_taylor(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER :: nint_taylor

!!$ Code starts here

      IF(Real_args_warn.AND.ABS(AIMAG(self)).GT.Real_args_tol) THEN
        nint_taylor = HUGE(0)
        RETURN
      ENDIF

      nint_taylor = NINT(REAL(self%drv(1,0)))

   END FUNCTION nint_taylor

!!$****************************************************************
!!$ Maximum of an array of taylors
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   PURE FUNCTION maxval_taylors(self,mask)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      LOGICAL, INTENT(IN),OPTIONAL :: mask(:)
      TYPE(taylor) :: maxval_taylors

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      IF(Real_args_warn.AND.ANY(ABS(AIMAG(self)).GT.Real_args_tol)) THEN
        maxval_taylors%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      maxval_taylors = -huge(0._dp)
      IF(.NOT.PRESENT(mask)) THEN
         DO i=1,SIZE(self)
            maxval_taylors = max(maxval_taylors,self(i))
         ENDDO
      ELSE
         DO i=1,SIZE(self)
            IF(mask(i)) maxval_taylors = max(maxval_taylors,self(i))
         ENDDO
      ENDIF

   END FUNCTION maxval_taylors

!!$****************************************************************
!!$ Minimum of an array of taylors
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   PURE FUNCTION minval_taylors(self,mask)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      LOGICAL, INTENT(IN),OPTIONAL :: mask(:)
      TYPE(taylor) :: minval_taylors

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      IF(Real_args_warn.AND.ANY(ABS(AIMAG(self)).GT.Real_args_tol)) THEN
        minval_taylors%drv(:,:) = IEEE_VALUE(0._dp,IEEE_QUIET_NAN)
        RETURN
      ENDIF

      minval_taylors = huge(0._dp)
      IF(.NOT.PRESENT(mask)) THEN
         DO i=1,SIZE(self)
            minval_taylors = min(minval_taylors,self(i))
         ENDDO
      ELSE
         DO i=1,SIZE(self)
            IF(mask(i)) minval_taylors = min(minval_taylors,self(i))
         ENDDO
      ENDIF

   END FUNCTION minval_taylors

!!$****************************************************************
!!$ Location of maximum of an array of taylors
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   PURE FUNCTION maxloc_taylors(self,mask)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      LOGICAL, INTENT(IN),OPTIONAL :: mask(:)
      INTEGER :: maxloc_taylors(1)

!!$ Code starts here

      IF(Real_args_warn.AND.ANY(ABS(AIMAG(self)).GT.Real_args_tol)) THEN
        maxloc_taylors = -HUGE(0)
        RETURN
      ENDIF

      IF(.NOT.PRESENT(mask)) THEN
         maxloc_taylors = maxloc(REAL(self%drv(1,0),dr_kind))
      ELSE
         maxloc_taylors = maxloc(REAL(self%drv(1,0),dr_kind),mask)
      ENDIF

   END FUNCTION maxloc_taylors

!!$****************************************************************
!!$ Location of minimum of an array of taylors
!!$****************************************************************

!!$ WARNING -- implicitly takes real part of argument

   PURE FUNCTION minloc_taylors(self,mask)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      LOGICAL, INTENT(IN),OPTIONAL :: mask(:)
      INTEGER :: minloc_taylors(1)

!!$ Code starts here

      IF(Real_args_warn.AND.ANY(ABS(AIMAG(self)).GT.Real_args_tol)) THEN
        minloc_taylors = -HUGE(0)
        RETURN
      ENDIF

      IF(.NOT.PRESENT(mask)) THEN
         minloc_taylors = minloc(REAL(self%drv(1,0),dr_kind))
      ELSE
         minloc_taylors = minloc(REAL(self%drv(1,0),dr_kind),mask)
      ENDIF

   END FUNCTION minloc_taylors

!!$****************************************************************
!!$ Sum of an array of taylors
!!$****************************************************************

   PURE FUNCTION sum_taylors(self,mask)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      LOGICAL, INTENT(IN),OPTIONAL :: mask(:)
      TYPE(taylor) :: sum_taylors

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      sum_taylors = 0
      IF(.NOT.PRESENT(mask)) THEN
         DO i=1,SIZE(self)
            sum_taylors = sum_taylors + self(i)
         ENDDO
      ELSE
         DO i=1,SIZE(self)
            IF(mask(i)) sum_taylors = sum_taylors + self(i)
         ENDDO
      ENDIF

   END FUNCTION sum_taylors

!!$****************************************************************
!!$ Product of an array of taylors
!!$****************************************************************

   PURE FUNCTION product_taylors(self,masked)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      LOGICAL, INTENT(IN),OPTIONAL :: masked(:)
      TYPE(taylor) :: product_taylors

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      product_taylors = 1
      IF(.NOT.PRESENT(masked)) THEN
         DO i=1,SIZE(self)
            product_taylors = product_taylors * self(i)
         ENDDO
      ELSE
         DO i=1,SIZE(self)
            IF(masked(i)) product_taylors = product_taylors * self(i)
         ENDDO
      ENDIF

   END FUNCTION product_taylors

!!$****************************************************************
!!$ Dot product of two vectors of taylors
!!$****************************************************************

   PURE FUNCTION dot_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:),other(:)
      TYPE(taylor) :: dot_taylors

!!$ Code starts here

      dot_taylors = SUM(CONJG(self)*other)

   END FUNCTION dot_taylors

!!$****************************************************************
!!$ Dot product of a vector of taylors and one of REALs
!!$****************************************************************

   PURE FUNCTION dot_taylors_reals(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      REAL(kind=dr_kind),INTENT(IN) :: other(:)
      TYPE(taylor) :: dot_taylors_reals

!!$ Code starts here

      dot_taylors_reals = SUM(self*other)

   END FUNCTION dot_taylors_reals

!!$****************************************************************
!!$ Dot product of a vector of taylors and one of REALs
!!$****************************************************************

   PURE FUNCTION dot_reals_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:)
      REAL(kind=dr_kind),INTENT(IN) :: self(:)
      TYPE(taylor) :: dot_reals_taylors

!!$ Code starts here

      dot_reals_taylors = SUM(self*other)

   END FUNCTION dot_reals_taylors

!!$****************************************************************
!!$ Dot product of a vector of taylors and one of default REALs
!!$****************************************************************

   PURE FUNCTION dot_taylors_sreals(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      REAL(kind=sr_kind),INTENT(IN) :: other(:)
      TYPE(taylor) :: dot_taylors_sreals

!!$ Code starts here

      dot_taylors_sreals = SUM(self*other)

   END FUNCTION dot_taylors_sreals

!!$****************************************************************
!!$ Dot product of a vector of taylors and one of default REALs
!!$****************************************************************

   PURE FUNCTION dot_sreals_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:)
      REAL(kind=sr_kind),INTENT(IN) :: self(:)
      TYPE(taylor) :: dot_sreals_taylors

!!$ Code starts here

      dot_sreals_taylors = SUM(self*other)

   END FUNCTION dot_sreals_taylors

!!$****************************************************************
!!$ Dot product of a vector of taylors and one of COMPLEXes
!!$****************************************************************

   PURE FUNCTION dot_taylors_complexes(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      COMPLEX(kind=dc_kind),INTENT(IN) :: other(:)
      TYPE(taylor) :: dot_taylors_complexes

!!$ Code starts here

      dot_taylors_complexes = SUM(CONJG(self)*other)

   END FUNCTION dot_taylors_complexes

!!$****************************************************************
!!$ Dot product of a vector of taylors and one of COMPLEXes
!!$****************************************************************

   PURE FUNCTION dot_complexes_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:)
      COMPLEX(kind=dc_kind),INTENT(IN) :: self(:)
      TYPE(taylor) :: dot_complexes_taylors

!!$ Code starts here

      dot_complexes_taylors = SUM(CONJG(self)*other)

   END FUNCTION dot_complexes_taylors

!!$****************************************************************
!!$ Dot product of a vector of taylors and one of default COMPLEXes
!!$****************************************************************

   PURE FUNCTION dot_taylors_scomplexes(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:)
      COMPLEX(kind=sc_kind),INTENT(IN) :: other(:)
      TYPE(taylor) :: dot_taylors_scomplexes

!!$ Code starts here

      dot_taylors_scomplexes = SUM(CONJG(self)*other)

   END FUNCTION dot_taylors_scomplexes

!!$****************************************************************
!!$ Dot product of a vector of taylors and one of default COMPLEXes
!!$****************************************************************

   PURE FUNCTION dot_scomplexes_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:)
      COMPLEX(kind=sc_kind),INTENT(IN) :: self(:)
      TYPE(taylor) :: dot_scomplexes_taylors

!!$ Code starts here

      dot_scomplexes_taylors = SUM(CONJG(self)*other)

   END FUNCTION dot_scomplexes_taylors

!!$****************************************************************
!!$ Matrix product of two taylor matrices
!!$****************************************************************

   PURE FUNCTION matmul_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:), other(:,:)
      TYPE(taylor) :: matmul_taylors(size(self,1),size(other,2))

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      FORALL(i=1:size(self,1), j=1:size(self,2)) &
           matmul_taylors(i,j) = SUM(self(i,:)*other(:,j))

   END FUNCTION matmul_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION matvec_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:), other(:)
      TYPE(taylor) :: matvec_taylors(size(self,1))

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      FORALL(i=1:size(self,1)) &
           matvec_taylors(i) = SUM(self(i,:)*other(:))

   END FUNCTION matvec_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION vecmat_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:), other(:,:)
      TYPE(taylor) :: vecmat_taylors(size(other,2))

!!$ Internal variables

      INTEGER :: j

!!$ Code starts here

      FORALL(j=1:size(other,2)) &
           vecmat_taylors(j) = SUM(self(:)*other(:,j))

   END FUNCTION vecmat_taylors

!!$****************************************************************
!!$ Matrix product of a taylor matrix and a REAL matrix
!!$****************************************************************

   PURE FUNCTION matmul_taylors_reals(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:)
      REAL(kind=dr_kind), INTENT(IN) :: other(:,:)
      TYPE(taylor) :: matmul_taylors_reals(size(self,1),size(other,2))

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      FORALL(i=1:size(self,1), j=1:size(self,2)) &
           matmul_taylors_reals(i,j) = SUM(self(i,:)*other(:,j))

   END FUNCTION matmul_taylors_reals

!!$****************************************************************
!!$ Matrix product of a taylor matrix and a REAL matrix
!!$****************************************************************

   PURE FUNCTION matmul_reals_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:,:)
      REAL(kind=dr_kind), INTENT(IN) :: self(:,:)
      TYPE(taylor) :: matmul_reals_taylors(size(self,1),size(other,2))

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      FORALL(i=1:size(self,1), j=1:size(self,2)) &
           matmul_reals_taylors(i,j) = SUM(self(i,:)*other(:,j))

   END FUNCTION matmul_reals_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and REAL vector
!!$****************************************************************

   PURE FUNCTION matvec_taylors_reals(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:)
      REAL(kind=dr_kind), INTENT(IN) :: other(:)
      TYPE(taylor) :: matvec_taylors_reals(size(self,1))

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      FORALL(i=1:size(self,1)) &
           matvec_taylors_reals(i) = SUM(self(i,:)*other(:))

   END FUNCTION matvec_taylors_reals

!!$****************************************************************
!!$ Matrix product of REAL matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION matvec_reals_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:)
      REAL(kind=dr_kind), INTENT(IN) :: self(:,:)
      TYPE(taylor) :: matvec_reals_taylors(size(self,1))

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      FORALL(i=1:size(self,1)) &
           matvec_reals_taylors(i) = SUM(self(i,:)*other(:))

   END FUNCTION matvec_reals_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and REAL vector
!!$****************************************************************

   PURE FUNCTION vecmat_reals_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: self(:)
      TYPE(taylor), INTENT(IN) :: other(:,:)
      TYPE(taylor) :: vecmat_reals_taylors(size(other,2))

!!$ Internal variables

      INTEGER :: j

!!$ Code starts here

      FORALL(j=1:size(other,2)) &
           vecmat_reals_taylors(j) = SUM(self(:)*other(:,j))

   END FUNCTION vecmat_reals_taylors

!!$****************************************************************
!!$ Matrix product of REAL matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION vecmat_taylors_reals(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=dr_kind), INTENT(IN) :: other(:,:)
      TYPE(taylor), INTENT(IN) :: self(:)
      TYPE(taylor) :: vecmat_taylors_reals(size(other,2))

!!$ Internal variables

      INTEGER :: j

!!$ Code starts here

      FORALL(j=1:size(other,2)) &
           vecmat_taylors_reals(j) = SUM(self(:)*other(:,j))

   END FUNCTION vecmat_taylors_reals

!!$****************************************************************
!!$ Matrix product of a taylor matrix and a default REAL matrix
!!$****************************************************************

   PURE FUNCTION matmul_taylors_sreals(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:)
      REAL(kind=sr_kind), INTENT(IN) :: other(:,:)
      TYPE(taylor) :: matmul_taylors_sreals(size(self,1),size(other,2))

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      FORALL(i=1:size(self,1), j=1:size(self,2)) &
           matmul_taylors_sreals(i,j) = SUM(self(i,:)*other(:,j))

   END FUNCTION matmul_taylors_sreals

!!$****************************************************************
!!$ Matrix product of a taylor matrix and a default REAL matrix
!!$****************************************************************

   PURE FUNCTION matmul_sreals_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:,:)
      REAL(kind=sr_kind), INTENT(IN) :: self(:,:)
      TYPE(taylor) :: matmul_sreals_taylors(size(self,1),size(other,2))

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      FORALL(i=1:size(self,1), j=1:size(self,2)) &
           matmul_sreals_taylors(i,j) = SUM(self(i,:)*other(:,j))

   END FUNCTION matmul_sreals_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and default REAL vector
!!$****************************************************************

   PURE FUNCTION matvec_taylors_sreals(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:)
      REAL(kind=sr_kind), INTENT(IN) :: other(:)
      TYPE(taylor) :: matvec_taylors_sreals(size(self,1))

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      FORALL(i=1:size(self,1)) &
           matvec_taylors_sreals(i) = SUM(self(i,:)*other(:))

   END FUNCTION matvec_taylors_sreals

!!$****************************************************************
!!$ Matrix product of default REAL matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION matvec_sreals_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:)
      REAL(kind=sr_kind), INTENT(IN) :: self(:,:)
      TYPE(taylor) :: matvec_sreals_taylors(size(self,1))

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      FORALL(i=1:size(self,1)) &
           matvec_sreals_taylors(i) = SUM(self(i,:)*other(:))

   END FUNCTION matvec_sreals_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and default REAL vector
!!$****************************************************************

   PURE FUNCTION vecmat_sreals_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: self(:)
      TYPE(taylor), INTENT(IN) :: other(:,:)
      TYPE(taylor) :: vecmat_sreals_taylors(size(other,2))

!!$ Internal variables

      INTEGER :: j

!!$ Code starts here

      FORALL(j=1:size(other,2)) &
           vecmat_sreals_taylors(j) = SUM(self(:)*other(:,j))

   END FUNCTION vecmat_sreals_taylors

!!$****************************************************************
!!$ Matrix product of default REAL matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION vecmat_taylors_sreals(self,other)

      IMPLICIT NONE

!!$ Passed variables

      REAL(kind=sr_kind), INTENT(IN) :: other(:,:)
      TYPE(taylor), INTENT(IN) :: self(:)
      TYPE(taylor) :: vecmat_taylors_sreals(size(other,2))

!!$ Internal variables

      INTEGER :: j

!!$ Code starts here

      FORALL(j=1:size(other,2)) &
           vecmat_taylors_sreals(j) = SUM(self(:)*other(:,j))

   END FUNCTION vecmat_taylors_sreals

!!$****************************************************************
!!$ Matrix product of a taylor matrix and a COMPLEX matrix
!!$****************************************************************

   PURE FUNCTION matmul_taylors_complexes(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:)
      COMPLEX(kind=dc_kind), INTENT(IN) :: other(:,:)
      TYPE(taylor) :: matmul_taylors_complexes(size(self,1),size(other,2))

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      FORALL(i=1:size(self,1), j=1:size(self,2)) &
           matmul_taylors_complexes(i,j) = SUM(self(i,:)*other(:,j))

   END FUNCTION matmul_taylors_complexes

!!$****************************************************************
!!$ Matrix product of a taylor matrix and a COMPLEX matrix
!!$****************************************************************

   PURE FUNCTION matmul_complexes_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:,:)
      COMPLEX(kind=dc_kind), INTENT(IN) :: self(:,:)
      TYPE(taylor) :: matmul_complexes_taylors(size(self,1),size(other,2))

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      FORALL(i=1:size(self,1), j=1:size(self,2)) &
           matmul_complexes_taylors(i,j) = SUM(self(i,:)*other(:,j))

   END FUNCTION matmul_complexes_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and COMPLEX vector
!!$****************************************************************

   PURE FUNCTION matvec_taylors_complexes(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:)
      COMPLEX(kind=dc_kind), INTENT(IN) :: other(:)
      TYPE(taylor) :: matvec_taylors_complexes(size(self,1))

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      FORALL(i=1:size(self,1)) &
           matvec_taylors_complexes(i) = SUM(self(i,:)*other(:))

   END FUNCTION matvec_taylors_complexes

!!$****************************************************************
!!$ Matrix product of COMPLEX matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION matvec_complexes_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:)
      COMPLEX(kind=dc_kind), INTENT(IN) :: self(:,:)
      TYPE(taylor) :: matvec_complexes_taylors(size(self,1))

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      FORALL(i=1:size(self,1)) &
           matvec_complexes_taylors(i) = SUM(self(i,:)*other(:))

   END FUNCTION matvec_complexes_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and COMPLEX vector
!!$****************************************************************

   PURE FUNCTION vecmat_complexes_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN) :: self(:)
      TYPE(taylor), INTENT(IN) :: other(:,:)
      TYPE(taylor) :: vecmat_complexes_taylors(size(other,2))

!!$ Internal variables

      INTEGER :: j

!!$ Code starts here

      FORALL(j=1:size(other,2)) &
           vecmat_complexes_taylors(j) = SUM(self(:)*other(:,j))

   END FUNCTION vecmat_complexes_taylors

!!$****************************************************************
!!$ Matrix product of COMPLEX matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION vecmat_taylors_complexes(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=dc_kind), INTENT(IN) :: other(:,:)
      TYPE(taylor), INTENT(IN) :: self(:)
      TYPE(taylor) :: vecmat_taylors_complexes(size(other,2))

!!$ Internal variables

      INTEGER :: j

!!$ Code starts here

      FORALL(j=1:size(other,2)) &
           vecmat_taylors_complexes(j) = SUM(self(:)*other(:,j))

   END FUNCTION vecmat_taylors_complexes

!!$****************************************************************
!!$ Matrix product of a taylor matrix and a default COMPLEX matrix
!!$****************************************************************

   PURE FUNCTION matmul_taylors_scomplexes(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:)
      COMPLEX(kind=sc_kind), INTENT(IN) :: other(:,:)
      TYPE(taylor) :: matmul_taylors_scomplexes(size(self,1),size(other,2))

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      FORALL(i=1:size(self,1), j=1:size(self,2)) &
           matmul_taylors_scomplexes(i,j) = SUM(self(i,:)*other(:,j))

   END FUNCTION matmul_taylors_scomplexes

!!$****************************************************************
!!$ Matrix product of a taylor matrix and a default COMPLEX matrix
!!$****************************************************************

   PURE FUNCTION matmul_scomplexes_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:,:)
      COMPLEX(kind=sc_kind), INTENT(IN) :: self(:,:)
      TYPE(taylor) :: matmul_scomplexes_taylors(size(self,1),size(other,2))

!!$ Internal variables

      INTEGER :: i,j

!!$ Code starts here

      FORALL(i=1:size(self,1), j=1:size(self,2)) &
           matmul_scomplexes_taylors(i,j) = SUM(self(i,:)*other(:,j))

   END FUNCTION matmul_scomplexes_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and default COMPLEX vector
!!$****************************************************************

   PURE FUNCTION matvec_taylors_scomplexes(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self(:,:)
      COMPLEX(kind=sc_kind), INTENT(IN) :: other(:)
      TYPE(taylor) :: matvec_taylors_scomplexes(size(self,1))

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      FORALL(i=1:size(self,1)) &
           matvec_taylors_scomplexes(i) = SUM(self(i,:)*other(:))

   END FUNCTION matvec_taylors_scomplexes

!!$****************************************************************
!!$ Matrix product of default COMPLEX matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION matvec_scomplexes_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: other(:)
      COMPLEX(kind=sc_kind), INTENT(IN) :: self(:,:)
      TYPE(taylor) :: matvec_scomplexes_taylors(size(self,1))

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      FORALL(i=1:size(self,1)) &
           matvec_scomplexes_taylors(i) = SUM(self(i,:)*other(:))

   END FUNCTION matvec_scomplexes_taylors

!!$****************************************************************
!!$ Matrix product of taylor matrix and default COMPLEX vector
!!$****************************************************************

   PURE FUNCTION vecmat_scomplexes_taylors(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=sc_kind), INTENT(IN) :: self(:)
      TYPE(taylor), INTENT(IN) :: other(:,:)
      TYPE(taylor) :: vecmat_scomplexes_taylors(size(other,2))

!!$ Internal variables

      INTEGER :: j

!!$ Code starts here

      FORALL(j=1:size(other,2)) &
           vecmat_scomplexes_taylors(j) = SUM(self(:)*other(:,j))

   END FUNCTION vecmat_scomplexes_taylors

!!$****************************************************************
!!$ Matrix product of default COMPLEX matrix and taylor vector
!!$****************************************************************

   PURE FUNCTION vecmat_taylors_scomplexes(self,other)

      IMPLICIT NONE

!!$ Passed variables

      COMPLEX(kind=sc_kind), INTENT(IN) :: other(:,:)
      TYPE(taylor), INTENT(IN) :: self(:)
      TYPE(taylor) :: vecmat_taylors_scomplexes(size(other,2))

!!$ Internal variables

      INTEGER :: j

!!$ Code starts here

      FORALL(j=1:size(other,2)) &
           vecmat_taylors_scomplexes(j) = SUM(self(:)*other(:,j))

   END FUNCTION vecmat_taylors_scomplexes


!!$****************************************************************
!!$ USER FUNCTIONS
!!$****************************************************************

!!$****************************************************************
!!$ Return an independent variable with COMPLEX value
!!$****************************************************************

   FUNCTION independent_complex(index,value)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: index
      COMPLEX(kind=dc_kind), INTENT(IN) :: value
      TYPE(taylor) :: independent_complex

!!$ Code starts here

      CALL fdb_generate(Taylor_order)

      independent_complex%drv(:,:) = 0._dp
      independent_complex%drv(:,0) = value
      independent_complex%drv(index,1) = 1._dp

   END FUNCTION independent_complex

!!$****************************************************************
!!$ Return an independent variable with default kind COMPLEX value
!!$****************************************************************

   FUNCTION independent_scomplex(index,value)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: index
      COMPLEX(kind=sc_kind), INTENT(IN) :: value
      TYPE(taylor) :: independent_scomplex

!!$ Code starts here

      CALL fdb_generate(Taylor_order)

      independent_scomplex%drv(:,:) = 0._dp
      independent_scomplex%drv(:,0) = value
      independent_scomplex%drv(index,1) = 1._dp

   END FUNCTION independent_scomplex

!!$****************************************************************
!!$ Return an independent variable with REAL value
!!$****************************************************************

   FUNCTION independent_real(index,value)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: index
      REAL(kind=dr_kind), INTENT(IN) :: value
      TYPE(taylor) :: independent_real

!!$ Code starts here

      CALL fdb_generate(Taylor_order)

      independent_real%drv(:,:) = 0._dp
      independent_real%drv(:,0) = value
      independent_real%drv(index,1) = 1._dp

   END FUNCTION independent_real

!!$****************************************************************
!!$ Return an independent variable with default kind REAL value
!!$****************************************************************

   FUNCTION independent_sreal(index,value)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: index
      REAL(kind=sr_kind), INTENT(IN) :: value
      TYPE(taylor) :: independent_sreal

!!$ Code starts here

      CALL fdb_generate(Taylor_order)

      independent_sreal%drv(:,:) = 0._dp
      independent_sreal%drv(:,0) = value
      independent_sreal%drv(index,1) = 1._dp

   END FUNCTION independent_sreal

!!$****************************************************************
!!$ Return an independent variable with INTEGER value
!!$****************************************************************

   FUNCTION independent_int(index,value)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(IN) :: index
      INTEGER, INTENT(IN) :: value
      TYPE(taylor) :: independent_int

!!$ Code starts here

      CALL fdb_generate(Taylor_order)

      independent_int%drv(:,:) = 0._dp
      independent_int%drv(:,0) = value
      independent_int%drv(index,1) = 1._dp

   END FUNCTION independent_int

!!$****************************************************************
!!$ Return the value of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION value(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=dc_kind) :: value

!!$ Code starts here

      value = self%drv(1,0)

   END FUNCTION value

!!$****************************************************************
!!$ Return the real part of the value of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION realvalue(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dc_kind) :: realvalue

!!$ Code starts here

      realvalue = REAL(self%drv(1,0))

   END FUNCTION realvalue

!!$****************************************************************
!!$ Return the imaginary part of the value of a taylor
!!$****************************************************************

   ELEMENTAL FUNCTION imagvalue(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      REAL(kind=dc_kind) :: imagvalue

!!$ Code starts here

      imagvalue = AIMAG(self%drv(1,0))

   END FUNCTION imagvalue

!!$****************************************************************
!!$ Return the vector of first derivatives
!!$****************************************************************

   PURE FUNCTION gradient(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=dc_kind) :: gradient(1:N_taylor_vars)

!!$ Code starts here

      gradient = self%drv(:,1)

   END FUNCTION gradient

!!$****************************************************************
!!$ Return the Laplacian
!!$****************************************************************

   ELEMENTAL FUNCTION laplacian(self)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      COMPLEX(kind=dc_kind) :: laplacian

!!$ Code starts here

      laplacian = SUM(self%drv(:,2))

   END FUNCTION laplacian

!!$****************************************************************
!!$ Return a specific derivative
!!$****************************************************************

   ELEMENTAL FUNCTION derivative(self,var,order)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: var
      INTEGER, INTENT(IN), OPTIONAL :: order
      COMPLEX(kind=dc_kind) :: derivative

!!$ Code starts here

      IF(PRESENT(order)) THEN
         derivative = self%drv(var,order)
      ELSE
         derivative = self%drv(var,1)
      ENDIF

   END FUNCTION derivative

!!$****************************************************************
!!$ Return the expansion in one variable
!!$****************************************************************

   PURE FUNCTION expansion(self,var)

      IMPLICIT NONE

!!$ Passed variables

      TYPE(taylor), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: var
      COMPLEX(kind=dc_kind), POINTER :: expansion(:)

!!$ Code starts here

      ALLOCATE(expansion(0:Taylor_order))

      expansion(:) = self%drv(var,0:Taylor_order)

   END FUNCTION expansion


!!$****************************************************************
!!$ AUXILIARY FUNCTIONS
!!$****************************************************************

!!$****************************************************************
!!$ The binomial coefficient (n over k)
!!$****************************************************************
!!$ choose(n,k) = n!/(k!(n-k)!)
!!$****************************************************************

   PURE FUNCTION choose(n,k)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(in) :: n,k
      INTEGER :: choose

!!$ Internal variables

      INTEGER :: i,perm

!!$ Code starts here

      choose = 0 ! Exclude pathological values of argument
      IF((n.LT.k).OR.(n.LT.0).OR.(k.LT.0)) RETURN
     
      choose = 1
      perm = 1

      DO i = 1,k
       choose = choose*(n-i+1)
       perm = perm*i
      ENDDO
      choose = choose/perm

   END FUNCTION choose

!!$****************************************************************
!!$ The factorial of n
!!$****************************************************************
!!$ fac(n,k) = n!
!!$****************************************************************

   PURE FUNCTION fac(n)

      IMPLICIT NONE

!!$ Passed variables

      INTEGER, INTENT(in) :: n
      INTEGER :: fac

!!$ Internal variables

      INTEGER :: i

!!$ Code starts here

      fac = 0 ! Exclude pathological values of argument
      IF(n.LT.0) RETURN
     
      fac = 1

      DO i = 1,n
       fac = fac*(n-i+1)
      ENDDO

   END FUNCTION fac

!!$****************************************************************
!!$ Compute numbers in Faa di Bruno's theorem
!!$****************************************************************
!!$ idx(:,1) are the derivative orders of the external derivatives,
!!$ idx(:,2:p+1) are the powers of the internal derivatives,
!!$ wgt(:) are the weights of the terms, i.e.:
!!$
!!$ [F(y(x))]^(p) = sum_j wgt(j) * F^(idx(j,1))(y(x)) * \
!!$                 (y^(1)(x))^idx(j,2) *...* (y^(p)(x))^idx(j,p+1)
!!$
!!$****************************************************************

   SUBROUTINE fdb_generate(maxp)

!!$ Passed variables

      INTEGER, INTENT(IN) :: maxp

!!$ Internal variables

      INTEGER :: k(1+Max_taylor_order),t(0:1+Max_taylor_order),m,i,nterms,p
      REAL :: factor(0:Max_taylor_order)

!!$ Code starts here

      nterms = have(maxp)
      
      IF(nterms.NE.0) GOTO 100

      FORALL(i=0:Max_taylor_order) factor(i) = fac(i)

      DO p=1,maxp
       nterms = 0
       t(0) = NINT(p/2.)
       t(1) = -1
       t(2:) = 0
      
       loop: DO
         t(1) = t(1) + 1
         DO i=1,p-1 
            IF(t(i).GT.MIN(t(i-1),p/2)) THEN
               t(1:i) = t(i+1) + 1 
               t(i+1) = t(i+1) + 1
            ENDIF
         ENDDO

         k(p) = t(p-1)
         k(2:p-1) = t(1:p-2) - t(2:p-1)
         k(1) = p - 2*t(1) - SUM(t(2:p-1))

         m = SUM(k(1:p))
         IF(ANY(k(1:p)<0)) CYCLE loop

         nterms = nterms + 1
         scratchi(nterms,1:p+1,p) = (/m,k(1:p)/)
         scratchr(nterms,p) = PRODUCT(1./factor(1:p)**k(1:p))*factor(p)/&
                               PRODUCT(factor(k(1:p))) 

         IF(k(p)==1) EXIT loop

       ENDDO loop

       have(p) = nterms

      ENDDO 

100   RETURN 

   END SUBROUTINE fdb_generate

!!$****************************************************************
!!$ Return precomputed numbers in Faa di Bruno's theorem
!!$****************************************************************

   PURE SUBROUTINE fdb(p,idx,wgt)
     
!!$ Passed variables
     
     INTEGER, INTENT(IN) :: p
     INTEGER, POINTER :: idx(:,:)
     INTEGER, POINTER :: wgt(:)
     
!!$ Code starts here
     
     ALLOCATE(idx(have(p),1+p))
     ALLOCATE(wgt(have(p)))
     
     idx(:,:) = scratchi(1:have(p),1:p+1,p)
     wgt(:) = scratchr(1:have(p),p)
     
   END SUBROUTINE fdb
   
 END MODULE cens_lib_ad
 
