MODULE MYDEFS
! Definition Module Src ESS 03/011/2010
  integer, parameter    :: dbl = selected_real_kind(15,307)
  integer, parameter    :: large = selected_int_kind(8)
  real (dbl), parameter :: PI = 3.14159265358979323_dbl
  real(dbl), parameter :: SQRT3=1.73205080756887729352744634151D0
  real (dbl), parameter :: TWOPI = 2.0_dbl*PI
  real (dbl), parameter :: ZERO = 0.0_dbl
  real (dbl), parameter :: ONE = 1.0_dbl
  complex(dbl), parameter :: ii = (zero,one)
END MODULE MYDEFS
