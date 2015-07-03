module precision
  !! Provides kind attributes for setting machine-compiler-independent precision for real numbers. Here, we define "single precision" to mean 32-bit precision, "double precision" to mean 64-bit precision, and "quadruple precision" to mean 128-bit precision. This ensures that precision is preserved regardless of the compiler or computer architecture.

  use iso_fortran_env
  use iso_c_binding

  implicit none
  
  integer, parameter :: sp = real32 !selected_real_kind(6, 37)
  !! Single precision: 32-bit real
  integer, parameter :: dp = real64 !selected_real_kind(15, 307)
  !! Double precision: 64-bit real
  integer, parameter :: qp = real128 !selected_real_kind(33, 4931)
  !! Quadruple precision: 128-bit real

  integer, parameter :: c_sp = c_float
  !! C type for single precision 
  integer, parameter :: c_dp = c_double
  !! C type for double precision 
  integer, parameter :: c_qp = c_long_double
  !! C type for quadruple precision

  integer, parameter :: kd = dp
  !! The precision of real types used in this program
  integer, parameter :: c_kd = c_dp
  !! The precision of C real types used in this program
  
contains 
  
end module precision
  
