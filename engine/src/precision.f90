module precision
  !! Provides kind attributes for setting machine-compiler-independent precision for real numbers. Here, we define "single precision" to mean 32-bit precision, "double precision" to mean 64-bit precision, and "quadruple precision" to mean 128-bit precision. This ensures that precision is preserved regardless of the compiler or computer architecture.

  implicit none
  
  integer, parameter :: sp = selected_real_kind(6, 37)
  !! Single precision: 32-bit real
  integer, parameter :: dp = selected_real_kind(15, 307)
  !! Double precision: 64-bit real
  integer, parameter :: qp = selected_real_kind(33, 4931)
  !! Quadruple precision: 128-bit real

  integer, parameter :: kd = dp
  !! The precision of real types used in this program
  
contains 
  
end module precision
  
