program main

  use precision
  use solver

  implicit none
  
  real(kind=kd), dimension(:,:), allocatable :: u, v, p
  integer :: nt, nx, ny, nit
  real(kind=kd) :: c, dt, dx, dy, rho, nu
  real(kind=kd), dimension(:), allocatable :: x, y

  integer :: i

  nx = 41; ny = 41
  nt = 1000; nit = 100
  c = 1

  dx = 2.0 / (nx-1)
  dy = 2.0 / (ny-1)

  allocate(x(nx))
  allocate(y(ny))
  
  x(:) = (/ (i, i=0,2,nx) /)
  y(:) = (/ (i, i=0,2,ny) /)

  rho = 1
  nu = 0.1
  dt = 0.001

  allocate(u(nx, ny))
  allocate(v(nx, ny))
  allocate(p(nx, ny))
  
  u(:,:) = 0
  v(:,:) = 0
  p(:,:) = 0

  call cavityFlow(nt, nit, u, v, dt, dx, dy, p, rho, nu)
  
end program main
