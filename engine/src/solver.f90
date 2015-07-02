module solver
  !! Contains the discretized 2D Navier-Stokes solver. The x-momentum equation is given by:
  !! $$ \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} = -\frac{1}{\rho} \frac{\partial p}{\partial x} + \nu \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) $$
  !! The y-momentum equation is similarly given by:
  !! $$ \frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} = -\frac{1}{\rho} \frac{\partial p}{\partial y} + \nu \left( \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right) $$
  !! Finally, the pressure equation is given by:
  !! $$ \frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2} = -\rho \left( \frac{\partial u}{\partial x} \frac{\partial u}{\partial x} + 2 \frac{\partial u}{\partial y} \frac{\partial v}{\partial x} + \frac{\partial v}{\partial y} \frac{\partial v}{\partial y} \right) $$

  use precision

  public :: cavityFlow
  private :: ppe, ppb
  
contains

  subroutine cavityFlow(nt, nit, u, v, dt, dx, dy, p, rho, nu)

    real(kind=kd), dimension(:,:), intent(in out) :: u, v, p
    real(kind=kd), intent(in) :: dt, dx, dy, rho, nu
    integer, intent(in) :: nt

    real(kind=kd), dimension(:,:), allocatable :: un, vn
    integer :: i, nx, ny, mx, my

    if (allocated(un)) deallocate(un)
    if (allocated(vn)) deallocate(vn)

    nx = size(u,1)
    ny = size(u,2)
    allocate(un(nx,ny))
    allocate(vn(nx,ny))

    mx = nx - 1
    my = ny - 1
        
    do i = 1, nt
       un(:,:) = u(:,:)
       vn(:,:) = v(:,:)

       p(:,:) = ppe(p, rho, dt, u, v, nx, ny, dx, dy, nit)
       
       u(2:mx,2:my) = un(2:mx,2:my) - &
            un(2:mx,2:my)*dt/dx*(un(2:mx,2:my)-un(2:mx,1:my-1)) - &
            vn(2:mx,2:my)*dt/dy*(un(2:mx,2:my)-un(1:mx-1,2:my)) - &
            dt/(2*rho*dx)*(p(2:mx,3:)-p(2:mx,1:my-1)) + &
            nu*(dt/dx**2*(un(2:mx,3:)-2*un(2:mx,2:my)+un(2:mx,1:my-1)) + &
            dt/dy**2*(un(3:,2:my)-2*un(2:mx,2:my)+un(1:mx-1,2:my)))

       v(2:mx,2:my) = vn(2:mx,2:my) - &
            un(2:mx,2:my)*dt/dx*(vn(2:mx,2:my)-vn(2:mx,1:my-1)) - &
            vn(2:mx,2:my)*dt/dy*(vn(2:mx,2:my)-vn(1:mx-1,2:my)) - &
            dt/(2*rho*dy)*(p(3:,2:my)-p(1:mx-1,2:my)) + &
            nu*(dt/dx**2*(vn(2:mx,3:)-2*vn(2:mx,2:my)+vn(2:mx,1:my-1)) + &
            (dt/dy**2*(vn(3:,2:my)-2*vn(2:mx,2:my)+vn(1:mx-1,2:my))))

       u(1,:) = 0; v(1,:) = 0
       u(:,1) = 1; v(:,1) = 0
       u(:,ny) = 0; v(:,ny) = 0
       u(nx,:) = 0; v(nx,:) = 0
    end do
    
  end subroutine cavityFlow
  
  function ppb(rho, dt, u, v, nx, ny, dx, dy)
    !! Solves a portion of the pressure-Poisson equation:
    !! $$ b = \left[ \frac{1}{\Delta t} \left( \frac{u_{i+1,j} - u_{i-1,j}}{2\Delta x} +
    !! \frac{v_{i,j+1} - v_{i,j-1}}{2\Delta y} \right) \right] $$

    integer, intent(in) :: nx, ny
    real(kind=kd), intent(in) :: rho, dt, dx, dy
    real(kind=kd), dimension(:,:), intent(in) :: u, v
    real(kind=kd), dimension(1:nx,1:ny) :: ppb
    integer :: mx, my

    mx = nx-1
    my = ny-1

    ppb(2:mx,2:my) = rho*(1/dt*((u(2:mx,3:)-u(2:mx,1:my-1))/(2*dx)+(v(3:,2:my)-v(1:mx-1,2:my))/(2*dy)) - &
         ((u(2:mx,3:)-u(2:mx,1:my-1))/(2*dx))**2 - &
         2*((u(3:,2:my)-u(1:mx-1,2:my))/(2*dy)*(v(2:mx,3:)-v(2:mx,1:my-1))/(2*dx)) - &
         ((v(3:,2:my)-v(1:mx-1,2:my))/(2*dy))**2)
    
  end function ppb

  function ppe(p, rho, dt, u, v, nx, ny, dx, dy, nit)
    !! Solves the pressure-Poisson equation:
    !! $$ p_{i,j}^n = \frac{\left(p_{i+1,j}^n + p_{i-1,j}^n \right) \Delta y^2 +
    !! \left( p_{i,j+1}^n + p_{i,j-1}^n \right) \Delta x^2}{2 \left(\Delta x^2 + \Delta y^2 \right)}
    !! - b $$
    !! where
    !! $$ b = \frac{1}{\Delta t} \left( \frac{u_{i+1,j} - u_{i-1,j}}{2\Delta x} +
    !! \frac{v_{i,j+1} - v_{i,j-1}}{2\Delta y} \right) - \left(\frac{u_{i+1,j} - u_{i-1,j}}{2\Delta x} \right)^2 - 2 \frac{u_{i,j+1}-u_{i,j-1}}{2\Delta y} \frac{v_{i+1,j}-v_{i-1,j}}{2\Delta x} - \left(\frac{v_{i,j+1}-v_{i,j-1}}{2\Delta y}\right)^2 $$

    integer, intent(in) :: nx, ny, nit
    real(kind=kd), intent(in) :: rho, dt, dx, dy
    real(kind=kd), dimension(:,:), intent(in) :: u, v
    real(kind=kd), dimension(:,:) :: p
    integer :: mx, my
    real(kind=kd), dimension(1:nx, 1:ny) :: b, pn, ppe

    mx = nx - 1
    my = ny - 1

    pn(:,:) = 0
    b(:,:) = ppb(rho, dt, u, v, nx, ny, dx, dy)
   
    do i = 1, nit
       pn(:,:) = p(:,:)
       p(2:mx,2:my) = ((pn(2:mx,3:)+pn(2:mx,1:my-1))*dy**2+(pn(3:,2:my)+pn(1:mx-1,2:my))*dx**2) / &
            (2*(dx**2+dy**2)) - &
            dx**2*dy**2/(2*(dx**2+dy**2))*b(2:mx,2:my)

       p(nx,:) = p(mx,:)
       p(1,:) = p(2,:)
       p(:,1) = p(:,2)
       p(:,ny) = 0
    end do

    ppe(:,:) = p(:,:)
    
  end function ppe
  
end module solver
