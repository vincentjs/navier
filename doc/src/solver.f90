module solver
  !! Contains the discretized 2D Navier-Stokes solver. The x-momentum equation is given by:
  !! $$ \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} = -\frac{1}{\rho} \frac{\partial p}{\partial x} + \nu \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) $$
  !! The y-momentum equation is similarly given by:
  !! $$ \frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} = -\frac{1}{\rho} \frac{\partial p}{\partial y} + \nu \left( \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right) $$
  !! Finally, the pressure equation is given by:
  !! $$ \frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2} = -\rho \left( \frac{\partial u}{\partial x} \frac{\partial u}{\partial x} + 2 \frac{\partial u}{\partial y} \frac{\partial v}{\partial x} + \frac{\partial v}{\partial y} \frac{\partial v}{\partial y} \right) $$

  use precision

  private :: ppb
  
contains

  pure real(kind=kd) function ppb(rho, dt, u, v, dx, dy) result(b)
    !! Solves a portion of the pressure-Poisson equation:
    !! $$ b = \left[ \frac{1}{\Delta t} \left( \frac{u_{i+1,j} - u_{i-1,j}}{2\Delta x} +
    !! \frac{v_{i,j+1} - v_{i,j-1}}{2\Delta y} \right) \right] $$

    real(kind=kd), intent(in) :: rho, dt, dx, dy
    real(kind=kd), dimension(-1:1,-1:1), intent(in) :: u, v

    b = ( (rho * dx**2 * dy**2) / (2 * (dx**2 + dy**2)) ) &
         * 1/dt * ( (u(1,0) - u(-1,0)) / (2*dx) + (v(0,1) - v(0,-1)) / (2*dy) ) &
         - ( (u(1,0) - u(-1,0)) / (2*dx) )**2 &
         - (-2 * ((u(0,1) - u(0,-1)) / (2*dy)) * ((v(1,0) - v(-1,0)) / (2*dx)) ) &
         - ( (v(0,1) - v(0,-1)) / (2*dy) )**2 
    
  end function ppb

  subroutine ppe(p, rho, dt, u, v, dx, dy)
    !! Solves the pressure-Poisson equation:
    !! $$ p_{i,j}^n = \frac{\left(p_{i+1,j}^n + p_{i-1,j}^n \right) \Delta y^2 +
    !! \left( p_{i,j+1}^n + p_{i,j-1}^n \right) \Delta x^2}{2 \left(\Delta x^2 + \Delta y^2 \right)}
    !! - b $$
    !! where
    !! $$ b = \frac{1}{\Delta t} \left( \frac{u_{i+1,j} - u_{i-1,j}}{2\Delta x} +
    !! \frac{v_{i,j+1} - v_{i,j-1}}{2\Delta y} \right) - \left(\frac{u_{i+1,j} - u_{i-1,j}}{2\Delta x} \right)^2 - 2 \frac{u_{i,j+1}-u_{i,j-1}}{2\Delta y} \frac{v_{i+1,j}-v_{i-1,j}}{2\Delta x} - \left(\frac{v_{i,j+1}-v_{i,j-1}}{2\Delta y}\right)^2 $$

    real(kind=kd), intent(in) :: rho, dt, dx, dy
    real(kind=kd), dimension(-1:1,-1:1), intent(in) :: u, v
    real(kind=kd), dimension(-1:1,-1:1), intent(in out) :: p
    real(kind=kd) :: b

    b = ppb(rho, dt, u, v, dx, dy)

    p(0,0) = ( (p(1,0) + p(-1,0)) * dy**2 + (p(0,1) + p(0,-1)) * dx**2 ) &
         / (2 * (dx**2 + dy**2)) &
         - b
    
  end subroutine ppe
  
  
  
end module solver
