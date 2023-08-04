subroutine ini_smooth( co_dash,co_bar)
   use comvar
   implicit none

   real    :: co_bar(nvar,0:Nx+1, 0:Ny+1, 0:Nz+1)
   real    :: co_dash(nvar,0:Nx+1, 0:Ny+1, 0:Nz+1)
   real    :: Rho(1:Nx, 1:Ny, 1:Nz), TDC
   integer :: i, j, k
   real    :: x, y, z, r, r_c
   
   ! Here co(1,:,:,:) = Density
   !    co(2:4,:,:,:) = Velocity 
   !      co(5,:,:,:) = Pressure
   
   ! In case of bar variables
   ! co(2, :, :, :) = Presser_Bar
   
      

   final_time   = 1.0
   
   Theta_Bar    = 300.0
   !Cp           = 1.0
   g            = 10.0
   P0           = 8.61 * 10E+04
   R0           = 8.31
   gamma        = 1.4
   Theta_dashc  = 0.5
   r_c          = 2.0

  

   xmin = -10.0
   xmax =  10.0
   
   ymin = -10.0
   ymax =  10.0
   
   zmin = 0.0
   zmax = 10.0
   
   
  

   dx = (xmax - xmin)/Nx
   dy = (ymax - ymin)/Ny
   dz = (zmax - zmin)/Nz
   
   co_dash(2, :, :, :) = 0.0  ! Velocity is zero
   co_dash(3, :, :, :) = 0.0
   co_dash(4, :, :, :) = 0.0
   
     ! Velocity is zero 
   co_bar(2, :, :, :) = 0.0
   co_bar(3, :, :, :) = 0.0
   co_bar(4, :, :, :) = 0.0
   
   co_bar(5, :, :, :) = Theta_Bar 
   
   
   
   
   
 

  

   do i=1,Nx
      do j = 1,Ny
         do k = 1, Nz
             x = xmin + (i-1)*dx + 0.5*dx
             y = ymin + (j-1)*dy + 0.5*dy
             z = zmin + (k-1)*dz + 0.5*dz
             
             co_bar(1, i, j, k) = P0/(Theta_Bar*R0)*(1.0 - (gamma - 1.0)/gamma*z*g/(R0*Theta_Bar))**(1.0/(gamma-1.0))
             
             r = sqrt(x**2 + y**2 + (z-2.0)**2)
             
             if ( r > r_c ) then
                co_dash(5, i, j, k) = 0.0
                Rho( i, j, k) = P0/((Theta_Bar) *R0)*(1.0 - (gamma - 1.0)/gamma*z*g/(R0*Theta_Bar))**(1.0/(gamma-1.0))
             else
                TDC    = 2.0*(1.0+ cos(M_PI*r_c/2.0)**2)
                co_dash(5, i, j, k) =  TDC
                Rho( i, j, k) = P0/((Theta_Bar+ TDC) *R0)*(1.0 - (gamma - 1.0)/gamma*z*g/(R0*Theta_Bar))**(1.0/(gamma-1.0))
             endif
            enddo
        enddo
    enddo    
             
    co_dash(1, 1:Nx, 1:Ny, 1:Nz)   = Rho(:, :, :) - co_bar(1, 1:Nx, 1:Ny, 1:Nz)         
     

    write(*,*)'Set the initial condition'

end subroutine ini_smooth
