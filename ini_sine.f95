subroutine ini_sine( co_dash,co_bar)
   use comvar
   implicit none

   real    :: co_bar(nvar,0:Nx+1, 0:Ny+1, 0:Nz+1)
   real    :: co_dash(nvar,0:Nx+1, 0:Ny+1, 0:Nz+1)
   !real    :: Rho(1:Nx, 1:Ny, 1:Nz), TDC
   integer :: i, j, k
   real    :: x, y, z, r_c
   
   ! Here co(1,:,:,:) = Density
   !    co(2:4,:,:,:) = Velocity 
   !      co(5,:,:,:) = Pressure
   
   ! In case of bar variables
   ! co(2, :, :, :) = Presser_Bar
   
      

   final_time   = 0.1
   
   Theta_Bar    = 0.0
   !Cp           = 1.0
   g            = 10.0
   P0           = 8.61 * 10E+04
   R0           = 8.31
   gamma        = 1.4
   Theta_dashc  = 0.5
   r_c          = 2.0

  

   xmin = 0.0
   xmax =  2.0*M_PI
   
   ymin = 0.0
   ymax =   2.0*M_PI
   
   zmin = 0.0
   zmax =  2.0*M_PI
   
   
  

   dx = (xmax - xmin)/Nx
   dy = (ymax - ymin)/Ny
   dz = (zmax - zmin)/Nz
   
   co_dash(2, :, :, :) = 1.0  ! Velocity is zero
   co_dash(3, :, :, :) = 1.0
   co_dash(4, :, :, :) = 1.0
   co_dash(5, :, :, :) = 1.0
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
             co_dash(1, i,j, k) = 1.0 + 0.2*sin(x+y+z)
             
            enddo
        enddo
    enddo    
             
    co_bar(1, 1:Nx, 1:Ny, 1:Nz)   = 0.0 !Rho(:, :, :) - co_bar(1, 1:Nx, 1:Ny, 1:Nz)         
     

    write(*,*)'Set the initial condition'

end subroutine ini_sine
