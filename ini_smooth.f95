subroutine ini_smooth(co_dash,co_bar)
   use comvar
   implicit none

   real    :: co_bar(nvar,0:Nx+1, 0:Ny+1, 0:Nz+1)
   real    :: co_dash(nvar,0:Nx+1, 0:Ny+1, 0:Nz+1)

   integer :: i, j, k
   real    :: x, y, z, r, r_c, Pi_e
   
   ! Here co(1,:,:,:) = Density
   !    co(2:4,:,:,:) = Velocity 
   !      co(5,:,:,:) = Pressure
   
      

   final_time   = 10.0
   
   Theta_Bar    = 300.0
   Cp           = 1.0
   g            = 9.8;
   P0           = 1.01 * 10E+05
   R0           = 8.31
   gamma        = 1.4
   Theta_dashc  = 0.5
   r_c          = 250

  

   xmin = 0.0
   xmax = 1000.0
   
   ymin = 0.0
   ymax = 1000.0
   
   zmin = 0.0
   zmax = 1000.0
   
   
  

   dx = (xmax - xmin)/Nx
   dy = (ymax - ymin)/Ny
   dz = (zmax - zmin)/Nz
   
   co_dash(2, :, :, :) = 0.0  ! Velocity is zero
   co_dash(3, :, :, :) = 0.0
   co_dash(4, :, :, :) = 0.0
   
   co_bar(2, :, :, :) = 0.0  ! Velocity is zero
   co_bar(3, :, :, :) = 0.0
   co_bar(4, :, :, :) = 0.0
   
   co_bar(5, :, :, :) = Theta_Bar 
   
 

  

   do i=1,Nx
      do j = 1,Ny
         do k = 1, Nz
             x = xmin + (i-1)*dx + 0.5*dx
             y = ymin + (j-1)*dy + 0.5*dy
             z = zmin + (k-1)*dz + 0.5*dz
             
             r = sqrt((x-500.0)**2 + (y-500.0)**2 + (z-350.0)**2)
             
             if ( r > r_c ) then
                co_dash(5, i, j, k) = 0.0
             else
                co_dash(5, i, j, k) =  Theta_dashc/2.0*(1+ cos(M_PI*r/r_c))
             endif
               
                
             Pi_e                = 1.0 - g*z/(Cp*co_bar(5, i, j, k))

             if (Pi_e .le. 0.0) then
                write (*,*) Pi_e, x, y, z
             endif
             
                
                co_bar(1, i, j, k)  = P0/(R*co_bar(5, i, j, k) )*(Pi_e)**(1.0/(gamma-1.0)) ! Rho_bar
                
                co_dash(1,i, j, k)  = -co_bar(1, i, j, k)*co_dash(5, i, j, k)/(co_bar(5, i, j, k))
                
                
            
            enddo
        enddo
    enddo    
             
              
     

    write(*,*)'Set the initial condition'

end subroutine ini_smooth
