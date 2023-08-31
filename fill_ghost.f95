subroutine fill_ghost(con)
   use comvar
   implicit none

   real :: con(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1), z
   integer :: k

   con(1  ,0,:,:) = con(1, 1, :, :)
   con(2  ,0,:,:) = con(2, 1, :, :)
   con(3  ,0,:,:) = con(3, 1, :, :)
   con(4  ,0,:,:) = con(4, 1, :, :)
   con(5  ,0,:,:) = con(5, 1, :, :)
   
   con(1  ,Nx+1,:,:) = con(1, Nx,:,:)
   con(2  ,Nx+1,:,:) = con(2, Nx,:,:)
   con(3  ,Nx+1,:,:) = con(3, Nx,:,:)
   con(4  ,Nx+1,:,:) = con(4, Nx,:,:)
   con(5  ,Nx+1,:,:) = con(5, Nx,:,:)

   
   con(1  ,:,0,:)   =  con(1,:,1,:)
   con(2  ,:,0,:)   =  con(2, :,1,:)
   con(3  ,:,0,:)   =  con(3, :,1,:)
   con(4  ,:,0,:)   =  con(4, :,1,:)
   con(5  ,:,0,:)   =  con(5, :,1,:)
   
   
   con(1  ,:,Ny+1,:)   = con(1,:,Ny,:)
   con(2  ,:,Ny+1,:)   = con(2, :,Ny,:)
   con(3  ,:,Ny+1,:)   = con(3, :,Ny,:)
   con(4,:,Ny+1,:)     = con(4, :, Ny,:)
   con(5,:,Ny+1,:)     = con(5, :, Ny,:)

  
   
   k = 0
   z = zmin + (k-1)*dz + 0.5*dz
   con(1,:,:,0)   =  P0/(Theta_Bar*R0)*(1.0 - (gamma - 1.0)/gamma*z*g/(R0*Theta_Bar))**(1.0/(gamma-1.0))
   con(2,:,:,0)   =  0.0
   con(3,:,:,0)   =  0.0
   con(4  ,:,:,0) =  0.0
   con(5  ,:,:,0) =  Theta_Bar 
   
   k = Nz+1
   z = zmin + (k-1)*dz + 0.5*dz
   con(1,:,:,Nz+1)   =  P0/(Theta_Bar*R0)*(1.0 - (gamma - 1.0)/gamma*z*g/(R0*Theta_Bar))**(1.0/(gamma-1.0))
   con(2,:,:,Nz+1)   =  0.0
   con(3,:,:,Nz+1)   =  0.0
   con(4  ,:,:,Nz+1) =  0.0
   con(5  ,:,:,Nz+1) =  Theta_Bar 
   
  
  
end subroutine fill_ghost
