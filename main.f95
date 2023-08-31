program main

   use comvar
   use mUMFPACK
   implicit none

   real, dimension(:), allocatable :: Con, Con1,Prim, Prim_Bar, res
   a =1
   Nx = 50*a
   Ny = 50*a
   Nz = 50*a
   
  
   itmax = 50000
   itsave= 50



   max_pit     = 10
   RESTOL      = 1.0e-6
   cfl         = 0.1


   
   SM  = 3*NZ-2

   ! file id for saving solution
   fileid_sol = 0
   fileid_omg = 0



   allocate( Con( nvar*(Nx+2)*(Ny+2)*(Nz+2) ) )
   allocate( Con1( nvar*(Nx+2)*(Ny+2)*(Nz+2) ) )
   allocate( Prim( nvar*(Nx+2)*(Ny+2)*(Nz+2)) )
   allocate( res( nvar*(Nx+2)*(Ny+2)*(Nz+2) ) )
   allocate( Prim_Bar( nvar*(Nx+2)*(Ny+2)*(Nz+2))) 

   call solveFVM(Con, Con1, Prim,Prim_Bar, res)

end program main
