program main

   use comvar
   use mUMFPACK
   implicit none

   real, dimension(:), allocatable :: Con, Prim, Prim_Bar, res

   Nx = 20
   Ny = 20
   Nz = 20
   
  
   itmax = 50000
   itsave= 1



   max_pit     = 10
   RESTOL      = 1.0e-6
   cfl         = 0.1


   Pre_Bar     = 10E+05
   SM  = 3*NZ-2

   ! file id for saving solution
   fileid_sol = 0
   fileid_omg = 0



   allocate( Con( nvar*(Nx+2)*(Ny+2)*(Nz+2) ) )
   allocate( Prim( nvar*(Nx+2)*(Ny+2)*(Nz+2)) )
   allocate( res( nvar*(Nx+2)*(Ny+2)*(Nz+2) ) )
   allocate( Prim_Bar( nvar*(Nx+2)*(Ny+2)*(Nz+2))) 

   call solveFVM(Con, Prim,Prim_Bar, res)

end program main
