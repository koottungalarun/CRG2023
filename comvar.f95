module comvar
   implicit none

   integer :: Nx, Ny, Nz, a
   real    :: xmin, xmax, ymin, ymax, zmin, zmax
   real    :: dx, dy, dz,  dt, dtp
   integer :: itmax
   integer :: itsave
   real    :: cfl
   real    :: final_time
   
   real    :: Theta_bar, Cp, g, P0, R0, gamma=1.4, Theta_dashc

  
   integer :: max_pit, SM
   real    :: RESTOL

  


   real :: M_PI = 4.0*atan(1.0)

   integer :: fileid_sol, fileid_omg

  




   
   integer :: nvar= 5, bnvar =3
   
end module comvar
