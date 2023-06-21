subroutine init_cond(co0, co1)
   use comvar
   implicit none
   
   real    :: co0(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)
   real    :: co1(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)

   
   !call init_cond_square(co1)
   !call init_cond_ring(co1)
   call ini_smooth(co0, co1)

end subroutine init_cond
