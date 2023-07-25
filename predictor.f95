subroutine predictor(Con,Con1, Prim_Bar, res)
   use comvar
   implicit none

   real :: Con(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: Con1(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)
   !real :: temp(nvar, 1:Nx, 1:Ny, 1:Nz)
   real :: Prim_Bar(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: res(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1), lambda
   integer :: i,j, k, l
  
   call compute_residual(Con, Prim_Bar, res)
   lambda = dt/(dx*dy*dz)
   !call saveprim(0.0, res)
   
   Con1(:,1:Nx,1:Ny,1:Nz) = Con(:,1:Nx,1:Ny,1:Nz)
   !Con(1:nvar,1:Nx,1:Ny,1:Nz) = Con1(1:nvar,1:Nx, 1:Ny,1:Nz)- lambda*res(1:nvar,1:Nx,1:Ny,1:Nz)
   do i = 1, Nx
      do j= 1, Ny
         do k= 1, Nz
             do l =1, nvar
              
              Con(l,i,j,k) = Con1(l,i, j,k)- lambda*res(l,i,j,k)
              
             enddo
         enddo
      enddo
   enddo
   
  
   

end subroutine predictor
