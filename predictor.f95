subroutine predictor(co, Prim_Bar, res)
   use comvar
   implicit none

   real :: co(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: co1(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: Prim_Bar(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: res(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1), lambda
   integer :: i,j, k, l
  
   call compute_residual(co, Prim_Bar, res)
   lambda = dt/(dx*dy*dz)
   !call saveprim(0.0, res)
   co1(:,:,:,:) = co(:,:,:,:)
   do i = 1, Nx
      do j= 1, Ny
         do k= 1, Nz
             do l =1, nvar
              
              co(l,i,j,k) = co1(l,i, j,k)- lambda*res(l,i,j,k)
              
             enddo
         enddo
      enddo
   enddo
   
  
   

end subroutine predictor
