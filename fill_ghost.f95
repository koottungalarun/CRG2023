subroutine fill_ghost(con)
   use comvar
   implicit none

   real :: con(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)

   integer :: i, j, k, l

   ! change to wall; reverse the sign of normal velocity

   do j= 0,Ny+1
     do k = 0, Nz+1
        do l= 1, nvar
          
          con(l, 0,  j, k)   = con(l, 1,  j, k)
          con(l, Nx+1, j, k) = con(l, Nx, j, k)
          
        enddo
      enddo
   enddo

    do i= 0,Nx+1
     do k = 0, Nz+1
        do l= 1, nvar
          
          con(l, i,  0, k)   = con(l, i,  1, k)
          con(l, i, Ny+1, k) = con(l, i, Ny, k)
          
        enddo
      enddo
   enddo

   do i= 0,Nx+1
     do j = 0, Ny+1
        do l= 1, nvar
          
          con(l, i,  j, 0)   = con(l, i,  j, 1)
          con(l, i, j, Nz+1) = con(l, i, j, Nz)
          
        enddo
      enddo
   enddo
  
end subroutine fill_ghost
