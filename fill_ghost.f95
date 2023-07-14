subroutine fill_ghost(con)
   use comvar
   implicit none

   real :: con(nvar, 0:Nx+1, 0:Ny+1, 0:Nz+1)

<<<<<<< Updated upstream
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
=======
  
   
   
   
   con(1:3,:,:,0) =  con(1:3,:,:,1)
   con(4  ,:,:,0) = -con(4, :,:,1)
   con(5  ,:,:,0) =  con(5, :, :,1)
   
   con(1:3,:,:,Nz+1) =  con(1:3,:,:,Nz)
   con(4  ,:,:,Nz+1) = -con(4, :,:,Nz)
   con(5  ,:,:,Nz+1) =  con(5, :, :,Nz)
   
   con(1  ,0,:,:) =  con(1, 1, :, :)
   con(2  ,0,:,:) = -con(2, 1, :, :)
   con(3:5,0,:,:) =  con(3:5, 1, :,:)
   
   con(1  ,Nx+1,:,:) = con(1, Nx,:,:)
   con(2  ,Nx+1,:,:) = -con(2, Nx,:,:)
   con(3:5,Nx+1,:,:) = con(3:5, Nx, :,:)
   
   con(1:2  ,:,0,:) =  con(1:2,:,1,:)
   con(3  ,:,0,:)   = -con(3, :,1,:)
   con(4:5,:,0,:)   =  con(4:5, :, 1,:)
   
   con(1:2  ,:,Ny+1,:) = con(1:2,:,Ny,:)
   con(3  ,:,Ny+1,:)   = -con(3, :,Ny,:)
   con(4:5,:,Ny+1,:)   = con(4:5, :, Ny,:)
   
>>>>>>> Stashed changes
  
end subroutine fill_ghost
