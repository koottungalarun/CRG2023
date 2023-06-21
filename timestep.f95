subroutine timestep(Con,Prim_Bar)
   use comvar
   implicit none
   real :: Prim_Bar(nvar, 0:nx+1, 0:Ny +1, 0:Nz+ 1)
   real :: Con(nvar, 0:nx+1, 0:Ny +1, 0:Nz+ 1)
   integer :: i, j, k
   real :: speed, eig, eig0(3)
   speed=0.0
   eig  = 0.0
    do i =1, Nx
       do j = 1, Ny
          do k = 1, Nz
             call max_eig(Con(:,i,j,k), Prim_Bar(:,i,j,k),1, 0, 0, eig0(1))
             call max_eig(Con(:,i,j,k), Prim_Bar(:,i,j,k),0, 1, 0, eig0(2))
             call max_eig(Con(:,i,j,k), Prim_Bar(:,i,j,k),0, 0, 1, eig0(3))
              eig = maxval(abs(eig0))
              speed = max(speed, eig)
          enddo
       enddo
   enddo
   if ( speed .le. 0.0) then
      dt =  10E-05
   else 
      dt = dx/speed
   endif
   !if(mu > 0.0)then
    !  dt = min(dt, dx**2/mu)
   !endif

   dt = cfl*dt

   write(*,*)'Time step =', dt


end subroutine timestep
