subroutine dash2con(Dash, Bar, Con)
   use comvar
   implicit none

   real :: Dash(nvar, 0:nx+1, 0:ny+1, 0:nz+1)
   real :: Con(nvar,  0:nx+1, 0:ny+1, 0:nz+1)
   real :: Bar(nvar,  0:nx+1, 0:ny+1, 0:nz+1)
   integer :: i, j, k
   
   
   

   do i= 0,nx+1
       do j = 1, ny+1
          do k = 1, nz +1
          
              Con(1, i, j, k) = Dash(1, i, j, k) + Bar(1, i, j, k)
              Con(2, i, j, k) = Dash(2, i, j, k) + Bar(2, i, j, k)
              Con(3, i, j, k) = Dash(3, i, j, k) + Bar(3, i, j, k)
              Con(4, i, j, k) = Dash(4, i, j, k) + Bar(4, i, j, k)
              Con(5, i, j, k) = Dash(5, i, j, k) + Bar(5, i, j, k)     
           enddo
       enddo
   enddo
   
end subroutine dash2con
