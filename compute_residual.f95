subroutine compute_residual(Con, Prim_Bar, res)

   use comvar

   implicit none

   real :: Prim_Bar(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real ::      Con(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real ::      res(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: xflux, yflux, zflux
   real :: PrimL(nvar), PrimR(nvar), ConL(nvar), ConR(nvar)
   real :: lam, lam1, Q(nvar)
   integer :: i, j, k
   
   res = 0.0
   lam = 0.0
   lam1= 0.0 
   
   do i =1, Nx
       do j = 1, Ny
          do k = 1, Nz
             call max_eig(Con(:,i,j,k),Prim_Bar(:, i, j, k), 1, 0, 0, lam1)
             lam = max(lam, abs(lam1))
          enddo
       enddo
   enddo
  
   do i=0,Nx
      do j = 1, Ny
         do k = 1, Nz          
            PrimL = Prim_Bar(:,   i, j, k)
            PrimR = Prim_Bar(:, i+1, j, k)
          
            ConL  = Con(:,   i, j, k)
            ConR  = Con(:, i+1, j, k)
           ! do l = 1, nvar
           ! write(*,*) PrimL(l), PrimR(l), ConL(l), ConR(l)
            !enddo
            !stop
            call numflux_x( ConL, ConR, PrimL, PrimR, lam, xflux)
         
        
            res(:,  i, j, k) = res(:,  i, j, k) + dy*dz*xflux
            res(:,i+1, j, k) = res(:,i+1, j, k) - dy*dz*xflux
         enddo
       enddo
   enddo
   
    do i =1, Nx
       do j = 1, Ny
          do k = 1, Nz
             call max_eig(Con(:,i,j,k),Prim_Bar(:, i, j, k), 0, 1, 0, lam1)
             lam = max(lam, abs(lam1))
          enddo
       enddo
   enddo
   
    do j=0,Ny
      do i = 1, Nx
         do k = 1, Nz          
            PrimL = Prim_Bar(:, i, j,  k)
            PrimR = Prim_Bar(:, i, j+1,k)
            ConL  = Con(:,   i, j, k)
            ConR  = Con(:, i, j+1, k)
            call numflux_y( ConL, ConR, PrimL, PrimR, lam, yflux)
         
        
            res(:,  i, j, k) = res(:,  i, j, k) + dx*dz*yflux
            res(:,i, j+1, k) = res(:,i, j+1, k) - dx*dz*yflux
         enddo
       enddo
   enddo
   
    do i =1, Nx
       do j = 1, Ny
          do k = 1, Nz
             call max_eig(Con(:,i,j,k),Prim_Bar(:, i, j, k), 0, 0, 1, lam1)
             lam = max(lam, abs(lam1))
          enddo
       enddo
   enddo
   
    do k=0,Nz
      do i = 1, Nx
         do j = 1, Ny          
            PrimL = Prim_Bar(:, i, j,  k)
            PrimR = Prim_Bar(:, i, j,k+1)
            ConL  = Con(:,   i, j, k)
            ConR  = Con(:, i, j, k+1)
            call numflux_z( ConL, ConR, PrimL, PrimR, lam, zflux)
         
            call source(ConL, Q)
            res(:,i, j,   k) = res(:,  i, j, k) + dx*dy*zflux + dx*dy*dz*Q
            res(:,i, j, k+1) = res(:,i, j, k+1) - dx*dy*zflux
         enddo
       enddo
   enddo
   
end subroutine compute_residual


subroutine source(co, Q)
    use comvar
    implicit none
    real :: co(nvar), Q(nvar)
    
    Q(1:3) = 0.0
    Q(4)   = -co(1)
    Q(5)   = 0.0
    return
end subroutine source
    

subroutine max_eig(co ,primb, ix,iy,iz, u)
     use comvar
     implicit none
     real ::  u
     real :: prim(nvar), primb(nvar), co(nvar)
     integer :: ix, iy, iz
     real ::  c, Pressure, Rho
    
    call cons2prims(co, prim, primb)
    
    Rho       = (prim(1)+primb(1))
    
    Pressure  = P0*(R0*Rho/P0)**gamma
  

    
    c  = sqrt(gamma*Pressure/Rho)  
     
    
     
     u = ix*abs(prim(2)) + iy*abs(prim(3)) + iz*abs(prim(4))
     
     u = u  + c
     
     return

endsubroutine max_eig
     
     





