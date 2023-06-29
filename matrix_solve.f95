subroutine matrix_solve( Con1, Con, Prim_Bar)
   use comvar
   implicit none

   real :: Con1(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: Con(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: Prim_Bar(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   
   integer :: i, j
   integer :: Ap(0:Nz)
   integer :: Ai(SM)
   real    :: Ax(SM)
   real    :: B(Nz), Sol(Nz)
   
   
   
   do i= 1, Nx
       do j = 1, Ny
         B =0.0
        Sol =0.0
         call compute_B(B, Con1(:,i,j,:), Prim_Bar(5,i,j,:))
         call local_mat(Ap, Ai, Ax,Prim_Bar(:,i,j,:))
         call solution(Ap, Ai, Ax, B, Sol)  
         Con(1, i,j,1:Nz ) = Sol(1:Nz)
       enddo
   enddo
   
end subroutine matrix_solve


subroutine compute_B(B, co0, co_prim)

    use comvar
    implicit none
    
    real :: B(Nz)
    real :: co0(nvar, 0:Nz+1)
    real :: co_prim(0:Nz+1)
    
    integer :: i
    
    do i = 1, Nz
       
       B(i) = co0(5,i) - co_prim(i) * dt/(2.0*dz)*(co0(4,i+1) - co0(4,i-1)) &
              + co_prim(i) * dt**2/(2.0*dz)*(co0(1, i+1) - co0(1,i-1))      &
              - dt**2/(2.0*dz)*(co0(5, i+1) - co0(5,i-1)) 
     
     enddo
     
end subroutine compute_B     


subroutine local_mat(Ap, Ai, Ax, con_bar)
    use comvar
    implicit none
    
    integer :: Ap(0:Nz)
    integer :: Ai(SM)
    real    :: Ax(SM), con_bar(nvar, 0:Nz+1)
    
   
    !real :: po(nvar,0:Nz+1)
    real :: c_sq
    
    integer :: i
    integer :: jl, jm, jr
    
    jl = 3
    jm = 4
    jr = 5
    
    Ap(0) = 0
    Ap(1) = 2
    
    Ai(1) = 0
    Ai(2) = 1
    
    c_sq  = gamma*con_bar(2, 1)/con_bar(1,1)
    Ax(1) = 1.0 + dt**2/dz + (dt/dz)**2* c_sq
    
    
    Ax(2) = -dt**2/dz*(c_sq/dz - 1.0/2.0)
    
    
    do i = 2, Nz-1
       
       Ap(i) = Ap(i-1) + 3
       
       Ai(jl) = MOD(i-2+Nz, Nz)
       Ai(jm) = MOD(i-1+Nz, Nz)
       Ai(jr) = MOD(i  +Nz, Nz)
       
       
       c_sq   = gamma*con_bar(2,i)/con_bar(1,i)
       
       
       Ax(jm) =  1.0 + 2.0* (dt/dz)**2* c_sq
       Ax(jm) =  -dt**2/dz*(c_sq/dz - 1.0/2.0)
       Ax(jl) =  -dt**2/dz*(c_sq/dz + 1.0/2.0)
       
       jl = jl + 3
       jm = jm + 3
       jr = jr + 3
       
    enddo
    
    Ap(Nz) = Ap(Nz-1) + 2
    
    Ai(3*Nz-3) = Nz - 1
    Ai(3*Nz-2) = Nz
    
    c_sq   = gamma*con_bar(2,Nz)/con_bar(1,Nz)
    
    Ax(3*Nz-3) = -dt**2/dz*(c_sq/dz + 1.0/2.0)
    Ax(3*Nz-2) = 1.0 + 2.0* (dt/dz)**2* c_sq
    
    
end subroutine 
       
       
    
    
    
        
     



