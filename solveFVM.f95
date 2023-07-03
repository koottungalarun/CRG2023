subroutine solveFVM(Con, Prim, Prim_Bar, res)

   use comvar

   implicit none

   real :: Con(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: Con1(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: Prim(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: res(nvar,   0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: Prim_Bar(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)


   integer :: it, i, j, k, l
   real    :: time, tmin, tmax
   real    :: resid
   logical :: tostop


   ! set initial condition
   call init_cond(Prim, Prim_Bar)
   call fill_ghost(Prim)
   call fill_ghost(Prim_Bar)
   ! converting primitive var to conservative var
   
   do i = 1, Nx
        do j = 1, Ny
           do k = 1, Nz
              
             call prim2cons(Prim(:,i,j,k) , Con(:,i,j,k),Prim_Bar(:,i,j,k))
             
           enddo
        enddo
   enddo
  
   call fill_ghost(Con)
   call saveprim(0.0, Con)
   
   time   = 0.0
   it     = 0

   do while(time < final_time .and. it < itmax)
      
      call timestep(Con, Prim_Bar)
      
      ! Exactly match final time
      tostop = .false.
      if(time + dt > final_time)then
         dt = final_time - time
         tostop = .true.
      endif
      
      Con1(:,:,:,:) = Con(:,:,:,:)
      !call saveprim(time, Con1)
      
      call predictor(Con1, Prim_Bar, res)
   !   call saveprim(time, Con1)
   !   call saveprim(time, Prim_Bar)
     
   !   stop
      call fill_ghost(Con1)
      
      call matrix_solve(Con1, Con, Prim_Bar)
      
       
      call fill_ghost(Con)
       
      call update_Rho(Con1, Con, Prim_Bar)
       
      call fill_ghost(Con)
      
      call update_Q3(Con1, Con,Prim_Bar)
       
      !call update_Q1(Con1, Con)
      Con(2,:, :, :) = Con1(2, :, :,:)
      !call update_Q2(Con1, Con)
      Con(3,:, :, :) = Con1(3, :, :,:)
     
      call fill_ghost(Con)
      it = it + 1
      
      
      do i = 1, Nx
         do j = 1, Ny
            do k = 1, Nz
              
             call cons2prims(Con(:,i,j,k), Prim(:,i,j,k) ,Prim_Bar(:,i,j,k))
             
            enddo
         enddo
      enddo
      call sol_residual(Con1, Con, resid)
      ! Compute min/max temperature
      tmin = +1.0e20
      tmax = -1.0e20
      do i=1,Nx
           do j = 1, Ny
               do k = 1, Nz
                  do l = 1, nvar
                      tmin = min(tmin, Con(l,i,j,k))
                      tmax = max(tmax, Con(l,i,j,k))
                  enddo
               enddo
           enddo
      enddo
 
      time = time + dt
      write(*,'(I6,F10.2,5E12.4)')it,time,tmin,tmax,resid

      if(mod(it,itsave)==0 .or. it==itmax .or. tostop)then
         call saveprim(time, Con)
      endif
     !stop 
   enddo ! time iteration loop


end subroutine solveFVM

subroutine update_Rho(Con1, Con, Prim_Bar)
    use comvar
    implicit none
    
    real :: Con1(nvar, 0: Nx+1, 0: Ny+1, 0: Nz+1)
    real :: Con(nvar, 0: Nx+1, 0: Ny+1, 0: Nz+1)
    real :: Prim_Bar(nvar, 0: Nx+1, 0: Ny+1, 0: Nz+1)
    
    integer :: i, j, k
    
    do i = 1, Nx
        do j = 1, Ny
           do k= 1, Nz
              
              Con(1, i,j,k) = Con1(1, i,j,k) + 1.0/Prim_Bar(5,i,j,k) &
                              *(Con(5,i,j,k) - Con1(5,i,j,k))
                              
           enddo
        enddo
    enddo
    
end subroutine update_Rho

subroutine update_Q3(Con1, Con,Prim_Bar)
    use comvar
    implicit none
    
    real :: Con1(nvar, 0: Nx+1, 0: Ny+1, 0: Nz+1)
    real :: Con(nvar, 0: Nx+1, 0: Ny+1, 0: Nz+1)
    real :: Prim_Bar(nvar, 0: Nx+1, 0: Ny+1, 0: Nz+1)
    real :: C_Bar(0: Nx+1, 0: Ny+1, 0: Nz+1)
    
    integer :: i, j, k
    
    C_Bar(:, :,:) = sqrt(gamma*Prim_Bar(2,:,:,:)/Prim_Bar(1,:,:,:))   
    
    do i = 1, Nx
       do j = 1, Ny
           do k =1, Nz
              
               Con(4, i,j,k) = Con1(4, i, j, k) - dt/(2.0*dz*Prim_bar(5,i,j,k)) &
                               *(C_Bar(i,j, k+1)**2*Con(5,i,j, k+1)- C_Bar(i,j, k-1)**2*Con(5,i,j, k-1)) &
                               - dt/(2.0*dz)*(Con(1,i,j, k+1)- Con(1,i,j, k-1))
                               
           enddo
       enddo
    enddo
    
    
end subroutine update_Q3
    

subroutine sol_residual(co0, co1, resid)
   use comvar

   implicit none

   real :: co0(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: co1(nvar,  0:Nx+1, 0:Ny+1, 0:Nz+1)
   real :: resid

   integer :: i, j, k, l

   resid = 0.0

   do i=1,Nx
        do j = 1, Ny
            do k = 1, Nz
                do l=1, nvar
                 resid = resid + (co1(l,i,j,k) - co0(l,i,j,k))**2
                enddo
            enddo
       enddo
   enddo

   resid = sqrt(resid/(Nx*Ny*Nz))

end subroutine sol_residual
