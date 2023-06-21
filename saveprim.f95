subroutine saveprim(t,Primal)
   use comvar
   implicit none

   real    :: t
   real    :: Primal(nvar, 0:Nx+1, 0:Ny+1, 0:nz+1)
  
   integer :: i, j, k
   real    :: x, y, z
   character(len=512) :: filename

   filename = 'sol'
   call getfilename(filename, fileid_sol)

   open(10,file=trim(filename))
   write(10,*)'TITLE = "vortex flow"'
   write(10,*)'VARIABLES = "x", "y", "z", "Density", "Velx", "Vely", "Pressure"'
   write(10,*)'ZONE STRANDID=1, SOLUTIONTIME=',t,', I=',Nx,', J=',Ny, 'K=', Nz &
              ,' DATAPACKING=POINT'
  do k = 1, Nz
     do j=1, Ny
        do i=1, Nx
            x = xmin + (i-1)*dx + 0.5*dx
            y = ymin + (j-1)*dy + 0.5*dy
            z = zmin + (k-1)*dz + 0.5*dz

            write(10,'(7E24.14)') x, y,z, Primal(1,i,j,k), Primal(2,i,j,k), &
                              Primal(3,i,j,k), Primal(5,i,j,k)
        enddo
      enddo
   enddo

   close(10)

end subroutine saveprim
