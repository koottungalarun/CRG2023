subroutine cons2prims(Con, Prim, Prim_Bar)
   use comvar
   implicit none

   real :: Prim(nvar)
   real :: Prim_Bar(nvar)
   real :: Con(nvar), Rho
 

   
             
              Prim(1) = Con(1)
              Rho     = (Prim(1) + Prim_Bar(1))
              Prim(2) = Con(2)/Rho
              Prim(3) = Con(3)/Rho
              Prim(4) = Con(4)/Rho
              Prim(5) = (Con(5)-Con(1)*Prim_Bar(5))/Rho
 
end subroutine cons2prims
