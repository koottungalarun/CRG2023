subroutine cons2prims(Con, Prim, Prim_Bar)
   use comvar
   implicit none

   real :: Prim(nvar)
   real :: Prim_Bar(nvar)
   real :: Con(nvar), s
 

   
             
              Prim(1) = Con(1)
              s       = (Prim(1) + Prim_Bar(1))
              Prim(2) = Con(2)/s
              Prim(3) = Con(3)/s
              Prim(4) = Con(4)/s
              Prim(5) = (Con(5)-Con(1)*Prim_Bar(5))/s  
 
end subroutine cons2prims
