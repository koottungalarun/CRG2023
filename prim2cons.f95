subroutine prim2cons(Prim, Con, Prim_Bar)
   use comvar
   implicit none

   real :: Prim(nvar)
   real :: Con(nvar)
   real :: Prim_Bar(nvar)

   Con(1) = Prim(1)
   Con(2) = Prim(2)*(Con(1) + Prim_Bar(1))
   Con(3) = Prim(3)*(Con(1) + Prim_Bar(1))
   Con(4) = Prim(4)*(Con(1) + Prim_Bar(1))
   Con(5) = (Con(1) + Prim_Bar(1))*Prim(5) + Prim(1)*Prim_Bar(5)
 
end subroutine prim2cons
