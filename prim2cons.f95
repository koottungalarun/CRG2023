subroutine prim2cons(Prim, Con, Prim_Bar)
   use comvar
   implicit none

   real :: Prim(nvar)
   real :: Con(nvar)
   real :: Prim_Bar(nvar), Rho

   Rho    = (Con(1) + Prim_Bar(1))             ! Density
   Con(1) = Prim(1)                            ! Desity_dash
   Con(2) = Prim(2)*Rho                        ! Rho*vx
   Con(3) = Prim(3)*Rho                        !Rho*vy
   Con(4) = Prim(4)*Rho                        ! Rho*vz
   
   ! Computation of (Rho*theta)_dash  
   Con(5) = Rho*Prim(5) + Prim(1)*Prim_Bar(5)
 
end subroutine prim2cons
