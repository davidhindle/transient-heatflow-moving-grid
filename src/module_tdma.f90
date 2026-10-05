!*******************************************************
!*    TDMA algorithm subroutine for F90                *
!*                                                     *
!*       D. Hindle, University of Göttingen            *
!*                                                     *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* Simone Sebben & B. Rabi Baliga (1995)               * 
!* SOME EXTENSIONS OF TRIDIAGONAL AND PENTADIAGONAL    *
!* MATRIX ALGORITHMS, Numerical Heat Transfer, Part B: * 
!* Fundamentals: An International Journal of           *
!* Computation and Methodology,28:3, 323-351,          * 
!* DOI:10.1080/10407799508928837                       *
!*                                                     *
!*                                                     * 
!*******************************************************
MODULE TDMA

CONTAINS

!  ***************************************************************
! Given an N x N matrix A, which is tridiagonal, 
! which is equivalent to a linear system
! au(i) + bu(i+1)  + cu(i-1) = d(i)
! where a = main diagonal, c is lr, b is ur diagonals and d is the load vector
! after rearranging for PDMA recursive algorithm 
! au(i) = d(i) - (bu(i+1) + cu(i-1))  hence, carry out -1*(b+c)
! store matrix entries in vectors a(N), b(N), c(N), d(N)
! where c(1), b(N), lie outside the matrix Aij and are always zero
! and use recurrence relationship 
! 
! forwards from i=2,N (i=1 first calculated directly)
! then set
! u(N) = q(N)  (these will be boundary conditions)
! implement backward recursion on
! u(i) = p(i) u(i+1) + q(i)
! from i=N-2,1
!  ***************************************************************

subroutine tri(n,a,b,c,d,u)

! a - main diagonal
! b - ur diagonals 
! c - lr diagonals
! p,q iteration variables
! d and u, vectors from system Au=d

! based on  Simone Sebben & B. Rabi Baliga (1995) SOME EXTENSIONS OF 
! TRIDIAGONAL AND PENTADIAGONAL MATRIX ALGORITHMS, Numerical Heat Transfer, Part B: 
! Fundamentals: An International Journal of Computation and Methodology,28:3, 323-351,
! DOI:10.1080/10407799508928837

! note, linear system is
! au(i) + bu(i+1) + cu(i-1)  = d(i)
! rearranging for TDMA recursive algorithm requires
! au(i) = d(i) - (bu(i+1) + cu(i-2))  
! hence, to correspond to Sebben's algorithm, carry out -1*(b+c)

implicit none



double precision :: p(n),q(n)
double precision :: a(n),b(n),c(n),d(n),u(n)
double precision :: bp(n),cp(n),dp(n),ep(n)
double precision :: edq(n),denom(n)
double precision, parameter :: zero=0.0d0
integer :: i,j,k,l,n

!initialise
edq=zero; ; denom=zero
p=zero; q=zero; u=zero

! rearrange linear system as described in Sebben - switches all values of b and c -ve including neumann bc term....
do i = 1,n
  bp(i)=-b(i)
  cp(i)=-c(i)
end do

!set known iteration variable values for i=1, i=2, (p,q)

p(1)= bp(1) / a(1)  ! = 0
q(1)= d(1) / a(1)  ! = 0

! carry out forward recursion
! i=2,n 

do i=2,n
  edq(i)=(d(i)+cp(i)*q(i-1))
  denom(i)=a(i)-cp(i)*p(i-1)
  p(i)= bp(i) / denom(i)
  q(i)= edq(i) / denom(i)
end do

! set end values of u (n) to start backward recursion
!print *, 'q(n) = ', q(n),p(n)
u(n)=q(n)

! backward recursion algorithm to solve for u (n-2,n-1,...,1)

do i=n-1,1,-1
  u(i)=q(i) + p(i)*u(i+1)
end do

return

END subroutine tri

!---------------------------------------------------------------------------------------------
! subroutine to multiply a pentadiagonal matrix and a vector (only elements of rows i=2, node-1)
! u = a*f

subroutine tdmm(n,a1,am1,ap1,f,u)

implicit none


double precision :: a1(n),am1(n),ap1(n),f(n),u(n)
integer :: i,j,k,l,n
u(1:n)=0.0d0

do i=2,n-1
  u(i) = am1(i)*f(i-1)+a1(i)*f(i)+ap1(i)*f(i+1)     
end do

return

END subroutine tdmm

END MODULE TDMA

! end of file pdma.f90

