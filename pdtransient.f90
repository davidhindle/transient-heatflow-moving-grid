! program testparab
! test solution of a parabolic 1d differential equation
! parabolic form of the equation
! du/dt - (D(u) u') =  q  ; u(0,t)=a(t) or u'(0,t)=aa(t), u(L,t)=b(t) or u'(L,t)=bb(t)
! boundary conditions - should work with either Dirichlet or Neumann
!--------------------------------------------------------------------------------------------
! PARAMETERS
! D = diffusion coefficient 
! inicond = initial value of function
! dx = grid spacing for numerical solution
! dt = time step 
! maxt = duration of model run
! -------------------------------------------------------------------------------------------
! METHOD
! 
! equation reduces to a finite difference approximation.
! which yields  u_r+1 = A(u) u_r+1 + q_r + u_r
! which is rearranged to (I - dt*A(u_r)) u_r+1 = u_r + dt*q_r 
! which is a linear system (linearised, non-linear system)
! A u = q
! solve with tridiagonal matrix algorithm (Thomas algorithm)
!
! -------------------------------------------------------------------------------------------
! INPUT VARIABLES
! kon - number of layers in the model
! layer specific properties given layer by layer  
! k1(kon) - thermal conductivity
! ao(kon) - radiogenic heat production
!  
! qo(i) - source term initial values, not currently variable with time
! ic - initial condition (could be a steady state solution to the problem)
! OUTPUT VARIABLES 
! tr, tr1 stored and new temperature solutions
!_____________________________________________________________________________________________

use TDMA
implicit none

double precision, parameter :: zero=0.0d0
integer, parameter :: dp = kind(1.d0)
integer i,j,l,ct,ierr,nsub,n,kon,nv,bcf,il,ir
integer outno,ll,tsteps,tsave,ios
integer, dimension(:), allocatable :: kdn
double precision, dimension(:), allocatable :: u,q,qo,inis,tr,tr1,tro,trj,trc,p1,p2,p3,D,diff,hf,x
double precision, dimension(:), allocatable :: hfo,hfc,hfj,dc,dj
double precision :: depth,b,c,v,uplift,totuplift,rhm,shm,aom,kom
double precision, dimension(:), allocatable :: ko,k1,ao,shc,rho,shco,rhoo
double precision dt,dx,dx2
double precision ::  bc1, bc2, kr, hc, rh
logical :: bckind, radiative
double precision :: inicond,maxt,ic

 character out1*11, out2*11, cnla*3              ! set character variables for storing sequential files

! read in parameter file, paramtrans.txt 
open (11,file='paramtrans.txt')
read (11,*) radiative  ! radiative heat affects thermal conductivity
read (11,*) bckind     ! boundary condition type - false = dirichlet
read (11,*) kr         ! temperature cut off for radiative heat effect in °K
read (11,*) dx       ! grid space metres
dx2=dx**2
read (11,*) dt       ! time step years
dt=dt*365*24*60*60 
read (11,*) maxt     ! max time years - end time for run
maxt=maxt*365*24*60*60
read (11,*) outno    ! number of times to save
read (11,*) bc1      ! boundary condition surface
read (11,*) bc2      ! bundary condition base (this is neumann if bckind = true)
read (11,*) v        ! vertical velocity, m/yr
read (11,*) kon      ! number of layers in lithosphere model
allocate (k1(kon),kdn(kon+1),ao(kon),rhoo(kon),shco(kon))
kdn(1)=1
i=1
do while (i.le.kon)
15 read (11,*) k1(i),ao(i),rhoo(i),shco(i),depth ! reads layer properties (thermal cond., radiogenic heat, density, spec heat cap., depth to base in metres)
kdn(i+1)=int(depth/dx)                           ! determines node number of layer base given grid spacing and depth 
print *, 'node number of base of layer ', i, ' = ',  kdn(i+1)
i=i+1
end do
read (11,*) kom,aom,rhm,shm
! finish parameter read

! allocate vectors
nsub=kdn(i)+1  ! determines value of total number of nodes needed for depths and grid spacing given in parameters. This is automatic
print *, 'nsub = ',nsub
allocate (u(nsub),tr(nsub),tr1(nsub),q(nsub),qo(nsub),inis(nsub),p1(nsub),p2(nsub),p3(nsub),D(nsub),ko(nsub),diff(nsub),hf(nsub))
allocate (x(nsub+1),tro(nsub),trj(nsub),trc(nsub),hfo(nsub),hfc(nsub),hfj(nsub),dc(nsub),dj(nsub))
allocate (shc(nsub),rho(nsub))

! open log file 
 open (15, file='log.dat')
! set initial layer properties across all nodes corresponding to the layers
do j=1,kon
print *, 'layer top node = ', kdn(j), 'layer ko = ', k1(j)
 do i=kdn(j),kdn(j+1)
 ko(i)=k1(j)
 shc(i)=shco(j) 
 rho(i)=rhoo(j)
 qo(i)=(dt/(rho(i)*shc(i)))*ao(j)    ! q +ve, because (I - A) u_r+1 = u_r + q_r, for steady state opposite 0 = Au + q, so Au = -q  
 D(i)=ko(i)
 x(i)=-dx*(i-1)
 write (15,*) i,x(i)
 print *, i
 end do
end do
x(i)=x(i-1)-dx
x(i+1)=x(i)-dx
write (15,*) i, x(i)
write (15,*) i+1, x(i+1)
write (15,*) '************************'
! initialise vectors of non-zero matrix diagonals
p1 = zero ; p2= zero ; p3=zero 


! initialise temperature to initial condition
open (40,file='transic.dat',status='old')
i=1
do
 read (40,*,IOSTAT=ios) tr(i)
 if (ios /= 0) exit  ! Check if the read operation encountered end of file or error
 i=i+1
 if (i > nsub) exit  ! Ensure you don't exceed the size of the tr array
end do

print *, i, ' lines read from transic', tr(i-1)


! set initial thermal conductivity based on initial condition
 call jaupart(nsub,radiative,ko,kr,tr,D)

!set timesteps and save intervals correctly
tsteps=nint(maxt/dt)
tsave=nint(maxt/(outno*dt)) 

! write log file for run
 
  print *, 'tsteps = ', tsteps
  print *, 'dt = ', dt, ' seconds'
  print *, 'maxt = ', maxt, ' seconds'
  print *, 'tsave = ', tsave, ' steps'
  do i = 1,nsub,10
   write (15,*) i,D(i), tr(i)
  end do
  write (15,*) '*************************'
  

!timestepping loop: initial condition set
uplift = 0.d0
totuplift = 0.d0
do ll=1,tsteps
! integrate x(i) in time with v(x(i))
! this calculates proportion of v at a point linear varying from 1 (x=0) to 0 (x=nsub)
do i=2,nsub
x(i) = x(i) - (v*dt*(x(1)-x(nsub))/x(i))
end do

! test where x(2) is - x(2) is the node stretching away from or approaching x(1)
if (x(2).gt.-1.d0) then
 do i=1,nsub-1
  ! tr(i) = tr(i+1)
  ! ko(i) = ko(i+1)
  ! qo(i) = qo(i+1)
  ! rho(i) = rho(i+1)
  ! shc(i) = shc(i+1)
  ! write (15,*) tr(i),ko(i),qo(i),rho(i),shc(i) 
   print *, 'stretch critical' 
 end do
end if
!print *, ll, 'uplift = ', uplift
!if (uplift.ge.dx) then 
!print *, ll, 'uplift = ', uplift
!write(15,*) 'time = ', ll*dt/(360*24*3600), 'total uplift = ', totuplift
!uplift = 0.d0
! ko(nsub)=kom
! rho(nsub)=rhm
! shc(nsub)=shm
! qo(nsub)=(dt/(rhm*shm))*aom 
! do i=1,nsub-1
!  tr(i) = tr(i+1)
!  ko(i) = ko(i+1)
!  qo(i) = qo(i+1)
!  rho(i) = rho(i+1)
!  shc(i) = shc(i+1)
!write (15,*) tr(i),ko(i),qo(i),rho(i),shc(i)  
! end do
!write (15,*) ko(nsub),qo(nsub),rho(nsub),shc(nsub)
!write (15,*) '**********************************' 
! call jaupart(nsub,radiative,ko,kr,tr,D) 
!end if

! set up matrix
! print *, 'call matrix'
 call matrix(D,nsub,nv,bckind,x,dt,bc1,bc2,ao,shc,rho,tr,qo,q,p1,p2,p3)
! print *, 'called matrix'
! solve
 call tri(nsub,p1,p2,p3,q,tr1) !n,a,b,c,q,u
! update coefficient (D(u) = D(tr1))
 call jaupart(nsub,radiative,ko,kr,tr1,D)
! set tr to tr1
 tr=tr1


! naming and storing of solution files  
  if (mod(ll,tsave).eq.0) then
    out1='sol1000.dat'
! calculate heat flow
    do i=2,nsub-1
      hf(i)=((D(i+1)+D(i)+D(i-1))/3.d0) * (tr1(i+1) - tr1(i-1))/(2.d0*dx)
    end do
    hf(1)=hf(2)
    hf(nsub)=hf(nsub-1) 
    write (cnla,'(i3.3)') int(ll/tsave)
    out1(5:7)=cnla(1:3)
    print *, out1    
    open(16,file=out1)
    do i=1,nsub
     write (16,*) x(i), tr(i), D(i),hf(i)   
    end do
    close (16)
!    open (14, file=out1)

   end if

 200 continue
    
end do    

 close(15)
end program
       

! *********************************************************************************************
! SUBROUTINES
! *********************************************************************************************
! subroutine Jaupart
! calculate temperature dependent thermal conductivity
! acccording to algorithm of Jaupart and Mareschal

! define D(T,z) for different, temperature dependent models and solve
!________________________________________________________________________________
! ---- this is the Jaupart algorithm --------------------------------------------
! ---- includes radiative heat contribution to thermal conductivity--------------
! _______________________________________________________________________________
! input variables
! temperature tr(i)
! initial thermal conductivity ko(i)
! nsub - total grid nodes
! radiative - logical parameter true of false, if there is radiative heat transfer or not
! output D(i) - thermal conductivity of node i 


subroutine jaupart(nsub,radiative,ko,kr,tr,D)

implicit none
integer :: nsub,i,j,k,l
double precision :: kr,D(nsub),tr(nsub),ko(nsub)
logical :: radiative

    do i=1,nsub
     D(i) =  2.26d0 - 618.251d0/tr(i) + ko(i)*(355.576/tr(i) - 0.30247)
     if (radiative .eqv. .true.) then
      if (tr(i).gt.kr) then
       D(i) = D(i) + (tr(i)**3.d0) * 0.37d-9
      end if
     end if 
    end do

  D(nsub)=D(nsub-1)

end subroutine

! ------------------------------------------------------------------------------------------------------

! subtroutine matrix
! set the matrix diagonals and boundary conditions for the tdma matrix solver

subroutine matrix(D,nsub,nv,bckind,x,dt,bc1,bc2,ao,shc,rho,tr,qo,q,p1,p2,p3)

implicit none
integer :: nsub,nv,i,j
double precision, parameter :: zero=0.0d0
double precision :: D(nsub),ao(nsub),q(nsub),qo(nsub),tr(nsub),p1(nsub),p2(nsub),p3(nsub),rho(nsub),shc(nsub),x(nsub+1)
double precision :: dt, bc1, bc2
logical :: bckind

p1 = zero ; p2= zero ; p3=zero 

! matrix diagonals
! p1 = principal (a), p2 = upper (b), p3 = lower (c)
! a,b,c correspond to form of tridiagonal solver in module_tdma.f90
!-----------------------------------------------------------------------------------------------
 do i=2,nsub-1

 if  (((x(i+1)-x(i-1))*(x(i+1)-x(i))).eq.zero) then
 print *, 'divide by zero p2'
 exit
 end if
  p2(i)=1.d0/((x(i+1)-x(i-1))*(x(i+1)-x(i))) *(D(i)+D(i+1)) ! upper = b, b(1) = part of bc, or zero  
  p2(i) = -dt/(rho(i)*shc(i))*p2(i)
 end do
!print *, 'p2 finished'
! p2(1) - part of boundary condition

 do i=2,nsub-1 

 if  (((x(i+1)-x(i-1))*(x(i+1)-x(i))).eq.zero) then
 print *, 'divide by zero p2'
 exit
 end if

  p1(i)=-1.d0/((x(i+1)-x(i-1))*(x(i+1)-x(i)))*(D(i+1)+D(i))  -1.d0/((x(i+1)-x(i-1))*(x(i)-x(i-1)))*(D(i)+D(i-1))  ! principal = a, a(1) = bc, a(N) = bc   
  p1(i) = 1.d0 - dt/(rho(i)*shc(i))*p1(i)
 end do
!print *, 'p1 finished'
! p1(1) part of boundary condition
! p1(nsub) part of boundary condition

 do i=2,nsub-1
  
  p3(i)=1.d0/((x(i+1)-x(i-1))*(x(i)-x(i-1)))*(D(i)+D(i-1))   ! lower = c, c(1)= 0, c(N) = part of bc or zero  
  p3(i) = -dt/(rho(i)*shc(i))*p3(i)
 end do
!print *, 'p3 finished' 
do i=2,nsub
 q(i)=qo(i)+tr(i)
end do

p1(nsub)=p1(nsub-1)
p2(nsub)=p2(nsub-1)

! ---------------------------------------------------------------------------------------------------
! boundary conditions: choice between top dirichlet and base either dirichlet or neumann
! neumann at base = +ve value = heat flow in at base
! bc temps in °C 
! ---------------------------------------------------------------------------------------------------
if (bckind .eqv. .false.) then

!! this sets dirichlet boundary conditions on both sides of model (bckind = false)
!!---------------------------------------------------------------------------------------------------

p1(1) = 1.d0 !a(1)
p1(nsub) = 1.d0 !a(n)
p2(1) = 0.d0 !b(1)
p3(nsub) = 0.d0 !c(n)

!! and
!! ***


q(1)=bc1 + 273.d0
q(nsub)=bc2 + 273.d0

!print *, 'dirichlet selected'

elseif (bckind .eqv. .true.) then

!! this sets basal neumann boundary condition (bckind = true)
!!-------------------------------------------------------------------------------------------------------

p1(1) = 1.d0 !a(1) dirichlet, z=0
p2(1) = 0.d0 !b(1) dirichlet, z=0
p1(nsub) = 1.d0/(x(nsub-1)-x(nsub)) !a(n) neumann, z=B, u(n) 
p3(nsub) = -1.d0/(x(nsub-1)-x(nsub)) !c(n) numeann, z=B, u(n-1)


!! and
!! ***

q(1)=bc1 + 273.d0
q(nsub)=bc2/D(nsub)

!print *, 'neumann selected'

end if


end subroutine

