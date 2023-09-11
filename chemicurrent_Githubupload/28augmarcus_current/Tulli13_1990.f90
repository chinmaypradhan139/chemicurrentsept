module mod_global

implicit none
integer::nsize,INFO,total_dimensions,Hi,Ne,nold,iflag,jflag,fflag 
  real*8 :: tau,planks_constant,kT,Vr,Energy_of_diabat,U0,dt,dtc,total_time
  complex*16 :: population
  real*8 :: time,rnd,omega,gh,dG,Band_width,RAND,rhos
  complex*16,dimension(:,:), allocatable :: c,c_old
  real*8,dimension(:,:), allocatable :: Energy_hamil,g,Identity,H,Grd,vdij
  real*8,dimension(:,:,:), allocatable :: Gradient_hamil,acw,old_Gradient_hamil
  real*8,dimension(:), allocatable :: Energy,acc1,E_levels,w_levels,knot_x
  real*8,dimension(:), allocatable ::pos,v,old_acc2,momentum,population_mat
  real*8, dimension(:), allocatable :: mass,pop_mat,e_metal,imp_mat,long_mat,popstate,chem_mat
  integer,dimension(:), allocatable :: lambda,old_lambda,diabstate
  complex*16, dimension(:), allocatable :: cwork
  integer, dimension(:,:), allocatable :: Hi_sp
  real*8, dimension(:,:), allocatable :: W_mat
  real*8, dimension(:,:), allocatable :: Ne_overlap
  real*8, dimension(:,:), allocatable :: new_force,delr,delp
  real*8, dimension(:,:), allocatable :: old_Energy_hamil  
  integer(kind=4) ::iseed
  integer :: ham_sp,hp,ntraj,locate,neu
  integer :: n_q,iter,traj_no
  integer :: hop_flg,collapse_flg 
  real*8, dimension(:), allocatable :: inpot
  integer, dimension(:),allocatable :: outputs
  real*8, dimension(:,:), allocatable :: tdata
  integer :: outdims,out_freq,pop_avg_time
  real*8 :: barrier
  integer, dimension(:), allocatable :: bins
  real*8, dimension(:,:), allocatable :: out_dist,avg_dist
  
contains
!.......................................................................
subroutine setup_initial_values1
implicit none
integer(kind=4) :: O,j
integer(kind=4), dimension(:),allocatable::x
real*8 :: wt
integer:: ut,xt,ip,kn,i
integer :: dims,parameters

parameters=17
outdims=7
dims=parameters+outdims

allocate(inpot(dims))

 open(25,file='fort.23')
 do ip=1,dims
 read(25,*) inpot(ip)
 enddo
 close(25)

allocate(outputs(outdims))

do j=1,outdims
    outputs(j)=int(inpot(parameters+j))
enddo
    
!write(11,*) inpot

Hi=inpot(1)
Ne=int(Hi/2)
tau=inpot(4)
Band_width=inpot(5)
ham_sp=int(inpot(6))
omega=inpot(7)
gh=inpot(8)
dG=inpot(9)
KT=inpot(10)
dtc=inpot(11)
dt=inpot(12)
total_dimensions=inpot(13)
wt=inpot(14)
total_time=wt*5000
!write(28,*) wt
rhos=real(Hi)/Band_width
Vr=sqrt(tau/(2*3.1416))
!barrier=0.022
barrier=0.d0
pop_avg_time=1000
out_freq=int(total_time/pop_avg_time)
  allocate(pos(total_dimensions))
  allocate(v(total_dimensions))
  allocate(momentum(total_dimensions))
  allocate(acc1(total_dimensions))
  allocate(old_acc2(total_dimensions))
  allocate(acw(Hi,Hi,total_dimensions))
  allocate(H(Hi,Hi))
  allocate(Energy_hamil(Hi,Hi)) 
  allocate(Energy(Hi))
  allocate(Gradient_hamil(Hi,Hi,total_dimensions))
  allocate(c(Ne,Hi))
  allocate(g(Ne,Hi))
  allocate(lambda(Ne))
  allocate(diabstate(Ne))
  allocate(Identity(Hi,Hi))
  allocate(Grd(Hi,Hi))
  allocate(vdij(Hi,Hi))
  allocate(population_mat(int(total_time/dtc)))
  allocate(E_levels(int(Hi/2)))
  allocate(knot_x(int(Hi/2)))
  allocate(w_levels(int(Hi/2)))
  allocate(mass(total_dimensions))
  allocate(tdata(int(total_time/dtc),outdims))
  allocate(Hi_sp(Ne**2+1,Ne))
  allocate(W_mat(Hi,Hi))
  allocate(Ne_overlap(Ne**2+1,Ne**2+1))
  allocate(new_force(Ne**2+1,total_dimensions))
  allocate(old_Energy_hamil(Hi,Hi))
  allocate(old_Gradient_hamil(Hi,Hi,total_dimensions))
  allocate(c_old(Ne,Hi))
  allocate(old_lambda(Ne))
  allocate(delr(Ne**2+1,total_dimensions))
  allocate(delp(Ne**2+1,total_dimensions))
  allocate(e_metal(Hi))
  allocate(bins(Hi))
  allocate(out_dist(out_freq,Hi))
  allocate(avg_dist(out_freq,Hi))

call random_seed(size=O)
  allocate(x(O))
   do j=1,O
   x(j)=j**6*iseed+23452132
   enddo
  call random_seed(put=x)
 open(2,file='rndom',status='new')
 write(2,*) iseed,x,O
 


open(167,file='raw_x.txt')
do kn=1,int(Hi/2)
read(167,*) knot_x(kn)
enddo


open(169,file='raw_w.txt')
do xt=1,int(Hi/2)
read(169,*)w_levels(xt)
enddo
close(169)

do i=1,int(Hi/2)
    e_metal(int(Hi/2)-i+1)=-band_width/2.d0*(1/2.d0+knot_x(i)/2.d0)
    e_metal(int(Hi/2)+i)=band_width/2.d0*(1/2.d0+knot_x(i)/2.d0)
!    Vc(nquant/2-i+1)=sqrt(gama_coup/(2*pi))*sqrt(band_width/2.d0*weights(i))
!    Vc(nquant/2+i)=sqrt(gama_coup/(2*pi))*sqrt(band_width/2.d0*weights(i))
enddo


mass(1)=2000


write(66,*) 1/(1+exp(dG/KT))

bins=0
end subroutine
!.........................................................
  
subroutine gaussian_random_number(rnd0)
!USE IFPORT
   !! generates gaussian distribution with center 0, sigma 1
   !! q0+sig*rnd gives center=q0, sigma=sig
   implicit none
   integer(kind=4) :: n,j,M,O,k
   real*8,intent(out)::rnd0
   real*8 rnd1,rnd2,pi
   pi=dacos(-1.d0)
   call random_number(rnd1)
   call random_number(rnd2)
   rnd0 = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!.............................................................
function fdist(x)
real*8 x,fdist

fdist=1/(1+exp(x/KT))


end function
!............................................................
subroutine setup_initial_values2
implicit none
integer :: i



call potential(ham_sp)


diabstate=0

do i=1,Ne
    diabstate(i)=i+1
enddo



end subroutine 
!........................................................................
subroutine Metropolis_Hastings 
integer :: i,j,k,l,binary(Hi),n
real*8 rnd6, start,finish,pr,nor
binary=0
do k=1,1000
    call CPU_time(start)
    do l=1,10000
        call random_number(rnd6)
        iloop : do i=1,Ne
             pr=0.d0
             nor=fdist(H(diabstate(i),diabstate(i)))

             do n=2,Hi-1
                if (.not.(any(n==diabstate))) then
                   nor=nor+fdist(H(n,n))*(1-fdist(H(diabstate(i),diabstate(i))))
                 end if
             enddo
            
             pr=fdist(H(diabstate(i),diabstate(i)))/nor
 !            write(135,*) pr,nor
             if (pr>rnd6) cycle
             jloop : do j=2,Hi-1
                if (.not.(any(j==diabstate))) then
                        
                        pr=pr+fdist(H(j,j))*(1-fdist(H(diabstate(i),diabstate(i))))/nor
                        if (pr>rnd6) then
                           diabstate(i)=j
  !                          write(135,*) diabstate(i),j,pr,rnd6
                           exit jloop
                        end if
                end if
            enddo jloop
            write(136,*) i,j,pr
        enddo iloop
      call CPU_time(finish)
      write(124,*)k,l,finish-start
    enddo
    binary=binary+binstate(diabstate)
enddo

!do i=2,Hi
!    write(623,*) H(i,i),binary(i)/real(k),fdist(H(i,i))
!enddo
        
end subroutine 
!.........................................................................
subroutine simple_fermi
integer :: binary(Hi),unfilled(Ne)
integer :: n,irnd1,i,irnd2,k,l
real*8 rnd7,rnd8,pr,rnd9,tot_e
integer :: forb,vorb

binary=0


n=1
do i=2,Hi
    if (.not.(any(i==diabstate))) then
        unfilled(n)=i
        n=n+1
    end if
enddo

do k=1,1
    do l=1,10000
        call random_number(rnd7)
        irnd1=ceiling(Ne*rnd7)
        forb=diabstate(irnd1)
        if (irnd1==1) cycle     
        call random_number(rnd8)
        irnd2=ceiling(Ne*rnd8)
        vorb=unfilled(irnd2)
        if (unfilled(irnd2)==1) cycle

        pr=dexp(-(H(vorb,vorb)-H(forb,forb))/KT)
        call random_number(rnd9)
        if (pr>rnd9) then
            diabstate(irnd1)=vorb
            unfilled(irnd2)=forb
        end if
    enddo

    binary=binary+binstate(diabstate)
enddo


do i=2,Hi
    write(623,*) H(i,i),binary(i)/real(k),fdist(H(i,i))
enddo

tot_e=0.d0
do i=2,Hi
    tot_e=tot_e+fdist(H(i,i))
enddo
    

write(673,*) tot_e
write(673,*)
end subroutine 
!.........................................................................
function ptrans(m,n)
real*8 ptrans
integer m,n

ptrans=fdist(H(n,n))*(1-fdist(H(m,m)))

end function 
!..........................................................................
subroutine diag_wrapper(matrix,nsize,eigen_values,eigen_vectors)
real*8, intent(inout) :: matrix(nsize,nsize)
real*8, dimension(nsize,nsize) :: mat
integer LWORK,nsize
real*8, allocatable :: WORK(:)
real*8, intent(out) :: eigen_vectors(nsize,nsize),eigen_values(nsize)
mat=matrix
LWORK=3*nsize-1
allocate(WORK(LWORK))
call dsyev('V','U',nsize,mat,nsize,eigen_values,WORK,LWORK,INFO)
eigen_vectors=mat

end subroutine

!..........................................................................
subroutine logm(mat,log_mat,n)
   !! http://arxiv.org/pdf/1203.6151v4.pdf
   implicit none
   integer,intent(in):: n
   real*8,intent(in):: mat(n,n)
   real*8,intent(out):: log_mat(n,n)
   integer i
   complex*16 T(n,n),en(n),vect(n,n)
   complex*16 dd(n,n)
    
   call schur(mat,T,n,en,vect,nold,cwork)
   dd=0.d0
   do i=1,n
     dd(i,i)=cdlog(t(i,i)/cdabs(t(i,i)))
   enddo
      
   log_mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine logm


!.............................................................................

subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
   !! Diaganalizing matrix using dsyevr. First m_values eigen values and
!eigenvectors computed.
   !! The module's common variables should contain:

   !! Initialize nold=0

   !! nold makes sure that everytime value of n changes, work and iwork
!are re-allocated for optimal performance.
   !! mat is destroyed after use.

   implicit none
   integer,intent(in) :: n
   integer,intent(inout) :: nold
   complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
   real*8,intent(in) :: mat(n,n)
   complex*16,intent(out) :: T(n,n)
   complex*16,allocatable,intent(inout):: cwork(:)
   real*8 rwork(n)
   complex*16 mat_c(n,n)

   integer lwork
   logical:: select
   logical bwork(n)
   integer sdim,info,AllocateStatus

   T=mat

   info=0
   sdim=0

   if(nold.ne.n .or. .not.allocated(cwork)) then
   !if(nold.ne.n) then
     lwork=-1
     if(allocated(cwork))deallocate(cwork)
     allocate(cwork(n))
     call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
     lwork=int(cwork(1))
     deallocate(cwork)
     allocate(cwork(lwork),STAT=AllocateStatus)
     if(allocatestatus.ne.0) write(6,*)"problem in schur, allocation"
     nold=n
   endif

   lwork=size(cwork)
   call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)

 !  write(57,*) mat
   if(info.ne.0) then
     write(6,*) "problem in scur",info
     stop
   endif

end subroutine schur
!...........................................................................
subroutine nonadiabaticvector(t,u)
implicit none
integer :: acwi
integer, intent(in) :: t,u

!write(94,*) lambda(t)

acw=0

do acwi=1,total_dimensions

if (t.ne.u) then
acw(lambda(t),u,acwi)=(sum(Energy_hamil(:,lambda(t))&
*matmul(Gradient_hamil(:,:,acwi),Energy_hamil(:,u))))/(Energy(lambda(t))-Energy(u))
else
acw(lambda(t),u,acwi)=0
endif

end do
end subroutine
!..........................................................................
subroutine force
implicit none
integer :: acc1i,z,u,q,i,j
real*8, dimension(total_dimensions) :: m_acc1
!call potential(1)
!do acc1i=1,total_dimensions

!do z=1,Ne
!m_acc1(acc1i,z)=-((sum((Energy_hamil(:,lambda(z)))*matmul(Gradient_hamil(:,:,acc1i),Energy_hamil(:,lambda(z)))))/mass(1))
!end do
!end do
!acc1=sum(m_acc1(:,1:Ne))-omega**2*pos(1)

m_acc1=0.d0
do j=1,Ne
  do i=1,total_dimensions
    m_acc1=m_acc1-sum(Energy_hamil(:,lambda(j))*matmul(Gradient_hamil(:,:,i),Energy_hamil(:,lambda(j))))/mass(i)
  enddo
enddo

acc1=m_acc1-omega**2*pos(1)


end subroutine
!..........................................................................
subroutine Rungekutta
implicit none
integer :: i,j,p,q
complex*16,dimension(Hi,Hi) :: Vad
complex*16,dimension(Ne,Hi) ::k1,k2,k3,k4
complex*16,dimension(Hi,Hi) ::M

do i=1,Hi
do j=1,Hi 
if (i.eq.j) then
Vad(i,j)=Energy(i)/(cmplx(0,1))
else
Vad(i,j)=0
end if
end do
end do

do p=1,Hi
do q=1,Hi
M(p,q)=Vad(p,q)-sum(v*acw(p,q,:))
enddo
enddo

k1=dt*matmul(c,transpose(M))
k2=dt*matmul((c+k1/2),transpose(M))
k3=dt*matmul((c+k2/2),transpose(M))
k4=dt*matmul((c+k3),transpose(M))
c=c+(k1+2*k2+2*k3+k4)/6
end subroutine
!.........................................................................

subroutine Rungefutta
implicit none
integer :: i,j,p,q,x,y
complex*16,dimension(Hi,Hi) :: Vad
complex*16,dimension(Ne,Hi) ::k1,k2,k3,k4
complex*16,dimension(Hi,Hi) ::M
!real*8,intent(in), dimension(Hi,Hi) :: dij
do i=1,Hi
do j=1,Hi
if (i.eq.j) then
Vad(i,j)=Energy(i)/(cmplx(0,1))
else
Vad(i,j)=0
end if
end do
end do


do p=1,Hi
do q=1,Hi
M(p,q)=Vad(p,q)-vdij(p,q)


enddo
enddo





k1=dt*matmul(c,transpose(M))
k2=dt*matmul((c+k1/2),transpose(M))
k3=dt*matmul((c+k2/2),transpose(M))
k4=dt*matmul((c+k3),transpose(M))

c=c+(k1+2*k2+2*k3+k4)/6
end subroutine



!...................................................................................
subroutine evolve_quantum
implicit none
integer, dimension(Ne) ::p
integer :: i,d,r,s
complex*16 :: det_S1,det_S2
complex*16, dimension(Ne,Ne) ::  S1,S2
complex*16 :: bdi,Adi
real :: Aii,gid,rnd,pr

pr=0.d0

call rungefutta

call random_number(rnd)
jloop : do i=1,Ne
  do d=1,Hi
  !  if(.not.(any(d==lambda))) then
    if (FINDLOC(lambda,d,1).eq.0) then
        p=lambda

       do r=1,Ne
         S1(:,r)=c(:,p(r))
       enddo

       call modulus(S1,Ne,det_S1)

       p(i)=d


       do s=1,Ne
         S2(:,s)=c(:,p(s))
       enddo

       call modulus(S2,Ne,det_S2)



       Adi=det_S2*conjg(det_S1)
       bdi=-2*real(conjg(Adi)*vdij(d,lambda(i)))
       Aii=det_S1*conjg(det_S1)
       gid=dt*real(bdi/Aii)
       if (gid<0) then
           gid=0.d0
       end if         
       pr=pr+gid
             

       if (pr>rnd) then
          write(22,*) time,'tim'
 !         call nonadiabaticvector(i,d)
          write(22,*) lambda
          write(22,*) lambda(i),d
          call hop1(i,d)
          write(22,*) lambda
          write(22,*) pr,rnd
        
        
          exit jloop
       
       endif
     end if
   enddo
enddo jloop




end subroutine
!.........................................................................
subroutine classical_evolution
implicit none
integer :: p,r,TT,yt,i,j,x,y
real*8,dimension(Hi) :: signature
real*8 :: KE,TE,EE,rnd1,Aii,tr_A
complex*16 :: c_diab(Ne,Hi),modcsq
real*8 :: Ax(Ne**2+1)
real*8 ::sig1(Ne**2+1),sig2(Ne**2+1)

call potential(ham_sp)
call force



TT=int(total_time/dtc)


do while(time.le.total_time)

if (outputs(4)==1) call impurity_pop

if (outputs(1)==1) call populations

call check_long_pop

!call A_mat_pop

call construct_state_space(Hi_sp)

call statepop

!call full_rho_diab

call reset_flags




!call construct_state_space(Hi_sp)

tr_A=0.d0
do i=1,Ne**2+1
   call calculate_sigma(Hi_sp(i,:),Aii,1)
   tr_A=tr_A+Aii
enddo


do i=1,Ne**2+1
   call calculate_sigma(Hi_sp(i,:),sig1(i),1)
enddo


!call changing_A_basis(sig1,sig2)




!call calculate_sigma(Hi_sp(1,:),Ax(1),1)
!write(94,*) time,Ax(1)
!Ax=0
!do i=1,Ne**2+1
!    call calculate_sigma(Hi_sp(i,:),Ax(i),1)
!    write(94,*) i,Ax(i)
!end do
!write(94,*)


write(90,*)time,tr_A


if (tr_A<0.90) then
   write(91,*) time,tr_A,'oh_no'
end if
!!...check state space..........
!do i=1,Ne**2+1
!   write(87,*) Hi_sp(i,:)
!enddo
!..............................

call save_old_state
call velocity_verlet2

call signt(signature)

do r=1,Hi
 if (signature(r)<0) then
    Energy_hamil(:,r)=-Energy_hamil(:,r)
 endif
enddo

call vdotd

n_q=int(dtc/dt)
do while(n_q>0) 

call evolve_quantum

n_q=n_q-1
enddo


if (int(inpot(15)).eq.1) then
    !write(232,*)'hi'
    call verlet_decoherence
    call check_collapse
end if

!if (mod(int(time),500).eq.0) then
!  call full_density_matrix

!end if


!! checking decoherence effect.......
!if (mod(int(time),50000).eq.0) then
!   c=0
!   do i=1,Ne 
!      c(i,lambda(i))=1
!   enddo
!end if
!!...................................

!! checking energy conservation......
!EE=Energy(lambda(1))
!do i=2,Ne
!    EE=EE+Energy(lambda(i))
!enddo  
!write(44,*) time,0.5*mass(1)*omega**2*pos(1)**2+0.5*sum(mass*v*v)+EE,pos

!! checking mod ci**2...............

!do i=1,Ne
!   modcsq=sum(c(i,1:Hi)*conjg(c(i,1:Hi)))
!   write(45,*)time,i,sum(c(i,1:Hi)*conjg(c(i,1:Hi)))
!   if ((real(modcsq)<0.99).or.(real(modcsq)>1.01)) then
!       write(45,*) time,i,modcsq
!   end if
!enddo
!..................................





!! printing coeff matrix...............
!write (217,*) time
!do x=1,Hi
! do y=1,Hi
!   write(217,'(f12.8$)')H(x,y)
! enddo
!    write(217,*)
!enddo
!......................................



if (outputs(5)==1) call check_long_pop

yt=int(time/dtc)
population_mat(yt)=real(population)



time=time+dtc
enddo







end subroutine
!................................................................................
subroutine hop1(t,u)
implicit none
integer, intent(in) :: t,u
real*8 :: para_v
real*8, dimension(:),allocatable :: perp_v
real*8 :: cond


!write(22,*) lambda(t),u,'dd'

allocate(perp_v(total_dimensions))

acw(lambda(t),u,1)=vdij(lambda(t),u)/v(1)

cond=((sum(v*acw(lambda(t),u,:))/norm2(acw(lambda(t),u,:)))**2+(2*Energy(lambda(t))/mass(1))-(2*Energy(u)/mass(1)))

write(22,*)time,cond,acw(lambda(t),u,:)


if (cond>0) then
    hp=hp+1
    para_v=sum(v*acw(lambda(t),u,:))/norm2(acw(lambda(t),u,:))
    perp_v=v-para_v*acw(lambda(t),u,:)/norm2(acw(lambda(t),u,:))
    para_v=((para_v)/abs(para_v))*sqrt(para_v**2+(2*Energy(lambda(t))/mass(1))-(2*Energy(u)/mass(1)))
    v=(para_v)*acw(lambda(t),u,:)/norm2(acw(lambda(t),u,:))+perp_v
    lambda(t)=u
    delr=0.d0
    delp=0.d0
!    call transversed_surfaces(hp)
!    v=-v
    write(122,*) lambda
    hop_flg=1
end if

!Note: (sum(v*acw(lambda(t),u,:))/norm(acw(lambda(t),u,:))) is the parallel component of velocity in the direction of nonadiabatic coupling vector
end subroutine
!................................................................................


subroutine hop(t,u)
implicit none
integer, intent(in) :: t,u
integer :: reciprocal_mass_loop,v_loop
real*8, dimension(total_dimensions) :: reciprocal_mass
real*8 :: a,b,gama,frustrated_condition




do reciprocal_mass_loop=1,total_dimensions
     reciprocal_mass(reciprocal_mass_loop)=1/mass(reciprocal_mass_loop)
end do



       a=0.5*sum(reciprocal_mass*acw(lambda(t),u,:))
       b=sum(v*acw(lambda(t),u,:))
     
   
       frustrated_condition=b**2+4*a*(Energy(lambda(t))-Energy(u))

       if (frustrated_condition>0) then
           if (b<0) then
           gama=(b+sqrt(frustrated_condition))/2*a
           else
           gama=(b-sqrt(frustrated_condition))/2*a
           end if
           do v_loop=1,total_dimensions
             v(v_loop)=v(v_loop)-gama*acw(lambda(t),u,v_loop)/mass(v_loop)
           end do
           lambda(t)=u
        end if




end subroutine hop
!.................................................................................................
subroutine velocity_verlet
real*8 :: delr(total_dimensions),delv(total_dimensions)
real*8 :: gama_dt,gamma_B,c0,c1,c2
gamma_B=2*omega
gama_dt=gamma_B*dtc
c0=dexp(-gama_dt)
c1=1.d0/gama_dt*(1.d0-c0)
c2=1.d0/gama_dt*(1.d0-c1)
call stochastic_force(delr,delv)
pos=pos+c1*dtc*v+c2*dtc*dtc*acc1+delr
old_acc2=acc1


call potential(ham_sp)
call force
v=c0*v+(c1-c2)*dtc*old_acc2+c2*dtc*acc1+delv





end subroutine
!................................................................................
subroutine vdotd

!variable for checking sign of the wave
!real*8, intent(out),dimension(Hi,Hi) :: dij
real*8, dimension(Hi,Hi) :: Ut,logarithm_Ut
integer :: i,x,y







Ut=matmul(transpose(old_Energy_hamil),Energy_hamil)



!W_mat=Ut

call orthoganalize(Ut,Hi)
W_mat=Ut

call logm(Ut,logarithm_Ut,Hi)
vdij=logarithm_Ut/dtc


call construct_Ne_overlap(Ne_overlap)


!! check N electron ovelap matrix
!do x=1,Ne**2+1
! do y=1,Ne**2+1
!   write(61,'(f12.8$)')Ne_overlap(x,y)
! enddo
! write(61,*)
!enddo

!!..................................




end subroutine
!......................................................................................
subroutine orthoganalize(mat,n)
  integer,intent(in)::n
  real*8,intent(inout)::mat(n,n)
  real*8 S_mat(n,n)

  S_mat=matmul(transpose(mat),mat)
  call inverse_squareroot(S_mat,n)
  mat=matmul(mat,S_mat)

end subroutine orthoganalize
!.......................................................................................
subroutine inverse_squareroot(mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(inout):: mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=1.d0/t(i,i)**0.5d0
  enddo

  mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine inverse_squareroot



!........................................................................................
subroutine signt(signature)


!variable for checking sign of the wave

integer :: i
real*8, intent(out),dimension(Hi) :: signature
do i=1,Hi
signature(i)=sum(old_Energy_hamil(:,i)*Energy_hamil(:,i))
enddo



end subroutine

!....................................................................................
subroutine velocity_verlet2
pos=pos+v*dtc+0.5*acc1*dtc*dtc
old_acc2=acc1
call potential(ham_sp)
call force
v=v+0.5*(acc1+old_acc2)*dtc
end subroutine

!...................................................................................
subroutine stochastic_force(delr,delv)
real*8, intent(out) :: delr(total_dimensions),delv(total_dimensions)
integer :: i
real*8 :: rnd1,rnd2,sig_r,sig_v,sig_rv,gdt,gamma_B
gamma_B=2*omega
gdt=gamma_B*dtc

do i=1,total_dimensions

sig_r=dtc*dsqrt(KT/mass(i)*1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
sig_v=dsqrt(KT/mass(i)*(1-dexp(-2*gdt)))
sig_rv=(dtc*KT/mass(i)*1.d0/gdt*(1-dexp(-gdt))**2)/(sig_r*sig_v) !correlation coefficient

call gaussian_random_number(rnd1)
call gaussian_random_number(rnd2)
delr(i)=sig_r*rnd1
delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
enddo

end subroutine stochastic_force
!....................................................................................
subroutine populations
integer :: a,i,j,k,Es
real*8, dimension(Hi,Hi,Ne) :: rho_a
real*8, dimension(Hi,Hi,Ne) :: rho_d
real*8, dimension(Ne) :: LUMO_population
real :: pchem,nex,sH, syspop(Hi)


do a=1,Ne
do i=1,Hi
do j=1,Hi
rho_a(i,j,a)=c(a,i)*conjg(c(a,j))
enddo
enddo
enddo

do k=1,Ne
    rho_d(:,:,k)=matmul(Energy_hamil,(matmul(rho_a(:,:,k),transpose(Energy_hamil))))    
enddo

population=0.d0
do k=1,Ne
   population=population+rho_d(1,1,k) 
enddo

tdata(int(time/dtc),1)=1-population

syspop=0.d0
do j=1,Hi
    syspop(j)=sum(rho_d(j,j,1:Ne))
enddo


do i=1,Hi
    if (i>barrier) then
        j=300+i
        write(j,*) time,syspop(i)
    end if
enddo

call calchemcurr(syspop,barrier,2)

end subroutine
!...................................................................................................................
subroutine gausscalchemcurr(pops,sH,sv)
real, intent(in) :: pops(Hi)
integer, intent(in) :: sv
real*8, intent(in) :: sH
integer :: i,k 
real*8 :: nex,pchem
real*8 :: weight(Hi)



do i=1,int(Hi/2)
    weight(int(Hi/2)-i+1)=w_levels(i)
    weight(int(Hi/2)+i)=w_levels(i)
enddo

!do i=1,Hi
!    write(129,*) i,weight(i)
!enddo


pchem=0.0
do i=2,Hi
    if ((H(i,i)>sH).and.i<Hi) then
        pchem=pchem+(pops(i)*gfac(H(i,i),sH)*weight(i-1))/(H(i+1,i+1)-H(i,i))
    end if
enddo

write(52,*) 'hi'
tdata(int(time/dtc),sv)=pchem*band_width/2
end subroutine 

!..................................................................................................................

subroutine calchemcurr(pops,sH,sv)
real, intent(in) :: pops(Hi)
integer, intent(in) :: sv
real*8, intent(in) :: sH
integer :: i,k 
real*8 :: nex,pchem

pchem=0.0
do i=2,Hi
    if (H(i,i)>sH) then
        pchem=pchem+pops(i)*gfac(H(i,i),sH)
    end if
enddo

tdata(int(time/dtc),sv)=pchem
end subroutine 



!....................................................................................................................
subroutine Mcurrents(state,lb,cs)
integer, intent(in) :: cs,lb
integer, intent(in) :: state(Ne)
integer i,nex,j
logical Alogic
real*8 :: pchem,Dx,lamb
real*8 gf,theta,thetac,dth

gf=0.d0
theta=0.d0
Dx=44.32
lamb=52.91
nex=0
pchem=0.d0


!write(355,*) traj_no,time,pchem,lb
!write(355,*) state
do i=1,Ne
    if (state(i)>lb) then
        nex=1
        thetac=acos(sqrt(H(lb,lb)/H(state(i),state(i))))
        dth=thetac/1000
        
        do j=0,1000
            theta=j*dth
            gf=gf+dth*sin(theta)*exp(-Dx/(lamb*cos(theta)))
        enddo
        pchem=pchem+nex*gf
        if (pchem>0.0) then
            write(355,*) traj_no,time,pchem
            write(355,*) state
        end if

        Alogic=ISNAN(pchem)
!        if (Alogic==.true.) then
         if (pchem>10) then
            write(455,*) traj_no,'trajno',lb,'lb',i,'i',state(i),'i'
            write(455,*) state
            write(455,*) time,'time',gf,'gf',nex,'nex',thetac,'thetac'
            write(455,*) H(lb,lb),'H(lb,lb)',H(state(i),state(i)),'H(state(i),state(i))'
            write(455,*) theta,'theta',sin(theta),'sin(theta)',cos(theta),'cos(theta)'
            write(455,*) pchem, 'pchem'

        end if
    end if
enddo


tdata(int(time/dtc),cs)=pchem

end subroutine
!................................................................................................................








!.................................................................................................................
subroutine modulus(matrix,n,determinant)
 IMPLICIT NONE
     complex*16, DIMENSION(n,n) :: matrix
     INTEGER, INTENT(IN) :: n
     complex*16 :: m, temp
     INTEGER :: i, j, k, l
     LOGICAL :: DetExists = .TRUE.
     complex*16,intent(out) :: determinant
     l = 1
     !Convert to upper triangular form
     
     DO k = 1, n-1
         IF (matrix(k,k) == 0) THEN
             DetExists = .FALSE.
             DO i = k+1, n
                 IF (matrix(i,k) /= 0) THEN
                     DO j = 1, n
                         temp = matrix(i,j)
                         matrix(i,j)= matrix(k,j)
                         matrix(k,j) = temp
                     END DO
                     DetExists = .TRUE.
                     l=-l
                     EXIT
                 ENDIF
             END DO
             IF (DetExists .EQV. .FALSE.) THEN
                 determinant= 0
                 return
             END IF
         ENDIF
         DO j = k+1, n
             m = matrix(j,k)/matrix(k,k)
             DO i = k+1, n
                 matrix(j,i) = matrix(j,i) - m*matrix(k,i)
             END DO
         END DO
     END DO

     !Calculate determinant by finding product of diagonal elements
     determinant= l
     DO i = 1, n
         determinant= determinant* matrix(i,i)
     END DO

END subroutine modulus
!..................................................................................
subroutine transversed_surfaces(ks)
!this subroutine draw the potential surface on which the dynamics occured by
!using the input as lambda
real*8 :: posx
real*8 :: PS
real*8 :: eig_val(Hi)
integer :: ip,ei,tp
integer, intent(in) :: ks
character :: pes



!open(ks,file= pes)
do ip=1,600
  posx=-30.0+0.1*(ip)
  call hamiltonian(posx,eig_val)
  do ei=1,Ne
     PS=U0+eig_val(lambda(ei))
  enddo
  write(ks,*) posx,PS
enddo
!close(ks)


end subroutine
!.....................................................................................
subroutine hamiltonian(xos,Es) !! Extra subroutine to plot energy surfaces
real*8 :: U1
real*8, intent(in) :: xos
real*8 ,intent(out) :: Es(Hi)
real*8 :: eig_ham(Hi,Hi)
real*8 :: ham(Hi,Hi)
integer :: i,j,x,y

    U0=0.5*mass(1)*(omega**2)*(xos)**2
    U1=0.5*mass(1)*(omega**2)*(xos-gh)**2+dG
    ham(1,1)=U1-U0
    do i=2,Hi
        if (i.le.int(Hi/2)) then
                ham(1,i)=sqrt(Band_width*w_levels(i))*Vr/2
                ham(i,1)=sqrt(Band_width*w_levels(i))*Vr/2
        else
                ham(1,i)=sqrt(Band_width*w_levels(i-int(Hi/2)))*Vr/2
                ham(i,1)=sqrt(Band_width*w_levels(i-int(Hi/2)))*Vr/2
        end if
        do j=2,Hi
                if (i.eq.j) then
                        if (i.le.(int(Hi/2)+1)) then
                                ham(i,j)=-(Band_width/2)*(0.5+0.5*knot_x(int(Hi/2)-i+2))
                        else
                                ham(i,j)=(Band_width/2)*(0.5+0.5*knot_x(i-int(Hi/2)-1))
                        end if
                else
                        ham(i,j)=0.0
                end if
        end do
   enddo


!do x=1,Hi
 !do y=1,Hi
 !  write(218,'(f12.8$)')ham(x,y)
 !enddo
 !write(218,*)
!enddo




nsize=Hi
call diag_wrapper(ham,nsize,Es,eig_ham)


end subroutine
!...............................................................................
subroutine construct_state_space(ex_space)
implicit none
integer :: i,j,k,l,m,n
integer, intent(out) :: ex_space(Ne**2+1,Ne)
integer :: ts1(Ne)

k=2
  do i=1,Ne
     do j=1,Hi
      ts1=lambda
          if(.not.(any(j==ts1)).and.(k.le.Ne**2+1)) then
   !     if ((FINDLOC(ts1,j,1).eq.0).and.(k.le.Ne**2+1)) then

             ex_space(1,:)=ts1
             ts1(i)=j
             ex_space(k,:)=ts1
             k=k+1
        end if
      enddo
   enddo

end subroutine

!......................................................................
subroutine position_difference(m,k,l,t)
integer, intent(in) :: m
integer , intent(out) :: k,l,t
integer :: diff(Ne),j

diff=Hi_sp(1,:)-Hi_sp(m,:)

do j=1,Ne
   if (diff(j).ne.0) then
      k=Hi_sp(1,j)
      l=Hi_sp(m,j)
      t=j
    endif
enddo


end subroutine

!.................................................................................
!subroutine construct_Ne_overlap(Ne_Wmat)
!implicit none
!integer :: i,j,k,l,m,n
!integer :: ts1(Ne)
!real*8, intent(out) :: Ne_Wmat(Ne**2+1,Ne**2+1)


!do i=1,Ne**2+1
!  do j=1,Ne**2+1
!      Ne_Wmat(i,j)=W_mat(Hi_sp(i,1),Hi_sp(j,1))
!      do m=2,Ne
!          Ne_Wmat(i,j)=Ne_Wmat(i,j)*W_mat(Hi_sp(i,m),Hi_sp(j,m))
!      end do
!  enddo
!enddo



!end subroutine
!....................................................................................

subroutine construct_Ne_overlap(Ne_Wmat)
implicit none
integer :: i,j,k,l,m,n
integer :: temp(Ne),lock,key,pre(Ne)
real*8, intent(out) :: Ne_Wmat(Ne**2+1,Ne**2+1)
real*8 :: store,prod


do i=1,Ne**2+1
  do j=1,Ne**2+1
      if (i.ne.j) then
      Ne_Wmat(i,j)=W_mat(Hi_sp(i,1),Hi_sp(j,1))
        do m=2,Ne
            Ne_Wmat(i,j)=Ne_Wmat(i,j)*W_mat(Hi_sp(i,m),Hi_sp(j,m))
        end do
      else
         lock=1
         store=0.d0
         pre=0
         do while (lock.le.Ne) 
            call throw_key(j,lock,key)
            if (lock.ne.key) then
                pre=Hi_sp(j,:)
                temp=Hi_sp(j,:)
                temp(lock)=pre(key)
                temp(key)=pre(lock)
                prod=W_mat(pre(1),temp(1))
                do m=2,Ne
                   prod=prod*W_mat(pre(m),temp(m))
                end do
                store=store+prod
                lock=lock+1
            else
                lock=lock+1
            end if
         enddo  
            Ne_Wmat(i,j)=W_mat(Hi_sp(i,1),Hi_sp(j,1))
            do m=2,Ne
               Ne_Wmat(i,j)=Ne_Wmat(i,j)*W_mat(Hi_sp(i,m),Hi_sp(j,m))
            end do
            Ne_Wmat(i,j)=Ne_Wmat(i,j)-store
     end if  
  enddo
enddo

end subroutine
!..................................................................................


!subroutine construct_Ne_overlap(Ne_Wmat)
!implicit none
!integer :: i,j,k,l,m,n
!integer :: temp(Ne),lock,key,pre(Ne)
!real*8, intent(out) :: Ne_Wmat(Ne**2+1,Ne**2+1)
!real*8 :: store,prod


!do i=1,Ne**2+1
!  do j=1,Ne**2+1
!         lock=1
!         store=0.d0
!         pre=0
!         do while (lock.le.Ne) 
!            call throw_key(j,lock,key)
!            if (lock.ne.key) then
!                pre=Hi_sp(j,:)
!                temp=Hi_sp(j,:)
!                temp(lock)=pre(key)
!                temp(key)=pre(lock)
!                prod=W_mat(pre(1),temp(1))
!                do m=2,Ne
!                   prod=prod*W_mat(pre(m),temp(m))
!                end do
!                store=store+prod
!                lock=lock+1
!            else
!                lock=lock+1
!            end if
!         enddo  
!            Ne_Wmat(i,j)=W_mat(Hi_sp(i,1),Hi_sp(j,1))
!            do m=2,Ne
!               Ne_Wmat(i,j)=Ne_Wmat(i,j)*W_mat(Hi_sp(i,m),Hi_sp(j,m))
!            end do
!            Ne_Wmat(i,j)=Ne_Wmat(i,j)-store
!            write(69,*) pos(1),store
!  enddo
!enddo




!end subroutine
!..................................................................................
subroutine throw_key(arr,l,k)
implicit none

integer, intent(in) :: arr,l
integer, intent(out) :: k
integer :: plamb(Ne),xi,xj


plamb=Hi_sp(arr,:)
xj=plamb(l)+1
xi=findloc(plamb,xj,1)

if (xi.ne.0) then
   k=xi
else
   k=l
end if

end subroutine

!....................................................................................
subroutine calculate_iforce(istate,iforce,old_new)
implicit none
integer, intent(in) :: istate(Ne)
real*8, intent(out) :: iforce(total_dimensions)
real*8 :: m_F(total_dimensions,Ne)
integer :: j,k
integer, intent(in) :: old_new


if (old_new.eq.1) then
  do j=1,total_dimensions
     do k=1,Ne
       m_F(j,k)=-sum(Energy_hamil(:,istate(k))*matmul(Gradient_hamil(:,:,j),Energy_hamil(:,istate(k))))
     enddo
  enddo
else
  do j=1,total_dimensions
    do k=1,Ne
      m_F(j,k)=-sum(old_Energy_hamil(:,istate(k))*matmul(old_Gradient_hamil(:,:,j),old_Energy_hamil(:,istate(k))))
    enddo
  enddo
endif




iforce=sum(m_F(:,1:Ne))-mass*omega**2*pos(1)





end subroutine
!......................................................................................
subroutine save_old_state


old_Energy_hamil=Energy_hamil
old_lambda=lambda
old_Gradient_hamil=Gradient_hamil
c_old=c



end subroutine




!........................................................................................
subroutine calculate_sigma(istate,Ajj,new_old)
implicit none
integer, intent(in) :: new_old
integer, intent(in) :: istate(Ne)
real*8, intent(out) :: Ajj
integer, dimension(Ne) ::p
integer :: i,d,r,s
complex*16 :: det_S1,det_S2
complex*16, dimension(Ne,Ne) ::  S1,S2

p=istate

if (new_old.eq.1) then
  do r=1,Ne
     S1(:,r)=c(:,p(r))
  enddo
else
  do r=1,Ne
     S1(:,r)=c_old(:,p(r))
  enddo
endif

call modulus(S1,Ne,det_S1)



Ajj=det_S1*conjg(det_S1)




end subroutine
!.....................................................................................

subroutine net_potential_energy(istate,P_energy)
integer, intent(in) :: istate(Ne)
real*8, intent(out) :: P_energy
integer :: i

P_energy=Energy(istate(1))

do i=1,Ne
   P_energy=P_energy+Energy(istate(i))
enddo




end subroutine

!...................................................................................
subroutine verlet_decoherence
implicit none
integer :: i,j,k,u,t
real*8 :: old_force(Ne**2+1,total_dimensions),delF(Ne**2+1,total_dimensions)
real*8 :: sig(Ne**2+1),acc_dec(Ne**2+1,total_dimensions)
real*8 :: temp_delr(Ne**2+1,total_dimensions),temp_delp(Ne**2+1,total_dimensions)

do i=1,Ne**2+1
   call calculate_iforce(Hi_sp(i,:),old_force(i,:),2)
enddo

do i=1,Ne**2+1
   call calculate_sigma(Hi_sp(i,:),sig(i),2)
enddo

do i=1,Ne**2+1
   delF(i,:)=old_force(i,:)-old_force(1,:)
   acc_dec(i,:)=delF(i,:)*sig(i)/mass
enddo


do i=1,Ne**2+1
    delr(i,:)=delr(i,:)+delp(i,:)/mass*dtc+0.5*acc_dec(i,:)*dtc**2
    delp(i,:)=delp(i,:)+0.5*mass*acc_dec(i,:)*dtc
enddo

do i=1,Ne**2+1
   call calculate_iforce(Hi_sp(i,:),new_force(i,:),1)
enddo

delF=0.d0
do j=1,Ne**2+1
   do k=1,Ne**2+1
      delF(j,:)=delF(j,:)+dabs(Ne_overlap(j,k)**2)*(new_force(k,:)-new_force(1,:))
   enddo
enddo

do i=1,Ne**2+1
   acc_dec(i,:)=delF(i,:)*sig(i)/mass
enddo


do i=1,Ne**2+1
    delp(i,:)=delp(i,:)+0.5*mass*acc_dec(i,:)*dtc
enddo


temp_delr=0.d0;temp_delp=0.d0
  do j=1,Ne**2+1
    do k=1,Ne**2+1
      temp_delr(j,:)=temp_delr(j,:)+dabs(Ne_overlap(k,j)**2)*delr(k,:)
      temp_delp(j,:)=temp_delp(j,:)+dabs(Ne_overlap(k,j)**2)*delp(k,:)
    enddo
  enddo
  delr=temp_delr
  delp=temp_delp





!do i=1,Ne**2+1
!   call position_difference(i,u,t)
!   write(122,*) i,u,t
!enddo



!do i=1,Ne**2+1
!   write(112,*) i,sig(i)
!   write(113,*) i,delF(i,:)
!   write(114,*) i,acc_dec(i,:)
!   write(120,*) i,delr(i,:)
!   write(121,*) i,delp(i,:)
!   write(31,*) i, old_force(i,:)
!enddo






end subroutine
!.....................................................................................

subroutine check_collapse
  implicit none
  real*8 :: V_k(Ne**2+1)
  real*8 :: gama_collapse,gama_reset,rnd
  complex*8 su1,c_dia(Ne,Hi)
  integer n,i,j,p,q,k,x,y
  real*8 :: clamb,cmod,N_ele1,N_ele2


do i=1,Ne**2+1
   call net_potential_energy(Hi_sp(i,:),V_k(i))
enddo


 do n=2,Ne**2+1
   gama_reset=sum((new_force(n,:)-new_force(1,:))*dble(delr(n,:)-delr(1,:)))/2
        !! CAUTION !! !! Assumes delr(n,n,:) is in direction of v(:) !!
   call position_difference(n,p,q,k)
   su1=abs((V_k(1)-V_k(n))*vdij(p,q)*sum((delr(n,:)-delr(1,:))*v))/sum(v*v)
   gama_collapse=gama_reset-2*abs(su1)
   gama_collapse=gama_collapse*dtc
   gama_reset=-gama_reset*dtc
   call random_number(rnd)

 
   if (gama_collapse>rnd) then
     c=0
     do i=1,Ne
  !    c(i,p)=(c(i,p)/abs(sqrt(c(i,p)*conjg(c(i,p)))))*abs(sqrt(c(i,p)*conjg(c(i,p))+c(i,q)*conjg(c(i,q))))
  !    c(i,q)=0
    !  write(72,*) time,sum(c(i,:)*conjg(c(i,:)))
         c(i,lambda(i))=1.d0
      enddo
      write(59,*) time,gama_collapse
    
   endif
           
          
     



         if(rnd<gama_collapse.or.rnd<gama_reset) then
   !         write(91,*) Hi_sp(1,:),'1'
   !         write(91,*) Hi_sp(n,:),'n'
           ! write(98,*) time,k,p,q
     
            write(99,*) time,gama_collapse,gama_reset
            do j=1,Ne**2+1
              delr(n,:)=0.d0;delr(j,:)=0.d0
              delp(n,:)=0.d0;delp(j,:)=0.d0
            enddo
            collapse_flg=1
          endif

enddo



end subroutine
!......................................................................................

subroutine impurity_pop
complex*8 :: rho(Hi,Hi,Ne),rho_d(Hi,Hi,Ne)
real*8, dimension(Ne) :: LUMO_population
integer :: a,i,j,k,Es
real :: syspop(Hi)



do a=1,Ne
   do i=1,Hi
     do j=1,Hi
        if ((i.eq.j).and.(i.eq.lambda(a))) then
           rho(i,j,a)=1
        else if (i.eq.j) then
           rho(i,j,a)=0
        else
           rho(i,j,a)=c(a,i)*conjg(c(a,j))
        end if
     enddo
   enddo
enddo

!do k=1,Ne
! do i=1,Hi
!   do j=1,Hi
!     write(218,'(f12.8$)')real(rho(i,j,k))
!   enddo
!   write(218,*)
! enddo
!  write(218,*)
!enddo



do k=1,Ne
    rho_d(:,:,k)=matmul(Energy_hamil,(matmul(rho(:,:,k),transpose(Energy_hamil))))    
enddo


population=0.d0
do k=1,Ne
   population=population+rho_d(1,1,k) 
enddo

!write(18,*)time,1-real(population)!,lambda
!write(18,*) time,real(population)

syspop=0.d0
do j=1,Hi
    syspop(j)=sum(rho_d(j,j,1:Ne))
enddo

!Es=int(Hi/2)+1

call calchemcurr(syspop,barrier,6)

tdata(int(time/dtc),4)=1-real(population)
end subroutine
!.....................................................................................
subroutine A_mat_pop
complex*16 :: c_dia(Ne,Hi),S1(Ne,Ne)
complex*16 :: det_S1,Aii,Tr_AA
real*8 :: pops,su,rnd3
Integer :: p(Ne),r,di_lambda(Ne),ts1(Ne)
integer :: ex_space(Ne**2+1,Ne),i,j,k
integer :: dia_state(Ne**2+1),in_flg
real*8 :: fullc

p=lambda

c_dia=matmul(c,transpose(Energy_hamil))

!fullc=0.d0
!do i=1,Ne
!   do j=1,Hi
!     fullc=fullc+c_dia(i,j)*conjg(c_dia(i,j))
!   enddo
!   if (fullc<0.98) then
!      write(344,*) time,i,fullc
!   end if
!   fullc=0.d0
!enddo

in_flg=1
if ((in_flg.eq.1).or.(hop_flg.eq.1).or.collapse_flg.eq.1) then
!....construct current diabatic state
  do j=1,Ne
     call random_number(rnd3)
     su=0.d0
    do i=1,Hi
      su=su+c_dia(j,i)*conjg(c_dia(j,i))
      if (rnd3<su) then
         di_lambda(j)=i
         exit 
      end if
    enddo
  enddo
  in_flg=0
end if
!.....................................

!write(900,*) di_lambda


!.....constructing Diabatic state space
k=2
  do i=1,Ne
     do j=1,Hi
      ts1=di_lambda
          if(.not.(any(j==ts1)).and.(k.le.Ne**2+1)) then
             ex_space(1,:)=ts1
             ts1(i)=j
             ex_space(k,:)=ts1
             k=k+1
        end if
      enddo
   enddo
!...........................................

!write(902,*)
!do i=1,Ne**2+1
!   write(902,*) i,ex_space(i,:)
!enddo
!write(902,*)




dia_state=FINDLOC(ex_space,value=1,dim=2)


!write(902,*)
!do i=1,Ne**2+1
!   write(902,*)i,dia_state(i)
!enddo
!write(902,*)


Tr_AA=0
do i=1,Ne**2+1
   call calculate_gen_A(ex_space(i,:),ex_space(i,:),c_dia,Aii)
   Tr_AA=Tr_AA+Aii
enddo
!write(901,*) time,Tr_AA





pops=0.d0
do i=1,Ne**2+1
   if (dia_state(i).ne.0) then
      call calculate_gen_A(ex_space(i,:),ex_space(i,:),c_dia,Aii)
      pops=pops+real(Aii)
      if ((real(Tr_AA)>1.01).or.(real(Tr_AA)<0.98)) then
         write(346,*) time,i,'trA',Tr_AA,'Ai',Aii
         write(346,*) ex_space(i,:)
         write(346,*)
      end if
   end if
enddo

!Tr_AA=0
!do i=1,Ne**2+1
!   call calculate_gen_A(ex_space(i,:),ex_space(i,:),c_dia,Aii)
!   Tr_AA=Tr_AA+Aii
!enddo
!write(901,*) time,Tr_AA




write(12,*) time,1-pops

end subroutine
!......................................................................................
subroutine all_e
integer :: i,j
real*8 :: t_el

t_el=c(1,1)*conjg(c(1,1))
do i=1,Ne
  do j=1,Hi
     t_el=t_el+c(i,j)*conjg(c(i,j))
  enddo
enddo

write(80,*) time,t_el


end subroutine
!....................................................................................
subroutine decor
integer :: i

c=0

do i=1,Ne
   c(i,lambda(i))=0
enddo

write(25,*) time
!write(14,*) c
write(13,*) lambda

end subroutine
!...................................................................................


subroutine potential(gq)
implicit none
integer :: i,x,y
real*8 h1,dh1(total_dimensions)
real*8 coup,Vc(Hi)
integer :: nquant
integer,intent(in) :: gq

if (mod(gq,2).eq.1) then  
  nquant=Hi  
  write(611,*) 'hi'
  H=0.d0
  Gradient_hamil=0.d0

  if (gq==3) then
      H(1,1)=dG
  else
      H(1,1)=0.5*mass(1)*omega**2*((pos(1)-gh)**2-pos(1)**2)+dG
  end if

  dh1(1) = mass(1)*omega**2*(-gh)
  coup = sqrt(tau/(2*3.1416))

  Gradient_hamil(1,1,:)=dh1

 
    write(123,*) int(nquant/2)
  do i=1,int(nquant/2)
    Vc(int(nquant/2)-i+1)=coup*sqrt(band_width*w_levels(i))/2.d0
    Vc(int(nquant/2)+i)=coup*sqrt(band_width*w_levels(i))/2.d0
  enddo


  do i=2,nquant
    H(i,i)=e_metal(i-1)
    H(i,1)=Vc(i-1)
    H(1,i)=Vc(i-1)
  enddo

  nsize=Hi
  !call diag_wrapper(H,nsize,Energy,Energy_hamil)

!   write(76,*) pos(1)
!   write(517,*) time
!   do x=1,Hi
!     do y=1,Hi
!       write(517,'(f12.8$)') H(x,y)
!     enddo
!     write(517,*)
!   enddo
else
  nquant=Hi  

  H=0.d0
  Gradient_hamil=0.d0

  if (gq==4) then 
      H(1,1)=dG
  else
      H(1,1)=0.5*mass(1)*omega**2*((pos(1)-gh)**2-pos(1)**2)+dG
  endif

  dh1(1) = mass(1)*omega**2*(-gh)
  coup = sqrt(tau/(2*3.1416*((Hi-1)/Band_width)))

  Gradient_hamil(1,1,:)=dh1

 



  do i=2,nquant
    H(i,i)=i*(Band_width/(Hi-2))-Band_width*(Hi+2)/(2*Hi-4)
!    H(i,i)=-Band_width/2.d0+(i-2)*Band_width/real(Hi-2)
    H(i,1)=coup
    H(1,i)=coup
  enddo

  nsize=Hi
!  call diag_wrapper(H,nsize,Energy,Energy_hamil)




end if
!   do x=1,Hi
!     do y=1,Hi
!       write(517,'(f12.8$)')Gradient_hamil(x,y,1)
!     enddo
!     write(517,*)
!   enddo


!do x=1,Hi
!  do y=1,Hi
!     write(517,'(f12.8$)') H(x,y)
!  enddo
!     write(517,*)
!enddo

end subroutine

!.......................................................................
subroutine draw_pes
integer :: i,j,k,l
real :: dx,TE

do i=1,Ne
   lambda(i)=i
enddo

call construct_state_space(Hi_sp)

do i=1,Ne**2+1
   write(27,*) i
   write(27,*) Hi_sp(i,:)
   write(27,*)
enddo



dx=60/1000.d0
do k=1,2*Ne
 pos=-40
 do i=1,1000
   pos=-40+i*dx
   call potential(1)
   TE=Energy(Hi_sp(k,1))
   do j=2,Ne
      TE=TE+Energy(Hi_sp(k,j))
   enddo
   l=100+k
   write(l,*) pos, 0.5*mass(1)*omega**2*pos**2+TE
 enddo
enddo

!k=0

!dx=60/1000.d0
!pos=-20
!do i=1,1000
!   pos=-20+i*dx
!   call potential(1)
!   do j=1,Ne
!      k=100+j
!      write(k,*) pos,Energy(lambda(j))
!   enddo
!enddo




end subroutine
!..................................................................
subroutine full_density_matrix
implicit none
integer :: i,j,kx,ky
complex*16 :: Aij(Ne**2+1,Ne**2+1),Ad_ij(Ne**2+1,Ne**2+1)
real*8 :: Si(Ne,Ne),Sj(Ne,Ne)
real*8 :: det_Si,det_Sj
complex*16 :: Tr_Ax
do i=1,Ne**2+1
   do j=1,Ne**2+1
       call calculate_As(Hi_sp(i,:),Hi_sp(j,:),Aij(i,j))
   end do
end do


Aij(1,1)=1

do i=2,Ne**2+1
   Aij(i,i)=0
enddo



Ad_ij=0

do kx=1,Ne**2+1
  do ky=1,Ne**2+1
    do i=1,Ne
      do j=1,Ne
        Si(i,j)=Energy_hamil(Hi_sp(kx,i),Hi_sp(kx,j))
      enddo
    enddo
    
    
    det_Si=FindDet(Si,Ne)
    
    do i=1,Ne
      do j=1,Ne
        Sj(i,j)=Energy_hamil(Hi_sp(ky,i),Hi_sp(ky,j))  
      enddo
    enddo
  
    Sj=transpose(Sj)

    det_Sj=FindDet(Sj,Ne)

   
    do i=1,Ne**2+1 
      do j=1,Ne**2+1
        Ad_ij(kx,ky)=Ad_ij(kx,ky)+det_Si*Aij(i,j)*det_Sj  
      end do
    end do
  enddo
enddo


write(101,*) time,1.d0-real(Ad_ij(1,1))

population=Ad_ij(1,1)

!pop_mat(int(time/dtc))=1-real(population)


Tr_Ax=0
do i=1,Ne**2+1
   Tr_Ax=Tr_Ax+Ad_ij(i,i)
enddo


write(887,*) time,Tr_Ax

end subroutine
!....................................................................
subroutine calculate_As(istate,jstate,Aij)
implicit none
integer, intent(in) :: istate(Ne),jstate(Ne)
complex*16, intent(out) :: Aij
complex*16 :: Ai(Ne,Ne), Aj(Ne,Ne)
complex*16 :: det_Ai,det_Aj
integer :: i,j

write(502,*) 'a'

do i=1,Ne
   Ai(:,i)=c(:,istate(i))
enddo

write(502,*) 'b'

call modulus(Ai,Ne,det_Ai)

write(502,*) 'c',det_Ai

do j=1,Ne
   Aj(:,j)=c(:,jstate(j))
enddo

write(502,*) 'd'

call modulus(Aj,Ne,det_Aj)

write(502,*) 'e'

Aij=det_Ai*conjg(det_Aj)

write(502,*) 'f',det_Aj

end subroutine
!....................................................................
subroutine check_mod
implicit none
real*8 :: test_mat(2,2),det_test

test_mat(1,1)=1.0
test_mat(1,2)=2.0
test_mat(2,1)=3.0
test_mat(2,2)=4.0

!call modulus(test_mat,2,det_test)

write(105,*) FindDet(test_mat,2)

end subroutine

!.......................................................................

REAL*8 FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    REAL*8, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL*8 :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
   
END FUNCTION FindDet
!..........................................................................
subroutine full_correct_density_matrix
implicit none
integer :: i,j,kx,ky,k
complex*16 :: Aij(Ne**2+1,Ne**2+1),Ad_ij(Ne**2+1,Ne**2+1)
real*8 :: Si(Ne,Ne),Sj(Ne,Ne)
real*8 :: det_Si,det_Sj,rnd4,su
real*8 :: U1_mat(Ne**2+1,Ne**2+1)
complex*16 :: Mmat(Ne,Hi),Tr_AA,Tr_Ai
integer :: diab_lamb(Ne),ts1(Ne)
integer :: ex_space(Ne**2+1,Ne)
integer :: dia_state(Ne**2+1)
real*8 :: pops


do kx=1,Ne**2+1
  do ky=1,Ne**2+1
  
    do i=1,Ne
      do j=1,Ne
        Si(i,j)=Energy_hamil(Hi_sp(kx,i),Hi_sp(kx,j))
      enddo
    enddo
    
    
    det_Si=FindDet(Si,Ne)
    
    do i=1,Ne
      do j=1,Ne
        Sj(i,j)=Energy_hamil(Hi_sp(ky,i),Hi_sp(ky,j))  
      enddo
    enddo
  
    Sj=transpose(Sj)

    det_Sj=FindDet(Sj,Ne)

    U1_mat(kx,ky)=det_Si*det_Sj   

   enddo
enddo

do i=1,Ne**2+1
   do j=1,Ne**2+1
       call calculate_As(Hi_sp(i,:),Hi_sp(j,:),Aij(i,j))
   end do
end do

Tr_Ai=0
do i=1,Ne**2+1
   Tr_Ai=Tr_Ai+Aij(i,i)
enddo

write(885,*) Tr_Ai


Aij(1,1)=1

do i=2,Ne**2+1
   Aij(i,i)=0
enddo

Ad_ij=matmul(U1_mat,transpose(Aij))


write(39,*) time,1.d0-real(Ad_ij(1,1))


Mmat=0.d0
do i=1,Ne
   Mmat(i,lambda(i))=1.d0
enddo

Mmat=matmul(Mmat,transpose(Energy_hamil))


do j=1,Ne
   call random_number(rnd4)
   su=0.d0
   do i=1,Hi
     su=su+Mmat(j,i)*conjg(Mmat(j,i))
     if (rnd4<su) then
        diab_lamb(j)=i
        exit 
     end if
   enddo
enddo



write(720,*) time,diab_lamb


!.....constructing Diabatic state space
k=2
  do i=1,Ne
     do j=1,Hi
      ts1=diab_lamb
          if(.not.(any(j==ts1)).and.(k.le.Ne**2+1)) then
             ex_space(1,:)=ts1
             ts1(i)=j
             ex_space(k,:)=ts1
             k=k+1
        end if
      enddo
   enddo
!...........................................
!write(721,*)
!do i=1,Ne**2+1
!   write(721,*) i,ex_space(i,:)
!enddo
!write(721,*)

!write(722,*)time,FINDLOC(ex_space,value=1,dim=2)

dia_state=FINDLOC(ex_space,value=1,dim=2)


pops=0.d0
do i=1,Ne**2+1
   if (dia_state(i).ne.0) then
      pops=pops+real(Ad_ij(i,i))
 !     write(774,*) time,i
 !     write(775,*) time, dia_state 
   end if
enddo


Tr_AA=0
do i=1,Ne**2+1
   Tr_AA=Tr_AA+Ad_ij(i,i)
enddo
write(884,*) time,Tr_AA

write(783,*) time
do i=1,Ne**2+1
   write(783,*)i,Ad_ij(i,i)
enddo
write(783,*)

write(19,*) time,1.d0-pops


end subroutine
!........................................................................

subroutine calculate_gen_A(istate,jstate,wave_fn,Aij)
implicit none
integer, intent(in) :: istate(Ne),jstate(Ne)
complex*16, intent(in) :: wave_fn(Ne,Hi)
complex*16, intent(out) :: Aij
complex*16 :: Ai(Ne,Ne), Aj(Ne,Ne)
complex*16 :: det_Ai,det_Aj
integer :: i,j

do i=1,Ne
   Ai(:,i)=wave_fn(:,istate(i))
enddo


call modulus(Ai,Ne,det_Ai)



do j=1,Ne
   Aj(:,j)=wave_fn(:,jstate(j))
enddo


call modulus(Aj,Ne,det_Aj)


Aij=det_Ai*conjg(det_Aj)


end subroutine


!.........................................................
subroutine full_rho_diab
implicit none
integer :: k,i,j
real*8 :: Sd(Ne,Ne)
complex*16 :: Ad_ij(Ne,Ne)
real*8 :: Tr_A,Aij(Ne,Ne)

write(500,*) 'M'

Ad_ij=0
Aij=0

write(500,*) 'N'

do i=1,Ne**2+1
   write(503,*) 'i',i
   write(503,*) Hi_sp(i,:)
enddo 


do i=1,Ne**2+1
   call calculate_sigma(Hi_sp(i,:),Aij(i,i),1)
enddo

write(505,*) time
do i=1,Ne**2+1
  write(505,*)i,Aij(i,i)
enddo

do k=1,Ne**2+1
  do i=1,Ne
    do j=1,Ne
      Sd(i,j)=Energy_hamil(Hi_sp(k,i),Hi_sp(k,j))
    enddo
  enddo
  write(500,*) 'q' 
  do i=1,Ne**2+1
     Ad_ij(k,k)=Ad_ij(k,k)+((FindDet(Sd,Ne))**2)*Aij(i,i)
  enddo
  write(500,*) 'r'
end do

write(500,*) 's'
Tr_A=0
do i=1,Ne*2+1
   Tr_A=Tr_A+real(Ad_ij(i,i))
enddo

write(500,*) 't'
write(999,*) time,Tr_A


end subroutine
!........................................................
subroutine changing_A_basis(siga,sigd)
real*8, intent(in) :: siga(Ne**2+1)
real*8 , intent(out) :: sigd(Ne**2+1)
integer :: i,j,k
real*8 :: Sd(Ne,Ne),Tr_A


do k=1,Ne**2+1
  do i=1,Ne
    do j=1,Ne
      Sd(i,j)=Energy_hamil(Hi_sp(k,i),Hi_sp(k,j))      
    enddo
  enddo
  write(909,*) time,k,FindDet(Sd,Ne)
  do i=1,Ne**2+1
     sigd(k)=sigd(k)+siga(i)*FindDet(Sd,Ne)**2
     write(998,*) k,i,siga(i),FindDet(Sd,Ne)**2
  enddo
end do

do j=1,Ne**2+1
   write(919,*) j,siga(j)*FindDet(Sd,Ne)**2
enddo


Tr_A=0
do i=1,Ne*2+1
   Tr_A=Tr_A+real(sigd(i))
enddo


write(999,*) time,Tr_A






end subroutine
!................................................
subroutine reset_flags

hop_flg=0

collapse_flg=0

end subroutine
!................................................
subroutine check_long_pop
implicit none
real*8 :: bar_loc,pop_l,pew
!bar_loc=6.864

pew=mass(1)*(omega**2)*gh
bar_loc=(gh/2)+dG/pew

if (pos(1)<bar_loc) then
   pop_l=0.d0
else
   pop_l=1.d0
end if

write(14,*) time,1.d0-pop_l
tdata(int(time/dtc),5)=1-real(pop_l)
end subroutine
!............................................................
subroutine pop_averaging(matr_pops,fil_no)!popl_m,fil_no)
implicit none
real*8, intent(in), dimension(int(total_time/dtc),ntraj) ::  matr_pops
real*8 :: popl_m(int(total_time/dtc))
integer :: i,j
integer, intent(in) :: fil_no


popl_m=0
do i=1,ntraj
   popl_m=popl_m+matr_pops(:,i)
enddo

do j=1,int(total_time/dtc)
       write(fil_no,*) j*dtc,popl_m(j)/real(ntraj)
enddo


end subroutine
!................................................
subroutine statepop
implicit none
real*8 :: dj,pop
integer :: i,loc,j
real*8 ::o1i(Hi)

do i=1,Hi
   o1i(i)=abs(Energy_hamil(1,i))
enddo
dj=maxval(o1i)
loc=maxloc(o1i,1)

if ((abs(dj)>0.7).and.(any(loc==lambda))) then
    pop=1.d0
!    write(76,*) dj,'dj',loc,'loc',lambda
    write(77,*) time,1.d0-pop,fermi(pos(1))!, 'loc', loc, 'lambda', lambda
    !write(77,*)
else
    pop=0.d0
    write(77,*) time,1.d0-pop,fermi(pos(1))!, 'loc', loc, 'lambda', lambda
    !write(77,*)
end if


end subroutine
!.......................................................................
  Function fermi(x)
!----------------------------------------
!Fermi Function 
!----------------------------------------
implicit none
real*8 :: fermi, x
real*8 :: Q0,Q1,dE

Q0=0.5*mass(1)*omega**2*x**2
Q1=0.5*mass(1)*omega**2*(x-gh)**2+dG
dE=Q1-Q0


fermi=1/(1+exp(dE/KT))

return
end

!................................................................
function gfac(epi,eps)
real*8 :: epi,eps, gfac
real*8 :: theta, thetac, Dx, lamb,dth
integer i

thetac=acos(eps/epi)
gfac=0.d0
theta=0.d0
dth=thetac/1000
Dx=44.32
lamb=52.91


do i=0,1000
    theta=i*dth
    gfac=gfac+dth*sin(theta)*exp(-Dx/(lamb*cos(theta)))
enddo


end
!...............................................................................
subroutine two_state_dynamics(outst)

integer, intent(inout) :: outst(Hi)
integer :: i,j
integer :: alamb(Ne)
real*8 :: exo,pr,rndt
integer :: temp(Hi)
real*8 :: kb,kf,kt,pops

if (outst(1)==1) then
    alamb=numstate(outst,0)
else
    alamb=numstate(outst,1)
end if

pr=0.d0
iloop : do i=1,Ne
! if (alamb(i)==6) then
    temp=outst
    if (temp(1)==1) then
        temp(1)=0 
        temp(alamb(i))=1
    else
        temp(1)=1 
        temp(alamb(i))=0
    end if
    exo=diab_energy(numstate(temp,1))-diab_energy(numstate(outst,1))
!    write(888,*) time
!    write(888,*)  numstate(temp,1)
!    write(888,*) numstate(outst,1)
    !if (outst(1)==0) then
    !    pr=pr+(1-fermi(exo))*(Mrate(exo,H(1,alamb(i))))*dtc
    !else
    !    pr=pr+fermi(exo)*(Mrate(exo,H(1,alamb(i))))*dtc
!    end if
!    pr=pr+(Mrate(exo,H(1,alamb(i))+Mrate(-exo,H(1,alamb(i)))))*dtc
    pr=pr+(Mrate(exo,H(1,alamb(i))))*dtc

    call random_number(rndt)
    !if (time<50000) write(666,*) time,pr   
    if (pr>rndt) then
        write(4001,*) traj_no,time
        write(4001,*) exo,H(1,alamb(i))
        outst=temp
        write(4001,*) outst
        exit iloop
    end if
!  endif
enddo iloop


end subroutine 
!...............................................................................
subroutine test_dynamics
implicit none
real*8:: kf,kb,kt,pops,exo,tim,exo2
integer :: temp(Hi),po

po=6
tim=0.d0
temp=binstate(diabstate)
temp(1)=1
temp(po)=0

write(565,*) binstate(diabstate)
write(565,*)
write(565,*)temp

exo=diab_energy(diabstate)-diab_energy(numstate(temp,1))
exo2=H(1,1)-H(po,po)
write(21,*) exo,exo2
do while (tim.le.total_time)
    kf=Mrate(exo,H(1,po))
    kb=Mrate(-exo,H(1,po))
    kt=kf+kb
    write(511,*) exo,H(1,po)
    pops=exp(-(kt*tim-log(kf)))/kt+(kb/kt)
    write(199,*) tim,pops
    tim=tim+dtc
enddo

end subroutine
!...............................................................................
subroutine simple_dynamics(outst,po)
integer, intent(inout) :: outst(Hi)
integer, intent(in) :: po
integer :: temp(Hi)
real*8 :: exo,pr,rndt

temp=outst
if (temp(1)==1) then
    temp(1)=0 
    temp(po)=1
else
    temp(1)=1 
    temp(po)=0
end if


exo=diab_energy(numstate(temp,1))-diab_energy(numstate(outst,1))

pr=pr+(Mrate(exo,H(1,po)))*dtc
call random_number(rndt)
if (pr>rndt) then
    write(4001,*) traj_no
    write(4001,*) time,pr,rndt
    write(4001,*) numstate(outst,1)
    outst=temp
end if

end subroutine 
!.............................................................................
subroutine state_dynamics
integer :: outstate(Hi),i,dims
integer :: Es,occu,LUMO,j,n

time=0.d0
Es=int(Hi/2)+1

do while (time.le.total_time)
    call hop_state(diabstate)
    outstate=binstate(diabstate)
    call check_occupied(diabstate,occu)
    if (occu==0) then
        LUMO=0
    else
        LUMO=1
    end if
    
    barrier=0.d0 !careful
    if ((jflag==0).or.fflag==0) then
!        if ((outputs(7)==1).and.(ham_sp==4)) call calchemcurr(real(outstate),barrier,7)
!        if ((outputs(7)==1).and.(ham_sp==3)) call gausscalchemcurr(real(outstate),barrier,7)

        if (outputs(7)==1) call calchemcurr(real(outstate),barrier,7)
        if (outputs(3)==1) tdata(int(time/dtc),3)=LUMO
    end if

    i=ceiling(time/pop_avg_time)
    out_dist(i,:)=outstate
    
    time=time+dtc

enddo

if (fflag==1) then
    if (jflag==1) then
        if (outstate(1)==0) then
             bins=bins+outstate
             neu=neu+1
             iflag=1
        else
             iflag=0
        end if
    end if
endif



end subroutine 
!................................................................................
subroutine Marcus(k01,k10,V_exo,coup)
real*8,intent(out) :: k01,k10
real*8, intent(in) :: V_exo,coup
real*8 longpop,V_reorg
real*8, parameter:: pi=22/7

V_reorg=0.5*mass(1)*omega**2*gh**2

k01=(2*pi*coup**2/sqrt(4*pi*V_reorg*kT))*exp(-(V_exo-V_reorg)**2/(4*V_reorg*KT))


k10=(2*pi*coup**2/sqrt(4*pi*V_reorg*kT))*exp(-(V_exo-V_reorg)**2/(4*V_reorg*KT))

longpop=exp(-V_exo/KT)

write(*,*) k01/k10,'marcus',1/longpop
write(*,*) k01,k10
end subroutine 
!..............................................................................
function Mrate(V_exo,coup)
real*8 :: pi
real*8 :: Mrate,V_exo,coup,V_reorg

pi=22/7

V_reorg=0.5*mass(1)*omega**2*gh**2
!V_reorg=0.00125

!Mrate=((2*pi*coup**2)/sqrt(4*pi*V_reorg*KT))*exp(-(V_exo+V_reorg)**2/(4*V_reorg*KT))
Mrate=2*pi*coup**2 * 1.d0/sqrt(4*pi*V_reorg*KT) * exp(-(V_reorg+V_exo)**2/(4*V_reorg*KT))

!Mrate=(2*pi*coup**2/sqrt(4*pi*V_reorg*KT))*exp(-(V_reorg+V_exo)**2/(4*V_reorg*KT))



end function
!..............................................................................
function diab_energy(istate)
real*8 :: diab_energy
integer :: i,istate(Ne)


diab_energy=0.d0
do i=1,Ne
    diab_energy=diab_energy+H(istate(i),istate(i))
enddo
    
end
!.............................................................................
function binstate(nstate)
integer :: nstate(Ne),binstate(Hi)
integer :: i


binstate=0
do i=1,Ne
    binstate(nstate(i))=1
enddo

end function
!......................................................................
function numstate(bstate,jflag)
integer :: numstate(Ne),bstate(Hi),jflag,i,k

k=1
do i=1,Hi
    if (bstate(i)==jflag) then
        numstate(k)=i
        k=k+1
    end if
enddo 


end function    
!.........................................................................

subroutine hop_state(state)
  implicit none
  integer i,state_new(Ne),occupied
  real*8 coup,exo,k_Marc,rnd,su,prob(Hi)
  real*8 lambda_B,pi,hbar
  integer, intent(inout) :: state(Ne)
    
  pi=22/7
  hbar=1.d0
  call check_occupied(state,occupied)
!  write(112,*) time,i
!  write(112,*) occupied
  lambda_B=0.5*mass(1)*omega**2*gh**2
  if(occupied==0) then !! Neutral
    su=0.d0
    call random_number(rnd)
    outerloop: do i=1,Ne
      state_new=state
      state_new(i)=1
      exo=diab_energy(state_new)-diab_energy(state)
      coup=H(1,state(i))
!      write(114,*) i,state(i),'part1'
      k_Marc=2*pi*coup**2/hbar * 1.d0/sqrt(4*pi*lambda_B*KT) * exp(-(lambda_B+exo)**2/(4*lambda_B*KT))
!      write(57,*) k_Marc
      su=su+k_Marc*dtc
      if(rnd<su) then
        state=state_new
        exit outerloop
      endif
    enddo outerloop
  else    !! Anionic
    su=0.d0
    call random_number(rnd)
    outerloop2: do i=1,Hi
      if(.not.(any(i==state))) then
        state_new=state
        state_new(occupied)=i
        exo=diab_energy(state_new)-diab_energy(state)
        coup=H(1,i)

!       write(114,*) i,state(i),'part2'
        k_Marc=2*pi*coup**2/hbar * 1.d0/sqrt(4*pi*lambda_B*KT) * exp(-(lambda_B+exo)**2/(4*lambda_B*KT))
        su=su+k_Marc*dtc
!        write(57,*) k_Marc
        if(rnd<su) then
          state=state_new
          exit outerloop2
        endif
      endif
    enddo outerloop2
  endif


end subroutine hop_state
!-----------------------------------------------------------------  

subroutine check_occupied(istate,occupied)
  implicit none
  integer, intent(in) :: istate(Ne)
  integer,intent(out)::occupied
  integer i

  occupied=0 !! Neutral state
  do i=1,Ne
    if(istate(i)==1) occupied=i !! Anionic
  enddo

end subroutine check_occupied
!.................................................................
function compliment(state)
integer state(Ne),compliment(Ne),i,j

j=1
do i=1,Hi
    if(.not.(any(i==state))) then
        compliment(j)=i
        j=j+1
    end if
enddo

end function 
!.............................................................
subroutine k_hopping(istate)
integer, intent(inout) :: istate(Ne)
integer alamb(Ne)
integer temp(Ne)
integer i,occupied,j
real*8 exo,pr,coup,rndr,V_k(Hi)

call check_occupied(istate,occupied)

if (occupied==0) alamb=istate
if (occupied.ne.0) alamb=compliment(istate)

!write(43,*) time
!write(43,*) istate
!write(43,*) compliment(istate)

pr=0.d0
iloop: do i=1,Ne

    if (occupied==0) then 
        coup=H(1,i)
        exo=H(1,1)-H(alamb(i),alamb(i))
    else
        coup=H(1,alamb(i))
        exo=H(alamb(i),alamb(i))-H(1,1)
    end if

    pr=pr+Mrate(exo,coup)*dtc
    call random_number(rndr)
    if (pr>rndr) then
        if (occupied==0) istate(i)=1
        if (occupied.ne.0) istate(occupied)=alamb(i)
        exit iloop
    end if
enddo iloop

end subroutine


!.............................................................
function scompliment(state)
integer state(Ne),scompliment(Hi-Ne),i,j

j=1
do i=1,Hi
    if(.not.(any(i==state))) then
        scompliment(j)=i
        j=j+1
    end if
enddo

end function 
!..............................................................


subroutine M_hopping(istate)
integer, intent(inout) :: istate(Ne)
integer alamb(Hi-Ne)
integer temp(Ne)
integer i,occupied,j
real*8 exo,pr,coup,rndr,V_k(Hi),Er

call check_occupied(istate,occupied)

if (occupied.ne.0) alamb=scompliment(istate)

!write(67,*) time
!rite(67,*) scompliment(istate)

!do j=1,Hi
!    write(122,*) H(1,j)
!enddo
!write(56,*) istate
pr=0.d0
iloop: do i=1,Ne
!   if (istate(i)==6) then
    if (occupied==0) then 
        coup=H(1,i)
        exo=H(1,1)-H(istate(i),istate(i))
    else
        coup=H(1,alamb(i))
        exo=H(alamb(i),alamb(i))-H(1,1)
    end if
    pr=pr+Mrate(exo,coup)*dtc
    call random_number(rndr)
   ! write(666,*) time,i,exo,coup,pr,occupied,rndr
   ! if (pr>rndr) then
   !     write(666,*) pr,'pr',rndr,'rndr'
   ! end if

    if (pr>rndr) then
        
!        write(53,*) traj_no,time
!        write(53,*) istate
        if (occupied==0) then
            istate(i)=1
            exit iloop
        else
            istate(occupied)=alamb(i)
            exit iloop
        end if
     end if
enddo iloop
!write(26,*) time,pr
end subroutine












end module














