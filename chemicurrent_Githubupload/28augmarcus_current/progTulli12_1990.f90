program tulli12
use mod_global
implicit none
integer :: i,j,n,total_trajectories,number_of_cores,TT,di,k
real*8 :: start,finish,nor
real*8, dimension(:,:,:), allocatable :: store

open(1, file ='input.txt')
read(1,*) number_of_cores,total_trajectories,iseed
close(1)

call setup_initial_values1
ntraj=int(total_trajectories/number_of_cores)


TT=int(total_time/dtc)

allocate(store(TT,ntraj,outdims))
neu=0
traj_no=1
!do traj_no=1,ntraj
avg_dist=0.d0
do while(traj_no.le.ntraj)
    call setup_initial_values2
    call CPU_TIME(start)
    fflag=0
    !call Metropolis_Hastings 
    call simple_fermi


    if (fflag==1) then
        jflag=1
!        call state_dynamics
        jflag=0
        if (iflag==0) cycle
    end if
    
    if (outputs(3)==1) call state_dynamics 

    do  k=1,outdims
        if (outputs(k)==1) store(:,traj_no,k)=tdata(:,k)
    enddo

    do i=1,out_freq
        avg_dist(i,:)=avg_dist(i,:)+out_dist(i,:)
    enddo
    
    

    call CPU_time(finish)
    write(119,*)finish-start
  !  write(689,*) ntraj,traj_no,neu
    traj_no=traj_no+1

enddo

do i=1,out_freq
    !write(432,*) "#", i*pop_avg_time
    do j=2,Hi
        write(432,*) H(j,j),avg_dist(i,j)/ntraj-fdist(H(j,j))
    enddo
!    write(432,*)
enddo

if (fflag==1) then
    do i=2,Hi
        nor=nor+fdist(H(i,i))
    enddo

    do i=2,Hi
        write(67,*) H(i,i),bins(i)/real(ntraj),fdist(H(i,i))
    enddo

end if


do i=2,Hi
    write(67,*) H(i,i),fdist(H(i,i))
enddo

!write(412,*) ntraj,neu

!if (fflag==1) ntraj=neu
k=0
do di=1,outdims
    k=100+di 
    if (outputs(di)==1) call pop_averaging(store(:,:,di),k)!
enddo

!call test_dynamics

open(345, file="ended",status='new')
    write(345,*) "ended"
close(345)

end program
!..............................................................................


