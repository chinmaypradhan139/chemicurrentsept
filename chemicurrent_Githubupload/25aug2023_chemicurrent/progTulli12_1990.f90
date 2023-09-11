program tulli12
use mod_global
implicit none
integer :: i,j,n,total_trajectories,number_of_cores,TT,di,k,l
real*8 :: start,finish
real*8, dimension(:,:,:), allocatable :: store


open(1, file ='input.txt')
read(1,*) number_of_cores,total_trajectories,iseed
close(1)

call setup_initial_values1
ntraj=int(total_trajectories/number_of_cores)


TT=int(total_time/dtc)

allocate(store(TT,ntraj,outdims))

avg_dist=0.d0

do traj_no=1,ntraj
    call setup_initial_values2

    !call draw_pes
    call CPU_TIME(start)

    call classical_evolution 
    write(111,*) traj_no


    do  k=1,outdims
        if (outputs(k)==1) store(:,traj_no,k)=tdata(:,k)
    enddo
    
    do n=1,2
        do i=1,out_freq
            avg_dist(i,:,n)=avg_dist(i,:,n)+out_dist(i,:,n)
        enddo
    enddo

    call CPU_time(finish)
    write(119,*)finish-start
enddo


!do i=1,out_freq
    !write(432,*) "#", i*pop_avg_time
!  do n=1,2
!    l=432+n    
!    do j=2,Hi
!        write(l,*) H(j,j),avg_dist(i,j,n)/ntraj-fdist(H(j,j))
!    enddo
!  enddo
!    write(432,*)
!enddo


k=0
do di=1,outdims
    k=100+di 
    if (outputs(di)==1) call pop_averaging(store(:,:,di),k)!
enddo



open(345, file="ended",status='new')
    write(345,*) "ended"
close(345)

end program
!..............................................................................


