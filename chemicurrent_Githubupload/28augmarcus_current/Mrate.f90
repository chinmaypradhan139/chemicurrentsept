program Mrate
real*8 :: coup,V_reorg,KT,exo(5),rate(5),pi,rate2(5)
integer i

coup=2.659d-4
KT=0.00095
V_reorg=0.00126
exo(1)=1.62d-02
exo(2)=1.12d-02
exo(3)=6.2d-03
exo(4)=1.2d-03
exo(5)=-0.0038

rate=0.d0
rate2=0.d0
pi=22/7
do i=1,5
    rate(i)=2*pi*coup**2 * 1.d0/dsqrt(4*pi*V_reorg*KT) * dexp(-(V_reorg+exo(i))**2/(4*V_reorg*KT))
    rate2(i)=2*pi*coup**2 * 1.d0/sqrt(4*pi*V_reorg*KT) * exp(-(V_reorg+exo(i))**2/(4*V_reorg*KT))
enddo


write(*,*) rate
write(*,*) rate2



end program Mrate
