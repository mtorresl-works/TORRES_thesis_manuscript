!! program for Monte Carlo (annealing) to determine twist angle of a double-twisted membrane


program siman


implicit none
 
integer :: i,ii,j,oo,icycle,nr,hnr,im1,nsample,nrej
real(8), dimension (:), allocatable ::  rr,pp,ppt
real(8) :: delr,ll,qq,ran,dee,ee,ee0,temp,temp0,disp,pold,hpi,dpp,acra
real(8) :: radius,xx,k3,lamt
real(8), parameter :: pi = 3.14159265358979323846264338327950288419
character(2) :: sys
character(60) :: dummy

! read in parameters

open (21,file = "./inp",action= "read")

read(21,*) nr,dummy ! # grid points
read(21,*) xx,dummy ! rod aspect ratio
read(21,*) lamt,dummy ! twist penetration length in units rod diameter
read(21,*) qq,dummy ! chiral strength via cholesteric pitch in units 1/rod diameter
read(21,*) radius,dummy ! total membrane radius in units rod diameter
read(21,*) temp0,dummy ! initial temperature

close(21)

! open log file

open (23,file = "./log",action= "write", status = "replace")
close(23)


! open  profile file

open (26,file = "./prof",action= "write", status = "replace")
close(26)


allocate(rr(nr),pp(nr),ppt(nr))

hpi = 0.5D0*pi



! elastic constants (bend-twist elastic ratio)

k3 = 0.25D0




! define r array

delr = radius/(nr*1.D0)

rr = (/ (i * delr, i=1,nr) /) 


! define linear starting array for pp ~ q*rr 

qq = 0.01D0

pp = (/ (i * qq * delr, i=1,nr) /) 


! initialize stepsize

disp = 0.2D0


! initialize energy and temperature

call energy(rr,pp,ee0)


CALL RANDOM_SEED()


! SEQUENCE OF MC CYCLES 


do  icycle = 1, 1000000


! quenching procedure

temp  = temp0*exp(-icycle*1.D0/50000.D0)


! monitor angular profile after nr MC cycles
 

if (  (mod(icycle,500) .eq. 0.D0)) then

open (22,file = "./2test",action= "write", status = "replace")

write(22,*) '# temperature : ', temp 


do i=1,nr

write(22,*) rr(i)," ",pp(i)


enddo

close(22)


open (23,file = "./log",action= "write", status = "old", position = "append")


write(23,*) icycle," ",temp," ",ee0


close(23)


endif




nrej = 0


! START MONTE CARLO CYCLE

do i=1,nr


! select random "particle"

CALL RANDOM_NUMBER(ran)

oo = int(nr*ran) + 1


! random displacement 

CALL RANDOM_NUMBER(ran)

pold = pp(oo)

pp(oo) = pp(oo) + (ran - 0.5D0)*disp


! calculate new energy

call energy(rr,pp,ee)


dee = ee - ee0


if (dee .le. 0.D0) then

! accept move

ee0 = ee

else

! accept with Boltzmann probability with temperature


CALL RANDOM_NUMBER(ran)

if (ran .le. exp(-dee/temp)) then

! accept

ee0  = ee

else

! reject

pp(oo)  = pold

nrej = nrej + 1

endif 

endif 

! END OF MC CYCLE


enddo

! END OF SAMPLING



enddo


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! calculate free energy Eq. 5.18

subroutine energy(rrr,ppp,ene)

implicit none
 
integer :: ii,jj,isel
real(8) :: sum,free,ndn,dp,sinp,cosp,hh,dhh
real(8), dimension(nr), intent(in) :: rrr,ppp
real(8), intent(out) :: ene
real(8), parameter :: pi = 3.14159265358979323846264338327950288419

sum = 0.D0

do ii=1,nr-1

sinp=sin(ppp(ii))
cosp=cos(ppp(ii))

! define gradient

dp=(ppp(ii+1) - ppp(ii))/delr


! local membrane height

hh = (xx/2.D0)*cos(ppp(ii))

dhh =(xx/2.D0)*(cos(ppp(ii+1)) - cos(ppp(ii)))/delr

ndn = dp + (1.D0/rrr(ii))*sinp*cosp

! twist-bend elasticity plus depletion

free = 2.D0*pi*rrr(ii)*delr*((ndn+ qq)**2 + k3*((sinp**4)/(rrr(ii)*rrr(ii))) + (1.D0/(lamt**2))*sqrt(1.D0 + (dhh**2)) )


! one-constant approximation (Kang, Lubensky)

! q0 = -kt/k2

! sfel = 2.D0*pi*delr*k2*(rrr(ii)*dp*dp+sin(2*ppp(ii))*dp+(1.D0/rrr(ii))*(sinp**2)-2.D0*q0*rrr(ii)*dp-q0*sin(2.D0*ppp(ii)))*hh


sum = sum + free

enddo


! add surface term

ene = free + 2.D0*pi*(1.D0/(lamt**2))*rrr(nr)*hh(rn)


end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program siman



