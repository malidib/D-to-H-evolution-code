Program CrankNicolson
  ! Mohamad Ali-Dib & Jean-Marc Petit 19/01/2015
  ! Crank-Nicolson method to solve the D/H ratio equation. Fully working.

  integer, save :: nb,timecounter,iii,taille,outstep,adv
  integer, parameter :: outunit=44

  double precision, DIMENSION(:), ALLOCATABLE :: p,r,&
       sigmag,viscg,csg,alphag,tc,hg,rcm,&
       kt,kt2,pr,kappa,afrac,dkrs,krs,a1,b1,&
       sol,a2,aa,bb,cc,dd,c,d1,b2,d2,omega,dsigma,e,vrd,&
       g1,g2,nsr,dnsr,insource,afrac2,vrd2
  double precision, DIMENSION(:,:), ALLOCATABLE :: frac
  double precision :: deltat,timeyears,fluxout,lambda,&
       dummy,deltaxm,fluxin,deltaxp,rsnow,msol,mdot,mstel,tsnow,&
       midot,tff,alpham
  character(len=70) :: fn
  
  ! Arrays size
  taille = 190
  ! Allocate arrays sizes
  Allocate (p(taille),r(taille),omega(taille),&
       sigmag(taille),viscg(taille),csg(taille),alphag(taille),&
       tc(taille),hg(taille),rcm(taille),afrac2(taille),&
       kt(taille),pr(taille),kappa(taille),afrac(taille),&
       dkrs(taille),krs(taille),a1(taille),b1(taille),c(taille),d1(taille),&
       sol(taille),a2(taille),b2(taille),d2(taille),kt2(taille),vrd2(taille),&
       aa(taille),bb(taille),cc(taille),dd(taille),dsigma(taille),e(taille),&
       vrd(taille),g1(taille),g2(taille),nsr(taille),dnsr(taille),insource(taille))
  Allocate (frac(800000,taille))

  ! Opening Output file
  open(unit=11, file='outputCN.txt', action ='write')
  open(unit=111, file='inputfrac.txt', action ='read')
  open(unit=1111, file='outputmain.txt', action ='write')

  ! -------------------
  
  ! Snowline evolution related constants
  
  msol=1
  mstel=1
  tsnow=160
  midot=1e-6
  tff=1e5
  alpham=1e-2


  ! adv is advection, on=1, off=0
  adv=1

  ! Lambda = 0.5 for semi implicit, 0 for explicit and 1 for fully implicit scheme
  lambda=1

  ! Flux at the outer and inner boundaries 
  fluxout=0
  fluxin=0
  ! ---------

  iii=1
  timecounter=1
  nb=1

  ! Initial fractionation value, uncomment to read (a non uniform) initial frac from file
  do i=1,taille
!  read(111,*)dummy,dummy,dummy,frac(1,i)
  enddo
  
  frac(1,:)=35
  ! ----------

  ! Deltat = time step in seconds. Stability in explicit scheme a=kappa*deltat/deltax2 must be less than 0.5
  deltat=5e8
  ! -------

  ! Outstep = time between outputs
  outstep = 100
  ! -------

  ! Reading input and changing units to cgs
  open(unit=1, file='input.dat')
  read(1,*)
  do i=1,taille
     read(1,*)r(i),sigmag(i),omega(i),viscg(i),tc(i),dummy,csg(i),dummy,dummy,p(i)
!     omega(i)=omega(i)/(3.15e7)
     rcm(i)=r(i)*1.49E13
     tc(i)=tc(i)
     csg(i)=csg(i)*100
     Hg(i)=Hg(i)*1.49E13

     ! kt(i) is the chemical equilibrium constant and afrac(i) the fractionation at equilibrium
     afrac(i)=exp(0.04-13.17/sqrt(tc(i))+583.75/tc(i))

     if (tc(i).le.1000) then
     kt(i)=(0.005374*.1/1.01325e5)*exp(-5.17d3/tc(i))*273/tc(i)
     else
     kt(i)=(0.005374*.1/1.01325e5)*exp(-5.17d3/1000.)*273/1000.     
     endif
     
     ! Inward drifting ices source
     insource=0
     insource(34)=15*0/(1e3*3600*365)
     
!     kt2(i)=(2e-22)*exp(5.17e3/tc(i))
!     afrac2(i)=0.9505867*exp(-tc(i)/427.1239)+1.39824*exp(-tc(i)/131.968)+&
!     0.5666*exp(-tc(i)/50.3673)+1.0490

     ! alphag is the turbulence parameter
!     alphag(44:taille)=2e-4
     ! pr is Prandtl number
     pr(i)=0.5

!     old kappa function of alpha 
!     kappa(i)=alphag(i)*(csg(i)**2)/(pr(i)*omega(i))

     !smoothing over a few grid points to avoid divergence
!     viscg(49)=0.32e13/5
!     viscg(48)=viscg(49)/5
!     viscg(47)=viscg(48)/5
!     viscg(46)=viscg(47)/5
!     viscg(45)=viscg(46)/5

     
     kappa(i)=viscg(i)/pr(i)
     krs(i)=kappa(i)*sigmag(i)*rcm(i)
     alphag(i)=viscg(i)*omega(i)/(csg(i)**2)
     nsr(i)=viscg(i)*sigmag(i)*sqrt(rcm(i))
  enddo
  
   ! Boundary conditions
     bb(1)=-1
     cc(1)=1
     aa(taille)=-1
     bb(taille)=1
     dd(1)=fluxin*(rcm(2)-rcm(1))
     dd(taille)=fluxout*(rcm(taille)-rcm(taille-1))
     
  ! Main time loop
   do while (iii.le.799999)
   
     ! Snowline position calculations 
     mdot=midot*exp(-(timecounter*timeyears)/tff)
     if (mdot.lt.1e-8) then
!     mdot=1e-8
!     rsnow=1
      exit
      
!     rsnow=1.24*((mstel/msol)**0.333)*((mdot/(1e-8*msol))**0.444)*((tsnow/170)**(-1.1133))&
!            *((alpham/0.01)**(-0.222))     
     rsnow=5.73*((mstel/msol)**0.333)*((mdot/(1e-8*msol))**0.222)*((tsnow/170)**(-0.8))

     else
     ! Dead Zone solution
     rsnow=5.73*((mstel/msol)**0.333)*((mdot/(1e-8*msol))**0.222)*((tsnow/170)**(-0.8))
     
     ! Classical solution
     !rsnow=1.24*((mstel/msol)**0.333)*((mdot/(1e-8*msol))**0.444)*((tsnow/170)**(-1.1133))&
     !       *((alpham/0.01)**(-0.222))
     endif
     
     
     
     do jjj=1,taille
     if (rsnow.le.r(jjj)) then
!     taille=jjj
     exit
     endif
     enddo
     
 

     ! Main spatial loop
     do jj = 2, taille-1
        deltaxp =  rcm(jj+1)-rcm(jj)
        deltaxm =  rcm(jj)-rcm(jj-1)

        dkrs(jj)=(krs(jj+1)-krs(jj-1))/(deltaxp+deltaxm)
        dsigma(jj)=(sigmag(jj+1)-sigmag(jj-1))/(deltaxp+deltaxm)
        dnsr(jj)=(nsr(jj+1)-nsr(jj-1))/(deltaxp+deltaxm)

        a1(jj)=kappa(jj)*deltat*lambda/(0.25*(deltaxm+deltaxp)**2)
        b1(jj)=(dkrs(jj)*deltat/(sigmag(jj)*rcm(jj)))*(lambda/(deltaxm+deltaxp))
        d1(jj)=-kt(jj)*p(jj)*deltat*lambda

        a2(jj)=kappa(jj)*deltat*(1-lambda)/(0.25*(deltaxm+deltaxp)**2)
        b2(jj)=(dkrs(jj)*deltat/(sigmag(jj)*rcm(jj)))*((1-lambda)/(deltaxm+deltaxp))
        d2(jj)=-kt(jj)*p(jj)*deltat*(1-lambda)
        c(jj)=(kt(jj)*p(jj)*afrac(jj)+insource(jj))*deltat
      
        vrd(jj)=(-3/(sigmag(jj)*(sqrt(rcm(jj)))))*dnsr(jj)        
        vrd2(jj)=-3*(kappa(jj)*pr(jj))/(2*rcm(jj))
!        vrd(jj)=0
        e(jj)=2*kappa(jj)*dsigma(jj)/sigmag(jj)
        g1(jj)=adv*(vrd(jj)-e(jj))*lambda*deltat/(deltaxm+deltaxp)  
        g2(jj)=adv*(vrd(jj)-e(jj))*(1-lambda)*deltat/(deltaxm+deltaxp)  

        ! aa,bb and cc are respectively the matrix lower, main and upper diagonal arrays 
        cc(jj)=-b1(jj)-a1(jj)+g1(jj)
        bb(jj)=1+2*a1(jj)-d1(jj)
        aa(jj)=-a1(jj)+b1(jj)-g1(jj)
        
        dd(jj)=(a2(jj)+b2(jj)-g2(jj))*frac(iii,jj+1) &
             + (1-2*a2(jj)+d2(jj))*frac(iii,jj) &
             + (-b2(jj)+a2(jj)+g2(jj))*frac(iii,jj-1)+c(jj)
        ! dd is the right hand term in the equation M(aa,..)X=dd.

     enddo

             
     call solve_tridiag(aa,bb,cc,dd,sol,taille)
     timeyears=iii*deltat/(3600*24*365)

     if (timecounter*timeyears.ge.nb) then
        ! build filename -- i.dat
        write(fn,fmt='(i0,a)') nb, '.res'

        ! open it with a fixed unit number
        open(unit=outunit,file=fn, form='formatted')

        do jj=1,taille
        if ((jj.ne.43).and.(jj.ne.44).and.(jj.ne.45)) then
           write(11,*)timecounter*timeyears,r(jj),tc(jj),sol(jj),vrd(jj),vrd2(jj)

           write(outunit, *) timecounter*timeyears,r(jj),tc(jj),sol(jj)
        endif  
!           if ((timecounter*timeyears).lt.400000) then   
           if (mdot.gt.2e-8) then        
           if (jj.eq.taille) then
           write(1111,*)timecounter*timeyears,r(jj),sol(jj)
           endif 
           elseif (mdot.le.2e-8) then
!           elseif ((timecounter*timeyears).ge.400000.and.(timecounter*timeyears).lt.500000) then
           write(1111,*)timecounter*timeyears,r(jj),sol(jj)
           endif
            
        enddo
        ! close it 
        close(outunit)

        if (nb .gt. 10*outstep) outstep = 10*outstep
        nb=nb+outstep

        print*,timecounter*timeyears,timecounter,timeyears,iii
     endif

     if (iii.le.799998) then

        do jj=1,taille
           frac(iii+1,jj)=sol(jj)
        enddo
     else
        do jj=1,taille
           frac(1,jj)=sol(jj)
        enddo
        iii=0
        timecounter=timecounter+1
     endif

     iii=iii+1
     if (nb .gt. 1000000) exit

  enddo
!  call system ("bash movedocs.sh")

!  call system (gnuplot gnuplotinput.plt")
  
!  call system (evince test.ps")

End Program


Subroutine solve_tridiag(a,b,c,d,x,n)
! Thomas algorithm (simplified form of Gaussian elimination
! that can be used to solve tridiagonal systems of equations)

           implicit none
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations

        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: a,b,c,d
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer i

! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do
!print*, x
return
end subroutine solve_tridiag
