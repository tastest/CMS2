c BAYES - program to check equivalency of Bayes and frequentist 
c         upper limit calculations
c
c                              John Conway   Nov. 2004

      implicit none

      integer nev  ! number of observed events
      real*8  sac  ! relative uncertainty in acceptance
      real*8  xbg  ! expected background
      real*8  sbg  ! uncertainty in background
      real*8  cl   ! desired CL
      real*8  smax ! upper limit of integration
      real*8  prec ! integration step size
      real*8  plim ! upper limit on Poisson process
      real*8  perr ! error on upper limit on Poisson process

      logical done
      character*1 yes

      print '('' ------------------------------------------------- '')'
      print '('' Bayes version 2 (Nov. 2004)                       '')'
      print '('' ------------------------------------------------- '')'
      print '('' This program, bayes.f, performs a calculation of  '')'
      print '('' the upper limit on a Poisson process with         '')'
      print '('' background, incorporating the effects of          '')' 
      print '('' uncertainty on the background and  signal         '')'
      print '('' acceptance.  The program assumes a uniform prior  '')'
      print '('' pdf in the signal rate, and Gaussian priors for   '')'
      print '('' the acceptance and expected background.  The      '')'
      print '('' Gaussian priors are truncated at zero.  Formally, '')'
      print '('' therefore, the integral of the posterior density  '')'
      print '('' with respect to the signal rate diverges          '')'
      print '('' logarithmically.  This is a small effect          '')'
      print '('' (requiring an integral out to very large signal   '')'
      print '('' rates to be observed) if the acceptance           '')'
      print '('' uncertainty is small, less than about 15-20%.  If '')'
      print '('' your acceptance uncertainty is larger than this   '')'
      print '('' your resulting limit may be sensitive to the upper'')'
      print '(''  limit of integration.  In this case, you have a  '')'
      print '('' problem, and you should consider an alternative   '')'
      print '('' approach.                                         '')'
      print '(''                                                   '')'
      print '('' If you have any questions about these issues      '')'
      print '('' please contact John Conway at conway@fnal.gov.    '')'
      print '('' ------------------------------------------------- '')'
      print '('' Be careful to enter a real number (with a decimal '')'
      print '('' point) for the values below requesting it.        '')'
      print '('' ------------------------------------------------- '')'
    
      print '($,'' enter number of observed events (integer): '')'
      read(*,'(i4)') nev

      print '($,'' enter relative error on acceptance (real): '')'
      read(*,'(f10.4)') sac

      if(sac.gt.0.15) then
        print '('' Your error is quite large - please contact  '')'
        print '('' John Conway at conway@fnal.gov !            '')'
      endif

      print '($,'' enter expected number of background events '',
     .        ''(real): '')'
      read(*,'(f10.4)') xbg

      print '($,'' enter error on number of background events '',
     .        ''(real): '')'
      read(*,'(f10.4)') sbg

      print '($,'' enter desired confidence level (real): '')'
      read(*,'(f10.4)') cl

      smax = 10*float(nev)
      if(smax.lt.10.0) smax = 10.0

      done = .false.
      do while(.not.done)
        print '($,'' enter integration upper limit (real): '')'
        read(*,'(f10.4)') smax
        if(smax.lt.10.0) then
          print '(/,'' Value of upper limit too small! '')'
        elseif(smax.lt.float(nev)+5.0*sqrt(float(nev))) then
          print '(/,'' Value of upper limit too small! '')'
        else
          done = .true.
        endif
      enddo

      print '($,'' enter integration step size '',
     .          ''(recommend 0.02 or less):'')'
      read(*,'(f10.4)') prec

      print '(/)'
      print '('' ----------------------------------------- '')'
      print '('' observed events                '',i10)',nev
      print '('' relative error on acceptance   '',f10.3)',sac
      print '('' expected background            '',f10.3)',xbg
      print '('' absolute error on background   '',f10.3)',sbg
      print '('' desired confidence level       '',f10.2)',cl
      print '('' integration upper limit        '',f10.2)',smax
      print '('' integration step size          '',f10.4)',prec
      print '('' ----------------------------------------- '')'
      print '(/)'
      print '($,'' Are the above correct? '')'
      read(*,'(a)') yes
      if(yes(1:1).ne.'y') stop

      call baylim(nev,sac,xbg,sbg,cl,prec,smax,plim,perr)

       print '(//,'' limit: less than'',
     .         f10.3,'' signal events '',//)',plim

      end


      subroutine baylim(nev,sac,xbg,sbg,cl,prec,smax,plim,perr)

c calculate Bayesian upper limit on Poisson process

      implicit none

c arguments

      integer nev       ! number of events
      real*8    sac     ! error on acceptance
      real*8    xbg     ! expected background
      real*8    sbg     ! error on background
      real*8    cl      ! desired CL
      real*8    smax    ! upper limit of integration
      real*8    prec    ! integration step size
      real*8    plim    ! Poisson upper limit  (returned)
      real*8    perr    ! error on upper limit (returned)

      real*8    xev,xevmax,dxev

      logical done

      real*8 blike,xlike

      integer i,nlist,icl
      real*8 blist(10000),xlist(10000),bint,bsum,bcl

c scan likelihood L(nev|xev,...)

      xevmax = smax
      dxev = prec

      nlist = 0
      bsum = 0.

      done = .false.

      xev = dxev/2.

      do while(.not.done)

        xlike = blike(nev,sac,xbg,sbg,xev)

        nlist = nlist + 1
        xlist(nlist) = xev
        blist(nlist) = xlike

        bsum = bsum + xlike

        print '($,i5,2g12.5,a)',nev,xev,xlike,char(13)

        xev = xev + dxev
        if(blist(nlist)/blist(1).lt.1.0e-6) done = .true.
        if(xev.gt.xevmax) done = .true.

      enddo

c find place just below threshold
      bint = 0.
      icl = 0
      do i=1,nlist
        if(bint         .le.cl*bsum.and.
     .     bint+blist(i).gt.cl*bsum) then
          icl = i
          bcl = bint
        endif
        bint = bint + blist(i)
      enddo

c interpolate linearly
      plim = xlist(icl) + (xlist(icl+1)-xlist(icl))*
     .                     (cl*bsum-bcl)/blist(icl+1)

      perr = 0.

      return
      end


      real*8 function blike(nev,sac,xbg,sbg,xev)

c return likelihood to observe nev events given expected

      implicit none

      integer nev     ! number of events observed
      real*8    xev   ! number of events expected
      real*8    sac   ! error on acceptance
      real*8    xbg   ! expected background
      real*8    sbg   ! error on background

      integer i,nmax
      parameter(nmax=100000)

      real*4  a,b
      real*8  yyev,ssev,yyex,xxev,ssac,xxbg,ssbg,yybg,xint,xxx
      real*8  stirling

      logical done

c perform double gaussian integral by Monte Carlo

      xxev = nev
      xxbg = xbg
      ssbg = sbg
      ssac = sac
      xxev = xev

      ssev = xxev*ssac

      xint = 0.

      done = .false.

      do i=1,nmax

c pick expected background and signal from Gaussian

        yybg = -1.
        yyev = -1.

        do while(yybg.lt.0..or.
     .           yyev.lt.0.)
          call rannor(a,b)
          yybg = xxbg + a*ssbg
          yyev = xxev*(1.0+b*ssac)
        enddo

c total expected         
        yyex = yybg + yyev

c value of integrand

        xxx = yyex**nev * dexp(-yyex) / stirling(nev)  ! Poisson term

        xint = xint + xxx

      enddo

      blike = xint/float(nmax)

      return
      end


      real*8 function stirling(n)

c Stirling's approximation for n! with first correction term

      real*8 x,fact(0:9)
      parameter (root2pi=2.506628)
      data fact/1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880./

      if(n.ge.0.and.n.le.9) then
        stirling=fact(n)
        return
      endif
      if(n.gt.0) then
        x=n
        stirling=root2pi*sqrt(x)*x**n*exp(-x)*(1.+1./12./x)
      else
        stirling=1.
      endif
      return
      end

