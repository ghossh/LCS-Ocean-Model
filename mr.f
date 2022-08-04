!  This the modified lcs_code, 
!  program numerical solutions part 2. indsgeos.for 
!  this program solves equations of motion of a continuously stratified
!  flow by time-stepping in a staggered grid system. the damper is added
!  in this code. 
!this version uses u_choice,v_choice,p_choice
!  the initial conditions for the time-stepping are u=v=w=p=rho=0.
! this code you will use to read realistic wind data. the program will read zonal wind field, meriodional wind filed data from input nc file,it will calculate wind stress from that values.all the units are in cgs.

   	include 'netcdf.inc' 
        
         
       include "./header_files/dimension.h"
        

      parameter (imxp=imax+1)
      logical storit,prntit
      parameter(nmodes=25)
      dimension sigma(nmodes),sigp(imax,jmax)

      common/amode/xpsi2(nmodes),cn(nmodes),zn(nmodes)

      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     1 rhon(imax,jmax),dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     5  fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
       

      common/realbd/istarr(jmax),indarr(jmax)    ! included

      dimension tauxx(imax,jmax) 
      dimension tauyy(imax,jmax)
      common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
      common/pcom/rmsk(imax,jmax)
      common/zcom/zetun(imax,jmax),zetvn(imax,jmax)
      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/msm, ijk, iii, idate3, iday2,i11
      common/cntrl8/dt
      real*8 rtime,dt

      logical wclos,sclos,eclos,nclos

c  start is true if restarting from previous run, otherwise false
      logical start

c  nstart is the day number of restart (it is only meaningful with start=true)
c  ndbeg and ndend are start and end times in days.
c  nbeg and nend are the corresponding time step numbers.
c  ndinc is the number of days to be run.
c  ndprin is the print interval in days.
c  ndstor is the storing interval in days.
c  ismth is the number of time steps between time smoothing  
c  logicals wclos,sclos,eclos,nclos set which boundaries are closed
      character*80 filplt, filresi,filreso
      character*8 tim,dat*9
      character*15 outfil
      character*60  maskfile  
      character*60  initpre
      character*80 maskin
      character*80 windfile
      character*2 text,text2
      character*6 dir, ext
      character*8 dir2
        character*10 dir3
        character*18 ext2
        character*17 ext3
        character*30 res_fil1
        character*29 res_fil2

      
      integer status, outid, ncid, length, kap, iap, jap,i,j
      integer iopen,idate
      integer londim, latdim, timedim
      integer lonid, latid, timeid, ntime,nt
      character*20 torigin
      real time1,time
      real delta,slat,elat,slon,elon  
      integer splons,splone,splats,splate,sbc_apply
      integer dplons,dplone,dplats,dplate,damp_equ,damp_u_south
      dimension istar2(jmax),indar2(jmax) ! needed for special int pts near 
      real step1,n_val,time_step                                                         ! india
      integer  midlat,midl(1000)
      real yo,fou,fo,fyo,save_it
      integer ddj
      integer start_mode,end_mode,start_program

	 include "./header_files/param.h"
        include "./header_files/input_files.h"
! For restart, please use start=.true., otherwise false----------------------
        if (start_program.eq.1) then 
	 start=.true.
        else
        start=.false.
        endif 
       
!-------------------------------------------------------------------------------
	 wclos=.false.
        sclos=.false.
        eclos=.false.
        nclos=.false.

	 if (sbc_apply.eq.1)then
        call lonlat(splons,splone,splats,splate)
	 endif
        step1=dt/wind_dif_time
        n_val=(86400.00*step1)/dt
        time_step =86400.00/(dt*n_val)

c        parameters for the damper damp(imax,jmax)
        if (damp_equ.eq.1)then      
        call  lonlat_damp(dplons,dplone,dplats,dplate)
         ii1=dplons
         ii2=dplone
         jj1=dplats
         jj4=dplate
         ddi=ii2-ii1
         ddj=int(2./delta)
         jj2=jj1+ddj
         jj3=jj4-ddj
c         write(*,*)jj3,'666'
         endif

c      data    zn/0.98,0.96,0.91,0.83,0.74,0.65,0.54,0.45,0.37,0.28/
c  ------------------------------------------------------------------------

      length=imax*jmax*12*4

      idate1=15 ! to write on (idate +1)th of every month

      data ndmon/31,28,31,30,31,30,31,31,30,31,30,31/
      call origin(torigin,iday,month,iyear)
!      data torigin / '01-jan-2000 00:00:00' / ! 20 char including spaces

c parameters for the model

      imaxm=imax-1
      jmaxm=jmax-1
      ijmax=imax*jmax
      ijmax3=ijmax*3  
      pi=4.0*abs(atan(1.))
      g=980.
      rearth=6370.e5
      dypdeg=pi/180.*rearth
      fday=2.*pi/86400.
c      beta=2.*fday/rearth
 
c  -----------------------
c  run parameters
c  -----------------------
                  
      istart=1
      iend=imax
      iskip=5
      jstart=1
      jend=jmax
      jskip=-5   
      ndbeg=0
      write(text2,109) start_mode
109   format(i2.2)
      if (start) then
      dir2='./input/'
        ext2='restart.txt'
        res_fil1=dir2//text2//ext2
       open(unit=24,file=res_fil1,status='old')
       read(24,*) ndbeg,kp,k,km,nsmth1
        ks=km
        km=k
        k=kp
        kp=ks
      endif
        if (start_program.eq.1) then 
        ndinc=ndrun+ndbeg
        else
        ndinc=ndrun
        endif  

!       open(unit=24,file='restart_read.txt',status='old')
!       read(24,*) ndbeg,kp,k,km,nsmth1
!        ks=km
!        km=k
!        k=kp
!        kp=ks
!      endif
      istore=0
c  set starting month
!      month=1
      iday_kount = month   !  counter for output
      iday1=month
!      iyear=2000
      oyear=iyear
      iday=1
      iday2=iday
      iday3=iday
       nn=month
!       ndyr=365
       if (iyear .ne. 0) then
        if (mod(iyear,4) .eq. 0) then
         ndmon(2)=29
         ndyr=366
        else
         ndmon(2)=28
         ndyr=365 
        endif
       endif
        
c  -----------------------
c  ocean parameters  


c  set basin and resolution parameters
       
      ymin=(slat+0.5*delta)*dypdeg
      ymax=(elat-0.5*delta)*dypdeg
      xmin=(slon+0.5*delta)*dypdeg
      xmax=(elon-0.5*delta)*dypdeg
      xlen=xmax-xmin
      dx=(xmax-xmin)/(imax-2)
      dy=(ymax-ymin)/(jmax-2)
      recdx=1./dx
      recdx2=recdx*recdx
      recdy=1./dy
      recdy2=recdy*recdy

      print 400,dx,dy,dypdeg,dt,g,
     2 ymin,ymax,xlen,xmin,xmax
400   format(/,1x,'model run parameters:',/,
     1       1x,'dx    = ',e10.4,/,1x,'dy    = ',e10.4,/
     1       1x,'dypdeg  = ',e10.4,/
     1       1x,'dt    = ',e10.4,/,1x,'g= ',e10.4,/
     1       1x,'ymin  = ',e10.4,/,1x,'ymax  = ',e10.4,/
     1       1x,'xlen  = ',e10.4,/,1x,'xmin  = ',e10.4,/
     1       1x,'xmax  = ',e10.4,/ )
c  ---------------------------------------------------------------------

c      calculation of midlatitude
       yo=midlat*pi/180
c       fyo=(midlat+delta)*dypdeg
       fyo=midlat*dypdeg
       beta=2.*fday*cos(yo)/rearth
       fo=2.*fday*sin(yo)
c       write(*,*)dypdeg,'midlat'
        write(71,*)beta,'888888'
c define coriolis force at f-plane
       if (f_plane.eq.1)then
       do 12 j=1,jmax
       do 12 i=1,imax
       fv(i,j)=fo
       fu(i,j)=fo
12    continue
	 else
c  define coriolis force at midlat plane
      do 11 j=1,jmax
      y=ymin+(j-1.)*dy
      do 11 i=1,imax
      fv(i,j)=fo+beta*(y-fyo)
      fu(i,j)=fo+beta*(y-fyo-dy/2.)
11    continue
      endif
     
c  define dampu array for the southern boundary damper
      x1=xlen-300.e5
      x2=xlen-150.e5
c  define damper along open southern boundary on u field only
      y1=ymin+300.e5
      y2=ymin+150.e5
      do 513 j=1,jmax
      yu=ymin+(j-1.5)*dy
      do 513 i=1,imax
      if(yu.gt.y1) then
        dampu(i,j)=0.
      else if(yu.lt.y2) then
        dampu(i,j)=1.
      else
        dampu(i,j)=(yu-y1)/(y2-y1)
      endif
513   continue


        call landmask(maskfile)

	 open(11,file='./in/cn.dat',status='old')
        open(12,file='./in/zn.dat',status='old')
        open(13,file='./in/xpsi2.dat',status='old')
         
!        open(22,file='restart_write.txt')

       do n=1,nmodes
       read (11,*) cn(n)
c       write(*,*) cn(n),'ccccccccc'
       enddo
       close(11)
        do n=1,nmodes
       read (13,*) xpsi2(n)
c       write(*,*) xpsi2(n),'xxxxxxxxxx'
       enddo
       close(13)
c       write(45,*) n,xpsi2(n)
       do 24 n=1,nmodes
        read(12,*)zn(n)
c        write(*,*) zn(n),'zzzzzzzzzz'
24     continue

       do 34 n=1,nmodes
       sigma(n)=cn(n)/(1.5*dx)
34     continue

c	nm = num_mode

       do 3500 nm=start_mode,end_mode
        
        if (damp_equ.eq.1)then
       do 44 i=1,imax
       do 44 j=1,jmax
       if(i.lt.ii1) damp(i,j)=0.0
       if(i.ge.ii1.and.j.gt.jj4) damp(i,j)=0.0
       if(i.ge.ii1.and.j.lt.jj1) damp(i,j)=0.0

        if(i.gt.ii2.and.j.ge.jj2.and.j.le.jj3)damp(i,j)=sigma(nm)
        if(i.gt.ii2.and.j.ge.jj1.and.j.lt.jj2)then
        damp(i,j)=sigma(nm)*((j-jj1)/ddj)

        elseif(i.gt.ii2.and.j.gt.jj3.and.j.le.jj4)then
        damp(i,j)=sigma(nm)*((jj4-j)/ddj)
        endif

        if(i.ge.ii1.and.i.le.ii2)go to 54
        go to 44
54      continue
        if(j.ge.jj1.and.j.le.jj2)then
        damp(i,j)=sigma(nm)*((i-ii1)/ddi)*((j-jj1)/ddj)
         elseif(j.gt.jj2.and.j.le.jj3)then
          damp(i,j)=sigma(nm)*((i-ii1)/ddi)
          elseif(j.ge.jj3.and.j.le.jj4)then
        damp(i,j)= sigma(nm)*((i-ii1)/ddi)*((jj4-j)/ddj)
         endif
c         write(21,*)i,j,damp(i,j),dampu(i,j)
44        continue
          else
         do j=1,jmax
         do i=1,imax
          damp(i,j)=0.0
c        write(22,*)i,j,damp(i,j),dampu(i,j)
         enddo
        enddo
        endif

        if (damp_u_south.eq.0) then
          do j=1,jmax
          do i=1,imax
          dampu(i,j)=0.0
          enddo
          enddo
          endif  
        nmm=nm
      write(text,106) nmm
106   format(i2.2)

      dir3='./restart/'
        ext3='restart.txt'
        res_fil2=dir3//text//ext3
        open(unit=22,file=res_fil2,status='unknown')

      dir='./out/'
      ext='out.nc'
      outfil=dir//text//ext
      nt=0
      idate=idate1
      iyear=oyear
      nn=month
      iday_kount=month
      iday=iday3
      ndyr=365
      if (iyear .ne. 0) then
       if (mod(iyear,4) .eq. 0) then
        ndmon(2)=29
        ndyr=366
       else
        ndmon(2)=28
        ndyr=365
       endif
       endif

c       initialize arrays
        tdt=dt
c        kp=3
c        k= 2
c        km=1
c        nsmth1=0

	if (start) then
	 call readrar(1)
        tdt=2*dt
	else
         kp=3
         k= 2
         km=1
         nsmth1=0
         call setrar(un,ijmax3,0.0)
         call setrar(vn,ijmax3,0.0)
         call setrar(pn,ijmax3,0.0)
	endif


c       c  define variables depending on run parameters
!      ndend=ndbeg+ndinc
      
      ndend=ndinc
      ndprin=100000
      ndstep=int(86400./dt +.5)
      nbeg=ndbeg
      nend=ndend*ndstep
      nprint=ndprin*ndstep
      nstore=ndstor*ndstep

	write(*,*) ndend, nend,ndstep,dt

c  ------------------------------------------------------------------------
c  print out print and store model parameters
      print 432, ndbeg,ndend,ndinc,ndprin,ndstor,ndstep,nprint,nstore
432   format(/,1x,'time/print/store parameters:'/,
     1           1x,'ndbeg = ',i5/,
     1           1x,'ndend = ',i5/,
     1           1x,'ndinc = ',i5/,
     1           1x,'ndprin= ',i5/,
     1           1x,'ndstor= ',i5/,
     1           1x,'ndstep= ',i5/,
     1           1x,'nprint= ',i10/,
     1           1x,'nstore= ',i5/)

c  begin time-stepping 
      npstp=0
      nsstp=0
      nsmth=nsmth1
      nstep=0
      kday=ndbeg
      ndsave=kday
      rtime=kday

      iopen=1
      call output(iopen,outfil,nt,torigin)


3000  continue        !  <================================== time-stepping loop
c  all time stepping loops are called from subroutine ocean

        call ocean(rtime,ndinc,splons,splone,splats,splate,sbc_apply,
     +    n_val,time_step)  
c  increment time parameters
      nstep=nstep+1
      npstp=npstp+1
      nsstp=nsstp+1
      nsmth=nsmth+1
      rtime=rtime+dt/86400.
      if(mod(nstep,ndstep).eq.0) then
        nstep=0
        kday=kday+1
        ndsave=kday
        rtime=kday
        iday=iday+1
        i11=0
      endif
      i11=i11+1
      if (i11 .ge. time_step ) then
        i11=0
        elseif (mod(rtime,1.0).eq.0.0) then
        i11=0
        endif

       if(iday.gt.ndmon(nn))then
2221    continue
        iday=iday-ndmon(nn)
        nn=nn+1
        if(nn.eq.13)then
         nn=1
         iyear=iyear+1
         if(mod(iyear,4).eq.0)then
          ndmon(2)=29
          ndyr=366
         else
          ndmon(2)=28
          ndyr=365
         endif
        endif
        if(iday.gt.ndmon(nn)) go to 2221
       endif

c  smooth fields, sm2dt for 2dt noise
      if(mod(nsmth,ismth).eq.0.0) then
        call sm2dt
        nsmth=0
      endif
      nsmth1=nsmth  !store for restart file


c  check for printing and/or storing fields
      prntit=.false.
      storit=.false.
c      storit=.true.

      if(mod(npstp,nprint).eq.0.0) then                       !print fields
        npstp=0 
        prntit=.true.
      endif  

       if (ndsave .eq. 0 .and. mod(nstep,ndstep) .eq. 1.0)then
        if (iday2 .lt. idate1) then
         idate=idate1-iday2+1
        else
         idate=idate1+ndmon(month)-iday2+1
         iday_kount=iday_kount+1
        endif
       endif
	if(save_it .ne. 0.0)then
         if((mod(rtime,save_it).eq.0.0).and.(rtime.gt.0.0)) then
         storit =.true.
	 endif
	else
	storit =.false.
        endif

        idate=idate+ndmon(iday_kount)
        iday_kount=iday_kount+1
        
        if (iday_kount .gt. 12 ) iday_kount=1
c        endif

      day=ndsave
        if(prntit) then     ! printing
        if(nmm.ne.1) go to 5186
        print 5196,day
5196    format('1','the taux-field after ',f10.1,' days ')
        print 31960,day
31960   format('0','the tauy-field after ',f10.1,' days ')
5186    continue
        print 439, day,nmm
439     format(1x,'the un field after ',f6.1,' days for mode',i4)

       endif

       if(storit) then      ! store the fields
       print 5660, day,nmm
5660    format(1x,'storing  fields for day ',f6.0,5x,'mode ',i4)
       iopen=2
       nt=nt+1
        call output(iopen, outfil,nt,torigin)
        endif

8000  continue

c	if(ndsave.eq.3600) goto 9000

      if(.not.(kday.eq.ndend.or.nstep.eq.nend))then
        tdt=2.*dt
        ks=km
        km=k
        k=kp
        kp=ks
        goto 3000    ! <================================== do another time step
      else
        goto 9000    ! <================================== done
      endif
c  -------------------------------------------------------------------------

9000  continue       ! program completed 

c  store last day output for restart
      print 5665, day,nmm
5665  format(1x,'saving restart fields for day ',f6.0,5x,'mode',i4)
     
c       write(23) ndsave,nsmth1,kp,k,km,un,vn,pn
	call restart(1,torigin)
        write(22,*) ndsave,kp,k,km,nsmth1
     
3500  continue

      stop
      end

! Refining rmask in included

        subroutine landmask(maskfile)

	include 'netcdf.inc'
	include "./header_files/dimension.h"

c  !      parameter (imax=ikmax-1)
c  !      parameter (jmax=jkmax-1)

        common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1  month,istore,nsmth1,ndsave,ndend,ndstep,
     2  ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3  wclos,sclos,eclos,nclos

        common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1  wn(imax,jmax,3),pn(imax,jmax,3),
     1  rhon(imax,jmax),dampu(imax,jmax),damp(imax,jmax),
     4  taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     5  fv(imax,jmax),tau1y(imax,jmax,2),
     5  rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)

        common/realbd/istarr(jmax),indarr(jmax)    ! included
        common/cntrl8/dt

        dimension tauxx(imax,jmax)
        dimension tauyy(imax,jmax)
        common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
        common/pcom/rmsk(imax,jmax)
        common/zcom/zetun(imax,jmax),zetvn(imax,jmax)
        character*60  maskfile
	real*8 dt
        
        call readmask(maskfile)
!DEFINE ISTARR(J), WEST-MOST INTERIOR POINT INDEX

        do j=1,jmax
         istarr(j)=1
         indarr(j)=imax
        enddo

       do J=2,JMAXM
        do I=1,IMAXM
         if(rmask(i,j).eq.0. .and. rmask(i+1,j).eq.1) then
          istarr(J)=I+1
          goto 111
         endif
        enddo
111     continue
       enddo

        istarr(1)=istarr(2)
        istarr(JMAX)=istarr(JMAXM)

!  DEFINE INDARR(J), EAST-MOST INTERIOR POINT INDEX
        do J=2,JMAXM
        do I=IMAX,2,-1
                if(rmask(i,j).eq.0. .and. rmask(i-1,j).eq.1) then
                indarr(J)=I-1
                goto 222
                endif
        enddo
 222    continue
        enddo

        indarr(1)=indarr(2)
        indarr(JMAX)=indarr(JMAXM)
 
        return
        end
 
       subroutine setrar(a,len,value)


         include "./header_files/dimension.h"

      dimension a(len)
      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      do 1 ll=1,len
    1 a(ll)=value
      return
      end

       subroutine ocean(rtime,ndinc, splons,splone,splats,splate,
     +  sbc_apply,n_val,time_step)

        include "./header_files/dimension.h"

!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
      parameter(nmodes=25)

      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/amode/xpsi2(nmodes),cn(nmodes),zn(nmodes)

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     2  rhon(imax,jmax),
     3 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     4 fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
     

      common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
      common/zcom/zetun(imax,jmax),zetvn(imax,jmax)
      logical wclos,sclos,eclos,nclos

      common/realbd/istarr(jmax),indarr(jmax)
      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      common/cntrl8/dt
c       common/specb/splons,splone,splats,splate
	integer splons,splone,splats,splate,sbc_apply
	real n_val
	real*8 rtime,dt

       include "./header_files/param.h"
      
!    calculate interpolated forcings for this time step
      call getobs(rtime,ndinc,n_val,time_step)

!  calculate un,vn,wn,pn,rhon fields for each mode
      call modecoef
       
!  apply boundary conditions on rectangular edges for un,vn,pn fields
      call bndry2(1)
       call ungeneral
      call vngeneral
	if (sbc_apply==1)then
       call unspecial(splons,splone,splats,splate)
      call vnspecial (splons,splone,splats,splate)
	endif
      call pnbound
 
      return
      end

      subroutine getobs(rtime,ndinc,n_val,time_step)

c  this routine interpolates monthly mean fields to obtain the required
c  forcing functions at each time step.

        include "./header_files/dimension.h"


!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
      parameter(nmodes=25)
      logical wclos,sclos,eclos,nclos

      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     2 rhon(imax,jmax),
     3 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     5  fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
     

      common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
      common/realbd/istarr(jmax),indarr(jmax)
      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/msm, ijk, iii, idate3, iday2,i11
      common/cntrl8/dt
      integer msm, m1, m2, imon,m3,flag
      real n_val,time_step
      real*8 dt,rtime
      
	character*80 windfile 
        character*23 maskfile
        character*60  initpre
        include "./header_files/input_files.h"
        include "./header_files/param.h"
      

c  determine months, m1 and m2, to use for interpolation
       
       if (ndsave.eq.0) then
        msm=0
        iii=month
        ijk=month
        idate3=ndmon(month)
       endif

2001   continue

        if (ndsave .eq. idate3 ) then
        msm=msm+ndmon(iii)
        iii=iii+1
        ijk=ijk+1
         if (iii.eq.13) iii=1
        idate3=idate3+ndmon(iii)
        endif
	m3=0
        m1=int(n_val*(rtime)+1)
	if(rtime.gt.0.0)then
	m3=int(n_val*(rtime-dt/86400)+1)

	endif
        m2=m1+1
    
!	write(44,*) 'call',ndsave,m1,m2,m3,rtime
       	
       call readwind(windfile,m1,m2,m3,rtime)
	
2002    continue

          fac=time_step
c             fac=15

      do 100 j=1,jmax
      do 100 i=1,imax

c  wind stress

       txsav(i,j)=(taux(i,j,1)*(fac-i11) +i11*taux(i,j,2))/fac
       tysav(i,j)=(tau1y(i,j,1)*(fac-i11) +i11*tau1y(i,j,2))/fac

!	if ((i .eq. 446) .and. (j .eq. 950)) then
!         write(42,*) ndsave, taux(i,j,1),taux(i,j,2),
!     #  tau1y(i,j,1),tau1y(i,j,2)
!	endif
      

100   continue
 
c        write(*,34)m1,m2,rtime,tau1y(1,1,1),tau1y(1,1,2),tysav(1,1),i11
c       write(*,34)m1,m2,rtime,taux(1,1,1),taux(1,1,2),txsav(1,1),i11
34      format(i4,1x,i4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,i4)

      return
      end

    
      subroutine modecoef
c  this subroutine time steps the equatoins for un,vn, pn and therefore
c  get wn and rhon 

        include "./header_files/dimension.h"


!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
      parameter(nmodes=25)

      common/amode/xpsi2(nmodes),cn(nmodes),zn(nmodes)

      logical wclos,sclos,eclos,nclos
      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     2 rhon(imax,jmax),
     3 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     5  fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
     

      common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
      common/zcom/zetun(imax,jmax),zetvn(imax,jmax)
      common/cntrl8/dt
       integer uchi(4,16000), uchj(4,16000)
       integer uch(4)
	real*8 dt

      common/pcom/rmsk(imax,jmax)
      common/realbd/istarr(jmax),indarr(jmax)
      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      timkwd=1.*86400.
     
!  calculate zeta's

       do 100 j=2,jmaxm
      jm=j-1
      jp=j+1
      is=istarr(j)
      ie=indarr(j)
      do 100 i=is,ie
      im=i-1
      ip=i+1
      zetun(i,j)=
     1  (un(im,j,km)+un(ip,j,km)-2.*un(i,j,km))*recdx2
     1  + ((un(i,jp,km)-un(i,j,km))-
     1  (un(i,j,km)-un(i,jm,km)))*recdy2


      zetvn(i,j)=
     2  ((vn(ip,j,km)-vn(i,j,km))-
     2  (vn(i,j,km)-vn(im,j,km)))*recdx2
     2  + (vn(i,jp,km)+vn(i,jm,km)-2.*vn(i,j,km))*recdy2	

100   continue

c  calculate interior fields
      do 9000 j=2,jmaxm
      jm=j-1 
      jp=j+1
      is=istarr(j)
      ie=indarr(j)


! Calculate un field

      do 9001 i=is,ie
      im=i-1
      ip=i+1


        if(rmask(i,j).eq.0.0) goto 9001

c  time step fields
	un(i,j,kp)=un(i,j,km)+ tdt*rmask(i,j)*(
     1   txsav(i,j)*zn(nmm)/xpsi2(nmm)+ visch*zetun(i,j)
     1   -recdx*(pn(ip,j,k)
     2   -pn(i,j,k))-a*un(i,j,km)/cn(nmm)**2+
     3   0.25* fu(i,j)*(vn(i,j,k)+vn(ip,j,k)+vn(i,jm,k)+
     1   vn(ip,jm,k))-(1./timkwd)*dampu(i,j)*un(i,j,km)
     1 -damp(i,j)*un(i,j,km))
9001  continue

! Calculate vn field

      do 9002 i=is,ie
      im=i-1
      ip=i+1

        if(rmask(i,j).eq.0.0) goto 9002

          vn(i,j,kp)=vn(i,j,km)+ tdt*rmask(i,j)*(
     1 tysav(i,j)*zn(nmm)/xpsi2(nmm)+visch*zetvn(i,j)
     2 -0.25*fv(i,j)*(un(im,jp,k)+un(i,jp,k)+un(im,j,k)+
     1 un(i,j,k))
     3 - recdy*(pn(i,jp,k)-pn(i,j,k))-a*vn(i,j,km)/cn(nmm)**2)

!       if ((i .eq. 446) .and. (j .eq. 950)) then
!	write(32,112) ndsave,kp,k,km,vn(i,j,kp),vn(i,j,k),vn(i,j,km),
!     #   tdt,tysav(i,j),fv(i,j),un(i,j,kp),un(i,j,k),pn(i,j,km),
!     #  visch
!       endif
9002  continue
112   format(4i4,1x,10f12.4)


! Calculate pn field

      do 9003 i=is,ie
      im=i-1
      ip=i+1
       if(rmask(i,j).eq.0.0) goto 9003
      pn(i,j,kp)=pn(i,j,km)+tdt*rmask(i,j)*(cn(nmm)**2*(
     1 -recdx*(un(i,j,k)-un(im,j,k))
     2 -recdy*(vn(i,j,k)-vn(i,jm,k)))-a*pn(i,j,k)/cn(nmm)**2
     1 -damp(i,j)*pn(i,j,km))
9003  continue
      do 9004 i=is,ie
      im=i-1
      ip=i+1
      rhon(i,j)=-pn(i,j,kp)/g
      wn(i,j,kp)=
     1 -recdx*(un(i,j,kp)-un(im,j,kp))
     2 -recdy*(vn(i,j,kp)-vn(i,jm,kp))

9004  continue

9000  continue

      return
      end

      subroutine bndry2(ibc)
c  this routine calls the subroutines that set the boundary conditions
c  for the field variables along the outer (rectangular) boundaries of
c  the basin

        include "./header_files/dimension.h"


!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
      parameter(nmodes=25)
      logical wclos,sclos,eclos,nclos


      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     1 rhon(imax,jmax),
     1 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     5  fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
   

      common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
      common/realbd/istarr(jmax),indarr(jmax)

      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      common/cntrl8/dt

	real*8 dt
 
      call nbnd(ibc)
      call sbnd(ibc)
      call ebnd(ibc)
      call wbnd(ibc)

      return
      end

      subroutine nbnd(ibc)
c  sets northern boundary conditions

        include "./header_files/dimension.h"


!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
 
      parameter(nmodes=25)
      logical wclos,sclos,eclos,nclos

      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     1 rhon(imax,jmax),
     3 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax), 
     4 fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
     

      common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
      common/realbd/istarr(jmax),indarr(jmax)

      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      common/cntrl8/dt

	real*8 dt
c  ----------------------------------
c  set northern boundary conditions
c  ----------------------------------
        do 100 i=1,imax
        pn(i,jmax,kp)=pn(i,jmaxm,kp)
        if(nclos) then
          un(i,jmax,kp) =0.-un(i,jmaxm,kp)
          vn(i,jmaxm,kp)=0.
        else
          un(i,jmax,kp)=un(i,jmaxm,kp)
          vn(i,jmax,kp)=vn(i,jmax-2,kp)
        endif
100     continue
       return
      end

      subroutine sbnd(ibc)

        include "./header_files/dimension.h"

c  sets southern boundary conditions
!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
      parameter(nmodes=25)

      logical wclos,sclos,eclos,nclos

      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/amode/xpsi2(nmodes),cn(nmodes),zn(nmodes)

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     1 rhon(imax,jmax),
     3 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     5  fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
     

      common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
      common/realbd/istarr(jmax),indarr(jmax)

      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      common/cntrl8/dt
        
      dimension zeta(imax,2) 
	real*8 dt
c  --------------------------------
c  set southern boundary conditions 
c  --------------------------------

        do 100 i=1,imax
        pn(i,1,kp)=pn(i,2,kp)
        if(sclos) then    !closed
           un(i,1,kp)=0.-un(i,2,kp)
c          vn initially set to zero and never re-calculated
c          for  closed southern boundary
       else               !open
c         fortran code for v below
          un(i,1,kp)=un(i,2,kp)
        endif
100     continue        

      if(.not.sclos.and.ibc.eq.1) then ! calculate v for open s boundary

        j=1
        jp=2
        jpp=jp+1
        jm=jp

        do 200 i=2,imaxm
        im=i-1
        ip=i+1
        zeta(i,j) =(vn(ip,j,km)+vn(im,j,km)-2.*vn(i,j,km))*
     1    recdx2
     1        +(vn(i,jp,km)+vn(i,jm,km)-2.*vn(i,j,km))*recdy2
        zeta(i,jp)=(vn(ip,jp,km)+vn(im,jp,km)-2.*vn(i,jp,km))*
     1    recdx2
     1        +(vn(i,jpp,km)+vn(i,j,km)-2.*vn(i,jp,km))*recdy2

200     continue

        if(wclos)then
          zeta(1,j)=0.-zeta(2,j)
        else
          zeta(1,j)=zeta(2,j)
        endif 
        if(eclos)then
          zeta(imax,j)=0.-zeta(imaxm,j)
        else
          zeta(imax,j)=zeta(imaxm,j)
        endif 

        do 300 i=2,imaxm
        im=i-1
        ip=i+1

      vn(i,j,kp)=vn(i,j,km)+ tdt*rmask(i,j)*(
     1 tysav(i,j)*zn(nmm)/xpsi2(nmm)+visch*zeta(i,j)
     2 -0.25*fv(i,j)*(un(im,jp,k)+un(i,jp,k)+un(im,j,k)+
     3  un(i,j,k))
     3 - recdy*(pn(i,jp,k)-pn(i,j,k))-a*vn(i,j,km)/cn(nmm)**2)
300     continue
      endif
      return
      end

      subroutine ebnd(ibc)
c  sets eastern boundary conditions


        include "./header_files/dimension.h"

!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
      parameter(nmodes=25)

      logical wclos,sclos,eclos,nclos

      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     1 rhon(imax,jmax),
     3 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax), 
     5 fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
     

      common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
      common/realbd/istarr(jmax),indarr(jmax)

      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii,idate3, iday2,i11
      common/cntrl8/dt

	real*8 dt

c  -------------------------------
c  set eastern boundary conditions
c  -------------------------------

        do 100 j=1,jmax
        pn(imax,j,kp)=pn(imaxm,j,kp)
        if(eclos) then !   closed
          un(imax,j,kp)=0. 
          un(imaxm,j,kp)=0.
          vn(imax,j,kp)=0.-vn(imaxm,j,kp)
        else  !    open
          un(imax,j,kp)=un(imax-2,j,kp)
          un(imax,j,kp)=un(imaxm,j,kp)
        endif
100     continue

      return
      end

      subroutine wbnd(ibc)
c  sets western boundary conditions

        include "./header_files/dimension.h"


!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
      parameter(nmodes=25)

      logical wclos,sclos,eclos,nclos


      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/amode/xpsi2(nmodes),cn(nmodes),zn(nmodes)

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     1 rhon(imax,jmax),
     3 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax), 
     5 fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
      

      common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
      common/realbd/istarr(jmax),indarr(jmax)

      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      common/cntrl8/dt
      dimension zeta(2,jmax)
	real*8 dt
c  --------------------------------
c  set western boundary conditions
c  --------------------------------
        do 100 j=1,jmax
        pn(1,j,kp)=pn(2,j,kp)
        if(wclos) then ! closed
          un(1,j,kp)=0.
          vn(1,j,kp)=0.-vn(2,j,kp)
        else   !   open
c         fortran code for un below
          vn(1,j,kp)=vn(2,j,kp)
        endif
100     continue

        if (.not.wclos) then
        i=1
        ip=2
        ipp=ip+1
        im=ip

        do 120 j=2,jmaxm
        jp=j+1
        jm=j-1
        zeta(i,j)=(un(ip,j,km)+un(im,j,km)-2.*un(i,j,km))*
     1         recdx2
     1         +(un(i,jp,km)+un(i,jm,km)-2.*un(i,j,km))*recdy2
        zeta(ip,j)=(un(ipp,j,km)+un(i,j,km)-2.*un(ip,j,km))*
     1         recdx2
     1        +(un(ip,jp,km)+un(ip,jm,km)-2.*un(ip,j,km))*recdy2
120     continue

        if(nclos)then
          zeta(i,jmax)=0.-zeta(i,jmaxm)
        else
          zeta(i,jmax)=zeta(i,jmaxm)
        endif     
        if(sclos)then
          zeta(i,1)=0.-zeta(i,2)
        else
          zeta(i,1)=zeta(i,2)
        endif     
130     continue
        do 160 j=2,jmaxm
        jp=j+1
        jm=j-1
        un(i,j,kp)=un(i,j,km)+ tdt*rmask(i,j)*(
     1   txsav(i,j)*zn(nmm)/xpsi2(nmm)+ visch*zeta(i,j)
     1    -recdx*(pn(ip,j,k)
     2   -pn(i,j,k))-a*un(i,j,km)/cn(nmm)**2+
     3  0.25* fu(i,j)*(vn(i,j,k)+vn(ip,j,k)+vn(i,jm,k)+
     4  vn(ip,jm,k)))

160     continue
      endif
      return
      end



      subroutine sm2dt
c  this routine averages fields every nsmth time steps to eliminate
c  time-step splitting

        include "./header_files/dimension.h"


!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
      parameter(nmodes=25)

      logical wclos,sclos,eclos,nclos

      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

      common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     1 rhon(imax,jmax),
     3 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax), 
     5 fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
      

      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      common/cntrl8/dt

	real*8 dt
 
      do 100 j=1,jmax
      do 100 i=1,imax
      un(i,j,kp)=.5*(un(i,j,k)+un(i,j,kp))
      un(i,j,k )=.5*(un(i,j,k)+un(i,j,km))
      vn(i,j,kp)=.5*(vn(i,j,k)+vn(i,j,kp))
      vn(i,j,k )=.5*(vn(i,j,k)+vn(i,j,km))
      pn(i,j,kp)=.5*(pn(i,j,k)+pn(i,j,kp))
      pn(i,j,k )=.5*(pn(i,j,k)+pn(i,j,km))
100   continue
5     continue

      return
      end

       subroutine readmask(maskfile)

       include 'netcdf.inc' 

        include "./header_files/dimension.h"


c  !      parameter (imax=ikmax-1)
c  !      parameter (jmax=jkmax-1)


        common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1  month,istore,nsmth1,ndsave,ndend,ndstep,
     2  ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3  wclos,sclos,eclos,nclos

        common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1  wn(imax,jmax,3),pn(imax,jmax,3),
     1  rhon(imax,jmax),dampu(imax,jmax),damp(imax,jmax),
     4  taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     5  fv(imax,jmax),tau1y(imax,jmax,2),
     5  rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
     

        common/realbd/istarr(jmax),indarr(jmax)    ! included

        dimension tauxx(imax,jmax)
        common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
        common/pcom/rmsk(imax,jmax)
        common/zcom/zetun(imax,jmax),zetvn(imax,jmax)
        common/cntrl8/dt
        character*60  maskfile
       integer ncid_nmask,jjndex_ncep(3),xid_ncep,yid_ncep,zid_ncep
      real xx_ncep,yy_ncep,zz_ncep, phiy,phix,mode,step,tau0
	real*8 dt

       include "./header_files/param.h"

       status = nf_open(maskfile,nf_nowrite,ncid_nmask)
      if (status.ne.nf_noerr) stop 'error with nf_open(maskfile)'


        status = nf_inq_varid (ncid_nmask,'mask',xid_ncep)
      if (status.ne.nf_noerr) stop 'error: nf_inq_varid: rmask'


      do jj=1,jmax
                jjndex_ncep(2) =jj
                do ii=1,imax
                jjndex_ncep(1) = ii

      status = nf_get_var1_real(ncid_nmask,xid_ncep,jjndex_ncep,xx_ncep)
       if (status.ne.nf_noerr) stop 'error: nf_get_var1_u'

      rmask(ii,jj)=xx_ncep

           enddo
           enddo
       status=nf_close(ncid_nmask)
        if (status.ne.nf_noerr) stop 'error with close'

       return
        end

