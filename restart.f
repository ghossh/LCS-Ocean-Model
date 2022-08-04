          subroutine restart(iopen,torigin)


        include 'netcdf.inc'
        include "./header_files/dimension.h"


      integer status, londim, latdim, timedim
      integer  lonid, latid, timeid,mid
      integer  undim(3), i,j,k11,ncdi1,unid1,vnid1,pnid1
      integer counta(3),starta(3)

      real time, x, y
      character*20  torigin
      character*40 resfil
      real delta,slat,elat,slon,elon

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


       dimension  var(imax,jmax),var1(imax,jmax),var2(imax,jmax)
       real   misg, valid_range(2),uu,vv,pp
       common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
       common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
       common/cntrl4/ncid1, unid1, vnid1,pnid1
       common/cntrl8/ dt
	real*8 dt
        character*10 dir
        character*10 ext
        character*28 res_fil
        character*2 text

      include "./header_files/param.h"
      if (start_program.eq.1) then
        ndinc=ndrun+ndbeg
        else
        ndinc=ndrun
        endif

       write(text,106) nmm
106     format(i2.2)
                dir='./restart/'
        ext='restart.nc'
        res_fil=dir//text//ext
          
!      status = nf_create('restart_write.nc',nf_share,ncid1)
!      if (status.ne.nf_noerr) stop 'error with nf_open(outfil)'

!	torigin='01-JAN-2000 00:00:00'

      iret=nf_create(res_fil,OR(NF_CLOBBER,NF_64BIT_OFFSET),ncid1)
      if (iret.ne.nf_noerr) stop 'error with nf_open(outfil)'

      status = nf_def_dim(ncid1,'lon',imax,londim)
      if (status.ne.nf_noerr) stop 'error with lon'
    
      status = nf_def_dim(ncid1,'lat',jmax,latdim)
      if (status.ne.nf_noerr) stop 'error with lat'

      status = nf_def_dim(ncid1,'time',ncunlim,timedim)
      if (status.ne.nf_noerr) stop 'error with time'

      status=nf_def_var(ncid1,'lon',nf_float,1,londim,lonid)
      if (status.ne.nf_noerr) stop 'error with lon1'
     
      status=nf_put_att_text(ncid1,lonid,'units',11,'degree east')
      if (status.ne.nf_noerr) stop 'error with lonunit'

      status=nf_def_var(ncid1,'lat',nf_float,1,latdim,latid)
      if (status.ne.nf_noerr) stop 'error with lat1'
   
      status=nf_put_att_text(ncid1,latid,'units',12,'degree north')
      if (status.ne.nf_noerr) stop 'error with latunit'

      status=nf_def_var(ncid1,'time',nf_float,1,timedim,timeid)
      if (status.ne.nf_noerr) stop 'error with time1'

      undim(1)=londim
      undim(2)=latdim
      undim(3)=timedim

      status=nf_def_var(ncid1,'un',nf_float,3,undim,unid1)
      if (status.ne.nf_noerr) stop 'error with u_restart'

      status=nf_def_var(ncid1,'vn',nf_float,3,undim,vnid1)
      if (status.ne.nf_noerr) stop 'error with v'

      status=nf_def_var(ncid1,'pn',nf_float,3,undim,pnid1)
      if (status.ne.nf_noerr) stop 'error with p'
      
      status=nf_put_att_text(ncid1,unid1,'long_name',1,'U')
      if (status.ne.nf_noerr) stop 'error with un1'

      status=nf_put_att_text(ncid1,vnid1,'long_name',1,'V')
      if (status.ne.nf_noerr) stop 'error with vn1'

      status=nf_put_att_text(ncid1,pnid1,'long_name',1,'P')
      if (status.ne.nf_noerr) stop 'error with pn1'
      
       status=nf_put_att_text(ncid1,unid1,'units',4,'cm/s')
      if (status.ne.nf_noerr) stop 'error with u unit'
      
       status=nf_put_att_text(ncid1,vnid1,'units',4,'cm/s')
      if (status.ne.nf_noerr) stop 'error with v unit'
    
      status=nf_put_att_text(ncid1,pnid1,'units',9,'dyne/cm^2')
      if (status.ne.nf_noerr) stop 'error with p unit'

      status=nf_put_att_text(ncid1,unid1,'FORTRAN_format',
     +    4,'f15.8')
      if (status.ne.nf_noerr) stop 'error with p unit'

      status=nf_put_att_text(ncid1,vnid1,'FORTRAN_format',
     +    4,'f15.8')
      if (status.ne.nf_noerr) stop 'error with p unit'

      status=nf_put_att_text(ncid1,pnid1,'FORTRAN_format',
     +    4,'f15.8')
      if (status.ne.nf_noerr) stop 'error with p unit'

      status=nf_put_att_text(ncid1,timeid,'long_name',4,'time')
      if (status.ne.nf_noerr) stop 'error with time2'

      status=nf_put_att_text(ncid1,timeid,'units',3,'day')
      if (status.ne.nf_noerr) stop 'error with time3'
        
      status=nf_put_att_text(ncid1,timeid,'time_origin',20,torigin)
      if (status.ne.nf_noerr) stop 'error with time3'

      status=nf_enddef(ncid1)
      if (status.ne.nf_noerr) stop 'error with endef'

         do i= 1,imax
          x = slon+delta*float(i-1)
          status= nf_put_var1_real(ncid1,lonid,i,x)
          if (status.ne.nf_noerr) stop 'error with lon2'
         enddo

         do j= 1,jmax
          y = slat+delta*float(j-1)
          status= nf_put_var1_real(ncid1,latid,j,y)
          if (status.ne.nf_noerr) stop 'error with lat2'
         enddo

	 do it=1,3
         time=float(it)
        status=nf_put_var1_real(ncid1,timeid,it,time)
        if (status.ne.nf_noerr) stop 'error with time21'
 	 enddo

        misg   = 0.0
        valid_range(1) = -1.000E+06
        valid_range(2) = 1.000E+06

        status=nf_inq_varid(ncid1,'un',unid1)
        if (status.ne.nf_noerr) stop 'error with un'

        status=nf_inq_varid(ncid1,'vn',vnid1)
        if (status.ne.nf_noerr) stop 'error with vn'

        status=nf_inq_varid(ncid1,'pn',pnid1)
        if (status.ne.nf_noerr) stop 'error with pn'

	do it=1,3
	starta(1)=1
        starta(2)=1
        starta(3)=it

        counta(3)=1
        counta(1)=imax
        counta(2)=jmax

        do i=1,imax
         do j=1,jmax
	var(i,j) =un(i,j,it)
        var1(i,j)=vn(i,j,it)
        var2(i,j)=pn(i,j,it)
	 enddo
	 enddo

        status=nf_put_vara_real(ncid1,unid1,starta,counta,var)
	if (iret.ne.nf_noerr) stop 'error with un3'
        status=nf_put_vara_real(ncid1,vnid1,starta,counta,var1)
	if (iret.ne.nf_noerr) stop 'error with un3'
        status=nf_put_vara_real(ncid1,pnid1,starta,counta,var2)
	if (iret.ne.nf_noerr) stop 'error with un3'


          enddo

       status=nf_close(ncid1)
       if (status.ne.nf_noerr) stop 'error with close'

      return
      end

