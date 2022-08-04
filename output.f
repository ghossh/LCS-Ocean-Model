          subroutine output(iopen, outfil,nt,torigin)


        include 'netcdf.inc'
        include "./header_files/dimension.h"


      integer ncid, londim, latdim, timedim
      integer lonid, latid, timeid,nt,iret
      integer  ntime, undim(3), unid,vnid,pnid,i,j
      integer starta(3),counta(3)
      real time, x,y
      real delta,slat,elat,slon,elon
      character*20 torigin
      character*15 outfil

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

      common/cntrl8/dt
       dimension  var(imax,jmax),var1(imax,jmax),var2(imax,jmax)
       real   misg, valid_range(2)
	real*8 dt
      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      common/cntrl5/ncid, unid, vnid,pnid

      include "./header_files/param.h"
       if (start_program.eq.1) then
        ndinc=ndrun+ndbeg
        else
        ndinc=ndrun
        endif

       
      if (iopen .eq. 1) then

        write(6,*)outfil

       data ntorig / 20 /
!      status = nf_create(outfil,nf_share,ncid)

      iret=nf_create(outfil,OR(NF_CLOBBER,NF_64BIT_OFFSET),ncid)
      if (iret.ne.nf_noerr) stop 'error with nf_open(outfil)'

      iret = nf_def_dim(ncid,'lon',imax,londim)
      if (iret.ne.nf_noerr) stop 'error with lon'

      iret = nf_def_dim(ncid,'lat',jmax,latdim)
      if (iret.ne.nf_noerr) stop 'error with lat'

      iret = nf_def_dim(ncid,'time',ncunlim,timedim)
      if (iret.ne.nf_noerr) stop 'error with time'

      iret=nf_def_var(ncid,'lon',nf_float,1,londim,lonid)
      if (iret.ne.nf_noerr) stop 'error with lon1'
     
      iret=nf_put_att_text(ncid,lonid,'units',11,'degree east')
      if (iret.ne.nf_noerr) stop 'error with lonunit'

      iret=nf_def_var(ncid,'lat',nf_float,1,latdim,latid)
      if (iret.ne.nf_noerr) stop 'error with lat1'
   
      iret=nf_put_att_text(ncid,latid,'units',12,'degree north')
      if (iret.ne.nf_noerr) stop 'error with latunit'

      iret=nf_def_var(ncid,'time',nf_float,1,timedim,timeid)
      if (iret.ne.nf_noerr) stop 'error with time1'

      undim(1)=londim
      undim(2)=latdim
      undim(3)=timedim

      iret=nf_def_var(ncid,'un',nf_float,3,undim,unid)
      if (iret.ne.nf_noerr) stop 'error with u_output'

      iret=nf_def_var(ncid,'vn',nf_float,3,undim,vnid)
      if (iret.ne.nf_noerr) stop 'error with v'

      iret=nf_def_var(ncid,'pn',nf_float,3,undim,pnid)
      if (iret.ne.nf_noerr) stop 'error with p'
      
      iret=nf_put_att_text(ncid,unid,'long_name',1,'U')
      if (iret.ne.nf_noerr) stop 'error with un1'

      iret=nf_put_att_text(ncid,vnid,'long_name',1,'V')
      if (iret.ne.nf_noerr) stop 'error with vn1'

      iret=nf_put_att_text(ncid,pnid,'long_name',1,'P')
      if (iret.ne.nf_noerr) stop 'error with pn1'
      
       iret=nf_put_att_text(ncid,unid,'units',4,'cm/s')
      if (iret.ne.nf_noerr) stop 'error with u unit'
      
       iret=nf_put_att_text(ncid,vnid,'units',4,'cm/s')
      if (iret.ne.nf_noerr) stop 'error with v unit'
    
      iret=nf_put_att_text(ncid,pnid,'units',9,'dyne/cm^2')
      if (iret.ne.nf_noerr) stop 'error with p unit'

      iret=nf_put_att_text(ncid,unid,'FORTRAN_format',
     + 5,'f15.8')
      if (iret.ne.nf_noerr) stop 'error with u unit'

      iret=nf_put_att_text(ncid,vnid,'FORTRAN_format',
     + 5,'f15.8')
      if (iret.ne.nf_noerr) stop 'error with v unit'

      iret=nf_put_att_text(ncid,pnid,'FORTRAN_format',
     + 5,'f15.8')
      if (iret.ne.nf_noerr) stop 'error with p unit'

      iret=nf_put_att_text(ncid,timeid,'long_name',4,'time')
      if (iret.ne.nf_noerr) stop 'error with time2'

      iret=nf_put_att_text(ncid,timeid,'units',3,'day')
      if (iret.ne.nf_noerr) stop 'error with time4'
        
      iret=nf_put_att_text(ncid,timeid,'time_origin',ntorig,torigin)
      if (iret.ne.nf_noerr) stop 'error with time2'

      iret=nf_enddef(ncid)
      if (iret.ne.nf_noerr) stop 'error with endef'

       do i= 1,imax
        x = slon+delta*float(i-1)
       iret= nf_put_var1_real(ncid,lonid,i,x)
       if (iret.ne.nf_noerr) stop 'error with lon2'
       enddo

       do j= 1,jmax
        y = slat+delta*float(j-1)
        iret= nf_put_var1_real(ncid,latid,j,y)
        if (iret.ne.nf_noerr) stop 'error with lat2'
       enddo

      endif

      if (iopen .eq. 2) then

      iret = nf_open(outfil,nf_write,ncid)
      if (iret.ne.nf_noerr) stop 'error with open'

      ntime=nt
      time=float(ndsave)
c      write(*,*) time, ndsave
         misg   = 9.969E+36
         valid_range(1) = -1.000E+06
         valid_range(2) = 1.000E+06

      iret=nf_inq_dimid(ncid,'time',timeid)
      if (iret.ne.nf_noerr)  stop 'error with time'

      iret=nf_inq_varid(ncid,'time',timeid)
      if (iret.ne.nf_noerr) stop 'error with time1'

      iret=nf_inq_varid(ncid,'un',unid)
      if (iret.ne.nf_noerr) stop 'error with un'

      iret=nf_inq_varid(ncid,'vn',vnid)
      if (iret.ne.nf_noerr) stop 'error with vn'

      iret=nf_inq_varid(ncid,'pn',pnid)
      if (iret.ne.nf_noerr) stop 'error with pn'

      iret=nf_put_var1_real(ncid,timeid,ntime,time)
      if (iret.ne.nf_noerr) stop 'error with time2'
  
        starta(1)=1
        starta(2)=1
        starta(3)=ntime
        counta(1)=imax 
        counta(2)=jmax
        counta(3)=1
       do j=1,jmax
        do i=1,imax
         if (rmask(i,j) .eq. 0.0) then
          var(i,j) =misg 
          var1(i,j)=misg 
          var2(i,j)=misg 
          elseif (rmask(i,j) .gt. 0.5) then
          if (i .gt. 1) then
           var(i,j) = 0.5*(un(i-1,j,kp)+un(i,j,kp))
          else
           var(i,j)=2.0*un(1,j,kp)-0.5*(un(1,j,kp)+un(2,j,kp))
          endif
          if (j .gt. 1) then
           var1(i,j) = 0.5*(vn(i,j-1,kp)+vn(i,j,kp))
          else
           var1(i,j)=2.0*vn(i,1,kp)-0.5*(vn(i,1,kp)+vn(i,2,kp))
          endif
          var2(i,j)=pn(i,j,kp) 
	  endif
          enddo
          enddo
         iret=nf_put_vara_real(ncid,unid,starta,counta,var)
         if (iret.ne.nf_noerr) stop 'error with un2'
         iret=nf_put_vara_real(ncid,vnid,starta,counta,var1)
         if (iret.ne.nf_noerr) stop 'error with vn2'
         iret=nf_put_vara_real(ncid,pnid,starta,counta,var2)
         if (iret.ne.nf_noerr) stop 'error with pn2'
           endif
       iret=nf_close(ncid)
       if (iret.ne.nf_noerr) stop 'error with close'

      return
      end

