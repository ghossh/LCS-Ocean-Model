         subroutine readwind(windfile,m1,m2,m3,rtime)

      character*80 windfile

       include 'netcdf.inc'
        include "./header_files/dimension.h"

      parameter (nmodes=1)
      integer kk,kk1,tt
      real uwnd(imax+1,jmax+1,2),vwnd(imax+1,jmax+1,2)
      integer iindex_ncep(3),starta(3),counta(3)
      real uu(imax+1,jmax+1,2),vv(imax+1,jmax+1,2),wwnd(imax+1,jmax+1,2)
      real u10(imax+1,jmax+1),v10(imax+1,jmax+1)
        
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
     

      common/cntrl3/msm, ijk, iii, idate3, iday2,i11
      common/cntrl2/ncid_ncep, uid_ncep, vid_ncep
      integer ii,jj,m11,m21,m31,uid_ncep,vid_ncep
      common/cntrl8/dt

	 real*8 rtime,rtime1,dt

       include "./header_files/param.h"
	 rtime1=rtime-ndbeg
        m11=m1-ndbeg
	 m21=m2-ndbeg
        m31=m3-ndbeg
        if (start_program.eq.1) then
        ndinc=ndrun+ndbeg
        else
        ndinc=ndrun
        endif

c      write(6,*)rtime1,rtime,ndbeg

c      write(6,*) windfile
      if (rtime1 .eq. 0.0) then

       status = nf_open(windfile,nf_nowrite,ncid_ncep)
      if (status.ne.nf_noerr) stop 'error with nf_open(windfile)'

      status = nf_inq_varid (ncid_ncep,'U10',uid_ncep)
      if (status.ne.nf_noerr) stop 'error: nf_inq_varid: u'

      status = nf_inq_varid (ncid_ncep,'V10',vid_ncep)
      if (status.ne.nf_noerr) stop 'error:nf_inq_varid: v'

      endif

        if(rtime1.eq.0.0)then  ! rtime if statement start

        kk=m11

        do tt=1,2

        starta(1)=1
        starta(2)=1
        starta(3)=kk+tt-1

        counta(1)=imax+1
        counta(2)=jmax+1
        counta(3)=1

       status = nf_get_vara_real(ncid_ncep,uid_ncep,starta,counta,u10)
       if (status.ne.nf_noerr) stop 'error: nf_get_var_u'
       status = nf_get_vara_real(ncid_ncep,vid_ncep,starta,counta,v10)
	  if (status.ne.nf_noerr) stop 'error with v'

	  do ii=1,imax+1
          do jj=1,jmax+1
          uu(ii,jj,tt)=u10(ii,jj)*100
          vv(ii,jj,tt)=v10(ii,jj)*100
          wwnd(ii,jj,tt)=sqrt((uu(ii,jj,tt)*uu(ii,jj,tt) +
     #    (vv(ii,jj,tt)*vv(ii,jj,tt))))
          enddo
       enddo

        do ii=1,imax
          do jj=1,jmax
        taux(ii,jj,tt)=(0.001175*0.0015*uu(ii,jj,tt)*wwnd(ii,jj,tt) +
     &  0.001175*0.0015*uu(ii+1,jj,tt)*wwnd(ii+1,jj,tt))/2.0
        tau1y(ii,jj,tt)=(0.001175*0.0015*vv(ii,jj,tt)*wwnd(ii,jj,tt) +
     &  0.001175*0.0015*vv(ii,jj+1,tt)*wwnd(ii,jj+1,tt))/2.0
          enddo
        enddo

        enddo
       write(6,*)u10(10,10)

!	write(42,22) 'f',ndsave,m1,m2,m3,m11,m21,m31,taux(8,24,1),
!     #  tau1y(8,24,1),taux(8,24,2),tau1y(8,24,2)

22	 format(a1,7i4,4f10.4)
24	 format(7i4,4f10.4)


         elseif(m11.gt.m31)then   ! rtime if statement else

        do  jj=1,jmax
        do  ii=1,imax

        taux(ii,jj,1)=taux(ii,jj,2)
        tau1y(ii,jj,1)=tau1y(ii,jj,2)
        enddo
        enddo

        kk1=m21

        starta(1)=1
        starta(2)=1
        starta(3)=kk1

        counta(1)=imax+1
        counta(2)=jmax+1
        counta(3)=1

        status = nf_get_vara_real(ncid_ncep,uid_ncep,starta,counta,u10)
	   if (status.ne.nf_noerr) stop 'error with u_read'
        status = nf_get_vara_real(ncid_ncep,vid_ncep,starta,counta,v10)
	   if (status.ne.nf_noerr) stop 'error with v_read'
      
   
	  do ii=1,imax+1
         do jj=1,jmax+1
         uu(ii,jj,1)=u10(ii,jj)*100
         vv(ii,jj,1)=v10(ii,jj)*100
         wwnd(ii,jj,1)=sqrt((uu(ii,jj,1)*uu(ii,jj,1) +
     #   (vv(ii,jj,1)*vv(ii,jj,1))))
         enddo
       enddo

        do ii=1,imax
         do jj=1,jmax
        taux(ii,jj,2)=(0.001175*0.0015*uu(ii,jj,1)*wwnd(ii,jj,1) +
     &  0.001175*0.0015*uu(ii+1,jj,1)*wwnd(ii+1,jj,1))/2.0
        tau1y(ii,jj,2)=(0.001175*0.0015*vv(ii,jj,1)*wwnd(ii,jj,1) +
     &  0.001175*0.0015*vv(ii,jj+1,1)*wwnd(ii,jj+1,1))/2.0
         enddo
        enddo
!	write(42,22) 's',ndsave,m1,m2,m3,m11,m21,m31,taux(8,24,1),
!     #  tau1y(8,24,1),taux(8,24,2),tau1y(8,24,2)
        endif  ! rtime if statement ends


      if (rtime.eq.ndinc)then
      status=nf_close(ncid_ncep)
       if (status.ne.nf_noerr) stop 'error with close'
      endif
      return
      end
