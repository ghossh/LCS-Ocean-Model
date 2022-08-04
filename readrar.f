         subroutine readrar(iopen)

       include 'netcdf.inc'
        include "./header_files/dimension.h"

         parameter (nmodes=25)
         integer kk,kk1,kk11,kk2
         integer iindex_ncep(3)
         real uu,vv,phiy,phix,mode,step,tau0,pp

      common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1 month,istore,nsmth1,ndsave,ndend,ndstep,
     2 ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3 wclos,sclos,eclos,nclos

       common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1 wn(imax,jmax,3),pn(imax,jmax,3),
     2 rhon(imax,jmax),
     3 dampu(imax,jmax),damp(imax,jmax),
     4 taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     5 fv(imax,jmax),tau1y(imax,jmax,2),
     5 rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
     

      common/cntrl3/msm, ijk, iii, idate3, iday2,i11
      common/cntrl6/ncid2,uid2,vid2,pid2
      common/cntrl8/dt
	integer ii,jj,uid2,vid2,pid2
        real*8 dt
        character*8 dir
        character*10 ext
        character*20 res_fil
        character*2 text

       include "./header_files/param.h"
       write(text,106) nmm
106     format(i2.2)
        dir='./input/'
        ext='restart.nc'
        res_fil=dir//text//ext

       status = nf_open(res_fil,nf_nowrite,ncid2)
      if (status.ne.nf_noerr) stop 'error with nf_open(restartfile)'

      status = nf_inq_varid (ncid2,'un',uid2)
      if (status.ne.nf_noerr) stop 'error: nf_inq_varid: u'

      status = nf_inq_varid (ncid2,'vn',vid2)
      if (status.ne.nf_noerr) stop 'error:nf_inq_varid: v'

      status = nf_inq_varid (ncid2,'pn',pid2)
      if (status.ne.nf_noerr) stop 'error:nf_inq_varid: v'
      
       do it=1,3
         iindex_ncep(3) =it
       do jj=1,jmax
         iindex_ncep(2) =jj
         do ii=1,imax
          iindex_ncep(1) = ii

       status = nf_get_var1_real(ncid2,uid2,iindex_ncep,uu)
       if (status.ne.nf_noerr) stop 'error: nf_get_var1_u2'

       status = nf_get_var1_real(ncid2,vid2,iindex_ncep,vv)
       if (status.ne.nf_noerr) stop 'error: nf_get_var1_v'

       status = nf_get_var1_real(ncid2,pid2,iindex_ncep,pp)
       if (status.ne.nf_noerr) stop 'error: nf_get_var1_v'

       un(ii,jj,it)=uu
       vn(ii,jj,it)=vv
       pn(ii,jj,it)=pp

        enddo
        enddo
        enddo

       status=nf_close(ncid2)
       if (status.ne.nf_noerr) stop 'error with close'

	
      return
      end
