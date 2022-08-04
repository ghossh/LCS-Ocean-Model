        subroutine vngeneral

        include "./header_files/dimension.h"

        parameter(nmodes=25)

        common/cntrl/dx,dy,a,g,visch,tdt,xlen,
     1  month,istore,nsmth1,ndsave,ndend,ndstep,
     2  ndbeg,nend,nprint,nstore,ismth,island,kp,k,km,nmm,
     3  wclos,sclos,eclos,nclos


        common/amode/xpsi2(nmodes),cn(nmodes),zn(nmodes)

        common/var/un(imax,jmax,3),vn(imax,jmax,3),
     1  wn(imax,jmax,3),pn(imax,jmax,3),
     2  rhon(imax,jmax),
     3  dampu(imax,jmax),damp(imax,jmax),
     4  taux(imax,jmax,2),tauy(imax,jmax,2),fu(imax,jmax),
     4  fv(imax,jmax),tau1y(imax,jmax,2),
     5  rmask(imax,jmax),txsav(imax,jmax),tysav(imax,jmax)
     

        common/com1/imaxm,jmaxm,ijmax,recdx,recdx2,recdy,recdy2
        common/zcom/zetun(imax,jmax),zetvn(imax,jmax)
        logical wclos,sclos,eclos,nclos
        common/pcom/rmsk(imax,jmax)
        common/realbd/istarr(jmax),indarr(jmax)
        common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
        common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
        common/cntrl8/dt
        integer  vchoice(imax,jmax), vchi(4,16000), vchj(4,16000)
        integer vch(4)
	real*8 dt

         do j=2,jmax-1
        do i=2,imax-1
        vchoice(i,j)=0
        if(rmask(i,j).eq.1.and.rmask(i,j+1).eq.0) vchoice(i,j)=1
c        write(131,*)'first',i,j,vchoice(i,j)
        if(rmask(i,j).eq.0.and.rmask(i,j+1).eq.1) vchoice(i,j)=2
        if(rmask(i,j).eq.0.and.rmask(i+1,j).eq.1.and.
     *  rmask(i,j+1).eq.0) vchoice(i,j)=3
        if(rmask(i,j).eq.0.and.rmask(i-1,j).eq.1.and.
     *  rmask(i,j+1).eq.0)  vchoice(i,j)=4
        enddo
        enddo

 
        do ll=1,4
        vch(ll)=0
        enddo
c==================================================================
c put the general boundary condition on solution

c=====================================================================
        do j=2, jmax-1
         do i=2, imax-1
          if (vchoice(i,j) .eq. 1) then
           vch(1) =vch(1)+1
           vchi(1,vch(1))=i
           vchj(1,vch(1))=j
          endif
          if (vchoice(i,j) .eq. 2) then
           vch(2) =vch(2)+1
           vchi(2,vch(2))=i
           vchj(2,vch(2))=j
          endif
          if (vchoice(i,j) .eq. 3) then
           vch(3) =vch(3)+1
           vchi(3,vch(3))=i
           vchj(3,vch(3))=j
          endif
          if (vchoice(i,j) .eq. 4) then
           vch(4) =vch(4)+1
           vchi(4,vch(4))=i
           vchj(4,vch(4))=j
          endif
          enddo
          enddo

c       southern box bound
         if (sclos)  then
         else
          do i=1,imax
           vch(2) =vch(2)+1
           vchi(2,vch(2))=i
           vchj(2,vch(2))=1
           enddo
           endif

c       northern box boun
         if (nclos)  then
         else
           do i=1, imax
           vch(1) =vch(1)+1
           vchi(1,vch(1))=i
           vchj(1,vch(1))=jmax
           enddo
          endif
c       eastern box boun
         if (eclos)  then
         else
         do j=1,jmax
           vch(4) =vch(4)+1
           vchi(4,vch(4))=imax
           vchj(4,vch(4))=j
           enddo
          endif

c       western box boun
         if (wclos)  then
         else
         do j=1,jmax
           vch(3) =vch(3)+1
           vchi(3,vch(3))=1
           vchj(3,vch(3))=j
           enddo
           endif

        do ll=1,4
        if (vch(ll).gt. 16000 ) then
        write(*,*) 'increase dimension of vch '
        stop
        endif
        enddo
        do ii=1,vch(1)
        i=vchi(1,ii)
        j=vchj(1,ii)
        vn(i,j,kp)=0
        enddo
        do ii=1,vch(2)
        i=vchi(2,ii)
        j=vchj(2,ii)
        vn(i,j,kp)=0
        enddo
        do ii=1,vch(3)
        i=vchi(3,ii)
        j=vchj(3,ii)
        vn(i,j,kp)=-vn(i+1,j,kp)
        enddo
        do ii=1,vch(4)
        i=vchi(4,ii)
        j=vchj(4,ii)
        vn(i,j,kp)=-vn(i-1,j,kp)
        enddo
       
        return
        end
                                          
