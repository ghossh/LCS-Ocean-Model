          subroutine pnbound

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

          integer hchoice(imax,jmax),hchi(8,16000), hchj(8,16000)
        integer hch(8)
        real dist,dx2,dy2
	real*8 dt

        dist=dx*dx+dy*dy
        dx2=dx*dx
        dy2=dy*dy
          do j=2,jmax-1
        do i=2,imax-1
        hchoice(i,j)=0
        if(rmask(i,j).eq.0.and.rmask(i-1,j).eq.1) hchoice(i,j)=1
        if(rmask(i,j).eq.0.and.rmask(i+1,j).eq.1) hchoice(i,j)=2
        if(rmask(i,j).eq.0.and.rmask(i,j-1).eq.1) hchoice(i,j)=3
        if(rmask(i,j).eq.0.and.rmask(i,j+1).eq.1) hchoice(i,j)=4
        if(rmask(i,j).eq.0.and.rmask(i,j+1).eq.1.and.
     *   rmask(i-1,j).eq.1) hchoice(i,j)=5
        if(rmask(i,j).eq.0.and.rmask(i+1,j).eq.1.and.
     *   rmask(i,j-1).eq.1) hchoice(i,j)=6
        if(rmask(i,j).eq.0.and.rmask(i-1,j).eq.1.and.
     *   rmask(i,j-1).eq.1) hchoice(i,j)=7
        if(rmask(i,j).eq.0.and.rmask(i+1,j).eq.1.and.
     *   rmask(i,j+1).eq.1) hchoice(i,j)=8
        enddo
        enddo

        do ll=1,8
        hch(ll)=0
        enddo

        do j=2, jmax-1
         do i=2, imax-1
          if (hchoice(i,j) .eq. 1) then
           hch(1) =hch(1)+1
           hchi(1,hch(1))=i
           hchj(1,hch(1))=j
          endif
          if (hchoice(i,j) .eq. 2) then
           hch(2) =hch(2)+1
           hchi(2,hch(2))=i
           hchj(2,hch(2))=j
          endif
          if (hchoice(i,j) .eq. 3) then
           hch(3) =hch(3)+1
           hchi(3,hch(3))=i
           hchj(3,hch(3))=j
          endif
          if (hchoice(i,j) .eq. 4) then
           hch(4) =hch(4)+1
           hchi(4,hch(4))=i
           hchj(4,hch(4))=j
          endif
          if (hchoice(i,j) .eq. 5) then
           hch(5) =hch(5)+1
           hchi(5,hch(5))=i
           hchj(5,hch(5))=j
          endif
          if (hchoice(i,j) .eq. 6) then
           hch(6) =hch(6)+1
           hchi(6,hch(6))=i
           hchj(6,hch(6))=j
          endif
          if (hchoice(i,j) .eq. 7) then
           hch(7) =hch(7)+1
           hchi(7,hch(7))=i
           hchj(7,hch(7))=j
          endif
          if (hchoice(i,j) .eq. 8) then
           hch(8) =hch(8)+1
           hchi(8,hch(8))=i
           hchj(8,hch(8))=j
          endif
        enddo
        enddo

        if (sclos)  then
         else
        do i=1,imax
         hch(4) =hch(4)+1
         hchi(4,hch(4))=i
         hchj(4,hch(4))=1
         enddo
         endif

c       northern box bound
         if (nclos)  then
         else
           do i=1, imax
           hch(3) =hch(3)+1
           hchi(3,hch(3))=i
           hchj(3,hch(3))=jmax
           enddo
          endif

c       eastern box boun
         if (eclos)  then
         else
         do j=1,jmax
           hch(1) =hch(1)+1
           hchi(1,hch(1))=imax
           hchj(1,hch(1))=j
           enddo
          endif
c       c       western box boun
         if (wclos)  then
         else
         do j=1,jmax
           hch(2) =hch(2)+1
           hchi(2,hch(2))=1
           hchj(2,hch(2))=j
           enddo
          endif
        do ll=1,8
        if (hch(ll).gt. 16000 ) then
        write(*,*) 'increase dimension of hch '
        stop
        endif
        enddo

c         write(*,*) dx2,'*******'
! Put the boundary condition on pn


        do ii=1,hch(1)
        i=hchi(1,ii)
        j=hchj(1,ii)
        pn(i,j,kp)=pn(i-1,j,kp)
        enddo
        do ii=1,hch(2)
        i=hchi(2,ii)
        j=hchj(2,ii)
        pn(i,j,kp)=pn(i+1,j,kp)
        enddo
        do ii=1,hch(3)
        i=hchi(3,ii)
        j=hchj(3,ii)
        pn(i,j,kp)=pn(i,j-1,kp)
        enddo
        do ii=1,hch(4)
        i=hchi(4,ii)
        j=hchj(4,ii)
        pn(i,j,kp)=pn(i,j+1,kp)
        enddo
        do ii=1,hch(5)
        i=hchi(5,ii)
        j=hchj(5,ii)
        pn(i,j,kp)=((pn(i,j+1,kp)*dx2)+(pn(i-1,j,kp)*dy2))/dist
        enddo
        do ii=1,hch(6)
        i=hchi(6,ii)
        j=hchj(6,ii)
        pn(i,j,kp)=((pn(i,j-1,kp)*dx2)+(pn(i+1,j,kp)*dy2))/dist
        enddo
        do ii=1,hch(7)
         i=hchi(7,ii)
        j=hchj(7,ii)
        pn(i,j,kp)=((pn(i,j-1,kp)*dx2)+(pn(i-1,j,kp)*dy2))/dist
        enddo
        do ii=1,hch(8)
        i=hchi(8,ii)
        j=hchj(8,ii)
        pn(i,j,kp)=((pn(i,j+1,kp)*dx2)+(pn(i+1,j,kp)*dy2))/dist
        enddo

c        write(*,*)hchoice(1,1),'####'
        return
        end

                                                                                          
