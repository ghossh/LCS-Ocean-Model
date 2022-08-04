         subroutine ungeneral
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

      common/pcom/rmsk(imax,jmax)
      common/realbd/istarr(jmax),indarr(jmax)
      common/cntrl1/nchoice,ndmon(12),ndyr,iyear,oyear,iday
      common/cntrl3/ msm, ijk, iii, idate3, iday2,i11
      common/cntrl8/dt
       
        integer uchoice(imax,jmax), uchi(4,16000), uchj(4,16000)
        integer uch(4)
	real*8 dt
        
        do j=2,jmax-1
        do i=2,imax-1
        uchoice(i,j)=0
        if(rmask(i,j).eq.1.and.rmask(i+1,j).eq.0) uchoice(i,j)=1
        if(rmask(i,j).eq.0.and.rmask(i+1,j).eq.1) uchoice(i,j)=2
        if(rmask(i,j).eq.0.and.rmask(i,j+1).eq.1.and.
     *   rmask(i+1,j).eq.0)uchoice(i,j)=3
        if(rmask(i,j).eq.0.and.rmask(i,j-1).eq.1.and.
     *   rmask(i+1,j).eq.0) uchoice(i,j)=4
c         write(25,*) rmask(i,j),i,j
        enddo
        enddo


         do ll=1,4
        uch(ll)=0
        enddo
c-----------------------------------------------------------
c put general bounday condition

c-------------------------------------------------------------------------
         do j=2, jmax-1
         do i=2, imax-1
          if (uchoice(i,j) .eq. 1) then
           uch(1) =uch(1)+1
           uchi(1,uch(1))=i
           uchj(1,uch(1))=j
          endif
          if (uchoice(i,j) .eq. 2) then
           uch(2) =uch(2)+1
           uchi(2,uch(2))=i
           uchj(2,uch(2))=j
          endif
          if (uchoice(i,j) .eq. 3) then
           uch(3) =uch(3)+1
           uchi(3,uch(3))=i
           uchj(3,uch(3))=j
          endif
          if (uchoice(i,j) .eq. 4) then
           uch(4) =uch(4)+1
           uchi(4,uch(4))=i
           uchj(4,uch(4))=j
          endif
          enddo
          enddo

c       southern box bound
         if (sclos)  then
         else
          do i=1,imax
           uch(3) =uch(3)+1
           uchi(3,uch(3))=i
           uchj(3,uch(3))=1
          enddo
          endif

c       northern box boundary
         if (nclos)  then
         else
           do i=1, imax
           uch(4)=uch(4)+1
           uchi(4,uch(4))=i
           uchj(4,uch(4))=jmax
           enddo
           endif
c           eastern box boun
         if (eclos)  then
         else
         do j=1,jmax
           uch(1) =uch(1)+1
           uchi(1,uch(1))=imax
           uchj(1,uch(1))=j
           enddo
          endif

c       western box boun
         if (wclos)  then
         else
         do j=1,jmax
           uch(2) =uch(2)+1
           uchi(2,uch(2))=1
           uchj(2,uch(2))=j
           enddo
          endif

        do ll=1,4
        if (uch(ll) .gt. 16000 ) then
        write(*,*) 'increase dimension of uch '
        stop
        endif
        enddo

        do ii=1,uch(1)
        i=uchi(1,ii)
        j=uchj(1,ii)
        un(i,j,kp)=0
c         write(139,*) 'first', i,j,uchoice(i,j),un(i,j,kp)
        enddo
        do ii=1,uch(2)
        i=uchi(2,ii)
        j=uchj(2,ii)
        un(i,j,kp)=0
        enddo
        do ii=1,uch(3)
        i=uchi(3,ii)
        j=uchj(3,ii)
        un(i,j,kp)=-un(i,j+1,kp)
c         write(139,*) 'first', i,j,uchoice(i,j),un(i,j,1)
        enddo
        do ii=1,uch(4)
        i=uchi(4,ii)
        j=uchj(4,ii)
        un(i,j,kp)=-un(i,j-1,kp)
        enddo
                 
        return
        end
      
