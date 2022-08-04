         subroutine unspecial( splons,splone,splats,splate)
        include "./header_files/dimension.h"


!      parameter (imax=ikmax-1)
!      parameter (jmax=jkmax-1)
      parameter(nmodes=25)

      common/amode/xpsi2(nmodes),cn(nmodes),zn(nmodes)

      logical wclos,sclos,eclos,nclos
      common/cntrl/dx,dy,a,g,visch,dt,tdt,xlen,
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
      integer splons,splone,splats,splate
c      real splons1,splone1,splats1,splate1
        integer uchoice(imax,jmax),  uchsi(4,16000), uchsj(4,16000)
        integer uchs(4)
        real capg
        do j=2,jmax-1
        do i=2,imax-1
        uchoice(i,j)=0
        if(rmask(i,j).eq.1.and.rmask(i+1,j).eq.0) uchoice(i,j)=1
        if(rmask(i,j).eq.0.and.rmask(i+1,j).eq.1) uchoice(i,j)=2
        if(rmask(i,j).eq.0.and.rmask(i,j+1).eq.1.and.
     *   rmask(i+1,j).eq.0)uchoice(i,j)=3
        if(rmask(i,j).eq.0.and.rmask(i,j-1).eq.1.and.
     *   rmask(i+1,j).eq.0) uchoice(i,j)=4
        enddo
        enddo

         do ll=1,4
        uchs(ll)=0
        enddo
         

c------------------------------------------
c       put the special boundary condition

c-------------------------------------------------------
        do j=splats,splate
         do i=splons,splone
c           do  j=178,245
c           do   i=235,260
          if (uchoice(i,j) .eq. 1) then
           uchs(1) =uchs(1)+1
           uchsi(1,uchs(1))=i
           uchsj(1,uchs(1))=j

c           write(131,*)'first',i,j,uchoice(i,j),uchsi(1,uchs(1))
          endif
          if (uchoice(i,j) .eq. 2) then
           uchs(2) =uchs(2)+1
           uchsi(2,uchs(2))=i
           uchsj(2,uchs(2))=j
c          write(131,*)'first',i,j,uchoice(i,j),uchsi(2,uchs(2))
          endif
          if (uchoice(i,j) .eq. 3) then
           uchs(3) =uchs(3)+1
           uchsi(3,uchs(3))=i
           uchsj(3,uchs(3))=j
          endif
          if (uchoice(i,j) .eq. 4) then
           uchs(4) =uchs(4)+1
           uchsi(4,uchs(4))=i
           uchsj(4,uchs(4))=j
          endif
          enddo
          enddo
c--------------------

        do ll=1,4
        if (uchs(ll).gt. 16000 ) then
        write(*,*) 'increase dimension of uch '
        stop
        endif
        enddo

 
         do ii=1,uchs(1)
        i=uchsi(1,ii)
        j=uchsj(1,ii)
        capg=.25*(tysav(i,j)+tysav(i,j-1)
     #   +tysav(i+1,j)+tysav(i+1,j-1))
     .  /xpsi2(nmm)
         un(i,j,kp)=capg/fu(i,j)
c         write(137,*)'second',i,j,capg,fu(i,j),un(i,j,kp)
         enddo

         do  ii=1,uchs(2)
        i=uchsi(2,ii)
        j=uchsj(2,ii)
        capg=.25*(tysav(i,j)+tysav(i+1,j)
     #   +tysav(i,j-1)+tysav(i+1,j-1))
     .   /xpsi2(nmm)
         un(i,j,kp)=capg/fu(i,j)
       write(141,*)'second',i,j,capg,fu(i,j), un(i,j,kp)
c         write(141,*)'second',i,j,xpsi2(1)  
         enddo

         do ii=1,uchs(3)
        i=uchsi(3,ii)
        j=uchsj(3,ii)
        un(i,j,kp)=-un(i,j+1,kp)
        enddo

        do ii=1,uchs(4)
        i=uchsi(4,ii)
        j=uchsj(4,ii)
        un(i,j,kp)=-un(i,j-1,kp)
        enddo



        return
        end


