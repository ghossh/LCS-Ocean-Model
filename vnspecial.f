         
         subroutine vnspecial(splons,splone,splats,splate)
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

        integer vchoice(imax,jmax), vchsi(4,16000), vchsj(4,16000)
        integer vchs(4)
        real capf
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
        vchs(ll)=0
        enddo
c------------------------------------------
c       put the special boundary condition

c-------------------------------------------------------
           do j=splats,splate
           do i=splons,splone
  
c          do  j=178,245
c           do   i=235,260
          if (vchoice(i,j) .eq. 1) then
           vchs(1) =vchs(1)+1
           vchsi(1,vchs(1))=i
           vchsj(1,vchs(1))=j

c           write(131,*)'first',i,j,uchoice(i,j),uchi(1,uch(1))
          endif
          if (vchoice(i,j) .eq. 2) then
           vchs(2) =vchs(2)+1
           vchsi(2,vchs(2))=i
           vchsj(2,vchs(2))=j
          endif
          if (vchoice(i,j) .eq. 3) then
           vchs(3) =vchs(3)+1
           vchsi(3,vchs(3))=i
           vchsj(3,vchs(3))=j
          endif
          if (vchoice(i,j) .eq. 4) then
           vchs(4) =vchs(4)+1
           vchsi(4,vchs(4))=i
           vchsj(4,vchs(4))=j
          endif
          enddo
          enddo
 
c-----------------------
c         southern box bound
c         if (sclos)  then
c         else
c          do i=1,imax
c           vchs(2) =vchs(2)+1
c           vchsi(2,vchs(2))=i
c           vchsj(2,vchs(2))=1
c           enddo
c           endif

c       northern box boun
c         if (nclos)  then
c         else
c           do i=1, imax
c           vchs(1) =vchs(1)+1
c           vchsi(1,vchs(1))=i
c           vchsj(1,vchs(1))=jmax
c           enddo
c          endif
c       eastern box boun
c         if (eclos)  then
c         else
c         do j=1,jmax
c           vchs(4) =vchs(4)+1
c           vchsi(4,vchs(4))=imax
c           vchsj(4,vchs(4))=j
c           enddo
c          endif

c       western box boun
c         if (wclos)  then
c         else
c         do j=1,jmax
c           vchs(3) =vchs(3)+1
c           vchsi(3,vchs(3))=1
c           vchsj(3,vchs(3))=j
c           enddo
c           endif

        do ll=1,4
        if (vchs(ll).gt. 16000 ) then
        write(*,*) 'increase dimension of vch '
        stop
        endif
        enddo

c=========
       do ii=1,vchs(1)
        i=vchsi(1,ii)
        j=vchsj(1,ii)
        capf=.25*(txsav(i,j)+txsav(i,j-1)
     #   +txsav(i+1,j)+txsav(i+1,j-1))
     .   /xpsi2(nmm)
         vn(i,j,kp)=-capf/fv(i,j)
c        write(137,*)'second',i,j,capf,fv(i,j),vn(i,j,kp)
         enddo

         do  ii=1,vchs(2)
        i=vchsi(2,ii)
        j=vchsj(2,ii)
        capf=.25*(txsav(i,j)+txsav(i+1,j)
     #   +txsav(i,j-1)+txsav(i+1,j-1))
     .  / xpsi2(nmm)
         vn(i,j,kp)=-capf/fv(i,j)
         write(142,*)'second',i,j,capf,fv(i,j), vn(i,j,kp) 
         enddo

         do ii=1,vchs(3)
        i=vchsi(3,ii)
        j=vchsj(3,ii)
        vn(i,j,kp)=-vn(i,j+1,kp)
        enddo

        do ii=1,vchs(4)
        i=vchsi(4,ii)
        j=vchsj(4,ii)
        vn(i,j,kp)=-vn(i,j-1,kp)
        enddo



        return
        end

