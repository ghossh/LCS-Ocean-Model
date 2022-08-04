! this is fortran program which will give i,j value for particular lat,lon
!this will help to apply special boundary condition and damp effect to any basin
!to run
!ifort -o exe   lonlat.f
!if you having different lon,lat value you replace 1,2
!and add each term`
              subroutine lonlat(splons,splone,splats,splate)

       include "./header_files/dimension.h"
c        logical wclos,sclos,eclos,nclos
c         logical start

        integer divi
         integer latdd
       real delta,slat,elat,slon,elon,slonp,slatp
       real spclons,spclone,spclats,spclate  
       integer lon(imax),slonf,dlon1,lon1(1),lon2(1),dlon2
       integer lat(jmax),deltaf,slatf,dlat1,lat1(1),dlat2,lat2(1)
       integer splons,splone,splats,splate
         include "./header_files/param.h" 
c)c  define i for particular longitude
       do i=1,imax
       do j=1,jmax
       deltaf=int(100*delta)
       lon(i)=int((i-1)*deltaf)
       slonp=slon+0.5*delta
       slonf=int(100*slonp)
       dlon1=int(100*spclons)
         divi1=((dlon1-slonf)/deltaf)
         londd1= deltaf*divi1
         dlon2=int(100*spclone)
         divi2=((dlon2-slonf)/deltaf)
         londd2= deltaf*divi2

        lat(j)=int((j-1)*deltaf)
        slatp=slat+0.5*delta 
       slatf=int(100*slatp)
       dlat1=int(100*spclats)
         divj1=((dlat1-slatf)/deltaf)
         latdd1= deltaf*divj1
         dlat2=int(100*spclate)
         divj2=((dlat2-slatf)/deltaf)
         latdd2= deltaf*divj2

         if (lon(i).eq.londd1 .and.lat(j).eq.latdd1) then
        lon1(1)=i
        lat1(1)=j
         splons=lon1(1)
         splats=lat1(1)
         elseif  (lon(i).eq.londd2.and.lat(j).eq.latdd2) then
        lon2(1)=i
        lat2(1)=j
        splone=lon2(1)
        splate=lat2(1)
         endif 
         enddo
         enddo
c        write(*,*) splons,splats,'3333333333333'
c        write(*,*) splone,splate,'qqqqqqqqqq'

        return
c       stop
       end
