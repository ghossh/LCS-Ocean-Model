! this is fortran program which will give i,j value for particular lat,lon
!this will help to apply special boundary condition and damp effect to any basin
!to run
!ifort -o exe   lonlat.f
!if you having different lon,lat value you replace 1,2
!and add each term`
              subroutine lonlat_damp(dplons,dplone,dplats,dplate)

       include "./header_files/dimension.h"
c        logical wclos,sclos,eclos,nclos
c         logical start

        integer divi
         integer latdd
       real delta,slat,elat,slon,elon,slonp,slatp
       real damplons,damplone,damplats,damplate  
       integer lon(imax),slonf,dlon3,lon3(1),lon4(1),dlon4
       integer lat(jmax),deltaf,slatf,dlat3,lat3(1),dlat4,lat4(1)
       integer dplons,dplone,dplats,dplate
         include "./header_files/param.h" 
c)c  define i for particular longitude
       do i=1,imax
       do j=1,jmax
       deltaf=int(100*delta)
       lon(i)=int((i-1)*deltaf)
       slonp=slon+0.5*delta
       slonf=int(100*slonp)
       dlon3=int(100*damplons)
         divi3=((dlon3-slonf)/deltaf)
         londd3= deltaf*divi3
         dlon4=int(100*damplone)
         divi4=((dlon4-slonf)/deltaf)
         londd4= deltaf*divi4

        lat(j)=int((j-1)*deltaf)
        slatp=slat+0.5*delta 
       slatf=int(100*slatp)
       dlat3=int(100*damplats)
         divj3=((dlat3-slatf)/deltaf)
         latdd3= deltaf*divj3
         dlat4=int(100*damplate)
         divj4=((dlat4-slatf)/deltaf)
         latdd4= deltaf*divj4

         if (lon(i).eq.londd3 .and.lat(j).eq.latdd3) then
        lon3(1)=i
        lat3(1)=j
         dplons=lon3(1)
         dplats=lat3(1)
         elseif  (lon(i).eq.londd4.and.lat(j).eq.latdd4) then
        lon4(1)=i
        lat4(1)=j
        dplone=lon4(1)
        dplate=lat4(1)
         endif 
         enddo
         enddo

        return
c       stop
       end
