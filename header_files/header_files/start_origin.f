! Please define start time origin

         subroutine origin(torigin,iday,month,iyear) 
         character*20 torigin, x_wind,y_wind 
         
         torigin='01-JAN-1991 00:00:00'  ! starting date of model
         iday=1                 ! start day of origin
         month=1                ! start month of origin
         iyear=2000             ! start year. Plase change time origin in the mr.f by changing torigin. This is for leap year calculation.

 
         return
          end
