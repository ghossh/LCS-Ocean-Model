c this is the main file where you will change all the parameters  according to your experiment.
c wind_dif_time is the time difference between input wind data 

        start_program=0        ! start_program=1 for restart only. Otherwise 0
	DELTA=0.1             ! resolution 
         
        slat=-29.8           !  start latitude
        elat=29.7            ! end latitude
         
        slon=30.2            !  start longitude
        elon=119.7           ! end longitude

	ndrun=1826               ! Number of days to be run. 
c 1826
        start_mode=1           ! Start mode
        end_mode=10             ! End mode

        dt=86400./96.          ! Model time step

	save_it= 1.0           ! (unit=days) save_it= 0.0 will save it in each time step 

      	visch=5.e6             ! Horizontal mixing
        a=0.00013              ! mode mixing coefficients


        wind_dif_time = 86400.00 ! Time frequency of wind. This is fixed for daily forcing. 
        midlat=0.0             ! Midlatitude
	ismth=41               ! Number of time steps between time smoothing
        f_plane=0  ! Choice of Coriolis force. if fplane=true  otherwise it will consider b-plane (put f_plane = 0)

c set parameter for damping
c 1 means condition is true, 0 means condition is false
        damp_u_south=1 ! if damp =true then damp_u_south=1
        damp_equ=0 ! if damp =true then damp_equ=1
        damplons=85.0
        damplone=98.0
        damplats=-5.0
        damplate=5.0

c Set parameter for spacial boundary condition
        sbc_apply=0 ! if sbc =true then sbc_apply=1
        spclons=82
        spclone=87
        spclats=6.5
        spclate=20
 

      
                
