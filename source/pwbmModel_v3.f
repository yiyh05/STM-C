
      program pwbm
      implicit none

      integer narray,ncells,inum,ndim
      parameter (narray=1206129,ndim=46)
      integer nrowvals(narray),ncolvals(narray)
      integer length1,name_length,length,lengthout
      integer i,j,k,icell,iyr,init,iter,nspinup
      integer id,iindx,jindx,cellid
      integer ifrstyr,ilstyr,ivegcov,iveg,iveg1
      integer nyear,idcomp,ivar 
      integer icell_fst, icell_lst,numgrid
      integer isoilclass,site_ID,iwrite
      integer iannout(20),idailyout(20)
      integer monthlyout(20),inewveg(narray)
      integer sandpercent,claypercent,dzpr_ilayer
      real*8 r_lat,r_lon
      real*8 tempvals,snowdvals,smsurfvals,smrtzvals,rosnowvals
      real*8 tminvals,dswrfvals,ndvivals
      real*8 smsurfvals1,smsurfvals2
      real*8 smrtzvals1,smrtzvals2
      real*8 tmean,tair,swc,swc1,swc2,swc_scalar
      real*8 rootzswc,dswrf,tmin,ndvi
      real*8 fieldcapacity,soil_depth
      real*8 soilporosity,wiltingpoint
      real*8 bulkdensity,thermconddry
      real*8 heatcapdry,heat_capdry
      real*8 msoil_porosity,sand_frac
      real*8 therm_conddry,therm_condpeat
      real*8 snow_thcnd,snow_heatc
      real*8 thcndpeat_dry_frozen,thcndpeat_sat_frozen
      real*8 thcndpeat_dry_thaw,thcndpeat_sat_thaw
      real*8 peatdensity,hcapacitypeatdry
      real*8 dzpr,varjunk,rmissing
      parameter(rmissing=-999.0d0)
      dimension tempvals(ndim,narray),snowdvals(ndim,narray)
      dimension smsurfvals(ndim,narray),smrtzvals(ndim,narray)
      dimension rosnowvals(ndim,narray),tminvals(ndim,narray)
      dimension dswrfvals(ndim,narray),ndvivals(ndim,narray)
      dimension therm_conddry(narray),heat_capdry(narray)
      dimension sand_frac(narray),msoil_porosity(narray) 
      character header*60,dataname*150,fname*100
      character domain_file*100,mineralsoilfile*100,meantempfile*100
      character*100 soiltextfile,peatfracfile,socpoolfile
      character*100 sitefile,cellidfile
      character*100 temp_path,snow_path,outpath
      character*100 rosnow_name,swc_path
      character*100 swc_surf_name,swc_rtz_name
      character*100 tmin_name,dswrf_name,ndvi_name
      character*100 tmin_path,dswrf_path,ndvi_path
      character*100 site_name,temp_name,snowd_name
      character*20  yearstr
      character  homedir*100,str1*120
      logical ONECELL_SWITCH
      data ONECELL_SWITCH /.false./

C/* FTDJN 
      INCLUDE 'pwbm.soil_temperature.h'
      integer isoil_freeze,idptz_rock,nlayer
      parameter (nlayer=15)
c      real*8  soilt_out(narray,46,23)         !soil temperature outputs 
      real*8  soilt_out(210000,46,15)
      real*8  soil_temperature(narray,sx:lx+1)!soil temperature
      real*8  soil_liqfraction(narray,sx:lx+1)!Liquid water fraction in soil
      real*8  soil_liqfraction_prev(narray,sx:lx+1) !Liquid water fraction in soil at the previous time step
      integer soil_phchange(narray,sx:lx+1) !Layers where the phase change occurs
      real*8  soil_saturation0(narray,sx:lx+1) !Fraction of pore space filled with liquid water at the previous time step
      real*8  soil_saturation1(narray,sx:lx+1) !Fraction of pore space filled with liquid water at the current time step

      real*8  soil_x(narray,sx:lx)             ! dynamic snow/soil depth for each layer.  

      real*8  solid_soil_thcnd(narray,lx)     ! thermal conductivity of soil solids for each layer                      
      real*8  dry_soil_thcnd(narray,lx)       !Thermal conductivity of dry soil
      real*8  dry_soil_heatc(narray,lx)       !Heat capacity of dry soil
      real*8  soil_freezing_pd(lx)            !Global array
      real*8  soil_cleaness(lx)     

      real*8  soil_timestep(narray)           !Automatically updated timestep 
      integer soil_gc(narray)                 !Automatically updated counter of converged iterations

      integer soil_frozen(narray)             !Check if the soil column is entirely frozen
      integer frozen_length(narray,nlayer)    ! the frozen length for each soil layer                  
      real*8  soil_front                      !Location of the freezing front closest to the surface
      real*8  soil_ALD(narray),ALD_out(narray) !The active layer depth
      real*8  soilthaw(narray,ndim)             ! soil thaw depth at each time step
      integer isoil_ALD(narray)               !XREF's node of the active layer
      real*8  soil_time(narray)               !Time in soil. For testing

      !SNOW	        
c      real*8  soil_frsnow(narray)             !Fresh snow density, according to M.Sturm
      integer soil_sn(narray)                 !Number of the first snow node in the grid
      real*8  soil_snowheight
      real*8  soil_snowthcnd,soil_snowheatc    !Snow thermal properties
      real*8  sfreezing_pd,scleaness           !Just for initialization
     
      real*8  dzmm(lx)

      !INFILTRATION into the frozen soil
c      real*8 SS0, SS1, max_dw, rel_dw
 
C*/ FTDJN

c  NEW VARIABLES FOR VERSION 3. Can set these in an include file (?)
      real*8 soilnode(sx:lx)
      real*8 water4layer(sx:lx,narray),watertoday(sx:lx)
      real*8 ice4layer(sx:lx,narray),icetoday(sx:lx)
      real*8 thermcondlayer(lx),peatfraction(1:lx)
      real*8 peatfractionarray(lx,narray)
      real*8 porositycombined(lx),porosityarray(lx,narray)
      real*8 porosity_4layer(lx),capacity_4layer(sx:lx)
      real*8 ro_soft_snow

      data soilnode/-0.696D0,-0.402D0,-0.292D0,-0.164D0,-0.070D0,
     $    0.000D0,0.01D0,0.03D0,0.08D0,0.13D0,0.23D0,0.33D0,0.45D0,
     $  0.55D0,0.70D0,1.05D0,1.40D0,1.75D0,2.25D0,2.75D0,3.25D0,4.75D0, 
     $  7.25D0,14.0D0,22.0D0,31.0D0,40.0D0,50.0D0,60.00D0/

C/* FTDJN

C   variable or constants used for carbon model
      INCLUDE 'pwbm.carbon.h'

      integer  ipool,npool,ilayer
      parameter (npool=7)  ! only account for SOC dynamics in the top 15 layers (~3m)        
      real*8  gpp,npp,Rh,ann_npp
      real*8  totsoc,totsoc_prev
      real*8  flitter_above(ndim),flitter_below(ndim)
      real*8  AD(npool),D0_array(narray)  ! diffusion rates 
      real*8  GPP_array(narray,ndim),NEE_array(narray,ndim)
      real*8  NPP_array(narray,ndim),Rh_array(narray,ndim)
      real*8  kscalar4layer(ndim,nlayer),socc_pool(narray,nlayer,npool)
c      real*8  Rh_out(narray,ndim,nlayer),Rh4layer(ndim,nlayer)
      real*8  Rh_out(210000,ndim,nlayer),Rh4layer(ndim,nlayer) 
   
C/* FTDJN 
      do icell=1,narray
        call initialize_soil(soil_temperature(icell,:),
     $   soil_liqfraction(icell,:),soil_timestep(icell),soil_gc(icell),
     $   soil_sn(icell),soil_time(icell),soil_front,soilnode)

        soil_phchange(icell,:)=0
        soil_ALD(icell)=0.0d0
        isoil_ALD(icell)=-1
      enddo
C*/ FTDJN

C*/ YYH
      sfreezing_pd = -0.03d0
      scleaness = -0.625d0
      do i=1,lx
        soil_freezing_pd(i)=sfreezing_pd
        soil_cleaness(i) = scleaness
      enddo
C*/ YYH
 
C OPEN CONTROL FILE AND READ RECORDS
      call getenv('HOME', homedir)
      inum = index(homedir, ' ') - 1
      OPEN(11,FILE='pwbm-stm-c.initialize_v1.1.txt')
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      domain_file = str1(1:length1) ! lat, lon, vegetation type
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      meantempfile = str1(1:length1)
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      mineralsoilfile = str1(1:length1) ! read thermal properties and porosity of mineral soils
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1 
      soiltextfile = str1(1:length1)    ! soil texture file    
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      peatfracfile = str1(1:length1)  ! peat fraction for each soil layer
      READ (11,'(a100)') str1 

!     input dataset name 
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      temp_name = str1(1:length1)   ! 8-day MODIS LST dataset name
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      snowd_name = str1(1:length1)  ! 8-day snow depth dataset name
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      rosnow_name = str1(1:length1) ! 8-day snow density dataset name
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      swc_surf_name = str1(1:length1) ! 8-day surface sm dataset name
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      swc_rtz_name = str1(1:length1) ! 8-day rootzone sm dataset name
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      ndvi_name = str1(1:length1)   ! 8-day NDVI dataset name
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      tmin_name = str1(1:length1)   ! 8-day Tmin dataset name
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      dswrf_name = str1(1:length1)  ! 8-day downward solar radiation            
      READ (11, *)

!     input file path
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      temp_path = str1(1:length1) ! 8-day MODIS LST path
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      snow_path = str1(1:length1) ! 8-day snow depth path
      READ (11,'(a100)') str1
      length1 = index (str1, ' ') - 1
      swc_path = str1(1:length1) ! 8-day surface soil moisture path
      READ (11,'(a100)') str1
      length1 = index(str1, ' ') - 1
      ndvi_path = str1(1:length1)    ! 8-day ndvi data path
      READ (11,'(a100)') str1
      length1 = index(str1, ' ') - 1
      dswrf_path = str1(1:length1)   ! SWGDN data path
      READ (11,'(a100)') str1
      lengthout = index (str1, ' ') - 1
      outpath = str1(1:lengthout)    ! output file path 

      READ (11,*)                      ! read blank line in init file                 
      READ (11,*) thcndpeat_dry_thaw 
      READ (11,*) thcndpeat_sat_thaw
      READ (11,*) thcndpeat_dry_frozen
      READ (11,*) thcndpeat_sat_frozen
      READ (11,*) hcapacitypeatdry
      READ (11,*) idptz_rock
      READ (11,*) nspinup
      READ (11,*)                      ! read blank line in itit file

      READ (11,*) icell_fst             ! the startinng pixel and ending pixel
      READ (11,*) icell_lst
      READ (11,*) 
      print*, icell_fst, icell_lst

      READ (11,*) ifrstyr
      READ (11,*) ilstyr
      READ (11,*) idailyout(1)         ! this is actually the 8-day composite period 
      READ (11,*) iannout(1)    
      CLOSE(11)  

!    reading the domain file, define the domain/sites/basin to run the model
      open(unit=14,file=domain_file)
      open(unit=17,file=meantempfile)
c  Read drainage cell information from table.  Indicies for
c  latitude and longitude are read in from file.  j is for lat, i for
c  lon. They are not used in model.  Fill arrays for basin info.
c  initialize arrays with vegetation
c  NLCD land cover type: ocean(0),open water(1),Ice/snow(2),developed(3),
c  barren_land(4),forest(5),shrub(6),grassland(7),cultivated(8),wetlands(9)                                          
      read(14,*) header                 ! skipping header/label
      read(17,*) header 
      do icell=1, narray
         read(14,16,end=99) id,iindx,jindx,r_lat,r_lon,iveg
         read(17,*) id, tmean
16       format (i10,2i6,2f10.2,i5)
         nrowvals(icell) = iindx  ! 1-based
         ncolvals(icell) = jindx 
         inewveg(icell) = iveg

         if(tmean .lt. 0.0d0) then
             D0_array(icell) = 5.0d0 ! cm2/yr, diffusion rate for permafrost areas
         else
             D0_array(icell) = 2.0d0
         endif                

         do i=sx, lx+1
            soil_temperature(icell,i) = tmean
         enddo
      enddo
      close(17)
      close(14) 
 99   ncells = icell - 1
      print*, ncells, 'Arctic Basin records read'

      open(62,file=peatfracfile)   
      do icell=1, ncells
         read(62,*) id,i,j,(peatfraction(i),i=1,10)
     
         do i=1, lx
           if(i .le. 10) then
               peatfractionarray(i,icell)=peatfraction(i)*0.01d0
           else
               peatfractionarray(i,icell)=0.0d0 
           endif
         enddo
      enddo
      close(62)

c  read soil texture and soil peat fraction file 
      open(61,file=mineralsoilfile)! porosity, thermal conductivity and heat capacity of mineral soils 
      open(63,file=soiltextfile)                    

      read(61,*) header
      read(63,*) header
      do icell = 1, ncells
         read(61,*) id,dzpr,dzpr_ilayer,soilporosity,
     $              thermconddry,heatcapdry
         read(63,*) id,sandpercent,claypercent

         if(sandpercent .lt. 0 .or. sandpercent .gt. 100) then
           print*, sandpercent
           stop
         endif

c recalculate the thermal conductivity of soil solids
         thermconddry = 7.7d0**(sandpercent*0.01d0)*2.5d0**
     $    (1.0d0-sandpercent*0.01d0)
         if(thermconddry .lt. 2.5d0 .or. thermconddry .gt. 7.7d0)
     $    then
            print*, thermconddry, 'out of range'
            stop
         endif
         
c  initialize soil water and ice profiles 
         do i=sx,lx+1
            soil_saturation0(icell,i) = 1.0d0
            soil_saturation1(icell,i) = 1.0d0
         enddo

c  Look up soil characteristics for given soil class.  Convert units for soil 
c  variables.  Also set some array values.
c         call lookup_v3(isoilclass,lx,peatfraction,soilporosity,
c     $     fieldcapacity,wiltingpoint,bulkdensity,thermconddry,
c     $     heatcapdry,sfreezing_pd,scleaness,porositycombined,
c     $     fcapacitycombined,hk_sat,psi_sat,beta_combined)

!  Set the soil porosity array by soil layer and grid cell.  
         do i=1,lx
            if(peatfractionarray(i,icell) .lt. 0.0d0 .or. 
     $        peatfractionarray(i,icell) .gt. 1.0d0) then 
               print*, i, peatfractionarray(i,icell)
               stop
            endif
            porositycombined(i)=(1.0d0-peatfractionarray(i,icell))
     $       *soilporosity + peatfractionarray(i,icell)*0.8d0
            if(porositycombined(i) .lt. 0.0d0 .or. 
     $        porositycombined(i) .gt. 0.8d0) then
               print*, i, porositycombined(i)
               stop
            endif
            porosityarray(i,icell) = porositycombined(i) ! unit: m3/m3
c            if(i .ge. idptz_rock) porosityarray(i,icell) = 0.05d0
c            print*, i, peatfractionarray(i,icell), porositycombined(i)
         enddo   

        msoil_porosity(icell) =  soilporosity 
        therm_conddry(icell) = thermconddry  ! save to array by grid cell
        heat_capdry(icell) = heatcapdry
        sand_frac(icell)   = sandpercent * 0.01d0

C/* FTDJN
        call compute_dzmm(soilnode*1.0d3, lx, dzmm) ! dzmm: unit: mm
 
        do i=sx,0
c   this is used to calculate the unfrozen water fraction below freezing point
           call soilsat(soil_temperature(icell,i),
     &       soil_liqfraction(icell,i),
     &        scleaness,sfreezing_pd)
           soil_saturation0(icell,i) = 0.0d0
           soil_saturation1(icell,i) = 0.0d0
           soil_liqfraction(icell,i) = 1.0d0
        enddo
c  initialize the "above surface" layer
        water4layer(sx:0,icell) = 0.0d0
        ice4layer(sx:0,icell) = 0.0d0

        do i=2,lx+1
          call soilsat(soil_temperature(icell,i),
     &       soil_liqfraction(icell,i),
     &       scleaness,sfreezing_pd)
          if(soil_liqfraction(icell,i) < CR_LIQ_FRACTION) then
             icetoday(i-1) = porosityarray(i-1,icell)*dzmm(i-1)*1.0d0 !0.9d0         
             watertoday(i-1) = 0.0d0
             soil_saturation0(icell,i) = icetoday(i-1)/
     &          (porosityarray(i-1,icell)*dzmm(i-1))
             soil_saturation1(icell,i)=soil_saturation0(icell,i)
          else
             icetoday(i-1) = 0.0d0
             watertoday(i-1) = porosityarray(i-1,icell)*dzmm(i-1)*1.0d0!0.9d0
             soil_saturation0(icell,i)=watertoday(i-1)/
     &        (porosityarray(i-1,icell)*dzmm(i-1))
             soil_saturation1(icell,i)=soil_saturation0(icell,i)
          endif
        enddo
C*/ FTDJN

        do i = 1, lx
           if(watertoday(i) .lt. 0.0d0) stop 'layer water < 0'
           if(icetoday(i) .lt. 0.0d0) stop 'layer ice < 0'
           water4layer(i,icell) = watertoday(i)    ! as input from spinup
           ice4layer(i,icell) = icetoday(i)        ! as input from spinup
            
c  Make sure water+ice of layer does not exceed capacity of layer.
           varjunk = dzmm(i) * porositycombined(i) ! junkvar is dummy capacity
           if(water4layer(i,icell)+ice4layer(i,icell).gt.varjunk) then
              if(soil_liqfraction(icell,i) < CR_LIQ_FRACTION) then
                 ice4layer(i,icell) = varjunk * 1.0d0 !0.9d0
                 water4layer(i,icell) = 0.0d0
              else
                 ice4layer(i,icell) = 0.0d0
                 water4layer(i,icell) = varjunk * 1.0d0 !0.9d0
              endif
           endif
  
c          print*, icell, i, water4layer(i,icell), ice4layer(i,icell)
        enddo
      enddo

c  ****************************************************************************
c  ****************** start spin-up loop **************************************
c  ****************************************************************************
      do init = 1, nspinup+1   ! a "spinup" loop around the time loop
                               ! usually won't use in this time series program      
        if(init .lt. nspinup+1) then
           nyear = 1
           iwrite = 0
        else
           iwrite = 1
           nyear = ilstyr - ifrstyr + 1

           print*, '************  Beginning transient run  ***********'
        endif

c       initialization
        if(init .eq. nspinup) then         

          do icell=1,ncells
          do i=1,nlayer
          do ipool=1,npool
             socc_pool(icell,i,ipool) = rmissing
          enddo
          enddo
          enddo
 
        endif

c  ****************************************************************************
c  ************************** start of year loop ******************************                              
c  ****************************************************************************
        do k = 1, nyear    

          if(init .ge. nspinup) then  
            do icell=1,ncells
              ALD_out(icell) = -999.0d0
              do idcomp=1,ndim
                 GPP_array(icell,idcomp) = rmissing
                 NPP_array(icell,idcomp) = rmissing
                 Rh_array(icell,idcomp)  = rmissing
                 NEE_array(icell,idcomp) = rmissing
                 soilthaw(icell,idcomp)  = rmissing

                 if(icell .le. 210000) then
                    do ilayer=1, 15
                       soilt_out(icell,idcomp,ilayer) = rmissing
                       Rh_out(icell,idcomp,ilayer) = rmissing
                    enddo
                 endif
              enddo
       
              do ilayer=1,nlayer
                 frozen_length(icell,ilayer) = 0
              enddo
            enddo             
          endif
          
          iyr = ifrstyr + k - 1 
          if(init.eq.1 .or. init.eq.nspinup+1) then
             print*, 'Reading surface meteorology file'
c  8-day composite period, including LST, snow depth, surface and rootzonec     
             call read_data_ann(temp_path,temp_name, iyr,
     $            nrowvals,ncolvals,tempvals)
             call read_data_ann(snow_path,snowd_name,iyr,
     $            nrowvals,ncolvals,snowdvals)
             call read_data_ann(snow_path,rosnow_name,iyr,
     $            nrowvals,ncolvals,rosnowvals)
             call read_data_ann(swc_path,swc_surf_name,1900,
     $            nrowvals,ncolvals,smsurfvals)
             call read_data_ann(swc_path,swc_rtz_name,1900,
     $            nrowvals,ncolvals,smrtzvals)
             call read_data_ann(temp_path,tmin_name,iyr,
     $            nrowvals,ncolvals,tminvals)
             call read_data_ann(dswrf_path,dswrf_name,iyr,
     $            nrowvals,ncolvals,dswrfvals)
             call read_data_ann(ndvi_path,ndvi_name,iyr,
     $            nrowvals,ncolvals,ndvivals)   
          endif

c  ****************************************************************************
c  ************************ begin gridcell loop ********** *********************
c  ****************************************************************************            
          do icell = icell_fst, icell_lst
c  ignore the areas covered by water and icesnow 
c  land cover type: water(0),open water(1),Ice/snow(2),Developed(3)
c  BarrenLand(4),Forest(5),Shrub(6),Grassland(7),Cultivated(8),wetlands(9)
             ivegcov = inewveg(icell)
             if(ivegcov .le. 2 .or. ivegcov .gt. 9) cycle
c  initialization 
             if(init .ge. nspinup) then
               do idcomp=1,ndim
                 do i=1,nlayer
                   kscalar4layer(idcomp,i) = rmissing
                   Rh4layer(idcomp,i)      = rmissing
                 enddo
               enddo
             endif
  
c  ****************************************************************************
c  ************************* begin day loop ******* ***************************
c  **************************************************************************** 
          do idcomp = 1, ndim  

             tair = tempvals(idcomp,icell)
             swc1 = smsurfvals(idcomp,icell)  ! surface soil wetness m3/m3
             swc2 = smrtzvals(idcomp,icell)
             soil_snowheight = snowdvals(idcomp,icell) ! snow depth, m
             ro_soft_snow = rosnowvals(idcomp,icell) ! snow density, kg/m3

c data anomalies 
             if(soil_snowheight .lt. 0.0d0) then
                stop 'snow depth out of range'
             endif
             if(ro_soft_snow .lt. 0.0d0) then
                stop 'snow density < 0.0'
             endif
             if(ro_soft_snow .gt. 600.0d0) then
                ro_soft_snow = 600.0d0
             endif
c fresh falling snow             
             if((soil_snowheight .gt. 0.0d0) .and. 
     $        (ro_soft_snow .eq. 0.0d0)) then
                ro_soft_snow = 50.0d0 ! fresh falling snow density
             endif 

             if(tair .gt. 150.0d0) tair = tair - 273.15d0
             if(swc1 .lt. 0.0d0 .or. swc1 .gt. 1.0d0) then 
                print *, 'swc1 out of range'
                stop
             endif
             if(swc2 .lt. 0.0d0 .or. swc2 .gt. 1.0d0) then
                print *, 'swc2 out of range'
                stop
             endif
c adjust the snow density -- this is for 8-day, should it be fine.
c if there is snow existence when temperature is high, the simulated
c soil temperature seems weird. 
             if(tair .gt. 2.0d0) then
                soil_snowheight = 0.0d0
                ro_soft_snow = 0.0d0
             endif

             do i=1, lx
                peatfraction(i) = peatfractionarray(i,icell)
                porosity_4layer(i) = porosityarray(i,icell) ! units:cm3/cm3
                capacity_4layer(i) = porosity_4layer(i)*dzmm(i)! units:cm3/cm3
             enddo      

c initialize the soil moisture profile -- need a more sophyisticated scheme
             do i=0,7
               if(i .le. 3) then 
                 swc = swc1 
               else 
                 swc = swc2
               endif
               if(soil_liqfraction(icell,i+1)<CR_LIQ_FRACTION) then
                 ice4layer(i,icell)=swc*porosityarray(i,icell)*dzmm(i)
                 water4layer(i,icell)=0.0d0
               else
                 water4layer(i,icell)=swc*porosityarray(i,icell)*dzmm(i)
                 ice4layer(i,icell)=0.0d0  
               endif
             enddo 
           
             do i=8,lx
               swc = swc2 
               if(soil_liqfraction(icell,i+1)<CR_LIQ_FRACTION) then
                 ice4layer(i,icell)=swc*porosityarray(i,icell)*dzmm(i)
                 water4layer(i,icell)=0.0d0
               else
                 water4layer(i,icell)=swc*porosityarray(i,icell)*dzmm(i)
                 ice4layer(i,icell)=0.0d0
               endif
             enddo          

c   this is only for dry soils and will be updated in the soil temperature routines.
c   the original thermal conductivity for dry soils seems too high, so I revised it.     
             therm_condpeat = 0.25d0  ! W/m/K   
             call conductivity_v3(lx,peatfraction,msoil_porosity(icell),
     $         porosity_4layer,therm_conddry(icell),therm_condpeat,
     $         dry_soil_thcnd(icell,:),solid_soil_thcnd(icell,:),
     $         heat_capdry(icell),hcapacitypeatdry,
     $         dry_soil_heatc(icell,:))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  Dima: At this point variable water4layer(0,icell) should have all surface
c         available water on days when the top soil layer ice4layer(1,icell)
c         has ice.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C/* FTDJN
c   snow parameters needed: snow density, thermal conductivity and heat capacity
c   both snow thermal conductivity and heat capacity were based on snow density 

!   ****  snow thermal conductivity (W m-1 K-1) ****
!            SNOW_THCND = 1.38D-1 - (1.01D0 * (ro_soft_snow/1.0D3)) +  ! eq 4
!     $                 (3.233D0 * ((ro_soft_snow/1.0D3) ** 2.0D0) )
             SNOW_THCND = 2.2D0 * (ro_soft_snow/1.0d3)**1.88
c     snow heat capacity 
c            SNOW_HEATC = 2.05114D3 * ro_soft_snow * 1.0e-6 ! units (MJ m3 K-1)
             SNOW_HEATC = CI * ro_soft_snow/(0.9167d3)  ! units (MJ m3 K-1) 
             soil_snowthcnd = SNOW_THCND
             soil_snowheatc = SNOW_HEATC 

             do i=2,lx+1
               if(soil_liqfraction(icell,i)<CR_LIQ_FRACTION) then
                 soil_saturation0(icell,i)=ice4layer(i-1,icell)/
     $             capacity_4layer(i-1) !Soil saturation at the previous time step
                 soil_saturation1(icell,i)=soil_saturation0(icell,i) 
               else
                 soil_saturation0(icell,i)=water4layer(i-1,icell)/
     &             capacity_4layer(i-1)  
                 soil_saturation1(icell,i)=soil_saturation0(icell,i)
               endif
               ! test the code
c               soil_saturation0(icell,i)=float(i-1)*0.1d0
c               if(soil_saturation0(icell,i) .gt. 1.0d0) then
c                  soil_saturation0(icell,i) = 1.0d0
c               endif
c               soil_saturation1(icell,i)=soil_saturation0(icell,i)
             enddo

             if(init .eq. nspinup) then 
               iwrite = 0
             else
               iwrite = 0
             endif
    
             call update_properties(soil_x(icell,:),soil_sn(icell), 
     &         soil_snowheight,soil_snowthcnd,soil_snowheatc,
     &         peatfraction,thcndpeat_dry_thaw,thcndpeat_sat_thaw,
     &         thcndpeat_dry_frozen,thcndpeat_sat_frozen,
     &         porosityarray(:,icell),dry_soil_thcnd(icell,:),
     &         solid_soil_thcnd(icell,:),dry_soil_heatc(icell,:),
     &         soil_liqfraction(icell,:),
     &         soil_saturation0(icell,:),soil_saturation1(icell,:),
c     &        soil_freezing_pd(icell,:),soil_cleaness(icell,:))
     &         soil_freezing_pd(:), soil_cleaness(:),iwrite)
             soil_time(icell)=soil_time(icell)+1.0D0 ![days]	!Used for testing purposes, shows how many days are computed

             SURFACE_TEMPERATURE_BEGIN=                      ![C]	!Initial snow/soil surface temperature at time="TIME"
     &        soil_temperature(icell,soil_sn(icell))
             SURFACE_TEMPERATURE_END  = tair                 ![C]	!Final snow/soil surface temperature   at time="TIME_END"

             soil_liqfraction_prev(icell,:)=soil_liqfraction(icell,:)
             call stemperature(soil_temperature(icell,:),
     &             soil_liqfraction(icell,:),
     &             soil_timestep(icell),
     &             soil_sn(icell),soil_gc(icell))

c check phase change 
             if(idcomp .eq. ndim) then
               if(maxval(soil_phchange(icell,:)) .eq. 1) then
                  i=2 
                  do while(I.LE.LX)
                    if(soil_phchange(icell,i).eq.1) then
                       i=i+1
                    else
                       isoil_ALD(icell)=i-1
                       i=i+lx
                    endif
                  enddo
                  
                  i=isoil_ALD(icell)
                  if(soil_frozen(icell).eq.1) then ! if the soil at at least one depth is completely frozen
                     soil_ALD(icell)=(XREF(i-1)+XREF(i))*0.5d0  ! Should calculate ALT differently                       
                  else  ! otherwise, it is seasonal freezing
c                     soil_ALD(icell)=-(XREF(i-1)+XREF(i))*0.5d0
                     soil_ALD(icell)=-999.0d0
                  endif
               else
                  isoil_ALD(icell)=-1
                  soil_ALD(icell)=-999.0d0
               endif
   
               soil_frozen(icell) = 0
               soil_phchange(icell,:)=0
             endif

             do i=2,lx+1
               if((soil_liqfraction_prev(icell,i)-CR_LIQ_FRACTION)*
     &          (soil_liqfraction(icell,i)-CR_LIQ_FRACTION)<0.0d0) then
                 soil_phchange(icell,i)=1
               endif

               if(soil_liqfraction(icell,i) < CR_LIQ_FRACTION) then
                  frozen_length(icell,i-1) = frozen_length(icell,i-1)+1
               endif
             enddo

c check whether the entire soil is frozen
             call check_soilfrozen(soil_liqfraction(icell,:),
     $         isoil_freeze)
             soil_frozen(icell)=max(soil_frozen(icell),isoil_freeze)

             if(init.eq.nspinup+1) then
             ! output
               do i=1, 15 
                 soilt_out(icell-icell_fst+1,idcomp,i)=
     $              soil_temperature(icell,i+1)
               enddo               

               if(soil_frozen(icell) .eq. 1) then
                 call check_soilthaw(soil_temperature(icell,:),
     $            soil_phchange(icell,:),soilnode,soilthaw(icell,idcomp)
     $            ,15)
               endif
             endif

             do i=2,lx+1
               if((soil_liqfraction(icell,i) < CR_LIQ_FRACTION).and. 
     &          (water4layer(i-1,icell) > 0.0d0) ) then
                  ice4layer(i-1,icell) = water4layer(i-1,icell)
                  water4layer(i-1,icell)=0.0d0  
               endif
               if((soil_liqfraction(icell,i) > CR_LIQ_FRACTION).and. 
     &          (ice4layer(i-1,icell) > 0.0d0) ) then
                  water4layer(i-1,icell) = ice4layer(i-1,icell)
                  ice4layer(i-1,icell)=0.0d0
               endif
             enddo
C*/ FTDJN

c   first calculate soil decomposition scalar for different depths
c   check the two variables: "water4layer" & "ice4layer"; they did not seem right to me
c   July 3rd, 2017: unfrozen water content is important; make sure it is right here.
             call kdecomp_scalar(dzmm,porosity_4layer,
     $        water4layer(:,icell),ice4layer(:,icell),
     $        soil_liqfraction(icell,:),soil_temperature(icell,:),
     $        kscalar4layer(idcomp,:),nlayer)  ! only calculate for soil depth 0-3m
           enddo   ! end of day loop
c *****************************************************************************
c ********************end of day loop ***************************************** 
c *****************************************************************************     

           ALD_out(icell) = maxval(soilthaw(icell,:)) 
c ****************************************************************************
c ***************************** run the SOC spin-up here *********************
c ****************************************************************************
           if(init .ge. nspinup) then

c             print*, icell, 'running cell '
c run the gpp module at the end of STM spin-up period
             call gpp_model(ivegcov,ndvivals(:,icell),
     $        tminvals(:,icell),dswrfvals(:,icell),smrtzvals(:,icell),
     $        GPP_array(icell,:),NPP_array(icell,:))
 
             ann_npp = 0.0d0
             do idcomp = 1, ndim
               if(NPP_array(icell,idcomp) .ge. 0.0d0) then
                 ann_npp = ann_npp + NPP_array(icell,idcomp)
               endif
             enddo ! loop for each temporal period
             if(ann_npp .lt. 0.0d0) stop 'annual npp < 0'

c A dynamic litterfall allocation scheme to account for litterfall seasonality
             call litter_allocate(ivegcov,ndvivals(:,icell),
     $        flitter_above,flitter_below)

c start SOC spinup process
             if(init .eq. nspinup) then
c first speed up the SOC spinup process through using accelerating rates
               do ipool=1, npool
                 AD(ipool) = 1.0d0
               enddo
               AD(5) = 5.0d0
               AD(6) = 100.0d0
               AD(7) = 1.0d0
c initialize the SOC pools
               do i=1,nlayer
               do ipool=1,npool
                  socc_pool(icell,i,ipool)=0.0d0
               enddo
               enddo
c use the total SOC content as a control variable
               call calSOCContent(socc_pool(icell,:,:),totsoc)

               iter = 1                  
               totsoc_prev = totsoc + 5.0d0
               do while((abs(totsoc-totsoc_prev) > 0.5d0) .and. 
     $             (iter < 3000))         
                 totsoc_prev = totsoc
                 call soc_decomp(ivegcov,AD,D0_array(icell),
     $             sand_frac(icell),ann_npp,
     $             flitter_above,flitter_below,kscalar4layer,
     $             socc_pool(icell,:,:),Rh4layer)
           
                 call calSOCContent(socc_pool(icell,:,:),totsoc)

                 iter = iter + 1
               end do ! end of while loop 

c scale the SOC pools
               do i=1,nlayer
               do ipool=1, npool
                 socc_pool(icell,i,ipool)=socc_pool(icell,i,ipool)
     $             *AD(ipool)
               enddo
               enddo
                
               call calSOCContent(socc_pool(icell,:,:),totsoc)

c then do the normal spin-up process 
               do ipool=1, npool
                 AD(ipool) = 1.0d0
               enddo

               totsoc_prev = totsoc + 1.0d0
               do while((abs(totsoc-totsoc_prev) > 0.5d0) .and. 
     $             (iter < 5000))
                 totsoc_prev = totsoc

                 call soc_decomp(ivegcov,AD,D0_array(icell),
     $             sand_frac(icell),ann_npp,
     $             flitter_above,flitter_below,kscalar4layer,
     $             socc_pool(icell,:,:),Rh4layer)

                 call calSOCContent(socc_pool(icell,:,:),totsoc)

                 iter = iter + 1
               enddo

             endif ! finish the SOC spinup
                
c followed by a transit run
             call soc_decomp(ivegcov,AD,D0_array(icell),
     $         sand_frac(icell),ann_npp,
     $         flitter_above,flitter_below,kscalar4layer,
     $         socc_pool(icell,:,:),Rh4layer)

c calculate the carbon flux 
             do id=1, ndim
                do i=1, nlayer
                  Rh_out(icell-icell_fst+1,id,i) = Rh4layer(id,i)
                enddo

                Rh = 0.0d0
                do i=1, nlayer
                  if(i .eq. 1) then
                    Rh = Rh + soilnode(i)*Rh4layer(id,i)
                  else
                    Rh = Rh+(soilnode(i)-soilnode(i-1))
     $               *(Rh4layer(id,i)+Rh4layer(id,i-1))/2.0d0
                  endif
                enddo

                if(Rh < 0.0d0) stop 'Rh < 0.'
                Rh_array(icell,id)  = Rh
                NEE_array(icell,id) = Rh - NPP_array(icell,id)
             enddo       

         endif  ! end fo SOC module

c *************************************************************************
c ******************** end of the SOC spin-up *****************************
c *************************************************************************
       enddo ! end of gridcell loop 
c *************************************************************************
c ******************** end of gridcell loop ******************************* 
c ************************************************************************* 
       if(init .eq. nspinup+1) then
          write(yearstr,'(i4)') iyr
          fname=outpath(1:lengthout)//'ALT.'//trim(yearstr)//'.txt'
          call write_ann_subset(fname,icell_fst,icell_lst,ALD_out) 

          numgrid = icell_lst - icell_fst + 1
          fname=outpath(1:lengthout)//'GPP.'//trim(yearstr)//'.bin'
          call write_ann_subset_v1(fname,icell_fst,icell_lst,numgrid,
     $        ndim,GPP_array)

          fname=outpath(1:lengthout)//'NEE.'//trim(yearstr)//'.bin'
          call write_ann_subset_v1(fname,icell_fst,icell_lst,numgrid,
     $        ndim,NEE_array)

          fname=outpath(1:lengthout)//'Rh.'//trim(yearstr)//'.bin'
          call write_ann_subset_v1(fname,icell_fst,icell_lst,numgrid,
     $        ndim,Rh_array)

          fname=outpath(1:lengthout)//'soilt.'//trim(yearstr)//
     $        '.0-1m.bin'
          call write_soil_subset(fname,icell_fst,icell_lst,numgrid,
     $        10,ndim,soilt_out(1:numgrid,:,1:10))

          fname=outpath(1:lengthout)//'Rh4layer.'//trim(yearstr)//
     $        '.0-1m.bin'
          call write_soil_subset(fname,icell_fst,icell_lst,numgrid,
     $        10,ndim,Rh_out(1:numgrid,:,1:10))

          fname=outpath(1:lengthout)//'SOC4layer.'//trim(yearstr)//
     $        '.0-1m.bin'
          call write_soil_subset(fname,icell_fst,icell_lst,numgrid,
     $        npool,10,socc_pool(icell_fst:icell_lst,1:10,:))
        endif

        write(*,110) 'Year ', iyr,  ' processed'
 110    format(a,i4,a)

      enddo ! end of year loop   
c *************************************************************************
     
      if (init .ge. 2) print*, init, ' times through spinup loop'
      enddo     
c *******************************************************************
c ******************** end of spinup loop ***************************
c *******************************************************************

c Sum the variable storing each basins' total (over time) discharge.
c Then compare this value to the total export for Arctic drainage
c (sum2ocean_alltime) and the non_Arctic drainages (sum_intrnl)
c yyh: the variables below actually were not calculated...

 600  print*, '******************************************************'
      print*, '********** Simulation sucessfully completed **********'
      print*, '******************************************************'

      stop
      end

