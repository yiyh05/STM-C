! I changed the land cover type, to be consistent with the NLCD land cover type
!
      integer nvegclass      ! number of vegetation class for carbon simulation 
      integer ntotclass      ! number of total land cover types used by pwbm      

! Oct 29, 2014: finally I found I have to add a woody carbon pool!
!     for carbon simulation, I include Ice, and polar desert also (set LUE very low); otherwise, there is no SOC for those areas.
!     NLCD land cover type: open water(1),Ice/snow(2),Develped(3),Barren Land(4),
!                           Forest(5),Shrub(6),Grassland(7),Cultivated(8),wetland(9) 
      PARAMETER(nvegclass       = 9,  ! Forest(5),Shrub(6),Grassland(7),Cultivated(8),Wetlands(9)
     &          ntotclass       = 10) ! water(0),open water (1),perennial Ice/Snow(2),Developed(3),Barren Land(4),...


 
