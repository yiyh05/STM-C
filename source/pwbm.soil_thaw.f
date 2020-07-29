
      subroutine check_soilthaw(stemperature,phchange,
     $    soilnode,zsoilthaw,nlayer)
c     *********************************************************
c     ****************** added by yyh *************************
c     *********************************************************
c     subroutine calculates the soil thaw depth for each day 
c     this code may not apply to the permafrost areas...
c     note: the definition of permafrost is: having soil remaining below 0degC for two or more years...
c     so in this case, I was not able to tell whether this is seasonal thaw above permafrost or not.  
c     inputs:   stemperautre  -  soil temperature at each node depth
c               soilnode      -  soil node depth 
c     outputs:  zsoilthaw      -  soil thaw depth

c     note that I only check the first 15 layers 

      implicit none
      INCLUDE  'pwbm.soil_temperature.h'

      integer phchange(sx:lx+1)
      real*8  stemperature(sx:lx+1)
      real*8  soilnode(sx:lx)     
      real*8  zsoilthaw,tsoil,tsoil1,tsoil2
      integer i,isoilnode,nlayer,isoilALD

c whether there is phase change
      if(maxval(phchange) .eq. 1) then
        i=2
        do while(i .le.nlayer+1)
           if(phchange(i) .eq. 1) then
             i=i+1
           else
             isoilALD = i-1
             i=i+lx
           endif
        enddo
      else
        isoilALD = -1
        zsoilthaw = -999.0d0
        return
      endif
 
      if(isoilALD .gt. nlayer) then
        zsoilthaw = -999.0d0
        return
      endif
  
c first look for the soil node in the upper layer 
      tsoil1 = stemperature(isoilALD) ! the upper layer
      tsoil = stemperature(isoilALD+1)
      tsoil2 = stemperature(isoilALD+2) ! the lower layer
      if(tsoil1*tsoil .le. 0.0d0) then
        zsoilthaw = soilnode(isoilALD-1)+(soilnode(isoilALD)-
     $    soilnode(isoilALD-1))*abs(tsoil1)/(abs(tsoil)+abs(tsoil1))

      else if(tsoil*tsoil2 .le. 0.0d0) then ! then look for the soil node in the lower layer
        zsoilthaw = soilnode(isoilALD)+(soilnode(isoilALD+1)
     $    -soilnode(isoilALD))*abs(tsoil)/(abs(tsoil2)+abs(tsoil))
      else
        zsoilthaw = -999.0d0
      endif           

      return
      end
