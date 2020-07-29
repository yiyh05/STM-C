
       subroutine conductivity_v3(isoillayers,peat_fraction,
     $  soilporosity,porositycombined,thermcondmineral,thermcondpeat,
     $  dry_soil_thcnd,soilsolid_thcnd,
     $  heatcapmineral,heatcappeat,dry_soil_heatc)

      implicit none
      integer ilayer,isoillayers,ibaselayer
      real*8  thermcondmineral,thermcondpeat,thermcondair
      real*8  heatcapmineral,heatcappeat
      real*8  soilporosity,mineraldensity,dry_mineral_thcnd
      real*8  porositycombined(isoillayers),peat_fraction(isoillayers)
      real*8  dry_soil_thcnd(isoillayers),soilsolid_thcnd(isoillayers)
      real*8  heatcaplayer(isoillayers),dry_soil_heatc(isoillayers) 

c  bulk density and thermal conductivity of the mineral soil
      mineraldensity    = 2700.0d0*(1.0d0 - soilporosity)
      dry_mineral_thcnd = (0.135d0*mineraldensity + 64.70d0)/
     $   (2700.0d0 - 0.947d0*mineraldensity)

      do ilayer = 1, isoillayers

c  thermal conductivity of the dry soil (including peat) 
         dry_soil_thcnd(ilayer) = 0.05d0 * peat_fraction(ilayer) +
     $     dry_mineral_thcnd * (1.0d0 - peat_fraction(ilayer))  

c  thermal conductivity of the mineral soil 
         soilsolid_thcnd(ilayer) = (1.0d0 - peat_fraction(ilayer))*
     $     thermcondmineral + thermcondpeat * peat_fraction(ilayer)

         heatcaplayer(ilayer) = 
     $     heatcappeat * peat_fraction(ilayer) +
     $     heatcapmineral * (1.0d0 - peat_fraction(ilayer))    
          
         dry_soil_heatc(ilayer) = 
     $     heatcaplayer(ilayer) * (1-porositycombined(ilayer)) 

c  I do not think I should treat the bedrock separately --- I will change the "porosity" here.          
c        write(111,*) ilayer,peat_fraction(ilayer),
c     $     dry_soil_thcnd(ilayer),soilsolid_thcnd(ilayer),
c     $     dry_soil_heatc(ilayer)
      enddo    ! end of loop over soil layers   


      end


