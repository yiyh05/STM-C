      subroutine litter_allocate(iveg,ndviarray,frac_above,frac_below)
c     I revised this to work for 8-day composite period as well 
c     subroutine distributes annual litterfall into daily time step using 
c     a simple dynamic allocation scheme based on land cover type and satellite 
c     NDVI times series, follow Randerson et al. (1996), but simplified
c     I assume 1) a constant litterfall throughout the year, representing constant turnover
c     of wood and fine roots, with the proportion dependent on land cover type (woody or herbaceous)
c     2) the rest of the litterfall falls during the later growing-season (senescence of leaves and fine roots),
c     and 3) throughout the growing season (turnover of fine roots during active growing season)
c     somehow, I think the seasonality of environmental factors has much big impact on the seasonality of carbon fluxes;
c     rather than the seasonality of litterfall.... 
c     however, all those ideas should be tested and refined later.. 
c
c     Inputs:  iveg        --- land cover type (mostly distinguish woody vs herbaceous plants)
c              ndviarray    --- satellite NDVI time series 
c     Outputs: frac_above  --- daily aboveground litterfall (unitless, percentage)
c              frac_below  --- daily belowground litterfall (unitless, percentage)                           
      
      implicit none
      INCLUDE 'pwbm.carbon.h'
      
      integer iveg,num,idcomp1,id,id1,ndim
      parameter(ndim=46)
      real*8  f1,f2,f3,f4,fleaf,froot,fmax,fsum,sum_var
      real*8  ndviarray(ndim),frac_above(ndim),frac_below(ndim)  ! data type should be consistent with the main code.
      real*8  fconst(nvegclass)
      real*8  leaf_ratio(nvegclass),root_ratio(nvegclass)
      real*8  flit_con(ndim),flit_var1(ndim),flit_var2(ndim)

c NLCD land cover type: open water(1),ice/snow(2),developed(3),barren_land(4),Forest(5),
c Shrub(6), Grassland(7), Cultivated(8), Wetlands(9)                                    
c it should also account for a constant proportion of leaf and fine root turn-over throughout the year.. 
      data fconst/0.0D0,0.0D0,0.0D0,0.0D0,0.5D0,0.3D0,0.15D0,0.15D0,
     $  0.15D0/
c leaf:root ratio represents the ratio of aboveground to belowground litterfall
      data leaf_ratio/0.0D0,0.0D0,0.5D0,0.5D0,0.4D0,0.4D0,0.4D0,0.5D0,
     $  0.4D0/ ! this is actually consistent with Yi et al. 2015
      data root_ratio/0.0D0,0.0D0,0.5D0,0.5D0,0.6D0,0.6D0,0.6D0,0.5D0,
     $  0.6D0/
                                  
c initialization 
      do idcomp1=1,ndim
         frac_above(idcomp1) = 0.0d0
         frac_below(idcomp1) = 0.0d0
      enddo
    
      fleaf = leaf_ratio(iveg)
      froot = 1.0d0 - fleaf

c (1) the constant portion of the litterfall throughout the year (woody components)
      do idcomp1=1,ndim
         flit_con(idcomp1) = fconst(iveg)/float(ndim)
c  initilization
         flit_var1(idcomp1) = 0.0d0 !this is the variable litterfall fraction due to leaf senescence
         flit_var2(idcomp1) = 0.0d0 !this represents the variable litterfall fraction occuring during active growing season (e.g. fine roots)
      enddo
c (2) find the variable leaf/fine root litterfall fraction due to leaf senescence, 
c     which is defined by the NDVI time series 
      sum_var = 0.0d0
      do idcomp1=1,ndim
c   here I use four-month window to get the changes in the leaves (LAI/NDVI)
         f1=0.0d0
         num=0
         do id=-3,-1
            id1 = idcomp1+id
            if(id1 .lt. 1) cycle
            if(ndviarray(id1) .lt. -1.0d0) cycle
            f1 = f1 + ndviarray(id1)
            num = num + 1
         enddo         
         if(num .eq. 3) then 
            f1=f1/float(num)
         else 
            f1=0.0d0
         endif

         f2 = 0.0d0
         num = 0
         do id=-6,-4
            id1 = idcomp1 + id 
            if(id1 .lt. 1) cycle
            if(ndviarray(id1) .lt. -1.0d0) cycle
            f2 = f2 + ndviarray(id1)
            num = num + 1     
         enddo
         if(num .eq. 3) then
            f2=f2/float(num)
         else
            f2=0.0d0
         endif
 
         f3 = 0.0d0
         num = 0
         do id=1, 3
            id1 = idcomp1+id
            if(id1 .gt. ndim) cycle
            if(ndviarray(id1) .lt. -1.0d0) cycle
            f3 = f3 + ndviarray(id1)
            num = num + 1
         enddo
         if(num .eq. 3) then
            f3=f3/float(num)
         else
            f3=0.0d0
         endif

         f4 = 0.0d0
         num = 0
         do id=4, 6
            id1 = idcomp1+id
            if(id1 .gt. ndim) cycle
            if(ndviarray(id1) .lt. -1.0d0) cycle
            f4 = f4 + ndviarray(id1)
            num = num + 1
         enddo
         if(num .eq. 3) then
            f4=f4/float(num)
         else
            f4=0.0d0
         endif
c          print *,idoy,f2,f1,f3,f4
         flit_var1(idcomp1) = (0.5*f2+f1-f3-f4*0.5) 
         if(flit_var1(idcomp1) .gt. 0.0d0) then
            sum_var = sum_var + flit_var1(idcomp1)
         else
            flit_var1(idcomp1) = 0.0d0
         endif
      enddo
      do idcomp1=1,ndim
         if(sum_var .gt. 0.0d0) then
            flit_var1(idcomp1)=(1.0d0-fconst(iveg))*
     $        flit_var1(idcomp1)/sum_var
         else  ! otherwise, no seasonality in litterfall 
            flit_var1(idcomp1)=(1.0d0-fconst(iveg))/float(ndim)
         endif
      enddo     

c (3) fine root turnover during the active growing season
      fmax = 0.0d0
      do idcomp1=1,ndim
         f1 = 0.0d0
         num = 0
         do id=-2,2
            id1 = idcomp1+id
            if((id1 .lt. 1) .or. (id1 .gt. ndim)) cycle
            if(ndviarray(id1) .lt. 0.0d0) cycle
            f1 = f1 + ndviarray(id1)
            num = num + 1
         enddo
         if(num .lt. 5) cycle
         flit_var2(idcomp1) = f1/float(num)
         if(flit_var2(idcomp1) .gt. fmax) then
           fmax=flit_var2(idcomp1)
         endif
      enddo
      sum_var = 0.0d0       
      do idcomp1=1, ndim
         if(fmax .gt. 0.0d0) then 
            flit_var2(idcomp1)=flit_var2(idcomp1)/fmax
            sum_var = sum_var + flit_var2(idcomp1)
         endif
      enddo

      do idcomp1=1, ndim
         if(sum_var .gt. 0.0d0) then 
            flit_var2(idcomp1)=(1.0d0-fconst(iveg))*flit_var2(idcomp1)
     $       /sum_var
         else
            flit_var2(idcomp1)=(1.0d0-fconst(iveg))/float(ndim) ! Otherwise, no seasonality 
         endif
      enddo

c finally, the daily percentage of annual litterfall
      fsum = 0.0d0
      do idcomp1=1, ndim
         frac_above(idcomp1)=(flit_con(idcomp1)+flit_var1(idcomp1))
     $     * fleaf
         frac_below(idcomp1)=(flit_con(idcomp1)+flit_var1(idcomp1)
     $     * 0.5+flit_var2(idcomp1)*0.5)*froot
         fsum = fsum+frac_above(idcomp1)+frac_below(idcomp1)
c         write(108,*) idcomp1,ndviarray(idcomp1),flit_con(idcomp1),
c     $     flit_var1(idcomp1),flit_var2(idcomp1),
c     $     frac_above(idcomp1),frac_below(idcomp1)
      enddo
      if(abs(fsum-1.0d0) .gt. 0.05d0) then
         print*,fsum
         stop 'something wrong with litterfall allocation'
      endif 
c make sure the total litterfall percentage is 1.0     
      do idcomp1=1, ndim
         frac_above(idcomp1)=frac_above(idcomp1)*(1.0d0+
     $      0.5*(1.0d0-fsum)/float(ndim))
         frac_below(idcomp1)=frac_below(idcomp1)*(1.0d0+
     $      0.5*(1.0d0-fsum)/float(ndim)) 
      enddo
      return
      end


      subroutine soc_decomp(iveg,AD,D0,sandfrac,ann_npp,
     $  frac_above,frac_below,kscalar_array,soccpool,Rh_4layer)
c     subroutine calculates the decomposition of the SOC pool, and updates
c     each carbon pool and also calculate the heteotrophic respiration and 
c     NEE fluxes
c     currently, I use the SOC structure (3 litter pools, and 4 soil carbon pools) from
c     BGC and CLM-CN; later on, need to incorporate the original TCF structure?!  
c     note that I add another deep soil carbon pool to account for deep soil carbon below active layer depth 
c     for permafrost area, and below root depth for non-permafrost areas
c
c     Inputs:  iveg       ---- vegetation type
c              AD            ---- accelerating rates for each carbon pool      
c              ann_npp       ---- annual NPP (gC/m2)
c              frac_above    ---- daily aboveground litterfall (gC/m2/day)
c              frac_below    ---- daily belowground litterfall (gC/m2/day) 
c              kscalar_array ---- the decomposition scalar rates 
c
c     Outputs: Rh_4layer --- heterotrophic respiration fluxes (gC/m2/d, or gC/m2/month)
c              soccpool --- update the SOC profile

      implicit none
      INCLUDE 'pwbm.carbon.h'
 
      integer iveg,i,id,ipool,iday,ilayer
      integer npool,ndim,nlayer,numdays
      parameter(npool=7,ndim=46,nlayer=15)
      real*8  ann_npp,Fs,Fd,Fs_daily,Fd_daily,Fs_wood,Fd_wood
      real*8  influx, fsom1_co2,fsom3_co2
      real*8  lit1_pool,lit2_pool,lit3_pool,cwd_pool
      real*8  k_lit1,k_lit2,k_lit3,k_som1,k_som2,k_som3,k_cwd
      real*8  Rh_lit1,Rh_lit2,Rh_lit3,Rh_som1,Rh_som2,Rh_som3,Rh_cwd
      real*8  litterfall,zr,zs,totsoc
      real*8  frac_above(ndim),frac_below(ndim)
      real*8  fwood_cellu,fwood_lignin,sandfrac
      real*8  D0, AD(npool),znode(nlayer),dz(nlayer)
      real*8  delt_t,delt_z2_up, delt_z2_down, alpha
      real*8  soccpool(nlayer,npool),kscalar_array(ndim,nlayer) 
      real*8  Rh_4layer(ndim,nlayer),Rh_arry(nlayer)
      real*8  D_daily(nlayer),XTR(nlayer)
      real*8  ATR(npool,nlayer),BTR(npool,nlayer)
      real*8  CTR(npool,nlayer),DTR(npool,nlayer)

      real*8 fLabile(nvegclass),fCellulose(nvegclass) ! litter quality, differ between grassy and woody components
      real*8 zr_depth(nvegclass),zs_depth(nvegclass)   ! e folding depth (m) of belowground and aboveground litterfall input  
      real*8 fwoody(nvegclass),fLignin(nvegclass)

c  NLCD land cover type: open water(1),Ice/snow(2),Developed(3),Barren_Land(4),
c     Forest(5),Shrub(6),Grassland(7),Cultivated(8),Wetlands(9)                                    
c  also include parameters for 'ice' and 'polar desert',set same as tundra
      data fLabile/0.0D0,0.0D0,0.68D0,0.68D0,0.31D0,0.495D0,
     $       0.68D0,0.68D0,0.68D0/      
      data fCellulose/0.0D0,0.0D0,0.23D0,0.23D0,0.45D0,0.34D0,
     $       0.23D0,0.23D0,0.23D0/
      data fLignin/0.0D0,0.0D0,0.09D0,0.09D0,0.24D0,0.165D0,
     $       0.09D0,0.09D0,0.09D0/
      data fwoody/0.0D0,0.0D0,0.0D0,0.0D0,0.2D0,0.15D0,0.1D0,
     $       0.05D0,0.05D0/
c      data zr_depth/0.0D0,0.0D0,0.25D0,0.25D0,0.25D0,0.25D0,
c     $       0.45D0,0.25D0,0.25D0/  ! Jackson et al. 1996
! this is from Jackson et al. 1996; boreal forest e_folding depth=0.17m
      data zr_depth/0.0D0,0.0D0,0.25D0,0.25D0,0.20D0,0.15D0,
     $       0.11D0,0.20D0,0.11D0/
      data zs_depth/0.0D0,0.0D0,0.1D0,0.1D0,0.1D0,0.1D0,0.1D0,
     $       0.1D0,0.1D0/ ! 8-10cm 
      data znode/0.01D0,0.03D0,0.08D0,0.13D0,0.23D0,0.33D0,0.45D0,
     $       0.55D0,0.70D0,1.05D0,1.40D0,1.75D0,2.25D0,2.75D0,3.25D0/
 
c the decomposition base rate, day-1
      k_lit1 = 0.040303d0
      k_lit2 = 0.01064d0
      k_lit3 = 0.01064d0
      k_cwd  = 0.000659d0
      k_som1 = 0.015647d0
      k_som2 = 0.000443d0
      k_som3 = 0.00001d0     

      fsom3_co2 = 1.0d0

      sandfrac = 0.34d0 
      fsom1_co2 = 0.85d0-0.68d0*(1.0d0-sandfrac)
c      print*,fsom1_co2
c accelerate the decomposition rates
      k_som2 = k_som2 * AD(5)
      k_som3 = k_som3 * AD(6)
      k_cwd  = k_cwd  * AD(7)

c other constants related to litterfall 
      fwood_cellu = 0.7d0
      fwood_lignin = 0.3d0
      zr = zr_depth(iveg)
      zs = zs_depth(iveg)

c constants to solve the SOC profile
      delt_t = 1.0d0  ! in days

c accelerating the SOC decomposition and diffusion rates
      do ilayer=1, nlayer
        if(znode(ilayer) .le. 1.0d0) then 
           D_daily(ilayer) = D0
        else if(znode(ilayer) .gt. 3.0d0) then
           D_daily(ilayer)=0.
        else
           D_daily(ilayer)=D0*(1.0d0-(znode(ilayer)-1.0d0)/2.0d0)
        endif
        if(D_daily(ilayer) .lt. 0.0d0 .or. D_daily(ilayer) .gt. D0) 
     $    stop 'D_daily out of range'
c        print*, ilayer, D_daily(ilayer)
      enddo 

c calculate the depth of each soil layer
      do i=nlayer,2,-1
         dz(i) = znode(i)-znode(i-1)
c         print*, i, dz(i)
      enddo
      dz(1) = znode(2)-znode(1)
c      print*, 1, dz(1)

c make sure there is no missing data 
      if(ann_npp .lt. 0.0d0) stop 'ann_npp < 0'

      do i=1, nlayer
      do ipool=1, npool
         if(soccpool(i,ipool)<0.0d0) stop 'soc<0'
      enddo
      enddo

c track the litter flow 
      litterfall = 0.0d0
c start the decomposition process
      do id=1, ndim

        if(id .lt. ndim) then
           numdays = 8
        else
           numdays = 5
        endif

        do ilayer=1, nlayer
           Rh_4layer(id,ilayer)=0.0d0
        enddo

c litterfall
        Fs_daily=ann_npp*frac_above(id)/float(numdays)
        Fd_daily=ann_npp*frac_below(id)/float(numdays)
c        Fs_daily = ann_npp*0.4d0/365.0d0
c        Fd_daily = ann_npp*0.6d0/365.0d0

        litterfall=litterfall+(Fs_daily+Fd_daily)*float(numdays)

        do iday=1, numdays

c initialization
          do ipool=1, npool
          do ilayer=1,nlayer
            ATR(ipool,ilayer) = -999.0d0
            BTR(ipool,ilayer) = -999.0d0
            CTR(ipool,ilayer) = -999.0d0
            DTR(ipool,ilayer) = -999.0d0
          enddo
          enddo                 

c parameterization for each layer 
          do ilayer=1, nlayer
            Rh_arry(ilayer) = -999.0d0

c aboveground and belowground litterfall, gC/m3
            Fd = Fd_daily/zr * exp(-1.0d0*znode(ilayer)/zr)
            Fs = Fs_daily/zs * exp(-1.0d0*znode(ilayer)/zs)
            Fd_wood = Fd * fwoody(iveg)
            Fs_wood = Fs * fwoody(iveg)

c the CWD pool
            cwd_pool = soccpool(ilayer,7)
            cwd_pool = cwd_pool + Fd_wood + Fs_wood
            Rh_cwd   = cwd_pool * k_cwd * kscalar_array(id,ilayer)
            cwd_pool = cwd_pool - Rh_cwd

c the other three litter pools
            lit1_pool = soccpool(ilayer,1)
            lit2_pool = soccpool(ilayer,2)
            lit3_pool = soccpool(ilayer,3)
            lit1_pool = lit1_pool+(Fs+Fd-Fs_wood-Fd_wood)*fLabile(iveg)
            lit2_pool = lit2_pool+(Fs+Fd-Fs_wood-Fd_wood)
     $        * fCellulose(iveg) + Rh_cwd*fwood_cellu
            lit3_pool = lit3_pool+(Fs+Fd-Fs_wood-Fd_wood)*fLignin(iveg)
     $        + Rh_cwd*fwood_lignin
 
            Rh_lit1 = lit1_pool*k_lit1*kscalar_array(id,ilayer)
            Rh_lit2 = lit2_pool*k_lit2*kscalar_array(id,ilayer)
            Rh_lit3 = lit3_pool*k_lit3*kscalar_array(id,ilayer)

            lit1_pool = lit1_pool - Rh_lit1
            lit2_pool = lit2_pool - Rh_lit2
            lit3_pool = lit3_pool - Rh_lit3

c update the litter pools 
            soccpool(ilayer,1) = lit1_pool
            soccpool(ilayer,2) = lit2_pool
            soccpool(ilayer,3) = lit3_pool
            soccpool(ilayer,7) = cwd_pool

c the respiration from SOM pools
            Rh_som1=soccpool(ilayer,4)*k_som1*kscalar_array(id,ilayer)
            Rh_som2=soccpool(ilayer,5)*k_som2*kscalar_array(id,ilayer)
            Rh_som3=soccpool(ilayer,6)*k_som3*kscalar_array(id,ilayer)

c the total soil respiration at this depth --- need introduce soil texture here          
            Rh_arry(ilayer) = Rh_lit1*0.55d0+(Rh_lit2+Rh_lit3)*0.5d0
     $       +fsom1_co2*Rh_som1+0.55d0*Rh_som2+fsom3_co2*Rh_som3
 
c update the parameters for the three SOM pools
            do ipool=4,6

              if(ilayer .eq. 1) then  ! surface layer 
                delt_z2_up = 2.0d0/(dz(ilayer+1)*
     $            (dz(ilayer)+dz(ilayer+1)))

                alpha = D_daily(ilayer)*AD(ipool)*delt_t*delt_z2_up/
     $            (365.0d0*10000.0d0)

                ATR(ipool,ilayer) = 0.0d0
                BTR(ipool,ilayer) = 1.0d0+2.0d0*alpha
                CTR(ipool,ilayer) = -2.0d0*alpha
              else if(ilayer .eq. nlayer) then  ! bottom layer 
                delt_z2_down = 2.0d0/(dz(ilayer)*
     $            (dz(ilayer)+dz(ilayer))) ! should it be 'ilayer+1'?
                
                alpha = D_daily(ilayer)*AD(ipool)*delt_t*delt_z2_down/
     $            (365.0d0*10000.0d0)

                ATR(ipool,ilayer) = -2.0d0*alpha
                BTR(ipool,ilayer) = 1.0d0+2.0d0*alpha
                CTR(ipool,ilayer) = 0.0d0
              else  ! intermediate layer 
                delt_z2_up = 2.0d0/(dz(ilayer+1)*
     $            (dz(ilayer)+dz(ilayer+1)))
                delt_z2_down = 2.0d0/(dz(ilayer)*
     $            (dz(ilayer)+dz(ilayer+1)))
           
                ATR(ipool,ilayer)=-1.0d0*D_daily(ilayer)*AD(ipool)
     $            *delt_t*delt_z2_down/(365.d0*10000.0d0)
                BTR(ipool,ilayer)=1.0d0+D_daily(ilayer)*AD(ipool)
     $            *delt_t*(delt_z2_up+delt_z2_down)/(365.0d0*10000.0d0)
                CTR(ipool,ilayer)=-1.0d0*D_daily(ilayer)*AD(ipool)
     $            *delt_t*delt_z2_up/(365.0d0*10000.0d0)
              endif

c  the carbon inputs/outputs for the three SOM pools 
             if(ipool .eq. 4) then  ! the active SOM pool
                influx = 0.45d0*Rh_lit1+0.50d0*Rh_lit2
     $            +0.42d0*Rh_som2+(1.0d0-fsom3_co2)*Rh_som3
                DTR(ipool,ilayer)=soccpool(ilayer,ipool)+
     $            (influx-Rh_som1)*delt_t

                if(DTR(ipool,ilayer) < 0.0d0) then
                  print*,ipool,ilayer,DTR(ipool,ilayer) 
                endif
 
             else if(ipool .eq. 5) then ! the slow SOM pool
                influx = 0.5d0*Rh_lit3+Rh_som1*(1.0d0-0.004d0-fsom1_co2)
                DTR(ipool,ilayer)=soccpool(ilayer,ipool)+
     $            (influx-Rh_som2)*delt_t
 
                if(DTR(ipool,ilayer) < 0.0d0) then
                  print*,ipool,ilayer,DTR(ipool,ilayer)
                endif

             else  ! the passive SOM pool
                influx = Rh_som1*0.004d0+Rh_som2*0.03d0
                DTR(ipool,ilayer)=soccpool(ilayer,ipool)+
     $            (influx-Rh_som3)*delt_t

                if(DTR(ipool,ilayer) < 0.0d0) then
                  print*,ipool,ilayer,DTR(ipool,ilayer)
                endif

             endif
          enddo ! loop for each SOM pool  
c         write(109,*) ilayer,kscalar_array(id,ilayer),Fd,Fs,Rh_lit1,
c     $     Rh_lit2,Rh_lit3,Rh_som1,Rh_som2,Rh_som3,
c     $     (soccpool(ilayer,ipool),ipool=1,7) 
        enddo  ! loop for each soil layer
        
c solve the tridiagonal matrix and update the SOM pools 
        do ipool=4,6
          do ilayer=1, nlayer
            XTR(ilayer) = -999.0d0
          enddo

          call TDMA(nlayer,ATR(ipool,:),BTR(ipool,:),CTR(ipool,:),
     $      DTR(ipool,:),XTR)

          do ilayer=1, nlayer
            soccpool(ilayer,ipool) = XTR(ilayer)
            if(XTR(ilayer)<0.0d0) stop 'XTR<0'
          enddo
        enddo

c update the C flux
        do ilayer=1, nlayer
          if(Rh_arry(ilayer) < 0.0d0) stop 'Rh < 0.'

          Rh_4layer(id,ilayer)=Rh_4layer(id,ilayer)+Rh_arry(ilayer)
        enddo 
      
      enddo ! loop for each day
      enddo !loop for each time step      

      if(abs(litterfall-ann_npp)>1.0d0) stop 'double check litterfall'
      
      return
      end      

      subroutine TDMA(nlayer,AR,BR,CR,DR,XR)
c     use the Thomas algorithm to solve the tridiagonal algorithm TDMA)
      
      implicit none
      
      integer i, nlayer
      real*8 AR(nlayer),BR(nlayer),CR(nlayer),DR(nlayer),XR(nlayer)
      real*8 m
 
      do i=1,nlayer
       if((abs(BR(i)).le.abs(AR(i))+abs(CR(i)))
     & .or.(DR(i) .lt. 0.0d0))  then
          print*, i,AR(i),BR(i),CR(i),DR(i)
          stop 'cannot find solutions for the tridiagonal matrix '
        endif
      enddo  

      do i=2,nlayer
        m = AR(i)/BR(i-1)
        BR(i) = BR(i) - m*CR(i-1)
        DR(i) = DR(i) - m*DR(i-1)
      enddo   

      XR(nlayer) = DR(nlayer)/BR(nlayer)
     
      do i=nlayer-1,1,-1
        XR(i) = (DR(i)-CR(i)*XR(i+1))/BR(i)
      enddo

      do i=1,nlayer
        if(XR(i)<0.0d0) stop 'XR<0!'
      enddo

      return
      end


      subroutine kdecomp_scalar(dzmm,porosity4layer,wat4_layer,
     $  ice4_layer,fliq_frac,soil_temp,kscalar_4layer,nlayer)

c     June 22, 2017: I revised this code, to make it simple

c     subroutine calculates the decomposition scalar from soil temperature,
c     soil moisture, and layer depth for three different soil zones (surface,
c     root-zone, deep soil)
c
c     Inputs:   wat4_layer     --- is this the total water content or ice content?
c               ice4_layer     --- ice content for each layer
c               soil_temp      --- soil temperature for each layer (deg C) 
c               dzmm           --- thickness of each soil layer (mm)
c               porosity4layer --- soil porosity (cm3/cm3)
c               fliq_frac      --- liquid water fraction   
c               nlayer         --- the number of soil layers accounting C dynamics       
c     Outputs:  
c               kscalar_4layer --- decomposotion rate scalar for each layer, range from 0-1
     
      implicit none
c      INCLUDE 'pwbm.carbon.h'
      INCLUDE 'pwbm.soil_temperature.h'

      integer  i, nlayer
      real*8   a1,a2,b1,b2,b3
      real*8   alpha1,alpha2
      real*8   topt,sm_opt,e_depth  
      real*8   tsoil,swc,tmult,wmult,zmult,kscalar
      real*8   kscalar_4layer(15),dzmm(lx)
      real*8   wat4_layer(sx:lx), ice4_layer(sx:lx),capacity4layer(lx)! mm 
      real*8   porosity4layer(lx) ! unit: cm3/cm3
      real*8   soil_temp(sx:lx+1),fliq_frac(sx:lx+1)   ! note that soil temperature is a bit different
      real*8   znode(15)
      parameter (e_depth=1.0d0)  ! e-folding depth for SOC decomposition
      parameter (a1=308.56d0,topt=25.0d0) ! should I add these to the header file?!
c      parameter (sm_opt=0.60d0,alpha1=2.77778,alpha2=1.0d0) !empirical parameters!!!
      parameter (sm_opt=0.30d0,alpha1=11.1111d0,alpha2=1.0d0) 
      data znode/0.01d0,0.03d0,0.08d0,0.13d0,0.23d0,0.33d0,0.45d0,
     $       0.55d0,0.70d0,1.05d0,1.40d0,1.70d0,2.25d0,2.75d0,3.25d0/

c constants for temperature constraint 
      a2 = 66.02d0+topt-20.0d0

      do i=1,nlayer
c note that root fraction will approach 0, but never = 0 (follow exponential decay curve)
         tsoil = soil_temp(i+1)  ! deg C, soil temperature array index range from 2 to lx+1..

         if(fliq_frac(i+1) < CR_LIQ_FRACTION) then
c            swc = ice4_layer(i)*fliq_frac(i+1)
            swc = ice4_layer(i)
         else
            swc = wat4_layer(i)  ! mm
         endif
      
         if((tsoil .le. -999.0d0) .or. (swc .lt. 0.0d0)) then
            print *, i, tsoil, swc
            stop 'check soil temperature and moisture data '
         endif 

c should I use field capacity or porosity?
         capacity4layer(i) = dzmm(i) * porosity4layer(i)! porosity4layer: cm3/cm3; capacity4layer: mm
c convert swc (mm) to saturation raction (%)
c here I think I should use field capacity, instead of saturation point...
         if(swc .gt. capacity4layer(i)) then
            swc = 1.0d0 ! this should not occur,just in case
         else 
            swc = swc/capacity4layer(i)
         endif

c the temperature constraint on decomposition
         if(tsoil .ge. topt) then
            tmult = 1.0d0
         else if(tsoil .lt. -10.0d0) then
            tmult = 0.0d0
         else
            tmult = exp(a1*(1.0d0/a2-1.0d0/(a2+tsoil-topt)))
         endif

         if(tmult.lt.0.0d0.or.tmult.gt.1.0d0) stop 'unreasonable tmult'

         if(swc .le. sm_opt) then
             wmult = 1.0d0 - alpha1*(swc-sm_opt)*(swc-sm_opt)
         else
c             wmult = 1.0d0 - alpha2*(swc-sm_opt)*(swc-sm_opt)
             wmult = 1.0d0
         endif

c impose additional constraints on decomposition (unfrozen water
content)
        wmult = wmult * fliq_frac(i+1)

c use different functions for unfrozen water content---double check this
c         if(fliq_frac(i+1) < CR_LIQ_FRACTIOn) then
c             wmult = fliq_frac(i+1)
c         else
c             wmult = 1.0d0
c         endif
c          if(tsoil < 0.0d0) then
c             wmult =((0.1d0-tsoil)/0.01d0)**(-0.5d0)/(10.0d0**(-0.5d0))
c          else
c             wmult = 1.0d0
c          endif

c         if(i .eq.2) write(113,*) i,tsoil,tmult,fliq_frac(i+1),swc,wmult
           
c         print*, tsoil,tmult,fliq_frac(i+1),wmult

         if(wmult .lt.0.0d0 .and. wmult.gt.-0.001d0) then
             wmult = 0.0d0
         endif

         if(wmult.lt.0.0d0.or.wmult.gt.1.0d0) then
             print *, dzmm(i),wat4_layer(i),ice4_layer(i),
     $          fliq_frac(i+1),swc, wmult
             stop 'unreasonable wmult'
         endif
  
c the depth related reduction in decomposition rates 
         zmult = exp(-1.0d0*znode(i)/e_depth)  !0.5m: e folding depth; need to adjust this parameter   
          
         kscalar = tmult*wmult       
c         kscalar = tmult*wmult*zmult

         if((kscalar .lt. 0.0d0) .or. (kscalar .gt. 1.0d0)) then
            stop 'error with calculating environmental constraint'
         endif

         kscalar_4layer(i) = kscalar
      enddo ! loop for each soil layer
      
      return 
      end

      subroutine calSOCContent(soccpool,totsoc)
c     Subroutine calculate the total SOC content for the whole profile (0-3m)
 
c     Inputs:  soild --- the node depth for each soil layer
c              soccpool  --- the soc concentration(kg/m3) for each layer

c     Outputs: totsoc --- total soc content (kg/m2) for layer 0-3m

      implicit none
      include 'pwbm.carbon.h'
      include 'pwbm.soil_temperature.h' 

      integer  i,ipool,npool,nlayer
      parameter(npool=7,nlayer=15)
      real*8   totsoc, socc
      real*8   soccpool(nlayer,npool)
      real*8   znode(1:nlayer)  ! only include the soil layers from 0 to 3m 
      data znode/0.01d0,0.03d0,0.08d0,0.13d0,0.23d0,0.33d0,0.45d0,
     $      0.55d0,0.70d0,1.05d0,1.40d0,1.75d0,2.25d0,2.75d0,3.25d0/

      totsoc = 0.0d0
      do ipool=1, npool
         totsoc = totsoc + (soccpool(1,ipool)*znode(1)) ! the first layer      
         do i=1, nlayer-1
           socc = (soccpool(i,ipool) + soccpool(i+1,ipool))/2.0d0
           if(socc < 0.0d0) stop 'socc < 0.'
           totsoc = totsoc + socc*(znode(i+1)-znode(i))
         enddo
      enddo     

      return
      end         

      subroutine gpp_model(iveg,ndviarray,tminarray,dswrfarray,
     $   rtzswcarray,gpparray,npparray)
c     Subroutine calculates daily gpp from meteorology and ndvi inputs 
c     note even for the desert and ice, I still include them in the calculation, but 
c     set LUE very low --- otherwise, there would be no soc at all...
c     revised on June, 2017 ---- run the model at 8-day time step 

c     Inputs:  tmin -- 8-day minimum temperature (deg C)
c              rootzswc -- 8-day soil saturation degree for the root zone (0-1)
c              dswrf  -- 8-day downward solar radiation (MJ/m2/day)
c              ndvi -- 8-day NDVI/FPAR (unitless)    
c              iveg -- vegetation type

c     Outputs: gpp --- 8-day gpp, unit: gC/m2/day
c              npp --- 8-day npp, unit: gC/m2/day ! npp should be defined at annual, but anyway...                
       
      implicit none
      include 'pwbm.carbon.h'

      integer   id,iveg,ndays,ndim
      parameter(ndim=46)
      real*8    ndvi,tmin,dswrf,rootzswc,swc_opt,lue
      real*8    tmult,wmult,rmissing
      real*8    b1,b2,b3,sm_opt
      real*8    gpp,npp
      real*8    ndviarray(ndim),tminarray(ndim)
      real*8    dswrfarray(ndim),rtzswcarray(ndim)
      real*8    gpparray(ndim),npparray(ndim) 
      parameter (rmissing=-9999.0d0)  ! there should be no missing data though
      parameter (sm_opt=0.60d0,b1=9.90d0,b2=6.13d0,b3=2.50d0)

      real*8  Tmin_min(nvegclass)   ! minimum Tmin constraint (deg C)
      real*8  Tmin_max(nvegclass)   ! maximum Tmin constraint (deg C)
      real*8  VPD_min(nvegclass)    ! minimum VPD constraint  (pa)
      real*8  VPD_max(nvegclass)    ! maximum VPD constraint  (pa)
      real*8  RootSM_min(nvegclass) ! minimum Rootzone soil moisture constraint (%)
      real*8  RootSM_max(nvegclass) ! maximum Rootzone soil moisture constraint (%)
 
      real*8  LUEmax(nvegclass)    ! maximum light use efficiency, kgC/(MJ/m2)
      real*8  CUE(nvegclass)       ! carbon use efficiency (unitless)
  
c    NLCD land cover type: 1: open water; 2: perennial Ice/Snow
c    3: developed; 4: barren land; 5: forest; 6: shrub; 7: grassland;
c    8: cultivated; 9: wetlands
c    the new lookup table for the gpp model
c    here, I treat the "grassland" as "tundra"; "shrub" as "Forest-tundra"
      data LUEmax/0.0D0, 0.0D0, 0.0001D0, 0.0001D0, 0.001055D0,
     $  0.00078D0, 0.0005D0, 0.0012D0, 0.001D0/
      data Tmin_min/-8.0D0, -8.0D0, -8.0D0, -8.0D0, -8.0D0, 
     $  -8.0D0, -8.0D0, -8.0D0, -8.0D0/
      data Tmin_max/12.0D0, 12.0D0, 12.0D0, 12.0D0,8.31D0, 
     $  8.5D0, 8.8D0, 12.02D0, 12.02D0/
      data VPD_min/500.0D0, 500.0D0, 500.0D0, 500.0D0, 500.0D0, 
     $  500.0D0, 500.0D0, 752.0D0, 752.0D0/
      data VPD_max/5000.0D0, 5000.0D0, 5000.0D0, 5000.0D0, 4000.0D0,
     $  4000.0D0, 4455.0D0, 5500.0D0, 5500.0D0/
      data RootSM_min/0.15D0,0.15D0,0.15D0,0.15D0,0.15D0,0.15D0,0.15D0,
     $  0.15D0, 0.15D0/
      data RootSM_max/0.75D0,0.75D0,0.75D0,0.75D0,0.75D0,0.75D0,0.75D0,
     $  0.75D0, 0.75D0/
      data CUE/0.60D0, 0.60D0, 0.60D0, 0.60D0, 0.55D0, 0.575D0, 0.6D0,
     $  0.60D0, 0.60D0/  
 
c    if not vegetated lands (still calculate GPP/NPP for the land cover
c    type #3 or #4)
      if(iveg .le. 2) then
         print *, iveg
         stop
      endif

      do id=1, ndim

         ndvi = ndviarray(id)
         tmin = tminarray(id)
         rootzswc = rtzswcarray(id)
         dswrf = dswrfarray(id)

c    in case there is any missing data in the inputs --- this should not occur though 
         if( tmin.le.rmissing .or. rootzswc.lt.0.0d0 
     $    .or. rootzswc .gt. 1.0d0 .or. dswrf .lt. 0.0d0  
     $    .or. ndvi .le. -1.0d0.or. ndvi .gt. 1.0d0) then
            stop 'double check the gpp model inputs'
         endif
c     sometimes ndvi <0. (generally water), set ndvi = 0.01  
         if(ndvi < 0.)  ndvi = 0.01  ! just to assign a small value 
         dswrf = dswrf * 0.0864      ! from W/m2 to MJ/m2/day

c     the temperature constraint
         if(tmin .lt. Tmin_min(iveg) ) then
           tmult = 0.0d0
         else if(tmin .gt. Tmin_max(iveg)) then
           tmult = 1.0d0
         else 
           tmult=(tmin-Tmin_min(iveg))/(Tmin_max(iveg)-Tmin_min(iveg))
         endif

c    the rootzone soil moisture constraint(wmult) --- there is large uncertainty regarding this
         swc_opt = 0.50d0
         if(rootzswc .gt. swc_opt) then
           wmult = 1.0d0
         else
           wmult = (1.0d0+rootzswc)/(1.0d0+swc_opt)
         endif

c     the light use efficiency
c      lue = LUEmax(iveg) * tmult * min(wmult, vmult)
         lue = LUEmax(iveg) * tmult * wmult ! how much impact will be from adding a root-zone soil moisture

c     the daily GPP (gC/m2/8day)
         if(id .lt. ndim) then
           ndays = 8.0d0
         else 
           ndays = 5.0d0
         endif
         gpp = ndays * lue * dswrf * 0.45d0 * ndvi * 1000.0d0 ! gC/m2/8day
c     according to npp definition, this should be defined at annual scale..but here
c     just for the sake of simplification
         npp = gpp * CUE(iveg) 
         gpparray(id) = gpp
         npparray(id) = npp
      enddo

      return
      end

      



