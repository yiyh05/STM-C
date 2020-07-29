      INTEGER LX             !Number of soil layer
      INTEGER SX             !Maximal number of snow layer

      REAL*8 CW              ![MJ/m/C]       !Heat capacity of liqiud water
      REAL*8 CI              ![MJ/m/C]       !Heat capacity of ice
      REAL*8 CA              ![MJ/m/C]       !Heat capacity of dry air
      REAL*8 CSNOW           ![MJ/m/C]       !Heat capacity of snow, assumed to be constant

      REAL*8 KW              ![W/m/C]       !thermal conductivity of liqiud water
      REAL*8 KI              ![W/m/C]       !thermal conductivity of ice
      REAL*8 KA              ![W/m/C]       !thermal conductivity of dry air

      REAL*8 LA              ![MJ/m3]          !Latent heat of fusion
      REAL*8 TSCALE          ![sec/10^6]       !Time scale
      REAL*8 CR_LIQ_FRACTION ![%]              !Critical liquid water fraction under which the soil is assumed to be thawed


      INTEGER NUM_MAX_ITER   ![5-9]            !Maximal number of allowed iterations before decreasing the time step
      REAL*8  MAX_DT         ![days]           !Maximal time step
      REAL*8  MAX_DELTA      ![C]              !Accuracy at each time step
      REAL*8  DLT            ![C]              !Regularization of computing a numerical derivative of the ethalpy

      REAL*8  XD             ![?]              !Regularization of unfrozen water content
      REAL*8  XS             ![?]              !Regularization of unfrozen water content
      REAL*8  XS_2                             !To reduce a lot of division by XS*XS

       
      REAL*8 LKWKI                             !Natural logarithm of KW/KI
      
       
      PARAMETER(LX              = 23, !145
     &          SX              =-5, !-55
     &
     &          NUM_MAX_ITER    = 7,
     &          MAX_DT          = 1.0D0,
     &          DLT             = 1.0D-10,
     &          MAX_DELTA       = 1.0D-3,
     &
     &          XD              = 3.0D-2,
     &          XS              = 2.0D1,
     &          XS_2            = 1.0D0/(XS*XS),
     &          CR_LIQ_FRACTION = 0.99D0,
     &
     &   
     &          CW              = 4.181D0,
     &          CI              = 2.114D0,
     &          CA              = 1.003D-3,
!     &          CSNOW           = 0.84D0,
     &   
     &          KA              = 0.025D0,
     &          KW              = 0.6D0, ! the original value is 0.56
     &          KI              = 2.2D0,
     &          LKWKI           =-1.36827585561721230D0, !=DLOG(KW/KI)
     &
     &          LA              = 333.2D0,
     &          TSCALE          = 1.0D6/691200.0D0) ! 8-day 
!     &          TSCALE         = 1.0D6/86400.0D0) ! daily
       
       

      !Every thing below is just to compute temperature fast
      REAL*8 SURFACE_TEMPERATURE_BEGIN
      REAL*8 SURFACE_TEMPERATURE_END


      REAL*8 XREF(SX:LX), DX_1(SX:LX), DIAG(SX:LX)
      COMMON /SOIL_GEOMETRY/   DX_1, DIAG, XREF

      REAL*8 LAYER_W0(SX:LX+1),LAYER_W1(SX:LX+1)
      REAL*8 LAYER_LF(SX:LX+1),LAYER_LT(SX:LX+1),LAYER_K(SX:LX+1)
      REAL*8 LAYER_CF(SX:LX+1),LAYER_CT(SX:LX+1),LAYER_C(SX:LX+1)

      REAL*8 LAYER_B(SX:LX+1)
      REAL*8 LAYER_Tp(SX:LX+1),LAYER_Tp_1(SX:LX+1)

      REAL*8 LAYER_S(SX:LX+1), LAYER_A(SX:LX+1)
      REAL*8 LAYER_L(SX:LX+1), LAYER_I(SX:LX+1)


      COMMON /SOIL_ITERATIONS/ LAYER_W0,LAYER_W1,LAYER_B,
     &                         LAYER_Tp,LAYER_Tp_1,
     &                         LAYER_LF,LAYER_LT,LAYER_K,
     &                         LAYER_CF,LAYER_CT,LAYER_C,
     &                         LAYER_S,LAYER_A,LAYER_L,LAYER_I,
     &                         SURFACE_TEMPERATURE_BEGIN,
     &                         SURFACE_TEMPERATURE_END
