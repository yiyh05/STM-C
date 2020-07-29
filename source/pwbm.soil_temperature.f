      SUBROUTINE STEMPERATURE(SOIL_TEMPERATURE,LIQUID_FRACTION,
     &                        TIME_STEP,SN,GC)
      IMPLICIT NONE
      INCLUDE 'pwbm.soil_temperature.h'
         
      REAL*8 SOIL_TEMPERATURE(SX:LX+1)        !Soil temperature       (C)   [In/Out]
      REAL*8 LIQUID_FRACTION(SX:LX+1)         !Liquid water fraction  (%)   [In/Out]
      REAL*8 TIME                             !Start of calculations  (days)[In]
      REAL*8 TIME_END                         !End of calculations    (days)[In]
      REAL*8 TIME_STEP                        !Current time step             (days)[In/Out]
      INTEGER SN                              !Number of the first node in the grid
      INTEGER GC                              !Counter of converged iterations


      REAL*8 OLD_TIME_STEP,DT_TEMP

      time = 0.0d0
      time_end = 1.0d0
      CALL STEMPERATURE_ATOMIC(SOIL_TEMPERATURE,LIQUID_FRACTION,
     &                         TIME_STEP,TIME,TIME_END,SN,GC)
      OLD_TIME_STEP=TIME_STEP
      DT_TEMP=TIME_END-TIME
      DO WHILE(1.0D-5.LT.DT_TEMP)
        CALL STEMPERATURE_ATOMIC(SOIL_TEMPERATURE,LIQUID_FRACTION,
     &                           DT_TEMP,TIME,TIME_END,SN,GC)
        DT_TEMP=TIME_END-TIME
      ENDDO
      TIME_STEP=OLD_TIME_STEP

      END 






      SUBROUTINE STEMPERATURE_ATOMIC(SOIL_TEMPERATURE,
     &                               LIQUID_FRACTION,
     &                               DT,TIME,TIME_END,SN,GC)
      IMPLICIT NONE
      INCLUDE 'pwbm.soil_temperature.h'

      INTEGER SN,GC
      REAL*8 SOIL_TEMPERATURE(SX:LX+1)
      REAL*8 LIQUID_FRACTION(SX:LX+1)

      INTEGER IT,I
      INTEGER FLAG

      REAL*8 TIME,TIME_END,DT
      REAL*8 T(SX:LX+2)
      REAL*8 Wc, Wp, dTemperature, Capp_water, Capp

      REAL*8 DELTA
      REAL*8 C1,B1,D1,R1,L1
      REAL*8 TL(SX:LX+1),   TR(SX:LX+1)
      REAL*8 ATR(SX:LX+1),  BTR(SX:LX+1),  CTR(SX:LX+1)
      REAL*8 ATRNR(SX:LX+1),BTRNR(SX:LX+1),CTRNR(SX:LX+1)
      REAL*8 DTR(SX:LX+1),  DTRNR(SX:LX+1)

      REAL*8 C(SX:LX+1), DC_DT(SX:LX+1)                       !Heat capacity
      REAL*8 H(SX:LX+1), DH_DT(SX:LX+1), D2H_D2T(SX:LX+1)     !Enthalpy due to water
      REAL*8 K(SX:LX+1), DK_DT(SX:LX+1)                       !Thermal conductivity
      REAL*8 Y(SX:LX+1), DY_DT(SX:LX+1), D2Y_D2T(SX:LX+1)     !Liquid water fraction

      REAL*8 SURFACE_TEMPERTURE



      DO WHILE(TIME+DT.LE.TIME_END)
        DO I=SN,LX+1
          T(I)=SOIL_TEMPERATURE(I)
        ENDDO
      
        FLAG=1
        IT=0
        DO WHILE (FLAG.EQ.1)
          DO I=SN,0
            CALL SNOW_TH_PROPERTIES(T(I),Y(I),DY_DT(I),D2Y_D2T(I), 
     &                              H(I),DH_DT(I),D2H_D2T(I),
     &                              C(I),DC_DT(I), 
     &                              K(I),DK_DT(I),
     &               LAYER_W1(I),LAYER_CF(I),LAYER_CT(I),LAYER_LF(I),
     &               LAYER_LT(I),LAYER_B(I),LAYER_Tp(I),LAYER_Tp_1(I))
          ENDDO
          DO I=2,LX+1
            CALL SOIL_TH_PROPERTIES(T(I),Y(I),DY_DT(I),D2Y_D2T(I), 
     &                              H(I),DH_DT(I),D2H_D2T(I),
     &                              C(I),DC_DT(I), 
     &                              K(I),DK_DT(I),
     &               LAYER_W1(I),LAYER_CF(I),LAYER_CT(I),LAYER_LF(I),
     &               LAYER_LT(I),LAYER_B(I),LAYER_Tp(I),LAYER_Tp_1(I))

c            print*, I, LAYER_LF(I), LAYER_LT(I)
c            if(I .eq. 3) then
c            if(T(I) .lt. 0.0d0) then
c             write(103,*) T(I),Y(I),C(I),K(I),LAYER_CF(I),LAYER_LF(I)
c            else
c             write(104,*) T(I),Y(I),C(I),K(I),LAYER_CT(I),LAYER_LT(I)
c            endif
c            endif

c            if(I .eq. 10) then
c            if(T(I) .lt. 0.0d0) then
c             write(105,*) T(I),Y(I),C(I),K(I),LAYER_CF(I),LAYER_LF(I)
c            else
c             write(106,*) T(I),Y(I),C(I),K(I),LAYER_CT(I),LAYER_LT(I)
c            endif
c            endif
          ENDDO

          !Snow model
          DO I=SN,-1
            TL(I)=(K(I)+K(I+1))*DT*DX_1(I)*0.5D0

            BTR(I)  = TL(I)
            CTR(I)  =-TL(I)
            ATR(I+1)=-TL(I)

            BTRNR(I)  = TL(I)
            CTRNR(I)  =-TL(I)
            ATRNR(I+1)=-TL(I)
          ENDDO
          BTR(0)=0.0D0
          BTRNR(0)=0.0D0

          !Flux
          CTR(0)= 1.0D0
          CTRNR(0)= 1.0D0
          !Mortar
          BTR(1)= 0.0D0
          BTRNR(1)= 0.0D0
          ATR(1)= 1.0D0
          ATRNR(1)= 1.0D0
          CTR(1)=-1.0D0
          CTRNR(1)=-1.0D0
          DTR(1)= 0.0D0      
          !Flux
          ATR(2)=-1.0D0; ATRNR(2)=-1.0D0

          ATR(SN)=0.0D0
          ATRNR(SN)=0.0D0
          !Soil
          DO I=2,LX
            TL(I) = (K(I)+K(I+1))*DT*DX_1(I-1)*0.5D0
            TR(I) = (T(I)-T(I+1))*DT*DX_1(I-1)*0.5D0

            BTR(I)  = TL(I)
            CTR(I)  =-TL(I)
            ATR(I+1)=-TL(I)

            BTRNR(I)  = TL(I)+TR(I)*DK_DT(I)
            CTRNR(I)  =-TL(I)+TR(I)*DK_DT(I+1)
            ATRNR(I+1)=-TL(I)-TR(I)*DK_DT(I)
          ENDDO
          BTR(LX+1)=0.0D0
          BTRNR(LX+1)=0.0D0
          CTR(LX+1)=0.0D0
          CTRNR(LX+1)=0.0D0

          !Snow
          DO I=SN,-1
            BTR(I+1)  =BTR(I+1)  +TL(I)
            BTRNR(I+1)=BTRNR(I+1)+TL(I)
          ENDDO
          !Soil
          DO I=2,LX
            BTR(I+1)  =BTR(I+1)  +TL(I)
            BTRNR(I+1)=BTRNR(I+1)+TL(I)-TR(I)*DK_DT(I+1)
          ENDDO

          !Snow            
          DO I=SN,0
            BTR(I)=BTR(I)+C(I)*DIAG(I)
            BTRNR(I)=BTRNR(I)+C(I)*DIAG(I)
            DTR(I)=C(I)*DIAG(I)*SOIL_TEMPERATURE(I)
          ENDDO
          !Soil
          DO I=2,LX+1
            dTemperature=T(I)-SOIL_TEMPERATURE(I)
            IF(ABS(dTemperature)<dlt) THEN                        !Second term in the apparent heat capacity due to freezing water
              Capp_water=DH_DT(I)                                 !if temperature in the soil layer is the same we approximate it analytically
              BTRNR(I)=BTRNR(I)+D2H_D2T(I)*DIAG(I-1)*dTemperature
            ELSE
              Wp=LAYER_W0(I)*LIQUID_FRACTION(I)*La*TSCALE          !Enthalpy due to unfrozen water at the previous step
              Wc=H(I)                                             !Enthalpy due to unfrozen water at the current step
              Capp_water=(Wc-Wp)/dTemperature                     !Otherwise we use an approximation
              BTRNR(I)=BTRNR(I)- (Capp_water-DH_DT(I))*DIAG(I-1)
            ENDIF
            Capp=C(I)+Capp_water                                  !Apparent heat capacity

            BTR(I)=BTR(I)+Capp*DIAG(I-1)
            BTRNR(I)=BTRNR(I)+Capp*DIAG(I-1)
     &              +DC_DT(I)*DIAG(I-1)*(T(I)-SOIL_TEMPERATURE(I))

            DTR(I)=Capp*DIAG(I-1)*SOIL_TEMPERATURE(I)
          ENDDO

          DO I=2,LX+1
              DTR(I)=DTR(I)+(LAYER_W1(I)-LAYER_W0(I))*DT*DIAG(I-1)
     &                                             *La*TSCALE
          ENDDO

          R1=SURFACE_TEMPERTURE(TIME+DT)
                  
          IF(IT<3) THEN  ! yyh: what does this mean?
            C1=CTR(SN)
            B1=BTR(SN)
            D1=DTR(SN) 
            CTR(SN)=0.0D0
            BTR(SN)=1.0D0
            DTR(SN)=R1
                  
            CALL TRIDIAG(SN,SX,LX+1,ATR,BTR,CTR,DTR)
            L1=D1-(B1*DTR(SN)+C1*DTR(SN+1))
                  
            DELTA=1.0D0
!            DELTA=0.0D0
!            DO I=SN,LX+1
!              delta=DMAX1(DABS(DTR(I)-T(I)),DELTA)
!            ENDDO

            DO I=SN,LX+1
              T(I)=DTR(I)
            ENDDO
            T(LX+2)=L1
          ELSE
            DTRNR(SN)=BTR(SN)*T(SN)+CTR(SN)*T(SN+1)+T(LX+2)-DTR(SN)
            DO I=SN+1,LX
              DTRNR(I)=ATR(I)*T(I-1)+BTR(I)*T(I)
     &                +CTR(I)*T(I+1)-DTR(I)
            ENDDO
            DTRNR(LX+1)=ATR(LX+1)*T(LX)+BTR(LX+1)*T(LX+1)-DTR(LX+1)

            C1=CTRNR(SN)
            B1=BTRNR(SN)
            D1=DTRNR(SN)
            CTRNR(SN)=0.0D0 
            BTRNR(SN) =1.0D0
            DTRNR(SN) =T(SN)-R1
            CALL TRIDIAG(SN,SX,LX+1,ATRNR,BTRNR,CTRNR,DTRNR)
            L1=D1-(B1*DTRNR(SN)+C1*DTRNR(SN+1))

            DELTA=0.0D0
            DO I=SN,LX
              T(I)=T(I)-DTRNR(I)
              delta=DMAX1(DABS(DTRNR(I)),DELTA)
            ENDDO
              T(LX+2)=T(LX+2)-L1
          ENDIF

C          WRITE(*,'(F10.6,F16.13 ,I3,F20.16,I4)') TIME+DT,DELTA,IT,DT,GC
          IF(DELTA.LT.MAX_DELTA) THEN
            !PRINT *,time,IT,DT
            TIME=TIME+DT
            GC=GC+1
            IF(GC.EQ.5) THEN
              DT=DT*2.0D0
              GC=0
            ENDIF
            DT=DMIN1(DT,MAX_DT)
            IT=0
            FLAG=0
          ELSE
            IT=IT+1
          ENDIF
          IF(IT.EQ.NUM_MAX_ITER) THEN
            DT=DT/2.0D0
            DO I=SN,LX+1
              T(I)=SOIL_TEMPERATURE(I)
            ENDDO
            GC=0
            IT=0
          ENDIF
        ENDDO


        DO I=SN,0
          SOIL_TEMPERATURE(I)=T(I)
          CALL SNOW_TH_PROPERTIES(SOIL_TEMPERATURE(I),
     &                            LIQUID_FRACTION(I),
     &                            DY_DT(I),D2Y_D2T(I), 
     &                            H(I),DH_DT(I),D2H_D2T(I),
     &                            C(I),DC_DT(I), 
     &                            K(I),DK_DT(I),
     &             LAYER_W1(I), LAYER_CF(I),LAYER_CT(I),LAYER_LF(I),
     &             LAYER_LT(I),LAYER_B(I), LAYER_Tp(I),LAYER_Tp_1(I))
          LAYER_L(I)=LAYER_W1(I)*LIQUID_FRACTION(I)
          LAYER_I(I)=LAYER_W1(I)*(1.0D0-LIQUID_FRACTION(I))
          LAYER_K(I)=K(I)
          LAYER_C(I)=C(I)/TSCALE
        ENDDO
        SOIL_TEMPERATURE(1)=T(1)
          DO I=2,LX+1
          SOIL_TEMPERATURE(I)=T(I)
          CALL SOIL_TH_PROPERTIES(SOIL_TEMPERATURE(I),
     &                            LIQUID_FRACTION(I),
     &                            DY_DT(I),D2Y_D2T(I), 
     &                            H(I),DH_DT(I),D2H_D2T(I),
     &                            C(I),DC_DT(I), 
     &                            K(I),DK_DT(I),
     &             LAYER_W1(I),LAYER_CF(I),LAYER_CT(I),LAYER_LF(I),
     &             LAYER_LT(I),LAYER_B(I),LAYER_Tp(I),LAYER_Tp_1(I))
          LAYER_L(I)=LAYER_W1(I)*LIQUID_FRACTION(I)
          LAYER_I(I)=LAYER_W1(I)*(1.0D0-LIQUID_FRACTION(I))
          LAYER_K(I)=K(I)
          LAYER_C(I)=C(I)/TSCALE
        ENDDO
        DO I=SX,SN-1
          SOIL_TEMPERATURE(I)=SOIL_TEMPERATURE(SN)
          LIQUID_FRACTION(I)=LIQUID_FRACTION(SN)
        ENDDO
      ENDDO
      END SUBROUTINE STEMPERATURE_ATOMIC

      



      
      SUBROUTINE SNOW_TH_PROPERTIES(T,Y,DY_DT,D2Y_D2T, 
     &                              H,DH_DT,D2H_D2T,
     &                              C,DC_DT, 
     &                              K,DK_DT,
     &                              W,CF,CT,LF,LT,B,Tp,Tp_1)
      IMPLICIT NONE
      INCLUDE 'pwbm.soil_temperature.h'

      REAL*8 T,W,LT,LF,CF,CT,B,Tp,Tp_1
      REAL*8      Y,DY_DT,D2Y_D2T
      REAL*8      H,DH_DT,D2H_D2T
      REAL*8      C,DC_DT 
      REAL*8      K,DK_DT

      C       =CF      *TSCALE
      K       =LF

      DK_DT   =0.0D0
      DC_DT   =0.0D0

      Y       =1.0D0
      DY_DT   =0.0D0
      D2Y_D2T =0.0D0 
      H       =0.0D0
      DH_DT   =0.0D0
      D2H_D2T =0.0D0

      END SUBROUTINE SNOW_TH_PROPERTIES


      SUBROUTINE SOIL_TH_PROPERTIES(T,Y,DY_DT,D2Y_D2T, 
     &                              H,DH_DT,D2H_D2T,
     &                              C,DC_DT, 
     &                              K,DK_DT,
     &                              W,CF,CT,LF,LT,B,Tp,Tp_1)
      IMPLICIT NONE
      INCLUDE 'pwbm.soil_temperature.h'

      REAL*8 T,W,LT,LF,CF,CT,B,Tp,Tp_1

      REAL*8      Y,DY_DT,D2Y_D2T
      REAL*8      H,DH_DT,D2H_D2T
      REAL*8      C,DC_DT 
      REAL*8      K,DK_DT


      REAL*8 F1, DF1, DDF1, KWKI
      REAL*8 E2, Q1,I1,I2,Z1,TS_2,Tp_2,TS


      TS=T-XD-Tp
      IF(TS.LT.0.0D0) THEN
        TS_2=1.0D0/(TS*TS)
        Tp_2=Tp_1*Tp_1
      
        Q1=TS_2*XS_2
        E2=DEXP(-Q1)
        F1=TS*E2+Tp
        DF1=E2*(1.0D0+2.0D0*Q1)
        DDF1=DF1*(2.0D0*Q1/TS)+E2*(-4.0D0*Q1/TS)

        Z1=DLOG(F1*Tp_1)
        Y         =DEXP(B*Z1) !(F1/Tp)**B
        I1=Y*Tp/F1
        I2=I1*Tp/F1
              
        DY_DT     =B*Tp_1*dF1*I1
        D2Y_D2T   =B*Tp_1*ddF1*I1+B*(B-1.0D0)*Tp_2*dF1*dF1*I2


        C        =(CF+W*Y*(CW-CI))      *TSCALE
        DC_DT    =W*DY_DT*(CW-CI)       *TSCALE
        H        =W*Y*La                *TSCALE
        DH_DT    =W*DY_DT*La            *TSCALE      
        D2H_D2T  =W*D2Y_D2T*La          *TSCALE      

        KWKI     =DEXP(Y*W*LKWKI)
        K        =LF*KWKI
c        K = LF*KWKI*2.0d0
c        write(103,*) T, Y, KWKI, K, C
        DK_DT    =LF*KWKI*W*LKWKI*DY_DT
c        DK_DT = LF*KWKI*W*LKWKI*DY_DT*2.0d0       
      ELSE
        Y        =1.0D0
        DY_DT    =0.0D0
        D2Y_D2T  =0.0D0

        C        =CT                    *TSCALE
        DC_DT    =0.0D0
        H        =W*La                  *TSCALE
        DH_DT    =0.0D0
        D2H_D2T  =0.0D0

        K         = LT
        DK_DT     = 0.0D0
c        write(104,*) T, K, C
      ENDIF
      END SUBROUTINE SOIL_TH_PROPERTIES


      SUBROUTINE SOILSAT(T,Y,B,Tp)
      IMPLICIT NONE
      INCLUDE 'pwbm.soil_temperature.h'

      REAL*8 Y,T,B,Tp,Tp_1
      REAL*8 F1,E2,Q1,Z1,TS_2,TS
      
      Tp_1=1.0D0/Tp
      TS=T-XD-Tp
      IF(TS.LT.0.0D0) THEN
        TS_2=1.0D0/(TS*TS)
      
        Q1=TS_2*XS_2
        E2=DEXP(-Q1)
        F1=TS*E2+Tp
        Z1=DLOG(F1*Tp_1)
        Y =DEXP(B*Z1) !(F1/Tp)**B
      ELSE
        Y=1.0D0
      ENDIF
      END SUBROUTINE SOILSAT




      !Solution of a tridiagonal matrix
      SUBROUTINE TRIDIAG(M,N0,N1,A,B,C,D)
      IMPLICIT NONE
      INTEGER M,N0,N1,I,K
      REAL*8 A(N0:N1),B(N0:N1),C(N0:N1),D(N0:N1), XM

      DO K=M+1,N1
        XM=A(K)/B(K-1)
        B(K)=B(K)-XM*C(K-1)
        D(K)=D(K)-XM*D(K-1)
      ENDDO

      D(N1)=D(N1)/B(N1)
      DO I=M+1,N1
        K=N1+M-I
        D(K)=(D(K)-C(K)*D(K+1))/B(K)
      ENDDO
      END SUBROUTINE TRIDIAG


      !Subroutine to compute an upper boundary condition
      REAL*8 FUNCTION SURFACE_TEMPERTURE(TIME)
      IMPLICIT NONE
      INCLUDE 'pwbm.soil_temperature.h'
      REAL*8 TIME

      SURFACE_TEMPERTURE=SURFACE_TEMPERATURE_BEGIN+
     &   (SURFACE_TEMPERATURE_END-SURFACE_TEMPERATURE_BEGIN)*TIME
C      SURFACE_TEMPERTURE=1.0d0
      END FUNCTION SURFACE_TEMPERTURE


      





      SUBROUTINE UPDATE_PROPERTIES(X, SN,
     &                  SNOW_HEIGHT,SNOW_THCND,SNOW_HEATC,
     &                  PEAT_FRAC,LT_PEAT_DRY,LT_PEAT_SAT,
     &                  LF_PEAT_DRY,LF_PEAT_SAT,
     &                  SOIL_POROSITY,THCND_SOIL_DRY, 
     &                  THCND_SOLIDS,HEATC_SOIL_DRY,
     &                  LIQUID_FRACTION,SATURATION0,SATURATION1,
     &                  SOIL_FREEZING_PD,SOIL_CLEANESS,IWRITE_FLAG)
      IMPLICIT NONE
      INCLUDE 'pwbm.soil_temperature.h'

      REAL*8 SNOW_HEIGHT                               !Snow height
      REAL*8 SNOW_THCND, SNOW_HEATC                    !Snow thermal conductivity and heat capacity
      REAL*8 SOIL_POROSITY(2:LX+1)                      !Soil porosity
      REAL*8 THCND_SOIL_SAT,KE
      REAL*8 PEAT_FRAC(2:LX+1)
      REAL*8 LT_PEAT_DRY,LT_PEAT_SAT       ! thermal conductivity of unfrozen peat (dry & saturation)          
      REAL*8 LF_PEAT_DRY,LF_PEAT_SAT       ! thermal conductivity of frozen peat (dry & saturation)                 

      REAL*8 THCND_SOLIDS(2:LX+1)     ! Thermal conductivity of soil solids
      REAL*8 HEATC_SOIL_DRY(2:LX+1),THCND_SOIL_DRY(2:LX+1)   !Thermal conductivity of dry peat and soil
      REAL*8 X(SX:LX)
      INTEGER SN
      INTEGER I,J
      REAL*8 DX
      INTEGER SOIL_RLAYER                              !Number of layer in a rootzone
      INTEGER IWRITE_FLAG 
      
      REAL*8 LIQUID_FRACTION(SX:LX+1)
      REAL*8 SATURATION0(SX:LX+1),SATURATION1(SX:LX+1)

      REAL*8 SOIL_CLEANESS(2:LX+1), SOIL_FREEZING_PD(2:LX+1) 

      !Provided the snow height, 
      !this code analyzes the grid, and inserts a new node representing the snow surface
      IF(1.0D-5.LT.SNOW_HEIGHT) THEN
        IF(DABS(XREF(SX))<SNOW_HEIGHT) THEN
          X(SX)=-SNOW_HEIGHT
          SN=SX
        ELSE
          I=SX
          DO WHILE (I<0)
            IF((DABS(XREF(I+1))<SNOW_HEIGHT).AND.
     &               (SNOW_HEIGHT<=DABS(XREF(I))+1.0d-10)) THEN
              IF((DABS(DABS(XREF(I+1))-SNOW_HEIGHT)<
     &                 DABS(XREF(I+1)-XREF(I))*0.5D0).AND.(I<-1)) THEN
                X(I+1)=-SNOW_HEIGHT
                SN=I+1
              ELSE
                X(I)=-SNOW_HEIGHT
                SN=I
              ENDIF
              I=1
            ELSE
              I=I+1
            ENDIF
          ENDDO
        ENDIF
      ELSE
        SN=2
      ENDIF

      DO J=SX,SN-1
        X(J)=XREF(J)
      ENDDO
      DO J=SN+1,2
        X(J)=XREF(J)
      ENDDO
      DO J=2,LX
        X(J)=XREF(J)
      ENDDO
      !Computes distances between each node in the grid in order to speed up consequent computations
      DO I=SN,1
        IF (I.EQ.SN) THEN
          DX=(X(SN+1)-X(SN))/2.0D0
        ELSEIF (I.EQ.LX) THEN
          DX=(X(LX)-X(LX-1))/2.0D0
        ELSEIF((I.LT.LX).AND.(SX.LT.I)) THEN
          DX=(X(I)-X(I-1))/2.0D0+(X(I+1)-X(I))/2.0D0
        ENDIF
        DIAG(I)=DX
        IF(I<LX) DX_1(I)=1.0D0/(X(I+1)-X(I))
      ENDDO

      !Update the snow thermal properties
      DO I=SN,0
        LAYER_CF(I)=SNOW_HEATC                                !Heat capacity of snow/ice
        LAYER_LF(I)=SNOW_THCND                                !thermal conductivity of snow
        LAYER_S(I)=0.0D0
        LAYER_A(I)=0.0D0                                      
        LAYER_W0(I)=0.0D0                                      !Water content=0,No phase change in snow
        LAYER_W1(I)=0.0D0                                      !Water content=0,No phase change in snow
        LAYER_CT(I)=LAYER_CF(I)                               !Heat capacity is the same
        LAYER_B(I)=-2.0D0                                     !Parameterization of the unfrozen water content
        LAYER_Tp(I)=-0.03D0                                   !Freezing point depression
        LAYER_Tp_1(I)=1.0D0/LAYER_Tp(I)                       !Need it to speed up computations
        LAYER_LT(I)=LAYER_LF(I)                               !Thermal conductivity of frozen soil
      ENDDO

      !Update thermal properties at each soil layer
       DO I=1,LX
        LAYER_S(I+1)=1.0D0-SOIL_POROSITY(I+1)                    !Volumetric content of soil particles
        LAYER_A(I+1)=SOIL_POROSITY(I+1)*(1.0D0-SATURATION1(I+1)) !Volumetric air content
        
        LAYER_W0(I+1)=SOIL_POROSITY(I+1)*SATURATION0(I+1)        !Volumetric water content at the previous time step
        LAYER_W1(I+1)=SOIL_POROSITY(I+1)*SATURATION1(I+1)        !Volumetric water content at the previous time step
        
       ! treat organic and mineral soil separately
        if(PEAT_FRAC(I+1) .ge. 0.5d0) then
           ! thawed peat soil
           LAYER_LT(I+1) = LT_PEAT_DRY+(LT_PEAT_SAT-LT_PEAT_DRY)*
     &       SATURATION1(I+1)**2.0d0
!            LAYER_LT(I+1) = LT_PEAT_DRY+(LT_PEAT_SAT-LT_PEAT_DRY)*
!    &       SATURATION1(I+1)
           if(LAYER_LT(I+1) .lt. LT_PEAT_DRY .or. LAYER_LT(I+1) .gt.
     &       LT_PEAT_SAT) then
             print*, LAYER_LT(I+1)
             stop
           endif 
           ! frozen peat soil
           LAYER_LF(I+1) = LF_PEAT_DRY*(LF_PEAT_SAT/LF_PEAT_DRY)
     &       **SATURATION1(I+1)
           if(LAYER_LF(I+1) .lt. LF_PEAT_DRY .or. LAYER_LF(I+1) .gt. 
     &       LF_PEAT_SAT) then
             print*, LAYER_LF(I+1)
             stop
           endif
        else
           !thawed mineral soil
           THCND_SOIL_SAT=THCND_SOLIDS(I+1)**(1.0-SOIL_POROSITY(I+1))*
     &      (KW**SOIL_POROSITY(I+1))

          if(SATURATION1(I+1) .lt. 0.1) then
             KE = 0.0d0
          else 
             KE = log10(SATURATION1(I+1))+1.0
          endif

           if(KE .lt. 0.0 .or. KE .gt. 1.0d0) then
              print*, KE
              stop
           endif

           LAYER_LT(I+1)=THCND_SOIL_DRY(I+1)*(1.0-KE)+
     &       KE*THCND_SOIL_SAT  ! thermal conductivity for thawed soil

           ! frozen mineral soil
           THCND_SOIL_SAT=THCND_SOLIDS(I+1)**(1.0-SOIL_POROSITY(I+1))*
     &       (KI**SOIL_POROSITY(I+1))
           LAYER_LF(I+1)=THCND_SOIL_DRY(I+1)*(1.0-SATURATION1(I+1))+
     &       THCND_SOIL_SAT*SATURATION1(I+1)
        endif     

        ! heat capacity 
        LAYER_CT(I+1)=HEATC_SOIL_DRY(I+1)+CW*LAYER_W1(I+1)!heat capacity of thawed peat           
        LAYER_CF(I+1)=HEATC_SOIL_DRY(I+1)+CI*LAYER_W1(I+1)  ! changed by YYH

        !for the base rock (>3.25 m)
        if(I .gt. 15) then
           LAYER_LT(I+1) = 1.5d0
           LAYER_LF(I+1) = 2.0d0
           LAYER_CT(I+1) = 2.5d0
           LAYER_CF(I+1) = 2.2d0
        endif

c        if(iwrite_flag .eq. 1) then
c           print*, I,SOIL_POROSITY(I+1),SATURATION1(I+1),LAYER_CT(I+1),
c     &       LAYER_CF(I+1),LAYER_LT(I+1),LAYER_LF(I+1)
c        endif

        LAYER_B(I+1)=SOIL_CLEANESS(I+1)                        !Parameterization of the unfrozen water content
        LAYER_Tp(I+1)=SOIL_FREEZING_PD(I+1)                    !Freezing point depression

        !For testing agains the so-called classicall Stefan solution uncomment the following five lines
        !The soil is considered to be homogeneous
        LAYER_Tp_1(I+1)=1.0D0/LAYER_Tp(I+1)               !Need it to speed up computations
c        print*, I, LAYER_LT(I), LAYER_LF(I)                
      ENDDO

      !Just for the sake of pretty looking graphics
      DO I=SX,SN-1
        LAYER_CF(I)=-1.0D0                                    !Heat capacity of ice
        LAYER_LF(I)=-1.0D0                                    !Thermal conductivity of thawed soil
        LAYER_K(I)=-1.0D0
        LAYER_C(I)=-1.0D0
                                                              !!! All other parameters are fictitious and introduced to keep computing algorithm efficient
        LAYER_S(I)=0.0D0
        LAYER_A(I)=0.0D0
        LAYER_W0(I)=0.0D0                                     !Porosity=0,No phase change
        LAYER_W1(I)=0.0D0                                     !Porosity=0,No phase change
        LAYER_CT(I)=LAYER_CF(I)                               !Heat capacity of ice
        LAYER_B(I)=-0.7D0                                     !Parameterization of the unfrozen water content
        LAYER_Tp(I)=-0.03D0                                   !Freezing point depression
        LAYER_Tp_1(I)=1.0D0/LAYER_Tp(I)                       !Need it to speed up computations
        LAYER_LT(I)=LAYER_LF(I)                               !Thermal conductivity of frozen soil
      ENDDO

      END SUBROUTINE UPDATE_PROPERTIES

      SUBROUTINE CHECK_SOILFROZEN(LIQ_FRACTION,FLG)
      Implicit none
      INCLUDE 'pwbm.soil_temperature.h'
      REAL*8 LIQ_FRACTION(SX:LX+1)
      INTEGER I, FLG
     
      FLG=1

      I=2
      DO WHILE(I.LE.LX)
        IF(LIQ_FRACTION(I)>CR_LIQ_FRACTION)THEN
          FLG=0
          I=LX+1
        ELSE
          I=I+1
        ENDIF
      ENDDO
      
      END SUBROUTINE CHECK_SOILFROZEN


      SUBROUTINE INITIALIZE_SOIL(SOIL_TEMPERATURE,LIQUID_FRACTION,
     &                           SOIL_DT,GC,SN,
     &                           SOIL_TIME,SOIL_FRONT,SOIL_NODES)
      IMPLICIT NONE
      INCLUDE 'pwbm.soil_temperature.h'
      REAL*8  SOIL_TEMPERATURE(SX:LX+1)                       !Temperature
      REAL*8  LIQUID_FRACTION(SX:LX+1)                        !Liquid water fraction 
      REAL*8  SOIL_DT                                         !Automatically updated timestep 
      REAL*8  SNOW_AGE                                        !Snow age
      INTEGER GC                                              !Automatically updated variable
      INTEGER SN                                              !Number of snow layers + 1
      REAL*8  SOIL_TIME                                       !Time in soil, for testing 
      REAL*8  SOIL_FRONT                                      !Location of the freezing front closest to the surface


      INTEGER I
      REAL*8 X(SX:LX),DX,SOIL_NODES(SX:LX)                                     !Local variables


      DO I=SX,-1
        XREF(I)=-0.2D0*(DEXP(0.3D0*DBLE(ABS(I)))-1.0D0) 
      ENDDO
      XREF(0)=0.0D0
      XREF(1)  =  0.00D0
      XREF(2:lx) = soil_nodes(2:lx)

      DO I=SX,LX
        X(I)=XREF(I)
      ENDDO

      DO I=SX,LX
        IF (I.EQ.SX) THEN
          DX=(X(SX+1)-X(SX))/2.0D0
        ELSEIF (I.EQ.LX) THEN
          DX=(X(LX)-X(LX-1))/2.0D0
        ELSEIF((I.LT.LX).AND.(SX.LT.I)) THEN
          DX=(X(I)-X(I-1))/2.0D0+(X(I+1)-X(I))/2.0D0
          ENDIF
        DIAG(I)=DX
        IF(I<LX) DX_1(I)=1.0D0/(X(I+1)-X(I))
      ENDDO

      DO I=SX,LX+1
        SOIL_TEMPERATURE(I)=1.0D0
        LIQUID_FRACTION(I) =1.0D0
      ENDDO

      SOIL_TIME =0.0D0
      SOIL_FRONT=0.0D0
      
      GC=0
      SOIL_DT=MAX_DT
      SN=0
      END SUBROUTINE INITIALIZE_SOIL




