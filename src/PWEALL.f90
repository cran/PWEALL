
SUBROUTINE XPWE(NT,NR,T,R,TCHANGE,OUTR)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwe_'::XPWE
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR
   REAL(8),INTENT(IN)::T(NT),R(NR),TCHANGE(NR)
   REAL(8),INTENT(OUT)::OUTR(NT,5)

   INTEGER::I,J,J1,TINDEX(NR),AI
   REAL(8)::TPLUS(NR),TMAX,RT(NR)

   TMAX=DMAX1(MAXVAL(T),MAXVAL(TCHANGE))+1.0
   TPLUS=0.0
   TPLUS(NR)=TMAX
   IF (NR>1) THEN
      TPLUS(1:(NR-1))=TCHANGE(2:NR)
   END IF

   DO I=1,NR,1
      TINDEX(I)=I
   END DO

   RT=R*(TPLUS-TCHANGE)
   OUTR=0.0
   OUTR(:,5)=1.0
   DO I=1,NT,1
      J=COUNT(TCHANGE<=T(I))
      IF (J>=1) THEN
         J1=J-1
         OUTR(I,1)=R(J)
         IF (J==1) THEN
            OUTR(I,2)=R(1)*T(I)
         ELSE
            OUTR(I,2)=SUM(RT(1:J1))+R(J)*(T(I)-TCHANGE(J))
         END IF
         OUTR(I,5)=DEXP(-OUTR(I,2))
         OUTR(I,4)=1.0-OUTR(I,5)
         OUTR(I,3)=OUTR(I,1)*OUTR(I,5)
      END IF
   END DO
   RETURN
END SUBROUTINE XPWE

SUBROUTINE XQPWE(NP,NR,P,R,TCHANGE,OUTR)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xqpwe_' ::XQPWE
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NP,NR
   REAL(8),INTENT(IN)::P(NP),R(NR),TCHANGE(NR)
   REAL(8),INTENT(OUT)::OUTR(NP)

   INTEGER::I,J,J1
   REAL(8)::SS(0:(NR-1)),AP(NP)

   SS=0.0
   IF (NR>=2) THEN
      DO I=1,NR-1,1
         SS(I)=SS(I-1)+R(I)*(TCHANGE(I+1)-TCHANGE(I))
      END DO
   END IF

   AP=-DLOG(1.0-P)
   DO I=1,NP,1
      J=COUNT(SS<=AP(I))
      IF (J==0) THEN
         OUTR(I)=0.0
      ELSE IF (J>=1) THEN
         J1=J-1
         OUTR(I)=TCHANGE(J)+(AP(I)-SS(J1))/R(J)
      END IF
   END DO
   RETURN
END SUBROUTINE XQPWE

SUBROUTINE XPWU(NT,NU,T,U,UT,OUTR)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwu_' ::XPWU
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NU
   REAL(8),INTENT(IN)::T(NT),U(NU),UT(NU)
   REAL(8),INTENT(OUT)::OUTR(NT)

   INTEGER::I,J,J1
   REAL(8)::UT0(0:NU),RUS(0:NU)

   UT0(0)=0.0;UT0(1:NU)=UT(:)
   RUS(0)=0.0
   DO I=1,NU,1
      RUS(I)=RUS(I-1)+U(I)*(UT0(I)-UT0(I-1))
   END DO

   OUTR=0.0
   DO I=1,NT,1
     J=COUNT(UT0<T(I))
     J1=J-1
     IF (J>=1 .AND. J<=NU) THEN
        OUTR(I)=RUS(J1)+U(J)*(T(I)-UT0(J1))
     ELSE IF (J>NU) THEN
        OUTR(I)=1.0
     END IF
   END DO
   RETURN
END SUBROUTINE XPWU

SUBROUTINE XQPWU(NP,NU,P,U,UT,OUTR)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xqpwu_' ::XQPWU
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NP,NU
   REAL(8),INTENT(IN)::P(NP),U(NU),UT(NU)
   REAL(8),INTENT(OUT)::OUTR(NP)

   INTEGER::I,J,J1
   REAL(8)::UT0(0:NU),RUS(0:NU)

   UT0(0)=0.0;UT0(1:NU)=UT(:)
   RUS(0)=0.0
   DO I=1,NU,1
      RUS(I)=RUS(I-1)+U(I)*(UT0(I)-UT0(I-1))
   END DO

   DO I=1,NP,1
      J=COUNT(RUS<=P(I))
      IF (J==0) THEN
         OUTR(I)=0.0
      ELSE IF (J>=1 .AND. J<=NU) THEN
         J1=J-1
         OUTR(I)=UT0(J1)+(P(I)-RUS(J1))/U(J)
      ELSE IF (J>NU) THEN
         OUTR(I)=UT(NU)
      END IF
   END DO
   RETURN
END SUBROUTINE XQPWU

SUBROUTINE XSPF(NX,X,EPS,FX)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xspf_' ::XSPF
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NX
   REAL(8),INTENT(IN)::X(NX),EPS
   REAL(8),INTENT(OUT)::FX(NX,3)

   INTEGER::I
   REAL(8)::EX

   DO I=1,NX,1
      EX=0.0
      IF (DABS(X(I))<=EPS) THEN
         FX(I,1)=1.0-X(I)/2.0+X(I)**2/6.0-X(I)**3/24.0+X(I)**4/120.0
         FX(I,2)=0.5-X(I)/3.0+X(I)**2/8.0-X(I)**3/30.0+X(I)**4/144.0
         FX(I,3)=1.0/3.0-X(I)/4.0+X(I)**2/10.0-X(I)**3/36.0+X(I)**4/168.0
      ELSE
         EX=DEXP(-X(I))
         FX(I,1)=(1.0-EX)/X(I)
         FX(I,2)=(1.0-EX-X(I)*EX)/X(I)**2
         FX(I,3)=(2.0*(1.0-EX-X(I)*EX)-X(I)**2*EX)/X(I)**3
      END IF
   END DO
   RETURN
END SUBROUTINE XSPF

!Modified 11/15/2017

SUBROUTINE XPWEFVPLUS(NT,NR,T,R1,R2,R3,R4,R5,R6,TCHANGE,TYPE,RP2,EPS,FX)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwefvplus_' ::XPWEFVPLUS
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR,TYPE
   REAL(8),INTENT(IN)::T(NT),R1(NR),R2(NR),R3(NR),R4(NR),R5(NR),R6(NR),TCHANGE(NR),RP2,EPS
   REAL(8),INTENT(OUT)::FX(NT,3)
   REAL(8)::TEM1(NT,3),TEM2(NT,3),PIR2(NR),PIR4(NR),RONE(NR)
   IF (TYPE==1) THEN
      CALL XPWEFV4(NT,NR,T,R2,R2+R6,R3,R1+R3-R2,TCHANGE,EPS,FX)
   ELSE IF (TYPE==2) THEN
      CALL XPWEFV4TYPE2(NT,NR,T,R1+R3,R2,R3,R6,TCHANGE,EPS,FX)
   ELSE IF (TYPE==3) THEN
      PIR2=RP2*R2
      PIR4=(1.0-RP2)*R4
      RONE=1.0
      CALL XPWEFV6(NT,NR,T,R1+R3-PIR4,PIR2,R3,PIR4+R6,RONE,PIR2,TCHANGE,EPS,TEM1)
      CALL XPWEFV6(NT,NR,T,R1+R3-PIR4,PIR2,R3,PIR4+R6,PIR4,RONE,TCHANGE,EPS,TEM2) 
      FX=TEM1+TEM2
   ELSE IF (TYPE==4) THEN
      CALL XPWEFV4TYPE2(NT,NR,T,R1+R3,R2,R3,R6,TCHANGE,EPS,TEM1) 
      CALL XPWEFV4(NT,NR,T,R4,R4+R6,R3,R1+R3-R4,TCHANGE,EPS,TEM2)
      FX=RP2*TEM1+(1.0-RP2)*TEM2                                 
   ELSE IF (TYPE==5) THEN
      CALL XPWEFV4TYPE2(NT,NR,T,R1+R3,R2,R3,R6,TCHANGE,EPS,TEM1)
      CALL XPWEFV4TYPE2(NT,NR,T,R1+R3,R4,R3,R6,TCHANGE,EPS,TEM2)
      FX=RP2*TEM1+(1.0-RP2)*TEM2                                 
   ELSE
      CALL XPWEFV4(NT,NR,T,R2,R2+R6,R3,R1+R3-R2,TCHANGE,EPS,FX)
   END IF
   RETURN
END SUBROUTINE XPWEFVPLUS

SUBROUTINE XPWECXPWU(NT,NR,NU,T,TAUR,U,UT,R1,R2,R3,R4,R5,RC,TCHANGE,TYPE,RP2,EPS,DU,DUPRIME)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwecxpwu_' ::XPWECXPWU
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR,NU,TYPE
   REAL(8),INTENT(IN)::T(NT),TAUR,U(NU),UT(NU),R1(NR),R2(NR),R3(NR),R4(NR),R5(NR),RC(NR),TCHANGE(NR),RP2,EPS
   REAL(8),INTENT(OUT)::DU(NT),DUPRIME(NT)

   INTEGER::I,N1
   REAL(8)::UTMINUS(NU),ATEMP,XT(NU+1)
   REAL(8),ALLOCATABLE::TEMP(:),TM(:,:),FX4(:,:),FX2(:,:),A(:),A0(:),A1(:)

   UTMINUS=0.0
   IF (NU>1) THEN
      UTMINUS(2:NU)=UT(1:(NU-1))
   END IF
   ATEMP=SUM(U*(UT-UTMINUS))
   DU=0.0;DUPRIME=0.0
   N1=NU+1
   IF (DABS(ATEMP-1.0)<=0.00001) THEN
      XT(1)=0.0
      XT(2:N1)=UT
      ALLOCATE(TEMP(N1),TM(N1,2),FX4(N1,3),FX2(N1,3),A(N1),A0(N1),A1(N1))
      TM(:,1)=0.0
      DO I=1,NT,1
         TM(:,2)=T(I)-XT
         TEMP=MAXVAL(TM,DIM=2)
         CALL XPWEFVPLUS(N1,NR,TEMP,R1,R2,R3,R4,R5,RC,TCHANGE,TYPE,RP2,EPS,FX4)
         CALL XPWEFV2(N1,NR,TEMP,R1,R1+R3+RC,TCHANGE,EPS,FX2)
         A=FX4(:,1)+FX2(:,1);A0=TEMP*A;A1=FX4(:,2)+FX2(:,2)
         DU(I)=SUM(U*(A0(1:NU)-A0(2:N1)))-SUM(U*(A1(1:NU)-A1(2:N1)))
         DUPRIME(I)=SUM(U*(A(1:NU)-A(2:N1)))
      END DO
      DEALLOCATE(TEMP,TM,FX4,FX2,A,A0,A1)
   END IF
   RETURN
END SUBROUTINE XPWECXPWU

SUBROUTINE XPWECXPWUFORVAR(TFIX,NT,NR,NU,T,TAUR,U,UT,R1,R2,R3,R4,R5,RC,TCHANGE,TYPE,RP2,EPS,F0,F1)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwecxpwuforvar_' ::XPWECXPWUFORVAR
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR,NU,TYPE
   REAL(8),INTENT(IN)::TFIX,T(NT),TAUR,U(NU),UT(NU),R1(NR),R2(NR),R3(NR),R4(NR),R5(NR),RC(NR),TCHANGE(NR),RP2,EPS
   REAL(8),INTENT(OUT)::F0(NT),F1(NT)

   INTEGER::I,N1,NTEMP
   REAL(8)::UTMINUS(NU),ATEMP,XT(NU+1),TXMINUS(NU),TX(NU)
   REAL(8),ALLOCATABLE::TEMP(:),TEMP1(:,:),TEMP0(:,:),FX4(:,:),FX2(:,:),A0(:),A1(:),A2(:),A0T(:),A1T(:),A2T(:),A0X(:),A1X(:),A2X(:)

   UTMINUS=0.0
   IF (NU>1) THEN
      UTMINUS(2:NU)=UT(1:(NU-1))
   END IF
   TXMINUS=TFIX-UTMINUS
   TX=TFIX-UT
   ATEMP=SUM(U*(UT-UTMINUS))
   F0=0.0;F1=0.0
   N1=NU+1
   IF (DABS(ATEMP-1.0)<=0.00001) THEN
      XT(1)=0.0
      XT(2:N1)=UT
      NTEMP=N1+NT
      ALLOCATE(TEMP(NTEMP),FX4(NTEMP,3),FX2(NTEMP,3),A0(NTEMP),A1(NTEMP),A2(NTEMP))
      ALLOCATE(A0T(NT),A1T(NT),A2T(NT),A0X(N1),A1X(N1),A2X(N1))
      ALLOCATE(TEMP1(NU,4),TEMP0(NU,4))
      TEMP(1:N1)=TFIX-XT
      TEMP((N1+1):NTEMP)=T
      !CALL XPWEFV4(NTEMP,NR,TEMP,R2,R2+RC,R3,R1+R3-R2,TCHANGE,EPS,FX4)
      CALL XPWEFVPLUS(NTEMP,NR,TEMP,R1,R2,R3,R4,R5,RC,TCHANGE,TYPE,RP2,EPS,FX4)
      CALL XPWEFV2(NTEMP,NR,TEMP,R1,R1+R3+RC,TCHANGE,EPS,FX2)
      A0=FX4(:,1)+FX2(:,1);A1=FX4(:,2)+FX2(:,2);A2=FX4(:,3)+FX2(:,3)
      A0T=A0((N1+1):NTEMP);A1T=A1((N1+1):NTEMP);A2T=A2((N1+1):NTEMP)
      A0X=A0(1:N1);A1X=A1(1:N1);A2X=A2(1:N1)

      DO I=1,NT,1
         TEMP1(:,1)=(TXMINUS*(A1T(I)-A1X(2:N1))-A2T(I)+A2X(2:N1))*U
         TEMP1(:,2)=(TXMINUS*(A1X(1:NU)-A1X(2:N1))-A2X(1:NU)+A2X(2:N1))*U
         TEMP1(:,3)=A1X(2:N1)*U*(UT-UTMINUS)
         TEMP1(:,4)=A1T(I)*U*(UT-UTMINUS)

         TEMP0(:,1)=(TXMINUS*(A0T(I)-A0X(2:N1))-A1T(I)+A1X(2:N1))*U
         TEMP0(:,2)=(TXMINUS*(A0X(1:NU)-A0X(2:N1))-A1X(1:NU)+A1X(2:N1))*U
         TEMP0(:,3)=A0X(2:N1)*U*(UT-UTMINUS)
         TEMP0(:,4)=A0T(I)*U*(UT-UTMINUS)

         F1(I)=SUM(TEMP1(:,1),TXMINUS>T(I) .AND. T(I)>=TX .AND. T(I)>0.0)
         F1(I)=F1(I)+SUM(TEMP1(:,2),T(I)>=TXMINUS .AND. TXMINUS>0.0)
         F1(I)=F1(I)+SUM(TEMP1(:,3),T(I)>TX .AND. TX>=0.0)
         F1(I)=F1(I)+SUM(TEMP1(:,4),TX>=T(I))

         F0(I)=SUM(TEMP0(:,1),TXMINUS>T(I) .AND. T(I)>=TX .AND. T(I)>0.0)
         F0(I)=F0(I)+SUM(TEMP0(:,2),T(I)>=TXMINUS .AND. TXMINUS>0.0)
         F0(I)=F0(I)+SUM(TEMP0(:,3),T(I)>TX .AND. TX>=0.0)
         F0(I)=F0(I)+SUM(TEMP0(:,4),TX>=T(I))
      END DO

      DEALLOCATE(TEMP,FX4,FX2,A0,A1,A2)
      DEALLOCATE(A0T,A1T,A2T,A0X,A1X,A2X,TEMP1,TEMP0)
   END IF
   RETURN
END SUBROUTINE XPWECXPWUFORVAR

SUBROUTINE XPWEFV2(NT,NR,T,R1,R2,TCHANGE,EPS,FX)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwefv2_' ::XPWEFV2
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR
   REAL(8),INTENT(IN)::T(NT),R1(NR),R2(NR),TCHANGE(NR),EPS
   REAL(8),INTENT(OUT)::FX(NT,3)

   INTEGER::I,N1,J
   REAL(8)::TMAX,TPLUS(NR),Y(NR),PR(NR,5),PX(NR,3),A(NR),A1,X0(NR),X1(NR),X2(NR),AX(1,3),Y1,TEMP(1)
   FX=0.0
   N1=COUNT(T>0.0)
   IF (N1>0) THEN
      TMAX=DMAX1(MAXVAL(T),MAXVAL(TCHANGE))+1.0
      TPLUS=0.0
      TPLUS(NR)=TMAX
      IF (NR>1) THEN
         TPLUS(1:(NR-1))=TCHANGE(2:NR)
      END IF
      Y=TPLUS-TCHANGE
      CALL XPWE(NR,NR,TCHANGE,R2,TCHANGE,PR)
      CALL XSPF(NR,R2*Y,EPS,PX)
      A=R1*PR(:,5)
      X0=A*Y*PX(:,1)
      X1=A*Y*(Y*PX(:,2)+TCHANGE*PX(:,1))
      X2=A*Y*(Y**2*PX(:,3)+2.0*Y*TCHANGE*PX(:,2)+TCHANGE**2*PX(:,1))
      DO I=1,NT,1
         IF (T(I)>0.0) THEN
            J=COUNT(TPLUS<=T(I))
            IF (J==0) THEN
               TEMP=R2(1)*T(I)
               A1=R1(1)*T(I)
               CALL XSPF(1,TEMP,EPS,AX)
               FX(I,1)=A1*AX(1,1)
               FX(I,2)=A1*T(I)*AX(1,2)
               FX(I,3)=A1*T(I)**2*AX(1,3)
            ELSE
               Y1=T(I)-TPLUS(J)
               TEMP=R2(J+1)*Y1
               A1=A(J+1)*Y1
               CALL XSPF(1,TEMP,EPS,AX)
               FX(I,1)=SUM(X0(1:J))+A1*AX(1,1)
               FX(I,2)=SUM(X1(1:J))+A1*(Y1*AX(1,2)+TCHANGE(J+1)*AX(1,1))
               FX(I,3)=SUM(X2(1:J))+A1*(Y1**2*AX(1,3)+2.0*Y1*TCHANGE(J+1)*AX(1,2)+TCHANGE(J+1)**2*AX(1,1))
            END IF
         END IF
      END DO
   END IF
   RETURN
END SUBROUTINE XPWEFV2

SUBROUTINE XPWEFV4(NT,NR,T, R1,R2,R3,R4,TCHANGE,EPS,FX)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwefv4_' ::XPWEFV4
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR
   REAL(8),INTENT(IN)::T(NT),R1(NR),R2(NR),R3(NR),R4(NR),TCHANGE(NR),EPS
   REAL(8),INTENT(OUT)::FX(NT,3)

   INTEGER::I,N1,J
   REAL(8)::TMAX,TPLUS(NR),Y(NR),PR2(NR,5),PR4(NR,5),PX2(NR,3),PX4(NR,3),A(NR),B(NR),AX2(1,3),AX4(1,3),Y1,TV2(1),TV4(1)
   REAL(8),DIMENSION(NR)::TEMP,SUMA,TEM,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,A01,A02,A11,A12,A21,A22
   REAL(8)::TEM1,TEMP11,TEMP21,TEMP31,TEMP41,TEMP51
   FX=0.0
   N1=COUNT(T>0.0)
   IF (N1>0) THEN
      TMAX=DMAX1(MAXVAL(T),MAXVAL(TCHANGE))+1.0
      TPLUS=0.0
      TPLUS(NR)=TMAX
      IF (NR>1) THEN
         TPLUS(1:(NR-1))=TCHANGE(2:NR)
      END IF
      Y=TPLUS-TCHANGE

      CALL XPWE(NR,NR,TCHANGE,R2,TCHANGE,PR2)  !a2=PR2(:,5)
      CALL XPWE(NR,NR,TCHANGE,R4,TCHANGE,PR4)  !a4=PR4(:,5)
      CALL XSPF(NR,R2*Y,EPS,PX2)               !b2$fxi=PX2(:,i)
      CALL XSPF(NR,R4*Y,EPS,PX4)               !b4$fxi=PX4(:,i)

      TEMP=R3*PR4(:,5)*Y*PX4(:,1)
      SUMA=0.0
      IF (NR>=2) THEN
         DO J=2,NR,1
            SUMA(J)=SUM(TEMP(1:(J-1)))
         END DO
      END IF

      TEM=DEXP(-R2*Y)*PX4(:,1)
      TEMP1=Y*(PX2(:,1)-TEM)/(R2+R4)
      TEMP2=Y**2*(PX2(:,2)-TEM)/(R2+R4)+TEMP1/(R2+R4)
      TEMP3=Y**2*PX2(:,2)+TCHANGE*Y*PX2(:,1)
      TEMP4=Y**3*PX2(:,3)+2.0*TCHANGE*Y**2*PX2(:,2)+TCHANGE**2*Y*PX2(:,1)
      TEMP5=Y**3*(PX2(:,3)-TEM)/(R2+R4)+2.0*TEMP2/(R2+R4)

      A=R1*PR2(:,5);B=R3*PR4(:,5)
      A01=A*Y*PX2(:,1)*SUMA
      A02=A*R3*PR4(:,5)*TEMP1
      A11=A*TEMP3*SUMA
      A12=A*R3*PR4(:,5)*(TEMP2+TCHANGE*TEMP1)
      A21=A*TEMP4*SUMA
      A22=A*R3*PR4(:,5)*(TEMP5+2.0*TCHANGE*TEMP2+TCHANGE**2*TEMP1)

      DO I=1,NT,1
         IF (T(I)>0.0) THEN
            J=COUNT(TPLUS<=T(I))
            IF (J==0) THEN
               TV2=R2(1)*T(I);TV4=R4(1)*T(I)
               CALL XSPF(1,TV2,EPS,AX2)
               CALL XSPF(1,TV4,EPS,AX4)
               TEM1=DEXP(-R2(1)*T(I))*AX4(1,1)
               TEMP11=T(I)*(AX2(1,1)-TEM1)/(R2(1)+R4(1))
               TEMP21=T(I)**2*(AX2(1,2)-TEM1)/(R2(1)+R4(1))+TEMP11/(R2(1)+R4(1))
               TEMP51=T(I)**3*(AX2(1,3)-TEM1)/(R2(1)+R4(1))+2.0*TEMP21/(R2(1)+R4(1))
               FX(I,1)=R1(1)*R3(1)*TEMP11
               FX(I,2)=R1(1)*R3(1)*TEMP21
               FX(I,3)=R1(1)*R3(1)*TEMP51
            ELSE
               Y1=T(I)-TPLUS(J)
               TV2=R2(J+1)*Y1;TV4=R4(J+1)*Y1
               CALL XSPF(1,TV2,EPS,AX2)
               CALL XSPF(1,TV4,EPS,AX4)
               TEM1=DEXP(-R2(J+1)*Y1)*AX4(1,1)
               TEMP11=Y1*(AX2(1,1)-TEM1)/(R2(J+1)+R4(J+1))
               TEMP21=Y1**2*(AX2(1,2)-TEM1)/(R2(J+1)+R4(J+1))+TEMP11/(R2(J+1)+R4(J+1))
               TEMP31=Y1**2*AX2(1,2)+TCHANGE(J+1)*Y1*AX2(1,1)
               TEMP41=Y1**3*AX2(1,3)+2.0*TCHANGE(J+1)*Y1**2*AX2(1,2)+TCHANGE(J+1)**2*Y1*AX2(1,1)
               TEMP51=Y1**3*(AX2(1,3)-TEM1)/(R2(J+1)+R4(J+1))+2.0*TEMP21/(R2(J+1)+R4(J+1))

               FX(I,1)=SUM(A01(1:J)+A02(1:J))+A(J+1)*Y1*AX2(1,1)*SUMA(J+1)
               FX(I,1)=FX(I,1)+A(J+1)*B(J+1)*TEMP11
               FX(I,2)=SUM(A11(1:J)+A12(1:J))+A(J+1)*TEMP31*SUMA(J+1)
               FX(I,2)=FX(I,2)+A(J+1)*B(J+1)*(TEMP21+TCHANGE(J+1)*TEMP11)
               FX(I,3)=SUM(A21(1:J)+A22(1:J))+A(J+1)*TEMP41*SUMA(J+1)
               FX(I,3)=FX(I,3)+A(J+1)*B(J+1)*(TEMP51+2.0*TCHANGE(J+1)*TEMP21+TCHANGE(J+1)**2*TEMP11)
            END IF
         END IF
      END DO
   END IF
   RETURN
END SUBROUTINE XPWEFV4

SUBROUTINE XPWEDILDPWU(NT,NR,NU,T,TAUR,U,UT,R1,R2,R3,RC,TCHANGE,EPS,DU,DUPRIME)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwedildpwu_' ::XPWEDILDPWU
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR,NU
   REAL(8),INTENT(IN)::T(NT),TAUR,U(NU),UT(NU),R1(NR),R2(NR),R3(NR),RC(NR),TCHANGE(NR),EPS
   REAL(8),INTENT(OUT)::DU(NT),DUPRIME(NT)

   INTEGER::I,N1
   REAL(8)::UTMINUS(NU),ATEMP,XT(NU+1)
   REAL(8),ALLOCATABLE::TEMP(:),TM(:,:),FX4(:,:),FX2(:,:),A(:),A0(:),A1(:)

   UTMINUS=0.0
   IF (NU>1) THEN
      UTMINUS(2:NU)=UT(1:(NU-1))
   END IF
   ATEMP=SUM(U*(UT-UTMINUS))
   DU=0.0;DUPRIME=0.0
   N1=NU+1
   IF (DABS(ATEMP-1.0)<=0.00001) THEN
      XT(1)=0.0
      XT(2:N1)=UT
      ALLOCATE(TEMP(N1),TM(N1,2),FX4(N1,3),FX2(N1,3),A(N1),A0(N1),A1(N1))
      TM(:,1)=0.0
      DO I=1,NT,1
         TM(:,2)=T(I)-XT
         TEMP=MAXVAL(TM,DIM=2)
         CALL XPWEFV4(N1,NR,TEMP,R2,R2+RC,R3,R1+R3-R2,TCHANGE,EPS,FX4)
         CALL XPWEFV2(N1,NR,TEMP,R1,R1+R3+RC,TCHANGE,EPS,FX2)
         A=FX4(:,1)+FX2(:,1);A0=TEMP*A;A1=FX4(:,2)+FX2(:,2)
         DU(I)=SUM(U*(A0(1:NU)-A0(2:N1)))-SUM(U*(A1(1:NU)-A1(2:N1)))
         DUPRIME(I)=SUM(U*(A(1:NU)-A(2:N1)))
      END DO
      DEALLOCATE(TEMP,TM,FX4,FX2,A,A0,A1)
   END IF
   RETURN
END SUBROUTINE XPWEDILDPWU


SUBROUTINE XPWEDILDPWUFORVAR(TFIX,NT,NR,NU,T,TAUR,U,UT,R1,R2,R3,RC,TCHANGE,EPS,F0,F1)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwedildpwuforvar_' ::XPWEDILDPWUFORVAR
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR,NU
   REAL(8),INTENT(IN)::TFIX,T(NT),TAUR,U(NU),UT(NU),R1(NR),R2(NR),R3(NR),RC(NR),TCHANGE(NR),EPS
   REAL(8),INTENT(OUT)::F0(NT),F1(NT)

   INTEGER::I,N1,NTEMP
   REAL(8)::UTMINUS(NU),ATEMP,XT(NU+1),TXMINUS(NU),TX(NU)
   REAL(8),ALLOCATABLE::TEMP(:),TEMP1(:,:),TEMP0(:,:),FX4(:,:),FX2(:,:),A0(:),A1(:),A2(:),A0T(:),A1T(:),A2T(:),A0X(:),A1X(:),A2X(:)

   UTMINUS=0.0
   IF (NU>1) THEN
      UTMINUS(2:NU)=UT(1:(NU-1))
   END IF
   TXMINUS=TFIX-UTMINUS
   TX=TFIX-UT
   ATEMP=SUM(U*(UT-UTMINUS))
   F0=0.0;F1=0.0
   N1=NU+1
   IF (DABS(ATEMP-1.0)<=0.00001) THEN
      XT(1)=0.0
      XT(2:N1)=UT
      NTEMP=N1+NT
      ALLOCATE(TEMP(NTEMP),FX4(NTEMP,3),FX2(NTEMP,3),A0(NTEMP),A1(NTEMP),A2(NTEMP))
      ALLOCATE(A0T(NT),A1T(NT),A2T(NT),A0X(N1),A1X(N1),A2X(N1))
      ALLOCATE(TEMP1(NU,4),TEMP0(NU,4))
      TEMP(1:N1)=TFIX-XT
      TEMP((N1+1):NTEMP)=T
      CALL XPWEFV4(NTEMP,NR,TEMP,R2,R2+RC,R3,R1+R3-R2,TCHANGE,EPS,FX4)
      CALL XPWEFV2(NTEMP,NR,TEMP,R1,R1+R3+RC,TCHANGE,EPS,FX2)
      A0=FX4(:,1)+FX2(:,1);A1=FX4(:,2)+FX2(:,2);A2=FX4(:,3)+FX2(:,3)
      A0T=A0((N1+1):NTEMP);A1T=A1((N1+1):NTEMP);A2T=A2((N1+1):NTEMP)
      A0X=A0(1:N1);A1X=A1(1:N1);A2X=A2(1:N1)

      DO I=1,NT,1
         TEMP1(:,1)=(TXMINUS*(A1T(I)-A1X(2:N1))-A2T(I)+A2X(2:N1))*U
         TEMP1(:,2)=(TXMINUS*(A1X(1:NU)-A1X(2:N1))-A2X(1:NU)+A2X(2:N1))*U
         TEMP1(:,3)=A1X(2:N1)*U*(UT-UTMINUS)
         TEMP1(:,4)=A1T(I)*U*(UT-UTMINUS)

         TEMP0(:,1)=(TXMINUS*(A0T(I)-A0X(2:N1))-A1T(I)+A1X(2:N1))*U
         TEMP0(:,2)=(TXMINUS*(A0X(1:NU)-A0X(2:N1))-A1X(1:NU)+A1X(2:N1))*U
         TEMP0(:,3)=A0X(2:N1)*U*(UT-UTMINUS)
         TEMP0(:,4)=A0T(I)*U*(UT-UTMINUS)

         F1(I)=SUM(TEMP1(:,1),TXMINUS>T(I) .AND. T(I)>=TX .AND. TX>0.0)
         F1(I)=F1(I)+SUM(TEMP1(:,2),T(I)>=TXMINUS)
         F1(I)=F1(I)+SUM(TEMP1(:,3),T(I)>TX)
         F1(I)=F1(I)+SUM(TEMP1(:,4),TX>=T(I))

         F0(I)=SUM(TEMP0(:,1),TXMINUS>T(I) .AND. T(I)>=TX .AND. TX>0.0)
         F0(I)=F0(I)+SUM(TEMP0(:,2),T(I)>=TXMINUS)
         F0(I)=F0(I)+SUM(TEMP0(:,3),T(I)>TX)
         F0(I)=F0(I)+SUM(TEMP0(:,4),TX>=T(I))
      END DO

      DEALLOCATE(TEMP,FX4,FX2,A0,A1,A2)
      DEALLOCATE(A0T,A1T,A2T,A0X,A1X,A2X,TEMP1,TEMP0)
   END IF
   RETURN
END SUBROUTINE XPWEDILDPWUFORVAR

SUBROUTINE LDYN(N,Y,YT)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ldyn_' ::LDYN
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N
   REAL(8),INTENT(IN)::Y(N)
   REAL(8),INTENT(OUT)::YT(N)

   INTEGER::I,J
   REAL(8),ALLOCATABLE,DIMENSION(:)::T

   YT=0.0
   DO I=1,N,1
      ALLOCATE(T(I))
      T=(/(REAL(J,8),J=1,I)/)
      YT(1:I)=YT(1:I)+T
      DEALLOCATE(T)
   END DO
   RETURN
END SUBROUTINE LDYN

SUBROUTINE XPWEFV4TYPE2(NT,NR,TIN,R1,R2,R3,R4,TCHANGE,EPS,FX)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwefv4type2_' ::XPWEFV4TYPE2
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR
   REAL(8),INTENT(IN)::TIN(NT),R1(NR),R2(NR),R3(NR),R4(NR),TCHANGE(NR),EPS
   REAL(8),INTENT(OUT)::FX(NT,3)

   INTEGER::N1,NR2,I,J,K,POS(NT),INDT(NT),J1,JM,JM1,SINDJ3,SINDJ2,SINDJ1,SINDJ0
   REAL(8)::TT(0:(2*NR-1))

   REAL(8)::TMAX,XX
   REAL(8),DIMENSION(0:NR,5)::PR1,PR2,PR4

   REAL(8)::TX20(1,3),TX40(1,3)

   REAL(8)::L2(1),L4(1),L24,ATEMP,ATEMP1,ATEMP2,ATEMP3
   REAL(8),DIMENSION(3)::TA,TB,TC

   REAL(8),ALLOCATABLE,DIMENSION(:,:)::FXX,XLA,XLB,ALLG,TX2,TX4
   REAL(8),ALLOCATABLE,DIMENSION(:)::T,BTEMP,BTEMP1,BTEMP2,TJ
   INTEGER,ALLOCATABLE,DIMENSION(:)::INDTS,INDJ,IND,INDS

   FX=0.0;POS=0
   INDT=(/(I,I=1,NT)/)
   WHERE (TIN>0.0)
      POS=1
   END WHERE
   N1=COUNT(POS==1)

   IF (N1>0) THEN
      ALLOCATE(FXX(N1,3),T(N1),ALLG(N1,3),INDTS(N1),IND(N1),INDJ(N1))
      INDTS=PACK(INDT,POS==1)
      IND=(/(I,I=1,N1)/)
      T=PACK(TIN,POS==1)
      FXX=0.0
      ALLG=0.0
      IF (NR==1) THEN
         ALLOCATE(TX2(N1,3),TX4(N1,3),XLA(N1,3),BTEMP(N1))
         L2=R2(1)+R4(1);L4=R1(1)-R2(1)
         L24=R1(1)+R4(1)
         CALL XSPF(N1,L2(1)*T,EPS,TX2)
         CALL XSPF(N1,L4(1)*T,EPS,TX4)
         BTEMP=DEXP(-L2(1)*T)*TX4(:,1)
         XLA(:,1)=T/L24*(TX2(:,1)-BTEMP)
         XLA(:,2)=T**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
         XLA(:,3)=T**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24
         XLA=R2(1)*R3(1)*XLA
         FXX=XLA
         DEALLOCATE(TX2,TX4,XLA,BTEMP)
      END IF

      IF (NR>1) THEN
         NR2=2*NR-1
         XX=TCHANGE(2)-TCHANGE(1)
         TMAX=DMAX1(MAXVAL(T),XX*REAL(NR2-1,8))+1.0
         TT(0:(NR-1))=TCHANGE(:)
         DO I=NR,NR2-1,1
            TT(I)=REAL(I,8)*XX
         END DO
         TT(0)=-0.00000001
         TT(NR2)=TMAX

         PR1=1.0;PR2=1.0;PR4=1.0
         CALL XPWE(NR,NR,TT(1:NR),R1,TCHANGE,PR1(1:NR,:))  !a1=PR1(:,5)
         CALL XPWE(NR,NR,TT(1:NR),R2,TCHANGE,PR2(1:NR,:))  !a2=PR2(:,5)
         CALL XPWE(NR,NR,TT(1:NR),R4,TCHANGE,PR4(1:NR,:))  !a4=PR4(:,5)

         ! G2jk
         DO J=2,NR,1
            J1=J-1;INDJ=0
            WHERE (T>=TT(J))
               INDJ=2
            ELSEWHERE (T>=TT(J1))
               INDJ=1
            END WHERE
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               DO K=1,J1,1
                  L2=R2(K)+R4(J);L4=R1(J-K)-R2(K)
                  L24=R1(J-K)+R4(J)
                  CALL XSPF(1,L4*XX,EPS,TX40)
                  ATEMP=PR1(J-K-1,5)*PR2(K,5)*PR4(J1,5)*R2(K)*R3(J-K)

                  IF (SINDJ2>0) THEN
                     ALLOCATE(INDS(SINDJ2))
		     INDS=PACK(IND,INDJ==2)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP1=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP1)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP1)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP1)+2.0*TA(2)/L24

                     TB(1)=TX20(1,1)
                     TB(2)=XX*TX20(1,2)
                     TB(3)=XX**2*TX20(1,3)
                     TB=TX40(1,1)*XX**2*TB

                     TC=TB-TA
                     ALLG(INDS,1)=TC(1)
                     ALLG(INDS,2)=TC(2)+TT(J1)*TC(1)
                     ALLG(INDS,3)=TC(3)+2.0*TT(J1)*TC(2)+TT(J1)**2*TC(1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		     DEALLOCATE(INDS)
                  END IF

                  IF (SINDJ1>0) THEN
		     ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		     INDS=PACK(IND,INDJ==1)
                     TJ=T(INDS)-TT(J1)
                     CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)
                     CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     XLA(:,1)=TX40(1,1)*XX*TJ*TX2(:,1)-XLA(:,1)
                     XLA(:,2)=TX40(1,1)*XX*TJ**2*TX2(:,2)-XLA(:,2)
                     XLA(:,3)=TX40(1,1)*XX*TJ**3*TX2(:,3)-XLA(:,3)

                     ALLG(INDS,1)=XLA(:,1)
                     ALLG(INDS,2)=XLA(:,2)+TT(J1)*XLA(:,1)
                     ALLG(INDS,3)=XLA(:,3)+2.0*TT(J1)*XLA(:,2)+TT(J1)**2*XLA(:,1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		     DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF
               END DO
            END IF
         END DO

         ! G1jk K=1,NR-1,J=K,NR-1
         DO J=1,NR-1,1
            J1=J-1;INDJ=0
            WHERE (T>=TT(J))
               INDJ=2
            ELSEWHERE (T>=TT(J1))
               INDJ=1
            END WHERE
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               DO K=1,J,1
                  L2=R2(K)+R4(J);L4=R1(J-K+1)-R2(K)
                  L24=R1(J-K+1)+R4(J)
                  ATEMP=PR1(J-K,5)*PR2(K-1,5)*PR4(J1,5)*R2(K)*R3(J-K+1)

                  IF (SINDJ2>0) THEN
		     ALLOCATE(INDS(SINDJ2))
		     INDS=PACK(IND,INDJ==2)
                     CALL XSPF(1,L4*XX,EPS,TX40)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP1=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP1)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP1)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP1)+2.0*TA(2)/L24

                     ALLG(INDS,1)=TA(1)
		     ALLG(INDS,2)=TA(2)+TT(J1)*TA(1)
		     ALLG(INDS,3)=TA(3)+2.0*TT(J1)*TA(2)+TT(J1)**2*TA(1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		     DEALLOCATE(INDS)
                  END IF

                  IF (SINDJ1>0) THEN
   	             ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		     INDS=PACK(IND,INDJ==1)
                     TJ=T(INDS)-TT(J1)
                     CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)
                     CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     ALLG(INDS,1)=XLA(:,1)
                     ALLG(INDS,2)=XLA(:,2)+TT(J1)*XLA(:,1)
                     ALLG(INDS,3)=XLA(:,3)+2.0*TT(J1)*XLA(:,2)+TT(J1)**2*XLA(:,1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		     DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF
               END DO
            END IF
         END DO


         ! G1jk K=1,NR-1,J=NR,NR+K-2 SET J=J-NR+1
        IF (NR>2) THEN
         DO J=1,NR-2,1
            JM=J+NR-1;JM1=JM-1;INDJ=0
            WHERE (T>=TT(JM+1))
               INDJ=3
            ELSEWHERE (T>=TT(JM))
               INDJ=2
            ELSEWHERE (T>=TT(JM1))
               INDJ=1
            END WHERE
            SINDJ3=COUNT(INDJ==3)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J+1)
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J) AND TT(J+1)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               DO K=J+1,NR-1,1
                  L2=R2(K)+R4(NR);L4=R1(J+NR-K)-R2(K)
                  L24=R1(J+NR-K)+R4(NR)
                  ATEMP=PR1(JM-K,5)*PR2(K-1,5)*PR4(NR-1,5)*R2(K)*R3(JM-K+1)
                  ATEMP1=DEXP(-R4(NR)*REAL(J-1,8)*XX)
                  ATEMP2=DEXP(-R4(NR)*REAL(J,8)*XX-R2(K)*XX)

                  IF (SINDJ3>0) THEN
		     ALLOCATE(INDS(SINDJ3))
		     INDS=PACK(IND,INDJ==3)
                     CALL XSPF(1,L4*XX,EPS,TX40)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP3=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP3)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP3)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP3)+2.0*TA(2)/L24

                     ALLG(INDS,1)=ATEMP1*TA(1)-ATEMP2*TA(1)+ATEMP2*TX40(1,1)*XX**2*TX20(1,1)

                     ALLG(INDS,2)=ATEMP1*(TA(2)+TT(JM1)*TA(1))
                     ALLG(INDS,2)=ALLG(INDS,2)-ATEMP2*(TA(2)+TT(JM)*TA(1))
                     ALLG(INDS,2)=ALLG(INDS,2)+ATEMP2*TX40(1,1)*XX**2*(XX*TX20(1,2)+TT(JM)*TX20(1,1))

                     ALLG(INDS,3)=ATEMP1*(TA(3)+2.0*TT(JM1)*TA(2)+TT(JM1)**2*TA(1))
                     ALLG(INDS,3)=ALLG(INDS,3)-ATEMP2*(TA(3)+2.0*TT(JM)*TA(2)+TT(JM)**2*TA(1))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX**2*(XX**2*TX20(1,3))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX**2*(2.0*TT(JM)*XX*TX20(1,2))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX**2*(TT(JM)**2*TX20(1,1))

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		     DEALLOCATE(INDS)
                  END IF

                  IF (SINDJ2>0) THEN
		     ALLOCATE(INDS(SINDJ2),TJ(SINDJ2),TX2(SINDJ2,3),TX4(SINDJ2,3),BTEMP(SINDJ2),XLA(SINDJ2,3))
		     INDS=PACK(IND,INDJ==2)
                     TJ=T(INDS)-TT(JM)
                     CALL XSPF(1,L4*XX,EPS,TX40)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP3=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP3)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP3)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP3)+2.0*TA(2)/L24

                     CALL XSPF(SINDJ2,L2(1)*TJ,EPS,TX2)
                     CALL XSPF(SINDJ2,L4(1)*TJ,EPS,TX4)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)
                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     ALLG(INDS,1)=ATEMP1*TA(1)-ATEMP2*XLA(:,1)+ATEMP2*TX40(1,1)*XX*TJ*TX2(:,1)
                     ALLG(INDS,2)=ATEMP1*(TA(2)+TT(JM1)*TA(1))
                     ALLG(INDS,2)=ALLG(INDS,2)-ATEMP2*(XLA(:,2)+TT(JM)*XLA(:,1))
                     ALLG(INDS,2)=ALLG(INDS,2)+ATEMP2*TX40(1,1)*XX*TJ*(TJ*TX2(:,2)+TT(JM)*TX2(:,1))

                     ALLG(INDS,3)=ATEMP1*(TA(3)+2.0*TT(JM1)*TA(2)+TT(JM1)**2*TA(1))
                     ALLG(INDS,3)=ALLG(INDS,3)-ATEMP2*(XLA(:,3)+2.0*TT(JM)*XLA(:,2)+TT(JM)**2*XLA(:,1))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(TJ**2*TX2(:,3))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(2.0*TT(JM)*TJ*TX2(:,2))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(TT(JM)**2*TX2(:,1))

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		     DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF

                  IF (SINDJ1>0) THEN
   		     ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		     INDS=PACK(IND,INDJ==1)
                     TJ=T(INDS)-TT(JM1)
                     CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)               !b2$fxi=PX2(:,i)
                     CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)               !b4$fxi=PX4(:,i)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     ALLG(INDS,1)=XLA(:,1)
                     ALLG(INDS,2)=XLA(:,2)+TT(JM1)*XLA(:,1)
                     ALLG(INDS,3)=XLA(:,3)+2.0*TT(JM1)*XLA(:,2)+TT(JM1)**2*XLA(:,1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ATEMP1*ALLG(INDS,:)
		     DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF
               END DO
            END IF
         END DO
        END IF

         ! G1jk K=NR,J=NR,NR+K-2 SET J=J-NR+1
         DO J=1,NR-1,1
            JM=J+NR-1;JM1=JM-1;INDJ=0
            WHERE (T>=TT(JM))
               INDJ=2
            ELSEWHERE (T>=TT(JM1))
               INDJ=1
            END WHERE
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               DO K=NR,NR,1
                  L2=R2(K)+R4(NR);L4=R1(J+NR-K)-R2(K)
                  L24=R1(J+NR-K)+R4(NR)
                  ATEMP=PR1(JM-K,5)*PR2(K-1,5)*PR4(NR-1,5)*R2(K)*R3(JM-K+1)
                  ATEMP1=DEXP(-R4(NR)*REAL(J-1,8)*XX)
                  ATEMP2=DEXP(-R4(NR)*REAL(J,8)*XX-R2(K)*XX)

                  IF (SINDJ2>0) THEN
		     ALLOCATE(INDS(SINDJ2),TJ(SINDJ2),TX2(SINDJ2,3))
		     INDS=PACK(IND,INDJ==2)
                     TJ=T(INDS)-TT(JM)
                     CALL XSPF(1,L4*XX,EPS,TX40)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP3=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP3)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP3)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP3)+2.0*TA(2)/L24

                     CALL XSPF(SINDJ2,L2(1)*TJ,EPS,TX2)

                     ALLG(INDS,1)=ATEMP1*TA(1)+ATEMP2*TX40(1,1)*XX*TJ*TX2(:,1)
                     ALLG(INDS,2)=ATEMP1*(TA(2)+TT(JM1)*TA(1))
                     ALLG(INDS,2)=ALLG(INDS,2)+ATEMP2*TX40(1,1)*XX*TJ*(TJ*TX2(:,2)+TT(JM)*TX2(:,1))

                     ALLG(INDS,3)=ATEMP1*(TA(3)+2.0*TT(JM1)*TA(2)+TT(JM1)**2*TA(1))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(TJ**2*TX2(:,3))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(2.0*TT(JM)*TJ*TX2(:,2))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(TT(JM)**2*TX2(:,1))

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		     DEALLOCATE(INDS,TJ,TX2)
                  END IF

                  IF (SINDJ1>0) THEN
   		     ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		     INDS=PACK(IND,INDJ==1)
                     TJ=T(INDS)-TT(JM1)
                     CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)               !b2$fxi=PX2(:,i)
                     CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)               !b4$fxi=PX4(:,i)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     ALLG(INDS,1)=XLA(:,1)
                     ALLG(INDS,2)=XLA(:,2)+TT(JM1)*XLA(:,1)
                     ALLG(INDS,3)=XLA(:,3)+2.0*TT(JM1)*XLA(:,2)+TT(JM1)**2*XLA(:,1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ATEMP1*ALLG(INDS,:)
		     DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF
               END DO
            END IF
         END DO


         ! G1jk K=1,NR,J=NR+K-1 SET J=J-K+1 SO J=NR
         DO K=1,NR,1
            JM=NR+K-1;JM1=JM-1;INDJ=0
            WHERE (T>=TT(JM))
               INDJ=2
            ELSEWHERE (T>=TT(JM1))
               INDJ=1
            END WHERE
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               L2=R2(K)+R4(NR);L4=R1(NR)-R2(K)
               L24=R1(NR)+R4(NR)
               ATEMP=PR1(NR-1,5)*PR2(K-1,5)*PR4(NR-1,5)*R2(K)*R3(NR)
               ATEMP1=DEXP(-R4(NR)*REAL(K-1,8)*XX)
               ATEMP2=DEXP(-R4(NR)*REAL(K,8)*XX-R2(K)*XX)
               IF (SINDJ2>0) THEN
		  ALLOCATE(INDS(SINDJ2),TJ(SINDJ2),TX2(SINDJ2,3),TX4(SINDJ2,3),BTEMP(SINDJ2),XLA(SINDJ2,3),XLB(SINDJ2,3))
		  INDS=PACK(IND,INDJ==2)
                  TJ=T(INDS)-TT(JM1)
                  CALL XSPF(SINDJ2,L4(1)*TJ,EPS,TX4)
                  CALL XSPF(SINDJ2,L2(1)*TJ,EPS,TX2)
                  BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                  XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                  XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                  XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                  TJ=T(INDS)-TT(JM)
                  CALL XSPF(SINDJ2,L4(1)*TJ,EPS,TX4)
                  CALL XSPF(SINDJ2,L2(1)*TJ,EPS,TX2)
                  BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                  XLB(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                  XLB(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLB(:,1)/L24
                  XLB(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLB(:,2)/L24

                  ALLG(INDS,1)=ATEMP1*XLA(:,1)-ATEMP2*XLB(:,1)
                  ALLG(INDS,2)=ATEMP1*(XLA(:,2)+TT(JM1)*XLA(:,1))-ATEMP2*(XLB(:,2)+TT(JM)*XLB(:,1))
                  ALLG(INDS,3)=ATEMP1*(XLA(:,3)+2.0*TT(JM1)*XLA(:,2)+TT(JM1)**2*XLA(:,1))
                  ALLG(INDS,3)=ALLG(INDS,3)-ATEMP2*(XLB(:,3)+2.0*TT(JM)*XLB(:,2)+TT(JM)**2*XLB(:,1))

                  FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
                  DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA,XLB)
               END IF

               IF (SINDJ1>0) THEN
	          ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		  INDS=PACK(IND,INDJ==1)
                  TJ=T(INDS)-TT(JM1)
                  CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)
                  CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)
                  BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                  XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                  XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                  XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                  ALLG(INDS,1)=XLA(:,1)
                  ALLG(INDS,2)=(XLA(:,2)+TT(JM1)*XLA(:,1))
                  ALLG(INDS,3)=(XLA(:,3)+2.0*TT(JM1)*XLA(:,2)+TT(JM1)**2*XLA(:,1))

                  FXX(INDS,:)=FXX(INDS,:)+ATEMP*ATEMP1*ALLG(INDS,:)
		  DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
               END IF
            END IF
         END DO
      END IF
      FX(INDTS,1)=FXX(:,1)
      FX(INDTS,2)=FXX(:,2)
      FX(INDTS,3)=FXX(:,3)
      DEALLOCATE(FXX,T,ALLG,INDTS,IND,INDJ)
   END IF
   RETURN
END SUBROUTINE XPWEFV4TYPE2

SUBROUTINE XFOURHR(NT,NR,TIN,R1,R2,R3,R4,TCHANGE,EPS,FX)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xfourhr_' ::XFOURHR
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR
   REAL(8),INTENT(IN)::TIN(NT),R1(NR),R2(NR),R3(NR),R4(NR),TCHANGE(NR),EPS
   REAL(8),INTENT(OUT)::FX(NT)


   INTEGER::J,K,N1,POS(NT),INDT(NT),SINDJ1,SINDJ2,SINDJ3,SINDJ4
   REAL(8)::TMAX,TPLUS(NR),CTEMP,B,A,A1,A2,B12,A12
   REAL(8),DIMENSION(NR,5)::PR2,PR4

   REAL(8),ALLOCATABLE,DIMENSION(:,:)::TX
   REAL(8),ALLOCATABLE,DIMENSION(:)::FXX,T,TJ,X2TEMP,X4TEMP,XX,XX1
   INTEGER,ALLOCATABLE,DIMENSION(:)::INDTS,INDJ,IND,INDS

   FX=0.0;POS=0
   INDT=(/(J,J=1,NT)/)
   WHERE (TIN>0.0)
      POS=1
   END WHERE
   N1=COUNT(POS==1)

   IF (N1>0) THEN
      ALLOCATE(FXX(N1),T(N1),INDTS(N1),IND(N1),INDJ(N1))
      INDTS=PACK(INDT,POS==1)
      IND=(/(J,J=1,N1)/)
      T=PACK(TIN,POS==1)
      FXX=0.0
 
      TMAX=DMAX1(MAXVAL(T),MAXVAL(TCHANGE))+0.1
      TPLUS=0.0
      TPLUS(NR)=TMAX
      IF (NR>1) THEN
         TPLUS(1:(NR-1))=TCHANGE(2:NR)
      END IF
      
      CALL XPWE(NR,NR,TCHANGE,R2,TCHANGE,PR2)
      CALL XPWE(NR,NR,TCHANGE,R4,TCHANGE,PR4)

      DO J=1,NR,1
         DO K=1,NR,1
            B=TPLUS(J)+TPLUS(K)
            A=TCHANGE(J)+TCHANGE(K)
            A1=TCHANGE(J)+TPLUS(K)
            A2=TPLUS(J)+TCHANGE(K)
            B12=DMAX1(A1,A2)
            A12=DMIN1(A1,A2)
            CTEMP=R1(J)*PR2(J,5)*R3(K)*PR4(K,5)
            WHERE (T<A12 .AND. T>A) 
                INDJ=1
            ELSEWHERE (T<A2 .AND. T>=A1)
                INDJ=2
            ELSEWHERE (T<A1 .AND. T>=A2)
                INDJ=3
            ELSEWHERE (T<B .AND. T>=B12)
                INDJ=4
            ELSEWHERE
                INDJ=0
            END WHERE
            SINDJ1=COUNT(INDJ==1)
            SINDJ2=COUNT(INDJ==2)
            SINDJ3=COUNT(INDJ==3)
            SINDJ4=COUNT(INDJ==4)
            IF (SINDJ1>0) THEN
               ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),X4TEMP(SINDJ1),XX(SINDJ1),XX1(SINDJ1),TX(SINDJ1,3))
               INDS=PACK(IND,INDJ==1)
               TJ=PACK(T,INDJ==1)
               X4TEMP=TJ-A
               XX=TJ-A
               XX1=(R2(J)-R4(K))*XX
               CALL XSPF(SINDJ1,XX1,EPS,TX)
               FXX(INDS)=FXX(INDS)+CTEMP*DEXP(-R4(K)*X4TEMP)*XX*TX(:,1)
               DEALLOCATE(INDS,TJ,X4TEMP,XX,XX1,TX)
            END IF


            IF (SINDJ2>0) THEN
               ALLOCATE(INDS(SINDJ2),TJ(SINDJ2),X2TEMP(SINDJ2),X4TEMP(1),XX(1),XX1(1),TX(1,3))
               INDS=PACK(IND,INDJ==2)
               TJ=PACK(T,INDJ==2)
               X2TEMP=TJ-A1
               X4TEMP=TPLUS(K)-TCHANGE(K)
               XX=TPLUS(K)-TCHANGE(K)
               XX1=(R2(J)-R4(K))*XX
               CALL XSPF(1,XX1,EPS,TX)
               FXX(INDS)=FXX(INDS)+CTEMP*DEXP(-R2(J)*X2TEMP-R4(K)*X4TEMP(1))*XX(1)*TX(1,1)
               DEALLOCATE(INDS,TJ,X2TEMP,X4TEMP,XX,XX1,TX)
            END IF

            IF (SINDJ3>0) THEN
               ALLOCATE(INDS(SINDJ3),TJ(SINDJ3),X4TEMP(SINDJ3),XX(1),XX1(1),TX(1,3))
               INDS=PACK(IND,INDJ==3)
               TJ=PACK(T,INDJ==3)
               X4TEMP=TJ-A
               XX=TPLUS(J)-TCHANGE(J)
               XX1=(R2(J)-R4(K))*XX
               CALL XSPF(1,XX1,EPS,TX)
               FXX(INDS)=FXX(INDS)+CTEMP*DEXP(-R4(K)*X4TEMP)*XX(1)*TX(1,1)
               DEALLOCATE(INDS,TJ,X4TEMP,XX,XX1,TX)
            END IF


            IF (SINDJ4>0) THEN
               ALLOCATE(INDS(SINDJ4),TJ(SINDJ4),X2TEMP(SINDJ4),X4TEMP(1),XX(SINDJ4),XX1(SINDJ4),TX(SINDJ4,3))
               INDS=PACK(IND,INDJ==4)
               TJ=PACK(T,INDJ==4)
               X2TEMP=TJ-A1
               X4TEMP=TPLUS(K)-TCHANGE(K)
               XX=B-TJ
               XX1=(R2(J)-R4(K))*XX
               CALL XSPF(SINDJ4,XX1,EPS,TX)
               FXX(INDS)=FXX(INDS)+CTEMP*DEXP(-R2(J)*X2TEMP-R4(K)*X4TEMP(1))*XX*TX(:,1)
               DEALLOCATE(INDS,TJ,X2TEMP,X4TEMP,XX,XX1,TX)
            END IF
         END DO
      END DO
      FX(INDTS)=FXX
      DEALLOCATE(FXX,T,INDTS,IND,INDJ)
   END IF
   RETURN
END SUBROUTINE XFOURHR


SUBROUTINE XRMSTH(N,Y,D,TFIX,NT,TE,EPS,RMST,VRMST,VADD)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xrmsth_' ::XRMSTH
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N,NT
   REAL(8),INTENT(IN)::Y(N),TFIX,TE(NT),EPS
   INTEGER,INTENT(IN)::D(N)
   REAL(8),INTENT(OUT)::RMST,VRMST,VADD

   INTEGER::I,KK
   REAL(8)::XLAM(NT),XD(NT),XS(NT+1),XR(NT),XK(NT)
   REAL(8)::XSTEMP,TEMP(NT),TEMP1(NT)

   DO I=1,NT,1
      XR(I)=REAL(COUNT(Y>=TE(I)),8)/REAL(N,8)
      XD(I)=REAL(COUNT(DABS(Y-TE(I))<EPS .AND. D==1),8)/REAL(N,8)
      XLAM(I)=XD(I)/XR(I)
   END DO
   XS(1)=1.0
   DO I=2,NT+1,1
      XS(I)=XS(I-1)*(1.0-XLAM(I-1))
   END DO
   KK=COUNT(TE<=TFIX)
   XSTEMP=XS(KK+1)
   TEMP=TE*XS(1:NT)*XLAM
   RMST=TFIX*XSTEMP+SUM(TEMP,TE<=TFIX)
   
   DO I=1,NT,1
      XK(I)=TE(I)*XS(I+1)+SUM(TEMP(1:I))
   END DO
   
   !This one seems to underestimate the var
   !TEMP1=(RMST-XK)*(1.0-XLAM)*XLAM/XR !
   TEMP1=(RMST-XK)*XLAM/XR
   TEMP=(RMST-XK)*TEMP1
   VRMST=SUM(TEMP,TE<=TFIX)
   VADD=SUM(TEMP1,TE<=TFIX)
   RETURN
END SUBROUTINE XRMSTH

SUBROUTINE XWLRCAL(N,NT,TE,TFIX,DD1,DD0,R1,R0,NW,WEIGHTS,EPS,XTEST,XVTEST,XLR,XLCOR)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xwlrcal_' ::XWLRCAL
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N,NT,NW
   REAL(8),INTENT(IN)::TE(NT),TFIX,DD1(NT),DD0(NT),R1(NT),R0(NT),WEIGHTS(NT,NW),EPS
   REAL(8),INTENT(OUT)::XTEST(NW),XVTEST(NW),XLR(NW),XLCOR(NW,NW)

   INTEGER::I,J
   REAL(8),ALLOCATABLE,DIMENSION(:)::R,DD,TEMP,TEMP1,TEM,TEM1
   REAL(8)::ETEMP

   ETEMP=1.0/REAL(N,8)
   ALLOCATE(R(NT),DD(NT),TEMP(NT),TEMP1(NT),TEM(NT),TEM1(NT))
   R=R1+R0;DD=DD1+DD0
   TEM=DD1-DD*R1/R
   !TEM1=DD*(R1*R0/R**2)*(R-DD)/(R-1.0/REAL(N,8))
   DO I=1,NT,1
      TEM1(I)=DD(I)*(R1(I)*R0(I)/R(I)**2)
      IF (DABS(R(I)-ETEMP)>EPS) THEN
         TEM1(I)=TEM1(I)*(R(I)-DD(I))/(R(I)-ETEMP)
      END IF
   END DO
   
   
   DO I=1,NW,1
      TEMP=TEM*WEIGHTS(:,I)
      TEMP1=TEM1*WEIGHTS(:,I)**2
      XTEST(I)=SUM(TEMP,TE<=TFIX)
      XVTEST(I)=SUM(TEMP1,TE<=TFIX)
      XLR(I)=DSQRT(REAL(N,8))*XTEST(I)/DSQRT(XVTEST(I))
   END DO
   
   XLCOR=1.0
   DO I=1,NW,1
      DO J=I,NW,1
         IF (J>I) THEN
            TEMP1=TEM1*WEIGHTS(:,I)*WEIGHTS(:,J)
            XLCOR(I,J)=SUM(TEMP1,TE<=TFIX)/DSQRT(XVTEST(I)*XVTEST(J))
            XLCOR(J,I)=XLCOR(I,J)
         END IF
      END DO
   END DO
   DEALLOCATE(R,DD,TEMP,TEMP1,TEM,TEM1)
   RETURN
END SUBROUTINE XWLRCAL

!Added 3/19/2018
SUBROUTINE XWCOXCAL(N,NT,TE,TFIX,DD1,DD0,R1,R0,NW,WEIGHTS,EPS,XTOL,MAXITER,XBETA,XVBETA,XLCOR)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xwcoxcal_' ::XWCOXCAL
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N,NT,NW,MAXITER
   REAL(8),INTENT(IN)::TE(NT),TFIX,DD1(NT),DD0(NT),R1(NT),R0(NT),WEIGHTS(NT,NW),EPS,XTOL
   REAL(8),INTENT(OUT)::XBETA(NW),XVBETA(NW),XLCOR(NW,NW)

   INTEGER::I,J,K,RI
   REAL(8),ALLOCATABLE,DIMENSION(:)::R,DD,TEMP,TEMP1,TEM,TEM1,TX
   REAL(8)::ETEMP,EP,B1,B0,XTEM,X1TEM

   ETEMP=1.0/REAL(N,8)
   ALLOCATE(R(NT),DD(NT),TEMP(NT),TEMP1(NT),TEM(NT),TEM1(NT),TX(NT))
   DD=DD1+DD0
   
   DO K=1,NW,1
      RI=0;EP=1.0;B1=0.0;B0=0.0
      DO WHILE (RI<MAXITER .AND. EP>XTOL)
         TX=DEXP(B0*WEIGHTS(:,K))
         R=TX*R1+R0
         TEM=(DD1-DD*TX*R1/R)*WEIGHTS(:,K)
         TEM1=DD*(TX*R1*R0/R**2)*WEIGHTS(:,K)**2
         XTEM=SUM(TEM,TE<=TFIX)
         X1TEM=SUM(TEM1,TE<=TFIX)
         B1=B0+XTEM/X1TEM
         RI=RI+1;EP=DABS(B1-B0)
         B0=B1
      END DO
      XBETA(K)=B1;XVBETA(K)=1.0/X1TEM
   END DO
   XLCOR=1.0
   
   DEALLOCATE(R,DD,TEMP,TEMP1,TEM,TEM1,TX)
   RETURN
END SUBROUTINE XWCOXCAL



SUBROUTINE XSURVFUNC(N,Y,D,NT,TE,EPS,NTR,TR,RR,ST)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xsurvfunc_' ::XSURVFUNC
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N,NT,NTR
   REAL(8),INTENT(IN)::Y(N),TE(NT),EPS,TR(NTR)
   INTEGER,INTENT(IN)::D(N)
   REAL(8),INTENT(OUT)::RR(NTR),ST(NTR)

   INTEGER::I,KK
   REAL(8)::XLAM(NT),XD(NT),XS(NT+1),XR(NT)
   !REAL(8)::XSTEMP,TEMP(NT),TEMP1(NT)

   DO I=1,NT,1
      XR(I)=REAL(COUNT(Y>=TE(I)),8)/REAL(N,8)
      XD(I)=REAL(COUNT(DABS(Y-TE(I))<EPS .AND. D==1),8)/REAL(N,8)
      XLAM(I)=1.0-XD(I)/XR(I)
   END DO
   XS(1)=1.0
   DO I=2,NT+1,1
      XS(I)=XS(I-1)*(1.0-XLAM(I-1))
   END DO

   
   DO I=1,NTR,1
      RR(I)=REAL(COUNT(Y>=TR(I)),8)
      KK=COUNT(TE<TR(I))
      ST(I)=XS(KK+1)
   END DO
   RR=RR/REAL(N,8)
   RETURN
END SUBROUTINE XSURVFUNC

SUBROUTINE XWLRUTIL(N,Y,D,Z,NT,TE,EPS,MFUNC)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xwlrutil_' ::XWLRUTIL
   IMPLICIT NONE
   INTEGER,PARAMETER::NTOUT=21
   INTEGER,INTENT(IN)::N,NT
   REAL(8),INTENT(IN)::Y(N),TE(NT),EPS
   INTEGER,INTENT(IN)::D(N),Z(N)
   REAL(8),INTENT(OUT)::MFUNC(NT,NTOUT)
   
   REAL(8),DIMENSION(NT)::DD1,DD0,DD,R1,R0,R,S1,S0,S,CLAM1,CLAM0,CLAM,LAM1,LAM0,LAM
   ! PETO-PETO WEIGHTS
   REAL(8),DIMENSION(NT)::PETOLAM1,PETOLAM0,PETOLAM,PETOS1,PETOS0,PETOS,TS1,TS0,TS
   INTEGER::I,XR1(NT),XR0(NT),XR(NT),XD1(NT),XD0(NT),XD(NT)

   LAM1=0.0;LAM0=0.0;LAM=0.0
   DO I=1,NT,1
      XR1(I)=COUNT(Y>=TE(I) .AND. Z==1)
      XR0(I)=COUNT(Y>=TE(I) .AND. Z==0)
      XR(I)=XR1(I)+XR0(I)
      XD1(I)=COUNT(DABS(Y-TE(I))<EPS .AND. D==1 .AND. Z==1)
      XD0(I)=COUNT(DABS(Y-TE(I))<EPS .AND. D==1 .AND. Z==0)
      XD(I)=XD1(I)+XD0(I)
      PETOLAM1(I)=REAL(XD1(I),8)/REAL(XR1(I)+1,8)
      PETOLAM0(I)=REAL(XD0(I),8)/REAL(XR0(I)+1,8)
      PETOLAM(I)=REAL(XD(I),8)/REAL(XR(I)+1,8)
      IF (XR1(I)>0) THEN
         LAM1(I)=REAL(XD1(I),8)/REAL(XR1(I),8)
      END IF
      IF (XR0(I)>0) THEN
         LAM0(I)=REAL(XD0(I),8)/REAL(XR0(I),8)
      END IF
      IF (XR(I)>0) THEN
         LAM(I)=REAL(XD(I),8)/REAL(XR(I),8)
      END IF
   END DO
   DD=REAL(XD,8)/REAL(N,8);DD1=REAL(XD1,8)/REAL(N,8);DD0=REAL(XD0,8)/REAL(N,8)
   R=REAL(XR,8)/REAL(N,8);R1=REAL(XR1,8)/REAL(N,8);R0=REAL(XR0,8)/REAL(N,8)
   PETOS(1)=1.0-PETOLAM(1);PETOS1(1)=1.0-PETOLAM1(1);PETOS0(1)=1.0-PETOLAM0(1)
   TS(1)=1.0-LAM(1);TS1(1)=1.0-LAM1(1);TS0(1)=1.0-LAM0(1)
   S(1)=1.0;S1(1)=1.0;S0(1)=1.0
   CLAM(1)=0.0;CLAM1(1)=0.0;CLAM0(1)=0.0
   DO I=2,NT,1
      S(I)=S(I-1)*(1.0-LAM(I-1))
      S1(I)=S1(I-1)*(1.0-LAM1(I-1))
      S0(I)=S0(I-1)*(1.0-LAM0(I-1))
      CLAM(I)=CLAM(I-1)+LAM(I-1)
      CLAM1(I)=CLAM1(I-1)+LAM1(I-1)
      CLAM0(I)=CLAM0(I-1)+LAM0(I-1)
      PETOS(I)=PETOS(I-1)*(1.0-PETOLAM(I))
      PETOS1(I)=PETOS1(I-1)*(1.0-PETOLAM1(I))
      PETOS0(I)=PETOS0(I-1)*(1.0-PETOLAM0(I))
      TS(I)=TS(I-1)*(1.0-LAM(I))
      TS1(I)=TS1(I-1)*(1.0-LAM1(I))
      TS0(I)=TS0(I-1)*(1.0-LAM0(I))
   END DO
   !REAL(8),DIMENSION(NT)::DD1,DD0,DD,R1,R0,R,S1,S0,S,CLAM1,CLAM0,CLAM,LAM1,LAM0,LAM
   MFUNC(:,1)=DD1;MFUNC(:,2)=DD0;MFUNC(:,3)=DD;
   MFUNC(:,4)=R1;MFUNC(:,5)=R0;MFUNC(:,6)=R;
   MFUNC(:,7)=S1;MFUNC(:,8)=S0;MFUNC(:,9)=S;
   MFUNC(:,10)=CLAM1;MFUNC(:,11)=CLAM0;MFUNC(:,12)=CLAM;
   MFUNC(:,13)=LAM1;MFUNC(:,14)=LAM0;MFUNC(:,15)=LAM;
   MFUNC(:,16)=PETOS1;MFUNC(:,17)=PETOS0;MFUNC(:,18)=PETOS;
   MFUNC(:,19)=TS1;MFUNC(:,20)=TS0;MFUNC(:,21)=TS;
   RETURN
END SUBROUTINE XWLRUTIL

!Added 11/15/2017
SUBROUTINE XPWEFV6(NT,NR,TIN,R1,R2,R3,R4,R5,R6,TCHANGE,EPS,FX)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'xpwefv6_' ::XPWEFV6
   IMPLICIT NONE
   INTEGER,INTENT(IN)::NT,NR
   REAL(8),INTENT(IN)::TIN(NT),R1(NR),R2(NR),R3(NR),R4(NR),R5(NR),R6(NR),TCHANGE(NR),EPS
   REAL(8),INTENT(OUT)::FX(NT,3)

   INTEGER::N1,NR2,I,J,K,POS(NT),INDT(NT),J1,JM,JM1,SINDJ3,SINDJ2,SINDJ1,SINDJ0
   REAL(8)::TT(0:(2*NR-1))

   REAL(8)::TMAX,XX
   REAL(8),DIMENSION(0:NR,5)::PR1,PR2,PR4

   REAL(8)::TX20(1,3),TX40(1,3)

   REAL(8)::L2(1),L4(1),L24,ATEMP,ATEMP1,ATEMP2,ATEMP3
   REAL(8),DIMENSION(3)::TA,TB,TC

   REAL(8),ALLOCATABLE,DIMENSION(:,:)::FXX,XLA,XLB,ALLG,TX2,TX4
   REAL(8),ALLOCATABLE,DIMENSION(:)::T,BTEMP,BTEMP1,BTEMP2,TJ
   INTEGER,ALLOCATABLE,DIMENSION(:)::INDTS,INDJ,IND,INDS

   FX=0.0;POS=0
   INDT=(/(I,I=1,NT)/)
   WHERE (TIN>0.0)
      POS=1
   END WHERE
   N1=COUNT(POS==1)

   IF (N1>0) THEN
      ALLOCATE(FXX(N1,3),T(N1),ALLG(N1,3),INDTS(N1),IND(N1),INDJ(N1))
      INDTS=PACK(INDT,POS==1)
      IND=(/(I,I=1,N1)/)
      T=PACK(TIN,POS==1)
      FXX=0.0
      ALLG=0.0
      IF (NR==1) THEN
         ALLOCATE(TX2(N1,3),TX4(N1,3),XLA(N1,3),BTEMP(N1))
         L2=R2(1)+R4(1);L4=R1(1)-R2(1)
         L24=R1(1)+R4(1)
         CALL XSPF(N1,L2(1)*T,EPS,TX2)
         CALL XSPF(N1,L4(1)*T,EPS,TX4)
         BTEMP=DEXP(-L2(1)*T)*TX4(:,1)
         XLA(:,1)=T/L24*(TX2(:,1)-BTEMP)
         XLA(:,2)=T**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
         XLA(:,3)=T**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24
         XLA=R2(1)*R3(1)*XLA
         FXX=XLA
         DEALLOCATE(TX2,TX4,XLA,BTEMP)
      END IF

      IF (NR>1) THEN
         NR2=2*NR-1
         XX=TCHANGE(2)-TCHANGE(1)
         TMAX=DMAX1(MAXVAL(T),XX*REAL(NR2-1,8))+1.0
         TT(0:(NR-1))=TCHANGE(:)
         DO I=NR,NR2-1,1
            TT(I)=REAL(I,8)*XX
         END DO
         TT(0)=-0.00000001
         TT(NR2)=TMAX

         PR1=1.0;PR2=1.0;PR4=1.0
         CALL XPWE(NR,NR,TT(1:NR),R1,TCHANGE,PR1(1:NR,:))  !a1=PR1(:,5)
         CALL XPWE(NR,NR,TT(1:NR),R2,TCHANGE,PR2(1:NR,:))  !a2=PR2(:,5)
         CALL XPWE(NR,NR,TT(1:NR),R4,TCHANGE,PR4(1:NR,:))  !a4=PR4(:,5)

         ! G2jk
         DO J=2,NR,1
            J1=J-1;INDJ=0
            WHERE (T>=TT(J))
               INDJ=2
            ELSEWHERE (T>=TT(J1))
               INDJ=1
            END WHERE
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               DO K=1,J1,1
                  L2=R2(K)+R4(J);L4=R1(J-K)-R2(K)
                  L24=R1(J-K)+R4(J)
                  CALL XSPF(1,L4*XX,EPS,TX40)
                  !ATEMP=PR1(J-K-1,5)*PR2(K,5)*PR4(J1,5)*R2(K)*R3(J-K)  !CHANGED
                  ATEMP=PR1(J-K-1,5)*PR2(K,5)*PR4(J1,5)*R3(J-K)*R5(J)*R6(K)  

                  IF (SINDJ2>0) THEN
                     ALLOCATE(INDS(SINDJ2))
		                 INDS=PACK(IND,INDJ==2)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP1=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP1)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP1)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP1)+2.0*TA(2)/L24

                     TB(1)=TX20(1,1)
                     TB(2)=XX*TX20(1,2)
                     TB(3)=XX**2*TX20(1,3)
                     TB=TX40(1,1)*XX**2*TB

                     TC=TB-TA
                     ALLG(INDS,1)=TC(1)
                     ALLG(INDS,2)=TC(2)+TT(J1)*TC(1)
                     ALLG(INDS,3)=TC(3)+2.0*TT(J1)*TC(2)+TT(J1)**2*TC(1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		                 DEALLOCATE(INDS)
                  END IF

                  IF (SINDJ1>0) THEN
		                 ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		                 INDS=PACK(IND,INDJ==1)
                     TJ=T(INDS)-TT(J1)
                     CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)
                     CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     XLA(:,1)=TX40(1,1)*XX*TJ*TX2(:,1)-XLA(:,1)
                     XLA(:,2)=TX40(1,1)*XX*TJ**2*TX2(:,2)-XLA(:,2)
                     XLA(:,3)=TX40(1,1)*XX*TJ**3*TX2(:,3)-XLA(:,3)

                     ALLG(INDS,1)=XLA(:,1)
                     ALLG(INDS,2)=XLA(:,2)+TT(J1)*XLA(:,1)
                     ALLG(INDS,3)=XLA(:,3)+2.0*TT(J1)*XLA(:,2)+TT(J1)**2*XLA(:,1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		                 DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF
               END DO
            END IF
         END DO

         ! G1jk K=1,NR-1,J=K,NR-1
         DO J=1,NR-1,1
            J1=J-1;INDJ=0
            WHERE (T>=TT(J))
               INDJ=2
            ELSEWHERE (T>=TT(J1))
               INDJ=1
            END WHERE
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               DO K=1,J,1
                  L2=R2(K)+R4(J);L4=R1(J-K+1)-R2(K)
                  L24=R1(J-K+1)+R4(J)
                  !ATEMP=PR1(J-K,5)*PR2(K-1,5)*PR4(J1,5)*R2(K)*R3(J-K+1) !CHANGED
                  ATEMP=PR1(J-K,5)*PR2(K-1,5)*PR4(J1,5)*R3(J-K+1)*R5(J)*R6(K) 

                  IF (SINDJ2>0) THEN
		                 ALLOCATE(INDS(SINDJ2))
		                 INDS=PACK(IND,INDJ==2)
                     CALL XSPF(1,L4*XX,EPS,TX40)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP1=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP1)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP1)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP1)+2.0*TA(2)/L24

                     ALLG(INDS,1)=TA(1)
		                 ALLG(INDS,2)=TA(2)+TT(J1)*TA(1)
		                 ALLG(INDS,3)=TA(3)+2.0*TT(J1)*TA(2)+TT(J1)**2*TA(1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		                 DEALLOCATE(INDS)
                  END IF

                  IF (SINDJ1>0) THEN
   	                 ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		                 INDS=PACK(IND,INDJ==1)
                     TJ=T(INDS)-TT(J1)
                     CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)
                     CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     ALLG(INDS,1)=XLA(:,1)
                     ALLG(INDS,2)=XLA(:,2)+TT(J1)*XLA(:,1)
                     ALLG(INDS,3)=XLA(:,3)+2.0*TT(J1)*XLA(:,2)+TT(J1)**2*XLA(:,1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		                 DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF
               END DO
            END IF
         END DO


         ! G1jk K=1,NR-1,J=NR,NR+K-2 SET J=J-NR+1
        IF (NR>2) THEN
         DO J=1,NR-2,1
            JM=J+NR-1;JM1=JM-1;INDJ=0
            WHERE (T>=TT(JM+1))
               INDJ=3
            ELSEWHERE (T>=TT(JM))
               INDJ=2
            ELSEWHERE (T>=TT(JM1))
               INDJ=1
            END WHERE
            SINDJ3=COUNT(INDJ==3)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J+1)
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J) AND TT(J+1)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               DO K=J+1,NR-1,1
                  L2=R2(K)+R4(NR);L4=R1(J+NR-K)-R2(K)
                  L24=R1(J+NR-K)+R4(NR)
                  !ATEMP=PR1(JM-K,5)*PR2(K-1,5)*PR4(NR-1,5)*R2(K)*R3(JM-K+1) !CHANGED
                  ATEMP=PR1(JM-K,5)*PR2(K-1,5)*PR4(NR-1,5)*R3(JM-K+1)*R5(NR)*R6(K)
                  
                  ATEMP1=DEXP(-R4(NR)*REAL(J-1,8)*XX)
                  ATEMP2=DEXP(-R4(NR)*REAL(J,8)*XX-R2(K)*XX)

                  IF (SINDJ3>0) THEN
		                 ALLOCATE(INDS(SINDJ3))
		                 INDS=PACK(IND,INDJ==3)
                     CALL XSPF(1,L4*XX,EPS,TX40)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP3=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP3)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP3)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP3)+2.0*TA(2)/L24

                     ALLG(INDS,1)=ATEMP1*TA(1)-ATEMP2*TA(1)+ATEMP2*TX40(1,1)*XX**2*TX20(1,1)

                     ALLG(INDS,2)=ATEMP1*(TA(2)+TT(JM1)*TA(1))
                     ALLG(INDS,2)=ALLG(INDS,2)-ATEMP2*(TA(2)+TT(JM)*TA(1))
                     ALLG(INDS,2)=ALLG(INDS,2)+ATEMP2*TX40(1,1)*XX**2*(XX*TX20(1,2)+TT(JM)*TX20(1,1))

                     ALLG(INDS,3)=ATEMP1*(TA(3)+2.0*TT(JM1)*TA(2)+TT(JM1)**2*TA(1))
                     ALLG(INDS,3)=ALLG(INDS,3)-ATEMP2*(TA(3)+2.0*TT(JM)*TA(2)+TT(JM)**2*TA(1))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX**2*(XX**2*TX20(1,3))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX**2*(2.0*TT(JM)*XX*TX20(1,2))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX**2*(TT(JM)**2*TX20(1,1))

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		                 DEALLOCATE(INDS)
                  END IF

                  IF (SINDJ2>0) THEN
		                 ALLOCATE(INDS(SINDJ2),TJ(SINDJ2),TX2(SINDJ2,3),TX4(SINDJ2,3),BTEMP(SINDJ2),XLA(SINDJ2,3))
		                 INDS=PACK(IND,INDJ==2)
                     TJ=T(INDS)-TT(JM)
                     CALL XSPF(1,L4*XX,EPS,TX40)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP3=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP3)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP3)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP3)+2.0*TA(2)/L24

                     CALL XSPF(SINDJ2,L2(1)*TJ,EPS,TX2)
                     CALL XSPF(SINDJ2,L4(1)*TJ,EPS,TX4)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)
                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     ALLG(INDS,1)=ATEMP1*TA(1)-ATEMP2*XLA(:,1)+ATEMP2*TX40(1,1)*XX*TJ*TX2(:,1)
                     ALLG(INDS,2)=ATEMP1*(TA(2)+TT(JM1)*TA(1))
                     ALLG(INDS,2)=ALLG(INDS,2)-ATEMP2*(XLA(:,2)+TT(JM)*XLA(:,1))
                     ALLG(INDS,2)=ALLG(INDS,2)+ATEMP2*TX40(1,1)*XX*TJ*(TJ*TX2(:,2)+TT(JM)*TX2(:,1))

                     ALLG(INDS,3)=ATEMP1*(TA(3)+2.0*TT(JM1)*TA(2)+TT(JM1)**2*TA(1))
                     ALLG(INDS,3)=ALLG(INDS,3)-ATEMP2*(XLA(:,3)+2.0*TT(JM)*XLA(:,2)+TT(JM)**2*XLA(:,1))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(TJ**2*TX2(:,3))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(2.0*TT(JM)*TJ*TX2(:,2))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(TT(JM)**2*TX2(:,1))

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		                 DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF

                  IF (SINDJ1>0) THEN
   		               ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		                 INDS=PACK(IND,INDJ==1)
                     TJ=T(INDS)-TT(JM1)
                     CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)               !b2$fxi=PX2(:,i)
                     CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)               !b4$fxi=PX4(:,i)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     ALLG(INDS,1)=XLA(:,1)
                     ALLG(INDS,2)=XLA(:,2)+TT(JM1)*XLA(:,1)
                     ALLG(INDS,3)=XLA(:,3)+2.0*TT(JM1)*XLA(:,2)+TT(JM1)**2*XLA(:,1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ATEMP1*ALLG(INDS,:)
		                 DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF
               END DO
            END IF
         END DO
        END IF

         ! G1jk K=NR,J=NR,NR+K-2 SET J=J-NR+1
         DO J=1,NR-1,1
            JM=J+NR-1;JM1=JM-1;INDJ=0
            WHERE (T>=TT(JM))
               INDJ=2
            ELSEWHERE (T>=TT(JM1))
               INDJ=1
            END WHERE
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               DO K=NR,NR,1
                  L2=R2(K)+R4(NR);L4=R1(J+NR-K)-R2(K)
                  L24=R1(J+NR-K)+R4(NR)
                  !ATEMP=PR1(JM-K,5)*PR2(K-1,5)*PR4(NR-1,5)*R2(K)*R3(JM-K+1) !CHANGED
                  ATEMP=PR1(JM-K,5)*PR2(K-1,5)*PR4(NR-1,5)*R3(JM-K+1)*R5(NR)*R6(K)  
                  
                  ATEMP1=DEXP(-R4(NR)*REAL(J-1,8)*XX)
                  ATEMP2=DEXP(-R4(NR)*REAL(J,8)*XX-R2(K)*XX)

                  IF (SINDJ2>0) THEN
		                 ALLOCATE(INDS(SINDJ2),TJ(SINDJ2),TX2(SINDJ2,3))
		                 INDS=PACK(IND,INDJ==2)
                     TJ=T(INDS)-TT(JM)
                     CALL XSPF(1,L4*XX,EPS,TX40)
                     CALL XSPF(1,L2*XX,EPS,TX20)
                     ATEMP3=DEXP(-L2(1)*XX)*TX40(1,1)
                     TA(1)=XX/L24*(TX20(1,1)-ATEMP3)
                     TA(2)=XX**2/L24*(TX20(1,2)-ATEMP3)+TA(1)/L24
                     TA(3)=XX**3/L24*(TX20(1,3)-ATEMP3)+2.0*TA(2)/L24

                     CALL XSPF(SINDJ2,L2(1)*TJ,EPS,TX2)

                     ALLG(INDS,1)=ATEMP1*TA(1)+ATEMP2*TX40(1,1)*XX*TJ*TX2(:,1)
                     ALLG(INDS,2)=ATEMP1*(TA(2)+TT(JM1)*TA(1))
                     ALLG(INDS,2)=ALLG(INDS,2)+ATEMP2*TX40(1,1)*XX*TJ*(TJ*TX2(:,2)+TT(JM)*TX2(:,1))

                     ALLG(INDS,3)=ATEMP1*(TA(3)+2.0*TT(JM1)*TA(2)+TT(JM1)**2*TA(1))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(TJ**2*TX2(:,3))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(2.0*TT(JM)*TJ*TX2(:,2))
                     ALLG(INDS,3)=ALLG(INDS,3)+ATEMP2*TX40(1,1)*XX*TJ*(TT(JM)**2*TX2(:,1))

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
		                 DEALLOCATE(INDS,TJ,TX2)
                  END IF

                  IF (SINDJ1>0) THEN
   		               ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		                 INDS=PACK(IND,INDJ==1)
                     TJ=T(INDS)-TT(JM1)
                     CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)               !b2$fxi=PX2(:,i)
                     CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)               !b4$fxi=PX4(:,i)
                     BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                     XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                     XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                     XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                     ALLG(INDS,1)=XLA(:,1)
                     ALLG(INDS,2)=XLA(:,2)+TT(JM1)*XLA(:,1)
                     ALLG(INDS,3)=XLA(:,3)+2.0*TT(JM1)*XLA(:,2)+TT(JM1)**2*XLA(:,1)

                     FXX(INDS,:)=FXX(INDS,:)+ATEMP*ATEMP1*ALLG(INDS,:)
		                 DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
                  END IF
               END DO
            END IF
         END DO


         ! G1jk K=1,NR,J=NR+K-1 SET J=J-K+1 SO J=NR
         DO K=1,NR,1
            JM=NR+K-1;JM1=JM-1;INDJ=0
            WHERE (T>=TT(JM))
               INDJ=2
            ELSEWHERE (T>=TT(JM1))
               INDJ=1
            END WHERE
            SINDJ2=COUNT(INDJ==2)  !HOW MANY ELEMENTS IN T THAT ARE GREATER THAN TT(J)
            SINDJ1=COUNT(INDJ==1)  !HOW MANY ELEMENTS IN T THAT ARE BETWEEN TT(J1) AND TT(J)
            SINDJ0=COUNT(INDJ==0)  !HOW MANY ELEMENTS IN T THAT ARE LESS THAN TT(J1)
            IF (SINDJ0<N1) THEN
               L2=R2(K)+R4(NR);L4=R1(NR)-R2(K)
               L24=R1(NR)+R4(NR)
               !ATEMP=PR1(NR-1,5)*PR2(K-1,5)*PR4(NR-1,5)*R2(K)*R3(NR) !CHANGED
               ATEMP=PR1(NR-1,5)*PR2(K-1,5)*PR4(NR-1,5)*R3(NR)*R5(NR)*R6(K)
               
               ATEMP1=DEXP(-R4(NR)*REAL(K-1,8)*XX)
               ATEMP2=DEXP(-R4(NR)*REAL(K,8)*XX-R2(K)*XX)
               IF (SINDJ2>0) THEN
		              ALLOCATE(INDS(SINDJ2),TJ(SINDJ2),TX2(SINDJ2,3),TX4(SINDJ2,3),BTEMP(SINDJ2),XLA(SINDJ2,3),XLB(SINDJ2,3))
		              INDS=PACK(IND,INDJ==2)
                  TJ=T(INDS)-TT(JM1)
                  CALL XSPF(SINDJ2,L4(1)*TJ,EPS,TX4)
                  CALL XSPF(SINDJ2,L2(1)*TJ,EPS,TX2)
                  BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                  XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                  XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                  XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                  TJ=T(INDS)-TT(JM)
                  CALL XSPF(SINDJ2,L4(1)*TJ,EPS,TX4)
                  CALL XSPF(SINDJ2,L2(1)*TJ,EPS,TX2)
                  BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                  XLB(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                  XLB(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLB(:,1)/L24
                  XLB(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLB(:,2)/L24

                  ALLG(INDS,1)=ATEMP1*XLA(:,1)-ATEMP2*XLB(:,1)
                  ALLG(INDS,2)=ATEMP1*(XLA(:,2)+TT(JM1)*XLA(:,1))-ATEMP2*(XLB(:,2)+TT(JM)*XLB(:,1))
                  ALLG(INDS,3)=ATEMP1*(XLA(:,3)+2.0*TT(JM1)*XLA(:,2)+TT(JM1)**2*XLA(:,1))
                  ALLG(INDS,3)=ALLG(INDS,3)-ATEMP2*(XLB(:,3)+2.0*TT(JM)*XLB(:,2)+TT(JM)**2*XLB(:,1))

                  FXX(INDS,:)=FXX(INDS,:)+ATEMP*ALLG(INDS,:)
                  DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA,XLB)
               END IF

               IF (SINDJ1>0) THEN
	                ALLOCATE(INDS(SINDJ1),TJ(SINDJ1),TX2(SINDJ1,3),TX4(SINDJ1,3),BTEMP(SINDJ1),XLA(SINDJ1,3))
		              INDS=PACK(IND,INDJ==1)
                  TJ=T(INDS)-TT(JM1)
                  CALL XSPF(SINDJ1,L4(1)*TJ,EPS,TX4)
                  CALL XSPF(SINDJ1,L2(1)*TJ,EPS,TX2)
                  BTEMP=DEXP(-L2(1)*TJ)*TX4(:,1)

                  XLA(:,1)=TJ/L24*(TX2(:,1)-BTEMP)
                  XLA(:,2)=TJ**2/L24*(TX2(:,2)-BTEMP)+XLA(:,1)/L24
                  XLA(:,3)=TJ**3/L24*(TX2(:,3)-BTEMP)+2.0*XLA(:,2)/L24

                  ALLG(INDS,1)=XLA(:,1)
                  ALLG(INDS,2)=(XLA(:,2)+TT(JM1)*XLA(:,1))
                  ALLG(INDS,3)=(XLA(:,3)+2.0*TT(JM1)*XLA(:,2)+TT(JM1)**2*XLA(:,1))

                  FXX(INDS,:)=FXX(INDS,:)+ATEMP*ATEMP1*ALLG(INDS,:)
		              DEALLOCATE(INDS,TJ,TX2,TX4,BTEMP,XLA)
               END IF
            END IF
         END DO
      END IF
      FX(INDTS,1)=FXX(:,1)
      FX(INDTS,2)=FXX(:,2)
      FX(INDTS,3)=FXX(:,3)
      DEALLOCATE(FXX,T,ALLG,INDTS,IND,INDJ)
   END IF
   RETURN
END SUBROUTINE XPWEFV6

