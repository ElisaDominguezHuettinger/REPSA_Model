C----------------------------------------------------------------------
C----------------------------------------------------------------------
C   smf : Stable manifold of saddle-focus
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C       Translation table:
C
C                U[1] <--> x
C                U[2] <--> y
C                U[3] <--> z
C	             p[1] <--> alpha1
C                p[2] <--> alpha2
C                p[3] <--> alpha3
C                p[4] <--> gamma1
C                p[5] <--> gamma2
C                p(6) <--> mu1
C	             p[7] <--> mu2
C                p[8] <--> beta1
C                p[9] <--> beta2
C                p[10] <--> beta3
C                p[12] <--> TS
C                p(13) <--> tau
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
C     ---------- ----
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(NDIM),PAR(*),F(NDIM)
C
       CALL FFFF(3,U,ICP,PAR,IJAC,F,DUMMY)
       PERIOD=PAR(11)
       DO 1 I=1,NDIM
         F(I)=PERIOD*F(I)
 1     CONTINUE
C
      RETURN
      END
C
      SUBROUTINE FFFF(NDM,U,ICP,PAR,IJAC,F,DFDU)
C     ---------- ----
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(NDM),PAR(*),F(NDM),DFDU(NDM,NDM),DFDP(NDM,NDM)
C
       x= U(1)
       y= U(2)
       z= U(3)
       a1= PAR(1)
       a2= PAR(2)
       a3= PAR(3)
       g1= PAR(4)
       g2= PAR(5)
       m1= PAR(6)
       m2= PAR(7)
       b1= PAR(8)
       b2= PAR(9)
       b3= PAR(10)
       TS= PAR(12)
       tau= PAR(13)
       SL= PAR(12)-(U(1)+U(2)+U(3))
       d1=(1+PAR(4)*U(3))**(-1)
       d2=(1+PAR(5)*U(3))**(-1)
C
       F(1)= a1*x*(PAR(12)-(U(1)+U(2)+U(3)))*(1+tau*y)*(1+g1*z)**(-1)
     1       -x*(b1+PAR(6)*z)
       F(2)= a2*y*(PAR(12)-(U(1)+U(2)+U(3)))*(1+g2*z)**(-1)
     1       -y*(b2+PAR(7)*z)
       F(3)= a3*z*(PAR(12)-(U(1)+U(2)+U(3)))-b3*z
C
      IF(IJAC.EQ.0)RETURN
C
        DFDU(1,1)= -b1-m1*z-a1*(1+tau*y)*(2*x+y+z-TS)*d1
        DFDU(1,2)= -a1*x*(1+tau*(x+2*y+z-TS))*d1
        DFDU(1,3)= -m1*x+a1*x*(1+tau*y)*(g1*(x+y-TS)-1)*d1**2
C
        DFDU(2,1)= -a2*y*d2
        DFDU(2,2)= -b2-m2*z+a2*(TS-x-2*y-z)*d2
        DFDU(2,3)= -m2*y+a2*y*(g2*(x+y-TS)-1)*d2**2
C
        DFDU(3,1)= -a3*z
        DFDU(3,2)= -a3*z
        DFDU(3,3)= -b3+a3*(TS-x-y-2*z)
C
      IF(IJAC.EQ.1)RETURN
C
        DFDP(1,1)= x*(1+tau*y)*(TS-x-y-z)*d1
        DFDP(2,1)= 0.d0
        DFDP(3,1)= 0.d0
C
        DFDP(1,2)= 0.d0
        DFDP(2,2)= y*(TS-x-y-z)*d2
        DFDP(3,2)=0.d0
C
        DFDP(1,3)= 0.d0
        DFDP(2,3)= 0.d0
        DFDP(3,3)= (TS-x-y-z)*z
C 
        DFDP(1,4)= -a1*x*(1+tau*y)*(TS-x-y-z)*z*d1**2
        DFDP(2,4)= 0.d0
        DFDP(3,4)= 0.d0
C
        DFDP(1,5)= 0.d0
        DFDP(2,5)= -a2*y*(TS-x-y-z)*z*d2**2
        DFDP(3,5)= 0.d0
C
        DFDP(1,6)= -x*z
        DFDP(2,6)= 0.d0
        DFDP(3,6)= 0.d0
C
        DFDP(1,7)= 0.d0
        DFDP(2,7)= -y*z
        DFDP(3,7)= 0.d0
C
        DFDP(1,8)= -x
        DFDP(2,8)= 0.d0
        DFDP(3,8)= 0.d0
C
        DFDP(1,9)= 0.d0
        DFDP(2,9)= -y
        DFDP(3,9)= 0.d0
C
        DFDP(1,10)= 0.d0
        DFDP(2,10)= 0.d0
        DFDP(3,10)= -z
C
        DFDP(1,12)= 0.d0
        DFDP(2,12)= -a2*y*(TS-x-y-z)*z*d2**2
        DFDP(3,12)= 0.d0
C
        DFDP(1,13)= a1*x*(1+tau*y)*d1
        DFDP(2,13)= a2*y*d2
        DFDP(3,13)= a3*z
C
      RETURN
      END
C
      SUBROUTINE STPNT(NDIM,U,PAR,T)
C     ---------- -----
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FRST/ IFRST
      DIMENSION U(NDIM),PAR(*)
C
      IF(IFRST.NE.1)THEN
        IFRST=1
C       Set the parameter values
        PERIOD=0.0001
        PAR(1) = 0.0011    ! a1
        PAR(2) = 0.00137        ! a2
        PAR(3) = 0.001    ! a3
        PAR(4) = 0.01      ! g1
        PAR(5) = 0.01       ! g2
        PAR(6) = 0.01    ! m1
        PAR(7) = 0.01        ! m2
        PAR(8) = 0.29    ! b1
        PAR(9) = 0.07      ! b2
        PAR(10) = 0.1       ! b3
        PAR(12) = 146.9      ! TS
        PAR(13) = 0.0995       ! tau
        PAR(11)= PERIOD
        PAR(14) = 0.001           ! epsilon-radio
        eps=PAR(14)
        PAR(15) = 0.3d0           ! theta
C        PAR(16) = (1+g1*z)**(-1)
        PAR(17) = 3.71187E+02
        PI=4*DATAN(1.0d0)
        P1 =  0.0670547d0
        P2 = 0.0670547
        lu = 0.0373385
        vssx= 0.999752               ! vs x
        vssy= -0.00856335              ! vs y
        vssz= 0.0205421                  ! vs z
        vsx= 0.0d0                 ! vss x
        vsy= 0.992335                 ! vss y
        vsz= 0.123573                 ! vss z
        a = 17.207372980603584
        b = 23.75989193173618
        c = 5.932735087660235
        a1= PAR(1)
        a2= PAR(2)
        a3= PAR(3)
        g1= PAR(4)
        g2= PAR(5)
        m1= PAR(6)
        m2= PAR(7)
        b1= PAR(8)
        b2= PAR(9)
        b3= PAR(10)
        TS= PAR(12)
        tau= PAR(13)
        qx= 0.0790
        qy= 15.4667
        qz= 8.9326
       ENDIF
C     Singularity coordinates:
        a = 17.207372980603584
        b = 23.75989193173618
        c = 5.932735087660235
C
C       Initial conditions for 1st run:
C
        U(1) = a+eps*(cos(2*PI*PAR(15))*vsx
     1       + sin(2*PI*PAR(15))*vssx)
        U(2) = b+eps*(cos(2*PI*PAR(15))*vsy
     1       + sin(2*PI*PAR(15))*vssy)
        U(3) = c+eps*(cos(2*PI*PAR(15))*vsz
     1       + sin(2*PI*PAR(15))*vssz)
C
      RETURN
      END
C
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
C     ---------- ----
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PAR(*),ICP(*),U0(NDIM),U1(NDIM),FB(NBC)
C Local
      PARAMETER (NDM=3)
C
        PI=4*DATAN(1.0d0)
        P1 = 0.0670547d0
        P2 = 0.0670547
        lu = 0.0373385
        vssx= 0.999752               ! vs x
        vssy= -0.00856335              ! vs y
        vssz= 0.0205421                  ! vs z
        vsx= 0.0d0                 ! vss x
        vsy= 0.992335                 ! vss y
        vsz= 0.123573                 ! vss z
        eps=PAR(14)
       a = 17.207372980603584
       b = 23.75989193173618
       c = 5.932735087660235
        qx= 0.0790
        qy= 15.4667
        qz= 8.9326
C
       FB(1)=U1(1)-(a+eps*(cos(2*PI*PAR(15))*vsx
     1      +sin(2*PI*PAR(15))*vssx))
       FB(2)=U1(2)-(b+eps*(cos(2*PI*PAR(15))*vsy
     1      +sin(2*PI*PAR(15))*vssy))
       FB(3)=U1(3)-(c+eps*(cos(2*PI*PAR(15))*vsz
     1      +sin(2*PI*PAR(15))*vssz))
       FB(4)=(qx-U0(1))**2+(qy-U0(2))**2+(qz-U0(3))**2
C
C
C
      RETURN
      END
C
      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
C     ---------- ----
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(NDIM),UOLD(NDIM),UDOT(NDIM),UPOLD(NDIM)
      DIMENSION FI(NINT),ICP(*),PAR(*)
C Local
      PARAMETER (NDM=3)
      DIMENSION F(NDM),F0(NDM),DFDU(NDM,NDM)
C
      CALL FFFF(NDM,U   ,ICP,PAR,1,F ,DFDU)
      CALL FFFF(NDM,UOLD,ICP,PAR,0,F0,DUMMY)
C
      FI(1)=PAR(11)*SQRT(F(1)*F(1)+F(2)*F(2)+F(3)*F(3)) - PAR(16)
C
      RETURN
      END
C
      SUBROUTINE FOPT
      RETURN
      END
C 
      SUBROUTINE PVLS
      RETURN 
      END 
