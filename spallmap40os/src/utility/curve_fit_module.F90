MODULE CURVE_FIT_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for curve fit
    !    linear curve fit used for obtaining earlier exfoliation time
    !
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: LINEAR_FIT

 CONTAINS

    SUBROUTINE LINEAR_FIT(X,Y,NDATA,SIG,MWT,A,B,SIGA,SIGB,CHI2,Q)
 ! linear regression
 ! sig is standard deviation
 ! y=a+b x 
 ! use mwt = 0; i.e. assume that standard deviations are not available
 ! 

    implicit none

    integer, intent(IN)                 :: ndata, mwt
    real, dimension(ndata), intent(IN)  :: x, y, sig
    real, intent(OUT)  :: a, b, siga, sigb, chi2, q

    real :: sx, sy, st2, ss, wt, SXOSS, t,SIGDAT
    integer :: i

      SX=0.
      SY=0.
      ST2=0.
      B=0.
      IF(MWT.NE.0) THEN
        SS=0.
        DO 11 I=1,NDATA
          WT=1./(SIG(I)**2)
          SS=SS+WT
          SX=SX+X(I)*WT
          SY=SY+Y(I)*WT
11      CONTINUE
      ELSE
        DO 12 I=1,NDATA
          SX=SX+X(I)
          SY=SY+Y(I)
12      CONTINUE
        SS=FLOAT(NDATA)
      ENDIF
      SXOSS=SX/SS
      IF(MWT.NE.0) THEN
        DO 13 I=1,NDATA
          T=(X(I)-SXOSS)/SIG(I)
          ST2=ST2+T*T
          B=B+T*Y(I)/SIG(I)
13      CONTINUE
      ELSE
        DO 14 I=1,NDATA
          T=X(I)-SXOSS
          ST2=ST2+T*T
          B=B+T*Y(I)
14      CONTINUE
      ENDIF
      B=B/ST2
      A=(SY-SX*B)/SS
      SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)
      SIGB=SQRT(1./ST2)
      CHI2=0.
      IF(MWT.EQ.0) THEN
        DO 15 I=1,NDATA
          CHI2=CHI2+(Y(I)-A-B*X(I))**2
15      CONTINUE
        Q=1.
        SIGDAT=SQRT(CHI2/(NDATA-2))
        SIGA=SIGA*SIGDAT
        SIGB=SIGB*SIGDAT
      ELSE
        DO 16 I=1,NDATA
          CHI2=CHI2+((Y(I)-A-B*X(I))/SIG(I))**2
16      CONTINUE
        ! Q=GAMMQ(0.5*(NDATA-2),0.5*CHI2)
      ENDIF
      RETURN
      END SUBROUTINE LINEAR_FIT

END MODULE CURVE_FIT_MODULE
