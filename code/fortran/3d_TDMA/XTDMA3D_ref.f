!              ********* ABHIJIT MUKHERJEE *******************
      SUBROUTINE XTDMA3D
!             
                IMPLICIT NONE
      INCLUDE "VAR3D.DEC"
!
!*******LOCAL VARIABLES*****************************************
!
C             *** COEFFICIENTS IN DISCRETIZATION EQUATION ***
                REAL(8) AE(ID,JD,KD),AW(ID,JD,KD),AN(ID,JD,KD)
                REAL(8) AS(ID,JD,KD),AT(ID,JD,KD),AB(ID,JD,KD)
                REAL(8) AP(ID,JD,KD),BB(ID,JD,KD),PHI(-1:ID+2,-1:JD+2,-1:KD+2)
                REAL(8) RP0, RPK, RPR,GAMMAP
                CHARACTER VAR
                COMMON /TDMA/AE,AW,AN,AS,AT,AB,AP,BB,PHI
                6 ,RP0,RPK,RPR,GAMMAP,VAR
!
                INTEGER I, J, K
      REAL(8) PX(ID), QX(ID), DEX(ID)
                REAL(8) AAX(ID), BBX(ID), CCX(ID), DDX(ID)
                REAL(8) PHICORRI(ID)
!
!*******X-LINES************************************************
!
                DO I=1, LEND
                                AAX(I)=0.
                                BBX(I)=0.
                                CCX(I)=0.
                                DDX(I)=0.
                                PX(I)=0.
                                QX(I)=0.
                                DEX(I)=0.
                                PHICORRI(I)=0.
                ENDDO
!
!              ******** X *******
                DO  K=1, NEND
                                DO  J=1, MEND
                                                DO  I=1, LEND
!
                                                AAX(I)=AAX(I)+AP(I,J,K)-AN(I,J,K)-AS(I,J,K)
     6                                                                         -AT(I,J,K)-AB(I,J,K)
                                                BBX(I)=BBX(I)+AE(I,J,K)
                                                CCX(I)=CCX(I)+AW(I,J,K)
                                DDX(I)=DDX(I)+AE(I,J,K)*PHI(I+1,J,K)+AW(I,J,K)*PHI(I-1,J,K)+
     6                                         AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)+
     6                                         AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)+
     6                                         BB(I,J,K)-AP(I,J,K)*PHI(I,J,K)
                                ENDDO
                                ENDDO
                                ENDDO
!
                                PX(1)=BBX(1)/AAX(1)
                                QX(1)=DDX(1)/AAX(1)
!
                                DO I=2, LEND
                                                PX(I)=BBX(I)/(AAX(I)-CCX(I)*PX(I-1))
                                                QX(I)=(DDX(I)+CCX(I)*QX(I-1))/(AAX(I)-CCX(I)*PX(I-1))
                                ENDDO
                                                               PHICORRI(LEND)=QX(LEND)
!
                                DO I=LEND-1, 1, -1
                                                PHICORRI(I)=PX(I)*PHICORRI(I+1)+QX(I)
                                ENDDO
!                             
                DO  K=1, NEND
                                                DO  J=1, MEND
                                                                DO  I=1, LEND
                                                PHI(I,J,K)=PHI(I,J,K)+PHICORRI(I)
                                ENDDO
                                ENDDO
                                ENDDO
!
!              &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                DO 550 K=1, NEND
                                DO 540 J=1, MEND
                                DO 520 I=1, LEND
                                                                IF (I==1) THEN
                                                                                PX(1)=AE(I,J,K)/AP(I,J,K)
                                                                ELSE
                                                                               PX(I)=AE(I,J,K)
     &                                                                              /(AP(I,J,K)-AW(I,J,K)*PX(I-1))
                                                                END IF
!                               
                                DEX(I)=BB(I,J,K)
     &                                             +AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)
     &                                  +AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)
!                                             
                                                                IF (I==1) THEN
                                                                                QX(I)=DEX(I)/AP(I,J,K)
                                                                ELSE
                                  QX(I)=(DEX(I)+AW(I,J,K)*QX(I-1))
     &                                                                              /(AP(I,J,K)-AW(I,J,K)*PX(I-1))
!
                                                                END IF
520                                         CONTINUE
!
                                                PHI(LEND,J,K)=QX(LEND)
!
 
                                                DO 530 I=LEND-1, 1, -1
                                                                PHI(I,J,K)=PX(I)*PHI(I+1,J,K)
     &                          +QX(I)
530                                         CONTINUE
!
540                         CONTINUE
!
550         CONTINUE
!
!              **************************** SUBROUTINE END ****************************   
                RETURN
                END