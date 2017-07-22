!	********* ABHIJIT MUKHERJEE *******************
      SUBROUTINE ZTDMA3D
!	
	IMPLICIT NONE
      INCLUDE "VAR3D.DEC" 
!
!*******LOCAL VARIABLES*****************************************
!
C	*** COEFFICIENTS IN DISCRETIZATION EQUATION ***
	REAL(8) AE(ID,JD,KD),AW(ID,JD,KD),AN(ID,JD,KD)
	REAL(8) AS(ID,JD,KD),AT(ID,JD,KD),AB(ID,JD,KD)
	REAL(8) AP(ID,JD,KD),BB(ID,JD,KD),PHI(-1:ID+2,-1:JD+2,-1:KD+2)
	REAL(8) RP0, RPK, RPR,GAMMAP 
	CHARACTER VAR
	COMMON /TDMA/AE,AW,AN,AS,AT,AB,AP,BB,PHI
	6 ,RP0,RPK,RPR,GAMMAP,VAR
!
!
	INTEGER I, J, K
      REAL(8) PZ(KD), QZ(KD), DEZ(KD)
	REAL(8) AAZ(KD), BBZ(KD), CCZ(KD), DDZ(KD)
	REAL(8) PHICORRK(KD)
!
!*******Z-LINES********************************************************
!
	DO K=1, NEND
		AAZ(K)=0.
		BBZ(K)=0.
		CCZ(K)=0.
		DDZ(K)=0.
		PZ(K)=0.
		QZ(K)=0.
		DEZ(K)=0.
		PHICORRK(K)=0.
	ENDDO
!
	!	***********Z*******************
	DO  K=1, NEND
		DO  J=1, MEND
			DO  I=1, LEND
!
			AAZ(K)=AAZ(K)+AP(I,J,K)-AE(I,J,K)-AW(I,J,K)
     6					-AN(I,J,K)-AS(I,J,K)
			BBZ(K)=BBZ(K)+AT(I,J,K)
			CCZ(K)=CCZ(K)+AB(I,J,K)
		DDZ(K)=DDZ(K)+AE(I,J,K)*PHI(I+1,J,K)+AW(I,J,K)*PHI(I-1,J,K)+
     6			AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)+
     6			AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)+
     6			BB(I,J,K)-AP(I,J,K)*PHI(I,J,K)
		ENDDO
		ENDDO
		ENDDO
!
		PZ(1)=BBZ(1)/AAZ(1)
		QZ(1)=DDZ(1)/AAZ(1)
!
		DO K=2, NEND
			PZ(K)=BBZ(K)/(AAZ(K)-CCZ(K)*PZ(K-1))
			QZ(K)=(DDZ(K)+CCZ(K)*QZ(K-1))/(AAZ(K)-CCZ(K)*PZ(K-1))
		ENDDO
			 	PHICORRK(NEND)=QZ(NEND)
!
		DO K=NEND-1, 1, -1
			PHICORRK(K)=PZ(K)*PHICORRK(K+1)+QZ(K)
		ENDDO
!	
	DO  K=1, NEND
			DO  J=1, MEND
				DO  I=1, LEND
			PHI(I,J,K)=PHI(I,J,K)+PHICORRK(K)
		ENDDO
		ENDDO
		ENDDO
!	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        DO 590 I=1, LEND
                DO 580 J=1, MEND
                        DO 560 K=1, NEND
				IF (K==1) THEN
					PZ(K)=AT(I,J,K)/AP(I,J,K)
				ELSE 
                          PZ(K)=AT(I,J,K)
     &                          	 /(AP(I,J,K)-AB(I,J,K)*PZ(K-1))
				ENDIF
!                                
                    DEZ(K)=BB(I,J,K)
     &				+AN(I,J,K)*PHI(I,J+1,K)+AE(I,J,K)*PHI(I+1,J,K)
     &                  +AW(I,J,K)*PHI(I-1,J,K)+AS(I,J,K)*PHI(I,J-1,K)
!                                
				IF (K==1) THEN
					QZ(K)=DEZ(K)/AP(I,J,K)
				ELSE
                           QZ(K)=(DEZ(K)+AB(I,J,K)*QZ(K-1))
     &                      	      /(AP(I,J,K)-AB(I,J,K)*PZ(K-1))
				END IF
560                     CONTINUE
!
                        PHI(I,J,NEND)=QZ(NEND)
!
                        DO 570 K=NEND-1, 1, -1
!				
            		PHI(I,J,K)=PZ(K)*PHI(I,J,K+1)
     &              			 +QZ(K)

570                     CONTINUE
!
580             CONTINUE
!
590     CONTINUE
!
!	**************************** SUBROUTINE END ****************************	
	RETURN
	END

