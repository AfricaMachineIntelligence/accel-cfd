!	********* ABHIJIT MUKHERJEE *******************
      SUBROUTINE YTDMA3D
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
      REAL(8) PY(JD), QY(JD), DEY(JD)
	REAL(8) AAY(JD), BBY(JD), CCY(JD), DDY(JD)
	REAL(8) PHICORRJ(JD)
!
! ******Y-LINES********************************************************
!
	DO J=1, MEND
		AAY(J)=0.
		BBY(J)=0.
		CCY(J)=0.
		DDY(J)=0.
		PY(J)=0.
		QY(J)=0.
		DEY(J)=0.
		PHICORRJ(J)=0.
	ENDDO
!
	!	******** Y *******
	DO  K=1, NEND
		DO  J=1, MEND
			DO  I=1, LEND
!
			AAY(J)=AAY(J)+AP(I,J,K)-AE(I,J,K)-AW(I,J,K)
     6					-AT(I,J,K)-AB(I,J,K)
			BBY(J)=BBY(J)+AN(I,J,K)
			CCY(J)=CCY(J)+AS(I,J,K)
		DDY(J)=DDY(J)+AE(I,J,K)*PHI(I+1,J,K)+AW(I,J,K)*PHI(I-1,J,K)+
     6			AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)+
     6			AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)+
     6			BB(I,J,K)-AP(I,J,K)*PHI(I,J,K)
		ENDDO
		ENDDO
		ENDDO
!
		PY(1)=BBY(1)/AAY(1)
		QY(1)=DDY(1)/AAY(1)
!
		DO J=2, MEND
			PY(J)=BBY(J)/(AAY(J)-CCY(J)*PY(J-1))
			QY(J)=(DDY(J)+CCY(J)*QY(J-1))/(AAY(J)-CCY(J)*PY(J-1))
		ENDDO
			 	PHICORRJ(MEND)=QY(MEND)
!
		DO J=MEND-1, 1, -1
			PHICORRJ(J)=PY(J)*PHICORRJ(J+1)+QY(J)
		ENDDO
!
	DO  K=1, NEND
			DO  J=1, MEND
				DO  I=1, LEND
			PHI(I,J,K)=PHI(I,J,K)+PHICORRJ(J)
 		ENDDO
		ENDDO
		ENDDO
!
!	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	
	DO 515 K=1, NEND
		DO 510 I=1, LEND
			DO 500 J=1, MEND
				IF (J==1)THEN
					PY(J)=AN(I,J,K)/AP(I,J,K)
				ELSE
					PY(J)=AN(I,J,K)
     &					/(AP(I,J,K)-AS(I,J,K)*PY(J-1))
				END IF
!                        	
				DEY(J)=BB(I,J,K)
     &				+AE(I,J,K)*PHI(I+1,J,K)+AW(I,J,K)*PHI(I-1,J,K)
     &					+AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)
!	
				IF (J==1) THEN
					QY(J)= DEY(J)/AP(I,J,K)
				ELSE
					QY(J)=(DEY(J)+AS(I,J,K)*QY(J-1))
     &                      /(AP(I,J,K)-AS(I,J,K)*PY(J-1))
				END IF
500     		CONTINUE
!
			PHI(I,MEND,K)=QY(MEND)
!
			DO 505 J=MEND-1, 1, -1
				PHI(I,J,K)=PY(J)*PHI(I,J+1,K)
     &	  	        	   	   +QY(J)
505			CONTINUE
!
510     	CONTINUE
!
515	CONTINUE
!
!	**************************** SUBROUTINE END ****************************	
	RETURN
	END

