      PROGRAM hiperbolica
      IMPLICIT NONE
      REAL*8 L,H,HT,ALPHA,R,U,A,PI,X,TI,T1,TF,T,T2,M,MAT2,F_IN,F_FRON
      INTEGER I,J,K,N1,N2,P
      character*100 z
      DIMENSION A(100,100),U(100,100),TI(100),M(100,100),MAT2(100,100)
      
      PI=ACOS(-1.0)
      !Condiciones de la cuerda
      !TF=0.5   !Tiempo final !*******************************
      !L=1   !Longitud de la cuerda   ! Ejemplo 1
      !N1=8   !intervalos x
      !N2=8   !intervalos t !******************************

      TF=10   !Tiempo final !*******************************
      L=1.5   !Longitud de la cuerda   ! Ejemplo 2
      N1=80   !intervalos x
      N2=80   !intervalos t !******************************

      H=L/N1  !delta de x
      HT=TF/N2  !delta de t
      ALPHA=1  !constante 1/velocidad

      R=ALPHA*HT/H !Constante

      !Comentar la constante R para el ejemplo 1
      R=0.982395  !Ejemplo 2

      !N1=INT(L/H)
      !N2=INT(TF/HT)

      OPEN(5,FILE='HIPER.TXT')
      WRITE(5,*)N1,N2

      !LAS MATRICES SE HACEN 0
      DO I=1,N2+1
         DO J=1,N1+1
            U(I,J)=0
         END DO
      END DO
      
      DO I=2,N1
         DO J=2,N1
            A(I,J)=0
         END DO
      END DO
      
      !CONDICIONES INICIALES
      X=0
      I=1
      DO WHILE(X.LE.L)
         U(1,I) = F_IN(X) !Funci¢n de condiciones iniciale
         WRITE(*,*)"I=",I,"u=",U(1,I)
         X=X+H
         I=I+1
      END DO
      !N1=I
      !Condiciones en la frontera
      T=0
      DO I=1,N2+1
         U(I,1)=F_FRON(T)
         U(I,N1+1)=0
         WRITE(*,*) U(I,N1+1)
         T=T+HT
      END DO
       WRITE(*,*)" "
      DO I=1,N2+1
         WRITE(*,*)(U(I,J),J=1,N1+1)
      END DO

      !Se realiza el primer paso, resolver la primera matriz de cofactores
      X=2*H
      I=2
      DO WHILE(X.LE.L+H/2)
         DO J=I-1,I+1
            IF(I.NE.2.OR.J.NE.1)THEN
            IF(J.EQ.I-1.OR.J.EQ.I+1)THEN
                A(I,J)=-R**2/2
            ELSE
                A(I,J)=(2.0+R**2)
            END IF
            END IF
         END DO
         !Se calcula el t‚rmino independiente
         TI(I)=(R**2/2)*(U(1,I-1)-2*U(1,I)+U(1,I+1))+2*U(1,I)
         WRITE(5,*)TI(I)
         A(I,N1+1)=TI(I)
         X=X+H
         I=I+1
      END DO

      WRITE(5,*)" "
      WRITE(5,*)"   "
      WRITE(5,*)"MATRIZ DE COFACTORES"
      DO I=1,N1
         WRITE(5,*)(A(I,J),J=1,N1+1)
      END DO

      CALL GAUSS(A,N1+1)
      !CALL LU(A,N1,TI)   !CALCULAR LA MATRIZ POR LU

      DO I=2,N1+1
         U(2,I)=A(I,N1+1) !COMENTAR SI CALCULA LA MATRIZ POR LU
         !U(2,I)=TI(I)  !DESCOMENTAR SI DESEA RESOLVER EL PROBLEMA POR LU
      END DO

      !SIGUIENTES MATRICES DE COFACTORES
      K=2
      T=2*HT
      DO WHILE(T.LE.TF)
         DO I=2,N1
            TI(I)=0
            DO J=2,N1+1
               A(I,J)=0
            END DO
         END DO
         X=2*H
         I=2
         DO WHILE(X.LE.L+H/2)
            DO J=I-1,I+1
               IF(I.NE.2.OR.J.NE.1)THEN
               IF(J.EQ.I-1.OR.J.EQ.I+1)THEN
                  A(I,J)=-R**2/4
               ELSE
                  A(I,J)=(1.0+R**2/2)
               END IF
               END IF
            END DO
            !T‚rmino independiente
            T1=0.25*(U(K-1,I-1)-2*U(K-1,I)+U(K-1,I+1))
            T2=(R**2)*(T1+.5*(U(K,I-1)-2*U(K,I)+U(K,I+1)))
            TI(I)=T2+2*U(K,I)-U(K-1,I)
            A(I,N1+1)=TI(I)
            X=X+H
            I=I+1
         END DO
         WRITE(5,*)" "
         DO I=1,N1
            WRITE(5,*)(A(I,J),J=1,N1+1)
         END DO
         CALL GAUSS(A,N1+1)    !RESOLVER LA MATRIZ POR GAUSS
         !CALL LU(A,N1,TI)    !RESOLVER LA MATRIZ POR LU
         DO I=2,N1+1
            U(K+1,I)=A(I,N1+1) !COMENTAR SI RESUELVE LA MATRIZ POR LU
            !U(K+1,I)=TI(I)   !DESCOMENTAR SI RESUELVE LA MATRIZ POR LU
         END DO
         T=T+HT
         K=K+1
      END DO

      WRITE(5,*)" "
      DO I=1,N2+1
         WRITE(5,101)(U(I,J),J=1,N1+1)
      END DO
      !Se llena la matriz de la que se va a graficar
      X=0
      DO I=1,N2+2
         DO J=1,N1+2
            IF(I.EQ.1)THEN
               M(J,I)=X
               X=X+H
            ELSE
               M(J,I)=U(I-1,J)
            END IF
         END DO
      END DO
      
      OPEN(6,FILE="implicito.txt")
      WRITE(5,*)" "
      DO I=1,N1+1
         WRITE(5,102)(M(I,J),J=1,N2+2)
         WRITE(6,*)(M(I,J),J=1,N2+2)
      END DO
      CLOSE(5)
      CLOSE(6)
      
      !Se llena un archivo llamado fotograma para hacer el gif
      OPEN(8,FILE="fotograma.dat")
      DO I=1,N2
         X=0
         DO J=1,N1+1
            WRITE(8,*)X,U(I,J)
            X=X+H
         END DO
         WRITE(8,*)"  "
         WRITE(8,*)"  "
      END DO
      CLOSE(8)
      
      !Script de gnuplot para graficar
      OPEN(7,FILE="implicito.plt")
      WRITE(7,*)'set title "Metodo implicito" font ",12"'
      write(7,*)"set grid"
      write(7,*)"set zeroaxis"
      write(7,*)'set xrange[0:',l,']'
      !write(7,*)'set yrange[-1.2:0]'
      write(7,*)'set xlabel "x" font ",12" '
      write(7,*)'set ylabel "u(x,t)" font ",12" '
      write(7,*)'plot "implicito.txt" using 1:2 w lp'
      do i=3,N2+1
         write(7,*)'replot "implicito.txt" using 1:',i,' w lp'
      end do
      close(7)
      CALL SYSTEM('implicito.plt')
      
      !Script de gnuplot para hacer la animaci¢n de la cuerda
      OPEN(9,FILE="implicit_animation.plt")
      WRITE(9,*)'stats "fotograma.dat" name "A"'
      WRITE(9,*)'set xrange [A_min_x:',L,']'
      WRITE(9,*)'set yrange [A_min_y:A_max_y]'
      WRITE(9,*)'set border 3'
      WRITE(9,*)'set xlabel "X" font ",12"'
      WRITE(9,*)'set ylabel "U(x,t)" font ",12"'
      WRITE(9,*)'set term gif animate delay 20'
      WRITE(9,*)'set output "cuerda_implicit.gif"'
      z='do for[i=0:int(A_blocks-1)]{plot "fotograma.dat" index i w lp}'
      WRITE(9,*)z
      CLOSE(9)
      CALL SYSTEM('implicit_animation.plt')
  101 FORMAT(16F9.5)
  102 FORMAT(10F9.5)

      PAUSE
      END PROGRAM
!**************************************************
C Funci¢n de condiciones iniciales
      FUNCTION F_IN(X)
      REAL*8 F_IN,X,PI
         PI=ACOS(-1.0)

         !F_IN=-SIN(PI*X)    !EJEMPLO 1

         F_IN=0.0

         !IF(X.LE.1) THEN         !**************************
          !    F_IN=0.05*X
         !ELSE IF(X.GT.1) THEN        !EJEMPLO 2
          !    F_IN=-0.1*X+0.15
         !END IF                 !*****************************
      END FUNCTION
      
C Funci¢n de condiciones en la frontera
      FUNCTION F_FRON(T)
      REAL*8 F_FRON,T,PI
      PI=ACOS(-1.0)
      IF(T.LT.1) THEN
      F_FRON=SIN(T*PI) !Ejemplo en el que las condiciones en la frontera son diferentes de 0
      ELSE
      F_FRON=00 !Ejemplos 1 y 2 tienen condiciones en la frontera igual a 0
      END IF
      END FUNCTION
      
!! * ***********************************************************
      !RESOLVER LA MATRIZ   por gauss
      subroutine gauss(mat1,B)
      integer a,b
      real*8, dimension (100,101):: mat1,mat2,mat3
      a=B-1


       l=2
      do while(l.le.(a-1))
      if (mat1(l,l).eq.0)then
      do while (mat1(l,l).eq.0)
      do i=l+1,a
      do j=1,b
      mat1(l,j)=mat1(l,j)+mat1(i,j)
      end do
      end do
      end do
      end if

      r=1/(mat1(l,l))
      do j=1,b
      mat1(l,j)= mat1(l,j)*r
      end do


      do i=l+1,a
      f=mat1(i,l)*mat1(l,l)
      do j=1,b
      mat1(i,j)=mat1(i,j)-mat1(l,j)*f
      end do
      end do
      l=l+1
      end do

      r=1/(mat1(l,l))
      do j=2,b
      mat1(l,j)= mat1(l,j)*r
      end do

      do while(l.ge.2)
      i=l-1
      do while(i.ge.1)
      f=mat1(i,l)*mat1(l,l)


       j=b
      do while(j.ge.l)
      mat1(i,j)=mat1(i,j)-mat1(l,j)*f

      j=j-1
      end do
      i=i-1

      end do
      l=l-1

      end do

C      write(*,*) "los nuevos valores son"
      do n=2,a
      write(5,*) n,"U=",mat1(n,b)
      end do
      RETURN
      end subroutine
!!***************************************************
C Subrutina para el m‚todo LU
      SUbroutine LU (A,N,R)
        implicit none
        iNTEGER n,i,J,k
        REAL*8   s1,s2,s3,s4,s5
        Real*8,  dimension (100,100):: A,L,U
        Real*8,  dimension (100):: R,T,S
      DO I=1,n
         DO J=1,n
            L(I,J)=0
            U(I,J)=0
        END DO
        T(I)=0
        S(I)=0
      END DO

      !Se realiza las operaciones para encontrar las matrices de la solucion
       L(2,2)=A(2,2)
       U(2,2)=1
       do k=2,n-1
          U(2,k+1)=A(2,k+1)/L(2,2)
          do I=3,K
             s1=0
             do J=2,I-1
                s1=s1+L(I,J)*U(J,k+1)
             end do
          U(I,k+1)=(A(I,k+1)-s1)/L(I,I)
          end do
          L(k+1,2)=A(k+1,2)
          do I=3,k
             s2=0
             do J=2,I-1
                s2=s2+U(J,I)*L(k+1,J)
             end do
             L(k+1,I)=A(k+1,I)-s2
          end do
          U(k+1,k+1)=1
          s3=0
          do I=2,k
             s3=s3+L(k+1,I)*U(I,k+1)
          end do
          L(k+1,k+1)=A(k+1,k+1)-s3
      end do
      !se encuentra las soluciones intermedias
       T(2)= R(2)/L(2,2)
       do I=3,n
          s4=0
          do J=2,I-1
             s4=s4+(L(I,J)*T(J))
          end do
          T(I)=(R(I)-s4)/L(I,I)
       end do
       !se encuentran la soluciones finales
       S(n)=T(n)/U(n,n)
       do I=n-1,2,-1
          s5=0
          do J=I+1,n
             s5=s5+U(I,J)*S(J)
          end do
          S(I)=(T(I)-s5)/U(I,I)
      end do
      write(5,*)"Las soluciones son las siguente"
      do I=2,n
      write(5,*)"x",I,"=",S(I)
      R(I)=S(I)
      end do
      RETURN
      end subroutine

