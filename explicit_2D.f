      program derivadasParciales
      implicitnone
!      es una ecuaciondiferencial que tiene los valores de frontera de
!      y(t,x) donde y(0,x)=0 y y(t,0)=0
      !se realizara el procedimiento utilizando la logica de las difer
!     encias finitas
      real*8:: h,ht,L,A,FinX,FinT,x,e,c,r,pendiente,Fmult,Fsum,t,y,corX
      real*8:: corY,corZ,c2,r2,max,W,PI,V,Sum_error,po
      integer::j,i,k,intervalos,or,nx,p,ny,u,s,NODX,NODY,cont
      character*100::z
      real*8,dimension(300,2000,2000)::mat
      real*8,dimension(1000,1)::error
      write(*,*)"ingrese el numero de intervalos de tiempo"
      read(*,*)intervalos


      !para poder evealuar necesitamos las ocndiciones inicales
      !suponiendo que tenemos una cuerda de longitudd
c     para llenar las condiciones iniciales de las cuales partiremos
c     prociedemos a iniciar las primeras condiciones de la matriz
      write(*,*)"Programa encargado de resolver EDP hiperbolicas"
      write(*,*)" de segundo orden, de dos variables espaciales "
      write(*,*)"con 1 nodo de vibracion"
      or=2

      write(*,*) "ingrese el ancho en x"
      read(*,*) A
      write(*,*) "ingrese el largo en y"
      read(*,*) L
      write(*,*)"ingrese el valor de los incrementos en x y y"
      read(*,*)h
      ht=h/100.0
      pendiente=(ht/h)**or
      
      NODX=3
      NODY=3
      PI=ACOS(-1.0)
      V=2000.0
      W=sqrt(V*((NODX**2)*(PI**2)/(A**2)+(NODY**2)*(PI**2)/(L**2)))
      Write(*,*)"w es ", w

      PRINT*,PENDIENTE
c     el numero de intervalos en lso que se partira la cuerda es grande
c     dependiendo del tama¤o de la misma para tener ma¤or presision

      r= A/h+1
      nx=int(r)
      r= L/h+1
      ny=int(r)

c     vamos a hacer I como el contador de T, J oara X y u para Y
      do i=1,intervalos
      do j=1,nx
      do u=1,ny
      mat(i,j,u)=0.0
      end do
      end do
      end do
      
c     empezamos nuestras condiciones en la frontera----------------------
      !condiciones para los valores de forntera de y
      do i=1,intervalos
      do j=1,nx
      mat(i,j,1)=0.0
      mat(i,j,ny)=0.0
      end do
      end do
      !condiciones para los valores de frontera de x
      do i=1,intervalos
      do u=1,ny
      mat(i,1,u)=0.0
      mat(i,nx,u)=0.0
      end do
      end do
      !condiciones para los valores de frontera para t=0
      do j=1,nx
      do u=1,ny
      mat(1,j,u)=0.0
      mat(1,j,u)=0.0
      end do
      end do
      !Condiciones iniciales
      corX=A/2
      corY=L/2
      corZ=0.003
      max=corZ
      t=0.0
      x=0.0+h
      do j=2,nx-1
      y=0.0+h
      do u=2,ny-1
      mat(1,j,u)=Fsum(x,y,t,NODX,NODY,W,A,L)
      y=y+h
      end do
      x=x+h
      end do
     
      
!!!!!!! aqui va a estar el corazon del programa el corazon de pollo
      write(*,*) pendiente
      c=2000.0*pendiente
      c2=2000.0*pendiente
!_---------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
      k=1
      t=0.0+ht
      do i=1,intervalos-1
      x=0.0+h
      do j=2,nx-1
      y=0.0+h
      do u=2,ny-1
      
      if (i.eq.1)then
      r=c*(mat(i,j+1,u)+mat(i,j-1,u))/2.0
      r2=c2*(mat(i,j,u-1)+mat(i,j,u+1))/2.0
      mat(i+1,j,u)=r+(1-c)*mat(i,j,u)
      !mat(i+1,j,u)=
      else
      mat(i+1,j,u)=Fmult(i,j,u,mat,c,c2)
      if (mat(i+1,j,u).gt.max)then
      max=mat(i+1,j,u)
      !go to 45
      end if
      end if
      po=Fsum(x,y,t,NODX,NODY,W)
      if (po.ne.0.0)then
      
      sum_error=sum_error+abs((mat(i+1,j,u)-po)/po)*100.0
      cont=cont+1
      end if
      y=y+h

      end do
      x=x+h
      end do
      Error(i,1)=sum_error/cont
      sum_error=0.0
      cont=0
      k=k+1
      t=t+ht
      end do
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------

   45  s=nx*ny
      k=1
      
      open(16,file="coor.xyz")

      t=0
      do i=1,intervalos
      write(16,*)"  ",s
      x=0.0
      do j=1,nx
      y=0.0
      do u=1,ny


      write(16,14)k,x,y,mat(i,j,u),22,k

      y=y+h
      k=k+1
      if (k.gt.s)then
      k=1
      end if
      end do
      x=x+h
      end do
      t=t+ht
      end do

      close(16)
      
      
      
      
      open(18,file="archivo.txt")
      x=0
      do j=1,nx
      y=0
      do u=1,ny
      write(18,*)x,y,mat(1,j,u)
      y=y+h
      end do
      x=x+h
      end do
      close(18)
      !----------------------------------------------------------
      !_----------------------------------------------------------
      !-----------------------------------------------------------
      open(30,file="plotCoordenadas.plt")
      write(30,*)"set terminal gif animate delay 2"
      write(30,*)"set output ""cuerda.gif"""
      write(30,*)"#set terminal x11"

      write(30,*)"set ylabel ""Y"""
      write(30,*)"set xlabel ""X"""
      write(30,*)"set zlabel ""Z"""
      write(30,*)"set xrange [0:",A,"]"
      write(30,*)"set yrange [0:",L,"]"
      write(30,*)"set zrange [",(-1.7*max),":",(1.7*max),"]"
      write(30,*)"set title ""COMPORTAMIENDO DE ONDA"""
      write(30,*)"set style data points"
      write(30,*)"set surface"
      write(30,*)"stats ""fotograma.txt"" name ""A"""


      write(30,*)"file1=""fotograma.txt"""
      z="for [i=0:int(A_blocks-1)]{splot file1 index i w p}"
      write(30,*)"do ",z

      close(30)

!-------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      open(43,file="ploteo.plt")
      write(43,*)"#set terminal x11"
      write(43,*)"set ylabel ""Y"""
      write(43,*)"set xlabel ""X"""
      write(43,*)"set zlabel ""Z"""
      write(43,*)"set xrange [0:",A,"]"
      write(43,*)"set yrange [0:",L,"]"
      write(43,*)"set zrange [",(-1.7*max),":",(1.7*max),"]"
      write(43,*)"set title ""COMPORTAMIENDO DE ONDA"""
      write(43,*)"set style data points"
      write(43,*)"set surface"
      write(43,*)"set termoption dashed"

!#Archivo gif de salida
      write(43,*)"set isosamples 51"
      write(43,*)"stats ""fotograma.txt"" name ""A"""


      write(43,*)"file1=""fotograma.txt"""
      z="for [i=0:int(A_blocks-1)]{splot file1 index i w p}"

      write(43,*)"do ",z

      close(43)
      
      
      open(85,file="fotograma.txt")
      do i=1,intervalos

      x=0.0
      do j=1,nx
      y=0.0
      do u=1,ny


      write(85,*)x,y,mat(i,j,u)

      y=y+h
      end do
      x=x+h
      end do
      write(85,*)"                                                 "
      write(85,*)"                                                 "
      
      end do
      close(85)
      
      
      !------------------------------------------------------------
      open(150,file="archivoPy.txt")
      do i=1,intervalos

      x=0.0
      do j=1,nx
      y=0.0
      do u=1,ny


      write(150,*)x,y,mat(i,j,u)

      y=y+h
      end do
      x=x+h
      end do

      end do
      close(150)


      !-----------------------------------------
      open(12,file="Error2.txt")
      t=ht
      do i=1,intervalos-1
      
      write(12,*)t,Error(i,1)
      t=t+ht
      end do
      close(12)
      
      call system('plotCoordenadas.plt')
      
  14     format(3x,I4,2x,'O',3x,f10.5,2x,f10.5,2x,f10.5,4x,I3,3x,I4)
      pause
      end program
      
      
      function Fmult(i,j,u,mat,c,c2)
          real*8,dimension(300,2000,2000):: mat
          real*8::c,r,c2,Fmult
          integer::i,j,u
          r=2.0*(1-c-c2)*mat(i,j,u)+c*(mat(i,j+1,u)+mat(i,j-1,u))
          Fmult=r+c2*(mat(i,j,u+1)+mat(i,j,u-1))-mat(i-1,j,u)
          end
          
          
      function Fsum(x,y,t,NODX,NODY,w,A,L)
      real*8::x,y,t,pi,w,A,L,Fsum
      integer::NODX,NODY
      pi=acos(-1.0)
      Fsum=0.003*(sin(NODX*pi*x/A)*sin(NODY*pi*y/L))*cos(w*t)
      END

      
      subroutine Piramide(mat,nx,ny,A,L,corX,corY,corZ,h,matr)
      real*8,dimension(300,2000,2000)::mat
      integer::i,j,u,nx,ny
      real*8::a1,b1,c1,d1,a2,a3,a4,b2,b3,b4,c2,c3,c4,d2,d3,d4,Y3,h
      real*8::corX,corY,corZ,A,L,x,y
      a1=0.0
      b1=corZ*A
      c1=corY*A*(-1)
      d1=(a1*corX+b1*corY+c1*corZ)*(-1)
      a2=corZ*L*(-1)
      b2=0.0
      c2=(corX-A)*L
      d2=(a2*corX+b2*corY+c2*corZ)*(-1)
      a3=0.0
      b3=corZ*A
      c3=(-1)*A*(corY-L)
      d3=(a3*corX+b3*corY+c3*corZ)*(-1)
      a4=(corZ*L)*(-1)
      b4=0.0
      c4=corX*L
      d4=(a4*corx+b4*corY+c4*corZ)*(-1)
      y=h
      do u=2,ny-1
      x=h
      do j=2,nx-1
      Y3=(y-L)*(A-corX)/(L-corY)+A
      if ((corX*y/corY).le.x.and.(((corX-A)/corY)*y+A).ge.x)then
      mat(1,j,u)=(a1*x+b1*y+d1)*(-1)/c1
      else if ((corX*y/corY).ge.x.and.((y-L)*corX/(corY-L)).ge.x)then
      mat(1,j,u)=(a4*x+b4*y+d4)*(-1)/c4
      else if (((y-L)*corX/(corY-L)).le.x.and.Y3.ge.x)then
      mat(1,j,u)=(a3*x+b3*y+d3)*(-1)/c3
      else
      mat(1,j,u)=(a2*x+b2*y+d2)*(-1)/c2
      end if
      x=x+h
      end do
      y=y+h
      end do
      
      end