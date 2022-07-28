      subroutine main(t, h, e, aw, s, dC, dL, em, dZ0)

      implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
      dimension z1(15),betam(15)  
      dimension s(2)
      dimension z2(15),z3(15),x3(15),qwork(400)
	  pi =  3.141593d0	      
	  zero =(0.d0,0.d0)   ! Нуль
         zi = (0.d0,1.d0)	  ! Мнимая единица
         nn = 15             ! Количество вершин
         w2 = aw/2        ! Половина ширины полосы
! Координаты вершин
      z1(1) =  0
      z1(2) = -w2
      z1(3) = -w2+zi*t
      z1(4) =  w2+zi*t
      z1(5) =  w2
      z1(6) =  0
      z1(7) = -zi*h
      z1(8) =  100000
      z1(9) =  w2+s(2)
      z1(10)=  w2+s(2)+zi*t 
      z1(11)=  100000
      z1(12)= -w2-s(1)+zi*t
      z1(13)= -w2-s(1)
      z1(14)= -100000
      z1(15)= -zi*h
! Углы в долях pi -1
      betam(1) = -0.5
      betam(2) =  0.5
      betam(3) =  0.5
      betam(4) =  0.5
      betam(5) =  0.5
      betam(6) = -0.5
      betam(7) = -0.5
      betam(8) =  -1
      betam(9) =  0.5
      betam(10) = 0.5
      betam(11) = -2
      betam(12) = 0.5
      betam(13) = 0.5
      betam(14) = -1
      betam(15) =-0.5
      z10 = dcmplx(0, h)			! Задан. конф.центра z10 в исх. многоугол. z1
! Compute nodes and weights for parameter problem:
       nptsq =12								   ! nptsq =12
       call qinit(nn,betam,nptsq,qwork)
! Solve parameter problem:
!   (initial guess must be given to avoid accidental exact solution)
       iguess = 1
! Начальное приближение точек для круга
       do 1 k = 1,nn
	z2(k) = exp(dcmplx(0.d0, k-4.d0))   
1	continue
! Точность
       tol = 1.d-8	! 1.d-14
       call scsolv(iprint,iguess,tol,errest,nn,c,z2,z10,z1,betam,nptsq,qwork)
!	z1(k) - вершины исх.многоугольника (в оригинале были обозначены как w(k) )
!	z2(k) - круговая каноническая область (в оригинале обозначена  как z(k) )
!	z3(k) - верхняя полуплоскость
! Дробно-линейное отображение круга z2(k) на верх.полуплоскость (окружность -> веществ.ось х3)
	z30 = zi;				      ! Задание конформ.центра z30 на верхней полуплоскости z3
	z301=-zi;		
	do 2 k = 1, nn-1
	z3(k)=(z2(k)*z301-z30)/(z2(k)-1.) ! Отображающая функция z3
	x3(k)=dreal(z3(k));			    ! выделение веществен.части (т.е. отбрасывание мнимой)
2	continue	  
!--------------------------------------------------------------------------
! Расчёт ёмкости всей структуры целиком [Книга 2007, с.47]
	a=x3(6)-x3(1)
	b=x3(7)-x3(6)
	ak=dsqrt(b/(a+b))
	bC=1./aKK1(ak)
!--------------------------------------------------------------------------
! Расчёт ёмкости в диэлектрической подложке (substrate) [Книга 2007, с.47]
	x12=dexp(pi*s(1)/h)
	x13=dexp(pi*(s(1)+aw)/h)
	x14=dexp(pi*(s(1)+aw+s(2))/h)
	a=x13-x12
	b1=x12-1.
	b2=x14-x13
	ak=dsqrt(b1*b2/((a+b1)*(a+b2)))
       ddC=1./aKK1(ak)
!--------------------------------------------------------------------------
       dC = 8.854*(bC + (e-1.)*ddC)
       dL=0.4*pi/bC  
       dZ0=dsqrt(1.d6*dL/dC)
       em = (dC/8.854)/bC
!--------------------------------------------------------------------------
    end 