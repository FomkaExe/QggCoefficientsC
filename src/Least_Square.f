      program GaussMethod
       implicit double precision (a-h, o-z)
      parameter(K = 14,K1 = 15, N= 9)
      
	dimension   X(K,N), Y(K,N),NI(K),D(K1),Z(K1), C(K1,K1)
	dimension   YC(K,N),dyc(K,N), sumdy(K)
	dimension   dY(K,N),p(K,N)
	
	open(11, file = 'scaleTa_IZ2.dat')! otkryvaem vhodnoy fail
 	open(21, file = 'outdata_Ta_IZ2.dat')
     
	
	Do 4 i = 1,K1
	D(i) = 0.0
	Do 2 j = 1,K1
	C(i,j) = 0.0
2     continue	
 4    continue
	read(11,*) K0
      IF(K0 <= K) THEN
	read(11,*) (NI(i), i=1,K0)
	Do 10 i = 1,K0
	read(11,*)
	NII = NI(i)
	Do 5 j = 1,NII
	read(11,*) iz, in, y(i,j), dy(i,j)
	x(i,j) = float(in)
	p(i,j) =0.15817/ dy(i,j)
	p(i,j) = 1
5     continue	
10    continue	
      ELSE
      GO TO 999
	END IF
	k01 = k0+1
	Do 20 i = 2,K01
	NII = NI(i-1)
	Do 15 j = 1,NII
	C(1,1) = C(1,1) +p(i-1,j)*X(i-1,j)**2
	C(1,i) = C(1,i) +p(i-1,j)* X(i-1,j)
	C(i,1) = C(i,1) + p(i-1,j)*X(i-1,j)
	D(1) = D(1) +p(i-1,j)* x(i-1,j)*y(i-1,j)
	D(i) = D(i) + p(i-1,j)*y(i-1,j)
	C(i,i) = C(i,i)+p(i-1,j)
15    continue	
20    continue	
      call Gauss_method( k0,k01, c, d,z)
      write(21, *) 'a=', z(1)
	Do 40 ik = 2,K01
      write(21, *) 'b(',ik-1,')=', z(ik)
40    continue	

ccc ------------------------------
      a = z(1)
	Do 60 ik = 1,K0
 	b = z(ik+1)
 	sumdy(ik) =0.0 
 	Do 50 in = 1,NI(ik)
      YC( ik,in) = a*X( ik,in)+b
      dyc(ik,in) = YC(ik,in) - Y(ik,in)
      if(abs(YC( ik,in)).gt.1.0) then
 	sumdy(ik) =sumdy(ik)+dyc(ik,in)**2/abs(yc(ik,in))
 	end if
50    continue	
      write(21,*) 'sigma', ik,NI(ik), sumdy(ik)
60    continue	
ccc ------------------------------

999    continue	
      stop
      end
      
      subroutine  Gauss_method(k,k1, c1, d1,b)
      implicit double precision (a-h, o-z)
      dimension c1(k1,K1), d1(k1),b(k1)
      dimension cnew(k1,K1), q(k1,k1)
      dimension  indD(k1),cold(K1)
	
	Do 2 i = 1,K1
	indD(i) = i
2     continue
	Do 5 i = 1,K1
 	Do 4 j = 1,K1
	Cnew(i,j) = 0.0
4     continue
5     continue
ccc privedenie matricy koefficientov k treugol'noy 
 	Do 100 J = 1,K
 	I0 = J
 	iflag = -1
 	J1 = J+1 
 	cmax = C1(J,J)
 	Do 10 is = J1,K1
	If(abs(cmax).lt.abs(C1(is,J))) then
	I0 = Is
	cmax = C1(is,J)
	iflag = 1
	end if 
10    continue	
      if (iflag.eq.1) then
 	Do 20 jc = 1,K1
	cold(jc) = C1(J,jc) 
	c1(J,jc) = c1(I0,jc)
	c1(i0,jc) = cold(jc)
20    continue	
	
	dold = d1(J)
	d1(J) = d1(I0)
	d1(I0) = dold
	indDold = indD(J) 
	indD(j) = indD(I0)
	indD(i0) = indDold  
      end if
 	Do 30 i = J1,K1
 	q(i,j) = C1(i,j)/C1(j,j)
 	cnew(i,j)= q(i,j)
	D1(i) = D1(i) -q(i,j)*D1(j)
 	Do 25 j2 = j1,K1
	C1(i,j2) = C1(i,j2) -q(i,j)*C1(j,j2)
25     continue	 
30    continue	 
 
100    continue
ccc nahozhdenie rezul'tatov
      B(K1) = D1(K1)/C1(K1,K1) 
	Do 200 i = 1,K
      IM = K-I+1
      aux = D1(IM)
      	do 110 j = IM, K
      aux= AUX-C1(im, j+1)*b(j+1)
 110    continue
      b(im) = aux/c1(im,im)	
200    continue	 
      return
      end

