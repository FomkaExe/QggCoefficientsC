      program GaussMethod
      implicit double precision (a-h, o-z)
      parameter (K = 7, K1 = 8, N = 11)

      dimension X(K,N), Y(K,N), NI(K), D(K1), Z(K1), C(K1,K1)
      dimension YC(K,N), dyc(K,N), sumdy(K)

      open(11, file = 'input/OTa_qgg.dat')
      open(21, file = 'output/OTa_qgg_out.dat')

      do 4 i = 1, K1
         D(i) = 0.0d0
         Z(i) = 0.0d0
         do 2 j = 1, K1
            C(i,j) = 0.0d0
 2       continue
 4    continue

c input data
      read(11,*) K0
      write(*,*) 'K0 = ', K0

      if (K0 .le. K) then
         read(11,*) (NI(i), i = 1, K0)
         write(*,*) 'NI: ', (NI(i), i = 1, K0)

         do 10 i = 1, K0
            read(11,*)
            NII = NI(i)
            do 5 j = 1, NII
               read(11,*) iz, in, x(i,j), y(i,j)
               write(*,*) 'i=', i, ' j=', j, ' iz=', iz, ' in=', in
               write(*,*) 'x=', x(i,j), ' y=', y(i,j)
 5          continue
 10      continue
      else
         go to 999
      end if

c Least squares matrix
      do 20 i = 2, K1
         NII = NI(i-1)
         C(i,i) = dfloat(NII)
         do 15 j = 1, NII
            C(1,1) = C(1,1) + X(i-1,j)**2
            C(1,i) = C(1,i) + X(i-1,j)
            C(i,1) = C(i,1) + X(i-1,j)
            D(1)   = D(1)   + X(i-1,j)*Y(i-1,j)
            D(i)   = D(i)   + Y(i-1,j)
 15      continue
 20   continue

      write(*,*) 'Matrix C:'
      do 30 i = 1, K1
         write(*,*) (C(i,j), j = 1, K1)
 30   continue

      write(*,*) 'Vector D:'
      write(*,*) (D(i), i = 1, K1)

      call Gauss_method(K1, C, D, Z)

      write(21,*) 'a=', z(1)
      do 40 ik = 2, K1
         write(21,*) 'b(', ik-1, ')=', z(ik)
 40   continue

      a = z(1)

c solution accuracy
      do 60 ik = 1, K0
         b = z(ik+1)
         sumdy(ik) = 0.0d0
         do 50 in = 1, NI(ik)
            YC(ik,in)  = a * X(ik,in) + b
            dyc(ik,in) = YC(ik,in) - Y(ik,in)
            if (abs(YC(ik,in)) .gt. 1.0d0) then
               sumdy(ik) = sumdy(ik) +
     &                     dyc(ik,in)**2 / abs(YC(ik,in))
            end if
 50      continue
         write(21,*) 'sigma', ik, NI(ik), sumdy(ik)
 60   continue

 999  continue
      stop
      end


      subroutine Gauss_method(nsys, a, rhs, sol)
      implicit double precision (a-h, o-z)
      parameter (K1MAX = 8)
      dimension a(K1MAX,K1MAX), rhs(K1MAX), sol(K1MAX)

c forward elimination with partial pivoting
      do 100 j = 1, nsys-1
         ipiv = j
         pmax = abs(a(j,j))

         do 10 i = j+1, nsys
            if (abs(a(i,j)) .gt. pmax) then
               pmax = abs(a(i,j))
               ipiv = i
            end if
 10      continue

         if (pmax .lt. 1.0d-20) then
            write(*,*) 'Error: singular matrix at step ', j
            stop
         end if

         if (ipiv .ne. j) then
            do 20 m = 1, nsys
               t = a(j,m)
               a(j,m) = a(ipiv,m)
               a(ipiv,m) = t
 20         continue

            t = rhs(j)
            rhs(j) = rhs(ipiv)
            rhs(ipiv) = t
         end if

         do 30 i = j+1, nsys
            factor = a(i,j) / a(j,j)
            a(i,j) = 0.0d0
            do 25 m = j+1, nsys
               a(i,m) = a(i,m) - factor * a(j,m)
 25         continue
            rhs(i) = rhs(i) - factor * rhs(j)
 30      continue

 100  continue

      if (abs(a(nsys,nsys)) .lt. 1.0d-20) then
         write(*,*) 'Error: zero last diagonal ', a(nsys,nsys)
         stop
      end if

c back substitution
      sol(nsys) = rhs(nsys) / a(nsys,nsys)

      do 200 i = nsys-1, 1, -1
         s = rhs(i)
         do 110 j = i+1, nsys
            s = s - a(i,j) * sol(j)
 110     continue
         if (abs(a(i,i)) .lt. 1.0d-20) then
            write(*,*) 'Error: zero diagonal in back substitution'
            write(*,*) 'i=', i, ' a(i,i)=', a(i,i)
            stop
         end if
         sol(i) = s / a(i,i)
 200  continue

      return
      end
