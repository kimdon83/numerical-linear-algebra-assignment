include "prec.f90"

program exer210

use prec
implicit none
integer(mpi)::n,i, k,j, temp, i_temp, j_temp
real(mp), allocatable :: A(:,:) , x_plot(:) , b(:), L(:,:) , U(:,:), y(:) , z(:) ,error(:,:)
real(mp)::h, denominator, i_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input

print *,"input denominator of h. (h=1/denominator)"
read (*,*) denominator
print *, "input value h=",1/denominator

h=1/denominator
n=1/h+1._mp

allocate(A(1:n,1:n),L(1:n,1:n), U(1:n,1:n) , x_plot(1:n) , b(1:n) ,y(1:n) , z(1:n) ,error(1:n,1:n))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector
A=0._mp
L=0._mp
U=0._mp
y=0._mp
z=0._mp
 
do i=1,n
 !   i_real=i
    x_plot(i)=(i)/(n-1)
    b(i)=exp(sin(x_plot(i)))*h*h
enddo
A(1,1)=b(1)
A(n,n)=b(n)

do i=2,n-1
   A(i,i)=2
   A(i,i-1)=-1
   A(i,i+1)=-1
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   doolittle

do i=1,n
    L(i,i)=1
enddo
do j=1,n
    do k=j,n            
        U(j,k)=(A(j,k)-dot_product(L(j,1:j-1),U(1:j-1,k)))
    enddo
    do k=j+1,n
        L(k,j)=(A(k,j)-dot_product(L(k,1:j-1),U(1:j-1,j)))/U(j,j)
    enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test the result of doolittle algorithm
!print *, "A"
!print *, A
!print *, "L"
!print *, L
!print *, "U"
!print *, U
!  do i=1,n
!     do j=1,n
!        error(i,j)=dot_product(L(i,1:n),U(1:n,j))
!     enddo
!  enddo
!  error=error-A
!print *, "error"
!print *, error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Forward substitution
y(1)=b(1)
do i=2,n
    y(i)=(b(i)-dot_product(L(i,1:i-1),y(1:i-1)))
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Backward substitution
z(n)=y(n)/U(n,n)
do i_temp=1,n-1
    i=n-i_temp
    z(i)=(y(i)-dot_product(U(i,i+1:n),z(i+1:n)))/U(i,i)
enddo

print *, "u"
print *, z

end program exer210
