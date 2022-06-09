include "prec.f90"

program exer209

use prec
implicit none
integer(mpi)::n, i  ,k,j  , temp, i_temp, j_temp
real(mp), allocatable :: A(:,:), b(:) ,L(:,:) , U(:,:), y(:) , x(:)  ,error(:,:)  !, x_plot(:) , ,  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input
n=100._mp

allocate(A(1:n,1:n) , b(1:n) ,L(1:n,1:n), U(1:n,1:n)   ,y(1:n) , x(1:n)   ,error(1:n,1:n)        )!, , x_plot(1:n)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector
A=0._mp
b=0._mp
L=0._mp
U=0._mp
y=0._mp
x=0._mp

A(1,1)=2
A(1,2)=-1
A(n,n-1)=-1
A(n,n)=2
b(1)=1

do i=2,n-1
   A(i,i)=2
   A(i,i-1)=-1
   A(i,i+1)=-1
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   crout

do i=1,n
    U(i,i)=1
enddo
do j=1,n
    do k=j,n
        L(k,j)=(A(k,j)-dot_product(L(k,1:j-1),U(1:j-1,j)))
    enddo
    do k=j+1,n            
        U(j,k)=(A(j,k)-dot_product(L(j,1:j-1),U(1:j-1,k)))/L(j,j)
    enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test the result of crout algorithm
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
! enddo
!  error=error-A
!print *, "error"
!print *, error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Forward substitution
y(1)=b(1)/L(1,1)
do i=2,n
    y(i)=(b(i)-dot_product(L(i,1:i-1),y(1:i-1)))/L(i,i)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Backward substitution
x(n)=y(n)
do i_temp=1,n-1
    i=n-i_temp
    x(i)=y(i)-dot_product(U(i,i+1:n),x(i+1:n))
enddo

print *, "answer x"
print *, x

end program exer209
