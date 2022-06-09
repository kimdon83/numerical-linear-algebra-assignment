program exer211

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n,i, k,j, temp, i_temp, j_temp, step
real(mp), allocatable :: A(:,:) , x_plot(:) , b(:), L(:,:) , U(:,:), x_step(:) ,error(:,:), A_inverse(:,:)
real(mp)::h, denominator, i_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input

n=10._mp

allocate(A(1:n,1:n),L(1:n,1:n), U(1:n,1:n) , x_plot(1:n) , b(1:n) , x_step(1:n) ,error(1:n,1:n) , A_inverse(1:n,1:n))

do step=1,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector

b=0._mp
b(step)=1

A=0._mp
A(1,1)=2
A(1,2)=-1
A(n,n)=2
A(n,n-1)=-1
do i=2,n-1
   A(i,i)=2
   A(i,i-1)=-1
   A(i,i+1)=-1
enddo

L=0._mp
U=0._mp
!y=0._mp
x_step=0._mp

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

do i=2,n
    b(i)=(b(i)-dot_product(L(i,1:i-1),b(1:i-1)))
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Backward substitution

x_step(n)=b(n)/U(n,n)
do i_temp=1,n-1
    i=n-i_temp
    x_step(i)=(b(i)-dot_product(U(i,i+1:n),x_step(i+1:n)))/U(i,i)
enddo

A_inverse(1:n,step)=x_step

enddo

print *, "A_inverse"

do i=1,n
    print *, A_inverse(1:n,i)
enddo

end program exer211
