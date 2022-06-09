program exer_218

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n,j, i, pl,k,p , j_temp,i_temp 
real(mp), allocatable :: A(:,:)  ,b(:) ,G(:,:),upper(:,:),y(:),u_ans(:) 
real(mp)  temp, temp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input

n=100

allocate(A(1:n,1:n)  ,b(1:n) ,G(1:n,1:n),  upper(1:n,1:n), y(1:n), u_ans(1:n) ) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector

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

b=0._mp
b(1)=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cholesky
do j=1,n
    temp=0
    temp2=0
    do p=1,j-1
        temp=temp+A(j,p)*A(j,p)
    enddo
    A(j,j)=sqrt(A(j,j)-temp)
    do k=j+1,n
        do pl=1,j-1
            temp2=temp2+A(k,pl)*A(j,pl)
        enddo
        A(k,j)=(A(k,j)-temp2)/A(j,j)
    enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! show the result of code and test the result
do i=1,n
    do j=1,i
        G(i,j)=A(i,j)
        upper(j,i)=G(i,j)
    enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! foward substitution
y(1)=b(1)
do i=2,n
    temp=0
    do j=1,i-1
        temp=temp+G(i,j)*y(j)
    enddo
    y(i)=(b(i)-temp)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Backward substitution
u_ans(n)=y(n)/upper(n,n)
do i_temp=1,n-1
    i=n-i_temp
    temp=0;
    do j_temp=1,n-i
        j=n+1-j_temp
        temp=temp+upper(i,j)*u_ans(j)
    enddo
    u_ans(i)=(y(i)-temp)/upper(i,i)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *,u_ans

!read (*,*)
end program exer_218
