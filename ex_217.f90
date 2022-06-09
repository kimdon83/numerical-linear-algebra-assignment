program exer_217

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n, domain,i  ,p,j,k,pl, j_temp,i_temp
real(mp), allocatable :: A(:,:) ,b(:), r(:), d(:), L(:,:), Diag(:,:), y(:),u_ans(:), upper(:,:) 
real(mp) h, temp, temp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input

h=0.01_mp
domain=1/h
n=domain-1

allocate(A(1:n,1:n),b(1:n),r(1:n), d(1:n), L(1:n,1:n),Diag(1:n,1:n) , y(1:n), u_ans(1:n), upper(1:n,1:n) )  

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
b=h*h
b(n)=b(n)+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LDLt reduced version
do j=1,n
    temp=0
    temp2=0
  
    do p=1,j-1
        r(p)=d(p)*A(j,p)
        temp=temp+r(p)*A(j,p)
    enddo
    d(j)=A(j,j)-temp
    do k=j+1,n
        do pl=1,j-1
            temp2=temp2+A(k,pl)*r(pl)
        enddo
        A(k,j)=(A(k,j)-temp2)/d(j)
    enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! show the result of code and test the result
L=0._mp
do i=1,n
    L(i,i)=1
enddo
Diag=0._mp
do i=1,n
    Diag(i,i)=1
enddo

do i=2,n
    do j=1,i-1
        L(i,j)=A(i,j)
    enddo
enddo
!do i=1,n
!    Diag(i,i)=d(i)
!enddo

!U=D*L';
do i=1,n
    do j=i,n
        upper(i,j)=d(i)*L(j,i)
    enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! foward substitution
y(1)=b(1)
do i=2,n
    temp=0
    do j=1,i-1
        temp=temp+L(i,j)*y(j)
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

print *,0, u_ans, 1

!read (*,*)
end program exer_217
