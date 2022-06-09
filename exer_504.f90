program exer_504

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::N,i,j,k, len_v
real(mp), allocatable :: A(:,:),b(:), A_save(:,:),A_temp(:,:),b_temp(:)
real(mp)  h, rjk, rjj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N=100
h=1/N
N=N-1
allocate(A(1:N,1:N),b(1:N),A_save(1:n,1:n),A_temp(1:n,1:n) ,b_temp(1:n) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector
A=0._mp
b=0._mp
A(1,1)=1
do i=2,N
    A(i,i)=2    
    A(i,i-1)=-1    
    A(i-1,i)=-1
enddo
b(1)=1

A_save = A

do j = 1,N
    rjj = sqrt(dot_product( A(1:N,j),A(1:N,j) ))
    A(1:N,j) = A(1:N,j)/rjj
    if ( j<N ) then
        do k = j+1,N
            rjk = dot_product(A(1:N,j),A(1:N,k))
            A(1:N,k) = A(1:N,k) - A(1:N,j)*rjk
        enddo
    endif
enddo

do i=1,n
    do j=1,n
        A_temp(i,j)=A(j,i)
    enddo
enddo
A=A_temp
do i=1,n
    b_temp(i)=dot_product(A(i,1:n),b(1:n))
enddo
b = b_temp
do i=1,n
    do j=1,n
        A_temp(i,j)=dot_product(A(i,1:n),A_save(1:n,j))
    enddo
enddo
A=A_temp
A = A*A_save

b(n)=b(n)/A(n,n);
do j = n-1,1-1
    b(j) = b(j) - dot_product(A(j,j+1:n),b(j+1:n))
    b(j) = b(j)/A(j,j)
enddo
!read (*,*) 

end program exer_504
