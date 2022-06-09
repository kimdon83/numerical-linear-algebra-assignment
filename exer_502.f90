program exer_502

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n,i,j,k
real(mp), allocatable :: A(:,:),b(:), v(:)
real(mp)  h,beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n=100
h=1/n
n=n-1
allocate(A(1:n,1:n),b(1:n)  )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector
A=0._mp
b=0._mp
A(1,1)=1
do i=2,n
    A(i,i)=2    
    A(i,i-1)=-1    
    A(i-1,i)=-1
enddo
b(1)=1

do j=1,n-1
    beta = sqrt(dot_product(A(j:n,j),A(j:n,j))) * A(j,j) / abs(A(j,j))
    allocate(v(1:n-j+1))
    v=A(j:n,j)
    v(1)=v(1)+beta
    A(j,j)=-beta
    A(j+1:n,j)=0
    do k = j+1,n
        A(j:n,k) = A(j:n,k) - 2*dot_product(v,A(j:n,k))/dot_product(v,v)*v
    enddo
    b(j:n) = b(j:n) - 2*dot_product(v,b(j:n))/dot_product(v,v) *v
    deallocate(v)
enddo

b(n)=b(n)/A(n,n);
do j = n-1,1-1
    b(j) = b(j) - dot_product(A(j,j+1:n),b(j+1:n))
    b(j) = b(j)/A(j,j)
enddo
!read (*,*) 

end program exer_502
