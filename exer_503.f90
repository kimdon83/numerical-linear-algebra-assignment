program exer_503

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::N,i,j,k, len_v
real(mp), allocatable :: A(:,:),b(:), v(:),w(:)
real(mp)  h,beta,c,s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N=100
h=1/N
N=N-1
allocate(A(1:N,1:N),b(1:N)  )
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


do j = 1,N-1
    len_v = N -j +1
    allocate(v(1:len_v),w(1:len_v))
    do k = j+1,N
        c = A(j,j) / sqrt( A(j,j)**2 + A(k,j)**2 )
        s = A(k,j) / sqrt( A(j,j)**2 + A(k,j)**2 )

        v(1:len_v) = A(j,j:N)
        w(1:len_v) = A(k,j:N)
        A(j,j:N) = c*v + s*w
        A(k,j:N) = -s*v + c*w
        
        v(1) = b(j)
        w(1) = b(k)
        b(j) = c*v(1) + s*w(1)
        b(k) = -s*v(1) + c*w(1)
    enddo
    deallocate(v,w)
enddo

b(n)=b(n)/A(n,n);
do j = n-1,1-1
    b(j) = b(j) - dot_product(A(j,j+1:n),b(j+1:n))
    b(j) = b(j)/A(j,j)
enddo
!read (*,*) 

end program exer_503
