program exer_215

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))


integer(mpi)::n ,i, j, k, p 
real(mp), allocatable :: A(:,:) , A_save(:,:), r(:), w(:), d(:) ,L(:,:),Diag(:,:),Mt(:,:), A_test(:,:), A_temp(:,:)
real(mp) temp, temp1, temp2
!integer(mpi), allocatable::perm(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input

n=100._mp

allocate(A(1:n,1:n) ,A_save(1:n,1:n), r(1:n), w(1:n) ,d(1:n), L(1:n,1:n))
allocate( Diag(1:n,1:n),  Mt(1:n,1:n) , A_test(1:n,1:n), A_temp(1:n,1:n) )     
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
A_save=A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LDMt algorithm
do j=1,n
    temp=0._mp
    temp1=0._mp
    temp2=0._mp
    do p=1,j-1
        r(p)=d(p)*A(p,j)
        w(p)=A(j,p)*d(p)
        temp=temp+w(p)*A(p,j)
    enddo
    d(j)=A(j,j)-temp;
    do k=j+1,n
        do p=1,j-1
            temp1=temp1+A(k,p)*r(p)
        enddo
        A(k,j)=(A(k,j)-temp1)/d(j)
        do p=1,j-1
            temp2=temp2+A(p,k)*w(p)
        enddo
        A(j,k)=(A(j,k)-temp2)/d(j)
    enddo
enddo

L=0._mp
Diag=0._mp
Mt=0._mp

do i=1,n
    L(i,i)=1
    Mt(i,i)=1
    Diag(i,i)=d(i)
enddo

do i=2,n
    do j=1,i-1
        L(i,j)=A(i,j)
        Mt(j,i)=A(i,j)
    enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test the result

!print *, A
!print *, "L"
!print *, L
!print *, "D"
!print *, Diag
!print *, "Mt"
!print *, Mt

do i=1,n
    do j=1,n
        A_temp(i,j)=dot_product(L(i,1:n),Diag(1:n,j))
    enddo
enddo
do i=1,n
    do j=1,n
        A_test(i,j)=dot_product(A_temp(i,1:n),Mt(1:n,j))
    enddo
enddo

A_test=A_test-A_save
temp=sum(A_test)
print *, "sum of all element of differrence between A_save and L*U"
print *, temp

!read (*,*)

end program exer_215
