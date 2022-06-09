include "prec.f90"

program exercies2_6
  use prec
  implicit none

  integer(mpi)::k,n,   i,j
  real(mp), allocatable ::A(:,:),A_u(:,:), delta(:,:), G(:,:), G_t(:,:),L(:,:),eye(:,:)

  n=3

  allocate(A(1:n,1:n),A_u(1:n,1:n), delta(1:n,1:n), G(1:n,1:n), G_t(1:n,1:n),L(1:n,1:n),eye(1:n,1:n))
 
A(1,1)=2._mp
A(1,2)=-1._mp
A(1,3)=0._mp
A(2,1)=-1._mp
A(2,2)=2._mp
A(2,3)=-1._mp
A(3,1)=0._mp
A(3,2)=-1._mp
A(3,3)=2._mp

A_u=A

L=0._mp

do i=1,n
   L(i,i)=1._mp
enddo

eye=0._mp

do i=1,n
   eye(i,i)=1._mp
enddo

!G=0._mp
!k=2.mp
do k=2,n
      call build_G(A_u,k-1,n,G)
      call Matrix_multi2(G,A_u,n)
      G_t=2*eye-G
!call build_G_t(G,n)
      call Matrix_multi1(L,G_t,n)
enddo

  do i=1,n
     do j=1,n
        delta(i,j)=dot_product(L(i,1:n),A_u(1:n,j))
     enddo
  enddo
  delta=delta-A

!!L*A_u-A!!
n=3
print *, "A"
do i=1,n
    print *, (A(i,j), j=1,3) 
end do

print *, "L"
do i=1,n
    print *, (L(i,j), j=1,3) 
end do


print *, "A_u"
do i=1,n
    print *, (A_u(i,j), j=1,3) 
end do

print *, "delta"
do i=1,n
    print *, (delta(i,j), j=1,3) 
end do

end program exercies2_6
