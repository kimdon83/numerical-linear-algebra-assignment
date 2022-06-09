
subroutine Matrix_multi1(L,G_t,n)
  use prec
  implicit none

  integer(mpi)::i,j,n
  real(mp)::m(1:n,1:n),G_t(1:n,1:n)
  real(mp),intent(inout)::L(1:n,1:n)

  m=L
  do i=1,n
     do j=1,n
        L(i,j)=dot_product(m(i,1:n),G_t(1:n,j))
     enddo
  enddo
end subroutine Matrix_multi1
