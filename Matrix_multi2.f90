
subroutine Matrix_multi2(G,A_u,n)
  use prec
  implicit none

  integer(mpi)::i,j,n
  real(mp)::m(1:n,1:n),G(1:n,1:n)
  real(mp), intent(inout)::A_u(1:n,1:n)
  
  m=A_u
  do i=1,n
     do j=1,n
        A_u(i,j)=dot_product(G(i,1:n),m(1:n,j))
     enddo
  enddo
end subroutine Matrix_multi2
