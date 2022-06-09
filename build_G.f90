subroutine build_G(A,k,n,G)
  use prec
  implicit none

  integer(mpi)::i,k,n
  real(mp) :: A(1:n,1:n)
  real(mp), intent(inout):: G(1:n,1:n)

  G=0._mp
  do i=1,n
     G(i,i)=1._mp
  enddo
  
  do i=k+1,n
      G(i,k)=-A(i,k)/A(k,k)
  enddo

end subroutine build_G
