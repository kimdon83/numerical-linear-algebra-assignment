program exer_223

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n,i,j 
real(mp), allocatable :: b(:), ad(:), al(:), au(:)
real(mp) h, denominator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input
print *,"input denominator of h. (h=1/denominator)"
read (*,*) denominator          !16, 32, 64, 128
print *, "input value h=",1/denominator

h=1/denominator
n=((1/h))-1._mp
!h=1/128._mp 

!n=(1/h)-1

allocate (b(1:n), ad(1:n), al(1:n-1), au(1:n-1) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector

b=0._mp
do i=1,n
    b(i)=exp(sin(h*i))*h*h 
enddo
b(1)=b(1)+1 
b(n)=b(n)+1 

ad=2._mp
al=-1._mp
au=-1._mp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Thomas algorithm
ad(1)=ad(1) 
al(1)=al(1)/ad(1) 
do j=2,n-1
    ad(j)=ad(j)-al(j-1)*au(j-1) 
    al(j)=al(j)/ad(j) 
enddo
ad(n)=ad(n)-al(n-1)*au(n-1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! foward substitution
b(1)=b(1) 
do j=2,n
    b(j)=(b(j)-al(j-1)*b(j-1))/ad(j) 
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Backward substitution
b(n)=b(n)/ad(n) 
do j=n-1,1,-1
    b(j)=(b(j)-au(j)*b(j+1))/ad(j) 
enddo

print *, "b(solution)"
print *, b

!read (*,*)
end program exer_223