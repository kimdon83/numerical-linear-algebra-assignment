program exer_214

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n ,i, jsave, j , js, k ,l !,i, k,j, temp, i_temp, j_temp, step
real(mp), allocatable :: A(:,:) , A_save(:,:) , swap(:)
real(mp) m
integer(mpi), allocatable::perm(:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input

n=10._mp

allocate(A(1:n,1:n) ,A_save(1:n,1:n) , swap(1:n), perm(1:n))     

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

do j=1,n
    jsave=maxloc(abs(A(j:n,j)),1)
    js=jsave+j-1
    perm(j)=js
    swap(1:n)=A(j,1:n)
    A(j,1:n)=A(js,1:n)
    A(js,1:n)=swap

    m=A(j,j)
    A(j,j)=1/A(j,j)
    A(j,1:j-1)=-A(j,j)*A(j,1:j-1)
    A(j,j+1:n)=-A(j,j)*A(j,j+1:n)

    A(1:j-1,j)=A(j,j)*A(1:j-1,j)
    A(j+1:n,j)=A(j,j)*A(j+1:n,j)

    do k=1,n
        if(k/=j) then
            do l = 1,n
                if (l/=j) A(k,l)=A(k,l)+A(k,j)*m*A(j,l)
            end do
        end if
    end do
enddo

do j=n-1,1,-1
    if(j/=perm(j)) then
        swap(1:n)=A(1:n,j)
        A(1:n,j)=A(1:n,perm(j))
        A(1:n,perm(j))=swap(1:n)
    end if
enddo


!print *, "A_inverse"

!do i=1,n
!    print *, A_inverse(1:n,i)
!enddo

print *,"inverse of A"
do i=1,n
    print *, A(1:n,i)
enddo

!print *, A_save 
print *, "perm"

print *, perm

end program exer_214
