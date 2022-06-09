program exer_227

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n, mesh,index1,i,j,bl,bu,k,p,pl !,p,j,k,pl, j_temp,i_temp
real(mp), allocatable :: A(:,:) ,b(:), A_save(:,:) ,b_save(:), G(:,:), G_trans(:,:) 
real(mp) h, h_square, temp,test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input

h=1/100._mp    !h=0.01 
mesh=(1/h)-1   
n=mesh*mesh !n=100

allocate(A(1:n,1:n),b(1:n),G(1:n,1:n), G_trans(1:n,1:n) )  !,A_save(1:n,1:n),b_save(1:n)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector

A=0._mp
b=0._mp
h_square=h*h 
b=h_square 

do i=1,mesh
    do j=1,mesh
        index1=mesh*(i-1)+j 
        A(index1,index1)=4 
        if (j/=1) then
            A(index1,index1-1)=-1 
            A(index1-1,index1)=-1 
        endif
        if (i/=mesh)then
            A(index1,index1+mesh)=-1 
            A(index1+mesh,index1)=-1 
        endif
        
        b(index1)=2*(3.141592*3.141592)*sin(i*h)*sin(j*h)*h*h
        
    enddo
enddo
!A_save=A 
!b_save=b 

bl=mesh 
bu=mesh 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cholesky
do j=1,n
    temp=0 
    
    if (j/=1) then
        do p=1,j-1
            temp=temp+(A(j,p)*A(j,p)) 
        enddo
    endif
    A(j,j)=sqrt(A(j,j)-temp) 
    if (j/=n) then
        do k=j+1,n
            temp=0
            do pl=1,j-1
                temp=temp+A(k,pl)*A(j,pl) 
            enddo
            A(k,j)=(A(k,j)-temp)/A(j,j) 
        enddo
    endif
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test the result
!G=0._mp
!G_trans=0._mp
!do i=1,n
!    do j=1,i
!        G(i,j)=A(i,j)
!        G_trans(j,i)=A(i,j)
!    enddo
!enddo
!test=0._mp
!do i=1,n
!    do j=1,n
!    test=test+dot_product(G(i,1:n),G_trans(1:n,j))-A_save(i,j)
!    test=test+dot_product(G(i,1:n),G(j,1:n))
!    test=test+G(i,j)-G_trans(j,i)
!    enddo
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! foward substitution
b(1)=b(1)/A(1,1)
do j=2,n
    k=max(1,j-bl)
    temp=0
    do i=k,j-1
        temp=temp+A(j,i)*b(i)
    enddo
    b(i)=(b(i)-temp)/A(j,j)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Backward substitution
b(n)=b(n)/A(n,n)
do j=n-1,1,-1
    temp=0 
    k=min(n,j+bu)
    do i=j+1,k
        temp=temp+A(i,j)*b(i)
    enddo
    b(j)=(b(j)-temp)/A(j,j)
enddo

!print *, "A_save"
!print *, A_save
!print *, "b_save"
!print *, b_save
!print *, "A"
!print *, A
!print *, "G"
!print *, G
!print *, "G_trans"
!print *, G_trans
!print *, "test"
!print *, test
print *, b

print *, b
!read (*,*) 
open(1,file='b_cholesky.dat')
write(1,*) b

close(1)

!read (*,*) 
end program exer_227
