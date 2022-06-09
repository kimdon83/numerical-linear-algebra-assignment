program exer_221

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n_sub, n,j, i, bl,bu,ml, mu  ,       pl,k,p , j_temp,i_temp ,index1
real(mp), allocatable :: A(:,:) ,b(:), A_save(:,:) ,b_save(:),b_plot(:,:) !,G(:,:),upper(:,:),y(:),u_ans(:) 
real(mp)  h, temp, temp2,h_square

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input
h=1/100._mp
n_sub=1/h-1
n=n_sub*n_sub
h_square=h*h
allocate(A(1:n,1:n)  ,b(1:n),A_save(1:n,1:n)  ,b_save(1:n),b_plot(1:n_sub+2,1:n_sub+2)  )  !,G(1:n,1:n),  upper(1:n,1:n), y(1:n), u_ans(1:n) ) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector

A=0._mp
b=0._mp
b=h_square
do i=1,n_sub
    do j=1,n_sub
        index1=n_sub*(i-1)+j 
        A(index1,index1)=4 
        if (j/=1)then
            A(index1,index1-1)=-1 
            A(index1-1,index1)=-1 
        endif
        if (i/=n_sub)then
            A(index1,index1+n_sub)=-1 
            A(index1+n_sub,index1)=-1 
        endif
        
        if (j==1)then
            b(index1)=b(index1)+1 
        endif
        if (j==n_sub) then
            b(index1)=b(index1)+1 
        endif
        if (i==1)then
            b(index1)=b(index1)+1 
        endif
        if (i==n_sub)then
            b(index1)=b(index1)+1 
        endif
       
    enddo
enddo
A_save=A 
b_save=b 

bl=n_sub 
bu=n_sub 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LU_decomposition w/o pivoting with bandwitdh
do j=1,n-1
    ml=min(j+bl,n)
    do i=j+1,ml
        A(i,j)=A(i,j)/A(j,j)
    enddo
    do k=j+1,ml
        mu=min(j+bu,n)
        do i=j+1,n
            A(k,i)=A(k,i)-A(k,j)*A(j,i)
        enddo
    enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=2,n
    k=max(1,j-bl)
    temp=0
    do i=k,j-1
        temp=temp+A(j,i)*b(i)
    enddo
    b(j)=b(j)-temp
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   back substitution
b(n)=b(n)/A(n,n)
do j=n-1,-1,1
    temp=0
    k=min(j+bu,n)
    do i=j+1,k
        temp=temp+A(j,i)*b(i)
   enddo
   b(j)=(b(j)-temp)/A(j,j)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make answer to matrix form to view at excel
b_plot=0._mp
do i=1,n_sub
    do j=1,n_sub
        index1=n_sub*(i-1)+j
        b_plot(i+1,j+1)=b(index1)
    enddo
enddo
do i=1,n_sub+2
    b_plot(n_sub+2,i)=1
    b_plot(i,n_sub+2)=1
    b_plot(1,i)=1
    b_plot(i,1)=1
enddo

print *, b_plot
!read (*,*) i
open(1,file='b_plot.dat')
!do i=1,n_sub+2
        write(1,*) b_plot !(i,:) 
!enddo
close(1)

end program exer_221
