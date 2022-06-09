program exer_4_14

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n,j,k, i ,step         ! bl,bu,ml, mu  ,       pl,k,p , j_temp,i_temp ,index1
real(mp), allocatable :: A(:,:) , A_save(:,:) ,v(:,:),v_cap(:,:),lamda(:),x(:,:),v_save(:),temp(:,:),temp1(:) ! ,b_save(:),b_plot(:,:),v_tmp(:) !,G(:,:),u_ans(:) 
real(mp)  check,lamda_save 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n=10
n=n-1
allocate(A(1:n,1:n)  ,A_save(1:n,1:n) ,v(1:n,1:n) , lamda(1:n),v_cap(1:n,1:n),x(1:n,1:n),v_save(1:n),temp(1:n,1:n) ,temp1(1:n) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector
A=0._mp
A(1,1)=1
do i=2,n
    A(i,i)=2
    A(i,i-1)=-1
    A(i-1,i)=-1
enddo
A=A/((n+1)*(n+1))

A_save=A

step=0
check=1
lamda=1._mp
v=1._mp
x=0._mp
do i=1,n
    x(1,i)=1
enddo

!power method
do i=1,n
    if (i/=1) then
        do j=1,n
            do k=1,n
                temp(j,k)=v_cap(j,i-1)*v_cap(k,i-1)
            enddo
        enddo
        A=A-lamda(i-1)*temp

        check=1
    endif
    do step=1,150
        v_save=v(1:n,i)

        do j=1,n
             temp1(j)=dot_product(A(j,1:n),v(1:n,i))
        enddo

        v(1:n,i)=temp1(1:n)/lamda(i)
        lamda_save=lamda(i)

        lamda(i)=lamda_save*dot_product(v(1:n,i),v(1:n,i))/dot_product(v(1:n,i),v_save)
        check=abs((lamda(i)-lamda_save)/lamda(i));
        v_cap(:,i)=v(1:n,i)/sqrt(dot_product(v(1:n,i),v(1:n,i)))
    end do
enddo

print *,A_save
print * 
print *, check

print * ,"eigenvalue"
print *, lamda
print * , "eigenvector"

do i=1,n
    print *, v_cap(1:n,i)
    print *
enddo
read (*,*) 

end program exer_4_14
