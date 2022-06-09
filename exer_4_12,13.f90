program exer_4_12

implicit none

integer, parameter::mpi=4
integer, parameter::mp=selected_real_kind(p=15,r=307)
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer(mpi)::n,j,k, i ,step         ! bl,bu,ml, mu  ,       pl,k,p , j_temp,i_temp ,index1
real(mp), allocatable :: A(:,:) , A_save(:,:) ,v(:,:),v_cap(:,:),lamda(:),x(:,:),v_save(:),temp(:,:),temp1(:) ! ,b_save(:),b_plot(:,:),v_tmp(:) !,G(:,:),u_ans(:) 
real(mp)  check,lamda_save 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n=3
allocate(A(1:n,1:n)  ,A_save(1:n,1:n) ,v(1:n,1:n) , lamda(1:n),v_cap(1:n,1:n),x(1:n,1:n),v_save(1:n),temp(1:n,1:n) ,temp1(1:n) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initinalize matrix and vector
A=reshape((/-261,-530,-800,209,422,631,-49,-98,-144/),(/N,N/))
A_save=A
step=0
check=1
lamda=1._mp
v=1._mp
x=reshape((/1,0,0,1,0,0,1,0,0/),(/N,N/))
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

print * 
print *, check

print * ,"eigenvalue"
print *, lamda
print * , "eigenvector"
print *, v_cap(1:n,1)
print *
print *, v_cap(1:n,2)
print *
print *, v_cap(1:n,3)
read (*,*) 

end program exer_4_12
