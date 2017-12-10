
! aa为原矩阵，b为存放aa的逆矩阵，n为矩阵aa的维数
subroutine nizhen(aa,b,n)
  use dat1
  integer n,i,j,k
  real*8,dimension(1:2,1:2)::aa,b,a
  a=aa
  do i=1,n
    b(i,i)=1
  enddo
  do i=1,n
    b(i,:)=b(i,:)/a(i,i)
    a(i,i:n)=a(i,i:n)/a(i,i)
   do j=i+1,n
     do k=1,n
   b(j,k)=b(j,k)-b(i,k)*a(j,i)
  enddo
  a(j,i:n)=a(j,i:n)-a(i,i:n)*a(j,i)
 enddo
  enddo
  do i=n,1,-1
    do j=i-1,1,-1
   do k=1,n
     b(j,k)=b(j,k)-b(i,k)*a(j,i)
    enddo
  enddo
   enddo
   10 format (2f10.5)
   write(*,10)(B(i,:),i=1,2)
end
