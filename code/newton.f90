!香蕉函数 的牛顿法 实现
module dat1
    implicit none
    real,dimension(:,:),allocatable,save::Hess,ivHess
    real,dimension(:),allocatable,save::Gra
    real::fx,mgra
    real,dimension(1:2)::X1,X0
endmodule

program main
    use dat1
    implicit none
    integer::i,j,k
    allocate (Hess(1:2,1:2))
    allocate(Gra(1:2))!Gradient,
    X0=[-1.2,1.0]!start
    fx=100*(x0(2)-x0(1)**2)**2+(1-x0(1))**2
    print*," The location of start point :    (",X0(1),",",X0(2),")"
    print*,"The original value of fx",fx
    do i=1,10000,1
         call G(X0,Gra,mGra,Hess)
         write(*,'("Hess :",1x)');write(*,10)(Hess(j,:),j=1,2)
         call nizhen(2)
          if (mGra<0.0001) exit
          X1=X0-matmul(ivHess,Gra)
          X0=X1
          fx=100*(x0(2)-x0(1)**2)**2+(1-x0(1))**2
          !fx=-fx
          10 format(2f10.4)
          write (*,'("Location",$)');write(*,10)(X0);
          write(*,'("Gradient",$)');write(*,10)(gra);
          !write(*,'("Muo     ",$)');write(*,10)(mgra);
          print*,"The value of the function";write(*,10)(fx)
          print*,"Finished one circulation "
          k=i
     enddo
     print*,"The value of the function",fx
     write(*,'("Gradient",$)');write(*,10)(gra);
     print*,k," step"
endprogram

subroutine G(x,Gx,mg,he)
!该步通过输入坐标X，计算梯度Gx，模mg, Hessian metrix (He)
    use dat1
    implicit none
    real::mg
    real,dimension (1:2)::x,Gx
    real,dimension(1:2,1:2)::He
    He(1,1)=-400*x(1)+1200*x(1)**2+2! Input the Hessian metrix
    He(1,2)=-400*x(1)
    He(2,1)=He(1,2)
    He(2,2)=200
    !He=-He!transfer the minus
    Gx(1)=-400*x(1)*(x(2)-x(1)**2)-2+2*x(1)
    Gx(2)=200*(x(2)-x(1)**2)
    !Gx=-Gx
    mg=(Gx(1)**2+Gx(2)**2)**(0.5)
endsubroutine

! aa为原矩阵，b为存放aa的逆矩阵，n为矩阵aa的维数
subroutine nizhen(n)
  use dat1
  implicit none
  integer n,i,j,k
  real,dimension(1:2,1:2)::b,a
  a=Hess
  b=0
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
   ivhess=b
   write(*,'("The inversemetrix is :",1x)')
   write(*,10)(ivhess(i,:),i=1,2)
end
