!最速下降法，对f(x,y)=x^2+4y^2 的实现
module dat1
    implicit none
    real*8,dimension(:),allocatable,save::Gra,Sk,X1,X0
    real*8,dimension(:,:),allocatable,save::Hess
    real*8::Alpha,mGra,fx
    endmodule

program main
    use dat1
    implicit none
    integer::i,pau
    allocate (Hess(1:2,1:2))
    allocate(Gra(1:2),Sk(1:2),X1(1:2),X0(1:2))!Gradient,
    X0=[2,2]
    call G(X0,Gra,mGra,Sk,Hess,Alpha)
    do i=1,10,1
          if (mGra<0.01) exit
          X1=X0+Alpha*Sk
          call G(X1,Gra,mGra,Sk,Hess,Alpha)
          X0=X1
          fx=x0(1)**2+4*x0(2)**2
          10 format(5f10.5)
          write (*,'("Location",$)');write(*,10)(X0);
          write(*,'("Gradient",$)');write(*,10)(gra);
          write(*,'("Muo     ",$)');write(*,10)(mgra);
          write(*,'("Direction",$)');write(*,10)(sk)
          print*,"The value of the function";write(*,10)(fx)
          print*,"Finished one circulation "
          write(*,*)
     enddo
endprogram

subroutine G(x,y,c,s,He,A)
    use dat1
    implicit none
    real*8::C,A
    real*8,dimension (1:2)::x,y,s
    real*8,dimension(1:2,1:2)::He
    y(1)=2*x(1)
    y(2)=8*x(2)
    c=(y(1)**2+y(2)**2)**(0.5)
    S=-1.0/c*(y)
    He(1,1)=2
    He(1,2)=0
    He(2,1)=0
    He(2,2)=8
    A=-1.0*(dot_product(y,S)/(dot_product(Matmul(S,He),S)))
    return
endsubroutine

