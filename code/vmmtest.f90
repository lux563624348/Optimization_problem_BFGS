!加入了Hessian矩阵的一维搜索。

module dat1
implicit none
real,dimension(:,:),allocatable,save::Hess0,Hess1,BHess,Chess
real,dimension(:),allocatable,save::Gra0,Gra1,sk
real::fx,mgra
real,dimension(1:2)::pX1,pX0,deltaX
endmodule

program main
use dat1
implicit none
real::alpha
integer::i,j,k
allocate (Hess0(1:2,1:2),Hess1(1:2,1:2),BHess(1:2,1:2),CHess(1:2,1:2))
allocate(Gra0(1:2),Gra1(1:2),Sk(1:2))!Gradient

px0=[-1.2,1.0]!start

print*," The location of start point :    (",pX0(1),",",pX0(2),")"
print*,"The original value of fx",fx
call setup(px0,gra0,Hess0,alpha,Sk)
px1=px0+alpha*Sk

do i=1,1000,1
    call dfx(Px1,Gra1,mgra)
    if (mGra<0.0001) exit

    call G(pX0,pX1,gra0,Gra1,hess0,hess1)
    call sam1(px1,gra1,hess1,alpha,sk,mgra)
    px0=px1
    gra0=gra1
    hess0=hess1
    px1=px0+alpha*sk
    k=i
    10 format(2f10.5)
    fx=100*(px1(2)-px1(1)**2)**2+(1-px1(1))**2
    write (*,'("Location",$)');write(*,10)(pX1);
    write(*,'("Muo     ",$)');write(*,10)(mgra);
    print*,"The value of the function";write(*,10)(fx)
    print*,"Finished one circulation "
enddo
print*,k," step"
endprogram

subroutine sam1(x,Gx,he,A,S,mg)
!该步为最优梯度法部分
use dat1
implicit none
real::mg,A
real,dimension (1:2)::x,S,Gx
real,dimension(1:2,1:2)::He,he0
call dfx(x,Gx,mg)
S=-(Matmul(He,Gx))
call Hessian(X,He0)! 用Hessi矩阵代替
A=-(dot_product(Gx,S)/(dot_product(Matmul(S,He0),S)))
endsubroutine

subroutine G(x0,x1,Gx0,gx1,he0,he1)
!该步通过输入坐标X，计算梯度Gx，模mg, Hessian metrix (He)
use dat1
implicit none
real::A,B,C,mg
integer ::j
real,dimension (1:2)::x0,x1,Gx0,Gx1,dx,dG,Z
real,dimension(1:2,1:2)::He0,He1,HeB,HeC
call dfx(x1,gx1,mg)
Dx=x1-x0
Dg=gx1-gx0
Z=Matmul(He0,dG)
B=1.0/(dot_product(dx,dG))
C=1.0/(dot_product(z,dG))
call Xcheng(Dx,HeB)
HeB=HeB*B
call Xcheng(Z,HeC)
HeC=HeC*C
He1=He0+HeB-HeC
!10 format(2f10.4)
!write(*,'("Hess:")');Write(*,10)(He1(j,:),j=1,2)
endsubroutine

subroutine Xcheng(x,XX)
implicit none
real,dimension(1:2)::x
real,dimension(1:2,1:2)::XX
XX(1,1)=x(1)**2
XX(1,2)=x(1)*x(2)
xx(2,1)=xx(1,2)
XX(2,2)=x(2)**2
endsubroutine

subroutine dfx(x,df,m)
    implicit none
    real::m
    real,dimension(1:2)::x,df
    df(1)=-400*x(1)*(x(2)-x(1)**2)-2+2*x(1)
    df(2)=200*(x(2)-x(1)**2)
    m=(df(1)**2+df(2)**2)**(0.5)
    endsubroutine

subroutine setup(x,Gx,he,A,S)
!该步为最优梯度法部分
use dat1
implicit none
real::mg,A
real,dimension (1:2)::x,S,Gx
real,dimension(1:2,1:2)::He,He1
He(1,1)=1! Input the Hessian metrix
He(1,2)=0
He(2,1)=0
He(2,2)=1
call dfx(x,Gx,mg)
S=-(Matmul(He,Gx))
call Hessian(x,He1)! 用Hessi矩阵代替
A=-(dot_product(Gx,S)/(dot_product(Matmul(S,He1),S)))
endsubroutine

subroutine Hessian(x,he)
implicit none
real,dimension(1:2)::X
real,dimension(1:2,1:2)::he
He(1,1)=-400*x(1)+1200*x(1)**2+2! Input the Hessian metrix
He(1,2)=-400*x(1)
He(2,1)=He(1,2)
He(2,2)=200
endsubroutine

