! 黄金分割搜索，实验函数fx=2*sin(x)-(x**2)/10
module dat
implicit none
real*8::x1,x2,xu,xl,d,dc,fx1,fx2
endmodule

program main
    use dat
implicit none
integer::i
print *,"Welcome !"
print *,"Input the range of the domain!"
print*,"Left"
read *,xl
print*,"Right"
read*,xu
    dc=(5**0.5-1)/2
        d=dc*(xu-xl)
        x1=xl+d
        x2=xu-d
        call cal(x1,fx1)
        call cal(x2,fx2)
do i=1,8,1
      print *,"The result is ",x1,x2
    x1=xl+d
        x2=xu-d
        call cal(x1,fx1)
        call cal(x2,fx2)
     if (fx1>fx2) then
        xl=x2
        x2=x1
        d=dc*(xu-xl)
        x1=xu-d
    else
        xu=x1
        x1=x2
        d=dc*(xu-xl)
        x2=xl+d
        endif

enddo
print *,"The result is ",fx1,fx2
        print *,"x1,x2",x1,x2
endprogram

subroutine cal(x,fx)
    use dat
    implicit none
    real*8::x,fx
    fx=2*sin(x)-(x**2)/10
    return
    endsubroutine

