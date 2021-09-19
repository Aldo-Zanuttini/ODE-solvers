import matplotlib.pyplot as plt
# FUNCTION RK
# Solves the IVP y'=f(x,y) using a Runge-Kutta method (either RK2 or RK4)
# INPUTS
# y0               initial condition (y(0)=y0)
# interval         interval to solve upon (a vector [a,b])
# h                step size
# fun              RHS of differential equation y'=f(x,y), lambda function! So write 
# K                order of RK method desired (2 or 4)
# w1               constant in y(n+1)=y(n)+w1*k1+w2*k2
# w2               constant in y(n+1)=y(n)+w1*k1+w2*k2
# alpha            constant in k2=h*f(x(n)+alpha*h,y(n)+beta*h)
# beta             constant in k2=h*f(x(n)+alpha*h,y(n)+beta*h)
# plot             optional: if ==True then the function plots results, default==False
# OUTPUTS
# y_vals           list containing solution for y'=f(x,y)
# x_vals           list containing the gridpoints
def RK(y0,interval,h,fun,order,w1,w2,alpha,beta,plot=False):
    a=interval[0]
    b=interval[1]
    N=(b-a)/h
    N=int(N)
    y=y0
    x=a
    x_vals=[]
    x_vals.append(x)
    y_vals=[]
    y_vals.append(y)
    if(order==2):
        for i in range(1,N):
            k1=h*fun(x,y)
            k2=h*fun(x+alpha*h,y+beta*k1)
            y=y+w1*k1+w2*k2
            x=x+h
            x_vals.append(x)
            y_vals.append(y)
    elif(order==4):
        for i in range(1,N):
            k1=h*fun(x,y)
            k2=h*fun(x+h/2,y+k1/2)
            k3=h*fun(x+h/2,y+k2/2)
            k4=h*fun(x+h,y+k3)
            x=x+h
            y=y+(k1+k2*2+k3*2+k4)/6
            x_vals.append(x)
            y_vals.append(y)
    else:
        print("Desired order of approximation not available: choose order=2 or order=4.")
    if plot==True:
        plt.plot(x_vals,y_vals)
        plt.xlabel('x')
        plt.ylabel('y')
    return x_vals, y_vals