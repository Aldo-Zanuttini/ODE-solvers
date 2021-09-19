import matplotlib.pyplot as plt
import numpy as np
# FUNCTION BetaMethod
# Solves the IVP y'=f(x,y) using the beta method
# INPUTS
# y0               initial condition (y(0)=y0)
# interval         interval to solve upon (a vector [a,b])
# h                step size
# f                RHS of differential equation y'=f(x,y), lambda function! So write 
# df               derivative of f with respect to y
# beta             defines the method (0 for forward Euler, 1 for backward Euler, 0.5 for Crank-Nicolson)
# plot             optional: if ==True then the function plots results, default==False
# OUTPUTS
# y_vals           list containing solution for y'=f(x,y)
# x_vals           list containing the gridpoints
def BetaMethod(f,df,interval,y0,beta,h,plot=False):
    y=y0
    x=interval[0]
    N=int((interval[1]-interval[0])/h)
    x_vals=[]
    y_vals=[]
    x_vals.append(x)
    y_vals.append(y)
    if type(y0)==float or type(y0)==int:
        Jacobian_dimensions=1
    else:
        Jacobian_dimensions=len(y0)
    if beta==0:
        for i in range(1,N):
            y=y+h*f(x,y)
            x=x+h
            y_vals.append(y)
            x_vals.append(x)
    else:
        for i in range(1,N):
            F=lambda u: u-y-h*f(x+beta*h,(1-beta)*y+beta*u)
            dF=lambda u: np.identity(Jacobian_dimensions)-h*beta*np.array(df(x+beta*h,(1-beta)*y+beta*u))
            error=1
            y_temp=y
            while error<0.001:
                y_temp_new=y_temp-F(y_temp)/dF(y_temp)
                error=np.absolute(y_temp_new-y_temp)
                y_temp=y_temp_new
            y=y+h*f(x+beta*h,y*(1-beta)+y_temp*beta)
            x=x+h
            y_vals.append(y)
            x_vals.append(x)
    if plot==True:
        plt.plot(x_vals,y_vals)
        plt.xlabel("x")
        plt.ylabel("y")
    return x_vals,y_vals