# Sample runs of ODE solvers
import numpy as np
"""RK2"""
#import RK_methods as RK
#x,y=RK.RK(1,[0,20],0.02,lambda x, y: y*np.sin(x),2,0.5,0.5,1,1,True)
"""BetaMethod"""
import beta_method as bm
x,y=bm.BetaMethod(lambda x,y: y*(1-y),lambda x,y: 1-2*y,[0,20],0.1,0.5,0.2,True)