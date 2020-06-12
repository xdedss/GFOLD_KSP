# -*- coding: utf-8 -*-
# GFOLD_static_p3p4


min_=min
from cvxpy import *
import cvxpy_codegen as cpg
from time import time
import numpy as np
import sys

import GFOLD_params as params

''' As defined in the paper...

 PROBLEM 3: Minimum Landing Error (tf roughly solved)
 MINIMIZE : norm of landing error vector
 SUBJ TO  :
            0) initial conditions satisfied (position, velocity)
            1) final conditions satisfied (altitude, velocity)
            2) dynamics always satisfied
            3) x stays in cone at all times
            4) relaxed convexified mass and thrust constraints
            5) thrust pointing constraint
            6) sub-surface flight constraint

 PROBLEM 4: Minimum Fuel Use
 MAXIMIZE : landing mass, opt variables are dynamical and
 SUBJ TO  :
            0) same constraints as p1, plus:
            1) landing point must be equal or better than that found by p1

'''


def GFOLD_gen(N, pmark, fname): # PRIMARY GFOLD SOLVER

    if pmark=='p3':
        program = 3
    elif pmark=='p4':
        program = 4

    g0 = 9.80665
    
    x0 = Parameter(6, 1, name="x0")
    z0_term_inv = Parameter(1, N, name="z0_term_inv")
    z0_term_log = Parameter(1, N, name="z0_term_log")
    g = Parameter(3, 1, name="g_vec")
    sparse_params = Parameter(10, 1, name="sparse_params")
    alpha_dt, G_max, V_max, y_gs_cot, p_cs_cos, m_wet_log, r1, r2, tf_, straight_fac = sparse_params
    
    dt = tf_ * (1/N)  # Integration dt
    
    print('N = ',N)

    x =Variable(6,N,name='var_x') # state vector (3position,3velocity)
    u =Variable(3,N,name='var_u') # u = Tc/mass because Tc[:,n]/m[n] is not allowed by DCP
    z= Variable(1,N,name='var_z')  # z = ln(mass)
    s= Variable(1,N,name='var_s') # thrust slack parameter
    

    con = []  # CONSTRAINTS LIST
    
    con += [x[0:3:1,0]  == x0[0:3]] # initial pos
    con += [x[3:6,0]  == x0[3:6]] # initial vel
    con += [x[3:6,N-1]== np.array([0,0,0])] # don't forget to slow down, buddy!

    con += [s[0,N-1] == 0] # thrust at the end must be zero
    con += [u[:,0] == s[0,0]*np.array([1,0,0])] # thrust direction starts straight
    con += [u[:,N-1] == s[0,N-1]*np.array([1,0,0])] # and ends straight
    con += [z[0,0] == m_wet_log] # convexified (7)

    if program==3:
        con += [x[0,N-1] == 0]

    elif program==4:
        con += [x[0:3,N-1] == np.array([0,0,0])] # force landing point equal to O

    for n in range(0,N-1):

        con += [x[3:6,n+1] == x[3:6,n] + (dt*0.5)*((u[:,n]+g[:,0]) + (u[:,n+1]+g[:,0]))]
        con += [x[0:3,n+1] == x[0:3,n] + (dt*0.5)*(x[3:6,n+1]+x[3:6,n])]

        # glideslope cone
        con += [ norm( (x[0:3,n])[1:3] ) - y_gs_cot*(x[0,n])  <= 0 ]
            
        con += [ norm(x[3:6,n]) <= V_max ] # velocity
        #con += [norm(u[:,n+1]-u[:,n]) <= dt*T_max/m_dry * 3]
        con += [z[0,n+1] == z[0,n] - (alpha_dt*0.5)*(s[0,n] + s[0,n+1])] # mass decreases
        con += [norm(u[:,n]) <= s[0,n]] # limit thrust magnitude & also therefore, mass

        # Thrust pointing constraint
        con += [ u[0,n] >= p_cs_cos*s[0,n]  ]

        if n > 0:
            
            #z0_term = m_wet - alpha * r2 * (n) * dt  # see ref [2], eq 34,35,36
            
            #z0 = log(z0_term)
            
            z0 = z0_term_log[0,n]
            
            mu_1 = r1*(z0_term_inv[0,n])
            mu_2 = r2*(z0_term_inv[0,n])
            
            #更正一处原项目与论文不符之处
            # 示意图：https://www.desmos.com/calculator/wtcfgnepe1
            con += [s[0,n] >= mu_1 * (1 - (z[0,n] - z0) + (z[0,n] - z0)**2 *0.5)] # lower thrust bound
            con += [s[0,n] <= mu_2 * (1 - (z[0,n] - z0))] # upper thrust bound


    #con += [x[0,0:N-1] >= 0] # no
    
    if program == 3:
        print('-----------------------------')
        #objective=Minimize(norm(x[0:3,N-1]-rf))
        expression = 0
#        for i in range(N):
#            expression += norm(x[4:6,i])*(1) # - rf[0:3,0]
#        expression *= straight_fac
        #增加一些正则项，希望i越大，r越接近于原点
        #这样导致p3生成的轨迹会在尽量短的时间内到达目标并一直停留在目标
        #这样可以求出到达目标需要的最短时间
        for i in range(N):
            expression += norm(x[0:3,i])*(i/N) # - rf[0:3,0]
        objective=Minimize(expression)
        problem=Problem(objective,con)
        print('generating p3')
        cpg.codegen(problem, fname)
        #obj_opt=problem.solve(solver=ECOS,verbose=True,feastol=5e-20)#solver=SCS,max_iters=5000,verbose=True,use_indirect=False)
        #print(x.value)
        print('-----------------------------')
        #print(z.value)
    elif program == 4:
        print('-----------------------------')
        #objective=Maximize(z[0,N-1])
        expression = 0
        #增加一些正则项，希望在接近目标时优先消除水平速度
        for i in range(N):
            expression += norm(x[4:6,i])*(i/N) # - rf[0:3,0]
        expression *= straight_fac
        expression += -z[0,N-1] * N
        objective=Minimize(expression)
        problem=Problem(objective,con)
        print('generating p4')
        cpg.codegen(problem, fname)
        #obj_opt=problem.solve(solver=ECOS,verbose=True)#solver=SCS,max_iters=5000,verbose=True,use_indirect=False,warm_start=True) # OK to warm start b/c p1 gave us a decent answer probably
        print('-----------------------------')

if __name__ == '__main__':
    
    if 'p3' in sys.argv[1:]:
        print ('--------------------------------P3------------')
        GFOLD_gen(params.N3, 'p3', 'GFOLD_p3')
        
    if 'p4' in sys.argv[1:]:
        print ('--------------------------------P4------------')
        GFOLD_gen(params.N4, 'p4', 'GFOLD_p4')
    
    
