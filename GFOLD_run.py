
import numpy as np
import time
import sys

import GFOLD_params as params

class solver:
    def __init__(self, v_data=None):
        if v_data != None:
            self.set_params(v_data)
    
    def set_params(self, v_data):
        self.Isp_inv = 1 / v_data['Isp']
        self.alpha = 1 / 9.80665 / v_data['Isp']
        self.G_max = v_data['G_max']
        self.V_max = v_data['V_max']
        self.y_gs_cot = 1 / np.tan(v_data['y_gs'])
        self.p_cs_cos = np.cos(v_data['p_cs'])
        self.m_wet = v_data['m_wet']
        self.m_wet_log = np.log(v_data['m_wet'])
        self.r1 = v_data['T_max'] * v_data['throt'][0]
        self.r2 = v_data['T_max'] * v_data['throt'][1]
        self.tf_ = v_data['tf']
        self.straight_fac = v_data['straight_fac']
        
        self.x0 = v_data['x0']
        self.g = v_data['g']
    
    def pack_data(self, N):
        dt = self.tf_ / N
        alpha_dt = self.alpha * dt
        t = np.linspace(0, (N-1) * dt, N)
        z0_term = self.m_wet - self.alpha * self.r2 * t
        z0_term_inv = (1 / z0_term).reshape(1, N)
        z0_term_log = np.log(z0_term).reshape(1, N)
        x0 = self.x0.reshape(6,1)
        g = self.g.reshape(3,1)
        sparse_params = np.array((alpha_dt, self.G_max, self.V_max, self.y_gs_cot, self.p_cs_cos, self.m_wet_log, self.r1, self.r2, self.tf_, self.straight_fac))
        sparse_params=sparse_params.reshape(len(sparse_params),1)
        return (x0,z0_term_inv,z0_term_log,g,sparse_params)
    
    def run_p3(self):
        import gfold_solver_p3 as solver3
        (x0,z0_term_inv,z0_term_log,g,sparse_params) = self.pack_data(params.N3)
        res = solver3.cg_solve(x0 = x0, g_vec = g, z0_term_log = z0_term_log, z0_term_inv = z0_term_inv, 
            sparse_params=sparse_params)
        
        if res[1]['status'] == 'optimal':
            tf_m = self.tf_
            x = res[0]['var_x']
            for i in range(x.shape[1]):
                if (np.linalg.norm(x[0:3,i]) + np.linalg.norm(x[3:6,i])) < 0.1:
                    tf_m = i / x.shape[1] * self.tf_
                    break
            return tf_m
        else:
            print(res)
            return self.tf_ #None
        
    def run_p4(self):
        import gfold_solver_p4 as solver4
        (x0,z0_term_inv,z0_term_log,g,sparse_params) = self.pack_data(params.N4)
        res = solver4.cg_solve(x0 = x0, g_vec = g, z0_term_log = z0_term_log, z0_term_inv = z0_term_inv, 
            sparse_params=sparse_params)
        
        if res[1]['status'] == 'optimal':
            m = np.exp(res[0]['var_z'])
            return (self.tf_,res[0]['var_x'],res[0]['var_u'],m,res[0]['var_s'],res[0]['var_z'])
            # (tf,x,u,m,s,z)
        else:
            print(res)
            m = np.exp(res[0]['var_z'])
            return (self.tf_,res[0]['var_x'],res[0]['var_u'],m,res[0]['var_s'],res[0]['var_z'])#None

    def solve(self):
        print("------solve_generated-------")
        start = time.time()
        
        tf_m = self.run_p3()
        if tf_m == None:
            print('p3 failed')
            return None
        self.tf_ = tf_m + 0.1 * self.straight_fac
#        tf_m = self.run_p3()
#        if tf_m == None:
#            print('p3- failed')
#            return None
#        self.tf_ = tf_m
        print('tf_m:' + str(tf_m))
        res = self.run_p4()
        if res == None:
            print('p4 failed')
            return None
        print("------solved in %fs-------" % (time.time() - start))
        return res

    def solve_direct(self, N3 = params.N3, N4 = params.N4):
        print("------solve_direct-------")
        import GFOLD_direct_exec as solver_direct
        start = time.time()
        packed_data = self.pack_data(N3)
        (obj_opt,x,u,m,s,z) = solver_direct.GFOLD_direct(N3, 'p3', packed_data)
        if obj_opt == None:
            print('p3 failed')
            return None
        tf_m = self.tf_
        for i in range(x.shape[1]):
            if (np.linalg.norm(x[0:3,i]) + np.linalg.norm(x[3:6,i])) < 0.1:
                tf_m = i / x.shape[1] * self.tf_
                break
        print('tf_m:' + str(tf_m))
        self.tf_ = tf_m + 0.1 * self.straight_fac
        packed_data = self.pack_data(N4)
        (obj_opt,x,u,m,s,z) = solver_direct.GFOLD_direct(N4, 'p4', packed_data)
        if obj_opt == None:
            print('p4 failed')
            return None
        print("------solved in %fs-------" % (time.time() - start))
        return (tf_m,x,u,m,s,z)
        

if __name__ == '__main__':
    from EvilPlotting import *
#    test_vessel = {
#        'Isp' : 203.94,
#        'G_max' : 3,
#        'V_max' : 90,
#        'y_gs' : np.radians(30),
#        'p_cs' : np.radians(45),
#        'm_wet' : (2)*1e3 + (0.3)*1e3,
#        'T_max' : 24000,
#        'throt' : [0.2,0.8],
#        #'x0' : np.array([2400, 450, -330,  -10,  -40,   10]),
#        'x0' : np.array([1400, 450, -330,  -20,  40,   40]),
#        'g' : np.array([-3.71,0,0]),
#        'tf' : 80,
#        'straight_fac' : 1,
#    }
    test_vessel = {
        'Isp': 256.83880615234375, 
        'G_max': 100, 
        'V_max': 150, 
        'y_gs': 0.7853981633974483, 
        'p_cs': 1.335176877775662, 
        'm_wet': 5904.27880859375, 
        'T_max': 172557.234375, 
        'throt': [0.1, 0.8], 
        'x0': np.array([  230.6261789 ,  -259.6486983 ,   107.25770515, -1424.00991313, -153.36160242,    25.85659298]), 
        'g': np.array([-9.807,  0.   ,  0.   ]), 
        'tf': 30, 
        'straight_fac': 5
    }
    if 'direct' in sys.argv[1:]:
        print('solving test vessel directly')
        (tf,x,u,m,s,z) = solver(test_vessel).solve_direct()
    else:
        print('solving test vessel using generated code')
        (tf,x,u,m,s,z) = solver(test_vessel).solve()
    plot_run3D(tf,x,u,m,s,z, test_vessel)


