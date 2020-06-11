

import krpc
import time
import math
import numpy as np
import numpy.linalg as npl

from GFOLD_run import solver

from threading import Thread

print('solver imported')

params = {}
with open('params.txt', 'r', encoding='utf-8') as f:
    for line in f:
        pair = line.split('#')[0].split('=')
        if len(pair) == 2:
            key = pair[0].strip()
            value = eval(pair[1])
            params[key] = value

def clamp(num, maxnum, minnum):
    if num > maxnum:
        return maxnum
    elif num < minnum:
        return minnum
    return num

def clamp_mag(vec, maxmag):
    mag = npl.norm(vec)
    if mag > maxmag:
        return vec / mag * maxmag
    return vec

def lerp(vec1, vec2, t):
    return t * vec2 + (1-t) * vec1

def sgn(f):
    if f > 0:
        return 1
    elif f < 0:
        return -1
    return 0

def v2(f1, f2):
    return np.array((f1, f2))

def v3(f1, f2, f3):
    return np.array((f1, f2, f3))

def v(vec):
    return np.array(vec)

def normalize(vec):
    return vec / npl.norm(vec)

def q(axis, angle):
    (x, y, z) = axis
    s = math.sin(angle / 2)
    c = math.cos(angle / 2)
    axis = v3(x+.0, y+.0, z+.0)
    axis /= npl.norm(axis)
    return (s * axis[0], s * axis[1], s * axis[2], c)

def rotation_mat(q):
    (x, y, z, w) = q
    return np.mat([
    [1-2*y**2-2*z**2, 2*x*y+2*w*z, 2*x*z-2*w*y],
    [2*x*y-2*w*z, 1-2*x**2-2*z**2, 2*y*z+2*w*x],
    [2*x*z+2*w*y, 2*y*z-2*w*x, 1-2*x**2-2*y**2]
    ])

def transform(vec, mat):
    res = np.mat(vec) * mat
    return v3(res[0,0], res[0,1], res[0,2])

def angle_around_axis(v1, v2, axis):
    axis = normalize(axis)
    v1 = normalize(np.cross(v1, axis))
    v2 = normalize(np.cross(v2, axis))
    direction = sgn(np.dot(np.cross(v1, v2), axis))
    return direction * math.acos(np.dot(v1, v2))

class FreeInertialControl:
    # 假设：几乎无干扰
    def __init__(self):
        self.unified_authority = 1
        self.redundancy = 0.1
        self.kp = 1
        self.kd = 0
    
    def update(self, error, velocity):
        acc = self.unified_authority * (1 - self.redundancy)
        target_vel = math.sqrt(2 * max(abs(error) - 0.5 * deg2rad, 0) * acc) * sgn(-error)
        self.result = (target_vel - velocity) * self.kp - velocity * self.kd
        return self.result

class PID:
    def __init__(self):
        self.ep = True
        self.ei = True
        self.ed = True
        self.kp = 1
        self.ki = 0
        self.kd = 1
        self.sd = 0
        self.diff = 0
        self.integral = 0
        self.integral_limit = 1
        self.error_prev = 0
        self.first = True
        self.second = True
        self.dumpf = None
    
    def update(self, error, dt):
        if self.first:
            self.first = False
            self.error_prev = error
        elif self.second:
            self.second = False
            self.diff = (error - self.error_prev) / dt
        
        self.integral += error * dt * self.ki
        self.integral = clamp(self.integral, self.integral_limit, -self.integral_limit)
        self.diff = lerp(self.diff, (error - self.error_prev) / dt, 1-self.sd)
        p = -error * self.kp
        i = -self.integral
        d = -self.diff * self.kd
        self.result = p * (1 if self.ep else 0) + i * (1 if self.ei else 0) + d * (1 if self.ed else 0)
        
        self.p_prev = p
        self.i_prev = i
        self.d_prev = d
        self.error_prev = error
        return self.result

deg2rad = math.pi / 180
rad2deg = 180 / math.pi
g0 = params['g0']

conn = krpc.connect(name='PID')
space_center = conn.space_center
vessel = space_center.active_vessel
flight = vessel.flight()
body = vessel.orbit.body

delta_time = 0.01

#target
target_lat = params['target_lat'] * deg2rad
target_lon = params['target_lon'] * deg2rad
target_height = params['target_height']
target_axis = target_height + body.surface_height(target_lat * rad2deg, target_lon * rad2deg) + body.equatorial_radius
target_body_pos = np.array((math.cos(target_lon) * math.cos(target_lat), math.sin(target_lat), math.sin(target_lon) * math.cos(target_lat))) * target_axis


#limit
max_tilt = params['max_tilt'] * deg2rad
max_tilt_sin = math.sin(max_tilt)
min_tilt_cos = math.cos(max_tilt)
max_tilt_tan = math.tan(max_tilt)

throttle_limit = params['throttle_limit']
throttle_limit_ctrl = params['throttle_limit_ctrl']

#生成求解器所需的输入
def vessel_profile1(vessel_d):
    est_time = 0.5
    #print('vel ' + str(vessel_d['vel']))
    #print('acceleration ' + str(vessel_d['acceleration']))
    #print('error ' + str(vessel_d['error']))
    vel_est = vessel_d['vel'] + vessel_d['acceleration'] * est_time
    pos_est = vessel_d['error'] + vessel_d['vel'] * est_time + 0.5 * vessel_d['acceleration'] * est_time * est_time
    return {
        'Isp' : vessel_d['specific_impulse'],
        'G_max' : params['G_max'],
        'V_max' : params['V_max'],
        'y_gs' : params['y_gs'] * deg2rad,
        'p_cs' : max_tilt * 0.85,
        'm_wet' : vessel_d['mass'],
        'T_max' : vessel_d['max_thrust'],
        'throt' : throttle_limit,
        'x0' : np.array([pos_est[0], pos_est[1], pos_est[2], vel_est[0], vel_est[1], vel_est[2]]),
        'g' : np.array([-g0,0,0]),
        'tf' : params['tf'],
        'straight_fac' : params['straight_fac'],
    }

#rotation 
ctrl_x_rot = PID()
ctrl_x_rot.kp = params['ctrl_x_rot.kp']
ctrl_x_rot.kd = params['ctrl_x_rot.kd']
#ctrl_x_rot.redundancy = 0.1
ctrl_y_avel_kp = params['ctrl_y_avel_kp']
ctrl_z_rot = PID()
ctrl_z_rot.kp = params['ctrl_z_rot.kp']
ctrl_z_rot.kd = params['ctrl_z_rot.kd']
#ctrl_z_rot.redundancy = 0.1
#测量值
#torque = v3(3.66e+04, 5000, 3.66e+04)
#torque_k = v3(8.2e+04-3.66e+04, 0, 8.2e+04-3.66e+04)
#torque = v3(10300.000011920929, 10300.000011920929, 10300.000011920929)
#torque_k = v3(15183.20083618,  10772.2761631,   15183.24184418)
#print(vessel.available_torque)
# k
k_x = params['k_x']
k_v = params['k_v']

# final
final_throttle = params['final_throttle']
final_kp = params['final_kp']

# time
game_delta_time = 0.02
game_prev_time = space_center.ut
start_time = time.time()

# references
ref_local = vessel.reference_frame 
ref_surface = vessel.surface_reference_frame #地面参考系
ref_body = body.reference_frame
ref_target_temp = space_center.ReferenceFrame.create_relative(ref_body, position=target_body_pos)
ref_target = space_center.ReferenceFrame.create_hybrid(ref_target_temp, rotation=ref_surface, velocity=ref_target_temp)

prev_vel = v(vessel.velocity(ref_target))
N = 80
gfold_path = None
n_i = -1
error = v(vessel.position(ref_target))
debug_lines = params['debug_lines']
if debug_lines:
    lines = [conn.drawing.add_line((0,0,0),(0,0,0), ref_target) for i in range(N-1)]
    directions = [conn.drawing.add_line((0,0,0), (1,0,0), ref_target) for i in range(N)]
    target_line = conn.drawing.add_line((0,0,0),(1,0,0),ref_target)
    target_line.color = (0,0,1)
    target2_line = conn.drawing.add_line((0,0,0),(1,0,0),ref_target)
    target2_line.color = (0,0,1)
    head_line = conn.drawing.add_line((0,0,0),(1,0,0),ref_target)
    head_line.color = (0,1,1)
    for line in directions:
        line.color = (0,1,0)

nav_mode = 'none'
frcount = 0

def update_lines(x, u):
    m_u = vessel.max_thrust / vessel.mass
    for i in range(N-1):
        lines[i].start = x[0:3, i]
        lines[i].end = x[0:3, i+1]
    for i in range(N):
        directions[i].start = x[0:3, i]
        directions[i].end = x[0:3, i] + u[:, i] * 5 / m_u

#找到路径上最近的点的索引值
def find_nearest_index(x, r, v):
    nearest_mag = npl.norm(x[0:3, 0] - r)
    nearest_i = 0
    for i in range(x.shape[1]):
        mag = npl.norm(x[0:3, i] - r) # + npl.norm(x[3:6, i] - v) * 0.2
        if mag < nearest_mag:
            nearest_mag = mag
            nearest_i = i
    v = x[3:6, nearest_i]
    v_norm = npl.norm(v)
    v_dir = v / v_norm
    frac = clamp(np.dot(r - x[0:3, nearest_i], v_dir) / (tf / N * v_norm), 0.5, -0.5)
    return nearest_i + frac

#路径上采样
def sample_index(index):
    #if index >= N-1:
    if index >= N-1:
        #return (x[0:3, N-1], x[3:6, N-1], u[:, N-1])
        return (v3(0,0,0), v3(0,0,0), v3( 9.807,0,0))
    elif index <= 0:
        i = 0
        frac = index
    else:
        i = math.floor(index)
        frac = index - i
    x_i_s = lerp(x[:, i], x[:, i+1], frac)
    u_i_s = lerp(u[:, i], u[:, i+1], frac)
    if index < 0:
        u_i_s = u[:, 1].copy()
        #print('u1234 ' + str(u[:, 0:4]))
        #print('u_i_s ' + str(u[:, 1]))
    return (x_i_s[0:3].copy(), x_i_s[3:6].copy(), u_i_s.copy())

def conic_clamp(target_a, min_mag, max_mag, max_tilt):
    a_mag = npl.norm(target_a)
    hor_dir = v3(0, target_a[1], target_a[2])
    hor_dir /= npl.norm(hor_dir)
    #target_direction = target_a / a_mag
    a_hor = npl.norm(target_a[1:3])
    a_ver = target_a[0]
    
    if (a_hor < min_mag * math.sin(max_tilt)):
        a_ver_min = math.sqrt(min_mag**2 - a_hor**2)
    else:
        a_ver_min = math.cos(max_tilt) * min_mag
    
    if (a_hor < max_mag * math.sin(max_tilt)):
        a_ver_max = math.sqrt(max_mag**2 - a_hor**2)
    else:
        a_ver_max = math.cos(max_tilt) * max_mag
    
    a_ver = clamp(a_ver, a_ver_max, a_ver_min)
    
    a_hor = min(a_hor, a_ver * math.tan(max_tilt))
    
    return hor_dir * a_hor + v3(a_ver, 0, 0)

#求解最优路径
def solve_gfold(v_data):
    global gfold_path, n_i, nav_mode
    if params['direct']:
        gfold_path = solver(v_data).solve_direct
    else:
        gfold_path = solver(v_data).solve()
    n_i = -100
    if gfold_path != None:
        (tf,x,u,m,s,z) = gfold_path
        if debug_lines:
            update_lines(x, u)
        #print('x0: ' + str(v_data['x0']))
        #print('tf: ' + str(tf))
        print('gfold')
        nav_mode = 'gfold'
        print(v_data)

while True:
    time.sleep(delta_time)
    real_time = time.time() - start_time
    ut = space_center.ut
    game_delta_time = ut - game_prev_time
    if game_delta_time < 0.01: #意味着游戏中还没有经过一个物理帧，所以不进行计算
        continue
        
# 取得一些之后要用的数据
    vessel_d = {}
    error = vessel_d['error'] = v(vessel.position(ref_target)) # 目标系里的偏差
    avel = vessel_d['avel'] = v(vessel.angular_velocity(ref_surface)) # 地面系下角速度（等于目标系角速度
    vel = vessel_d['vel'] = v(vessel.velocity(ref_target)) # 地面速度
    
    rotation_local2srf = rotation_mat(v(vessel.rotation(ref_surface))) # 机体系到地面系旋转矩阵
    rotation_srf2local = npl.inv(rotation_local2srf) # 地面系到机体系旋转矩阵
    moment_of_inertia_local = v(vessel.moment_of_inertia) # 转动惯量
    mass = vessel_d['mass'] = vessel.mass
    max_thrust = vessel_d['max_thrust'] = vessel.max_thrust
    acceleration = vessel_d['acceleration'] = (vel - prev_vel) / game_delta_time
    
    vessel_d['specific_impulse'] = vessel.specific_impulse
    #print(game_delta_time)
    
    
        
    if nav_mode == 'gfold': #跟随gfold路径
        (tf,x,u,m,s,z) = gfold_path
        #n_i = max(n_i + game_delta_time * 0.2 * N/tf, find_nearest_index(x, error, vel))
        n_i = max(n_i - game_delta_time * 0.2 * N/tf, find_nearest_index(x, error, vel))
        print(n_i)
#        x_i = (x[0:3, n_i] + x[0:3, min(n_i+1 ,N-1)]) / 2
#        v_i = (x[3:6, n_i] + x[3:6, min(n_i+1 ,N-1)]) / 2
#        u_i = (u[:, n_i] + u[:, min(n_i+1 ,N-1)]) / 2
        (x_i, v_i, u_i) = sample_index(n_i)
        (x_i_, v_i_, u_i_) = sample_index(n_i + min(1.5 * N/tf, npl.norm(vel) / 50 * N/tf))
        
        #v_i_dir = v_i / npl.norm(v_i)
#        target_a = u_i_
#        target_v = v_i_
#        target_x = x_i # + np.dot((error - x_i), v_i_dir) * v_i_dir
#        #print(n_i, target_a, target_v, target_x)
#        target_a += (target_v - vel) * k_v + (target_x - error) * k_x
        
        target_a = u_i + (v_i - vel) * k_v + (x_i - error) * k_x
        target_a_ = u_i_ + (v_i_ - vel) * k_v + (x_i - error) * k_x
        
        if debug_lines:
            target_line.start = error
            target_line.end = (x_i[0], x_i[1], x_i[2])
            target2_line.start = error
            target2_line.end = (x_i_[0], x_i_[1], x_i_[2])
        
        max_throttle_ctrl = throttle_limit_ctrl[1] * (max_thrust / mass)
        min_throttle_ctrl = throttle_limit_ctrl[0] * (max_thrust / mass)
        target_a = conic_clamp(target_a, min_throttle_ctrl, max_throttle_ctrl, max_tilt)
        target_a_ = conic_clamp(target_a_, min_throttle_ctrl, max_throttle_ctrl, max_tilt)
        if n_i < 0 :
            target_a = np.array([g0, 0, 0]) + u_i
            
#        target_direction[0] = max(target_direction[0], (1.0/max_tilt_tan) * npl.norm(target_direction[1:3]))
#        #target_direction[1:3] *= (1/npl.norm(target_direction[1:3])) * math.sqrt(1-target_direction[0]**2)
        
        target_direction = target_a_ / npl.norm(target_a_)
        target_throttle = npl.norm(target_a) / (max_thrust / mass)
        #print(target_a)
        if debug_lines:
            head_line.start = error
            head_line.end = error + target_direction * 8
        
        if n_i > 0:
            vessel.control.throttle = target_throttle
        
        if (N - n_i) * tf / N < 10:
            vessel.control.gear = True
        if npl.norm(error[1:3]) < params['final_radius'] and npl.norm(error[0]) < params['final_height']:
            vessel.control.gear = True
            print('final')
            nav_mode = 'final'
    elif nav_mode == 'final': #最终降落阶段，直接pid
        max_acc = throttle_limit_ctrl[1] * (max_thrust / mass) - g0
        max_acc_low = throttle_limit_ctrl[1] * final_throttle * (max_thrust / mass) - g0
        est_h = error[0] - vel[0]**2 / (2 * max_acc)
        est_h_low = error[0] - vel[0]**2 / (2 * max_acc_low)
        est_h_center = (est_h + est_h_low) / 2
        
        vessel.control.throttle = clamp(lerp(throttle_limit_ctrl[1] * final_throttle, throttle_limit_ctrl[1], -est_h_low / (est_h - est_h_low) * (1+final_kp)), throttle_limit_ctrl[1], throttle_limit_ctrl[0])
        
        error_hor = v3(0, error[1], error[2])
        vel_hor = v3(0, vel[1], vel[2])
        ctrl_hor = -error_hor * 0.03 - vel_hor * 0.06
        target_direction = ctrl_hor + v3(1, 0, 0)
        target_direction /= npl.norm(target_direction)
        target_direction = conic_clamp(target_direction, 1, 1, max_tilt)
    elif nav_mode == 'simple':
        target_a = (-vel) * k_v*3 + (-error) * k_x*0 + v3(g0, 0, 0)
        a_mag = npl.norm(target_a)
        target_throttle = a_mag / (max_thrust / mass)
        target_direction = target_a / a_mag
        target_direction[0] = max(target_direction[0], min_tilt_cos)
        target_direction[1:3] *= (1/npl.norm(target_direction[1:3])) * math.sqrt(1-target_direction[0]**2)
        vessel.control.throttle = target_throttle
    else:
        target_direction = -vel
        vessel.control.throttle = 0
        
    if error[0] < params['start_altitude']:
        #frcount -= 1
        if frcount <= 0:
            frcount = 1
#            v_data = {
#                'Isp' : vessel.specific_impulse,
#                'G_max' : 5,
#                'V_max' : 180,
#                'y_gs' : np.radians(30),
#                'p_cs' : max_tilt,
#                'm_wet' : mass,
#                'T_max' : max_thrust,
#                'throt' : [0.1, 1.0],
#                'x0' : np.array([pos_est[0], pos_est[1], pos_est[2], vel_est[0], vel_est[1], vel_est[2]]),
#                'g' : np.array([-g0,0,0]),
#                'tf' : 20,
#            }
            v_data = vessel_profile1(vessel_d)
            Thread(target=solve_gfold, args=(v_data,)).start()
        
# 变换到机体坐标系计算姿态控制，以下xyz均指机体系
    target_direction_local = transform(target_direction, rotation_srf2local) # 机体系的目标姿态的机体y轴指向
    avel_local = transform(avel, rotation_srf2local)  # 机体系角速度
    # 三个轴方向能提供的最大角加速度
    #authority_local = (torque + torque_k * vessel.control.throttle) / moment_of_inertia_local 
    #authority_local = np.abs(v(vessel.available_torque[0])) / moment_of_inertia_local 
    #ctrl_x_rot.unified_authority = authority_local[0]
    #ctrl_z_rot.unified_authority = authority_local[2]
    # pid控制，roll直接消除角速度
    #control_pitch = -clamp(ctrl_x_rot.update(angle_around_axis(target_direction_local, v3(0, 1, 0), v3(1, 0, 0)), avel_local[0]), 1, -1)
    #control_yaw = -clamp(ctrl_z_rot.update(angle_around_axis(target_direction_local, v3(0, 1, 0), v3(0, 0, 1)), avel_local[2]), 1, -1)
    control_pitch = -clamp(ctrl_x_rot.update(angle_around_axis(target_direction_local, v3(0, 1, 0), v3(1, 0, 0)), game_delta_time), 1, -1)
    control_yaw = -clamp(ctrl_z_rot.update(angle_around_axis(target_direction_local, v3(0, 1, 0), v3(0, 0, 1)), game_delta_time), 1, -1)
    control_roll = clamp(avel_local[1] * ctrl_y_avel_kp, 1, -1)
    
    vessel.control.pitch = control_pitch
    vessel.control.yaw = control_yaw
    vessel.control.roll = control_roll
    
# 终止条件
    if (npl.norm(error[1:3]) < 3 and npl.norm(error[0]) < 1 and npl.norm(vel[1:3]) < 0.3 and npl.norm(vel[0]) < 0.5 and npl.norm(avel) < 0.2) or (vel[0] > 0 and npl.norm(error[0]) < 1):
        vessel.control.throttle = 0
        break
    
    prev_vel = vel
    game_prev_time = ut