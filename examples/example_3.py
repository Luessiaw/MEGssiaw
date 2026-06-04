'''固定源，由测量值进行源成像'''
from MEGssiaw import *

# 真实源, 用于产生信号
head_radius = 0.08
rp0 = np.array([0,0,head_radius])
rp0s = np.array([rp0,rp0]) # 两个方向
np0s = np.array([unit_x,unit_y])
# 源强度
Q0 = np.array([0,100e-9]) # 单位:Am

# 用于成像的格点。在半径为 8 cm 的上半球面上均匀划分。
grid_length = 0.005
grid_point_num = int(8*head_radius**2/grid_length**2)
print("Number of grid points: {0:d}".format(grid_point_num))

# 源格点通道
rps = np.array(fibonacci_half_sphere(grid_point_num))*head_radius
nps = np.zeros((grid_point_num*2,3))
for i,rp in enumerate(rps):
    e1,e2,e3 = getSphericalUnitVector(rp)
    nps[i,:] = e2
    nps[i+grid_point_num,:] = e3
rps = np.concatenate([rps,rps])

# 探头通道，标量磁力仪，测量方向沿z轴
N = 256
sensor_radius = 0.11
rss = np.array(fibonacci_half_sphere(N))*sensor_radius
nss = np.array([unit_z,]*N)

# 用于计算测量值的导联场
L0 = computeLeadFieldMatrix(rp0s,np0s,rss,nss)
B0 = L0 @ Q0
B0 += np.random.randn()*100e-15 # 添加噪声

# 用于源定位的导联场
L = computeLeadFieldMatrix(rps,nps,rss,nss)
W = computeInverseOperator(L,1e-5)
Q = W @ B0

rps = rps[:grid_point_num]
powers = Q[:grid_point_num]**2 + Q[grid_point_num:]**2
locRes = computeLocRes(rps,powers,0.5)

print("True Source Location: {0} cm".format(rp0*1e2))
print("Localized Source: {0} cm".format(locRes*1e2))
print("Localization Error: {0:.2f} mm.".format(np.linalg.norm(rp0-locRes)*1e3),)

print()