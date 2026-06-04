from MEGssiaw import *

'''固定位置的源和单个探头通道，计算测量值'''

# 探头通道。
# 考虑探头位于头顶位置
rs = np.array([0,0,0.11])
# 1. 矢量磁力仪, 三轴
rss_v = np.array([rs,]*3)
# 2. 标量磁力仪, 测量方向沿特定方向。
rss_s = rs[None,:]

# 探头通道的方向
# 矢量磁力仪
nss_v = np.array([unit_x,unit_y,unit_z])
# 标量磁力仪
theta = np.deg2rad(30)
unit = unit_x*np.cos(theta) + unit_y*np.sin(theta)
nss_s = unit[None,:]

# 源通道, 用于产生信号
rp = np.array([0,0,0.08])
rps = np.array([rp,rp]) # 两个通道，因为源有两个方向
nps = np.array([unit_x,unit_y])
# 源强度
Q = np.array([0, 10e-9]) # 单位:Am

# 导联场, 用于产生信号
# 矢量磁力仪
L0_v = computeLeadFieldMatrix(rps,nps,rss_v,nss_v)
# 标量磁力仪
L0_s = computeLeadFieldMatrix(rps,nps,rss_s,nss_s)

# 信号
B_v = L0_v @ Q
B_s = L0_s @ Q

print(B_v)
print(B_s)

print("Done")