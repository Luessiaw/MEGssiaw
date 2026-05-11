from MEGssiaw import *
import matplotlib.pyplot as plt
import matplotlib as mpl

'''使用 L Curve 方法，搜索最佳参数'''

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
N = 128
sensor_radius = 0.11
rss = np.array(fibonacci_half_sphere(N))*sensor_radius
nss = np.array([unit_z,]*N)

# 用于计算测量值的导联场
L0 = computeLeadFieldMatrix(rp0s,np0s,rss,nss)
B0 = L0 @ Q0
B = B0 + np.random.randn()*100e-15 # 添加噪声

L = computeLeadFieldMatrix(rps,nps,rss,nss)
params_exp = np.arange(-5,1,0.5)
params = np.array([10**float(p) for p in params_exp])

normQ, normResidual = computeLCurve(L,B,params)


print("L Curve Result:")
print("lambda:   "+" \t".join([f"{p:.2e}" for p in params]))
print("norm Q:   "+" \t".join([f"{p:.2e}" for p in normQ]))
print("norm Res: "+" \t".join([f"{p:.2e}" for p in normResidual]))


fig, axs = plt.subplots(1,2,figsize=(10,4),constrained_layout=True)

norm = mpl.colors.Normalize(
    vmin=np.min(params_exp),
    vmax=np.max(params_exp)
)
cmap = plt.cm.viridis
sc1 = axs[0].scatter(
    params, normQ,
    c=params_exp,
    cmap=cmap,
    norm=norm,
    s=20
)
axs[0].scatter(
    params, normResidual,
    c=params_exp,
    cmap=cmap,
    norm=norm,
    s=20,
    marker='x'
)

axs[0].set_xlabel("lambda")
axs[0].set_ylabel("Norm")
axs[0].set_title("Norm v.s. lambda")
axs[0].grid()
axs[0].legend(["|Q|","|LQ-B|"])

sc2 = axs[1].scatter(
    normResidual, normQ,
    c=params_exp,
    cmap=cmap,
    norm=norm,
    s=20
)
axs[1].set_xlabel("|LQ-B|")
axs[1].set_ylabel("|Q|")
axs[1].set_title("solution v.s. residual")
axs[1].grid()

cbar = fig.colorbar(
    sc2,          # 任意一个 scatter 即可
    ax=axs,       # 作用于所有子图
    location='right'
)
cbar.set_label("lg(lambda)")

for ax in axs:
    ax.set_xscale('log')
    ax.set_yscale('log')

fig.savefig("examples/example_4_L_curve.png",dpi=300)
print("Done.")
