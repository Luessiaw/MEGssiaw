from MEGssiaw import *
from MEGssiaw.visualise import plt

'''固定位置的源，在探头球面上产生的测量值分布，并绘制散点图'''

# 源通道, 用于产生信号
rp = np.array([0,0,0.08])
rps = rp[None,:] # 只考虑一个方向
nps = unit_y[None,:]
# 源强度
Q = np.array([100e-9,]) # 单位:Am

fig,axs = plt.subplots(1,4,sharex=True,sharey=True,figsize=(20,5),constrained_layout=True )
vmax = 4
vmin = -4

# 探头通道。
N = 10000
rss = np.array(fibonacci_sphere(N*2)[:N])*0.11
# 1. 测量方向沿径向
nss = []
for rs in rss:
    ns = rs/np.linalg.norm(rs)
    nss.append(ns)
nss = np.array(nss)
for i,unit in enumerate(["radial",unit_x,unit_y,unit_z]):
    if i:
        # 2. 测量方向沿 x,y,z 轴
        nss = np.array([unit,]*N)
    L = computeLeadFieldMatrix(rps,nps,rss,nss)
    B = L @ Q
    print("Max B: {0:.3f}\t Min B: {1:.3f}".format(np.max(B)*1e12,np.min(B)*1e12))

    ax = axs[i]
    sc = ax.scatter(rss[:,0]*1e2,rss[:,1]*1e2,s=20,c=B*1e12,cmap="viridis",vmin=vmin,vmax=vmax)
    ax.set_xlabel("x (cm)")
    ax.set_xlim([-12,12])
    ax.set_ylim([-12,12])
    ax.set_aspect("equal")
    ax.set_title("B_{0:s}".format(["r","x","y","z"][i]))

axs[0].set_ylabel("y (cm)")
cbar = fig.colorbar(
    sc,          # 任意一个 scatter 返回值即可
    ax=axs,      # 绑定所有子图
    location="right",
    shrink=0.9
)
cbar.set_label("Measured Value (pT)")

fig.savefig("examples/figs/example_1-MeasuredValue.png",dpi=300)
print("Done")
