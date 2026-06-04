from .exceptions import *
from .tools import *
from .maths import *
import logging
logger = logging.getLogger(__name__)

import numpy as np

def computeLeadFieldMatrix(sourcePositions:np.ndarray,
                           sourceOrientations:np.ndarray,
                           sensorPositions:np.ndarray,
                           sensorOrientations:np.ndarray):
    '''此函数用于构建前向模型矩阵。
    sourcePositions: (m,3) 数组，表示源空间中各源通道的位置。
    sourceOrientations: (m,3) 数组，表示各源通道的方向。
    sensorPositions: (n,3) 数组，表示探头通道空间中各探头通道的位置。
    sensorOrientations: (n,3) 数组，表示各探头通道的方向。
    注意，所有参数均以通道为单位，而非以源或探头位置。即若一个源或探头处包含多个通道，则视为多个源或探头。
    return:
    L: 前向矩阵，(n,m) 数组
    '''

    sp = sourcePositions.shape
    so = sourceOrientations.shape
    if (sp != so):
        msg = "源通道的位置与方向数组形状不相同。源通道位置的形状为 %s, 而源通道方向的形状为 %s." % (sp,so)
        logging.warning(msg)
        raise MEGssiawError(msg)
    if sp[1] != 3:
        msg = "源通道数组的第2维度长度应为3, 而当前为 %d." % sp[1]
        logging.warning(msg)
        raise MEGssiawError(msg)

    sp = sensorPositions.shape
    so = sensorOrientations.shape
    if (sp != so):
        msg = "探头通道的位置与方向数组形状不相同。探头通道位置的形状为 %s, 而探头通道方向的形状为 %s." % (sp,so)
        logging.warning(msg)
        raise MEGssiawError(msg)
    if sp[1] != 3:
        msg = "探头通道数组的第2维度长度应为3, 而当前为 %d." % sp[1]
        logging.warning(msg)
        raise MEGssiawError(msg)

    numOfSource = sourcePositions.shape[0]
    numOfSensor = sensorPositions.shape[0]
    
    nsx = sensorOrientations[:,0] # n^s_x
    nsy = sensorOrientations[:,1] # n^s_y
    nsz = sensorOrientations[:,2] # n^s_z

    rsx = sensorPositions[:,0] # r^s_x
    rsy = sensorPositions[:,1] # r^s_y
    rsz = sensorPositions[:,2] # r^s_z

    npx = sourceOrientations[:,0] # n^p_x
    npy = sourceOrientations[:,1] # n^p_y
    npz = sourceOrientations[:,2] # n^p_z

    rpx = sourcePositions[:,0] # r^p_x
    rpy = sourcePositions[:,1] # r^p_y
    rpz = sourcePositions[:,2] # r^p_z

    shape = (numOfSensor,numOfSource)
    nsx = np.broadcast_to(nsx[:,None],shape)
    nsy = np.broadcast_to(nsy[:,None],shape)
    nsz = np.broadcast_to(nsz[:,None],shape)
    
    rsx = np.broadcast_to(rsx[:,None],shape)
    rsy = np.broadcast_to(rsy[:,None],shape)
    rsz = np.broadcast_to(rsz[:,None],shape)
    sr = np.sqrt(rsx**2+rsy**2+rsz**2)
    
    npx = np.broadcast_to(npx[None:,],shape)
    npy = np.broadcast_to(npy[None:,],shape)
    npz = np.broadcast_to(npz[None:,],shape)
    
    rpx = np.broadcast_to(rpx[None:,],shape)
    rpy = np.broadcast_to(rpy[None:,],shape)
    rpz = np.broadcast_to(rpz[None:,],shape)
    
    Rx = rsx - rpx
    Ry = rsy - rpy
    Rz = rsz - rpz
    sR = np.sqrt(Rx**2+Ry**2+Rz**2)
    
    np_rp_x = npy*rpz - npz*rpy # (n^p_j\times r^p_j)_x
    np_rp_y = npz*rpx - npx*rpz # (n^p_j\times r^p_j)_y
    np_rp_z = npx*rpy - npy*rpx # (n^p_j\times r^p_j)_z
    
    rdR = rsx*Rx + rsy*Ry + rsz*Rz # r^s_i \cdot R_ij
    F = sR*(sr*sR+rdR)
    
    L1x = np_rp_x/F
    L1y = np_rp_y/F
    L1z = np_rp_z/F
    
    L2k = rsx*np_rp_x + rsy*np_rp_y + rsz*np_rp_z
    L2k /= F**2
    L2x = -((sR**2/sr + sR)*rsx + (2*sr + sR + rdR/sR)*Rx)*L2k
    L2y = -((sR**2/sr + sR)*rsy + (2*sr + sR + rdR/sR)*Ry)*L2k
    L2z = -((sR**2/sr + sR)*rsz + (2*sr + sR + rdR/sR)*Rz)*L2k
    
    Lx = L1x + L2x
    Ly = L1y + L2y
    Lz = L1z + L2z
    
    L = nsx*Lx + nsy*Ly + nsz*Lz
    L *= k_mu
    
    return L

def computeInverseOperator(L:np.ndarray,regular_param=1e-5,C:np.ndarray=None,CQ:np.ndarray=None):
    '''由导联场矩阵计算逆算符。
    L: 导联场矩阵，(N,M) 数组
    C: 探头通道的协方差，(N,N) 数组
    CQ: 源通道的协方差，(M,M) 数组
    regular_param: 正则化参数
    
    return 
    W: 逆矩阵，(M,N) 数组
    '''
    N,M = L.shape
    if C is None:
        C = np.eye(N,N)
    if CQ is None:
        CQ = np.eye(M,M)

    if C.shape[0] != C.shape[1]:
        raise MEGssiawError("矩阵 C 不是方阵。当前 C 的形状：", C.shape)
    if CQ.shape[0] != CQ.shape[1]:
        raise MEGssiawError("矩阵 CQ 不是方阵。当前 CQ 的形状：", CQ.shape)
    if C.shape[0] != N:
        raise MEGssiawError("协方差矩阵 C 的大小与 L 的形状不匹配。当前 C 的维度为：{0:d}, 应为 {1:d}".format(C.shape[0],N))
    if CQ.shape[0] != M:
        raise MEGssiawError("协方差矩阵 CQ 的大小与 L 的形状不匹配。当前 CQ 的维度为：{0:d}, 应为 {1:d}".format(CQ.shape[0],M))
    
    LT = L.transpose()
    Llam = L @ CQ @ LT + regular_param * C
    Linv = np.linalg.inv(Llam)
    W = CQ @ LT @ Linv
    return W

def computeLocRes(rps:np.ndarray,powers:np.ndarray,threshold=0.5):
    '''评估定位精度及弥散度。目前只针对单个偶极子使用。'''        
    powers /= np.max(powers)
    powers[powers<threshold] = 0
    powers = powers/np.sum(powers)
    locPos = np.einsum("ij,i->j",rps,powers)

    return locPos

def computeLCurve(L:np.ndarray,B:np.ndarray,params:np.ndarray=None,C:np.ndarray=None,CQ:np.ndarray=None):
    '''L curve 方法搜索最佳正则化参数。
    L: 导联场矩阵，(N,M) 数组
    B: 测量值，(N,) 数组
    params: 待搜索的正则化参数，(k,) 数组
    C: 探头通道的协方差，(N,N) 数组
    CQ: 源通道的协方差，(M,M) 数组
    
    return: (normQ, normResidual)
    normQ: 与 param 对应的 Q 的模，(k,) 数组
    normResidual: 与 param 对应的残差的模，(k,) 数组
    '''
    if params is None:
        params = np.array([10**p for p in range(-5,2)])

    normQ = np.zeros(params.shape)
    normResidual = np.zeros(params.shape)
    for k,param in enumerate(params):
        W = computeInverseOperator(L,param,C,CQ)
        Q = W @ B
        residual = L @ Q - B
        normQ[k] = np.linalg.norm(Q)
        normResidual[k] = np.linalg.norm(residual)
    
    return normQ,normResidual
