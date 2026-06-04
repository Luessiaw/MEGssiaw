import numpy as np

unit_x = np.array([1,0,0])
unit_y = np.array([0,1,0])
unit_z = np.array([0,0,1])
origin = np.array([0,0,0])

k_mu = 1e-7

def fibonacci_sphere(M=100):
    '''产生单位球面上均匀分布的 M 个点
    返回 (cos phi sin theta, sin phi sin theta, cos theta)'''
    gold_ratio = (np.sqrt(5)+1)/2
    
    points = []
    for i in range(M):
        a = i/gold_ratio
        a = a - int(a)
        phi = 2*np.pi*a
        z = 1 - (2*i+1)/M
        theta = np.arccos(z)
        p = (np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),z)
        p = np.array(p)
        # p = np.reshape(p,(1,3))
        points.append(p)

    return points


def fibonacci_half_sphere(M=100):
    '''产生单位上半球面上均匀分布的 M 个点
    返回 (cos phi sin theta, sin phi sin theta, cos theta)'''
    gold_ratio = (np.sqrt(5)+1)/2
    
    points = []
    for i in range(M):
        a = i/gold_ratio
        a = a - int(a)
        phi = 2*np.pi*a
        z = 1 - (2*i+1)/M/2
        theta = np.arccos(z)
        p = (np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),z)
        p = np.array(p)
        # p = np.reshape(p,(1,3))
        points.append(p)

    return points

def getSphericalUnitVector(r):
    '''在r处的三个单位向量'''
    e1 = r/np.linalg.norm(r)
    e2 = np.cross(unit_z,e1)
    if np.linalg.norm(e2) > 1e-5:
        e2 = e2/np.linalg.norm(e2)
        e3 = np.cross(e1,e2)
        e3 = e3/np.linalg.norm(e3)
    else:
        e2 = unit_x
        e3 = unit_y
    return e1,e2,e3