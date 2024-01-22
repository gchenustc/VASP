"""
使用方法
1. 更改 main('a.vasp', 'b.vasp', switch='up')
    输出一系列文件夹，文件夹名字对应第二个晶体 ('b.vasp') 向上(a轴)移动的距离
    其中 a.vasp 和 b.vasp 对应需要合并的两个晶体文件名
2. 得到 b.vasp 向上移动的最佳距离 x 时进行下一步
    更改 main('a.vasp', 'b.vasp', x, switch='forward')
    输出一系列文件夹，文件夹名字对应第二个晶体 ('b.vasp') 向前(b轴)移动的距离
3. 得到 b.vasp 向前移动的最佳距离 y 时进行下一步
    更改 main('a.vasp', 'b.vasp', x, y, switch='center')
    输出一系列文件夹，文件夹名字对应中心晶界的距离
4. 得到中心晶界的最佳距离 z1 时进行下一步
    更改 main('a.vasp', 'b.vasp', x, y, z1, switch='border')
    输出一系列文件夹，文件夹名字对应边界晶界的距离
注意:
1. 脚本只适用于单种原子的晶体
2. a.vasp 和 b.vasp 是 Cartesian 坐标
3. 脚本中定义的上下是沿着a轴上下，前后是沿着b轴前后。上和前都是这一轴的正方向
"""

import numpy as np
import os
from math import sqrt,pow
import shutil

def rmOverlap(arra,arrb,const,distance=1.0):
    """
    移除arrb上与arra重合的原子，distance定义的是当两个原子距离小于这个值时认为它们重合
    注意，arra 与 arrb 需要是笛卡尔坐标
    N 原子半径是 0.80，N-N 的距离是1.45
    """
    repeatIndex=[]
    
    a = const[0,0]
    b = const[1,1]
    c = const[2,2]
    #print(c)
    arra = arra.copy()
    arrb = arrb.copy()
    for arr in [arra,arrb]:
        for row in range(len(arr)):
            while arr[row,0] >= a:
                arr[row,0] -= a
            while arr[row,0] < 0:
                arr[row,0] += a
            while arr[row,1] >= b:
                arr[row,1] -= b
            while arr[row,1] < 0:
                arr[row,1] += b
            while arr[row,2] >= c:
                arr[row,2] -= c
            while arr[row,2] < 0:
                arr[row,2] += c

    # 考虑周期性边界的重合
    arra_1 = arra.copy()
    arra_1[:,0] -= a
    arra_2 = arra.copy()
    arra_2[:,0] += a
    arra_3 = arra.copy()
    arra_3[:,1] -= b
    arra_4 = arra.copy()
    arra_4[:,1] += b
    
    arra_5 = arra.copy()
    arra_5[:,2] += c
    
    arra_6 = arra_5.copy()
    arra_6[:,0] -= a
    arra_7 = arra_5.copy()
    arra_7[:,0] += a
    arra_8 = arra_5.copy()
    arra_8[:,1] -= b
    arra_9 = arra_5.copy()
    arra_9[:,1] += b
    #print(arrb[:3],arrb_[:3])

    
    for i in range(arrb.shape[0]):  
        for j in range(arra.shape[0]):
            for arra_part in [arra,arra_1,arra_2,arra_3,arra_4,arra_5,arra_6,arra_7,arra_8,arra_9]:
                dist = sqrt(np.sum(np.square((arrb[i] - arra_part[j]))))
                #dist1 = sqrt(np.sum(np.square((arrb_1[i] - arra[j]))))
           
                if dist < distance:
                    #print(dist,dist1)
                    #print(arra[i],arrb[i],arrb_[i])                
                    repeatIndex.append(i)
                    break
                    
    repeatIndex = np.unique(np.array(repeatIndex))
    if repeatIndex.tolist():
        arrb = np.delete(arrb,repeatIndex,axis=0)  # 删除对应行，返回的是一个拷贝
    return arrb
    
    
def retPoscar(poscar, const, arra, arrb):
    """
    输出指定的 POSCAR
    poscar 是任意一个晶体的 POSCAR 路径，比如 '.\a.vasp'
    arra 指的是第一个晶体中原子的坐标，arrb 指的是第二个晶体中原子的坐标。
    """
    
    #const = np.loadtxt(poscar,dtype=np.float64,skiprows=2,max_rows=3)
    
    const_new = const.copy()
    #const_new[2,2] = const[2,2] * 2 + c_dist_center + c_dist_border
    
    atomsNum = 2 * arra.shape[0] - (arra.shape[0] - arrb.shape[0])
    #print(atomsNum)
    
    with open(poscar,'r') as f:
        head = f.readlines()
        #atomsNum_merge = 2 * int(head[6]) - atoms_dec
        head_copy = head[:] # bak
        head = head[:8]
        head[2] = '    '.join(str(const_new[0].tolist()).strip('[]').split(',')) + '\n'
        head[3] = '    '.join(str(const_new[1].tolist()).strip('[]').split(',')) + '\n'
        head[4] = '    '.join(str(const_new[2].tolist()).strip('[]').split(',')) + '\n'
        head[6] = str(atomsNum) + '\n'
        with open('POSCAR','w') as fw:
            fw.writelines(head)
            np.savetxt(fw,arra,fmt='%.8f')
            np.savetxt(fw,arrb,fmt='%.8f')
     
     
def main(poscara, poscarb, upmove=0, forwardmove=0, centermove=0, bordermove=0, step=0.2, maxM=1.6, maxM_left=None, switch='up'):
    """
    upmove 是定下第二个晶体向上移动的距离
    centermove 是定下中心晶界的距离
    bordermove 是定下边界晶界的距离
    step 是调整距离时每次移动的步长
    maxM 是向右移动的最大距离
    maxM_left 指的是向左搜索的最大范围，可以不赋值，默认为 maxM 的一半
    """
    
    # 确保中心和边界的晶界不同时移动
    assert not centermove or not bordermove
    
    arra_ori = np.loadtxt(poscara,dtype=np.float64,skiprows=8)
    arrb_ori = np.loadtxt(poscarb,dtype=np.float64,skiprows=8)
    
    const_ori = np.loadtxt(poscara,dtype=np.float64,skiprows=2,max_rows=3)
    a_ori = const_ori[0,0]  # a轴长度
    b_ori = const_ori[1,1]  # b轴长度
    c_ori = const_ori[2,2]   # 原始的c轴长度
    
    if switch == 'up':
        # j,k = 0,0 # j, k 分别是中心晶界和周期性边界晶界的距离
        for i in np.arange(0, a_ori, step):
            const = const_ori.copy()
            const[2,2] = 2*c_ori
            
            arrb = arrb_ori.copy()
            
            arrb[:,0] += i
            arrb[:,2] += c_ori
            arrb = rmOverlap(arra_ori,arrb,const)
            
            retPoscar(poscara,const,arra_ori,arrb)
            
            dirname = f'{i:.1f}'
#             print(dirname)
            os.mkdir(dirname)
            shutil.move('POSCAR',dirname+os.sep+'POSCAR')
            
    elif switch == 'forward':
        arrb_move = arrb_ori.copy()
        arrb_move[:,0] += upmove

        for i1 in np.arange(0, b_ori, step):
            const = const_ori.copy()
            const[2,2] = 2*c_ori
            
            arrb = arrb_move.copy()
            
            arrb[:,1] += i1
            arrb[:,2] += c_ori
            
            arrb = rmOverlap(arra_ori,arrb,const)
            
            retPoscar(poscara,const,arra_ori,arrb)
            
            dirname = f'{i1:.1f}'
            os.mkdir(dirname)
            shutil.move('POSCAR',dirname+os.sep+'POSCAR')

    elif switch == 'center':
        arrb_move = arrb_ori.copy()
        arrb_move[:,0] += upmove
        arrb_move[:,1] += forwardmove
        # arrb_upmove = rmOverlap(arra_ori,arrb_upmove)
        
        # 如果没有对 maxM_left 做赋值的化，默认取 maxM 的 -1/2
        # maxM_left 是负值
        if maxM_left is None:
            maxM_left = -float('%.2f' % (maxM/2)) 
 
        for j in np.arange(maxM_left, maxM, step):
            const = const_ori.copy()
            const[2,2] = 2*c_ori+j
            
            arrb = arrb_move.copy()
            arrb[:,2] = arrb[:,2] + c_ori + j
            
            arrb = rmOverlap(arra_ori,arrb,const)

            retPoscar(poscara, const, arra_ori, arrb)
            
            dirname = f'{j:.1f}'
            os.mkdir(dirname)
            shutil.move('POSCAR',dirname+os.sep+'POSCAR')
    
    elif switch == 'border':
        arrb_move = arrb_ori.copy()
        arrb_move[:,0] += upmove
        arrb_move[:,1] += forwardmove
        #arrb_upmove = rmOverlap(arra_ori,arrb_upmove)
        
        if maxM_left is None:
            maxM_left = -float('%.2f' % (maxM/2))
            
        for k in np.arange(maxM_left, maxM, step):
            const = const_ori.copy()
            const[2,2] = 2*c_ori+centermove+k
            
            arrb = arrb_move.copy()
            arrb[:,2] = arrb[:,2] + c_ori + centermove

            arrb = rmOverlap(arra_ori,arrb,const)
            
            retPoscar(poscara, const, arra_ori, arrb)
            
            dirname = f'{k:.1f}'
            os.mkdir(dirname)
            shutil.move('POSCAR',dirname+os.sep+'POSCAR')

if __name__ == "__main__":
    main('a.vasp', 'b.vasp', 0, 0, switch='center', maxM=0.9,maxM_left=-0.8)
