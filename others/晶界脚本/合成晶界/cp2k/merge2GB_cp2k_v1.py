"""
使用方法
0. 更改 26 行的 OUTFILE，确定脚本输出的文件名， 更改 93 行的 def（rmOverlap(arra,arrb,const,distance=1): 更改判定两原子重合的最小距离

1. 更改 main('a.cp2k', 'b.cp2k', switch='center')
    输出一系列文件夹，文件夹名字对应中心晶界的距离
    其中 a.vasp 和 b.vasp 对应需要合并的两个晶体文件名
2. 得到中心晶界的最佳距离 z1 时进行下一步
    更改 main('a.cp2k', 'b.cp2k', x, switch='border')
    输出一系列文件夹，文件夹名字对应边界晶界的距离
3. 得到 b.vasp 向前移动的最佳距离 y 时进行下一步
    更改 main('a.cp2k', 'b.cp2k', x, y, switch='forward')
    输出一系列文件夹，文件夹名字对应沿着b轴正方向移动的距离
注意:
1. 坐标文件需要时笛卡尔坐标
2. 脚本中定义的上下是沿着a轴上下，前后是沿着b轴前后。上和前都是这一轴的正方向
"""

import numpy as np
import os
from math import sqrt,pow
import shutil

OUTFILE = 'subsys.inc'

def getLocation(posFile):
    """
    定位 cp2k 输入文件的 cell 坐标 和 原子坐标 的起始和结束的行数
    返回 (胞坐标起始行, 胞坐标结束行, 原子坐标起始行，原子坐标结束行)
    注意以 0 行作为第一行
    """
    with open(posFile,'r') as f:
        head = f.readlines()
        head_copy = head[:] # bak

        begin_line = None # &subsys 所在行，以0为第一行, 下面相同
        end_line = None  # &end subsys 所在行
        cell_begin = None # 胞坐标开始的行数
        cell_end = None # 胞坐标结束的行数
        coord_begin = None # 胞坐标开始的行数
        coord_end = None # 胞坐标结束的行数

        for i in range(len(head)):
            if i > 0:
                last_content = head[i-1].strip()
            else:
                last_content = ""
            content = head[i].strip()

            #print(last_content[1:].lower())
            if content.upper() == '&SUBSYS':
                begin_line = i

            if content and content[0] != '&' and last_content and last_content[1:].upper() == 'CELL':
                cell_begin = i
                cell_end = i + 2

            if content and content[0] != '&' and last_content and last_content[1:].upper() == 'COORD':
                coord_begin = i

            if coord_begin and not coord_end and content and content[0] == '&':
                coord_end = i-1

            if 'SUBSYS' in content.upper():
                end_line = i
        return cell_begin,cell_end,coord_begin,coord_end
        

def cp2k2vasp(arr):
    """
    将 cp2k 格式的原子坐标转形式转换成 vasp 格式
    如 [[N 1 2 3],[N 4 5 6]] --> [[1 2 3],[4 5 6]]
    返回值 (元素矩阵，坐标矩阵)
    """
    arr_element, arr = arr.copy()[:,0][:,np.newaxis], arr.copy()[:,[1,2,3]]
    arr = arr.astype('f4')
    #print(arr_element)
    return arr_element,arr

def vasp2cp2k(arr_element,arr):
    """
    将 vasp 格式的原子坐标转形式转换成 cp2k 格式
    如 [[1 2 3],[4 5 6]] --> [[N 1 2 3],[N 4 5 6]]
    要求 arr_element.shape = [3,1] 而不是 [3,]
    """
    arr_element = arr_element.copy()
    arr = arr.copy()
    arr_new = np.hstack([arr_element,arr])
    return arr_new

def rmOverlap(arra,arrb,const,distance=1):
    """
    移除arrb上与arra重合的原子，distance定义的是当两个原子距离小于这个值时认为它们重合
    注意，arra 与 arrb 需要是笛卡尔坐标，arrb 是 元组 arrb = (arrb_element,arrb)
    const 指的是胞的晶格常数，用来判断周期性
    N 原子半径是 0.80，N-N 的距离是1.45
    """
    repeatIndex=[]
    arrb_element,arrb = arrb[0].copy(), arrb[1].copy()
    arra = arra.copy()
    
    a = const[0,0]
    b = const[1,1]
    c = const[2,2]
    
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
        arrb_element = np.delete(arrb_element,repeatIndex,axis=0)
    return arrb_element,arrb
    
def retPoscar(inp, const, arra, arrb):
    """
    输出指定的 cp2k 坐标文件
    inp 是任意一个晶体的 POSCAR 路径，比如 '.\a.vasp'
    const 指的是晶体原始的晶格常数
    arra 指的是第一个晶体中元素及其坐标，arrb 指的是第二个晶体的元素及其坐标。
    """
    const_new = const.copy()
    #const_new[2,2] = const[2,2] * 2 + c_dist_center + c_dist_border
    
    #atomsNum = 2 * arra.shape[0] - (arra.shape[0] - arrb.shape[0])
    #print(atomsNum)
    
    with open(inp,'r') as f:
        file = f.readlines()
        cell_begin, cell_end, coord_begin, coord_end = getLocation(inp)

        file_copy = file[:] # bak
        
        head = file[:coord_begin]
        head[cell_begin] = 'A    ' + '    '.join(str(const_new[0].tolist()).strip('[]').split(',')) + '\n'
        head[cell_begin+1] = 'B    ' + '    '.join(str(const_new[1].tolist()).strip('[]').split(',')) + '\n'
        head[cell_begin+2] = 'C    ' + '    '.join(str(const_new[2].tolist()).strip('[]').split(',')) + '\n'
        
        with open(OUTFILE,'w') as fw:
            fw.writelines(head)
            np.savetxt(fw,arra,fmt='%s')
            np.savetxt(fw,arrb,fmt='%s')
            fw.writelines(file[coord_end+1:])

def main(inpa, inpb, centermove=0, bordermove=0, step=0.2, maxM=1.6, maxM_left=None, switch='center'):
    """
    upmove 是定下第二个晶体向上移动的距离
    centermove 是定下中心晶界的距离
    bordermove 是定下边界晶界的距离
    step 是调整距离时每次移动的步长
    maxM 是向右移动的最大距离
    maxM_left 指的是向左搜索的最大范围，可以不赋值，默认为 maxM 的一半
    """
    
    cell_begin_a, cell_end_a, coord_begin_a, coord_end_a = getLocation(inpa)
    cell_begin_b, cell_end_b, coord_begin_b, coord_end_b = getLocation(inpb)
    
    arra_ori = np.loadtxt(inpa,dtype='O', skiprows=coord_begin_a, max_rows=coord_end_a-coord_begin_a+1)
    arrb_ori = np.loadtxt(inpb,dtype='O', skiprows=coord_begin_b, max_rows=coord_end_b-coord_begin_b+1)
    
    arra_element_ori,arra_ori = cp2k2vasp(arra_ori)
    arrb_element_ori,arrb_ori = cp2k2vasp(arrb_ori)
    #print(len(arrb_ori))
    #print(arrb_ori)
    #print(arra_ori)
    
    const_ori = np.loadtxt(inpa, dtype='O', skiprows=cell_begin_a, max_rows=3)[:,1:4].astype('f8')
    a_ori = const_ori[0,0]  # a轴长度
    b_ori = const_ori[1,1]  # b轴长度
    c_ori = const_ori[2,2]   # 原始的c轴长度
    #print(const)

    if switch == 'center':
        
        # 如果没有对 maxM_left 做赋值的化，默认取 maxM 的 -1/2
        # maxM_left 是负值
        if maxM_left is None:
            maxM_left = -float('%.2f' % (maxM/2)) 
 
        for j in np.arange(maxM_left, maxM, step):
            const = const_ori.copy()
            const[2,2] = 2*c_ori+j
            
            arrb_element = arrb_element_ori.copy()
            arrb = arrb_ori.copy()
            
            arrb[:,2] = arrb[:,2] + c_ori + j
            
            arrb_element, arrb = rmOverlap(arra_ori,(arrb_element,arrb),const)

            arra_vasp = vasp2cp2k(arra_element_ori,arra_ori)
            arrb_vasp = vasp2cp2k(arrb_element,arrb)
            retPoscar(inpa, const, arra_vasp, arrb_vasp)
            
            dirname = f'{j:.1f}'
            os.mkdir(dirname)
            shutil.move(OUTFILE,dirname+os.sep+OUTFILE)
    
    elif switch == 'border':
        
        if maxM_left is None:
            maxM_left = -float('%.2f' % (maxM/2))
            
        for k in np.arange(maxM_left, maxM, step):
            const = const_ori.copy()
            const[2,2] = 2*c_ori+centermove+k
            
            arrb_element = arrb_element_ori.copy()
            arrb = arrb_ori.copy()
            
            arrb[:,2] = arrb[:,2] + c_ori + centermove

            arrb_element, arrb = rmOverlap(arra_ori,(arrb_element,arrb),const)

            arra_vasp = vasp2cp2k(arra_element_ori,arra_ori)
            arrb_vasp = vasp2cp2k(arrb_element,arrb)
            retPoscar(inpa, const, arra_vasp, arrb_vasp)
            
            dirname = f'{k:.1f}'
            os.mkdir(dirname)
            shutil.move(OUTFILE,dirname+os.sep+OUTFILE)
            
            
    elif switch == 'forward':

        for i1 in np.arange(0, b_ori, step):
            const = const_ori.copy()
            const[2,2] = 2*c_ori + centermove + bordermove
            
            arrb_element = arrb_element_ori.copy()
            arrb = arrb_ori.copy()
            
            arrb[:,1] += i1
            arrb[:,2] = arrb[:,2] + c_ori + centermove
            
            arrb_element, arrb = rmOverlap(arra_ori,(arrb_element,arrb),const)

            arra_vasp = vasp2cp2k(arra_element_ori,arra_ori)
            arrb_vasp = vasp2cp2k(arrb_element,arrb)
            retPoscar(inpa, const, arra_vasp, arrb_vasp)
            
            dirname = f'{i1:.1f}'
            os.mkdir(dirname)
            shutil.move(OUTFILE,dirname+os.sep+OUTFILE)
            
if __name__ == '__main__':
    main('a.cp2k','b.cp2k',10, 10, switch='forward')
