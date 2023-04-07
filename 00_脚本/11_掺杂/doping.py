from pymatgen.core import  Structure as Sr
import general_fun as gf 
from ase.io import write,read
from itertools import combinations
import os
import numpy as np
import random
import time
import logging
"""
简介
随机掺杂，需要提供模板结构，并更改 “tem_stru_path” 这个变量为结构文件名，其他参数请看末尾
在“old_strus_dir_name”中的子文件夹中的结构为排除的结构，脚本生成的结构会排除与其等价的结构。可以不提供该文件夹，但此变量必须存在
比如old_strus_dir_name=old_strus
```shell
tree old_strus
```
>>>
old_strus/a/a.vasp
old_strus/a/b.vasp
old_strus/b/a.vasp
old_strus/b/b.vasp
<<<
不要将排除结构直接放在old_strus下


使用方法：
与general_fun放在同一个文件夹下，直接运行该文件
python doping.py
"""
logging.basicConfig(level=logging.INFO, filename="log.txt", filemode="a",
                            format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')

def get_doped_index_groups(tem_stru, doped_element, n_doped):
    """返回掺杂原子在tem_stru中的结构索引，例如 n_doped=3, 返回值为[(3,12,34),(1,63,1)...]，列表里的每一个元组都是一种掺杂可能性"""
    doped_elements_index = [] # 被取代元素所在的index
    for index,stru in enumerate(tem_stru):
        if stru.species_string == doped_element:
            doped_elements_index.append(index)
    random.shuffle(doped_elements_index)
    n_doped_total = len(doped_elements_index) # 被取代元素一共的数量
    assert n_doped_total >= n_doped # 被取代元素在晶胞中的数量一定小于等于要被取代的数量
    doped_index_groups = combinations(doped_elements_index, n_doped)
    return doped_index_groups

def doping(tem_stru, doping_element, doped_index):
    """获得对应的一个掺杂结构"""
    stru_doped = tem_stru.copy()
    for index in doped_index:
        stru_doped[index]=doping_element
    stru_doped.sort()
    return stru_doped

def duplicated_check(stru, strus_list):
    """相同结构检测，return False means no duplicated strus"""
    if len(strus_list) == 0:
        return False
    
    for stru_ in strus_list:
        if gf.symmetryCheck(stru_,stru):
            return True
    return False
    #return np.any(list(map(lambda x: gf.symmetryCheck(x,stru) ,strus_list)))
    
def get_old_sturs_path(old_strus_dir_name):
    """输出的结构不会与old_strus_dir_name文件夹中的一级子文件夹中的结构等价，获得old_strus_dir_name中的一级子文件夹内的结构路径，返回一个列表"""
    if not os.path.exists(old_strus_dir_name):
        return []
    old_strus_dirs_inner_list = list(os.walk(old_strus_dir_name))[0][1]
    existed_strus_path_dirs = [os.path.join(old_strus_dir_name,i) for i in old_strus_dirs_inner_list]
    existed_strus_path_list_files = [list(os.walk(i))[0][2] for i in existed_strus_path_dirs]
    old_strus_list=[]
    for dir_,files in zip(existed_strus_path_dirs,existed_strus_path_list_files):
        for file in files:
            old_strus_list.append(os.path.join(dir_,file))
    return old_strus_list

def get_all_doping_strus(tem_stru, doped_element, doping_element, n_doped, n_out, index_start, prefix, old_strus_dir_name):
    """
    获得所有掺杂结构，会剔除等价结构和old_strus_dir_name中的结构
    tem_stru：模板结构，根据此结构掺杂
    doping_element：掺杂的元素
    n_doped：掺杂的数量
    n_out: 输出的结构数量
    index_start：输出结构的文件名开始的索引号
    prefix：输出结构的文件名的前缀
    old_strus_dir_name：输出的结构不会与old_strus_dir_name文件夹中的结构等价，只识别该文件夹内的一级子文件夹
    """
    assert n_out > 0
    doped_index_groups =  get_doped_index_groups(tem_stru, doped_element, n_doped)
    #print(doped_index_groups)
    #random.shuffle(doped_index_groups)  # 打乱列表
    
    old_strus_list = get_old_sturs_path(old_strus_dir_name)
    existed_strus_list = list(map(lambda x:Sr.from_file(x),old_strus_list)) # 转换成 pymatgen 结构
    #print(existed_strus_list)

    strus_doped_list = []
    n = 0   # n是需要的结构的数量
    n_total = 0  # n_total是检测的结构的总量
    for doped_index in doped_index_groups:
        n_total += 1 
        logging.info(f"checking no.{n_total} stru.")
        stru_doped = doping(tem_stru_pymatgen, doping_element, doped_index)
        if not duplicated_check(stru_doped,existed_strus_list+strus_doped_list):
            n += 1
            logging.info("\t\t\t\t\t... This is the required structure ...")
            logging.info(f"\t\t\t\t\t... the index(from 0) respectively are {doped_index} ...")
            strus_doped_list.append(stru_doped)
            stru_doped.to(filename=f"{prefix}_{n_doped}{doping_element}_no.{n+index_start-1}.vasp",fmt="POSCAR")
            if n>=n_out:
                break
        else:
            logging.info("\t\t\t\t\t... repeated ...")
    return strus_doped_list
    


if __name__ == "__main__":
    # 开始时间
    starttime = time.time()

    # 已创建的结构的目录（新创建的结构不包含该目录的结构或者等价的结构）
    old_strus_dir_name = "old_strus"
    # 模板结构的名称
    tem_stru_path = "BPN_30GPa_3x1x2.vasp" 
    # 被掺杂元素和掺杂元素
    doped_element = "N"
    doping_element = "P"
    # 取代的数量
    n_doped = 8
    # 输出结构的数量
    n_out = 100
    # 输出结构的索引从数字几开始
    index_start=1
    # 输出结构的前缀
    prefix="BPN_30GPa_3x1x2"
    # 读取模板结构-ase,pymatgen
    tem_stru_ase = read(tem_stru_path)
    tem_stru_pymatgen = Sr.from_file(tem_stru_path)
    # 获得掺杂的结构 - 并输出
    strus_doped_list = get_all_doping_strus(tem_stru_pymatgen, doped_element, doping_element, n_doped, n_out, index_start, prefix, old_strus_dir_name)

    # 结束时间
    endtime = time.time()
    runtime = endtime-starttime
    logging.info("\nEnd of calculation.")
    logging.info("Program was running for %.2f seconds." % runtime)
