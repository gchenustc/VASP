from adsorpFun import *
from vasp_set import *
import time

if __name__ == "__main__":
    startTime = time.time()
    # 模板 - 拓展结构的固定方式和模板结构一样
    atoms_tem = read("100b-8vac-stable.vasp")
    # 吸附的位置
    adsorb_pos_list = [(7.135891290,5.028322562,10.794521724),
                       (7.135891290,1.225822666,10.794521724), 
                       (3.333391054,5.028322562,10.794521724), 
                       (3.333391054,1.225822666,10.794521724),
                       (5.542600611,7.411263812,11.532419419),
                       (5.542600611,3.608764029,11.532419419),
                       (1.740100602,7.411263812,11.532419419),
                       (1.740100602,3.608764029,11.532419419)]
    # 数据库路径
    db_path = "100b-8vac-H.db"

    # init_data(atoms_tem=atoms_tem, adsorb_element="H", adsorb_pos_list=adsorb_pos_list, db_path=db_path)
    # scf(vasp_scf, db_path, num_list=[1,10,10,10,10,10,10,10,1], adsorb_info={"H":1}) # num_list为吸附0，1，2...个结构的计算数量
    # sr(vasp_relax, db_path, num_list=[0,0,0,0,0,0,0,0,0], adsorb_info={"H":1}, unconvergence_selected=False, only_unconvergence=False)
    # scf_id(vasp_scf, db_path, id=1)
    # sr_id(vasp_relax, db_path, id=1)
    # freeze_id(db_path, id=1)
    # unfreeze_id(db_path, id=1)
    # update_label_decomposed(db_path, id_list=[1,2,3], value_list=[True,True,True]) # 将id=1，2和3的结构的decomposed标签更改为True，其他id的该标签值不变
    # view_freeze_stru(db_path)
    # view_decomposed_stru(db_path)   # 查看已经结构优化的结构中分解的结构，首先要使用update_label_decomposed函数更新一下
    # view_undecomposed_stru(db_path)   # 查看已经结构优化的结构中的未分解的结构
    # view_ori_stru(db_path)
    # view_scf_strus(db_path, energy_sort_way="adsorb_e", single_adsorb_element_energy=-3.38794825) # -3.38794825为vasp计算的H2的单个H的能量，不包含未收敛的结构
    # view_relax_strus(db_path, energy_sort_way="total", single_adsorb_element_energy=-3.38794825) # 不包含未收敛的结构
    # view_relax_stru_total(db_path) # 查看所有的不包含冻结的其他所有结构优化的结构
    # view_strus(db_path, "H=2")
    # view_strus_id(db_path, id_list=[1,2])
    # rmscf(db_path, id=1)
    # rmrelax(db_path, id=87, rmrange="current") # rmrange: "all" or "current"
    # write_id(db_path, id=[1,2], format="vasp")

    endTime = time.time()
    runtime = endTime - startTime
    logging.info("=============== --- ==============")
    logging.info("End of calculation.")
    logging.info("Program was running for %d min %d s." % (runtime/60, runtime%60))
    logging.info("=============== end ==============\n")