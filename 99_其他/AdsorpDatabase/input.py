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

    init_data(atoms_tem=atoms_tem, adsorb_element="H", adsorb_pos_list=adsorb_pos_list, db_path=db_path)
    scf(vasp_scf, db_path, num_list=[1,100,100,100,100,100,100,100,100]) # num_list为吸附0，1，2...个结构的计算数量
    # sr(vasp_relax, db_path, num_list=[0,2,2,2,2])
    # scf_id(vasp_scf, db_path, id=1)
    # sr_id(vasp_relax, db_path, id=2)
    # freeze_id(db_path, 2)
    # unfreeze_id(db_path, 2)
    # view_freeze_stru(db_path)
    # view_ori_stru(db_path)
    # view_scf_stru(db_path, sort_way="adsorb_e", single_adsorb_element_energy=-3.38794825) # -3.38794825为vasp计算的H2的单个H的能量
    # view_relax_stru(db_path)
    # view_strus(db_path, id_list=[3,12,13,1])
    # view_stru("id=12")
    # write_id(db_path, id=[7], format="vasp")
    # rmscf(db_path, id=1)
    # rmrelax(db_path, id=1)
    endTime = time.time()
    runtime = endTime - startTime
    logging.info("=============== --- ==============")
    logging.info("End of calculation.")
    logging.info("Program was running for %d min %d s." % (runtime/60, runtime%60))
    logging.info("=============== end ==============")