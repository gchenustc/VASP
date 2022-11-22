import itertools as it
from typing import Sequence
from typing import Union 
import logging
import os
import time
import ase
import ase.db
import numpy as np
from ase import Atom, db
from ase.build import add_adsorbate
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
from ase.io import read, write
from ase.optimize import BFGS
from ase.visualize import view
from ase.utils.structure_comparator import SymmetryEquivalenceCheck

"""
描述：
给定一个slab模型，给定吸附位和吸附物质，遍历所有的吸附情况，并将这些结构纳入到数据库中，可以有选择性地对这些结构进行单点或者结构计算
"""
# 初始化日志模块
logging.basicConfig(level=logging.INFO, filename="log.txt", filemode="a",
                    format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')


def copyfile(path_1, path_2):
    """
    将路径path_1下的所有目录及文件复制到路径path_2下
    path_1: 待复制的目录或者文件路径
    path_2: 目标路径
    """
    if os.path.isdir(path_1):  # path_1是目录

        list_1 = os.listdir(path_1)
        if not list_1:  # 复制目录，仅仅是复制空目录
            os.mkdir(path_2)
        else:
            os.mkdir(path_2)  # 先复制最上层空目录
            for i in list_1:
                path_r = os.path.join(path_1, i)  # 下层目录或文件的绝对路径
                path_w = os.path.join(path_2, i)  # 目标目录或文件的绝对路径
                if os.path.isfile(path_r):  # 是文件则直接进行复制
                    with open(path_r, 'rb') as rstream:
                        container = rstream.read()
                        with open(path_w, 'wb') as wstream:
                            wstream.write(container)
                else:  # 是目录则调用本函数
                    copyfile(path_r, path_w)

    else:  # path_1是文件
        with open(path_1, 'rb') as rstream:  # 是文件直接复制文件
            container = rstream.read()
            file_name = os.path.basename(path_1)
            path_2 = os.path.join(path_2, file_name)
            with open(path_2, 'wb') as wstream:
                wstream.write(container)


def removedir(dir):
    dir = dir.replace('\\', '/')
    if(os.path.isdir(dir)):
        for p in os.listdir(dir):
            removedir(os.path.join(dir, p))
        if(os.path.exists(dir)):
            os.rmdir(dir)
    else:
        if(os.path.exists(dir)):
            os.remove(dir)


def _duplicateCheck(db, atoms, check_range):
    """
    check_range: 传入db.select()中的参数，在其中进行重复检测
    """
    # 检测database是否为空，空则返回True
    if len(db) == 0:
        return False

    if check_range:
        atoms_to_be_check_list = list(
            map(lambda x: x.toatoms(), list(db.select(check_range))))
    else:
        atoms_to_be_check_list = list(
            map(lambda x: x.toatoms(), list(db.select())))

    return any(len(atoms.numbers) == len(dp_atoms_each.numbers) and np.all(atoms.numbers == dp_atoms_each.numbers) and np.all(np.abs(atoms.positions - dp_atoms_each.positions) < 0.001) for dp_atoms_each in atoms_to_be_check_list)


def _symmetryCheck(db, atoms, check_range=None):
    """
    check_range: 传入db.select()中的参数，在其中进行重复检测
    """
    if len(db) == 0:
        return False

    comp = SymmetryEquivalenceCheck()

    if check_range:
        atoms_to_be_check_list = list(
            map(lambda x: x.toatoms(), list(db.select(check_range))))
    else:
        atoms_to_be_check_list = list(
            map(lambda x: x.toatoms(), list(db.select())))

    return any(comp.compare(ats, atoms) for ats in atoms_to_be_check_list)


def vaspMove(db_row, mode):
    if os.path.exists(f"./calc_record/id_{db_row.id}_{mode}"):
        removedir(f"./calc_record/id_{db_row.id}_{mode}")
    copyfile("workdir", f"./calc_record/id_{db_row.id}_{mode}")
    removedir("workdir")


def init_data(atoms_tem, adsorb_element, adsorb_pos_list, db_path):
    """
    给定模板结构和吸附位，使用该函数找到所有的吸附结构
    注意，如果模板结构的POSCAR固定了，以后所有的衍生结构都按照此固定方式进行结构优化
    """
    logging.info("--------------- 开始初始化数据库 ---------------")
    n_strus = 0
    n_strus_total = list(map(lambda x: 0, adsorb_pos_list.copy()))
    n_strus_write = n_strus_total.copy()
    n_strus_repeat = n_strus_total.copy()
    mydb = db.connect(db_path)

    # 写入模板结构
    if not _symmetryCheck(mydb, atoms_tem, check_range="H=0"):
        mydb.write(atoms_tem, scf=False, relaxed=False, ori_stru_id=0)
        n_strus += 1

    for index, n in enumerate(range(len(adsorb_pos_list))):  # 1 or 2 or 3 ...
        n_adsorb = index+1
        # (pos1, pos2, pos3 ....)
        for pos_list in it.combinations(adsorb_pos_list, n+1):
            atoms = atoms_tem.copy()
            for pos in pos_list:
                atom = Atom(symbol=adsorb_element, position=pos, tag=1)
                atoms.append(atom)

            # n_adsorb如果指定，则只在吸附n_adsorb个的结构中检查
            if not _symmetryCheck(mydb, atoms, check_range=f"H={n_adsorb}"):
                mydb.write(atoms, scf=False, relaxed=False, ori_stru_id=0)
                n_strus += 1
                n_strus_write[n] += 1
            else:
                n_strus_repeat[n] += 1
            n_strus_total[n] += 1

    logging.info(f"有1,2,3...个吸附物的结构一共找到{str(n_strus_total).strip('[]')}个")
    logging.info(f"因为等价而剔除的结构分别分{str(n_strus_repeat).strip('[]')}个")
    logging.info(f"写入数据库的结构分别分{str(n_strus_write).strip('[]')}个")
    logging.info(f"一共写入{str(n_strus).strip('[]')}个数据。")
    logging.info("--------------- 结束初始化数据库 ---------------\n")


def rmscf(db_path, id):
    """
    删除init结构的scf标签，删除标签后不会抹除上一次计算的能量力等信息，但是之后允许对该结构进行新的计算，在view_scf_stru也不会对该结构进行展示
    """
    logging.info("--------------- remove scf result ---------------")
    mydb = db.connect(db_path)
    row = mydb.get(id=id)
    if row.scf:
        mydb.update(id=id, scf=False)
        logging.info("成功删除scf结果\n")
    else:
        logging.info("该id没有进行scf计算\n")


def rmrelax(db_path, id, rmrange="all"):
    """
    删除init结构的relax的计算结果
    rmrange: all or current; all为删除所有realx结果，rmrange是从当前id开始删除
    """
    assert rmrange in ["all", "current"]
    logging.info("--------------- remove relax result ---------------")

    mydb = db.connect(db_path)
    row = mydb.get(id=id)

    if _isfreeze(db_path, row):
        logging.info(f"id={id}的结构的init结构(id={_get_pre_ori_id(db_path, row)})已经被冻结，请先解冻\n")
        return

    if row.ori_stru_id == 0:
        rmrange = "all"

    if rmrange == "all":
        pre_ori_id = _get_pre_ori_id(db_path, row)
        mydb.update(id=pre_ori_id, relaxed=False)
    elif rmrange == "current":
        pre_ori_id = row.ori_stru_id
        if mydb.get(id=pre_ori_id).ori_stru_id == 0:
            mydb.update(id=pre_ori_id, relaxed=False)
    pre_ori_row = mydb.get(id=pre_ori_id)

    relax_chain_id = []
    current_row = pre_ori_row
    while 1:
        try:
            current_relax_target_id = current_row.relax_target_id
            if current_relax_target_id <= 0:
                break
            relax_chain_id.append(current_relax_target_id)
            current_row = mydb.get(id=current_relax_target_id)
        except Exception: 
            break

    if not relax_chain_id:
        logging.info(f"该结构id={id}没有进行relax计算\n")
        return

    for id in relax_chain_id:
        del mydb[id]
        logging.info(f"成功删除init结构(id={pre_ori_id})的relax结构(id={id})")
    logging.info("")

    mydb.update(id=pre_ori_id, relax_target_id=-1)


def scf(calc, db_path, num_list, adsorb_info={"H": 1}):
    """
    adsorb_kind: 如果吸附的是CO2，则 adsorb_kind = {"C":1, "O":2}
    num_list: 指定进行自洽计算地结构个数，列表索引的第一个是不吸附的结构的选取数量，第二个是吸附一个group的结构的选取数量，以此类推。
    """
    database = db.connect(db_path)
    if isinstance(calc, Vasp):
        # 将vasp计算的结果扔进calc_record保存，所以首先要创建这个文件夹
        if not (os.path.exists("calc_record") and os.path.isdir("calc_record")):
            os.mkdir("calc_record")
        calc.directory = ("workdir")
    logging.info("--------------- 开始批量 scf ---------------")
    logging.info(
        f"0,1,2,3...个吸附物结构的scf数量: {str(num_list).strip('[]')}, 一共需要计算{sum(num_list)}个结构")
    n_clacs = 0

    for ad_num, str_num in enumerate(num_list):

        if str_num == 0:
            continue

        logging.info(f">>> 开始对吸附了{ad_num}个的结构进行scf <<<")

        adsorb_kinds = adsorb_info.keys()  # ["H", ...]
        adsorb_info_str = ",".join(
            map(lambda x: f"{x}={ad_num*adsorb_info[x]}", adsorb_kinds))  # "C=1,O=2"
        # ori_stru_id=0 是选取init中的结构
        row_list = list(database.select(
            f"scf=False, ori_stru_id=0, {adsorb_info_str}"))

        if len(row_list) < str_num:
            str_num_old = str_num
            str_num = len(row_list)
            logging.info(
                f">>> 需要计算的结构有{str_num_old}个，但是数据库中最多只能选择{str_num}个 <<<")

        # 选取str_num个结构进行scf
        for row_groups in it.combinations(row_list, str_num):
            for row in row_groups:
                atoms = row.toatoms()
                atoms.calc = calc
                logging.info(f"开始进行第{n_clacs+1}个结构的计算, 该结构的id={row.id}")
                try:
                    atoms.get_potential_energy()
                    if isinstance(calc, Vasp):
                        # 计算完成后取得电子步
                        nelm_actual = int(
                            os.popen(r"grep -E 'Iteration[ ]+1' ./workdir/OUTCAR  | wc -l").read())
                        # 判断收敛
                        if nelm_actual == int(calc.parameters["nelm"]):
                            database.update(
                                id=row.id, atoms=atoms, scf=True, converg=False)
                        else:
                            database.update(
                                id=row.id, atoms=atoms, scf=True, converg=True)
                    else:
                        database.update(id=row.id, atoms=atoms, scf=True, converg=True)
                    logging.info(f"id={row.id}的结构scf成功")
                    n_clacs += 1
                except Exception:
                    logging.info(f"id={row.id}的结构scf失败，请查看计算日志")
                    time.sleep(2)
                if isinstance(calc, Vasp):
                    vaspMove(row, mode="scf")
                # 更新
            break  # 这一行一定要有，否则it.combinations会遍历所有可能的情况
    logging.info(f"一共计算了{n_clacs}个结构")
    logging.info("--------------- 结束批量 scf ---------------\n")


# num_list:[2,2...] 吸附0个，1个...的结构都取两个
def sr(calc, db_path, num_list, adsorb_info={"H": 1}, unconvergence_selected=False, only_unconvergence=True, bfgs=True):
    """
    adsorb_kind: 如果吸附的是CO2，则 adsorb_kind = {"C":1, "O":2}
    num_list: 指定进行自洽计算地结构个数，列表索引的第一个是不吸附的结构的选取数量，第二个是吸附一个group的结构的选取数量，以此类推。
    unconvergence_selected: 选中包含未收敛的结构
    only_unconvergence: 如果unconvergence_selected=True，此参数有效，只选中未收敛的结构
    """
    database = db.connect(db_path)
    if isinstance(calc, Vasp):
        bfgs = False
        # 将vasp计算的结果扔进calc_record保存，所以首先要创建这个文件夹
        if not (os.path.exists("calc_record") and os.path.isdir("calc_record")):
            os.mkdir("calc_record")
        calc.directory = ("workdir")

    logging.info("--------------- 开始批量 relax ---------------")
    logging.info(f">>> 是否包含对未收敛的结构进行计算: {unconvergence_selected} <<<")
    logging.info(
        f"0,1,2,3...个吸附物结构的relax数量: {str(num_list).strip('[]')}, 一共需要计算{sum(num_list)}个结构")
    n_clacs = 0

    for ad_num, str_num in enumerate(num_list):

        if str_num == 0:
            continue

        logging.info(f">>> 开始对吸附了{ad_num}个的结构进行relax <<<")

        adsorb_kinds = adsorb_info.keys()  # ["H", ...]
        adsorb_info_str = ",".join(
            map(lambda x: f"{x}={ad_num*adsorb_info[x]}", adsorb_kinds))  # "C=1,O=2"
        # ori_stru_id=0 是选取init中的结构
        row_list = list(database.select(
            f"relaxed=False, ori_stru_id=0, {adsorb_info_str}"))
        if unconvergence_selected:
            row_list_extent = list(database.select(
            f"relaxed_converg=False, {adsorb_info_str}"))
            # 去除冻结和分解的情况
            update_label_decomposed(db_path,id_list=[],value_list=[]) # 给没有decomposed标签的row打上标签，值为False
            rm_list = [row for row in row_list_extent if _isfreeze(db_path, row) or row.decomposed == True]

            for dropped in rm_list:
                row_list_extent.remove(dropped)

            # 去除非末端结构的情况
            rm_list = []
            for row in row_list_extent:
                try:
                    select_relax_target_id = row.relax_target_id
                except Exception:
                    select_relax_target_id = -1
                if select_relax_target_id > 0:
                    rm_list.append(row)

            for dropped in rm_list:
                row_list_extent.remove(dropped)

            if only_unconvergence:
                row_list = row_list_extent
            else:
                row_list.extend(row_list_extent)

        if len(row_list) < str_num:
            str_num_old = str_num
            str_num = len(row_list)
            logging.info(
                f">>> 需要计算的结构有{str_num_old}个，但是数据库中最多只能选择{str_num}个 <<<")

        for row_groups in it.combinations(row_list, str_num):
            for row in row_groups:
                atoms = row.toatoms()
                atoms.calc = calc
                logging.info(f"开始进行第{n_clacs+1}个结构的计算, 该结构的id={row.id}")
                if bfgs:
                    bfgs = BFGS(atoms)
                    bfgs.run(fmax=1e-1)
                try:
                    atoms.get_potential_energy()
                    if isinstance(calc, Vasp):
                        # 计算完成后取得离子步
                        nsw_actual = int(
                            os.popen(r"grep -E 'F=' ./workdir/OSZICAR | wc -l").read())
                        # 判断收敛
                        if nsw_actual == int(calc.parameters["nsw"]):
                            # 另外写入scf后的atoms
                            _id = database.write(
                                atoms, relaxed=True, ori_stru_id=row.id, relaxed_converg=False)
                        else:
                            _id = database.write(
                                atoms, relaxed=True, ori_stru_id=row.id, relaxed_converg=True)
                    else:
                        _id = database.write(
                            atoms, relaxed=True, ori_stru_id=row.id, relaxed_converg=True)
                    # 更新原始的atoms
                    database.update(
                        id=row.id, relaxed=True, relax_target_id=_id)
                    logging.info(f"id={row.id}的结构relax成功")
                    n_clacs += 1
                except Exception:
                    logging.info(f"id={row.id}的结构relax失败，请查看计算日志")
                    time.sleep(2)
                if isinstance(calc, Vasp):
                    vaspMove(row, mode="relax")

            break  # 这一行一定要有，否则it.combinations会遍历所有可能的情况
    
    update_label_decomposed(db_path,id_list=[],value_list=[]) # 给没有decomposed标签的row打上标签，值为False
    logging.info(f"一共计算了{n_clacs}个结构")
    logging.info("--------------- 结束批量 relax ---------------\n")


def scf_id(calc, db_path, id):
    """
    对指定id的结构进行scf计算
    """
    logging.info("--------------- scf for single stru ---------------")
    if isinstance(calc, Vasp):
        # 将vasp计算的结果扔进calc_record保存，所以首先要创建这个文件夹
        if not (os.path.exists("calc_record") and os.path.isdir("calc_record")):
            os.mkdir("calc_record")
        calc.directory = ("workdir")

    mydb = db.connect(db_path)
    select = mydb.get(id=id)

    # 检查init结构是否冻结
    pre_ori_id = _get_pre_ori_id(db_path, select)
    if _isfreeze(db_path, select):
        logging.info(f"id={select.id}的结构的init结构(id={pre_ori_id})已经被冻结,请先解冻\n")
        return

    try:
        bool_scf = select.scf
        if not bool_scf:
            logging.info(f"id={select.id}开始进行scf")
            atoms = select.toatoms()
            atoms.calc = calc
            try:
                atoms.get_potential_energy()
                if isinstance(calc, Vasp):
                    # 计算完成后取得电子步
                    nelm_actual = int(
                        os.popen(r"grep -E 'Iteration[ ]+1' ./workdir/OUTCAR  | wc -l").read())
                    # 判断收敛
                    if nelm_actual == int(calc.parameters["nelm"]):
                        mydb.update(id=select.id, atoms=atoms,
                                    scf=True, converg=True)
                    else:
                        mydb.update(id=select.id, atoms=atoms,
                                    scf=True, converg=False)
                else:
                    mydb.update(id=select.id, atoms=atoms, scf=True, converg=True)
                logging.info(f"id={select.id}的结构scf成功\n")
            except Exception:
                logging.info(f"id={select.id}的结构scf失败，请查看计算日志\n")
            if isinstance(calc, Vasp):
                vaspMove(select, mode="scf")
        else:
            logging.info("该结构已经进行过scf\n")
    except AttributeError:  # 没有scf标签的是relax的结构，对relax的结构进行scf
        logging.info(f"id={select.id}开始进行scf")
        atoms = select.toatoms()
        atoms.calc = calc
        try:
            atoms.get_potential_energy()
            if isinstance(calc, Vasp):
                # 计算完成后取得电子步
                nelm_actual = int(
                    os.popen(r"grep -E 'Iteration[ ]+1' ./workdir/OUTCAR  | wc -l").read())
                # 判断收敛
                if nelm_actual == int(calc.parameters["nelm"]):
                    mydb.update(id=select.id, atoms=atoms,
                                scf=True, converg=True)
                else:
                    mydb.update(id=select.id, atoms=atoms,
                                scf=True, converg=False)
            else:
                mydb.update(id=select.id, atoms=atoms, scf=True)
            logging.info(f"id={select.id}的结构scf成功\n")
        except Exception:
            logging.info(f"id={select.id}的结构scf失败，请查看计算日志\n")
        if isinstance(calc, Vasp):
            vaspMove(select, mode="scf")


def sr_id(calc, db_path, id, bfgs=True):
    """对指定id的结构进行relax计算"""
    logging.info("--------------- relax for single stru ---------------")
    if isinstance(calc, Vasp):
        bfgs = False
        # 将vasp计算的结果扔进calc_record保存，所以首先要创建这个文件夹
        if not (os.path.exists("calc_record") and os.path.isdir("calc_record")):
            os.mkdir("calc_record")
        calc.directory = ("workdir")

    mydb = db.connect(db_path)
    select = mydb.get(id=id)

    pre_ori_id = _get_pre_ori_id(db_path, select)
    if _isfreeze(db_path, select):
        logging.info(f"id={select.id}的结构的init结构(id={pre_ori_id})已经被冻结,请先解冻\n")
        return

    try:
        select_relax_target_id = select.relax_target_id
    except Exception:
        select_relax_target_id = -1

    if (not select.relaxed) or (select.ori_stru_id > 0 and select.relaxed_converg == False and select_relax_target_id<=0):
        logging.info(f"id={select.id}开始进行relax")
        atoms = select.toatoms()
        atoms.calc = calc
        if bfgs:
            bfgs = BFGS(atoms)
            bfgs.run(fmax=1e-1)
        try:
            atoms.get_potential_energy()
            if isinstance(calc, Vasp):
                # 计算完成后取得离子步
                nsw_actual = int(
                    os.popen(r"grep -E 'F=' ./workdir/OSZICAR | wc -l").read())
                # 判断收敛
                if nsw_actual == int(calc.parameters["nsw"]):
                    # 另外写入scf后的atoms
                    _id = mydb.write(
                        atoms, relaxed=True, ori_stru_id=select.id, relaxed_converg=False)
                else:
                    _id = mydb.write(
                        atoms, relaxed=True, ori_stru_id=select.id, relaxed_converg=True)
            else:
                _id = mydb.write(
                    atoms, relaxed=True, ori_stru_id=select.id, relaxed_converg=True)
            # 更新原始的atoms
            mydb.update(
                id=select.id, relaxed=True, relax_target_id=_id)
            logging.info(f"id={select.id}的结构relax成功\n")
        except Exception:
            logging.info(f"id={select.id}的结构relax失败，请查看计算日志\n")
        if isinstance(calc, Vasp):
            vaspMove(select, mode="relax")
    else:
        logging.info("该结构已经进行过relax或者已经收敛\n")
    update_label_decomposed(db_path,id_list=[],value_list=[]) # 给没有decomposed标签的row打上标签，值为False


def _get_pre_ori_id(db_path, row):
    """返回原始init结构的id"""
    mydb = db.connect(db_path)
    try:
        pre_row = mydb[row.ori_stru_id]
        while pre_row.ori_stru_id != [0, -1]:
            row = mydb[row.ori_stru_id]
            pre_row = mydb[row.ori_stru_id]
        return row.ori_stru_id
    except Exception:
        return row.id


def _isfreeze(db_path, row):
    """追踪init结构，如果它冻结了，返回True，否则返回False"""
    mydb = db.connect(db_path)
    pre_ori_id = _get_pre_ori_id(db_path, row)
    pre_ori_row = mydb.get(id=pre_ori_id)
    judge = pre_ori_row.ori_stru_id
    if judge == 0:
        return False
    elif judge == -1:
        return True
    else:
        logging.warning(f"出现一个错误，请检查源文件\n")


def freeze_id(db_path, id):
    """
    对指定id的结构的init结构进行冻结，冻结的结构将不进行scf和relax的计算，也不进行展示，标记为ori_stru_id=-1。
    """
    logging.info("--------------- freeze by id ---------------")
    mydb = db.connect(db_path)
    row = mydb.get(id=id)
    pre_ori_id = _get_pre_ori_id(db_path, row)
    isfreeze = _isfreeze(db_path, row)
    if isfreeze:
        logging.warning(f"id={id}的init结构(id={pre_ori_id})已经被冻结\n")
    else:
        mydb.update(id=pre_ori_id, ori_stru_id=-1)
        logging.info(f"id={id}的init结构(id={pre_ori_id})冻结成功\n")


def unfreeze_id(db_path, id):
    """
    对指定id的结构的init结构进行解冻
    """
    logging.info("--------------- unfreeze by id ---------------")
    mydb = db.connect(db_path)
    row = mydb.get(id=id)
    pre_ori_id = _get_pre_ori_id(db_path, row)
    isfreeze = _isfreeze(db_path, row)
    if isfreeze:
        mydb.update(id=pre_ori_id, ori_stru_id=0)
        logging.info(f"id={id}的init结构解冻成功(id={pre_ori_id})\n")
    else:
        logging.warning(f"id={id}的init结构(id={pre_ori_id})没被冻结\n")


def write_id(db_path, id, format="vasp"):
    """对指定id的结构进行输出，传入的id可以是单个(数字)也可以是多个(列表)"""
    logging.info("--------------- write by id ---------------")
    if isinstance(id, int):
        atoms = db.connect(db_path).get(id=id).toatoms()
        write(f"id_{id}.vasp", images=atoms, format=format)
        logging.info(f"将id={id}的结构以{format}的形式写入成功!\n")
    elif isinstance(id, list):
        for _id in id:
            atoms = db.connect(db_path).get(id=_id).toatoms()
            write(f"id_{_id}.vasp", images=atoms, format=format)
        logging.info(f"将id={id}的结构以{format}的形式写入成功!\n")
    else:
        logging.error("输入有误\n")


def _view_strus(db_path, *args, **kwargs):
    """
    传入db.select()中的参数，查看所选结构，并返回bd列表和对应的原子列表
    """
    # atoms_list = read(f"{db_path}@ori_stru_id=0") # 加上@查看所有结构
    mydb = db.connect(db_path)
    selected_list = list(mydb.select(*args, **kwargs))
    atoms = list(map(lambda x: x.toatoms(), selected_list))
    if not atoms:
        logging.info("没有合适的结构\n")
    else:
        view(atoms)
    return selected_list, atoms


def view_strus(db_path, *args, **kwargs):
    """
    传入db.select()中的参数，查看所选结构，并打印相关信息
    """
    logging.info("--------------- view selected strus ---------------")
    selected_list = _view_strus(db_path, *args, **kwargs)[0]
    # 打印信息
    logging.info("No.  id  is_relaxed")
    for index, row in enumerate(selected_list):
        logging.info("%-5d%-5d%-10s" %
                     (index+1, row.id, row.relaxed))
    logging.info("")


def view_ori_stru(db_path):
    """查看init的除了已经冻结的所有结构"""
    logging.info("--------------- view ori strus ---------------")
    view_strus(db_path, ori_stru_id=0)


def view_strus_id(db_path, id_list):
    """传入id列表(可以查看多个结构)查看结构"""
    mydb = db.connect(db_path)
    selected_list = [mydb.get(id) for id in id_list]
    view(map(lambda x: x.toatoms(), selected_list))
    # 打印信息
    logging.info("No.  id  is_scf  is_relaxed")
    for index, row in enumerate(selected_list):
        logging.info("%-5d%-5d%-8s%-8s" %
                     (index+1, row.id, row.scf, row.relaxed))
    logging.info("")


def view_freeze_stru(db_path):
    logging.info("--------------- view freeze stru ---------------")
    view_strus(db_path, ori_stru_id=-1)


def view_decomposed_stru(db_path):
    logging.info("--------------- view decomposed stru ---------------")
    view_strus(db_path, "decomposed=True, ori_stru_id>0")


def view_undecomposed_stru(db_path):
    logging.info("--------------- view undecomposed stru ---------------")
    view_strus(db_path, "decomposed=False, ori_stru_id>0")


def _view_strus_sorted_by_energy(db_path, select, energy_sort_way="total", single_adsorb_element_energy=None):
    """
    查看指定条件的所有结构。
    结构出现的顺序可以按照总能量排序(sort_way="total")，也可以按照吸附能(sort_way="adsorb_e")，按照吸附能排序需要给出\
    single_adsorb_element_energy，即吸附物进行独立计算得到的单个原子的能量，EMT()计算H2得到单个H的能量为0.5352706311
    """
    mydb = db.connect(db_path)

    if energy_sort_way == "total":
        selected_list = list(mydb.select(select))
        if not selected_list:
            logging.info("没有符合条件的结构")
            return

        row_sorted = sorted(selected_list, key=lambda x: x.energy)

        pre_origin_strus_id = [_get_pre_ori_id(db_path, row) for row in row_sorted]

        atoms = list(map(lambda x: x.toatoms(), row_sorted))

        # 打印
        logging.info("No.   id     total_energy    ori_stru_id    pre_ori_stru_id")  # ori_stru_id是上一个结构的id， pre_ori_stru_id是init结构id
        for index, row in enumerate(row_sorted):
            logging.info("%-6d%-6d%8.4f%15s%17s" % (index+1, row.id, row.energy, row.ori_stru_id, pre_origin_strus_id[index]))

    if energy_sort_way == "adsorb_e":
        try:
            tem_energy = mydb.get(id=1).energy  # 获得模板结构的能量
        except AttributeError:
            logging.info("需要先对模板结构（未吸附物质）进行单点计算")
            return

        n_atoms_tem = len(list(mydb.select(id=1))[0].toatoms().numbers) # 模板结构的原子数
        selected_list = list(mydb.select(f"{select}, natoms!={n_atoms_tem}"))
        if not selected_list:
            logging.info("没有找到相关结构")
            return

        for select in selected_list:
            n_adsorbs = len(select.toatoms().numbers) - \
                len(list(mydb.select(id=1))[0].toatoms().numbers)  # 吸附物的数量
            select.ave_energy = (select.energy - tem_energy -
                                 n_adsorbs * single_adsorb_element_energy) / n_adsorbs
        row_sorted = sorted(selected_list, key=lambda x: x.ave_energy)

        pre_origin_strus_id = [_get_pre_ori_id(db_path, row) for row in row_sorted]

        logging.info("No.   id    adsorb_energy    ori_stru_id    pre_ori_stru_id")
        # 打印排序后的信息 - 方便从数据库中找结构
        for index, row in enumerate(row_sorted):
            logging.info("%-6d%-6d%8.4f%15s%17s" % (index+1, row.id, row.ave_energy, row.ori_stru_id, pre_origin_strus_id[index]))
        atoms = list(map(lambda x: x.toatoms(), row_sorted))
    logging.info("")
    view(atoms)


def view_relax_stru_total(db_path):
    """
    查看结构优化的不包含已经冻结的所有结构
    """
    logging.info("--------------- view relax stru ---------------")
    mydb = db.connect(db_path)
    selected = list(mydb.select("relaxed=True, ori_stru_id>0"))

    # pre_origin_strus_id
    pre_origin_strus_id = [_get_pre_ori_id(db_path, row) for row in selected]

    # 排除已经冻结的结构
    rm_list = [row for row in selected if _isfreeze(db_path, row)]
    for dropped in rm_list:
        selected.remove(dropped)

    atoms = list(map(lambda x: x.toatoms(), selected))

    # 打印排序后的信息 - 方便从数据库中找结构
    logging.info("No.  id     energy    ori_stru_id    pre_ori_stru_id    max_force")
    for index, row in enumerate(selected):
        logging.info("%-5d%-5d%8.4f%9s%16s%19.4f" %
                     (index+1, row.id, row.energy, row.ori_stru_id, pre_origin_strus_id[index], row.fmax))
    logging.info("")
    view(atoms)


def view_scf_strus(db_path, select="scf=True, ori_stru_id=0, converg=True", energy_sort_way="total", single_adsorb_element_energy=None):
    logging.info("--------------- view scf stru (only converged) ---------------")
    _view_strus_sorted_by_energy(db_path, select, energy_sort_way, single_adsorb_element_energy)
    

def view_relax_strus(db_path, select="relaxed=True, ori_stru_id>0, relaxed_converg=True", energy_sort_way="total", single_adsorb_element_energy=None):
    """
    不包含为收敛的结构
    """
    logging.info("--------------- view relax stru (only converged)---------------")
    _view_strus_sorted_by_energy(db_path, select, energy_sort_way, single_adsorb_element_energy)


def update_label_decomposed(db_path: str, id_list: Sequence[int], value_list: Sequence[bool]) -> None:
    """
    更新id_list的decomposed标签改为True，同时，给没有decomposed标签的结构添加该标签，值为False
    """
    # 没有分解的结构添加decompose=False标签
    mydb = db.connect(db_path)
    decomposed = list(map(lambda x: x.id ,mydb.select("decomposed=True")))
    for row in mydb.select():
        if row.id not in decomposed and row.ori_stru_id > 0:
            mydb.update(id = row.id, decomposed=False)
    
    for id,value in zip(id_list, value_list):
        if mydb.get(id=id).ori_stru_id in [0,-1]:
            logging.info(f"该结构id={id}没有经过结构优化，无法添加标签")
        else:
            mydb.update(id = id, decomposed=value)
    logging.info("")