``` toc
min_depth: 1
```
# dpgen 概述
## dpgen 任务的基本流程
-   init：调用第一性原理计算软件 (如 vasp，cp2k) 产生初始数据
-   run：调用 DeePMD-kit、动力学模拟软件 Lammps、第一性原理计算软件 (如 vasp，cp2k，gaussian)，实现同步学习策略，循环迭代自动扩充数据集并提升模型质量
-   test：应用 DP 模型、第一性原理方法及经验力场，调用软件计算参考体系的多种性质

---
## dpgen init
### 帮助
```shell
dpgen -h
dpgen sub_command -h
```

### 输入文件
``` shell
# VASP文件
vasp_env.sh    # vasp计算所需要的环境
POSCAR         # 初始结构
INCAR_md       # 分子动力学的INCAR
INCAR_rlx      # 结构优化的INCAR
POTCAR_N

# dpgen文件
machine.json   # 调用作业系统的配置
param.json     # 参数设置
```

### 基本流程
0. 调用 dpgen                                                       —— machine. json
1. stage1. 产生初始平衡结构 (可跳过)           —— param. json  POSCAR  INCAR_rlx  POTCAR_N
cg-N 0GPa 下的结构优化
2. stage2. 产生 MD 起始构象                           —— param. json
![[../../00_Attachment/dpgen-scale.png]]
3. stage3. 执行 AIMD                                         —— param. json  "stage2 产生的 POSCAR" INCAR_md  POTCAR_N
4. stage4. 将 OUTCAR 打包成数据集     

>[! note]  
> 一次 init_bulk 任务只对一个相进行处理  
> 只进行一个温度下的 AIMD —— 其他温度的 AIMD 用 dpgen run 进行扩充

### 提交任务
``` shell
nohup dpgen init_bulk param.json machine.json > log.out &
```

- 查看提交的任务
``` shell
# 同一个shell下
jobs

# 新的shell下
ps -ef | grep dpgen
```

### 配置文件
#### `param.json`
``` json
{
    "stages" : [1, 2, 3, 4],                         // 指定 init_bulk 的执行步骤。
    "elements": ["N"],                               // 给定了POSCAR就不需要指定原子类型
    "super_cell": [2, 3, 2],                         // 对初始的POSCAR进行扩胞 
    "cell_type": "",                                 // 给定了POSCAR就不需要指定胞类型
    "_latt": "",                                     // 给定了POSCAR就不需要指定晶格常数
    "from_poscar": true,                             // 为true则指定初始POSCAR, 否则自动产生
    "from_poscar_path": "POSCAR",                    // POSCAR所在的路径
    "potcars":  ["./POTCAR_N"],                      // POTCAR所在的路径
    "relax_incar": "./INCAR_rlx",
    "md_incar" : "./INCAR_md",                                        
    "skip_relax": true,                              // 是否对POSCAR进行结构优化
    "scale": [1.00, 0.990],     // 指定 stages 2 放缩 stages 1 结构超胞盒子大小的标量倍率(可执行多个放缩)
    "pert_numb": 5,                                 // 每个放缩结构进行扰动的数目
    "pert_box": 0.03,                                // 晶格常数扰动的程度
    "pert_atom": 0.01,                               // 原子扰动的程度
    "md_nstep" : 30,                                 // 指定 stage 3 vasp中AIMD模拟的步数注意：当赋值与 md_incar 指定文件中 NSW (VASP控制MD步数的关键词)值不同时， stages 3 将使用 NSW 值
    "coll_ndata": 5000,                              // 指定 stages 4 收集 stages 3 数据的最大数量
    "type_map" : ["N"],                              // 指定 stages 4 整理训练数据时，按 deepmd 文件格式记录 元素种类 的顺序注意: 顺序与 elements potcars 及 POSCAR 文件等处指定的元素顺序一致
    "_comment": "that's all"
}
```

^0aa169

#### `machine.json`
^9c1a6d

> [! note]  
> dpgen init 阶段不用 train 和 model_devi，只要设置好 fp (调用 VASP 等计算软件的过程) 就行
``` json
{
  "api_version": "1.0",
  "deepmd_version": "2.1.0",
  "train" :[
    {
      "command": "dp",
      "machine": {
        "batch_type": "Slurm",                         // The batch job system type.
        "machine_type": "Slurm",                       // The connection used to remote machine.
        "context_type": "local",
        "local_root" : "./",                           // The dir where the tasks and relating files locate.
        "remote_root" : "..../init/work/train",              // 执行运算过程的文件夹(临时)，当运算结束后会把里面的文件提取出来
        "_remote_profile": {
            "_hostname": "202.127.207.132",
            "_username": "gchen",
            "_port": "22"
        }
      },
      "resources": {
        "number_node": 1,
        "cpu_per_node": 16,
        "gpu_per_node": 0,
        "queue_name": "hfacnormal01",
        "_custom_flags" : [
            "#SBATCH --partition=hfacnormal01"
            ],
        "_exclude_list": [],
        "_module_list": [],
        "source_list": [
            "..../init/dp_env.sh"
            ],
        "_time_limit": "24:00:00",
        "_mem_limit": "24:00:00",
        "group_size": 1                                  // 一个组的任务个数，一个组只会提交一个作业，组内的任务依次执行
      }
    }
  ],
  "model_devi":
    [{
      "command": "lmp -i input.lammps -v restart 0",
      "machine": {
        "batch_type": "Slurm",
        "machine_type": "Slurm",
        "context_type": "local",
        "local_root" : "./",
        "remote_root" : "..../init/work/model_devi"
      },
      "resources": {
        "number_node": 1,
        "cpu_per_node": 4,
        "gpu_per_node": 0,
        "queue_name": "hfacnormal01",
        "group_size": 1,
        "source_list": [
            "..../init/dp_env.sh"
            ]
      }
    }
  ],
  "fp":
    [{
      "command": "srun vasp_std",
      "machine": {
        "batch_type": "Slurm",
        "machine_type": "Slurm",
        "context_type": "local",
        "local_root" : "./",
        "remote_root" : "..../init/work/fp"
      },
      "resources": {
        "number_node": 1,
        "cpu_per_node": 64,
        "gpu_per_node": 0,
        "queue_name": "hfacnormal01",
        "group_size": 10, 
        "source_list": [
            "..../init/vasp_env.sh"
            ]
      }
    }
  ]
}
```

### 目录结构
``` text
.
├── log.out                               # 详细的任务输出
├── dpdispatcher.log                      # 分配任务和提示任务结束的日志
├── dpgen.log                             # 展示当前所在的stage

├── INCAR_md
├── INCAR_rlx
└── POTCAR_N
├── machine.json
├── param.json

├── cgN.02x03x02
  ├── param.json
  ├── 00.place_ele
  │   ├── 备份文件等
  │   ├── INCAR                           # INCAR_rlx
  │   ├── POSCAR                          # 扩胞后的结构
  │   ├── POSCAR.copid                    # 未扩胞的结构(原始)备份
  │   ├── POTCAR
  │   ├── sys-0096                        # stage1 后结构优化的结构拿到这里
  │      ├── CONTCAR            
  │      ├── INCAR -> ../INCAR
  │      ├── OUTCAR
  │      ├── POSCAR
  │      ├── POTCAR -> ../POTCAR
  ├── 01.scale_pert                        
  │   └── sys-0096
  │       └── scale-1.000
  │           │   POSCAR
  │           ├── 000000
  │           │   └── POSCAR
  │           ├── 000000
  │           │   └── POSCAR
  │           ├── 000001
  │           │   └── POSCAR
  │           ├── 000002
  │           │   └── POSCAR
  │           ├── 000003
  │           │   └── POSCAR
  │           ├── 000004
  │               └── POSCAR
  │       └── scale-0.990
  │           │   POSCAR
  │           ├── 000000
  │           │   └── POSCAR
  │           ├── 000000
  │           │   └── POSCAR
  │           ├── 000001
  │           │   └── POSCAR
  │           ├── 000002
  │           │   └── POSCAR
  │           ├── 000003
  │           │   └── POSCAR
  │           ├── 000004
  │           │   └── POSCAR
  ├── 02.md
    ├── 备份文件等
    ├── INCAR
    ├── POTCAR
    ├── sys-0096
      ├── deepmd
      │   ├── box.raw
      │   ├── coord.raw
      │   ├── energy.raw
      │   ├── force.raw
      │   ├── set.000
      │   │   ├── box.npy
      │   │   ├── coord.npy
      │   │   ├── energy.npy
      │   │   ├── force.npy
      │   │   └── virial.npy
      │   ├── type_map.raw
      │   ├── type.raw
      │   └── virial.raw
      └── scale-1.000
          ├── 000000
          │   ├── CONTCAR
          │   ├── INCAR -> ../../../INCAR
          │   ├── OUTCAR
          │   ├── POSCAR
          │   ├── POTCAR -> ../../../POTCAR
          ├── 000001
          │   └── ...
          ├── 000002
          │   └── ...
          ├── 000003
          │   └── ...
          └── 000004
              └── ...
```
- dpgen. log
``` shell
2022-10-20 21:45:19,890 - INFO : # working dir 1_cgN. POSCAR.02x03x02
2022-10-20 21:45:19,923 - INFO : Elements are N
2022-10-20 21:45:19,923 - INFO : Current stage is 1, relax
2022-10-20 21:45:20,162 - INFO : [96]
2022-10-20 21:58:26,222 - INFO : Current stage is 2, perturb and scale
2022-10-20 21:59:14,619 - INFO : Current stage is 3, run a short md
2022-10-21 20:51:27,745 - INFO : Current stage is 4, collect data
```
---

## dpgen run
### 输入文件
``` shell
# VASP文件
vasp_env.sh    # vasp环境
INCAR_scf      # scf计算所需的INCAR
POTCAR_N

# dp 文件
dp_env.sh      # deepmd-kit环境
machine.json   # 调用作业系统的配置
param.json     # 参数设置

# dpgen init 过程产生的文件
..../init/cgN.POSCAR.02x03x02
..../init/cgN.POSCAR.02x03x02/02.md/sys-0096/deepmd   # init过程产生的训练集
..../init/cgN.POSCAR.02x03x02/01.scale_pert/sys-0096/scale-*/0000**/POSCAR  # init中经过scale的POSCAR
```

### 基本流程
- 00. train : 
> 使用 init 准备的初始训练数据与之前迭代积累的训练数据，调用 DeePMD-kit 训练多个 (默认 4 个) 模型。模型间的唯一区别来自于初始化神经网络时使用不同的随机数种子。

- 01. model_devi : 
> 调用 LAMMPS 使用 00. train 的 **1** 个模型进行 MD 模拟。对于任一 MD 中 snapshot，模型间 (默认 4 个) 预测偏差越大，意味着当前模型系综对该 snapshot 构象的精度越低
> 通过引入模型偏差作为误差判据并设定上下限, 挑选出有希望有效改进模型对 PES 整体预测精度的 snapshot 构象，作为准备加入训练数据集的候选构象。
![[../../00_Attachment/未标题-1.png]]

- 02. fp : 
> 调用 VASP 对 01. model_devi 选取的候选构象进行第一性原理定标 (单点计算)，并调用 dpdata 收集整理所得数据加入到训练数据集中。

> 详细原理，见 [2019_DPGEN_ZhangYZ.pdf](https://dptechnology.feishu.cn/wiki/wikcn5bwMcYXWRQsIs7YWq4e3Nf)

#### 步骤分解
|迭代序号|各迭代的阶段序号|进程|
|---|---|---|
|0|0|make_train|
|0|1|run_train|
|0|2|post_train|
|0|3|make_model_dev|
|0|4|run_model_devi|
|0|5|post_model_devi|
|0|6|make_fp|
|0|7|run_fp|
|0|8|post_fp|
|1|0|make_train|
|1|1|run_train|
|……|……|……|

> 00 `make_train` ：生成训练需要的输入文件，如 input. json；  
> 01 `run_train` ：依据机器配置上传输入文件并执行训练任务；  
> 02 `post_train` ：收集整理输出文件，并分析训练任务的结果。
> ...

### 提交任务
```shell
nohup dpgen run param.json machine.json > log.out &
```

- 查看提交的任务
``` shell
# 同一个shell下
jobs

# 新的shell下
ps -ef | grep dpgen
```

### 配置文件
#### `param.json`
``` json
{
	"_comment": "***** basis set *****",
    "type_map": [
        "N"                                   //指定原子类型。体系有多少种原子，就指定多少种,注意和POSCAR和POTCAR中的原子顺序一样
    ],
    "mass_map": [
		14                                   //type_map中对应元素的 相对原子质量
    ],
    "init_data_prefix": "..../init/",        //为00.train 训练指定 初始数据目录的前缀
    "init_data_sys": [                       //为00.train 指定训练数据集
        "cgN.POSCAR.02x03x02/02.md/sys-0096/deepmd"
    ],
    "init_batch_size": [
        2
    ],
    "sys_configs_prefix": "..../init/",       //01.model_devi 的初始结构 目录前缀
    "sys_configs": [                          //01.model_devi MD 的起始构象。可以使用绝对路径或相对路径，支持通配符。同一个中括号里的结构属于同一个sys
         [
            "cgN.POSCAR.01x01x01/01.scale_pert/sys-0096/scale-1.000/00000[0-2]/POSCAR"
            ],
         [
            "cgN.POSCAR.01x01x01/01.scale_pert/sys-0096/scale-1.000/00000[3-4]/POSCAR",
            "cgN.POSCAR.01x01x01/01.scale_pert/sys-0096/scale-0.990/00000*/POSCAR"
            ]
    ],
	
    "_comment": "***** train *****",
    "numb_models": 4,
    "default_training_param": {
        "model": {
            "type_map": [
                "N"
            ],
            "descriptor": {
                "type": "se_e2_a",
                "sel": [
					96
                ],
                "rcut_smth": 1,
                "rcut": 6.0,
                "neuron": [
                    25,
                    50,
                    100
                ],
                "resnet_dt": false,
                "axis_neuron": 12,
                "seed": 1
            },
            "fitting_net": {
                "neuron": [
                    240,
                    240,
                    240 
                ],
                "resnet_dt": true,
                "seed": 1
            }
        },
        "learning_rate": {
            "type": "exp",
            "start_lr": 0.001,
            "decay_steps": 4000
        },
        "loss": {
            "start_pref_e": 0.02,
            "limit_pref_e": 2,
            "start_pref_f": 1000,
            "limit_pref_f": 1,
            "start_pref_v": 0.0,
            "limit_pref_v": 0.0
        },
        "training": {
            "_set_prefix": "",
            "numb_steps": 800000,
            "disp_file": "lcurve.out",
            "disp_freq": 2000,
            "save_freq": 2000,
            "save_ckpt": "model.ckpt",
            "disp_training": true,
            "time_training": true,
            "profiling": false,
            "profiling_file": "timeline.json",
            "_comment": "that's all"
        }
    },
	"_comment": "***** model_devi *****",
    "model_devi_dt": 0.001,                         //指定 本次迭代中 01.model_devi MD 的时间步长，单位: ps
    "model_devi_skip": 0,
    "model_devi_f_trust_lo": 0.10,                  //指定 01.model_devi 筛选构象的模型偏差力判据下界。
    "model_devi_f_trust_hi": 0.35,
    "model_devi_e_trust_lo": 10000000000.0,
    "model_devi_e_trust_hi": 10000000000.0,
    "model_devi_clean_traj": false,                 //在lammps进行完分子动力学之后是否清除结构文件
    "model_devi_jobs": [
       {"sys_idx": [0],"temps": [100],"press": [1.0], "trj_freq": 5, "nsteps":   1000, "ensemble": "nvt", "_idx":     "00" },
       {"sys_idx": [1],"temps": [100],"press": [1.0], "trj_freq": 5, "nsteps":   1000, "ensemble": "nvt", "_idx":     "01" },
       {"sys_idx": [0],"temps": [100],"press": [1.0], "trj_freq": 10, "nsteps":  3000, "ensemble": "nvt", "_idx":     "02" },
       {"sys_idx": [1],"temps": [100],"press": [1.0], "trj_freq": 10, "nsteps":  3000, "ensemble": "nvt", "_idx":     "03" },
       {"sys_idx": [0],"temps": [200, 300, 400, 500],"press": [1.0], "trj_freq": 5, "nsteps":   1000, "ensemble": "nvt", "_idx":     "04" },
       {"sys_idx": [1],"temps": [200, 300, 400, 500],"press": [1.0], "trj_freq": 5, "nsteps":   1000, "ensemble": "nvt", "_idx":     "05" },
       {"sys_idx": [0],"temps": [200, 300, 400, 500],"press": [1.0], "trj_freq": 10, "nsteps":  3000, "ensemble": "nvt", "_idx":     "06" },
       {"sys_idx": [1],"temps": [200, 300, 400, 500],"press": [1.0], "trj_freq": 10, "nsteps":  3000, "ensemble": "nvt", "_idx":     "07" },
       {"sys_idx": [0],"temps": [200, 300, 400, 500],"press": [1.0], "trj_freq": 20, "nsteps":  5000, "ensemble": "nvt", "_idx":     "08" },
       {"sys_idx": [1],"temps": [200, 300, 400, 500],"press": [1.0], "trj_freq": 20, "nsteps":  5000, "ensemble": "nvt", "_idx":     "09" },
       {"sys_idx": [0],"temps": [200, 300, 400, 500],"press": [1.0], "trj_freq": 40, "nsteps":  8000, "ensemble": "nvt", "_idx":     "10" },
       {"sys_idx": [1],"temps": [200, 300, 400, 500],"press": [1.0], "trj_freq": 40, "nsteps":  8000, "ensemble": "nvt", "_idx":     "11" }
    ],
	"_comment": "***** fp *****",
    "fp_style": "vasp",
    "shuffle_poscar": false,
    "fp_task_max": 130,                           //每一轮迭代选取的最大轨迹数
    "fp_task_min": 8,                             //每一轮迭代选取的最小轨迹数
    "fp_accurate_threshold": 0.999,                //每一轮fp的accurate比例，大于此值，fp_task_max=0
    "fp_accurate_soft_threshold": 0.997,           //accurate在此值和fp_accurate_threshold之间时，fp_task_max则线性衰减
    "_ratio_failed": 0.10,
    "fp_pp_path": "./",
    "fp_pp_files": [
        "POTCAR_N"
    ],
    "fp_incar": "./INCAR_scf"
}
```
`params`参数详解:
|关键词|数据结构类型|例子|描述|
|---|---|---|---|
|type_map|List of String|["Li","Si"]|指定原子类型。体系有多少种原子，就指定多少种|
|mass_map|List of Float|[6.9, 28.1]|type_map 中对应元素的相对原子质量|
|**读入用于训练的初始数据（AIMD 能量，力）**| | | |
|init_data_prefix|String|“~/data/1_LiSi_Example”|为迭代 0 中的 00. train 训练指定初始数据目录的前缀|
|init_data_sys|List of String|[“LiSi. POSCAR.01x01x01/02. md/sys-0016-0016/deepmd”]|为 00. train 指定初始训练数据所在目录，可以使用绝对路径或相对路径。（例如，可以对应 stage 4 中 deepmd 中的数据）。|
|init_batch_size （deepkit 新版本已隐藏该参数）|List of Integer|[1]|每个值均为指定 00. train 过程使用 init_data_sys 初始数据进行训练时，其对应体系的 batch_size。batch size 是指单个 batch（对应 DPMD-kit 训练设置）所包含随机选取构象的帧（snapshot, 一个构象）数量。注：a. 指定值 (batch size) 的个数与 init_data_sys 列表中的字符串（元素）个数要一一对应，且元素顺序也需一一对应。 b.设定 init_batch_size (以及 sys_batch_size ) 的一种推荐原则是：使指定值乘以单帧结构的原子数大于 32。如果设定为 auto (默认), 那 batch size 会自动指定为 32 除以单帧结构的原子数。这个参数理比较复杂，现在已隐藏。|
|**读入 init 过程后获得的初始结构**| | | |
|sys_configs_prefix|String|“~/data/1_LiSi_Example”|指定进行 01. model_devi MD 模拟的起始构象所在目录的前缀|
|sys_configs|List of list of String|[ ["LiSi. POSCAR.01x01x01/02. md/sys-0016-0016/scale-1.000/00000 [0-9]/POSCAR"],[ "LiSi. POSCAR.01x01x01/02. md/sys-0016-0016/scale-1.000/00001 [0-9]/POSCAR"] ]|指定 01. model_devi MD 的起始构象。可以使用绝对路径或相对路径，支持通配符。|
|sys_batch_size （deepkit 新版本已隐藏该参数）|List of Integer|[1, 1]|每个值均为指定 01. model_devi 过程使用 sys_configs 起始进行 MD 样本空间探索时，其对应体系的 batch_size。batch size 是指单个 batch（对应 DPMD-kit 训练设置）所包含随机选取构象的帧数量。注：可设置为 auto , batch size 将指定为 32 除以单帧结构的原子数|
|**train 过程参数设置 (DP 参数)** | | | |
|numb_models|Integer|4|指定 00. train 训练的模型数量|
|default_training_param|Dictionaries|{..."sel": [60, 60],"rcut_smth": 2.0,"rcut": 6.0,"neuron": [25, 50, 100],...}|指定迭代中 00. train 使用 DeePMD-kit 进行训练时的训练参数。相关设置可以参考：`https://github. com/deepmodeling/deepmd-kit`. 注：a.其中的 batch_size 值会被 init_batch_size 和 sys_batch_size 值覆盖;b.一般推荐令: stop_batch = 200 decay_steps|
|# exploration 过程的 model_devi 参数设置 (lammps MD 参数)| | | |
|---|---|---|---|
|model_devi_dt|Float|0.001|指定本次迭代中 01. model_devi MD 的时间步长，单位: ps|
|model_devi_skip|Integer|0|指定 01. model_devi 筛选构象时，跳过每个 MD 轨迹起始帧的数量|
|model_devi_f_trust_lo|Float|0.15|指定 01. model_devi 筛选构象的模型偏差力判据下界。注：力偏差上界和下界是非常重要的两个参数，在训练过程中需要随时根据训练结果调整，如何分析，见 3.4.。详情见文献 DOI：https://doi. org/10.1016/j.cpc.2020.107206|
|model_devi_f_trust_hi|Float|0.35|指定 01. model_devi 筛选构象的模型偏差力判据上界|
|model_devi_e_trust_lo|Float|1e10|指定 01. model_devi 筛选构象的模型偏差能量判据下界。注：建议设置一个较高的数字，因为力偏差可以提供更精确的信息。特殊情况下，如能量最小化可能需要这个|
|model_devi_e_trust_hi|Float|1e10|指定 01. model_devi 筛选构象的模型偏差能量判据上界|
|model_devi_clean_traj|Boolean|false|指定 01. model_devi 是否清除输出文件中记载 MD 轨迹的 traj 文件夹 (以节省空间)。|
|model_devi_jobs|List of Dict|`[{"sys_idx": [0],"temps": [200, 300, 400, 500],"press":[1.0, 100.0, 1000.0, 10000.0, 50000.0],"trj_freq": 5,"nsteps": 1000,"ensembles": "npt","_idx": "00"},...]`|设置 01. model_devi 探索样本空间的 MD 参数。探索是指调用 MD 模拟软件（lammps）进行 DPMD 计算；参数是指包含进行 MD 模拟的温度、压强、步长等。注：每一个 dict {}内的参数对应于一次迭代，即有多个迭代则设置多次 dict{}。其中， model_devi_jobs 中每一个 dict {}的序号与迭代的序号严格对应，即第 n 个 dict {}对应第 n 次迭代。|
|model_devi_jobs["sys_idx"]|List of Integer|[0]|指定 01. model_devi 本次迭代选 sys_configs 中的哪些构象作为 MD 的起始构象。注：选取的构象（帧）的序号与 sys_configs 列表中元素序号（分配构象）需要一一对应。例如，上面参照 param. json 中的 POSCAR 分为 2 个元素—[0]和[1], 那该参数调用这些构象的序号也需要设置为[0], [1], 分为 2 个 dict{}参数; 或者设置为[0, 1]写在一个 dict{}中。|
|model_devi_jobs["temps"]|List of Integer|[200, 250, 300, 400, 500]|指定 01. model_devi 本次迭代 MD 的温度 (K)，可在一次迭代中设置多个温度值。|
|model_devi_jobs["press"]|List of Integer|[1.0, 10.0, 100.0, 1000.0, 10000.0]|指定 01. model_devi 本次迭代 MD 的压力 (Bar)，可在一次迭代中设置多个压力值。|
|model_devi_jobs["trj_freq"]|Integer|20|指定 01. model_devi 本次迭代存储 MD 轨迹中 snapshot 的频率（存储的 snapshot 构象可供筛选）。|
|model_devi_jobs["nsteps"]|Integer|1000|指定 01. model_devi 本次迭代 MD 的步数|
|model_devi_jobs["ensembles"]|String|“npt”|指定 01. model_devi 本次迭代 MD 选取的系综。选项包括 “npt” 和 “nvt”。|
|model_devi_jobs["idx "]|String|“00”|指定 01. model_devi 本次迭代的迭代序号，注：与 model_devi_jobs 的每一个 dict {}的序号一致。|
|model_devi_jobs["taut"]|Float|0.1|恒温器的耦合周期 (ps)|
|model_devi_jobs["taup"]|Float|0.5|恒压器的耦合周期 (ps)|
|# labeling 过程的 first principle (fp) 参数设置 ( vasp/gaussian/cpk2 单点计算参数)| | | |
|---|---|---|---|
|fp_style == VASP| | | |
|fp_style|String|“VASP”|指定 02. fp 调用的第一性原理计算软件。选项：目前包括 "vasp", "pwscf", "siesta","gaussian" 和 "cp2k"。|
|fp_skip_bad_box|String|"length_ratio: 5; height_ratio: 5"| |
|fp_accurate_threshold|Float|0.999|指定准确率，如果准确比率大于这个数字，则不执行 02. fp 中 fp 计算，即 fp_task_max = 0|
|fp_accurate_soft_threshold|Float|0.997|如果准确率在这个指定值与 fp_accurate_threshold 之间，那么 fp_task_max 线性衰减为零|
|shuffle_poscar|Boolean|false| |
|fp_task_max|Integer|50|指定每次迭代 02. fp 中计算选取的最大构象数量，即进行第一性原理单点计算的构象个数上限。注：a.这些构象来自 DPMD 探索样本空间后整个轨迹中的候选构象。轨迹构象成为候选构象的条件：其力偏差在设置好的力偏差上下限之间；b. 当候选构象数量大于指定值时，将随机筛选指定值个构象进行计算。|
|fp_task_min|Integer|3|指定每次迭代 02. fp 中计算选取的最小构象数量，即进行第一性原理单点计算的构象个数下限。注：当候选构象数量小于指定值时，将忽略此次筛选的候选构象，也意味着模型训练精度收敛，可结束该训练过程。|
|ratio_failed|Float|0.05|指定 02. fp 中计算失败的比率值|
|fp_pp_path|String|"~/data/1_LiSi_Example"|指定 02. fp vasp 计算所需赝势文件的存储路径|
|fp_pp_files|List of String|[ "POTCAR_Li","POTCAR_Si"]|指定 02. fp vasp 计算所需的赝势文件。原子类型的顺序（元素）需与 type_map 中记录的原子类型（元素）顺序一致。|
|fp_incar|String|"./INCAR_LiSi"|指定 02. fp vasp 计算所需的 INCAR 输入文件，INCAR 中必须指定 KSPACING。|
|fp_style == Gaussian| | | |
|use_clusters|Boolean|false|指定 02. fp 是否对局域团簇而非整个系统进行计算，需要 DeePMD-kit 版本 1. x。|
|cluster_cutoff|Float|3.5|指定 02. fp 所计算团簇的截断半径，需要 use_clusters ：true|
|fp_params|Dict| |指定 02. fp Gaussian 计算的参数。|
|fp_params["keywords"]|String or list|"mn15/6-31g** nosymm scf (maxcyc=512)"|指定 02. fp Gaussian 计算的关键词。|
|fp_params["multiplicity"]|Integer or String|1|指定 02. fp Gaussian 计算的自旋多重度。注：可以设置为 auto ：自旋多重度将被自动指定。还可以设置为 frag ： "fragment=N"片段组合波函数方法将会被使用。|
|fp_params["nproc"]|Integer|4|指定 02. fp Gaussian 计算使用的核心数量。|
|fp_style == siesta| | | |
|use_clusters|Boolean|false|指定 02. fp 是否对局域团簇而非整个系统进行计算，需要 DeePMD-kit 版本 1. x。|
|cluster_cutoff|Float|3.5|指定 02. fp 所计算团簇的截断半径，需要 use_clusters true 。|
|fp_params|Dict| |指定 02. fp siesta 计算的参数。|
|fp_params["ecut"]|Integer|300|指定 02. fp siesta 计算中平面波对格点的截断。|
|fp_params["ediff"]|Float|1e-4|指定 02. fp siesta 计算中密度矩阵的容差。|
|fp_params["kspacing"]|Float|0.4|指定 02. fp siesta 计算中对布里渊区的采样因子。|
|fp_params["mixingweight"]|Float|0.05|指定 02. fp siesta 计算中 scf 对前一步 scf 所得密度矩阵的继承比例 (线性混合)?|
|fp_params["NumberPulay"]|Integer|5|指定 02. fp siesta 计算对 Pulay convergence accelerator 的设置.|
|fp_style == cp2k| | | |
|fp_params|Dict| |指定 02. fp cp2k 计算的参数。请参考`http://manual. cp2k. org`。只有 kind section 部分是必须设置|

### `machine.json`

> [[dpgen指南#^9c1a6d|点击这里]]

### 目录结构
- `tree .`
```text
├── log.out
├── dpdispatcher.log
├── dpgen.log                   # 了解任务执行的详细信息，比如01.model_devi阶段结束后的候选结构是failed, accurate还是candidate
├── record.dpgen                # 监控当前任务的执行阶段

├── INCAR_scf
├── POTCAR_N
└── vasp_env.sh
├── dp_env.sh

├── iter.000000                 # 目录
├── iter.000001                 # 目录
├── iter.000002                 # 目录
...
├── iter.0000..                 # 目录

├── machine.json
├── param.json
```

- record.dpgen
```text
 0 0
 0 1
 0 2
 0 3
 0 4
 0 5
 0 6
 0 7
 0 8
 1 0
 …….
```

- `ls ./iter.000000/`
```text
00.train  01.model_devi  02.fp
```

- `tree ./iter.000000/00.train`
```shell
├── 000
│   ├── checkpoint
│   ├── frozen_model.pb
│   ├── input.json
│   ├── lcurve.out
│   ├── model.ckpt.data-00000-of-00001
│   ├── model.ckpt.index
│   ├── model.ckpt.meta
│   └── train.log
├── 001
│   ├── checkpoint
│   ├── frozen_model.pb
│   ├── input.json
│   ├── lcurve.out
│   ├── model.ckpt.data-00000-of-00001
│   ├── model.ckpt.index
│   ├── model.ckpt.meta
│   └── train.log
├── 002
│   ├── checkpoint
│   ├── frozen_model.pb
│   ├── input.json
│   ├── lcurve.out
│   ├── model.ckpt.data-00000-of-00001
│   ├── model.ckpt.index
│   ├── model.ckpt.meta
│   └── train.log
├── 003
│   ├── checkpoint
│   ├── frozen_model.pb
│   ├── input.json
│   ├── lcurve.out
│   ├── model.ckpt.data-00000-of-00001
│   ├── model.ckpt.index
│   ├── model.ckpt.meta
│   └── train.log
├── data.init -> ..../init
├── data.iters                                  # 前几轮循环结束后的到的训练集，用于迭代循环训练
├── graph.000.pb -> 000/frozen_model.pb
├── graph.001.pb -> 001/frozen_model.pb
├── graph.002.pb -> 002/frozen_model.pb
└── graph.003.pb -> 003/frozen_model.pb
```

- `tree 01.model.devi`
```shell
├── 01.model_devi
│  ├── confs
│  │  ├── 000.0000.lmp                  # 转换成LAMMPS 输入格式起始结构
│  │  ├── 000.0000.poscar -> ..../cgN.02x03x02/02.md/sys-0096/scale-1.000/000000/POSCAR
│  │  ├── 000.0001.lmp
│  │  ├── 000.0001.poscar -> ..../cgN.02x03x02/02.md/sys-0096/scale-1.000/000001/POSCAR
│  │  …….
│  │  …….
│  ├── cur_job.json                      # param.json 中 model_devi_jobs 指定的当前迭代任务的MD参数信息：
│  ├── graph.000.pb -> ..../iter.000000/00.train/graph.000.pb
│  ├── graph.001.pb -> ..../iter.000000/00.train/graph.001.pb
│  ├── graph.002.pb -> ..../iter.000000/00.train/graph.002.pb
│  ├── graph.003.pb -> ..../iter.000000/00.train/graph.003.pb
│  ├── task.000.000000                   # 000.0000.lmp（第一个结构）的分子动力学
│  │  ├── conf.lmp -> ../confs/000.0000.lmp  # conf.lmp 输入文件：当前任务中， LAMMPS 输入格式的起始构象。
│  │  ├── input.lammps                       # 当前任务由 DP-GEN 自动生成的 LAMMPS 输入脚本。
│  │  ├── job.json
│  │  ├── model_devi.log
│  │  ├── model_devi.out                     # 记录 DPMD 采样构象的能量和力的模型偏差。
│  │  └── traj                               # 根据 param.json 中 "model_devi_f_trust_lo"：0.10 及 "model_devi_f_trust_hi"：0.35 ， DP-GEN 会对保存的构象通过力误差0.10和0.35 eV/Angstrom进行判定和筛选（**这两个参数最重要的参数之一**），将力误差在0.15和0.35 eV/Angstrom之间的探索构象选为“候选构象”，发送到02.fp 文件中。
│  │  ├── 0.lammpstrj
│  │  ├── 20.lammpstrj
│  │  ├  …….
│  ├── task.000.000001                    # 000.0001.lmp（第二个结构）的分子动力学
│  ├── task.000.000002
│  └── ………………..
```

- `tree ./iter.000000/02.fp`
``` shell
─ 02.fp
 ├── candidate.shuffled.000.out           # 记录被筛选为候选构象的结构在 01.model_devi 中所对应的任务号和帧数。
 ├── rest_accurate.shuffled.000.out
 ├── rest_failed.shuffled.000.out
 ├── data.000                             # 候选构象的能量、力、维里等性质的第一性原理计算完成后， DP-GEN 会收集整理这些计算结果并转换为 DeePMD-kit 训练所需数据的格式，存储在该文件夹中。这些数据会与初始数据一起被加入扩大的“数据集”中，在下一个迭代iter.0000(x+1)的 00.train 中，开始新的DP势函数训练。
 │  ├── box.raw
 │  ├── coord.raw
 │  ├── energy.raw
 │  ├── force.raw
 │  ├── set.000
 │  │  ├── box.npy
 │  │  ├── coord.npy
 │  │  ├── energy.npy
 │  │  ├── force.npy
 │  │  └── virial.npy
 │  ├── type_map.raw
 │  ├── type.raw
 │  └── virial.raw
 ├── POTCAR.000                             #  由 param.json 中 fp_pp_path 和 fp_pp_files 指定的 VASP 计算所需两种元素赝势文件合并而成（注意相应位置指定的元素顺序），所有第一性原理任务都使用该赝势文件。
 ├── task.000.000000
 │  ├── conf.dump -> ....iter.000000/01.model_devi/task.000.000227/traj/100.lammpstrj
 │  ├── INCAR
 │  ├── job.json -> ....iter.000000/01.model_devi/task.000.000227/job.json
 │  ├── KPOINTS
 │  ├── OUTCAR
 │  ├── POSCAR
 │  ├── POTCAR -> ../POTCAR.000
│  └── vasprun.xml
├── task.000.000001
├── task.000.000002
├── task.000.000003
``` shell

```shell
cat candidate.shuffled.000.out | grep task.000.000002
>>>
iter.000000/01.model_devi/task.000.000002 60
iter.000000/01.model_devi/task.000.000002 100
iter.000000/01.model_devi/task.000.000002 20
iter.000000/01.model_devi/task.000.000002 80
iter.000000/01.model_devi/task.000.000002 40
<<<
```

-   在下一迭代iter.0000(x+1)的 00.train 中，新的数据将被扩大的“数据集”, 查看和对比 DP-GEN 为iter.000000 和 iter.000001 生成的训练输入脚本看到：
```shell
grep -A 3 "systems" ./iter.000000/00.train/000/input.json
 "systems": [
 "../data.init/cgN.POSCAR.02x03x03/02.md/sys-0096/deepmd"
 ],
 
grep -A 3 "systems" ./iter.000001/00.train/000/input.json
 "systems": [
 "../data.init/LiSi.POSCAR.02x03x02/02.md/sys-0096/deepmd",
 "../data.iters/iter.000000/02.fp/data.000"
 ],
```
>[!attention]
> 注：假设最后一个迭代 iter.000005 （对应6迭代）的 候选构象数 小于fp_task_min 指定值时, 那该迭代实际上并没有进行新的 00.train ，而是拷贝了（软连接）iter.000004的 00.train 文件夹。
