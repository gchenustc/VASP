""" 
给定旋转轴 和 重合阵点比例 Sigma, 获得 A 结构和 B 结构的晶格矢量
"""
import logging
import numpy as np
from sympy import Max, Symbol, solve, atan, cos, sin
from math import sqrt, pow
import itertools
import sys
import logging
import time


savedStdout = sys.stdout  # 保存标准输出流
printGB = open("OUTPUT.txt", "w")
sys.stdout = printGB

# 第一步，创建一个logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)  # Log等级总开关
# 第三步，再创建一个handler，用于输出到控制台
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)   # 输出到console的log等级的开关
# 第四步，定义handler的输出格式
formatter = logging.Formatter(
    "%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s")
ch.setFormatter(formatter)
# 第五步，将logger添加到handler里面
logger.addHandler(ch)


def vectorize(array):
    """ 单位化向量 array """
    return array / np.linalg.norm(array)

def isgcd(array):
    """求列表是否存在公因数"""

    gcd = 1
    # 只有整数可以被公因数整除, 为什么不直接array.astype(int)？是为了防止小于1的浮点数的误差
    x, y, z = list(map(lambda x:abs(x), sorted(array.round().astype('i4'))))
    
    if z == 0:
        for i in range(1, y + 1):
            if x % i == 0 and y % i == 0:
                gcd = i
    else:
        for i in range(1, z + 1):
            if x % i == 0 and y % i == 0 and z % i == 0:
                gcd = i

    if gcd == 1:
        return False
    return True


def isInteger(x, accuracy=1E-4):
    """ 
    判断 x 是否是整数 ，不用int(x) == x 判断， 
    因为计算机会产生浮点误差，比如就算x是整数,x=1.0000004，int(x) != x
    """
    if -1 < x < 1 and abs(x) < accuracy:
        return True

    if abs(1 - round(x)/x) < accuracy:
        return True
    return False

def isSameMatrix(lst, matrix):
    """
    判断目标矩阵 matrix 是否与之前求得的矩阵等价，之前得到的矩阵都储存在 lst 列表中
    """
    matrix = np.abs(matrix[:,np.argsort(np.abs(matrix)[0])])
    for each_matrix in lst:
        if np.all(each_matrix == matrix):
            return True


def isDestiMatrix(R):
    """
    判断矩阵是否合法，这里的R不是旋转矩阵，而是最终得到的目标矩阵
    1. 矢量间互相垂直
    2. 某一轴内的元素无公因数
    3. 所有元素都为整数
    """
    # 判断是否为整数
    for line in R:
        for num in line:
            if not isInteger(num):
                return False

    b1, b2, b3 = R.round()

    # 判断矢量元素是否为0
    if np.all(b1 == 0) or np.all(b2 == 0) or np.all(b3 == 0):
        return False

    # 如果存在公因数返回False
    if isgcd(b1) or isgcd(b2) or isgcd(b3):
        return False
    
    # 判断垂直
    if np.inner(b1, b2) or np.inner(b1, b3) or np.inner(b2, b3):
        return False
    # print(b1,b2,b3)

    return True
# isDestVector 的 debug
# arr =np.array([[-1,1,-1],[1,1.00000033,1],[2,0,-1]])
# print(isDestiVetor(arr))


def getTheta(Sigma, b, angle):
    """
    由 Sigma 和 旋转轴 b 获得旋转角 theta
    """

    xy_container = []
    theta_container = []

    N = np.inner(b, b)
    # print(N)

    for Sigma in [Sigma * i for i in [1,2,4,8,16]]:
        x = Symbol('x')
        y = Symbol('y')
        # print(Sigma)

        # 求 y 的最大值
        y_max_mulsolve = solve([N * y**2 - Sigma], [y])
        # print(y_max_mulsolve)
        for each_y in y_max_mulsolve:
            each_y = eval(str(each_y[0]))
            if each_y > 0:
                y_max = int(each_y)
        # print(y_max)

        # 解重合点阵坐标 (x,y)
        solved_value = solve([x**2 + N * y**2 - Sigma], [x, y])
        # print(solved_value)
        for each_solve, y in itertools.product(solved_value, range(y_max + 1)):
            x = eval(str(each_solve[0]))
            # print(y,x)
            # if x > 0:
            # print("*")
            #print(abs(1 - round(x)/x))
            # 为了防止 x = 0.99997，int(x) = 0 (此时 x = 1)，所以用 round(x)，四舍五入取整
            if round(x) >= 1 and isInteger(x):
                xy_container.append((round(x), y))
                # print((int(x),y))
        # print(xy_container)

        # 解旋转角度 theta
        if xy_container:
            for each_xy in xy_container:
                # print(each_xy)
                theta_rad = 2 * atan(each_xy[1]/each_xy[0] * sqrt(N))
                # print(theta_rad)
                theta_rad = float('%0.8f' % theta_rad)
                theta_deg = np.rad2deg(theta_rad)
                # 对于立方晶系来说，角度大于90度的情况被重复考虑，比如 顺时针旋转 100度 就等于顺时针旋转 80度

                if 0 < theta_deg <= angle:
                    theta_container.append((theta_rad, theta_deg))
    # 去重
    theta_container = list(set(theta_container))
    if not theta_container:
        return False
    return theta_container


def getR(b, each_theta):
    """
    由旋转轴 b 和旋转角度 theta 求旋转矩阵 R
    """
    b_vect = vectorize(b)
    theta_rad, theta_deg = each_theta
    c = cos(theta_rad)
    s = sin(theta_rad)
    a = 1 - cos(theta_rad)
    Rx = np.array([a * b_vect[0] ** 2 + c, a * b_vect[1] * b_vect[0] +
                    s * b_vect[2], a * b_vect[2] * b_vect[0] - s * b_vect[1]])
    Ry = np.array([a * b_vect[0] * b_vect[1] - s * b_vect[2], a *
                    b_vect[1] ** 2 + c, a * b_vect[2] * b_vect[1] + s * b_vect[0]])
    Rz = np.array([a * b_vect[0] * b_vect[2] + s * b_vect[1], a *
                    b_vect[1] * b_vect[2] - s * b_vect[0], a * b_vect[2] ** 2 + c])
    return np.array([Rx, Ry, Rz]).T


def getM(a,b,c,R,MatrixContainer,flag=0):
    """
    以元组的形式返回，元组内存放目标矩阵，即互成一定角度的一对矩阵
    """

    initMatrix = np.array([a, b, c])
    initMatrix = initMatrix.astype(np.int8)
    # if flag == 1:
    #     print(initMatrix)

    if not isDestiMatrix(initMatrix):
        return
    # 现在的 initVector 是一个合法的向量
    # 计算目标矩阵
    # 这里得到的dtype是object, 要转换成为float，为什么不转换为int？因为还没有对矩阵进行整数判断

    endMatrix = np.dot(initMatrix, R)
    endMatrix = endMatrix.astype(np.float32)
    # if flag == 1:
    #     print('----')
    #     print(initMatrix)
    #     print(R)
    #     print(endMatrix)
    #     print()
        
    if not isDestiMatrix(endMatrix):
        return
    # 现在的 destVector 是一个合法的向量
    # 已经是合法向量了，所以矩阵内都是整数，但是原来的数字诸如 0.99999998, 2.00000012, 0.9999999992不好看, 将之转换为int
    # print("------------------mark1---------------")

    # 检测是否有重复矩阵    
    if isSameMatrix(MatrixContainer, initMatrix) or isSameMatrix(MatrixContainer, endMatrix):
        return

    # 把合法矩阵按照第一行(a矢量)进行递增排序，再取绝对值，添加到容器中
    # 以后得到的矩阵再按照同样的方法和MatrixContainer中的矩阵对比，如果相同视为一个矩阵
    MatrixContainer.append(np.abs(initMatrix[:,np.argsort(np.abs(initMatrix)[0])]))
    MatrixContainer.append(np.abs(endMatrix[:,np.argsort(np.abs(endMatrix)[0])]))
    initM = initMatrix[[2,1,0]]
    endM = endMatrix[[2,1,0]]
    MatrixContainer.append(np.abs(initM[:,np.argsort(np.abs(initM)[0])]))
    MatrixContainer.append(np.abs(endM[:,np.argsort(np.abs(endM)[0])]))

    endMatrix = endMatrix.astype(np.int8)

    return (initMatrix, endMatrix)


def getMatrix(b, R, minnum1=-10, maxnum1=10, minnum2=None, maxnum2=None, minnum3=None, maxnum3=None):
    """
    以循环的形式逐个猜测目标矩阵元素的值，调用 getM 函数求目标矩阵
    """
    
    if minnum2 is None:
        minnum2 = minnum1
    if maxnum2 is None:
        maxnum2 = maxnum1
    if minnum3 is None:
        minnum3 = minnum2
    if maxnum3 is None:
        maxnum3 = maxnum2
    max_M = max(list(map(lambda x:abs(x), [minnum1, minnum2, minnum3, maxnum1, maxnum2, maxnum3])))
    
    matrixContainer = []
    # 下面列表是用来去重的，得到的符合条件的矩阵和下面列表的矩阵对比，如果是等价矩阵则不添加
    MatrixContainer = []

    for i, j, k in itertools.product(range(minnum1, maxnum1+1), range(minnum2, maxnum2+1), range(minnum3, maxnum3+1)):
        # logging.info('The first Matrix circle --> %d,%d,%d, if three number all achieve %d, the current circle is over' % (i,j,k,maxnum1))

        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        solved = solve([x*i+y*j+z*k, x*b[0]+y*b[1]+z*b[2]], [x,y])
        # print(solved)

        if len(solved) == 2:
            for z in range(-max_M,max_M+1):
                names = locals()
                x_ = solved[x] # 这里 x = z符号相关的变量
                y_ = solved[y] # 同理这里也是
                # print(x_,y_)
                x_ = eval(str(x_))
                names['y_'] = eval(str(y_))
                # print(z)
                # print(names['x_'],names['y_'],names['z'])

                
                if names['x_'] > max_M or names['x_'] < -max_M or names['y_'] > max_M or names['y_'] < -max_M:
                    continue

                c = np.array([i, j, k])
                # 如果存在公因数, 则 continue
                a = np.array([names['x_'], names['y_'], names['z']])

                ret = getM(a,b,c,R.T,MatrixContainer)
                
                if not ret:
                    continue
                
                matrixContainer.append(ret)
                
        elif len(solved) == 1:
            if x in list(solved.keys()):
                # print(solved)
                for y, z in itertools.product(range(minnum2, maxnum2+1), range(minnum3, maxnum3+1)):
                    names = locals()
                    x_ = solved[x] # 这里 x = z符号相关的变量
                    names['x_'] = eval(str(x_))

                    c = np.array([i, j, k])
                    a = np.array([names['x_'],y,z])


                    ret = getM(a,b,c,R.T,MatrixContainer)
                    
                    if not ret:
                        continue
                    
                    matrixContainer.append(ret)

            else:
                # print(solved)
                for x, z in itertools.product(range(minnum1, maxnum1+1), range(minnum3, maxnum3+1)):
                    names = locals()
                    y_ = solved[y] # 这里 x = z符号相关的变量
                    names['y'] = eval(str(x_))

                    c = np.array([i, j, k])
                    a = np.array([x,names['y'],z])

                    ret = getM(a,b,c,R.T,MatrixContainer)
                    
                    if not ret:
                        continue
                    
                    matrixContainer.append(ret)

        else:
            for i in range(3):
                for value in range(-max_M, max_M+1):

                    c = np.array([i, j, k])
                    a = np.array([0,0,0])
                    a[i] =  value

                    ret = getM(a,b,c,R.T,MatrixContainer)
                    
                    if not ret:
                        continue
                    
                    matrixContainer.append(ret)

    return matrixContainer


def printf(matrixPir):
    """
    打印目标矩阵A信息
    """

    # 逆向输出(从小到大输出)
    matrixPir = list(reversed(matrixPir))
    lenth = len(matrixPir)
    print("The total number of matrix pair is %d " % lenth)
    print('The next output are InitMatrix and Transformed Matrix respectively:\n')

    for count, each in enumerate(matrixPir, start=1):
        # print(each)
        print('------')
        print("The no.%d matrix pair is:" % count)
        print(each[0], each[1],sep='\n\n')
        print('------')
        print()


if __name__ == "__main__":

    sys.stdout = savedStdout
    b_s = []
    while True:
        arr = input("input the index of the axis you want to search: (eg: 1 0 0) ")
        num1,num2,num3 = list(map(lambda x: int(x), arr.split()))
        b_s.append(np.array([num1,num2,num3]))
        if input("continue? (y/n) ") == 'n':
            break
    # 晶格矢量初始的最大取值
    M = int(input("input the max value of Matrix: "))
    # Sigma 的范围
    Max_Sigma = int(input("input the max value of sigma: "))
    sigma_s = range(3, Max_Sigma+1, 2)
    # 是否考虑大于90°的晶界
    if input("if consider the anger bigger than 90°, input 'y' else input 'n' ") == 'y':
        angle = 181
    else:
        angle = 91
    print("The setting is complete, see OUTPUT.txt for the result")
    #print("The Matrix printed may not correspond to the angle(theta) above, pleate note the structure, it may be the 180-theta, if in this condition, reverse the a and c axis. \
    #In other words, if reverse the a and c axis, theta = 180 - theta")
    sys.stdout = printGB

    # 打印开始时间
    t1 = int(time.time())
    time_local = time.localtime(t1)
    dt = time.strftime("%Y-%m-%d %H:%M:%S", time_local)
    print('---')
    print("start time:", dt)
    print("Max Maxtrix Number: %d, Max Sigma: %d" % (M, list(sigma_s)[-1]))
    print('---')
    print('\n\n')

    for b, sigma in itertools.product(b_s, sigma_s):
        theta_container = getTheta(sigma, b, angle)
        if theta_container:
            for each_theta in theta_container:

                R = getR(b, each_theta)

                # 打印1
                logging.info('rotation vector: %s, Sigma: %d, theta: %s' % (b, sigma, each_theta[1]))
                logging.info('waiting...')
                print('rotate vector:', b, '\nSigma:', sigma, '\ntheta_rad:', each_theta[0], '\ntheta_deg:', each_theta[1])
                print('Rotation Matrix:')
                print(R.astype(np.float32))
                sys.stdout.flush()

                destMatrix = getMatrix(b, R, -M, M)

                if not destMatrix:
                    print('exist matrixPir, but vector value is too small, try larger')
                    print("-"*60)
                    print()
                    sys.stdout.flush()
                    continue

                # 打印2
                print()
                printf(destMatrix)
                print("-"*60)
                print()
                sys.stdout.flush()

    print('---')
    t2 = int(time.time())
    time_local = time.localtime(t2)
    dt = time.strftime("%Y-%m-%d %H:%M:%S", time_local)
    print("end time:", dt)
    print("used time: %d min" % ((t2-t1)/60.0))
    print("---")