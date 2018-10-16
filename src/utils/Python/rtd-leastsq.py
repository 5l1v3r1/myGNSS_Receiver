from nav import *
import obs
from rec import *
import numpy as np
import matplotlib.pyplot as plt
from coord import *


def aallinb(a,b):
    for aa in a:
        if aa not in b:
            return False
    return True


PLS = ['G24','G27','G25','G32','G29']
ob1 = obs.Obs('ling1.obs')
ob2 = obs.Obs('ling2.obs')
# 参考站坐标
pos_r1 = np.array(REC15)
# 待求点坐标
pos_r2 = np.array(REC15)

# n = len(PLS)
# 预处理
list_epo = []   # 共有历元
list_dat1 = []  # 处理后的数据(参考站)
list_dat2 = []
for o1 in ob1.obs_data:
    for o2 in ob2.obs_data:
        if o1.epoch == o2.epoch and aallinb(PLS,o1.sat_list) and aallinb(PLS,o2.sat_list):
            list_epo.append(o1.epoch)
            dict_dat1 = {}
            for s in o1.sat_list:
                if s in PLS:
                    dict_dat1[s] = {}
            for i in range(o1.sat_num):
                if o1.sat_list[i] in PLS:
                    dict_dat1[o1.sat_list[i]]['C1'] = o1.data[i].c1
                    dict_dat1[o1.sat_list[i]]['L1'] = o1.data[i].l1
                    dict_dat1[o1.sat_list[i]]['D1'] = o1.data[i].d1
                    dict_dat1[o1.sat_list[i]]['S1'] = o1.data[i].s1
            list_dat1.append(dict_dat1)

            dict_dat2 = {}
            for s in o2.sat_list:
                if s in PLS:
                    dict_dat2[s] = {}
            for i in range(o2.sat_num):
                if o2.sat_list[i] in PLS:
                    dict_dat2[o2.sat_list[i]]['C1'] = o2.data[i].c1
                    dict_dat2[o2.sat_list[i]]['L1'] = o2.data[i].l1
                    dict_dat2[o2.sat_list[i]]['D1'] = o2.data[i].d1
                    dict_dat2[o2.sat_list[i]]['S1'] = o2.data[i].s1
            list_dat2.append(dict_dat2)
            break

# 参考星
ref_sat = max(list_dat1[0].items(), key=lambda x: x[1]['S1'])[0]

# 非参考星列表
no_ref_sat = list(list_dat1[0].keys())
no_ref_sat.remove(ref_sat)

# 画图用
xyz = []
enu = []
ccc = []

[xr, yr, zr] = np.array(pos_r2)+np.array([10, -7, 9])

for i in range(0,len(list_epo)):
    while True:
        Llist = []
        Alist = []
		cc = []
        # 参考星
        pos_s1 = np.array(SAT_POS[ref_sat])
        r_n1 = np.linalg.norm(pos_s1 - pos_r1)
        r_n11 = np.linalg.norm(pos_s1 - [xr, yr, zr])
        dc1 = list_dat2[i][ref_sat]['C1']-list_dat1[i][ref_sat]['C1']
        l0 = -(pos_s1[0] - xr) / r_n11
        m0 = -(pos_s1[1] - yr) / r_n11
        n0 = -(pos_s1[2] - zr) / r_n11

        # 非参考星
        for k in no_ref_sat:
            pos_s2 = np.array(SAT_POS[k])
            r_n2 = np.linalg.norm(pos_s2 - pos_r1)
            r_n22 = np.linalg.norm(pos_s2 - [xr, yr, zr])
            dc2 = list_dat2[i][k]['C1'] - list_dat1[i][k]['C1']
            l1 = -(pos_s2[0] - xr) / r_n22
            m1 = -(pos_s2[1] - yr) / r_n22
            n1 = -(pos_s2[2] - zr) / r_n22

            Llist.append([dc2-dc1+r_n2-r_n1-r_n22+r_n11])
            cc.append(dc2-dc1+r_n2-r_n1-r_n22+r_n11)
            Alist.append([l1-l0,m1-m0,n1-n0])

        L_Mat = np.mat(Llist)
        A_Mat = np.mat(Alist)
        # print(A_Mat.shape)
        X_Mat = (A_Mat.T * A_Mat).I * A_Mat.T * L_Mat
        xr, yr, zr = X_Mat[0, 0] + xr, X_Mat[1, 0] + yr, X_Mat[2, 0] + zr

        if sqrt(X_Mat[0, 0] ** 2 + X_Mat[1, 0] ** 2 + X_Mat[2, 0] ** 2) < 1E-4:
            print(xr, yr, zr)
            break
    xyz.append([xr, yr, zr]-pos_r2)
    enu.append(xyz2enu(pos_r2,[xr, yr, zr])[0:3])
	ccc.append(cc)
xyz = np.array(xyz)
enu = np.array(enu)
ccc = np.array(ccc)
