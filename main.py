# encoding=UTF-8

import numpy as np
import pandas as pd
import time


# ======计算所用函数=======
# 追负荷函数
def h_mo(ex, h1c, qx1):
    if ex < 10:
        a = h1_v(v1(float(h1c)) + (qx1 * (60 * 60)) / (10 ** 8))
        if a > 1880:
            a = 1880

        return a

    a = h1_v(v1(float(h1c)) - (2024.4 * (60 * 60)) / (10 ** 8))
    x_bujin = 10000  # 设置搜索精度值
    for i in range(int((h1c + 2) * x_bujin), int((h1c - 2) * x_bujin), -1):
        bi = i / x_bujin
        if abs(n12(h1c, bi, qx1) - ex) < 0.08 * ex:  # 这里可以改变试算出力精度
            a = bi
            if a > 1880:
                a = 1880
            break

    return a


# 定义出力计算函数
def n12(h1x, h2x, qx):
    n = 0
    h_pingjun = (float(h1x) + float(h2x)) / 2
    q_xiayou = ((v1(float(h1x)) - v1(float(h2x))) * (10 ** 8) / (60 * 60)) + float(qx)
    if q_xiayou > 14100:  # 防止下泄流量过大导致拟合函数失真
        q_xiayou = 14100
    h_xiayou = h2(q_xiayou)
    if q_xiayou < 0:
        n = float('-inf')
    elif q_xiayou > 2024.4:
        n = 8.5 * (h_pingjun - h_xiayou) * 2024.4 / 1000
        if n < 0:  # 最小出力限制条件
            n = float('-inf')
        elif n > 3600:  # 最大出力限制条件(装机容量)
            n = 3600
    else:
        n = 8.5 * (h_pingjun - h_xiayou) * q_xiayou / 1000
        if n < 0:  # 最小出力限制条件
            n = float('-inf')
        elif n > 3600:  # 最大出力限制条件(装机容量)
            n = 3600
    return n


# 根据上游库容计算上有水位 计算结果是m  精确替换完毕
def h1_v(vx):
    h = -1.68443652e-06 * (float(vx) ** 4) + 4.77965039e-04 * (float(vx) ** 3) - 5.72681255e-02 * (float(vx) ** 2) + \
        4.60415170e+00 * float(vx) + 1.70524490e+03
    return h


# 根据上游水位计算上游库容 计算结果库容是亿立方米  精确替换完毕
def v1(hx):
    v = -4.07925391e-09 * (float(hx) ** 4) + 3.31546219e-05 * (float(hx) ** 3) - 9.76341457e-02 * (float(hx) ** 2) + \
        1.24804462e+02 * float(hx) - 5.88202330e+04
    return v


# 根据下泄流量计算下游平均水位  精确替换完毕
def h2(qx):
    h = -1.17856933e-15 * (float(qx) ** 4) + 3.91564987e-11 * (float(qx) ** 3) - 4.86880318e-07 * (float(qx) ** 2) + \
        4.20518306e-03 * float(qx) + 1.63312600e+03
    return h


# 二分法确定电网负荷线函数
def binary_search(i_part, h1_start_part, h1_end_part, fh_DianWang1_part, q_Ru1_part, feng_chuli_yc1_part, guang_chuli_yc1_part):
    low = 0
    high = 3600
    mid = 1800
    h1x_YuCeJiHua = np.zeros((25, 1))  # 记录水位变化过程
    nx_zong = np.zeros((24, 1))  # 局部变量，记录电网负荷值

    while abs(h1x_YuCeJiHua[24, 0] - h1_end_part) > 0.07:
        fhx_min = fh_DianWang1_part.min()
        nx_zong = np.array(fh_DianWang1_part) - fhx_min + mid
        nx_zong = np.minimum(nx_zong, 3600)

        # 确定水电出力
        nx_shui = np.maximum(nx_zong - feng_chuli_yc1_part - guang_chuli_yc1_part, 0)

        # 按照当前水电计划出力后，计算水库末水位
        h1x_YuCeJiHua[0, 0] = h1_start_part
        for j in range(24):
            h1x_YuCeJiHua[j + 1, 0] = h_mo(nx_shui[j, 0], h1x_YuCeJiHua[j, 0], q_Ru1_part[j, 0])

        if mid > 3599 and h1x_YuCeJiHua[24, 0] > h1_end_part:
            print("按照通道满发")
            break
        else:
            if h1x_YuCeJiHua[24, 0] < h1_end_part:
                high = mid
                mid = (low + high) / 2
            else:
                low = mid
                mid = (low + high) / 2

        if high - low < 0.005:
            print("渐进波动收敛")
            break

    global h1_YuCeJiHua_res
    h1_YuCeJiHua_res[:, i_part] = h1x_YuCeJiHua[:, 0]  # 记录水电预测计划（水位）
    global n_zong_YuCe_res
    n_zong_YuCe_res[:, i_part] = nx_zong[:, 0]  # 记录系统总出力计划

    return mid


# 水电追负荷出力计算函数
def waterE_persue(i_part, h1_start_part, fh_MinJiHua_part, fh_DianWang1_part, q_Ru1_part, feng_chuli_sc1_part, guang_chuli_sc1_part):
    h1x_ShiJi = np.zeros((25, 1))  # 记录水位变化过程
    nx_zong = np.zeros((24, 1))  # 局部变量，记录电网负荷值

    fhx_min = fh_DianWang1_part.min()
    nx_zong = np.array(fh_DianWang1_part) - fhx_min + fh_MinJiHua_part
    nx_zong = np.minimum(nx_zong, 3600)

    # 确定水电出力
    nx_shui = np.maximum(nx_zong - feng_chuli_sc1_part - guang_chuli_sc1_part, 0)

    # 按照当前水电出力，计算水库末水位
    h1x_ShiJi[0, 0] = h1_start_part
    for j in range(24):
        h1x_ShiJi[j + 1, 0] = h_mo(nx_shui[j, 0], h1x_ShiJi[j, 0], q_Ru1_part[j, 0])

    global h1_ShiJi_res
    h1_ShiJi_res[:, i_part] = h1x_ShiJi[:, 0]
    global n_shui_ShiJi_res
    n_shui_ShiJi_res[:, i_part] = nx_shui[:, 0]

    return h1x_ShiJi[24, 0]


# ======导入数据========
# 导入水电站数据
pd_data_shuidianzhan1 = pd.read_excel("水电站数据.xlsx", sheet_name="2020")  # 导入数据表
pd_array_Q = np.array(pd_data_shuidianzhan1.loc[:, "入库流量"])
pd_array_H1 = np.array(pd_data_shuidianzhan1.loc[:, "上游水位"])
pd_array_N = np.array(pd_data_shuidianzhan1.loc[:, "出力"])
# 导入风光数据
pd_data_fengguang1 = pd.read_excel("风光数据.xlsx", sheet_name="2020")  # 导入数据表
pd_array_feng_sc = np.array(pd_data_fengguang1.loc[:, "风电实测"])
pd_array_guang_sc = np.array(pd_data_fengguang1.loc[:, "光伏实测"])
pd_array_feng_yc = np.array(pd_data_fengguang1.loc[:, "风电预测"])
pd_array_guang_yc = np.array(pd_data_fengguang1.loc[:, "光伏预测"])
# 用于装载计算数据的7个矩阵均为一维数组


# ======传入计算日具体数据=======
h1_start = pd_array_H1[0]
h1_end = pd_array_H1[24]
fh_DianWang1 = np.array(pd_array_N[0:24]).reshape(24, 1)
q_Ru1 = np.array(pd_array_Q[0:24]).reshape(24, 1)
feng_chuli_yc1 = np.array(pd_array_feng_yc[0:24]).reshape(24, 1)
guang_chuli_yc1 = np.array(pd_array_guang_yc[0:24]).reshape(24, 1)

feng_chuli_sc1 = np.array(pd_array_feng_sc[0:24]).reshape(24, 1)
guang_chuli_sc1 = np.array(pd_array_guang_sc[0:24]).reshape(24, 1)


# 计算过程记录
h1_YuCeJiHua_res = np.zeros((25, 365))  # 记录预测的计划（水位）
n_zong_YuCe_res = np.zeros((24, 365))  # 记录根据预测指定的电网负荷计划
h1_ShiJi_res = np.zeros((25, 365))  # 记录实际的水电计划（水位）
n_shui_ShiJi_res = np.zeros((24, 365))  # 记录实际的水电出力

fh_MinJiHua = np.zeros((1, 365))

# ======计算第一日======
fh_MinJiHua[0, 0] = binary_search(0, h1_start, h1_end, fh_DianWang1, q_Ru1, feng_chuli_yc1, guang_chuli_yc1)
h1_365_mo = waterE_persue(0, h1_start, fh_MinJiHua[0, 0], fh_DianWang1, q_Ru1, feng_chuli_sc1, guang_chuli_sc1)
print("第一日计算完成")
print("末水位为:", h1_365_mo)


# ======计算全年======
for i in range(364):
    h1_start = h1_ShiJi_res[24, i]
    h1_end = pd_array_H1[24 * (i+2)]
    fh_DianWang1 = np.array(pd_array_N[(24 * (i+1)):(24 * (i+2))]).reshape(24, 1)
    q_Ru1 = np.array(pd_array_Q[(24 * (i+1)):(24 * (i+2))]).reshape(24, 1)
    feng_chuli_yc1 = np.array(pd_array_feng_yc[(24 * (i+1)):(24 * (i+2))]).reshape(24, 1)
    guang_chuli_yc1 = np.array(pd_array_guang_yc[(24 * (i+1)):(24 * (i+2))]).reshape(24, 1)

    feng_chuli_sc1 = np.array(pd_array_feng_sc[(24 * (i+1)):(24 * (i+2))]).reshape(24, 1)
    guang_chuli_sc1 = np.array(pd_array_guang_sc[(24 * (i+1)):(24 * (i+2))]).reshape(24, 1)
    print("导入第", (i+2), "天数据")

    # 计算第i+1日
    fh_MinJiHua[0, i+1] = binary_search(i+1, h1_start, h1_end, fh_DianWang1, q_Ru1, feng_chuli_yc1, guang_chuli_yc1)
    h1_365_mo = waterE_persue(i+1, h1_start, fh_MinJiHua[0, i+1], fh_DianWang1, q_Ru1, feng_chuli_sc1, guang_chuli_sc1)
    print("第", (i+2), "日计算完成")
    print("末水位为:", h1_365_mo)


# 将计算结果存到excel中
pd_result_h1_YuCe = pd.DataFrame(h1_YuCeJiHua_res)
pd_result_h1_ShiJi = pd.DataFrame(h1_ShiJi_res)
pd_result_n_zongYuCe = pd.DataFrame(n_zong_YuCe_res)
pd_result_n_shuiShiJi = pd.DataFrame(n_shui_ShiJi_res)

writer = pd.ExcelWriter("output.xlsx")
pd_result_h1_YuCe.to_excel(writer, sheet_name="上游水位-计划", index=False)
pd_result_h1_ShiJi.to_excel(writer, sheet_name="上游水位-实际", index=False)
pd_result_n_zongYuCe.to_excel(writer, sheet_name="计划负荷线", index=False)
pd_result_n_shuiShiJi.to_excel(writer, sheet_name="水电出力-实际", index=False)
writer.close()
print("结果存储完成")
