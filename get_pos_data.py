import os
import numpy as np
import pandas as pd

import cooler

#import imputation function
from imputation.imputation import *

sv_list = pd.read_csv("./ref_data/sv_list.csv")
k562_list = sv_list[sv_list['Cell Line']=="K562"]

#只选择chrom1=chrom2的
k562_list_intra = k562_list[k562_list['chrom1']==k562_list['chrom2']]
k562_list_intra.index = range(len(k562_list_intra))


#列出raw_data/K562/cooler/文件夹下所有的文件
file_list = os.listdir("./raw_data/scihic/K562/cooler/")
mcool_list = [i for i in file_list if i.endswith(".mcool")]#读取所有的mcool文件

#333*36
#rwr config
resolution_rwr=100000
logscale=False
pad=1
std=1
rp=0.5 #restart probability to balance th information between global and local network structures
tol=0.01 #是什么意思
window_size=500000000
step_size=10000000
output_dist=500000000
min_cutoff=0
n_iters=20

def get_pos_submatrix_intra(clr,window,label_list,impute=False,resolution=100000
                            ,logscale=False,pad=1,std=1,rp=0.5,tol=0.01,window_size=500000000
                            ,step_size=10000000,output_dist=500000000,min_cutoff=0,n_iters=20):
    
    #找到有label对应的位置
    #创建空np
    pos_data_list = []
    pos_label_list = []
    # pos_data_list = np.array([])
    # pos_label_list = np.array([])
    
    
    #遍历每一个sv的breakpoint
    for i in range(len(label_list)):
        #chrom1=chrom2
        chrom1 = label_list['chrom1'].iat[i]
        chrom2 = label_list['chrom2'].iat[i]
        
        #取对应chromosome1的bin list
        bin_table = clr.bins().fetch(chrom1)[:]
        bin_table.index = range(len(bin_table))
        #matrix的shape
        max_length= clr.matrix(balance=False).fetch(chrom1).shape[0]

        #得到breakpoint
        breakpoint1 = label_list['breakpoint1'].iat[i]
        breakpoint2 = label_list['breakpoint2'].iat[i]
        #得到标签
        string = label_list['strands'].iat[i]

        #找到breakpoint的pos对应的bin的index，用来取submatrix
        bin1 = bin_table[(bin_table['chrom']==chrom1) &(bin_table['start'] < breakpoint1 ) & (breakpoint1 < bin_table['end'] )]
        bin2 = bin_table[(bin_table['chrom']==chrom2) &(bin_table['start'] < breakpoint2 ) & (breakpoint2 < bin_table['end'] )]
        #中心点的坐标
        x_c = bin1.index[0]#得到的是索引
        y_c = bin2.index[0]

        #起止点的坐标
        if window % 2 == 0:
            x1 = x_c - int((window-1)/2)
            x2 = x_c + int((window-1)/2) + 1
        
            y1 = y_c - int((window-1)/2)
            y2 = y_c + int((window-1)/2) + 1
        else:
            x1 = x_c - int((window-1)/2)
            x2 = x_c + int((window-1)/2)
            
            y1 = y_c - int((window-1)/2)
            y2 = y_c + int((window-1)/2)
  
        #特殊情况下breakpoint不在submatrix的中心
        #如果x1或y1小于0，就取0
        if x1 < 0:
            x1 = 0
            x2 = window-1
        if y1 < 0:
            y1 = 0
            y2 = window-1
        #如果x2或y2大于区域的长度，就取区域的长度
        if x2 > max_length:
            x2 = max_length-1
            x1 = max_length-window
        if y2 > max_length:
            y2 = max_length-1
            y1 = max_length-window
        
        if impute:
            #impute by rwr
            matrix = imputation_rwr(clr,chrom1,resolution,logscale,pad,std,rp,tol,window_size
                    ,step_size,output_dist,min_cutoff,n_iters)
            matrix = np.array(matrix)
   
        else:
            matrix = clr.matrix(balance=False)[:]
 
        submatrix = matrix[x1:x2+1, y1:y2+1]
        print("shape",submatrix.shape)
        if (submatrix.shape[0] == window):
            
            pos_data_list.append(submatrix)
            pos_label_list.append(string)

   
    return pos_data_list,pos_label_list


#这个暂时先不动
resolution = 100000
window = 32
cool_dir = "raw_data/scihic/K562/cooler/"


#创建空np
# pos_data_list = np.array([])
# pos_label_list = np.array([])
pos_data_list = []
pos_label_list = []

#读取所有的mcool文件
i = 1
for mc in mcool_list:
    print(i,mc)
    i = i+1
    clr = cooler.Cooler(cool_dir+mc+"::/resolutions/"+str(resolution))
    sc_pos_data_list,sc_pos_label_list = get_pos_submatrix_intra(clr,window,k562_list_intra
                                                                ,False,resolution_rwr,
                                                                logscale,pad,std,rp,tol,
                                                                window_size,step_size,
                                                                output_dist,min_cutoff,n_iters) 

    pos_data_list.append(sc_pos_data_list)
    pos_label_list.append(sc_pos_label_list)

#fetch data from list    
pos_data_list = [item for sublist in pos_data_list for item in sublist]
pos_label_list = [item for sublist in pos_label_list for item in sublist]
pos_data_array = np.dstack(pos_data_list)
pos_data_array = np.rollaxis(pos_data_array,-1)
print(pos_data_array.shape)
pos_label_list = np.array(pos_label_list)
print(pos_label_list.shape)
#保存

np.save("input_data/pos_data_w32_noimpute.npy",pos_data_array)
np.save("input_data/pos_label_w32_noimpute.npy",pos_label_list)


