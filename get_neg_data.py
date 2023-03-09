import os

import numpy as np
import pandas as pd

import cooler

#import imputation function
from imputation.imputation import *
def get_neg_submatrix_intra(clr,window,label_list,impute=False,resolution=100000
                            ,logscale=False,pad=1,std=1,rp=0.5,tol=0.01,window_size=500000000
                            ,step_size=10000000,output_dist=500000000,min_cutoff=0,n_iters=20):
    
    #找到有label对应的位置
    neg_data_list = []
    neg_label_list = []
    
    #遍历每一个sv的breakpoint
    for i in range(len(label_list)):
        #chrom1=chrom2
        chrom1 = label_list['chrom1'].iat[i]
        chrom2 = label_list['chrom2'].iat[i]
        #只取chromosome1的bin
        bin_table = clr.bins().fetch(chrom1)[:]
        #需要重新设置index
        bin_table.index = range(len(bin_table))
        
        breakpoint1 = label_list['breakpoint1'].iat[i]
        breakpoint2 = label_list['breakpoint2'].iat[i]
        
        #找到对应的bin
        bin1 = bin_table[(bin_table['chrom']==chrom1) &(bin_table['start'] < breakpoint1 ) & (breakpoint1 < bin_table['end'] )]
        bin2 = bin_table[(bin_table['chrom']==chrom2) &(bin_table['start'] < breakpoint2 ) & (breakpoint2 < bin_table['end'] )]

        #breakpoint的坐标
        x_c = bin1.index[0]#得到的是索引
        y_c = bin2.index[0]
        
        x1 = x_c - int((window-1)/2)
        x2 = x_c + int((window-1)/2)
        
        y1 = y_c - int((window-1)/2)
        y2 = y_c + int((window-1)/2)
        
        
        #如果x1或y1小于0，就取0
        if x1 < 0:
            x1 = 0
            x2 = window-1
        if y1 < 0:
            y1 = 0
            y2 = window-1
        #如果x2或y2大于区域的长度，就取区域的长度
        max_length= clr.matrix(balance=False).fetch(chrom1).shape[0]

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
            matrix = np.array(matrix) #将matrix转成array
            
        else:
            #raw matrix
            matrix = clr.matrix(balance=False)[:]
        #fetch submatrix
        submatrix = matrix[x1:x2+1, y1:y2+1]

        
        if (submatrix.shape[0] == window):
            neg_data_list.append(submatrix)
            neg_label_list.append("cor")

    return neg_data_list,neg_label_list

#随机选择其中的n个位点，取
def get_neg_loop_matrix(clr,window,loop_list,n_sample,impute,resolution_rwr,
                        logscale,pad,std,rp,tol,
                        window_size,step_size,
                        output_dist,min_cutoff,n_iters):
    
    neg_data_list = []
    neg_label_list = []
    #在loop_list中随机取n行
    sample_loop_list = loop_list.sample(n=n_sample)
    for i in range(len(sample_loop_list)):
        chrom1 = sample_loop_list['chr1'].iat[i]
        chrom2 = sample_loop_list['chr2'].iat[i]

        bin_table = clr.bins().fetch(chrom1)[:]
        bin_table.index = range(len(bin_table))

        pos1 = (sample_loop_list['x1'].iat[i] + sample_loop_list['x2'].iat[i])/2
        pos2 = (sample_loop_list['y1'].iat[i] + sample_loop_list['y2'].iat[i])/2

        #找到对应的bin

        bin1 = bin_table[(bin_table['start']<=pos1) & (bin_table['end']>=pos1)]
        bin2 = bin_table[(bin_table['start']<=pos2) & (bin_table['end']>=pos2)]

        #index用来fetch矩阵
        x_c = bin1.index[0]
        y_c = bin2.index[0]

        x1 = x_c - int((window-1)/2)
        x2 = x_c + int((window-1)/2)

        
        y1 = y_c - int((window-1)/2)
        y2 = y_c + int((window-1)/2)

        if x1 < 0:
            x1 = 0
            x2 = window-1
        if y1 < 0:
            y1 = 0
            y2 = window-1
        #如果x2或y2大于区域的长度，就取区域的长度
        max_length= clr.matrix(balance=False).fetch(chrom1).shape[0]

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
            #将matrix转成array
            matrix = np.array(matrix)
            
        else:
            #raw matrix
            matrix = clr.matrix(balance=False)[:]
        #fetch submatrix
        submatrix = matrix[x1:x2+1, y1:y2+1]

        if(submatrix.shape[0]==window and submatrix.shape[1]==window):
            neg_data_list.append(submatrix)
            neg_label_list.append("loop")
    
    return neg_data_list,neg_label_list


#cancer cell line的
def get_neg_cancer(clr,window,label_list,impute=False,resolution=100000
                            ,logscale=False,pad=1,std=1,rp=0.5,tol=0.01,window_size=500000000
                            ,step_size=10000000,output_dist=500000000,min_cutoff=0,n_iters=20):
    
    #找到有label对应的位置
    neg_data_list = []
    neg_label_list = []
    
    #遍历每一个sv的breakpoint
    for i in range(len(label_list)):
        #chrom1=chrom2
        chrom1 = label_list['chrom1'].iat[i]
        chrom2 = label_list['chrom2'].iat[i]
        #只取chromosome1的bin
        bin_table = clr.bins().fetch(chrom1)[:]
        #需要重新设置index
        bin_table.index = range(len(bin_table))
        
        breakpoint1 = label_list['breakpoint1'].iat[i]
        breakpoint2 = label_list['breakpoint2'].iat[i]
        
        #找到对应的bin
        bin1 = bin_table[(bin_table['chrom']==chrom1) &(bin_table['start'] < breakpoint1 ) & (breakpoint1 < bin_table['end'] )]
        bin2 = bin_table[(bin_table['chrom']==chrom2) &(bin_table['start'] < breakpoint2 ) & (breakpoint2 < bin_table['end'] )]

        #breakpoint的坐标
        x_c = bin1.index[0]#得到的是索引
        y_c = bin2.index[0]

        x_c1 = x_c - window
        y_c1 = y_c - window
        x_c2 = x_c + window
        y_c2 = y_c + window
        x_c3 = x_c - window
        y_c3 = y_c + window
        x_c4 = x_c + window
        y_c4 = y_c - window

        x_c_list = [x_c1,x_c2,x_c3,x_c4]
        y_c_list = [y_c1,y_c2,y_c3,y_c4]

        for i in range(len(x_c_list)):
            x_c = x_c_list[i]
            y_c = y_c_list[i]

            #新的x1x2y1y2
            x1 = x_c - int((window-1)/2)
            x2 = x_c + int((window-1)/2)

            y1 = y_c - int((window-1)/2)
            y2 = y_c + int((window-1)/2)


            #如果x1或y1小于0，就取0
            if x1 < 0:
                x1 = 0
                x2 = window-1
            if y1 < 0:
                y1 = 0
                y2 = window-1
            #如果x2或y2大于区域的长度，就取区域的长度
            max_length= clr.matrix(balance=False).fetch(chrom1).shape[0]

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
                #raw matrix
                matrix = clr.matrix(balance=False)[:]
            #fetch submatrix
            submatrix = matrix[x1:x2+1, y1:y2+1]
            #如果形状正确
            if(submatrix.shape[0]==window and submatrix.shape[1]==window):
                neg_data_list.append(submatrix)
                neg_label_list.append("cancer_no")
        
    return neg_data_list,neg_label_list


#sv list
sv_list = pd.read_csv("./ref_data/sv_list.csv")
k562_list = sv_list[sv_list['Cell Line']=="K562"]
k562_list_intra = k562_list[k562_list['chrom1']==k562_list['chrom2']]
k562_list_intra.index = range(len(k562_list_intra))

#基本配置
resolution = 100000
window = 21


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


gm12878_cool_dir = "./raw_data/scihic/GM12878/cooler/"
gm12878_file_list = os.listdir("./raw_data/scihic/GM12878/cooler/")
gm12878_mcool_list = [i for i in gm12878_file_list if i.endswith(".mcool")]

k562_cool_dir = "./raw_data/scihic/K562/cooler/"
k562_file_list = os.listdir("./raw_data/scihic/K562/cooler/")
k562_mcool_list = [i for i in k562_file_list if i.endswith(".mcool")]




#读取所有的mcool文件
neg_data_list = []
neg_label_list = []

#sv对应的区域
for mc in gm12878_mcool_list:
    clr = cooler.Cooler(gm12878_cool_dir+mc+"::/resolutions/"+str(resolution))
    sc_neg_data_list,sc_neg_label_list = get_neg_submatrix_intra(clr,window,k562_list_intra
                                                                ,False,resolution_rwr,
                                                                logscale,pad,std,rp,tol,
                                                                window_size,step_size,
                                                                output_dist,min_cutoff,n_iters) 
    
    neg_data_list.append(sc_neg_data_list)
    neg_label_list.append(sc_neg_label_list)

#fetch data from list    
neg_data_list = [item for sublist in neg_data_list for item in sublist]
neg_label_list = [item for sublist in neg_label_list for item in sublist]
neg_data_array = np.dstack(neg_data_list)
neg_data_array = np.rollaxis(neg_data_array,-1)
neg_label_list = np.array(neg_label_list)

#保存
np.save("input_data/neg_data_cor_w21_noimpute.npy",neg_data_array)
np.save("input_data/neg_label_cor_w21_noimpute.npy",neg_label_list)


#loop的
loop_path = "/share/home/mliu/sc_sv/raw_data/GM12878/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt.gz"
gm12878_loop_list = pd.read_csv(loop_path,sep="\t",compression='gzip')[['chr1','x1','x2','chr2','y1','y2']]
gm12878_loop_list['chr1'] = gm12878_loop_list['chr1'].apply(lambda x: "chr"+str(x))
gm12878_loop_list['chr2'] = gm12878_loop_list['chr2'].apply(lambda x: "chr"+str(x))




#读取所有的mcool文件
neg_data_loop_list = []
neg_label_loop_list = []

#按细胞做处理
n_sample = 30
for mc in gm12878_mcool_list:
    clr = cooler.Cooler(gm12878_cool_dir+mc+"::/resolutions/"+str(resolution))
    
    #取对应位置的
    sc_neg_data_list,sc_neg_label_list = get_neg_loop_matrix(clr,window,gm12878_loop_list,
                                                                n_sample
                                                                ,False
                                                                ,resolution_rwr,
                                                                logscale,pad,std,rp,tol,
                                                                window_size,step_size,
                                                                output_dist,min_cutoff,n_iters) 
    
    neg_data_loop_list.append(sc_neg_data_list)
    neg_label_loop_list.append(sc_neg_label_list)

#fetch data from list    
neg_data_loop_list = [item for sublist in neg_data_loop_list for item in sublist]
neg_label_loop_list = [item for sublist in neg_label_loop_list for item in sublist]
neg_data_loop_array = np.dstack(neg_data_loop_list)
neg_data_loop_array = np.rollaxis(neg_data_loop_array,-1)
neg_label_loop_list = np.array(neg_label_loop_list)

#保存
np.save("input_data/neg_data_loop_w21_noimpute.npy",neg_data_loop_array)
np.save("input_data/neg_label_loop_w21_noimpute.npy",neg_label_loop_list)


neg_data_cancer_list = []
neg_label_cancer_list = []

for mc in k562_mcool_list:
    clr = cooler.Cooler(k562_cool_dir+mc+"::/resolutions/"+str(resolution))
    #取对应位置的
    sc_neg_data_list,sc_neg_label_list = get_neg_cancer(clr,window,k562_list_intra
                                                                ,False,resolution_rwr,
                                                                logscale,pad,std,rp,tol,
                                                                window_size,step_size,
                                                                output_dist,min_cutoff,n_iters)
    
    neg_data_cancer_list.append(sc_neg_data_list)
    neg_label_cancer_list.append(sc_neg_label_list)

#fetch data from list    
neg_data_cancer_list = [item for sublist in neg_data_cancer_list for item in sublist]
neg_label_cancer_list = [item for sublist in neg_label_cancer_list for item in sublist]
neg_data_cancer_array = np.dstack(neg_data_cancer_list)
neg_data_cancer_array = np.rollaxis(neg_data_cancer_array,-1)
neg_label_cancer_list = np.array(neg_label_cancer_list)

#保存
np.save("input_data/neg_data_cancer_w21_noimpute.npy",neg_data_cancer_list)
np.save("input_data/neg_label_cancer_w21_noimpute.npy",neg_label_cancer_list)