# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 13:44:55 2019

@author: lin
"""

import pandas as pd


with open(r'D:\学习\毕设\数据\data\HELA\SRR3589958.T.txt') as f1:
    
    # 读入ENSG数据，储存于列表
    data = pd.read_table(r'D:\学习\毕设\数据\data\HELA\GC0.5 ENSG.txt', encoding = 'UTF-8', engine = 'python')
    seq_engs = data['EnsemblGeneID'].tolist()
    
    # 初始化变量
    i = 0
    expr = []
    engs_expr = {}
    
    # 根据ENSG号寻找TPM，写入列表中
    for i in range(len(seq_engs)):
        for line in f1:
            if seq_engs[i] in line:
                expr.append(float(line.strip().split(sep = '\t')[3]))
        
        # TPM从大到小排序后写入字典中
        expr.sort(reverse = True)
        engs_expr.setdefault(seq_engs[i], expr)
        
        # 初始化变量，文件指针回到0
        expr = []
        f1.seek(0)
    
    print(engs_expr)
    
    # 字典转为数据框
    writing = pd.DataFrame(pd.Series(engs_expr))
    
    # 将数据框写入文件中
    writing.to_csv(r'D:\学习\毕设\数据\data\HELA\GC0.5 TPM SRR3589958.csv', sep = ',', na_rep = 'NULL', header = True)
        
