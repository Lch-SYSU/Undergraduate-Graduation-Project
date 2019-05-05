# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 17:57:34 2019

@author: lin
"""

with open(r'D:\学习\毕设\数据\data\cds大于等于300bp.txt') as f1:

    all_seq_list = []
    
    for line in f1:
        if line.startswith('>'):
            next
        elif not line.startswith('>'):
            all_seq_list.append(line.strip())
    
    
    for i in range(len(all_seq_list)):
        print(isinstance(len(all_seq_list[i])/3, int))
    