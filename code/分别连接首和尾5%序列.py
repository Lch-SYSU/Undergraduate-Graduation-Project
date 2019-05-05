# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:10:09 2019

@author: lin
"""

with open(r'D:\学习\毕设\数据\hk_gene seq\ENC极端低 序列.txt') as f1:
    
    all_seq_list = []
    
    for line in f1:
        if line.startswith('>') or not(line):
            next
        elif not line.startswith('>'):
            all_seq_list.append(line.strip())
    all_seq = ''.join(all_seq_list)
    
    f2 = open(r'D:\学习\毕设\数据\hk_gene seq\ENC极端低序列 融合为一条.txt', 'a+')
    write_seq = '>ENC extreme low' + '\n' + all_seq
    f2.write(write_seq)
    f2.close()
    
    
with open(r'D:\学习\毕设\数据\hk_gene seq\ENC极端高 序列.txt') as f3:
    
    all_seq_list = []
    
    for line in f3:
        if line.startswith('>') or not(line):
            next
        elif not line.startswith('>'):
            all_seq_list.append(line.strip())
    all_seq = ''.join(all_seq_list)
    
    f4 = open(r'D:\学习\毕设\数据\hk_gene seq\ENC极端高序列 融合为一条.txt', 'a+')
    write_seq = '>ENC extreme high' + '\n' + all_seq
    f4.write(write_seq)
    f4.close()