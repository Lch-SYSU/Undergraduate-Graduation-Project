# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:30:26 2019

@author: lin
"""

with open(r'D:\学习\毕设\数据\hk_gene seq\ENC小于等于35 CDS.txt') as f1:
    
    all_seq_list = []
    
    for line in f1:
        if line.startswith('>'):
            next
        elif not line.startswith('>'):
            all_seq_list.append(line.strip())
    all_seq = ''.join(all_seq_list)
    
    f2 = open(r'D:\学习\毕设\数据\hk_gene seq\ENC小于等于35 CDS融合为一条.txt', 'a+')
    write_seq = '>all hk gene(>=300bp)' + '\n' + all_seq
    f2.write(write_seq)
    f2.close()