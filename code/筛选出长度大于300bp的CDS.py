# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 11:45:48 2018

@author: lin
"""

with open(r'D:\学习\毕设\数据\hk_gene seq\hk cds.txt') as f1:
    seq = ''
    cds_len = []
    
    for line in f1:
        
        if line.startswith('>') and seq != '':
            len_seq = len(seq)
        
        # 只将序列长度大于300的写入新文件
            if len_seq >= 300:
                f2 = open(r'D:\学习\毕设\数据\hk_gene seq\cds大于等于300bp.txt', 'a+')
                write_seq = title + '\n' + seq + '\n'
                f2.write(write_seq)
            
                cds_len.append(len_seq)
            seq = ''
            
        # 生成注释行
        if line.startswith('>l'):
            title = '>' + line[5:line.index('_c')]

        # 连接纯序列                       
        else:
            seq = seq + line.strip()
            
    f2.close()
    print(cds_len)
    
    
            