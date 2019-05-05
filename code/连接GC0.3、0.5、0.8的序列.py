# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:02:24 2019

@author: lin
"""
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:10:09 2019

@author: lin
"""

with open(r'D:\学习\毕设\数据\data\HELA\gc0.3 cds.txt') as f1:
    
    all_seq_list = []
    
    for line in f1:
        if line.startswith('>') or not(line):
            next
        elif not line.startswith('>'):
            all_seq_list.append(line.strip())
    all_seq = ''.join(all_seq_list)
    
    f2 = open(r'D:\学习\毕设\数据\data\HELA\gc0.3 0.5 0.8 cds.txt', 'a+')
    write_seq = '>gc0.3 cds' + '\n' + all_seq + '\n'
    f2.write(write_seq)
    f2.flush()
    
    
with open(r'D:\学习\毕设\数据\data\HELA\gc0.5 cds.txt') as f3:
    
    all_seq_list = []
    
    for line in f3:
        if line.startswith('>') or not(line):
            next
        elif not line.startswith('>'):
            all_seq_list.append(line.strip())
    all_seq = ''.join(all_seq_list)
    
    write_seq = '>gc0.5 cds' + '\n' + all_seq + '\n'
    f2.write(write_seq)
    f2.flush()
    
     
with open(r'D:\学习\毕设\数据\data\HELA\gc0.8 cds.txt') as f5:
    
    all_seq_list = []
    
    for line in f5:
        if line.startswith('>') or not(line):
            next
        elif not line.startswith('>'):
            all_seq_list.append(line.strip())
    all_seq = ''.join(all_seq_list)
    
    write_seq = '>gc0.8 cds' + '\n' + all_seq + '\n'
    f2.write(write_seq)
    f2.close()