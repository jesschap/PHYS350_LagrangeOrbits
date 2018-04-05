# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 01:41:27 2018

@author: aaron
"""

with open('collision_log_2R_EARTH_Jupiter_True.txt','r') as file:
    data = file.readlines()
    
    for line in data:
        velocities = line.split(':')
    print(velocities)
    

