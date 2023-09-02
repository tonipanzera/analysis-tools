# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

class rebinning:
    
    def __init__(self, data, old, new):
        
        self.data = data
        self.old = old
        self.new = new
        self.new_data = np.zeros(len(self.new)-1)
        self.list0 = []
        self.list1 = []
        
    
    def clean(self):
        
        if self.old[0]<self.new[0]:
            self.old = np.delete(self.old, 0)
            self.data = np.delete(self.data, 0)
        if self.new[-1]<self.old[-1]:
            self.old = np.delete(self.old, -1)
            self.data = np.delete(self.data, -1)
        
        for i in self.new:
            
            diff = i-self.old
            pos = np.sum(diff > 0)
            self.list0.append(pos)
        
        for i in range(len(self.list0)-1):
            self.list1.append(self.list0[i+1]-self.list0[i])
            
        
        if self.list1[0]==0:
            self.new = np.delete(self.new, 0)
            self.list1 = np.delete(self.list1, 0)
            self.new_data = np.delete(self.new_data, 0)
        
        if self.list1[-1]==0:
            self.new = np.delete(self.new, -1)
            self.list1 = np.delete(self.list1, -1)
            self.new_data = np.delete(self.new_data, -1)
            
    def rebin(self):
        
        zero_count = 0
        
        for i in range(len(self.list1)):
            
            start = self.list1[i]
            
            if i==0:
                
                for j in range(start-1):
                    self.new_data[i] +=self.data[i+j]
                    
                old_bin_size = self.old[i+start]-self.old[i+start-1]
                frac = (self.new[i+1]-self.old[i+start-1])/old_bin_size
                self.new_data[i] += frac*self.data[i+start-1]
            
            elif i==len(self.list1)-1:
                
                for j in range(start-1):
                    self.new_data[i] += self.data[i+j-zero_count]
                
                old_bin_size = self.old[i-zero_count]-self.old[i-1-zero_count]
                frac = (self.old[i-zero_count]-self.new[i])/old_bin_size
                self.new_data[i] += frac*self.data[i-1-zero_count]
            
            else:
                
                if start==0:
                    
                    old_bin_size = self.old[i+start-zero_count]-self.old[i+start-zero_count-1]
                    new_bin_size = self.new[i+1]-self.new[i]
                    frac = new_bin_size/old_bin_size
                    self.new_data[i] += frac*self.data[i+start-1-zero_count]
                
                else:
                    
                    for j in range(start-1):
                        self.new_data[i] += self.data[i+j-zero_count]
                    
                    old_bin_size_left = self.old[i-zero_count]-self.old[i-1-zero_count]
                    frac_left = (self.old[i-zero_count]-self.new[i])/old_bin_size_left
                    new_left = frac_left*self.data[i-1-zero_count]
                    
                    old_bin_size_right = self.old[i+start-zero_count]-self.old[i+start-zero_count-1]
                    frac_right = (self.new[i+1]-self.old[i-zero_count+start-1])/old_bin_size_right
                    new_right = frac_right*self.data[i-zero_count+start-1]
                    
                    self.new_data[i] += new_left + new_right
            
            if start==0:
                zero_count += 1
            elif start==1:
                zero_count = zero_count
            elif start==2:
                zero_count -= 1
            elif start==3:
                zero_count -= 2
            elif start==4:
                zero_count -= 3
            elif start==5:
                zero_count -= 4
            
        
        return self.new, self.new_data
            
            
        
        
        
        
        