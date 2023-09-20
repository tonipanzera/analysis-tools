import numpy as np
from scipy.optimize import curve_fit
import functions

class Rebin:
    
    def __init__(self, data, error, old, new):
        
        self.data = data
        self.error = error
        self.old = old
        self.new = new
        self.new_data = np.zeros(len(self.new)-1)
        self.new_error = np.zeros(len(self.new)-1)
        self.list0 = []
        self.list1 = []
        
    
    def clean(self):
        list_temp = []
        for i in range(len(self.old)):
            
            #list_temp = []
            
            if self.old[i]<self.new[0]:
                #print(i)
                list_temp.append(i)
            if self.new[-1]<self.old[-1-i]:
                list_temp.append(-1-i)
        #print(list_temp)
        self.data = np.delete(self.data, list_temp)
        self.error = np.delete(self.error, list_temp)
        self.old = np.delete(self.old, list_temp)
        
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
            self.new_error = np.delete(self.new_error, 0)
        
        if self.list1[-1]==0:
            self.new = np.delete(self.new, -1)
            self.list1 = np.delete(self.list1, -1)
            self.new_data = np.delete(self.new_data, -1)
            self.new_error = np.delete(self.new_error, -1)
        
        #for i in range(len(self.data)):
            
            #if self.data[i] < 0:
                
                #self.data[i]=0
            
    def rebinData(self):
        
        zero_count = 0
        
        for i in range(len(self.list1)):
            
            start = self.list1[i]
            error_list = []
            
            new_bin = self.new[i+1]-self.new[i]
            
            if i==0:
                
                
                for j in range(start-1):
                    old_bin = self.old[i+j+1]-self.old[i+j]
                    self.new_data[i] +=self.data[i+j]*old_bin/new_bin
                    error_list.append(((old_bin/new_bin)*self.error[i+j])**2)
                   
                old_bin_size = self.old[i+start]-self.old[i+start-1]
                frac = (self.new[i+1]-self.old[i+start-1])/old_bin_size
                
                error_list.append((frac*self.error[i+start-1]*(old_bin_size/new_bin))**2)
                error_list = np.array(error_list)
                
                self.new_data[i] += frac*self.data[i+start-1]*(old_bin_size/new_bin)
                self.new_error[i] = np.sqrt(np.sum(error_list))
            
            
            elif i==len(self.list1)-1:
                
                for j in range(start-1):
                    old_bin = self.old[i+j-zero_count]-self.old[i+j-zero_count-1]
                    self.new_data[i] += self.data[i+j-zero_count]*old_bin/new_bin
                    error_list.append((self.error[i+j-zero_count]*(old_bin/new_bin))**2)
                
                old_bin_size = self.old[i-zero_count]-self.old[i-zero_count-1]
                frac = (self.old[i-zero_count]-self.new[i])/old_bin_size
                
                error_list.append((frac*self.error[i-1-zero_count]*(old_bin_size/new_bin))**2)
                
                error_list = np.array(error_list)
                
                self.new_data[i] += frac*self.data[i-zero_count-1]*(old_bin_size/new_bin)
                self.new_error[i] = np.sqrt(np.sum(error_list))
                
            else:
                
                if start==0:
                    
                    old_bin_size = self.old[i+start-zero_count]-self.old[i+start-zero_count-1]
                    new_bin_size = self.new[i+1]-self.new[i]
                    frac = new_bin_size/old_bin_size
                    
                    self.new_data[i] += frac*self.data[i+start-1-zero_count]
                    self.new_error[i] += np.sqrt((frac*self.error[i+start-1-zero_count])**2)
                
                else:
                    
                    for j in range(start-1):
                        old_bin = self.old[i+j-zero_count]-self.old[i+j-zero_count-1]
                        self.new_data[i] += self.data[i+j-zero_count]*old_bin/new_bin
                        error_list.append((self.error[i+j-zero_count]*(old_bin/new_bin))**2)
                    
                    old_bin_size_left = self.old[i-zero_count]-self.old[i-zero_count-1]
                    #new_bin_left = self.new[i]
                    #print(old_bin_size_left)
                    frac_left = (self.old[i-zero_count]-self.new[i])/old_bin_size_left
                    #print(zero_count)
                    new_left = frac_left*self.data[i-zero_count-1]*(old_bin_size_left/new_bin)
                    error_list.append((frac_left*self.error[i-1-zero_count]*(old_bin_size_left/new_bin))**2)
                    
                    old_bin_size_right = self.old[i+start-zero_count]-self.old[i+start-zero_count-1]
                    frac_right = (self.new[i+1]-self.old[i-zero_count+start-1])/old_bin_size_right
                    new_right = frac_right*self.data[i-zero_count+start-1]*(old_bin_size_right/new_bin)
                    error_list.append((frac_right*self.error[i-zero_count+start-1]*(old_bin_size_right/new_bin))**2)
                    
                    error_list = np.array(error_list)
                    print(new_left,new_right)
                    self.new_data[i] += new_left + new_right
                    self.new_error[i] = np.sqrt(np.sum(error_list))
            #print(zero_count)
            zero_count += -start+1
            
        return self.new[0:len(self.new_data)], self.new_data, self.new_error
    

class FitContinuum:
    
    def __init__(self, x, y, error, x_range_left, x_range_right):
        
        self.x = x
        self.y = y
        self.error = error
        self.x_range_left = x_range_left
        self.x_range_right = x_range_right
    
    def fit_continuum(self):
        
        #Grab the index starts and ends of left and right regions
        start_idx_left = functions.closest_value_idx(self.x, self.x_range_left[0])
        end_idx_left = functions.closest_value_idx(self.x, self.x_range_left[1])
        start_idx_right = functions.closest_value_idx(self.x, self.x_range_right[0])
        end_idx_right = functions.closest_value_idx(self.x, self.x_range_right[1])
        
        #print(self.y)
        #Isolate the x, y, and errors in question
        x_left = self.x[start_idx_left:end_idx_left]
        x_right = self.x[start_idx_right:end_idx_right]
        y_left = self.y[start_idx_left:end_idx_left]
        y_right = self.y[start_idx_right:end_idx_right]
        err_left = self.error[start_idx_left:end_idx_left]
        err_right = self.error[start_idx_right:end_idx_right]
        #print(y_left)
        #print(y_right)
        
        #Concatenate these to have one array for each
        x_tot = np.concatenate((x_left, x_right))
        y_tot = np.concatenate((y_left, y_right))
        err_tot = np.concatenate((err_left, err_right))
        
        #Calculate uncertainty in left + right regions: sigma**2 = sum_i sigma_i**2 / sqrt(N)
        sigma_c = np.sqrt(np.sum(err_tot**2)/len(err_tot))
            
        popt, pcov = curve_fit(functions.linear_func, x_tot, y_tot)
            
        a = popt[0]
        b = popt[1]
        
        return a, b, sigma_c


class Integrate:
    
    '''
    Class to integrate flux in velocity space.
    '''
    
    
    def __init__(self, x, y, error, left_end, right_end, a, b, sigma_c, rest_wave):
        
        self.x = x
        self.y = y
        self.error = error
        self.left_end = left_end
        self.right_end = right_end
        self.a = a
        self.b = b
        self.sigma_c = sigma_c
        self.rest_wave = rest_wave
    
    def integrate(self):
        
        c = 3*10**5
        
        #Isolate indices of start and end of line region
        start = functions.closest_value_idx(self.x, self.left_end)
        end = functions.closest_value_idx(self.x, self.right_end)
        
        #Cut each array accordingly
        cut_y = self.y[start:end]
        cut_x = self.x[start:end]
        cut_err = self.error[start:end]
        
        #Create an array for continuum fit
        xdata = np.linspace(self.left_end, self.right_end, len(cut_y))
        
        #Create continuum fit
        ydata = self.a*xdata + self.b
        
        #Subtract continuum
        cut_y = cut_y - ydata
        
        
        
        y_array = np.zeros(len(cut_y))
        sigma_array = np.zeros(len(cut_y))
        
        #print('Before integration, fractional error is: '+str(cut_err/np.abs(cut_y)))
        #print('Uncertainty in continuum is: '+str(self.sigma_c))
        
        
        
        for i in range(len(cut_y)):
            
            deltaX = 5
            
            y_array[i] = cut_y[i]*deltaX*(self.rest_wave/c)
            sigma_array[i] = cut_err[i]*deltaX*(self.rest_wave/c)
            
        #print('After integration, fractional error is: '+str(sigma_array/np.abs(y_array)))
        
        #Calculate uncertainty in y continuum: sigma Fc = sqrt(sigma**2)*deltaLambda
        sigma_Fc = np.sqrt(self.sigma_c**2)  
        
        #Calculate uncertainty in y
        #sigma_m = np.sqrt(np.sum(sigma_array**2))
        
        #Calculate uncertainty in y line: sigma Fm = sqrt(sigma**2)*deltaLambda
        sigma_Fm = np.sqrt(np.sum(sigma_array**2))#*5*(self.rest_wave/c)
        
        
        #Add up continuum-subtracted y values
        integrated = np.sum(y_array)
        
        sigma = np.sqrt(sigma_Fc**2 + sigma_Fm**2)
        
        print('Sigma line = '+str(sigma_Fm))
        print('Sigma continuum = '+str(sigma_Fc))
        print('Sigma total = '+str(sigma))
        
        return integrated, sigma

class WeHateH2:
    
    def __init__(self, x, y, error, start, end):
        
        self.x = x
        self.y = y
        self.error = error
        self.start = start
        self.end = end
        
    def remove_H2(self):
        
        start = functions.closest_value_idx(self.x, self.start)
        end = functions.closest_value_idx(self.x, self.end)
        
        cut_y = self.y[start:end]
        cut_err = self.error[start:end]
        
        x = np.array([self.y[start], self.y[end]])
        
        y = ((self.y[end]-self.y[start])/(self.x[end]-self.x[start]))*x
        
        y_middle = np.linspace(self.y[start], self.y[end], len(cut_y))
        
        x_new = np.linspace(self.x[start], self.x[end], len(y_middle))
        
        y_new = np.concatenate((self.y[:start], y_middle, self.y[end:]))
        
        return y_middle, x_new, y_new
        
        
            
            
        
        
        
        
        
            
            
        
        
        
        
        