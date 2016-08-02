# -*- coding: utf-8 -*-
"""
Created on Wed May 18 12:49:01 2016

@author: part0075
"""
import numpy as np


def mesh_overlay(layers):
    dz = 1    
    min(layers[0])


def fit_mesh(layers, lay_min = 0.1, lay_max = 100.0):
    mesh = np.zeros(len(layers[0])/2 + 1)
    zone = np.zeros(len(layers[0])/2)
    shifted = np.full((len(layers[0])/2), False, dtype=bool)
    processed = np.full((len(layers[0])/2), False, dtype=bool)    
    index_adjuster = 0
    shift = 0
    
    for index, layer in enumerate(layers[0]):
        # Check thickness of layer is suitable
        if index%2 == 1:
            continue
        #End if
        
        if layers[1][index] == False:
            shifted[index/2] == True
            continue
        #end if
            
        if index/2 != 0:
            if shifted[index/2-1] == True:
                shift = 0                
                for i in range(index/2 - 1, -1, -1):
                    if shifted[i] == True:
                         shift += 1                                           
                    #End if
                 #End for
             #end if
        #end if
                         
        if layers[0][index+1]-layer < lay_min:
            # Check first that it is not top layer and if it is assign minimum thickness and set zone to inactive (i.e. -1)
            if index == len(layers[0])-2:
                mesh[index/2+index_adjuster-shift] = layer
                mesh[index/2+1+index_adjuster-shift] = layer+x_min
                zone[index/2+index_adjuster-shift] = -1
            else:
                shifted[index/2+index_adjuster] = True    
                if layers[1][index+2] == False:
                    shifted[index/2+1] = True
                    continue
                #end if
                mesh[index/2+index_adjuster-shift] = layers[0][index+2]
                mesh[index/2+1+index_adjuster-shift] = layers[0][index+3]
                zone[index/2+index_adjuster-shift] = index/2 + 2
            #End if    
        elif layers[0][index+1]-layer > lay_max:        
            mesh[index/2+index_adjuster-shift] = layer
            new_layers = int((layers[0][index+1]-layer)/x_max) + 1        
            new_thickness = (layers[0][index+1]-layer)/float(new_layers)        
            old_index_adjuster = index_adjuster
            for i in range(new_layers):
                if i < max(range(new_layers)):
                    mesh = np.append(mesh, np.zeros(1), axis=0)        
                    zone = np.append(zone, np.zeros(1), axis=0)
                    shifted = np.append(shifted, np.full(1, False, dtype=bool))                    
                    index_adjuster += 1    
                #End if
                mesh[index/2+1+i+old_index_adjuster-shift] = layer + (i+1) * new_thickness
                zone[index/2+i+old_index_adjuster-shift] = index/2 + 1
            #end for
        else:    
            if shifted[index/2-1] == True:
                shift = 1                
                for i in range(index/2-1 - 1, 0, -1):
                    if shifted[i] == True:
                         shift += 1                                           
                    #End if
                 #End for
                mesh[index/2+index_adjuster-shift] = layer
                mesh[index/2+1+index_adjuster-shift] = layers[0][index+1]
                zone[index/2+index_adjuster-shift] = index/2 + 1
            else:                         
                mesh[index/2+index_adjuster] = layer
                mesh[index/2+1+index_adjuster] = layers[0][index+1]
                zone[index/2+index_adjuster] = index/2 + 1
            #End if
        #End if
    return mesh, zone

if __name__ == "__main__":

    # Cases:
        
    x_min = 0.1
    x_max = 4.0
    
    a_bot = 5.95
    a_top = 6
    b_bot = 6
    b_top = 12
    c_bot = 12
    c_top = 12.01
    
    layers = [[a_bot, a_top, b_bot, b_top, c_bot, c_top], [True, True, False, False, True, True]]
    
    mesh, zone = fit_mesh(layers, lay_min = x_min, lay_max = x_max)