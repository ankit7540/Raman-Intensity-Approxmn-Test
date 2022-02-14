# -*- coding: utf-8 -*-
import numpy as np
######################################################

def P_branch(inp, vi, vf, Jmax):
    
    
    nJ=inp.shape[1]
    
    if (Jmax > nJ):
        print('\t Error : JMax exceeds avaibale data. Quitting.')
        exit()
    
    
    out=np.zeros( Jmax )  
    
    count=0
    
    
    for i in range(Jmax, 0, -1):
        #print(i, inp[vf,i-1]- inp[vi,i])
        
        out[count] = inp[vf,i-1]- inp[vi,i]
        count=count+1
    
    #print(' ------------- ')

    return clean_zero(out)

##########################################################

def R_branch(inp, vi, vf, Jmax):
    
    nJ=inp.shape[1]
    
    if (Jmax > nJ):
        print('\t Error : JMax exceeds avaibale data. Quitting.')
        exit()
    
    
    out=np.zeros(Jmax   )  
    
    count=0
    
    
    for i in range(0, Jmax, +1):
        
        t=inp[vf,i+1]- inp[vi,i]
        #print(i, t)
        out[count] = t
        count=count+1
    

    
    return clean_zero(out)


##########################################

def Q_branch(inp, vi, vf, Jmax):
    
    nJ=inp.shape[1]
    
    if (Jmax > nJ):
        print('\t Error : JMax exceeds avaibale data. Quitting.')
        exit()
    
    
    out=np.zeros(Jmax   )  
    
    count=0
    
    
    for i in range(1, Jmax+1, +1):
        
        t=inp[vf,i]- inp[vi,i]
        #print(i, t)
        out[count] = t
        count=count+1
    
    return clean_zero(out)


##########################################

def O_branch(inp, vi, vf, Jmax):
    
    nJ=inp.shape[1]
    
    if (Jmax > nJ):
        print('\t Error : JMax exceeds avaibale data. Quitting.')
        exit()
    
    
    out=np.zeros(Jmax-1   )  
    
    count=0
    
    
    for i in range(Jmax, 1, -1):
        t=inp[vf,i-2]- inp[vi,i]
        #print(i, t)
        
        out[count] = t
        count=count+1
    
    #print(' ------------- ')
    

    
    return clean_zero(out)


##########################################
def S_branch(inp, vi, vf, Jmax):
    
    nJ=inp.shape[1]
    
    if (Jmax > nJ):
        print('\t Error : JMax exceeds avaibale data. Quitting.')
        exit()
    
    
    out=np.zeros(Jmax-1   )  
    
    count=0
    
    
    for i in range(0, Jmax-1, +1):
        
        t=inp[vf,i+2]- inp[vi,i]
        #print(i, t)
        out[count] = t
        count=count+1
    

    
    return clean_zero(out)


##########################################

def clean_zero(inp):
    index= np.where(inp == 0)[0]
    data = np.delete(inp, index)
    return data

##########################################

def clean_zero_2D(inp):   
    return inp[~np.all(inp == 0, axis=1)]    

##########################################
