import pandas as pd
import numpy as np
import os
#
from scipy.interpolate import interp2d
#
#
def get_coor(data_file, patch_name):
    '''
    data_file: this is the path to the file generate by running the
               command writeCellCentres
    patch_name: this is the name of the patch for which to return 
                coordinates, for example 'inlet'
    '''
    with open(data_file,mode='r') as f_handle:
        #
        # read lines into a list, splitting on new lines
        #
        lines = f_handle.read().splitlines()
    #
    data_strt = 0
    data_end = 0
    #
    # find the data ranges in terms of line numbers
    #
    for (i,line) in enumerate(lines):
        if patch_name in line:
            data_strt = i+6
            data_end = int(lines[i+4])+data_strt
            break
            
    ccx = np.zeros( data_end - data_strt )
            
    for i in range(data_strt, data_end):
        ccx[i-data_strt] = float(lines[i])
        
    return ccx

def return_boundary_velocities(coor_file_base, patch_name, data):
    '''
    coor_file_base: this is the path to the directory where the
                    files generate by running the command
                    writeCellCentres are located
    patch_name: this is the name of the patch
    data: is the pandas data frame generated from the .dat file
    '''
    f_name = os.path.join(coor_file_base,'Cx')
    ccx = get_coor(f_name, patch_name)
    #
    f_name = os.path.join(coor_file_base,'Cy')
    ccy = get_coor(f_name, patch_name)
    #
    f_name = os.path.join(coor_file_base,'Cz')
    ccz = get_coor(f_name, patch_name)
    #
    f_vx = interp2d(np.array(data['z'].values),np.array(data['y'].values), data['Vx'].values)
    f_vy = interp2d(np.array(data['z'].values),np.array(data['y'].values), data['Vy'].values)
    f_vz = interp2d(np.array(data['z'].values),np.array(data['y'].values), data['Vz'].values)
    #
    #
    lines = []
    #
    for x,y,z in zip(ccx,ccy,ccz):
        line = '('
        vx = f_vx(z*1000,y*1000)
        vy = f_vy(z*1000,y*1000)
        vz = f_vz(z*1000,y*1000)
        line = line + '%.8f %.8f %.8f'%(vx,vy,vz) + ')\n'
        lines.append(line)
        
    return lines

def return_boundary_k(coor_file_base, patch_name, data):
    '''
    coor_file_base: this is the path to the directory where the
                    files generate by running the command
                    writeCellCentres are located
    patch_name: this is the name of the patch
    data: is the pandas data frame generated from the .dat file
    '''
    f_name = os.path.join(coor_file_base,'Cx')
    ccx = get_coor(f_name, patch_name)
    #
    f_name = os.path.join(coor_file_base,'Cy')
    ccy = get_coor(f_name, patch_name)
    #
    f_name = os.path.join(coor_file_base,'Cz')
    ccz = get_coor(f_name, patch_name)
    #
    f_k = interp2d(np.array(data['z'].values),np.array(data['y'].values), data['k'].values)
    # 
    k_vals = []
    #
    for x,y,z in zip(ccx,ccy,ccz):
        line = '('
        k = f_k(z*1000,y*1000)
        k_vals.append(k)
        
    return k_vals    

#==============================================================
# Specify where the path to the .dat file you want to use
# here as the data_file variable
#===============================================================
data_file = '/home/josh/Documents/academic/nuen689_cfd/project/ecamp_docs/Vel_Inlet_X-50_u06ms.dat'
#
# This is copied from the .dat file header
VARIABLES = ("x", "y", "z", "Vx", "U_x95", "Vy", "U_y95", "Vz", "U_z95", "RMSVx", "RMSVxUpper", "RMSVxLower", "RMSVy", "RMSVyUpper", "RMSVyLower", "RMSVz", "RMSVzUpper", "RMSVzLower", "k", "kLOWER", "kUPPER")
with open(data_file, mode='r') as f_handle:
    #reading the .dat file into a pandas data frame called data
    data = pd.read_csv(f_handle, sep="\t", skiprows=5, names=VARIABLES, header=None, usecols=np.arange(0,21,1))
#
# the .dat file contains data for both the top and bottom face. This is
# split into two different data frames for interpolation purposes
#
slice_top = []
slice_bottom = []
stride = 13
length = 26
strt = 0
while True:
    slice_top.append( data[strt:strt+stride] )
    slice_bottom.append( data[strt+stride:strt+stride*2] )
    strt += stride*2
    if strt >= stride*length:
        break
df_top = pd.concat( slice_top )
df_bottom = pd.concat( slice_bottom )
#
#==============================================================
# 1) Specify the inlet patch name as patch_name. Note there is 
#    a top and bottom data frame, df_top and df_bottom respectivley,
#    so the patch should be the correct one
# 2) Specify the paths to folder that contains the three 
#    coordinate files generated after running the
#    writeCellCentres openFoam command
#==============================================================
#
# inlet1 is the top face
#
patch_name = 'inlet1'
base_path = '/home/josh/Documents/academic/nuen689_cfd/project/base_bouyant/diff_t/0'
lines = return_boundary_velocities(base_path, patch_name, df_top)


with open(patch_name+'_U.out',mode='w') as f_out:
    f_out.write('value nonuniform List<vector>\n')
    f_out.write(str( len(lines) )+'\n')
    f_out.write('(\n')
    
    for line in lines:
        f_out.write(line)
    
    f_out.write(');\n')
#
# inllet2 is the bottom face
#
patch_name = 'inlet2'
lines = return_boundary_velocities(base_path, patch_name, df_bottom)
#
#
with open(patch_name+'_U.out',mode='w') as f_out:
    f_out.write('value nonuniform List<vector>\n')
    f_out.write(str( len(lines) )+'\n')
    f_out.write('(\n')
        
    for line in lines:
        f_out.write(line)
        
    f_out.write(');\n')
#
#
patch_name = 'inlet1'
k_vals = return_boundary_k(base_path, patch_name, df_top)
#
L = 0.0035
cmu = 0.09
Dh = 2*0.025*0.05/(0.05+0.025)
TINY = 1e-10
SMALL = 1e-7
#
with open(patch_name+'_omega.out',mode='w') as omega_out, \
     open(patch_name+'_k.out',mode='w') as k_out, \
     open(patch_name+'_eps.out',mode='w') as eps_out:
    omega_out.write('value nonuniform List<scalar>\n')
    omega_out.write(str( len(k_vals) )+'\n')
    omega_out.write('(\n')
    k_out.write('value nonuniform List<scalar>\n')
    k_out.write(str( len(k_vals) )+'\n')
    k_out.write('(\n')  
    eps_out.write('value nonuniform List<scalar>\n')
    eps_out.write(str( len(k_vals) )+'\n')
    eps_out.write('(\n')       
    for k in k_vals:
        #
        if k < TINY:
            k = SMALL
        #            
        omega = np.sqrt(k)/Dh/0.07/cmu**(1/4)
        line = '%.8f'%omega+'\n'
        omega_out.write(line)
        #
        line = '%.8f'%k+'\n'
        k_out.write(line)
        #
        eps = 0.09**(3./4)*k**(3./2)/Dh/0.07
        line = '%.8f'%eps+'\n'
        eps_out.write(line)
        
    omega_out.write(');\n')
    k_out.write(');\n')
    eps_out.write(');\n')
#
#
patch_name = 'inlet2'
k_vals = return_boundary_k(base_path, patch_name, df_bottom)
#
#
with open(patch_name+'_omega.out',mode='w') as omega_out, \
     open(patch_name+'_k.out',mode='w') as k_out, \
     open(patch_name+'_eps.out',mode='w') as eps_out:
    omega_out.write('value nonuniform List<scalar>\n')
    omega_out.write(str( len(k_vals) )+'\n')
    omega_out.write('(\n')
    k_out.write('value nonuniform List<scalar>\n')
    k_out.write(str( len(k_vals) )+'\n')
    k_out.write('(\n')  
    eps_out.write('value nonuniform List<scalar>\n')
    eps_out.write(str( len(k_vals) )+'\n')
    eps_out.write('(\n')       
    for k in k_vals:
        #
        if k < TINY:
            k = SMALL
        #            
        omega = np.sqrt(k)/Dh/0.07/cmu**(1/4)
        line = '%.8f'%omega+'\n'
        omega_out.write(line)
        #
        line = '%.8f'%k+'\n'
        k_out.write(line)
        #
        eps = 0.09**(3./4)*k**(3./2)/Dh/0.07
        line = '%.8f'%eps+'\n'
        eps_out.write(line)
        
    omega_out.write(');\n')
    k_out.write(');\n')
    eps_out.write(');\n')
