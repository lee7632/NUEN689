import importlib
import lhs_engine as lhs
importlib.reload(lhs)
import os
from shutil import rmtree, copyfile
#
#
def write_velocity(i_run, inlet_num, temp_dir, out_dir):
    
    fname = 'inlet%i_U.out'%inlet_num
    
    with open(join(temp_dir,fname),'r') as u_handle:
        template_U_lines = u_handle.readlines()
    #
    u_out = [None]*len(template_U_lines)
    #
    u_out[:3] = template_U_lines[:3]
    u_out[-1] = template_U_lines[-1]
    #
    for i_line in range(3, len(template_U_lines) -1 ):
        data = template_U_lines[i_line].split()
        
        vx = float(data[0][1:])
        vy = float(data[1])
        vz = float(data[2][:-1])
        
        scale_x = lhs.param_dict['Vx_scale_%i'%inlet_num]['values'][i_run]
        scale_y = lhs.param_dict['Vy_scale_%i'%inlet_num]['values'][i_run]
        scale_z = lhs.param_dict['Vz_scale_%i'%inlet_num]['values'][i_run]
        
        vx *= scale_x
        vy *= scale_y
        vz *= scale_z
        
        u_out[i_line] = '('+str(vx)+' '+str(vy) + ' ' + str(vz) +')\n'
        
        
    with open(join(out_dir,fname),'w') as fout:
        for line in u_out:
            fout.write(line)
#
def write_k(i_run, inlet_num, temp_dir, out_dir):
    
    fname = 'inlet%i_k.out'%inlet_num
    
    with open(join(temp_dir,fname),'r') as u_handle:
        template_k_lines = u_handle.readlines()    
    
    #
    k_out = [None]*len(template_k_lines)
    #
    k_out[:3] = template_k_lines[:3]
    k_out[-1] = template_k_lines[-1]
    #    
    #
    for i_line in range(3, len(template_k_lines) -1 ):
        
        k = float(template_k_lines[i_line])
        
        scale_k = lhs.param_dict['k_scale_%i'%inlet_num]['values'][i_run]
        
        k *= scale_k
        
        k_out[i_line] = str(k)+'\n'
        
        
    with open(join(out_dir,fname),'w') as fout:
        for line in k_out:
            fout.write(line)    

#
def write_omega(i_run, inlet_num, temp_dir, out_dir):
    
    fname = 'inlet%i_omega.out'%inlet_num
    
    with open(join(temp_dir,fname),'r') as u_handle:
        template_om_lines = u_handle.readlines()    
    
    #
    Cmul_nom = 0.0035*0.09
    #
    om_out = [None]*len(template_om_lines)
    #
    om_out[:3] = template_om_lines[:3]
    om_out[-1] = template_om_lines[-1]
    #    
    #
    for i_line in range(3, len(template_om_lines) -1 ):
        
        om = float(template_om_lines[i_line])
        
        CmuL = lhs.param_dict['CmuL']['values'][i_run]
        
        om *= Cmul_nom/CmuL
        
        om_out[i_line] = str(om)+'\n'
        
        
    with open(join(out_dir,fname),'w') as fout:
        for line in om_out:
            fout.write(line)
    
#
def write_block_mesh(i_run, temp_dir, out_dir):
    
    with open(join(temp_dir,'blockMeshDict'),'r') as u_handle:
        block_mesh_lines = u_handle.readlines()
    
    sp_angle = lhs.param_dict['sp_angle']['values'][i_run]
    
    block_mesh_lines[18] = "Angle       "+str(sp_angle)+   \
       ";                              // angle between the mixing channels\n"
       
    inlet_L = lhs.param_dict['inlet_L']['values'][i_run]
       
    block_mesh_lines[22] = "TopL        "+str(inlet_L)+ \
        ";                             // inflow section\n"
    
    with open(join(out_dir,'blockMeshDict'),'w') as fout:
        for line in block_mesh_lines:
            fout.write(line)        
 
def write_nu(i_run, temp_dir, out_dir):
    with open(join(temp_dir,'transportProperties'),'r') as u_handle:
        prop_lines = u_handle.readlines()    
        
    nu = lhs.param_dict['nu']['values'][i_run]
    
    prop_lines[18] = "nuvar  "+str(nu)+";\n"
    
    with open(join(out_dir,'transportProperties'),'w') as fout:
        for line in prop_lines:
            fout.write(line)
  
def write_omega_wall(i_run, temp_dir, out_dir):
    with open(join(temp_dir,'omega'),'r') as u_handle:
        prop_lines = u_handle.readlines()    
        
    beta = lhs.param_dict['beta']['values'][i_run]
    nu = lhs.param_dict['nu']['values'][i_run]
    
    prop_lines[17] = "beta  "+str(beta)+";\n"
    
    prop_lines[18] = "nuvar  "+str(nu)+";\n"
    
    with open(join(out_dir,'omega'),'w') as fout:
        for line in prop_lines:
            fout.write(line)     
#
join = os.path.join
#
template_dir = r'/home/josh/Documents/academic/nuen689_cfd/project/LHS/LHS_template'
output_base_dir = '/home/josh/Documents/academic/nuen689_cfd/project/debug_lhs_output'
#
#
ic_folder = os.path.join(template_dir,"0")
#
runs = lhs.bins
#    
#
#
for i_run in range(runs):
    
    base_name = 'LHS_%i'%i_run
    
    output_dir = os.path.join(output_base_dir, base_name)
    
    if  not  os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        rmtree(output_dir)
        os.makedirs(output_dir)
    #
    # handle the files in the 0 folder
    #
    ic_bc_dir_out = join(output_dir,'0')
    os.makedirs(ic_bc_dir_out)
    #
    ic_bc_dir_tempalte = join(template_dir,'0')
    #
    # do velocity scaling
    #
    write_velocity(i_run, 1, ic_bc_dir_tempalte, ic_bc_dir_out)
    write_velocity(i_run, 2, ic_bc_dir_tempalte, ic_bc_dir_out)
    #
    # Do k scaling
    #
    write_k(i_run, 1, ic_bc_dir_tempalte, ic_bc_dir_out)
    write_k(i_run, 2, ic_bc_dir_tempalte, ic_bc_dir_out)
    #
    # scaling for Cmu*L term in omega calc
    #
    write_omega(i_run, 1, ic_bc_dir_tempalte, ic_bc_dir_out)
    write_omega(i_run, 2, ic_bc_dir_tempalte, ic_bc_dir_out)
    #
    # Account for the beta parameter in the omega file
    #
    write_omega_wall(i_run, ic_bc_dir_tempalte, ic_bc_dir_out)
    #
    # Copy over all other files from the 0 folder
    #
    copy_files = ['k', 'nut', 'p', 'U']
    #
    for f_name in copy_files:
        #
        source = join(ic_bc_dir_tempalte,f_name)
        destination = join(ic_bc_dir_out, f_name)
        #
        copyfile(source, destination)    
    #
    # Replace the Angle and entry length
    #
    sys_tempalte = join(template_dir,'system')
    sys_out_dir = join(output_dir,'system')
    os.makedirs(sys_out_dir)
    #
    #
    write_block_mesh(i_run, sys_tempalte, sys_out_dir)
    #
    # copy over all the other files in the system directory
    #
    copy_files = ['controlDict','decomposeParDict',
                  'fvSchemes', 'fvSolution']
    
    for f_name in copy_files:
        #
        source = join(sys_tempalte,f_name)
        destination = join(sys_out_dir, f_name)
        #
        copyfile(source, destination)
    
    #
    # Doing things in the constant directory
    #
    const_template = join(template_dir,'constant')
    #
    const_out_dir = join(output_dir,'constant')
    os.makedirs(const_out_dir)
    #
    # Copy over the turbulenceProperties file
    #
    copyfile(join(const_template,'turbulenceProperties'), 
             join(const_out_dir,'turbulenceProperties') )
    #
    # Input the uncertainy viscosity value
    #
    write_nu(i_run, const_template, const_out_dir)

slurm_script = join(output_base_dir,'run.slurm')

with open(slurm_script,'w') as f_handle:
    f_handle.write('#!/bin/bash\n')
    f_handle.write('#SBATCH --export=NONE                    #Do not propagate environment\n')
    f_handle.write('#SBATCH --get-user-env=L                 #Replicate login environment\n')
    f_handle.write('#SBATCH --job-name=LHS_NUEN_689\n')
    f_handle.write('#SBATCH --time=4:00:00\n')
    f_handle.write('#SBATCH --nodes=3                        #Request 3 node\n')
    f_handle.write('#SBATCH --ntasks=64\n')
    f_handle.write('#SBATCH --cpus-per-task=1\n')
    f_handle.write('#SBATCH --mem=128G                        #Request 8GB per node\n')
    f_handle.write('#SBATCH --output=out.%j                  #Send stdout/err to "Example2Out.[jobID]"\n')
    #
    #
    for i in range(runs):
        if i > 0 :
            f_handle.write('cd ../LHS_%i\n'%i)
        else:
            f_handle.write('cd LHS_0\n')
        f_handle.write('blockMesh\n')
        f_handle.write('decomposePar\n')
        f_handle.write('mpirun -n 64 simpleFoam -parallel\n')
        f_handle.write('reconstructPar\n')
        
    
keys = tuple( lhs.param_dict.keys() )
line = ''
for key in keys:
    line = line + key + ' '
#
line = line.strip()
print(line)
#
for i in range(runs):
    line = ''
    
    for key in keys:
        val = lhs.param_dict[key]['values'][i]
        line = line + str(val) + ' '
    
    line = line.strip()
    
    print(line)
        
    
    
