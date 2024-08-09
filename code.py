#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from PIL import Image
from matplotlib.image import imread


# In[2]:


# defining the projection matrices
pm0 = np.array([[776.649963, -298.408539, -32.048386, 993.1581875], 
                [132.852554, 120.885834, -759.210876, 1982.174000], 
                [0.744869, 0.662592, -0.078377, 4.629312012]])
pm1 = np.array([[431.503540, 586.251892, -137.094040, 1982.053375], 
                [23.799522, 1.964373, -657.832764, 1725.253500], 
                [-0.321776, 0.869462, -0.374826, 5.538025391]])
pm2 = np.array([[-153.607925, 722.067139, -127.204468, 2182.4950], 
                [141.564346, 74.195686, -637.070984, 1551.185125], 
                [-0.769772, 0.354474, -0.530847, 4.737782227]])
pm3 = np.array([[-823.909119, 55.557896, -82.577644, 2498.20825], 
                [-31.429972, 42.725830, -777.534546, 2083.363250], 
                [-0.484634, -0.807611, -0.335998, 4.934550781]])
pm4 = np.array([[-715.434998, -351.073730, -147.460815, 1978.534875], 
                [29.429260, -2.156084, -779.121704, 2028.892750], 
                [0.030776, -0.941587, -0.335361, 4.141203125]])
pm5 = np.array([[-417.221649, -700.318726, -27.361042, 1599.565000], 
                [111.925537, -169.101776, -752.020142, 1982.983750], 
                [0.542421, -0.837170, -0.070180, 3.929336426]])
pm6 = np.array([[94.934860, -668.213623, -331.895508, 769.8633125], 
                [-549.403137, -58.174614, -342.555359, 1286.971000], 
                [0.196630, -0.136065, -0.970991, 3.574729736]])
pm7 = np.array([[452.159027, -658.943909, -279.703522, 883.495000], 
                [-262.442566, 1.231108, -751.532349, 1884.149625], 
                [0.776201, 0.215114, -0.592653, 4.235517090]])


# In[3]:


# Reading the camera images
s0 = imread('silhouette0.pbm')
s1 = imread('silhouette1.pbm')
s2 = imread('silhouette2.pbm')
s3 = imread('silhouette3.pbm')
s4 = imread('silhouette4.pbm')
s5 = imread('silhouette5.pbm')
s6 = imread('silhouette6.pbm')
s7 = imread('silhouette07.pbm')


# In[4]:


# Initializing all the matrices to zeros
projection_matrix = np.zeros((3, 4, 8))
sil_matrix = np.zeros((582, 780, 8))
im = np.zeros((582, 780, 3, 8))


# In[5]:


# Defining the projection matrices
projection_matrix[:,:,0] = pm0
projection_matrix[:,:,1] = pm1
projection_matrix[:,:,2] = pm2
projection_matrix[:,:,3] = pm3
projection_matrix[:,:,4] = pm4
projection_matrix[:,:,5] = pm5
projection_matrix[:,:,6] = pm6
projection_matrix[:,:,7] = pm7


# In[6]:


# this step is necessary to convert silhouette image to gray scale so that there is no discrepency in the other
# steps. If we do it directly the array dimensions do not match
s0_gray = 0.2125 * s0[:,:,0] + 0.7154 * s0[:,:,1] + 0.0721 * s0[:,:,2]
s1_gray = 0.2125 * s1[:,:,0] + 0.7154 * s1[:,:,1] + 0.0721 * s1[:,:,2]
s2_gray = 0.2125 * s2[:,:,0] + 0.7154 * s2[:,:,1] + 0.0721 * s2[:,:,2]
s3_gray = 0.2125 * s3[:,:,0] + 0.7154 * s3[:,:,1] + 0.0721 * s3[:,:,2]
s4_gray = 0.2125 * s4[:,:,0] + 0.7154 * s4[:,:,1] + 0.0721 * s4[:,:,2]
s5_gray = 0.2125 * s5[:,:,0] + 0.7154 * s5[:,:,1] + 0.0721 * s5[:,:,2]
s6_gray = 0.2125 * s6[:,:,0] + 0.7154 * s6[:,:,1] + 0.0721 * s6[:,:,2]
s7_gray = 0.2125 * s7[:,:,0] + 0.7154 * s7[:,:,1] + 0.0721 * s7[:,:,2]

#Defining the silhouettes matrices
sil_matrix[:,:,0] = s0_gray
sil_matrix[:,:,1] = s1_gray
sil_matrix[:,:,2] = s2_gray
sil_matrix[:,:,3] = s3_gray
sil_matrix[:,:,4] = s4_gray
sil_matrix[:,:,5] = s5_gray
sil_matrix[:,:,6] = s6_gray
sil_matrix[:,:,7] = s7_gray


# In[7]:


# Defining the camera matrices
i0 = imread('cam0.png')
i1 = imread('cam1.png')
i2 = imread('cam2.png')
i3 = imread('cam3.png')
i4 = imread('cam4.png')
i5 = imread('cam5.png')
i6 = imread('cam6.png')
i7 = imread('cam7.png')


# In[8]:


im[:,:,:,0] = i0
im[:,:,:,1] = i1
im[:,:,:,2] = i2
im[:,:,:,3] = i3
im[:,:,:,4] = i4
im[:,:,:,5] = i5
im[:,:,:,6] = i6
im[:,:,:,7] = i7


# In[9]:


x_grid = 5
y_grid = 6
z_grid = 2.5


# In[10]:


volume = x_grid * y_grid * z_grid
no_of_voxels = 10000000
vox_s = round(volume**(1/3)/no_of_voxels**(1/3),3)


# In[11]:


true_voxs = 0
total_voxs = 0
vox_matrix = []
surface_vox_matrix = []
color_matrix = []
fl = 0
prev_vec = []


# In[12]:


for x in np.arange(-x_grid/2, x_grid/2, vox_s):
    for y in np.arange(-y_grid/2, y_grid/2, vox_s):
        for z in np.arange(0, z_grid, vox_s):
            # Initialize decision vector
            dec_v = np.zeros(8)
            # Increment total voxels counter
            total_voxs += 1
            # Create world coordinate vector
            w_col = np.array([x, y, z, 1.0]).T
            # Loop over projection matrices
            for i in range(8):
                # Project world coordinate vector onto image plane
                uv_cor = np.round(np.dot(projection_matrix[:,:,i], w_col) / np.dot(projection_matrix[:,:,i], w_col)[2])
                # Check if projected coordinates are within image bounds
                #if (1<=uv_cor[0]<=780) and (1<=uv_cor[1]<=582):
                if (0<=uv_cor[0]<=sil_matrix.shape[1]-1) and (0<=uv_cor[1]<=sil_matrix.shape[0]-1):
                    # Update decision vector with corresponding silhouette pixel value
                    dec_v[i] = sil_matrix[int(uv_cor[1]), int(uv_cor[0]), i]
                    # Check if all decision vector values are non-zero
                    if np.all(dec_v):
                        # Increment true voxels counter
                        true_voxs += 1
                        # Update voxel matrix with world coordinates
                        vox_matrix.append([x, y, z])
                        # Update color matrix with corresponding color value from image
                        r = im[int(uv_cor[1]), int(uv_cor[0]), 0, 7]
                        g = im[int(uv_cor[1]), int(uv_cor[0]), 1, 7]
                        b = im[int(uv_cor[1]), int(uv_cor[0]), 2, 7]
                        color_matrix.append([r, g, b])
                        # Check if voxel is part of surface
                        if fl == 0:
                            # Increment surface voxel counter
                            surface_vox_matrix.append([x, y, z])
                            fl += 1
                            prev_vec = [x, y, z]
                            continue
                            # Add surface voxels to matrix
                            if (prev_vec[0]==x) and (prev_vec[1]==y):
                                fl += 1
                            else:
                                if fl > 1:
                                    surface_vox_matrix.append(prev_vec)
                                    surface_vox_matrix.append([x, y, z])
                                    f1 = 1
                                else:
                                    # Reset surface voxel counter
                                    surface_vox_matrix.append([x, y, z])
                                    f1 = 1
                            fl = 1
                            # Update previous vector
                            prev_vec = [x, y, z]


# In[1]:


# % Computed the surface voxels
# if(size(surface_vox_matrix, 1) > 1):
#     x_s = surface_vox_matrix(:,1);
#     y_s = surface_vox_matrix(:,2);
#     z_s = surface_vox_matrix(:,3);
#     surface_vox_matrix = [x_s, y_s, z_s];
# end

# % Plot surface voxels
# figure;
# scatter3(surface_vox_matrix(:,1), surface_vox_matrix(:,2), surface_vox_matrix(:,3), 10, color_matrix, 'filled');
# title('Surface Voxels');
# xlabel('X');
# ylabel('Y');
# zlabel('Z');
# axis equal;


# In[14]:


import open3d as o3d


# In[15]:


# Convert color matrix to numpy array

color_matrix = np.array(color_matrix)

vox_mat_ptc = o3d.geometry.PointCloud()
vox_mat_ptc.points = o3d.utility.Vector3dVector(vox_matrix)
vox_mat_ptc.colors = o3d.utility.Vector3dVector(color_matrix.astype(float)/255)


# In[16]:


# Write point cloud to PLY file and display in figure 1

o3d.io.write_point_cloud('dancer_full_colored_2.ply', vox_mat_ptc, write_ascii=True)


# In[17]:


import os
os.system('meshlab dancer_full_colored_2.ply')


# In[18]:


# Write point cloud to PLY file and display in figure 1

o3d.io.write_point_cloud('dancer_full_colored_2.ply', vox_mat_ptc, write_ascii=True)
full_mod = o3d.io.read_point_cloud('dancer_full_colored_2.ply')
o3d.visualization.draw_geometries([full_mod])


# In[19]:


# Create point cloud from surface voxels

sur_vox_mat_ptc = o3d.geometry.PointCloud()
sur_vox_mat_ptc.points = o3d.utility.Vector3dVector(surface_vox_matrix)
sur_vox_mat_ptc.colors = vox_mat_ptc.colors


# In[20]:


# Write point cloud to PLY file and display in figure 2

o3d.io.write_point_cloud('dancer_surf_2.ply', sur_vox_mat_ptc, write_ascii=True)
surf_mod = o3d.io.read_point_cloud('dancer_surf_2.ply')
o3d.visualization.draw_geometries([surf_mod])

