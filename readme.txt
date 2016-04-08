A simple voxelizer that converts ".obj" meshes files to uniform voxels using conservative voxelization.
To Compile:
Use CMake or create a new visual studio solution that includes all source ".cpp" files.

To use
>vox.exe number_of_meshes mesh_files.obj -r resolution 3
or
>vox.exe number_of_meshes mesh_files.obj -l cell_size 3

number_of_meshes: integer specifying how many meshes you would like to voxelizer

mesh_files a list of ".obj" mesh files.

-r specifies the number of voxels to put along the longest axis of the object.

-l specifies the side length of a voxel.

The last number 3 specifies dimension.

For example, to voxelize the included example mesh file
>vox.exe 1 meshes/bunny.obj -r 64 3

Output:
voxel_out.obj a mesh file for quick preview of the voxelized mesh
vox_out.txt
First line is number f voxels in x,y,z axes.
The rest of file is a list of integers where "0" stands for empty voxel, "1", "2" stands for filled voxels.
The integers are listed in z,y,x order.
i.e., the first row of integers correspond to a row of voxels in (0,0,:).