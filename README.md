
# Finding Critical Points with Marching Tetrahedrons

In this project, we use marching tetrahedrons algorithm to find the critical points, which is only a rough search of the critical points, without considering the degradation.(It's not accurate. Give up)

Now the project changes to the way of speed interpolation point coordinates.
GPU is used to speed up.

### Note: 
1. Many similar points may be found. 
2. Only VTK files with tetrahedral "UNSTRUCTURED_GRID" are supported
##### Next, we will filter them to remove similar points.

### build:
    make clean
    make

### run:
    ./CriticalPoints vtk_file_path [1:Use Gpu]
