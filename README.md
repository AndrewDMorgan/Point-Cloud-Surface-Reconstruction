# Point-Cloud-Surface-Reconstruction
This is a program I've been working on for a few months. It reconstructs the surface of a point cloud. This is a continuation of my previous year's supercomputing project.

Currently, I'm working on translating the initial python version into c++ to get better performance. I am also working on implementing an Octree data structure for the point cloud to improve memory consumption along with performance. I'm utilizing the MATAR library to also get better memory usage and performance. I plan on running multiple of the steps in parallel utalizing MATAR.

The method I'm using for surface reconstruction relies on the generation of a signed distance field to represent the surface. Additional steps were implemented to improve the quality of the results.

## NM Supercomputing Challange 2025

More information is present in the presentation I've been working on for this project.
Presentation -> https://docs.google.com/presentation/d/1DU2jcye-E8gFEpKhPdf2-ofJ4YAFJ_p3borTgJbGals/edit?usp=sharing


Compile using c++ 17
g++ --std=c++17 main.cpp

