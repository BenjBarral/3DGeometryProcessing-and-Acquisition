# An implementation of the Laplacian matrix computation, with different applications.

![](media/giphy.gif)

An application that allows to align two point clouds of the same object in an interactive way, using [ICP](http://www-evasion.inrialpes.fr/people/Franck.Hetroy/Teaching/ProjetsImage/2007/Bib/besl_mckay-pami1992.pdf) [1].
Implementation of the Iterative Closest Point algorithm, and the [Point-To-Plane](https://www.comp.nus.edu.sg/~lowkl/publications/lowk_point-to-plane_icp_techrep.pdf) [2] variant.


## Usage
### Dependencies
Libraries needed : 
- [LibIGL](https://github.com/libigl/libigl) . Put the libigl source code folder two levels above this repository.
- [Spectra](https://spectralib.org/download.html) (for large scale eigen value computations). Put the Spectra source code folder one level above this repository, or change the ```include_directories``` and in the CMake to the corresponding paths.

### Building and running the code
Build using CMake and the CMakeLists.txt provided.
Run the ```Laplacian_processing_viewer``` target.
Change the line ```string mesh_file_name``` with the right filenames of your point clouds, in main.cpp.

### Using the interface
The libigl GUI menu has several options you can manipulate : 
- Change the "ICP mode" slider in order to perfom point-to-point ICP (0) or point-to-plane ICP (1).
- Click the "ICP step" button to perform one registration step of the algorithm.
- Click the "ICP registration" button to perform the algorithm until convergence.
- Change the values of "Max num iterations" and "Error threshold" to define the convergence stopping conditions.
- You can roughly pre-align the point clouds by manipulating the "Rotation slider" and "Translation slider".


## References
[1] Paul J. Besl et N.D. McKay. A Method for Registration of 3-D Shapes. IEEE Trans. on Pattern Analysis and Machine Intelligence.

[2] LOW, K.-L. 2004. Linear least-squares optimization for point-to- plane icp surface registration.
Tech. rep., Chapel Hill, University of North Carolina.
