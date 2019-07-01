# An implementation of the Laplace Beltrami computation, with different applications.

![](media/giphy.gif)

Implementation of [Laplacian Mesh Processing](https://people.eecs.berkeley.edu/~jrs/meshpapers/Sorkine.pdf).
Implementation of [Laplacien Mesh Spectral Reconstruction](https://members.loria.fr/Bruno.Levy/papers/Laplacian_SMI_2006.pdf).
Implementation of [Implicit Laplacian Smoothing](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.2.5347&rep=rep1&type=pdf).

Coded in C++, with the LibIGL library.

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
The libigl GUI menu I built has several options you can manipulate : 
- Play with the "Curvature - Spectral reconstruction - Smoothing" slider to change the visualization mode : 
0 : visualize the different curvatures - 1 : visualize spectral reconstruction - 2 : visualize mesh smoothing
- If in  "Spectral reconstruction" mode, change the values of "Number of eigen vectors" to set the number of eigen vectors used for reconstruction
- If in "Smoothing" mode, choose between "Explicit mode" and "Implicit mode", to perform explicit or implicit smoothing iterations.  
- If in "Smoothing" mode : change the value of "Lambda explicit/implicit smoothing" to change the value of the smoothing parameter. Change the "Number of smoothing iterations" using the +/- buttons.
- If in "Smoothing" mode : click "Smoothing iteration" to perform an iteration of the smoothing algorithm. Press "Reset" to go back to the initial state.
- You can add Gaussian noise to the vertex positions before smoothing, by changing the value of "Amount of noise", and setting the "Add noise" slider to 1.
