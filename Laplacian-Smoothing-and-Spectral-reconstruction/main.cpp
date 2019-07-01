//
//  main.cpp
//  
//
//  Created by Benjamin Barral on 05/02/2019.
//

#include <stdio.h>



#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <igl/writeOBJ.h>
#include <igl/file_exists.h>
#include <Eigen/Geometry>
#include <cmath>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Eigen/Core>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include<Eigen/IterativeLinearSolvers>

#include <iterator>
#include <random>

using namespace std;
using namespace Eigen;
using namespace Spectra;

string mesh_file_name = "../../../example_meshes/bunny.obj";
bool cotan_form = true;
bool curvatureModeLB = true;
bool viewSpectralReconstructionMode = false;
bool meshSmoothing = false;
bool explicitSmoothingMode = true;

bool normalizeColors = true;

bool cotFormula = true;

int maxKSpec = 3;
// Change this to one of the three modes to perform one of the three courseworks kinds of task
/* For all the other tasks, stick to TASK_1 and change the parameters in the UI :
 - For task 3 : change to 1 on the 'Noise mode' slider parameter
 - For task 4 : change to 1 on the 'Subsampling mode' slider parameter
 - For task 6 : change to 1 on the 'ICP mode' slider parameter in order to perform global point-to-plane 
*/

void computeMeanCurvatureDiscrete(const MatrixXd& V, const vector<vector<int>>& adjacencyVectors,
                                  VectorXd& H, MatrixXd& normals, const int& num_vertices)
{
    int valence;
    vector<int> adjacencyVector;
    Vector3d currVec,neighborVec;
    Vector3d laplacVec;
    Vector3d normal;
    double mean_curv;
    for (int i = 0; i < num_vertices; i++){
        laplacVec = VectorXd::Zero(3);
        currVec = V.row(i);
        adjacencyVector =adjacencyVectors.at(i);
        valence = adjacencyVector.size();
        for (auto const& j : adjacencyVector){
            neighborVec = V.row(j);
            laplacVec += (neighborVec - currVec);
        }
        laplacVec /= valence;
        mean_curv = laplacVec.norm() / 2.;
        normal = -laplacVec / (2. * mean_curv);
        //normal = -laplacVec / laplacVec.norm();
        //mean_curv = - laplacVec[0] / (2 * normal[0]);
        normals.row(i) = normal;
        H(i) = mean_curv;
    }
}

void computeLaplaceBeltramiMatrixUniform(const MatrixXd& V, const vector<vector<int>>& adjacencyVectors,
                            SparseMatrix<double>& LB, const int& num_vertices)
{
    int valence;
    vector<int> adjacencyVector;
    for (int i = 0; i<num_vertices; i++){
        adjacencyVector = adjacencyVectors.at(i);
        valence = adjacencyVector.size();
        for (auto const& j : adjacencyVector){
            LB.coeffRef(i, j) = 1.0 / double(valence);
            LB.coeffRef(i, i) = -1.0;
        }
    }
}

void computeLaplaceBeltramiMatrixCotan(const MatrixXd& V, const MatrixXi& F,
                                       const vector<vector<int>>& adjacencyVectors,
                                       const vector<vector<int>>& adjacencyFaceVectors,
                                         SparseMatrix<double>& C,
                                       SparseMatrix<double>& Minv,
                                       SparseMatrix<double>& M,
                                       SparseMatrix<double>& M_minus_half,
                                       const int& num_vertices, const int& num_faces)
{
    double area;
    vector<int> adjacencyVector, adjacencyFaceVector;
    vector <int> adjVert;
    Vector3d currVer, neighborVer1, neighborVer2;
    for (int i = 0; i<num_vertices; i++){
        adjacencyVector = adjacencyVectors.at(i);
        adjacencyFaceVector = adjacencyFaceVectors.at(i);
        area = 0.;
        currVer = V.row(i);
        // Compute area and angles
        for (auto const& f : adjacencyFaceVector){
            adjVert.clear();
            for (int k = 0; k < 3; k++){
                int v = F(f,k);
                if (v!=i){
                    adjVert.push_back(v);
                }
                //else cout << "cx,io" << endl;
            }
            if (adjVert.size() > 2) cout << "adVert size > 2 wshhhh" << endl;
            int i1 = adjVert.at(0);
            neighborVer1 = V.row(i1);
            int i2 = adjVert.at(1);
            neighborVer2 = V.row(i2);
            Vector3d edge_II1 = neighborVer1 - currVer;
            Vector3d edge_II2 = neighborVer2 - currVer;
            Vector3d edge_I1I2 = neighborVer2 - neighborVer1;
            area += ((edge_II1.cross(edge_II2)).norm() / 2.) / 3.;
            edge_II1.normalize();// /= edge_II1.norm();
            edge_II2.normalize(); // /= edge_II2.norm();
            edge_I1I2.normalize(); // /= edge_I1I2.norm();
            double angle1 = acos(edge_II2.dot(edge_I1I2));
            double angle2 = acos(edge_II1.dot(-edge_I1I2));
            double cot1 = tan(M_PI_2 - angle1); // 1. / tan(angle1);
            double cot2 = tan(M_PI_2 - angle2); // 1. / tan(angle2);
            //if (cot1 < 0 || cot2 < 0) cout << "Cot negative wsh" << endl;
            C.coeffRef(i, i1) += cotFormula? cot1 : angle1; //cot1;
            C.coeffRef(i, i2) += cotFormula? cot2 : angle2; //cot2;
            C.coeffRef(i, i) -= cotFormula? cot1 : angle1; //cot1;
            C.coeffRef(i, i) -= cotFormula? cot2 : angle2; //cot2;
        }
        Minv.coeffRef(i, i) = 1. / (2. * area);
        M.coeffRef(i, i) = 2. * area;
        M_minus_half.coeffRef(i,i) = sqrt(1. / (2. * area));
    }
}

void computeMeanCurvatureLaplaceBeltrami(const SparseMatrix<double>& LB, const MatrixXd& V,
                                         VectorXd& H, MatrixXd& normals, const int& num_vertices,
                                         const vector < vector <int> >& adjacencyList)
{
    MatrixXd laplacProduct = LB * V;
    Vector3d laplacVec, currVec;
    Vector3d normal;
    H = VectorXd::Zero(num_vertices);
    double mean_curv;
    // Compute the barycenter of the mesh
    Vector3d barycenter = Vector3d::Zero(3);
    for (int i = 0; i < num_vertices; i++){
        barycenter +=  V.row(i);
    }
    barycenter /= double(num_vertices);
    for (int i = 0; i < num_vertices; i++){
        laplacVec = laplacProduct.row(i);
        //cout << "Laplac vec = " << laplacVec << endl;
        mean_curv = laplacVec.norm() / 2.;
        normal = -laplacVec / (2. * mean_curv);
        currVec = V.row(i);
        
        normals.row(i) = normal;
        H(i) = mean_curv;
    }
}

vector < vector<int> >  computeVertexFaceAdjacency(const MatrixXd& V, const MatrixXi& F,
                                                   const int& num_faces, const int& num_vertices)
{
    vector < vector<int> > adjVertFaceList;
    adjVertFaceList.resize(num_vertices); // , vector< int> ());
    for (int f = 0; f < num_faces; f++){
        for (int k = 0; k < 3; k++){
            int i = F(f,k);
            adjVertFaceList[i].push_back(f);
        }
    }
    return adjVertFaceList;
}

void computeGaussCurvatureAngleDeficit(const MatrixXd& V, const MatrixXi& F, VectorXd& G,
                                       const int& num_vertices,
                                       const vector<vector<int>>& adjacencyFaceVectors,
                                       const bool& cotanMode)
{
    int numAdjFaces;
    double faceAngle;
    vector <int> adjFaces;
    vector <int> adjVert;
    Vector3d currVer, neighborVer1, neighborVer2;
    double area; // for area normalization
    for (int i = 0; i < num_vertices; i++){
        currVer = V.row(i);
        faceAngle = 0;
        adjFaces = adjacencyFaceVectors.at(i);
        numAdjFaces = adjFaces.size();
        if (cotanMode) area = 0;
        for (const auto f : adjFaces){
            adjVert.clear();
            for (int k = 0; k < 3; k++){
                int v = F(f,k);
                if (v!=i){
                    adjVert.push_back(v);
                }
            }
            neighborVer1 = V.row(adjVert.at(0));
            neighborVer2 = V.row(adjVert.at(1));
            Vector3d edge1 = neighborVer1 - currVer;
            Vector3d edge2 = neighborVer2 - currVer;
            if (cotanMode) area += ((edge1.cross(edge2)).norm() / 2.) / 3.;
            edge1/= edge1.norm();
            edge2/= edge2.norm();
            faceAngle += acos( edge1.dot(edge2) );
        }
        faceAngle = 2 * M_PI - faceAngle;
        if (cotanMode) faceAngle /= area;
        faceAngle;
        G(i) = faceAngle;
    }
}

void computeSpectralReconstruction(MatrixXd& VSpec, const MatrixXd& V, const int& num_vertices,
                                   SparseGenMatProd<double>& op, int& kSpec,
                                   MatrixXd& eigenVectors,
                                   const bool& computeEigenDecompo)
{
    VSpec = MatrixXd::Zero(num_vertices, 3);
    int kEffective = kSpec;
    if (computeEigenDecompo){
        GenEigsSolver< double, SMALLEST_MAGN, SparseGenMatProd<double> > eigs(&op, kSpec, min(20 * kSpec + 1, num_vertices));
        // Initialize and compute
        eigs.init();
        int nconv = eigs.compute();
        eigenVectors = eigs.eigenvectors().real();
        kEffective = eigenVectors.cols();
        cout << "num eigen values found : " << kEffective << endl;
        kSpec = kEffective;
    }
    for (int k = 0; k < kSpec; k++){
        //(es.eigenvalues()[k]).real();
        VectorXd eigenVecK = eigenVectors.col(k);
        for (int i = 0; i<3; i++){
            VSpec.col(i) += eigenVecK.dot(V.col(i)) * eigenVecK;
        }
    }
}

void computeSpectralReconstructionCotan(MatrixXd& VSpec, const MatrixXd& V,
                                        const SparseMatrix<double>& M_minus_half,
                                        const SparseMatrix<double>& M,
                                        const int& num_vertices,
                                   SparseSymMatProd<double>& op, int& kSpec,
                                   MatrixXd& eigenVectors,
                                   const bool& computeEigenDecompo)
{
    VSpec = MatrixXd::Zero(num_vertices, 3);
    int kEffective = kSpec;
    if (computeEigenDecompo){
        SymEigsSolver< double, SMALLEST_MAGN, SparseSymMatProd<double> > eigs(&op, kSpec, min(20 * kSpec + 1, num_vertices));
        // Initialize and compute
        eigs.init();
        int nconv = eigs.compute();
        eigenVectors = eigs.eigenvectors().real();
        kEffective = eigenVectors.cols();
        cout << "num eigen values found : " << kEffective << endl;
        kSpec = kEffective;
    }
    for (int k = 0; k < kSpec; k++){
        VectorXd eigenVecK = M_minus_half * eigenVectors.col(k);
        //eigenVecK /= eigenVecK.norm();
        for (int i = 0; i<3; i++){
            VSpec.col(i) += eigenVecK.dot(M * V.col(i)) * eigenVecK;
            VSpec.col(i) += eigenVecK.dot(M * V.col(i)) * eigenVecK;
        }
    }
}

SparseMatrix<double> computeSmoothMatrixExplicit(const int& num_vertices, const SparseMatrix<double>& LB,
                                                 const double& lambda)
{
    SparseMatrix<double> S(num_vertices,num_vertices);
    S.setIdentity();
    S = S + lambda * LB;
    return S;
}

SparseMatrix<double> computeSmoothMatrixImplicit(const int& num_vertices, const SparseMatrix<double>& LB,
                                                 const double& lambda)
{
    SparseMatrix<double> S(num_vertices,num_vertices);
    S.setIdentity();
    S = S - lambda * LB;
    return S;
}

void smoothMeshExplicit(MatrixXd& prev_V_smooth, MatrixXd& curr_V_smooth,
                const double& lambda_explicit,const SparseMatrix<double>& LB,
                const int& num_vertices, const int& numIterations)
{
    SparseMatrix<double> smooth_matrix = computeSmoothMatrixExplicit(num_vertices, LB, lambda_explicit);
    for (int i = 0; i < numIterations; i++){
        curr_V_smooth = smooth_matrix * prev_V_smooth;
        prev_V_smooth = curr_V_smooth;
    }
}

void smoothMeshExplicitUntilConvergence(MatrixXd& prev_V_smooth, MatrixXd& curr_V_smooth,
                        const double& lambda_explicit,const SparseMatrix<double>& LB,
                        const int& num_vertices, const double& epsilon, const int& maxIterations)
{
    bool converged = false;
    int cout = 0;
    SparseMatrix<double> smooth_matrix = computeSmoothMatrixExplicit(num_vertices, LB, lambda_explicit);
    double change = DBL_MAX;
    while (change > epsilon){
        curr_V_smooth = smooth_matrix * prev_V_smooth;
        
        prev_V_smooth = curr_V_smooth;
    }
}


void smoothMeshImplicit(MatrixXd& prev_V_smooth, MatrixXd& curr_V_smooth,
                        const double& lambda_implicit,const SparseMatrix<double>& LB,
                        const int& num_vertices, const int& numIterations)
{
    SparseMatrix<double> smooth_matrix_impl = computeSmoothMatrixImplicit(num_vertices, LB, lambda_implicit);
    BiCGSTAB< SparseMatrix<double> > solverImpl;
    solverImpl.compute(smooth_matrix_impl);
    for (int i = 0; i < numIterations; i++){
        curr_V_smooth = solverImpl.solve(prev_V_smooth);
        prev_V_smooth = curr_V_smooth;
    }
}

void computeBoundingBoxSize(const MatrixXd& V, const int& num_vertices, Vector3d& boundingSizes, double& diagSize)
{
    Vector3d maxPoint = DBL_MIN * Vector3d::Ones(3);
    Vector3d minPoint = DBL_MAX * Vector3d::Ones(3);
    for (int i = 0; i < num_vertices; i++){
        for (int j = 0; j < 3; j++){
            double coord = V(i,j);
            if (coord > maxPoint(j)) maxPoint(j) = coord;
            if (coord < minPoint(j)) minPoint(j) = coord;
        }
    }
    boundingSizes = maxPoint - minPoint;
    diagSize = boundingSizes.norm();
}

void addNoise(const MatrixXd& V0, MatrixXd& noisyV, const int& num_vertices, const double& noise_ratio, const Vector3d& boundingBoxSizes)
{
    noisyV = V0;
    std::default_random_engine generator;
    normal_distribution<double>* dist = new normal_distribution<double>[3];
    double sigma;
    for (int j = 0; j < 3; j++){
        sigma = noise_ratio * boundingBoxSizes[j];
        dist[j] = normal_distribution<double>(0,sigma);
    }
    //distX(generator);
    for (int i = 0; i < num_vertices; i++){
        for (int j = 0; j < 3; j++){
            noisyV(i,j) += dist[j](generator);
        }
    }
    delete dist;
}

int main(int argc, char *argv[])
{
    // Init the viewer
    igl::opengl::glfw::Viewer viewer;
    
    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    MatrixXd V;//,N,TC,FTC,FN;
    MatrixXi F;
    igl::readOBJ(mesh_file_name,V,F);
    //igl::readOFF(mesh_file_name,V,F);
    
    //viewer.data().set_vertices(V);
    viewer.data().set_mesh(V, F);
    
    int num_vertices = V.rows();
    int num_faces = F.rows();
    cout << "Num vertices : " << num_vertices << endl;
    cout << "Num faces : " << num_faces << endl;
    
    vector< vector< int > > adjacencyList;
    igl::adjacency_list(F, adjacencyList);
    vector< vector< int > > adjacencyFaceList = computeVertexFaceAdjacency(V,F, num_faces, num_vertices);
    
    // CURVATURES
    /*Mean curvature*/
    // Discrete version
    MatrixXd normals_discrete(num_vertices,3);
    VectorXd mean_curvatures_discrete(num_vertices);
    computeMeanCurvatureDiscrete(V, adjacencyList, mean_curvatures_discrete, normals_discrete, num_vertices);
    // Laplace Beltrami version
    //Uniform
    MatrixXd normals_LB_uniform(num_vertices,3);
    VectorXd mean_curvatures_LB_uniform(num_vertices);
    SparseMatrix<double> laplace_beltrami_uniform(num_vertices,num_vertices);
    computeLaplaceBeltramiMatrixUniform(V, adjacencyList, laplace_beltrami_uniform, num_vertices);
    computeMeanCurvatureLaplaceBeltrami(laplace_beltrami_uniform, V, mean_curvatures_LB_uniform, normals_LB_uniform, num_vertices, adjacencyList);
    //Cotan
    MatrixXd normals_LB_cotan(num_vertices,3);
    VectorXd mean_curvatures_LB_cotan(num_vertices);
    SparseMatrix<double> laplace_beltrami_cotan(num_vertices,num_vertices);
    SparseMatrix<double> laplace_beltrami_Minv(num_vertices,num_vertices),
    laplace_beltrami_M(num_vertices,num_vertices),
    M_minus_half(num_vertices,num_vertices),
    laplace_beltrami_C(num_vertices,num_vertices);
    computeLaplaceBeltramiMatrixCotan(V, F, adjacencyList, adjacencyFaceList, laplace_beltrami_C,
                                      laplace_beltrami_Minv, laplace_beltrami_M, M_minus_half,
                                      num_vertices, num_faces);
    
    
    // DEBUG : use IGL laplace
    SparseMatrix<double> laplace_beltrami_IGL,M_IGL, M_inv_IGL;
    igl::cotmatrix(V,F,laplace_beltrami_IGL);
    
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC, M_IGL);
    igl::invert_diag(M_IGL,M_inv_IGL);
    
    laplace_beltrami_cotan = laplace_beltrami_Minv * laplace_beltrami_C;
    computeMeanCurvatureLaplaceBeltrami(laplace_beltrami_cotan, V, mean_curvatures_LB_cotan, normals_LB_cotan, num_vertices, adjacencyList);
    
    
    /*Gauss curvature*/
    VectorXd gauss_curvatures_uniform(num_vertices);
    computeGaussCurvatureAngleDeficit(V, F, gauss_curvatures_uniform, num_vertices, adjacencyFaceList, false);
    
    VectorXd gauss_curvatures_cotan(num_vertices);
    computeGaussCurvatureAngleDeficit(V, F, gauss_curvatures_cotan, num_vertices, adjacencyFaceList, true);
    
    MatrixXd color_mean_curvatures_LB_uniform(num_vertices,3);
    igl::jet(mean_curvatures_LB_uniform, normalizeColors, color_mean_curvatures_LB_uniform);
    MatrixXd color_mean_curvatures_LB_cotan(num_vertices,3);
    igl::jet(mean_curvatures_LB_cotan, normalizeColors, color_mean_curvatures_LB_cotan);
    MatrixXd color_mean_curvatures_discrete(num_vertices,3);
    igl::jet(mean_curvatures_discrete, normalizeColors, color_mean_curvatures_discrete);
    
    MatrixXd color_gauss_curvatures_uniform(num_vertices,3);
    igl::jet(gauss_curvatures_uniform, normalizeColors, color_gauss_curvatures_uniform);
    MatrixXd color_gauss_curvatures_cotan(num_vertices,3);
    igl::jet(gauss_curvatures_cotan, normalizeColors, color_gauss_curvatures_cotan);
    
    MatrixXd normals = curvatureModeLB ? (cotan_form ? normals_LB_cotan: normals_LB_uniform)
    : normals_discrete;
    MatrixXd color_mean_curvatures = curvatureModeLB ? (cotan_form ? color_mean_curvatures_LB_cotan : color_mean_curvatures_LB_uniform)
    : color_mean_curvatures_discrete;
    MatrixXd color_gauss_curvatures = cotan_form ? color_gauss_curvatures_cotan : color_gauss_curvatures_uniform;
    
    MatrixXd C = color_mean_curvatures;
    //viewer.data().set_normals(normals);
    viewer.data().set_colors(C);
    
    // SPECTRAL RECONSTRUCTION
    int kSpec = maxKSpec;
    MatrixXd VSpecUniform = MatrixXd::Zero(num_vertices, 3);
    MatrixXd VSpecCotan = MatrixXd::Zero(num_vertices, 3);
    MatrixXd eigenVectorsCotan, eigenVectorsUniform;
    // Reconstruction from Uniform Laplace
    SparseMatrix<double> laplace_beltrami_cotan_sym(num_vertices,num_vertices);
    laplace_beltrami_cotan_sym = M_minus_half * laplace_beltrami_C * M_minus_half;
    SparseGenMatProd<double> opUniform(laplace_beltrami_uniform);
    SparseSymMatProd<double> opCotan(laplace_beltrami_cotan_sym);
    
    /*computeSpectralReconstructionCotan(VSpecCotan, V, M_minus_half, laplace_beltrami_M, num_vertices, opCotan, maxKSpec, eigenVectorsCotan, true);
    computeSpectralReconstruction(VSpecUniform, V, num_vertices, opUniform, maxKSpec, eigenVectorsUniform, true);
    viewer.data().clear();
    viewer.data().set_mesh(cotan_form ? VSpecCotan : VSpecUniform, F);*/
    
    // MESH SMOOTHING
    MatrixXd curr_V_smooth = V; //(num_vertices,3);
    MatrixXd prev_V_smooth = V;
    double lambda_explicit = 0.2;
    double lambda_implicit = 0.2;
    
    // NOISE
    double noise_ratio = 0.01;
    int noiseModeInt = 0;
    Vector3d boundingSizes;
    double diagSize;
    computeBoundingBoxSize(V, num_vertices, boundingSizes, diagSize);
    MatrixXd noisyV;
    addNoise(V, noisyV, num_vertices, noise_ratio, boundingSizes);
    
    // Menu UI parameters
    bool transformUserInput = false;
    bool updateShading = false;
    bool updateMesh = false;
    bool reset = false;
    int gaussianCurvatureMode = 0;
    int cotanMode = int(cotan_form);
    int laplaceMode = int(curvatureModeLB);
    int implicitMode =  int(!explicitSmoothingMode);
    int numIterations = 1;
    int viewMode = 2;
    
    menu.callback_draw_custom_window = [&]()
    {
        updateShading = false;
        updateMesh = false;
        transformUserInput = false;
        reset = false;
        bool changeLaplaceForm = false;
        bool addNoiseAndReset = false;
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(800, 300), ImGuiSetCond_FirstUseEver);
        ImGui::Begin(
                     "MyProperties", nullptr,
                     ImGuiWindowFlags_NoSavedSettings
                     );
        
        if(ImGui::SliderInt("Mean curvature - Gauss curvature", &gaussianCurvatureMode, 0, 1)){
            updateShading = true;
        }
        if(ImGui::SliderInt("Uniform mode - Cotan mode", &cotanMode, 0, 1)){
            if (cotanMode == 1) cotan_form = true;
            else cotan_form = false;
            updateShading = true;
            changeLaplaceForm = true;
        }
        if(ImGui::SliderInt("Discrete mode - Laplace mode", &laplaceMode, 0, 1)){
            if (laplaceMode == 1) curvatureModeLB = true;
            else curvatureModeLB = false;
            updateShading = true;
        }
        if (ImGui::SliderInt("Curvature - Spectral reconstruction - Smoothing", &viewMode, 0, 2)){
            if (viewMode == 0){
                updateShading = true;
                viewSpectralReconstructionMode = false;
                meshSmoothing = false;
            }
            else if (viewMode == 1){
                updateShading = false;
                viewSpectralReconstructionMode = true;
                updateMesh = true;
                meshSmoothing = false;
            }
            else{
                updateShading = false;
                viewSpectralReconstructionMode = false;
                updateMesh = true;
                meshSmoothing = true;
            }
        }
        bool spectralReconstruction = viewSpectralReconstructionMode &&
        (ImGui::SliderInt("Number of eigen vectors", &kSpec, 1, maxKSpec) || changeLaplaceForm);
        if(spectralReconstruction){
            if (cotan_form  || !(changeLaplaceForm))
                computeSpectralReconstructionCotan(VSpecCotan, V, M_minus_half, laplace_beltrami_M, num_vertices, opCotan, kSpec, eigenVectorsCotan, false);
            if (!cotan_form || !(changeLaplaceForm))
                computeSpectralReconstruction(VSpecUniform, V, num_vertices, opUniform, kSpec, eigenVectorsUniform, false);
            updateMesh = true;
        }
        if (ImGui::SliderInt("Explicit mode - Implicit mode", &implicitMode, 0, 1)){
            if (implicitMode == 1){
                explicitSmoothingMode = false;
            }
            else{
                explicitSmoothingMode = true;
            }
        }
        
        ImGui::InputInt("Number of smoothing iterations", &numIterations);
        if(ImGui::Button("Smoothing iteration")){
            updateMesh = true;
            if(meshSmoothing){
                if (explicitSmoothingMode){
                    smoothMeshExplicit(prev_V_smooth,curr_V_smooth, lambda_explicit,
                                       cotan_form? laplace_beltrami_cotan : laplace_beltrami_uniform, num_vertices, numIterations);
                }
                else{
                    smoothMeshImplicit(prev_V_smooth,curr_V_smooth, lambda_implicit,
                                       cotan_form? laplace_beltrami_cotan : laplace_beltrami_uniform, num_vertices, numIterations);
                }
            }
        }
        
        ImGui::InputDouble("Lambda explicit smoothing", &lambda_explicit);
        ImGui::InputDouble("Lambda implicit smoothing", &lambda_implicit);
        if(ImGui::Button("Reset")){
            updateMesh = true;
            if (noiseModeInt == 0){
            curr_V_smooth = V;
            prev_V_smooth = V;
            }
            else addNoiseAndReset = true;
        }
        bool addNoiseTest = false;
        if (ImGui::SliderInt("Add noise", &noiseModeInt, 0, 1)){
            updateMesh = true;
            addNoiseAndReset = true;
            if (noiseModeInt == 0) noisyV = V;
            else addNoiseTest = true;
        }
        bool addNoiseTest2 = (ImGui::InputDouble("Amount of noise (ratio of the bounding box dimensions", &noise_ratio));
        addNoiseTest = addNoiseTest || addNoiseTest2;
        
        if (addNoiseTest){
            if (noise_ratio !=0) addNoise(V, noisyV, num_vertices, noise_ratio, boundingSizes);
            addNoiseAndReset = true;
            updateMesh = true;
        }
        
        if (updateShading && viewMode == 0){
            if (gaussianCurvatureMode){
                color_gauss_curvatures = cotan_form ? color_gauss_curvatures_cotan : color_gauss_curvatures_uniform;
                C = color_gauss_curvatures;
            }
            else{
                color_mean_curvatures = curvatureModeLB ? (cotan_form ? color_mean_curvatures_LB_cotan : color_mean_curvatures_LB_uniform)
                : color_mean_curvatures_discrete;
                C = color_mean_curvatures;
            }
            
            normals = curvatureModeLB ? (cotan_form ? normals_LB_cotan: normals_LB_uniform)
            : normals_discrete;
            
            viewer.data().clear();
            viewer.data().set_mesh(V, F);
            //viewer.data().set_normals(normals);
            viewer.data().set_colors(C);
        }
        if (addNoiseAndReset){
            curr_V_smooth = noisyV;
            prev_V_smooth = noisyV;
        }
        if (viewSpectralReconstructionMode && updateMesh){
            viewer.data().clear();
            viewer.data().set_mesh(cotan_form ? VSpecCotan : VSpecUniform, F);
        }
        if (meshSmoothing && updateMesh){
            viewer.data().clear();
            viewer.data().set_mesh(curr_V_smooth, F);
        }
 
        ImGui::End();
    };
    
    // Call GUI
    viewer.launch();
 
}


