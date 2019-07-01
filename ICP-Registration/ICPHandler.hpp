//
//  ICPHandler.hpp
//  3DA_assignment1_bin
//
//  Created by Benjamin Barral on 11/02/2019.
// Description : Class that manages the ICP algorithm for the registration of one point-cloud to another.

#ifndef ICPHandler_hpp
#define ICPHandler_hpp

#include <igl/readPLY.h>
#include <ANN/ANN.h>

using namespace Eigen;
using namespace std;

struct Subsampling{
    bool subsample_mode_;
    double subsample_rate_;
    Subsampling(){
        subsample_mode_ = false;
        subsample_rate_ = 0.0005;
    }
    Subsampling(const bool& subsampleMode, const float& subsampleRate){
        subsample_mode_ = subsampleMode;
        subsample_rate_ = subsampleRate;
    }
};

enum ICPMode{
    POINT_TO_POINT = 0,
    POINT_TO_PLANE
};

class ICPHandler{
private:
    MatrixXd verticesA_, verticesB_; // Register verticesB_ to verticesA_
    int nA_,nB_; // Number of points
    Subsampling subsampling_;
    vector<int> icp_points_indices; // The indices of the points in verticesB_ used for ICP (c.f. subsampling)
    vector<Vector3d> icp_pointsA_, icp_pointsB_; // Correspondences (nearest neighbor);
    int num_ICP_points_;
    ANNkd_tree* kd_tree_;
    Vector3d barycenterA_, barycenterB_;
    vector<Vector3d> recentered_pointsA_, recentered_pointsB_;
    Matrix4d rigid_transform_;
    double error_; // Square root of Mean-square registration error
    vector<double> errors_in_steps_; // to study behaviors of algorithm
    double convergence_threshold_;
    int max_iterations_;
    int num_iterations_;
    // For Point-to-plane ICP
    MatrixXd normalsA_;
    vector<Vector3d> icp_normalsA_;
    ICPMode icp_mode_;
    MatrixXd matrix_A_;
    VectorXd vector_b_;
    
    void SubsampleUniformly();
    void ComputeKDTree();
    void ComputeCorrespondences();
    void ComputeBarycentersAndRecenterPoints();
    void ComputeRigidTransformation();
    void ComputeLinearApproximationMatrices(); // for point to plane
    void TransformPoints();
    void ComputeError();
    void PrintStepError() const;
    void PrintConvergenceState() const;
    
public:
    ICPHandler();
    ICPHandler(const MatrixXd& verticesA, const MatrixXd& verticesB, const Subsampling& subsampling,
               const double& convergenceThreshold, const ICPMode& icpMode);
    //ICPHandler& operator=(const ICPHandler& other);
    void set_verticesA(const MatrixXd& verticesA);
    void set_verticesB(const MatrixXd& verticesB);
    void set_subsampling(const Subsampling& subsampling);
    void set_convergence_threshold(const double& convergenceThreshold);
    void set_max_iterations(const int& maxIterations);
    void set_icp_mode(const ICPMode& icpMode);
    MatrixXd verticesB() const;
    MatrixXd normalsA() const;
    vector<double> errors_in_steps() const;
    void ComputeIcpStep();
    void Reset();
    void ComputeIcpRegistration();
    void ComputeNormals();
};

#endif /* ICPHandler_hpp */
