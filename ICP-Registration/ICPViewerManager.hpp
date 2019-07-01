//
//  ICPViewerManager.hpp
//  3DA_assignment1_bin
//
//  Created by Benjamin Barral on 11/02/2019.
//

#ifndef ICPViewerManager_hpp
#define ICPViewerManager_hpp


#include <igl/opengl/glfw/Viewer.h>
#include <stdio.h>
#include "ICPHandler.hpp"
#include "PointCloudData.hpp"
#include "ICPHandler.hpp"
#include <math.h>

using namespace std;
using namespace Eigen;

struct ManualRigidTransformationParams{
    // Handles rigid transformation from user input
    int rot_ind_ = 0, trans_ind_ = 0;
    double rotation_angles_[3] = {0,0,0};
    double translation_[3] = {0,0,0};
    ManualRigidTransformationParams(){
    }
    ManualRigidTransformationParams(const double* angles, const double* trans){
        for (int i = 0; i<3; i++){
            rotation_angles_[i] = angles[i];
            translation_[i] = trans[i];
        }
    }
    void Reset(){
        double rotation_angles_reset[3] = {0,0,0};
        double translation_reset[3] = {0,0,0};
        rot_ind_ = 0;
        trans_ind_ = 0;
        for (int i = 0; i<3; i++){
            rotation_angles_[i] = rotation_angles_reset[i];
            translation_[i] = translation_reset[i];
        }
    }
    void UpdateAngle(const double& angleValue){
        rotation_angles_[rot_ind_] = angleValue * M_PI / 180;
    }
    void UpdateTranslation(const double& transValue){
        translation_[trans_ind_] = transValue;
    }
};

struct NoiseParams{
    bool noise_mode_ = false;
    float sigma_rate_ = 0;
    NoiseParams(){
        
    }
    NoiseParams(const bool& noiseMode, const float& sigmaRate){
        noise_mode_ = noiseMode;
        sigma_rate_ = sigmaRate;
    }
};

struct ICPParams{
    double convergence_threshold_ = 0.0008;
    int max_num_iterations_ = 120;
    Subsampling subsampling_;
    ICPMode icp_mode_;
    ICPParams(){
        
    }
    ICPParams(const double& convergenceThreshold, const Subsampling& subsampling, const ICPMode& icpMode){
        convergence_threshold_ = convergenceThreshold;
        subsampling_ = subsampling;
        icp_mode_ = icpMode;
    }
};

class ICPViewerManager{
private:
    ICPHandler icp_handler_;
    PointCloudData* point_cloud_data_array_;
    Transform <double , 3, Affine > manual_rigid_transform_matrix_;
    int num_point_clouds_;
    int *point_cloud_registration_correspondent_indices; // for task 5 : global registration
    
public:
    ManualRigidTransformationParams manual_rigid_transform_params;
    ICPParams icp_params;
    NoiseParams noise_params;
    int index_pcA_,index_pcB_; // indices of the too point clouds used for local registration in case of multiple point clouds (global registration)
    float point_size_;
    
    ICPViewerManager();
    ICPViewerManager(const MatrixXd* pointClouds, const int& numPC, const ICPMode& icpMode, int* pcCorrespIndices,
                     const bool& randColors);
    void ComputeManualRigidTransformFromParameters();
    void TransformFromParameters();
    void ChangeColorsRandomly(igl::opengl::glfw::Viewer& viewer);
    void AddNoise();
    void DeleteNoise();
    void SetIcpSubsampling();
    void ShowNormals(igl::opengl::glfw::Viewer& viewer);
    void SetIcpMode();
    void SetConvergenceThreshold();
    void SetMaxIterations();
    void ConfigureIcpPair(igl::opengl::glfw::Viewer& viewer);
    void PerformIcpStep();
    void PerformIcpRegistration();
    void InitDisplay(igl::opengl::glfw::Viewer& viewer);
    void UpdateDisplay(igl::opengl::glfw::Viewer& viewer);
    
};

#endif /* ICPViewerManager_hpp */
