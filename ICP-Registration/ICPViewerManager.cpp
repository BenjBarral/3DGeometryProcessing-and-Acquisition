//
//  ICPViewerManager.cpp
//  3DA_assignment1_bin
//
//  Created by Benjamin Barral on 11/02/2019.
//

#include "ICPViewerManager.hpp"

static const Vector3d colorReferencePointCloud = Vector3d::Constant(3, 1.);
static const Vector3d colorRegisteringPointCloud = Vector3d::Constant(3, 0.);
static const Vector3d initialGlobalOffset = Vector3d::Constant(3,0.2);

ICPViewerManager::ICPViewerManager() {
    
}

ICPViewerManager::ICPViewerManager(const MatrixXd* pointClouds, const int& numPC, const ICPMode& icpMode,
                                   int* pcCorrIndices, const bool& randColors) {
    num_point_clouds_ = numPC;
    point_cloud_data_array_ = new PointCloudData[numPC];
    index_pcA_ = 0;
    index_pcB_ = 1;
    point_size_ = 1.;
    point_cloud_data_array_[index_pcA_] = randColors ? PointCloudData(pointClouds[index_pcA_]) : PointCloudData(pointClouds[index_pcA_], colorReferencePointCloud);
    point_cloud_data_array_[index_pcB_] = randColors ? PointCloudData(pointClouds[index_pcB_]) : PointCloudData(pointClouds[index_pcB_], colorRegisteringPointCloud);
    
    Transform <double , 3, Affine > initial_global_transform = Transform <double , 3, Affine >:: Identity();
    initial_global_transform.translate(initialGlobalOffset);
    for (int i=2; i<numPC; i++){
        point_cloud_data_array_[i] = PointCloudData(pointClouds[i]);
        point_cloud_data_array_[i].TransformVerticesUser(initial_global_transform);
    }
    icp_handler_ = ICPHandler(pointClouds[index_pcA_], pointClouds[index_pcB_], icp_params.subsampling_, icp_params.convergence_threshold_, icpMode);
    point_cloud_registration_correspondent_indices = pcCorrIndices;
    ComputeManualRigidTransformFromParameters();
    TransformFromParameters();
}

void ICPViewerManager::ComputeManualRigidTransformFromParameters() { 
    manual_rigid_transform_matrix_ = Transform <double , 3, Affine >:: Identity();
    manual_rigid_transform_matrix_.rotate( AngleAxisd(manual_rigid_transform_params.rotation_angles_[0], Vector3d::UnitX() ) );
    manual_rigid_transform_matrix_.rotate( AngleAxisd(manual_rigid_transform_params.rotation_angles_[1], Vector3d::UnitY() ) );
    manual_rigid_transform_matrix_.rotate( AngleAxisd(manual_rigid_transform_params.rotation_angles_[2], Vector3d::UnitZ() ) );
    manual_rigid_transform_matrix_.translate( Vector3d(manual_rigid_transform_params.translation_[0],
                                                       manual_rigid_transform_params.translation_[1],
                                                       manual_rigid_transform_params.translation_[2]) );
}

void ICPViewerManager::TransformFromParameters() {
    ComputeManualRigidTransformFromParameters();
    point_cloud_data_array_[index_pcB_].TransformVerticesUser(manual_rigid_transform_matrix_);
    icp_handler_.set_verticesB(point_cloud_data_array_[index_pcB_].updated_vertices());
}

void ICPViewerManager::ChangeColorsRandomly(igl::opengl::glfw::Viewer& viewer) {
    for (int i = 0; i< num_point_clouds_; i++){
        point_cloud_data_array_[i].ChangeColorRandomly();
        viewer.data_list[i].clear();
        viewer.data_list[i].add_points(point_cloud_data_array_[i].updated_vertices(), point_cloud_data_array_[i].colors());
    }
}

void ICPViewerManager::AddNoise(){
    if (noise_params.sigma_rate_ !=0) point_cloud_data_array_[index_pcA_].AddNoise(noise_params.sigma_rate_);
    icp_handler_.set_verticesA(point_cloud_data_array_[index_pcA_].updated_vertices());
}

void ICPViewerManager::DeleteNoise() {
    point_cloud_data_array_[index_pcA_].ResetVertices();
    icp_handler_.set_verticesA(point_cloud_data_array_[index_pcA_].updated_vertices());
}

void ICPViewerManager::SetIcpSubsampling() {
    icp_handler_.set_subsampling(icp_params.subsampling_);
}

void ICPViewerManager::SetConvergenceThreshold(){
    icp_handler_.set_convergence_threshold(icp_params.convergence_threshold_);
}

void ICPViewerManager::SetMaxIterations(){
    icp_handler_.set_max_iterations(icp_params.max_num_iterations_);
}

void ICPViewerManager::SetIcpMode() {
    icp_handler_.set_icp_mode(icp_params.icp_mode_);
}

void ICPViewerManager::ShowNormals(igl::opengl::glfw::Viewer& viewer) {
    if (icp_params.icp_mode_ == POINT_TO_PLANE){
        viewer.data_list[index_pcA_].clear();
        MatrixXd normalColors = MatrixXd::Constant(icp_handler_.normalsA().rows(), 3, 0.5) + 0.5 * icp_handler_.normalsA();
        viewer.data_list[index_pcA_].add_points(point_cloud_data_array_[index_pcA_].updated_vertices(), normalColors);
    }
}

void ICPViewerManager::ConfigureIcpPair(igl::opengl::glfw::Viewer& viewer){
    // Change the color of the current registering point cloud to white for validation
    point_cloud_data_array_[index_pcB_].set_colors(colorReferencePointCloud);
    viewer.data_list[index_pcB_].clear();
    viewer.data_list[index_pcB_].add_points(point_cloud_data_array_[index_pcB_].updated_vertices(), point_cloud_data_array_[index_pcB_].colors());
    // Change the registering point cloud index
    index_pcB_ = index_pcB_ + 1;
    index_pcA_ = point_cloud_registration_correspondent_indices[index_pcB_];
    
    // Move back the points of the new registering point cloud to their initial position, and change their colour to black
    point_cloud_data_array_[index_pcB_].set_colors(colorRegisteringPointCloud);
    point_cloud_data_array_[index_pcB_].ResetVertices();
    icp_handler_.set_verticesA(point_cloud_data_array_[index_pcA_].updated_vertices());
    icp_handler_.set_verticesB(point_cloud_data_array_[index_pcB_].updated_vertices());
}

void ICPViewerManager::PerformIcpStep() { 
    icp_handler_.ComputeIcpStep();
    point_cloud_data_array_[index_pcB_].set_updated_vertices(icp_handler_.verticesB());
}

void ICPViewerManager::PerformIcpRegistration() { 
    icp_handler_.ComputeIcpRegistration();
    point_cloud_data_array_[index_pcB_].set_updated_vertices(icp_handler_.verticesB());
}

void ICPViewerManager::InitDisplay(igl::opengl::glfw::Viewer& viewer) {
    viewer.data().point_size = point_size_;
    viewer.data().add_points(point_cloud_data_array_[0].original_vertices(), point_cloud_data_array_[0].colors());
    for (int i = 1; i < num_point_clouds_; i++){
        viewer.append_mesh();
        viewer.data().point_size = point_size_;
        viewer.data().add_points(point_cloud_data_array_[i].updated_vertices(), point_cloud_data_array_[i].colors());
    }
    if(num_point_clouds_ > 2){
        cout << "GLOBAL REGISTRATION : " << endl << "Align Point cloud 1 with 0." << endl;
    }
}

void ICPViewerManager::UpdateDisplay(igl::opengl::glfw::Viewer& viewer) {
    /*if(icp_params.icp_mode_ == POINT_TO_PLANE){
        viewer.data_list[0].clear();
        viewer.data_list[0].set_vertices(point_cloud_dataA_.updated_vertices());
        viewer.data().set_colors(point_cloud_dataA_.colors());
        viewer.data().set_normals(icp_handler_.normalsA());
    }*/
    //viewer.data_list[index_pcA_].clear();
    //viewer.data_list[index_pcA_].add_points(point_cloud_data_array_[index_pcA_].updated_vertices(), point_cloud_data_array_[index_pcA_].colors());
    viewer.data_list[index_pcB_].clear();
    viewer.data_list[index_pcB_].point_size = point_size_;
    viewer.data_list[index_pcB_].add_points(point_cloud_data_array_[index_pcB_].updated_vertices(), point_cloud_data_array_[index_pcB_].colors());
    viewer.data_list[index_pcA_].clear();
    viewer.data_list[index_pcA_].point_size = point_size_;
    viewer.data_list[index_pcA_].add_points(point_cloud_data_array_[index_pcA_].updated_vertices(), point_cloud_data_array_[index_pcA_].colors());
}



