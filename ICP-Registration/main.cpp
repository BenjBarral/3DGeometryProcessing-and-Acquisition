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
#include <ANN/ANN.h>
#include "ICPViewerManager.hpp"

#include <iterator>
#include <random>

using namespace std;
using namespace Eigen;

enum TaskId{
    TASK_1 = 0,
    TASK_2,
    TASK_5
};

const TaskId task_Id = TASK_1;
// Change this to one of the three modes to perform one of the three courseworks kinds of task
/* For all the other tasks, stick to TASK_1 and change the parameters in the UI :
 - For task 3 : change to 1 on the 'Noise mode' slider parameter
 - For task 4 : change to 1 on the 'Subsampling mode' slider parameter
 - For task 6 : change to 1 on the 'ICP mode' slider parameter in order to perform global point-to-plane 
 */


int main(int argc, char *argv[])
{
    // Init the viewer
    igl::opengl::glfw::Viewer viewer;
    
    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    // READ THE DATA AND PREPARE THE OBJECTS, BASED ON THE TASK SELECTED
    const int num_data_total = 2;
    string filenames[num_data_total] =
    {"../../data_examples/bun000.ply", //0
        "../../data_examples/bun045.ply"}; /*, //1
        "../../data/bun090.ply", //2
        "../../data/bun180.ply", //3
        "../../data/bun270.ply", //4
        "../../data/bun315.ply", //5
        "../../data/chin.ply", //6
        "../../data/ear_back.ply", //7
        "../../data/top2.ply", //8
        "../../data/top3.ply" //9
    };*/
    
    // Order for PC-to-PC registration :
    // 1->0 - 2->1 - 3->2 - 5->0 - 4->5 - 6->5 - 7->5 8->7 9->8
    const int pc_order_indices[num_data_total] =    {0, 1};//, 2, 3, 5, 4, 6, 7, 8, 9};
    int pc_correspondent_indices[num_data_total] = {-1, 0};//, 1, 2, 0, 4, 4, 4, 7, 8 };
    
    int num_point_clouds;
    MatrixXd* pointClouds;
    MatrixXd F;
    if (task_Id == TASK_1){
        // Normal ICP with too point clouds
        num_point_clouds = 2;
        int ind_first_pc = 0;
        int ind_second_pc = 1;
        pointClouds = new MatrixXd[num_point_clouds];
        igl::readPLY(filenames[ind_first_pc], pointClouds[0], F);
        igl::readPLY(filenames[ind_second_pc], pointClouds[1], F);
    }
    else if (task_Id == TASK_2){
        // Same point cloud with rotation
        num_point_clouds = 2;
        int ind_first_pc = 0;
        pointClouds = new MatrixXd[num_point_clouds];
        igl::readPLY(filenames[ind_first_pc], pointClouds[0], F);
        igl::readPLY(filenames[ind_first_pc], pointClouds[1], F);

    }
    else {
        // Global registration
        num_point_clouds = num_data_total;
        pointClouds = new MatrixXd[num_point_clouds];
        // Put the point cloud data in the specific order described above
        for (int i = 0; i<num_point_clouds; i++){
            igl::readPLY(filenames[pc_order_indices[i]], pointClouds[i], F);
        }
    }
    
    // ICP parameters
    ICPMode icpMode = POINT_TO_POINT;
    
    ICPViewerManager icp_viewer_manager(pointClouds, num_point_clouds, icpMode, pc_correspondent_indices,
                                        task_Id == TASK_1 || task_Id == TASK_2);

    // Menu UI parameters
    int noiseModeSlider = 0;
    int samplingModeSlider = 0;
    int icpModeSlider = 0;
    float userRotateValueX,userRotateValueY,userRotateValueZ = 0;
    float userTranslateValueX, userTranslateValueY, userTranslateValueZ = 0;
    bool transformUserInput = false;
    bool updateView = false;
    bool reset = false;
    
    // ADD THE POINTS TO THE VIEWER
    icp_viewer_manager.InitDisplay(viewer);
    
    menu.callback_draw_custom_window = [&]()
    {
        updateView = false;
        transformUserInput = false;
        reset = false;
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(500, 300), ImGuiSetCond_FirstUseEver);
        ImGui::Begin(
                     "MyProperties", nullptr,
                     ImGuiWindowFlags_NoSavedSettings
                     );
        
        // PICK THE POINT CLOUD TO REGISTER
        if(ImGui::Button("Global registration step validation (task 5)")){
            icp_viewer_manager.ConfigureIcpPair(viewer);
            reset = true;
            updateView = true;
        }
        
        // TRANSFORM FROM USER INPUT
        if (ImGui::SliderFloat("Rotation slider (X axis)", &userRotateValueX, -180, 180) ){
            transformUserInput = true;
            icp_viewer_manager.manual_rigid_transform_params.rot_ind_ = 0;
            icp_viewer_manager.manual_rigid_transform_params.UpdateAngle(userRotateValueX);
        }
        if (ImGui::SliderFloat("Rotation slider (Y axis)", &userRotateValueY, -180, 180) ){
        //if (ImGui::InputFloat("Rotation slider (Y axis)", &userRotateValueY)){
            transformUserInput = true;
            icp_viewer_manager.manual_rigid_transform_params.rot_ind_ = 1;
            icp_viewer_manager.manual_rigid_transform_params.UpdateAngle(userRotateValueY);
        }
        if (ImGui::SliderFloat("Rotation slider (Z axis)", &userRotateValueZ, -180, 180) ){
            transformUserInput = true;
            icp_viewer_manager.manual_rigid_transform_params.rot_ind_ = 2;
            icp_viewer_manager.manual_rigid_transform_params.UpdateAngle(userRotateValueZ);
        }
        if (ImGui::SliderFloat("Translation slider (X axis)", &userTranslateValueX, -0.2, 0.2) ){
            transformUserInput = true;
            icp_viewer_manager.manual_rigid_transform_params.trans_ind_ = 0;
            icp_viewer_manager.manual_rigid_transform_params.UpdateTranslation(userTranslateValueX);
        }
        if (ImGui::SliderFloat("Translation slider (Y axis)", &userTranslateValueY, -0.2, 0.2) ){
            transformUserInput = true;
            icp_viewer_manager.manual_rigid_transform_params.trans_ind_ = 1;
            icp_viewer_manager.manual_rigid_transform_params.UpdateTranslation(userTranslateValueY);
        }
        if (ImGui::SliderFloat("Translation slider (Z axis)", &userTranslateValueZ, -0.2, 0.2) ){
            transformUserInput = true;
            icp_viewer_manager.manual_rigid_transform_params.trans_ind_ = 2;
            icp_viewer_manager.manual_rigid_transform_params.UpdateTranslation(userTranslateValueZ);
        }
        if (ImGui::Button("Reset", ImVec2(-1, 0))){
            reset = true;
            transformUserInput = true;
        }
        if (reset){
            icp_viewer_manager.manual_rigid_transform_params.Reset();
            userRotateValueX = 0;
            userRotateValueY = 0;
            userRotateValueZ = 0;
            userTranslateValueX = 0;
            userTranslateValueY = 0;
            userTranslateValueZ = 0;
        }
        if (transformUserInput){
            icp_viewer_manager.TransformFromParameters();
            updateView = true;
        }
        // CHANGE COLOR RANDOMLY
        if (ImGui::Button("Change colors", ImVec2(-1, 0))){
            icp_viewer_manager.ChangeColorsRandomly(viewer);
            updateView = true;
        }
        // CHANGE POINT SIZE
        if (ImGui::SliderFloat("Point size", &(icp_viewer_manager.point_size_), .5, 3.)) updateView = true;
        // SHOW NORMALS IN CASE OF POINT-TO-PLANE
        if(ImGui::Button("Show normals")){
            icp_viewer_manager.ShowNormals(viewer);
        }
        
        // ADD NOISE
        if(ImGui::SliderInt("Noise mode (task 3)", &noiseModeSlider, 0, 1)){
            icp_viewer_manager.noise_params.noise_mode_ = bool(noiseModeSlider);
            if(!icp_viewer_manager.noise_params.noise_mode_){
                icp_viewer_manager.DeleteNoise();
                updateView = true;
            }
        }
        if (ImGui::SliderFloat("Noise sigma rate (task 3)", &(icp_viewer_manager.noise_params.sigma_rate_),0,0.2)){
            if(icp_viewer_manager.noise_params.noise_mode_){
                icp_viewer_manager.AddNoise();
            }
            updateView = true;
        }
        
        // ICP
        // ICP parameters
        if(ImGui::SliderInt("Subsampling mode (task 4)", &samplingModeSlider, 0, 1)){
            icp_viewer_manager.icp_params.subsampling_.subsample_mode_ = bool(samplingModeSlider);
            if (!icp_viewer_manager.icp_params.subsampling_.subsample_mode_) icp_viewer_manager.SetIcpSubsampling();
        }
        if(ImGui::InputDouble("Sumbsampling rate (task 4)", &(icp_viewer_manager.icp_params.subsampling_.subsample_rate_))){
            icp_viewer_manager.SetIcpSubsampling();
        }
        if(ImGui::InputDouble("Error threshold", &(icp_viewer_manager.icp_params.convergence_threshold_))){
            icp_viewer_manager.SetConvergenceThreshold();
        }
        if(ImGui::InputInt("Max num iterations", &(icp_viewer_manager.icp_params.max_num_iterations_))){
            icp_viewer_manager.SetMaxIterations();
        }
        if(ImGui::SliderInt("ICP mode (task 6)", &icpModeSlider, 0, 1)){
            icp_viewer_manager.icp_params.icp_mode_ = ICPMode(icpModeSlider);
            icp_viewer_manager.SetIcpMode();
            //updateView = true;
        }
        
        // Run ICP
        if (ImGui::Button("ICP Step")){
            icp_viewer_manager.PerformIcpStep();
            updateView = true;
        }
        if (ImGui::Button("ICP Registration")){
            icp_viewer_manager.PerformIcpRegistration();
            updateView = true;
        }
        
        // UPDATE THE POINT CLOUDS
        if (updateView){
            icp_viewer_manager.UpdateDisplay(viewer);
        }
 
        ImGui::End();
    };
    
    // Call GUI
    viewer.launch();
    
    annClose();
    // readPLY
    // add_point
    
}


