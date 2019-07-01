//
//  PointCloudData.cpp
//  3DA_assignment1_bin
//
//  Created by Benjamin Barral on 11/02/2019.
//

#include "PointCloudData.hpp"
#include <random>

using namespace std;

PointCloudData::PointCloudData() {
    num_points_ = 0;
}

PointCloudData::PointCloudData(const MatrixXd &originalVertices) {
    original_vertices_ = originalVertices;
    updated_vertices_ = originalVertices;
    ChangeColorRandomly();
    num_points_ = originalVertices.rows();
}

PointCloudData::PointCloudData(const MatrixXd &originalVertices, const Vector3d& color) {
    original_vertices_ = originalVertices;
    updated_vertices_ = originalVertices;
    for (int i = 0; i<3; i++){
        colors_(0,i) = color(i);
    }
    num_points_ = originalVertices.rows();
}

MatrixXd PointCloudData::original_vertices() const { 
    return original_vertices_;
}

void PointCloudData::set_original_vertices(const MatrixXd &originalVertices) { 
    original_vertices_ = originalVertices;
}

MatrixXd PointCloudData::updated_vertices() const { 
    return updated_vertices_;
}

void PointCloudData::ResetVertices() {
    updated_vertices_ = original_vertices_;
}

void PointCloudData::ChangeColorRandomly() { 
    for (int i=0; i<3; i++){
        colors_(0,i) = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    }
}

void PointCloudData::TransformVerticesUser(const Transform <double , 3, Affine > &transformMatrix) {
    Matrix<double, 4, Dynamic> verticesBHomogeneous = (original_vertices_.rowwise().homogeneous()).transpose();
    MatrixXd temp = (transformMatrix * verticesBHomogeneous).transpose();
    updated_vertices_ = temp.rowwise().hnormalized();
}

void PointCloudData::UpdateVertices(const MatrixXd &transformMatrix) { 
    Matrix<double, 4, Dynamic> verticesBHomogeneous = (updated_vertices_.rowwise().homogeneous()).transpose();
    MatrixXd temp = (transformMatrix * verticesBHomogeneous).transpose();
    updated_vertices_ = temp.rowwise().hnormalized();;
}

Matrix<double, 1, 3> PointCloudData::colors() const { 
    return colors_;
}

void PointCloudData::set_updated_vertices(const MatrixXd &updatedVertices) { 
    updated_vertices_ = updatedVertices;
}

void PointCloudData::set_colors(const Vector3d &color){
    for (int i=0; i<3; i++){
        colors_(0,i) = color(i);
    }
}

void PointCloudData::ComputeBoundingBox(){
    xMin = __DBL_MAX__;
    yMin = __DBL_MAX__;
    zMin = __DBL_MAX__;
    xMax = __DBL_MIN__;
    yMax = __DBL_MIN__;
    zMax = __DBL_MIN__;
    for (int i = 0; i < num_points_; i++){
        Vector3d vec = original_vertices_.row(i);
        double x = vec(0);
        double y = vec(1);
        double z = vec(2);
        if (x > xMax){
            xMax = x;
        }
        if (x < xMin){
            xMin = x;
        }
        if (y > yMax){
            yMax = y;
        }
        if (y < yMin){
            yMin = y;
        }
        if (z > zMax){
            zMax = z;
        }
        if (z < zMin){
            zMin = z;
        }
    }
}

void PointCloudData::AddNoise(const double &sigmaRate) {
    ComputeBoundingBox();
    std::default_random_engine generator;
    double sigmaX = sigmaRate * (xMax - xMin);
    double sigmaY = sigmaRate * (yMax - yMin);
    double sigmaZ = sigmaRate * (zMax - zMin);
    normal_distribution<double> distX(0, sigmaX);
    normal_distribution<double> distY(0, sigmaY);
    normal_distribution<double> distZ(0, sigmaZ);
    for (int i =0; i<num_points_; i++){
        updated_vertices_(i,0) = original_vertices_(i,0) + distX(generator);
        updated_vertices_(i,1) = original_vertices_(i,1) + distY(generator);
        updated_vertices_(i,2) = original_vertices_(i,2) + distZ(generator);
    }
}









