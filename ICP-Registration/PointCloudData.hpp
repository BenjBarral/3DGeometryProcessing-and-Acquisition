//
//  PointCloudData.hpp
//  3DA_assignment1_bin
//
//  Created by Benjamin Barral on 11/02/2019.
//

#ifndef PointCloudData_hpp
#define PointCloudData_hpp

#include <stdio.h>
#include <igl/readPLY.h>
using namespace Eigen;

class PointCloudData{
private:
    MatrixXd original_vertices_,updated_vertices_;
    Matrix<double,1,3> colors_;
    int num_points_;
    double xMin,xMax,yMin,yMax,zMin,zMax; // for the bounding box computation
public:
    PointCloudData();
    PointCloudData(const MatrixXd& originalVertices);
    PointCloudData(const MatrixXd &originalVertices, const Vector3d& color);
    MatrixXd original_vertices() const;
    void set_original_vertices(const MatrixXd& originalVertices);
    void set_updated_vertices(const MatrixXd& updatedVertices);
    void set_colors(const Vector3d& color);
    MatrixXd updated_vertices() const;
    Matrix<double,1,3> colors() const;
    void ResetVertices();
    void ChangeColorRandomly();
    void ComputeBoundingBox();
    void AddNoise(const double& sigmaRate);
    void TransformVerticesUser(const Transform <double , 3, Affine >& transformMatrix);
    void UpdateVertices(const MatrixXd& transformMatrix);
    
};

#endif /* PointCloudData_hpp */
