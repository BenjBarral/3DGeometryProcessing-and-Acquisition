//
//  ICPHandler.cpp
//  3DA_assignment1_bin
//
//  Created by Benjamin Barral on 11/02/2019.
//

#include "ICPHandler.hpp"
#include <Eigen/QR>



static const double MAX_DIST = 0.0028; // For step 2 : point selection
static const int K_NORMAL_ESTIMATION = 7;

ICPHandler::ICPHandler() {
    num_ICP_points_ = 0;
    convergence_threshold_ = 0;
    error_ = 0;
    max_iterations_ = 120;
    num_iterations_ = 0;
    nA_ = 0;
    nB_ = 0;
    icp_mode_ = POINT_TO_POINT;
}

ICPHandler::ICPHandler(const MatrixXd &verticesA, const MatrixXd &verticesB, const Subsampling &subsampling, const double &convergenceThreshold, const ICPMode& icpMode) {
    verticesA_ = verticesA;
    nA_ = verticesA.rows();
    verticesB_ = verticesB;
    nB_ = verticesB.rows();
    ComputeKDTree();
    subsampling_ = subsampling;
    convergence_threshold_ = convergenceThreshold;
    max_iterations_ = 120;
    SubsampleUniformly(); // Based on the value of subsampling_ : subsample or use the whole point cloud
    error_ = 0;
    num_iterations_ = 0;
    rigid_transform_ = Matrix4d::Zero(4, 4);
    icp_mode_ = icpMode;
    if (icp_mode_ == POINT_TO_PLANE){
        ComputeNormals();
    }
}

void ICPHandler::set_verticesA(const MatrixXd& verticesA) {
    verticesA_ = verticesA;
    nA_ = verticesA.rows();
    ComputeKDTree();
    if (icp_mode_ == POINT_TO_PLANE){
        ComputeNormals();
    }
}

void ICPHandler::set_verticesB(const MatrixXd& verticesB) {
    verticesB_ = verticesB;
    nB_ = verticesB.rows();
    SubsampleUniformly();
}

void ICPHandler::set_subsampling(const Subsampling& subsampling) {
    subsampling_ = subsampling;
    SubsampleUniformly();
}

void ICPHandler::set_convergence_threshold(const double& convergenceThreshold) {
    convergence_threshold_ = convergenceThreshold;
}

void ICPHandler::set_max_iterations(const int &maxIterations){
    max_iterations_ = maxIterations;
}

void ICPHandler::set_icp_mode(const ICPMode &icpMode){
    icp_mode_ = icpMode;
    if (icp_mode_ == POINT_TO_PLANE){
        ComputeNormals();
    }
}

MatrixXd ICPHandler::verticesB() const {
    return verticesB_;
}

MatrixXd ICPHandler::normalsA() const{
    return normalsA_;
}

vector<double> ICPHandler::errors_in_steps() const {
    return errors_in_steps_;
}

void ICPHandler::SubsampleUniformly() {
    icp_points_indices.clear();
    int d = (subsampling_.subsample_mode_ && subsampling_.subsample_rate_ !=0) ? int(1.0/subsampling_.subsample_rate_) : 1;
    num_ICP_points_ = 0;
    for (int i = 0; i<nB_; i+=d){
        num_ICP_points_ +=1;
        icp_points_indices.push_back(i);
    }
}

void ICPHandler::ComputeKDTree() {
    int dim = 3;
    ANNpointArray annPointArray = new ANNpoint[nA_]; // annAllocPts(nA_, dim);
    for (int i = 0; i<nA_; i++){
        annPointArray[i] = new ANNcoord[dim]; // annAllocPt(dim);
        for (int j = 0; j<3; j++){
            annPointArray[i][j] = verticesA_(i,j);
        }
    }
    kd_tree_ = new ANNkd_tree(annPointArray,nA_,dim);
}

void ICPHandler::ComputeCorrespondences() { 
    icp_pointsA_.clear();
    icp_pointsB_.clear();
    icp_normalsA_.clear();
    num_ICP_points_ = 0;
    const int k = 1;
    ANNidxArray nnIdx = new ANNidx[k];
    ANNdistArray dists = new ANNdist[k];
    ANNpoint queryPt = new ANNcoord[3];
    for (int i = 0; i<icp_points_indices.size(); i++){
        int ind = icp_points_indices.at(i);
        for (int j = 0; j<3; j++){
            queryPt[j] = verticesB_(ind,j);
        }
        kd_tree_->annkSearch(queryPt, k, nnIdx, dists, 0.0);
        int corrInd = nnIdx[0];
        if (dists[0] <= MAX_DIST * MAX_DIST){
            icp_pointsB_.push_back(verticesB_.row(ind));
            icp_pointsA_.push_back(verticesA_.row(corrInd));
            if (icp_mode_ == POINT_TO_PLANE){
                icp_normalsA_.push_back(normalsA_.row(corrInd));
            }
            num_ICP_points_ +=1;
        }
    }
    delete[] nnIdx;
    delete[] dists;
    delete[] queryPt;
}

void ICPHandler::ComputeBarycentersAndRecenterPoints() {
    barycenterA_ = Vector3d::Zero(3);
    barycenterB_ = Vector3d::Zero(3);
    recentered_pointsA_.clear();
    recentered_pointsB_.clear();
    for (int i = 0; i<num_ICP_points_; i++){
        barycenterA_ += icp_pointsA_.at(i);
        barycenterB_ += icp_pointsB_.at(i);
    }
    barycenterA_ /= num_ICP_points_;
    barycenterB_ /= num_ICP_points_;
    for (int i = 0; i<num_ICP_points_; i++){
        recentered_pointsA_.push_back(icp_pointsA_.at(i) - barycenterA_);
        recentered_pointsB_.push_back(icp_pointsB_.at(i) - barycenterB_);
    }
}

void ICPHandler::ComputeLinearApproximationMatrices() {
    matrix_A_ = MatrixXd::Zero(num_ICP_points_, 6);
    vector_b_ = VectorXd::Zero(num_ICP_points_);
    Vector3d normal, vecA, vecB, crossProd;
    for (int i = 0; i < num_ICP_points_; i++){
        normal = icp_normalsA_.at(i);
        vecA = icp_pointsA_.at(i);
        vecB = icp_pointsB_.at(i);
        vector_b_(i) = normal.dot(vecA - vecB);
        crossProd = vecB.cross(normal);
        for (int j = 0; j < 3; j++){
            matrix_A_(i,j) = crossProd(j);
            matrix_A_(i,j+3) = normal(j);
        }
    }
}

void ICPHandler::ComputeRigidTransformation() {
    if (icp_mode_ == POINT_TO_POINT){
        Matrix3d A = Matrix3d::Zero(3,3);
        for (int i = 0; i<num_ICP_points_; i++){
            A += recentered_pointsB_.at(i) * (recentered_pointsA_.at(i)).transpose();
        }
        JacobiSVD<Matrix3d> svd(A, ComputeFullU | ComputeFullV);
        Matrix3d U = svd.matrixU();
        Matrix3d V = svd.matrixV();
        Matrix3d R = V * U.transpose();
        Vector3d t = barycenterA_ - R * barycenterB_;
        rigid_transform_.block(0, 0, 3, 3) = R;
        rigid_transform_.col(3) = t.homogeneous();
    }
    else{
        VectorXd sol = matrix_A_.colPivHouseholderQr().solve(vector_b_);
        Transform <double , 3, Affine > tempTransform = Transform <double , 3, Affine >:: Identity();
        tempTransform.rotate( AngleAxisd(sol[0], Vector3d::UnitX() ) );
        tempTransform.rotate( AngleAxisd(sol[1], Vector3d::UnitY() ) );
        tempTransform.rotate( AngleAxisd(sol[2], Vector3d::UnitZ() ) );
        tempTransform.translate( Vector3d(sol[3], sol[4], sol[5]));
        rigid_transform_ = tempTransform.matrix();
    }
}

void ICPHandler::TransformPoints() {
    Matrix<double, 4, Dynamic> verticesBHomogeneous = (verticesB_.rowwise().homogeneous()).transpose();
    MatrixXd temp = (rigid_transform_ * verticesBHomogeneous).transpose();
    verticesB_ = temp.rowwise().hnormalized();
}

void ICPHandler::ComputeError() {
    error_ = 0;
    Vector3d errorVec, vecA, vecB, normal;
    Vector4d vecBHom;
    for (int i = 0; i<num_ICP_points_; i++){
        vecA = icp_pointsA_.at(i);
        vecBHom = (icp_pointsB_.at(i)).homogeneous() ;
        vecB = (rigid_transform_ * vecBHom).hnormalized();
        errorVec = vecA - vecB;
        if (icp_mode_ == POINT_TO_POINT) error_ += errorVec.squaredNorm();
        else {
            normal = icp_normalsA_.at(i);
            error_ += pow(errorVec.dot(normal) , 2);
        }
    }
    error_ = sqrt(error_/num_ICP_points_);
}

void ICPHandler::PrintStepError() const {
    cout << "Error at step " << num_iterations_ << " : " << error_ << " (square root of mean squared distances)" << endl;
}

void ICPHandler::ComputeIcpStep() {
    ComputeCorrespondences();
    if (icp_mode_ == POINT_TO_POINT) ComputeBarycentersAndRecenterPoints();
    else ComputeLinearApproximationMatrices();
    ComputeRigidTransformation();
    TransformPoints();
    ComputeError();
}


void ICPHandler::Reset() {
    error_ = MAXFLOAT;
    num_iterations_ = 0;
    errors_in_steps().clear();
}

void ICPHandler::PrintConvergenceState() const {
    if (num_iterations_ == max_iterations_) cout << "No convergence under " << max_iterations_ << " iterations." << " Final error : " << error_ << endl;
    else cout << "Convergence of ICP after " << num_iterations_ << " iterations." << " Final error : " << error_ <<
        "Number of ICP points : " << num_ICP_points_ << endl;
}

void ICPHandler::ComputeIcpRegistration() {
    Reset();
    double convThreshold = (icp_mode_ == POINT_TO_POINT)?convergence_threshold_ : convergence_threshold_ * 0.5;
    while(error_ > convThreshold && num_iterations_ < max_iterations_){
        ComputeIcpStep();
        num_iterations_ +=1;
        errors_in_steps().push_back(error_);
        PrintStepError();
    }
    PrintConvergenceState();
}

void ICPHandler::ComputeNormals(){
    normalsA_ = MatrixXd::Zero(nA_,3);
    ANNidxArray nnIdx = new ANNidx[K_NORMAL_ESTIMATION];
    ANNdistArray dists = new ANNdist[K_NORMAL_ESTIMATION];
    ANNpoint queryPt = new ANNcoord[3];
    Vector3d point,neighborPt;
    Matrix3d covariance;
    Vector3d normal;
    for (int i = 0; i<nA_;i++){
        point = verticesA_.row(i);
        for (int j = 0; j<3; j++){
            queryPt[j] = point(j);
        }
        kd_tree_->annkSearch(queryPt, K_NORMAL_ESTIMATION, nnIdx, dists, 0.0);
        covariance = Matrix3d::Zero(3,3);
        for (int q = 1; q < K_NORMAL_ESTIMATION; q++){
            neighborPt = verticesA_.row(nnIdx[q]);
            covariance += (neighborPt - point)*((neighborPt - point).transpose());
        }
        covariance /= (K_NORMAL_ESTIMATION - 1);
        EigenSolver<MatrixXd> es(covariance);
        int minEigenValueIndex = - 1;
        double minNorm = DBL_MAX;
        double maxNorm = DBL_MIN;
        for (int r = 0; r < 3; r++){
            double a = (es.eigenvalues()[r]).real();
            double b = (es.eigenvalues()[r]).imag();
            double eigenNorm = a * a + b * b;
            if (eigenNorm < minNorm){
                minEigenValueIndex = r;
                minNorm = eigenNorm;
            }
            /*if (eigenNorm > maxNorm){
                minEigenValueIndex = r;
                maxNorm = eigenNorm;
            }*/
        }
        normal = es.eigenvectors().col(minEigenValueIndex).real();
        normalsA_.row(i) = normal;
    }
    
}





