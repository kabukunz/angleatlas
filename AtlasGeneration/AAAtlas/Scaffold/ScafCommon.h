#pragma once
/*some code from libigl*/
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <queue>
#include <iostream>
#include <iterator>
using namespace std;
void adjacency_matrix(const Eigen::MatrixXi& F, Eigen::SparseMatrix<int>& A);
void components(const Eigen::SparseMatrix<int>& A, Eigen::MatrixXi& C, Eigen::MatrixXi& counts);
void doublearea(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXd& dblA);
void remove_unreferenced(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, Eigen::MatrixXd& PV, Eigen::MatrixXi& PF, Eigen::VectorXi& VI);

void boundary_loop(const Eigen::MatrixXi &F_ref, std::vector<std::vector<int>>& boundaryloop);

void readobj(const std::string & filename, Eigen::MatrixXd& V, Eigen::MatrixXd& UV_V, Eigen::MatrixXi& F, Eigen::MatrixXi& UV_F);
