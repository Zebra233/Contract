//
// Created by Zebra on 2023/6/2.
//
#include "../include/user.h"
#include <iostream>

User::User(double t, double beta, double gamma, int t_type, int beta_type, int gamma_type, Edge edge, double distribution): t(t), beta(beta), gamma(gamma), t_type(t_type), beta_type(beta_type), gamma_type(gamma_type), edge(edge), distribution(distribution) {
    this->MRS = - this->beta + this->gamma * this->t * (edge.price_edge - edge.price_cloud);
}

const double User::utility(ContractItem contractItem) const {
    return this->beta * contractItem.x - this->gamma * this->t * contractItem.x * (this->edge.price_edge - this->edge.price_cloud) + contractItem.r;
}

void User::printUser() {
    // std::cout << "t: " << this->t  << " beta: " << this->beta << " gamma: " << this->beta << " t_type: " << " MRS: " << this->MRS << " distribution: " << this->distribution << std::endl;
    cout << "user type-(" << this->t_type << ", " << this->beta_type << ", " << this->gamma_type << ") " << "MRS " << this->MRS << " 分布 " << this->distribution << endl;
}