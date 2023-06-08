//
// Created by Zebra on 2023/6/2.
//
#include "../include/user.h"


User::User(double t, double beta, double gamma, int t_type, int beta_type, int gamma_type, Edge edge): t(t), beta(beta), gamma(gamma), t_type(t_type), beta_type(beta_type), gamma_type(gamma_type), edge(edge) {
    this->MRS = this->beta - this->gamma * this->t;
}

const double User::utility(ContractItem contractItem) const {
    return V * this->t -this->beta * contractItem.x + this->gamma * this->t * contractItem.x + contractItem.r - ((this->t - contractItem.x) * this->edge.price_edge + contractItem.x * this->edge.price_cloud);
}