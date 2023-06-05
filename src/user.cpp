//
// Created by Zebra on 2023/6/2.
//
#include "../include/user.h"


User::User(double t, double beta, double gamma, Edge edge, ContractItem contractItem): t(t), beta(beta), gamma(gamma), edge(edge), contractItem(contractItem) {

}
User::User(double t, double beta, double gamma, Edge edge): t(t), beta(beta), gamma(gamma), edge(edge) {

};

double User::utility(ContractItem contractItem) {
    return -this->beta * contractItem.x - this->gamma * this->t * contractItem.x + contractItem.r - ((this->t - contractItem.x) * this->edge.price_edge + contractItem.x * this->edge.price_cloud);
}