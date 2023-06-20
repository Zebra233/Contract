//
// Created by Zebra on 2023/6/2.
//

#ifndef CONTRACT_USER_H
#define CONTRACT_USER_H

#include "setting.h"
#include "contract.h"
#include "edge.h"

class User {
public:
    User(double t, double beta, double gamma, int t_type, int beta_type, int gamma_type, Edge edge, double distribution);
    double t;
    double beta;
    double gamma;
    int t_type;
    int beta_type;
    int gamma_type;
    double distribution; // 该用户的分布
    // 属于哪个 edge
    Edge edge;
    double MRS;

    const double utility(ContractItem contractItem) const;

    void chooseContractItem(Contract contract);

    bool operator < (const User &user) const {
        if (this->t_type < user.t_type) {
            return true;
        } else if (this->t_type == user.t_type) {
            if (this->beta_type < user.beta_type) {
                return true;
            } else if (this->beta_type == user.beta_type) {
                if (this->gamma_type < user.gamma_type) {
                    return true;
                }
            }
        }
        return false;
    }

    void printUser();
};

#endif //CONTRACT_USER_H
