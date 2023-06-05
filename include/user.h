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
    User(double t, double beta, double gamma, Edge edge, ContractItem contractItem);
    User(double t, double beta, double gamma, Edge edge);
    double t;
    double beta;
    double gamma;
    ContractItem contractItem;
    // 属于哪个 edge
    Edge edge;

    double utility(ContractItem contractItem);

    void chooseContractItem(Contract contract);
};

#endif //CONTRACT_USER_H
