//
// Created by Zebra on 2023/6/2.
//

#ifndef CONTRACT_EDGE_H
#define CONTRACT_EDGE_H

#include "setting.h"
#include "contract.h"

class Edge {
public:
    Edge(double price_edge, double price_cloud);

    double price_edge;
    double price_cloud;
    Contract contract;

};

#endif //CONTRACT_EDGE_H
