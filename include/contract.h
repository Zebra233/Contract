//
// Created by Zebra on 2023/6/2.
//

#ifndef CONTRACT_CONTRACT_H
#define CONTRACT_CONTRACT_H

#include <vector>
#include <map>
#include <tuple>
#include <string>
#include "setting.h"

using namespace std;

class ContractItem {
public:
    ContractItem(double x, double r);

    double x; // 卸载到云到比例
    double r; // 补贴

};

class Contract {
public:
    map<tuple<int, int, int>, ContractItem> contractMap;
};

#endif //CONTRACT_CONTRACT_H
