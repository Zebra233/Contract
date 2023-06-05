#include <iostream>
#include "include/user.h"
#include "include/contract.h"
#include "include/setting.h"
#include "include/edge.h"


using namespace std;


int main() {
    ContractItem contractItem1(0.5, 0.5);
    ContractItem contractItem2(0.6, 0.6);
    ContractItem contractItem3(0.7, 0.7);

    Contract contract;

    contract.contractMap[make_tuple(1, 1, 1)] = contractItem1;
    contract.contractMap[make_tuple(2, 2, 2)] = contractItem2;
    contract.contractMap[make_tuple(3, 3, 3)] = contractItem3;

    Edge edge(0.5, 0.5);

    User user1(1, 1, 1, edge);
    User user2(2, 2, 2, edge);
    User user3(3, 3, 3, edge);



    cout << "user type-(" << user1.t << ", " << user1.beta << ", " << user1.gamma << ")  contractItem" << "1 " <<  user1.utility(contractItem1) << endl;

    return 0;
}
