#include <iostream>
#include "include/user.h"
#include "include/contract.h"
#include "include/setting.h"
#include "include/edge.h"


using namespace std;

vector<User> users;

// double 表示该类型用户分布
map<User, double> userMap;

// 用户虚拟类型排序
vector<User> userMRS;

// 伪函数用于用户 MRS 排序
class userMRSSort
{
public:
    bool operator() (const User &user1, const User &user2) const
    {
        return user1.MRS < user2.MRS;
    }
};

double getEdgeUtility(Edge edge, ContractItem contractItem) {
    int sum = 0;
    for (auto user: users) {
        sum += (user.t - contractItem.x) * edge.price_edge - contractItem.r - (user.t - contractItem.x);
    }
    return sum;
}

// 按照分布生成用户
// t ~ U(10, 30)
// beta ~ U(1, 2)
// gamma ~ U(0.75, 1.5)

//void userInit(const Edge &edge) {
//    for (int i = 0; i < NumOfUser; i++) {
//        double t = rand() % 20 + 10.0;
//        double beta = (rand() % 10) / 10.0 + 1.0;
//        double gamma = (rand() % 75) / 100.0 + 0.75;
//        User user(t, beta, gamma, edge, 0, 0, 0);
//        users.push_back(user);
//    }
//}

// 直接按照 type 来生成用户
void userInitByType222(const Edge &edge) {
    double t[] = {10, 30};
    double beta[] = {1, 2};
    double gamma[] = {0.75, 1.5};

    // 均匀分布 4种用户类型 12% 另外4种用户类型 13%
    for (int i = 0; i < 4; i++) {
        User user(t[0], beta[i / 2], gamma[i % 2], 0, i / 2, i % 2, edge);
        userMap.insert(make_pair(user, 0.12));
        userMRS.push_back(user);
    }
    for (int i = 4; i < 8; i++) {
        User user(t[1], beta[(i - 4) / 2], gamma[i % 2], 1, (i - 4) / 2, i % 2, edge);
        userMap.insert(make_pair(user, 0.13));
        userMRS.push_back(user);
    }
}


// 按照 MRS 将用户排序
int userSortByMrs() {
    sort(userMRS.begin(), userMRS.end(), userMRSSort());
}

// 求解得到最优契约



int main() {
    ContractItem contractItem1(0.5, 0.5);
    ContractItem contractItem2(0.6, 0.6);
    ContractItem contractItem3(0.7, 0.7);

    Contract contract;

    contract.contractMap.insert(make_pair(make_tuple(1, 1, 1), contractItem1));
    contract.contractMap.insert(make_pair(make_tuple(2, 2, 2), contractItem2));
    contract.contractMap.insert(make_pair(make_tuple(3, 3, 3), contractItem3));

    Edge edge(18, 15);

    userInitByType222(edge);

//    for (auto i : users) {
//        for (auto j : contract.contractMap) {
//            cout << "user type-(" << i.t << ", " << i.beta << ", " << i.gamma << ") 选择契约项 contractItem " << get<0>(j.first) << " " << get<1>(j.first) << " " << get<2>(j.first) << " "
//                    << " x " << j.second.x << " r " << j.second.r << " Utility " <<  i.utility(j.second) << endl;
//        }
//    }

    for (auto i : userMap) {
        for (auto j : contract.contractMap) {
            cout << "user type-(" << i.first.t << ", " << i.first.beta << ", " << i.first.gamma << ") " << "MRS " << i.first.MRS << " 分布 " << i.second << " 选择契约项 contractItem " << get<0>(j.first) << " " << get<1>(j.first) << " " << get<2>(j.first) << " "
                 << " x " << j.second.x << " r " << j.second.r << " Utility " <<  i.first.utility(j.second) << endl;
        }
    }

    userSortByMrs();

    // 输出排序后的用户
    for (auto i : userMRS) {
        cout << "user type-(" << i.t << ", " << i.beta << ", " << i.gamma << ") " << "MRS " << i.MRS << endl;
    }



    return 0;
}


//#include <ilcplex/ilocplex.h>
//
//int main() {
//    // 创建Cplex环境
//    IloEnv env;
//
//    try {
//        // 创建Cplex求解器
//        IloModel model(env);
//        IloCplex cplex(model);
//
//        // 创建决策变量
//        IloNumVar x(env, 0.0, IloInfinity, ILOFLOAT, "x");
//        IloNumVar y(env, 0.0, IloInfinity, ILOFLOAT, "y");
//
//        // 添加目标函数和约束条件
//        model.add(IloMaximize(env, 2*x + 3*y));
//        model.add(3*x + 2*y <= 12);
//        model.add(2*x + y <= 8);
//
//        // 求解优化问题
//        cplex.solve();
//
//        // 打印结果
//        if (cplex.getStatus() == IloAlgorithm::Optimal) {
//            std::cout << "Optimal solution found!" << std::endl;
//            std::cout << "Objective value: " << cplex.getObjValue() << std::endl;
//            std::cout << "x = " << cplex.getValue(x) << std::endl;
//            std::cout << "y = " << cplex.getValue(y) << std::endl;
//        } else {
//            std::cout << "No solution found!" << std::endl;
//        }
//    } catch (IloException& ex) {
//        std::cerr << "Error: " << ex << std::endl;
//    }
//
//    // 释放资源
//    env.end();
//
//    return 0;
//}
