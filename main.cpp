#include <iostream>
#include "include/user.h"
#include "include/contract.h"
#include "include/setting.h"
#include "include/edge.h"
#include <ilcplex/ilocplex.h>

using namespace std;

vector<User> users;

// double 表示该类型用户分布
map<User, double> userMap;

// 用户虚拟类型排序
vector<User> userMRS;

// 契约项
vector<ContractItem> contractItems;

// 用户 type-i 选择 contract i 的效用
// 用于验证 IC 约束
vector<double> UserIConIUtil;

// 伪函数用于用户 MRS 排序
class userMRSSort
{
public:
    bool operator() (const User &user1, const User &user2) const
    {
        return user1.MRS < user2.MRS;
    }
};

// 得到边缘效用
double getEdgeUtility(Edge edge) {
    int sum = 0;
    for (int i = 0; i < userMRS.size(); i++) {
        sum += NumOfUser * userMRS[i].distribution * (userMRS[i].t * contractItems[i].x * edge.price_edge - contractItems[i].r - (userMRS[i].t * contractItems[i].x * userMRS[i].t * contractItems[i].x));
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
    // 能得到不符合 IC 约束的解
//    double t[] = {8, 12};
//    double beta[] = {3, 8};
//    double gamma[] = {0.6, 0.8};
// 不符合 IC 3个
//    double t[] = {8, 16};
//    double beta[] = {3, 8};
//    double gamma[] = {0.6, 0.8};
// 不符合 IC 3个
    double t[] = {8, 16};
    double beta[] = {3, 8};
    double gamma[] = {0.6, 0.88};


    // 均匀分布 4种用户类型 12% 另外4种用户类型 13%
    for (int i = 0; i < 4; i++) {
        User user(t[0], beta[i / 2], gamma[i % 2], 0, i / 2, i % 2, edge, 0.12);
        userMap.insert(make_pair(user, 0.12));
        userMRS.push_back(user);
    }
    for (int i = 4; i < 8; i++) {
        User user(t[1], beta[(i - 4) / 2], gamma[i % 2], 1, (i - 4) / 2, i % 2, edge, 0.13);
        userMap.insert(make_pair(user, 0.13));
        userMRS.push_back(user);
    }
}

// 按照 MRS 将用户排序
void userSortByMrs() {
    sort(userMRS.begin(), userMRS.end(), userMRSSort());
}

// 契约项初始化
void contractItemInit() {
    for (int i = 0; i < K * S * M; i++) {
        ContractItem contractItem(0.0, 0.0);
        contractItems.push_back(contractItem);
    }
}

// 计算目标函数中的分布和
double Q1ToN(int n) {
    double sum = 0.0;
    for (int i = 0; i <= n; i++) {
        sum += userMRS[i].distribution;
    }
    return sum;
}

// Y(x, Gamma)
double Y(double x, int i) {
    return userMRS[i].beta * x - userMRS[i].gamma * userMRS[i].t * x * (userMRS[i].edge.price_edge - userMRS[i].edge.price_cloud);
}

// 通过 x 得到 r
double getRFromX(double x, int i) {
    if (i == K * S * M - 1) {
        return - Y(x ,i);
    } else {
        return contractItems[i + 1].r + Y(contractItems[i + 1].x, i) - Y(x, i);
    }
//    return sqrt(abs(userMRS[i].MRS)) * x;
}

// old 求解优化问题得到 x_i
double getXBySolveopt(int i) {
    // 创建Cplex环境
    IloEnv env;

    try {
        // 创建Cplex求解器
        IloModel model(env);
        IloCplex cplex(model);

        // 定义变量 x_i
        IloNumVar x(env, 0, userMRS[i].t, ILOFLOAT, "x");

        // 松弛了单调性条件的约束

        // 设置目标函数
        IloExpr objectiveExpr(env);
        IloExpr income(env);
        IloExpr cost(env);
        IloExpr Yi(env);

        income = (userMRS[i].t - x) * userMRS[i].edge.price_edge;
        //cost =  coefficient * IloLog(userMRS[0].t - vars[0]);
        cost =  1.0 / (1.0 - (double)CostItemCoefficient) * IloPower(userMRS[i].t - x, 1.0 - (double)CostItemCoefficient);
        Yi =  (V * userMRS[i].t - userMRS[i].beta * x + userMRS[i].gamma * userMRS[i].t * x -
               ((userMRS[i].t - x) * userMRS[i].edge.price_edge +
                x * userMRS[i].edge.price_cloud)) * Q1ToN(i);

        if (i != 0) {
            IloExpr Yi_1(env);
            Yi_1 = (V * userMRS[i - 1].t - userMRS[i - 1].beta * x + userMRS[i - 1].gamma * userMRS[i - 1].t * x -
                    ((userMRS[i - 1].t - x) * userMRS[i - 1].edge.price_edge +
                    x * userMRS[i - 1].edge.price_cloud)) * Q1ToN(i - 1);
            objectiveExpr =
                    userMRS[i].distribution * (income  - cost) -
                    Yi_1 +
                    Yi;
        } else {
            objectiveExpr =
                    userMRS[i].distribution * (income  - cost) +
                    Yi;
        }

        IloObjective objective = IloMaximize(env, objectiveExpr);
        model.add(objective);

        // 求解优化问题
        cplex.solve();

        // 打印结果
        if (cplex.getStatus() == IloAlgorithm::Optimal) {
            // 获取结果
            IloNumArray solution(env);
            double objectiveValue = cplex.getObjValue();
            IloNum xResult = cplex.getValue(x);

            cout << "---------- opt result ----------" << endl;
            std::cout << "x: " << xResult << std::endl;
            std::cout << "Objective Value: " << objectiveValue << std::endl;

            contractItems[i].x = (double)xResult;
            contractItems[i].r = getRFromX(xResult, i);
            cout << "--------------------------------" << endl;

            cout << "契约项: " << contractItems[i].x << " " << contractItems[i].r << endl;
            std::cout << "用户效用: " << userMRS[i].utility(contractItems[i]) << std::endl;
            std::cout << "边缘效用: " << (userMRS[i].t - xResult) * userMRS[i].edge.price_edge - contractItems[i].r -  1 / (1 - CostItemCoefficient) * pow(userMRS[i].t - xResult, 1 - CostItemCoefficient) << std::endl;
            std::cout << "边缘效用2: " << (userMRS[i].t - 5) * userMRS[i].edge.price_edge - contractItems[i].r -  1 / (1 - CostItemCoefficient) * pow(userMRS[i].t - 5, 1 - CostItemCoefficient) << std::endl;
        } else {
            std::cout << "No solution found!" << std::endl;
            exit(-1);
        }
    } catch (IloException& ex) {
        std::cerr << "Error: " << ex << std::endl;
    }

    // 释放资源
    env.end();
}

// 修改模型，改为卸载到边缘的比例，求解优化问题得到 x_i
// start 表示这一段递增自序列的开始是哪，
double getXBySolveopt2(int start, int i, double lb, double ub) {
    // 创建Cplex环境
    IloEnv env;

    try {
        // 创建Cplex求解器
        IloModel model(env);
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());  // 将日志输出重定向到空流，禁止显示求解信息

        // 定义变量 x_i
        IloNumVar x(env, lb, ub, ILOFLOAT, "x");

        // 松弛了单调性条件的约束

        // 设置目标函数
        IloExpr objectiveExpr(env);
        IloExpr income(env);
        IloExpr cost(env);
        IloExpr Yi(env);
        IloExpr Yi_1(env);

        income = userMRS[i].t * x * userMRS[i].edge.price_edge;
        cost =  userMRS[i].t * x * userMRS[i].t * x;
        Yi =  (userMRS[i].beta * x -
                userMRS[i].gamma * userMRS[i].t * x * (userMRS[i].edge.price_edge - userMRS[i].edge.price_cloud))
                * Q1ToN(i);

        if (i != start) {
            Yi_1 = (userMRS[i - 1].beta * x -
                    userMRS[i - 1].gamma * userMRS[i - 1].t * x * (userMRS[i - 1].edge.price_edge - userMRS[i - 1].edge.price_cloud))
                            * Q1ToN(i - 1);
            objectiveExpr =
                    userMRS[i].distribution * (income  - cost) -
                    Yi_1 +
                    Yi;
        } else {
            objectiveExpr =
                    userMRS[i].distribution * (income  - cost) +
                    Yi;
        }

        IloObjective objective = IloMaximize(env, objectiveExpr);
        model.add(objective);

        // 求解优化问题
        cplex.solve();

        // 打印结果
        if (cplex.getStatus() == IloAlgorithm::Optimal) {
            // 获取结果
            IloNumArray solution(env);
            double objectiveValue = cplex.getObjValue();
            IloNum xResult = cplex.getValue(x);

//            cout << "---------- opt result ----------" << endl;
//            std::cout << "x: " << xResult << std::endl;
//            std::cout << "Objective Value: " << objectiveValue << std::endl;

            contractItems[i].x = (double)xResult;
            //contractItems[i].r = getRFromX(xResult, i);
            //cout << "--------------------------------" << endl;

            //cout << "契约项: " << contractItems[i].x << " " << contractItems[i].r << endl;
            //std::cout << "用户效用: " << userMRS[i].utility(contractItems[i]) << std::endl;
            //std::cout << "边缘效用: " << (userMRS[i].t - xResult) * userMRS[i].edge.price_edge - contractItems[i].r -  1 / (1 - CostItemCoefficient) * pow(userMRS[i].t - xResult, 1 - CostItemCoefficient) << std::endl;
            //std::cout << "边缘效用2: " << (userMRS[i].t - 5) * userMRS[i].edge.price_edge - contractItems[i].r -  1 / (1 - CostItemCoefficient) * pow(userMRS[i].t - 5, 1 - CostItemCoefficient) << std::endl;
        } else {
            std::cout << "No solution found!" << std::endl;
            exit(-1);
        }
    } catch (IloException& ex) {
        std::cerr << "Error: " << ex << std::endl;
    }

    // 释放资源
    env.end();
}

// 使用最原始的问题模型解，不分离单个的 x
void solveAll() {
    Edge edge(10, 5);

    userInitByType222(edge);

    userSortByMrs();

    contractItemInit();

    for (auto i : userMRS) {
        i.printUser();
    }

    // 创建Cplex环境
    IloEnv env;

    try {
        // 创建Cplex求解器
        IloModel model(env);
        IloCplex cplex(model);
#pragma region 创建变量
        int numVariables = 8;  // 变量数量

        IloNumVarArray variables(env, numVariables);  // 创建变量数组

        for (int i = 0; i < numVariables; i++) {
            variables[i] = IloNumVar(env, 0, userMRS[i].t);  // 设置每个变量的取值范围
        }

        // 将变量数组添加到模型中
        model.add(variables);
#pragma endregion


        // 单调性约束
        // 添加单调递增约束
        for (int i = 0; i < numVariables - 1; i++) {
            model.add(variables[i] >= variables[i + 1]);
        }
        // IR 约束
        model.add(V * userMRS[K*S*M - 1].t - userMRS[K*S*M - 1].beta * variables[numVariables - 1] + userMRS[K*S*M - 1].gamma * userMRS[K*S*M - 1].t * variables[numVariables - 1] +
                          sqrt(abs(userMRS[K*S*M - 1].MRS)) * variables[numVariables - 1] -
                          ((userMRS[K*S*M - 1].t - variables[numVariables - 1]) * userMRS[K*S*M - 1].edge.price_edge + variables[numVariables - 1] * userMRS[K*S*M - 1].edge.price_cloud)
                            >= 0 );

        // IC 约束
        for (int i = 0; i < K*S*M - 1; i++) {
            model.add(sqrt(abs(userMRS[i + 1].MRS)) * variables[i + 1] -
            // Y(x_i, Gamma_{i+1})
            (V * userMRS[i + 1].t - userMRS[i + 1].beta * variables[i] +
            userMRS[i + 1].gamma * userMRS[i + 1].t * variables[i] -
            ((userMRS[i + 1].t - variables[i]) * userMRS[i + 1].edge.price_edge +
            variables[i] * userMRS[i + 1].edge.price_cloud)) +
            // Y(x_{i+1}, Gamma_{i+1})
            (V * userMRS[i + 1].t - userMRS[i + 1].beta * variables[i + 1] +
            userMRS[i + 1].gamma * userMRS[i + 1].t * variables[i + 1] -
            ((userMRS[i + 1].t - variables[i + 1]) * userMRS[i + 1].edge.price_edge +
            variables[i + 1] * userMRS[i + 1].edge.price_cloud))

            >= sqrt(abs(userMRS[i].MRS)) * variables[i]);

            model.add(sqrt(abs(userMRS[i].MRS)) * variables[i] >=
                    sqrt(abs(userMRS[i + 1].MRS)) * variables[i + 1] +
                    // Y(x_{i+1}, Gamma_i)
                    (V * userMRS[i].t - userMRS[i].beta * variables[i + 1]) +
                    userMRS[i].gamma * userMRS[i].t * variables[i + 1] -
                    ((userMRS[i].t - variables[i + 1]) * userMRS[i].edge.price_edge +
                    variables[i + 1] * userMRS[i].edge.price_cloud) -
                    // Y(x_i, Gamma_i)
                    (V * userMRS[i].t - userMRS[i].beta * variables[i]) +
                    userMRS[i].gamma * userMRS[i].t * variables[i] -
                    ((userMRS[i].t - variables[i]) * userMRS[i].edge.price_edge +
                    variables[i] * userMRS[i].edge.price_cloud)
            );
        }


#pragma region 设置目标函数
        // 设置目标函数
        IloExpr objectiveExpr(env);
        IloExpr income(env);
        IloExpr cost(env);
        IloExpr Yi(env);
        IloExpr Yi_1(env);
//        for (int i = 0; i < numVariables; i++) {
//            income = (userMRS[i].t - variables[i]) * userMRS[i].edge.price_edge;
//            double coefficient = 10;  // 对数项的系数
////        cost =  coefficient * IloLog(userMRS[0].t - x);
//            cost =  coefficient * (userMRS[i].t - variables[i]);
//            Yi =  (V * userMRS[i].t - userMRS[i].beta * variables[i] + userMRS[i].gamma * userMRS[i].t * variables[i] -
//                   ((userMRS[i].t - variables[i]) * userMRS[i].edge.price_edge +
//                           variables[i] * userMRS[i].edge.price_cloud)) * Q1ToN(i);
//            if (i != 0) {
//                Yi_1 = (V * userMRS[i - 1].t - userMRS[i - 1].beta * variables[i] +
//                        userMRS[i - 1].gamma * userMRS[i - 1].t * variables[i] -
//                        ((userMRS[i - 1].t - variables[i]) * userMRS[i - 1].edge.price_edge +
//                         variables[i] * userMRS[i - 1].edge.price_cloud)) * Q1ToN(i - 1);
//
//                objectiveExpr += userMRS[i].distribution * (income - cost) - Yi_1 + Yi;
//            } else {
//                objectiveExpr += userMRS[i].distribution * (income - cost) + Yi;
//            }
//        }
//        objectiveExpr *= K*S*M;
        for (int i = 0; i < K*S*M; i++) {
            objectiveExpr += K*S*M * userMRS[i].distribution * (
                                                                       (userMRS[i].t - variables[i]) * userMRS[i].edge.price_edge -
                                                                               sqrt(abs(userMRS[i].MRS)) * variables[i]
                                                                               - 10 * (userMRS[i].t - variables[i])
                    );
        }

        IloObjective objective = IloMaximize(env, objectiveExpr);
#pragma endregion
        model.add(objective);

        // 求解优化问题
        cplex.solve();

#pragma region 打印结果
        // 打印结果
        if (cplex.getStatus() == IloAlgorithm::Optimal) {
            // 获取结果
            IloNumArray solution(env);
            double objectiveValue = cplex.getObjValue();
            for (int i = 0; i < 8; i++) {
                IloNum xResult = cplex.getValue(variables[i]);

                cout << "---------- opt result ----------" << endl;
                std::cout << "x: " << xResult << std::endl;
                std::cout << "Objective Value: " << objectiveValue << std::endl;

                contractItems[i].x = (double)xResult;
                contractItems[i].r = getRFromX(xResult, i);
                cout << "--------------------------------" << endl;

                cout << "契约项: " << contractItems[i].x << " " << contractItems[i].r << endl;
                std::cout << "用户效用: " << userMRS[i].utility(contractItems[i]) << std::endl;
                std::cout << "边缘效用: " << (userMRS[i].t - xResult) * userMRS[i].edge.price_edge - contractItems[i].r -  1 / (1 - CostItemCoefficient) * pow(userMRS[i].t - xResult, 1 - CostItemCoefficient) << std::endl;
                //std::cout << "边缘效用2: " << (userMRS[i].t - 5) * userMRS[i].edge.price_edge - contractItems[i].r -  1 / (1 - CostItemCoefficient) * pow(userMRS[i].t - 5, 1 - CostItemCoefficient) << std::endl;
            }

        } else {
            std::cout << "No solution found!" << std::endl;
        }
    } catch (IloException& ex) {
        std::cerr << "Error: " << ex << std::endl;
    }
#pragma endregion

    // 释放资源
    env.end();
}

double solveAll2() {
    // 创建Cplex环境
    IloEnv env;

    try {
        // 创建Cplex求解器
        IloModel model(env);
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());  // 将日志输出重定向到空流，禁止显示求解信息

#pragma region 创建变量
        int numVariables = K * S * M;  // 变量数量

        IloNumVarArray variables(env, numVariables);  // 创建变量数组

        for (int i = 0; i < numVariables; i++) {
            variables[i] = IloNumVar(env, 0, 1);  // 设置每个变量的取值范围
        }

        // 将变量数组添加到模型中
        model.add(variables);
#pragma endregion

#pragma region 约束条件
        // 单调性约束
        // 添加单调递增约束
        for (int i = 0; i < numVariables - 1; i++) {
            model.add(variables[i] >= variables[i + 1]);
        }


#pragma endregion

        // 设置目标函数
        IloExpr objectiveExpr(env);
        IloExpr income(env);
        IloExpr cost(env);
        IloExpr Yi(env);
        IloExpr Yi_1(env);

        for (int i = 0; i < K*S*M; i++) {
            Yi =  userMRS[i].beta * variables[i] -
                    userMRS[i].gamma * userMRS[i].t * variables[i] * (userMRS[i].edge.price_edge - userMRS[i].edge.price_cloud);
            if (i > 0) {
                Yi_1 = userMRS[i - 1].beta * variables[i] -
                   userMRS[i - 1].gamma * userMRS[i - 1].t * variables[i] * (userMRS[i - 1].edge.price_edge - userMRS[i - 1].edge.price_cloud);
                objectiveExpr += K*S*M * (
                    // 收入 - 支出
                    userMRS[i].distribution * (userMRS[i].t * variables[i] * userMRS[i].edge.price_edge -
                    - userMRS[i].t * variables[i] * userMRS[i].t * variables[i]) -
                    // -Y(x_i, Gamma_{i - 1}) Q
                    Yi_1 * Q1ToN(i - 1) +
                    // + Yi Q
                    Yi * Q1ToN(i)
                    );
            } else {
                objectiveExpr += K*S*M * (
                        // 收入 - 支出
                        userMRS[i].distribution * (userMRS[i].t * variables[i] * userMRS[i].edge.price_edge -
                                                   - userMRS[i].t * variables[i] * userMRS[i].t * variables[i]) +
                        // + Yi Q
                        Yi * Q1ToN(i)
                );
            }
        }

        IloObjective objective = IloMaximize(env, objectiveExpr);
        model.add(objective);

        // 求解优化问题
        cplex.solve();

        // 打印结果
        if (cplex.getStatus() == IloAlgorithm::Optimal) {
            // 获取结果
            IloNumArray solution(env);
            double objectiveValue = cplex.getObjValue();


//            cout << "---------- opt result ----------" << endl;
//            std::cout << "x: " << xResult << std::endl;
//            std::cout << "Objective Value: " << objectiveValue << std::endl;
            for (int i = 0; i < K * S * M; i++) {
                IloNum xResult = cplex.getValue(variables[i]);
                contractItems[i].x = (double) xResult;
                contractItems[i].r = getRFromX(xResult, i);
            }
            //cout << "--------------------------------" << endl;

            //cout << "契约项: " << contractItems[i].x << " " << contractItems[i].r << endl;
            //std::cout << "用户效用: " << userMRS[i].utility(contractItems[i]) << std::endl;
            //std::cout << "边缘效用: " << (userMRS[i].t - xResult) * userMRS[i].edge.price_edge - contractItems[i].r -  1 / (1 - CostItemCoefficient) * pow(userMRS[i].t - xResult, 1 - CostItemCoefficient) << std::endl;
            //std::cout << "边缘效用2: " << (userMRS[i].t - 5) * userMRS[i].edge.price_edge - contractItems[i].r -  1 / (1 - CostItemCoefficient) * pow(userMRS[i].t - 5, 1 - CostItemCoefficient) << std::endl;
        } else {
            std::cout << "No solution found!" << std::endl;
            exit(-1);
        }
    } catch (IloException& ex) {
        std::cerr << "Error: " << ex << std::endl;
    }

    // 释放资源
    env.end();
}

// 检测契约项是否是单调的，返回不单调的部分的开头
// 然后从开头直接重新计算
// 如果返回 K*S*M，说明没有问题，子序列全部符合单调性条件
int monTest() {
    // 严格递减就++
    vector<pair<int, int>> temp;
    int i = 1;
    while (i < K*S*M && contractItems[i - 1].x > contractItems[i].x) {
        i++;
    }
    // 非递减
    return i;
}




int main() {
    // test();
//    Edge edge(8, 6);
    Edge edge(8, 6);

    userInitByType222(edge);
    userSortByMrs();
    contractItemInit();

    int start = 0;

    getXBySolveopt2(start, K * S * M - 1, 0, 1);
    for (int i = K * S * M - 2; i >= 0; i--) {
//        getXBySolveopt2(start, i, contractItems[i + 1].x, 1);
        getXBySolveopt2(start, i, 0, 1);
    }

    //solveAll2();

    int monIndex = 0;
    monIndex = monTest();

    while (monIndex != K*S*M) {
        for (int i = K * S * M - 1; i >= monIndex; i--) {
            getXBySolveopt2(monIndex, i, 0, contractItems[monIndex - 1].x - contractItems[monIndex - 1].x / 10);
        }
//        for (int i = 0; i < 8; i++) {
//            std::cout << "契约项: " << contractItems[i].x << " " << contractItems[i].r << endl;
//            std::cout << endl;
//        }
        monIndex = monTest();
    }

    for (int i = K * S * M - 1; i >= 0; i--) {
        contractItems[i].r = getRFromX(contractItems[i].x, i);
    }


//    solveAll();

    std::cout << "----------- 契约结果 -----------" << endl;
    for (int i = 0; i < 8; i++) {
        userMRS[i].printUser();
        std::cout << "契约项: " << contractItems[i].x << " " << contractItems[i].r << endl;
        std::cout << "用户效用: " << userMRS[i].utility(contractItems[i]) << std::endl;
        UserIConIUtil.push_back(userMRS[i].utility(contractItems[i]));
        std::cout << endl;
    }

    double edgeUtility = 0;
    for (int i = 0; i < 8; i++) {
        edgeUtility += NumOfUser * userMRS[i].distribution *
                (userMRS[i].t * contractItems[i].x * userMRS[i].edge.price_edge - contractItems[i].r - userMRS[i].t * contractItems[i].x * userMRS[i].t * contractItems[i].x);
    }
    std::cout << "边缘效用: " << edgeUtility << std::endl;



    // 测试 IC 约束
    for (int userI = 0; userI < K*S*M; userI++) {
        for (int conI = 0; conI < K*S*M; conI++) {
            if (userI != conI && userMRS[userI].utility(contractItems[conI]) >= UserIConIUtil[userI]) {
                std::cout << "不符合 IC 约束：用户 " << userI << " 契约项 " << conI << endl;
                std::cout << "原效用 " << UserIConIUtil[userI] << " 新效用 " << userMRS[userI].utility(contractItems[conI]) << endl;
                std::cout << endl;
            }
        }
    }


    return 0;
}
