//
// Created by Zebra on 2023/6/2.
//

#ifndef CONTRACT_SETTING_H
#define CONTRACT_SETTING_H
//
// @brief: 参数设置
// @birth: created by zebra on 2023-05-26
// @version: V1
// 用户数量
#define NumOfUser 100 // 用户数量可以指定，但是具体的用户类型信息应该按照分布随机生成

// 定义用户私有信息类型的数量
#define K 2 // 需要卸载的总的任务量
#define S 2 // 延迟敏感度类型数量
#define M 2 // 价格敏感度类型数量

// 定义 EDGE 成本项系数
#define CostItemCoefficient 0.6

#define V 10 // 完成单位任务量得到的收益

#endif //CONTRACT_SETTING_H
