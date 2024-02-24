#include <iostream>
#include <vector>

struct Node {
    int id;
    int degree;
};

int main() {
    const int totalNodes = 500;
    const int minDegree = 1000;

    // 定义节点向量
    std::vector<Node> nodes(totalNodes);

    // 初始化节点，假设每个节点都有相同的自由度
    for (int i = 0; i < totalNodes; ++i) {
        nodes[i].id = i + 1;
        nodes[i].degree = minDegree / 10; // 假设每个节点的自由度是 100
    }

    // 计算总的自由度
    int totalDegrees = minDegree * totalNodes / 10;

    // 调整自由度使每个节点的总和大于 1000
    for (int i = 0; i < totalNodes; ++i) {
        while (nodes[i].degree * 10 <= minDegree) {
            for (int j = 0; j < totalNodes; ++j) {
                if (i != j && nodes[j].degree > minDegree / 10) {
                    nodes[i].degree++;
                    nodes[j].degree--;
                    break;
                }
            }
        }
    }

    // 输出每个节点的ID和自由度
    for (int i = 0; i < totalNodes; ++i) {
        std::cout << "Node " << nodes[i].id << ": Degree " << nodes[i].degree * 10 << std::endl;
    }

    // 验证总的自由度是否大于 1000
    int sum = 0;
    for (int i = 0; i < totalNodes; ++i) {
        sum += nodes[i].degree * 10;
    }
    std::cout << "Total degrees: " << sum << std::endl;

    return 0;
}
