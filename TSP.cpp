#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

using namespace std;

const int INF = numeric_limits<int>::max();

struct Edge {
    int id, u, v, weight;
};

int solveTSP(int n, int m, const vector<Edge>& edges, int startNode) {
    vector<vector<int>> adjMatrix(n + 1, vector<int>(n + 1, INF));
    vector<int> edgeRoute;
    int totalCost = 0;

    vector<vector<int>> edgeIds(n + 1, vector<int>(n + 1, -1));
    for (const auto& e : edges) {
        adjMatrix[e.u][e.v] = e.weight;
        adjMatrix[e.v][e.u] = e.weight;
        edgeIds[e.u][e.v] = e.id;
        edgeIds[e.v][e.u] = e.id;
    }

    vector<bool> visited(n + 1, false);
    visited[startNode] = true;

    int currentNode = startNode;

    for (int count = 1; count < n; ++count) {
        int nextNode = -1;
        int minCost = INF;

        for (int i = 1; i <= n; ++i) {
            if (!visited[i] && adjMatrix[currentNode][i] < minCost) {
                minCost = adjMatrix[currentNode][i];
                nextNode = i;
            }
        }

        if (nextNode == -1) break;  

        totalCost += minCost;
        edgeRoute.push_back(edgeIds[currentNode][nextNode]);  
        visited[nextNode] = true;
        currentNode = nextNode;
    }


    totalCost += adjMatrix[currentNode][startNode];
    edgeRoute.push_back(edgeIds[currentNode][startNode]); 

    reverse(edgeRoute.begin(), edgeRoute.end());

    cout << "Cost: " << totalCost << endl;
    cout << "Route: ";
    for (size_t i = 0; i < edgeRoute.size(); ++i) {
        cout << edgeRoute[i];
        if (i != edgeRoute.size() - 1) cout << ", ";
    }
    cout << endl;

    return totalCost;
}

int main() {
    int n, m, startNode;
    cin >> n; 
    cin >> m;
    vector<Edge> edges(m);

    for (int i = 0; i < m; ++i) {
        cin >> edges[i].id >> edges[i].u >> edges[i].v >> edges[i].weight;
    }

    cin >> startNode;
    solveTSP(n, m, edges, startNode);

    return 0;
}
