# Graph Theory W7 Assignment - Group 1

Members:

| Name  | NRP |
| ------------- | ------------- |
| Muiz Surya Fata     | 5025231005  |
| Alfa Radithya Fanany   | 5025231008  |
| Muhammad Iqbal Shafarel    | 5025231080 |
| Ali Ridho   | 5025231162|

# Code Explanation

## A. Travelling Salesman Problem (TSP)

### Headers

```
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

using namespace std;

```

These are the headers we used in this code:
* `#include <iostream>` = Standard input output operations
* `#include <vector>` = For dynamic arrays
* `#include <limits>` = Define numerical limit
* `#include <algorithm>` = using various algorithm (in this case = reverse)
* `using namespace std;` = removing the `std::` prefix

### Constant and Struct Definitions

```
const int INF = numeric_limits<int>::max();

struct Edge {
    int id, u, v, weight;
};


```
* `INF` = represents a very large number to initialize the adj matrix
* `struct Edge` = A structure to hold the info about an edge in the graph including the edge identifier `id`, one end point and the other `u` and `v`, and the cost of the edge `weight`.

### TSP Solver Function

```
int solveTSP(int n, int m, const vector<Edge>& edges, int startNode) {

    vector<vector<int>> adjMatrix(n + 1, vector<int>(n + 1, INF));
    vector<int> edgeRoute;
    int totalCost = 0;

    vector<vector<int>> edgeIds(n + 1, vector<int>(n + 1, -1));


```

This part initializes the Data Structure for this code

* `solveTSP` = A function declaration for solving the problem
* `adjMatrix` = a 2D vector that represents the adjacency matrix of the graph
* `edgeRoute` = store the sequences of `id` in order
* `totalCost` = accumulate the total cost

### Build the Matrix

```
    for (const auto& e : edges) {
        adjMatrix[e.u][e.v] = e.weight;
        adjMatrix[e.v][e.u] = e.weight;
        edgeIds[e.u][e.v] = e.id;
        edgeIds[e.v][e.u] = e.id;
    }

```
In this part, the loop iterates over each edge in the `edges` vector. It will set the weight of both direction in `adjMatrix`, and it will set the edge ID for both directions.

### TSP Logic

```
    vector<bool> visited(n + 1, false);
    visited[startNode] = true;

    int currentNode = startNode;
```
In this part, a boolean vector will track which node have been visited and mark the starting node as visited. It will also keep track of the current node with `currentNode`.

### Find the Next Node

```
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
```

In this part, the steps are:
* `nextNode` and `minCost` are initialized to keep track the next node to visit and also the minimum cost to go there.
* The inner loop check each node to find the unvisited node with the smallest edge cost
* if there is no unvisited node, the loop breaks

The accumulated minimum cost will be added to `totalCost`, and the ID of the edge will be added to `edgeRoute`.

### Completing the Route

```
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
```

In this part, the cost to return to starting node is added to the total cost and the ID of the returning edge is added to the edge route.

after that part, the `totalCost` is printed, and the `edgeRoute` is reversed to display the correct order of the traversal


### Main Function

```
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
```
In the main function, it handles the input such as:
* nodes `n`, edges `m`, and the edges themselves to `edges` vector
* the starting node

after reading the inputs, it will call the `solveTSP` function to calculate the cost and the route

### Result

after we run the code, the image below is the output of the code:

![image](https://github.com/user-attachments/assets/269b959f-530f-4471-a8bd-a38669c91466)


## B. Chinese Postman Problem (CPP)

## C. Knights Tours
