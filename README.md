# Graph Theory W7 Assignment - Group 1

Group 1 Members:

| Name  | NRP |
| ------------- | ------------- |
| Muiz Surya Fata     | 5025231005  |
| Alfa Radithya Fanany   | 5025231008  |
| Muhammad Iqbal Shafarel    | 5025231080 |
| Ali Ridho   | 5025231162|

This repository will give the explanation for several Graph problems such as Travelling Salesman Problem (TSP), Chinese Postman Problem (CPP), and Knight's Tour Problem. The detailed explanations will be explained below, and the related source code is attached in this repository.


# Code Explanation

## 1. Travelling Salesman Problem (TSP)

### Solution Overview

The solution for this code implements a greedy algorithm for solving the Traveling Salesman Problem (TSP). It builds a graph based on input edges and calculates the route that visits each node exactly once. 

### A. Headers

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

### B. Constant and Struct Definitions

```
const int INF = numeric_limits<int>::max();

struct Edge {
    int id, u, v, weight;
};


```
* `INF` = represents a very large number to initialize the adj matrix
* `struct Edge` = A structure to hold the info about an edge in the graph including the edge identifier `id`, one end point and the other `u` and `v`, and the cost of the edge `weight`.

### C. TSP Solver Function

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

### D. Build the Matrix

```
    for (const auto& e : edges) {
        if (adjMatrix[e.u][e.v] > e.weight) {
            adjMatrix[e.u][e.v] = e.weight;
            adjMatrix[e.v][e.u] = e.weight;
            edgeIds[e.u][e.v] = e.id;
            edgeIds[e.v][e.u] = e.id;
        }
    }

```
In this part, the loop iterates over each edge in the `edges` vector. This loop will update the matrix so that it will only choose the path with the cheapest cost for the `adjMatrix`, and it will set the edge ID for both directions.

### E. TSP Logic

```
    vector<bool> visited(n + 1, false);
    visited[startNode] = true;

    int currentNode = startNode;
```
In this part, a boolean vector will track which node have been visited and mark the starting node as visited. It will also keep track of the current node with `currentNode`.

### F. Find the Next Node

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

### G. Completing the Route

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


### H. Main Function

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

### I. Result

after we run the code, the image below is the output of the code:

![image](https://github.com/user-attachments/assets/269b959f-530f-4471-a8bd-a38669c91466)


## 2. Chinese Postman Problem (CPP)
### A. Header
```
#include <bits/stdc++.h>
using namespace std;
#define MAX_NODES 16
#define INF LLONG_MAX
```
- Define for max interval of nodes (16).
- `INF` to make sure for max interval of the path cost.


### B. Initialize variable
```
int n, m; 
long long totalCost = 0; 
long long adj[MAX_NODES][MAX_NODES]; 
int degree[MAX_NODES]; 
vector<int> oddDegreeNodes; 
vector<int> bestPath; 
```
- Using two-dimensional  adjency matrix for the graph.
- Initialize `oddpath` and `bestPath` with vector due to graph with matrix type.

### C. Initializing Graph
```
void initializeGraph()
{
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                adj[i][j] = 0; 
            } 
            else{
                adj[i][j] = INF; 
            }
        }
        degree[i] = 0; 
    }
    cout << "Graph initialized.\n";
}
```
- Function Declaration: void initializeGraph(): Declares a function that returns no value.
Outer Loop (Node Iteration):
- `for (int i = 0; i < n; i++)`: Iterates over each node in the graph, where n is the total number of nodes.
- `for (int j = 0; j < n; j++)`: Iterates through all nodes again for each node i.
Set Distances:
- if `(i == j):` Checks if the current node (i) is the same as the node being compared (j).
- `adj[i][j] = 0;`: Sets the distance from a node to itself to 0.
- `else: If the nodes are different:
adj[i][j] = INF;`: Sets the distance between different nodes to infinity (indicating no direct connection).
- `cout << "Graph initialized.\n";`: indicating the graph has been successfully initialized.

### D. Read the edges
```
void readEdges()
{
    for (int i = 0; i < m; i++){
        int edgeId, x, y;
        long long cost; 
        cin >> edgeId >> x >> y >> cost; 
        x--; y--; 
        adj[x][y] = min(adj[x][y], cost); 
        adj[y][x] = min(adj[y][x], cost);
        totalCost += cost; 
        degree[x]++;
        degree[y]++;
        cout << "Edge added: " << edgeId << " from " << x+1 << " to " << y+1 << " with cost " << cost << "\n";
    }
    cout << "Total cost of edges: " << totalCost << "\n";
}
```
- `void readEdges():` Declares a function that returns no value.
- `for (int i = 0; i < m; i++):` Iterates m times, where m is the total number of edges to be read.
- `int edgeId, x, y;`: Declares variables to store the edge identifier and the two nodes that the edge connects.
- `long long cost;`: Declares a variable to store the cost of the edge.
- `cin >> edgeId >> x >> y >> cost;`: Reads the edge identifier, two nodes (x and y), and the cost from standard input.
- `x--; y--;`: Decrements x and y by 1 to convert the input from 1-based indexing (common in user input) to 0-based indexing (used in programming).
- `adj[x][y] = min(adj[x][y], cost);`: Updates the adjacency matrix at position [x][y] with the minimum cost between the existing cost and the new cost.
- `adj[y][x] = min(adj[y][x], cost);`: Similarly updates the adjacency matrix at position `[y][x]` to ensure the graph remains undirected.
- `totalCost += cost;`: Adds the cost of the current edge to the totalCost variable, which keeps track of the total cost of all edges.
- `degree[x]++;`: Increases the degree of node x by 1 to indicate that an edge is connected to it.
- `degree[y]++;`: Increases the degree of node y by 1.
Output Edge Information:
- `cout << "Edge added: " << edgeId << " from " << x+1 << " to " << y+1 << " with cost " << cost << "\n";`: Outputs showing the edge ID, the nodes it connects (converted back to 1-based indexing), and the cost of the edge.

### E. Using floydWarshall Algortihm to find the minimum cost
```
    for (int k = 0; k < n; k++){
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (dist[i][j] > dist[i][k] + dist[k][j]){
                    dist[i][j] = dist[i][k] + dist[k][j];
                    cout << "Updated distance from " << i+1 << " to " << j+1 << " via " << k+1 << " to " << dist[i][j] << "\n";
                }
            }
        }
    }
```
- The algorithm works by iteratively improving the path distances between pairs of nodes through intermediate nodes.
- For every possible intermediate node, it checks whether the path from node `i`to node `j` can be shortened by going through that intermediate node `k`. If the path through `k` is shorter, the algorithm updates the distance.

### F. Identify Odd Degree
```
void findOddDegreeNodes()
{
    for (int i = 0; i < n; i++){
        if (degree[i] % 2 == 1){
            oddDegreeNodes.push_back(i); 
            cout << "Odd degree node found: " << i+1 << "\n";
        }
    }
}
```
- Identify nodes that have an odd degree in the graph. In the context of graph theory, a node with an odd degree is important when determining Eulerian paths (paths that traverse every edge exactly once).

### G. Finding Mini Cost Function
```
long long findMinCost(const vector<bool> &mask, long long dist[MAX_NODES][MAX_NODES], vector<int> &path) 
{
    bool allNodesRemoved = true;
    for (bool inMask : mask) {
        if (inMask){
            allNodesRemoved = false; 
            break; 
        }
    }
    if (allNodesRemoved) return 0; 

    long long minCost = INF; 
    int firstNode = -1;

    for (int i = 0; i < n; i++){
        if (mask[i]){ 
            firstNode = i; 
            break; 
        }
    }

    for (int secondNode = firstNode + 1; secondNode < n; secondNode++){
        if (mask[secondNode]){ 
            vector<bool> newMask = mask;
            newMask[firstNode] = false; 
            newMask[secondNode] = false; 

            long long newCost = dist[firstNode][secondNode] + 
                                findMinCost(newMask, dist, path);

            if (newCost < minCost){
                minCost = newCost; 
                path.push_back(firstNode); 
                path.push_back(secondNode); 
                cout << "Pairing nodes: " << firstNode + 1 << " and " << secondNode + 1 << " with cost: " << newCost << "\n";
            }
        }
    }

    return minCost; 
}
```
- This loop identifies the first node that is still available (not paired). It updates firstNode to the index of this node.
- This nested loop attempts to pair firstNode with every other available node `(secondNode)` that comes after it in the list.
- A copy of vector is created as `newMask`, and both `firstNode` and `secondNode` are marked as unavailable by setting their corresponding values to false.
- The cost for this particular pairing is calculated by adding the distance between `firstNode` and `secondNode` `(dist[firstNode][secondNode])` to the cost returned by a recursive call to findMinCost with the updated mask (newMask).
- If the calculated `newCost` is less than the current `minCost`, it updates `minCost` and adds `firstNode` and `secondNode` to the path vector. A message is printed to indicate the nodes that have been paired and the associated cost.
- 
### H. Main Function
```
int main() 
{
    cin >> n >> m; 
    initializeGraph(); 
    readEdges(); 

    long long dist[MAX_NODES][MAX_NODES];
    floydWarshall(dist); 

    findOddDegreeNodes(); 

    vector<bool> mask(n, false);
    for (int node : oddDegreeNodes){
        mask[node] = true; 
    }

    int startingPoint;
    cin >> startingPoint; 
    startingPoint--; 

    vector<int> path;

    long long minCost = totalCost+findMinCost(mask, dist, path); 

    cout << "Minimum Cost: " << minCost << "\n";
    return 0; 
}
```
- `initializeGraph();`
Calls the initializeGraph function
- `readEdges();`calls the readEdges function to read the edge data and update the adjacency matrix with edge costs 
- `long long dist[MAX_NODES][MAX_NODES];`
Declares a 2D array dist to store shortest distances between nodes.
`floydWarshall(dist);`calls the *floydWarshall* function to compute the shortest paths between all pairs of nodes and populate the dist array.
- `findOddDegreeNodes();`
Calls the findOddDegreeNodes function to identify and store nodes with odd degrees.
- `vector<bool> mask(n, false);`
Initializes a boolean vector mask of size n, where all values are set to false.`for (int node : oddDegreeNodes){ mask[node] = true; }`Iterates through the oddDegreeNodes vector and sets the corresponding indices in the mask vector to true, marking these nodes as available for pairing.
- `int startingPoint;`
Declares an integer variable to store the starting node for traversal.`cin >> startingPoint;`
Reads the starting point from the user.`startingPoint--;`
Decreases the value of startingPoint by 1 
- `vector<int> path;`
Declares a vector path to store the sequence of paired nodes during the computation.
- `long long minCost = totalCost + findMinCost(mask, dist, path);`
Calls the findMinCost function to calculate the minimum cost of pairing the odd-degree nodes, adding it to the totalCost of edges. The result is stored in minCost.
- `cout << "Minimum Cost: " << minCost << "\n";`
Outputs the calculated minimum cost to the console.

### I. Results
![image](https://github.com/user-attachments/assets/d81ac070-a82d-4ca6-b3ac-dd9b342d3e33)

## 3. Knight's Tour

### Solution Overview

This solution aims to solve the Knight's Tour problem using a backtracking approach with a predefined set of moves in an anticlockwise pattern. The goal is to visit every square on the chessboard exactly once using the knight's legal moves.

### A. Valid Move

```
def is_valid(x, y, n, m, visited):
    return 0 <= x < n and 0 <= y < m and not visited[x][y]

```
This part is a function to check if a move is valid and within the boundary of the board.

The parameters are as follows:
* `x, y` = the coordinate for next move
* `n` = rows
* `m` = columns
* `visited` = A matrix to track if a cell is already visited

### B. Create Graph

```
def create_custom_graph(n, m):
    custom_moves = [(2, -1), (1, -2), (-1, -2), (-2, -1), (-2, 1), (-1, 2), (1, 2), (2, 1)]
    graph = {}

    for i in range(n):
        for j in range(m):
            neighbors = []
            for move in custom_moves:
                new_x, new_y = i + move[0], j + move[1]
                if 0 <= new_x < n and 0 <= new_y < m:
                    neighbors.append((new_x, new_y))
            graph[(i, j)] = neighbors

    return graph
```

In this part, the function builds a graph where each position is a node, and the neighbors are the possible moves. For each cell on the board, the function checks if moving to a neighboring cell using any of the knight's moves is valid. If valid, that neighboring cell is added to the node’s neighbors.

The purpose is to precompute all valid moves from each square and organize them into a graph structure, allowing efficient traversal during backtracking.

### C. Knight Tour Graph

```
def knight_tour_graph(n, m, start_x, start_y, graph, visited, path, move_count):
    if move_count == n * m:  
        return True

    for neighbor in graph[(start_x, start_y)]:
        new_x, new_y = neighbor
        if is_valid(new_x, new_y, n, m, visited):

            visited[new_x][new_y] = True
            path.append((new_x, new_y))

            if knight_tour_graph(n, m, new_x, new_y, graph, visited, path, move_count + 1):
                return True

            visited[new_x][new_y] = False
            path.pop()

    return False
```
This part is a recursive backtracking function that attempts to find a path where the knight visits every square on the board exactly once. It tries each valid move from the current square (from the graph) in the anticlockwise order.

For each valid move, it will:
* Marks the square as visited.
* Adds it to the path.
* Recursively calls itself to continue the tour from the new square.
* If a solution is found (all squares are visited), the function returns True. Otherwise, it backtracks by unmarking the square and trying other possible moves.

The purpose is to explore all possible knight paths, making decisions recursively and backtracking when needed until a solution is found.

### D. Solve Knight Tour

```
def solve_knight_tour_custom(n, m, start_x, start_y):
    graph = create_custom_graph(n, m)

    visited = [[False for _ in range(m)] for _ in range(n)]
    path = [(start_x, start_y)]
    visited[start_x][start_y] = True

    if knight_tour_graph(n, m, start_x, start_y, graph, visited, path, 1):

        for coord in path:
            print(coord)
    else:
        print("No solution exists")
```

In this part, it initialize the process by creating a graph and setting up the visited array.
it will call `knight_tour_graph` to start the tour from the given starting position `(start_x, start_y)`.
If the knight successfully visits all squares, it prints the resulting path. Otherwise, it prints `"No solution exists"`.

The purpose is to be the driver function that sets up the chessboard, the knight's starting position, and initiates the recursive backtracking process.

### E. Main Function

```
n, m = 5, 5  
start_x, start_y = 2, 2  

solve_knight_tour_custom(n, m, start_x, start_y)
```

this part acts as a main function. It sets the dimension of the board to `5 x 5` and the starting position to `2 x 2`. It will then call `solve_knight_tour_custom` to calculate the result of this code

### F. Result

after we run the code, the following image is the output result:
-Once the knight visits all squares, the solve_knight_tour_custom function prints each coordinate from the path, one by one, resulting in the knight’s tour being displayed in the order it was executed.

![image](https://github.com/user-attachments/assets/741a3b22-e619-454f-85fc-de7cb6b17d56)

