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

### Solution Overview

This solution aims to solve the Knight's Tour problem using a backtracking approach with a predefined set of moves in an anticlockwise pattern. The goal is to visit every square on the chessboard exactly once using the knight's legal moves.

### Valid Move

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

### Create Graph

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

### Knight Tour Graph

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

### Solve Knight Tour

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

### Main Function

```
n, m = 5, 5  
start_x, start_y = 2, 2  

solve_knight_tour_custom(n, m, start_x, start_y)
```

this part acts as a main function. It sets the dimension of the board to `5 x 5` and the starting position to `2 x 2`. It will then call `solve_knight_tour_custom` to calculate the result of this code

### Result

after we run the code, the following image is the output result:

![image](https://github.com/user-attachments/assets/741a3b22-e619-454f-85fc-de7cb6b17d56)

