#include <bits/stdc++.h>
using namespace std;
#define MAX_NODES 16
#define INF LLONG_MAX

int n, m; 
long long totalCost = 0; 
long long adj[MAX_NODES][MAX_NODES]; 
int degree[MAX_NODES]; 
vector<int> oddDegreeNodes; 
vector<int> bestPath; 

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
        cout << "Edge added: " << edgeId << " from " << 
        x+1 << " to " << y+1 << " with cost " << cost << "\n";
    }
    cout << "Total cost of edges: " << totalCost << "\n";
}

void floydWarshall(long long dist[MAX_NODES][MAX_NODES])
{
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            dist[i][j] = adj[i][j];
        }
    }
    
    cout << "Initial distance matrix:\n";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (dist[i][j] == INF){
                cout << "INF ";
            } 
            else{
                cout << dist[i][j] << " "; // Output directly
            }

        }
        cout << "\n";
    }

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

    cout << "Final distance matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (dist[i][j] == INF){
                cout << "INF ";
            } 
            else {
                cout << dist[i][j] << " "; // Output directly
            }

        }
        cout << "\n";
    }
}

void findOddDegreeNodes()
{
    for (int i = 0; i < n; i++){
        if (degree[i] % 2 == 1){
            oddDegreeNodes.push_back(i); 
            cout << "Odd degree node found: " << i+1 << "\n";
        }
    }
}

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

            long long newCost = dist[firstNode][secondNode] + findMinCost(newMask, dist, path);

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
    vector<int> fullRoute;
    return 0; 
}
