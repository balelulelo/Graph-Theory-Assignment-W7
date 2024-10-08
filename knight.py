#Ali Ridho 5025231162
def is_valid(x, y, n, m, visited):
    return 0 <= x < n and 0 <= y < m and not visited[x][y]

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

n = int(input())
m = int(input())
start_x = int(input())
start_y = int(input())

solve_knight_tour_custom(n, m, start_x, start_y)
