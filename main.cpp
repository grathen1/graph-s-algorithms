#include <SFML/Graphics.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <cmath>
#include <limits>
#include <functional>
#include <set>

int state = 0;

class Graph {
private:
    std::vector<std::vector<bool>> connections;
    std::vector<bool> activeVertices;
    std::vector<std::vector<int>> weights;
    std::vector<std::vector<bool>> activeEdges;
    std::vector<std::vector<int>> weightConnections;
    float centerX, centerY;
    float radius;

public:
    Graph() : centerX(460), centerY(400), radius(300) {}

    bool read_from_file(const std::string& filename, bool weighted = true) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Unable to open file: " << filename << std::endl;
            return false;
        }

        int from, to, weight;
        int max_node = 0;

        while (file >> from >> to >> weight) {
            max_node = std::max(max_node, std::max(from, to));
        }

        // Resize vectors to accommodate the maximum index found
        connections.resize(max_node + 1, std::vector<bool>(max_node + 1, false));
        weights.resize(max_node + 1, std::vector<int>(max_node + 1, 0));  // Store weights here
        activeEdges.resize(max_node + 1, std::vector<bool>(max_node + 1, false));
        activeVertices.resize(max_node + 1, false);

        // Reset file stream to the beginning to read the data again
        file.clear();
        file.seekg(0, std::ios::beg);

        while (file >> from >> to >> weight) {
            connections[from][to] = true;
            connections[to][from] = true; // For undirected graphs
            weights[from][to] = weight;
            weights[to][from] = weight; // Assuming the graph is undirected
            activeVertices[from] = true;
            activeVertices[to] = true;
        }
        file.close();
        return true;
    }

    void draw_edge(size_t from, size_t to, sf::RenderWindow &window, sf::Color color, const sf::Font &font) const {
        if (!connections[from][to]) return; // Only draw if there's an edge

        float angleFrom = 2 * M_PI * from / connections.size();
        float xFrom = centerX + radius * cos(angleFrom);
        float yFrom = centerY + radius * sin(angleFrom);

        float angleTo = 2 * M_PI * to / connections.size();
        float xTo = centerX + radius * cos(angleTo);
        float yTo = centerY + radius * sin(angleTo);

        float midX = (xFrom + xTo) / 2;
        float midY = (yFrom + yTo) / 2;

        sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(xFrom, yFrom), color),
                sf::Vertex(sf::Vector2f(xTo, yTo), color)
        };
        window.draw(line, 2, sf::Lines);

//        if (weights[from][to] != 0) {
//            sf::Text weightText;
//            weightText.setFont(font);
//            weightText.setString(std::to_string(weights[from][to]));
//            weightText.setCharacterSize(15);
//            weightText.setFillColor(sf::Color::Blue);
//            weightText.setPosition(midX, midY);
//            window.draw(weightText);
//        }
    }

    void display_vertex(int vertex, sf::RenderWindow &window, sf::Font &font, sf::Color color) const {
        if (!activeVertices[vertex]) return;  // Only display active vertices

        float angleStep = 2 * M_PI / connections.size();
        float angle = vertex * angleStep;

        float x = centerX + radius * cos(angle) - 10;
        float y = centerY + radius * sin(angle) - 10;

        sf::CircleShape shape(10);
        shape.setFillColor(color);
        shape.setOutlineColor(sf::Color::Black);
        shape.setOutlineThickness(1);
        shape.setPosition(x, y);

        sf::Text text(std::to_string(vertex), font, 15);
        text.setFillColor(sf::Color::Black);
        text.setPosition(x + 10 - text.getLocalBounds().width / 2, y + 10 - text.getLocalBounds().height);

        window.draw(shape);
        window.draw(text);
    }

    void display(sf::RenderWindow &window, sf::Font &font) {
        for (size_t i = 0; i < connections.size(); ++i) {
            if (!activeVertices[i]) continue;
            for (size_t j = 0; j < connections[i].size(); ++j) {
                if (connections[i][j]) {
                    draw_edge(i, j, window, activeEdges[i][j] ? sf::Color::Red : sf::Color::Black, font);
                }
            }
        }
        for (size_t i = 0; i < connections.size(); ++i) {
            display_vertex(i, window, font, sf::Color::Yellow);
        }
    }

    void bfs(int start_vertex, sf::RenderWindow &window, sf::Font &font) {
        std::vector<bool> visited(connections.size(), false);
        std::queue<int> queue;
        visited[start_vertex] = true;
        queue.push(start_vertex);

        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();

            window.clear(sf::Color::White);
            display(window, font);
            display_vertex(current, window, font, sf::Color::Red);
            window.display();
            sf::sleep(sf::milliseconds(500));

            for (int i = 0; i < connections[current].size(); ++i) {
                if (connections[current][i] && !visited[i]) {
                    visited[i] = true;
                    queue.push(i);
                    activeEdges[current][i] = true;
                    activeEdges[i][current] = true;

                    window.clear(sf::Color::White);
                    display(window, font);
                    display_vertex(current, window, font, sf::Color::Red);
//                    display_vertex(i, window, font, sf::Color::Green);
                    window.display();
                    sf::sleep(sf::milliseconds(500));  // Pause to observe the connection
                }
            }
        }

        window.clear(sf::Color::White);
        display(window, font);
        window.display();
    }

    void dfs(int start_point, sf::RenderWindow &window, sf::Font &font) {
        if (connections.empty() || start_point >= connections.size()) {
            std::cout << "Invalid start point or empty graph." << std::endl;
            return;
        }

        std::vector<bool> visited(connections.size(), false);
        std::stack<int> stack;
        visited[start_point] = true;
        stack.push(start_point);

        while (!stack.empty()) {
            int current = stack.top();
            stack.pop();

            // Ensure that every step is visually updated
            window.clear(sf::Color::White);
            display(window, font);
            display_vertex(current, window, font, sf::Color::Red);
            window.display();
            sf::sleep(sf::milliseconds(500));  // Visual delay for better understanding

            for (int i = 0; i < connections[current].size(); i++) {
                if (connections[current][i] && !visited[i]) {
                    visited[i] = true;
                    stack.push(i);

                    activeEdges[current][i] = true;
                    activeEdges[i][current] = true;

                    window.clear(sf::Color::White);
                    display(window, font);
                    display_vertex(current, window, font, sf::Color::Red);
                    display_vertex(i, window, font, sf::Color::Green);
                    window.display();
                    sf::sleep(sf::milliseconds(500));
                }
            }
        }

        window.clear(sf::Color::White);
        display(window, font);
        window.display();
    }

    void resetVisualization() {
        for (size_t i = 0; i < activeEdges.size(); i++) {
            for (size_t j = 0; j < activeEdges[i].size(); j++) {
                activeEdges[i][j] = false;  // Reset all edges to inactive
            }
        }
    }

    void add_edge(int from, int to, int weight) {
        connections[from][to] = true;
        connections[to][from] = true;
        weights[from][to] = weight;
        weights[to][from] = weight;  // Ensure weight is symmetrical for undirected graph
        activeVertices[from] = true;
        activeVertices[to] = true;
    }

    int node_count() const {
        return connections.size();
    }

    bool is_connected(int from, int to) const {
        return connections[from][to] != 0;
    }

    int weight(int from, int to) const {
        return connections[from][to];
    }

    void prims_min_tree(sf::RenderWindow &window, sf::Font &font, int start_point = 0) {
        int n = node_count();
        if (n == 0) return;

        std::vector<bool> visited(n, false);
        std::priority_queue<std::pair<int, std::pair<int, int>>,
                std::vector<std::pair<int, std::pair<int, int>>>,
                std::greater<std::pair<int, std::pair<int, int>>>> pq;

        for (int i = 0; i < n; i++) {
            if (connections[start_point][i]) {
                pq.push({weights[start_point][i], {start_point, i}});
            }
        }

        visited[start_point] = true;
        int mst_cost = 0;

        window.clear(sf::Color::White);
        display(window, font);
        display_vertex(start_point, window, font, sf::Color::Green);
        window.display();
        sf::sleep(sf::milliseconds(500));

        while (!pq.empty()) {
            auto edge = pq.top();
            pq.pop();
            int w = edge.first;
            int u = edge.second.first;
            int v = edge.second.second;

            if (visited[v]) continue;

            visited[v] = true;
            mst_cost += w;
            activeEdges[u][v] = true;
            activeEdges[v][u] = true;

            window.clear(sf::Color::White);
            display(window, font);
            for (int node = 0; node < n; node++) {
                if (visited[node]) {
                    display_vertex(node, window, font, sf::Color::Green);  // Highlight the visited nodes
                } else {
                    display_vertex(node, window, font, sf::Color::Yellow);  // Non-visited nodes
                }
            }
            window.display();
            sf::sleep(sf::milliseconds(500));

            for (int i = 0; i < n; i++) {
                if (!visited[i] && connections[v][i]) {
                    pq.push({weights[v][i], {v, i}});
                }
            }
        }

        std::cout << "Total weight of the MST: " << mst_cost << std::endl;
    }

    void dijkstra_path(sf::RenderWindow &window, sf::Font &font, int from, int to) {
        int n = node_count();
        std::vector<int> distances(n, std::numeric_limits<int>::max());
        std::vector<int> previous(n, -1);
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

        distances[from] = 0;
        pq.push({0, from});

        resetVisualization();  // Reset visualization, especially the edges

        while (!pq.empty() && window.isOpen()) {
            int current_distance = pq.top().first;
            int current_vertex = pq.top().second;
            pq.pop();

            if (current_vertex == to) break;

            for (int v = 0; v < n; v++) {
                if (connections[current_vertex][v]) {
                    int weight = weights[current_vertex][v];
                    int distanceThroughU = current_distance + weight;

                    if (distanceThroughU < distances[v]) {
                        distances[v] = distanceThroughU;
                        previous[v] = current_vertex;
                        pq.push({distanceThroughU, v});
                    }
                }
            }

            // Visualize the current step
            window.clear(sf::Color::White);
            display(window, font);  // Display all vertices and edges with current active status
            display_vertex(current_vertex, window, font, sf::Color::Red);
            window.display();
            sf::sleep(sf::milliseconds(500)); // Pause for visualization effect
        }

        // Highlight only the edges in the final shortest path
        std::vector<int> path;
        for (int at = to; at != -1; at = previous[at]) {
            if (previous[at] != -1) {
                activeEdges[at][previous[at]] = true;
                activeEdges[previous[at]][at] = true;
            }
            path.push_back(at);
        }
        std::reverse(path.begin(), path.end());

        // Final visualization of the path
        window.clear(sf::Color::White);
        display(window, font);  // Redraw all vertices and now correctly highlight the path
        for (int node : path) {
            display_vertex(node, window, font, sf::Color::Green);  // Highlight the path nodes
        }
        window.display();
        sf::sleep(sf::milliseconds(1000)); // Final visualization pause

        std::cout << "Path: ";
        for (int node : path) {
            std::cout << node << " ";
        }
        std::cout << "\nCost: " << distances[to] << std::endl;
    }


    int heuristic(int from, int to) const {
        return std::abs(from - to);
    }

    void astar_path(sf::RenderWindow &window, sf::Font &font, int from, int to) {
        int n = node_count();
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> open_set;
        std::vector<int> g_costs(n, std::numeric_limits<int>::max());
        std::vector<bool> closed_set(n, false);
        std::vector<int> prev(n, -1);

        g_costs[from] = 0;
        open_set.push({heuristic(from, to), from});

        // Reset visualization at the start
        resetVisualization();

        while (!open_set.empty()) {
            int current = open_set.top().second;
            open_set.pop();

            if (closed_set[current]) continue;
            closed_set[current] = true;

            // Visualization
            window.clear(sf::Color::White);
            display(window, font);
            display_vertex(current, window, font, sf::Color::Red);
            window.display();
            sf::sleep(sf::milliseconds(500)); // Shorter pause for quicker visualization

            if (current == to) {
                break;
            }

            for (int neighbor = 0; neighbor < n; ++neighbor) {
                if (!is_connected(current, neighbor) || closed_set[neighbor]) continue;

                int tentative_g_cost = g_costs[current] + weights[current][neighbor];
                if (tentative_g_cost < g_costs[neighbor]) {
                    prev[neighbor] = current;
                    g_costs[neighbor] = tentative_g_cost;
                    int f_cost = tentative_g_cost + heuristic(neighbor, to);
                    open_set.push({f_cost, neighbor});
                }
            }
        }

        // Extract the path if it exists and visualize it
        if (g_costs[to] != std::numeric_limits<int>::max()) {
            std::vector<int> path;
            for (int at = to; at != -1; at = prev[at]) {
                path.push_back(at);
            }
            std::reverse(path.begin(), path.end());

            // Redraw the graph with the final path
            window.clear(sf::Color::White);
            resetVisualization();
            for (size_t i = 0; i < path.size(); i++) {
                int node = path[i];
                display_vertex(node, window, font, sf::Color::Green);
                if (i > 0) {
                    activeEdges[path[i - 1]][node] = true;
                    activeEdges[node][path[i - 1]] = true; // For undirected graphs
                }
            }
            display(window, font); // Display the final state of the graph with the path highlighted
            window.display();
            sf::sleep(sf::milliseconds(1000)); // Pause to show the final path

            // Output the path and total cost
            std::cout << "Path: ";
            for (int node : path) {
                std::cout << node << " ";
            }
            std::cout << "\nTotal Cost: " << g_costs[to] << std::endl;
        } else {
            std::cout << "No path exists." << std::endl;
        }
    }

    // Place this function within your Graph class
    int fordFulkerson(sf::RenderWindow &window, sf::Font &font, int source, int sink) {
        std::vector<std::vector<int>> residualGraph(weights); // Copy weights into residual capacities
        std::vector<int> parent(connections.size(), -1);

        int max_flow = 0; // There is no flow initially

        while (_bfs(window, font, source, sink, parent, residualGraph)) {
            int path_flow = std::numeric_limits<int>::max();
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                path_flow = std::min(path_flow, residualGraph[u][v]);
            }

            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                residualGraph[u][v] -= path_flow;
                residualGraph[v][u] += path_flow;

                // Visualization of flow along the path
                display(window, font);
                draw_edge(u, v, window, sf::Color::Cyan, font); // Draw path in Cyan
                window.display();
                sf::sleep(sf::milliseconds(500)); // Pause for visualization effect
            }
            max_flow += path_flow;
        }
        std::cout << "The maximum possible flow is " << max_flow << std::endl;

        return max_flow;
    }

    // Utility function for Ford-Fulkerson algorithm, BFS to find path from source to sink
    bool _bfs(sf::RenderWindow &window, sf::Font &font, int s, int t, std::vector<int>& parent, const std::vector<std::vector<int>>& residualGraph) {
        std::vector<bool> visited(parent.size(), false);
        std::queue<int> queue;
        queue.push(s);
        visited[s] = true;
        parent[s] = -1;

        while (!queue.empty()) {
            int u = queue.front();
            queue.pop();

            for (int v = 0; v < connections.size(); v++) {
                if (!visited[v] && residualGraph[u][v] > 0) {
                    // Visualization for the BFS search path
                    display(window, font);
                    draw_edge(u, v, window, sf::Color::Green, font); // Draw path in Green
                    window.display();
                    sf::sleep(sf::milliseconds(500)); // Pause for visualization effect

                    queue.push(v);
                    parent[v] = u;
                    visited[v] = true;
                    if (v == t) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    int _dijkstra_step(std::vector<int>& distances, std::vector<int>& parent, std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>& q, std::set<int>& checked_nodes) {
        int buff = q.top().second;
        q.pop();
        checked_nodes.insert(buff);
        for (int i = 0; i < node_count(); i++) {
            if (is_connected(buff, i)) {
                int weight = weights[buff][i];
                if (distances[i] > distances[buff] + weight) {
                    distances[i] = distances[buff] + weight;
                    parent[i] = buff;
                    q.push({distances[i], i});
                }
            }
        }
        return buff;
    }

    struct bidirect_result {
        std::vector<int> path;
        std::vector<int> checked_nodes;

        bidirect_result() : path(), checked_nodes() {}
        bidirect_result(std::vector<int> path, std::vector<int> checked_nodes) :
                path(path), checked_nodes(checked_nodes) {}
    };

    bidirect_result bidirect_dijkstra_path(int from, int to, sf::RenderWindow &window, sf::Font &font) {
        int n = node_count();
        if (n == 0) return {};

        std::set<int> checked_nodes;
        std::vector<int> distances_from(n, std::numeric_limits<int>::max());
        std::vector<int> distances_to(n, std::numeric_limits<int>::max());
        std::vector<int> parent_from(n, -1);
        std::vector<int> parent_to(n, -1);
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> q_from;
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> q_to;

        distances_from[from] = 0;
        q_from.push({0, from});
        distances_to[to] = 0;
        q_to.push({0, to});
        int middle = -1;

        while (!q_from.empty() && !q_to.empty()) {
            int buff = q_from.top().second;
            q_from.pop();
            checked_nodes.insert(buff);
            if (distances_to[buff] != std::numeric_limits<int>::max()) {
                middle = buff;
                break;
            }

            for (int i = 0; i < n; i++) {
                if (is_connected(buff, i)) {
                    int weight = weights[buff][i];
                    if (distances_from[i] > distances_from[buff] + weight) {
                        distances_from[i] = distances_from[buff] + weight;
                        parent_from[i] = buff;
                        q_from.push({distances_from[i], i});

                        // Visualization
                        window.clear(sf::Color::White);
                        display(window, font);  // Redraw all vertices and edges with current active status
                        display_vertex(buff, window, font, sf::Color::Red);  // Highlight the current node
                        display_vertex(i, window, font, sf::Color::Green);  // Highlight the neighbor node
                        window.display();
                        sf::sleep(sf::milliseconds(500)); // Pause for visualization effect
                    }
                }
            }

            buff = q_to.top().second;
            q_to.pop();
            checked_nodes.insert(buff);
            if (distances_from[buff] != std::numeric_limits<int>::max()) {
                middle = buff;
                break;
            }

            for (int i = 0; i < n; i++) {
                if (is_connected(buff, i)) {
                    int weight = weights[buff][i];
                    if (distances_to[i] > distances_to[buff] + weight) {
                        distances_to[i] = distances_to[buff] + weight;
                        parent_to[i] = buff;
                        q_to.push({distances_to[i], i});

                        // Visualization
                        window.clear(sf::Color::White);
                        display(window, font);  // Redraw all vertices and edges with current active status
                        display_vertex(buff, window, font, sf::Color::Red);  // Highlight the current node
                        display_vertex(i, window, font, sf::Color::Green);  // Highlight the neighbor node
                        window.display();
                        sf::sleep(sf::milliseconds(500)); // Pause for visualization effect
                    }
                }
            }
        }

        if (middle == -1) return {};

        std::vector<int> ret;
        int cur = middle;
        while (cur != to) {
            ret.push_back(cur);
            cur = parent_to[cur];
        }
        ret.push_back(to);

        std::reverse(ret.begin(), ret.end());
        cur = parent_from[middle];
        while (cur != from) {
            ret.push_back(cur);
            cur = parent_from[cur];
        }
        ret.push_back(from);

        std::reverse(ret.begin(), ret.end());

        window.clear(sf::Color::White);
        display(window, font);  // Redraw all vertices and edges with current active status
        for (size_t i = 0; i < ret.size() - 1; i++) {
            int node1 = ret[i];
            int node2 = ret[i + 1];
            activeEdges[node1][node2] = true;  // Mark edges in the final path as active
            activeEdges[node2][node1] = true;  // For undirected graphs
        }
        window.display();
        sf::sleep(sf::milliseconds(1000));

        // Print out the path and total distance
        int total_distance = 0;
        for (int i = 0; i < ret.size() - 1; i++) {
            int node1 = ret[i];
            int node2 = ret[i + 1];
            total_distance += weights[node1][node2];
        }
        std::cout << "Total Distance: " << total_distance << std::endl;

        return {ret, {checked_nodes.begin(), checked_nodes.end()}};
    }
    bidirect_result bidirect_astar_path(sf::RenderWindow &window, sf::Font &font, int from, int to) {
        const int n = node_count();
        if (n == 0) return {};
        if (to == -1) to = n - 1;

        std::set<int> checked_nodes;
        std::vector<int> distances_from(n, std::numeric_limits<int>::max()), distances_to(n, std::numeric_limits<int>::max()),
                parent_from(n, -1), parent_to(n, -1);
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> q_from, q_to;

        distances_from[from] = 0;
        q_from.push({0, from});
        distances_to[to] = 0;
        q_to.push({0, to});
        int middle = -1;
        int buff;
        while (!q_from.empty() && !q_to.empty()) {
            buff = _dijkstra_step(distances_from, parent_from, q_from, checked_nodes);
            if (distances_to[buff] != std::numeric_limits<int>::max()) {
                middle = buff;
                break;
            }

            buff = _dijkstra_step(distances_to, parent_to, q_to, checked_nodes);
            if (distances_from[buff] != std::numeric_limits<int>::max()) {
                middle = buff;
                break;
            }
        }

        if (middle == -1) return {};
        int prev_d = std::numeric_limits<int>::max();
        std::pair<int, int> connection{-1, -1};
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (distances_to[j] == std::numeric_limits<int>::max() || distances_from[i] == std::numeric_limits<int>::max()) continue;
                if (!is_connected(i, j)) continue;
                int cur_d = distances_to[j] + distances_from[i] + weight(i, j) + heuristic(i, to) - heuristic(j, to);
                if (prev_d == std::numeric_limits<int>::max() || prev_d > cur_d) {
                    connection = {i, j};
                    prev_d = cur_d;
                }
            }
        }

        if (prev_d == std::numeric_limits<int>::max()) return {};

        std::vector<int> ret;
        int cur = connection.second;
        while (cur != to) {
            ret.push_back(cur);
            cur = parent_to[cur];
        }
        ret.push_back(to);

        std::reverse(ret.begin(), ret.end());
        cur = connection.first;
        while (cur != from) {
            ret.push_back(cur);
            cur = parent_from[cur];
        }
        ret.push_back(from);

        std::reverse(ret.begin(), ret.end());

        // Calculate total distance
        int total_distance = 0;
        for (int i = 0; i < ret.size() - 1; i++) {
            int node1 = ret[i];
            int node2 = ret[i + 1];
            total_distance += weights[node1][node2];
        }

        // Visualization
        window.clear(sf::Color::White);
        for (size_t i = 0; i < ret.size() - 1; i++) {
            int node1 = ret[i];
            int node2 = ret[i + 1];
            activeEdges[node1][node2] = true;  // Mark edges in the final path as active
            activeEdges[node2][node1] = true;  // For undirected graphs
        }
        display(window, font);  // Redraw all vertices and edges with current active status
        for (int node : checked_nodes) {
            display_vertex(node, window, font, sf::Color::Green);  // Mark checked nodes in green
        }
        window.display();
        sf::sleep(sf::milliseconds(1000));

        std::cout << "Total Distance: " << total_distance << std::endl;

        return {ret, {checked_nodes.begin(), checked_nodes.end()}};
    }
};

void menu(sf::RenderWindow& window, sf::Font& font, bool is_main){
    window.clear(sf::Color::White);
    if (is_main) {
        sf::Text name;
        name.setFont(font);
        name.setString("Graph`s visualizer");
        name.setCharacterSize(50);
        name.setFillColor(sf::Color::Black);
        name.setPosition(10 + (870 - name.getLocalBounds().width) / 2,
                         10 + (500 - name.getLocalBounds().height) / 2 - 10);
        window.draw(name);

        sf::Text author;
        author.setFont(font);
        author.setString("made by Rita Kulesh");
        author.setCharacterSize(15);
        author.setFillColor(sf::Color::Black);
        author.setPosition(10 + (170 - author.getLocalBounds().width) / 2,
                           780 + (10 - author.getLocalBounds().height) / 2 - 10);
        window.draw(author);

        sf::RectangleShape bfsButt(sf::Vector2f(200, 50));
        bfsButt.setFillColor(sf::Color::Magenta);
        bfsButt.setPosition(30, 300);
        window.draw(bfsButt);

        sf::Text bfsText;
        bfsText.setFont(font);
        bfsText.setString("BFS");
        bfsText.setCharacterSize(25);
        bfsText.setFillColor(sf::Color::White);
        bfsText.setPosition(30 + (200 - bfsText.getLocalBounds().width) / 2 - 10,
                            300 + (60 - bfsText.getLocalBounds().height) / 2 - 10);
        window.draw(bfsText);

        sf::RectangleShape dijButt(sf::Vector2f(200, 50));
        dijButt.setFillColor(sf::Color::Magenta);
        dijButt.setPosition(30, 370);
        window.draw(dijButt);

        sf::Text dijText;
        dijText.setFont(font);
        dijText.setString("Dijkstra");
        dijText.setCharacterSize(25);
        dijText.setFillColor(sf::Color::White);
        dijText.setPosition(30 + (200 - bfsText.getLocalBounds().width) / 2 - 10,
                            370 + (60 - bfsText.getLocalBounds().height) / 2 - 10);
        window.draw(dijText);

        sf::RectangleShape ddButt(sf::Vector2f(200, 50));
        ddButt.setFillColor(sf::Color::Magenta);
        ddButt.setPosition(30, 440);
        window.draw(ddButt);

        sf::Text ddText;
        ddText.setFont(font);
        ddText.setString("Bidirectional Dijkstra");
        ddText.setCharacterSize(25);
        ddText.setFillColor(sf::Color::White);
        ddText.setPosition(30 + (200 - ddText.getLocalBounds().width) / 2,
                            440 + (60 - ddText.getLocalBounds().height) / 2 - 10);
        window.draw(ddText);

        sf::RectangleShape astarButt(sf::Vector2f(200, 50));
        astarButt.setFillColor(sf::Color::Magenta);
        astarButt.setPosition(250, 370);
        window.draw(astarButt);

        sf::Text astarText;
        astarText.setFont(font);
        astarText.setString("A*");
        astarText.setCharacterSize(25);
        astarText.setFillColor(sf::Color::White);
        astarText.setPosition(250 + (200 - bfsText.getLocalBounds().width) / 2 + 7,
                            370 + (60 - bfsText.getLocalBounds().height) / 2 - 10);
        window.draw(astarText);

        sf::RectangleShape dualaButt(sf::Vector2f(200, 50));
        dualaButt.setFillColor(sf::Color::Magenta);
        dualaButt.setPosition(250, 440);
        window.draw(dualaButt);

        sf::Text dualastarText;
        dualastarText.setFont(font);
        dualastarText.setString("Bidirectional A*");
        dualastarText.setCharacterSize(25);
        dualastarText.setFillColor(sf::Color::White);
        dualastarText.setPosition(250 + (200 - dualastarText.getLocalBounds().width) / 2 + 7,
                              440 + (60 - dualastarText.getLocalBounds().height) / 2 - 10);
        window.draw(dualastarText);

        sf::RectangleShape ffButt(sf::Vector2f(200, 50));
        ffButt.setFillColor(sf::Color::Magenta);
        ffButt.setPosition(470, 370);
        window.draw(ffButt);

        sf::Text ffText;
        ffText.setFont(font);
        ffText.setString("Ford-Fulkerson");
        ffText.setCharacterSize(25);
        ffText.setFillColor(sf::Color::White);
        ffText.setPosition(470 + (200 - bfsText.getLocalBounds().width) / 2 -40,
                              370 + (60 - bfsText.getLocalBounds().height) / 2 - 10);
        window.draw(ffText);

        sf::RectangleShape primButt(sf::Vector2f(200, 50));
        primButt.setFillColor(sf::Color::Magenta);
        primButt.setPosition(470, 300);
        window.draw(primButt);

        sf::Text primText;
        primText.setFont(font);
        primText.setString("Prim");
        primText.setCharacterSize(25);
        primText.setFillColor(sf::Color::White);
        primText.setPosition(470 + (200 - primText.getLocalBounds().width) / 2,
                            300 + (60 - primText.getLocalBounds().height) / 2 - 10);
        window.draw(primText);

        sf::RectangleShape dfsButt(sf::Vector2f(200, 50));
        dfsButt.setFillColor(sf::Color::Magenta);
        dfsButt.setPosition(250, 300);
        window.draw(dfsButt);

        sf::Text dfsText;
        dfsText.setFont(font);
        dfsText.setString("DFS");
        dfsText.setCharacterSize(25);
        dfsText.setFillColor(sf::Color::White);
        dfsText.setPosition(250 + (200 - dfsText.getLocalBounds().width) / 2,
                            300 + (60 - dfsText.getLocalBounds().height) / 2 - 10);
        window.draw(dfsText);

        sf::RectangleShape exitButt(sf::Vector2f(200, 50));
        exitButt.setFillColor(sf::Color::Red);
        exitButt.setPosition(690, 300);
        window.draw(exitButt);

        sf::Text exitText;
        exitText.setFont(font);
        exitText.setString("Exit");
        exitText.setCharacterSize(25);
        exitText.setFillColor(sf::Color::White);
        exitText.setPosition(690 + (200 - exitText.getLocalBounds().width) / 2,
                             300 + (60 - exitText.getLocalBounds().height) / 2 - 10);
        window.draw(exitText);
    }
    if(!is_main) {
        sf::RectangleShape startButt(sf::Vector2f(200, 50));
        startButt.setFillColor(sf::Color::Green);
        startButt.setPosition(10, 10);
        window.draw(startButt);

        sf::Text startText;
        startText.setFont(font);
        startText.setString("Start");
        startText.setCharacterSize(25);
        startText.setFillColor(sf::Color::Black);
        startText.setPosition(10 + (200 - startText.getLocalBounds().width) / 2,
                              15 + (50 - startText.getLocalBounds().height) / 2 - 10);
        window.draw(startText);

        sf::RectangleShape menuButt(sf::Vector2f(200, 50));
        menuButt.setFillColor(sf::Color::Red);
        menuButt.setPosition(10, 70);
        window.draw(menuButt);

        sf::Text menuText;
        menuText.setFont(font);
        menuText.setString("Menu");
        menuText.setCharacterSize(25);
        menuText.setFillColor(sf::Color::Black);
        menuText.setPosition(10 + (200 - menuText.getLocalBounds().width) / 2,
                             75 + (50 - menuText.getLocalBounds().height) / 2 - 10);
        window.draw(menuText);
    }
}

int main() {
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;

    sf::RenderWindow window(sf::VideoMode(920, 800), "Graph Visualization", sf::Style::Default, settings);
    sf::Font font;
    if (!font.loadFromFile("C:/Users/ACER/downloads/arial_light.ttf")) {
        std::cerr << "Failed to load font" << std::endl;
        return -1;
    }

    Graph graph;
    if (!graph.read_from_file("C:/Users/ACER/downloads/graph_data_wei.txt", true)) {
        std::cerr << "Failed to read points from file or file is empty.\n";
        return -1;
    }

    int start_vertex = 1;
    int end_vertex = 53;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    if (state == 0 && event.mouseButton.y >= 300 && event.mouseButton.y <= 350) {
                        if (event.mouseButton.x >= 30 && event.mouseButton.x <= 230) {
                            state = 1;
                        } else if (event.mouseButton.x >= 250 && event.mouseButton.x <= 450) {
                            state = 2;
                        } else if (event.mouseButton.x >= 470 && event.mouseButton.x <= 670) {
                            state = 3;
                        }else if (event.mouseButton.x >= 690 && event.mouseButton.x <= 890) {
                            window.close();
                        }
                    } if(state == 0 && event.mouseButton.y >= 370 && event.mouseButton.y <= 420){
                        if (event.mouseButton.x >= 30 && event.mouseButton.x <= 230) {
                            state = 4;
                        } else if (event.mouseButton.x >= 250 && event.mouseButton.x <= 450) {
                            state = 5;
                        }else if (event.mouseButton.x >= 470 && event.mouseButton.x <= 670) {
                            state = 6;
                        }
                    } if(state == 0 && event.mouseButton.y >= 440 && event.mouseButton.y <= 490){
                        if (event.mouseButton.x >= 30 && event.mouseButton.x <= 230) {
                            state = 7;
                        } else if (event.mouseButton.x >= 250 && event.mouseButton.x <= 450) {
                            state = 8;
                        }
                    }else if ((state == 1 || state == 2 || state == 3 || state == 4 || state == 5 || state == 6 || state == 7 || state == 8) && event.mouseButton.x >= 10 && event.mouseButton.x <= 210) {
                        if (event.mouseButton.y >= 10 && event.mouseButton.y <= 60) {
                            if (state == 1){
                                graph.bfs(start_vertex, window, font); // Perform BFS
                            } else if (state == 2) {
                                graph.dfs(start_vertex, window, font); // Perform DFS
                            } else if (state == 3){
                                graph.prims_min_tree(window, font, start_vertex);
                            }else if (state == 4){
                                graph.dijkstra_path(window, font, start_vertex, end_vertex);
                            }else if (state == 5){
                                graph.astar_path(window, font, start_vertex, end_vertex);
                            }else if (state == 6){
                                graph.fordFulkerson(window, font, start_vertex, end_vertex);
                            }else if (state == 7){
                                graph.bidirect_dijkstra_path(start_vertex, end_vertex, window, font);
                            }else if (state == 8){
                                graph.bidirect_astar_path(window, font, start_vertex, end_vertex);
                            }
                            window.display();
                        }
                        if (event.mouseButton.y >= 70 && event.mouseButton.y <= 120) {
                            graph.resetVisualization();
                            state = 0;
                        }
                    }
                }
            }
        }

        window.clear(sf::Color::White);

        if (state == 0) {
            menu(window, font, true);
        } else {
            menu(window, font, false);
            graph.display(window, font);
        }

        window.display();
    }

    return 0;
}