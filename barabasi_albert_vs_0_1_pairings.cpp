#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <vector>


struct Endpoint {
    double x{0.0};
    size_t chord_id{0};
    bool is_right{false};

    Endpoint() = default;
    Endpoint(double x_, size_t chord_id_, bool is_right_) : x(x_), chord_id(chord_id_), is_right(is_right_) {}
};

struct GraphEdge {
    size_t u{0};
    size_t v{0};

    GraphEdge() = default;
    GraphEdge(size_t u_, size_t v_) : u(u_), v(v_) {}
};

struct Counter {
    size_t sum{0};
    size_t cnt{0};

    Counter() = default;
    Counter(size_t sum_, size_t cnt_) : sum(sum_), cnt(cnt_) {}
};


std::vector<GraphEdge> generate_ba_pairing_01(const size_t n, const size_t m, std::mt19937& gen) {
    std::uniform_real_distribution<> dis(0.0, 1.0);

    const size_t N = m * n;
    std::vector<Endpoint> points;
    for(size_t i = 0; i < N; i++) {
        const double r1 = dis(gen);
        const double r2 = dis(gen);

        if(r1 < r2) {
            points.emplace_back(r1, i, false);
            points.emplace_back(r2, i, true);
        }
        else {
            points.emplace_back(r2, i, false);
            points.emplace_back(r1, i, true);
        }
    }

    std::sort(points.begin(), points.end(), [](const Endpoint& a, const Endpoint& b) {
        return a.x < b.x;
    });

    std::vector<GraphEdge> edges;
    edges.reserve(N);

    std::vector<ssize_t> pair_first_vertex(N, -1);
    size_t current_vertex_id = 0;

    for(const auto& p : points) {
        if (pair_first_vertex[p.chord_id] == -1) {
            pair_first_vertex[p.chord_id] = current_vertex_id;
        } 
        else {
            const size_t u = static_cast<size_t>(pair_first_vertex[p.chord_id]) / m;
            const size_t v = current_vertex_id / m;

            if (u < v) {
                edges.emplace_back(u, v);
            }
            else {
                edges.emplace_back(v, u);
            }
        }

        if (p.is_right) {
            current_vertex_id++;
        }
    }

    return edges;
}

std::vector<GraphEdge> generate_ba_batagelj_urlika_brandes(const size_t n, const size_t m, std::mt19937& gen) {
    const size_t N = m * n;
    std::vector<size_t> points(2 * N);

    for(size_t i = 0; i < N; i++) {    
        size_t v = i / m;
        points[2 * i] = v;

        std::uniform_int_distribution<> dis(0, 2 * i);
        const size_t r = dis(gen);

        points[2 * i + 1] = points[r];
    }

    std::vector<GraphEdge> edges;
    edges.reserve(N);

    for(size_t i = 0; i < N; i++) {
        edges.emplace_back(points[2 * i], points[2 * i + 1]);
    }

    return edges;
}

std::map<size_t, size_t> generate_histogram(const std::vector<GraphEdge>& edges, const size_t n) {
    std::vector<size_t> deg(n, 0);
    for(const auto& edge: edges) {
        deg[edge.u]++;
        deg[edge.v]++;
    }

    std::map<size_t, size_t> histogram;
    for (size_t i = 0; i < n; i++) {
        histogram[deg[i]]++;
    }

    return histogram;
}

size_t count_loops(const std::vector<GraphEdge>& edges) {
    size_t loops_cnt = 0;
    for(const auto& edge: edges) {
        if(edge.u == edge.v) {
            loops_cnt++;
        }
    }

    return loops_cnt;
}

size_t count_multiple_edges(const std::vector<GraphEdge>& edges) {
    std::map<std::pair<size_t, size_t>, size_t> edge_counts;

    for (const auto& edge : edges) {
        edge_counts[{edge.u, edge.v}]++;
    }

    size_t multiple_edge_cnt = 0;
    for (const auto& kv : edge_counts) {
        if (kv.second > 1) multiple_edge_cnt += (kv.second - 1);
    }

    return multiple_edge_cnt;
}

double avg(const Counter& counter) {
    return static_cast<double>(counter.sum) / static_cast<double>(counter.cnt);
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <n> <m> <reps>\n";
        return 1;
    }

    const size_t n = std::stoi(argv[1]);
    const size_t m = std::stoi(argv[2]);
    const size_t reps = std::stoi(argv[3]);

    if (n == 0 || m == 0) {
        std::cerr << "n and m must be positive integers\n";
        return 1;
    }
    if (m > n) {
        std::cerr << "m must be <= n\n";
        return 1;
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    std::map<size_t, Counter> histogram_pairing_01;
    std::map<size_t, Counter> histogram_batagelj_urlika_brandes;
    Counter loops_count_pairing_01;
    Counter loops_count_batagelj_urlika_brandes;
    Counter multiple_edges_count_pairing_01;
    Counter multiple_edges_batagelj_urlika_brandes;

    for(size_t rep = 0; rep < reps; rep++) {
        const std::vector<GraphEdge> edges_pairing_01 = generate_ba_pairing_01(n, m, gen);
        const std::vector<GraphEdge> edges_batagelj_urlika_brandes = generate_ba_batagelj_urlika_brandes(n, m, gen);

        const std::map<size_t, size_t> hist_pairing_01 = generate_histogram(edges_pairing_01, n);
        const std::map<size_t, size_t> hist_batagelj_urlika_brandes = generate_histogram(edges_batagelj_urlika_brandes, n);

        for (const auto& kv : hist_pairing_01) {
            const size_t degree = kv.first;
            const size_t count = kv.second;
            
            // if(count == 0) {
            //     continue;
            // }

            auto& dc = histogram_pairing_01[degree];
            dc.sum += count;
            dc.cnt++;
        }

        for (const auto& kv : hist_batagelj_urlika_brandes) {
            const size_t degree = kv.first;
            const size_t count = kv.second;

            // if(count == 0) {
            //     continue;
            // }

            auto& dc = histogram_batagelj_urlika_brandes[degree];
            dc.sum += count;
            dc.cnt++;
        }

        loops_count_pairing_01.sum += count_loops(edges_pairing_01);
        loops_count_pairing_01.cnt++;

        loops_count_batagelj_urlika_brandes.sum += count_loops(edges_batagelj_urlika_brandes);
        loops_count_batagelj_urlika_brandes.cnt++;

        multiple_edges_count_pairing_01.sum += count_multiple_edges(edges_pairing_01);
        multiple_edges_count_pairing_01.cnt++;

        multiple_edges_batagelj_urlika_brandes.sum += count_multiple_edges(edges_batagelj_urlika_brandes);
        multiple_edges_batagelj_urlika_brandes.cnt++;
    }

    std::ofstream file1("result_pairing_01_" + std::to_string(m) + "_" + std::to_string(n) + "_" + std::to_string(reps) + ".txt");
    file1 << n << " " << m << std::endl;
    file1 << avg(loops_count_pairing_01) << " " << avg(multiple_edges_count_pairing_01) << std::endl;
    for(const auto& kv: histogram_pairing_01) {
        file1 << kv.first << " " << avg(kv.second) << std::endl;
    }

    std::ofstream file2("result_batagelj_urlika_brandes_" + std::to_string(m) + "_" + std::to_string(n) + "_" + std::to_string(reps) + ".txt");
    file2 << n << " " << m << std::endl;
    file2 << avg(loops_count_batagelj_urlika_brandes) << " " << avg(multiple_edges_batagelj_urlika_brandes) << std::endl;
    for(const auto& kv: histogram_batagelj_urlika_brandes) {
        file2 << kv.first << " " << avg(kv.second) << std::endl;
    }
}