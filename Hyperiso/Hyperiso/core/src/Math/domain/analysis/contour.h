#ifndef __CONTOUR_H__
#define __CONTOUR_H__

#include <vector>
#include <array>
#include <set>
#include <map>
#include <stack>
#include <optional>
#include <memory>
#include "functions.h"

using Point = std::pair<double, double>;
using Segment = std::pair<Point, Point>;
using Index = std::pair<int, int>;
using Path = std::vector<Point>;

inline const std::vector<std::set<Index>> MASK_LOOKUP {
    {},
    {{3, 0}},
    {{0, 1}},
    {{3, 1}},
    {{1, 2}},
    {{3, 2}, {0, 1}},
    {{0, 2}},
    {{3, 2}},
    {{2, 3}},
    {{0, 2}},
    {{0, 3}, {1, 2}},
    {{1, 2}},
    {{1, 3}},
    {{0, 1}},
    {{3, 0}},
    {}
};

enum class CellStatus {
    INTERNAL,
    EMPTY_LEAF,
    ACTIVE_LEAF
};

struct Cell {
    static inline const std::vector<Index> EDGES {{0, 1}, {1, 2}, {2, 3}, {3, 0}}; 

    std::array<Index, 4> vertices;
    std::array<double, 4> values;
    std::size_t depth;
    CellStatus status;
    std::optional<std::array<std::unique_ptr<Cell>, 4>> children;
};

struct MSGraph {
    std::map<Index, Point> vertices;
    std::map<Index, std::vector<Index>> adjacency;
};

class MarchingSquaresExtractor {
public:
    MarchingSquaresExtractor(RealValuedForm f, std::array<double, 4> bounds, size_t max_depth = 7);

    std::set<Path> find_iso_contour(double lvl);

private:
    std::unique_ptr<Cell> build(int i, int j, int depth, int force_refine = 2);

    Index point_to_idx(Point xy) const;
    Point idx_to_point(Index ij) const;
    std::array<Index, 4> cell_vertices(Index ij, std::size_t depth) const;
    double vertex_value(Index ij);
    Cell* locate_leaf(Point xy) const;
    Point edge_point(Cell* cell, int edge);
    std::set<std::pair<int, int>> get_crossed_edges(Cell* cell);
    bool point_on_boundary(Point xy, double tol = 1e-8) const;

    void set_iso_value(double lvl);
    void update_active_cells(Cell* cell);
    std::vector<Segment> compute_segments();
    MSGraph build_adjacency(std::vector<Segment> segments) const;
    std::set<Path> extract_paths(MSGraph graph) const;
    std::set<Point> detect_dangling_ends(std::set<Path> paths) const;
    void refine_leaf(Cell *leaf);
    std::set<Path> correct_topology(std::set<Path> raw_paths);

private:
    double xmin, xmax, ymin, ymax;
    std::size_t max_depth;
    int N;
    std::map<Index, double> vertices;
    std::set<Cell*> active_cells;
    std::unique_ptr<Cell> root;
    RealValuedForm f;
    double k;
};

namespace std {
    template <>
    struct hash<std::pair<int, int>> {
        std::size_t operator()(const Index& p) const noexcept {
            std::size_t h = static_cast<std::size_t>(p.first);
            h ^= (static_cast<std::size_t>(p.second) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
            return h;
        }
    };
}

#endif // __CONTOUR_H__
