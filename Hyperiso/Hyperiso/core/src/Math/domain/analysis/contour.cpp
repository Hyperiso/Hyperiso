#include "contour.h"

MarchingSquaresExtractor::MarchingSquaresExtractor(const ScalarField2D& f, std::array<double, 4> bounds, size_t max_depth) :
    f(f), max_depth(max_depth)
{
    this->xmin = bounds[0];
    this->xmax = bounds[1];
    this->ymin = bounds[2];
    this->ymax = bounds[3];
    this->N = 1 << this->max_depth;
    this->root = nullptr;
}

std::set<Path> MarchingSquaresExtractor::find_iso_contour(double lvl) {
    this->set_iso_value(lvl);
    MSGraph first_pass_graph = this->build_adjacency(this->compute_segments());
    std::set<Path> first_pass_paths = this->extract_paths(first_pass_graph);
    std::set<Path> second_pass_paths = this->correct_topology(first_pass_paths);
    return second_pass_paths;
}

std::unique_ptr<Cell> MarchingSquaresExtractor::build(int i, int j, int depth, int force_refine) {
    std::array<Index, 4> verts = this->cell_vertices({i, j}, depth);
    std::array<double, 4> vals;
    for (size_t i = 0; i < 4; i++)
        vals[i] = this->vertex_value(verts[i]);

    std::cout << "Vertex values: ";
    for (auto &&v : vals) {
        std::cout << v << ",";
    }
    std::cout << std::endl;
    

    auto positive = [] (double a) { return a > 0; };
    bool topology = !std::all_of(vals.begin(), vals.end(), positive) && std::any_of(vals.begin(), vals.end(), positive);

    std::cout << std::boolalpha << "Topology = " << topology << std::endl;

    if (!topology && depth == this->max_depth)
        return std::make_unique<Cell>(verts, vals, depth, CellStatus::EMPTY_LEAF);

    bool refine = topology || depth < force_refine;
    if (!refine && depth < this->max_depth - 1) {
        int I0 = verts[0].first, I1 = verts[1].first, J0 = verts[0].second, J1 = verts[2].second;
        int Im = (I0 + I1) / 2, Jm = (J0 + J1) / 2;
        double fc = this->vertex_value({Im, Jm});
        bool inside_loop = fc * vals[0] < 0;

        if (!inside_loop) {
            std::set<Index> edge_midpoints {{Im, J0}, {I1, Jm}, {Im, J1}, {I0, Jm}};
            auto opposite = [this, vals] (Index ijm) { return this->vertex_value(ijm) * vals[0] < 0; };
            inside_loop = std::any_of(edge_midpoints.begin(), edge_midpoints.end(), opposite);

            if (!inside_loop) {
                int Iq1 = (3 * I0 + I1) / 4, Iq2 = (I0 + 3 * I1) / 4;
                int Jq1 = (3 * J0 + J1) / 4, Jq2 = (J0 + 3 * J1) / 4;
                std::set<Index> quarter_points {{Iq1, Jq1}, {Iq2, Jq1}, {Iq2, Jq2}, {Iq1, Jq2}};
                inside_loop = std::any_of(quarter_points.begin(), quarter_points.end(), opposite);
            }
        }

        refine = refine || inside_loop;
    }

    if (!refine || depth == this->max_depth) {
        std::cout << (topology ? "Active" : "Empty") << " leaf" << std::endl;
        return std::make_unique<Cell>(verts, vals, depth, topology ? CellStatus::ACTIVE_LEAF : CellStatus::EMPTY_LEAF);
    }
        
    std::array<std::unique_ptr<Cell>, 4> children;
    std::array<Index, 4> deltas {{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
    for (size_t k = 0; k < 4; k++)
        children[k] = this->build(2 * i + deltas[k].first, 2 * j + deltas[k].second, depth + 1, force_refine);

    return std::make_unique<Cell>(verts, vals, depth, CellStatus::INTERNAL, std::move(children));
}

Index MarchingSquaresExtractor::point_to_idx(Point xy) const {
    double u = (xy.first - this->xmin) / (this->xmax - this->xmin);
    double v = (xy.second - this->ymin) / (this->ymax - this->ymin);

    int I = static_cast<int>(u * this->N);
    int J = static_cast<int>(v * this->N); 

    return {I, J};
}

Point MarchingSquaresExtractor::idx_to_point(Index ij) const {
    double x = this->xmin + (this->xmax - this->xmin) * ij.first / this->N;
    double y = this->ymin + (this->ymax - this->ymin) * ij.second / this->N;

    return {x, y};
}

std::array<Index, 4> MarchingSquaresExtractor::cell_vertices(Index ij, std::size_t depth) const {
    int scale = 1 << (this->max_depth - depth);

    int I0 = ij.first * scale;
    int I1 = (ij.first + 1) * scale;
    int J0 = ij.second * scale;
    int J1 = (ij.second + 1) * scale;

    return {{{I0, J0}, {I1, J0}, {I1, J1}, {I0, J1}}};
}

double MarchingSquaresExtractor::vertex_value(Index ij) {
    if (!this->vertices.contains(ij)) {
        Point xy = idx_to_point(ij);
        this->vertices.emplace(ij, this->f(xy.first, xy.second) - this->k);
    }

    return this->vertices.at(ij);
}

Cell* MarchingSquaresExtractor::locate_leaf(Point xy) const {
    Cell* cell = this->root.get();
    Index IJ = point_to_idx(xy);

    while (cell->children.has_value()) {
        auto verts = cell->vertices;
        int Im = (verts[0].first + verts[2].first) / 2;
        int Jm = (verts[0].second + verts[2].second) / 2;

        if (IJ.first < Im)
            cell = IJ.second < Jm ? cell->children.value()[0].get() : cell->children.value()[3].get();
        else
            cell = IJ.second < Jm ? cell->children.value()[1].get() : cell->children.value()[2].get();
    }

    return cell;
}

Point MarchingSquaresExtractor::edge_point(Cell *cell, int edge) {
    auto verts = cell->vertices;
    auto vals = cell->values;

    Index ab = Cell::EDGES[edge];
    Index IJa = verts[ab.first];
    Index IJb = verts[ab.second];

    double fa = vals[ab.first];
    double fb = vals[ab.second];
    double t = -fa / (fb - fa);

    Point xy_a = idx_to_point(IJa);
    Point xy_b = idx_to_point(IJb);

    return {
        xy_a.first + t * (xy_b.first - xy_a.first),
        xy_a.second + t * (xy_b.second - xy_a.second)
    };
}

std::set<std::pair<int, int>> MarchingSquaresExtractor::get_crossed_edges(Cell *cell) {
    int mask = 0;
    for (size_t i = 0; i < 4; i++)
        if (cell->values[i] > 0) mask |= (1 << i);
    
    if (mask == 5 || mask == 10) {
        auto verts = cell->vertices;
        int Im = (verts[0].first + verts[2].first) / 2;
        int Jm = (verts[0].second + verts[2].second) / 2;
        double fc = vertex_value({Im, Jm});
        mask = fc * cell->values[0] > 0 ? 5 : 10;
    }   

    return MASK_LOOKUP[mask];
}

bool MarchingSquaresExtractor::point_on_boundary(Point xy, double tol) const {
    double x = xy.first, y = xy.second;
    return (
        std::abs(x - this->xmin) < tol ||
        std::abs(x - this->xmax) < tol ||
        std::abs(y - this->ymin) < tol ||
        std::abs(y - this->ymax) < tol
    );
}

void MarchingSquaresExtractor::set_iso_value(double lvl) {
    this->k = lvl;
    this->vertices.clear();
    this->root = this->build(0, 0, 0);
    this->active_cells.clear();
    this->update_active_cells(this->root.get());

    std::cout << "Active cells: " << this->active_cells.size() << std::endl;
}

void MarchingSquaresExtractor::update_active_cells(Cell *cell) {
    if (cell->status == CellStatus::EMPTY_LEAF)
        return;

    if (cell->status == CellStatus::ACTIVE_LEAF) {
        this->active_cells.emplace(cell);
        return;
    }

    for (std::unique_ptr<Cell>& child : cell->children.value()) {
        this->update_active_cells(child.get());
    }
}

std::vector<Segment> MarchingSquaresExtractor::compute_segments() {
    std::vector<Segment> segments;

    for (Cell* cell : this->active_cells) {
        auto crossed_edges = this->get_crossed_edges(cell);
        for (auto edge_pair : crossed_edges) {
            Point p0 = this->edge_point(cell, edge_pair.first);
            Point p1 = this->edge_point(cell, edge_pair.second);
            segments.emplace_back(p0, p1);
        }
    }

    return segments;
}

MSGraph MarchingSquaresExtractor::build_adjacency(std::vector<Segment> segments) const {
    double coord_max = std::max({std::abs(this->xmin), std::abs(this->xmax), std::abs(this->ymin), std::abs(this->ymax)});
    double epsilon = 2 * coord_max / INT32_MAX;
    auto quantize = [epsilon] (Point p) -> std::pair<int, int> {
        return {static_cast<int>(p.first / epsilon), static_cast<int>(p.second / epsilon)};
    };

    std::map<Index, Point> vertices;
    std::map<Index, std::vector<Index>> adjacency;

    for (auto& s : segments) {
        Index k0 = quantize(s.first);
        Index k1 = quantize(s.second);
        vertices.emplace(k0, s.first);
        vertices.emplace(k1, s.second);

        if (adjacency.contains(k0))
            adjacency.at(k0).emplace_back(k1);
        else 
            adjacency.emplace(k0, std::vector {k1});

        if (adjacency.contains(k1))
            adjacency.at(k1).emplace_back(k0);
        else 
            adjacency.emplace(k1, std::vector {k0});
    }

    return MSGraph {std::move(vertices), std::move(adjacency)};
}

std::set<Path> MarchingSquaresExtractor::extract_paths(MSGraph graph) const {
    auto unused = get_keys(graph.adjacency);
    std::set<Path> paths;
    std::set<Index> path_ends;

    while (!unused.empty()) {
        auto start = *unused.begin();
        unused.erase(start);

        if (graph.adjacency.at(start).size() == 2) {
            std::stack<Index> stack;
            stack.push(start);
            std::set<Index> seen = {start};
            Index open_start;
            bool has_open_start {false};
            while (!stack.empty()) {
                Index u = stack.top();                
                stack.pop();

                if (graph.adjacency.at(u).size() == 1) {
                    open_start = u;
                    has_open_start = true;
                    break;
                }

                for (auto v : graph.adjacency.at(u)) {
                    if (!seen.contains(v)) {
                        seen.emplace(v);
                        stack.push(v);
                    }
                }
            }

            if (has_open_start)
                start = open_start;

            std::vector<Index> path = {start};
            Index prev;
            Index cur = start;
            
            if (path_ends.contains(start)) continue;

            while (true) {
                auto neigh = graph.adjacency.at(cur);
                Index next;
                bool has_next {false};
                for (Index v : neigh) {
                    if (v != prev) {
                        next = v;
                        has_next = true;
                        break;
                    }
                }
                if (!has_next) break;

                path.emplace_back(next);
                unused.erase(next);
                prev = cur;
                cur = next;

                if (cur == start) break;
            }

            path_ends.emplace(start);
            path_ends.emplace(cur);

            Path point_path;
            for (Index idx : path)
                point_path.emplace_back(graph.vertices.at(idx));
            
            paths.emplace(point_path);
        }
    }

    return paths;
}

std::set<Point> MarchingSquaresExtractor::detect_dangling_ends(std::set<Path> paths) const {
    std::set<Point> dangling;
    for (Path p : paths) {
        if (p.front() != p.back() && !(point_on_boundary(p.front()) && point_on_boundary(p.back()))) {
            dangling.emplace(p.front());
            dangling.emplace(p.back());
        }
    }
    return dangling;
}

void MarchingSquaresExtractor::refine_leaf(Cell *leaf) {
    if (leaf->status == CellStatus::INTERNAL || leaf->depth == this->max_depth)
        return;

    size_t d = leaf->depth;
    int scale = 1 << (this->max_depth - d);
    Index IJ0 = leaf->vertices[0];
    double i = IJ0.first / scale, j = IJ0.second / scale;
    std::unique_ptr<Cell> refined = this->build(i, j, d, std::min(d + 2, this->max_depth));
    *leaf = std::move(*refined);
}

std::set<Path> MarchingSquaresExtractor::correct_topology(std::set<Path> raw_paths) {
    std::set<Point> dangling = this->detect_dangling_ends(raw_paths);
    std::set<Cell*> to_refine;

    for (Point p : dangling)
        to_refine.emplace(this->locate_leaf(p));

    for (Cell* cell : to_refine)
        this->refine_leaf(cell);

    this->update_active_cells(this->root.get());
    MSGraph updated_graph = this->build_adjacency(this->compute_segments());
    return this->extract_paths(updated_graph);
}
