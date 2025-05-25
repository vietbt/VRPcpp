#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cfloat>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>
#define CHAR_LEN 100
#define EPS 1e-9
#define INF 0x3f3f3f3f
using namespace std;
const int OUTLIER = -2;
const int NOTCLASSIFIED = -1;
#define TERMINATION 25000 * ACTUAL_PROBLEM_SIZE
#define STOP_CNT 25000 * ACTUAL_PROBLEM_SIZE - 5
struct floydWarshall {
    int nV;
    vector<vector<int>> dist;
    vector<vector<int>> next;
};
struct node {
    int id;
    double x;
    double y;
};
struct solution {
    int *tour;
    int id;
    double tour_length;
    int steps;
};
struct point {
    int id;
    int neighbours_count = 0;
    int cluster_ID;
    vector<int> adjacentPoints;
};
class EVRPEnv {
    public:
        struct node *node_list;
        int *cust_demand;
        bool *charging_station;
        double **distances;
        int problem_size;
        double energy_consumption;
        int DEPOT;
        int NUM_OF_CUSTOMERS;
        int ACTUAL_PROBLEM_SIZE;
        double OPTIMUM;
        int NUM_OF_STATIONS;
        int BATTERY_CAPACITY;
        int MAX_CAPACITY;
        vector<int> CUSTOMERS;
        vector<int> AFSs;
        map<int, int> afsIdMap;
        floydWarshall* fw;
        double evals;
        double current_best;
        solution *best_sol;
        double euclidean_distance(int i, int j) {
            double xd = node_list[i].x - node_list[j].x;
            double yd = node_list[i].y - node_list[j].y;
            double r = sqrt(xd * xd + yd * yd);
            return r;
        }
        void compute_distances(bool round_int) {
            for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++)
                for (int j = 0; j < ACTUAL_PROBLEM_SIZE; j++) {
                    double d = euclidean_distance(i, j);
                    if (round_int) {
                        d = (double) (int) (d+0.5);
                    }
                    distances[i][j] = d;
                }
                    
        }
        double **generate_2D_matrix_double(int n, int m) {
            double **matrix;
            matrix = new double *[n];
            for (int i = 0; i < n; i++) matrix[i] = new double[m];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                    matrix[i][j] = 0.0;
            return matrix;
        }
        void init_evals() { evals = 0; }
        void init_current_best() { current_best = INT_MAX; }
        void initialize_heuristic() {
            best_sol = new solution;
            best_sol->tour = new int[4 * NUM_OF_CUSTOMERS];
            best_sol->id = 1;
            best_sol->steps = 0;
            best_sol->tour_length = INT_MAX;
        }
        bool is_charging_station(int node) {
            bool flag = false;
            if (charging_station[node] == true) flag = true; 
            else flag = false;
            return flag;
        }
        double get_energy_consumption(int from, int to) {
            return energy_consumption * distances[from][to];
        }
        void initFloydWarshall(int nV) {
            fw = new floydWarshall;
            fw->nV = nV;
            fw->next = vector<vector<int>>(nV, vector<int>(nV, -1));
            fw->dist = vector<vector<int>>(nV, vector<int>(nV, INF));
            for (int i = 0; i < nV; i++) {
                for (int j = i; j < nV; j++) {
                    if (i == j) {
                        fw->dist[i][j] = 0;
                        fw->next[i][j] = j;
                    } else {
                        int start = AFSs[i];
                        int goal = AFSs[j];
                        double consumption = get_energy_consumption(start, goal);
                        if (consumption <= BATTERY_CAPACITY) {
                            fw->dist[i][j] = consumption;
                            fw->dist[j][i] = consumption;
                            fw->next[i][j] = j;
                            fw->next[j][i] = i;
                        }
                    }
                }
            }
        }
        void floydWarshall2(vector<vector<int>> &graph) {
            fw->next = vector<vector<int>>(fw->nV, vector<int>(fw->nV, -1));
            fw->nV = graph.size();
            for (int i = 0; i < fw->nV; i++)
                for (int j = 0; j < fw->nV; j++)
                    fw->next[i][j] = j;
        }
        void planPaths() {
            int i, j, k;
            for (k = 0; k < fw->nV; k++)
                for (i = 0; i < fw->nV; i++)
                    for (j = 0; j < fw->nV; j++)
                        if (fw->dist[i][k] + fw->dist[k][j] < fw->dist[i][j]) {
                            fw->dist[i][j] = fw->dist[i][k] + fw->dist[k][j];
                            fw->next[i][j] = fw->next[i][k];
                        }
        }
        vector<int> getPath(int u, int v, bool afsIds) {
            vector<int> path;
            if (fw->next[u][v] == -1) return path;
            path.push_back(u);
            while (u != v) {
                u = fw->next[u][v];
                path.push_back(u);
            }
            if (afsIds) for (auto &n : path) n = AFSs[n];
            return path;
        }
        void initMyStructures() {
            AFSs.clear();
            CUSTOMERS.clear();
            int afsId = 0;
            for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++) {
                if (is_charging_station(i)) {
                    AFSs.push_back(i);
                    afsIdMap.insert(pair<int, int>(i, afsId++));
                } else CUSTOMERS.push_back(i);
            };
            initFloydWarshall(AFSs.size());
            planPaths();
        }
        double get_energy_per_unit() { return energy_consumption; }
        double get_distance(int from, int to) {
            evals += (1.0 / ACTUAL_PROBLEM_SIZE);
            return distances[from][to];
        }
        void density_connected(int current_p, int cluster, int min_pt, vector<point> &points) {
            points.at(current_p).cluster_ID = cluster;
            if (points.at(current_p).neighbours_count < min_pt) return;
            for (auto &next : points.at(current_p).adjacentPoints) {
                if (points.at(next).cluster_ID != NOTCLASSIFIED && points.at(next).cluster_ID != OUTLIER)
                    continue;
                density_connected(next, cluster, min_pt, points);
            }
        }
        vector<vector<vector<int>>> dbca() {
            double reach = BATTERY_CAPACITY / get_energy_per_unit();
            vector<double> epss;
            vector<int> min_pts;
            epss.push_back(reach / (2));
            epss.push_back(reach / (3));
            for (int i = 2; i < 6; i++) epss.push_back(reach / (i * 2));
            epss.push_back(reach / (15));
            epss.push_back(reach / (20));
            for (int i = 2; i < 6; i++) min_pts.push_back(i);
            vector<vector<vector<int>>> cluster_sets;
            for (int x = 0; x < epss.size(); x++) {
                double eps = epss.at(x);
                point* points_base = new point[ACTUAL_PROBLEM_SIZE];
                for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++) {
                    points_base[i].cluster_ID = NOTCLASSIFIED;
                    for (int j = 0; j < ACTUAL_PROBLEM_SIZE; j++) {
                        if (i == j) continue;
                        if (get_distance(i, j) <= eps) {
                            points_base[i].neighbours_count++;
                            points_base[i].adjacentPoints.push_back(j);
                        }
                    }
                }
                for (int y = 0; y < min_pts.size(); y++) {
                    int min_pt = min_pts.at(y);
                    int cluster_idx = -1;
                    vector<point> points;
                    vector<vector<int>> cluster_set;
                    for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++) points.push_back(points_base[i]);
                    for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++) {
                        if (points.at(i).cluster_ID != NOTCLASSIFIED) continue;
                        if (points.at(i).neighbours_count >= min_pt) {
                            cluster_idx++;
                            density_connected(i, cluster_idx, min_pt, points);
                        } else points.at(i).cluster_ID = OUTLIER;
                    }
                    for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++) {
                        if (points.at(i).cluster_ID != OUTLIER) continue;
                        double minDist = DBL_MAX;
                        int min_node_id = -2;
                        for (int j = 0; j < ACTUAL_PROBLEM_SIZE; j++) {
                            if (i == j) continue;
                            if (points.at(j).cluster_ID == OUTLIER) continue;
                            double dist = get_distance(i, j);
                            if (dist < minDist) {
                                minDist = dist;
                                min_node_id = j;
                            }
                        }
                        if (min_node_id == -2) {
                            cluster_idx++;
                            points.at(i).cluster_ID = cluster_idx;
                        } else points.at(i).cluster_ID = points.at(min_node_id).cluster_ID;
                    }
                    cluster_set.resize(cluster_idx + 1);
                    for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++)
                        if (points.at(i).cluster_ID != OUTLIER)
                            cluster_set[points.at(i).cluster_ID].push_back(i);
                    for (int i = 0; i < cluster_set.size(); i++) {
                        bool has_customer = false;
                        for (auto &next : cluster_set.at(i))
                            if (!is_charging_station(next)) has_customer = true;
                        if (!has_customer) {
                            cluster_set.erase(cluster_set.begin() + i);
                        }
                    }
                    for (int i = 0; i < cluster_set.size(); i++) {
                        bool has_depo = false;
                        for (auto &next : cluster_set.at(i)) if (next == DEPOT) has_depo = true;
                        if (!has_depo) cluster_set.at(i).push_back(DEPOT);
                    }
                    bool has_duplicate = false;
                    for (auto &c : cluster_sets) if (c == cluster_set) has_duplicate = true;
                    if (has_duplicate) continue;
                    cluster_sets.push_back(cluster_set);
                }
            }
            return cluster_sets;
        }
        int get_customer_demand(int customer) { return cust_demand[customer]; }
        vector<int> clarke_wright(bool capacitated, bool clusters, vector<int> node_list) {
            set<int> unusedCustomers;
            if (clusters) {
                for (int node : node_list)
                    if (!is_charging_station(node))
                        unusedCustomers.insert(node);
            } else {
                for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++)
                    if (!is_charging_station(i))
                        unusedCustomers.insert(i);
            }
            vector<vector<int>> subtours;
            vector<int> subtours_length;
            while (unusedCustomers.size() > 0) {
                vector<int> subtour;
                int remaining_capacity = MAX_CAPACITY;
                double length = 0;
                auto maxDist = 0;
                int furthest;
                for (int cand : unusedCustomers) {
                    double dist = get_distance(0, cand);
                    if (dist > maxDist) {
                        maxDist = dist;
                        furthest = cand;
                    }
                }
                subtour.push_back(furthest);
                remaining_capacity -= get_customer_demand(furthest);
                unusedCustomers.erase(furthest);
                double dist_front = maxDist;
                double dist_back = maxDist;
                bool enough_capacity = true;
                while (enough_capacity) {
                    enough_capacity = false;
                    auto distImprovement = DBL_MAX;
                    int closest;
                    int neigbour;
                    int front = subtour.front();
                    int back = subtour.back();
                    bool at_front = false;
                    double dist_to_depo;
                    for (int cand : unusedCustomers) {
                        if ((get_customer_demand(cand) <= remaining_capacity) || !capacitated) {
                            enough_capacity = true;
                            double dist_candidate_depo = get_distance(0, cand);
                            double dist_saved_front = get_distance(front, cand) - dist_candidate_depo - dist_front;
                            double dist_saved_back = get_distance(back, cand) - dist_candidate_depo - dist_back;
                            if (dist_saved_front < distImprovement) {
                                at_front = true;
                                distImprovement = dist_saved_front;
                                closest = cand;
                                dist_to_depo = dist_candidate_depo;
                            }
                            if (dist_saved_back < distImprovement) {
                                at_front = false;
                                distImprovement = dist_saved_back;
                                closest = cand;
                                dist_to_depo = dist_candidate_depo;
                            }
                        }
                    }
                    if (!enough_capacity) break;
                    if (at_front) {
                        dist_front = dist_to_depo;
                        subtour.insert(subtour.begin(), closest);
                    } else {
                        dist_back = dist_to_depo;
                        subtour.push_back(closest);
                    }
                    remaining_capacity -= get_customer_demand(closest);
                    unusedCustomers.erase(closest);
                    if (unusedCustomers.size() == 0) {
                        enough_capacity = false;
                        break;
                    }
                }
                subtours.push_back(subtour);
            }
            vector<int> tmp;
            if (capacitated) {
                tmp.push_back(0);
                for (auto tour : subtours) {
                    for (int customer : tour)
                        tmp.push_back(customer);
                    tmp.push_back(0);
                }
            } else tmp = subtours.at(0);
            return tmp;
        }
        int getRemainingLoad(vector<int> evrpTour) {
            int load = 0;
            for (auto node : evrpTour) {
                if (node == 0) load = MAX_CAPACITY;
                else load -= get_customer_demand(node);
            }
            return load;
        }
        double getRemainingBattery(vector<int> evrpTour) {
            double battery = 0;
            for (int i = 0; i < evrpTour.size(); i++) {
                int cur = evrpTour[i];
                if (i > 0) {
                    int prev = evrpTour[i - 1];
                    battery -= get_energy_consumption(prev, cur);
                }
                if (is_charging_station(cur)) battery = BATTERY_CAPACITY;
            }
            return battery;
        }
        int getClosestAFS(int node) {
            auto minDist = DBL_MAX;
            int closest;
            for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++) {
                if (is_charging_station(i) && i != node) {
                    double dist = get_distance(node, i);
                    if (dist < minDist) {
                        minDist = dist;
                        closest = i;
                    }
                }
            }
            return closest;
        }
        bool addAndCheckLastN(int node, bool reset = false) {
            static int lastN[4];
            static int nextLastNIndex = 0;
            static int counter = 0;
            if (reset) {
                for (int i = 0; i < 4; i++) lastN[i] = -1;
                counter = 0;
            }
            counter++;
            bool check = true;
            if (counter > 4) {
                check = false;
                for (int i = 0; i < 2; i++) {
                    int index1 = i % 4;
                    int index2 = (i + 2) % 4;
                    if (lastN[index1] != lastN[index2]) {
                        check = true;
                        break;
                    }
                }
            }
            if (check) {
                lastN[nextLastNIndex] = node;
                nextLastNIndex = (nextLastNIndex + 1) % 4;
                return true;
            } else {
                return false;
            }
        }
        bool get_to_depot_possibly_through_afss(vector<int> &evrpTour) {
            bool canReachDepot = getRemainingBattery(evrpTour) - get_energy_consumption(evrpTour.back(), 0) >= 0;
            if (canReachDepot) evrpTour.push_back(0);
            else {
                int closestAFS = getClosestAFS(evrpTour.back());
                vector<int> afsSubpath = getPath(afsIdMap[closestAFS], afsIdMap[0], true);
                evrpTour.insert(evrpTour.end(), afsSubpath.begin(), afsSubpath.end());
            }
            return addAndCheckLastN(0);
        }
        int getReachableAFSClosestToGoal(int cur, int goal, double battery) {
            auto minDist = DBL_MAX;
            int closest = -1;
            for (int i = 0; i < ACTUAL_PROBLEM_SIZE; i++) {
                if (is_charging_station(i) && i != cur && battery >= get_energy_consumption(cur, i)) {
                    double dist = get_distance(i, goal);
                    if (dist < minDist) {
                        minDist = dist;
                        closest = i;
                    }
                }
            }
            return closest;
        }
        vector<int> tsp2evrp_zga_relaxed(vector<int> tspTour) {
            vector<int> evrpTour;
            int nextId = 0;
            evrpTour.push_back(0);
            while (nextId != tspTour.size()) {
                if (getRemainingLoad(evrpTour) < get_customer_demand(tspTour[nextId])) {
                    bool check = get_to_depot_possibly_through_afss(evrpTour);
                    if (!check) break;
                } else {
                    int closestAFSToGoal = getClosestAFS(tspTour[nextId]);
                    double remainingBattery = getRemainingBattery(evrpTour);
                    double energyToNext = get_energy_consumption(evrpTour.back(), tspTour[nextId]);
                    double nextToAFS = get_energy_consumption(tspTour[nextId], closestAFSToGoal);
                    if (remainingBattery - energyToNext >= nextToAFS) {
                        evrpTour.push_back(tspTour[nextId]);
                        nextId++;
                        if (!addAndCheckLastN(nextId - 1)) break;
                    } else {
                        int closestAFS = getReachableAFSClosestToGoal(evrpTour.back(), tspTour[nextId], getRemainingBattery(evrpTour));
                        bool canReach = getRemainingBattery(evrpTour) - get_energy_consumption(evrpTour.back(), closestAFS) >= 0;
                        assert(canReach);
                        evrpTour.push_back(closestAFS);
                        if (!addAndCheckLastN(closestAFS)) break;
                    }
                }
            }
            get_to_depot_possibly_through_afss(evrpTour);
            return evrpTour;
        }
        double fitness_evaluation(vector<int> &tour) {
            double tour_length = 0;
            for (int i = 0; i < tour.size() - 1; i++)
                tour_length += distances[tour[i]][tour[i + 1]];
            if (tour_length < current_best) current_best = tour_length;
            evals++;
            return tour_length;
        }
        vector<int> init_from_dbca() {
            auto cluster_sets = dbca();
            vector<vector<int>> evrp_tours;
            for (int x = 0; x < cluster_sets.size(); x++) {
                auto cluster_set = cluster_sets.at(x);
                vector<int> evrp_tour;
                for (int y = 0; y < cluster_set.size(); y++) {
                    auto cluster = cluster_set.at(y);
                    vector<int> initTour;
                    vector<int> tmp_evrp_tour;
                    initTour = clarke_wright(true, true, cluster);
                    tmp_evrp_tour = tsp2evrp_zga_relaxed(initTour);
                    evrp_tour.insert(evrp_tour.end(), tmp_evrp_tour.begin(), tmp_evrp_tour.end());
                }
                int last = -1;
                for (int i = 0; i < evrp_tour.size(); i++) {
                    if (last == 0 && evrp_tour.at(i) == 0) {
                        evrp_tour.erase(evrp_tour.begin() + i);
                        i--;
                    }
                    last = evrp_tour.at(i);
                }
                evrp_tours.push_back(evrp_tour);
            }
            double minEval = DBL_MAX;
            vector<int> bestTour;
            for (auto &tour : evrp_tours) {
                double eval = fitness_evaluation(tour);
                if (eval < minEval) {
                    minEval = eval;
                    bestTour = tour;
                }
            }
            return bestTour;
        }
        double get_evals() { return evals; }
        double triangle_adition(int afs, int b, int c) {
            return get_distance(b, afs) + get_distance(afs, c) - get_distance(b, c);
        }
        vector<int> tsp2evrp_zga_mini(vector<int> tspTour) {
            vector<int> evrpTour;
            int nextId = 0;
            addAndCheckLastN(-1, true);
            evrpTour.push_back(0);
            while (nextId != tspTour.size()) {
                int closestAFSToGoal = getClosestAFS(tspTour[nextId]);
                double remainingBattery = getRemainingBattery(evrpTour);
                double energyToNext = get_energy_consumption(evrpTour.back(), tspTour[nextId]);
                double nextToAFS = get_energy_consumption(tspTour[nextId], closestAFSToGoal);
                if (remainingBattery - energyToNext >= nextToAFS) {
                    evrpTour.push_back(tspTour[nextId]);
                    nextId++;
                    if (!addAndCheckLastN(tspTour[nextId - 1]))
                        break;
                } else {
                    int closestAFS = getReachableAFSClosestToGoal(evrpTour.back(), tspTour[nextId], getRemainingBattery(evrpTour));
                    bool canReach = getRemainingBattery(evrpTour) - get_energy_consumption(evrpTour.back(), closestAFS) >= 0;
                    assert(canReach);
                    evrpTour.push_back(closestAFS);
                    if (!addAndCheckLastN(closestAFS)) break;
                }
            }
            get_to_depot_possibly_through_afss(evrpTour);
            return evrpTour;
        }
        bool lost_customers(vector<int> test_tour, int orig_number) {
            int new_number = 0;
            for (int i = 0; i < test_tour.size(); i++) {
                int node = test_tour.at(i);
                if (!is_charging_station(node)) new_number++;
            }
            if (orig_number != new_number) return true;
            return false;
        }
        double triangle_adition2(int afs, int b, int c, double dist_bc) {
            return get_distance(b, afs) + get_distance(afs, c) - dist_bc;
        }
        double AFSrealoc_one_afs(vector<int> &tour, int original_afs_at) {
            int afs = tour.at(original_afs_at);
            int b = tour.at(original_afs_at - 1);
            int c = tour.at(original_afs_at + 1);
            double additional_dist_original = triangle_adition(afs, b, c);
            vector<int> r_tour;
            int n_of_cust_orig = 0;
            for (int i = (tour.size() - 1); i >= 0; i--) {
                int node = tour.at(i);
                if (!is_charging_station(node)) {
                    r_tour.push_back(node);
                    n_of_cust_orig++;
                }
            }
            vector<int> f_tour = r_tour;
            reverse(f_tour.begin(), f_tour.end());
            f_tour = tsp2evrp_zga_mini(f_tour);
            if (lost_customers(f_tour, n_of_cust_orig)) {
            }
            int front_afs_at = 0;
            int n_of_front_afs = 0;
            for (int i = 1; i < f_tour.size() - 1; i++) {
                int node = f_tour.at(i);
                if (is_charging_station(node)) {
                    n_of_front_afs++;
                    front_afs_at = i;
                }
            }
            if (n_of_front_afs > 1) return 0;
            else if (n_of_front_afs < 1) {
                double improvement = additional_dist_original;
                if (improvement > EPS) tour = f_tour;
                return improvement;
            }
            r_tour = tsp2evrp_zga_mini(r_tour);
            if (lost_customers(r_tour, n_of_cust_orig)){
            }
            int back_afs_at = 0;
            int n_of_back_afs = 0;
            for (int i = 1; i < r_tour.size() - 1; i++) {
                int node = r_tour.at(i);
                if (is_charging_station(node)) {
                    n_of_back_afs++;
                    back_afs_at = i;
                }
            }
            if (n_of_back_afs > 1) return 0;
            else if (n_of_back_afs < 1) {
                double improvement = additional_dist_original;
                if (improvement > EPS) tour = r_tour;
                return improvement;
            }
            afs = f_tour.at(front_afs_at);
            b = f_tour.at(front_afs_at - 1);
            c = f_tour.at(front_afs_at + 1);
            double additional_dist_front_afs = triangle_adition(afs, b, c);
            afs = r_tour.at(back_afs_at);
            b = r_tour.at(back_afs_at - 1);
            c = r_tour.at(back_afs_at + 1);
            double additional_dist_back_afs = triangle_adition(afs, b, c);
            int back_pos_in_orig = tour.size() - back_afs_at - 1;
            if (front_afs_at == back_pos_in_orig && r_tour.at(back_afs_at) == tour.at(front_afs_at))
                return 0;
            else if (abs(front_afs_at - back_pos_in_orig) == 1) {
                if (additional_dist_original <= additional_dist_front_afs && additional_dist_original <= additional_dist_back_afs)
                    return 0;
                else if (additional_dist_front_afs <= additional_dist_back_afs) {
                    double improvement = additional_dist_original - additional_dist_front_afs;
                    if (improvement > EPS) tour = f_tour;
                    return improvement;
                } else {
                    reverse(r_tour.begin(), r_tour.end());
                    double improvement = additional_dist_original - additional_dist_back_afs;
                    if (improvement > EPS) tour = r_tour;
                    return improvement;
                }
            }
            int front = back_pos_in_orig;
            int back = front_afs_at;
            double minAdd = 0;
            bool front_smaller = false;
            bool back_smaller = false;
            bool original_smaller = false;
            double improvement = -1;
            if (additional_dist_original <= additional_dist_front_afs && additional_dist_original <= additional_dist_back_afs) {
                minAdd = additional_dist_original;
                original_smaller = true;
            } else if (additional_dist_front_afs <= additional_dist_back_afs) {
                front_smaller = true;
                minAdd = additional_dist_front_afs;
            } else {
                back_smaller = true;
                minAdd = additional_dist_back_afs;
            }
            vector<int> clean_tour;
            for (int i = 0; i < f_tour.size(); i++) {
                int node = f_tour.at(i);
                if (i != back) clean_tour.push_back(node);
            }
            int best_index = -1;
            int best_afs = -1;
            for (int i = front; i < back - 1; i++) {
                double bi = clean_tour.at(i);
                double ci = clean_tour.at(i + 1);
                double dist_bc = get_distance(bi, ci);
                for (int j = 0; j < ACTUAL_PROBLEM_SIZE; j++) {
                    if (is_charging_station(j)) {
                        int p_afs = j;
                        double addi = triangle_adition2(p_afs, bi, ci, dist_bc);
                        if (addi < minAdd) {
                            minAdd = addi;
                            best_index = i;
                            best_afs = p_afs;
                        }
                    }
                }
            }
            improvement = additional_dist_original - minAdd;
            if (improvement <= EPS) return improvement;
            if (best_index == -1) {
                if (original_smaller) return 0;
                if (front_smaller) {
                    improvement = additional_dist_original - additional_dist_front_afs;
                    if (improvement > EPS) tour = f_tour;
                    return improvement;
                }
                if (back_smaller) {
                    improvement = additional_dist_original - additional_dist_back_afs;
                    if (improvement > EPS) tour = r_tour;
                    return improvement;
                }
            }
            clean_tour.insert(clean_tour.begin() + best_index + 1, best_afs);
            if (improvement > EPS) {
                if (lost_customers(clean_tour, n_of_cust_orig)){
                }
                tour = clean_tour;
            }
            return improvement;
        }
        double AFSrealoc_more(vector<int> &tour) {
            double tour_dist = 0;
            double f_tour_dist = 0;
            double r_tour_dist = 0;
            vector<int> r_tour;
            int num_c_orig = 0;
            for (int i = (tour.size() - 1); i >= 0; i--) {
                int node = tour.at(i);
                if (!is_charging_station(node)) {
                    r_tour.push_back(node);
                    num_c_orig++;
                }
            }
            vector<int> f_tour = r_tour;
            reverse(f_tour.begin(), f_tour.end());
            f_tour = tsp2evrp_zga_mini(f_tour);
            r_tour = tsp2evrp_zga_mini(r_tour);
            for (int i = 0; i < tour.size() - 1; i++) {
                int n1 = tour.at(i);
                int n2 = tour.at(i + 1);
                tour_dist += get_distance(n1, n2);
            }
            for (int i = 0; i < f_tour.size() - 1; i++) {
                int n1 = f_tour.at(i);
                int n2 = f_tour.at(i + 1);
                f_tour_dist += get_distance(n1, n2);
            }
            for (int i = 0; i < r_tour.size() - 1; i++) {
                int n1 = r_tour.at(i);
                int n2 = r_tour.at(i + 1);
                r_tour_dist += get_distance(n1, n2);
            }
            int num_c_new = 0;
            if (tour_dist <= f_tour_dist && tour_dist <= r_tour_dist) return 0;
            if (f_tour_dist <= r_tour_dist) {
                for (int i = 0; i < f_tour.size(); i++) {
                    int node = f_tour.at(i);
                    if (!is_charging_station(node)) num_c_new++;
                }
                if (num_c_orig != num_c_new) return 0;
                double improvement = tour_dist - f_tour_dist;
                if (improvement > EPS) tour = f_tour;
                return improvement;
            } else {
                for (int i = 0; i < r_tour.size(); i++) {
                    int node = r_tour.at(i);
                    if (!is_charging_station(node)) num_c_new++;
                }
                if (num_c_orig != num_c_new) return 0;
                reverse(r_tour.begin(), r_tour.end());
                double improvement = tour_dist - r_tour_dist;
                if (improvement > EPS) tour = r_tour;
                return improvement;
            }
        }
        bool AFSrealoc_common(vector<int> &tour, bool firstImprove, bool more_than_one) {
            double evals_start = get_evals();
            bool eval_ok = true;
            bool did_it_improve = false;
            vector<vector<int>> subtours;
            vector<int> sub_tour;
            bool first = true;
            for (int i = 0; i < tour.size(); i++) {
                int node = tour.at(i);
                if (node == 0 && !first) {
                    sub_tour.push_back(node);
                    subtours.push_back(sub_tour);
                    sub_tour.clear();
                }
                sub_tour.push_back(node);
                first = false;
            }
            vector<int> sub_tour_afs_n;
            vector<int> sub_tour_afs_at;
            int n_of_one_AFS_subtours = 0;
            for (int j = 0; j < subtours.size(); j++) {
                auto st = subtours.at(j);
                int n_of_afs = 0;
                int afs_at = 0;
                for (int i = 1; i < st.size() - 1; i++) {
                    int node = st.at(i);
                    if (is_charging_station(node)) {
                        n_of_afs++;
                        afs_at = i;
                    }
                }
                sub_tour_afs_n.push_back(n_of_afs);
                sub_tour_afs_at.push_back(afs_at);
                if (n_of_afs == 1) n_of_one_AFS_subtours++;
            }
            bool first_update = false;
            double total_improvement = 0;
            vector<vector<int>> subtours_new;
            for (int j = 0; j < subtours.size(); j++) {
                if (get_evals() > TERMINATION - 25) eval_ok = false;
                auto st = subtours.at(j);
                auto st_orig = subtours.at(j);
                int n_of_afs = sub_tour_afs_n.at(j);
                double improvement = 0;
                if (n_of_afs == 1 && !more_than_one && eval_ok) {
                    if (!first_update || !firstImprove) {
                        int afs_at = sub_tour_afs_at.at(j);
                        improvement = AFSrealoc_one_afs(st, afs_at);
                        if (improvement > EPS) first_update = true;
                    }
                } else if (n_of_afs > 1 && more_than_one && eval_ok) {
                    if (!first_update || !firstImprove) {
                        int afs_at = sub_tour_afs_at.at(j);
                        improvement = AFSrealoc_more(st);
                        if (improvement > EPS) first_update = true;
                    }
                }
                if (improvement > EPS) {
                    total_improvement += improvement;
                    subtours_new.push_back(st);
                } else subtours_new.push_back(st_orig);
            }
            vector<int> evrp_tour;
            for (auto &st : subtours_new)
                evrp_tour.insert(evrp_tour.end(), st.begin(), st.end());
            int last = -1;
            for (int i = 0; i < evrp_tour.size(); i++) {
                if (last == 0 && evrp_tour.at(i) == 0) {
                    evrp_tour.erase(evrp_tour.begin() + i);
                    i--;
                }
                last = evrp_tour.at(i);
            }
            if (total_improvement > EPS) {
                tour = evrp_tour;
                did_it_improve = true;
            }
            return did_it_improve;
        }
        bool AFSrealoc_one(vector<int> &tour, bool firstImprove) {
            return AFSrealoc_common(tour, firstImprove, false);
        }
        bool AFSrealoc_more_than_one(vector<int> &tour, bool firstImprove) {
            return AFSrealoc_common(tour, firstImprove, true);
        }
        double twoOptCostUpdate(vector<int> &tour, int i, int j) {
            double cut1 = get_distance(tour[i - 1], tour[i]);
            double cut2 = get_distance(tour[j], tour[j + 1]);
            double removed = cut1 + cut2;
            double add1 = get_distance(tour[i - 1], tour[j]);
            double add2 = get_distance(tour[i], tour[j + 1]);
            double added = add1 + add2;
            return removed - added;
        }
        void twoOptMove(vector<int> &tour, int i, int j) {
            reverse(tour.begin() + i, tour.begin() + j + 1);
        }
        bool isValidTour(vector<int> tour) {
            int load = 0;
            double battery = 0;
            int custCnt = 0;
            for (int i = 0; i < tour.size(); i++) {
                int cur = tour[i];
                if (cur == 0) load = MAX_CAPACITY;
                else {
                    load -= get_customer_demand(cur);
                    if (load < 0) return false;
                }
                if (i > 0) {
                    int prev = tour[i - 1];
                    battery -= get_energy_consumption(prev, cur);
                    if (battery < 0) return false;
                }
                if (is_charging_station(cur)) battery = BATTERY_CAPACITY;
                else custCnt++;
            }
            if ((tour[0] != 0) || (tour[tour.size() - 1] != 0)) return false;
            if (custCnt < NUM_OF_CUSTOMERS) return false;
            return true;
        }
        bool twoOpt(vector<int> &tour, bool firstImprove) {
            vector<int> tmpTour;
            vector<int> bestTour;
            tmpTour.reserve(tour.size());
            bestTour.reserve(tour.size());
            bool improved = false;
            double impBest = EPS;
            for (int i = 1; i < tour.size() - 1; i++) {
                if (get_evals() > STOP_CNT || (firstImprove && improved)) break;
                for (int j = i + 1; j < tour.size() - 1; j++) {
                    if (get_evals() > STOP_CNT || (firstImprove && improved)) break;
                    auto impCurr = twoOptCostUpdate(tour, i, j);
                    if (impCurr > impBest) {
                        tmpTour = tour;
                        twoOptMove(tmpTour, i, j);
                        if (isValidTour(tmpTour)) {
                            impBest = impCurr;
                            bestTour = tmpTour;
                            improved = true;
                        }
                    }
                }
            }
            if (improved) tour = bestTour;
            return improved;
        }
        double twoStringCostUpdate(vector<int> &tour, int i, int j, int X, int Y) {
            auto last = tour.size() - 1;
            double cut1 = get_distance(tour[i], tour[i + 1]);
            double cut2 = X != 0 ? get_distance(tour[i + X], tour[i + X + 1]) : 0;
            double cut3 = (j != last && j != i + X) ? get_distance(tour[j], tour[j + 1]) : 0;
            double cut4 = (Y != 0 && j + Y != last) ? get_distance(tour[j + Y], tour[j + Y + 1]) : 0;
            double removed = cut1 + cut2 + cut3 + cut4;
            double add1 = Y != 0 ? get_distance(tour[i], tour[j + 1]) : get_distance(tour[i], tour[i + X + 1]);
            double add2 = 0;
            if (Y != 0)
                add2 = j != i + X ? get_distance(tour[j + Y], tour[i + X + 1]) : get_distance(tour[j + Y], tour[i + 1]);
            double add3 = 0;
            if (j != i + X) {
                if (X != 0) add3 = get_distance(tour[j], tour[i + 1]);
                else if (j != last && j + Y != last) add3 = get_distance(tour[j], tour[j + Y + 1]);
            }
            double add4 = (X != 0 && j + Y != last) ? get_distance(tour[i + X], tour[j + Y + 1]) : 0;
            double added = add1 + add2 + add3 + add4;
            return removed - added;
        }
        bool twoStringMove(vector<int> &tour, int i, int j, int X, int Y) {
            if ((X == 0 || Y == 0) && i == j) return false;
            if (j < i) return twoStringMove(tour, j, i, Y, X);
            else {
                vector<int> auxIPlus1ToIPlusX(X);
                vector<int> auxIPlusXPlus1ToJ(j - i - X);
                vector<int> auxJPlus1ToJPlusY(Y);
                for (int k = 0; k < auxIPlus1ToIPlusX.size(); ++k)
                    auxIPlus1ToIPlusX[k] = tour[i + 1 + k];
                for (int k = 0; k < auxIPlusXPlus1ToJ.size(); ++k)
                    auxIPlusXPlus1ToJ[k] = tour[i + X + 1 + k];
                for (int k = 0; k < auxJPlus1ToJPlusY.size(); ++k)
                    auxJPlus1ToJPlusY[k] = tour[j + 1 + k];
                for (int k = 0; k < auxJPlus1ToJPlusY.size(); ++k)
                    tour[i + 1 + k] = auxJPlus1ToJPlusY[k];
                for (int k = 0; k < auxIPlusXPlus1ToJ.size(); ++k)
                    tour[i + Y + 1 + k] = auxIPlusXPlus1ToJ[k];
                for (int k = 0; k < auxIPlus1ToIPlusX.size(); ++k)
                    tour[j + Y - X + 1 + k] = auxIPlus1ToIPlusX[k];
                return true;
            }
        }
        bool twoString(vector<int> &tour, int X, int Y, bool firstImprove) {
            vector<int> tmpTour;
            vector<int> bestTour;
            tmpTour.reserve(tour.size());
            bestTour.reserve(tour.size());
            double impBest = EPS;
            bool improved = false;
            for (int i = 0; i < tour.size() - X; ++i) {
                if (get_evals() > STOP_CNT || (firstImprove && improved)) break;
                for (int j = (X == Y) ? i + X : 0; j < tour.size() - Y; ++j) {
                    if (get_evals() > STOP_CNT || (firstImprove && improved)) break;
                    bool cond1 = i + X == j && (X == 0 || Y == 0);
                    bool cond2 = j + Y == i && (X == 0 || Y == 0);
                    if (!cond1 && !cond2 && ((j - i >= X) || (i - j >= Y))) {
                        auto impCurr = j < i ? twoStringCostUpdate(tour, j, i, Y, X) : twoStringCostUpdate(tour, i, j, X, Y);
                        if (impCurr > impBest) {
                            tmpTour = tour;
                            twoStringMove(tmpTour, i, j, X, Y);
                            if (isValidTour(tmpTour)) {
                                impBest = impCurr;
                                bestTour = tmpTour;
                                improved = true;
                            }
                        }
                    }
                }
            }
            if (improved) tour = bestTour;
            return improved;
        }
        inline bool onePoint(vector<int> &tour, bool firstImprove) {
            return twoString(tour, 0, 1, firstImprove);
        }
        inline bool twoPoint(vector<int> &tour, bool firstImprove) {
            return twoString(tour, 1, 1, firstImprove);
        }
        inline bool threePoint(vector<int> &tour, bool firstImprove) {
            return twoString(tour, 1, 2, firstImprove);
        }
        void mergeAFSs(vector<int> &tour) {
            int i = 0;
            while (i != tour.size() - 1)
                if (tour[i] == tour[i + 1]) tour.erase(tour.begin() + i + 1); else i++;
        }
        void generalizedDoubleBridge(vector<int> &tour, unsigned p) {
            auto n = tour.size();
            p = (p + 1 > n) ? n - 1 : p;
            default_random_engine generator(rand());
            auto distribution = uniform_int_distribution<int>(0, n - 2);
            vector<bool> throwsBool(n, false);
            for (unsigned i = 0; i < p; ++i) {
                int throwI;
                do throwI = distribution(generator);
                while (throwsBool[throwI]);
                throwsBool[throwI] = true;
            }
            vector<vector<int>> pathStrings{};
            for (unsigned i = 0; i < n; ++i) {
                vector<int> pathString{};
                while (!throwsBool[i] && i < n - 1) pathString.emplace_back(i++);
                pathString.emplace_back(i);
                pathStrings.push_back(pathString);
            }
            vector<int> newPositions(pathStrings[0]);
            vector<int> roulette = vector<int>(pathStrings.size());
            for (unsigned i = 0; i < roulette.size(); ++i) roulette[i] = i;
            shuffle(roulette.begin() + 1, roulette.end(), generator);
            for (unsigned k = 1; k < pathStrings.size(); ++k) {
                auto &nextString = pathStrings[roulette[k]];
                auto side = uniform_int_distribution<unsigned>(0, 1)(generator);
                if (side == 0) for (int i : nextString) newPositions.push_back(i);
                else for (unsigned i = 0; i < nextString.size(); ++i)
                        newPositions.push_back(nextString[nextString.size() - i - 1]);
            }
            vector<int> newPath(n);
            for (unsigned j = 0; j < n; ++j) newPath[j] = tour[newPositions[j]];
            newPath = tsp2evrp_zga_relaxed(newPath);
            mergeAFSs(newPath);
            tour = newPath;
        }
        double fitness_evaluation(int *routes, int size) {
            int i;
            double tour_length = 0.0;
            for (i = 0; i < size - 1; i++) tour_length += distances[routes[i]][routes[i + 1]];
            if (tour_length < current_best) current_best = tour_length;
            evals++;
            return tour_length;
        }
        int getRandomAmongAvailable(unsigned availableCount, const vector<bool> &available) {
            default_random_engine generator(rand());
            uniform_int_distribution<int> distribution(0, availableCount - 1);
            auto r = distribution(generator);
            for (unsigned i = 0; i < available.size(); ++i)
                if (available[i]) {
                    if (r == 0) return i;
                    --r;
                }
            return -1;
        }
        void rvnd(vector<int> &tour) {
            int i = 0;
            auto nlSize = static_cast<unsigned>(6);
            vector<bool> available(6, true);
            while (nlSize > 0) {
                auto neighborhoodIdx = getRandomAmongAvailable(nlSize, available);
                bool is_valid;
                if (neighborhoodIdx == 0) is_valid = AFSrealoc_one(tour, false);
                else if (neighborhoodIdx == 1) is_valid = AFSrealoc_more_than_one(tour, false);
                else if (neighborhoodIdx == 2) is_valid = twoOpt(tour, false);
                else if (neighborhoodIdx == 3) is_valid = onePoint(tour, false);
                else if (neighborhoodIdx == 4) is_valid = twoPoint(tour, false);
                else if (neighborhoodIdx == 5) is_valid = threePoint(tour, false);
                if (is_valid) {
                    nlSize = static_cast<unsigned>(6);
                    fill(available.begin(), available.end(), true);
                    mergeAFSs(tour);
                    double a = get_evals();
                    double b = TERMINATION;
                } else {
                    available[neighborhoodIdx] = false;
                    --nlSize;
                }
            }
        }
        vector<int> ms_vns() {
            int vns_restarts = ACTUAL_PROBLEM_SIZE * 0.35;
            int vns_cnt = 0;
            vector<int> very_best;
            double very_best_score;
            while (get_evals() < STOP_CNT) {
                auto best = init_from_dbca();
                double best_score = fitness_evaluation(best);
                if (very_best.size() == 0) {
                    very_best = best;
                    very_best_score = best_score;
                }
                while (vns_cnt < vns_restarts && get_evals() < STOP_CNT) {
                    auto current = best;
                    generalizedDoubleBridge(current, 2);
                    rvnd(current);
                    double current_score = fitness_evaluation(current);
                    double a = get_evals();
                    double b = TERMINATION;
                    if (current_score < best_score) {
                        vns_cnt = 0;
                        best = current;
                        best_score = current_score;
                    } else vns_cnt++;
                }
                vns_cnt = 0;
                if (best_score < very_best_score) {
                    very_best = best;
                    very_best_score = best_score;
                }
            }
            return very_best;
        }
        void run_heuristic() {
            initMyStructures();
            auto evrpTour = ms_vns();
            for (int i = 0; i < evrpTour.size(); i++) best_sol->tour[i] = evrpTour[i];
            best_sol->steps = evrpTour.size();
            best_sol->tour_length = fitness_evaluation(evrpTour);
        }
};
namespace py = pybind11;
class EVRP {
    public:
        EVRPEnv evrp;
        vector<int> solution;
        EVRP(int dimension, int stations, int capacity, int energy_capacity, double energy_consumption, py::list nodes, int seed, bool isRoundingInteger) {
            srand(seed);
            evrp.ACTUAL_PROBLEM_SIZE = dimension + stations;
            evrp.NUM_OF_STATIONS = stations;
            evrp.MAX_CAPACITY = capacity;
            evrp.BATTERY_CAPACITY = energy_capacity;
            evrp.energy_consumption = energy_consumption;
            evrp.NUM_OF_CUSTOMERS = dimension - 1;
            evrp.node_list = new node[evrp.ACTUAL_PROBLEM_SIZE];
            evrp.cust_demand = new int[evrp.ACTUAL_PROBLEM_SIZE];
            evrp.charging_station = new bool[evrp.ACTUAL_PROBLEM_SIZE];
            int node_id = 0;
            for (py::handle node : nodes) {
                double x = node.attr("x").cast<double>();
                double y = node.attr("y").cast<double>();
                int demand = node.attr("demand").cast<int>();
                bool is_depot = node.attr("is_depot").cast<bool>();
                add_node(node_id++, x, y, demand, is_depot);
            }
            evrp.distances = evrp.generate_2D_matrix_double(evrp.ACTUAL_PROBLEM_SIZE, evrp.ACTUAL_PROBLEM_SIZE);
            int i = 0;
            for (py::handle node_i : nodes) {
                int j = 0;
                for (py::handle node_j : nodes) {
                    double d = node_i.attr("distance_to")(node_j).cast<double>();
                    evrp.distances[i][j] = d;
                    j++;
                }
                i++;
            }
            // evrp.compute_distances(isRoundingInteger);
            evrp.init_evals();
            evrp.init_current_best();
            evrp.initialize_heuristic();
            evrp.initMyStructures();
        }
        void init(int dimension, int stations, int capacity, int energy_capacity, double energy_consumption) { 
            evrp.ACTUAL_PROBLEM_SIZE = dimension + stations;
            evrp.NUM_OF_STATIONS = stations;
            evrp.MAX_CAPACITY = capacity;
            evrp.BATTERY_CAPACITY = energy_capacity;
            evrp.energy_consumption = energy_consumption;
            evrp.NUM_OF_CUSTOMERS = dimension - 1;
            evrp.node_list = new node[evrp.ACTUAL_PROBLEM_SIZE];
            evrp.cust_demand = new int[evrp.ACTUAL_PROBLEM_SIZE];
            evrp.charging_station = new bool[evrp.ACTUAL_PROBLEM_SIZE];
        }
        void add_node(int node_id, double x, double y, int demand, bool is_depot) {
            if (is_depot) evrp.DEPOT = node_id;
            evrp.cust_demand[node_id] = demand;
            evrp.charging_station[node_id] = demand == 0 || is_depot;
            evrp.node_list[node_id].id = node_id;
            evrp.node_list[node_id].x = x;
            evrp.node_list[node_id].y = y;
        }
        py::list init_solution() {
            solution = evrp.init_from_dbca();
            return convert_solution(solution);
        }
        py::list step(py::list py_solution) {
            solution = py_solution.cast<vector<int>>();
            return sub_step();
        }
        py::list sub_step() {
            solution = evrp.tsp2evrp_zga_relaxed(solution);
            evrp.mergeAFSs(solution);
            evrp.rvnd(solution);
            return convert_solution(solution);
        }
        py::list convert_solution(vector<int> solution) {
            py::list py_solution = py::cast(solution);
            return py_solution;
		}
        py::list get_best_solution() {
            return convert_solution(solution);
		}
        py::list get_offspring() {
            return convert_solution(solution);
		}
};
PYBIND11_MODULE(evrp_cpp, m) {
    m.doc() = R"pbdoc(Pybind11 example plugin)pbdoc";
    py::class_<EVRP>(m, "EVRP")
        .def(py::init<int, int, int, int, double, py::list &, int, bool>())
        .def("init_solution", &EVRP::init_solution)
        .def("step", &EVRP::step)
        .def("sub_step", &EVRP::sub_step)
        .def("get_best_solution", &EVRP::get_best_solution)
        .def("get_offspring", &EVRP::get_offspring)
        ;
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
