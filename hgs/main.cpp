#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#include <algorithm>
#include "Genetic.h"
#include "commandline.h"
#include "LocalSearch.h"
#include "Split.h"
using namespace std;

namespace py = pybind11;

Params::Params(int capacity, py::list nodes, int nveh, bool isRoundingInteger) {
	int dim = (int) nodes.size();
	cli = vector<Client>(dim);
	nbClients = dim - 1;
	totalDemand = 0;
	vehicleCapacity = capacity;
	durationLimit = 1.e30;
	isRoundingInteger = isRoundingInteger;
	maxDemand = 0;

	maxDist = 0.;
	timeCost = vector < vector< double > >(dim, vector <double>(dim));

	int node_id = 0;
	for (py::handle node : nodes) {
		double x = node.attr("x").cast<double>();
		double y = node.attr("y").cast<double>();
		int demand = node.attr("demand").cast<int>();
		cli[node_id].custNum = node_id;
		cli[node_id].coordX = x;
		cli[node_id].coordY = y;
		cli[node_id].demand = demand;
		cli[node_id].serviceDuration = 0.;
		double angle = atan2(cli[node_id].coordY - cli[0].coordY, cli[node_id].coordX - cli[0].coordX);
		cli[node_id].polarAngle = CircleSector::positive_mod((int)(32768.*angle / PI));
		totalDemand += demand;
		if (cli[node_id].demand > maxDemand) maxDemand = cli[node_id].demand;
		node_id++;
	}
	
	int i = 0;
	for (py::handle node_i : nodes) {
		int j = 0;
		for (py::handle node_j : nodes) {
			double d = node_i.attr("distance_to")(node_j).cast<double>();
			if (d > maxDist) maxDist = d;
			timeCost[i][j] = d;
			j++;
		}
		i++;
	}

	if (nveh <= 0) {
		nbVehicles = (int) ceil(1.3*totalDemand/vehicleCapacity) + 3;
	} else {
		nbVehicles = nveh;
	}

	correlatedVertices = vector < vector < int > >(dim);
	vector < set < int > > setCorrelatedVertices = vector < set <int> >(dim);
	vector < pair <double, int> > orderProximity;
	for (int i = 1; i < dim; i++)
	{
		orderProximity.clear();
		for (int j = 1; j < dim; j++)
			if (i != j) orderProximity.push_back(pair <double, int>(timeCost[i][j], j));
		sort(orderProximity.begin(), orderProximity.end());
		for (int j = 0; j < min<int>(nbGranular, nbClients - 1); j++)
		{
			setCorrelatedVertices[i].insert(orderProximity[j].second);
			setCorrelatedVertices[orderProximity[j].second].insert(i);
		}
	}

	for (int i = 1; i < dim; i++)
		for (int x : setCorrelatedVertices[i])
			correlatedVertices[i].push_back(x);
	
	penaltyDuration = 1;
	penaltyCapacity = max<double>(0.1, min<double>(1000.,maxDist / maxDemand));
}

class CVRP {
    public:
		Params* params;
		Split* split;
		LocalSearch* localSearch;
		Population* population;
		Genetic* solver;
		
        CVRP(int capacity, py::list nodes, int seed, int nveh, bool isRoundingInteger) {
			srand(seed);
			params = new Params(capacity, nodes, nveh, isRoundingInteger);
			split = new Split(params);
			localSearch = new LocalSearch(params);
			population = new Population(params, split, localSearch);
			solver = new Genetic(params, split, population, localSearch);
        }
		py::list init_solution() {
            return convert_solution(solver->population->bestSolutionOverall.chromR);
        }
		void read_solution(py::list py_solution) {
			vector<int> solution = py_solution.cast<vector<int>>();
			int n_chromT = (int) solver->offspring->chromT.size();
			int n_chromR = (int) solver->offspring->chromR.size();
			for (int i = 0; i < n_chromR; i++)
				solver->offspring->chromR[i].clear();
			int posT = 0;
			int posR = 0;
			for (int i = 0; i < (int) solution.size(); i++)
			{
				int id = solution[i];
				if (id != 0) {
					solver->offspring->chromT[posT] = id;
					solver->offspring->chromR[posR].push_back(id);
					if (posT < n_chromT-1) posT++;
				} else {
					if (posR < n_chromR-1) posR++;
				}
			}
		}
        py::list step(py::list py_solution) {
			read_solution(py_solution);
			solver->run(64, 4);
			return convert_solution(solver->offspring->chromR);
        }
        void sub_step() {
            solver->run(1, 4);
        }
        py::list convert_solution(vector < vector <int> >  sol) {
            vector<int> solution;
			solution.push_back(0);
			for (int i = 0; i < (int)sol.size(); i++)
			{
				if (sol[i].size() == 0) {
					continue;
				}
				for (int j = 0; j < (int)sol[i].size(); j++)
				{
					solution.push_back(sol[i][j]);
				}
				solution.push_back(0);
			}
            py::list py_solution = py::cast(solution);
            return py_solution;
		}
        py::list get_best_solution() {
            return convert_solution(solver->population->bestSolutionRestart.chromR);
		}
		py::list get_offspring() {
            return convert_solution(solver->offspring->chromR);
		}
};

PYBIND11_MODULE(cvrp_cpp, m) {
    m.doc() = R"pbdoc(Pybind11 example plugin)pbdoc";
    py::class_<CVRP>(m, "CVRP")
        .def(py::init<int, py::list &, int, int, bool>())
        .def("init_solution", &CVRP::init_solution)
        .def("step", &CVRP::step)
        .def("sub_step", &CVRP::sub_step)
        .def("get_best_solution", &CVRP::get_best_solution)
        .def("get_offspring", &CVRP::get_offspring)
        ;
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
