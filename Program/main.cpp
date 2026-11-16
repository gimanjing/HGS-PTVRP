#include "Genetic.h"
#include "commandline.h"
#include "LocalSearch.h"
#include "Split.h"
#include "InstanceCVRPLIB.h"
#include <sstream>
#include <iomanip>
#include <fstream>
using namespace std;

// Appends one line describing the run to the results file
void writeResults(const CommandLine & commandline, const Params & params, Genetic & solver)
{
	// Instance name and the folder that contains it (the instance set)
	std::string path = commandline.pathInstance;
	for (char & c : path) if (c == '\\') c = '/';
	size_t slash = path.find_last_of('/');
	std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
	std::string folder = (slash == std::string::npos) ? "" : path.substr(0, slash);
	size_t slash2 = folder.find_last_of('/');
	if (slash2 != std::string::npos) folder = folder.substr(slash2 + 1);
	if (name.size() > 4 && name.substr(name.size() - 4) == ".vrp") name = name.substr(0, name.size() - 4);

	const Individual * best = solver.population.getBestFound();
	std::ostringstream line;
	line << std::setprecision(12);
	line << folder << "," << name << "," << params.nbClients << "," << params.vehicleCapacity << "," << params.durationLimit << "," << params.ap.seed << ",";
	if (best != NULL)
	{
		// Largest share of the horizon used by a template, and largest number of executions
		double maxUse = 0.; int maxTrips = 0;
		for (const std::vector<int> & route : best->chromR)
		{
			if (route.empty()) continue;
			double load = 0., time = params.timeCost[0][route[0]];
			for (size_t i = 0; i < route.size(); i++)
			{
				load += params.cli[route[i]].demand;
				time += params.cli[route[i]].serviceDuration;
				time += (i + 1 < route.size()) ? params.timeCost[route[i]][route[i + 1]] : params.timeCost[route[i]][0];
			}
			int m = params.nbTrips(load);
			if (m > maxTrips) maxTrips = m;
			if ((double)m * time / params.durationLimit > maxUse) maxUse = (double)m * time / params.durationLimit;
		}
		double timeBest = (double)(solver.population.getSearchProgress().back().first - params.startTime) / (double)CLOCKS_PER_SEC;
		line << "1," << best->eval.penalizedCost << "," << best->eval.nbRoutes << "," << best->eval.nbTrips << "," << maxTrips << "," << maxUse << "," << timeBest << ",";
	}
	else line << "0,,,,,,,";
	line << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << "," << solver.nbIterTotal << "\n";

	// Header if the file is new or empty. The line is sent in one write, so that parallel runs can share a file
	bool needHeader = true;
	{ std::ifstream check(commandline.pathResults); needHeader = !check.good() || check.peek() == std::ifstream::traits_type::eof(); }
	std::ofstream out(commandline.pathResults, std::ios::app);
	if (!out.is_open()) throw std::string("Impossible to open results file: " + commandline.pathResults);
	std::string text = line.str();
	if (needHeader) text = "set,instance,n,capacity,horizon,seed,feasible,cost,templates,trips,max_trips,max_horizon_use,time_to_best_s,time_s,iterations\n" + text;
	out << text;
	out.flush();
}

int main(int argc, char *argv[])
{
	try
	{
		// Reading the arguments of the program
		CommandLine commandline(argc, argv);

		// Print all algorithm parameter values
		if (commandline.verbose) print_algorithm_parameters(commandline.ap);

		// Reading the data file and initializing some data structures
		if (commandline.verbose) std::cout << "----- READING INSTANCE: " << commandline.pathInstance << std::endl;
		InstanceCVRPLIB cvrp(commandline.pathInstance, commandline.isRoundingInteger);

		Params params(cvrp.x_coords,cvrp.y_coords,cvrp.dist_mtx,cvrp.service_time,cvrp.demands,
			          cvrp.vehicleCapacity,cvrp.durationLimit,commandline.nbVeh,cvrp.isDurationConstraint,commandline.verbose,commandline.ap);

		// Running HGS
		Genetic solver(params);
		solver.run();
		
		if (commandline.pathResults != "") writeResults(commandline, params, solver);

		// Exporting the best solution
		if (solver.population.getBestFound() != NULL)
		{
			if (params.verbose) std::cout << "----- WRITING BEST SOLUTION IN : " << commandline.pathSolution << std::endl;
			solver.population.exportCVRPLibFormat(*solver.population.getBestFound(),commandline.pathSolution);
			solver.population.exportSearchProgress(commandline.pathSolution + ".PG.csv", commandline.pathInstance);
		}
	}
	catch (const string& e) { std::cout << "EXCEPTION | " << e << std::endl; }
	catch (const std::exception& e) { std::cout << "EXCEPTION | " << e.what() << std::endl; }
	return 0;
}
