#include "Split.h" 

void Split::generalSplit(Individual & indiv, int nbMaxVehicles)
{
	// The number of templates is bounded by the fleet size.
	// The original HGS additionally caps the Split of an offspring at the number of routes of its parent
	// (nbMaxVehicles). That cap is not applied here : in the PT-VRP a fleet-limited Split needs one layered
	// pass per vehicle, and in the tests the capped variant was slower without giving better solutions.
	(void)nbMaxVehicles;
	maxVehicles = std::max<int>(1, params.nbVehicles);

	const int n = params.nbClients;
	for (int i = 1; i <= n; i++)
	{
		cliSplit[i].demand = params.cli[indiv.chromT[i - 1]].demand;
		cliSplit[i].serviceTime = params.cli[indiv.chromT[i - 1]].serviceDuration;
		cliSplit[i].d0_x = params.timeCost[0][indiv.chromT[i - 1]];
		cliSplit[i].dx_0 = params.timeCost[indiv.chromT[i - 1]][0];
		if (i < n) cliSplit[i].dnext = params.timeCost[indiv.chromT[i - 1]][indiv.chromT[i]];
		else cliSplit[i].dnext = -1.e30;
		sumLoad[i] = sumLoad[i - 1] + cliSplit[i].demand;
		sumService[i] = sumService[i - 1] + cliSplit[i].serviceTime;
		sumDistance[i] = sumDistance[i - 1] + cliSplit[i - 1].dnext;
	}

	// Separable form of the distance and duration of the template i+1..j
	for (int i = 0; i < n; i++)
	{
		A[i]  = cliSplit[i + 1].d0_x - sumDistance[i + 1];
		At[i] = A[i] - sumService[i];
	}
	for (int j = 1; j <= n; j++)
	{
		B[j]  = sumDistance[j] + cliSplit[j].dx_0;
		Bt[j] = B[j] + sumService[j];
		BoutT[j] = sumDistance[j] + sumService[j]; // path out : same prefix WITHOUT the leg home. Non-decreasing in j, unconditionally.
	}

	// Unlimited fleet first, then limited fleet, then limited fleet with a penalized horizon
	int nbRoutes = splitSimple(indiv);
	if (nbRoutes == 0 || nbRoutes > maxVehicles)
		if (splitLF(indiv) == 0)
			splitLFSoft(indiv);

	// Build up the rest of the Individual structure
	indiv.evaluateCompleteCost(params);
}

void Split::layeredPass(const std::vector <double> & src, std::vector <double> & dst, std::vector <int> & prd)
{
	const int n = params.nbClients;
	const double horizon = params.durationLimit;
	const bool hasHorizon = params.isDurationConstraint;

	for (size_t k = 0; k < queues.size(); k++)
	{
		firstLoadLE[k] = 0; firstTimeLE[k] = 0; nextAdd[k] = 0; sortedByKey[k] = 1;
		queues[k].clear();
	}

	// Path-out frontier : a start dropped here can never fit again, whatever the distances
	int frontierPathOut = 0;

	for (int j = 1; j <= n; j++)
	{
		// Sound : both factors are non-decreasing in j, so once this fires for a start it fires for ever
		if (hasHorizon)
			while (frontierPathOut < j && (At[frontierPathOut] + BoutT[j]) * (double)trips(frontierPathOut, j) > horizon + 1.e-9)
				frontierPathOut++;

		if (frontierPathOut >= j)
			continue; // not even the singleton {j} can fit, on a bound that ignores the leg home

		int kHi = trips(frontierPathOut, j);
		if (kHi + 1 > (int)queues.size())
		{
			int grown = 2 * kHi + 2;
			firstLoadLE.resize(grown + 1, 0);
			firstTimeLE.resize(grown + 1, 0);
			nextAdd.resize(grown, 0);
			sortedByKey.resize(grown, 1);
			queues.resize(grown);
		}

		for (int k = 1; k <= kHi; k++)
		{
			// Starts of layer k at column j : those with m(i,j) == k, a window [loadLower, loadUpper)
			while (firstLoadLE[k] < j && trips(firstLoadLE[k], j) > k)
				firstLoadLE[k]++;

			int lower = firstLoadLE[k];
			int upper = (k == 1) ? j : firstLoadLE[k - 1];
			if (lower >= upper)
				continue;

			if (hasHorizon)
			{
				// Sound per-layer bound : At[i] > T_H/k - BoutT[j] means i cannot fit at k trips now, and BoutT
				// only grows while the trip count only grows, so it cannot fit at k or above ever again
				double timeBound = horizon / (double)k - BoutT[j];
				while (firstTimeLE[k] < j && At[firstTimeLE[k]] > timeBound + 1.e-9)
					firstTimeLE[k]++;
				if (firstTimeLE[k] > lower) lower = firstTimeLE[k];
				if (frontierPathOut > lower) lower = frontierPathOut;
				if (lower >= upper)
					continue;
			}

			std::deque <int> & dq = queues[k];
			while (!dq.empty() && dq.front() < lower)
				dq.pop_front();

			if (nextAdd[k] < lower)
				nextAdd[k] = lower;

			while (nextAdd[k] < upper)
			{
				int i = nextAdd[k];
				nextAdd[k]++;
				if (src[i] > 1.e29)
					continue;
				double ki = src[i] + (double)k * A[i];
				// Evict only on JOINT dominance : the newcomer must be at least as cheap AND at least as feasible.
				// At[new] <= At[old] means "wherever old fits, new fits", at every j.
				// Without a horizon every start fits, and the eviction on cost alone is exact.
				while (!dq.empty() && src[dq.back()] + (double)k * A[dq.back()] >= ki)
				{
					if (hasHorizon && At[i] > At[dq.back()] + 1.e-9)
					{
						// cheaper but harder to fit : keep the incumbent. The layer is no longer sorted by key,
						// so its query has to scan
						sortedByKey[k] = 0;
						nbBlockedPops++;
						break;
					}
					dq.pop_back();
				}
				dq.push_back(i);
			}

			// Feasibility itself is tested exactly, on the true duration. Nothing is popped on time.
			int front = -1;
			double bestKey = 1.e30;
			if (!hasHorizon)
			{
				if (!dq.empty()) { front = dq.front(); bestKey = src[front] + (double)k * A[front]; }
			}
			else if (sortedByKey[k])
			{
				// The first feasible entry is the cheapest one
				for (size_t q = 0; q < dq.size(); q++)
				{
					if ((At[dq[q]] + Bt[j]) * (double)k <= horizon + 1.e-9)
					{
						front = dq[q]; bestKey = src[front] + (double)k * A[front];
						break;
					}
				}
			}
			else
			{
				// The guard has blocked a pop in this layer : cheapest feasible entry anywhere in the deque
				for (size_t q = 0; q < dq.size(); q++)
				{
					if ((At[dq[q]] + Bt[j]) * (double)k > horizon + 1.e-9) continue;
					double kk = src[dq[q]] + (double)k * A[dq[q]];
					if (kk < bestKey) { bestKey = kk; front = dq[q]; }
				}
			}

			if (front < 0)
				continue;

			double cand = bestKey + (double)k * B[j];
			if (cand < dst[j])
			{
				dst[j] = cand;
				prd[j] = front;
			}
		}
	}
}

int Split::splitSimple(Individual & indiv)
{
	nbCallsSimple++;
	const int n = params.nbClients;

	potential[0][0] = 0;
	for (int i = 1; i <= n; i++)
		potential[0][i] = 1.e30;

	layeredPass(potential[0], potential[0], pred[0]);

	if (potential[0][n] > 1.e29)
		return 0;

	int nbRoutes = 0;
	for (int end = n; end != 0; end = pred[0][end])
		nbRoutes++;

	if (nbRoutes > maxVehicles)
		return nbRoutes;

	// Filling the chromR structure
	for (int k = params.nbVehicles - 1; k >= nbRoutes; k--)
		indiv.chromR[k].clear();

	int end = n;
	for (int k = nbRoutes - 1; k >= 0; k--)
	{
		indiv.chromR[k].clear();
		int begin = pred[0][end];
		for (int ii = begin; ii < end; ii++)
			indiv.chromR[k].push_back(indiv.chromT[ii]);
		end = begin;
	}

	return nbRoutes;
}

// Split for problems with limited fleet
int Split::splitLF(Individual & indiv)
{
	nbCallsLF++;
	const int n = params.nbClients;

	potential[0][0] = 0;
	for (int i = 1; i <= n; i++)
		potential[0][i] = 1.e30;
	for (int k = 1; k <= maxVehicles; k++)
		for (int i = 0; i <= n; i++)
			potential[k][i] = 1.e30;

	for (int k = 0; k < maxVehicles; k++)
		layeredPass(potential[k], potential[k + 1], pred[k + 1]);

	// It could be cheaper to use a smaller number of vehicles
	double minCost = 1.e30;
	int nbRoutes = 0;
	for (int k = 1; k <= maxVehicles; k++)
		if (potential[k][n] < minCost)
			{minCost = potential[k][n]; nbRoutes = k;}

	if (minCost > 1.e29)
		return 0;

	fillRoutesLF(indiv, nbRoutes);
	return nbRoutes;
}

// Last resort : no split fits the horizon with maxVehicles templates
int Split::splitLFSoft(Individual & indiv)
{
	nbCallsLFSoft++;
	const int n = params.nbClients;

	for (int k = 0; k <= maxVehicles; k++)
		for (int i = 0; i <= n; i++)
			potential[k][i] = 1.e30;
	potential[0][0] = 0;

	// Simple Split using Bellman's algorithm in topological order, horizon excess penalized
	for (int k = 0; k < maxVehicles; k++)
	{
		for (int i = k; i < n; i++)
		{
			if (potential[k][i] > 1.e29) continue;
			for (int j = i + 1; j <= n; j++)
			{
				double m = (double)trips(i, j);
				double cost = m * (A[i] + B[j])
					+ params.penaltyDuration * std::max<double>(m * (At[i] + Bt[j]) - params.durationLimit, 0.);
				if (potential[k][i] + cost < potential[k + 1][j])
				{
					potential[k + 1][j] = potential[k][i] + cost;
					pred[k + 1][j] = i;
				}
			}
		}
	}

	double minCost = 1.e30;
	int nbRoutes = 0;
	for (int k = 1; k <= maxVehicles; k++)
		if (potential[k][n] < minCost)
			{minCost = potential[k][n]; nbRoutes = k;}

	if (minCost > 1.e29)
		throw std::string("ERROR : no Split solution has been propagated until the last node");

	fillRoutesLF(indiv, nbRoutes);
	return nbRoutes;
}

void Split::fillRoutesLF(Individual & indiv, int nbRoutes)
{
	for (int k = params.nbVehicles - 1; k >= nbRoutes; k--)
		indiv.chromR[k].clear();

	int end = params.nbClients;
	for (int k = nbRoutes - 1; k >= 0; k--)
	{
		indiv.chromR[k].clear();
		int begin = pred[k + 1][end];
		for (int ii = begin; ii < end; ii++)
			indiv.chromR[k].push_back(indiv.chromT[ii]);
		end = begin;
	}
}

Split::Split(const Params & params): params(params)
{
	const int n = params.nbClients;
	cliSplit = std::vector <ClientSplit>(n + 1);
	sumDistance = std::vector <double>(n + 1, 0.);
	sumLoad = std::vector <double>(n + 1, 0.);
	sumService = std::vector <double>(n + 1, 0.);
	A = std::vector <double>(n + 1, 0.);
	At = std::vector <double>(n + 1, 0.);
	B = std::vector <double>(n + 1, 0.);
	Bt = std::vector <double>(n + 1, 0.);
	BoutT = std::vector <double>(n + 1, 0.);
	potential = std::vector < std::vector <double> >(params.nbVehicles + 1, std::vector <double>(n + 1, 1.e30));
	pred = std::vector < std::vector <int> >(params.nbVehicles + 1, std::vector <int>(n + 1, 0));
	firstLoadLE = std::vector <int>(3, 0);
	firstTimeLE = std::vector <int>(3, 0);
	nextAdd = std::vector <int>(2, 0);
	sortedByKey = std::vector <char>(2, 1);
	queues = std::vector < std::deque <int> >(2);
}
