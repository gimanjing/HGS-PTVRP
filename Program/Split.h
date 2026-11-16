/*MIT License

Copyright(c) 2020 Thibaut Vidal

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#ifndef SPLIT_H
#define SPLIT_H

#include "Params.h"
#include "Individual.h"
#include <deque>

// SPLIT FOR THE PERIODIC-TEMPLATE VRP (PT-VRP)
//
// A route (template) sigma is executed m(sigma) = ceil(q(sigma)/Q) times, so
//
//     cost(sigma)     = m(sigma) * d(sigma)
//     duration(sigma) = m(sigma) * tau(sigma)  <=  T_H      (T_H = params.durationLimit, when the instance has one)
//
// The multiplier is a step function of the load, which breaks the dominance rule of the O(n) linear Split.
// The decoder used here is the layered Split of the split-linear repository (Split_Layered_PTVRP_safe) :
// one deque per trip count k, inside which k is constant and k*d(i,j) = k*A[i] + k*B[j] separates again.
// It carries the two repairs of that repository :
//
//   FRONT (pruning)  : a start is dropped on the PATH OUT (route without the leg home), which never shrinks
//                      as j grows, so the drop is permanent without assuming the triangle inequality.
//   BACK (eviction)  : an older start is evicted only on JOINT dominance -- the newcomer must be at least as
//                      cheap AND at least as easy to fit (At[new] <= At[old]).
//
// The horizon is a hard constraint inside the layered passes. Only if no split fits the horizon with the
// allowed number of templates does the code fall back to a Bellman pass with a penalized horizon, so that
// the search always receives a solution (exactly as the original HGS Split tolerates infeasibility).

struct ClientSplit
{
	double demand;
	double serviceTime;
	double d0_x;
	double dx_0;
	double dnext;
	ClientSplit() : demand(0.), serviceTime(0.), d0_x(0.), dx_0(0.), dnext(0.) {};
};

class Split
{

 private:

 // Problem parameters
 const Params & params ;
 int maxVehicles ;

 /* Auxiliary data structures */
 std::vector < ClientSplit > cliSplit;
 std::vector < std::vector < double > > potential;  // potential[k][i] : best cost to serve the first i clients with k templates (k = 0 is also used by the unlimited-fleet pass)
 std::vector < std::vector < int > > pred;          // Index of the predecessor in an optimal path
 std::vector <double> sumDistance; // sumDistance[i] for i > 1 contains the sum of distances : sum_{k=1}^{i-1} d_{k,k+1}
 std::vector <double> sumLoad;     // sumLoad[i] for i >= 1 contains the sum of loads : sum_{k=1}^{i} q_k
 std::vector <double> sumService;  // sumService[i] for i >= 1 contains the sum of service time : sum_{k=1}^{i} s_k

 // d(i,j) = A[i] + B[j],  tau(i,j) = At[i] + Bt[j],  pathout_time(i,j) = At[i] + BoutT[j]
 // for the template visiting the clients i+1..j of the giant tour
 std::vector <double> A, B, At, Bt, BoutT;

 // Work structures of a layered pass, indexed by the trip count k
 std::vector <int> firstLoadLE, firstTimeLE, nextAdd;
 std::vector <char> sortedByKey;
 std::vector < std::deque <int> > queues;

 // m(i,j) : number of executions of the template visiting the clients i+1..j
 inline int trips(int i, int j) const
 {
	 return params.nbTrips(sumLoad[j] - sumLoad[i]);
 }

 // One layered pass : dst[j] = min over the feasible starts i < j of src[i] + m(i,j) * d(i,j).
 // Called with src == dst for the unlimited fleet, and with consecutive fleet layers otherwise.
 void layeredPass(const std::vector <double> & src, std::vector <double> & dst, std::vector <int> & prd);

 // Split for unlimited fleet, returns the number of templates (0 if no split fits the horizon)
 int splitSimple(Individual & indiv);

 // Split for limited fleet, hard horizon. Returns 0 if no split fits the horizon with maxVehicles templates
 int splitLF(Individual & indiv);

 // Last resort : Bellman split for limited fleet with a penalized horizon. Always returns a solution
 int splitLFSoft(Individual & indiv);

 // Fills chromR from the predecessors of the fleet-limited potentials
 void fillRoutesLF(Individual & indiv, int nbRoutes);

public:

  // Counters, to measure how often each decoder is used
  long nbCallsSimple = 0, nbCallsLF = 0, nbCallsLFSoft = 0;

  // Evictions refused by the joint-dominance guard (cheaper start, but harder to fit in the horizon)
  long nbBlockedPops = 0;

  // General Split function (tests the unlimited fleet, and only if it does not produce a solution within the fleet size, runs the Split algorithm for limited fleet).
  // nbMaxVehicles is kept for compatibility with the original interface and is not used, see Split.cpp
  void generalSplit(Individual & indiv, int nbMaxVehicles);

  // Constructor
  Split(const Params & params);

};
#endif
