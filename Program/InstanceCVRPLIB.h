//
// Created by chkwon on 3/22/22.
//

#ifndef INSTANCECVRPLIB_H
#define INSTANCECVRPLIB_H
#include<string>
#include<vector>

class InstanceCVRPLIB
{
public:
	std::vector<double> x_coords;
	std::vector<double> y_coords;
	std::vector< std::vector<double> > dist_mtx;
	std::vector<double> service_time;
	std::vector<double> demands;
	double durationLimit = 86400.;							// PT-VRP horizon T_H : 24 hours by default, replaced by MAX_ROUTE or HORIZON if the file has one (distances double as travel times)
	double vehicleCapacity = 1.e30;							// Capacity limit
	bool isDurationConstraint = true; 						// Indicates if the problem includes duration constraints
	int nbClients ;											// Number of clients (excluding the depot)

	InstanceCVRPLIB(std::string pathToInstance, bool isRoundingInteger);
};


#endif //INSTANCECVRPLIB_H
