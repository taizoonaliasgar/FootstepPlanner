#ifndef SHARED_DATA
#define SHARED_DATA

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Sparse"
#include "global_loco_opts.h"

#include "mutex"
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#define SET_DATA 1
#define GET_DATA 0
#define HL_DATA 1
#define LL_DATA 0
#define SIM_DATA 2
boost::mutex mtx;

struct sharedData
{
	// Provided by HL
	Eigen::Matrix<double,  12, 1> comDes = Eigen::Matrix<double, 12, 1>::Zero();
	Eigen::Matrix<double,  17, 1> fDes   = Eigen::Matrix<double, 17, 1>::Zero();
	int MPC_cnt = 0;
	
	int ind[4] = {1};
	// Provided by LL
	size_t gait = STAND;
	int control_Tick = 0;

	double domLen = 200;
	double phaseVar = 0;
	
	double q[18] = {0};
	double dq[18] = {0};
	double rotMatrixDouble[9] = {0};
	
	Eigen::Matrix<double, 12, 1> QPforce = Eigen::Matrix<double, 12, 1>::Zero();
	Eigen::VectorXd tau = Eigen::MatrixXd::Zero(18,1);
	Eigen::Matrix<double, 3, 4> toePos = Eigen::MatrixXd::Zero(3,4);

	int solvetime = 0;

	int ind_LL[4] = {1,1,1,1};
	Eigen::Matrix<double, 3, 4> toe_prev = Eigen::MatrixXd::Zero(3,4);
	
};

sharedData data;

void updateData(int setget, int highlow, sharedData *newData){
	// set=1,  get=0
	// high=1, low=0
	// boost::lock_guard<boost::mutex> guard(mtx);
	if(setget==SET_DATA){
		
		if(highlow==SIM_DATA){ // set high level data
			memcpy(data.q,newData->q,18*sizeof(double));
			memcpy(data.dq,newData->dq,18*sizeof(double));
			memcpy(data.rotMatrixDouble,newData->rotMatrixDouble,9*sizeof(double));
			data.control_Tick = newData->control_Tick;
			
		
		}else if(highlow==HL_DATA){ // set low level data
			Eigen::Matrix<double, 17, 1> fDes_temp = data.fDes;
			Eigen::Matrix<double, 12, 1> comDes_temp = data.comDes;
			memcpy(&data,newData,sizeof(sharedData));
			data.fDes = fDes_temp;
			data.comDes = comDes_temp;
			memcpy(data.ind,newData->ind,4*sizeof(int));
			data.solvetime = newData->solvetime;
		
		}else{
			data.tau = newData->tau;
			data.QPforce = newData->QPforce;
			data.toePos = newData->toePos;
			data.toe_prev = newData->toe_prev;
			memcpy(data.ind_LL,newData->ind_LL,4*sizeof(int));
		}
	}else{
		
		if(highlow==SIM_DATA){ // get data for high level
			newData->tau = data.tau;
			memcpy(newData->ind_LL,data.ind_LL,4*sizeof(int));
			newData->toePos = data.toePos;
			newData->toe_prev = data.toe_prev;
		
		}else if(highlow==HL_DATA){ // get data for low level
			memcpy(newData->q,data.q,18*sizeof(double)); 
			memcpy(newData->dq,data.dq,18*sizeof(double));
			newData->control_Tick = data.control_Tick;
			newData->QPforce = data.QPforce;
			newData->toePos = data.toePos;
		
		}else{
			Eigen::Matrix<double, 17, 1> fDes_temp = newData->fDes;
			Eigen::Matrix<double, 12, 1> comDes_temp = newData->comDes;
			memcpy(newData,&data,sizeof(sharedData));
			newData->fDes = fDes_temp;
			newData->comDes = comDes_temp;
			
			memcpy(newData->q,data.q,18*sizeof(double));
			memcpy(newData->dq,data.dq,18*sizeof(double));
			memcpy(newData->rotMatrixDouble,data.rotMatrixDouble,9*sizeof(double));
			
			memcpy(newData->ind,data.ind,4*sizeof(int));
			newData->control_Tick = data.control_Tick;

			newData->solvetime = data.solvetime;

		}

	}
};


#endif
