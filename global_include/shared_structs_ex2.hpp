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
#define IMU_DATA 3
boost::mutex mtx;

struct sharedData
{
	// Provided by HL
	Eigen::Matrix<double,  12, 1> comDes = Eigen::Matrix<double, 12, 1>::Zero();
	Eigen::Matrix<double,  17, 1> fDes   = Eigen::Matrix<double, 17, 1>::Zero();
	int MPC_cnt = 0;
	
	int ind[5] = {1};
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

	double att_euler[3] = {0};
	double comp_angular_rate[3] = {0};
	Eigen::Matrix<double, 3, 1> IMUacc = Eigen::MatrixXd::Zero(3,1); // scaled acc for IMU data
};

sharedData data;

void updateData(int setget, int highlow, sharedData *newData){
	// set=1,  get=0
	// high=1, low=0
	boost::lock_guard<boost::mutex> guard(mtx);
	if(setget==SET_DATA){
		
		if(highlow==SIM_DATA){ // set high level data
			memcpy(data.q,newData->q,18*sizeof(double));
			memcpy(data.dq,newData->dq,18*sizeof(double));
			memcpy(data.rotMatrixDouble,newData->rotMatrixDouble,9*sizeof(double));
			data.control_Tick = newData->control_Tick;
			
		
		}else if(highlow==HL_DATA){ // set low level data
			data.fDes = newData->fDes;
			data.comDes = newData->comDes;
			memcpy(data.ind,newData->ind,5*sizeof(int));
			data.solvetime = newData->solvetime;
		
		}else if(highlow==LL_DATA){
			data.tau = newData->tau;
			data.QPforce = newData->QPforce;
			data.toePos = newData->toePos;
			data.toe_prev = newData->toe_prev;
			memcpy(data.ind_LL,newData->ind_LL,4*sizeof(int));
		
		}else{
			memcpy(data.att_euler,newData->att_euler,3*sizeof(double));
			memcpy(data.comp_angular_rate,newData->comp_angular_rate,3*sizeof(double));
		}
	}else{
		
		if(highlow==SIM_DATA){ // get data for high level
			newData->tau = data.tau;
			memcpy(newData->ind_LL,data.ind_LL,4*sizeof(int));
			newData->toePos = data.toePos;
			newData->toe_prev = data.toe_prev;
			memcpy(newData->att_euler,data.att_euler,3*sizeof(double));
			memcpy(newData->comp_angular_rate,data.comp_angular_rate,3*sizeof(double));
		
		}else if(highlow==HL_DATA){ // get data for low level
			memcpy(newData->q,data.q,18*sizeof(double)); 
			memcpy(newData->dq,data.dq,18*sizeof(double));
			newData->control_Tick = data.control_Tick;
			newData->QPforce = data.QPforce;
			newData->toePos = data.toePos;
		
		}else{
			
			newData->fDes = data.fDes;
			newData->comDes = data.comDes;
			
			memcpy(newData->q,data.q,18*sizeof(double));
			memcpy(newData->dq,data.dq,18*sizeof(double));
			memcpy(newData->rotMatrixDouble,data.rotMatrixDouble,9*sizeof(double));
			
			memcpy(newData->ind,data.ind,5*sizeof(int));
			newData->control_Tick = data.control_Tick;

			newData->solvetime = data.solvetime;

		}

	}
};

void updateDataExp(int setget, int highlow, sharedData *newData){
	// set=1,  get=0
	// high=1, low=0
	boost::lock_guard<boost::mutex> guard(mtx);
	if(setget==SET_DATA){
		
		if(highlow==HL_DATA){ // set low level data
			data.fDes = newData->fDes;
			data.comDes = newData->comDes;
			memcpy(data.ind,newData->ind,5*sizeof(int));
			data.solvetime = newData->solvetime;
		
		}else if(highlow==LL_DATA){
			
			data.control_Tick = newData->control_Tick;
			memcpy(data.q,newData->q,18*sizeof(double));
			memcpy(data.dq,newData->dq,18*sizeof(double));
			data.QPforce = newData->QPforce;
			data.toePos = newData->toePos;
			
		}else{
			memcpy(data.att_euler,newData->att_euler,3*sizeof(double));
			memcpy(data.comp_angular_rate,newData->comp_angular_rate,3*sizeof(double));
			// memcpy(data.IMUacc,newData->IMUacc,3*sizeof(double));
			data.IMUacc = newData->IMUacc;
		}
	}else{
		
		if(highlow==HL_DATA){ // get data for low level
			memcpy(newData->q,data.q,18*sizeof(double)); 
			memcpy(newData->dq,data.dq,18*sizeof(double));
			newData->control_Tick = data.control_Tick;
			newData->QPforce = data.QPforce;
			newData->toePos = data.toePos;
		
		}else{
			
			newData->fDes = data.fDes;
			newData->comDes = data.comDes;
			memcpy(newData->ind,data.ind,5*sizeof(int));
			newData->solvetime = data.solvetime;
			memcpy(newData->att_euler,data.att_euler,3*sizeof(double));
			memcpy(newData->comp_angular_rate,data.comp_angular_rate,3*sizeof(double));
			// memcpy(newData->IMUacc,data.IMUacc,3*sizeof(double));
			newData->IMUacc = data.IMUacc;

		}

	}
};
#endif
