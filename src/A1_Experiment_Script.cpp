//
// Authror: Randy Fawcett on 12/2021.
//
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "raisim/OgreVis.hpp"
//#include "randyImguiPanel.hpp"
#include "raisimBasicImguiPanel.hpp"
#include "raisimKeyboardCallback.hpp" // THIS IS WHERE cmd.vel AND cmd.pose ARE IMPLEMENTED (as globals)
#include "raisim/RaisimServer.hpp"
#include "helper.hpp"
//#include "helper2.hpp"
#include "Filters.h"
#include <iostream>
#include <filesystem>
#include <fstream>

#include "LocoWrapper.hpp"
#include "SRBNMPC.hpp"
#include "A1_Dynamics_full.h"

using namespace std;
// namespace fs = std::filesystem;

FiltStruct_f* filt = (FiltStruct_f*)malloc(sizeof(FiltStruct_f));

void comparingfootv(raisim::Vec<3> FRfoot,raisim::Vec<3> FLfoot,raisim::Vec<3> RRfoot,raisim::Vec<3> RLfoot,Eigen::Matrix<double, 3, 4> foot_vel){
    std::cout << FRfoot[0]-foot_vel(0,0) << "\t" << FRfoot[1]-foot_vel(1,0) << "\t" << FRfoot[2]-foot_vel(2,0) << "\t" <<
                    FLfoot[0]-foot_vel(0,1) << "\t" << FLfoot[1]-foot_vel(1,1) << "\t" << FLfoot[2]-foot_vel(2,1) << "\t" <<
                        RRfoot[0]-foot_vel(0,2) << "\t" << RRfoot[1]-foot_vel(1,2) << "\t" << RRfoot[2]-foot_vel(2,2) << "\t" <<
                            RLfoot[0]-foot_vel(0,3) << "\t" << RLfoot[1]-foot_vel(1,3) << "\t" << RLfoot[2]-foot_vel(2,3) << std::endl;

}

double distx(size_t controlTick, double amp){
    
    return amp*sin((controlTick-8000)/1000.0/M_PI);
    //return dist;
};

void setupCallback() {
    raisim::OgreVis *vis = raisim::OgreVis::get();

    /// light
    vis->getLight()->setDiffuseColour(1, 1, 1);
    vis->getLight()->setCastShadows(false);
    Ogre::Vector3 lightdir(-3,3,-0.5); // Light shines on ROBOTS top/front/right side
    // Ogre::Vector3 lightdir(-3,-3,-0.5); // Light shines on ROBOTS top/front/left side
    lightdir.normalise();
    vis->getLightNode()->setDirection({lightdir});
    vis->setCameraSpeed(300);

    vis->addResourceDirectory(raisim::loadResource("material"));
    vis->loadMaterialFile("myMaterials.material");

    vis->addResourceDirectory(vis->getResourceDir() + "/material/skybox/violentdays");
    vis->loadMaterialFile("violentdays.material");

    /// shdow setting
    vis->getSceneManager()->setShadowTechnique(Ogre::SHADOWTYPE_TEXTURE_ADDITIVE);
    vis->getSceneManager()->setShadowTextureSettings(2048, 3);

    /// scale related settings!! Please adapt it depending on your map size
    // beyond this distance, shadow disappears
    vis->getSceneManager()->setShadowFarDistance(10);
    // size of contact points and contact forces
    vis->setContactVisObjectSize(0.03, 0.6);
    // speed of camera motion in freelook mode
    vis->getCameraMan()->setTopSpeed(5);
}

void kinestimator(double q[18], double dq[18], const int* contact, Eigen::Matrix<double,3,3> R, int robotdown){
    
    //float numContact = (weightedCon[0]+weightedCon[1]+weightedCon[2]+weightedCon[3]);

    float numContact = contact[0]+contact[1]+contact[2]+contact[3];
	// ================================== //
	// ========= Kin Estimator ========== //
	// ================================== //

	// toe pos
	double fr_toe[3], fl_toe[3], rl_toe[3], rr_toe[3];
	static double COM[3]= {0,0,0};
    double Jfr_toe[54], Jfl_toe[54], Jrl_toe[54], Jrr_toe[54];
	double COM_vel[3] = {0,0,0};
	

	q[0] = 0; q[1] = 0; q[2] = 0;
    if(robotdown){
	    FK_FR_toe(fr_toe, q); FK_FL_toe(fl_toe, q);
	    FK_RR_toe(rr_toe, q); FK_RL_toe(rl_toe, q);
        J_FR_toe(Jfr_toe, q); J_FL_toe(Jfl_toe, q);
	    J_RR_toe(Jrr_toe, q); J_RL_toe(Jrl_toe, q);
    }else{
        FK_FR_toe_u(fr_toe, q); FK_FL_toe_u(fl_toe, q);
	    FK_RR_toe_u(rr_toe, q); FK_RL_toe_u(rl_toe, q);
        J_FR_toe_u(Jfr_toe, q); J_FL_toe_u(Jfl_toe, q);
	    J_RR_toe_u(Jrr_toe, q); J_RL_toe_u(Jrl_toe, q);
    }
	
	// update change in com pos
	static double fr_prev[3] = {fr_toe[0],fr_toe[1],fr_toe[2]};
	static double fl_prev[3] = {fl_toe[0],fl_toe[1],fl_toe[2]};
	static double rr_prev[3] = {rr_toe[0],rr_toe[1],rr_toe[2]};
	static double rl_prev[3] = {rl_toe[0],rl_toe[1],rl_toe[2]};
	
    double deltaPos[2] = {0.0};
    for(int i=0; i<2; ++i){
        deltaPos[i] -= (fr_toe[i]-fr_prev[i])*contact[0];
        deltaPos[i] -= (fl_toe[i]-fl_prev[i])*contact[1];
        deltaPos[i] -= (rr_toe[i]-rr_prev[i])*contact[2];
        deltaPos[i] -= (rl_toe[i]-rl_prev[i])*contact[3];
        deltaPos[i] /= numContact;
    }    
    
	COM[0] += deltaPos[0];
	COM[1] += deltaPos[1];
    COM[2]  = -1.0*(fr_toe[2]*contact[0]+fl_toe[2]*contact[1]+rr_toe[2]*contact[2]+rl_toe[2]*contact[3])/numContact;
	
	for(int i=0; i<3; ++i){
		fr_prev[i] = fr_toe[i]; fl_prev[i] = fl_toe[i];
		rr_prev[i] = rr_toe[i]; rl_prev[i] = rl_toe[i];		
	}
	
	numContact = (contact[0]+contact[1])*robotdown + contact[2]+contact[3];
	Eigen::Matrix<double,3,1> dq_temp = {dq[3],dq[4],dq[5]};
	toWorld(&dq[3],dq_temp,R);
	for (int i = 3; i < 18; ++i){
		COM_vel[0] -= (Jfr_toe[3*i+0]*contact[0]*robotdown + Jfl_toe[3*i+0]*contact[1]*robotdown + Jrr_toe[3*i+0]*contact[2] + Jrl_toe[3*i+0]*contact[3])*dq[i];
	 	COM_vel[1] -= (Jfr_toe[3*i+1]*contact[0]*robotdown + Jfl_toe[3*i+1]*contact[1]*robotdown + Jrr_toe[3*i+1]*contact[2] + Jrl_toe[3*i+1]*contact[3])*dq[i];
	 	COM_vel[2] -= (Jfr_toe[3*i+2]*contact[0]*robotdown + Jfl_toe[3*i+2]*contact[1]*robotdown + Jrr_toe[3*i+2]*contact[2] + Jrl_toe[3*i+2]*contact[3])*dq[i];
	}
	COM_vel[0] /= numContact;
	COM_vel[1] /= numContact;
	COM_vel[2] /= numContact;
	
	dq_temp = {dq[3],dq[4],dq[5]};
	toBody(&dq[3],dq_temp,R);

	// Set results
	q[0] = COM[0]; q[1] = COM[1]; q[2] = COM[2]+0.02;
	dq[0] = COM_vel[0]; dq[1] = COM_vel[1]; dq[2] = COM_vel[2];

}

void plannerNMPC(size_t controlTick, LocoWrapper *loco_obj, SRBNMPC* loco_plan, casadi::Function solver, Eigen::Matrix<double,16,1> q0, Eigen::Matrix<double, 3, 4> foot_position, Eigen::Matrix<double, 12, 1> lastQPforce){//} casadi::Function solver, int duration_data) {
    std::map<std::string, casadi::DM> arg, res;

    int controlMPC = std::floor(controlTick/10); 
    //std::cout << "controlMPC:" << controlMPC << std::endl;
    //casadi::DM X_prev = loco_plan->getprevioussol_ll(q0,foot_position,controlMPC);//casadi::DM::zeros(NFS*(HORIZ+1)+NFI*HORIZ,1); 
    casadi::DM X_prev = loco_plan->getprevioussol_fullsim(q0,foot_position,lastQPforce,controlMPC);//casadi::DM::zeros(NFS*(HORIZ+1)+NFI*HORIZ,1);
    if(controlMPC==2740){
        q0.block(12,0,4,1) << foot_position(0,0),foot_position(0,1),foot_position(0,2),foot_position(0,3);//0.15,0.15,-0.1,-0.1;
    }else{
        q0(12) = double(X_prev(12));
        q0(13) = double(X_prev(13));
        q0(14) = double(X_prev(14));
        q0(15) = double(X_prev(15));
    }
            
    // if(controlMPC<2740){
    //     casadi::DM p = loco_plan->motionPlannerN3(q0,controlMPC);
    // }else{
        casadi::DM p = loco_plan->motionPlannerN(q0,controlMPC);
    // }
    if(controlMPC<2740){
        p = loco_plan->motionPlannerN3(q0,controlMPC);
        //std::cout << "controlMPC:" << controlMPC << std::endl;
    }
    loco_plan->setpreviousp(p);

    arg["lbx"] = loco_plan->lowerboundx(p, controlMPC);
    arg["ubx"] =  loco_plan->upperboundx(p);
    arg["lbg"] =  loco_plan->lowerboundg();
    arg["ubg"] =  loco_plan->upperboundg();
    arg["x0"] = X_prev;
    arg["p"] = p;
            
    auto start = std::chrono::high_resolution_clock::now();
    res = solver(arg);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    int duration_data = static_cast<int>(duration.count());
    //std::cout << "NMPC Solve Time: " << duration_data << "ms" << std::endl;

    loco_plan->setprevioussol(res.at("x"));
            
    Eigen::Matrix<double, 33, 1> opt_HLMPC_state = loco_plan->getNMPCsol2(controlMPC);
            
    loco_obj->setoptNLstate(opt_HLMPC_state);
    loco_obj->setcontactconfig(controlMPC);
    loco_plan->mpcdataLog(q0, opt_HLMPC_state.block(16,0,12,1), controlMPC, Eigen::Matrix<double, 12, 1>::Zero());
            
}

void controller(std::vector<raisim::ArticulatedSystem *> A1, LocoWrapper *loco_obj, SRBNMPC* loco_plan, casadi::Function solver, size_t controlTick, raisim::Contact contactInstance, std::ofstream &file_est){
    /////////////////////////////////////////////////////////////////////
    //////////////////////////// INITIALIZE
    /////////////////////////////////////////////////////////////////////
    static size_t loco_kind = TROT;                        // Gait pattern to use
    size_t pose_kind = POSE_CMD;                   // Pose type to use (if loco_kind is set to POSE)
    size_t settling = 0.2*ctrlHz;                   // Settling down
    size_t duration = 1.8*ctrlHz;                   // Stand up 
    size_t loco_start = settling + duration;        // Start the locomotion pattern

    size_t shifttime = 0.75*ctrlHz;
    size_t movetime = 0.3*ctrlHz;
    size_t shifttime2 = 1.0*ctrlHz;
    size_t movetime2 = 0.3*ctrlHz;
    size_t movetime3 = 0.2*ctrlHz;

    double switchtime = 24;

    double *tau;
    //double tau[18] = {0};
    double jpos[18], jvel[18], jpos_est[18], jvel_est[18];
    
    Eigen::VectorXd jointTorqueFF = Eigen::MatrixXd::Zero(TOTAL_DOF,1);
    Eigen::VectorXd jointPosTotal = Eigen::MatrixXd::Zero(TOTAL_DOF+1,1);
    Eigen::VectorXd jointVelTotal = Eigen::MatrixXd::Zero(TOTAL_DOF,1);    
        
    raisim::Mat<3,3> rotMat;
    double rotMatrixDouble[9] = {1,0,0,0,1,0,0,0,1};
    Eigen::Matrix<double, 3, 1> eul;
    Eigen::Matrix<double, 4, 1> quat;

    //For Taizoon High level
    Eigen::Matrix<double, 12, 1> QP_Force;
    const int* foot_state;
    Eigen::Matrix<double, 3, 4> foot_position;
    Eigen::Matrix<double, 3, 4> hip_position;
    Eigen::Matrix<double, 33, 1> opt_HLMPC_state;

    Eigen::Matrix<double, 3, 1>  eul_state = Eigen::MatrixXd::Zero(3,1);
    Eigen::Matrix<double, 3, 1>  omega_state = Eigen::MatrixXd::Zero(3,1);

    static double eul_prev[3] = {0,0,0};
    static double rotMat_prev[9] = {1,0,0,0,1,0,0,0,1};
    /////////////////////////////////////////////////////////////////////
    //////////////////////////// UPDATE STATE
    /////////////////////////////////////////////////////////////////////
    A1.back()->getState(jointPosTotal, jointVelTotal);
    A1.back()->getBaseOrientation(rotMat);

    int robotdown = controlTick < switchtime*ctrlHz ? 1 : 0;
    
    if(!robotdown){

        auto imu = A1.back()->getSensorSet("imu_parent")->getSensor<raisim::InertialMeasurementUnit>("imu");
        auto imu_o = imu->getOrientation();   // Quaternion
        auto imu_w = imu->getAngularVelocity();  // Angular velocity in radians/s
        quat(0) = imu_o[0];
        quat(1) = imu_o[1];
        quat(2) = imu_o[2];
        quat(3) = imu_o[3];

        omega_state(0) = imu_w[0];
        omega_state(1) = imu_w[1];
        omega_state(2) = imu_w[2];

        quat_to_XYZ(quat,eul_state);
        Eigen::Matrix<double,3,3> rotIMU = Eigen::MatrixXd::Zero(3,3);
        quat_to_R(quat,rotIMU);
        for(size_t i=0;i<3;i++){
            for (size_t j = 0; j < 3; j++){
                rotMat[3*i+j] = rotIMU(j,i);
            }
        }
            
    }else{
    
        quat = jointPosTotal.block(3,0,4,1);
        quat_to_XYZ(quat,eul_state);
    }


    // if(std::isnan(eul_state(0)) || std::isnan(eul_state(1)) || std::isnan(eul_state(2))){
    //     eul_state(0) = eul_prev[0];
    //     eul_state(1) = eul_prev[1];
    //     eul_state(2) = eul_prev[2];
    //     for(size_t i=0;i<9;i++){
    //         rotMatrixDouble[i] = rotMat_prev[i];
    //     }
    //     std::cout << "=====================================" << std::endl;
    // }else{
    //     eul_prev[0] = eul_state(0);
    //     eul_prev[1] = eul_state(1);
    //     eul_prev[2] = eul_state(2);
    //     for(size_t i=0;i<9;i++){
    //         rotMatrixDouble[i] = rotMat[i];
    //         rotMat_prev[i] = rotMat[i];
    //     }
    // }

    for(size_t i=0;i<9;i++){
        rotMatrixDouble[i] = rotMat[i];
    }
    
    Eigen::Map< Eigen::Matrix<double, 3, 3> > rotE(rotMatrixDouble, 3, 3);
    
    if(robotdown){
        omega_state = rotE.transpose()*jointVelTotal.segment(3,3); // convert to body frame, like robot measurements
    }

    for(size_t i=0; i<3; ++i){
        jpos[i] = jointPosTotal(i);
        jvel[i] = jointVelTotal(i);
        jpos_est[i] = jointPosTotal(i);
        jvel_est[i] = jointVelTotal(i);
        jpos[i+3] = eul_state(i);
        jvel[i+3] = omega_state(i);
        jpos_est[i+3] = eul_state(i);
        jvel_est[i+3] = omega_state(i);
    }

    for(size_t i=6; i<18; ++i){
        jpos[i] = jointPosTotal(i+1);
        jvel[i] = jointVelTotal(i);
        jpos_est[i] = jointPosTotal(i+1);
        jvel_est[i] = jointVelTotal(i);
    }


    const int* contactMat = loco_obj->getConDes();
    // if(controlTick>2499){
    //     if(controlTick < switchtime*ctrlHz){
    //         loco_obj->getStateEstimatefull(jpos_est,jvel_est,contactMat,rotE,robotdown,false,controlTick);
    //     }else{
    //         loco_obj->getStateEstimatefull(jpos_est,jvel_est,contactMat,rotE,robotdown,true,controlTick);
    //     }
    // }else if(controlTick>0){
    //     kinestimator(jpos_est,jvel_est,contactMat,rotE,robotdown);
    // }

    float vel_temp[3] = {jvel_est[0],jvel_est[1],jvel_est[2]};
    // if(controlTick>=31000){
    //     discrete_butter_f(filt,vel_temp);
    // }
    file_est << controlTick << "," << jpos[0] << "," << jpos[1] << "," << jpos[2] << "," << jvel[0] << "," << jvel[1] << "," << jvel[2] << ","
         << jpos[3] << "," << jpos[4] << "," << jpos[5] << "," << jvel[3] << "," << jvel[4] << "," << jvel[5] << ","
         << jpos_est[0] << "," << jpos_est[1] << "," << jpos_est[2] << "," << jvel_est[0] << "," << jvel_est[1] << "," << jvel_est[2] << ","
         << jpos_est[3] << "," << jpos_est[4] << "," << jpos_est[5] << "," << jvel_est[3] << "," << jvel_est[4] << "," << jvel_est[5] << ","
        // << imu_eul(0) << "," << imu_eul(1) << "," << imu_eul(2) << "," << imu_omega(0) << "," << imu_omega(1) << "," << imu_omega(2) << ","
         << rotE(0,0) << "," << rotE(0,1) << "," << rotE(0,2) << "," 
         << rotE(1,0) << "," << rotE(1,1) << "," << rotE(1,2) << ","
         << rotE(2,0) << "," << rotE(2,1) << "," << rotE(2,2) << ","
         << vel_temp[0] << "," << vel_temp[1] << "," << vel_temp[2] << "\n";
            
    // if(controlTick>=31000){
    //     jvel_est[0] = vel_temp[0];
    //     jvel_est[1] = vel_temp[1];
    //     jvel_est[2] = vel_temp[2];
    // }
    
    // for (size_t i = 0; i < 18; i++)
    // {
    //     jpos_est[i] = jpos[i];
    //     jvel_est[i] = jvel[i];
    // }
    
    Eigen::Matrix<double,16,1> q0;
    q0.setZero(16,1);
    q0.block(0,0,3,1) << jpos_est[0],jpos_est[1],jpos_est[2];//= jointPosTotal.block(0,0,3,1);
    q0.block(3,0,3,1) << jvel_est[0],jvel_est[1],jvel_est[2];//= jointVelTotal.block(0,0,3,1);
    q0.block(6,0,3,1) << jpos_est[3],jpos_est[4],jpos_est[5];
    q0.block(9,0,3,1) << jvel_est[3],jvel_est[4],jvel_est[5];//.block(3,0,3,1);
    std::map<std::string, casadi::DM> arg, res;
    
    int force[4] = {0};
    for(auto &con: A1.back()->getContacts()){
        int conInd = con.getlocalBodyIndex();
        force[conInd/3-1] = 500;
        force[conInd/3-1] = con.getNormal().e().norm();
    }

    // float vel_temp[3] = {0,0,0};//{cmd.vel[0],cmd.vel[1],cmd.vel[2]};
    // float pose_temp[6] = {0,0,0,0,0,0};
    // float filt_vel_temp[3] = {jvel_est[0],jvel_est[1],jvel_est[2]};
    // discrete_butter_f(filt,filt_vel_temp);
    int duration_data =0;
    Eigen::Matrix<double,4,1> nextcon = Eigen::MatrixXd::Ones(4,1);
    int maxsteps = 14;
    int settlingsteps = 4;
    double rearweight = 6;
    

    Eigen::Matrix<double, 12, 1> lastQPforce = Eigen::MatrixXd::Zero(12,1);
    
    /////////////////////////////////////////////////////////////////////
    //////////////////////////// CONTROL
    /////////////////////////////////////////////////////////////////////
    // Update the desired torques
    if(controlTick < settling){ // Settle down
        // loco_obj->posSetup(jointPosTotal.head(3));
        double temp[18] = {0};
        tau = temp;
        loco_obj->initStandVars(jointPosTotal.block(0,0,3,1),jointPosTotal(5),(int)duration);
    }
    else if(controlTick >= settling & controlTick < loco_start){ // Start standing
        loco_obj->calcTau2(jpos_est,jvel_est,rotMatrixDouble,STAND,controlTick,loco_start);//,shifttime,movetime,shifttime2,movetime2,movetime3);
        tau = loco_obj->getTorque();

    }
    else if(controlTick >= loco_start & controlTick < loco_start+shifttime){
        
        nextcon(0) = 0;
        
        if(controlTick==loco_start){
            Eigen::Matrix<double, 4, 1> wfoot = rearweight*Eigen::MatrixXd::Ones(4,1);
            wfoot(0)=1;//nextcon(0);
            wfoot(1)=1;//nextcon(1);
            loco_obj->getshiftedCoM(wfoot);
            loco_obj->setshiftedCoM();
        }

        loco_obj->setswingContact(nextcon);
        loco_obj->calcTau2(jpos_est,jvel_est,rotMatrixDouble,STANDUP,controlTick,loco_start);//,shifttime,movetime,shifttime2,movetime2,movetime3);
        tau = loco_obj->getTorque();

    }else if(controlTick >= loco_start+shifttime & controlTick < switchtime*ctrlHz){// & controlTick < loco_start + shifttime){ // Start locomotion
        
        int stepind = std::floor((controlTick-loco_start-shifttime)/(shifttime+movetime));
        //loco_obj->tookfirststep();
        if(stepind<maxsteps){
            loco_obj->stepsonwall(stepind);
        
            if(stepind==1){
                loco_obj->tookfirststep();
            }
        
            if(stepind%2==0){
                nextcon(0) = 0;
                if(controlTick == loco_start+shifttime + stepind*(shifttime+movetime)){
                    loco_obj->incstep();
                }         
            }else{
                nextcon(1) = 0;
            }
         
            if(controlTick==loco_start + (stepind+1)*(movetime + shifttime)){
                Eigen::Matrix<double, 4, 1> wfoot = rearweight*Eigen::MatrixXd::Ones(4,1);//3
                wfoot(0)=1;
                wfoot(1)=1;
                
                loco_obj->getshiftedCoM(wfoot);
                loco_obj->setshiftedCoM();
                
            }
            
            loco_obj->setswingContact(nextcon);
            loco_obj->calcTau2(jpos_est,jvel_est,rotMatrixDouble,STANDUP,controlTick,loco_start);//,shifttime,movetime,shifttime2,movetime2,movetime3);
            tau = loco_obj->getTorque();
        
        
        }else{
  
            int stepind2 = std::floor((controlTick-loco_start-shifttime-(maxsteps)*(movetime + shifttime))/(shifttime2+movetime2));
            loco_obj->settlesteps(stepind2);
            if(stepind2==settlingsteps){
                loco_obj->gotfinalstate();
            }
            if(controlTick==loco_start + shifttime + (maxsteps)*(movetime + shifttime)){// + stepind*(shifttime+movetime+shifttime2)){
                loco_obj->stopclimbing();
                loco_obj->stepsonwall(maxsteps);
                loco_obj->setfinalCoM();
            }
            
            loco_obj->calcTau2(jpos_est,jvel_est,rotMatrixDouble,STANDUP,controlTick,loco_start);//,shifttime,movetime,shifttime2,movetime2,movetime3);
            tau = loco_obj->getTorque();
        }
        
    }else if(controlTick >= switchtime*ctrlHz){

        int stepind2 = std::floor((controlTick-switchtime*ctrlHz)/(shifttime2+movetime3));
        loco_obj->settlesteps(stepind2);

        if(controlTick == switchtime*ctrlHz+2*(shifttime2+movetime3)+shifttime2){
            loco_obj->startwalking();
        }

        if(controlTick==30000 || controlTick==31000){
            loco_obj->readytoreallywalk();
        }

        if(controlTick==32000){
            loco_plan->letsgo();
        }
        
        //Eigen::Matrix<double,12,1> q_est = Eigen::MatrixXd::Zero(12,1);
        if(controlTick == switchtime*ctrlHz){
            //q_est.block(0,0,3,1) = q0.block(0,0,3,1);
            loco_obj->readytowalk();
            lastQPforce = loco_obj->getpreviousQPforce();
        }//else{
         //   q_est = loco_obj->getStateEstimate(jpos,jointVelTotal,imu_eul,imu_omega);
        //}

        if(controlTick == switchtime*ctrlHz + stepind2*(shifttime2+movetime3)){
            //q_est.block(0,0,3,1) = q0.block(0,0,3,1);
            loco_obj->setfinalCoM2(stepind2);
        }

        // q_est.block(6,0,3,1) = imu_eul;
        // q_est(9) = imu_w(0);
        // q_est(10) = imu_w(1);
        // q_est(11) = imu_w(2);
        // for(size_t i=0;i<3;i++){
        //     for (size_t j = 0; j < 3; j++)
        //     {
        //         rotMatrixDouble[3*i+j] = rotIMU(j,i);
        //     }
        // }

        // for(size_t i=0; i<3; ++i){
        //     //jpos[i] = q_est(i);
        //     //jvel[i] = q_est(3+i);
        //     jpos[3+i] = q_est(6+i);
        //     jvel[3+i] = q_est(9+i);
        // }

        //Eigen::Matrix<double,16,1> q0;
        //q0.setZero(16,1);
        // q0.block(0,0,3,1) << jpos[0],jpos[1],jpos[2];//= jointPosTotal.block(0,0,3,1);
        // q0.block(3,0,3,1) << jvel[0],jvel[1],jvel[2];//= jointVelTotal.block(0,0,3,1);
        // q0.block(6,0,3,1) << jpos[3],jpos[4],jpos[5];
        // q0.block(9,0,3,1) << jvel[3],jvel[4],jvel[5];//.block(3,0,3,1);

        loco_obj->updatestate(jpos_est,jvel_est,rotMatrixDouble);
        foot_position = loco_obj->getfootposition();
        hip_position = loco_obj->gethipposition();
        int stancephase = loco_obj->stancecounter();
        if(controlTick%10==0){
            plannerNMPC(controlTick, loco_obj, loco_plan, solver, q0, foot_position, lastQPforce);
        }
        //loco_obj->setRaisimD(Dr);
        //loco_obj->setRaisimH(Hr);
        loco_obj->calcTau2(jpos_est,jvel_est,rotMatrixDouble,UPWALK,controlTick,loco_start);//,shifttime,movetime,shifttime2,movetime2,movetime3);
        
        tau = loco_obj->getTorque();
    }

    
    jointTorqueFF = Eigen::Map< Eigen::Matrix<double,18,1> >(tau,18);
    jointTorqueFF.block(0,0,6,1).setZero();
    
    // Set the desired torques
    A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
    A1.back()->setGeneralizedForce(jointTorqueFF);

    // Apply the disturbance force to the desired body
    // Assuming the body index is 0, adjust as needed
    // if(controlTick > 500 && controlTick < 1000){

    //     double amp = 600.0;
    //     //double dist = distx(controlTick, amp);
    //     Eigen::Vector3d disturbanceForce(0,amp,0); 
    //     A1.back()->setExternalForce(0, disturbanceForce);
    // }

};



int main(int argc, char *argv[]) {
    
    long int terrain_number = atoll(argv[3]); // Use atoll for long long int

    // ============================================================ //
    // =================== SETUP RAISIM/VISUALS =================== //
    // ============================================================ //
    /// create raisim world
    raisim::World::setActivationKey(raisim::loadResource("activation.raisim"));
    raisim::World world;
    world.setTimeStep(simfreq_raisim);
    
    raisim::OgreVis *vis = raisim::OgreVis::get();

    /// these method must be called before initApp
    vis->setWorld(&world);
    vis->setWindowSize(1792, 1200); // Should be evenly divisible by 16!!
    vis->setImguiSetupCallback(imguiSetupCallback); // These 2 lines make the interactable gui visible
    vis->setImguiRenderCallback(imguiRenderCallBack);
    vis->setKeyboardCallback(raisimKeyboardCallback);
    vis->setSetUpCallback(setupCallback);
    vis->setAntiAliasing(2);

    /// starts visualizer thread
    vis->initApp();
    
    /// create raisim objects
    raisim::TerrainProperties terrainProperties;
    terrainProperties.frequency = 0.0;
    terrainProperties.zScale = 0.0;
    terrainProperties.xSize = 300.0;
    terrainProperties.ySize = 300.0;
    terrainProperties.xSamples = 50;
    terrainProperties.ySamples = 50;
    terrainProperties.fractalOctaves = 0;
    terrainProperties.fractalLacunarity = 0.0;
    terrainProperties.fractalGain = 0.0;

    raisim::HeightMap *ground = world.addHeightMap(0.0, 0.0, terrainProperties);
    vis->createGraphicalObject(ground, "terrain", "checkerboard_blue");
    world.setDefaultMaterial(0.8, 0.0, 0.0); //surface friction could be 0.8 or 1.0
    vis->addVisualObject("extForceArrow", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP); 
    
    
    // // WEIGHT VISUALIZATION FOR LCSS PAPER, KEEP FOR NOW.
    // // create raisim objects
    // double scale = 1;
    // vis->loadMeshFile("/home/kavehakbarihamed/raisim/workspace/A1_LL_Exp/rsc/Dumbell_5lb.STL", "weight1", false);
    // raisim::VisualObject *weightVis1 = vis->addVisualObject("weight_vis1", "weight1", "purple", {scale, scale, scale}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    // vis->loadMeshFile("/home/kavehakbarihamed/raisim/workspace/A1_LL_Exp/rsc/Dumbell_5lb.STL", "weight2", false);
    // raisim::VisualObject *weightVis2 = vis->addVisualObject("weight_vis2", "weight2", "purple", {scale, scale, scale}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);

    // raisim::Box *box = world.addBox(0.24,0.08,0.15,4.54);
    // box->setPosition(-0.005,0,0.25);
    // // vis->createGraphicalObject(box,"Payload","purple");

    // raisim::Box *box = world.addBox(0.04,0.04,0.01,20);
    // box->setPosition(0.,0,0.0);
    // // box->setPosition(-0.1,0,0.25);
    // vis->createGraphicalObject(box,"Payload","purple"); 

    auto& list = vis->getVisualObjectList();

    // ============================================================ //
    // ======================= SETUP Robot ======================== //
    // ============================================================ //
    std::vector<raisim::ArticulatedSystem*> A1;
    // A1.push_back(world.addArticulatedSystem(raisim::loadResource("Go1/Go1.urdf"))); // WHEN USING Go1, BE SURE TO CHANGE CMAKE TO USE CORRECT DYNAMICS
    A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_modified_new.urdf")));
    //A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_modified_up_mod_sensored.urdf")));
    
    
    vis->createGraphicalObject(A1.back(), "A1");

    //Rear offset -0.1
    //A1.back()->setGeneralizedCoordinate({0, 0, 0.5, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                  -0.7337, 1.0175, -2.035, 0.7337, 1.0175, -2.035, 0.0, 2.4532, -1.1582, 0.0, 2.4532, -1.1582});
    
    A1.back()->setGeneralizedCoordinate({0, 0, 0.12, 1, 0, 0, 0,
                                        0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6});   
    //Cutting to the transition point 5.10208e-05
    //A1.back()->setGeneralizedCoordinate({-0.0552767,0,0.496916,0.8257,-0.0079,-0.5640,-0.0080,//-0.0119154,-1.19863,-0.0108896,
    //                                        -0.535698,0.818817,-2.2517,0.511951,0.790705,-2.24335,-0.00892291,1.93392,-1.36839,0.00366515,1.94833,-1.3735});
    //A1.back()->setGeneralizedCoordinate({-0.0262742,0.0139459,0.493275,0.8203,-0.1058,-0.5509,-0.1111,//-0.176766,-1.1874,-0.136993,
    //                                        -0.654018,0.886213,-2.23036,0.542235,0.870764,-2.55789,0.0609142,1.99266,-1.34162,0.059542,2.05143,-1.40771});
    A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
    A1.back()->setName("A1_Robot");
    
    raisim::Box *box_right = world.addBox(200.0, 0.2, 0.8, 1000000, "rubber");//terrainProperties);
    raisim::Box *box_left = world.addBox(200.0, 0.2, 0.8, 1000000, "rubber");

    box_right->setPosition(0,-0.32,0.4);
    box_left->setPosition(0,0.32,0.4);

    //vis->createGraphicalObject(box_right, "right_wall", "checkerboard_blue");
    vis->createGraphicalObject(box_left, "left_wall", "checkerboard_blue");
    
    A1.back()->getCollisionBody("FR_foot/0").setMaterial("wood");
    A1.back()->getCollisionBody("FL_foot/0").setMaterial("wood");

    world.setMaterialPairProp("wood", "rubber", 0.8, 0, 0);

    int roughterrain = 0;
    if(roughterrain){
       //long int randomSeed = std::time(nullptr);
       std::srand(terrain_number);
       //std::srand(163727240);
       //std::srand(1642619542);
       //std::cout << "RNG Seed: " << randomSeed << std::endl;
       bool GroundHeightVariation = true;
       if(GroundHeightVariation){
           int numBlk = 20;//150;
           double percent = 55;//70 // there will be a block x percent of the time
           int direction = 0; // 0 for x, 1 for y
           int fwd_bwd = 1; // 1 for forward, -1 for backward
           int maxHeight = 2;//4; // max height in centimeters
           double startPos = fwd_bwd*0.5;//fwd_bwd*2;
           double height, width, length;
           double left, right;
           const char * colors[3] ={"yellow","orange","red"};
           int j = 0;
           width = 0.14; length = 0.8;
           double mass = 0.9;
           for (int i=0; i<numBlk; i++){
               left = rand() % 100;
               right = rand() % 100;
               if (left>(100-percent)){
                   height = 1.0*(rand() % (maxHeight+1))/100.0;
                   j = (height<=(1.0*(maxHeight)/300.0)) ? 0 : (height<=(2.0*(maxHeight)/300.0)) ? 1 : 2;
                   if(direction == 0){
                       raisim::Box *box = world.addBox(width,length,height,mass);
                       box->setPosition(startPos,length/2,0.01);
                       vis->createGraphicalObject(box,"box"+std::to_string(i)+"left",colors[j]);
                   }else{
                       raisim::Box *box = world.addBox(length,width,height,mass);
                       box->setPosition(length/2,startPos,0.01);
                       vis->createGraphicalObject(box,"box"+std::to_string(i)+"left",colors[j]);
                   }
               }
               if (right>(100-percent)){
                   height = 1.0*(rand() % (maxHeight+1))/100.0;
                   j = (height<=(1.0*(maxHeight)/300.0)) ? 0 : (height<=(2.0*(maxHeight)/300.0)) ? 1 : 2;
                   if(direction == 0){
                       raisim::Box *box = world.addBox(width,length,height,mass);
                       box->setPosition(startPos,-length/2,0.01);
                       vis->createGraphicalObject(box,"box"+std::to_string(i)+"right",colors[j]);
                   }else{
                       raisim::Box *box = world.addBox(length,width,height,mass);
                       box->setPosition(-length/2,startPos,0.01);
                       vis->createGraphicalObject(box,"box"+std::to_string(i)+"right",colors[j]);
                   }
               }
               startPos+=fwd_bwd*width;
           }
           if(maxHeight ==2 && numBlk != 150){
               // stack blocks on top
               startPos = fwd_bwd*2;
               for (int i=numBlk+0; i<numBlk+numBlk; i++){
                   left = rand() % 100;
                   right = rand() % 100;
                   if (left>(100-percent)){
                       height = 1.0*(rand() % (maxHeight+1))/100.0;
                       j = (height<=(1.0*(maxHeight)/300.0)) ? 0 : (height<=(2.0*(maxHeight)/300.0)) ? 1 : 2;
                       if(direction == 0){
                           //raisim::Box *box = world.addBox(width,length,height,mass);
                           //box->setPosition(startPos,length/2,0.01);
                           raisim::Box *box = world.addBox(width*((rand()%200)*0.001+1),length*((rand()%300)*0.001+1),height,mass);
                           box->setPosition(startPos*((rand()%5)*0.1+1),length/2,0.1);
                           vis->createGraphicalObject(box,"box"+std::to_string(i)+"left",colors[j]);
                       }else{
                           raisim::Box *box = world.addBox(length,width,height,mass);
                           box->setPosition(length/2,startPos,0.01);
                           vis->createGraphicalObject(box,"box"+std::to_string(i)+"left",colors[j]);
                       }
                   }
                   if (right>(100-percent)){
                       height = 1.0*(rand() % (maxHeight+1))/100.0;
                       j = (height<=(1.0*(maxHeight)/300.0)) ? 0 : (height<=(2.0*(maxHeight)/300.0)) ? 1 : 2;
                       if(direction == 0){
                           //raisim::Box *box = world.addBox(width,length,height,mass);
                           //box->setPosition(startPos,-length/2,0.01);
                           raisim::Box *box = world.addBox(width*((rand()%200)*0.001+1),length*((rand()%300)*0.001+1),height,mass);
                           box->setPosition(startPos*((rand()%5)*0.1+1),-length/2,0.1);
                           vis->createGraphicalObject(box,"box"+std::to_string(i)+"right",colors[j]);
                       }else{
                           raisim::Box *box = world.addBox(length,width,height,mass);
                           box->setPosition(-length/2,startPos,0.01);
                           vis->createGraphicalObject(box,"box"+std::to_string(i)+"right",colors[j]);
                       }
                   }
                   startPos+=fwd_bwd*width;
               }
           }
       }
    }
    LocoWrapper* loco_obj1 = new LocoWrapper(argc,argv);
    //LocoWrapperwalk* loco_obj2 = new LocoWrapperwalk(argc,argv);
    //loco_obj1->setRFfalse();
    SRBNMPC* loco_plan = new SRBNMPC(argc,argv,1,0);
    loco_plan->generator();
    std::string file_name = "Go2_w0p2_ro_0p05_35b";//"take2_1";
    
    std::string prefix_code = std::filesystem::current_path().string() + "/";
    std::string prefix_lib = std::filesystem::current_path().string() + "/";

    // Create a new NLP solver instance from the compiled code
    std::string lib_name = prefix_lib + file_name + ".so";
    casadi::Dict opts = {{"ipopt.print_level", 1}, {"print_time", 0},{"ipopt.max_iter", 10},{"ipopt.acceptable_tol", 1e-2},{"ipopt.acceptable_obj_change_tol", 1e-2}};
    casadi::Function solver = casadi::nlpsol("solver", "ipopt", lib_name, opts);

    float a[3] = {1.00000000, -1.99555712, 0.99556697};
    float b[3] = {0.00000246, 0.00000492, 0.00000246};
    // float a[3] = {1.0, -1.14298050253990, 0.41280159809619};
    // float b[3] = {0.06745527388907, 0.13491054777814, 0.06745527388907};
    populate_filter_f(filt, a, b, 3, 3);


    // ============================================================ //
    // ================= VIEW AND RECORDING OPTIONS =============== //
    // ============================================================ //
    raisim::gui::showContacts = false;
    raisim::gui::showForces = false;
    raisim::gui::showCollision = false;
    raisim::gui::showBodies = true;

    std::string cameraview = "side";
    bool panX = true;                // Pan view with robot during walking (X direction)
    bool panY = false;                // Pan view with robot during walking (Y direction)
    bool record = true;             // Record?
    double startTime = 0*ctrlHz;    // Recording start time
    double simlength = 50000;//60000;//300*ctrlHz;   // Sim end time
    double fps = 30;            
    //std::string directory = "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/Oct10/";
    std::string directory = "../data25/Feb28/";
    // std::string filename = "Payload_Inplace";
    std::string filename = "fullsim";//"JacVCL_OWCL_rt55_3";
    // std::string filename = "inplace_sim";


    // ============================================================ //
    // ========================= VIEW SETUP ======================= //
    // ============================================================ //
    // NOTE: Pi is defined in /dynamics/dynamicsSupportFunctions.h
    // NOTE: This section still needs some work. 
    double Pi=3.14;
    if(cameraview == "iso"){
        // vis->getCameraMan()->getCamera()->setPosition(3, 3, 2);
        vis->getCameraMan()->getCamera()->setPosition(2, -1, 0.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(5*Pi/6-Pi/2));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else if(cameraview == "isoside"){
        vis->getCameraMan()->getCamera()->setPosition(1.1, -2, 0.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(4*Pi/6-Pi/2));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else if(cameraview == "side"){
        vis->getCameraMan()->getCamera()->setPosition(0, -2, 0.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(0));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else if(cameraview == "front"){
        vis->getCameraMan()->getCamera()->setPosition(2, 0, 0.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(Pi/2));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else if(cameraview == "top"){
        //vis->getCameraMan()->getCamera()->setPosition(1, -3, 2.5);
        //vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(1.0));
        vis->getCameraMan()->getCamera()->setPosition(1, 0, 2.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(Pi/2));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else{
        vis->getCameraMan()->getCamera()->setPosition(1, -3, 2.5);
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(1.0));
    }
    unsigned long mask = 0;
    if(raisim::gui::showBodies) mask |= raisim::OgreVis::RAISIM_OBJECT_GROUP;
    if(raisim::gui::showCollision) mask |= raisim::OgreVis::RAISIM_COLLISION_BODY_GROUP;
    if(raisim::gui::showContacts) mask |= raisim::OgreVis::RAISIM_CONTACT_POINT_GROUP;
    if(raisim::gui::showForces) mask |= raisim::OgreVis::RAISIM_CONTACT_FORCE_GROUP;
    vis->setVisibilityMask(mask);

    //if(panX) raisim::gui::panViewX = panX;
    //if(panY) raisim::gui::panViewY = panY;
    
    // ============================================================ //
    // ========================== RUN SIM ========================= //
    // ============================================================ //
    const std::string name = directory+filename+"_"+cameraview+".mp4";
    vis->setDesiredFPS(fps);
    long simcounter = 0;//22350;
    static bool added = false;

    raisim::Contact contactInstance;

    std::ofstream file_est("../data25/estimator60.csv");

    while (!vis->getRoot()->endRenderingQueued() && simcounter <= simlength){

        size_t dist_start = 40000*ctrlHz;               // Start the disturbance (if any)
        size_t dist_stop  = dist_start+200;             // Stop the disturbance (if any)

        // FOR LCSS PAPER, KEEP FOR NOW.
        // if (simcounter>1500 & !added){
        //     box->setPosition(-0.1,0,0.30);
        //     added = true;
        // }
        // raisim::Vec<3> pos = box->getComPosition();
        // weightVis1->offset = {pos[0]-0.245/2,pos[1]-.078/2,pos[2]-0.078};
        // weightVis2->offset = {pos[0]-0.245/2,pos[1]-.078/2,pos[2]-0.015};
        // // std::cout<<pos[0]<<"\t"<<pos[1]<<"\t"<<pos[2]<<std::endl;

        controller(A1,loco_obj1,loco_plan,solver,simcounter,contactInstance,file_est);
        world.integrate();        
        
        if (simcounter%15 == 0)
            vis->renderOneFrame();
        
        if (!vis->isRecording() & record & simcounter>=startTime)
            vis->startRecordingVideo(name);
        
        auto currentPos = vis->getCameraMan()->getCamera()->getPosition();
        //if (raisim::gui::panViewX){
            Eigen::VectorXd jointPosTotal(18 + 1);
            Eigen::VectorXd jointVelTotal(18);
            jointPosTotal.setZero();
            jointVelTotal.setZero();
            A1.back()->getState(jointPosTotal, jointVelTotal);
            if (cameraview=="front"){
                currentPos[0] = jointPosTotal(0)+2;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="side"){
                currentPos[0] = jointPosTotal(0);
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="iso"){
                currentPos[0] = jointPosTotal(0)+2;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="isoside"){
                currentPos[0] = jointPosTotal(0)+1.1;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            }
        //}
        /*if (raisim::gui::panViewY){
            Eigen::VectorXd jointPosTotal(18 + 1);
            Eigen::VectorXd jointVelTotal(18);
            jointPosTotal.setZero();
            jointVelTotal.setZero();
            A1.back()->getState(jointPosTotal, jointVelTotal);
            if (cameraview=="front"){
                currentPos[1] = jointPosTotal(1);
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="side"){
                currentPos[1] = jointPosTotal(1)-2;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="iso"){
                currentPos[1] = jointPosTotal(1)-1;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="isoside"){
                currentPos[1] = jointPosTotal(1)-2;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            }
        }*/
        std::cout << "simcounter" << "\t" << simcounter << std::endl;

        // if(abs(jointPosTotal(1))>0.04){
        //     std::cout << simcounter << "\t" << jointPosTotal(0) << "\t" << jointPosTotal(1) << "\t" << jointPosTotal(2) << "\t" << jointPosTotal(15) << "\t" << jointPosTotal(18) << "\t" << -1 << std::endl;
        //     break;
        // }else if(abs(jointPosTotal(2)-0.5)>0.05){
        //     std::cout << simcounter << "\t" << jointPosTotal(0) << "\t" << jointPosTotal(1) << "\t" << jointPosTotal(2) << "\t" << jointPosTotal(15) << "\t" << jointPosTotal(18) << "\t" << -2 << std::endl;
        //     break;
        // }else if(jointPosTotal(15) > 0 || jointPosTotal(18) > 0){
        //     std::cout << simcounter << "\t" << jointPosTotal(0) << "\t" << jointPosTotal(1) << "\t" << jointPosTotal(2) << "\t" << jointPosTotal(15) << "\t" << jointPosTotal(18) << "\t" << -3 << std::endl;
        //     break;
        // }else if(jointPosTotal(0)>7){
        //     std::cout << simcounter << "\t" << jointPosTotal(0) << "\t" << jointPosTotal(1) << "\t" << jointPosTotal(2) << "\t" << jointPosTotal(15) << "\t" << jointPosTotal(18) << "\t" << 1 << std::endl;
        //     break;
        // }
        
        simcounter++;
        
    }
    file_est.close();

    // End recording if still recording
    if (vis->isRecording())
        vis->stopRecordingVideoAndSave();

    /// terminate the app
    vis->closeApp();

    delete loco_obj1;
    //delete loco_obj2;
    clear_filter_f(filt);

    return 0;
}

