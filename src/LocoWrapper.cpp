//
// Authror: Randy Fawcett on 12/2021.
//
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "LocoWrapper.hpp"
#include "iostream"

LocoWrapper::LocoWrapper(int argc, char *argv[]) : Parameters(argc,argv){

//    std::string filename = "/media/kavehakbarihamed/Data/A1_RaiSim_Outputs/LCSS_2021/Payload_Trot_10cm.txt";
//    std::string filename = "/media/kavehakbarihamed/Data/A1_RaiSim_Outputs/nothing.txt";
    std::string filename = "../data25/nothing.csv";
//    std::string filename = ""; // empty string will produce no output file
    
    data = std::unique_ptr<DataLog>( new DataLog(filename) ); // make_unique DNE in c++11
    quad = new RobotModel();
    conEst = new ContactEst();
    LL = new LowLevelCtrl();
    VC = new VirtualConstraints();
    PP = new MotionPlanner();

    state = quad->getStatePointer();
    dyn = quad->getDynamicsPointer();
    kin = quad->getKinematicsPointer();
    con = conEst->getConInfoPointer();
    traj = PP->getTrajInfoPointer();
    vcon = VC->getVCPointer();
    ll = LL->getllPointer();
    
    contact_horizon.block(1,25,1,15) = Eigen::MatrixXd::Zero(1,15);//15
    contact_horizon.block(2,28,1,12) = Eigen::MatrixXd::Zero(1,12);
    contact_horizon.block(0,5,1,15) = Eigen::MatrixXd::Zero(1,15);//15
    contact_horizon.block(3,8,1,12) = Eigen::MatrixXd::Zero(1,12);
}

LocoWrapper::~LocoWrapper(){
    delete quad;
    delete conEst;
    delete LL;
    delete VC;
    delete PP;
}

void LocoWrapper::calcTau(const double q[18], const double dq[18], const double R[9], const int force[4], size_t gait, size_t ctrlTick, size_t duration){
    
    //phaseVar = getPhase(1.0*locoTick, 0.0, 1.0*traj->domLen);   // update phase variable
    quad->updateState(q,dq,R);                                  // update state
    

    // float footPos[4] = {0};  // DUMMY VARS
    // if (gait!=gaitTemp || (phaseVar>maxPhase && gait!=STAND) ){ 
    //     // Change domain immediately since gait changed
    //     conEst->forceDomChange();                                                       // force con->changeDomain=1 to plan properly
    //     // std::cout << "Time trigger: " << phaseVar << std::endl;
    //     phaseVar = 0;
    //     PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);     // plan trajectory
    //     conEst->updateConState(footPos,phaseVar,force);                                 // update contact detection
    //     locoTick = 0;
    // }else {
    //     // Wait for impact to change domain
    //     conEst->updateConState(footPos,phaseVar,force);                                 // impact detection
    //     if (con->changeDomain==1 && gait!=STAND){
    //         locoTick = 0;
    //         // std::cout << "Contact trigger: " << phaseVar << std::endl;
    //         phaseVar = 0;
    //     }
    //     PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params,  opt_HLstate, NLstep);     // plan trajectory
    // }

    float footPos[4] = {0};  // DUMMY VARS
    if (ctrlTick<1 || gait!=gaitTemp || (phaseVar>maxPhase && gait!=STAND) ){ 
        locoTick = 0;
        phaseVar = getPhase(1.0*locoTick, 0.0, 199);
        flphase = 0;
        rlphase = 0;
        //conEst->forceDomChange();
        //std::cout << ctrlTick << "\t" << locoTick << "\t" << phaseVar << "\t" << "Yes1" << std::endl;
        //traj->toeInit = kin->toePos;
            
    }else{
        phaseVar = getPhase(1.0*locoTick, 0.0, 199);
        flphase = getPhase(1.0*locoTick-49, 0.0, 145);
        rlphase = getPhase(1.0*locoTick-79, 0.0, 115);
        flphase = (flphase<0)?0:(flphase>1)?1:flphase;
        rlphase = (rlphase<0)?0:(rlphase>1)?1:rlphase;
        if(locoTick==50){//||locoTick==80){
            conEst->forceDomChange();
        }
        //std::cout << ctrlTick << "\t" << locoTick << "\t" << phaseVar << "\t" << "Yes2" << std::endl;           
    }
    
    PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);
    bool reachedwall = PP->getReachedWall();
    quad->updateSwingMatrices(con->ind,con->cnt);                                             // update the jacobians
    if(gait==STANDUP){
        VC->updateVirtualConstraintssetfoot(state, kin, traj, con, gait, phaseVar, &motion_params, ll, flphase,rlphase,reachedwall);    // update VC's
        if(reachedwall){
            LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
        }else{  
            LL->calcTorque(state, dyn, kin, vcon, con, &ll_params);    
        } 
    }else{
        VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
        LL->calcTorque(state, dyn, kin, vcon, con, &ll_params);
    }
                                    // run low level controller
    //data->writeData(state,vcon,traj,ll,ctrlTick,force);                                        // log relavent data
    plottingfoothd(vcon,con);
    data->writeData(state,vcon,con,traj,ll,kin,ctrlTick,force,opt_HLstate,locoTick,phaseVar,flphase,rlphase,0.0,0.0,NLstep,duration);
    locoTick += (ctrlHz)/LL_Hz;     // increment locoTick
    gaitTemp = gait;                // update the previous gait used
}

void LocoWrapper::calcTau2(const double q[18], const double dq[18], const double R[9], size_t gait, size_t ctrlTick, size_t solveduration){

    phaseVar = getPhase(1.0*locoTick, 0.0, 1.0*traj->domLen);   // update phase variable
    phaseVar = (phaseVar>1) ? 1 : ((phaseVar<0) ? 0 : phaseVar);

    quad->updateState(q,dq,R); 
    
    if(gait==STAND){
        //conEst->updateConState(footPos,phaseVar,force);
        //phaseVar = getPhase(1.0*locoTick, 0.0, 1.0*traj->domLen);
        
        PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);
        quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
        VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
        LL->calcTorque(state, dyn, kin, vcon, con, &ll_params);   
                                         // run low level controller                                       // log relavent data
    }else if(gait == UPWALK){

        
        if(readytowalkf){
        
            // if (gait!=gaitTemp || (phaseVar>maxPhase) || ctrlTick == switchtime*ctrlHz + 2*(shifttime2+movetime3) + shifttime2){ 
            if (HLMTphase == 0 || HLMTphase == 20 || ctrlTick == switchtime*ctrlHz + 2*(shifttime2+movetime3) + shifttime2){ 

                locoTick = 0;
                phaseVar = getPhase(1.0*locoTick, 0.0, 199);
                flphase = 0;
                rlphase = 0;
                //conEst->forceDomChange();
                z = Eigen::MatrixXd::Zero(12,1);
            
            }else{
                phaseVar = getPhase(1.0*locoTick, 0.0, 199);
                flphase = getPhase(1.0*locoTick-49, 0.0, 145);
                rlphase = getPhase(1.0*locoTick-79, 0.0, 115);
                flphase = (flphase<0)?0:(flphase>1)?1:flphase;
                rlphase = (rlphase<0)?0:(rlphase>1)?1:rlphase;
                // if(locoTick==50||locoTick==80){
                if (HLMTphase == 5 || HLMTphase == 8 || HLMTphase == 25 || HLMTphase == 28){
                    conEst->forceDomChange();
                }            
            }
   
            PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);  
            VC->updateVirtualConstraintswalk(state, kin, traj, con, gait, flphase,rlphase, &motion_params, ll);    // update VC's    
            VC->setDesiredForce(opt_HLstate.block(12,0,12,1));
            z.block(6,0,3*(4-con->cnt),0) += vcon->y.block(6,0,3*(4-con->cnt),0)/ctrlHz;         
            LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, HRai, z, Ki);                                     // run low level controller
        
        }else{
            
            if( ctrlTick == switchtime*ctrlHz+settlestep*(shifttime2+movetime3) || ctrlTick == switchtime*ctrlHz+settlestep*(shifttime2+movetime3)+shifttime2){
                locoTick = 0;
                phaseVar = 0;
                PP->setToeInit(kin);
                PP->setx0y0z0(state->q(0),state->q(1),state->q(2),state->q(4));
            }

            if(ctrlTick < switchtime*ctrlHz + settlestep*(shifttime2+movetime3)+shifttime2){
                PP->shiftCoM3(conEst,phaseVar,shifttime2,true);
                quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
                VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
                VC->setDesiredForce(opt_HLstate.block(12,0,12,1));
                LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
            }else{
                if(settlestep%2==0){
                    nextContact[0] = 0;
                    nextContact[1] = 1;
                }else{
                    nextContact[0] = 1;
                    nextContact[1] = 0;
                }
                conEst->setDesDomain(nextContact);
                PP->movefoot3(movetime3);
                quad->updateSwingMatrices(con->ind,con->cnt);
                VC->updateVirtualConstraintssetfoot(state, kin, traj, con, gait, phaseVar, &motion_params, ll, flphase,rlphase,true);
                VC->setDesiredForce(opt_HLstate.block(12,0,12,1));
                LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
            }

        }


    }else{

        if(stopclimb){

            if( ctrlTick == loco_start_e+shifttime+(wallstep)*(shifttime+movetime)+settlestep*(shifttime2+movetime2) || ctrlTick == loco_start_e+shifttime+(wallstep)*(shifttime+movetime)+settlestep*(shifttime2+movetime2)+shifttime2){
                locoTick = 0;
                phaseVar = 0;
                PP->setToeInit(kin);
                PP->setx0y0z0(state->q(0),state->q(1),state->q(2),state->q(4));
            }
            
            if(atfinalstate){
                PP->shiftCoM2(conEst,phaseVar,shifttime2,true);
                quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
                VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
                LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
            }else{
                
                if(ctrlTick < loco_start_e+shifttime+(wallstep)*(shifttime+movetime) + settlestep*(shifttime2+movetime2)+shifttime2){
                    PP->shiftCoM2(conEst,phaseVar,shifttime2,true);
                    quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
                    VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
                    LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
                }else{
                    if(settlestep%2==0){
                        nextContact[0] = 0;
                        nextContact[1] = 1;
                    }else{
                        nextContact[0] = 1;
                        nextContact[1] = 0;
                    }
                    conEst->setDesDomain(nextContact);
                    setxzsteplength(movetime2);
                    //PP->movefoot(movetime,wallstep);
                    quad->updateSwingMatrices(con->ind,con->cnt);
                    VC->updateVirtualConstraintssetfoot(state, kin, traj, con, gait, phaseVar, &motion_params, ll, flphase,rlphase,true);
                    LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
                }
            }

        }else{
        
            if(gait!=gaitTemp || ctrlTick == loco_start_e+shifttime || ctrlTick == loco_start_e+(wallstep+1)*(shifttime+movetime) || ctrlTick == loco_start_e+shifttime+(wallstep)*(shifttime+movetime)){
                locoTick = 0;
                phaseVar = 0;
                PP->setToeInit(kin);
                PP->setx0y0z0(state->q(0),state->q(1),state->q(2),state->q(4));
                if(ctrlTick == loco_start_e+shifttime){
                    setrearhippose();
                }
            }

            if(ctrlTick < loco_start_e+shifttime){
            
                PP->shiftCoM(conEst,phaseVar,shifttime);
                quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
                VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
                LL->calcTorque(state, dyn, kin, vcon, con, &ll_params);
        
            }else if(ctrlTick < loco_start_e+(wallstep+1)*(shifttime+movetime)){

                conEst->setDesDomain(nextContact);
                //if(wallstep<1){
                PP->movefoot(movetime,wallstep);//,wallstep);//(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);
                //}else{
                //PP->movefoot2(movetime,phaseVar);
                //}
                quad->updateSwingMatrices(con->ind,con->cnt);
                VC->updateVirtualConstraintssetfoot(state, kin, traj, con, gait, phaseVar, &motion_params, ll, flphase,rlphase,true);
                LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
        
            }else{

                PP->shiftCoM2(conEst,phaseVar,shifttime,false);
                quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
                VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
                LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);  
            }
        }
    }
    
    //PP->datalogger(ctrlTick);
    data->writeData(state,vcon,con,traj,ll,kin,ctrlTick,force_LL,opt_HLstate,locoTick,phaseVar,flphase,rlphase,0.0,0.0,NLstep,solveduration);
    locoTick += (ctrlHz)/LL_Hz;     // increment locoTick
    gaitTemp = gait;

}

void LocoWrapper::setcontactconfig(int controlMPC){
    
    if(controlMPC<2740){
        desired_contact[0] = 1;
        desired_contact[1] = 1;
        desired_contact[2] = 1;
        desired_contact[3] = 1;
    }else{
        for(int i=0;i<4;i++){
            desired_contact[i]=contact_horizon(i,controlMPC%40); 
        }
    }
    conEst->setDesDomain(desired_contact);
    quad->updateSwingMatrices(con->ind,con->cnt); 
}

void LocoWrapper::plottingfoothd(const VCInfo *vc, const ContactInfo *con){

    NLstep(0) = (1-con->ind[0])*vc->hd(7);
    NLstep(1) = (1-con->ind[1])*vc->hd(7);
    NLstep(2) = (1-con->ind[2])*vc->hd(11);
    NLstep(3) = (1-con->ind[0])*vc->hd(11);
}

void LocoWrapper::getshiftedCoM(Eigen::Matrix<double, 4, 1> footweight){

    Eigen::Matrix<double, 3, 4> footpos = kin->toePos;
    CoMnew.block(0,0,2,1) = footweight(0)*footpos.block(0,0,2,1) + footweight(1)*footpos.block(0,1,2,1) + footweight(2)*footpos.block(0,2,2,1) + footweight(3)*footpos.block(0,3,2,1);
    // x_new = newCoM(0);
    // y_new = newCoM(1);
    // z_new = newCoM(2);
    CoMnew.block(0,0,2,1) = CoMnew.block(0,0,2,1)/footweight.sum();
    CoMnew(2)=0.25;//CoMnew(2);
    CoMnew(3)= -(std::floor(wallstep/2)+1)*0.25;
}

void LocoWrapper::setrearhippose(){

    double hip_x = (kin->hipPos(0,2)+kin->hipPos(0,3))/2;
    double hip_y = (kin->hipPos(1,2)+kin->hipPos(1,3))/2;
    double hip_z = (kin->hipPos(2,2)+kin->hipPos(2,3))/2;
    double foot_x = (kin->toePos(0,2)+kin->toePos(0,3))/2;
    PP->setrearhip(hip_x,hip_y,hip_z,foot_x);
}

void LocoWrapper::setfinalCoM(){
    
    Eigen::Matrix<double,4,1> CoM_final = Eigen::MatrixXd::Zero(4,1);
    CoM_final(0) = (kin->toePos(0,2)+kin->toePos(0,3))/2+0.1;
    CoM_final(2) = 0.5;
    CoM_final(3) = -1.2;
    PP->setshiftedCoM(CoM_final);

};

void LocoWrapper::setfinalCoM2(){
    
    double pitch = ((0.3-0.1*(settlestep+1))>0) ? (0.3-0.1*(settlestep+1)) : 0;

    double rhipz = (kin->hipPos(2,2) + kin->hipPos(2,3))/2;
    Eigen::Matrix<double,4,1> CoM_final = Eigen::MatrixXd::Zero(4,1);
    CoM_final(0) = (kin->toePos(0,2)+kin->toePos(0,3))/2+0.1;// + 0.05; //(kin->hipPos(0,2) + kin->hipPos(0,3))/2 + 0.183*sin(pitch);//
    CoM_final(2) = 0.5;//rhipz + 0.183*cos(pitch)-0.05;//0.5
    CoM_final(3) = pitch;///(phaseVar>1) ? 1 : ((phaseVar<0) ? 0 : phaseVar)
    PP->setshiftedCoM(CoM_final);

};

void LocoWrapper::setxzsteplength(size_t movetime){
                    
    Eigen::Matrix<double, 4, 1> xzsteplenth = Eigen::MatrixXd::Zero(4,1);
    xzsteplenth(0) = 0.1;//kin->HipPos(0,0)+0.1-kin->ToePos(0,0);
    xzsteplenth(1) = 0;//kin->HipPos(0,1)+0.1-kin->ToePos(0,1);
    xzsteplenth(2) = 0.1;//kin->HipPos(0,2)-kin->ToePos(0,2);
    xzsteplenth(3) = 0;//kin->HipPos(0,3)-kin->ToePos(0,3);
    PP->movefoot2(movetime,xzsteplenth);
}

void LocoWrapper::setoptNLstate(Eigen::Matrix<double, 33, 1> HLopt){
    
    opt_HLstate.block(0,0,12,1) = HLopt.block(0,0,12,1);
    //opt_HLstate(2) = 0.35;
    opt_HLstate.block(12,0,12,1) = HLopt.block(16,0,12,1);
    NLstep = HLopt.block(28,0,5,1);
}

Eigen::Matrix<double, 12, 1> LocoWrapper::getStateEstimate(double jointPos[18], Eigen::VectorXd jointVelTotal, Eigen::Matrix<double, 3, 1> imu_eul, Eigen::Matrix<double, 3, 1> imu_omega){
    
    Eigen::Matrix<double, 12, 1> p_est =Eigen::MatrixXd::Zero(12,1);

    double jointPosIMU[18] = {0};
    for (size_t i = 0; i < 18; i++)
    {
        jointPosIMU[i] = jointPos[i];
    }
    jointPosIMU[3] = imu_eul(0);
    jointPosIMU[4] = imu_eul(1);
    jointPosIMU[5] = imu_eul(2);

    Eigen::VectorXd jointVelTotalIMU = jointVelTotal;
    jointVelTotalIMU.block(3,0,3,1) = imu_omega;
    
    Eigen::Matrix<double, 3, 4> stance_feet = Eigen::MatrixXd::Zero(3,4);
    for (size_t i = 0; i < 4; i++)
    {
        stance_feet.block(0,i,3,1) = con->ind[i]*kin->toePos.block(0,i,3,1); 
    }
    
    Eigen::Matrix<double, 3, 4> feet0CoM = quad->FootEstimator(jointPosIMU);

    Eigen::Matrix<double, 3, 4> pCoM_raw = Eigen::MatrixXd::Zero(3,4);
    for (size_t i = 0; i < 4; i++)
    {
            pCoM_raw.block(0,i,3,1) = stance_feet.block(0,i,3,1)-con->ind[i]*feet0CoM.block(0,i,3,1); 
    }

    p_est.block(0,0,3,1) = pCoM_raw.rowwise().sum()/con->cnt;
    
    Eigen::Matrix<double, 12, 18> JacobianFull = quad->JacobianEstimator(jointPosIMU);
    Eigen::Matrix<double, 3, 1> pdot_Raw = Eigen::MatrixXd::Zero(3,1);
    for(size_t i=2; i<4; i++){
        
        pdot_Raw -= con->ind[i]*JacobianFull.block(3*i,3,3,15)*jointVelTotalIMU.block(3,0,15,1);
        
    }
    p_est.block(3,0,3,1) = pdot_Raw/(con->ind[2]+con->ind[3]);
    
    return p_est;

}


void LocoWrapper::getStateEstimatefull(double q[18], double dq[18], const int* contact, Eigen::Matrix<double,3,3> R, int robotdown, bool dynswitch, size_t ctrlTick){
    
    float numContact = (contact[0]+contact[1])+rearfootweight*contact[2]+rearfootweight*contact[3];
    // if(ctrlTick>27399){
        // numContact = (contact[0]+contact[1])*robotdown+rearfootweight*contact[2]+rearfootweight*contact[3];
    // }
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
	double fr_prev[3] = {kin->toePos(0,0),kin->toePos(1,0),kin->toePos(2,0)};//{fr_toe[0],fr_toe[1],fr_toe[2]};
	double fl_prev[3] = {kin->toePos(0,1),kin->toePos(1,1),kin->toePos(2,1)};//{fl_toe[0],fl_toe[1],fl_toe[2]};
	double rr_prev[3] = {kin->toePos(0,2),kin->toePos(1,2),kin->toePos(2,2)};//{rr_toe[0],rr_toe[1],rr_toe[2]};
	double rl_prev[3] = {kin->toePos(0,3),kin->toePos(1,3),kin->toePos(2,3)};//{rl_toe[0],rl_toe[1],rl_toe[2]};
	
    double deltaPos[3] = {0.0,0.0,0.0};
    for(int i=0; i<3; ++i){
        // if(ctrlTick<27400){
            deltaPos[i] -= (fr_toe[i]-fr_prev[i])*contact[0];
            deltaPos[i] -= (fl_toe[i]-fl_prev[i])*contact[1];
        // }
        deltaPos[i] -= (rr_toe[i]-rr_prev[i])*contact[2]*rearfootweight;
        deltaPos[i] -= (rl_toe[i]-rl_prev[i])*contact[3]*rearfootweight;
        deltaPos[i] /= numContact;
    }    
    
	COM[0] = deltaPos[0];
	COM[1] = deltaPos[1];
    COM[2] = deltaPos[2];//-1.0*(fr_toe[2]*contact[0]+fl_toe[2]*contact[1]+rr_toe[2]*contact[2]+rl_toe[2]*contact[3])/numContact;
	
	// for(int i=0; i<3; ++i){
	// 	fr_prev[i] = fr_toe[i]; fl_prev[i] = fl_toe[i];
	// 	rr_prev[i] = rr_toe[i]; rl_prev[i] = rl_toe[i];		
	// }
	
	numContact = (contact[0]+contact[1])*robotdown + rearfootweight*contact[2]+rearfootweight*contact[3];
    //numContact = (contact[0]+contact[1]) + rearfootweight*contact[2]+rearfootweight*contact[3];
	Eigen::Matrix<double,3,1> dq_temp = {dq[3],dq[4],dq[5]};
	
    if(dynswitch){

        for (int i = 3; i < 18; ++i){
		    COM_vel[0] -= (Jfr_toe[3*i+0]*contact[0]*robotdown + Jfl_toe[3*i+0]*contact[1]*robotdown + Jrr_toe[3*i+0]*contact[2]*rearfootweight + Jrl_toe[3*i+0]*contact[3]*rearfootweight)*dq[i];
	 	    COM_vel[1] -= (Jfr_toe[3*i+1]*contact[0]*robotdown + Jfl_toe[3*i+1]*contact[1]*robotdown + Jrr_toe[3*i+1]*contact[2]*rearfootweight + Jrl_toe[3*i+1]*contact[3]*rearfootweight)*dq[i];
	 	    COM_vel[2] -= (Jfr_toe[3*i+2]*contact[0]*robotdown + Jfl_toe[3*i+2]*contact[1]*robotdown + Jrr_toe[3*i+2]*contact[2]*rearfootweight + Jrl_toe[3*i+2]*contact[3]*rearfootweight)*dq[i];
            // COM_vel[0] -= (Jfr_toe[3*i+0]*contact[0] + Jfl_toe[3*i+0]*contact[1] + Jrr_toe[3*i+0]*contact[2]*rearfootweight + Jrl_toe[3*i+0]*contact[3]*rearfootweight)*dq[i];
	 	    // COM_vel[1] -= (Jfr_toe[3*i+1]*contact[0] + Jfl_toe[3*i+1]*contact[1] + Jrr_toe[3*i+1]*contact[2]*rearfootweight + Jrl_toe[3*i+1]*contact[3]*rearfootweight)*dq[i];
	 	    // COM_vel[2] -= (Jfr_toe[3*i+2]*contact[0] + Jfl_toe[3*i+2]*contact[1] + Jrr_toe[3*i+2]*contact[2]*rearfootweight + Jrl_toe[3*i+2]*contact[3]*rearfootweight)*dq[i];

        }
	    COM_vel[0] /= numContact;
	    COM_vel[1] /= numContact;
	    COM_vel[2] /= numContact;

    }else{

        toWorld(&dq[3],dq_temp,R);
	    for (int i = 3; i < 18; ++i){
		    COM_vel[0] -= (Jfr_toe[3*i+0]*contact[0]*robotdown + Jfl_toe[3*i+0]*contact[1]*robotdown + Jrr_toe[3*i+0]*contact[2]*rearfootweight + Jrl_toe[3*i+0]*contact[3]*rearfootweight)*dq[i];
	 	    COM_vel[1] -= (Jfr_toe[3*i+1]*contact[0]*robotdown + Jfl_toe[3*i+1]*contact[1]*robotdown + Jrr_toe[3*i+1]*contact[2]*rearfootweight + Jrl_toe[3*i+1]*contact[3]*rearfootweight)*dq[i];
	 	    COM_vel[2] -= (Jfr_toe[3*i+2]*contact[0]*robotdown + Jfl_toe[3*i+2]*contact[1]*robotdown + Jrr_toe[3*i+2]*contact[2]*rearfootweight + Jrl_toe[3*i+2]*contact[3]*rearfootweight)*dq[i];
	    }
	    COM_vel[0] /= numContact;
	    COM_vel[1] /= numContact;
	    COM_vel[2] /= numContact;
	
	    dq_temp = {dq[3],dq[4],dq[5]};
	    toBody(&dq[3],dq_temp,R);
    }

	// Set results
	q[0] = COM[0]; q[1] = COM[1]; q[2] = COM[2];
    if(ctrlTick<27000){
	    dq[0] = COM_vel[0] > xdot_thresh ? xdot_thresh : (COM_vel[0] < -xdot_thresh ? -xdot_thresh : COM_vel[0]); 
        dq[1] = COM_vel[1] > yzdot_thresh ? yzdot_thresh : (COM_vel[1] < -yzdot_thresh ? -yzdot_thresh : COM_vel[1]);
        dq[2] = COM_vel[2] > yzdot_thresh ? yzdot_thresh : (COM_vel[2] < -yzdot_thresh ? -yzdot_thresh : COM_vel[2]); 
        //dq[2] = COM_vel[2];
    }else{
        dq[0] = COM_vel[0] > xdot_thresh2 ? xdot_thresh2 : (COM_vel[0] < -xdot_thresh2 ? -xdot_thresh2 : COM_vel[0]); 
        dq[1] = COM_vel[1] > yzdot_thresh2 ? yzdot_thresh2 : (COM_vel[1] < -yzdot_thresh2 ? -yzdot_thresh2 : COM_vel[1]);
        dq[2] = COM_vel[2] > yzdot_thresh2 ? yzdot_thresh2 : (COM_vel[2] < -yzdot_thresh2 ? -yzdot_thresh2 : COM_vel[2]); 
    }
}


void LocoWrapper::ExpWrapper(const double jpos_est[18], const double jvel_est[18], const double rotMatrixDouble[9], size_t control_Tick, size_t solveduration, 
                                int HLContactIndex[5], Eigen::Matrix<double, 12, 1> comDes, Eigen::Matrix<double, 17, 1> fDes){

    setoptNLstateExp(comDes,fDes);
    setHLphase(HLContactIndex[4]);
    if(control_Tick < loco_start_e+shifttime){
        nextcon_e(0) = 0;
        
        if(control_Tick==loco_start_e){
            Eigen::Matrix<double, 4, 1> wfoot = rearweight*Eigen::MatrixXd::Ones(4,1);
            wfoot(0)=1;wfoot(1)=1;
            getshiftedCoM(wfoot);setshiftedCoM();
        }
        setswingContact(nextcon_e);
        calcTau2(jpos_est,jvel_est,rotMatrixDouble,STANDUP,control_Tick,solveduration);

    }else if(control_Tick >= loco_start_e + shifttime & control_Tick < switchtime*ctrlHz){// & controlTick < loco_start + shifttime){ // Start locomotion
    
        stepind_e = std::floor((control_Tick-loco_start_e-shifttime)/(shifttime+movetime));
    
        if(stepind_e<maxsteps){
            stepsonwall(stepind_e);
    
            if(stepind_e==1){tookfirststep();}
    
            if(stepind_e%2==0){
                nextcon_e(0) = 0;nextcon_e(1) = 1;
                if(control_Tick == loco_start_e +shifttime + stepind_e*(shifttime+movetime)){incstep();}       
            }else{
                nextcon_e(0) = 1;nextcon_e(1) = 0;
            }
     
            if(control_Tick==loco_start_e + (stepind_e+1)*(movetime + shifttime)){
                Eigen::Matrix<double, 4, 1> wfoot = rearweight*Eigen::MatrixXd::Ones(4,1);//3
                wfoot(0)=2*nextcon_e(1);
                wfoot(1)=2*nextcon_e(0);  
                getshiftedCoM(wfoot);setshiftedCoM();
            }
            setswingContact(nextcon_e);
     
        }else{

            stepind2_e = std::floor((control_Tick-loco_start_e-shifttime-(maxsteps)*(movetime + shifttime))/(shifttime2+movetime2));
            settlesteps(stepind2_e);
            if(stepind2_e==settlingsteps){gotfinalstate();}

            if(control_Tick==loco_start_e + shifttime + (maxsteps)*(movetime + shifttime)){
                stopclimbing();
                stepsonwall(maxsteps);
                setfinalCoM();
            } 
        }
        calcTau2(jpos_est,jvel_est,rotMatrixDouble,STANDUP,control_Tick,solveduration);
        
    }else{

        stepind2_e = std::floor((control_Tick-switchtime*ctrlHz)/(shifttime2+movetime3));
        settlesteps(stepind2_e);

        // switch (control_Tick) {
        //     case switchtime*ctrlHz: readytowalk(); break;
        //     case switchtime*ctrlHz + stepind2_e*(shifttime2+movetime3): setfinalCoM2(); break;
        //     case switchtime*ctrlHz+2*(shifttime2+movetime3)+shifttime2: startwalking(); break; // This executes
        //     case 30000: readytoreallywalk(); break;
        //     case 31000: readytoreallywalk(); break;
        //     case 32000: letsgo(); break;
        // }

        if(control_Tick == switchtime*ctrlHz){readytowalk();}//break;}
        if(control_Tick == switchtime*ctrlHz + stepind2_e*(shifttime2+movetime3)){setfinalCoM2();}//break;}
        if(control_Tick == switchtime*ctrlHz+2*(shifttime2+movetime3)+shifttime2){startwalking();}//break;}
        if(control_Tick == (switchtime+6)*ctrlHz){readytoreallywalk();}//break;}
        if(control_Tick == (switchtime+7)*ctrlHz){readytoreallywalk();}//break;}

        //updatestate(jpos_est,jvel_est,rotMatrixDouble);
        if(readytowalkf){setcontactconfigExp(HLContactIndex);}
        calcTau2(jpos_est,jvel_est,rotMatrixDouble,UPWALK,control_Tick,solveduration);

    }
    
}

void LocoWrapper::setcontactconfigExp(int HLContactIndex[5]){
    
    desired_contact[0] = HLContactIndex[0];
    desired_contact[1] = HLContactIndex[1];
    desired_contact[2] = HLContactIndex[2];
    desired_contact[3] = HLContactIndex[3];
   
    conEst->setDesDomain(desired_contact);
    quad->updateSwingMatrices(con->ind,con->cnt); 
}

void LocoWrapper::setoptNLstateExp(Eigen::Matrix<double, 12, 1> comDes, Eigen::Matrix<double, 17, 1> fDes){
    
    opt_HLstate.block(0,0,12,1) = comDes;
    //opt_HLstate(2) = 0.35;
    opt_HLstate.block(12,0,12,1) = fDes.block(0,0,12,1);
    NLstep = fDes.block(12,0,5,1);
}


void LocoWrapper::calcTau2k(const double q[18], const double dq[18], const double R[9], size_t gait, size_t ctrlTick, size_t solveduration){

    phaseVar = getPhase(1.0*locoTick, 0.0, 1.0*traj->domLen);   // update phase variable
    phaseVar = (phaseVar>1) ? 1 : ((phaseVar<0) ? 0 : phaseVar);

    quad->updateState(q,dq,R); 
    
    if(gait==STAND){
        //conEst->updateConState(footPos,phaseVar,force);
        //phaseVar = getPhase(1.0*locoTick, 0.0, 1.0*traj->domLen);
        
        PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);
        quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
        VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
        LL->calcTorque(state, dyn, kin, vcon, con, &ll_params);   
                                         // run low level controller                                       // log relavent data
    }else if(gait == UPWALK){

        
        if(readytowalkf){
        
            // if (gait!=gaitTemp || (phaseVar>maxPhase) || ctrlTick == switchtime*ctrlHz + 2*(shifttime2+movetime3) + shifttime2){ 
            if (HLMTphase == 0 || HLMTphase == 20 || ctrlTick == switchtime*ctrlHz + 2*(shifttime2+movetime3) + shifttime2){ 

                locoTick = 0;
                phaseVar = getPhase(1.0*locoTick, 0.0, 199);
                flphase = 0;
                rlphase = 0;
                //conEst->forceDomChange();
                z = Eigen::MatrixXd::Zero(12,1);
            
            }else{
                phaseVar = getPhase(1.0*locoTick, 0.0, 199);
                flphase = getPhase(1.0*locoTick-49, 0.0, 145);
                rlphase = getPhase(1.0*locoTick-79, 0.0, 115);
                flphase = (flphase<0)?0:(flphase>1)?1:flphase;
                rlphase = (rlphase<0)?0:(rlphase>1)?1:rlphase;
                // if(locoTick==50||locoTick==80){
                if (HLMTphase == 5 || HLMTphase == 8 || HLMTphase == 25 || HLMTphase == 28){
                    conEst->forceDomChange();
                }            
            }
   
            PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);  
            VC->updateVirtualConstraintswalk(state, kin, traj, con, gait, flphase,rlphase, &motion_params, ll);    // update VC's    
            VC->setDesiredForce(opt_HLstate.block(12,0,12,1));
            z.block(6,0,3*(4-con->cnt),0) += vcon->y.block(6,0,3*(4-con->cnt),0)/ctrlHz;         
            LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, HRai, z, Ki);                                     // run low level controller
        
        }else{
            
            if( ctrlTick == switchtime*ctrlHz+settlestep*(shifttime2+movetime3) || ctrlTick == switchtime*ctrlHz+settlestep*(shifttime2+movetime3)+shifttime2){
                locoTick = 0;
                phaseVar = 0;
                PP->setToeInit(kin);
                PP->setx0y0z0(state->q(0),state->q(1),state->q(2),state->q(4));
            }

            if(ctrlTick < switchtime*ctrlHz + settlestep*(shifttime2+movetime3)+shifttime2){
                PP->shiftCoM3(conEst,phaseVar,shifttime2,true);
                quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
                VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
                VC->setDesiredForce(opt_HLstate.block(12,0,12,1));
                LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
            }else{
                if(settlestep%2==0){
                    nextContact[0] = 0;
                    nextContact[1] = 1;
                }else{
                    nextContact[0] = 1;
                    nextContact[1] = 0;
                }
                conEst->setDesDomain(nextContact);
                PP->movefoot3(movetime3);
                quad->updateSwingMatrices(con->ind,con->cnt);
                VC->updateVirtualConstraintssetfoot(state, kin, traj, con, gait, phaseVar, &motion_params, ll, flphase,rlphase,true);
                VC->setDesiredForce(opt_HLstate.block(12,0,12,1));
                LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
            }

        }


    }else{

        if(stopclimb){
            wallstep = wallstep-maxsteps;
        }
        if(ctrlTick == standoffset+(wallstep)*(currentshift+movetime)||ctrlTick == standoffset+(wallstep)*(currentshift+movetime)+currentshift){
            locoTick = 0;
            phaseVar = 0;
            PP->setToeInit(kin);
            PP->setx0y0z0(state->q(0),state->q(1),state->q(2),state->q(4));
            if(ctrlTick == loco_start_e+shifttime){
                setrearhippose();
            }
        }
        
        if(ctrlTick < standoffset+currentshift+(wallstep)*(currentshift+movetime)){
            PP->shiftCoMk(conEst,wallstep,phaseVar,currentshift,stopclimb);
            quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
            VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll); 
        }else{
            conEst->setDesDomain(nextContact);
            PP->movefootk(movetime,wallstep,stopclimb);//,wallstep);//(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);
            quad->updateSwingMatrices(con->ind,con->cnt);
            VC->updateVirtualConstraintssetfoot(state, kin, traj, con, gait, phaseVar, &motion_params, ll, flphase,rlphase,true);
        }

        if(ctrlTick < loco_start_e+shifttime){
            LL->calcTorque(state, dyn, kin, vcon, con, &ll_params);
        }else{
            LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, Hr, z, Ki);
        }
    }
                             
    data->writeData(state,vcon,con,traj,ll,kin,ctrlTick,force_LL,opt_HLstate,locoTick,phaseVar,flphase,rlphase,0.0,0.0,NLstep,solveduration);
    locoTick += (ctrlHz)/LL_Hz;     // increment locoTick
    gaitTemp = gait;

}


void LocoWrapper::ExpWrapperk(const double jpos_est[18], const double jvel_est[18], const double rotMatrixDouble[9], size_t control_Tick, size_t solveduration, 
    int HLContactIndex[5], Eigen::Matrix<double, 12, 1> comDes, Eigen::Matrix<double, 17, 1> fDes){

    setoptNLstateExp(comDes,fDes);
    setHLphase(HLContactIndex[4]);
    if(control_Tick < switchtime*ctrlHz){//loco_start_e+shifttime){
        
        stepind_e = std::floor((control_Tick-loco_start_e)/(shifttime+movetime));
        stepsonwall(stepind_e);
        
        if(stepind_e%2==0){nextcon_e(0) = 0;nextcon_e(1) = 1;}      
        else{nextcon_e(0) = 1;nextcon_e(1) = 0;}
        setswingContact(nextcon_e);

        if(stepind_e<maxsteps){
            standoffset = loco_start_e;
            currentshift = shifttime;
            
            if(stepind_e==1){tookfirststep();}
            
            if(control_Tick==loco_start_e + stepind_e*(shifttime+movetime)){
                Eigen::Matrix<double, 4, 1> wfoot = rearweight*Eigen::MatrixXd::Ones(4,1);
                wfoot(0)=2*nextcon_e(1);
                wfoot(1)=2*nextcon_e(0);
                getshiftedCoM(wfoot);setshiftedCoM();
            }
            
        }else{
            standoffset = loco_start_e + maxsteps*(shifttime+movetime);
            currentshift = shifttime2;
            if(control_Tick==standoffset){// + (maxsteps)*(movetime + shifttime)){
                stopclimbing();
                setfinalCoM();
            } 
        }
        calcTau2k(jpos_est,jvel_est,rotMatrixDouble,STANDUP,control_Tick,solveduration);

    }else{

        stepind2_e = std::floor((control_Tick-switchtime*ctrlHz)/(shifttime2+movetime3));
        settlesteps(stepind2_e);


        if(control_Tick == switchtime*ctrlHz){readytowalk();}//break;}
        if(control_Tick == switchtime*ctrlHz + stepind2_e*(shifttime2+movetime3)){setfinalCoM2();}//break;}
        if(control_Tick == switchtime*ctrlHz+2*(shifttime2+movetime3)+shifttime2){startwalking();}//break;}
        if(control_Tick == (switchtime+6)*ctrlHz){readytoreallywalk();}//break;}
        if(control_Tick == (switchtime+7)*ctrlHz){readytoreallywalk();}//break;}

        //updatestate(jpos_est,jvel_est,rotMatrixDouble);
        if(readytowalkf){setcontactconfigExp(HLContactIndex);}
        calcTau2k(jpos_est,jvel_est,rotMatrixDouble,UPWALK,control_Tick,solveduration);

    }

}

// if(control_Tick == switchtime*ctrlHz+2*(shifttime2+movetime3)+shifttime2){startwalking();}
// if(control_Tick==30000 || control_Tick==31000){readytoreallywalk();}
// if(control_Tick==32000){letsgo();}
// if(control_Tick == switchtime*ctrlHz){readytowalk();}//lastQPforce = getpreviousQPforce();
// if(control_Tick == switchtime*ctrlHz + stepind2_e*(shifttime2+movetime3)){setfinalCoM2();}