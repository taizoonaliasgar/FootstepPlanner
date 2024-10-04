//
// Authror: Randy Fawcett on 12/2021.
//
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "LocoWrapperwalk.hpp"
#include "iostream"

LocoWrapperwalk::LocoWrapperwalk(int argc, char *argv[]) : Parameters(argc,argv){

//    std::string filename = "/media/kavehakbarihamed/Data/A1_RaiSim_Outputs/LCSS_2021/Payload_Trot_10cm.txt";
//    std::string filename = "/media/kavehakbarihamed/Data/A1_RaiSim_Outputs/nothing.txt";
    std::string filename = "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/nothing.csv";
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
    
    locoTick = 0;
    maxPhase = 0.996;

    opt_HLstate = Eigen::MatrixXd::Zero(24,1);
    opt_HLstate(14) = 0.35;
    two_contact = 10;

    
    contact_horizon.block(1,25,1,15) = Eigen::MatrixXd::Zero(1,15);//15
    contact_horizon.block(2,28,1,12) = Eigen::MatrixXd::Zero(1,12);
    contact_horizon.block(0,5,1,15) = Eigen::MatrixXd::Zero(1,15);//15
    contact_horizon.block(3,8,1,12) = Eigen::MatrixXd::Zero(1,12);

   /*  contact_horizon.block(1,20-two_contact,1,two_contact) = Eigen::MatrixXd::Zero(1,two_contact);
    contact_horizon.block(2,20-two_contact,1,two_contact) = Eigen::MatrixXd::Zero(1,two_contact);
    contact_horizon.block(0,20+two_contact,1,two_contact) = Eigen::MatrixXd::Zero(1,two_contact);
    contact_horizon.block(3,20+two_contact,1,two_contact) = Eigen::MatrixXd::Zero(1,two_contact); */

    //contact_horizon.block(1,0,1,20) = Eigen::MatrixXd::Zero(1,20);
    //contact_horizon.block(2,0,1,20) = Eigen::MatrixXd::Zero(1,20);
    //contact_horizon.block(0,20,1,20) = Eigen::MatrixXd::Zero(1,20);
    //contact_horizon.block(3,20,1,20) = Eigen::MatrixXd::Zero(1,20);
    /* contact_horizon.block(1,20,1,10) = Eigen::MatrixXd::Zero(1,10);
    contact_horizon.block(2,20,1,10) = Eigen::MatrixXd::Zero(1,10);
    contact_horizon.block(0,30,1,10) = Eigen::MatrixXd::Zero(1,10);
    contact_horizon.block(3,30,1,10) = Eigen::MatrixXd::Zero(1,10); */

    //NLstep = Eigen::MatrixXd::Zero(4,1);
    
}

LocoWrapperwalk::~LocoWrapperwalk(){
    delete quad;
    delete conEst;
    delete LL;
    delete VC;
    delete PP;
}

void LocoWrapperwalk::calcTau(const double q[18], const double dq[18], const double R[9], const int force[4], size_t gait, size_t ctrlTick){//}, Eigen::Matrix<double, 24, 1>& opt_HLstate){
    
    quad->updateState(q,dq,R);      // update state
    float footPos[4] = {0};  // DUMMY VARS

    //std::cout<< "Step 1" << std::endl;
    if(gait == STAND){

        conEst->updateConState(footPos,phaseVar,force);
        PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);
        quad->updateSwingMatrices(con->ind,con->cnt);                                               // update the jacobian    
        VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);    // update VC's
        LL->calcTorque(state, dyn, kin, vcon, con, &ll_params);                                     // run low level controller                                       // log relavent data

    }else{
        
        float footPos[4] = {0};  // DUMMY VARS
        if (ctrlTick<1 || gait!=gaitTemp || (phaseVar>maxPhase && gait!=STAND) ){ 
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
            if(locoTick==50||locoTick==80){
                conEst->forceDomChange();
            }            
        }

        
        PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);  
        VC->updateVirtualConstraintswalk(state, kin, traj, con, gait, flphase,rlphase, &motion_params, ll);    // update VC's
        VC->setDesiredForce(opt_HLstate.block(12,0,12,1));
        
        z.block(6,0,3*(4-con->cnt),0) += vcon->y.block(6,0,3*(4-con->cnt),0)/ctrlHz;         

        LL->calcTorquewalk(state, dyn, kin, vcon, con, &ll_params, HRai, z, Ki);                                     // run low level controller
        //if(locoTick==0){
        //    conEst->forceDom0();
        //}
    }

    //Eigen::FullPivLU<Eigen::Matrix<double,12,18>> lu_decomp(kin->Jtoe);
    //auto rank = lu_decomp.rank();
    //std::cout << rank << std::endl;

    Ddiff = (dyn->D-DRai).norm()/dyn->D.norm();
    Hdiff = (dyn->H-HRai).norm()/dyn->H.norm();

    data->writeData(state,vcon,con,traj,ll,kin,ctrlTick,force,opt_HLstate,locoTick,phaseVar,flphase,rlphase,Ddiff,Hdiff,NLstep); 
    
    locoTick += (ctrlHz)/LL_Hz;     // increment locoTick
    gaitTemp = gait;
}

void LocoWrapperwalk::setoptNLstate(Eigen::Matrix<double, 33, 1> HLopt){
    
    opt_HLstate.block(0,0,12,1) = HLopt.block(0,0,12,1);
    //opt_HLstate(2) = 0.35;
    opt_HLstate.block(12,0,12,1) = HLopt.block(16,0,12,1);
    NLstep = HLopt.block(28,0,5,1);
}


void LocoWrapperwalk::setcontactconfig(int controlMPC){
    
    int u= 0;
    for(int i=0;i<4;i++){
        desired_contact[i]=contact_horizon(i,controlMPC%40);
        if(desired_contact[i]<1){
            Ki(6+3*u,6+3*u) = 1-desired_contact[i];
            u++;
        }
    }
    Ki = 0*Ki;
    
    conEst->setDesDomain(desired_contact);
    quad->updateSwingMatrices(con->ind,con->cnt); 
}