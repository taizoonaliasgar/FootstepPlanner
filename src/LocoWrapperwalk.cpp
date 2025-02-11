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
    //std::string filename = "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/nothing.csv";
    std::string filename = "../datalog/nothing.csv";
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

    CoMhistory.block(2,0,1,fitsample+1) = 0.5*Eigen::MatrixXd::Ones(1,fitsample+1);
    R_KF.block(3,3,3,3) = (Eigen::Matrix<double, 3, 3>() << 1.4892,  0.0431, -0.0366,
                                        0.0431, 0.6475, -0.0204,
                                        -0.0366, -0.0204, 0.2310).finished();
    Q_KF.block(3,3,3,3) = (Eigen::Matrix<double, 3, 3>() << 3.9631, 0.0695, 1.7835,
                                                0.0695, 5.4780, -0.0241,
                                                1.7835, -0.0241, 7.8273).finished();
    P_KF.block(3,3,3,3) = (Eigen::Matrix<double, 3, 3>() << 0.0066, 0.0015, -0.0058,
                                                                0.0015, 0.0110, -0.0006,
                                                                -0.0058, -0.0006, 0.0072).finished();
    P_KF.block(0,0,3,3) = 0.01*Eigen::MatrixXd::Identity(3,3);                                     
    x_est_prev(2) = 0.5;

    A_KF.block(0,3,3,3) = 0.001*Eigen::MatrixXd::Identity(3,3);
    B_KF.block(3,0,3,3) = 0.001*Eigen::MatrixXd::Identity(3,3);
}

LocoWrapperwalk::~LocoWrapperwalk(){
    delete quad;
    delete conEst;
    delete LL;
    delete VC;
    delete PP;
}

void LocoWrapperwalk::calcTau(const double q[18], const double dq[18], const double R[9], const int force[4], size_t gait, size_t ctrlTick, size_t duration){//}, Eigen::Matrix<double, 24, 1>& opt_HLstate){
    
    //std::cout << "LocoWrapperwalk" << "\t" << duration << std::endl;
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

        //std::cout << "Settin desired state" << std::endl;
        PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, opt_HLstate, NLstep);  
        
        if(RaibFlag){
            VC->updateVirtualConstraintswalkR(state, kin, traj, con, gait, flphase,rlphase, &motion_params, ll);    // update VC's   
        }else{
            //std::cout << "Settin VC" << std::endl;
            VC->updateVirtualConstraintswalk(state, kin, traj, con, gait, flphase,rlphase, &motion_params, ll);    // update VC's    
        }

        VC->setDesiredForce(opt_HLstate.block(12,0,12,1));
        
        z.block(6,0,3*(4-con->cnt),0) += vcon->y.block(6,0,3*(4-con->cnt),0)/ctrlHz;         
        //std::cout << "Getting torque" << std::endl;
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

    data->writeData(state,vcon,con,traj,ll,kin,ctrlTick,force,opt_HLstate,locoTick,phaseVar,flphase,rlphase,Ddiff,Hdiff,NLstep,duration); 
    
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

Eigen::Matrix<double, 12, 1> LocoWrapperwalk::getStateEstimate(double jointPos[18], Eigen::VectorXd jointVelTotal, Eigen::Matrix<double, 3, 1> imu_eul, Eigen::Matrix<double, 3, 1> imu_omega){
    
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
        //if(con->ind[i]==1){
            stance_feet.block(0,i,3,1) = con->ind[i]*kin->toePos.block(0,i,3,1); 
        //    std::cout << "stance_feet" << "\t" << i << "\t" << stance_feet.block(0,i,3,1).transpose() << std::endl;
        //}  
    }
    
    Eigen::Matrix<double, 3, 4> feet0CoM = quad->FootEstimator(jointPosIMU);
    //std::cout << "feet0CoM" << "\t" << feet0CoM << std::endl;

    Eigen::Matrix<double, 3, 4> pCoM_raw = Eigen::MatrixXd::Zero(3,4);
    for (size_t i = 0; i < 4; i++)
    {
        //if(con->ind[i]==1){
            pCoM_raw.block(0,i,3,1) = stance_feet.block(0,i,3,1)-con->ind[i]*feet0CoM.block(0,i,3,1); 
      //      std::cout << "pCoM_raw" << "\t" << i << "\t" << pCoM_raw.block(0,i,3,1).transpose() << std::endl;   
        //}
    }

    p_est.block(0,0,3,1) = pCoM_raw.rowwise().sum()/con->cnt;
    
    Eigen::Matrix<double, 12, 18> JacobianFull = quad->JacobianEstimator(jointPosIMU);
    Eigen::Matrix<double, 3, 1> pdot_Raw = Eigen::MatrixXd::Zero(3,1);
    for(size_t i=2; i<4; i++){
        
        //pdot_Raw.block(3*i,0,3,1) = -con->ind[i]*JacobianFull.block(3*i,3,3,15)*jointVelTotal.block(3,0,15,1);
        pdot_Raw -= con->ind[i]*JacobianFull.block(3*i,3,3,15)*jointVelTotalIMU.block(3,0,15,1);
        
    }
    p_est.block(3,0,3,1) = pdot_Raw/(con->ind[2]+con->ind[3]);

    // CoMhistory.block(0,0,3,fitsample) = CoMhistory.block(0,1,3,fitsample);
    // CoMhistory.block(0,fitsample,3,1) = p_est.block(0,0,3,1);
    // p_est.block(3,0,3,1) = getvEstimate();
    
    return p_est;

}

Eigen::Matrix<double, 3, 4> LocoWrapperwalk::getfootv(double jointPos[18],Eigen::VectorXd jointVelTotal){
    //std::cout << "getfootv" << std::endl;
    Eigen::Matrix<double, 12, 18> JacobianFull = quad->JacobianEstimator(jointPos);
    //std::cout << "JacobianFull" << "\t" << JacobianFull << std::endl;
    Eigen::Matrix<double, 3, 4> footv = Eigen::MatrixXd::Zero(3,4);
    for(size_t i=0; i<4; i++){
        footv.block(0,i,3,1) = JacobianFull.block(3*i,0,3,18)*jointVelTotal;
    }
    return footv;

}

Eigen::Matrix<double,3,1> LocoWrapperwalk::getvEstimate(){

    Eigen::Matrix<double,3,1> v_estimate = Eigen::MatrixXd::Zero(3,1);
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(fitsample + 1, fitorder + 1);
    
    X(0,0)=1;
    for (size_t i = 1; i <= fitsample; ++i){
        for (int j = 0; j <= fitorder; ++j){
            X(i, j) = std::pow(i*0.001, j);
        }
    }

    Eigen::MatrixXd Y = Eigen::VectorXd::Zero(fitsample + 1,1);
    
    for (size_t f_ind = 0; f_ind < 3; ++f_ind){

        Y = CoMhistory.block(f_ind,0,1,fitsample+1).transpose();    
        Eigen::VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * Y);
        
        for(int i=1; i<=fitorder; i++){
            v_estimate(f_ind) += i*coeffs(i)*std::pow((fitsample)*0.001, i-1);
        }
        
    }

    return v_estimate;
}

Eigen::Matrix<double, 6, 1> LocoWrapperwalk::VelKF(Eigen::Matrix<double, 6, 1> x_est, Eigen::Matrix<double, 3, 1> a_est){
    
    Eigen::Matrix<double, 6, 1> x_est_process = x_est_prev + B_KF*a_est;
    
    P_KF = A_KF*P_KF*A_KF.transpose() + Q_KF;
    K_KF = P_KF*(P_KF+R_KF).inverse();
    
    x_est_prev = x_est_process + K_KF*(x_est - x_est_process);
    P_KF = (Eigen::MatrixXd::Identity(6,6)-K_KF)*P_KF;

    return x_est_prev;
}