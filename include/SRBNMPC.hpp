#ifndef SINGLERIGIDBODY_NMPC
#define SINGLERIGIDBODY_NMPC

#include "global_loco_structs.hpp"
#include "Parameters.hpp"
#include "EigenUtils.hpp"
#include "Transforms.hpp"
#include "Bezier.h"
#include "MPC_global_params.hpp"
//#include "iswift_qp.h"
#include <casadi/casadi.hpp>

using DynInf = DynamicsInfo;
using KinInf = KinematicsInfo;
using ConInf = ContactInfo;
using MPCP = Settings::MPC_params;

#include "fstream"

const size_t FILE_CNT_MPC = 5;
const std::string FILE_NAMES_MIT[FILE_CNT_MPC] = {
    "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/state.txt",
    "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/inputs.txt",
    "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/desired_state.txt",
    "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/foot_position.txt",
    "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/contact_index.txt"
//     "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/footposition.txt",
//    // "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/forceinitial.txt",
//     "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/steplength.txt",
//     "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/foothold.txt",
//     "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/forcesol.txt"
};


class SRBNMPC : public Parameters
{
public:
    std::fstream file[FILE_CNT_MPC];
    SRBNMPC(int argc, char *argv[], int numRobots, int id);
    virtual ~SRBNMPC(){
        for(size_t i=0; i<FILE_CNT_MPC; i++){
            if(file[i].is_open()){
                file[i].close();
            }
        }
        std::cout << "file close done" << std::endl;   
    };

    void RunMPC(size_t tick, Eigen::Matrix<double,12,1> &q0, size_t gait);
    void UpdateCost();
    void UpdateConstraints(Eigen::Matrix<double,12,1> &q0);
    void GetLinearDynamics(Eigen::Matrix<double,12,1> &q0);
    void PlanTraj(size_t tick, Eigen::Matrix<double,NRS,1> &q, size_t gait);
    int getRobotID(){return agentID;};
    void setStandVars(double x, double y, double z);

    // Helper functions
    void ExtractYPR(const Eigen::Matrix<double,3,3> &R, Eigen::Matrix<double,3,1> &eul);
    void SkewMat(const Eigen::Matrix<double, 3, 4> &pf, Eigen::Matrix<double,3,12> &skew);
    Eigen::Matrix<double,3,1> getdisturbance(){return disturbforce;};
    inline void Rotate_z(Eigen::Matrix<double,3,3> &R, double ang){
        R(0,0) = cos(ang);  R(0,1) = -sin(ang);  R(0,2) = 0;
        R(1,0) = sin(ang);  R(1,1) =  cos(ang);  R(1,2) = 0;
        R(2,0) = 0;         R(2,1) =  0;         R(2,2) = 1;
    }

    //Taizoon functions
    void impactDetection(size_t tick, Eigen::Matrix<double,12,1> &q0, size_t gait);
    void impactDetectionTrot(size_t tick, Eigen::Matrix<double,12,1> &q0, size_t gait);
    void footstepplanner(Eigen::Matrix<double,12,1> &q0);
    void mpcdataLog(Eigen::Matrix<double,16,1> q0, Eigen::Matrix<double,12,1> force, size_t tick,Eigen::Matrix<double, 12, 1> p_foot);//(Eigen::Matrix<double,12,1> &q0, size_t tick);
    void motionPlannerwll(Eigen::Matrix<double,12,1> &q0);

    //Integration with low level
    void plannerforlowlevel(size_t tick, Eigen::Matrix<double,12,1> &q0, size_t gait, Eigen::Matrix<double, 3, 4> foot_position, Eigen::Matrix<double, 3, 4> hip_position);
    Eigen::Matrix<double, 24, 1> getHLMPCresults(){return opt_HLsol_;};

    //NMPC
    void generator();//(casadi::DM p);
    casadi::DM motionPlannerN(Eigen::Matrix<double,16,1> q0,size_t controlTick);
    casadi::SX UpdateCostN(casadi::SX x, casadi::SX x_des);
    casadi::SX UpdateConstraintsN(casadi::SX x, casadi::SX p);
    casadi::SX NonlinearDynamics(casadi::SX st,casadi::SX con, casadi::SX conc,casadi::SX conp1,casadi::SX conm1);
    //casadi::MX NonlinearDynamics(casadi::MX st, casadi::MX con);//, int i, casadi::MX p);
    casadi::SX GetTorque(casadi::SX st,casadi::SX con);
    casadi::SX inequalitycons(casadi::SX con);
    casadi::SX skewsym(casadi::SX v3);
    casadi::DM lowerboundx(casadi::DM p, int controlMPC);
    casadi::DM upperboundx(casadi::DM p);
    casadi::DM lowerboundg();
    casadi::DM upperboundg();
    //void setprevioussol(casadi::Matrix<double> sol0){  previous_sol = sol0; };
    void setprevioussol(casadi::DM sol0);
    casadi::DM getprevioussol(Eigen::Matrix<double,16,1> q0, size_t controlTick);
    Eigen::Matrix<double,12,1> getOptforce(int phase, int order);
    void getOptForceCoeff(int order);
    Eigen::Matrix<double,12,1> getFootPos();
    void writeMatrixToFile(const casadi::SX& matrix, const std::string& filename);
    void writeMatrixToFileDM(const casadi::DM& matrix, const std::string& filename);
    //Eigen::VectorXd polyfit(const std::vector<double>& x, const std::vector<double>& y, int order);
    
    Eigen::Matrix<double,12,HORIZ> arrangeOptforce();
    void getForce();
    Eigen::Matrix<double,12,1> returnForce(){return OPTforce;};
    void setpreviousp(casadi::DM p);
    void setsteplb(int controlMPC);

    //Fill simulation
    Eigen::Matrix<double,32,1> getNMPCsol(int controlMPC);
    Eigen::Matrix<double,33,1> getNMPCsol2(int controlMPC);
    casadi::DM getprevioussol_ll(Eigen::Matrix<double,16,1> q0, Eigen::Matrix<double,3,4> foothold, size_t controlTick);

private: 
    std::string filename;
    std::fstream fid;
    int success = 0;
    int agentID;
    int numAgents;

    Eigen::Matrix<double,3,3> Jinv;
    Eigen::Matrix<double,3,3> J;

    Eigen::Matrix<double,12,12> A, B;
    Eigen::Matrix<double,12,1> D;

    Eigen::Matrix<double,NRS*HORIZ,1> qdes = Eigen::MatrixXd::Zero(NRS*HORIZ,1);

    Eigen::Matrix<double,3,3> Rz;
    Eigen::Matrix<double,3,4> hip_pos;
    Eigen::Matrix<double,3,4> foot_pos;
    Eigen::Matrix<double,NRI,HORIZ> optforce;
    Eigen::Matrix<double,NRS,HORIZ> opttraj;

    // TOO LARGE FOR THE STACK. PUT IS ON THE HEAP
    // Eigen::Matrix<double,(NRS+NRI)*HORIZ,(NRS+NRI)*HORIZ> P_QP;
    // Eigen::Matrix<double,(NRS+NRI)*HORIZ,1> c_QP;
    // Eigen::Matrix<double,(NRS)*HORIZ,(NRS+NRI)*HORIZ> A_QP;
    // Eigen::Matrix<double,(NRS)*HORIZ,1> b_QP;
    // Eigen::Matrix<double,24*HORIZ,(NRS+NRI)*HORIZ> G_QP;
    // Eigen::Matrix<double,24*HORIZ,1> h_QP;

    Eigen::MatrixXd P_QP;
    Eigen::MatrixXd c_QP;
    Eigen::MatrixXd A_QP;
    Eigen::MatrixXd b_QP;
    Eigen::MatrixXd G_QP;
    Eigen::MatrixXd h_QP;

    
    // For trajectory planning
    // Not likely to be used when the low-level controller is involved
    double x0;
    double xf;
    double y0;
    double yf;
    double z0;
    double standTime;
    TrajInfo traj;
    size_t poseType = POSE_Z;

    // "Static" locomotion variables
    int phaseIdx = 2;
    int startTrot = 0;
    Eigen::Matrix<double, 3, 1> desVel = {0.5,0,0};
    Eigen::Matrix<double, 3, 1> desOmega = {0,0,0};
    Eigen::Matrix<double, 3, 1> desVelWrld = {0,0,0};
    double yawLock = 0;

    int stepcnt = 0;

    int con[4] = {1,0,0,1};

    //Taizoon params
    int con_prev[4] = {1,0,0,1};
    int locomotionTick = 0; // reset every domain
    size_t mpcTick = 0;
    int16_t stance_flag = 1;
    int16_t simcounter = 0;
    int16_t planindex = 0;
    int16_t M = HORIZ;
    double stand_height = 0.5;
    Eigen::Matrix<double, DOMAINSTEPS+HORIZ,1> vz_domain;
    double z_travel = 0;
    double fzmaxr = 250;
    double fzmaxf = 150;
    double mu= 0.6/sqrt(2); //mpc_params.mu_MPC;
    Eigen::Matrix<double,3,4> foot_current = Eigen::MatrixXd::Zero(3,4);
    Eigen::Matrix<double,12,12> Anext, Bnext;

    //Biped walk
    Eigen::Matrix<double,4,40> contact_horizon = Eigen::MatrixXd::Ones(4,40);
    Eigen::Matrix<double,3,3> Rstand;
    Eigen::Matrix<double,3,3> Jstand;
    Eigen::Matrix<double,3,3> Jinvstand;
    int16_t contactindex = 0;
    Eigen::Matrix<double,3,1> disturbforce = {0,0,0};
    Eigen::Matrix<double, 24, 1> opt_HLsol_;

    casadi::MX contact_sequence = casadi::MX::ones(4,40);//,{0});
    casadi::DM contact_sequence_dm = casadi::DM::ones(4,40);//,{0});
    casadi::DM gravityN = {0,0,9.81};
    casadi::DM Jstandcasadi = casadi::DM::zeros(3,3);
    casadi::DM Jinvcasadi = casadi::DM::zeros(3,3);
    casadi::DM Raibheur = 0;//0.5*Tstance*desVel(0);
    int RaibMult = 2;//8;//8;
    casadi::DM previous_sol = casadi::DM::zeros(NFS*(HORIZ+1)+NFI*HORIZ,1);
    double localvelocity = 0;
    casadi::DM x_dom_init = 0;

    //int forcefitorder = 4;
    Eigen::MatrixXd forceCoeff;
    //Eigen::Matrix<double, 12, 5> forceCoeff = Eigen::MatrixXd::Zero(12,5);
    Eigen::Matrix<double, HORIZ,1> forcefitx = Eigen::MatrixXd::Zero(HORIZ,1);
    Eigen::Matrix<double,12,1> OPTforce = Eigen::Matrix<double,12,1>::Zero();
    casadi::DM Raibstep = 0;//0.5*Tstance*desVel(0);
    casadi::DM vRaibstep = 0;//0.5*Tstance*desVel(0);
    int velsign = 1;

    public:
    casadi::DM previousp = casadi::DM::zeros(NFS*(HORIZ+1)+NFI*(HORIZ)+4*(HORIZ+2),1);
    casadi::DM optstephll = casadi::DM::zeros(4,1);
    casadi::DM front_off = 0.1;
    casadi::DM rear_off = -0.1; 
    casadi::DM pitch_ref = 0;
    int Raibflag = 0;
};

#endif