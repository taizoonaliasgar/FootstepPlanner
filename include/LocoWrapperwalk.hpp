#ifndef LOCO_WRAPPERWALK
#define LOCO_WRAPPERWALK

// #include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/global_include/global_loco_structs.hpp"
// #include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/params/Parameters.hpp"
// #include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/dynamics/RobotModel.hpp"
#include "../global_include/global_loco_structs.hpp"
#include "../params/Parameters.hpp"
#include "../LL_integrator/dynamics/RobotModel.hpp"
#ifdef USE_OSQP
#include "LowLevelCtrl_OSQP.hpp"
#else
//#include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/include/LowLevelCtrl.hpp"
#include "LowLevelCtrl.hpp"
#endif
#include "VirtualConstraints.hpp"
#include "ContactEst.hpp"
#include "MotionPlanner.hpp"
#include "DataLog.hpp"

#include <memory>

class LocoWrapperwalk : public Parameters
{
public:
    LocoWrapperwalk(int argc, char *argv[]);
    virtual ~LocoWrapperwalk();

    void calcTau(const double q[18], const double dq[18], const double R[9], const int force[4], size_t gait, size_t ctrlTick, size_t duration);//, Eigen::Matrix<double, 24, 1>& opt_HL_state);
    double* getTorque(){return LL->getTorque();};
    const int* getConDes(){return con->des;};
    Eigen::Matrix<double, 18, 1> getJointPosCmd(){return ll->q;};
    Eigen::Matrix<double, 18, 1> getJointVelCmd(){return ll->dq;};
    void initStandVars(Eigen::Matrix<double,3,1> com, double yaw, double standTime){ PP->updateStandVars(com,yaw,standTime);};
    void updateDesiredForce(Eigen::Matrix<double, 12, 1> fDes){VC->setDesiredForce(fDes);};
    void updateVel(const float vel[3]){PP->setVel(vel);};
    void updatePose(const float pose[6]){PP->setPose(pose);};
    void updatePoseType(size_t poseType_){PP->setPoseType(poseType_);};
    
    

    //Taizoon changes
    Eigen::Matrix<double, 12, 1> getQPForce(){return ll->QP_force;};
    Eigen::Matrix<double, 3, 4> getfootposition(){return kin->toePos;};
    Eigen::Matrix<double, 3, 4> gethipposition(){return kin->hipPos;};
    const int* getffootstate(){return con->ind;};
    void setoptHLstate(Eigen::Matrix<double, 24, 1> HLopt){opt_HLstate = HLopt;};
    int stanceflag(){return stance_phase_flag;};
    int stancecounter(){return walk_tick;};
    void updatestate(const double q[18], const double dq[18], const double R[9]){quad->updateState(q,dq,R);};
    void setRaisimD(Eigen::Matrix<double, 18, 18> DRaisim){DRai = DRaisim;};
    void setRaisimH(Eigen::Matrix<double, 18, 1> HRaisim){HRai = HRaisim;};

    //NMPC
    void setoptNLstate(Eigen::Matrix<double, 33, 1> HLopt);
    void setcontactconfig(int controlMPC);

    // Pointers to structs
    const StateInfo *state;
    const DynamicsInfo *dyn;
    const KinematicsInfo *kin;
    const ContactInfo *con;
    const TrajInfo *traj;
    const VCInfo *vcon;
    const LLInfo *ll;

private:
    // size_t newDom = 0;
    size_t locoTick;
    double phaseVar = 0;
    double maxPhase;// = 0.96;
    size_t gaitTemp = STAND;
    size_t forceDomainChange = 0;

    //HLstate
    Eigen::Matrix<double, 24, 1> opt_HLstate;

    //Flying trot
    size_t flightTickRL = 0;
    size_t flightTickLR = 0;
    double flphase = 0;
    double rlphase = 0;
    int stance_phase_flag = 1;
    size_t second_jump_flag = 0;

    // Pointers to class objects
    std::unique_ptr<DataLog> data;
    RobotModel *quad;
    LowLevelCtrl *LL;
    VirtualConstraints *VC;
    ContactEst *conEst;
    MotionPlanner *PP;

    //flying trot changes
    double phaseVarRL = 0;
    double phaseVarLR = 0;
    double stepLenRL[3] = {0};
    double stepLenLR[3] = {0};
    std::vector<int> nextContact;
    double takeoffheight = 0.25;

    double mass_ = 12.453;
    double g_ = 9.81;

    //Upright walk
    int walk_tick = 0;
    std::vector<int> desired_contact{1,1,1,1};
    Eigen::Matrix<double,4,40> contact_horizon = Eigen::MatrixXd::Ones(4,40);

    Eigen::Matrix<double,18,18> DRai;// = Eigen::MatrixXd::Zero(18,18);
    Eigen::Matrix<double,18,1> HRai;// = Eigen::MatrixXd::Zero(18,1);
    double Ddiff;
    double Hdiff;
    double two_contact;
    Eigen::Matrix<double,12,1> z = Eigen::MatrixXd::Zero(12,1);
    Eigen::Matrix<double,12,12> Ki = Eigen::MatrixXd::Zero(12,12);

    //NLState
    Eigen::Matrix<double, 5, 1> NLstep = Eigen::MatrixXd::Zero(5,1);

};

inline double getPhase(double time, double time_0, double time_f){
    return (1.0*time-1.0*time_0)/(1.0*time_f-1.0*time_0);
};

#endif
