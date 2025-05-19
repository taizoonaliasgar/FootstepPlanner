#ifndef LOCO_WRAPPER
#define LOCO_WRAPPER

#include "global_loco_structs.hpp"
#include "Parameters.hpp"
#include "RobotModel.hpp"
#ifdef USE_OSQP
#include "LowLevelCtrl_OSQP.hpp"
#else
#include "LowLevelCtrl.hpp"
#endif
#include "VirtualConstraints.hpp"
#include "ContactEst.hpp"
#include "MotionPlanner.hpp"
#include "DataLog.hpp"

#include <memory>

class LocoWrapper : public Parameters
{
public:
    LocoWrapper(int argc, char *argv[]);
    virtual ~LocoWrapper();

    void calcTau(const double q[18], const double dq[18], const double R[9], const int force[4], size_t gait, size_t ctrlTick, size_t duration);
    double* getTorque(){return LL->getTorque();};
    const int* getConDes(){return con->des;};
    Eigen::Matrix<double, 18, 1> getJointPosCmd(){return ll->q;};
    Eigen::Matrix<double, 18, 1> getJointVelCmd(){return ll->dq;};
    void initStandVars(Eigen::Matrix<double,3,1> com, double yaw, double standTime){ PP->updateStandVars(com,yaw,standTime);};
    void updateDesiredForce(Eigen::Matrix<double, 12, 1> fDes){VC->setDesiredForce(fDes);};
    void updateVel(const float vel[3]){PP->setVel(vel);};
    void updatePose(const float pose[6]){PP->setPose(pose);};
    void updatePoseType(size_t poseType_){PP->setPoseType(poseType_);};

    //Taizoon functions to make upright
    void calcTau2(const double q[18], const double dq[18], const double R[9], size_t gait, size_t ctrlTick,size_t solveduration);
    void setcontactconfig(int controlMPC);
    void plottingfoothd(const VCInfo *vc, const ContactInfo *con);
    //Eigen::Matrix<double, 3, 4> getfootposition(){return kin->toePos;};
    void getshiftedCoM(Eigen::Matrix<double, 4, 1> footweight);
    void setshiftedCoM(){PP->setshiftedCoM(CoMnew);};
    void setswingContact(Eigen::Matrix<double,4,1>nContact){nextContact[0] = nContact(0);nextContact[1] = nContact(1);nextContact[2] = nContact(2);nextContact[3] = nContact(3);};
    void tookfirststep(){LL->afterfirststep();};
    void stepsonwall(int stepind){wallstep = stepind;};
    void stopclimbing(){stopclimb = true;};
    void incstep(){PP->increasesteplenth();};
    void setrearhippose();
    void setfinalCoM();
    void settlesteps(int setsteps){settlestep = setsteps;};
    void gotfinalstate(){atfinalstate = true;};
    void setxzsteplength(size_t movetime);

    //NMPC upright walk
    int stancecounter(){return walk_tick;};
    void setoptNLstate(Eigen::Matrix<double, 33, 1> HLopt);
    void setRFfalse(){RaibFlag = false;};
    Eigen::Matrix<double, 3, 4> getfootposition(){return kin->toePos;};
    Eigen::Matrix<double, 3, 4> gethipposition(){return kin->hipPos;};
    void updatestate(const double q[18], const double dq[18], const double R[9]){quad->updateState(q,dq,R);};
    Eigen::Matrix<double, 12, 1> getStateEstimate(double jointPos[18], Eigen::VectorXd jointVelTotal, Eigen::Matrix<double, 3, 1> imu_eul, Eigen::Matrix<double, 3, 1> imu_omega);

    //Stitching together
    void readytowalk(){quad->fullyupright();LL->fullyupright();};
    void setfinalCoM2(int settlesteps);
    void startwalking(){readytowalkf = true;LL->fullyupright();LL->enablehipposcntrl();};
    Eigen::Matrix<double, 12, 1> getpreviousQPforce(){return ll->QP_force;};
    void readytoreallywalk(){LL->keepwalking();};
    void getStateEstimatefull(double q[18], double dq[18], const int* contact, Eigen::Matrix<double,3,3> R, int robotdown, bool dynswitch,size_t ctrlTick);

    //Hardware
    void ExpWrapper(const double jpos_est[18], const double jvel_est[18], const double rotMatrixDouble[9], size_t control_Tick, size_t solveduration, 
                                    int HLContactIndex[5], Eigen::Matrix<double, 12, 1> comDes, Eigen::Matrix<double, 17, 1> fDes);
    void setcontactconfigExp(int HLContactIndex[5]);
    void setoptNLstateExp(Eigen::Matrix<double, 12, 1> comDes, Eigen::Matrix<double, 17, 1> fDes);
    Eigen::Matrix<double, 3, 4> gettoe_prev(){return toe_prev;};
    void settoe_prev(){toe_prev = kin->toePos;};

    //Multi-threading
    void setHLphase(int HLContact5){HLMTphase = HLContact5;};
    //Managing the knee
    void ExpWrapperk(const double jpos_est[18], const double jvel_est[18], const double rotMatrixDouble[9], size_t control_Tick, size_t solveduration, 
        int HLContactIndex[5], Eigen::Matrix<double, 12, 1> comDes, Eigen::Matrix<double, 17, 1> fDes);
    void calcTau2k(const double q[18], const double dq[18], const double R[9], size_t gait, size_t ctrlTick, size_t solveduration);
    void setfinalCoM3(Eigen::Matrix<double, 4, 1> footweight, size_t setsteps);
    void tooksecondstep(){LL->aftersecondstep();};
    
    //Logging IMU
    void setIMUdata(double att_euler[3]){data->setIMUeuler(att_euler);};

    //Switching Kp
    void switchKp(){LL->switchKp();};

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
    size_t locoTick = 0;
    double phaseVar = 0;
    double maxPhase = 0.996;
    size_t gaitTemp = STAND;
    size_t forceDomainChange = 0;

    // Pointers to class objects
    std::unique_ptr<DataLog> data;
    RobotModel *quad;
    LowLevelCtrl *LL;
    VirtualConstraints *VC;
    ContactEst *conEst;
    MotionPlanner *PP;

    //Taizoon Changes to go upright
    Eigen::Matrix<double, 24, 1> opt_HLstate = Eigen::MatrixXd::Zero(24,1);
    Eigen::Matrix<double, 5, 1> NLstep = Eigen::MatrixXd::Zero(5,1);
    Eigen::Matrix<double,4,40> contact_horizon = Eigen::MatrixXd::Ones(4,40);
    std::vector<int> desired_contact{1,1,1,1};
    double flphase = 0.0;
    double rlphase = 0.0;
    Eigen::Matrix<double,18,1> Hr = Eigen::MatrixXd::Zero(18,1);
    Eigen::Matrix<double,12,1> z = Eigen::MatrixXd::Zero(12,1); 
    Eigen::Matrix<double,12,12> Ki = Eigen::MatrixXd::Zero(12,12);
    Eigen::Matrix<double,4,1> CoMnew = Eigen::MatrixXd::Zero(4,1);
    bool CoMshifted = false;
    std::vector<int> nextContact = {1,1,1,1};
    int wallstep = 0;
    bool stopclimb = false;
    int settlestep = 0;
    bool atfinalstate = false;

    //NMPC upright walk
    int walk_tick = 0;
    bool RaibFlag = true;
    Eigen::Matrix<double,18,18> DRai = Eigen::MatrixXd::Zero(18,18);
    Eigen::Matrix<double,18,1> HRai = Eigen::MatrixXd::Zero(18,1);

    bool readytowalkf = false;
    //double switchtime = 24;

    //Estimator stuff
    double yzdot_thresh = 0.3;
    double xdot_thresh = 0.3;
    int rearfootweight = 4;
    double yzdot_thresh2 = 0.8;
    double xdot_thresh2 = 0.5;
    
    //Hardware params
    size_t settling_e = 0.2*ctrlHz;                   // Settling down
    size_t duration_e = 1.8*ctrlHz;                   // Stand up 
    size_t loco_start_e = settling_e + duration_e;        // Start the locomotion pattern

    size_t shifttime = 1*ctrlHz;//1
    size_t movetime = 0.3*ctrlHz;
    size_t shifttime2 = 1.0*ctrlHz;
    size_t movetime2 = 0.3*ctrlHz;
    size_t movetime3 = 0.2*ctrlHz;
    size_t switchtime = 20;

    Eigen::Matrix<double,4,1> nextcon_e = Eigen::MatrixXd::Ones(4,1);
    int maxsteps = 7;//14;
    int settlingsteps = 6;//4;
    double rearweight = 6;
    int stepind_e = 0;
    int stepind2_e = 0;
    double *tau_LL;
    const int force_LL[4] = {0,0,0,0};
    Eigen::Matrix<double, 3, 4> toe_prev = Eigen::MatrixXd::Zero(3,4);

    //Multi-threading
    int HLMTphase = 0;

    int standoffset = 0;
    int currentshift = 0;

    //Settling down
    double x00 = 0;
    double z00 = 0;
    double setdx = 0;
    double setdz = 0;    
};

inline double getPhase(double time, double time_0, double time_f){
    return (1.0*time-1.0*time_0)/(1.0*time_f-1.0*time_0);
};

#endif
