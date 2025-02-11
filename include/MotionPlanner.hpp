#ifndef PLANNER_HPP
#define PLANNER_HPP

// #include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/global_include/global_loco_structs.hpp"
// #include "ContactEst.hpp"
// #include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/util_include/Transforms.hpp"
// #include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/util_include/Bezier.h"
// #include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/util_include/EigenUtils.hpp"

#include "../global_include/global_loco_structs.hpp"
#include "ContactEst.hpp"
#include "../util_include/Transforms.hpp"
#include "../util_include/Bezier.h"
#include "../util_include/EigenUtils.hpp"

#define MAX_SL_F_X 0.16     // Max forward step length (magnitude)
#define MAX_SL_R_X 0.16     // Max backward step length (magnitude)
#define MAX_SL_Y 0.12       // Max lateral step length (magnitude)  

using MP = Settings::Motion_params;

class MotionPlanner
{
public:
    MotionPlanner();
    virtual ~MotionPlanner(){};

    void planTraj(const StateInfo *state, const KinematicsInfo *kin, ContactEst *con_obj, size_t gait, double phase, size_t ctrlTick, MP *params, Eigen::Matrix<double, 24, 1>& opt_HLstate, Eigen::Matrix<double, 5, 1>& NLstep);
    void updateStandVars(const Eigen::Matrix<double,3,1> &com, double yaw, double timeToStand);
    void updateVel(Eigen::Matrix<double,3,1> &desVel, Eigen::Matrix<double, 3, 1> &desOmega, MP *params);
    void setStep_Raibert(const StateInfo *state, double domLenSec, const Eigen::Matrix<double,3,1> &desVel, std::vector<double> KP);
    void setVel(const float vel[3]){desVel_(0)=vel[0]; desVel_(1)=vel[1]; desOmega_(2)=vel[2];};
    void setPose(const float pose[6]){for(int i=0;i<6;++i){pose_(i)=pose[i];};};
    void setPoseType(size_t poseType_){poseType = poseType_;};
    const TrajInfo* getTrajInfoPointer(){return &traj;};
    void setToeInit(const KinematicsInfo *kin){traj.toeInit = kin->toePos;};
    

    //Taizoon changes
    void savesteplen(ContactEst *con_obj);
    void updatesteplen(ContactEst *con_obj);
    //NMPC footstep planner
    void setStep_NMPC(Eigen::Matrix<double,5,1> NLstep,double vdes, const StateInfo *state, MP * params, double phase);
    void setFoot(const KinematicsInfo *kin);
    bool getReachedWall(){return reachedWall;};
    void setshiftedCoM(Eigen::Matrix<double,4,1> CoMnew){xnew = CoMnew(0); ynew = CoMnew(1); znew = CoMnew(2); pitchnew = CoMnew(3);};
    void shiftCoM(ContactEst *con_obj, double phase, size_t shifttime);
    void movefoot(size_t movetime, size_t wallsteps);//, size_t wallstep);
    void movefoot2(size_t movetime, double phase);
    void setx0y0z0(double x0_, double y0_, double z0_, double p0_){x0 = x0_; y0 = y0_; z0 = z0_;p0 = p0_;};
    void datalogger(size_t ctrlTick){std::cout << ctrlTick << "\t" << x0 << "\t" << y0 << "\t" << z0 << "\t" << xnew << "\t" << ynew << "\t" << znew << std::endl;};
    void increasesteplenth(){if(upstep<-0.01){upstep = upstep+0.04;}else{upstep = -0.01;};};
    void setrearhip(double x, double y, double z){rhip_x = x; rhip_y = y; rhip_z = z;};
    void shiftCoM2(ContactEst *con_obj, double phase, size_t shifttime, bool maxsteps);
    void movefoot2(size_t movetime, Eigen::Matrix<double,4,1> xzsteps);
    void shiftCoM3(ContactEst *con_obj, double phase, size_t shifttime, bool maxsteps);
    void movefoot3(size_t movetime);

protected:
    inline void setStepLen(double x, double y, double z){
        traj.stepLen[0] = x; traj.stepLen[1] = y; traj.stepLen[2] = z;
    }

private:
    double x0;
    double xf;
    double y0;
    double yf;
    double z0;
    double standTime;
    TrajInfo traj;
    size_t poseType = POSE_CMD;

    Eigen::Matrix<double, 3, 1> desVel_;
    Eigen::Matrix<double, 3, 1> desOmega_;
    Eigen::Matrix<double, 6, 1> pose_;
    double yawOffset;

    double stepLenRL[3] = {0.0,0.0,0.0};
    double stepLenLR[3] = {0.0,0.0,0.0};
    bool reachedWall = false;
    double xnew = 0;
    double ynew = 0;
    double znew = 0;
    double pitchnew = 0;
    double upstep = -0.2;
    double rhip_x = 0;
    double rhip_y = 0;
    double rhip_z = 0;
    double p0 = 0;
    double minpitch = -1.2;
    
};

#endif
