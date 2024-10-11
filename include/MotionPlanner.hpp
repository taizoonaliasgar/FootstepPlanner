#ifndef PLANNER_HPP
#define PLANNER_HPP

#include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/global_include/global_loco_structs.hpp"
#include "ContactEst.hpp"
#include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/util_include/Transforms.hpp"
#include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/util_include/Bezier.h"
#include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/util_include/EigenUtils.hpp"

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
    

    //Taizoon changes
    void savesteplen(ContactEst *con_obj);
    void updatesteplen(ContactEst *con_obj);
    //NMPC footstep planner
    void setStep_NMPC(Eigen::Matrix<double,5,1> NLstep,double vdes, const StateInfo *state, MP * params, double phase);


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
};

#endif
