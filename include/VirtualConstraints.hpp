#ifndef VIRTUALCONSTRAINTS
#define VIRTUALCONSTRAINTS

//#include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/global_include/global_loco_structs.hpp"
//#include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/util_include/Bezier.h"
#include "../global_include/global_loco_structs.hpp"
#include "../util_include/Bezier.h"

using KinInf = KinematicsInfo;
using ConInf = ContactInfo;
using MP = Settings::Motion_params;

class VirtualConstraints{
public:
    VirtualConstraints();
    virtual ~VirtualConstraints(){};

    void updateVirtualConstraints(const StateInfo *state, const KinInf *kin, const TrajInfo *traj, const ConInf *con, size_t gait, double phaseVar, MP *params, const LLInfo *ll);
    void updateTime();
    void setDesiredForce(Eigen::Matrix<double, 12, 1> fDes){VC.fDes = fDes;};
    const VCInfo* getVCPointer(){return &VC;};

    void updateVirtualConstraintsflight(const StateInfo *state, const KinInf *kin, const TrajInfo *traj, const ConInf *con, size_t gait, 
                                            double phaseVar, double phaseVarLR, double stepLenRL[3], double stepLenLR[3], MP *params, const LLInfo *ll);

    void updateVirtualConstraintswalk(const StateInfo *state, const KinInf *kin, const TrajInfo *traj, const ConInf *con, size_t gait, 
                                            double frontphase, double rearphase, MP *params, const LLInfo *ll);                                        

private:
    VCInfo VC; // Needed by the low-level controller
    double h_sw = 0.08;
    Eigen::MatrixXd h0, dh0;
    double wall_y = 0.25;
    double two_contact = 20;
};




#endif
