#ifndef DATALOG_HPP
#define DATALOG_HPP

#include "/home/taizoon/raisimEnv/raisimWorkspace/A1_LL_Exp-Arch_Change/global_include/global_loco_structs.hpp"
#include "iostream"
#include "fstream"


class DataLog{
public:
    DataLog(std::string filename);
    virtual ~DataLog();
    void openFile();
    void writeData(const StateInfo *state, const VCInfo *vc, const ContactInfo *con, 
                    const TrajInfo *traj, const LLInfo *LL, const KinematicsInfo *kin, const size_t ctrlTick, 
                    const int force[4], Eigen::Matrix<double, 24, 1>& opt_HL_state, size_t locoTick, double locoPhase, double flightPhaseRL, double flightPhaseLR, double DDiff, double HDiff, Eigen::Matrix<double, 5, 1> NLstep, size_t duration);

private:
    std::fstream fid;
    Eigen::Matrix<double, 12, 1> y_, dy_, hd_, dhd_, ddhd_;

    int success = 0;

};



#endif
