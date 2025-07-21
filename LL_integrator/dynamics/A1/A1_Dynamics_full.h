#ifndef __A1_DYNAMICS_H__
#define __A1_DYNAMICS_H__

#define MAX(X,Y)  ((X) < (Y) ? (Y) : (X))
#define MIN(X,Y)  ((X) > (Y) ? (Y) : (X))

#include "math.h"
#include "dynamicsSupportFunctions.h"

// Dynamics
void D_mat(double *p_output1,const double *var1);
void C_vec(double *p_output1,const double *var1,const double *var2);
void B_mat(double *p_output1,const double *var1);
void G_vec(double *p_output1,const double *var1);

// Forward Kinematics
void FK_FL_toe(double *p_output1,const double *var1);
void FK_FR_toe(double *p_output1,const double *var1);
void FK_RL_toe(double *p_output1,const double *var1);
void FK_RR_toe(double *p_output1,const double *var1);

void FK_FL_hip(double *p_output1,const double *var1);
void FK_FR_hip(double *p_output1,const double *var1);
void FK_RL_hip(double *p_output1,const double *var1);
void FK_RR_hip(double *p_output1,const double *var1);

// Jacobians
void J_FL_toe(double *p_output1,const double *var1);
void J_FR_toe(double *p_output1,const double *var1);
void J_RL_toe(double *p_output1,const double *var1);
void J_RR_toe(double *p_output1,const double *var1);

void J_FL_hip(double *p_output1,const double *var1);
void J_FR_hip(double *p_output1,const double *var1);
void J_RL_hip(double *p_output1,const double *var1);
void J_RR_hip(double *p_output1,const double *var1);

// Jacobian Dot
void dJ_FL_toe(double *p_output1,const double *var1,const double *var2);
void dJ_FR_toe(double *p_output1,const double *var1,const double *var2);
void dJ_RL_toe(double *p_output1,const double *var1,const double *var2);
void dJ_RR_toe(double *p_output1,const double *var1,const double *var2);

void dJ_FL_hip(double *p_output1,const double *var1,const double *var2);
void dJ_FR_hip(double *p_output1,const double *var1,const double *var2);
void dJ_RL_hip(double *p_output1,const double *var1,const double *var2);
void dJ_RR_hip(double *p_output1,const double *var1,const double *var2);


//For the upright configuration
void D_mat_u(double *p_output1,const double *var1);
void C_vec_u(double *p_output1,const double *var1,const double *var2);
void B_mat_u(double *p_output1,const double *var1);
void G_vec_u(double *p_output1,const double *var1);

// Forward Kinematics
void FK_FL_toe_u(double *p_output1,const double *var1);
void FK_FR_toe_u(double *p_output1,const double *var1);
void FK_RL_toe_u(double *p_output1,const double *var1);
void FK_RR_toe_u(double *p_output1,const double *var1);

void FK_FL_hip_u(double *p_output1,const double *var1);
void FK_FR_hip_u(double *p_output1,const double *var1);
void FK_RL_hip_u(double *p_output1,const double *var1);
void FK_RR_hip_u(double *p_output1,const double *var1);

// Jacobians
void J_FL_toe_u(double *p_output1,const double *var1);
void J_FR_toe_u(double *p_output1,const double *var1);
void J_RL_toe_u(double *p_output1,const double *var1);
void J_RR_toe_u(double *p_output1,const double *var1);

void J_FL_hip_u(double *p_output1,const double *var1);
void J_FR_hip_u(double *p_output1,const double *var1);
void J_RL_hip_u(double *p_output1,const double *var1);
void J_RR_hip_u(double *p_output1,const double *var1);

// Jacobian Dot
void dJ_FL_toe_u(double *p_output1,const double *var1,const double *var2);
void dJ_FR_toe_u(double *p_output1,const double *var1,const double *var2);
void dJ_RL_toe_u(double *p_output1,const double *var1,const double *var2);
void dJ_RR_toe_u(double *p_output1,const double *var1,const double *var2);

void dJ_FL_hip_u(double *p_output1,const double *var1,const double *var2);
void dJ_FR_hip_u(double *p_output1,const double *var1,const double *var2);
void dJ_RL_hip_u(double *p_output1,const double *var1,const double *var2);
void dJ_RR_hip_u(double *p_output1,const double *var1,const double *var2);

// Hip correction
void FK_RR_hip2(double *p_output1,const double *var1);
void FK_RL_hip2(double *p_output1,const double *var1);
void FK_FR_hip2(double *p_output1,const double *var1);
void FK_FL_hip2(double *p_output1,const double *var1);

void J_RR_hip2(double *p_output1,const double *var1);
void J_RL_hip2(double *p_output1,const double *var1);
void J_FR_hip2(double *p_output1,const double *var1);
void J_FL_hip2(double *p_output1,const double *var1);

void dJ_RR_hip2(double *p_output1,const double *var1,const double *var2);
void dJ_RL_hip2(double *p_output1,const double *var1,const double *var2);
void dJ_FR_hip2(double *p_output1,const double *var1,const double *var2);
void dJ_FL_hip2(double *p_output1,const double *var1,const double *var2);


#endif
