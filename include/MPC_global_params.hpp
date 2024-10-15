#ifndef MPC_GLOBAL_PARAMS
#define MPC_GLOBAL_PARAMS


// #ifndef DFLOAT
//     using osqp_mat = Eigen::MatrixXd;
//     using osqp_arr = Eigen::ArrayXd;
// #else
//     using osqp_mat = Eigen::MatrixXf;
//     using osqp_arr = Eigen::ArrayXf;
// #endif

#define NRS 12
#define NRI 12
#define HORIZ 5
#define TINI 10

#define MPC_Hz 100
#define MPC_dt ( (double)1.0 / MPC_Hz )
#define RAND 8

#define MASS 12.453
// #define MASS 13.75
#define GRAV 9.81

#define Tstance 0.2
#define simdt 0.001
#define dt 0.005
#define DOMAINSTEPS 20

#define NFS 16
#define NFI 16

#define NFSR 12
#define NFIR 12

#endif