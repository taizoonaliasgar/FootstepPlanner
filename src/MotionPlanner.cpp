//
// Author: Randy Fawcett on 12/2021.
//
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "MotionPlanner.hpp"

MotionPlanner::MotionPlanner(){
    x0 = 0;
    y0 = 0;
    z0 = 0.05;
    
    traj.domLen = 1*ctrlHz;
    standTime = 1*ctrlHz;

    traj.comDes.setZero();
    traj.redDes.setZero();
    traj.toeInit.setZero();
    traj.toeFinal.setZero();
    traj.toeOffset[2] = Z_TOE_OFFSET;

    desVel_.setZero();
    desOmega_.setZero();
    pose_.setZero();

    yawOffset = 0;
}

void MotionPlanner::planTraj(const StateInfo *state, const KinematicsInfo *kin, ContactEst *con_obj, size_t gait, double phase, size_t ctrlTick, MP * params, Eigen::Matrix<double, 24, 1>& opt_HLstate, Eigen::Matrix<double, 5, 1>& NLstep){
    const ContactInfo* con = con_obj->getConInfoPointer();   

    static Eigen::Matrix<double, 3, 1> desVel = {0,0,0};
    static Eigen::Matrix<double, 3, 1> desOmega = {0,0,0};

    // ****************************************** //
    // add yaw as class variable and initialize it
    // during updateStandVars
    // ****************************************** //

    // static double yaw = state->q(5);
    static bool standTrigger = false;
    double normMax = 0.05;
    if(gait==STAND){
        double s = (phase>1) ? 1 : ((phase<0) ? 0 : phase);

        double xFinal = x0-0.04*0;
        double yFinal = y0-0.04*0;
        double zFinal = params->standHeight;
        double alpha_x[8] = { x0,x0,x0,
                 x0+(xFinal-x0)/4,
                 x0+3*(xFinal-x0)/4,
                 xFinal,xFinal,xFinal};
        double alpha_y[8] = {y0,y0,y0,
                 y0+(yFinal-y0)/4,
                 y0+3*(yFinal-y0)/4,
                 yFinal,yFinal,yFinal};
        double alpha_z[8] = {z0,z0,z0,
                 z0+(zFinal-z0)/4,
                 z0+3*(zFinal-z0)/4,
                 zFinal,zFinal,zFinal};

        double traj_x[3], traj_y[3], traj_z[3];
        calcBezierAll((int)8, alpha_x, s, traj_x);
        calcBezierAll((int)8, alpha_y, s, traj_y);
        calcBezierAll((int)8, alpha_z, s, traj_z);

        // traj.comDes -> pos, vel, theta, omega
        traj.comDes.block(0,0,3,1) << traj_x[0], traj_y[0], traj_z[0];
        traj.comDes.block(3,0,3,1) << traj_x[1], traj_y[1], traj_z[1];
        traj.comDes.block(6,0,3,1) << 0, 0, yawOffset;
        traj.comDes.block(9,0,3,1) << 0, 0, 0;

        traj.domLen = standTime;
        con_obj->setDesDomain({1,1,1,1});

        traj.toeInit = kin->toePos; 
        traj.toeFinal = kin->toePos;
        traj.redDes = traj.comDes;
        
    }else if(gait==POSE) {
        traj.comDes.block(3,0,3,1) << 0,0,0;
        traj.comDes.block(9,0,3,1) << 0,0,0;
        double t = 1.0*ctrlTick/(1.0*ctrlHz);
        static double t_init = t;
        if (poseType!=POSE_CMD){
            static Eigen::Matrix<double, 12, 1> lock = traj.comDes;
            traj.comDes = lock;
        }
        if(poseType==POSE_X){
			double freq = 0.8*MY_PI;
            double mag = 0.04;
            traj.comDes(0) += mag*sin(freq*(t-t_init));
            traj.comDes(3) += mag*freq*cos(freq*(t-t_init));
        }
        else if(poseType==POSE_Y){
            double freq = 0.8*MY_PI;
            double mag = 0.04;
            traj.comDes(1) += mag*sin(freq*(t-t_init));
            traj.comDes(4) += mag*freq*cos(freq*(t-t_init));
        }
        else if(poseType==POSE_Z){
            double freq = 0.8*MY_PI;
            double mag = 0.05;
            traj.comDes(2) += mag*cos(freq*(t-t_init))-mag;
            traj.comDes(5) += -1.0*mag*freq*sin(freq*(t-t_init));
        }
        else if(poseType==POSE_ROLL){
            double freq = 0.8*MY_PI;
			double mag = 0.3491;
            traj.comDes(6) += mag*sin(freq*(t-t_init));
            traj.comDes(9) += mag*freq*cos(freq*(t-t_init));
        }
        else if(poseType==POSE_PITCH){
            double freq = 0.8*MY_PI;
            double mag = 0.17453;
            traj.comDes(7) += mag*sin(freq*(t-t_init));
            traj.comDes(10) += mag*freq*cos(freq*(t-t_init));
        }
        else if(poseType==POSE_YAW){
            double freq = 0.8*MY_PI;
            double mag = 0.13963;
            traj.comDes(8) += mag*sin(freq*(t-t_init));
            traj.comDes(11) += mag*freq*cos(freq*(t-t_init));
        }
        else if(poseType==POSE_COMB){
        	double freq = 0.6*MY_PI;
        	double mag = 0.3491;
        	traj.comDes(7) += mag*sin(freq*(t-t_init));
        	traj.comDes(10) += mag*freq*cos(freq*(t-t_init));
        	
            static int triggerStart = 0;
            if (triggerStart || cos(freq*(t-t_init))<0){
                traj.comDes(8) += mag*cos(freq*(t-t_init));
                traj.comDes(11) += -mag*freq*sin(freq*(t-t_init));
                triggerStart = 1;
            }
        }
        else if(poseType==POSE_CMD){
            // pose_ = {delta x, delta y, delta z, roll, pitch, delta yaw}
            standTrigger = true;
            Eigen::Matrix<double, 3, 1> pose; 
            Eigen::Matrix<double, 4, 1> quatyaw;
            Eigen::Matrix<double, 4, 1> quatBody;
            Eigen::Matrix<double, 4, 1> quatWorld;
            XYZ_to_quat(0,0,yawOffset,quatyaw);
            XYZ_to_quat(pose_(3),pose_(4),pose_(5), quatBody);
            quatMult(quatyaw, quatBody, quatWorld);
            quat_to_XYZ(quatWorld,pose);
            traj.comDes(2) = params->standHeight+pose_(2);
            traj.comDes.block(6,0,3,1) = pose;
            traj.comDes.block(3,0,3,1) << 0,0,0;
            traj.comDes.block(9,0,3,1) << 0,0,0;
        }
    }else if(gait==TAP){
        setStepLen(0.0,0.0,0.0);
        if(con->changeDomain){ 
        	static double domLenSec = 1.5;
            con_obj->setDesDomain({1, 0, 1, 1});
            traj.domLen = domLenSec*ctrlHz;
            traj.toeInit = kin->toePos;
//            domLenSec -= (domLenSec>1.0) ? 0.25 : 0;
        }
    }else if(gait==INPLACE_WALK){
        if(con->changeDomain==1){
            traj.toeInit = kin->toePos;
            
            static int n = 0;
            Eigen::Matrix<int, 4, 4> doms;
            doms << 0,1,1,1,
                    1,0,1,1,
                    1,1,1,0,
                    1,1,0,1;
            n = (++n) % 4;
            con_obj->setDesDomain({doms(n,0),doms(n,1),doms(n,2),doms(n,3)});
            traj.domLen = 0.3*ctrlHz;
        }
    }else if(gait==INPLACE_TROT){
        if(con->changeDomain){
            if(con->des[0]==1){
                con_obj->setDesDomain({0, 1, 1, 0});}
            else{
                con_obj->setDesDomain({1, 0, 0, 1});}
            traj.domLen = 0.16*ctrlHz;
            traj.toeInit.block(0,0,2,4) = kin->hipPos.block(0,0,2,4);
            traj.toeInit.block(2,0,1,4) = kin->toePos.block(2,0,1,4);
        }
    }else if(gait==WALK){
        double dt = (1.0/LL_Hz);
        double domLenSec = 0.2;

        if(con->changeDomain==1){
            updateVel(desVel,desOmega,params);
            
            // ================================================ //
            // Update contact matrix and initial toe position(s)
            // ================================================ //
            static int n = 2;
            Eigen::Matrix<int, 4, 4> doms;
            doms << 1,1,0,1,
                    0,1,1,1,
                    1,1,1,0,
                    1,0,1,1;
            n = (++n) % 4;
            con_obj->setDesDomain({doms(n,0),doms(n,1),doms(n,2),doms(n,3)});
            traj.domLen = domLenSec*ctrlHz;
            traj.toeInit = kin->toePos;

            // ================================================ //
            // Marc Raibert foothold selection 
            // ================================================ //
            std::vector<double> KP = {0.01,0.0,0.0};
            setStep_Raibert(state,3*domLenSec,desVel,KP);
        }
    }else if(gait==TROT){
        double dt = (1.0/LL_Hz);
        double domLenSec = 0.2;

        //updateVel(desVel, desOmega, params);
        static int startTrot = 0;

        //if(startTrot<1500){
        //    opt_HLstate.block(15,0,3,1) << 0,0,0;
        //    startTrot+=1;
        //}

        if(con->changeDomain==1){

        	if(startTrot>5){
               //opt_HLstate.block(15,0,3,1) << 0,0,0;
               updateVel(desVel, desOmega, params);
            }
            startTrot+=1;
            // Eigen::Vector3d desVelBody = {state->R(0,0)*opt_HLstate(15,0)+state->R(1,0)*opt_HLstate(16,0)+state->R(2,0)*opt_HLstate(17,0),
            //                                  state->R(0,1)*opt_HLstate(15,0)+state->R(1,1)*opt_HLstate(16,0)+state->R(2,1)*opt_HLstate(17,0),
            //                                     state->R(0,2)*opt_HLstate(15,0)+state->R(1,2)*opt_HLstate(16,0)+state->R(2,2)*opt_HLstate(17,0)}; 
 
            
            // ================================================ //
            // Update contact matrix and initial toe position(s)
            // ================================================ //
            if(con->ind[0]==1){
                con_obj->setDesDomain({0, 1, 1, 0});
                //traj.toeInit.block(0,0,3,1) = kin->toePos.block(0,0,3,1);
                //traj.toeInit.block(0,3,3,1) = kin->toePos.block(0,3,3,1);   
            }else{
                con_obj->setDesDomain({1, 0, 0, 1});
                //traj.toeInit.block(0,1,3,2) = kin->toePos.block(0,1,3,2);
            }
            traj.domLen = domLenSec*ctrlHz;
            traj.toeInit = kin->toePos;

            // ================================================ //
            // Marc Raibert foothold selection (similar)
            // ================================================ //
            
            //std::vector<double> KP = {0.1596,0.1596,0};
            std::vector<double> KP = {0.04,0.02,0.0};
            setStep_Raibert(state,domLenSec,desVel,KP);
        }
   
    }else if(gait==FLY){

        traj.comDes.block(0,0,3,1) = state->q.block(0,0,3,1);
        traj.comDes.block(3,0,3,1) = state->dq.block(0,0,3,1);
        traj.comDes.block(6,0,3,1) = state->q.block(3,0,3,1);
        traj.comDes.block(9,0,3,1) = state->dq.block(3,0,3,1);
    
    }else if(gait==PACE){
        double dt = (1.0/LL_Hz);
        double domLenSec = 0.1;

        if(con->changeDomain==1){
            updateVel(desVel, desOmega, params);
            
            // ================================================ //
            // Update contact matrix and initial toe position(s)
            // ================================================ //
            if(con->des[0]==1){
                con_obj->setDesDomain({0, 1, 0, 1});}
            else{
                con_obj->setDesDomain({1, 0, 1, 0});}
            traj.domLen = domLenSec*ctrlHz;
            traj.toeInit = kin->toePos;

            // ================================================ //
            // Marc Raibert foothold selection 
            // ================================================ //
            std::vector<double> KP = {0.04,0.02,0.0};
            setStep_Raibert(state,domLenSec,desVel,KP);
        }

        Eigen::Matrix<double, 3, 1> desVelWorld = toWorld(desVel,state->R);
    }else if(gait==PRONK){
        double dt = (1.0/LL_Hz);
        double domLenSec = 0.20;
        traj.domLen = domLenSec*ctrlHz;

        static int startTrot = 0;

        if(con->changeDomain==1 ){
            updateVel(desVel, desOmega, params);
        	
            Eigen::Vector3d desVelBody = {state->R(0,0)*opt_HLstate(15,0)+state->R(1,0)*opt_HLstate(16,0)+state->R(2,0)*opt_HLstate(17,0),
                                             state->R(0,1)*opt_HLstate(15,0)+state->R(1,1)*opt_HLstate(16,0)+state->R(2,1)*opt_HLstate(17,0),
                                                state->R(0,2)*opt_HLstate(15,0)+state->R(1,2)*opt_HLstate(16,0)+state->R(2,2)*opt_HLstate(17,0)}; 
 
            
            traj.toeInit = kin->toePos;
            //std::vector<double> KP = {0.1596,0.1596,0};
            std::vector<double> KP = {0.04,0.02,0.0};
            setStep_Raibert(state,domLenSec,desVelBody,KP);
        }
   
    }else if(gait==UPWALK){
        double dt = (1.0/LL_Hz);
        double domLenSec = 0.30;
        traj.domLen = domLenSec*ctrlHz;

        static int startTrot = 0;

        // if(ctrlTick<1 || con->changeDomain==1){

        //     traj.toeInit = kin->toePos;

        // //     // updateVel(desVel, desOmega, params);
        	
        // //     // //Eigen::Vector3d desVelBody = {state->R(0,0)*opt_HLstate(15,0)+state->R(1,0)*opt_HLstate(16,0)+state->R(2,0)*opt_HLstate(17,0),
        // //     // //                                 state->R(0,1)*opt_HLstate(15,0)+state->R(1,1)*opt_HLstate(16,0)+state->R(2,1)*opt_HLstate(17,0),
        // //     // //                                    state->R(0,2)*opt_HLstate(15,0)+state->R(1,2)*opt_HLstate(16,0)+state->R(2,2)*opt_HLstate(17,0)}; 
        // //     // Eigen::Vector3d desVelBody = {state->R(0,0)*desVel(0)+state->R(1,0)*desVel(1)+state->R(2,0)*desVel(2),
        // //     //                                  state->R(0,1)*desVel(0)+state->R(1,1)*desVel(1)+state->R(2,1)*desVel(2),
        // //     //                                     state->R(0,2)*desVel(0)+state->R(1,2)*desVel(1)+state->R(2,2)*desVel(2)}; 
            
            
        // //     // std::vector<double> KP = {0.1596,0,0};//{0.1596,0,0};
        // //     // //std::vector<double> KP = {0.4,0.0,0.0};
        // //     // setStep_Raibert(state,domLenSec,desVel,KP);

        //     setStep_NMPC(NLstep,opt_HLstate(3),state,params,phase);
        //     con_obj->forceDom0();
            
        // }
        if(ctrlTick<1 || phase==0){
            traj.toeInit = kin->toePos;
        }
        if(con->changeDomain==1){
            setStep_NMPC(NLstep,opt_HLstate(3),state,params,phase);
            con_obj->forceDom0(); 
        }

        //if(ctrlTick%10==0){
        //    setStep_NMPC(NLstep,opt_HLstate(3),state,params);
        //}
    }else if(gait==STANDUP){
        //if(ctrlTick<1 || phase==0){
        //    traj.toeInit = kin->toePos;
        //}
        if(con->changeDomain==1){
            //if(phase<0.3){
                traj.toeInit = kin->toePos;
                setFoot(kin);
            //    traj.toeInit.block(0,0,3,2) = kin->toePos.block(0,0,3,2);
            //}else{
            //    traj.toeInit.block(0,2,3,2) = kin->toePos.block(0,2,3,2);
            //}
            con_obj->forceDom0(); 
        }
    }
    
    if(standTrigger && params->neverStopTrot!=1){
        con_obj->setDesDomain({1,1,1,1});
    }

    // Step length saturation
    // traj.stepLen[0] = (traj.stepLen[0]>MAX_SL_F_X) ? MAX_SL_F_X : (traj.stepLen[0]< -MAX_SL_R_X) ? -MAX_SL_R_X : traj.stepLen[0];
    // traj.stepLen[1] = (traj.stepLen[1]>MAX_SL_Y) ? MAX_SL_Y : (traj.stepLen[1]< -MAX_SL_Y) ? -MAX_SL_Y : traj.stepLen[1];
    double dt = (1.0/LL_Hz);
    if(gait!=STAND && gait!=POSE && gait!=TAP && gait!=FLY && gait!=UPWALK && gait!=STANDUP){
        
        
        Eigen::Matrix<double, 3, 1> desVelWorld = toWorld(desVel,state->R);
        Eigen::Matrix<double, 3, 1> desOmegaWorld = toWorld(desOmega,state->R);

        yawOffset = (desOmegaWorld(2)==0) ? yawOffset : state->q(5) + desOmegaWorld(2)*dt;

        Eigen::Matrix<double, 3, 1> pose;
        if(params->extCmd==0){
            pose.setZero();
        }else{
            // pose_ = {delta x, delta y, delta z, roll, pitch, delta yaw}
            Eigen::Matrix<double, 4, 1> quatyaw;
            Eigen::Matrix<double, 4, 1> quatBody;
            Eigen::Matrix<double, 4, 1> quatWorld;
            
            XYZ_to_quat(0,0,state->bodyPos(5),quatyaw);
            XYZ_to_quat(0,state->bodyPos(4),0, quatBody);//pose_(4)
            
            quatMult(quatyaw, quatBody, quatWorld);
            quat_to_XYZ(quatWorld,pose);
        }
        
        traj.comDes.block(0,0,3,1) << state->q.block(0,0,3,1) + desVelWorld*dt;
        traj.comDes(2) = params->standHeight;
        traj.comDes.block(3,0,3,1) << desVelWorld;
        //traj.comDes(4) = 0;
        //traj.comDes(5) = 0;
        traj.comDes(6) = pose(0);
        traj.comDes(7) = pose(1);
        traj.comDes(8) = yawOffset;
        //traj.comDes.block(6,0,3,1) = opt_HLstate.block(6,0,3,1);
        traj.comDes.block(9,0,3,1) = desOmegaWorld;      

    }else if(gait == UPWALK){
        // traj.comDes.block(0,0,3,1) << state->q.block(0,0,3,1);// + opt_HLstate.block(3,0,3,1)*dt;
        // traj.comDes(0) = xnew;
        // traj.comDes(1) = 0;
        traj.comDes(2) = 0.5;//params->standHeight;
        if(ctrlTick>=100000){
            traj.comDes.block(0,0,2,1) << state->q.block(0,0,2,1) + opt_HLstate.block(3,0,2,1)*dt;
            traj.comDes.block(3,0,3,1) = opt_HLstate.block(3,0,3,1);
            traj.comDes.block(6,0,3,1) = opt_HLstate.block(6,0,3,1);
            traj.comDes.block(9,0,3,1) = opt_HLstate.block(9,0,3,1);
        }else{
            traj.comDes.block(0,0,2,1) = state->q.block(0,0,2,1);//state->q(0);// + opt_HLstate.block(3,0,3,1)*dt;
            // traj.comDes(1) = 0; //state->q.block(0,0,2,1);// + opt_HLstate.block(3,0,3,1)*dt;
            traj.comDes.block(3,0,3,1) = Eigen::MatrixXd::Zero(3,1);//opt_HLstate.block(3,0,3,1);
            traj.comDes.block(6,0,3,1) = Eigen::MatrixXd::Zero(3,1);//opt_HLstate.block(6,0,3,1);
            traj.comDes.block(9,0,3,1) = Eigen::MatrixXd::Zero(3,1);//opt_HLstate.block(9,0,3,1);
            traj.comDes(7) = 0.35;
        }
        if(ctrlTick>=34000){ 
            traj.comDes(0)=x_fixed;
            traj.comDes(1)=0.0;
        }
    
    }else if(gait == STANDUP){
        //shiftCoM(state->q(0), state->q(1), state->q(2),0.0);
        traj.comDes.block(0,0,12,1) = Eigen::MatrixXd::Zero(12,1);
        //traj.comDes.block(0,0,3,1) = state->q.block(0,0,3,1);
        traj.comDes(2) = params->standHeight;
    }
}

void MotionPlanner::updateStandVars(const Eigen::Matrix<double,3,1> &com, double yaw, double timeToStand){
    x0 = com(0);
    y0 = com(1);
    z0 = com(2);
    traj.domLen = timeToStand;
    standTime = timeToStand;
    if(timeToStand<100){
        traj.domLen *=ctrlHz;
        standTime *= ctrlHz;
    }
    yawOffset = yaw;
};

void MotionPlanner::updateVel(Eigen::Matrix<double,3,1> &desVel, Eigen::Matrix<double, 3, 1> &desOmega, MP *params){
    if(params->extCmd==0){
        float rate = 0.1;
        int fwd_sgn = (params->fwdSpeed>0) ? 1 : (params->fwdSpeed<0) ? -1 : 0;
        int lat_sgn = (params->latSpeed>0) ? 1 : (params->latSpeed<0) ? -1 : 0;
        desVel(0) += (fwd_sgn*desVel(0)<fwd_sgn*params->fwdSpeed) ? rate*fwd_sgn : 0;
        desVel(1) += (lat_sgn*desVel(1)<lat_sgn*params->latSpeed) ? rate*lat_sgn : 0;

        rate = 0.1;
        int yaw_sgn = (params->yawSpeed>0) ? 1 : (params->yawSpeed<0) ? -1 : 0;
        desOmega(2) += (yaw_sgn*desOmega(2)<yaw_sgn*params->yawSpeed) ? rate*yaw_sgn : 0;
        pose_.setZero();
    }else{
        //desVel(0) = (std::abs(desVel_(0))<0.02) ? 0 : desVel_(0);
        //desVel(1) = (std::abs(desVel_(1))<0.02) ? 0 : desVel_(1);
        desVel(0) += 0.05;
        if(desVel(0)>=desVel_(0)){
            desVel(0)=desVel_(0);
        }
        desVel(1) = (std::abs(desVel_(1))<0.02) ? 0 : desVel_(1);
        desVel(2) = 0;
        desOmega(0) = 0;
        desOmega(1) = 0;
        desOmega(2) = (std::abs(desOmega_(2))<0.05) ? 0 : desOmega_(2);
    }
};

void MotionPlanner::setStep_Raibert(const StateInfo *state, double domLenSec, const Eigen::Matrix<double,3,1> &desVel, std::vector<double> KP){
    // ================================================ //
    // Marc Raibert foothold selection
    // ================================================ //
    setStepLen(0.0,0.0,0.0);
    Eigen::Matrix<double, 3, 1> stepLenTemp;
    stepLenTemp = toBody(state->comFiltered,state->R)-desVel;                   // Vel error
    stepLenTemp(0) *= KP[0]; stepLenTemp(1) *= KP[1]; stepLenTemp(2) *= KP[2];  // mult by KP
    stepLenTemp += (domLenSec*desVel)/2;                                        // Raibert Heuristic
    toWorld(traj.stepLen,stepLenTemp,state->R);                                 // set step length in world frame
}


void MotionPlanner::savesteplen(ContactEst *con_obj){
    const ContactInfo* con = con_obj->getConInfoPointer();

    if(con->des[0]==0){           
        for(int sl=0;sl<3;sl++){
            stepLenRL[sl] = traj.stepLen[sl];
        }        
    }else{           
        for(int sl=0;sl<3;sl++){
            stepLenLR[sl] = traj.stepLen[sl];
        }
    }
}

void MotionPlanner::updatesteplen(ContactEst *con_obj){

    const ContactInfo* con = con_obj->getConInfoPointer();
    if(con->des[0]==0){           
        setStepLen(stepLenRL[0], stepLenRL[1], stepLenRL[2]);      
    }else{           
        setStepLen(stepLenLR[0], stepLenLR[1], stepLenLR[2]);
    }

}

void MotionPlanner::setStep_NMPC(Eigen::Matrix<double,5,1> NLstep, double vdes, const StateInfo *state, MP * params, double phase){
    
    
    // double stepLenTemp = 0;//Eigen::MatrixXd::Zero(3,1);
    // stepLenTemp = sqrt(params->standHeight/9.81)*(state->dq(0)-vdes);//,0,3,1); 

    // traj.FRstepLen = 4*vdes*0.2/2 + 2*stepLenTemp;//(0); 
    // traj.RLstepLen = 4*vdes*0.2/2 + 2*stepLenTemp;//(0);   
    // traj.FLstepLen = 4*vdes*0.2/2 + 2*stepLenTemp;//(0); 
    // traj.RRstepLen = 4*vdes*0.2/2 + 2*stepLenTemp;//(0);
    if(phase<0.3){
        traj.FRstepLen(0) =  NLstep(4);  
        traj.FLstepLen(0) =  NLstep(4); 
    }else{
        traj.RLstepLen(0) =  NLstep(4); 
        traj.RRstepLen(0) =  NLstep(4);  
    }                                
                                    
}


void MotionPlanner::setFoot(const KinematicsInfo *kin){

    double wall_y = 0.2;

    double wall_thresh = 0.18;
    double lat_step = 0.01;

    if(kin->toePos(1,0)>-wall_thresh){
        traj.FRstepLen(1) = -lat_step;
        traj.FRstepLen(2) = 0;
    }else{
        reachedWall = true;
        traj.FRstepLen(0) = lat_step;
        traj.FRstepLen(1) = -wall_y-kin->toePos(1,0);
        traj.FRstepLen(2) = 0*lat_step;
    }

    if(kin->toePos(1,1)<wall_thresh){
        traj.FLstepLen(1) = lat_step;
        traj.FLstepLen(2) = 0;
    }else{
        reachedWall = true;
        traj.FLstepLen(0) = lat_step;
        traj.FLstepLen(1) = wall_y-kin->toePos(1,1);//0;//wall_thresh+lat_step-kin->toePos(1,0);
        traj.FLstepLen(2) = 0*lat_step;
    }

    //traj.FRstepLen =  0.0; 
    //traj.FLstepLen =  0.0;
    // traj.RLstepLen(1) =  0.0; 
    // traj.RRstepLen(1) =  0.0;
}

void MotionPlanner::shiftCoM(ContactEst *con_obj, double phase, size_t shifttime){
    double s = (phase>1) ? 1 : ((phase<0) ? 0 : phase);
    
    double alpha_x[8] = { x0,x0,x0,
                 x0+(xnew-x0)/4,
                 x0+3*(xnew-x0)/4,
                 xnew,xnew,xnew};
    double alpha_y[8] = {y0,y0,y0,
                 y0+(ynew-y0)/4,
                 y0+3*(ynew-y0)/4,
                 ynew,ynew,ynew};
    double alpha_z[8] = {z0,z0,z0,
                 z0+(znew-z0)/4,
                 z0+3*(znew-z0)/4,
                 znew,znew,znew};

    double alpha_p[8] = {p0,p0,p0,
                    p0+(pitchnew-p0)/4,
                    p0+3*(pitchnew-p0)/4,
                    pitchnew,pitchnew,pitchnew};

    double traj_x[3], traj_y[3], traj_z[3], traj_p[3];
    calcBezierAll((int)8, alpha_x, s, traj_x);
    calcBezierAll((int)8, alpha_y, s, traj_y);
    calcBezierAll((int)8, alpha_z, s, traj_z);
    calcBezierAll((int)8, alpha_p, s, traj_p);

    // traj.comDes -> pos, vel, theta, omega
    traj.comDes.block(0,0,3,1) << traj_x[0], traj_y[0], traj_z[0];
    traj.comDes.block(3,0,3,1) << traj_x[1], traj_y[1], traj_z[1];
    traj.comDes.block(6,0,3,1) << 0, traj_p[0], 0;
    traj.comDes.block(9,0,3,1) << 0, 0, 0;

    con_obj->setDesDomain({1,1,1,1});

    traj.domLen = shifttime-30;
    // traj.toeInit.setZero();
    // traj.toeFinal.setZero();
    // traj.toeOffset[2] = Z_TOE_OFFSET;

}

void MotionPlanner::movefoot(size_t movetime, size_t wallsteps){
    
    double pitch_imp = (pitchnew > minpitch) ? pitchnew : minpitch; 
    
    traj.comDes.block(0,0,3,1) << xnew,ynew,znew;
    traj.comDes.block(3,0,3,1) << 0, 0, 0;
    traj.comDes.block(6,0,3,1) << 0, pitch_imp, 0;
    traj.comDes.block(9,0,3,1) << 0, 0, 0;

    //con_obj->setDesDomain({1,1,1,1});
    double frontstep = std::floor(wallsteps/2)*0.05;
    frontstep = (frontstep<0.20) ? frontstep : 0.20;

    traj.domLen = movetime+10;
    traj.FRstepLen = {frontstep,0,upstep};//-0.15};
    traj.FLstepLen = {frontstep,0,upstep};//-0.15};
    traj.RLstepLen = {0,0,0};
    traj.RRstepLen = {0,0,0};
}


void MotionPlanner::movefoot2(size_t movetime, double phase){
    
    double s = (phase>1) ? 1 : ((phase<0) ? 0 : phase);
    
    double alpha_x[8] = { x0,x0,x0,
                 x0+(xnew-x0)/4,
                 x0+3*(xnew-x0)/4,
                 xnew,xnew,xnew};
    double alpha_y[8] = {y0,y0,y0,
                 y0+(ynew-y0)/4,
                 y0+3*(ynew-y0)/4,
                 ynew,ynew,ynew};
    double alpha_z[8] = {z0,z0,z0,
                 z0+(znew-z0)/4,
                 z0+3*(znew-z0)/4,
                 znew,znew,znew};

    double traj_x[3], traj_y[3], traj_z[3];
    calcBezierAll((int)8, alpha_x, s, traj_x);
    calcBezierAll((int)8, alpha_y, s, traj_y);
    calcBezierAll((int)8, alpha_z, s, traj_z);
    
    
    traj.comDes.block(0,0,3,1) << traj_x[0], traj_y[0], traj_z[0];
    traj.comDes.block(3,0,3,1) << traj_x[1], traj_y[1], traj_z[1];
    //traj.comDes.block(0,0,3,1) << xnew,ynew,znew;
    //traj.comDes.block(3,0,3,1) << 0, 0, 0;
    traj.comDes.block(6,0,3,1) << 0, 0, 0;
    traj.comDes.block(9,0,3,1) << 0, 0, 0;

    //con_obj->setDesDomain({1,1,1,1});

    traj.domLen = movetime+20;
    traj.FRstepLen = {0,0,-0.15};
    traj.FLstepLen = {0,0,-0.15};
    traj.RLstepLen = {0,0,0};
    traj.RRstepLen = {0,0,0};
    //std::cout << "Entered here" << std::endl;
}


void MotionPlanner::shiftCoM2(ContactEst *con_obj, double phase, size_t shifttime, bool maxsteps){
    double s = (phase>1) ? 1 : ((phase<0) ? 0 : phase);
    
    double pitch_imp = (pitchnew > minpitch) ? pitchnew : minpitch;
    
    //if(maxsteps){
    //    pitch_imp = minpitch;
    //}else{
        // xnew = rhip_x + 0.183*cos(pitch_imp);
        // znew = rhip_z - 0.183*sin(pitch_imp);
    //    pitch_imp = (pitchnew > minpitch) ? pitchnew : minpitch;
   // }

    if(!maxsteps){
        xnew = rhip_x + 0.183*cos(pitch_imp);
        znew = rhip_z - 0.183*sin(pitch_imp);
    }

    double alpha_x[8] = { x0,x0,x0,
                 x0+(xnew-x0)/4,
                 x0+3*(xnew-x0)/4,
                 xnew,xnew,xnew};
    double alpha_y[8] = {y0,y0,y0,
                 y0+(ynew-y0)/4,
                 y0+3*(ynew-y0)/4,
                 ynew,ynew,ynew};
    double alpha_z[8] = {z0,z0,z0,
                 z0+(znew-z0)/4,
                 z0+3*(znew-z0)/4,
                 znew,znew,znew};

    double alpha_p[8] = {p0,p0,p0,
                 p0+(pitch_imp-p0)/4,
                 p0+3*(pitch_imp-p0)/4,
                 pitch_imp,pitch_imp,pitch_imp};

    double traj_x[3], traj_y[3], traj_z[3], traj_p[3];
    calcBezierAll((int)8, alpha_x, s, traj_x);
    calcBezierAll((int)8, alpha_y, s, traj_y);
    calcBezierAll((int)8, alpha_z, s, traj_z);
    calcBezierAll((int)8, alpha_p, s, traj_p);

    // traj.comDes -> pos, vel, theta, omega
    traj.comDes.block(0,0,3,1) << traj_x[0], traj_y[0], traj_z[0];
    traj.comDes.block(3,0,3,1) << traj_x[1], traj_y[1], traj_z[1];
    traj.comDes.block(6,0,3,1) << 0, traj_p[0], 0;
    traj.comDes.block(9,0,3,1) << 0, 0, 0;

    con_obj->setDesDomain({1,1,1,1});

    traj.domLen = shifttime-30;
    // traj.toeInit.setZero();
    // traj.toeFinal.setZero();
    // traj.toeOffset[2] = Z_TOE_OFFSET;

}

void MotionPlanner::movefoot2(size_t movetime, Eigen::Matrix<double,2,1> xzsteps){
    
    traj.comDes.block(0,0,3,1) << xnew,ynew,znew;
    traj.comDes.block(3,0,3,1) << 0, 0, 0;
    traj.comDes.block(6,0,3,1) << 0, pitchnew, 0;
    traj.comDes.block(9,0,3,1) << 0, 0, 0;

    traj.domLen = movetime+10;
    traj.FRstepLen = {xzsteps(0),0,xzsteps(1)};//-0.15};
    traj.FLstepLen = {xzsteps(0),0,xzsteps(1)};//-0.15};
    traj.RLstepLen = {0,0,0};
    traj.RRstepLen = {0,0,0};
}


void MotionPlanner::shiftCoM3(ContactEst *con_obj, double phase, size_t shifttime, bool maxsteps){
    double s = (phase>1) ? 1 : ((phase<0) ? 0 : phase);
    
    double pitch_imp = pitchnew;

    double alpha_x[8] = { x0,x0,x0,
                 x0+(xnew-x0)/4,
                 x0+3*(xnew-x0)/4,
                 xnew,xnew,xnew};
    double alpha_y[8] = {y0,y0,y0,
                 y0+(ynew-y0)/4,
                 y0+3*(ynew-y0)/4,
                 ynew,ynew,ynew};
    double alpha_z[8] = {z0,z0,z0,
                 z0+(znew-z0)/4,
                 z0+3*(znew-z0)/4,
                 znew,znew,znew};

    double alpha_p[8] = {p0,p0,p0,
                 p0+(pitch_imp-p0)/4,
                 p0+3*(pitch_imp-p0)/4,
                 pitch_imp,pitch_imp,pitch_imp};

    double traj_x[3], traj_y[3], traj_z[3], traj_p[3];
    calcBezierAll((int)8, alpha_x, s, traj_x);
    calcBezierAll((int)8, alpha_y, s, traj_y);
    calcBezierAll((int)8, alpha_z, s, traj_z);
    calcBezierAll((int)8, alpha_p, s, traj_p);

    // traj.comDes -> pos, vel, theta, omega
    traj.comDes.block(0,0,3,1) << traj_x[0], traj_y[0], traj_z[0];
    traj.comDes.block(3,0,3,1) << traj_x[1], traj_y[1], traj_z[1];
    traj.comDes.block(6,0,3,1) << 0, traj_p[0], 0;
    traj.comDes.block(9,0,3,1) << 0, 0, 0;

    con_obj->setDesDomain({1,1,1,1});

    traj.domLen = shifttime;
    // traj.toeInit.setZero();
    // traj.toeFinal.setZero();
    // traj.toeOffset[2] = Z_TOE_OFFSET;

}

void MotionPlanner::movefoot3(size_t movetime){
    
    traj.comDes.block(0,0,3,1) << xnew,ynew,znew;
    traj.comDes.block(3,0,3,1) << 0, 0, 0;
    traj.comDes.block(6,0,3,1) << 0, pitchnew, 0;
    traj.comDes.block(9,0,3,1) << 0, 0, 0;

    traj.domLen = movetime;
    traj.FRstepLen = {0.15,-0.05,0};//-0.15};
    traj.FLstepLen = {0.15,-0.05,0};//-0.15};
    traj.RLstepLen = {0,0,0};
    traj.RRstepLen = {0,0,0};
}


void MotionPlanner::shiftCoMk(ContactEst *con_obj, int wallstep, double phase, size_t shifttime, bool maxsteps){
    double s = (phase>1) ? 1 : ((phase<0) ? 0 : phase);
    
    double pitch_imp = (pitchnew > minpitch) ? pitchnew : minpitch;
    
    //if(maxsteps){
    //    pitch_imp = minpitch;
    //}else{
        // xnew = rhip_x + 0.183*cos(pitch_imp);
        // znew = rhip_z - 0.183*sin(pitch_imp);
    //    pitch_imp = (pitchnew > minpitch) ? pitchnew : minpitch;
   // }

    if(wallstep>0 && !maxsteps){
        xnew = rhip_x + 0.183*cos(pitch_imp);
        znew = rhip_z - 0.183*sin(pitch_imp);
    }

    double alpha_x[8] = { x0,x0,x0,
                 x0+(xnew-x0)/4,
                 x0+3*(xnew-x0)/4,
                 xnew,xnew,xnew};
    double alpha_y[8] = {y0,y0,y0,
                 y0+(ynew-y0)/4,
                 y0+3*(ynew-y0)/4,
                 ynew,ynew,ynew};
    double alpha_z[8] = {z0,z0,z0,
                 z0+(znew-z0)/4,
                 z0+3*(znew-z0)/4,
                 znew,znew,znew};

    double alpha_p[8] = {p0,p0,p0,
                 p0+(pitch_imp-p0)/4,
                 p0+3*(pitch_imp-p0)/4,
                 pitch_imp,pitch_imp,pitch_imp};

    double traj_x[3], traj_y[3], traj_z[3], traj_p[3];
    calcBezierAll((int)8, alpha_x, s, traj_x);
    calcBezierAll((int)8, alpha_y, s, traj_y);
    calcBezierAll((int)8, alpha_z, s, traj_z);
    calcBezierAll((int)8, alpha_p, s, traj_p);

    // traj.comDes -> pos, vel, theta, omega
    traj.comDes.block(0,0,3,1) << traj_x[0], traj_y[0], traj_z[0];
    traj.comDes.block(3,0,3,1) << traj_x[1], traj_y[1], traj_z[1];
    traj.comDes.block(6,0,3,1) << 0, traj_p[0], 0;
    traj.comDes.block(9,0,3,1) << 0, 0, 0;

    con_obj->setDesDomain({1,1,1,1});

    traj.domLen = shifttime-30;
    // traj.toeInit.setZero();
    // traj.toeFinal.setZero();
    // traj.toeOffset[2] = Z_TOE_OFFSET;

}

void MotionPlanner::movefootk(size_t movetime, size_t wallsteps, bool maxsteps){
    
    double pitch_imp = (pitchnew > minpitch) ? pitchnew : minpitch; 
    
    traj.comDes.block(0,0,3,1) << xnew,ynew,znew;
    traj.comDes.block(3,0,3,1) << 0, 0, 0;
    traj.comDes.block(6,0,3,1) << 0, pitch_imp, 0;
    traj.comDes.block(9,0,3,1) << 0, 0, 0;

    double frontstep = std::floor(wallsteps/2)*0.05;
    double upstep2 = -0.2+std::floor(wallsteps/2)*0.04;

    if(maxsteps){
        frontstep = 0.1;
        upstep2 = 0;
    }else{    
        frontstep = (frontstep<0.20) ? frontstep : 0.2;
        upstep2 = (frontstep<-0.01) ? upstep2 : -0.01;
    }
    traj.domLen = movetime+10;
    traj.FRstepLen = {frontstep,0,upstep2};//-0.15};
    traj.FLstepLen = {frontstep,0,upstep2};//-0.15};
    traj.RLstepLen = {0,0,0};
    traj.RRstepLen = {0,0,0};
}