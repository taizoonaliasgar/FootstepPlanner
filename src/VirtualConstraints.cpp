//
// Authror: Randy Fawcett on 12/2021.
//
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "VirtualConstraints.hpp"

using VirtCon = VirtualConstraints;

VirtCon::VirtualConstraints(){
    VC.fDes.setZero();
}

void VirtCon::updateVirtualConstraints(const StateInfo *state, const KinInf *kin, const TrajInfo *traj, const ConInf *con, size_t gait, double phaseVar, MP *params, const LLInfo *ll){
    size_t outDim = 6+3*(4-con->cnt);
    size_t conDim = 3*con->cnt;

    h0.setZero(outDim,1);
    dh0.setZero(outDim,1);
    VC.H0.setZero(outDim,TOTAL_DOF);
    VC.dH0.setZero(outDim,1);
    VC.hd.setZero(outDim,1);
    VC.dhd.setZero(outDim,1);
    VC.ddhd.setZero(outDim,1);
    VC.y.setZero(outDim,1);
    VC.dy.setZero(outDim,1);
    VC.y_ST.setZero(conDim,1);
    VC.dy_ST.setZero(conDim,1);
    VC.hd_ST.setZero(conDim,1);
    VC.dhd_ST.setZero(conDim,1);
    VC.H0.block(0,0,6,6) = Eigen::MatrixXd::Identity(6,6);
    
    h0 = VC.H0*state->q;

    if (gait==STAND){
        VC.hd.block(0,0,3,1)   << traj->comDes.block(0,0,3,1);
        VC.dhd.block(0,0,3,1)  << traj->comDes.block(3,0,3,1);;
        VC.ddhd.block(0,0,3,1) << 0, 0, 0;
        VC.hd.block(3,0,3,1)   << traj->comDes.block(6,0,3,1);
        VC.dhd.block(3,0,3,1)  << 0, 0, 0;
        VC.ddhd.block(3,0,3,1) << 0, 0, 0;
    }
    else {
        VC.hd.block(0,0,3,1) = traj->comDes.block(0,0,3,1);
        VC.hd.block(3,0,3,1) = traj->comDes.block(6,0,3,1);
        VC.dhd.block(0,0,3,1) = traj->comDes.block(3,0,3,1);
        VC.dhd.block(3,0,3,1) = traj->comDes.block(9,0,3,1);
        VC.ddhd.block(0,0,3,1) << 0,0,0;
        VC.ddhd.block(3,0,3,1) << 0,0,0;

        Eigen::VectorXd hipAcc(3);
        Eigen::VectorXd hipVel(3);

        size_t cnts = 0;
        size_t cntc = 0;
        double to = traj->toeOffset[2];
        double ds = (1.0*ctrlHz)/200;//traj->domLen;
//        double phase = (phaseVar>1.0) ? 1.0 : phaseVar;
		double phase = phaseVar;
        //double dt = traj->domLen/(1.0*ctrlHz);
        double dt = 200/(1.0*ctrlHz);
        for(size_t i=0; i<4; i++){
            if(con->ind[i]==0){
                h0.block(6+cnts,0,3,1) = kin->toePos.block(0,i,3,1);
                VC.H0.block(6+cnts,0,3,TOTAL_DOF) = kin->Jtoe.block(3*i,0,3,TOTAL_DOF);
                VC.dH0.block(6+cnts,0,3,1) = kin->dJtoe.block(3*i,0,3,1);
                

                // Swing leg to follow time varying bezier
                hipVel = kin->Jhip.block(3*i,0,3,18)*state->dq;
                hipAcc = ( kin->Jhip.block(3*i,0,3,18)*ll->ddq + kin->dJhip.block(3*i,0,3,1) );
                // hipAcc.setZero();
                
                double tune=0;
                if (gait==PACE){
                    tune = (2*(i%2==0)-1)*0.04; // eventually this needs to be rotated by R (body to world)
                }
                
                int n = 4;
                double tx[3] = {0};
                double ax[4]{traj->toeInit(0,i), traj->toeInit(0,i), 
                            kin->hipPos(0,i)+traj->stepLen[0], kin->hipPos(0,i)+traj->stepLen[0]};
                double dax[4] {0, 0, hipVel(0), hipVel(0)};
                double ddax[4] {0, 0, hipAcc(0), hipAcc(0)};
                calcVaryingBezierAll(n,dt,ax,dax,ddax,phase,tx);
                
                double ty[3] = {0};
                double ay[4] = {traj->toeInit(1,i), traj->toeInit(1,i), 
                                kin->hipPos(1,i)+traj->stepLen[1]+tune, kin->hipPos(1,i)+traj->stepLen[1]+tune};
                double day[4] = {0, 0, hipVel(1), hipVel(1)};
                double dday[4] = {0, 0, hipAcc(1), hipAcc(1)};
                calcVaryingBezierAll(n,dt,ay,day,dday,phase,ty);

                n = 8;
                double tz[3] = {0};
                double az[8]= {traj->toeInit(2,i), traj->toeInit(2,i), params->swingHeight, params->swingHeight, 
                               params->swingHeight, to+0.005, to+0.005, to};
                calcBezierAll(n, az, phaseVar, tz);

                // Save foot traj
                VC.hd.block(6+cnts,0,3,1)   << tx[0], ty[0], tz[0];
                VC.dhd.block(6+cnts,0,3,1)  << tx[1], ty[1], tz[1]*ds;
                VC.ddhd.block(6+cnts,0,3,1) << tx[2], ty[2], tz[2]*ds*ds; // z scaling necessary!!

                cnts+=3;
            }
        }
    }
    dh0 = VC.H0*state->dq;

    VC.y = h0-VC.hd;
    VC.dy = dh0-VC.dhd;
}

void VirtCon::updateVirtualConstraintsflight(const StateInfo *state, const KinInf *kin, const TrajInfo *traj, const ConInf *con, size_t gait, 
                                                double phaseVarRL, double phaseVarLR, double stepLenRL[3], double stepLenLR[3], MP *params, const LLInfo *ll){

    size_t outDim = 12;//6+3*(4-con->cnt);
    size_t conDim = 0;//3*con->cnt;

    // Only initialize size if necessary
    h0.setZero(outDim,1);
    dh0.setZero(outDim,1);
    VC.H0.setZero(outDim,TOTAL_DOF);
    VC.dH0.setZero(outDim,1);
    VC.hd.setZero(outDim,1);
    VC.dhd.setZero(outDim,1);
    VC.ddhd.setZero(outDim,1);
    VC.y.setZero(outDim,1);
    VC.dy.setZero(outDim,1);
        
    Eigen::VectorXd hipAcc(3);
    Eigen::VectorXd hipVel(3);

    size_t cnts = 0;
    size_t cntc = 0;
    double to = traj->toeOffset[2];
    double ds = (1.0*ctrlHz)/200;
    double sstepLen[3];
    double phase;

    double dt = 200/(1.0*ctrlHz);
    for(size_t i=0; i<4; i++){

        h0.block(3*i,0,3,1) = kin->toePos.block(0,i,3,1);
        VC.H0.block(3*i,0,3,TOTAL_DOF) = kin->Jtoe.block(3*i,0,3,TOTAL_DOF);
        VC.dH0.block(3*i,0,3,1) = kin->dJtoe.block(3*i,0,3,1);
        
        // Swing leg to follow time varying bezier
        hipVel = kin->Jhip.block(3*i,0,3,18)*state->dq;
        hipAcc = ( kin->Jhip.block(3*i,0,3,18)*ll->ddq + kin->dJhip.block(3*i,0,3,1) );
        // hipAcc.setZero();
                
        int n = 4;
        if(i==0 || i==1){
            phase = phaseVarRL;
            sstepLen[0] = stepLenRL[0];
            sstepLen[1] = stepLenRL[1];
            sstepLen[2] = stepLenRL[2];
        
        }else{
            phase = phaseVarLR;
            sstepLen[0] = stepLenLR[0];
            sstepLen[1] = stepLenLR[1];
            sstepLen[2] = stepLenLR[2];
        }

        double tx[3] = {0};
        double ax[4]{traj->toeInit(0,i), traj->toeInit(0,i), 
                    kin->hipPos(0,i)+sstepLen[0], kin->hipPos(0,i)+sstepLen[0]};
        double dax[4] {0, 0, hipVel(0), hipVel(0)};
        double ddax[4] {0, 0, hipAcc(0), hipAcc(0)};
        calcVaryingBezierAll(n,dt,ax,dax,ddax,phase,tx);
                
        double ty[3] = {0};
        double ay[4] = {traj->toeInit(1,i), traj->toeInit(1,i), 
                        kin->hipPos(1,i)+sstepLen[1], kin->hipPos(1,i)+sstepLen[1]};
        double day[4] = {0, 0, hipVel(1), hipVel(1)};
        double dday[4] = {0, 0, hipAcc(1), hipAcc(1)};
        calcVaryingBezierAll(n,dt,ay,day,dday,phase,ty);

        n = 8;
        double tz[3] = {0};
        double az[8]= {traj->toeInit(2,i), traj->toeInit(2,i), params->swingHeight, params->swingHeight, 
                        params->swingHeight, to+0.005, to+0.005, to};
        calcBezierAll(n, az, phase, tz);

        // Save foot traj
        VC.hd.block(3*i,0,3,1)   << tx[0], ty[0], tz[0];
        VC.dhd.block(3*i,0,3,1)  << tx[1], ty[1], tz[1]*ds;
        VC.ddhd.block(3*i,0,3,1) << tx[2], ty[2], tz[2]*ds*ds; // z scaling necessary!!
    }

    dh0 = VC.H0*state->dq;

    VC.y = h0-VC.hd;
    VC.dy = dh0-VC.dhd;

}

void VirtCon::updateVirtualConstraintswalk(const StateInfo *state, const KinInf *kin, const TrajInfo *traj, const ConInf *con, size_t gait, double frontphase, double rearphase, MP *params, const LLInfo *ll){
    size_t outDim = 6+3*(4-con->cnt);
    size_t conDim = 3*con->cnt;

    h0.setZero(outDim,1);
    dh0.setZero(outDim,1);
    VC.H0.setZero(outDim,TOTAL_DOF);
    VC.dH0.setZero(outDim,1);
    VC.hd.setZero(outDim,1);
    VC.dhd.setZero(outDim,1);
    VC.ddhd.setZero(outDim,1);
    VC.y.setZero(outDim,1);
    VC.dy.setZero(outDim,1);
    VC.y_ST.setZero(conDim,1);
    VC.dy_ST.setZero(conDim,1);
    VC.hd_ST.setZero(conDim,1);
    VC.dhd_ST.setZero(conDim,1);
    VC.H0.block(0,0,6,6) = Eigen::MatrixXd::Identity(6,6);
    
    h0 = VC.H0*state->q;

    double step[3] = {0.0,0.0,0.0};//{traj->stepLen[0],traj->stepLen[1],traj->stepLen[2]};
    //double temp = step[1];

    //std::cout << "Swingheight" << params->swingHeight << std::endl;

    if (gait==STAND){
        VC.hd.block(0,0,3,1)   << traj->comDes.block(0,0,3,1);
        VC.dhd.block(0,0,3,1)  << traj->comDes.block(3,0,3,1);;
        VC.ddhd.block(0,0,3,1) << 0, 0, 0;
        VC.hd.block(3,0,3,1)   << traj->comDes.block(6,0,3,1);
        VC.dhd.block(3,0,3,1)  << 0, 0, 0;
        VC.ddhd.block(3,0,3,1) << 0, 0, 0;
    }
    else {
        VC.hd.block(0,0,3,1) = traj->comDes.block(0,0,3,1);
        VC.hd.block(3,0,3,1) = traj->comDes.block(6,0,3,1);
        VC.dhd.block(0,0,3,1) = traj->comDes.block(3,0,3,1);
        VC.dhd.block(3,0,3,1) = traj->comDes.block(9,0,3,1);
        VC.ddhd.block(0,0,3,1) << 0,0,0;
        VC.ddhd.block(3,0,3,1) << 0,0,0;

        Eigen::VectorXd hipAcc(3);
        Eigen::VectorXd hipVel(3);

        size_t cnts = 0;
        size_t cntc = 0;
        double to = traj->toeOffset[2];
        //double ds = (1.0*ctrlHz)/(two_contact*10-5);
        //double dt = (two_contact*10-5)/(1.0*ctrlHz);
        double ds = (1.0*ctrlHz)/(145);
        double dt = (145)/(1.0*ctrlHz);
        //double phase = (phaseVar>1.0) ? 1.0 : phaseVar;
		double phase = frontphase;
        //double dt = traj->domLen/(1.0*ctrlHz);
        
        for(size_t i=0; i<4; i++){
            if(con->ind[i]==0){
                
                switch (i)
                {
                case 0:
                    step[0] = traj->FRstepLen;
                    break;
                case 1:
                    step[0] = traj->FLstepLen;
                    break;
                case 2:
                    step[0] = traj->RRstepLen;
                    break;
                case 3:
                    step[0] = traj->RLstepLen;
                    break;
                }

                h0.block(6+cnts,0,3,1) = kin->toePos.block(0,i,3,1);
                VC.H0.block(6+cnts,0,3,TOTAL_DOF) = kin->Jtoe.block(3*i,0,3,TOTAL_DOF);
                VC.dH0.block(6+cnts,0,3,1) = kin->dJtoe.block(3*i,0,3,1);
                

                // Swing leg to follow time varying bezier
                hipVel = kin->Jhip.block(3*i,0,3,18)*state->dq;
                hipAcc = ( kin->Jhip.block(3*i,0,3,18)*ll->ddq + kin->dJhip.block(3*i,0,3,1) );
                // hipAcc.setZero();
                
                double tune=0;
                if (gait==PACE){
                    tune = (2*(i%2==0)-1)*0.04; // eventually this needs to be rotated by R (body to world)
                }

                double tx[3] = {0};
                double ty[3] = {0};
                double tz[3] = {0};
                int n = 4;
                
                if(i>1){
                    phase = rearphase;
                    double dt = (115)/(1.0*ctrlHz);
                    double ds = (1.0*ctrlHz)/(115);
                    //double dt = (two_contact*10-5)/(1.0*ctrlHz);
                    //double ds = (1.0*ctrlHz)/(two_contact*10-5);
                    //step[0] = step[0]/2;
                    //step[1] = step[1]/2;
                    //step[2] = step[2]/2;

                    //double ty[3] = {0};
                    double ay[4] = {traj->toeInit(1,i), traj->toeInit(1,i), 
                                    kin->hipPos(1,i)+step[1]+tune+pow(-1,i+1)*0.0838, kin->hipPos(1,i)+step[1]+tune+pow(-1,i+1)*0.0838};
                    double day[4] = {0, 0, hipVel(1), hipVel(1)};
                    double dday[4] = {0, 0, hipAcc(1), hipAcc(1)};
                    calcVaryingBezierAll(n,dt,ay,day,dday,phase,ty);

                    
                    //double ax[4]{traj->toeInit(0,i), traj->toeInit(0,i), 
                    //        kin->hipPos(0,i)+step[0], kin->hipPos(0,i)+step[0]};
                    double ax[4]{traj->toeInit(0,i), traj->toeInit(0,i)+step[0], 
                            traj->toeInit(0,i)+step[0], traj->toeInit(0,i)+step[0]};
                    double dax[4] {0, 0, hipVel(0), hipVel(0)};
                    double ddax[4] {0, 0, hipAcc(0), hipAcc(0)};
                    calcVaryingBezierAll(n,dt,ax,dax,ddax,phase,tx);

                    n = 8;
                    //double tz[3] = {0};
                    //double az[8]= {traj->toeInit(2,i), traj->toeInit(2,i), params->swingHeight, 2*params->swingHeight, 
                    //            params->swingHeight, to+0.005, to, to};
                    double az[8]= {traj->toeInit(2,i), traj->toeInit(2,i), params->swingHeight, 2*params->swingHeight, 
                                2*params->swingHeight, params->swingHeight, 0.019,0.019};//traj->toeInit(2,i),traj->toeInit(2,i)};//0.02, 0.02};//to+0.01, to};
                    calcBezierAll(n, az, phase, tz);

                }else{
                
                    n=4;
                    //double tz[3] = {0};
                    double az[4] = {traj->toeInit(2,i), traj->toeInit(2,i), 
                                    kin->hipPos(2,i)+step[2]+tune, kin->hipPos(2,i)+step[2]+tune};
                    double daz[4] = {0, 0, hipVel(2), hipVel(2)};
                    double ddaz[4] = {0, 0, hipAcc(2), hipAcc(2)};
                    calcVaryingBezierAll(n,dt,az,daz,ddaz,phase,tz);

                    //double ax[4]{traj->toeInit(0,i), traj->toeInit(0,i)+0.1, 
                    //            kin->hipPos(0,i)+2*step[0], kin->hipPos(0,i)+2*step[0]};
                    double ax[4]{traj->toeInit(0,i), traj->toeInit(0,i),//+step[0]/2, //+0.1
                                traj->toeInit(0,i)+step[0], traj->toeInit(0,i)+step[0]};
                    double dax[4] {0, 0, hipVel(0), hipVel(0)};
                    double ddax[4] {0, 0, hipAcc(0), hipAcc(0)};
                    calcVaryingBezierAll(n,dt,ax,dax,ddax,phase,tx);

                    n = 8;
                    //double ty[3] = {0};
                    //double ay[8]= {traj->toeInit(1,i), traj->toeInit(1,i), pow(-1,i+1)*wall_y+pow(-1,i)*params->swingHeight, pow(-1,i+1)*wall_y+pow(-1,i)*params->swingHeight, 
                    //            pow(-1,i+1)*wall_y+pow(-1,i)*params->swingHeight, pow(-1,i+1)*wall_y+pow(-1,i)*(to+0.005), pow(-1,i+1)*wall_y, pow(-1,i+1)*wall_y+pow(-1,i)*to};
                    double ay[8]= {traj->toeInit(1,i), traj->toeInit(1,i), pow(-1,i+1)*wall_y+pow(-1,i)*params->swingHeight, pow(-1,i+1)*wall_y+pow(-1,i)*params->swingHeight, 
                                pow(-1,i+1)*wall_y+pow(-1,i)*params->swingHeight, pow(-1,i+1)*wall_y+pow(-1,i)*params->swingHeight, pow(-1,i+1)*wall_y, pow(-1,i+1)*wall_y};
                    
                    calcBezierAll(n, ay, phase, ty);
                }
                
                // Save foot traj
                VC.hd.block(6+cnts,0,3,1)   << tx[0], ty[0], tz[0];
                if(i>1){
                    ds = (1.0*ctrlHz)/115;//(two_contact*10-5);
                    VC.dhd.block(6+cnts,0,3,1)  << tx[1], ty[1], tz[1]*ds;
                    VC.ddhd.block(6+cnts,0,3,1) << tx[2], ty[2], tz[2]*ds*ds; // z scaling necessary!!
                }else{
                    ds = (1.0*ctrlHz)/145;//(two_contact*10-5);
                    VC.dhd.block(6+cnts,0,3,1)  << tx[1], ty[1]*ds, tz[1];
                    VC.ddhd.block(6+cnts,0,3,1) << tx[2], ty[2]*ds*ds, tz[2]; // z scaling necessary!!
                }
                cnts+=3;
            }
        }
    }
    dh0 = VC.H0*state->dq;

    VC.y = h0-VC.hd;
    VC.dy = dh0-VC.dhd;
}
