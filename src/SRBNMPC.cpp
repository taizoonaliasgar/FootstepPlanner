#include "SRBNMPC.hpp"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <casadi/casadi.hpp>
#include <casadi/core/generic_type.hpp>
#include <casadi/core/dae_builder.hpp>
#include <casadi/core/optistack.hpp>



namespace fs = std::filesystem;

SRBNMPC::SRBNMPC(int argc, char *argv[], int numRobots, int id) : Parameters(argc,argv){
    // filename = ""; // empty string will produce no output file
    //if(!filename.empty()){
    //    fid = std::fstream(filename, std::ios::out);
    //}
    //success = (fid.is_open()) ? 1 : 0;
    for (size_t i=0; i< FILE_CNT_MPC; i++){
        file[i].open(FILE_NAMES_MIT[i].c_str(), std::ios_base::out);
    }



    J << 0.016840,0.000084,0.000598,0.000084,0.056579,0.000025,0.000598,0.000025,0.064714;
    Jinv << 59.402578,-0.087845,-0.548594,-0.087845,17.674526,-0.006053,-0.548594,-0.006053,15.457771;
    
    Rstand << 0,0,-1,  0,1,0,  1,0,0;
    Jstand = Rstand*J*Rstand.transpose();
    Jinvstand = Rstand*Jinv*Rstand.transpose();

    std::copy(Jstand.data(), Jstand.data() + Jstand.size(), Jstandcasadi.ptr());
    std::copy(Jinvstand.data(), Jinvstand.data() + Jinvstand.size(), Jinvcasadi.ptr());

    // contact_sequence(1,casadi::Slice(5,20)) = casadi::MX::zeros(1,15);
    // contact_sequence(2,casadi::Slice(8,20)) = casadi::MX::zeros(1,12);
    // contact_sequence(0,casadi::Slice(25,40)) = casadi::MX::zeros(1,15);
    // contact_sequence(3,casadi::Slice(28,40)) = casadi::MX::zeros(1,12);

    contact_sequence_dm(1,casadi::Slice(5,20)) = casadi::DM::zeros(1,15);
    contact_sequence_dm(2,casadi::Slice(8,20)) = casadi::DM::zeros(1,12);
    contact_sequence_dm(0,casadi::Slice(25,40)) = casadi::DM::zeros(1,15);
    contact_sequence_dm(3,casadi::Slice(28,40)) = casadi::DM::zeros(1,12);

    for (size_t i = 0; i < HORIZ; i++)
    {
        forcefitx(i) = i*10;
    }

    forcefitx = forcefitx/(HORIZ-1);
    forcefitx = forcefitx/10.0;
    //std::cout << forcefitx << std::endl;

    
}


void SRBNMPC::impactDetection(size_t tick, Eigen::Matrix<double,12,1> &q0, size_t gait){}

void SRBNMPC::generator(){
    
    casadi::SX p = casadi::SX::sym("p",(NFS+NFI)*HORIZ+NFS+4*(HORIZ+1),1);
    casadi::SX x = casadi::SX::sym("x", (NFS+NFI)*HORIZ+NFS);

    // Objective
    casadi::SX f = UpdateCostN(x,p);
    std::string basePath = "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/build/tmp/";
    writeMatrixToFile(f, basePath + "obj.txt");
    
    // Constraints
    casadi::SX g = UpdateConstraintsN(x,p);
    writeMatrixToFile(g, basePath + "const.txt");

    casadi::Dict opts;
    opts["ipopt.max_iter"] = 200;  // Replace Max_mpciter with its actual value
    opts["ipopt.print_level"] = 0;  // Can be changed to 0 or 3 based on verbosity required
    opts["print_time"] = 0;  // Disable printing solver time
    opts["ipopt.acceptable_tol"] = 1e-2;  // Tolerance for stopping criterion
    opts["ipopt.acceptable_obj_change_tol"] = 1e-2;  // Objective change tolerance for stopping

    casadi::Function solver = casadi::nlpsol("solver", "ipopt", {{"x", x}, {"f", f}, {"g", g}, {"p", p}}, opts);

    // file name
    std::string file_name = "upright_nlp_70";
    // code predix
    std::string prefix_code = fs::current_path().string() + "/";

    // Generate C code for the NLP functions
    solver.generate_dependencies(file_name + ".c");

    std::string prefix_lib = fs::current_path().string() + "/";

    // compile c code to a shared library
    std::string compile_command = "gcc -fPIC -shared -O3 " + 
        prefix_code + file_name + ".c -o " +
        prefix_lib + file_name + ".so";

    std::cout << compile_command << std::endl;
    int compile_flag = std::system(compile_command.c_str());
    casadi_assert(compile_flag==0, "Compilation failed");
    std::cout << "Compilation successed!" << std::endl;

}

casadi::DM SRBNMPC::motionPlannerN(Eigen::Matrix<double,16,1> q0, size_t controlTick){
    

    casadi::DM x_des = casadi::DM::zeros(NFS*(HORIZ+1)+NFI*(HORIZ)+4*(HORIZ+1)); 
    casadi::DM q0_dm = casadi::DM::zeros(16,1);
    for (int i=0; i<16; i++){
        q0_dm(i) = q0(i);
    }
    
    x_des(casadi::Slice(0,NFS)) = q0_dm(casadi::Slice(0,NFS));
    
    casadi::DM conp1=0;

    if(controlTick%40 == 0){
        if(localvelocity < desVel(0)){
            localvelocity = localvelocity + 0.05;
        }else{
            localvelocity = desVel(0);
        }
    }
    //Raibheur = 0.5*Tstance*abs(q0(3));

    if(controlTick%20 == 0){
        Raibheur = 0.5*Tstance*(q0(3));//localvelocity;
        x_dom_init = q0_dm(0)+3*Raibheur;//Tstance*localvelocity;
    }

    //Raibheur < 0 ? Raibheur = 0 : Raibheur = Raibheur; 
    
    for(size_t i= 0; i< HORIZ; i++){
        
        conp1 = remainder(controlTick+i,40);
        x_des((i+1)*NFS) = x_des(0) + (i+1)*localvelocity*MPC_dt, 
        x_des((i+1)*NFS+1) = 0;//x_des(1) + (i+1)*desVel(1)*MPC_dt; 
        x_des((i+1)*NFS+2) = stand_height;                                 
        x_des((i+1)*NFS+3) = localvelocity, 
        x_des((i+1)*NFS+4) = desVel(1);

        //x_des((i+1)*NFS+12) = x_des(12)+Raibheur;
        //x_des((i+1)*NFS+13) = x_des(13)+Raibheur;
        //x_des((i+1)*NFS+14) = x_des(14)+Raibheur;
        //x_des((i+1)*NFS+15) = x_des(15)+Raibheur;
        
        x_des((i+1)*NFS+12) = contact_sequence_dm(0,conp1)*x_des(12)+(1-contact_sequence_dm(0,conp1))*x_dom_init;//x_des(12)+Raibheur
        x_des((i+1)*NFS+13) = contact_sequence_dm(1,conp1)*x_des(13)+(1-contact_sequence_dm(1,conp1))*x_dom_init;//x_des(13)+Raibheur
        x_des((i+1)*NFS+14) = contact_sequence_dm(2,conp1)*x_des(14)+(1-contact_sequence_dm(2,conp1))*x_dom_init;//x_des(14)+Raibheur
        x_des((i+1)*NFS+15) = contact_sequence_dm(3,conp1)*x_des(15)+(1-contact_sequence_dm(3,conp1))*x_dom_init;//x_des(15)+Raibheur

        x_des((HORIZ+1)*NFS+i*NFI+12) = 2*(1-contact_sequence_dm(0,conp1))*Raibheur;//0.4
        x_des((HORIZ+1)*NFS+i*NFI+13) = 2*(1-contact_sequence_dm(1,conp1))*Raibheur;//0.4
        x_des((HORIZ+1)*NFS+i*NFI+14) = 2*(1-contact_sequence_dm(2,conp1))*Raibheur;//0.4
        x_des((HORIZ+1)*NFS+i*NFI+15) = 2*(1-contact_sequence_dm(3,conp1))*Raibheur;//0.4
        //x_des(casadi::Slice((HORIZ+1)*NFS+i*NFI+12,(HORIZ+1)*NFS+i*NFI+16)) = 2*(casadi::DM::ones(4,1)-contact_sequence_dm(casadi::Slice(0,4),conp1))*Raibheur;

        x_des(NFS*(HORIZ+1)+NFI*(HORIZ)+4*i) = contact_sequence_dm(0,conp1);
        x_des(NFS*(HORIZ+1)+NFI*(HORIZ)+4*i+1) = contact_sequence_dm(1,conp1);
        x_des(NFS*(HORIZ+1)+NFI*(HORIZ)+4*i+2) = contact_sequence_dm(2,conp1);
        x_des(NFS*(HORIZ+1)+NFI*(HORIZ)+4*i+3) = contact_sequence_dm(3,conp1);
    }

    
    conp1=remainder(controlTick+HORIZ,40);
    x_des(NFS*(HORIZ+1)+NFI*(HORIZ)+4*HORIZ) = contact_sequence_dm(0,conp1);
    x_des(NFS*(HORIZ+1)+NFI*(HORIZ)+4*HORIZ+1) = contact_sequence_dm(1,conp1);
    x_des(NFS*(HORIZ+1)+NFI*(HORIZ)+4*HORIZ+2) = contact_sequence_dm(2,conp1);
    x_des(NFS*(HORIZ+1)+NFI*(HORIZ)+4*HORIZ+3) = contact_sequence_dm(3,conp1);
    
    
    return x_des;
    
}


casadi::SX SRBNMPC::UpdateCostN(casadi::SX x, casadi::SX x_des){
    
    casadi::SX f = 0;
    
    Eigen::Matrix3d R_force = Eigen::Matrix3d::Zero();
    Eigen::Matrix<double, NFS, NFS> Q = Eigen::Matrix<double, NFS, NFS>::Zero();
    Eigen::Matrix<double, NFI-4, NFI-4> R = Eigen::Matrix<double, NFI-4, NFI-4>::Zero();
    Eigen::Matrix<double, 4, 4> Rv = Eigen::Matrix<double, 4, 4>::Zero();
    
    
    Q.setZero();
    Q.block(0,0,3,3).diagonal() <<  1e4,1e4,8e5;//mpc_params.qpx, mpc_params.qpy, mpc_params.qpz;
    Q.block(3,3,3,3).diagonal() <<  1e4,1e4,1e4;//mpc_params.qvx, mpc_params.qvy, mpc_params.qvz;
    Q.block(6,6,3,3).diagonal() <<  8e4,8e5,3e4;//mpc_params.qrr, mpc_params.qrp, mpc_params.qry;
    //Yaw:3e3
    Q.block(9,9,3,3).diagonal() <<  1,1,1;//mpc_params.qwr, mpc_params.qwp, mpc_params.qwy;
    Q.block(12,12,4,4).diagonal() <<  1e-2,1e-2,1e-2,1e-2;//1e3//1e5,1e5,1e5,1e5;
    //repdiag(Q,Q_rep,HORIZ+1);
    
    R_force.setZero();
    R_force.diagonal() << 0.01,0.01,0.01;//mpc_params.rx, mpc_params.ry, mpc_params.rz;
    for(int i=0;i<4;i++){
        if(i<2){
            R.block(3*i,3*i,3,3) = 0.1*R_force;
        }else{
            R.block(3*i,3*i,3,3) = R_force;
        }
        //R.block(12+i,12+i,1,1) << 0.01;
    }
    
    //R = 1e2*R;
    Rv.setZero();
    Rv.diagonal() << 0.01,0.01,0.01,0.01;//1e2,1e2,1e2,1e2;//0.01,0.01,0.01,0.01;
    //Rv = 1e5*Rv;//1e5

    casadi::DM Qc = casadi::DM::zeros(Q.rows(),Q.cols());
    std::copy(Q.data(), Q.data() + Q.size(), Qc.ptr());

    casadi::DM Rc = casadi::DM::zeros(R.rows(),R.cols());
    std::copy(R.data(), R.data() + R.size(), Rc.ptr());
    
    casadi::DM Rvc = casadi::DM::zeros(Rv.rows(),Rv.cols());
    std::copy(Rv.data(), Rv.data() + Rv.size(), Rvc.ptr());

    casadi::SX st,con;
    casadi::SX p_Raibheur, f_des;// = 0.5*Tstance*localvelocity*casadi::SX::ones(4,1);

    for(int i=0;i<HORIZ;i++){
        st = x(casadi::Slice(i*NFS,i*NFS+NFS));
        con = x(casadi::Slice(NFS*(HORIZ+1) + i*NFI,NFS*(HORIZ+1) + i*NFI+NFI));
        p_Raibheur = x_des(casadi::Slice(NFS*(HORIZ+1)+NFI*i+12,NFS*(HORIZ+1)+NFI*i+NFI));
        f_des = x_des(casadi::Slice(NFS*(HORIZ+1) + i*NFI,NFS*(HORIZ+1) + i*NFI+12));

        f= f+ mtimes(mtimes(transpose(st-x_des(casadi::Slice(i*NFS,i*NFS+NFS))),casadi::SX(Qc)),st-x_des(casadi::Slice(i*NFS,i*NFS+NFS)))
                + mtimes(mtimes(transpose(con(casadi::Slice(0,12))-f_des),casadi::SX(Rc)),con(casadi::Slice(0,12))-f_des) 
                    +mtimes(mtimes(transpose(con(casadi::Slice(12,16))-p_Raibheur),casadi::SX(Rvc)),con(casadi::Slice(12,16))-p_Raibheur);
    }
    

    return f;
}

casadi::SX SRBNMPC:: UpdateConstraintsN(casadi::SX x, casadi::SX p){
    // ====================================================================== //
    // ======================== Equality Constraints ======================== //
    // ====================================================================== //
    casadi::SX g_constraints = casadi::SX::zeros(NFS*(HORIZ+1)+NFI*HORIZ,1);
    
    casadi::SX st,con,f_value,st_next,st_next_euler,conp1;

    g_constraints(casadi::Slice(0,NFS)) = x(casadi::Slice(0,NFS))-p(casadi::Slice(0,NFS));
    

    for(int i=0; i<HORIZ; ++i){
        
        conp1 = p(casadi::Slice(NFS*(HORIZ+1)+HORIZ*NFI+4*(i+1),NFS*(HORIZ+1)+HORIZ*NFI+4*(i+2)));
        
        st = x(casadi::Slice(i*NFS,i*NFS+NFS));
         
        con = x(casadi::Slice(NFS*(HORIZ+1) + i*NFI,NFS*(HORIZ+1) + i*NFI+NFI));
        
        st_next = x(casadi::Slice((i+1)*NFS,(i+1)*NFS+NFS));
        
        f_value = NonlinearDynamics(st,con,conp1); 

        st_next_euler = st + (f_value);
        
        g_constraints(casadi::Slice((i+1)*NFS,(i+1)*NFS+NFS)) = st_next-st_next_euler;
        
        g_constraints(casadi::Slice((HORIZ+1)*NFS+i*NFI,(HORIZ+1)*NFS+i*NFI+NFI)) = inequalitycons(con);
        
    }

    return g_constraints;

}

casadi::SX SRBNMPC::NonlinearDynamics(casadi::SX st,casadi::SX con, casadi::SX conp1){

    casadi::SX rhs = casadi::SX::zeros(NFS,1);

    casadi::SX pose = st(casadi::Slice(0,3));

    casadi::SX velocity = st(casadi::Slice(3,6));

    casadi::SX phi = st(6);
    casadi::SX theta = st(7);
    casadi::SX psi = st(8);
    casadi::SX w = st(casadi::Slice(9,12));
    
    casadi::SX f1 = con(casadi::Slice(0,3));
    casadi::SX f2 = con(casadi::Slice(3,6));
    casadi::SX f3 = con(casadi::Slice(6,9));
    casadi::SX f4 = con(casadi::Slice(9,12));
    
    casadi::SX A = casadi::SX::zeros(3,3);
    A(0,0) = 1;
    A(0,1) = sin(phi)*tan(theta);
    A(0,2) = cos(phi)*tan(theta);
    A(1,0) = 0;
    A(1,1) = cos(phi);
    A(1,2) = -sin(phi);
    A(2,0) = 0;
    A(2,1) = sin(phi)/cos(theta);
    A(2,2) = cos(phi)/cos(theta);

    casadi::SX tau = GetTorque(st,con);

    rhs(casadi::Slice(0,3)) = st(casadi::Slice(3,6)) * MPC_dt;

    rhs(casadi::Slice(3,6)) = (f1+f2+f3+f4)*MPC_dt/MASS - casadi::SX(gravityN) * MPC_dt;

    rhs(casadi::Slice(6,9)) = mtimes(A,st(casadi::Slice(9,12)))*MPC_dt;

    casadi::SX Jw = mtimes(casadi::SX(Jstandcasadi),st(casadi::Slice(9,12)));
    casadi::SX what = skewsym(st(casadi::Slice(9,12)));
    casadi::SX wJw = mtimes(what,Jw);
    rhs(casadi::Slice(9,12)) = mtimes(casadi::SX(Jinvcasadi),tau-wJw) * MPC_dt;
    
    rhs(12) = con(12)*conp1(0);//contact(0);
    rhs(13) = con(13)*conp1(1);//contact(1);
    rhs(14) = con(14)*conp1(2);//contact(2);
    rhs(15) = con(15)*conp1(3);//contact(3);

    return rhs;
}

casadi::SX SRBNMPC::GetTorque(casadi::SX st,casadi::SX con){

    casadi::SX r1 = casadi::SX::zeros(3,1);
    casadi::SX r2 = casadi::SX::zeros(3,1);
    casadi::SX r3 = casadi::SX::zeros(3,1);
    casadi::SX r4 = casadi::SX::zeros(3,1);
    r1(0) = st(12)-st(0);
    r1(1) = -0.25-st(1);
    r1(2) = 0.533-st(2);

    r2(0) = st(13)-st(0);
    r2(1) = 0.25-st(1);
    r2(2) = 0.533-st(2);
    
    r3(0) = st(14)-st(0);
    r3(1) = -0.1321-st(1);
    r3(2) = 0-st(2);
    
    r4(0) = st(15)-st(0);
    r4(1) = 0.1321-st(1);
    r4(2) = 0-st(2);
    
    casadi::SX Mr1 = skewsym(r1);
    casadi::SX Mr2 = skewsym(r2);
    casadi::SX Mr3 = skewsym(r3);  
    casadi::SX Mr4 = skewsym(r4);
    
    return mtimes(Mr1,con(casadi::Slice(0,3)))+mtimes(Mr2,con(casadi::Slice(3,6)))+mtimes(Mr3,con(casadi::Slice(6,9)))+mtimes(Mr4,con(casadi::Slice(9,12)));

}

casadi::SX SRBNMPC::inequalitycons(casadi::SX con){

       casadi::SX g_inequality = casadi::SX::zeros(16,1);

       for(int leg=0;leg<4;leg++){

            if(leg<1){
                g_inequality(4*leg) = con(0) - mu*con(1);
                g_inequality(4*leg+1) = con(2) - mu*con(1);
                g_inequality(4*leg+2) = con(0) + mu*con(1);
                g_inequality(4*leg+3) = con(2) + mu*con(1);
            }else if(leg<2){
                g_inequality(4*leg) = con(3) + mu*con(4);
                g_inequality(4*leg+1) = con(5) + mu*con(4);
                g_inequality(4*leg+2) = con(3) - mu*con(4);
                g_inequality(4*leg+3) = con(5) - mu*con(4);
            }else{
                g_inequality(4*leg) = con(3*leg) - mu*con(3*leg+2);
                g_inequality(4*leg+1) = con(3*leg+1) - mu*con(3*leg+2);
                g_inequality(4*leg+2) = con(3*leg) + mu*con(3*leg+2);
                g_inequality(4*leg+3) = con(3*leg+1) + mu*con(3*leg+2);
            }

       }
       return g_inequality;

}

casadi::DM SRBNMPC::lowerboundx(casadi::DM p){
    
    casadi::DM lbx = -casadi::DM::inf(NFS*(HORIZ+1)+NFI*HORIZ,1);
    casadi::DM contact_index = casadi::DM::zeros(4,1); 

    for(int k=0 ; k<HORIZ ; k++){

        contact_index = p(casadi::Slice(NFS*(HORIZ+1)+HORIZ*NFI+4*(k),NFS*(HORIZ+1)+HORIZ*NFI+4*(k+1)));
        for(int leg=0 ; leg<4 ; leg++){
            if(leg<1){
                lbx(NFS*(HORIZ+1)+NFI*k+1) = 0;

            }else if(leg<2){
                lbx(NFS*(HORIZ+1)+NFI*(k)+4) = -contact_index(leg)*fzmaxf;
                
            }else{
                lbx(NFS*(HORIZ+1)+NFI*(k)+3*(leg)+2) = 0;
        
            }
        }
        //lbx(casadi::Slice(NFS*(HORIZ+1)+NFI*k+12,NFS*(HORIZ+1)+NFI*k+16)) = (casadi::DM::ones(4,1)-contact_index)*Raibheur;//casadi::DM::zeros(4,1);//-casadi::DM::ones(4,1)*Raibheur*10;
        lbx(NFS*(HORIZ+1)+NFI*k+12) = 0*(1-contact_index(0))*Raibheur;//0.2
        lbx(NFS*(HORIZ+1)+NFI*k+13) = 0*(1-contact_index(1))*Raibheur;
        lbx(NFS*(HORIZ+1)+NFI*k+14) = 0*(1-contact_index(2))*Raibheur;
        lbx(NFS*(HORIZ+1)+NFI*k+15) = 0*(1-contact_index(3))*Raibheur;
    }
    
    return lbx;

}

casadi::DM SRBNMPC::upperboundx(casadi::DM p){
    
    casadi::DM ubx = casadi::DM::inf(NFS*(HORIZ+1)+NFI*HORIZ,1);
    casadi::DM contact_index = casadi::DM::zeros(4,1);
    

    for(int k=0 ; k<HORIZ ; k++){
        
        contact_index = p(casadi::Slice(NFS*(HORIZ+1)+HORIZ*NFI+4*(k),NFS*(HORIZ+1)+HORIZ*NFI+4*(k+1)));
        
        for(int leg=0 ; leg<4 ; leg++){
            if(leg<1){
                ubx(NFS*(HORIZ+1)+NFI*(k)+1) = contact_index(leg)*fzmaxf;
            }else if(leg<2){
                ubx(NFS*(HORIZ+1)+NFI*(k)+4) = 0;
            }else{
                ubx(NFS*(HORIZ+1)+NFI*(k)+3*(leg)+2) = contact_index(leg)*fzmaxr;
            }
        }
        //ubx(casadi::Slice(NFS*(HORIZ+1)+NFI*(k)+12,NFS*(HORIZ+1)+NFI*(k)+16)) = (casadi::DM::ones(4,1)-contact_index)*RaibMult*Raibheur;
        // ubx(NFS*(HORIZ+1)+NFI*k+12) = RaibMult*(1-contact_index(0))*Raibheur;//0.6
        // ubx(NFS*(HORIZ+1)+NFI*k+13) = 0.6*(1-contact_index(1))*Raibheur;
        // ubx(NFS*(HORIZ+1)+NFI*k+14) = 0.6*(1-contact_index(2))*Raibheur;
        // ubx(NFS*(HORIZ+1)+NFI*k+15) = 0.6*(1-contact_index(3))*Raibheur;

        ubx(casadi::Slice(NFS*(HORIZ+1)+NFI*k+12,NFS*(HORIZ+1)+NFI*k+16)) = RaibMult*(casadi::DM::ones(4,1)-contact_index)*Raibheur;
    }
    
    return ubx;
}

casadi::DM SRBNMPC::lowerboundg(){
    
    casadi::DM lbg = casadi::DM::zeros(NFS*(HORIZ+1)+NFI*HORIZ,1);
    
    for(int k=0 ; k<HORIZ ; k++){
        
        for(int leg=0 ; leg<4 ; leg++){     
            lbg(casadi::Slice(NFS*(HORIZ+1)+NFI*k+4*leg,NFS*(HORIZ+1)+NFI*k+4*leg+2)) = -casadi::DM::inf();    
        }
    }
    
    return lbg;
}

casadi::DM SRBNMPC::upperboundg(){
    
    casadi::DM ubg = casadi::DM::zeros(NFS*(HORIZ+1)+NFI*HORIZ,1);

    for(int k=0 ; k<HORIZ ; k++){
        
        for(int leg=0 ; leg<4 ; leg++){     
            ubg(casadi::Slice(NFS*(HORIZ+1)+NFI*k+4*leg+2,NFS*(HORIZ+1)+NFI*k+4*leg+4)) = casadi::DM::inf();    
        }
    }
    return ubg;
}

casadi::SX SRBNMPC::skewsym(casadi::SX v3){

    casadi::SX M3 = casadi::SX::zeros(3, 3);
    M3(0,0) = 0;
    M3(0,1) = -v3(2);
    M3(0,2) = v3(1);
    
    M3(1,0) = v3(2);
    M3(1,1) = 0;
    M3(1,2) = -v3(0);
    
    M3(2,0) = -v3(1);
    M3(2,1) = v3(0);
    M3(2,2) = 0;

    return M3;

}


Eigen::Matrix<double,12,HORIZ> SRBNMPC::arrangeOptforce(){
    
    Eigen::Matrix<double,12,HORIZ> optforce = Eigen::Matrix<double,12,HORIZ>::Zero();
    std::vector<double> v = previous_sol(casadi::Slice(static_cast<casadi_int>(NFS*(HORIZ+1)), static_cast<casadi_int>(NFS*(HORIZ+1)+12))).get_elements();

    optforce.block(0,0,12,1) = Eigen::Map<Eigen::Matrix<double,12,1>>(v.data());
    for (size_t i = 1; i < HORIZ; i++)
    {
        v = previous_sol(casadi::Slice(static_cast<casadi_int>(NFS*(HORIZ+1)+NFI*i), static_cast<casadi_int>(NFS*(HORIZ+1)+NFI*i+12))).get_elements();
        optforce.block(0,i,12,1) = Eigen::Map<Eigen::Matrix<double,12,1>>(v.data());
    }

    //Eigen::Matrix<double,12,1> optforce_vector = Eigen::Map<Eigen::Matrix<double,12,1>>(v.data());
    
    return optforce;

}

Eigen::Matrix<double,12,1> SRBNMPC::getFootPos(){
    
    Eigen::Matrix<double,12,1> foot_pos;

    std::vector<double> v = previous_sol(casadi::Slice(12,16)).get_elements();
    Eigen::Matrix<double,4,1> foot_x = Eigen::Map<Eigen::Matrix<double,4,1>>(v.data());

    for(int leg=0;leg<4;leg++){
        foot_pos(3*leg) = foot_x(leg);
        if(leg<2){
            foot_pos(3*leg+1) = pow(-1,leg+1)*0.25;
            foot_pos(3*leg+2) = 0.533;
        }else{
            foot_pos(3*leg+1) = pow(-1,leg+1)*0.1321;
            foot_pos(3*leg+2) = 0;
        }
    }
    
    return foot_pos;

}

casadi::DM SRBNMPC::getprevioussol(Eigen::Matrix<double,16,1> q0, size_t controlTick){
        
    casadi::DM x0 = casadi::DM::zeros(NFS*(HORIZ+1)+NFI*HORIZ);
    casadi::DM q0_dm = casadi::DM::zeros(16,1);
    for (int i=0; i<12; i++){
        q0_dm(i) = q0(i);
    }

    x0(casadi::Slice(0,12)) = q0_dm(casadi::Slice(0,12));
    if(controlTick<1){
        x0(12) = 0.2;//0.02
        x0(13) = 0.2;
        x0(14) = -0.1;
        x0(15) = -0.1;
        casadi::DM sum_conseq;// = casadi::DM::sum1(contact_sequence_dm); 
        
        for(int i = 0; i<HORIZ; i++){
            x0(casadi::Slice(NFS*(i+1),NFS*(i+2))) = x0(casadi::Slice(0,NFS));
            sum_conseq = contact_sequence_dm(2,i)+contact_sequence_dm(3,i)+contact_sequence_dm(0,i)+contact_sequence_dm(1,i);
            
            for (size_t leg = 0; leg < 4; leg++)
            {
                if(leg>1){
                    x0(NFS*(HORIZ+1)+NFI*(i)+3*(leg+1)-1) = contact_sequence_dm(leg,i)*MASS*gravityN(2)/sum_conseq;
                }else{
                    x0(NFS*(HORIZ+1)+NFI*(i)+3*(leg+1)-1) = contact_sequence_dm(leg,i)*MASS*gravityN(2)/sum_conseq;//*mu*10;
                    x0(NFS*(HORIZ+1)+NFI*(i)+3*(leg+1)-2) = 2*pow(-1,leg)*x0(NFS*(HORIZ+1)+NFI*(i)+3*(leg+1)-1);
                    
                } 
            }
            
        }
    }else{
        x0(casadi::Slice(12,16)) = previous_sol(casadi::Slice(28,32));
        x0(casadi::Slice(NFS,NFS*(HORIZ))) = previous_sol(casadi::Slice(NFS*2,NFS*(HORIZ+1)));
        x0(casadi::Slice(NFS*HORIZ,NFS*(HORIZ+1))) = previous_sol(casadi::Slice(NFS*HORIZ,NFS*(HORIZ+1)));

        x0(casadi::Slice(NFS*(HORIZ+1),NFS*(HORIZ+1)+NFI*(HORIZ-1))) = previous_sol(casadi::Slice(NFS*(HORIZ+1)+NFI,NFS*(HORIZ+1)+NFI*(HORIZ)));
        x0(casadi::Slice(NFS*(HORIZ+1)+NFI*(HORIZ-1),NFS*(HORIZ+1)+NFI*(HORIZ))) = previous_sol(casadi::Slice(NFS*(HORIZ+1)+NFI*(HORIZ-1),NFS*(HORIZ+1)+NFI*(HORIZ)));

    }
    
    return x0;
}


void SRBNMPC::writeMatrixToFile(const casadi::SX& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Print matrix dimensions
        file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << std::endl;
        //Print matrix content preserving original shape
        for (int i = 0; i < matrix.size1(); ++i) {
            for (int j = 0; j < matrix.size2(); ++j) {
                file << matrix(i, j);
                if (j < matrix.size2() - 1) file << ", ";
            }
            file << std::endl;
        }
        file.close();

    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

void SRBNMPC::writeMatrixToFileDM(const casadi::DM& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Print matrix dimensions
        file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << std::endl;
        //Print matrix content preserving original shape
        for (int i = 0; i < matrix.size1(); ++i) {
            for (int j = 0; j < matrix.size2(); ++j) {
                file << matrix(i, j);
                if (j < matrix.size2() - 1) file << ", ";
            }
            file << std::endl;
        }
        file.close();

    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

Eigen::Matrix<double,24,1> SRBNMPC::getNMPCsol(){
    
    Eigen::Matrix<double,24,1> nmpc_sol;
    std::vector<double> v = previous_sol.get_elements();

    std::vector<double> optstate(v.begin() + NFS, v.begin() + 2*NFS);
    std::vector<double> optforce(v.begin() + NFS*(HORIZ+1), v.begin() + NFS*(HORIZ+1) + NFI);
    
    nmpc_sol.block(0,0,12,1) = Eigen::Map<Eigen::Matrix<double,12,1>>(optstate.data());
    nmpc_sol.block(0,12,12,1) = Eigen::Map<Eigen::Matrix<double,12,1>>(optforce.data());
    
    
    return nmpc_sol;

}

void SRBNMPC::mpcdataLog(Eigen::Matrix<double,16,1> q0, Eigen::Matrix<double,12,1> force, size_t tick){//casadi::DM X_prev,(Eigen::Matrix<double,12,1> &q0, size_t tick){
      
    file[0] << tick << "\t" << q0(0) << "\t" << q0(1) << "\t" 
                << q0(2) << "\t" << q0(3) << "\t" 
                << q0(4) << "\t" << q0(5) << "\t" 
                << q0(6) << "\t" << q0(7) << "\t" 
                << q0(8) << "\t" << q0(9) << "\t" 
                << q0(10) << "\t" << q0(11) << "\t"
                << q0(12) << "\t" << q0(13) << "\t" 
                << q0(14) << "\t" << q0(15) <<std::endl;

      
    file[1] << force(0) << "\t" << force(1) << "\t" 
                << force(2) << "\t" << force(3) << "\t" 
                << force(4) << "\t" << force(5) << "\t" 
                << force(6) << "\t" << force(7) << "\t" 
                << force(8) << "\t" << force(9) << "\t" 
                << force(10) << "\t" << force(11) << "\t"
                << previous_sol(NFS*(HORIZ+1)+12) << "\t" << previous_sol(NFS*(HORIZ+1)+13) << "\t" 
                << previous_sol(NFS*(HORIZ+1)+14) << "\t" << previous_sol(NFS*(HORIZ+1)+15) << std::endl;

    file[2] <<  previousp(16) << "\t" << previousp(17) << "\t" 
                << previousp(18) << "\t" << previousp(19) << "\t" 
                << previousp(20) << "\t" << previousp(21) << "\t" 
                << previousp(22) << "\t" << previousp(23) << "\t" 
                << previousp(24) << "\t" << previousp(25) << "\t" 
                << previousp(26) << "\t" << previousp(27) << "\t"
                << previousp(28) << "\t" << previousp(29) << "\t" 
                << previousp(30) << "\t" << previousp(31) << std::endl;
     
    file[3] << previousp(NFS*(HORIZ+1)+NFI*(HORIZ)) << "\t" << previousp(NFS*(HORIZ+1)+NFI*(HORIZ)+1) << "\t" 
                << previousp(NFS*(HORIZ+1)+NFI*(HORIZ)+2) << "\t" << previousp(NFS*(HORIZ+1)+NFI*(HORIZ)+3) << std::endl;
    
    file[4] << previousp(12) << "\t" << previousp(13) << "\t" << previousp(14) << "\t" << previousp(15) << std::endl;

    // file[5] <<  X_prev(NFS*(HORIZ+1))     << "\t" << X_prev(NFS*(HORIZ+1) + 1) << "\t" << 
    //             X_prev(NFS*(HORIZ+1) + 2) << "\t" << X_prev(NFS*(HORIZ+1) + 3) << "\t" << 
    //             X_prev(NFS*(HORIZ+1) + 4) << "\t" << X_prev(NFS*(HORIZ+1) + 5) << "\t" << 
    //             X_prev(NFS*(HORIZ+1) + 6) << "\t" << X_prev(NFS*(HORIZ+1) + 7) << "\t" << 
    //             X_prev(NFS*(HORIZ+1) + 8) << "\t" << X_prev(NFS*(HORIZ+1) + 9) << "\t" << 
    //             X_prev(NFS*(HORIZ+1) + 10) << "\t" << X_prev(NFS*(HORIZ+1) + 11) << "\t" << std::endl;

    
    /*file[6] << for(size_t i = 0; i < count; i++)
    {
        /* code */
    //}*/
    for(size_t i = 0; i < HORIZ; i++){
        for (size_t j = 0; j < 4; j++)
        {
            file[5] << previous_sol(NFS*(HORIZ+1)+NFI*i+12+j) << "\t";    /* code */
        }
        
        file[5] << std::endl;
    }

    for(size_t i = 0; i < HORIZ; i++){
        for (size_t j = 0; j < 4; j++)
        {
            file[6] << previous_sol(NFS*(i+1)+12+j) << "\t";    /* code */
        }
        
        file[6] << std::endl;
    }
    
                         
}



void SRBNMPC::getOptForceCoeff(int order){

    forceCoeff = Eigen::MatrixXd::Zero(12,order+1);
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(HORIZ, order + 1);
    Eigen::MatrixXd Y = Eigen::VectorXd::Zero(HORIZ,1);
    Eigen::MatrixXd optforce = arrangeOptforce();  
    //std::cout << "Optimum force: " << "\n" << optforce << std::endl;
    
    for (size_t f_ind = 0; f_ind < 12; ++f_ind){
        Y = optforce.block(f_ind,0,1,HORIZ).transpose();
    
        for (size_t i = 0; i < HORIZ; ++i) {
            for (int j = 0; j <= order; ++j) {
                X(i, j) = std::pow(forcefitx(i), j);
            }
        }
        
        Eigen::VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * Y);
        forceCoeff.block(f_ind,0,1,order+1) = coeffs.transpose();
        
    }
    //writeMatrixToFileDM(forceCoeff, "forceCoeff.txt");
    //std::cout << "Force Coefficients: " << "\n" << forceCoeff << std::endl;
}


Eigen::Matrix<double,12,1> SRBNMPC::getOptforce(int phase,int order){

    Eigen::Matrix<double,12,1> result = Eigen::MatrixXd::Zero(12,1);
    //std::cout << "Phase: " << phase << std::endl;
    double phased = double(phase)/(HORIZ-1);
    phased = phased/10.0;
    //std::cout << "Phased: " << phased << std::endl;
    for (size_t f_ind = 0; f_ind < 12; f_ind++)
    {
        for (int i = 0; i < order+1; ++i) {
            result(f_ind) = result(f_ind) + forceCoeff(f_ind,i) * std::pow(phased, i);
        }
    }
    //std::cout << result.transpose() << std::endl;
    
    
    return result;
}

void SRBNMPC::getForce(){
    
    //Eigen::Matrix<double,12,1> optforce = Eigen::Matrix<double,12,1>::Zero();
    std::vector<double> v = previous_sol(casadi::Slice(static_cast<casadi_int>(NFS*(HORIZ+1)), static_cast<casadi_int>(NFS*(HORIZ+1)+12))).get_elements();
    OPTforce = Eigen::Map<Eigen::Matrix<double,12,1>>(v.data());

    //return OPTforce;

}