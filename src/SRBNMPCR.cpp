#include "SRBNMPCR.hpp"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <casadi/casadi.hpp>
#include <casadi/core/generic_type.hpp>



namespace fs = std::filesystem;

SRBNMPCR::SRBNMPCR(int argc, char *argv[], int numRobots, int id) : Parameters(argc,argv){
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
    contact_sequence_dm = casadi::DM::ones(4,40);
    contact_sequence_dm(1,casadi::Slice(5,20)) = casadi::DM::zeros(1,15);
    contact_sequence_dm(2,casadi::Slice(8,20)) = casadi::DM::zeros(1,12);
    contact_sequence_dm(0,casadi::Slice(25,40)) = casadi::DM::zeros(1,15);
    contact_sequence_dm(3,casadi::Slice(28,40)) = casadi::DM::zeros(1,12);

    localvelocity = 0;
    for (size_t i = 0; i < HORIZ; i++)
    {
        forcefitx(i) = i*10;
    }

    forcefitx = forcefitx/(HORIZ-1);
    forcefitx = forcefitx/10.0;
    
}


void SRBNMPCR::impactDetection(size_t tick, Eigen::Matrix<double,12,1> &q0, size_t gait){}

void SRBNMPCR::mpcdataLog(Eigen::Matrix<double,17,1> q0, Eigen::Matrix<double,12,1> force, Eigen::Matrix<double, 12, 1> p_foot,size_t tick){//(Eigen::Matrix<double,12,1> &q0, size_t tick){
     
    file[0] << tick << "\t" << q0(0) << "\t" << q0(1) << "\t" 
                << q0(2) << "\t" << q0(3) << "\t" 
                << q0(4) << "\t" << q0(5) << "\t" 
                << q0(6) << "\t" << q0(7) << "\t" 
                << q0(8) << "\t" << q0(9) << "\t" 
                << q0(10) << "\t" << q0(11) <<std::endl;
                

      
    file[1] << force(0) << "\t" << force(1) << "\t" 
                << force(2) << "\t" << force(3) << "\t" 
                << force(4) << "\t" << force(5) << "\t" 
                << force(6) << "\t" << force(7) << "\t" 
                << force(8) << "\t" << force(9) << "\t" 
                << force(10) << "\t" << force(11) << std::endl;

    file[2] <<  p_prev(12) << "\t" << p_prev(13) << "\t" 
                << p_prev(14) << "\t" << p_prev(15) << "\t" 
                << p_prev(16) << "\t" << p_prev(17) << "\t" 
                << p_prev(18) << "\t" << p_prev(19) << "\t" 
                << p_prev(20) << "\t" << p_prev(21) << "\t" 
                << p_prev(22) << "\t" << p_prev(23) << std::endl;
     
    file[3] << p_prev(NFSR*(HORIZ+1)+NFIR*HORIZ) << "\t" << p_prev(NFSR*(HORIZ+1)+NFIR*HORIZ+1)  << "\t" 
                << p_prev(NFSR*(HORIZ+1)+NFIR*HORIZ+2) << "\t" << p_prev(NFSR*(HORIZ+1)+NFIR*HORIZ+3)  << std::endl; 

    file[4] << p_prev(NFSR*(HORIZ+1)+(NFIR+4)*HORIZ) << "\t" << p_prev(NFSR*(HORIZ+1)+(NFIR+4)*HORIZ+1) << "\t" 
                << p_prev(NFSR*(HORIZ+1)+(NFIR+4)*HORIZ+2) << "\t" << p_prev(NFSR*(HORIZ+1)+(NFIR+4)*HORIZ+3) << std::endl;

    file[5] << p_foot(0) << "\t" << p_foot(3) << "\t" 
                << p_foot(6) << "\t" << p_foot(9) << std::endl;

    // file[5] <<  X_prev(NFSR*(HORIZ+1))     << "\t" << X_prev(NFSR*(HORIZ+1) + 1) << "\t" << 
    //             X_prev(NFSR*(HORIZ+1) + 2) << "\t" << X_prev(NFSR*(HORIZ+1) + 3) << "\t" << 
    //             X_prev(NFSR*(HORIZ+1) + 4) << "\t" << X_prev(NFSR*(HORIZ+1) + 5) << "\t" << 
    //             X_prev(NFSR*(HORIZ+1) + 6) << "\t" << X_prev(NFSR*(HORIZ+1) + 7) << "\t" << 
    //             X_prev(NFSR*(HORIZ+1) + 8) << "\t" << X_prev(NFSR*(HORIZ+1) + 9) << "\t" << 
    //             X_prev(NFSR*(HORIZ+1) + 10) << "\t" << X_prev(NFSR*(HORIZ+1) + 11) << "\t" << std::endl;
               
}



void SRBNMPCR::generator(){
    

    //casadi::SX p = casadi::SX::sym("p",(NFSR+NFIR)*HORIZ+NFSR+1,1);
    casadi::SX p = casadi::SX::sym("p",(NFSR+NFIR)*HORIZ+NFSR+8*(HORIZ),1);
    //casadi::SX p = casadi::SX::sym("p",(NFSR+NFIR)*HORIZ+NFSR+4*(HORIZ),1);
    //casadi::SX p = casadi::SX::sym("p",(NFSR+NFIR)*HORIZ+NFSR,1);
    casadi::SX x = casadi::SX::sym("x", (NFSR+NFIR)*HORIZ+NFSR);
    

    // Objective
    casadi::SX f = UpdateCostN(x,p);
    //std::cout << "f"<<std::endl;
    std::string basePath = "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/build/tmp/";
    writeMatrixToFile(f, basePath + "obj.txt");
    
    // Constraints
    casadi::SX g = UpdateConstraintsN(x,p);
    //std::cout << "g"<<std::endl;
    writeMatrixToFile(g, basePath + "const.txt");

    // Create an NLP solver instance 
    casadi::Dict opts = {{"ipopt.print_level", 0}, {"print_time", 0},{"ipopt.max_iter", 200},{"ipopt.acceptable_tol", 1e-2},{"ipopt.acceptable_obj_change_tol", 1e-2}};
    //std::cout << "opts"<<std::endl;
    casadi::Function solver = casadi::nlpsol("solver", "ipopt", {{"x", x}, {"f", f}, {"g", g}, {"p", p}}, opts);
    //casadi::Function solver = casadi::nlpsol("solver", "ipopt", nlppp);
    //std::cout << "solver"<<std::endl;
    int num_variables = solver.n_out();
    //std::cout << num_variables << std::endl;

    // file name
    std::string file_name = "upright_nlp_9r";
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

casadi::DM SRBNMPCR::motionPlannerN(Eigen::Matrix<double,17,1> q0, size_t controlTick){
    
    //casadi::DM x_des = casadi::DM::zeros(NFSR*(HORIZ+1)+NFIR*(HORIZ));
    //casadi::DM x_des = casadi::DM::zeros(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*(HORIZ));
    casadi::DM x_des = casadi::DM::zeros(NFSR*(HORIZ+1)+NFIR*(HORIZ)+8*(HORIZ));
    casadi::DM q0_dm(q0.rows(), q0.cols());
    for (int i=0; i<17; i++){
        q0_dm(i) = q0(i);
    }
    //x_des(casadi::Slice(0,2)) = q0_dm(casadi::Slice(0,2));
    x_des(casadi::Slice(0,NFSR)) = q0_dm(casadi::Slice(0,NFSR));
    
    casadi::DM conp1, conp1_prev, contact_index;
    if(controlTick%40 == 0){
        if(localvelocity < desVel(0)){
            localvelocity = localvelocity + 0.05;
        }else{
            localvelocity = desVel(0);
        }
        Raibheur = 0.5*Tstance*(abs(q0(3))) + sqrt(stand_height/9.81)*(abs(q0(3))-localvelocity);
    }
    //Raibheur = 0.5*Tstance*localvelocity;//desVel(0);
    
    //contact_index = remainder(controlTick,40);
    for(size_t i= 0; i< HORIZ; i++){
        
        conp1 =(controlTick+i)%40;//remainder(contact_index+i,40);
        if(controlTick+i<1){
            conp1_prev = 0;//remainder(controlTick+i-1,40);
        }else{
            conp1_prev = (controlTick+i-1)%40;
        }

        x_des((i+1)*NFSR) = x_des(0) + (i+1)*localvelocity*MPC_dt; 
        x_des((i+1)*NFSR+1) = 0;//x_des(1) + (i+1)*desVel(1)*MPC_dt; 
        x_des((i+1)*NFSR+2) = stand_height;                                 
        x_des((i+1)*NFSR+3) = localvelocity, 
        x_des((i+1)*NFSR+4) = desVel(1);
        // x_des((i+1)*NFSR+12) = x_des(12)+Raibheur;
        // x_des((i+1)*NFSR+13) = x_des(13)+Raibheur;
        // x_des((i+1)*NFSR+14) = x_des(14)+Raibheur;
        // x_des((i+1)*NFSR+15) = x_des(15)+Raibheur;

        // x_des((HORIZ+1)*NFSR+i*NFIR+12) = Raibheur;
        // x_des((HORIZ+1)*NFSR+i*NFIR+13) = Raibheur;
        // x_des((HORIZ+1)*NFSR+i*NFIR+14) = Raibheur;
        // x_des((HORIZ+1)*NFSR+i*NFIR+15) = Raibheur;

        
        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*i) = contact_sequence_dm(0,conp1);
        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*i+1) = contact_sequence_dm(1,conp1);
        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*i+2) = contact_sequence_dm(2,conp1);
        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*i+3) = contact_sequence_dm(3,conp1);

        
        //x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i) = 0.2;//contact_sequence_dm(0,conp1);
        //x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+1) = 0.2;//contact_sequence_dm(1,conp1);
        //x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+2) = -0.1;//contact_sequence_dm(2,conp1);
        //x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+3) = -0.1;//contact_sequence_dm(3,conp1);
        
        //if(controlTick>0){
            for (int leg = 0;leg<4;leg++){
                if(contact_sequence_dm(leg,conp1).scalar() > contact_sequence_dm(leg,conp1_prev).scalar()){
                    if(i<1){
                        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+leg) = x_des(0)+Raibheur;//+0.2;//contact_sequence_dm(0,conp1);
                    }else{
                        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+leg) = x_des(NFSR*(i-1))+Raibheur;//-0.1;
                    }    
                }else{
                    if(i<1){
                        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+leg) = q0_dm(12+leg);
                    }else{
                        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+leg) = x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*(i-1)+leg);//x_des(NFSR*i);
                    }
                }
            }
        // }else{
        //     x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i) = 0.2;//contact_sequence_dm(0,conp1);
        //     x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+1) = 0.2;//contact_sequence_dm(1,conp1);
        //     x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+2) = -0.1;//contact_sequence_dm(2,conp1);
        //     x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+3) = -0.1;//contact_sequence_dm(3,conp1);

        // }
        
    }
    
    return x_des;
    
}


casadi::SX SRBNMPCR::UpdateCostN(casadi::SX x, casadi::SX x_des){
    
    casadi::SX f = 0;
    
    Eigen::Matrix3d R_force = Eigen::Matrix3d::Zero();
    Eigen::Matrix<double, NFSR, NFSR> Q = Eigen::Matrix<double, NFSR, NFSR>::Zero();
    Eigen::Matrix<double, NFIR, NFIR> R = Eigen::Matrix<double, NFIR, NFIR>::Zero();
    
    
    Q.setZero();
    Q.block(0,0,3,3).diagonal() <<  1e4,1e4,8e5;//mpc_params.qpx, mpc_params.qpy, mpc_params.qpz;
    Q.block(3,3,3,3).diagonal() <<  1e4,1e4,1e4;//mpc_params.qvx, mpc_params.qvy, mpc_params.qvz;
    Q.block(6,6,3,3).diagonal() <<  8e4,8e5,3e3;//mpc_params.qrr, mpc_params.qrp, mpc_params.qry;
    Q.block(9,9,3,3).diagonal() <<  1,1,1;//mpc_params.qwr, mpc_params.qwp, mpc_params.qwy;
    
    
    R_force.setZero();
    R_force.diagonal() << 0.01,0.01,0.01;//mpc_params.rx, mpc_params.ry, mpc_params.rz;
    for(int i=0;i<4;i++){
        R.block(3*i,3*i,3,3) = R_force;
    }
    

    casadi::DM Qc = casadi::DM::zeros(Q.rows(),Q.cols());
    std::copy(Q.data(), Q.data() + Q.size(), Qc.ptr());

    casadi::DM Rc = casadi::DM::zeros(R.rows(),R.cols());
    std::copy(R.data(), R.data() + R.size(), Rc.ptr());
    
    

    casadi::SX st = casadi::SX::zeros(NFSR,1);
    casadi::SX con = casadi::SX::zeros(NFIR,1);
    

    for(int i=0;i<HORIZ;i++){
        st = x(casadi::Slice(i*NFSR,i*NFSR+NFSR));
        con = x(casadi::Slice(NFSR*(HORIZ+1) + i*NFIR,NFSR*(HORIZ+1) + i*NFIR+NFIR));
        f= f+ mtimes(mtimes(transpose(st-x_des(casadi::Slice(i*NFSR,i*NFSR+NFSR))),casadi::SX(Qc)),st-x_des(casadi::Slice(i*NFSR,i*NFSR+NFSR)))
                + mtimes(mtimes(transpose(con),casadi::SX(Rc)),con);
    }

    st = x(casadi::Slice(HORIZ*NFSR,HORIZ*NFSR+NFSR));
    f = f+ mtimes(mtimes(transpose(st-x_des(casadi::Slice(HORIZ*NFSR,HORIZ*NFSR+NFSR))),casadi::SX(Qc)),st-x_des(casadi::Slice(HORIZ*NFSR,HORIZ*NFSR+NFSR)));
    return f;
}

casadi::SX SRBNMPCR::UpdateConstraintsN(casadi::SX x, casadi::SX p){
    // ====================================================================== //
    // ======================== Equality Constraints ======================== //
    // ====================================================================== //
    casadi::SX g_constraints = casadi::SX::zeros(NFSR*(HORIZ+1)+NFI*HORIZ,1);
    
    double mu = 0.8;//mpc_params.mu_MPC;
    
    casadi::SX contact_index, conind, foot;
    casadi::SX st,con,f_value,st_next,st_next_euler;

    g_constraints(casadi::Slice(0,NFSR)) = x(casadi::Slice(0,NFSR))-p(casadi::Slice(0,NFSR));
    

    for(int i=0; i<HORIZ; ++i){
        
        foot = p(casadi::Slice(NFSR*(HORIZ+1)+HORIZ*NFIR+4*HORIZ+4*(i),NFSR*(HORIZ+1)+HORIZ*NFIR+4*HORIZ+4*(i+1)));
        
        st = x(casadi::Slice(i*NFSR,i*NFSR+NFSR)); 
        con = x(casadi::Slice(NFSR*(HORIZ+1) + i*NFIR,NFSR*(HORIZ+1) + i*NFIR+NFIR));
        st_next = x(casadi::Slice((i+1)*NFSR,(i+1)*NFSR+NFSR));
        
        f_value = NonlinearDynamics(st,con,foot);
        st_next_euler = st + (f_value);
        
        g_constraints(casadi::Slice((i+1)*NFSR,(i+1)*NFSR+NFSR)) = st_next-st_next_euler; 
        g_constraints(casadi::Slice((HORIZ+1)*NFSR+i*NFI,(HORIZ+1)*NFSR+i*NFI+NFI)) = inequalitycons(con);
        
    }

    return g_constraints;

}

casadi::SX SRBNMPCR::NonlinearDynamics(casadi::SX st,casadi::SX con, casadi::SX conp1){//}, int i, casadi::MX p){

    casadi::SX rhs = casadi::SX::zeros(NFSR,1);

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
    

    casadi::SX tau = GetTorque(st,con,conp1);

    rhs(casadi::Slice(0,3)) = st(casadi::Slice(3,6)) * MPC_dt;
    rhs(casadi::Slice(3,6)) = (f1+f2+f3+f4)*MPC_dt/MASS - casadi::SX(gravityN) * MPC_dt;
    rhs(casadi::Slice(6,9)) = mtimes(A,st(casadi::Slice(9,12)))*MPC_dt;

    casadi::SX Jw = mtimes(casadi::SX(Jstandcasadi),st(casadi::Slice(9,12)));
    casadi::SX what = skewsym(st(casadi::Slice(9,12)));
    casadi::SX wJw = mtimes(what,Jw);
    rhs(casadi::Slice(9,12)) = mtimes(casadi::SX(Jinvcasadi),tau-wJw) * MPC_dt;

    return rhs;
}

casadi::SX SRBNMPCR::GetTorque(casadi::SX st,casadi::SX con, casadi::SX conp1){

    casadi::SX r1 = casadi::SX::zeros(3,1);
    casadi::SX r2 = casadi::SX::zeros(3,1);
    casadi::SX r3 = casadi::SX::zeros(3,1);
    casadi::SX r4 = casadi::SX::zeros(3,1);
    r1(0) = conp1(0)-st(0);
    r1(1) = -0.25-st(1);
    r1(2) = 0.533-st(2);

    r2(0) = conp1(1)-st(0);
    r2(1) = 0.25-st(1);
    r2(2) = 0.533-st(2);
    
    r3(0) = conp1(2)-st(0);
    r3(1) = -0.1321-st(1);
    r3(2) = 0-st(2);
    
    r4(0) = conp1(3)-st(0);
    r4(1) = 0.1321-st(1);
    r4(2) = 0-st(2);
    
    casadi::SX Mr1 = skewsym(r1);
    casadi::SX Mr2 = skewsym(r2);
    casadi::SX Mr3 = skewsym(r3);  
    casadi::SX Mr4 = skewsym(r4);
    
    return mtimes(Mr1,con(casadi::Slice(0,3)))+mtimes(Mr2,con(casadi::Slice(3,6)))+mtimes(Mr3,con(casadi::Slice(6,9)))+mtimes(Mr4,con(casadi::Slice(9,12)));

}

casadi::SX SRBNMPCR::inequalitycons(casadi::SX con){

       casadi::SX g_inequality = casadi::SX::zeros(NFI,1);

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

casadi::DM SRBNMPCR::lowerboundx(casadi::DM p){
    
    casadi::DM lbx = -casadi::DM::inf(NFSR*(HORIZ+1)+NFIR*HORIZ,1); 
    casadi::DM contact_index = casadi::DM::zeros(4,1); 
     
    for(int k=0 ; k<HORIZ ; k++){

        contact_index = p(casadi::Slice(NFSR*(HORIZ+1)+HORIZ*NFIR+4*(k),NFSR*(HORIZ+1)+HORIZ*NFIR+4*(k+1)));

        for(int leg=0 ; leg<4 ; leg++){
            if(leg<1){
                lbx(NFSR*(HORIZ+1)+NFIR*k+1) = 0;
            }else if(leg<2){
                lbx(NFSR*(HORIZ+1)+NFIR*(k)+4) = -contact_index(leg)*fzmaxf;   
            }else{
                lbx(NFSR*(HORIZ+1)+NFIR*(k)+3*(leg)+2) = 0;
            }
        }    
    }
    
    return lbx;

}

casadi::DM SRBNMPCR::upperboundx(casadi::DM p){
    
    casadi::DM ubx = casadi::DM::inf(NFSR*(HORIZ+1)+NFIR*HORIZ,1);
    casadi::DM contact_index = casadi::DM::zeros(4,1);

    for(int k=0 ; k<HORIZ ; k++){
        
        contact_index = p(casadi::Slice(NFSR*(HORIZ+1)+HORIZ*NFIR+4*(k),NFSR*(HORIZ+1)+HORIZ*NFIR+4*(k+1)));
        
        for(int leg=0 ; leg<4 ; leg++){
            if(leg<1){
                ubx(NFSR*(HORIZ+1)+NFIR*(k)+1) = contact_index(leg)*fzmaxf;
            }else if(leg<2){
                ubx(NFSR*(HORIZ+1)+NFIR*(k)+4) = 0;
            }else{
                ubx(NFSR*(HORIZ+1)+NFIR*(k)+3*(leg)+2) = contact_index(leg)*fzmaxr;
            }
        }
    }
    
    return ubx;
}

casadi::DM SRBNMPCR::lowerboundg(){
    
    casadi::DM lbg = casadi::DM::zeros(NFSR*(HORIZ+1)+NFI*HORIZ,1);
    
    for(int k=0 ; k<HORIZ ; k++){
        
        for(int leg=0 ; leg<4 ; leg++){     
            lbg(casadi::Slice(NFSR*(HORIZ+1)+NFI*k+4*leg,NFSR*(HORIZ+1)+NFI*k+4*leg+2)) = -casadi::DM::inf();    
        }
    }
    
    return lbg;
}

casadi::DM SRBNMPCR::upperboundg(){
    
    casadi::DM ubg = casadi::DM::zeros(NFSR*(HORIZ+1)+NFI*HORIZ,1);

    for(int k=0 ; k<HORIZ ; k++){
        
        for(int leg=0 ; leg<4 ; leg++){     
            ubg(casadi::Slice(NFSR*(HORIZ+1)+NFI*k+4*leg+2,NFSR*(HORIZ+1)+NFI*k+4*leg+4)) = casadi::DM::inf();    
        }
    }
    return ubg;
}

casadi::SX SRBNMPCR::skewsym(casadi::SX v3){

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


//Eigen::Matrix<double,12,1> SRBNMPCR::getOptforce(){
    
//    std::vector<double> v = previous_sol(casadi::Slice(NFSR*(HORIZ+1),NFSR*(HORIZ+1)+12)).get_elements();
//    Eigen::Matrix<double,12,1> optforce = Eigen::Map<Eigen::Matrix<double,12,1>>(v.data());
    
//    return optforce;

//}

Eigen::Matrix<double,12,1> SRBNMPCR::getFootPos(){
    
    Eigen::Matrix<double,12,1> foot_pos;

    for(int leg=0;leg<4;leg++){
        foot_pos(3*leg) = double(p_prev(NFSR*(HORIZ+1)+(NFIR+4)*(HORIZ)+leg));//.get_elements();
        if(leg<2){
            //foot_pos(3*leg) = 0.2;
            foot_pos(3*leg+1) = pow(-1,leg+1)*0.25;
            foot_pos(3*leg+2) = 0.533;
        }else{
            //foot_pos(3*leg) = -0.1;
            foot_pos(3*leg+1) = pow(-1,leg+1)*0.1321;
            foot_pos(3*leg+2) = 0;
        }
    }
    
    return foot_pos;

}

casadi::DM SRBNMPCR::getprevioussol(Eigen::Matrix<double,17,1> q0, size_t controlTick){
        
    casadi::DM x0 = casadi::DM::zeros(NFSR*(HORIZ+1)+NFIR*HORIZ,1);
    casadi::DM q0_dm(q0.rows(), q0.cols());
    for (int i=0; i<12; i++){
        q0_dm(i) = q0(i);
    }

    x0(casadi::Slice(0,12)) = q0_dm(casadi::Slice(0,12));
    if(controlTick<1){
        
        casadi::DM sum_conseq = casadi::DM::sum1(contact_sequence_dm); 
        
        for(int i = 0; i<HORIZ; i++){

            x0(casadi::Slice(NFSR*(i+1),NFSR*(i+2))) = x0(casadi::Slice(0,NFSR));
            
            for (size_t leg = 0; leg < 4; leg++)
            {
                x0(NFSR*(HORIZ+1)+NFIR*(i)+3*(leg+1)-1) = contact_sequence_dm(leg,i)*MASS*gravityN(2)/sum_conseq(i);
                if(leg<2){
                     x0(NFSR*(HORIZ+1)+NFIR*(i)+3*(leg+1)-2) = 2*pow(-1,leg)*x0(NFSR*(HORIZ+1)+NFIR*(i)+3*(leg+1)-1);//contact_sequence_dm(leg,i)*MASS*gravityN(2)/sum_conseq(i);
                }
                
            }
            
        }
    }else{
        x0(casadi::Slice(NFSR,NFSR*(HORIZ))) = previous_sol(casadi::Slice(NFSR*2,NFSR*(HORIZ+1)));
        x0(casadi::Slice(NFSR*HORIZ,NFSR*(HORIZ+1))) = previous_sol(casadi::Slice(NFSR*HORIZ,NFSR*(HORIZ+1)));

        x0(casadi::Slice(NFSR*(HORIZ+1),NFSR*(HORIZ+1)+NFIR*(HORIZ-1))) = previous_sol(casadi::Slice(NFSR*(HORIZ+1)+NFIR,NFSR*(HORIZ+1)+NFIR*(HORIZ)));
        x0(casadi::Slice(NFSR*(HORIZ+1)+NFIR*(HORIZ-1),NFSR*(HORIZ+1)+NFIR*(HORIZ))) = previous_sol(casadi::Slice(NFSR*(HORIZ+1)+NFIR*(HORIZ-1),NFSR*(HORIZ+1)+NFIR*(HORIZ)));

    }
    
    return x0;
}


void SRBNMPCR::writeMatrixToFile(const casadi::SX& matrix, const std::string& filename) {
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
        // file << matrix;
        // file << std::endl;
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

Eigen::Matrix<double,4,1> SRBNMPCR::getpreviousfoot(){
    
    Eigen::Matrix<double,4,1> p_prev_eig;
    for(int leg=0;leg<4;leg++){
        p_prev_eig(leg) = p_prev(NFSR*(HORIZ+1)+NFIR*HORIZ+4*HORIZ+leg).scalar();
    }
    
    return p_prev_eig;//(casadi::Slice((NFSR+NFIR)*HORIZ+NFSR+4*(HORIZ),(NFSR+NFIR)*HORIZ+NFSR+4*(HORIZ)+4));
};

void SRBNMPCR::getOptForceCoeff(int order){

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

Eigen::Matrix<double,12,1> SRBNMPCR::getOptforce(int phase,int order){

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

Eigen::Matrix<double,12,HORIZ> SRBNMPCR::arrangeOptforce(){
    
    Eigen::Matrix<double,12,HORIZ> optforce = Eigen::Matrix<double,12,HORIZ>::Zero();
    std::vector<double> v = previous_sol(casadi::Slice(static_cast<casadi_int>(NFSR*(HORIZ+1)), static_cast<casadi_int>(NFSR*(HORIZ+1)+12))).get_elements();

    optforce.block(0,0,12,1) = Eigen::Map<Eigen::Matrix<double,12,1>>(v.data());
    for (size_t i = 1; i < HORIZ; i++)
    {
        v = previous_sol(casadi::Slice(static_cast<casadi_int>(NFSR*(HORIZ+1)+NFIR*i), static_cast<casadi_int>(NFSR*(HORIZ+1)+NFIR*i+12))).get_elements();
        optforce.block(0,i,12,1) = Eigen::Map<Eigen::Matrix<double,12,1>>(v.data());
    }

    //Eigen::Matrix<double,12,1> optforce_vector = Eigen::Map<Eigen::Matrix<double,12,1>>(v.data());
    
    return optforce;

}

void SRBNMPCR::setpreviousp(casadi::DM p){ 
        //previousp = p;
        //std::cout << "p size: " << p.size1() << " " << p.size2() << std::endl;
        for (size_t i = 0; i < p.size1(); i++)
        {
            p_prev(i) = p(i);
          //  std::cout <<i <<std::endl;
        } 
};




