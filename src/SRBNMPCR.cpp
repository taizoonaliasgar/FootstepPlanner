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

    
    contact_sequence(1,casadi::Slice(5,20)) = casadi::MX::zeros(1,15);
    contact_sequence(2,casadi::Slice(8,20)) = casadi::MX::zeros(1,12);
    contact_sequence(0,casadi::Slice(25,40)) = casadi::MX::zeros(1,15);
    contact_sequence(3,casadi::Slice(28,40)) = casadi::MX::zeros(1,12);

    contact_sequence_dm(1,casadi::Slice(5,20)) = casadi::DM::zeros(1,15);
    contact_sequence_dm(2,casadi::Slice(8,20)) = casadi::DM::zeros(1,12);
    contact_sequence_dm(0,casadi::Slice(25,40)) = casadi::DM::zeros(1,15);
    contact_sequence_dm(3,casadi::Slice(28,40)) = casadi::DM::zeros(1,12);

    localvelocity = 0;
    
}


void SRBNMPCR::impactDetection(size_t tick, Eigen::Matrix<double,12,1> &q0, size_t gait){}

void SRBNMPCR::mpcdataLog(Eigen::Matrix<double,17,1> q0, casadi::DM p, casadi::DM X_prev, size_t tick){//(Eigen::Matrix<double,12,1> &q0, size_t tick){
    //size_t timeshift  = 0;
           
    //file[0] << std::fixed << std::setprecision(14) << std::endl;
    // file[0] << con[0] << "\t" << con[1] << "\t" << con[2] << "\t" << con[3] << "\t" << M+1 << "\t" << planindex-1 << "\t" << locomotionTick << "\t" << tick << std::endl;
    //casadi::DM contact_index = remainder(p(NFSR*(HORIZ+1)+NFIR*HORIZ),40);
    //file[1] << std::fixed << std::setprecision(14) << std::endl;   
    file[0] << tick << "\t" << q0(0) << "\t" << q0(1) << "\t" 
                << q0(2) << "\t" << q0(3) << "\t" 
                << q0(4) << "\t" << q0(5) << "\t" 
                << q0(6) << "\t" << q0(7) << "\t" 
                << q0(8) << "\t" << q0(9) << "\t" 
                << q0(10) << "\t" << q0(11) <<std::endl;//<< "\t"
                //<< q0(12) << "\t" << q0(13) << "\t" 
                //<< q0(14) << "\t" << q0(15) <<std::endl;

    //file[2] << std::fixed << std::setprecision(14) << std::endl;  
    file[1] << previous_sol(NFSR*(HORIZ+1)) << "\t" << previous_sol(NFSR*(HORIZ+1)+1) << "\t" 
                << previous_sol(NFSR*(HORIZ+1)+2) << "\t" << previous_sol(NFSR*(HORIZ+1)+3) << "\t" 
                << previous_sol(NFSR*(HORIZ+1)+4) << "\t" << previous_sol(NFSR*(HORIZ+1)+5) << "\t" 
                << previous_sol(NFSR*(HORIZ+1)+6) << "\t" << previous_sol(NFSR*(HORIZ+1)+7) << "\t" 
                << previous_sol(NFSR*(HORIZ+1)+8) << "\t" << previous_sol(NFSR*(HORIZ+1)+9) << "\t" 
                << previous_sol(NFSR*(HORIZ+1)+10) << "\t" << previous_sol(NFSR*(HORIZ+1)+11) << std::endl;

    file[2] <<  p(16-4) << "\t" << p(17-4) << "\t" 
                << p(18-4) << "\t" << p(19-4) << "\t" 
                << p(20-4) << "\t" << p(21-4) << "\t" 
                << p(22-4) << "\t" << p(23-4) << "\t" 
                << p(24-4) << "\t" << p(25-4) << "\t" 
                << p(26-4) << "\t" << p(27-4) << std::endl;
    //file[3] << std::fixed << std::setprecision(14) << std::endl;  
    file[3] << p(NFSR*(HORIZ+1)+NFIR*HORIZ) << "\t" << p(NFSR*(HORIZ+1)+NFIR*HORIZ+1)  << "\t" 
                << p(NFSR*(HORIZ+1)+NFIR*HORIZ+2) << "\t" << p(NFSR*(HORIZ+1)+NFIR*HORIZ+3)  << std::endl; 

    file[4] << p(NFSR*(HORIZ+1)+(NFIR+4)*HORIZ) << "\t" << p(NFSR*(HORIZ+1)+(NFIR+4)*HORIZ+1) << "\t" 
                << p(NFSR*(HORIZ+1)+(NFIR+4)*HORIZ+2) << "\t" << p(NFSR*(HORIZ+1)+(NFIR+4)*HORIZ+3) << std::endl;

    file[5] <<  X_prev(NFSR*(HORIZ+1))     << "\t" << X_prev(NFSR*(HORIZ+1) + 1) << "\t" << 
                X_prev(NFSR*(HORIZ+1) + 2) << "\t" << X_prev(NFSR*(HORIZ+1) + 3) << "\t" << 
                X_prev(NFSR*(HORIZ+1) + 4) << "\t" << X_prev(NFSR*(HORIZ+1) + 5) << "\t" << 
                X_prev(NFSR*(HORIZ+1) + 6) << "\t" << X_prev(NFSR*(HORIZ+1) + 7) << "\t" << 
                X_prev(NFSR*(HORIZ+1) + 8) << "\t" << X_prev(NFSR*(HORIZ+1) + 9) << "\t" << 
                X_prev(NFSR*(HORIZ+1) + 10) << "\t" << X_prev(NFSR*(HORIZ+1) + 11) << "\t" << std::endl;
    //            << contact_sequence_dm(2,contact_index) << "\t" << contact_sequence_dm(3,contact_index)  << "\t"
    //            << ubx(NFSR*(HORIZ+1)+1) << "\t" << ubx(NFSR*(HORIZ+1)+4) << "\t" << ubx(NFSR*(HORIZ+1)+8) << "\t" 
    //            << ubx(NFSR*(HORIZ+1)+11) << std::endl; 
    //             << q0(4) << "\t" << q0(5) << "\t" 
    //             << q0(6) << "\t" << q0(7) << "\t" 
    //             << q0(8) << "\t" << q0(9) << "\t" 
    //             << q0(10) << "\t" << q0(11) << "\t" <<std::endl;

    //file[4] << std::fixed << std::setprecision(14) << std::endl;  
    // file[4] << opttraj(0,0) << "\t" << opttraj(1,0) << "\t" 
    //             << opttraj(2,0) << "\t" << opttraj(3,0) << "\t" 
    //             << opttraj(4,0) << "\t" << opttraj(5,0) << "\t" 
    //             << opttraj(6,0) << "\t" << opttraj(7,0) << "\t" 
    //             << opttraj(8,0) <<  "\t" << opttraj(9,0) << "\t" 
    //             << opttraj(10,0) << "\t" << opttraj(11,0) << "\t" <<std::endl;
               
}



void SRBNMPCR::generator(){
    

    //casadi::SX p = casadi::SX::sym("p",(NFSR+NFIR)*HORIZ+NFSR+1,1);
    casadi::SX p = casadi::SX::sym("p",(NFSR+NFIR)*HORIZ+NFSR+8*(HORIZ),1);
    //casadi::SX p = casadi::SX::sym("p",(NFSR+NFIR)*HORIZ+NFSR+4*(HORIZ),1);
    //casadi::SX p = casadi::SX::sym("p",(NFSR+NFIR)*HORIZ+NFSR,1);
    casadi::SX x = casadi::SX::sym("x", (NFSR+NFIR)*HORIZ+NFSR);
    

    // Objective
    casadi::SX f = UpdateCostN(x,p);
    std::cout << "f"<<std::endl;
    std::string basePath = "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/build/tmp/";
    writeMatrixToFile(f, basePath + "obj.txt");
    
    // Constraints
    casadi::SX g = UpdateConstraintsN(x,p);
    std::cout << "g"<<std::endl;
    writeMatrixToFile(g, basePath + "const.txt");

    // Create an NLP solver instance
    casadi::Dict opts = {{"ipopt.print_level", 0}, {"print_time", 0},{"ipopt.max_iter", 100}};
    std::cout << "opts"<<std::endl;
    casadi::Function solver = casadi::nlpsol("solver", "ipopt", {{"x", x}, {"f", f}, {"g", g}, {"p", p}}, opts);
    //casadi::Function solver = casadi::nlpsol("solver", "ipopt", nlppp);
    std::cout << "solver"<<std::endl;
    int num_variables = solver.n_out();
    std::cout << num_variables << std::endl;

    // file name
    std::string file_name = "upright_nlp_2r";
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
    if(localvelocity < desVel(0)){
        localvelocity = localvelocity + 0.01;
    }else{
        localvelocity = desVel(0);
    }
    Raibheur = 0.5*Tstance*localvelocity;//desVel(0);
    
    contact_index = remainder(controlTick,40);
    for(size_t i= 0; i< HORIZ; i++){
        
        conp1 = remainder(contact_index+i,40);
        if(controlTick>0){
            conp1_prev = remainder(controlTick+i-1,40);
        }else{
            conp1_prev = 0;
        }

        x_des((i+1)*NFSR) = x_des(0) + (i+1)*localvelocity*MPC_dt; 
        x_des((i+1)*NFSR+1) = x_des(1) + (i+1)*desVel(1)*MPC_dt; 
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

        /*
        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i) = 0.2;//contact_sequence_dm(0,conp1);
        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+1) = 0.2;//contact_sequence_dm(1,conp1);
        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+2) = -0.1;//contact_sequence_dm(2,conp1);
        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+3) = -0.1;//contact_sequence_dm(3,conp1);
        */
        if(controlTick>1){
            for (int leg = 0;leg<4;leg++){
                if(contact_sequence_dm(leg,conp1).scalar() > contact_sequence_dm(leg,conp1_prev).scalar()){
                    if(leg<2){
                        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+leg) = x_des(NFSR*i)+2*Raibheur;//+0.2;//contact_sequence_dm(0,conp1);
                    }else{
                        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+leg) = x_des(NFSR*i)+Raibheur;//-0.1;
                    }    
                }else{
                    if(i<1){
                        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+leg) = q0_dm(12+leg);
                    }else{
                        x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+leg) = x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*(i-1)+leg);//x_des(NFSR*i);
                    }
                }
            }
        }else{
            x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i) = 0.2;//contact_sequence_dm(0,conp1);
            x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+1) = 0.2;//contact_sequence_dm(1,conp1);
            x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+2) = -0.1;//contact_sequence_dm(2,conp1);
            x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+4*i+3) = -0.1;//contact_sequence_dm(3,conp1);

        }
        
    }

    // x_des((HORIZ+1)*NFSR+HORIZ*NFIR) = controlTick;
    // conp1=remainder(contact_index+HORIZ,40);
    // x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ) = contact_sequence_dm(0,conp1);
    // x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+1) = contact_sequence_dm(1,conp1);
    // x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+2) = contact_sequence_dm(2,conp1);
    // x_des(NFSR*(HORIZ+1)+NFIR*(HORIZ)+4*HORIZ+3) = contact_sequence_dm(3,conp1);

    
    return x_des;
    
}


casadi::SX SRBNMPCR::UpdateCostN(casadi::SX x, casadi::SX x_des){
    
    casadi::SX f = 0;
    
    Eigen::Matrix3d R_force;
    Eigen::Matrix<double, NFSR, NFSR> Q;
    Eigen::Matrix<double, NFIR, NFIR> R;
    //Eigen::Matrix<double, 4, 4> Rv;
    
    //Eigen::MatrixXd Q_rep, R_rep, Q_cost;
    
    //Q_rep.setZero(NFSR*(HORIZ+1),NFSR*(HORIZ+1));
    //R_rep.setZero(NFIR*HORIZ,NFIR*HORIZ);
    
    Q.setZero();
    Q.block(0,0,3,3).diagonal() <<  1e4,1e4,8e5;//mpc_params.qpx, mpc_params.qpy, mpc_params.qpz;
    Q.block(3,3,3,3).diagonal() <<  1e4,1e4,1e4;//mpc_params.qvx, mpc_params.qvy, mpc_params.qvz;
    Q.block(6,6,3,3).diagonal() <<  8e4,8e5,3e3;//mpc_params.qrr, mpc_params.qrp, mpc_params.qry;
    Q.block(9,9,3,3).diagonal() <<  1,1,1;//mpc_params.qwr, mpc_params.qwp, mpc_params.qwy;
    //Q.block(12,12,4,4).diagonal() <<  1,1,1,1;
    //repdiag(Q,Q_rep,HORIZ+1);
    
    R_force.setZero();
    R_force.diagonal() << 0.01,0.01,0.01;//mpc_params.rx, mpc_params.ry, mpc_params.rz;
    for(int i=0;i<4;i++){
        R.block(3*i,3*i,3,3) = R_force;
        //R.block(12+i,12+i,1,1) << 0.01;
    }
    //Rv.setZero();
    //Rv.diagonal() << 0.01,0.01,0.01,0.01;

    casadi::DM Qc = casadi::DM::zeros(Q.rows(),Q.cols());
    std::copy(Q.data(), Q.data() + Q.size(), Qc.ptr());

    casadi::DM Rc = casadi::DM::zeros(R.rows(),R.cols());
    std::copy(R.data(), R.data() + R.size(), Rc.ptr());
    
    //casadi::DM Rvc = casadi::DM::zeros(Rv.rows(),Rv.cols());
    //std::copy(Rv.data(), Rv.data() + Rv.size(), Rvc.ptr());

    casadi::SX st,con;
    //casadi::SX Raibheurcol = 0.5*Tstance*desVel(0)*casadi::SX::ones(4,1);

    for(int i=0;i<HORIZ;i++){
        st = x(casadi::Slice(i*NFSR,i*NFSR+NFSR));
        con = x(casadi::Slice(NFSR*(HORIZ+1) + i*NFIR,NFSR*(HORIZ+1) + i*NFIR+NFIR));
        f= f+ mtimes(mtimes(transpose(st-x_des(casadi::Slice(i*NFSR,i*NFSR+NFSR))),casadi::SX(Qc)),st-x_des(casadi::Slice(i*NFSR,i*NFSR+NFSR)))
                + mtimes(mtimes(transpose(con),casadi::SX(Rc)),con);
                    //+mtimes(mtimes(transpose(con(casadi::Slice(12,16))-Raibheurcol),casadi::SX(Rvc)),con(casadi::Slice(12,16))-Raibheurcol);
    }
    
    // repdiag(R,R_rep,HORIZ);
    
    // Q_cost = Eigen::MatrixXd::Zero(NFSR*(HORIZ+1)+NFIR*(HORIZ),NFSR*(HORIZ+1)+NFIR*(HORIZ));
    
    // Q_cost.block(0,0,NFSR*(HORIZ+1),NFSR*(HORIZ+1)) = Q_rep;
    // Q_cost.block(NFSR*(HORIZ+1),NFSR*(HORIZ+1),NFIR*(HORIZ),NFIR*(HORIZ)) = R_rep;

    // casadi::DM Qcasadi = casadi::DM::zeros(Q_cost.rows(),Q_cost.cols());
    // std::copy(Q_cost.data(), Q_cost.data() + Q_cost.size(), Qcasadi.ptr());
    
    // casadi::SX xT = transpose(x)-transpose(x_des(casadi::Slice(0,NFSR*(HORIZ+1)+NFIR*(HORIZ))));
    // casadi::SX xTQ = mtimes(xT,casadi::SX(Qcasadi));
    // f = mtimes(xTQ, x-x_des(casadi::Slice(0,NFSR*(HORIZ+1)+NFIR*(HORIZ))));
    //f = mtimes(mtimes(transpose(x)-transpose(x_des),Qcasadi), x-x_des);

    return f;
}

casadi::SX SRBNMPCR::UpdateConstraintsN(casadi::SX x, casadi::SX p){
    // ====================================================================== //
    // ======================== Equality Constraints ======================== //
    // ====================================================================== //
    casadi::SX g_constraints = casadi::SX::zeros(NFSR*(HORIZ+1)+NFI*HORIZ,1);
    
    double mu = 0.8;//mpc_params.mu_MPC;
    
    casadi::SX contact_index, conind, conp1;
    //casadi::DM conp1;
    casadi::SX st,con,f_value,st_next,st_next_euler;
   
    //casadi::SX p_control = p(NFSR*(HORIZ+1)+NFIR*(HORIZ));

    g_constraints(casadi::Slice(0,NFSR)) = x(casadi::Slice(0,NFSR))-p(casadi::Slice(0,NFSR));
    

    for(int i=0; i<HORIZ; ++i){
        
        //contact_index = remainder((p_control+i),40);
        //conp1 = remainder((contact_index+1),40);
        //conp1 = casadi::DM(remainder((p_control+i+1),40));
        //std::cout << "conp1: " << conp1 << std::endl;

        conp1 = p(casadi::Slice(NFSR*(HORIZ+1)+HORIZ*NFIR+4*HORIZ+4*(i),NFSR*(HORIZ+1)+HORIZ*NFIR+4*HORIZ+4*(i+1)));//remainder((contact_index+1),40);
        //conp1 = casadi::SX::zeros(4,1);
        
        st = x(casadi::Slice(i*NFSR,i*NFSR+NFSR));
         
        con = x(casadi::Slice(NFSR*(HORIZ+1) + i*NFIR,NFSR*(HORIZ+1) + i*NFIR+NFIR));
        
        st_next = x(casadi::Slice((i+1)*NFSR,(i+1)*NFSR+NFSR));
        
        f_value = NonlinearDynamics(st,con,conp1); 
        //f_value = NonlinearDynamics(st,con);//,i,p);  

        st_next_euler = st + (f_value);
        
        g_constraints(casadi::Slice((i+1)*NFSR,(i+1)*NFSR+NFSR)) = st_next-st_next_euler;
        
        g_constraints(casadi::Slice((HORIZ+1)*NFSR+i*NFI,(HORIZ+1)*NFSR+i*NFI+NFI)) = inequalitycons(con);
        

    }

    return g_constraints;

}

casadi::SX SRBNMPCR::NonlinearDynamics(casadi::SX st,casadi::SX con, casadi::SX conp1){//}, int i, casadi::MX p){

    casadi::SX rhs = casadi::SX::zeros(NFSR,1);
    //std::cout << "rhs" << std::endl;
    // casadi::SX contact = casadi::SX::zeros(4,1);
    // contact(0) = contact_sequence_dm(0,conp1);
    // contact(1) = contact_sequence_dm(1,conp1);
    // contact(2) = contact_sequence_dm(2,conp1);
    // contact(3) = contact_sequence_dm(3,conp1);
    //std::cout << "contact" << std::endl;
    
    //contact = p(casadi::Slice(NFSR*(HORIZ+1)+NFIR*HORIZ+(i)*4,NFSR*(HORIZ+1)+NFIR*HORIZ+(i+1)*4));

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
    // A(0,0) = cos(psi)/cos(theta);
    // A(0,1) = sin(psi)/cos(theta);
    // A(0,2) = 0;
    // A(1,0) = -sin(psi);
    // A(1,1) =cos(psi);
    // A(1,2) =0;
    // A(2,0) =cos(psi)*tan(theta);
    // A(2,1) =sin(psi)*tan(theta);
    // A(2,2) = 1;

    casadi::SX tau = GetTorque(st,con,conp1);

    rhs(casadi::Slice(0,3)) = st(casadi::Slice(3,6)) * MPC_dt;

    rhs(casadi::Slice(3,6)) = (f1+f2+f3+f4)*MPC_dt/MASS - casadi::SX(gravityN) * MPC_dt;

    rhs(casadi::Slice(6,9)) = mtimes(A,st(casadi::Slice(9,12)))*MPC_dt;

    casadi::SX Jw = mtimes(casadi::SX(Jstandcasadi),st(casadi::Slice(9,12)));
    casadi::SX what = skewsym(st(casadi::Slice(9,12)));
    casadi::SX wJw = mtimes(what,Jw);
    rhs(casadi::Slice(9,12)) = mtimes(casadi::SX(Jinvcasadi),tau-wJw) * MPC_dt;//mtimes(skewsym(st(casadi::Slice(9,12)),mtimes(casadi::SX(Jstandcasadi),st(casadi::Slice(9,12))))))*MPC_dt;

    //rhs(12) = con(12)*conp1(0);//contact(0);
    //rhs(13) = con(13)*conp1(1);//contact(1);
    //rhs(14) = con(14)*conp1(2);//contact(2);
    //rhs(15) = con(15)*conp1(3);//contact(3);

    return rhs;
}

casadi::SX SRBNMPCR::GetTorque(casadi::SX st,casadi::SX con, casadi::SX conp1){

    casadi::SX r1 = casadi::SX::zeros(3,1);
    casadi::SX r2 = casadi::SX::zeros(3,1);
    casadi::SX r3 = casadi::SX::zeros(3,1);
    casadi::SX r4 = casadi::SX::zeros(3,1);
    r1(0) = conp1(0)-st(0);//conp1(0)
    r1(1) = -0.25-st(1);
    r1(2) = 0.533-st(2);

    r2(0) = conp1(1)-st(0);//conp1(1)
    r2(1) = 0.25-st(1);
    r2(2) = 0.533-st(2);
    
    r3(0) = conp1(2)-st(0);//conp1(2)
    r3(1) = -0.1321-st(1);
    r3(2) = 0-st(2);
    
    r4(0) = conp1(3)-st(0);//conp1(3)
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
    //casadi::DM p_control = p(NFSR*(HORIZ+1)+NFIR*(HORIZ));
     

    for(int k=0 ; k<HORIZ ; k++){
        //casadi::DM contact_index = remainder((p_control+k),40);
        casadi::DM contact_index = p(casadi::Slice(NFSR*(HORIZ+1)+HORIZ*NFIR+4*(k),NFSR*(HORIZ+1)+HORIZ*NFIR+4*(k+1)));
        for(int leg=0 ; leg<4 ; leg++){
            if(leg<1){
                lbx(NFSR*(HORIZ+1)+NFIR*k+1) = 0;

            }else if(leg<2){
                lbx(NFSR*(HORIZ+1)+NFIR*(k)+4) = -contact_index(leg)*fzmaxf;//contact_index(leg)*//-contact_sequence_dm(leg,contact_index)*fzmaxf;//
                
            }else{
                lbx(NFSR*(HORIZ+1)+NFIR*(k)+3*(leg)+2) = 0;
        
            }
        }
        //lbx(casadi::Slice(NFSR*(HORIZ+1)+NFIR*k+12,NFSR*(HORIZ+1)+NFIR*k+16)) = casadi::DM::zeros(4,1);
    }
    
    return lbx;

}

casadi::DM SRBNMPCR::upperboundx(casadi::DM p){
    
    casadi::DM ubx = casadi::DM::inf(NFSR*(HORIZ+1)+NFIR*HORIZ,1);
    //casadi::DM p_control = p(NFSR*(HORIZ+1)+NFIR*(HORIZ));

    for(int k=0 ; k<HORIZ ; k++){
        //casadi::DM contact_index = remainder((p_control+k),40);
        //casadi::DM concol = contact_sequence_dm(casadi::Slice(0,4),contact_index);
        casadi::DM contact_index = p(casadi::Slice(NFSR*(HORIZ+1)+HORIZ*NFIR+4*(k),NFSR*(HORIZ+1)+HORIZ*NFIR+4*(k+1)));
        
        for(int leg=0 ; leg<4 ; leg++){
            if(leg<1){
                ubx(NFSR*(HORIZ+1)+NFIR*(k)+1) = contact_index(leg)*fzmaxf;//contact_index(leg)*//contact_sequence_dm(leg,contact_index)*fzmaxf;//
            }else if(leg<2){
                ubx(NFSR*(HORIZ+1)+NFIR*(k)+4) = 0;
            }else{
                ubx(NFSR*(HORIZ+1)+NFIR*(k)+3*(leg)+2) = contact_index(leg)*fzmaxr;//contact_index(leg)*//contact_sequence_dm(leg,contact_index)*fzmaxr;//
            }
        }
        //ubx(casadi::Slice(NFSR*(HORIZ+1)+NFIR*(k)+12,NFSR*(HORIZ+1)+NFIR*(k)+16)) = (casadi::DM::ones(4,1)-contact_index)*RaibMult*Raibheur; // -concol Use the created casadi::Slice object
    }
    
    return ubx;
}

casadi::DM SRBNMPCR::lowerboundg(){
    
    casadi::DM lbg = casadi::DM::zeros(NFSR*(HORIZ+1)+NFI*HORIZ,1);
    
    //lbg(casadi::Slice(0,NFSR*(HORIZ+1)))=-0.001*casadi::DM::ones(NFSR*(HORIZ+1),1);
    for(int k=0 ; k<HORIZ ; k++){
        
        for(int leg=0 ; leg<4 ; leg++){     
            lbg(casadi::Slice(NFSR*(HORIZ+1)+NFI*k+4*leg,NFSR*(HORIZ+1)+NFI*k+4*leg+2)) = -casadi::DM::inf();    
        }
    }
    
    return lbg;
}

casadi::DM SRBNMPCR::upperboundg(){
    
    casadi::DM ubg = casadi::DM::zeros(NFSR*(HORIZ+1)+NFI*HORIZ,1);
    //ubg(casadi::Slice(0,NFSR*(HORIZ+1)))=0.001*casadi::DM::ones(NFSR*(HORIZ+1),1);

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


Eigen::Matrix<double,12,1> SRBNMPCR::getOptforce(){
    
    std::vector<double> v = previous_sol(casadi::Slice(NFSR*(HORIZ+1),NFSR*(HORIZ+1)+12)).get_elements();
    Eigen::Matrix<double,12,1> optforce = Eigen::Map<Eigen::Matrix<double,12,1>>(v.data());
    
    return optforce;

}

Eigen::Matrix<double,12,1> SRBNMPCR::getFootPos(casadi::DM pp){
    
    Eigen::Matrix<double,12,1> foot_pos;

    //std::vector<double> v = previous_sol(casadi::Slice(12,16)).get_elements();
    //Eigen::Matrix<double,4,1> foot_x = Eigen::Map<Eigen::Matrix<double,4,1>>(v.data());

    for(int leg=0;leg<4;leg++){
        foot_pos(3*leg) = double(pp(NFSR*(HORIZ+1)+(NFIR+4)*(HORIZ)+leg));//.get_elements();
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
        //x0(12) = 0.2;
        //x0(13) = 0.2;
        //x0(14) = -0.1;
        //x0(15) = -0.1;
        casadi::DM sum_conseq = casadi::DM::sum1(contact_sequence_dm); 
        
        for(int i = 0; i<HORIZ; i++){
            x0(casadi::Slice(NFSR*(i+1),NFSR*(i+2))) = x0(casadi::Slice(0,NFSR));
            
            for (size_t leg = 0; leg < 4; leg++)
            {
                // x0(NFSR*(HORIZ+1)+NFIR*(i)+3*(leg+1)-1) = contact_sequence_dm(leg,i)*MASS*gravityN(2)/sum_conseq(i);
                // if(leg<2){
                //     x0(NFSR*(HORIZ+1)+NFIR*(i)+3*(leg+1)-2) = 2*pow(-1,leg)*x0(NFSR*(HORIZ+1)+NFIR*(i)+3*(leg+1)-1);//contact_sequence_dm(leg,i)*MASS*gravityN(2)/sum_conseq(i);
                // }
                if(leg>1){
                   x0(NFSR*(HORIZ+1)+NFIR*(i)+3*(leg+1)-1) = 44; 
                }else{
                    x0(NFSR*(HORIZ+1)+NFIR*(i)+3*(leg+1)-2) = pow(-1,leg)*22;
                    x0(NFSR*(HORIZ+1)+NFIR*(i)+3*(leg+1)-1) = 17.5;
                }
            }
            
        }
    }else{
        //x0(casadi::Slice(12,16)) = previous_sol(casadi::Slice(28,32));
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









