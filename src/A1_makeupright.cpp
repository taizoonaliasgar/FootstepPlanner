//
// Authror: Randy Fawcett on 12/2021.
//
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "raisim/OgreVis.hpp"
//#include "randyImguiPanel.hpp"
#include "raisimBasicImguiPanel.hpp"
#include "raisimKeyboardCallback.hpp" // THIS IS WHERE cmd.vel AND cmd.pose ARE IMPLEMENTED (as globals)
#include "raisim/RaisimServer.hpp"
#include "helper.hpp"
//#include "helper2.hpp"
#include "Filters.h"
#include <iostream>
#include <filesystem>
#include <fstream>

//HDSRL header
//#include "/home/taizoon/raisimEnv/raisimWorkspace/flying_trot/include/LocoWrapperfly.hpp"
//#include "/home/taizoon/raisimEnv/raisimWorkspace/flying_trot/include/locomotion_planner.h"
//#include "/home/taizoon/raisimEnv/raisimWorkspace/flying_trot/include/LocoWrapperwalk.hpp"
#include "LocoWrapper.hpp"
#include "SRBNMPC.hpp"
using namespace std;
// namespace fs = std::filesystem;

FiltStruct_f* filt = (FiltStruct_f*)malloc(sizeof(FiltStruct_f));

void comparingfootv(raisim::Vec<3> FRfoot,raisim::Vec<3> FLfoot,raisim::Vec<3> RRfoot,raisim::Vec<3> RLfoot,Eigen::Matrix<double, 3, 4> foot_vel){
    std::cout << FRfoot[0]-foot_vel(0,0) << "\t" << FRfoot[1]-foot_vel(1,0) << "\t" << FRfoot[2]-foot_vel(2,0) << "\t" <<
                    FLfoot[0]-foot_vel(0,1) << "\t" << FLfoot[1]-foot_vel(1,1) << "\t" << FLfoot[2]-foot_vel(2,1) << "\t" <<
                        RRfoot[0]-foot_vel(0,2) << "\t" << RRfoot[1]-foot_vel(1,2) << "\t" << RRfoot[2]-foot_vel(2,2) << "\t" <<
                            RLfoot[0]-foot_vel(0,3) << "\t" << RLfoot[1]-foot_vel(1,3) << "\t" << RLfoot[2]-foot_vel(2,3) << std::endl;

}

double distx(size_t controlTick, double amp){
    
    return amp*sin((controlTick-8000)/1000.0/M_PI);
    //return dist;
};

void setupCallback() {
    raisim::OgreVis *vis = raisim::OgreVis::get();

    /// light
    vis->getLight()->setDiffuseColour(1, 1, 1);
    vis->getLight()->setCastShadows(false);
    Ogre::Vector3 lightdir(-3,3,-0.5); // Light shines on ROBOTS top/front/right side
    // Ogre::Vector3 lightdir(-3,-3,-0.5); // Light shines on ROBOTS top/front/left side
    lightdir.normalise();
    vis->getLightNode()->setDirection({lightdir});
    vis->setCameraSpeed(300);

    vis->addResourceDirectory(raisim::loadResource("material"));
    vis->loadMaterialFile("myMaterials.material");

    vis->addResourceDirectory(vis->getResourceDir() + "/material/skybox/violentdays");
    vis->loadMaterialFile("violentdays.material");

    /// shdow setting
    vis->getSceneManager()->setShadowTechnique(Ogre::SHADOWTYPE_TEXTURE_ADDITIVE);
    vis->getSceneManager()->setShadowTextureSettings(2048, 3);

    /// scale related settings!! Please adapt it depending on your map size
    // beyond this distance, shadow disappears
    vis->getSceneManager()->setShadowFarDistance(10);
    // size of contact points and contact forces
    vis->setContactVisObjectSize(0.03, 0.6);
    // speed of camera motion in freelook mode
    vis->getCameraMan()->setTopSpeed(5);
}

void controller(std::vector<raisim::ArticulatedSystem *> A1, LocoWrapper *loco_obj, SRBNMPC* loco_plan, casadi::Function solver, size_t controlTick, raisim::Contact contactInstance) {
    /////////////////////////////////////////////////////////////////////
    //////////////////////////// INITIALIZE
    /////////////////////////////////////////////////////////////////////
    static size_t loco_kind = TROT;                        // Gait pattern to use
    size_t pose_kind = POSE_CMD;                   // Pose type to use (if loco_kind is set to POSE)
    size_t settling = 0.2*ctrlHz;                   // Settling down
    size_t duration = 1.8*ctrlHz;                   // Stand up 
    size_t loco_start = settling + duration;        // Start the locomotion pattern

    size_t shifttime = 0.75*ctrlHz;
    size_t movetime = 0.3*ctrlHz;
    size_t shifttime2 = 0.6*ctrlHz;
    
    double *tau;
    //double tau[18] = {0};
    double jpos[18], jvel[18];
    
    Eigen::VectorXd jointTorqueFF = Eigen::MatrixXd::Zero(TOTAL_DOF,1);
    Eigen::VectorXd jointPosTotal = Eigen::MatrixXd::Zero(TOTAL_DOF+1,1);
    Eigen::VectorXd jointVelTotal = Eigen::MatrixXd::Zero(TOTAL_DOF,1);    
        
    raisim::Mat<3,3> rotMat;
    Eigen::Matrix<double, 3, 1> eul;
    Eigen::Matrix<double, 4, 1> quat;

    Eigen::Matrix<double, 3, 1> imu_eul = Eigen::MatrixXd::Zero(3,1);
    Eigen::Matrix<double, 4, 1> imu_quat = Eigen::MatrixXd::Zero(4,1);
    Eigen::Matrix<double, 3, 1> imu_omega = Eigen::MatrixXd::Zero(3,1);

    //For Taizoon High level
    Eigen::Matrix<double, 12, 1> QP_Force;
    const int* foot_state;
    Eigen::Matrix<double, 3, 4> foot_position;
    Eigen::Matrix<double, 3, 4> hip_position;
    Eigen::Matrix<double, 33, 1> opt_HLMPC_state;
    
    /////////////////////////////////////////////////////////////////////
    //////////////////////////// UPDATE STATE
    /////////////////////////////////////////////////////////////////////
    A1.back()->getState(jointPosTotal, jointVelTotal);
    A1.back()->getBaseOrientation(rotMat);
     
    auto imu = A1.back()->getSensorSet("imu_parent")->getSensor<raisim::InertialMeasurementUnit>("imu");
    raisim::Vec<3> linearAcceleration = imu->getLinearAcceleration();
    Eigen::Matrix<double,3,1> acc_bFrame = Eigen::MatrixXd::Zero(3,1);
    acc_bFrame(0) = linearAcceleration(0);
    acc_bFrame(1) = linearAcceleration(1);
    acc_bFrame(2) = linearAcceleration(2);

    auto imu_o = imu->getOrientation();   // Quaternion
    auto imu_w = imu->getAngularVelocity();  // Angular velocity in radians/s
    
    // std::cout << quat_base[0] << "\t" << quat_base[1] << "\t" << quat_base[2] << "\t" << quat_base[3] << "\t" <<
    //                 imu_o[0] << "\t"    << imu_o[1]     << "\t" << imu_o[2]     << "\t"  << imu_o[3] << std::endl;

    imu_quat(0) = imu_o[0];
    imu_quat(1) = imu_o[1];
    imu_quat(2) = imu_o[2];
    imu_quat(3) = imu_o[3];

    imu_omega(0) = imu_w[0];
    imu_omega(1) = imu_w[1];
    imu_omega(2) = imu_w[2];

    quat_to_XYZ(imu_quat,imu_eul);

    raisim::MatDyn D = A1.back()->getMassMatrix();
    raisim::VecDyn H = A1.back()->getNonlinearities({0,0,-9.81});

    Eigen::Matrix<double,18,18> Dr;
    Eigen::Matrix<double,18,1> Hr;
    

    for(int i=0;i<18;i++){
        for(int j=0;j<18;j++){
            Dr(i,j) = D(i,j);
        }
        Hr(i)=H(i,1);
    }
    
    double rotMatrixDouble[9];
    for(size_t i=0;i<9;i++){
        rotMatrixDouble[i] = rotMat[i];
    }
    
    
    Eigen::Map< Eigen::Matrix<double, 3, 3> > rotE(rotMatrixDouble, 3, 3);
    jointVelTotal.segment(3,3) = rotE.transpose()*jointVelTotal.segment(3,3); // convert to body frame, like robot measurements

    quat = jointPosTotal.block(3,0,4,1);
    quat_to_XYZ(quat,eul);
    
    for(size_t i=0; i<3; ++i){
        jpos[i] = jointPosTotal(i);
        jvel[i] = jointVelTotal(i);
        jpos[i+3] = eul(i);
        jvel[i+3] = jointVelTotal(i+3);
    }

    for(size_t i=6; i<18; ++i){
        jpos[i] = jointPosTotal(i+1);
        jvel[i] = jointVelTotal(i);
    }

    // std::cout << jpos[0] << "\t" << jpos[1] << "\t" << jpos[2] << "\t" << jvel[0] << "\t" << jvel[1] << "\t" << jvel[2] << "\t"
    //                 << q_est(0) << "\t" << q_est(1) << "\t" << q_est(2) << "\t" << q_est(3) << "\t" << q_est(4) << "\t" << q_est(5) << "\t"
    //                                              << trunk_acc(0) << "\t" << trunk_acc(1) << "\t" << trunk_acc(2) << "\t" 
    //                                                 << acc_wFrame(0) << "\t" << acc_wFrame(1) << "\t" << acc_wFrame(2) << "\t" 
    //                                                 << acc_bFrame(0) << "\t" << acc_bFrame(1) << "\t" << acc_bFrame(2) << "\t"
    //                                                 << jpos[3] << "\t" << jpos[4] << "\t" << jpos[5] << "\t" 
    //                                                 << jvel[3] << "\t" << jvel[4] << "\t" << jvel[5] << "\t"
    //                                                 << q_est(6) << "\t" << q_est(7) << "\t" << q_est(8) << "\t" 
    
    // std::cout << controlTick << "\t" << acc_bFrame(0) << "\t" << acc_bFrame(1) << "\t" << acc_bFrame(2) << "\t"
    //                                                << imu_eul[0] << "\t" << imu_eul[1] << "\t" << imu_eul[2] << "\t"
    //                                                 << eul(0) << "\t" << eul(1) << "\t" << eul(2) //<< "\t"
    //                                                 //<< q_est(6) << "\t" << q_est(7) << "\t" << q_est(8) << "\t" 
    //                                                 << std::endl; 

    Eigen::Matrix<double,16,1> q0;
    q0.setZero(16,1);
    q0.block(0,0,3,1) << jpos[0],jpos[1],jpos[2];//= jointPosTotal.block(0,0,3,1);
    q0.block(3,0,3,1) << jvel[0],jvel[1],jvel[2];//= jointVelTotal.block(0,0,3,1);
    q0.block(6,0,3,1) << jpos[3],jpos[4],jpos[5];
    q0.block(9,0,3,1) << jvel[3],jvel[4],jvel[5];//.block(3,0,3,1);
    std::map<std::string, casadi::DM> arg, res;
    
    int force[4] = {0};
    for(auto &con: A1.back()->getContacts()){
        int conInd = con.getlocalBodyIndex();
        force[conInd/3-1] = 500;
        force[conInd/3-1] = con.getNormal().e().norm();
    }

    float vel_temp[3] = {0,0,0};//{cmd.vel[0],cmd.vel[1],cmd.vel[2]};
    float pose_temp[6] = {0,0,0,0,0,0};
    float filt_vel_temp[3] = {jvel[0],jvel[1],jvel[2]};
    discrete_butter_f(filt,filt_vel_temp);
    int duration_data =0;
    Eigen::Matrix<double,4,1> nextcon = Eigen::MatrixXd::Ones(4,1);
    int maxsteps = 14;
    
    /////////////////////////////////////////////////////////////////////
    //////////////////////////// CONTROL
    /////////////////////////////////////////////////////////////////////
    // Update the desired torques
    if(controlTick < settling){ // Settle down
        // loco_obj->posSetup(jointPosTotal.head(3));
        double temp[18] = {0};
        tau = temp;
        loco_obj->initStandVars(jointPosTotal.block(0,0,3,1),jointPosTotal(5),(int)duration);
    }
    else if(controlTick >= settling & controlTick < loco_start){ // Start standing
        loco_obj->calcTau2(jpos,jvel,rotMatrixDouble,STAND,controlTick,loco_start);//,shifttime,movetime,shifttime2,0,0);
        tau = loco_obj->getTorque();

    }
    else if(controlTick >= loco_start & controlTick < loco_start+shifttime){
        
        nextcon(0) = 0;
        
        if(controlTick==loco_start){
            Eigen::Matrix<double, 4, 1> wfoot = 6*Eigen::MatrixXd::Ones(4,1);
            wfoot(0)=1;//nextcon(0);
            wfoot(1)=1;//nextcon(1);
            loco_obj->getshiftedCoM(wfoot);
            loco_obj->setshiftedCoM();
        }

        loco_obj->setswingContact(nextcon);
        loco_obj->calcTau2(jpos,jvel,rotMatrixDouble,STANDUP,controlTick,loco_start);//,shifttime,movetime,shifttime2,0,0);
        tau = loco_obj->getTorque();

    }else if(controlTick >= loco_start+shifttime){// & controlTick < loco_start + shifttime){ // Start locomotion
        
        int stepind = std::floor((controlTick-loco_start-shifttime)/(shifttime+movetime));
        
        if(stepind<maxsteps){
            loco_obj->stepsonwall(stepind);
        
            if(stepind==1){
                loco_obj->tookfirststep();
                //std::cout << "First step" << "\t" << controlTick << std::endl; 
            }
        
            if(stepind%2==0){
                nextcon(0) = 0;
                if(controlTick == loco_start+shifttime + stepind*(shifttime+movetime)){
                    loco_obj->incstep();
                }
                
            }else{
                nextcon(1) = 0;
            }

            
        
         
            if(controlTick==loco_start + (stepind+1)*(movetime + shifttime)){// + stepind*(shifttime+movetime+shifttime2)){
                Eigen::Matrix<double, 4, 1> wfoot = 6*Eigen::MatrixXd::Ones(4,1);//3
                wfoot(0)=1;//nextcon(1);
                wfoot(1)=1;//nextcon(0);
                //wfoot(2)=1+1*wfoot(0);
                //wfoot(3)=1+1*wfoot(1);
                loco_obj->getshiftedCoM(wfoot);
                loco_obj->setshiftedCoM();
                //std::cout << "controlTick:" << "\t" << controlTick << "\t" << "stepind" << "\t" << stepind << std::endl;
            }
            // else if(controlTick == loco_start + (stepind)*(movetime + shifttime) + shifttime && stepind>0){
            //     Eigen::Matrix<double, 4, 1> wfoot = 6*Eigen::MatrixXd::Ones(4,1);//3
            //     wfoot(0)=1;//nextcon(0);
            //     wfoot(1)=1;//nextcon(1);
            //     //wfoot(2)=2+2*wfoot(0);
            //     //wfoot(3)=2+2*wfoot(1);
            //     loco_obj->getshiftedCoM(wfoot);
            //     loco_obj->setshiftedCoM();
            // }

            loco_obj->setswingContact(nextcon);
            loco_obj->calcTau2(jpos,jvel,rotMatrixDouble,STANDUP,controlTick,loco_start);//,shifttime,movetime,shifttime2,0,0);
            tau = loco_obj->getTorque();
        
        
        }else{
  
            
            if(controlTick==loco_start + shifttime + (maxsteps)*(movetime + shifttime)){// + stepind*(shifttime+movetime+shifttime2)){
                loco_obj->stopclimbing();
                loco_obj->stepsonwall(maxsteps);
                Eigen::Matrix<double, 4, 1> wfoot = 6*Eigen::MatrixXd::Ones(4,1);
                wfoot(0)=1;
                wfoot(1)=1;
                //wfoot(2)=1+1*wfoot(0);
                //wfoot(3)=1+1*wfoot(1);
                loco_obj->getshiftedCoM(wfoot);
                loco_obj->setshiftedCoM();
                // nextcon(0) = 1;
                // nextcon(1) = 0;
                // loco_obj->setswingContact(nextcon);
                //std::cout << "controlTick:" << "\t" << controlTick << "\t" << "stepind" << "\t" << stepind << std::endl; 
            }
            
            loco_obj->calcTau2(jpos,jvel,rotMatrixDouble,STANDUP,controlTick,loco_start);//,shifttime,movetime,shifttime2,0,0);
            tau = loco_obj->getTorque();
        }

        // loco_obj->updateVel(vel_temp);
        // loco_obj->updatePose(pose_temp);
        // if(controlTick%10==0){
        //     int controlMPC = std::floor(controlTick/10);
        //     loco_obj->setcontactconfig(controlMPC);
        // }
        
    }

    
    jointTorqueFF = Eigen::Map< Eigen::Matrix<double,18,1> >(tau,18);
    jointTorqueFF.block(0,0,6,1).setZero();
    
    // Set the desired torques
    A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
    A1.back()->setGeneralizedForce(jointTorqueFF);

    // Apply the disturbance force to the desired body
    // Assuming the body index is 0, adjust as needed
    // if(controlTick > 500 && controlTick < 1000){

    //     double amp = 600.0;
    //     //double dist = distx(controlTick, amp);
    //     Eigen::Vector3d disturbanceForce(0,amp,0); 
    //     A1.back()->setExternalForce(0, disturbanceForce);
    // }

};



int main(int argc, char *argv[]) {
    
    //long int terrain_number = atoll(argv[3]); // Use atoll for long long int

    // ============================================================ //
    // =================== SETUP RAISIM/VISUALS =================== //
    // ============================================================ //
    /// create raisim world
    raisim::World::setActivationKey(raisim::loadResource("activation.raisim"));
    raisim::World world;
    world.setTimeStep(simfreq_raisim);
    
    raisim::OgreVis *vis = raisim::OgreVis::get();

    /// these method must be called before initApp
    vis->setWorld(&world);
    vis->setWindowSize(1792, 1200); // Should be evenly divisible by 16!!
    vis->setImguiSetupCallback(imguiSetupCallback); // These 2 lines make the interactable gui visible
    vis->setImguiRenderCallback(imguiRenderCallBack);
    vis->setKeyboardCallback(raisimKeyboardCallback);
    vis->setSetUpCallback(setupCallback);
    vis->setAntiAliasing(2);

    /// starts visualizer thread
    vis->initApp();
    
    /// create raisim objects
    raisim::TerrainProperties terrainProperties;
    terrainProperties.frequency = 0.0;
    terrainProperties.zScale = 0.0;
    terrainProperties.xSize = 300.0;
    terrainProperties.ySize = 300.0;
    terrainProperties.xSamples = 50;
    terrainProperties.ySamples = 50;
    terrainProperties.fractalOctaves = 0;
    terrainProperties.fractalLacunarity = 0.0;
    terrainProperties.fractalGain = 0.0;

    raisim::HeightMap *ground = world.addHeightMap(0.0, 0.0, terrainProperties);
    vis->createGraphicalObject(ground, "terrain", "checkerboard_blue");
    world.setDefaultMaterial(0.8, 0.0, 0.0); //surface friction could be 0.8 or 1.0
    vis->addVisualObject("extForceArrow", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP); 
    
    
    // // WEIGHT VISUALIZATION FOR LCSS PAPER, KEEP FOR NOW.
    // // create raisim objects
    // double scale = 1;
    // vis->loadMeshFile("/home/kavehakbarihamed/raisim/workspace/A1_LL_Exp/rsc/Dumbell_5lb.STL", "weight1", false);
    // raisim::VisualObject *weightVis1 = vis->addVisualObject("weight_vis1", "weight1", "purple", {scale, scale, scale}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    // vis->loadMeshFile("/home/kavehakbarihamed/raisim/workspace/A1_LL_Exp/rsc/Dumbell_5lb.STL", "weight2", false);
    // raisim::VisualObject *weightVis2 = vis->addVisualObject("weight_vis2", "weight2", "purple", {scale, scale, scale}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);

    // raisim::Box *box = world.addBox(0.24,0.08,0.15,4.54);
    // box->setPosition(-0.005,0,0.25);
    // // vis->createGraphicalObject(box,"Payload","purple");

    // raisim::Box *box = world.addBox(0.04,0.04,0.01,20);
    // box->setPosition(0.,0,0.0);
    // // box->setPosition(-0.1,0,0.25);
    // vis->createGraphicalObject(box,"Payload","purple"); 

    auto& list = vis->getVisualObjectList();

    // ============================================================ //
    // ======================= SETUP Robot ======================== //
    // ============================================================ //
    std::vector<raisim::ArticulatedSystem*> A1;
    // A1.push_back(world.addArticulatedSystem(raisim::loadResource("Go1/Go1.urdf"))); // WHEN USING Go1, BE SURE TO CHANGE CMAKE TO USE CORRECT DYNAMICS
    A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_modified_new.urdf")));
    //A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_modified_up_mod_sensored.urdf")));
    
    
    vis->createGraphicalObject(A1.back(), "A1");

    //Rear offset -0.1
    //A1.back()->setGeneralizedCoordinate({0, 0, 0.5, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                  -0.7337, 1.0175, -2.035, 0.7337, 1.0175, -2.035, 0.0, 2.4532, -1.1582, 0.0, 2.4532, -1.1582});
    A1.back()->setGeneralizedCoordinate({0, 0, 0.12, 1, 0, 0, 0,
                                        0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6});   
    A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
    A1.back()->setName("A1_Robot");
    
    raisim::Box *box_right = world.addBox(200.0, 0.2, 0.8, 1000000, "rubber");//terrainProperties);
    raisim::Box *box_left = world.addBox(200.0, 0.2, 0.8, 1000000, "rubber");

    box_right->setPosition(0,-0.32,0.4);
    box_left->setPosition(0,0.32,0.4);

    //vis->createGraphicalObject(box_right, "right_wall", "checkerboard_blue");
    vis->createGraphicalObject(box_left, "left_wall", "checkerboard_blue");
    
    A1.back()->getCollisionBody("FR_foot/0").setMaterial("wood");
    A1.back()->getCollisionBody("FL_foot/0").setMaterial("wood");

    world.setMaterialPairProp("wood", "rubber", 0.8, 0, 0);


    LocoWrapper* loco_obj = new LocoWrapper(argc,argv);
    //loco_obj->setRFfalse();
    SRBNMPC* loco_plan = new SRBNMPC(argc,argv,1,0);
    //loco_plan->generator();
    std::string file_name = "take2_1";
    
    std::string prefix_code = std::filesystem::current_path().string() + "/";
    std::string prefix_lib = std::filesystem::current_path().string() + "/";

    // Create a new NLP solver instance from the compiled code
    std::string lib_name = prefix_lib + file_name + ".so";
    casadi::Dict opts = {{"ipopt.print_level", 1}, {"print_time", 0},{"ipopt.max_iter", 10},{"ipopt.acceptable_tol", 1e-2},{"ipopt.acceptable_obj_change_tol", 1e-2}};
    casadi::Function solver = casadi::nlpsol("solver", "ipopt", lib_name, opts);

    float a[3] = {1.00000000, -1.99555712, 0.99556697};
    float b[3] = {0.00000246, 0.00000492, 0.00000246};
    // float a[3] = {1.0, -1.14298050253990, 0.41280159809619};
    // float b[3] = {0.06745527388907, 0.13491054777814, 0.06745527388907};
    populate_filter_f(filt, a, b, 3, 2);


    // ============================================================ //
    // ================= VIEW AND RECORDING OPTIONS =============== //
    // ============================================================ //
    raisim::gui::showContacts = false;
    raisim::gui::showForces = false;
    raisim::gui::showCollision = false;
    raisim::gui::showBodies = true;

    std::string cameraview = "side";
    bool panX = true;                // Pan view with robot during walking (X direction)
    bool panY = false;                // Pan view with robot during walking (Y direction)
    bool record = true;             // Record?
    double startTime = 0*ctrlHz;    // Recording start time
    double simlength = 20000;//60000;//300*ctrlHz;   // Sim end time
    double fps = 30;            
    //std::string directory = "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/Oct10/";
    std::string directory = "../data25/Jan17/";
    // std::string filename = "Payload_Inplace";
    std::string filename = "climb_11";//"JacVCL_OWCL_rt55_3";
    // std::string filename = "inplace_sim";


    // ============================================================ //
    // ========================= VIEW SETUP ======================= //
    // ============================================================ //
    // NOTE: Pi is defined in /dynamics/dynamicsSupportFunctions.h
    // NOTE: This section still needs some work. 
    double Pi=3.14;
    if(cameraview == "iso"){
        // vis->getCameraMan()->getCamera()->setPosition(3, 3, 2);
        vis->getCameraMan()->getCamera()->setPosition(2, -1, 0.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(5*Pi/6-Pi/2));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else if(cameraview == "isoside"){
        vis->getCameraMan()->getCamera()->setPosition(1.1, -2, 0.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(4*Pi/6-Pi/2));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else if(cameraview == "side"){
        vis->getCameraMan()->getCamera()->setPosition(0, -2, 0.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(0));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else if(cameraview == "front"){
        vis->getCameraMan()->getCamera()->setPosition(2, 0, 0.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(Pi/2));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else if(cameraview == "top"){
        //vis->getCameraMan()->getCamera()->setPosition(1, -3, 2.5);
        //vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(1.0));
        vis->getCameraMan()->getCamera()->setPosition(1, 0, 2.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(Pi/2));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
    }else{
        vis->getCameraMan()->getCamera()->setPosition(1, -3, 2.5);
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(1.0));
    }
    unsigned long mask = 0;
    if(raisim::gui::showBodies) mask |= raisim::OgreVis::RAISIM_OBJECT_GROUP;
    if(raisim::gui::showCollision) mask |= raisim::OgreVis::RAISIM_COLLISION_BODY_GROUP;
    if(raisim::gui::showContacts) mask |= raisim::OgreVis::RAISIM_CONTACT_POINT_GROUP;
    if(raisim::gui::showForces) mask |= raisim::OgreVis::RAISIM_CONTACT_FORCE_GROUP;
    vis->setVisibilityMask(mask);

    //if(panX) raisim::gui::panViewX = panX;
    //if(panY) raisim::gui::panViewY = panY;
    
    // ============================================================ //
    // ========================== RUN SIM ========================= //
    // ============================================================ //
    const std::string name = directory+filename+"_"+cameraview+".mp4";
    vis->setDesiredFPS(fps);
    long simcounter = 0;
    static bool added = false;

    raisim::Contact contactInstance;

    while (!vis->getRoot()->endRenderingQueued() && simcounter <= simlength){

        size_t dist_start = 40000*ctrlHz;               // Start the disturbance (if any)
        size_t dist_stop  = dist_start+200;             // Stop the disturbance (if any)

        // FOR LCSS PAPER, KEEP FOR NOW.
        // if (simcounter>1500 & !added){
        //     box->setPosition(-0.1,0,0.30);
        //     added = true;
        // }
        // raisim::Vec<3> pos = box->getComPosition();
        // weightVis1->offset = {pos[0]-0.245/2,pos[1]-.078/2,pos[2]-0.078};
        // weightVis2->offset = {pos[0]-0.245/2,pos[1]-.078/2,pos[2]-0.015};
        // // std::cout<<pos[0]<<"\t"<<pos[1]<<"\t"<<pos[2]<<std::endl;

        controller(A1,loco_obj,loco_plan,solver,simcounter,contactInstance);
        world.integrate();        
        
        if (simcounter%15 == 0)
            vis->renderOneFrame();
        
        if (!vis->isRecording() & record & simcounter>=startTime)
            vis->startRecordingVideo(name);
        
        auto currentPos = vis->getCameraMan()->getCamera()->getPosition();
        //if (raisim::gui::panViewX){
            Eigen::VectorXd jointPosTotal(18 + 1);
            Eigen::VectorXd jointVelTotal(18);
            jointPosTotal.setZero();
            jointVelTotal.setZero();
            A1.back()->getState(jointPosTotal, jointVelTotal);
            if (cameraview=="front"){
                currentPos[0] = jointPosTotal(0)+2;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="side"){
                currentPos[0] = jointPosTotal(0);
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="iso"){
                currentPos[0] = jointPosTotal(0)+2;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="isoside"){
                currentPos[0] = jointPosTotal(0)+1.1;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            }
        //}
        /*if (raisim::gui::panViewY){
            Eigen::VectorXd jointPosTotal(18 + 1);
            Eigen::VectorXd jointVelTotal(18);
            jointPosTotal.setZero();
            jointVelTotal.setZero();
            A1.back()->getState(jointPosTotal, jointVelTotal);
            if (cameraview=="front"){
                currentPos[1] = jointPosTotal(1);
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="side"){
                currentPos[1] = jointPosTotal(1)-2;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="iso"){
                currentPos[1] = jointPosTotal(1)-1;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            } else if(cameraview=="isoside"){
                currentPos[1] = jointPosTotal(1)-2;
                vis->getCameraMan()->getCamera()->setPosition(currentPos);
            }
        }*/
        //std::cout << "simcounter" << "\t" << simcounter << std::endl;

        // if(abs(jointPosTotal(1))>0.04){
        //     std::cout << simcounter << "\t" << jointPosTotal(0) << "\t" << jointPosTotal(1) << "\t" << jointPosTotal(2) << "\t" << jointPosTotal(15) << "\t" << jointPosTotal(18) << "\t" << -1 << std::endl;
        //     break;
        // }else if(abs(jointPosTotal(2)-0.5)>0.05){
        //     std::cout << simcounter << "\t" << jointPosTotal(0) << "\t" << jointPosTotal(1) << "\t" << jointPosTotal(2) << "\t" << jointPosTotal(15) << "\t" << jointPosTotal(18) << "\t" << -2 << std::endl;
        //     break;
        // }else if(jointPosTotal(15) > 0 || jointPosTotal(18) > 0){
        //     std::cout << simcounter << "\t" << jointPosTotal(0) << "\t" << jointPosTotal(1) << "\t" << jointPosTotal(2) << "\t" << jointPosTotal(15) << "\t" << jointPosTotal(18) << "\t" << -3 << std::endl;
        //     break;
        // }else if(jointPosTotal(0)>7){
        //     std::cout << simcounter << "\t" << jointPosTotal(0) << "\t" << jointPosTotal(1) << "\t" << jointPosTotal(2) << "\t" << jointPosTotal(15) << "\t" << jointPosTotal(18) << "\t" << 1 << std::endl;
        //     break;
        // }
        
        simcounter++;
        
    }

    // End recording if still recording
    if (vis->isRecording())
        vis->stopRecordingVideoAndSave();

    /// terminate the app
    vis->closeApp();

    delete loco_obj;
    clear_filter_f(filt);

    return 0;
}

