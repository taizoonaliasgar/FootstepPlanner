//
// Authror: Randy Fawcett on 12/2021.
//
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "raisim/OgreVis.hpp"
//#include "randyImguiPanel.hpp"
#include "raisimBasicImguiPanel.hpp"
#include "raisimKeyboardCallback.hpp" // THIS IS WHERE cmd.vel AND cmd.pose ARE IMPLEMENTED (as globals)
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
#include "LocoWrapperwalk.hpp"
#include "SRBNMPC.hpp"
using namespace std;
// namespace fs = std::filesystem;

FiltStruct_f* filt = (FiltStruct_f*)malloc(sizeof(FiltStruct_f));

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

void controller(std::vector<raisim::ArticulatedSystem *> A1, LocoWrapperwalk *loco_obj, SRBNMPC* loco_plan, casadi::Function solver, size_t controlTick, raisim::Contact contactInstance) {
    /////////////////////////////////////////////////////////////////////
    //////////////////////////// INITIALIZE
    /////////////////////////////////////////////////////////////////////
    static size_t loco_kind = TROT;                        // Gait pattern to use
    size_t pose_kind = POSE_CMD;                   // Pose type to use (if loco_kind is set to POSE)
    size_t settling = 0*0.2*ctrlHz;                   // Settling down
    size_t duration = 0*0.8*ctrlHz;                   // Stand up 
    size_t loco_start = settling + duration;        // Start the locomotion pattern

    
    double *tau;
    //double tau[18] = {0};
    double jpos[18], jvel[18];
    
    Eigen::VectorXd jointTorqueFF = Eigen::MatrixXd::Zero(TOTAL_DOF,1);
    Eigen::VectorXd jointPosTotal = Eigen::MatrixXd::Zero(TOTAL_DOF+1,1);
    Eigen::VectorXd jointVelTotal = Eigen::MatrixXd::Zero(TOTAL_DOF,1);    
        
    raisim::Mat<3,3> rotMat;
    Eigen::Matrix<double, 3, 1> eul;
    Eigen::Matrix<double, 4, 1> quat;

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
    
    raisim::MatDyn D = A1.back()->getMassMatrix();
    raisim::VecDyn H = A1.back()->getNonlinearities({0,0,-9.81});
    
    //raisim::Vec<3>& impulse = contactInstance.getImpulse();

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
    //eul(1) = eul(1)-1.57079632679;
    //std::cout << "quat" << "\t" << quat << std::endl;
    //std::cout << "eul" << "\t" << eul(0) << "\t" << eul(1) << "\t" << eul(2) << std::endl; 
    for(size_t i=0; i<3; ++i){
        jpos[i] = jointPosTotal(i);
        jvel[i] = jointVelTotal(i);
        jpos[i+3] = eul(i);
        jvel[i+3] = jointVelTotal(i+3);
    }

    for(size_t i=6; i<18; ++i){
        /*if((i-1)%3==0){
            jpos[i] = jointPosTotal(i+1);//-1.57079632679;
        }else{
            jpos[i] = jointPosTotal(i+1);//-1.57079632679;
        }*/
        jpos[i] = jointPosTotal(i+1);
        jvel[i] = jointVelTotal(i);
    }
    

    Eigen::Matrix<double,16,1> q0;
    q0.setZero(16,1);
    q0.block(0,0,3,1) << jpos[0],jpos[1],jpos[2];//= jointPosTotal.block(0,0,3,1);
    q0.block(3,0,3,1) << jvel[0],jvel[1],jvel[2];//= jointVelTotal.block(0,0,3,1);
    q0.block(6,0,3,1) = eul;
    q0.block(9,0,3,1) << jvel[3],jvel[4],jvel[5];//.block(3,0,3,1);
    std::map<std::string, casadi::DM> arg, res;

    //std::cout << q0 << std::endl;

    int force[4] = {0};
    for(auto &con: A1.back()->getContacts()){
        int conInd = con.getlocalBodyIndex();
        force[conInd/3-1] = 500;
        force[conInd/3-1] = con.getNormal().e().norm();
    }

    float filt_vel_temp[3] = {jvel[0],jvel[1],jvel[2]};
    discrete_butter_f(filt,filt_vel_temp);
    //printf("Fwd Vel: %0.3f || Lat Vel: %0.3f || Yaw Vel: %0.3f || Pitch: %0.3f\n",filt_vel_temp[0],filt_vel_temp[1],filt_vel_temp[2],jpos[4]);//cmd.pose[1]);


    /////////////////////////////////////////////////////////////////////
    //////////////////////////// CONTROL
    /////////////////////////////////////////////////////////////////////
    // Update the desired torques
    if(controlTick < settling){ // Settle down
        // loco_obj->posSetup(jointPosTotal.head(3));
        //double temp[18] = {0};
        //tau = temp;
        //loco_obj->initStandVars(jointPosTotal.block(0,0,3,1),jointPosTotal(5),(int)duration);
    }
    else if(controlTick >= settling & controlTick < loco_start){ // Start standing
        //loco_obj->calcTau(jpos,jvel,rotMatrixDouble,force,STAND,controlTick);
        //tau = loco_obj->getTorque();

    }
    else if(controlTick >= loco_start){ // Start locomotion
        
        //if(controlTick == loco_start){
        
        loco_obj->updatestate(jpos,jvel,rotMatrixDouble);
        int duration_data =0;
        //}
        foot_position = loco_obj->getfootposition();
        hip_position = loco_obj->gethipposition();
        int stancephase = loco_obj->stancecounter();
        //std::cout << foot_position.block(0,1,3,1).transpose() << std::endl;
        //std::cout << foot_position.block(0,0,3,1).transpose() << std::endl;
        if(controlTick%10==0){

            int controlMPC = std::floor(controlTick/10); 
            //std::cout << "controlMPC:" << controlMPC << std::endl;
            casadi::DM X_prev = loco_plan->getprevioussol_ll(q0,foot_position,controlMPC);//casadi::DM::zeros(NFS*(HORIZ+1)+NFI*HORIZ,1); 
            if(controlMPC<1){
                q0.block(12,0,4,1) << foot_position(0,0),foot_position(0,1),foot_position(0,2),foot_position(0,3);//0.15,0.15,-0.1,-0.1;
            }else{
                q0(12) = double(X_prev(12));
                q0(13) = double(X_prev(13));
                q0(14) = double(X_prev(14));
                q0(15) = double(X_prev(15));
            }
            
            casadi::DM p = loco_plan->motionPlannerN(q0,controlMPC);
            loco_plan->setpreviousp(p);

            arg["lbx"] = loco_plan->lowerboundx(p, controlMPC);
            arg["ubx"] =  loco_plan->upperboundx(p);
            arg["lbg"] =  loco_plan->lowerboundg();
            arg["ubg"] =  loco_plan->upperboundg();
            arg["x0"] = X_prev;
            arg["p"] = p;
            
            auto start = std::chrono::high_resolution_clock::now();
            res = solver(arg);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            duration_data = static_cast<int>(duration.count());
            //std::cout << "NMPC Solve Time: " << duration_data << "ms" << std::endl;

            loco_plan->setprevioussol(res.at("x"));
            
            opt_HLMPC_state = loco_plan->getNMPCsol2(controlMPC);
            
            loco_obj->setoptNLstate(opt_HLMPC_state);
            loco_obj->setcontactconfig(controlMPC);
            loco_plan->mpcdataLog(q0, opt_HLMPC_state.block(16,0,12,1), controlMPC, Eigen::Matrix<double, 12, 1>::Zero());
            

            //float pose_temp[3] = {opt_HLMPC_state(0),opt_HLMPC_state(1),opt_HLMPC_state(2)};
            //float vel_temp[3] = {opt_HLMPC_state(3),opt_HLMPC_state(4),opt_HLMPC_state(15)};
            //loco_obj->updateVel(vel_temp);
            //loco_obj->updatePose(pose_temp);
        }
        //float vel_temp[3] = {0.5,0,0};
        //loco_obj->updateVel(vel_temp);
        
        loco_obj->setRaisimD(Dr);
        loco_obj->setRaisimH(Hr);
        loco_obj->calcTau(jpos,jvel,rotMatrixDouble,force,UPWALK,controlTick,duration_data);
        
        tau = loco_obj->getTorque();
        
        //tau = tau*0;//Eigen::MatrixXd::Zero(18,1);// >(tau,18);

        //double temp[18] = {0};
        //tau = temp;
        //std::cout<< "Torque deploy" << std::endl;
        //for(int i=0;i<18;i++){
        //    tau[i] = 0;
        //}
    }

    //std::cout<< "Torque deploy 1.3" << std::endl;
    jointTorqueFF = Eigen::Map< Eigen::Matrix<double,18,1> >(tau,18);
    //for(int i=0;i<18;i++){
    //        jointTorqueFF(i) = tau[i];
    //}
    //std::cout<< "Torque deploy 1.5" << std::endl;

    jointTorqueFF.block(0,0,6,1).setZero();
    
    // Set the desired torques
    A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
    A1.back()->setGeneralizedForce(jointTorqueFF);

    // Example force in the x-direction

    // Apply the disturbance force to the desired body
    // Assuming the body index is 0, adjust as needed
    if(controlTick > 8000 && controlTick < 9000){

        double amp = 120.0;
        double dist = distx(controlTick, amp);
        Eigen::Vector3d disturbanceForce(dist, 0.0, 0.0); 
        A1.back()->setExternalForce(0, disturbanceForce);
    }
};



int main(int argc, char *argv[]) {
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
    //A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_modified_new.urdf")));
    //A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_modified_up_mod.urdf")));
    A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_modified_up_mod_2.urdf")));
    vis->createGraphicalObject(A1.back(), "A1");

    //A1.back()->setGeneralizedCoordinate({0, 0, 0.35, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                    -0.5804, 1.8973, -2.3609, 0.5804, 1.8973, -2.3609, 0.0, 3.0156, -2.2091, 0.0, 3.0156, -2.2091});

    //A1.back()->setGeneralizedCoordinate({0, 0, 0.45, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                    -0.4944, 0.8218, -1.6435, 0.4944, 0.8218, -1.6435, 0.0, 2.6951, -1.5379, 0.0, 2.6951, -1.5379});
    
    //Full Order Simulation
    //Rear offset -0.1
    A1.back()->setGeneralizedCoordinate({0, 0, 0.5, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
                                       -0.7337, 1.0175, -2.035, 0.7337, 1.0175, -2.035, 0.0, 2.4532, -1.1582, 0.0, 2.4532, -1.1582});
    //Rear offset -0.01
    //A1.back()->setGeneralizedCoordinate({0, 0, 0.5, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                    -0.7337, 1.0175, -2.035, 0.7337, 1.0175, -2.035, 0.0, 2.247, -1.2899, 0.0, 2.247, -1.2899});
    //Rear offset -0.04
    //A1.back()->setGeneralizedCoordinate({0, 0, 0.5, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                    -0.7337, 1.0175, -2.035, 0.7337, 1.0175, -2.035, 0.0, 2.3305, -1.2703, 0.0, 2.3305, -1.2703});
    //Rear offset -0.08
    //A1.back()->setGeneralizedCoordinate({0, 0, 0.5, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                    -0.7337, 1.0175, -2.035, 0.7337, 1.0175, -2.035, 0.0, 2.4195, -1.2068, 0.0, 2.4195, -1.2068});

    //Rear offset -0.06
    //A1.back()->setGeneralizedCoordinate({0, 0, 0.5, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                    -0.7337, 1.0175, -2.035, 0.7337, 1.0175, -2.035, 0.0, 2.3784, -1.244, 0.0, 2.3784, -1.244});
    //Rear offset -0.05
    //A1.back()->setGeneralizedCoordinate({0, 0, 0.5, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                     -0.7337, 1.0175, -2.035, 0.7337, 1.0175, -2.035, 0.0, 2.3553, -1.2585, 0.0, 2.3553, -1.2585});
    //Front offset 0.15
    //A1.back()->setGeneralizedCoordinate({0, 0, 0.5, 1,0,0,0,//0.9238795,0,0.3826834,0,//1, 0, 0, 0,
    //                                     -0.596, 0.9333, -1.8665, 0.596, 0.9333, -1.8665, 0.0, 2.4532, -1.1582, 0.0, 2.4532, -1.1582});
    
    A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
    A1.back()->setName("A1_Robot");
    //assignBodyColors(vis,"A1","collision","gray");

    raisim::Box *box_right = world.addBox(200.0, 0.2, 0.8, 1000000, "rubber");//terrainProperties);
    raisim::Box *box_left = world.addBox(200.0, 0.2, 0.8, 1000000, "rubber");

    box_right->setPosition(0,-0.37,0.4);
    box_left->setPosition(0,0.37,0.4);

    //vis->createGraphicalObject(box_right, "right_wall", "checkerboard_blue");
    vis->createGraphicalObject(box_left, "left_wall", "checkerboard_blue");
    
    A1.back()->getCollisionBody("FR_foot/0").setMaterial("wood");
    A1.back()->getCollisionBody("FL_foot/0").setMaterial("wood");

    world.setMaterialPairProp("wood", "rubber", 0.8, 0, 0);

    //auto& MMP = world.getMaterialPairProp(A1.back()->getCollisionBody("FR_foot/0").getMaterial(),
    //                                    ground->getCollisionObject().getMaterial());//box_right->getCollisionObject().getMaterial()"");//raisim::MaterialPairProperties
    //const Ogre::MaterialPtr& footmaterial = A1.back()->getCollisionBody("FR_FOOT").getMaterial();

    auto block1 = world.addBox(0.1, 0.46, 0.05, 1.0); // size (1x1x1) and mass (1.0)
    auto block2 = world.addBox(0.1, 0.46, 0.05, 1.0); // size (2x0.5x0.5) and mass (2.0)

    // Set the position of the blocks
    block1->setPosition(1.0, 0.0, 0.025); // Position (x=0, y=0, z=0.5)
    block2->setPosition(1.2, 0.0, 0.025); // Position (x=2, y=2, z=0.25)
    
    vis->createGraphicalObject(block1, "block1", "checker_yellow");
    vis->createGraphicalObject(block2, "block2", "checker_yellow");

    LocoWrapperwalk* loco_obj = new LocoWrapperwalk(argc,argv);
    
    SRBNMPC* loco_plan = new SRBNMPC(argc,argv,1,0);
    //loco_plan->generator();
    std::string file_name = "upright_h5_71";
    // code predix
    // std::string prefix_code = "/home/trec/WorkRaj/raisim_legged/FootstepPlanner/build/";//fs::current_path().string() + "/";
    std::string prefix_code = std::filesystem::current_path().string() + "/";
    // shared library prefix
    // std::string prefix_lib = "/home/trec/WorkRaj/raisim_legged/FootstepPlanner/build/";//fs::current_path().string() + "/";
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
    bool record = false;             // Record?
    double startTime = 0*ctrlHz;    // Recording start time
    double simlength = 60000;//60000;//300*ctrlHz;   // Sim end time
    double fps = 30;            
    //std::string directory = "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/datalog/Oct10/";
    std::string directory = "../datalog/Oct14/";
    // std::string filename = "Payload_Inplace";
    std::string filename = "upright_A1";
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
        std::cout << "simcounter" << "\t" << simcounter << std::endl;
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
