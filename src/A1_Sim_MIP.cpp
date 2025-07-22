//
// Authror: Basit M. Imran.
// Date : 2024-04-17
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "unitree_legged_sdk/unitree_legged_sdk.h"
#include "unitree_legged_sdk/unitree_joystick.h"
#include "raisim/OgreVis.hpp"
//#include "randyImguiPanel.hpp"
#include "raisimBasicImguiPanel.hpp"
#include "raisimKeyboardCallback.hpp"
#include "raisim/RaisimServer.hpp"

#include "helper.hpp"
#include "Filters.h"

#include "timer.h"

#include <fstream>
#include <iostream>
#include <filesystem>
//HDSRL header
#include "LocoWrapper.hpp"
#include "SRBNMPC.hpp"
#include "A1_Dynamics_full.h"

#include "shared_structs_ex2.hpp"
#include "Transforms.hpp"
//#include "OtherUtils.hpp"
//#include <yaml-cpp/yaml.h>
#include "stdio.h"

#include "mip/mip_all.hpp"
#include "microstrain/connections/serial/serial_connection.hpp"
#include "example_utils.hpp"
#include "fusedVelocityEstimate.hpp"

#include <chrono>
#include <thread>
#include <memory>

using namespace UNITREE_LEGGED_SDK;

sharedData HLData;
sharedData LLData;
sharedData SimData;
sharedData IMUData;

class ExternalComm
{
private:
    std::shared_ptr<microstrain::connections::SerialConnection> connection;
    const size_t parseBufferSize = 1024;
    std::unique_ptr<uint8_t[]> parseBuffer = std::make_unique<uint8_t[]>(parseBufferSize);

    // std::unique_ptr<uint8_t[]> parseBuffer = new uint8_t[parseBufferSize];
    

    // Device pointer
    std::unique_ptr<mip::Interface> device;

    float sensor_to_vehicle_rotation_euler[3] = {0.0, 0.0, 0.0};

    // Data stores
    mip::data_sensor::GpsTimestamp sensor_gps_time;
    mip::data_sensor::ScaledAccel sensor_accel;
    mip::data_sensor::ScaledGyro sensor_gyro;
    mip::data_sensor::CompEulerAngles sensor_comp_euler_angles;

    mip::data_filter::Timestamp filter_gps_time;
    mip::data_filter::Status filter_status;
    mip::data_filter::EulerAngles filter_euler_angles;
    mip::data_filter::CompAngularRate filter_comp_angular_rate;
    mip::data_filter::CompAccel filter_comp_accel;

    // Dispatch handlers
    mip::DispatchHandler sensor_data_handlers[4];
    mip::DispatchHandler filter_data_handlers[5];

    // State tracking
    bool filter_state_running = false;
    bool is_initialized = false;
    std::mutex update_mutex;

    bool SIMFlag = true;
    std::ofstream csvFile;//("contact_forces.csv");
    VelocityKalman3D velocity_filter;


public:
    ExternalComm(int argc, char* argv[]):velocity_filter(0.001, 1e-2, 2.7e-4)
    // : parseBufferSize(1024), 
        // parseBuffer(new uint8_t[parseBufferSize])
    {
		
        // const std::string port = "/dev/ttyACM0";//argv[2];
        // const uint32_t baudrate = 9600;//std::stoi(argv[3]);
        
        // std::cout << "Initializing IMU Sensor on port " << port << " at " << baudrate << " baud" << std::endl;
        
        // // Create serial connection
        // connection = std::make_shared<microstrain::connections::SerialConnection>(port, baudrate);
        
        // // Try to connect
        // if (!connection->connect()) {
        //     throw std::runtime_error("Failed to connect to " + port + " at " + std::to_string(baudrate) + " baud");
        // }
        
        // // Create device interface
        // device = std::make_unique<mip::Interface>(
        //     connection.get(),       // Connection pointer
        //     parseBuffer.get(),      // Parse buffer
        //     parseBufferSize,        // Parse buffer size
        //     1000,                   // Parse timeout (ms)
        //     1000                    // Base reply timeout (ms)
        // );
        
        // std::cout << "Connected to device on " << port << " at " << baudrate << " baud" << std::endl;

        double ad[3] = {1.0, -1.47548044359265, 0.58691950806119};
        double bd[3] = {0.02785976611714, 0.05571953223427, 0.02785976611714};
        populate_filter_d(jointfilter,ad,bd,3,12);

        // 0.75 Hz
        float af[3] = {1.00000000,-1.99333570,0.99335783};
        float bf[3] = {0.00000553,0.00001107,0.00000553};
        populate_filter_f(remotefilter, af, bf, 3, 2);

        // 2 Hz 
        float aa[3] = {1.00000000,-1.98222893,0.98238545};
        float ba[3] = {0.00003913,0.00007826,0.00003913};
        populate_filter_f(angfilter, aa, ba, 3, 2);

        csvFile.open("../data25/kneepose.csv");

        // StandDuration = 10000;
        // SettlingTime = 8000;

        /********************************************RASIM INIT********************************************** */
        // const int NUMBER_OF_SIMS = 1;
        // const float threshold = 0.4;
        // bool shared_data_backed_up = 0;
        //raisim::OgreVis *vis = raisim::OgreVis::get();		
        // std::cout << "Setup filter" << std::endl;	
        
    }	

	virtual ~ExternalComm(){
		clear_filter_d(jointfilter);
		clear_filter_f(remotefilter);
		clear_filter_f(angfilter);

		
        delete vis;
        delete ground;
        delete list;
        auto test = A1.back();
        A1.pop_back();
        delete test;
        csvFile.close();
    }




    //support functions
	void setupCallback();
	// void plotGRFs(std::map<std::string, raisim::VisualObject>* list, const std::vector<double>& GRF, const std::vector<double>& feet_vec, const std::vector<double>& contacts);
    void plotGRFs(std::map<std::string, raisim::VisualObject>* list, Eigen::Matrix<double,17,1> GRF, Eigen::Matrix<double, 3, 4> toePos, const int contacts[4]);
    void setupRaisim();

    // main thread execution functions
	void Calc();
	void HighLevel();
	void SimExec();//(std::ofstream &file_est);  

    void connectIMU();//(mip::Interface& device);
    void setupIMUfilter();//(mip::Interface& device);
    void configureIMU();//(mip::Interface& device);
    void getIMMUdata();//(mip::Interface& device);

    std::unique_ptr<LocoWrapper> loco_obj;
	std::unique_ptr<SRBNMPC> nmpc_obj;
	
	FiltStruct_d* jointfilter  = (FilterStructure_d*)malloc(sizeof(FilterStructure_d));
	FiltStruct_f* angfilter    = (FilterStructure_f*)malloc(sizeof(FilterStructure_f));
	FiltStruct_f* remotefilter = (FilterStructure_f*)malloc(sizeof(FilterStructure_f));

    float LLdt = 0.00100001f;
    float HLdt = 0.0100001f;
    long simcounter = 0;
    size_t settling = 0.2*ctrlHz;                   // Settling down
    size_t duration = 1.8*ctrlHz;                   // Stand up 
    size_t loco_start = settling + duration;        // Start the locomotion pattern
    double switchtime = 20;

    //Eigen::VectorXd jointTorqueFF = Eigen::MatrixXd::Zero(TOTAL_DOF,1);
    Eigen::VectorXd jointPosTotal = Eigen::MatrixXd::Zero(TOTAL_DOF+1,1); // +1 is for 4th Component of Quaternion 
    Eigen::VectorXd jointVelTotal = Eigen::MatrixXd::Zero(TOTAL_DOF,1);
    raisim::Mat<3,3> rotMat;

    //Raisim stuff
    bool setup_raisim = true;
    raisim::World world;
    raisim::OgreVis *vis = raisim::OgreVis::get();
    raisim::HeightMap *ground;
    std::map<std::string, raisim::VisualObject>* list = nullptr;
    std::string cameraview = "side";
    bool panX = true;                // Pan view with robot during walking (X direction)
    bool panY = false;                // Pan view with robot during walking (Y direction)
    bool record = false;            // Record?
    double fps = 30;            
    std::string directory = "../data25/May26/";
    std::string filename = "MTSim";
    std::string name = directory+filename+"_"+".mp4";
    
    double startTime = 0*ctrlHz;    // Recording start time
    double simlength = 50*ctrlHz;

    //Estimator
    int rearweight_est = 4;
    double yzdot_thresh = 0.3;
    double xdot_thresh = 0.3;
    double yzdot_thresh2 = 0.8;
    double xdot_thresh2 = 0.5;

    //Estimator
    void getthetadot(double q[18],double dq[18]);
    void kinestimatorrr(double q[18], double dq[18], int contact[4], Eigen::Matrix<double,3,3> R);
    void getStateEstimatefullll(double q[18], double dq[18], int contact[4], Eigen::Matrix<double,3,3> R, Eigen::Matrix<double,3,4> toes, int robotdown, size_t ctrlTick, Eigen::Vector3d acc_bFrame);
    //A1
    std::vector<raisim::ArticulatedSystem*> A1;
    
    float rollIMU = 0.0f, pitchIMU = 0.0f, yawIMU = 0.0f;
    float gyroXIMU = 0.0f, gyroYIMU = 0.0f, gyroZIMU = 0.0f;
    uint64_t timestampIMU = 0;

    struct ForceData {
        double normalForce;      // Normal force magnitude
        double frictionForce[2]; // Tangential friction force [x,y]
        double position[3];      // Contact position in world frame
        bool inContact;          // Contact state
    };

    ForceData footForces[4];
    void getDetailedContactForces();
    size_t calfIdx[4]= {0,1,2,3};
    raisim::Vec<3> calfPose = {0.0, 0.0, 0.0};
    raisim::Vec<3> localPosition = {0.0, 0.0, -0.2};
    // raisim::ArticulatedSystem::JointRef jointFR(static_cast<size_t>(0),A1.back());
    // Eigen::Matrix<double,3,3> Rpremult = Eigen::MatrixXd::Zero(3,3);
    // Eigen::Matrix<double,3,3> Rpostmult = Eigen::MatrixXd::Identity(3,3);
    // Rpremult(0,2)=-1;Rpremult(1,1)=-1;Rpremult(2,0)=-1;
    // Rpostmult(0,0)=1;Rpostmult(1,1)=-1;Rpostmult(2,2)=-1;
};


void ExternalComm::setupRaisim(){  
	
    raisim::World::setActivationKey(raisim::loadResource("activation.raisim"));
    world.setTimeStep(simfreq_raisim);

    /// these method must be called before initApp
    vis->setWorld(&world);
    vis->setWindowSize(1920, 1080); // Should be evenly divisible by 16!!
    vis->setImguiSetupCallback(imguiSetupCallback); // These 2 lines make the interactable gui visible
    vis->setImguiRenderCallback(imguiRenderCallBack);
    vis->setKeyboardCallback(raisimKeyboardCallback);
    vis->setSetUpCallback(std::bind(&ExternalComm::setupCallback, this));
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

    //raisim::HeightMap 
    ground = world.addHeightMap(0.0, 0.0, terrainProperties);
    vis->createGraphicalObject(ground, "terrain", "checkerboard_blue");
    world.setDefaultMaterial(0.8, 0.0, 0.0); //surface friction could be 0.8 or 1.0
    vis->addVisualObject("extForceArrow", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);

    // Foot force visualization arrows
    vis->addVisualObject("GRF1", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF2", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF3", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF4", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);

    auto& list = vis->getVisualObjectList();
    // std::vector<raisim::ArticulatedSystem*> A1;
    A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_modified_new.urdf")));   // with NMPC and LL 
    vis->createGraphicalObject(A1.back(), "A1");
    A1.back()->setName("A1_Robot");
    A1.back()->setGeneralizedCoordinate({0, 0, 0.12, 1 , 0, 0, 0,-0.2, Pi/3, -2.6, 0.2, Pi/3, -2.6, -0, Pi/3, -2.6, 0, Pi/3, -2.6});
    A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);

    raisim::Box *box_right = world.addBox(200.0, 0.2, 0.8, 1000000, "rubber");//terrainProperties);
    raisim::Box *box_left = world.addBox(200.0, 0.2, 0.8, 1000000, "rubber");

    box_right->setPosition(0,-0.32,0.4);
    box_left->setPosition(0,0.32,0.4);

    // vis->createGraphicalObject(box_right, "right_wall", "checkerboard_blue");
    vis->createGraphicalObject(box_left, "left_wall", "checkerboard_blue");
    
    A1.back()->getCollisionBody("FR_foot/0").setMaterial("wood");
    A1.back()->getCollisionBody("FL_foot/0").setMaterial("wood");

    world.setMaterialPairProp("wood", "rubber", 0.7, 0, 0);

    raisim::gui::showContacts = false;
    raisim::gui::showForces = false;
    raisim::gui::showCollision = false;
    raisim::gui::showBodies = true;        

    // ============================================================ //
    // ========================= VIEW SETUP ======================= //
    // ============================================================ //
    if(cameraview == "iso"){
        vis->getCameraMan()->getCamera()->setPosition(-1, -3, 0.5);
        vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(4.5*Pi/6-Pi/2));
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
        vis->getCameraMan()->getCamera()->setPosition(2, -1, 6);
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(0));
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
    vis->setDesiredFPS(fps);
    static bool added = false; 

    calfIdx[0] = A1.back()->getBodyIdx("FR_thigh");
    calfIdx[1] = A1.back()->getBodyIdx("FL_thigh");
    calfIdx[2] = A1.back()->getBodyIdx("RR_thigh");
    calfIdx[3] = A1.back()->getBodyIdx("RL_thigh");
    // jointFR = A1.back()->getJoint("FR_calf_joint");
}

void ExternalComm::setupCallback() {

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

void ExternalComm::getDetailedContactForces() {
    // Reset all force data
    for (int i = 0; i < 4; i++) {
        footForces[i].normalForce = 0.0;
        footForces[i].frictionForce[0] = 0.0;
        footForces[i].frictionForce[1] = 0.0;
        footForces[i].position[0] = 0.0;
        footForces[i].position[1] = 0.0;
        footForces[i].position[2] = 0.0;
        footForces[i].inContact = false;
    }
    
    // Get all contacts from Raisim
    for (auto &con: A1.back()->getContacts()) {
        // Get foot index (0-3 for FR, FL, RR, RL)
        int conInd = con.getlocalBodyIndex();
        int footIndex = conInd/3 - 1;
        
        // Ensure valid foot index
        if (footIndex >= 0 && footIndex < 4) {
            // Get contact position
            const auto& position = con.getPosition();
            footForces[footIndex].position[0] = position[0];
            footForces[footIndex].position[1] = position[1];
            footForces[footIndex].position[2] = position[2];
            
            // Get normal force (Z component)
            const auto& normal = con.getNormal();
            footForces[footIndex].normalForce = normal.e().norm();
            
            // Get friction forces (X,Y components)
            const auto& force = con.getImpulse();
            // Project force to get tangential components
            Eigen::Vector3d normalVec(normal.e()(0), normal.e()(1), normal.e()(2));
            normalVec.normalize();
            Eigen::Vector3d forceVec(force[0], force[1], force[2]);
            
            // Calculate friction force by subtracting normal component
            double normalComponent = forceVec.dot(normalVec);
            Eigen::Vector3d frictionVec = forceVec - normalComponent * normalVec;
            
            footForces[footIndex].frictionForce[0] = frictionVec(0);
            footForces[footIndex].frictionForce[1] = frictionVec(1);
            
            // Mark as in contact
            footForces[footIndex].inContact = true;
            
            // Print detailed force information
            printf("Foot %d: Normal=%.2f N, Friction=[%.2f, %.2f] N, Pos=[%.3f, %.3f, %.3f]\n", 
                   footIndex, 
                   footForces[footIndex].normalForce,
                   footForces[footIndex].frictionForce[0],
                   footForces[footIndex].frictionForce[1],
                   footForces[footIndex].position[0],
                   footForces[footIndex].position[1],
                   footForces[footIndex].position[2]);
        }
    }
}



//void ExternalComm::plotGRFs(std::map<std::string, raisim::VisualObject>* list, const std::vector<double>& GRF, const std::vector<double>& feet_vec, const std::vector<double>& contacts) {
void ExternalComm::plotGRFs(std::map<std::string, raisim::VisualObject>* list, Eigen::Matrix<double,17,1> GRF, Eigen::Matrix<double, 3, 4> toePos, const int contacts[4]) {

    // Ensure the vectors are of the correct size
    // if (GRF.size() < 12 || feet_vec.size() < 12 || contacts.size() < 4) {
    //     std::cerr << "Error: Input vectors are of incorrect size." << std::endl;
    //     return;
    // }

    for (int i = 0; i < 4; ++i) {
        raisim::Vec<3> dir;
        for (int j = 0; j < 3; ++j) {
            dir[j] = GRF(3 * i + j);
        }

        // Normalize the direction vector if it is not a zero vector
        double norm = dir.norm();
        if (norm > 1e-6) {  // Check if the vector is non-zero to avoid division by zero
            dir /= norm;
        } else {
            dir.setZero();
        }

        // Convert direction vector to a rotation matrix that aligns the z-axis with the direction vector
        raisim::Mat<3, 3> rot;
        if (contacts[i] == 1 && norm > 1e-6) {
            raisim::zaxisToRotMat(dir, rot);
        } else {
            rot.setIdentity();  // Set rotation to identity if no contact or zero norm

        }

        // Visual object key
        std::string objKey = "GRF" + std::to_string(i + 1);
        (*list)[objKey].offset = {toePos(0,i), toePos(1,i),toePos(2,i)};
        (*list)[objKey].scale = {0.2, 0.2, 0.005 * norm};  // Scaling based on the norm of the GRF vector

        // Set the rotation matrix
        (*list)[objKey].rotationOffset = rot;

    }
}

void ExternalComm::HighLevel(){

    //std::cout << "Inhighlevel" << std::endl;
    
    updateData(GET_DATA, HL_DATA, &HLData);
    if(HLData.control_Tick > switchtime*1000+2999){//} && HLData.control_Tick%10==0){ // Settle down
        auto start = std::chrono::high_resolution_clock::now();
        nmpc_obj->planner_MT(HLData.control_Tick, HLData.q, HLData.dq, HLData.toePos, HLData.QPforce);
        HLData.comDes= nmpc_obj->returncomDes();
        HLData.fDes= nmpc_obj->returnfDes();
        HLData.solvetime = nmpc_obj->returnSolveTime();
        int* indcon = nmpc_obj->returnConInd30(HLData.control_Tick);
        HLData.ind[0] = indcon[0];
        HLData.ind[1] = indcon[1];
        HLData.ind[2] = indcon[2];
        HLData.ind[3] = indcon[3];
        HLData.ind[4] = indcon[4];

        updateData(SET_DATA, HL_DATA, &HLData);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        // std::cout << duration.count() << "\t" << "Full high level time" << std::endl;
    }
    //std::cout << "Exitinghighlevel" << std::endl;
    
    
        
}

void ExternalComm::Calc(){

    // std::cout << "Inlowlevel" << std::endl;

    updateData(GET_DATA, LL_DATA, &LLData);

    if(LLData.control_Tick < settling){ // Settle down
        //double temp[18] = {0};
        //tau = temp;
        loco_obj->initStandVars(jointPosTotal.block(0,0,3,1),jointPosTotal(5),(int)duration);
        
    }else if(LLData.control_Tick >= settling & LLData.control_Tick < loco_start){ // Start standing
        
        loco_obj->calcTau2(LLData.q,LLData.dq,LLData.rotMatrixDouble,STAND,LLData.control_Tick,LLData.solvetime);  
    }else{

        loco_obj->ExpWrapper(LLData.q,LLData.dq,LLData.rotMatrixDouble,LLData.control_Tick,LLData.solvetime,LLData.ind,LLData.comDes,LLData.fDes);
    }

    LLData.tau = Eigen::Map<Eigen::VectorXd>(loco_obj->getTorque(),18);
    LLData.tau.block(0,0,6,1).setZero();
    LLData.toePos = loco_obj->getfootposition();
    LLData.toe_prev = loco_obj->gettoe_prev();
    LLData.QPforce = loco_obj->getpreviousQPforce();
    const int* ind_LL = loco_obj->getConDes();
    LLData.ind_LL[0] = ind_LL[0];
    LLData.ind_LL[1] = ind_LL[1];
    LLData.ind_LL[2] = ind_LL[2];
    LLData.ind_LL[3] = ind_LL[3];

    updateData(SET_DATA, LL_DATA, &LLData);
    loco_obj->settoe_prev();
    // std::cout << "Exitinglowlevel" << std::endl;
    
}



void ExternalComm::kinestimatorrr(double q[18], double dq[18], int contact[4], Eigen::Matrix<double,3,3> R){
    
    float numContact = contact[0]+contact[1]+contact[2]+contact[3];
	// ================================== //
	// ========= Kin Estimator ========== //
	// ================================== //

	// toe pos
	double fr_toe[3], fl_toe[3], rl_toe[3], rr_toe[3];
	static double COM[3]= {0,0,0};
    double Jfr_toe[54], Jfl_toe[54], Jrl_toe[54], Jrr_toe[54];
	double COM_vel[3] = {0,0,0};
	
	q[0] = 0; q[1] = 0; q[2] = 0;
    
    FK_FR_toe(fr_toe, q); FK_FL_toe(fl_toe, q);
    FK_RR_toe(rr_toe, q); FK_RL_toe(rl_toe, q);
    J_FR_toe(Jfr_toe, q); J_FL_toe(Jfl_toe, q);
    J_RR_toe(Jrr_toe, q); J_RL_toe(Jrl_toe, q);
    
	// update change in com pos
	static double fr_prev[3] = {fr_toe[0],fr_toe[1],fr_toe[2]};
	static double fl_prev[3] = {fl_toe[0],fl_toe[1],fl_toe[2]};
	static double rr_prev[3] = {rr_toe[0],rr_toe[1],rr_toe[2]};
	static double rl_prev[3] = {rl_toe[0],rl_toe[1],rl_toe[2]};
	
    double deltaPos[2] = {0.0};
    for(int i=0; i<2; ++i){
        deltaPos[i] -= (fr_toe[i]-fr_prev[i])*contact[0];
        deltaPos[i] -= (fl_toe[i]-fl_prev[i])*contact[1];
        deltaPos[i] -= (rr_toe[i]-rr_prev[i])*contact[2];
        deltaPos[i] -= (rl_toe[i]-rl_prev[i])*contact[3];
        deltaPos[i] /= numContact;
    }    
    
	COM[0] += deltaPos[0];
	COM[1] += deltaPos[1];
    COM[2]  = -1.0*(fr_toe[2]*contact[0]+fl_toe[2]*contact[1]+rr_toe[2]*contact[2]+rl_toe[2]*contact[3])/numContact;
	
	for(int i=0; i<3; ++i){
		fr_prev[i] = fr_toe[i]; fl_prev[i] = fl_toe[i];
		rr_prev[i] = rr_toe[i]; rl_prev[i] = rl_toe[i];		
	}
	
	numContact = (contact[0]+contact[1]) + contact[2]+contact[3];
	Eigen::Matrix<double,3,1> dq_temp = {dq[3],dq[4],dq[5]};
	toWorld(&dq[3],dq_temp,R);
	for (int i = 3; i < 18; ++i){
		COM_vel[0] -= (Jfr_toe[3*i+0]*contact[0] + Jfl_toe[3*i+0]*contact[1] + Jrr_toe[3*i+0]*contact[2] + Jrl_toe[3*i+0]*contact[3])*dq[i];
	 	COM_vel[1] -= (Jfr_toe[3*i+1]*contact[0] + Jfl_toe[3*i+1]*contact[1] + Jrr_toe[3*i+1]*contact[2] + Jrl_toe[3*i+1]*contact[3])*dq[i];
	 	COM_vel[2] -= (Jfr_toe[3*i+2]*contact[0] + Jfl_toe[3*i+2]*contact[1] + Jrr_toe[3*i+2]*contact[2] + Jrl_toe[3*i+2]*contact[3])*dq[i];
	}
	COM_vel[0] /= numContact;
	COM_vel[1] /= numContact;
	COM_vel[2] /= numContact;
	
	dq_temp = {dq[3],dq[4],dq[5]};
	toBody(&dq[3],dq_temp,R);

	// Set results
	q[0] = COM[0]; q[1] = COM[1]; q[2] = COM[2]+0.02;
	dq[0] = COM_vel[0]; dq[1] = COM_vel[1]; dq[2] = COM_vel[2];

}

void ExternalComm::getthetadot(double q[18],double dq[18]){

    Eigen::Matrix<double,3,3> A;
    double phi = q[3];
    double theta = q[4];
    
    A(0,0) = 1;     A(0,1) = sin(phi)*tan(theta);   A(0,2) = cos(phi)*tan(theta);
    A(1,0) = 0;     A(1,1) = cos(phi);              A(1,2) = -sin(phi);
    A(2,0) = 0;     A(2,1) = sin(phi)/cos(theta);   A(2,2) = cos(phi)/cos(theta);

    Eigen::Matrix<double,3,1> dq_temp = {dq[3],dq[4],dq[5]};
    Eigen::Matrix<double,3,1> thetadot = A*dq_temp;
    dq[3] = thetadot(0);
    dq[4] = thetadot(1);
    dq[5] = thetadot(2);
}

void ExternalComm::getStateEstimatefullll(double q[18], double dq[18], int contact[4], Eigen::Matrix<double,3,3> R, Eigen::Matrix<double,3,4> toes, int robotdown, size_t ctrlTick, Eigen::Vector3d acc_bFrame){
    
    float numContact = (contact[0]+contact[1])+rearweight_est*(contact[2]+contact[3]);
    // if(ctrlTick>27399){
        // numContact = (contact[0]+contact[1])*robotdown+rearfootweight*contact[2]+rearfootweight*contact[3];
    // }
	// ================================== //
	// ========= Kin Estimator ========== //
	// ================================== //

	// toe pos
	double fr_toe_e[3], fl_toe_e[3], rl_toe_e[3], rr_toe_e[3];
	double COM_e[3]= {0,0,0};
    double Jfr_toe_e[54], Jfl_toe_e[54], Jrl_toe_e[54], Jrr_toe_e[54];
	double COM_vel_e[3] = {0,0,0};
	

	q[0] = 0; q[1] = 0; q[2] = 0;
    if(robotdown){
	    FK_FR_toe(fr_toe_e, q); FK_FL_toe(fl_toe_e, q);
	    FK_RR_toe(rr_toe_e, q); FK_RL_toe(rl_toe_e, q);
        J_FR_toe(Jfr_toe_e, q); J_FL_toe(Jfl_toe_e, q);
	    J_RR_toe(Jrr_toe_e, q); J_RL_toe(Jrl_toe_e, q);
    }else{
        FK_FR_toe_u(fr_toe_e, q); FK_FL_toe_u(fl_toe_e, q);
	    FK_RR_toe_u(rr_toe_e, q); FK_RL_toe_u(rl_toe_e, q);
        J_FR_toe_u(Jfr_toe_e, q); J_FL_toe_u(Jfl_toe_e, q);
	    J_RR_toe_u(Jrr_toe_e, q); J_RL_toe_u(Jrl_toe_e, q);
    }
	
	// update change in com pos
	double fr_prev_e[3] = {toes(0,0),toes(1,0),toes(2,0)};//{fr_toe[0],fr_toe[1],fr_toe[2]};
	double fl_prev_e[3] = {toes(0,1),toes(1,1),toes(2,1)};//{fl_toe[0],fl_toe[1],fl_toe[2]};
	double rr_prev_e[3] = {toes(0,2),toes(1,2),toes(2,2)};//{rr_toe[0],rr_toe[1],rr_toe[2]};
	double rl_prev_e[3] = {toes(0,3),toes(1,3),toes(2,3)};//{rl_toe[0],rl_toe[1],rl_toe[2]};
	
    double deltaPos[3] = {0.0,0.0,0.0};
    for(int i=0; i<3; ++i){
        // if(ctrlTick<27400){
            deltaPos[i] -= (fr_toe_e[i]-fr_prev_e[i])*contact[0];
            deltaPos[i] -= (fl_toe_e[i]-fl_prev_e[i])*contact[1];
        // }
        deltaPos[i] -= (rr_toe_e[i]-rr_prev_e[i])*contact[2]*rearweight_est;
        deltaPos[i] -= (rl_toe_e[i]-rl_prev_e[i])*contact[3]*rearweight_est;
        deltaPos[i] /= numContact;
    }    
    
	COM_e[0] = deltaPos[0];
	COM_e[1] = deltaPos[1];
    COM_e[2] = deltaPos[2];//-1.0*(fr_toe[2]*contact[0]+fl_toe[2]*contact[1]+rr_toe[2]*contact[2]+rl_toe[2]*contact[3])/numContact;
	
	// for(int i=0; i<3; ++i){
	// 	fr_prev[i] = fr_toe[i]; fl_prev[i] = fl_toe[i];
	// 	rr_prev[i] = rr_toe[i]; rl_prev[i] = rl_toe[i];		
	// }
	// getthetadot(q,dq);
    int robotdown2=0;
	
	
    if(!robotdown){
    // if(ctrlTick>switchtime*1000-1){
            
        numContact = (contact[0]+contact[1])*robotdown2 + rearweight_est*(contact[2]+contact[3]);
        // Eigen::Matrix<double,3,1> dq_temp = {dq[3],dq[4],dq[5]};
        // getthetadot(q,dq);
        // Eigen::Matrix<double,3,1> dq_temp2 = {dq[3],dq[4],dq[5]};
        // toWorld(&dq[3],dq_temp,R);
        // std::cout << simcounter << "\t" << dq_temp[0] << "\t" << dq_temp[1] << "\t" << dq_temp[2] << "\t" << dq[3] << "\t" << dq[4] << "\t" << dq[5] << "\t"
        //                                             << dq_temp2[0] << "\t" << dq_temp2[1] << "\t" << dq_temp2[2] << std::endl;
        // dq[3] = dq_temp[0];
        // dq[4] = dq_temp[1];
        // dq[5] = dq_temp[2];
        for (int i = 3; i < 18; ++i){
		    COM_vel_e[0] -= (Jfr_toe_e[3*i+0]*contact[0]*robotdown2 + Jfl_toe_e[3*i+0]*contact[1]*robotdown2 + Jrr_toe_e[3*i+0]*contact[2]*rearweight_est + Jrl_toe_e[3*i+0]*contact[3]*rearweight_est)*dq[i];
	 	    COM_vel_e[1] -= (Jfr_toe_e[3*i+1]*contact[0]*robotdown2 + Jfl_toe_e[3*i+1]*contact[1]*robotdown2 + Jrr_toe_e[3*i+1]*contact[2]*rearweight_est + Jrl_toe_e[3*i+1]*contact[3]*rearweight_est)*dq[i];
	 	    COM_vel_e[2] -= (Jfr_toe_e[3*i+2]*contact[0]*robotdown2 + Jfl_toe_e[3*i+2]*contact[1]*robotdown2 + Jrr_toe_e[3*i+2]*contact[2]*rearweight_est + Jrl_toe_e[3*i+2]*contact[3]*rearweight_est)*dq[i];
        }
	    COM_vel_e[0] /= numContact;
	    COM_vel_e[1] /= numContact;
	    COM_vel_e[2] /= numContact;

        // std::cout << "accel_scaled" << "\t" << acc_bFrame.transpose() << std::endl;
        // std::cout << "a_world" << "\t" << (R*acc_bFrame).transpose() << std::endl;
        // Eigen::Vector3d a_world = R*acc_bFrame - Eigen::Vector3d(0,0,9.81);
        // dq_temp = {dq[3],dq[4],dq[5]};
	    // toBody(&dq[3],dq_temp,R);

        velocity_filter.step(simcounter,acc_bFrame, R, COM_vel_e);
        Eigen::Vector3d fused = velocity_filter.getVelocity();
        
        std::cout << simcounter << "\t" << dq[0] << "\t" << dq[1] << "\t" << dq[2] << "\t" << fused(0) << "\t" << fused(1) << "\t" << fused(2) << "\t"
                                            << COM_vel_e[0] << "\t" << COM_vel_e[1] << "\t" << COM_vel_e[2] << std::endl;
        COM_vel_e[0] = fused(0);
        COM_vel_e[1] = fused(1);
        COM_vel_e[2] = fused(2);

    }else{
        numContact = (contact[0]+contact[1])*robotdown + rearweight_est*(contact[2]+contact[3]);
        
        Eigen::Matrix<double,3,1> dq_temp = {dq[3],dq[4],dq[5]};
        toWorld(&dq[3],dq_temp,R);
        Eigen::Matrix<double,3,1> dq_temp2 = {dq[3],dq[4],dq[5]};
        // dq[3] = dq_temp(0);
        // dq[4] = dq_temp(1);
        // dq[5] = dq_temp(2);
        // getthetadot(q,dq);
        
        // std::cout << simcounter << "\t" << dq_temp[0] << "\t" << dq_temp[1] << "\t" << dq_temp[2] << "\t" << dq[3] << "\t" << dq[4] << "\t" << dq[5] << "\t"
        //                                             << dq_temp2[0] << "\t" << dq_temp2[1] << "\t" << dq_temp2[2] << std::endl;
                
	    for (int i = 3; i < 18; ++i){
		    COM_vel_e[0] -= (Jfr_toe_e[3*i+0]*contact[0]*robotdown + Jfl_toe_e[3*i+0]*contact[1]*robotdown + Jrr_toe_e[3*i+0]*contact[2]*rearweight_est + Jrl_toe_e[3*i+0]*contact[3]*rearweight_est)*dq[i];
	 	    COM_vel_e[1] -= (Jfr_toe_e[3*i+1]*contact[0]*robotdown + Jfl_toe_e[3*i+1]*contact[1]*robotdown + Jrr_toe_e[3*i+1]*contact[2]*rearweight_est + Jrl_toe_e[3*i+1]*contact[3]*rearweight_est)*dq[i];
	 	    COM_vel_e[2] -= (Jfr_toe_e[3*i+2]*contact[0]*robotdown + Jfl_toe_e[3*i+2]*contact[1]*robotdown + Jrr_toe_e[3*i+2]*contact[2]*rearweight_est + Jrl_toe_e[3*i+2]*contact[3]*rearweight_est)*dq[i];
	    }
	    COM_vel_e[0] /= numContact;
	    COM_vel_e[1] /= numContact;
	    COM_vel_e[2] /= numContact;
	
	    dq_temp = {dq[3],dq[4],dq[5]};
	    toBody(&dq[3],dq_temp,R);
        // dq[3] = dq_temp[0];
        // dq[4] = dq_temp[1];
        // dq[5] = dq_temp[2];
    }

	// Set results
	q[0] = COM_e[0]; q[1] = COM_e[1]; q[2] = COM_e[2];
    if(ctrlTick<27000){
	    dq[0] = COM_vel_e[0] > xdot_thresh ? xdot_thresh : (COM_vel_e[0] < -xdot_thresh ? -xdot_thresh : COM_vel_e[0]); 
        dq[1] = COM_vel_e[1] > yzdot_thresh ? yzdot_thresh : (COM_vel_e[1] < -yzdot_thresh ? -yzdot_thresh : COM_vel_e[1]);
        dq[2] = COM_vel_e[2] > yzdot_thresh ? yzdot_thresh : (COM_vel_e[2] < -yzdot_thresh ? -yzdot_thresh : COM_vel_e[2]); 
        //dq[2] = COM_vel[2];
    }else{
        dq[0] = COM_vel_e[0] > xdot_thresh2 ? xdot_thresh2 : (COM_vel_e[0] < -xdot_thresh2 ? -xdot_thresh2 : COM_vel_e[0]); 
        dq[1] = COM_vel_e[1] > yzdot_thresh2 ? yzdot_thresh2 : (COM_vel_e[1] < -yzdot_thresh2 ? -yzdot_thresh2 : COM_vel_e[1]);
        dq[2] = COM_vel_e[2] > yzdot_thresh2 ? yzdot_thresh2 : (COM_vel_e[2] < -yzdot_thresh2 ? -yzdot_thresh2 : COM_vel_e[2]); 
    }
}



void ExternalComm::SimExec(){//(std::ofstream &file_est){
 
    //std::cout << "InSimExec" << std::endl;
    if (setup_raisim){
        setupRaisim();
        setup_raisim = false; 
    }
    updateData(GET_DATA, SIM_DATA, &SimData);
    if(!vis->getRoot()->endRenderingQueued() && simcounter < simlength){
        
        A1.back()->setGeneralizedForce(SimData.tau);
        // std::vector<double> GRF(std::begin(HLData.fDes), std::end(HLData.fDes)-5);
		// std::vector<double> feet_vec(std::begin(SimData.toePos), std::end(SimData.toePos));
		// std::vector<double> contacts(std::begin(SimData.ind_LL), std::end(SimData.ind_LL));
        // if(GRF.size() >0) plotGRFs(list, GRF, feet_vec, contacts);
        //plotGRFs(list, HLData.fDes, SimData.toePos, SimData.ind_LL);
        world.integrate();        
        
        if (simcounter%60 == 0)
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
        
        // std::cout << "simcounter" << "\t" << simcounter << std::endl;
        simcounter++; 
        
    }
	else if (simcounter > simlength-1){
        if (vis->isRecording()){vis->stopRecordingVideoAndSave();}
        vis->closeApp();
    }

    //////////////////////////////////
    //      STATE ESTIMATION        //
    //////////////////////////////////

    double jpos[18],jpos_est[18], jvel[18],jvel_est[18];
    double rotMatrixDouble[9] = {1,0,0,0,1,0,0,0,1};

    Eigen::Matrix<double, 3, 1> eul = Eigen::MatrixXd::Zero(3,1);
    Eigen::Matrix<double, 4, 1> quat = Eigen::MatrixXd::Zero(4,1);
    Eigen::Matrix<double, 3, 1>  eul_state = Eigen::MatrixXd::Zero(3,1);
    Eigen::Matrix<double, 3, 1>  omega_state = Eigen::MatrixXd::Zero(3,1);
    Eigen::Matrix<double,3,1> acc_bFrame = Eigen::MatrixXd::Zero(3,1);

    A1.back()->getState(jointPosTotal, jointVelTotal);
    A1.back()->getBaseOrientation(rotMat);

    int robotdown = simcounter < switchtime*ctrlHz ? 1 : 0;
    
    if(!robotdown){
        Eigen::Matrix<double,3,3> rotIMU = Eigen::MatrixXd::Zero(3,3);
        if(SIMFlag){
            auto imu = A1.back()->getSensorSet("imu_parent")->getSensor<raisim::InertialMeasurementUnit>("imu");
            auto imu_o = imu->getOrientation();   // Quaternion
            auto imu_w = imu->getAngularVelocity();  // Angular velocity in radians/s
            raisim::Vec<3> linearAcceleration = imu->getLinearAcceleration();
            acc_bFrame(0) = linearAcceleration(0);
            acc_bFrame(1) = linearAcceleration(1);
            acc_bFrame(2) = linearAcceleration(2);
            quat(0) = imu_o[0];
            quat(1) = imu_o[1];
            quat(2) = imu_o[2];
            quat(3) = imu_o[3];
            quat_to_XYZ(quat,eul_state);
            quat_to_R(quat,rotIMU);
            omega_state(0) = imu_w[0];
            omega_state(1) = imu_w[1];
            omega_state(2) = imu_w[2];
            
        }else{
            eul_state(0) = SimData.att_euler[0];
            eul_state(1) = SimData.att_euler[1];
            eul_state(2) = SimData.att_euler[2];
            R_XYZ(eul_state,rotIMU);
            omega_state(0) = SimData.comp_angular_rate[0];
            omega_state(1) = SimData.comp_angular_rate[1];
            omega_state(2) = SimData.comp_angular_rate[2];       
        }
        for(size_t i=0;i<3;i++){
            for (size_t j = 0; j < 3; j++){
                rotMat[3*i+j] = rotIMU(j,i);
            }
        } 
        // }else{
            // for(size_t i=0;i<3;i++){
            //     for (size_t j = 0; j < 3; j++){
            //         rotIMU(j,i)=rotMat[3*i+j];
            //     }
            // } 
            // rotIMU(0,0) = rotMat(0,0)   ;   rotIMU(0,1) =  rotMat(0,1)  ;   rotIMU(0,2) =   rotMat(0,2);
            // rotMat(0,0) = rotMat(2,0)  ;   rotMat(0,1) =  rotMat(2,1) ;   rotMat(0,2) =   rotMat(2,2);
            // rotMat(2,0) = -rotIMU(0,0)   ;   rotMat(2,1) =  -rotIMU(0,1)  ;   rotMat(2,2) =   -rotIMU(0,2);
            // eul_state(0) = atan2(rotMat(1,2),rotMat(2,2));
		    // eul_state(1) = -asin(rotMat(0,2));
		    // eul_state(2) = atan2(rotMat(0,1),rotMat(0,0));
        // }  
               
    }else{
        quat = jointPosTotal.block(3,0,4,1);
        quat_to_XYZ(quat,eul_state);
    }

    // std::cout << rotMat(0,0) << "\t" << rotMat(0,1) << "\t" << rotMat(0,2) << std::endl;

    for(size_t i = 0 ; i < 9 ; i++){
        rotMatrixDouble[i] = rotMat[i];
    }
    
    Eigen::Map< Eigen::Matrix<double, 3, 3> > rotE(rotMatrixDouble, 3, 3);
    
    if(robotdown){
        omega_state = rotE.transpose()*jointVelTotal.segment(3,3); // convert to body frame, like robot measurements
        // omega_state = jointVelTotal.segment(3,3);
    }
    
    for(size_t i=0; i<3; ++i){
        jpos[i] = jointPosTotal(i);
        jvel[i] = jointVelTotal(i);
        jpos_est[i] = jointPosTotal(i);
        jvel_est[i] = jointVelTotal(i);
        jpos[i+3] = eul_state(i);
        jvel[i+3] = omega_state(i);
        jpos_est[i+3] = eul_state(i);
        jvel_est[i+3] = omega_state(i);
    }

    for(size_t i=6; i<18; ++i){
        jpos[i] = jointPosTotal(i+1);
        jvel[i] = jointVelTotal(i);
        jpos_est[i] = jointPosTotal(i+1);
        jvel_est[i] = jointVelTotal(i);
    }
    
    // int raisimForce[4] = {0};
    // for(auto &con: A1.back()->getContacts()){
    //     int conInd = con.getlocalBodyIndex();//SimData.ind_LL;
    //     raisimForce[conInd/3-1] = 500;
    //     raisimForce[conInd/3-1] = con.getNormal().e().norm();
    // }

    // Assuming raisim::World world and raisim::ArticulatedSystem* a1 are already created and set up
    // auto& contacts = world.getContactForce();  // get all contact points
    // For mapping contact forces by foot name
    // std::map<std::string, Eigen::Vector3d> footContactForces = {
    //     {"FR_foot", Eigen::Vector3d::Zero()},
    //     {"FL_foot", Eigen::Vector3d::Zero()},
    //     {"RR_foot", Eigen::Vector3d::Zero()},
    //     {"RL_foot", Eigen::Vector3d::Zero()}
    // };
    // Iterate through contacts
    // for (const auto& contact : A1.back()->getContacts()) {
    //     std::string contactBodyName = A1.back()->getBodyNames()[contact.getlocalBodyIndex()];

    //     if (footContactForces.find(contactBodyName) != footContactForces.end()) {
    //         Eigen::Vector3d force = Eigen::Vector3d(contact.getImpulse().e());  // in world frame
    //         footContactForces[contactBodyName] = force;

    //         csvFile << contactBodyName << " Contact Force: "
    //                 << "X: " << force.x() << ", "
    //                 << "Y: " << force.y() << ", "
    //                 << "Z: " << force.z() << std::endl;
    //     }
    // }
    csvFile << simcounter << ",";
    for(int i=0; i<4; ++i){
        A1.back()->getPosition(calfIdx[i], localPosition, calfPose);
        // jointFR.getPosition(calfPose);
        csvFile << calfPose[0] << "," << calfPose[1] << "," << calfPose[2] << ",";
    }
    csvFile << "\n";

    if(simcounter>2499){
        getStateEstimatefullll(jpos_est,jvel_est,SimData.ind_LL,rotE,SimData.toePos,robotdown,simcounter,acc_bFrame);
        
    }else if(simcounter>0){
        kinestimatorrr(jpos_est,jvel_est,SimData.ind_LL,rotE);
    }

    // file_est << simcounter << "," << jpos[0] << "," << jpos[1] << "," << jpos[2] << "," << jvel[0] << "," << jvel[1] << "," << jvel[2] << ","
    //      << jpos[3] << "," << jpos[4] << "," << jpos[5] << "," << jvel[3] << "," << jvel[4] << "," << jvel[5] << ","
    //      << jpos_est[0] << "," << jpos_est[1] << "," << jpos_est[2] << "," << jvel_est[0] << "," << jvel_est[1] << "," << jvel_est[2] << ","
    //      << jpos_est[3] << "," << jpos_est[4] << "," << jpos_est[5] << "," << jvel_est[3] << "," << jvel_est[4] << "," << jvel_est[5] << ","
        // << imu_eul(0) << "," << imu_eul(1) << "," << imu_eul(2) << "," << imu_omega(0) << "," << imu_omega(1) << "," << imu_omega(2) << ","
        //  << rotE(0,0) << "," << rotE(0,1) << "," << rotE(0,2) << "," 
        //  << rotE(1,0) << "," << rotE(1,1) << "," << rotE(1,2) << ","
        //  << rotE(2,0) << "," << rotE(2,1) << "," << rotE(2,2) << ","
        //  << vel_temp[0] << "," << vel_temp[1] << "," << vel_temp[2] 
        // << "\n";
   
    memcpy(SimData.q,jpos_est,18*sizeof(double));
    memcpy(SimData.dq,jvel_est,18*sizeof(double));
	memcpy(SimData.rotMatrixDouble,rotMatrixDouble,9*sizeof(double));
    SimData.control_Tick = simcounter;

	// Set Updated data for MPC/LL
	updateData(SET_DATA, SIM_DATA, &SimData);  
    
}

void ExternalComm::connectIMU(){//(mip::Interface& device){
    // if (!utils->device) {
    //     std::cerr << "Failed to create device" << std::endl;
    //     //return -1;
    // }
    
    // Ping the device to verify communication
    if (mip::commands_base::ping(*device) != mip::CmdResult::ACK_OK) {
        std::cerr << "ERROR: Could not ping the device!" << std::endl;
        //return -1;
    }
    
    std::cout << "Successfully pinged the device!" << std::endl;
    float gyro_bias[3] = {0, 0, 0};
    // Rest of your code...
    const uint32_t sampling_time = 2000; // The default is 15000 ms and longer sample times are recommended but shortened for convenience
    const mip::Timeout old_mip_sdk_timeout = device->baseReplyTimeout();
    printf("Capturing gyro bias. This will take %d seconds \n", sampling_time/1000);
    device->setBaseReplyTimeout(sampling_time * 2);

    if(mip::commands_3dm::captureGyroBias(*device, sampling_time, gyro_bias) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not capture gyro bias!");

    if(mip::commands_3dm::saveGyroBias(*device) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not save gyro bias!");
    
    const uint8_t device_selector = 3;
    const uint8_t enable_flag = 1;
    if(mip::commands_3dm::writeDatastreamControl(*device, device_selector, enable_flag) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not enable device data stream!");

    // Reset the timeout
    device->setBaseReplyTimeout(old_mip_sdk_timeout);

    printf("Gyro bias captured with sampling time: %d, and gyro bias captured as: %f %f %f.\n", sampling_time, gyro_bias[0], gyro_bias[1], gyro_bias[2]);

}

void ExternalComm::configureIMU(){//(mip::Interface& device){
    uint16_t sensor_base_rate;

    //Note: Querying the device base rate is only one way to calculate the descriptor decimation.
    //We could have also set it directly with information from the datasheet.

    if(mip::commands_3dm::imuGetBaseRate(*device, &sensor_base_rate) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not get sensor base rate format!");

    const uint16_t sensor_sample_rate = 1000; // Hz
    const uint16_t sensor_decimation = sensor_base_rate / sensor_sample_rate;

    std::array<mip::DescriptorRate, 4> sensor_descriptors = {{
        { mip::data_sensor::DATA_TIME_STAMP_GPS, sensor_decimation },
        { mip::data_sensor::DATA_ACCEL_SCALED,   sensor_decimation },
        { mip::data_sensor::DATA_GYRO_SCALED,    sensor_decimation },
        { mip::data_sensor::DATA_COMP_EULER_ANGLES, sensor_decimation },
    }};

    if(mip::commands_3dm::writeImuMessageFormat(*device, static_cast<uint8_t>(sensor_descriptors.size()), sensor_descriptors.data()) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not set sensor message format!");

}

void ExternalComm::setupIMUfilter(){//(mip::Interface& device){
    uint16_t filter_base_rate;

    if(mip::commands_3dm::filterGetBaseRate(*device, &filter_base_rate) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not get filter base rate format!");

    const uint16_t filter_sample_rate = 1000; // Hz
    const uint16_t filter_decimation = filter_base_rate / filter_sample_rate;

    std::array<mip::DescriptorRate, 5> filter_descriptors = {{
        { mip::data_filter::DATA_FILTER_TIMESTAMP, filter_decimation },
        { mip::data_filter::DATA_FILTER_STATUS,    filter_decimation },
        { mip::data_filter::DATA_ATT_EULER_ANGLES, filter_decimation },
        { mip::data_filter::DATA_COMPENSATED_ANGULAR_RATE, filter_decimation },
        { mip::data_filter::DATA_COMPENSATED_ACCELERATION, filter_decimation },
    }};

    if(mip::commands_3dm::writeFilterMessageFormat(*device, static_cast<uint8_t>(filter_descriptors.size()), filter_descriptors.data()) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not set filter message format!");

    if(mip::commands_filter::writeSensorToVehicleRotationEuler(*device, sensor_to_vehicle_rotation_euler[0], sensor_to_vehicle_rotation_euler[1], sensor_to_vehicle_rotation_euler[2]) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not set sensor-2-vehicle rotation!");

    if(mip::commands_filter::writeAutoInitControl(*device, 1) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not set filter autoinit control!");

    if(mip::commands_filter::reset(*device) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not reset the filter!");

    device->registerExtractor(sensor_data_handlers[0], &sensor_gps_time);
    device->registerExtractor(sensor_data_handlers[1], &sensor_accel);
    device->registerExtractor(sensor_data_handlers[2], &sensor_gyro);
    device->registerExtractor(sensor_data_handlers[3], &sensor_comp_euler_angles);

    //Filter Data
    device->registerExtractor(filter_data_handlers[0], &filter_gps_time);
    device->registerExtractor(filter_data_handlers[1], &filter_status);
    device->registerExtractor(filter_data_handlers[2], &filter_euler_angles);
    device->registerExtractor(filter_data_handlers[3], &filter_comp_angular_rate);
    device->registerExtractor(filter_data_handlers[4], &filter_comp_accel);

    if(mip::commands_base::resume(*device) != mip::CmdResult::ACK_OK)
        printf("ERROR: Could not resume the device!");

    mip::Timestamp prev_print_timestamp = getCurrentTimestamp();
    printf("Sensor is configured... waiting for filter to enter running mode.\n");

}

void ExternalComm::getIMMUdata(){//(mip::Interface& device){
    
    // if (!device) {
    //     printf("ERROR: Device pointer is null\n");
    //     return;
    // }
    
    device->update();
    //Check Filter State
    if((!this->filter_state_running) && ((this->filter_status.filter_state == mip::data_filter::FilterMode::GX5_RUN_SOLUTION_ERROR) || (this->filter_status.filter_state == mip::data_filter::FilterMode::GX5_RUN_SOLUTION_VALID)))
    {
        printf("NOTE: Filter has entered running mode.\n");
        this->filter_state_running = true;
    }
    //Once in running mode, print out data at 1 Hz
    if(this->filter_state_running)
    {
        auto now = std::chrono::system_clock::now();
        auto unix_timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
    
        // printf("Timestamp = %lld ms: TOW = %f: ATT_EULER = [%f %f %f]: COMP_ANG_RATE = [%f %f %f]\n",//: COMP_ACCEL = [%f %f %f]\n",
        //         unix_timestamp,
        //         this->filter_gps_time.tow, 
        //         this->filter_euler_angles.roll, 
        //         this->filter_euler_angles.pitch, 
        //         this->filter_euler_angles.yaw,
        //         // this->filter_comp_angular_rate.gyro[0], 
        //         // this->filter_comp_angular_rate.gyro[1], 
        //         // this->filter_comp_angular_rate.gyro[2]);
        //         this->sensor_comp_euler_angles.roll,
		// 		this->sensor_comp_euler_angles.pitch,
		// 		this->sensor_comp_euler_angles.yaw);
        IMUData.att_euler[0] = this->filter_euler_angles.roll;
        IMUData.att_euler[1] = this->filter_euler_angles.pitch;
        IMUData.att_euler[2] = this->filter_euler_angles.yaw;
        IMUData.comp_angular_rate[0] = this->filter_comp_angular_rate.gyro[0];
        IMUData.comp_angular_rate[1] = this->filter_comp_angular_rate.gyro[1];
        IMUData.comp_angular_rate[2] = this->filter_comp_angular_rate.gyro[2];        
        updateData(SET_DATA, IMU_DATA, &IMUData);                          
    }
    // }catch (const std::exception& e) {
    //     printf("ERROR in getIMMUdata: %s\n", e.what());
    // }
    // catch (...) {
    //     printf("Unknown ERROR in getIMMUdata\n");
    // }
}





int main(int argc, char *argv[]) {

    
    ExternalComm extComm(argc, argv);
    extComm.loco_obj = std::unique_ptr<LocoWrapper>(new LocoWrapper(argc, argv));
    extComm.nmpc_obj  = std::unique_ptr<SRBNMPC>(new SRBNMPC(argc,argv,1,0));
    
    // extComm.connectIMU();
    // extComm.configureIMU();
    // extComm.setupIMUfilter();
    
    int simIMU = 0;

    LoopFunc loop_calc("calc_loop", extComm.LLdt,1, boost::bind(&ExternalComm::Calc, &extComm));
	LoopFunc loop_mpc("mpc_loop", extComm.HLdt,2, boost::bind(&ExternalComm::HighLevel, &extComm));
	LoopFunc loop_sim("sim_loop", extComm.LLdt,3, boost::bind(&ExternalComm::SimExec, &extComm));
    // // LoopFunc loop_imu("imu_loop", extComm.LLdt,4, boost::bind(&ExternalComm::getIMMUdata, &extComm));
	
	loop_sim.start();
	sleep(1.0);
	loop_mpc.start();
	sleep(1.0);
	loop_calc.start();

    // sleep(1.0);
    // loop_imu.start();
    // // loop_flush.start();

    while(true)// (simIMU < 500000)
    {
        sleep(0.1);
        // extComm.getIMUread2();
        // simIMU++;
    }
    
    // std::ofstream file_est("../data25/estimatorMT13.csv");
    // while (true)
	// {
			
    //     // sleep(0.1);
    //     // extComm.getIMUread();
    //     extComm.SimExec();//(file_est);
    //     // std::cout << "SimExec" << std::endl;
    //     extComm.HighLevel();
    //     // std::cout << "HighLevel" << std::endl;
    //     extComm.Calc();
    //     // std::cout << "Calc" << std::endl;
    //     // sim_setup = false;

	// } 

    // file_est.close();

    
    return 0;
}


