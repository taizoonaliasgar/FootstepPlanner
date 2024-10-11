//
// Created by Jeeseop Kim on 02/08/2022.
// MIT License
//
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "raisim/OgreVis.hpp"
#include "raisimBasicImguiPanel.hpp"
#include "raisimKeyboardCallback.hpp"
#include "helper.hpp"
//for Unity
#include "raisim/RaisimServer.hpp"
#include "raisim/World.hpp"
#include <iostream>
#include <filesystem>
#include "math_define.h"

//HDSRL header
//#include "v60_controller/locomotion_controller.h"
//#include "locomotion_planner.h"
#include "SRBNMPCR.hpp"

//For multithreading
#include "pthread.h"

//header helps stop while loop with ctrl+c
#include "signal.h"

//Boost to std-C
//#include "chrono"
//#include "thread"

// using namespace std::chrono_literals;
using namespace std;
namespace fs = std::filesystem;

size_t signaled = 0;
void signalhandler(int param){
    signaled = 1;
}

void writeMatrixToFile(const casadi::DM& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Print matrix dimensions
        file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << std::endl;
        // Print matrix content preserving original shape
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

void setupCallback() {
    auto vis = raisim::OgreVis::get();

    /// light
    // light setting 1
    vis->getLight()->setDiffuseColour(1, 1, 1);
    vis->getLight()->setCastShadows(true);
    //Ogre::Vector3 lightdir(-3,-3,-0.5);
    Ogre::Vector3 lightdir(-1,-1,-1);
    lightdir.normalise();
    vis->getLightNode()->setDirection({lightdir});
    
    /* //light setting 2
    vis->getLight()->setDiffuseColour(1, 1, 1);
    vis->getLight()->setCastShadows(true);
    vis->setCameraSpeed(300);
    Ogre::ColourValue rgba(0.5, 0.5, 0.5, 1);
    vis->setAmbientLight(rgba);
    Ogre::Vector3 lightdir(0.1, 0, -1);
    Ogre::Vector3 lightdir2(0, 0, 3);
    lightdir.normalise();
    lightdir2.normalise();
    vis->getLightNode()->setDirection({lightdir});
    vis->addLight("Light2");
    vis->getLightNode("Light2")->setDirection({lightdir2});
    */

    vis->setCameraSpeed(300);
    
    /// load  textures
    vis->addResourceDirectory(vis->getResourceDir() + "/material/checkerboard");
    vis->loadMaterialFile("checkerboard.material");

    vis->addResourceDirectory(vis->getResourceDir() + "/material/skybox/violentdays");
    vis->loadMaterialFile("violentdays.material");

    //std::string currpath;
    //currpath = binaryPath.getDirectory();
    vis->addResourceDirectory("/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/material_collection");
    vis->loadMaterialFile("collection.material");
    
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

    //sky setting
    Ogre::Quaternion quat;
    quat.FromAngleAxis(Ogre::Radian(M_PI_2), {1,0,0});
    vis->getSceneManager()->setSkyBox(true,
                                      "black", //"white" //"skybox/violentdays"
                                      500,
                                      true,
                                      quat,
                                      Ogre::ResourceGroupManager::AUTODETECT_RESOURCE_GROUP_NAME);
}

// ================================================================ //
// ================== functions for convenience =================== //
// ================================================================ //
// void printAllBodyNames(raisim::OgreVis* vis){
//     // This function prints the names of all of the visual bodies being rendered in the
//     // Ogre environment
//     Ogre::SceneManager* scmg = vis->getSceneManager(); 
//     Ogre::SceneManager::MovableObjectIterator it = scmg->getMovableObjectIterator("Entity");
//     while(it.hasMoreElements()){
//         std::cout<<it.getNext()->getName()<<std::endl;
//     }
// }

void printAllBodyNames(raisim::OgreVis *vis){
    // This function prints the names of all of the visual bodies being rendered in the
    // Ogre environment

    Ogre::SceneManager *scmg = vis->getSceneManager(); 
    Ogre::SceneManager::MovableObjectMap obj = scmg->getMovableObjects("Entity");
    Ogre::SceneManager::MovableObjectMap::iterator it = obj.begin();
    while(it!=obj.cend()){
        std::cout<<(*it).second->getName()<<std::endl;
        it++;
    }
}

void raisim2EigenVec(raisim::Vec<3> a, Eigen::Vector3d& b){
    for(size_t i=0; i<3; i++){
        b(i) = a(i);
    }
}

void raisim2EigenMat(raisim::Mat<3,3> a, Eigen::Matrix<double,3,3>& b){
    b<< a(0,0),a(0,1),a(0,2),a(1,0),a(1,1),a(1,2),a(2,0),a(2,1),a(2,2); //row major
}

void eigen2RaisimVec(Eigen::Vector3d a, raisim::Vec<3>& b){
    for(size_t i=0; i<3; i++){
        b(i) = a(i);
    }
}

void visExtforce(std::vector<raisim::ArticulatedSystem*> srbx2, raisim::Vec<3> extforce, raisim::Vec<3>& scale, raisim::Mat<3,3>& rot){
    raisim::Vec<3> dir;
    dir = extforce;
    dir /= dir.norm();
    raisim::zaxisToRotMat(dir, rot);
    scale = {0.2, 0.2, 0.01*extforce.norm()};
}

void setExtforcewithvis(std::vector<raisim::ArticulatedSystem *> srbx2, size_t agnum, std::map<std::string, raisim::VisualObject>& objlist,
                 const std::string& objname, size_t bodyidx, Eigen::Vector3d bodypos, Eigen::MatrixXd bodyrot, Eigen::Vector3d offset, Eigen::Vector3d extforce_e){
    Eigen::Vector3d offset_inertial;
    offset_inertial = offset;//bodyrot * offset; //offset;// 
    offset_inertial -= bodypos;

    if(extforce_e.norm()<1e-4){
        extforce_e *= 0;
    }

    srbx2[agnum]->setExternalForce(bodyidx, {offset_inertial(0), offset_inertial(1), offset_inertial(2)}, {extforce_e(0), extforce_e(1), extforce_e(2)});
    //srbx2[agnum]->setExternalForce(bodyidx,raisim::ArticulatedSystem::Frame::WORLD_FRAME,{extforce_e(0), extforce_e(1), extforce_e(2)},raisim::ArticulatedSystem::Frame::WORLD_FRAME,{offset_inertial(0), offset_inertial(1), offset_inertial(2)});
    objlist[objname].offset = {offset_inertial(0) + bodypos(0), offset_inertial(1)+ bodypos(1), offset_inertial(2)+bodypos(2)};
    //objlist[objname].offset = {offset_inertial(0) + bodypos(0), offset_inertial(1)+ bodypos(1), 0};
    visExtforce(srbx2, {extforce_e(0),extforce_e(1),extforce_e(2)}, objlist[objname].scale, objlist[objname].rotationOffset);
}

void setdisturbforcewithvis(std::vector<raisim::ArticulatedSystem *> srbx2, size_t agnum, std::map<std::string, raisim::VisualObject>& objlist,
                 const std::string& objname, size_t bodyidx, Eigen::Vector3d bodypos, Eigen::MatrixXd bodyrot, Eigen::Vector3d offset, Eigen::Vector3d extforce_e){
    Eigen::Vector3d offset_inertial;
    offset_inertial = offset;//bodyrot * offset; //offset;// 
    //offset_inertial -= bodypos;
    //std::cout <<offset_inertial.transpose() << std::endl;

    srbx2[agnum]->setExternalForce(bodyidx, {offset_inertial(0), offset_inertial(1), offset_inertial(2)}, {extforce_e(0), extforce_e(1), extforce_e(2)});
    objlist[objname].offset = {offset_inertial(0) + bodypos(0), offset_inertial(1)+ bodypos(1), offset_inertial(2)+bodypos(2)};
    //objlist[objname].offset = {offset_inertial(0) + bodypos(0), offset_inertial(1)+ bodypos(1), 0};
    visExtforce(srbx2, {extforce_e(0),extforce_e(1),extforce_e(2)}, objlist[objname].scale, objlist[objname].rotationOffset);
}
// =============================================== //
// ================ planner ====================== //
// =============================================== //
void planner(std::vector<raisim::ArticulatedSystem*> srbx2, 
                size_t totalagentnumber, 
                size_t eachagentnumber,
                SRBNMPCR* mpc_obj, casadi::Function solver, casadi::DM contact_sequence){

    // ================================ //
    // ========= INITIALIZE =========== //
    // ================================ //
    

    static auto& list = raisim::OgreVis::get()->getVisualObjectList();
    static size_t controlTick = 0;
    

    size_t settling_temp  = 0*ctrlHz;   // PD control
    size_t settling = 0*ctrlHz;         // start torque control and settling down
    size_t duration = 0*ctrlHz;         // standing up
    size_t loco_start = settling + duration;
    //loco_pln->timeSetup(duration, loco_start);

    size_t loco_kind = TROT;        // kind of locomotion user want to execute
    size_t loco_status = STAND;     // current status of locomotion

    // =========================================================== //
    // ================ subscribing body states ================== //
    // =========================================================== //
    raisim::Vec<3> bodypos = {0,0,0};
    raisim::Mat<3,3> bodyrot = {0,0,0,0,0,0,0,0,0};
    raisim::Vec<3> bodyvel = {0,0,0};
    raisim::Vec<3> bodyomega = {0,0,0};
    
    Eigen::Vector3d bodypos_e;
    Eigen::Matrix<double,3,3> bodyrot_e;
    Eigen::Vector3d bodyvel_e;
    Eigen::Vector3d bodyomega_e;

    bodypos_e.setZero();
    bodyrot_e.setZero(3,3);
    bodyvel_e.setZero();
    bodyomega.setZero();

    Eigen::VectorXd jointPosTotal = Eigen::MatrixXd::Zero(7,1);
    Eigen::VectorXd jointVelTotal = Eigen::MatrixXd::Zero(6,1);
    double jpos[18], jvel[18];
    Eigen::Matrix<double, 3, 1> eul;
    Eigen::Matrix<double, 4, 1> quat;
    Eigen::Matrix<double, 4, 1> quatcomposition;

    if(eachagentnumber == 0){
        srbx2[eachagentnumber]->getPosition(0, bodypos);
        srbx2[eachagentnumber]->getOrientation(0, bodyrot);
        srbx2[eachagentnumber]->getVelocity(0, bodyvel);
        srbx2[eachagentnumber]->getAngularVelocity(0,bodyomega);

        srbx2[eachagentnumber]->getState(jointPosTotal, jointVelTotal);
        quat = jointPosTotal.block(3,0,4,1);

        quat_to_XYZ(quat,eul);
        
        for(size_t i=0; i<3; ++i){
            jpos[i] = jointPosTotal(i);
            jvel[i] = jointVelTotal(i);
            jpos[i+3] = eul(i);
            jvel[i+3] = jointVelTotal(i+3);
        }

    }

    raisim2EigenVec(bodypos, bodypos_e);
    raisim2EigenMat(bodyrot, bodyrot_e);
    raisim2EigenVec(bodyvel, bodyvel_e);
    raisim2EigenVec(bodyomega, bodyomega_e);


    // ================================================= //
    // ================ control input ================== //
    // ================================================= //
    Eigen::Matrix<double, 12, 1> forceFFvec;
    forceFFvec.setZero(12,1);

    Eigen::Matrix<double, 3, 1> disturbance;
    disturbance.setZero(3,1);

    Eigen::Matrix<double, 12, 1> p_foot;
    p_foot.setZero(12,1);

    Eigen::Matrix<double, 3, 4> p_foot_mpc;
    p_foot_mpc.setZero(3,4);
    Eigen::Matrix<double,17,1> q0;
    q0.setZero(17,1);
    // =========================================== //
    // ================ CONTROL ================== //
    // =========================================== //
    double xoffset = 0.183;//0.267/2;
    double yoffset = 0.1321;//0.194/2;
    double a1mass = 12.4530;
    Eigen::Vector3d extf(0, 0, a1mass*9.81*0.25);
    Eigen::Vector3d extf1(0, 0, a1mass*9.81*0.5);
    Eigen::Vector3d extf0(0, 0, 0);
    Eigen::Vector3d offset0( xoffset, -yoffset, 0); //right front (a1 rule)
    Eigen::Vector3d offset1( xoffset,  yoffset, 0); //left front (a1 rule)
    Eigen::Vector3d offset2(-xoffset, -yoffset, 0); //right hind (a1 rule)
    Eigen::Vector3d offset3(-xoffset,  yoffset, 0); //left hind (a1 rule)
   
    q0.block(0,0,3,1) = jointPosTotal.block(0,0,3,1);
    q0.block(3,0,3,1) = jointVelTotal.block(0,0,3,1);
    q0.block(6,0,3,1) = eul;
    q0.block(9,0,3,1) = jointVelTotal.block(3,0,3,1);
    q0(16,0) = controlTick;
    
    casadi::DM X_prev = casadi::DM::zeros(NFSR*(HORIZ+1)+NFIR*HORIZ,1); 
    std::map<std::string, casadi::DM> arg, res;
    
    int fit_order = 6;
    Eigen::Matrix<double,3,1> distforce = {50,0,0};
    //robot loadingq

    if(controlTick < settling_temp){
        if(eachagentnumber == 0){ 
            setExtforcewithvis(srbx2, 0, list, "extForceArrow0", 0, bodypos_e, bodyrot_e, offset0, extf);
            setExtforcewithvis(srbx2, 0, list, "extForceArrow1", 0, bodypos_e, bodyrot_e, offset1, extf);
            setExtforcewithvis(srbx2, 0, list, "extForceArrow2", 0, bodypos_e, bodyrot_e, offset2, extf);
            setExtforcewithvis(srbx2, 0, list, "extForceArrow3", 0, bodypos_e, bodyrot_e, offset3, extf);
        }

    }
    //locomotion
    else if(controlTick >= settling_temp){
        
        if(controlTick%10==0){
            
            int controlMPC = std::floor(controlTick/10);

            if(controlTick<1){
                q0.block(12,0,4,1) << 0.2,0.2,-0.1,-0.1;
            }else{
                //Eigen::Matrix<double, 4, 1> prevp = mpc_obj->getpreviousfoot();
                q0.block(12,0,4,1) = mpc_obj->getpreviousfoot();
            }
    
            X_prev = mpc_obj->getprevioussol(q0,controlMPC); 
    
            casadi::DM p = mpc_obj->motionPlannerN(q0,controlMPC);
    
            arg["lbx"] = mpc_obj->lowerboundx(p); 
            arg["ubx"] =  mpc_obj->upperboundx(p); 
            arg["lbg"] =  mpc_obj->lowerboundg(); 
            arg["ubg"] =  mpc_obj->upperboundg();
    
            arg["x0"] = X_prev;
            arg["p"] = p;

            std::string basePath = "/home/taizoon/raisimEnv/raisimWorkspace/footstep_planner/build/tmp/";
            writeMatrixToFile(arg["lbx"], basePath + "args_lbx.txt");
            writeMatrixToFile(arg["ubx"], basePath + "args_ubx.txt");
            writeMatrixToFile(arg["lbg"], basePath + "args_lbg.txt");
            writeMatrixToFile(arg["ubg"], basePath + "args_ubg.txt");
            writeMatrixToFile(arg["x0"], basePath + "args_x0.txt");
            writeMatrixToFile(arg["p"], basePath + "args_p.txt");

            casadi::DM ubxdata = mpc_obj->upperboundx(p); 
            res = solver(arg);
            mpc_obj->setprevioussol(res.at("x"));
            mpc_obj->getOptForceCoeff(fit_order);

            
            mpc_obj->setpreviousp(p);
        }

        forceFFvec = mpc_obj->getOptforce(controlTick%10,fit_order);
        p_foot = mpc_obj->getFootPos();
        
        setExtforcewithvis(srbx2,0, list, "extForceArrow0", 0, bodypos_e, bodyrot_e, p_foot.block(0,0,3,1), forceFFvec.block(0,0,3,1));
        setExtforcewithvis(srbx2,0, list, "extForceArrow1", 0, bodypos_e, bodyrot_e, p_foot.block(3,0,3,1), forceFFvec.block(3,0,3,1));
        setExtforcewithvis(srbx2,0, list, "extForceArrow2", 0, bodypos_e, bodyrot_e, p_foot.block(6,0,3,1), forceFFvec.block(6,0,3,1));
        setExtforcewithvis(srbx2,0, list, "extForceArrow3", 0, bodypos_e, bodyrot_e, p_foot.block(9,0,3,1), forceFFvec.block(9,0,3,1));
        if(controlTick>5000 && controlTick<5500){
            setExtforcewithvis(srbx2,0, list, "disturbForceArrow0", 0, bodypos_e, bodyrot_e, q0.block(0,0,3,1), distforce);
        }
           
        
        mpc_obj->mpcdataLog(q0,forceFFvec,p_foot,controlTick);
    }
    std:: cout<< "controlTick: " << controlTick << std::endl;
    controlTick++;
}

// ================================================= //
// ============= MainFunction Start ================ //
// ================================================= //

int main(int argc, char **argv) {

    //std::cout << std::setprecision(14);
    std::cout.setf(ios::fixed,ios::floatfield);
    std::cout.precision(14);
    raisim::World::setActivationKey(raisim::loadResource("../activation.raisim"));
    /// create raisim world
    raisim::World world;
    world.setTimeStep(simfreq_raisim);

    auto vis = raisim::OgreVis::get();
    auto& list = raisim::OgreVis::get()->getVisualObjectList();

    /// these method must be called before initApp
    vis->setWorld(&world);
    vis->setWindowSize(1800, 1200);
    vis->setImguiSetupCallback(imguiSetupCallback);
    vis->setImguiRenderCallback(imguiRenderCallBack);
    vis->setKeyboardCallback(raisimKeyboardCallback);
    vis->setSetUpCallback(setupCallback);
    vis->setAntiAliasing(2);

    /// starts visualizer thread
    vis->initApp();

    bool simpleground = true;
    
    if(simpleground == true){
        /// create raisim objects
        auto ground = world.addGround();
        ground->setName("checkerboard");
        /// create visualizer objects
        vis->createGraphicalObject(ground, 20, "floor", "checkerboard_green");
        
        // raisim::TerrainProperties terrainProperties;
        // terrainProperties.frequency = 0.0;
        // terrainProperties.zScale = 0.0;
        // terrainProperties.xSize = 300.0;
        // terrainProperties.ySize = 300.0;
        // terrainProperties.xSamples = 50;
        // terrainProperties.ySamples = 50;
        // terrainProperties.fractalOctaves = 0.0;
        // terrainProperties.fractalLacunarity = 0.0;
        // terrainProperties.fractalGain = 0.0;
        // raisim::HeightMap *ground = world.addHeightMap(0.0, 0.0, terrainProperties);

        // world.setDefaultMaterial(0.8, 0.0, 0.0); //surface friction could be 0.8 or 1.0
        // vis->createGraphicalObject(ground, "terrain", "checkerboard_blue");
    }

    if(simpleground == false){
        /// create raisim objects    
        raisim::TerrainProperties terrainProperties;
        terrainProperties.frequency = 0.0;
        terrainProperties.zScale = 0.0;
        terrainProperties.xSize = 50.0;
        terrainProperties.ySize = 50.0;
        terrainProperties.xSamples = 10;
        terrainProperties.ySamples = 10;
        terrainProperties.fractalOctaves = 0;
        terrainProperties.fractalLacunarity = 0.0;
        terrainProperties.fractalGain = 0.0;
        auto ground = world.addHeightMap(0.0, 0.0, terrainProperties);
        /// create visualizer objects
        vis->createGraphicalObject(ground, "terrain", "checkerboard_green");
    }

    world.setDefaultMaterial(0.9, 0.0, 0.0); //surface friction could be 0.8 or 1.0
    //world.setDefaultMaterial(1.0, 0.0, 0.0); //surface friction could be 0.8 or 1.0

    // ============================================= //
    // ======= External force Visualization ======== //
    // ============================================= //
    vis->addVisualObject("extForceArrow0", "arrowMesh", "blue", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("extForceArrow1", "arrowMesh", "blue", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("extForceArrow2", "arrowMesh", "blue", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("extForceArrow3", "arrowMesh", "blue", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("extForceArrow4", "arrowMesh", "blue", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("extForceArrow5", "arrowMesh", "blue", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("extForceArrow6", "arrowMesh", "blue", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("extForceArrow7", "arrowMesh", "blue", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("disturbForceArrow0", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("disturbForceArrow1", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);

    // ============================================= //
    // ============ SETUP ROBOT ==================== //
    // ============================================= //
    std::vector<raisim::ArticulatedSystem*> srbx2;

    Eigen::VectorXd jointPgain(6), jointDgain(6);
    jointPgain.setZero();
    jointDgain.setZero();
    jointPgain.tail(6).setConstant(0.0);
    jointDgain.tail(6).setConstant(0.0);

    // ============================================= //
    // === Number of agent and initial pos setup === //
    // ============================================= //
    const size_t totalagentnum = 1;
    bool constraint = true; // this sim is using constraint from raisim: it has some compliance which makes result less precise
    bool payload = false;
    bool collaboration = true;

    // =========================================== //
    Eigen::MatrixXd agentinitial(totalagentnum, 3); // x, y, yawangle
    agentinitial.setZero(totalagentnum, 3);

    // ========= one-step delay stack ============= //
    //especially for distributed control with communication delay
    bool communicationdelay = false;
    size_t delaytick = 1;
    Eigen::MatrixXd stateDelaystack(4*totalagentnum, delaytick+1);
    stateDelaystack.setZero(4*totalagentnum, delaytick+1);
    
    if(totalagentnum == 1){
        agentinitial(0,0) = 0.0;
        agentinitial(0,1) = 0.0;
        agentinitial(0,2) = 0.0; //agent0 initial yaw angle
    }

    if(totalagentnum == 2){
        agentinitial(0,0) = 0.0;
        agentinitial(0,1) = 0.0;
        agentinitial(0,2) = 0.0; //agent0 initial yaw angle
        
        agentinitial(1,0) = 0.0;//1/sqrt(2);//0.0;
        agentinitial(1,1) = -1.0;//1/sqrt(2);//1.0;
        agentinitial(1,2) = 0.0; //agent1 initial yaw angle

        for(size_t i=0; i<delaytick; i++){
            stateDelaystack(0,i)=agentinitial(0,0);
            stateDelaystack(2,i)=agentinitial(0,1);
            stateDelaystack(4,i)=agentinitial(1,0);
            stateDelaystack(6,i)=agentinitial(1,1);
        }
    }
    // =========================================== //


    // for(size_t i=0; i<totalagentnum; i++){
        int i=0;
        double qw, qx, qy, qz;
        qw = 1.0; qx = 0.0; qy = 0.0; qz = 0.0;
        toQuartonian(qw, qx, qy, qz, 0.0, 0.0, agentinitial(i,2));
        //qw = 1/sqrt(2); qx = 0.0; qy = -1/sqrt(2); qz = 0.0;
        //toQuartonian(qw, qx, qy, qz, 0.0, -1.5707963, agentinitial(i,2));

        std::cout << "qw:"<< qw <<"qx:"<< qx <<"qy:"<< qy <<"qz:"<< qz << std::endl;
        srbx2.push_back(world.addArticulatedSystem(raisim::loadResource("srb_a1_upright.urdf")));
        //srbx2.push_back(world.addArticulatedSystem(raisim::loadResource("togo/togo_manual_measure.urdf")));
        auto A1_vis = vis->createGraphicalObject(srbx2.back(), "SRB" + std::to_string(i));

        srbx2.back()->setGeneralizedCoordinate({agentinitial(i, 0), agentinitial(i,1), 0.35, qw, qx, qy, qz});

        srbx2.back()->setGeneralizedForce(Eigen::VectorXd::Zero(srbx2.back()->getDOF()));
        srbx2.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);//PD_PLUS_FEEDFORWARD_TORQUE);
        srbx2.back()->setPdGains(jointPgain, jointDgain);
        srbx2.back()->setName("SRB" + std::to_string(i));
    // }

    printAllBodyNames(vis);

    srbx2[0]->printOutBodyNamesInOrder();

    std::string cameraview;
    double scale;
    scale = 1;
    cameraview = "side";

    vis->select(A1_vis->at(0),false);
    //iso
    if(cameraview == "iso"){
        vis->getCameraMan()->getCamera()->setPosition(3*scale, 2.5, 3);
        vis->getCameraMan()->getCamera()->roll(Ogre::Radian(2.3));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(0.9));
    }
    //side
    else if(cameraview == "side"){
        vis->getCameraMan()->getCamera()->setPosition(0.5*scale, 4, 0.5);
        vis->getCameraMan()->getCamera()->roll(Ogre::Radian(3.14));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(1.57));
    }
    //sidetop
    else if(cameraview == "sidetop"){
        vis->getCameraMan()->getCamera()->setPosition(0.5*scale, 4, 3);
        vis->getCameraMan()->getCamera()->roll(Ogre::Radian(3.14));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(0.9));
    }
    //front
    else if(cameraview == "front"){
        vis->getCameraMan()->getCamera()->setPosition(4*scale, 0, 0.5);
        vis->getCameraMan()->getCamera()->roll(Ogre::Radian(1.57));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(1.57));
    }
    //top
    else if(cameraview == "top"){
        vis->getCameraMan()->getCamera()->setPosition(0.5*scale, 0, 4.*scale);
        vis->getCameraMan()->getCamera()->roll(Ogre::Radian(-1.57));
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(0.));
    }
    else{
        /// set camera
        vis->getCameraMan()->getCamera()->setPosition(1, -2-int(1), 1.5+1);
        vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(1.));
    }
    vis->getCameraMan()->setYawPitchDist(Ogre::Radian(-0.5*M_PI/5+M_PI), -Ogre::Radian(0.70*M_PI/2), 3);
    
    //-----Video stuff-------//
    bool record = true;             // Record?
    double startTime = 0*ctrlHz;    // Recording start time
    double simlength = 8*ctrlHz;   // Sim end time
    double fps = 30;            
    std::string directory = "/home/taizoon/raisimEnv/raisimWorkspace/foot_videos/";
    std::string filename = "upright_walk_5r";
    const std::string name = directory+filename+"_"+cameraview+".mp4";
    vis->setDesiredFPS(fps);
    
    //--------------end-----------//
    /// run the app
    //vis->run();

    //for Unity
    //raisim::RaisimServer server(&world);
    //server.launchServer();
    SRBNMPCR* mpc_obj = new SRBNMPCR(argc,argv,1,0);   
    mpc_obj->generator();

    std::string file_name = "upright_nlp_9r";

    // code predix
    std::string prefix_code = fs::current_path().string() + "/";
    // shared library prefix
    std::string prefix_lib = fs::current_path().string() + "/";

    // Create a new NLP solver instance from the compiled code
    std::string lib_name = prefix_lib + file_name + ".so";
    casadi::Dict opts = {{"ipopt.print_level", 0}, {"print_time", 0},{"ipopt.max_iter", 200},{"ipopt.acceptable_tol", 1e-2},{"ipopt.acceptable_obj_change_tol", 1e-2}};
    casadi::Function solver = casadi::nlpsol("solver", "ipopt", lib_name, opts);

    /// run the app using while loop
    signal(SIGINT, signalhandler);
    long simcounter = 0;
    casadi::DM contact_sequence_dm = casadi::DM::ones(4,40);
    contact_sequence_dm(1,casadi::Slice(5,20)) = casadi::DM::zeros(1,15);
    contact_sequence_dm(2,casadi::Slice(8,20)) = casadi::DM::zeros(1,12);
    contact_sequence_dm(0,casadi::Slice(25,40)) = casadi::DM::zeros(1,15);
    contact_sequence_dm(3,casadi::Slice(28,40)) = casadi::DM::zeros(1,12);

    while (!vis->getRoot()->endRenderingQueued() && simcounter <  3e7){

        planner(srbx2, 1, 0,mpc_obj,solver,contact_sequence_dm);// loco_pln0, 
        world.integrate();
        if (simcounter%2 == 0){
            vis->renderOneFrame();
        }

        if (!vis->isRecording() & record & simcounter>=startTime){
            vis->startRecordingVideo(name);
        }

        //for Unity
        //raisim::MSLEEP(2);
        //server.integrateWorldThreadSafe();
        
        simcounter++;
        if(signaled)
            break;
    }
    /// terminate
    if (vis->isRecording())
        vis->stopRecordingVideoAndSave();

    /// terminate the app
    vis->closeApp();
    //for Unity
    //server.killServer();

    //delete loco_pln0;

    //std::cout << "delete done" <<std::endl;

    return 0;
}