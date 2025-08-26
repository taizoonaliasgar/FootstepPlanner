//
// Author: Basit M. Imran.
// Date : 2024-04-17
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

// #include "unitree_legged_sdk/unitree_legged_sdk.h"
// #include "unitree_legged_sdk/unitree_joystick.h"

#include "timer.h"
#include "stdio.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <filesystem>

#include "Transforms.hpp"
#include "Filters.h"
// #include "OtherUtils.hpp"
// #include <yaml-cpp/yaml.h>
#include "helper.hpp"
//Planner
#include "LocoWrapper.hpp"
#include "SRBNMPC.hpp"
// #include "A1_Dynamics_full.h"
#include "Go2_fulldynamics.h"
#include "shared_structs_ex2.hpp"
//IMU Sensor
#include "mip/mip_all.hpp"
#include "microstrain/connections/serial/serial_connection.hpp"
#include "example_utils.hpp"
#include "fusedVelocityEstimate.hpp"

#include <chrono>
// #include <thread>
#include <memory>

#include <unitree/robot/channel/channel_publisher.hpp>
#include <unitree/robot/channel/channel_subscriber.hpp>
#include <unitree/idl/go2/LowState_.hpp>
#include <unitree/idl/go2/LowCmd_.hpp>
#include <unitree/common/time/time_tool.hpp>
#include <unitree/common/thread/thread.hpp>

#include <unitree/idl/ros2/String_.hpp>

#include <unitree/robot/b2/motion_switcher/motion_switcher_client.hpp>
#include "gamepad.hpp"


using namespace unitree::common;
using namespace unitree::robot;
using namespace unitree::robot::b2;

#define TOPIC_LOWCMD "rt/lowcmd"
#define TOPIC_LOWSTATE "rt/lowstate"
#define TOPIC_LIDAR "rt/utlidar/switch"

// using namespace UNITREE_LEGGED_SDK;
constexpr double PosStopF = (2.146E+9f);
constexpr double VelStopF = (16000.0f);

bool rough_terrain_en = 0;
bool obstacle_en = 0;

sharedData HLData;
sharedData LLData;
sharedData IMUData;

class ExternalComm {

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
	mip::data_sensor::CompOrientationMatrix sensor_comp_orientation_matrix;
	mip::data_sensor::CompQuaternion sensor_comp_quaternion;
	
	mip::data_filter::Timestamp filter_gps_time;
	mip::data_filter::Status filter_status;
	mip::data_filter::EulerAngles filter_euler_angles;
	mip::data_filter::CompAngularRate filter_comp_angular_rate;
	mip::data_filter::CompAccel filter_comp_accel;

	// Dispatch handlers
	mip::DispatchHandler sensor_data_handlers[6];
	mip::DispatchHandler filter_data_handlers[5];

	// State tracking
	bool filter_state_running = false;
	bool is_initialized = false;
	std::mutex update_mutex;
	std::ofstream csvFile;

	VelocityKalman3D velocity_filter;
	

public:
		// ExternalComm() : udpComp0(8082, "192.168.123.10", 8007, sizeof(LowCmd), sizeof(LowState)){
		ExternalComm() : velocity_filter(0.005, 0.0001, 0.01){//udpComp(LOWLEVEL), 

            fid = fopen("/home/taizoon/raisimEnv/Workspace/FootstepPlanner/stateData_1.csv", "w");
            
            const std::string port = "/dev/ttyACM0";//argv[2];
			const uint32_t baudrate = 9600;//std::stoi(argv[3]);
			
			std::cout << "Initializing IMU Sensor on port " << port << " at " << baudrate << " baud" << std::endl;
			
			// Create serial connection
			connection = std::make_shared<microstrain::connections::SerialConnection>(port, baudrate);
			
			// Try to connect
			if (!connection->connect()) {
				throw std::runtime_error("Failed to connect to " + port + " at " + std::to_string(baudrate) + " baud");
			}
			
			// Create device interface
			device = std::make_unique<mip::Interface>(
				connection.get(),       // Connection pointer
				parseBuffer.get(),      // Parse buffer
				parseBufferSize,        // Parse buffer size
				1000,                   // Parse timeout (ms)
				1000                    // Base reply timeout (ms)
			);
			
			std::cout << "Connected to device on " << port << " at " << baudrate << " baud" << std::endl;
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

			//4 Hz 
            double av[5] = { 1.000000000, -3.934325821, 5.805125421,-3.807232457, 0.936433243}; 
			double bv[5] = {0.000000024, 0.000000097, 0.000000145,0.000000097, 0.000000024}; 
            populate_filter_d(linearvelfilter, av, bv, 5, 3);

			csvFile.open("../data25/torqueEstimates.csv");
			
		}
	
		

	virtual ~ExternalComm(){
		clear_filter_d(jointfilter);
		clear_filter_f(remotefilter);
		clear_filter_f(angfilter);
		clear_filter_d(linearvelfilter);

		// nmpc_obj->logData();
        nmpc_obj.reset();
        loco_obj.reset();
		csvFile.close();
        
	}

	// main thread execution functions
	void Calc();
	void HighLevel();

    void Init();
	void InitLowCmd();
    void LowStateMessageHandler(const void* messages);
    void LowCmdWrite();
	int queryMotionStatus();
	void looper(ThreadPtr& thread, std::string name, float LLdt, int cpu_id, std::function<void()> func);
	std::string queryServiceName(std::string form,std::string name);

	void connectIMU();
    void setupIMUfilter();
    void configureIMU();
    void getIMMUdata();
	
	std::unique_ptr<LocoWrapper> loco_obj;
	std::unique_ptr<SRBNMPC> nmpc_obj;
    // std::unique_ptr<MovingAverageFilter> moving_avg_filter;
	// UDP udpComp;
    FILE *fid;

	long motiontime = 0;
	bool softFall = false;
	bool beginCommand = false;
	bool beginPose = false;
	bool setup = true;
	double q[18] = {0.0};
	double dq[18] = {0.0};
	double tauEst[12] = {0.0};
    double motorTemp[12] = {0.0};

	FiltStruct_d* jointfilter  = (FilterStructure_d*)malloc(sizeof(FilterStructure_d));
	FiltStruct_f* angfilter    = (FilterStructure_f*)malloc(sizeof(FilterStructure_f));
	FiltStruct_f* remotefilter = (FilterStructure_f*)malloc(sizeof(FilterStructure_f));
    FiltStruct_d* linearvelfilter = (FilterStructure_d*)malloc(sizeof(FilterStructure_d));

    // LowState state = {0};
	// LowCmd cmd = {0};
	xRockerBtnDataStruct remote;

    unitree_go::msg::dds_::LowCmd_ low_cmd{};      // default init
    unitree_go::msg::dds_::LowState_ low_state{};  // default init
    std_msgs::msg::dds_::String_ lidar_switch_msg;
    /*publisher subscriber*/
    ChannelPublisherPtr<unitree_go::msg::dds_::LowCmd_> lowcmd_publisher;
    ChannelSubscriberPtr<unitree_go::msg::dds_::LowState_> lowstate_subscriber;
    ChannelPublisherPtr<std_msgs::msg::dds_::String_> lidar_switch_pub;
    /*LowCmd write thread*/
    ThreadPtr imuThreadPtr;
    ThreadPtr mpcThreadPtr;
    ThreadPtr calcThreadPtr;
	MotionSwitcherClient msc;

	int stop = 0;
	float vel[3] = {0};
	double filtered_vel[3] = {0};
	float pose[6] = {0};
	float ang[2] = {0};
	Eigen::Vector3d scaled_acc = Eigen::MatrixXd::Zero(3,1);
	

	float LLdt = 0.001f;
 // Some parameters in controller0
	size_t loco_kind = STAND;                        // Gait pattern to use
    size_t settling = 0.2*ctrlHz;                   // Settling down
    size_t duration = 1.8*ctrlHz;                   // Stand up
    size_t loco_start = settling + duration;        // Start the locomotion pattern
    timer tset;
	double switchtime = 20;
	int robotdown = 1;
	
    Eigen::VectorXd jointTorqueFF = Eigen::MatrixXd::Zero(TOTAL_DOF,1);
	Eigen::VectorXd jointPosTotal = Eigen::MatrixXd::Zero(TOTAL_DOF+1,1); // +1 is for 4th Component of Quaternion 
    Eigen::VectorXd jointVelTotal = Eigen::MatrixXd::Zero(TOTAL_DOF,1);
    int conIndDes[4] = {1,1,1,1};

	double startTime = 0*ctrlHz;    // Recording start time
    double simlength = 50*ctrlHz;

    //Estimator
    int rearweight_est = 4;
    double yzdot_thresh = 0.3;
    double xdot_thresh = 0.3;
    double yzdot_thresh2 = 0.35;
    double xdot_thresh2 = 0.4;

	// Eigen::Matrix<double,3,1> eigen_eul = Eigen::MatrixXd::Zero(3,1);
	
	Eigen::Matrix<double,3,3> IMUframeoffset = (Eigen::Matrix3d() << 
																1, 0, 0,
																0, -1, 0,
																0, 0, -1).finished();

	Eigen::Matrix<double,3,3> IMUR = Eigen::MatrixXd::Identity(3,3);
	Eigen::Matrix<double,3,3> IMURotation = Eigen::MatrixXd::Identity(3,3);
	// Eigen::Matrix<double,3,1> imurot_eul = Eigen::MatrixXd::Zero(3,1);
	// Eigen::Matrix<double,3,1> imurot_eul2 = Eigen::MatrixXd::Zero(3,1);
	// Eigen::Matrix<double,3,1> imurot_eul3 = Eigen::MatrixXd::Zero(3,1);
	// Eigen::Matrix<double,4,1> imuquat = Eigen::MatrixXd::Zero(4,1);
	// Eigen::Matrix<double,3,1> imuquat_eul = Eigen::MatrixXd::Zero(3,1);

    //Estimator
	void getthetadot(double q[3],double dq[3]);
    void kinestimatorrr(double q[18], double dq[18], int contact[4], Eigen::Matrix<double,3,3> R);
    void getStateEstimatefullll(double q[18], double dq[18], int contact[4], Eigen::Matrix<double,3,3> R, Eigen::Matrix<double,3,4> toes, int robotdown, size_t ctrlTick, Eigen::Vector3d acc_wFrame);
  
};


uint32_t crc32_core(uint32_t* ptr, uint32_t len){

    unsigned int xbit = 0;
    unsigned int data = 0;
    unsigned int CRC32 = 0xFFFFFFFF;
    const unsigned int dwPolynomial = 0x04c11db7;

    for (unsigned int i = 0; i < len; i++)
    {
        xbit = 1 << 31;
        data = ptr[i];
        for (unsigned int bits = 0; bits < 32; bits++)
        {
            if (CRC32 & 0x80000000)
            {
                CRC32 <<= 1;
                CRC32 ^= dwPolynomial;
            }
            else
            {
                CRC32 <<= 1;
            }

            if (data & xbit)
                CRC32 ^= dwPolynomial;
            xbit >>= 1;
        }
    }

    return CRC32;
}

void ExternalComm::Init(){
    InitLowCmd();

    /*create publisher*/
    lowcmd_publisher.reset(new ChannelPublisher<unitree_go::msg::dds_::LowCmd_>(TOPIC_LOWCMD));
    lowcmd_publisher->InitChannel();

    /*create subscriber*/
    lowstate_subscriber.reset(new ChannelSubscriber<unitree_go::msg::dds_::LowState_>(TOPIC_LOWSTATE));
    lowstate_subscriber->InitChannel(std::bind(&ExternalComm::LowStateMessageHandler, this, std::placeholders::_1), 1);

    /*LIDAR publisher*/
    lidar_switch_pub.reset(new ChannelPublisher<std_msgs::msg::dds_::String_>(TOPIC_LIDAR));
    lidar_switch_pub->InitChannel();
    lidar_switch_msg.data() = "OFF";  // turn OFF LIDAR
    lidar_switch_pub->Write(lidar_switch_msg);

	/*init MotionSwitcherClient*/
    msc.SetTimeout(10.0f); 
    msc.Init();
    /*Shut down motion control-related service*/
    while(queryMotionStatus())
    {
        std::cout << "Try to deactivate the motion control-related service." << std::endl;
        int32_t ret = msc.ReleaseMode(); 
        if (ret == 0) {
            std::cout << "ReleaseMode succeeded." << std::endl;
        } else {
            std::cout << "ReleaseMode failed. Error code: " << ret << std::endl;
        }
        sleep(5);
    }

}

void ExternalComm::InitLowCmd(){
    low_cmd.head()[0] = 0xFE;
    low_cmd.head()[1] = 0xEF;
    low_cmd.level_flag() = 0xFF;
    low_cmd.gpio() = 0;

    for(int i=0; i<20; i++)
    {
        low_cmd.motor_cmd()[i].mode() = (0x01);   // motor switch to servo (PMSM) mode
        low_cmd.motor_cmd()[i].q() = (PosStopF);
        low_cmd.motor_cmd()[i].kp() = (0);
        low_cmd.motor_cmd()[i].dq() = (VelStopF);
        low_cmd.motor_cmd()[i].kd() = (0);
        low_cmd.motor_cmd()[i].tau() = (0);
    }
}

void ExternalComm::LowStateMessageHandler(const void* message){
    low_state = *(unitree_go::msg::dds_::LowState_*)message;
}

int ExternalComm::queryMotionStatus(){

    std::string robotForm,motionName;
    int motionStatus;
    int32_t ret = msc.CheckMode(robotForm,motionName);
    if (ret == 0) {
        std::cout << "CheckMode succeeded." << std::endl;
    } else {
        std::cout << "CheckMode failed. Error code: " << ret << std::endl;
    }
    if(motionName.empty())
    {
        std::cout << "The motion control-related service is deactivated." << std::endl;
        motionStatus = 0;
    }
    else
    {
        std::string serviceName = queryServiceName(robotForm,motionName);
        std::cout << "Service: "<< serviceName<< " is activate" << std::endl;
        motionStatus = 1;
    }
    return motionStatus;
}

std::string ExternalComm::queryServiceName(std::string form,std::string name){

    if(form == "0")
    {
        if(name == "normal" ) return "sport_mode"; 
        if(name == "ai" ) return "ai_sport"; 
        if(name == "advanced" ) return "advanced_sport"; 
    }
    else
    {
        if(name == "ai-w" ) return "wheeled_sport(go2W)"; 
        if(name == "normal-w" ) return "wheeled_sport(b2W)";
    }
    return "";
}

void ExternalComm::looper(ThreadPtr& thread, std::string name, float LLdt, int cpu_id, std::function<void()> func){
	/*loop publishing thread*/
    thread = CreateRecurrentThreadEx(name, cpu_id, 1000*LLdt, func);
	
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

void ExternalComm::getthetadot(double q[3],double dq[3]){

    Eigen::Matrix<double,3,3> A;
    double phi = q[1];
    double theta = q[2];
    
    A(0,0) = 1;     A(0,1) = sin(phi)*tan(theta);   A(0,2) = cos(phi)*tan(theta);
    A(1,0) = 0;     A(1,1) = cos(phi);              A(1,2) = -sin(phi);
    A(2,0) = 0;     A(2,1) = sin(phi)/cos(theta);   A(2,2) = cos(phi)/cos(theta);

    Eigen::Matrix<double,3,1> dq_temp = {dq[0],dq[1],dq[2]};
    Eigen::Matrix<double,3,1> thetadot = A*dq_temp;
    dq[0] = thetadot(0);
    dq[1] = thetadot(1);
    dq[2] = thetadot(2);
}

void ExternalComm::getStateEstimatefullll(double q[18], double dq[18], int contact[4], Eigen::Matrix<double,3,3> R, Eigen::Matrix<double,3,4> toes, int robotdown, size_t ctrlTick, Eigen::Vector3d acc_wFrame){
    
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
	
	numContact = (contact[0]+contact[1])*robotdown + rearweight_est*(contact[2]+contact[3]);
	
    if(!robotdown){

        for (int i = 3; i < 18; ++i){
		    COM_vel_e[0] -= (Jfr_toe_e[3*i+0]*contact[0]*robotdown + Jfl_toe_e[3*i+0]*contact[1]*robotdown + Jrr_toe_e[3*i+0]*contact[2]*rearweight_est + Jrl_toe_e[3*i+0]*contact[3]*rearweight_est)*dq[i];
	 	    COM_vel_e[1] -= (Jfr_toe_e[3*i+1]*contact[0]*robotdown + Jfl_toe_e[3*i+1]*contact[1]*robotdown + Jrr_toe_e[3*i+1]*contact[2]*rearweight_est + Jrl_toe_e[3*i+1]*contact[3]*rearweight_est)*dq[i];
	 	    COM_vel_e[2] -= (Jfr_toe_e[3*i+2]*contact[0]*robotdown + Jfl_toe_e[3*i+2]*contact[1]*robotdown + Jrr_toe_e[3*i+2]*contact[2]*rearweight_est + Jrl_toe_e[3*i+2]*contact[3]*rearweight_est)*dq[i];
        }
	    COM_vel_e[0] /= numContact;
	    COM_vel_e[1] /= numContact;
	    COM_vel_e[2] /= numContact;
		if(ctrlTick > 29999){
			velocity_filter.stepExp(acc_wFrame, COM_vel_e,ctrlTick);
		}

    }else{
        Eigen::Matrix<double,3,1> dq_temp = {dq[3],dq[4],dq[5]};
        toWorld(&dq[3],dq_temp,R);
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
    }

	// Set results
	q[0] = COM_e[0]; q[1] = COM_e[1]; q[2] = COM_e[2];
    if(ctrlTick<(switchtime+3)*ctrlHz){
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



void ExternalComm::HighLevel(){

    updateDataExp(GET_DATA, HL_DATA, &HLData);
    if(HLData.control_Tick > switchtime*1000+1999){//} && HLData.control_Tick%10==0){ // Settle down
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

        updateDataExp(SET_DATA, HL_DATA, &HLData);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        // std::cout << duration.count() << "\t" << "Full high level time" << std::endl;
    }
}

void ExternalComm::Calc(){

    updateDataExp(GET_DATA, LL_DATA, &LLData);
	robotdown = motiontime < switchtime*ctrlHz ? 1 : 0;

    // udpComp.Recv();
	// udpComp.GetRecv(state);
    
    // ===================================================== //
	// ============== Wireless Remote Stuff  =============== //
	// ===================================================== //
    memcpy(&remote, low_state.wireless_remote().data(), 40);
	// if((int)remote.btn.components.B!=0 && (int)remote.btn.components.L2!=0){
    if((int)remote.btn.components.R2!=0){
		stop = 1;
		low_cmd.reserve() = 0;
	}
	if((int)remote.btn.components.A!=0 && (int)remote.btn.components.L2!=0){
		softFall = true;
		// soft_fall_both = true;
	}
	if( ((int)remote.btn.components.A!=0) && ((int)remote.btn.components.L2 == 0) ){
		beginCommand = 1;
		low_cmd.reserve() = 1;
		// start_both = true;
	}
	// std::cout << "Setting up first command \n";
	vel[0] = 0.75f*remote.ly; // x vel
	vel[1] = -0.4f*remote.rx; // y vel
	ang[0] = 20.0f*3.14f/180.0f*remote.ry; // yaw vel
	ang[1] = -2.0f*remote.lx; // pitch pos
	discrete_butter_f(remotefilter,vel);
	discrete_butter_f(angfilter,ang);
	pose[4] = ang[0]; // set filtered pose
	vel[2] = ang[1];  // set filtered ang vel
	// std::cout << "Updating State \n";
	double *rotMat;
	Eigen::Matrix3d R;
	// ===================================================== //
	// ================= Update the State ================== //
	// ===================================================== //
	for (int i=0; i<12; ++i){
		q[i+6] = low_state.motor_state()[i].q();
		dq[i+6] = low_state.motor_state()[i].dq();
        tauEst[i] = low_state.motor_state()[i].tau_est();
		motorTemp[i] = low_state.motor_state()[i].temperature();
	}
	discrete_butter_d(jointfilter,&dq[6]);

	if(robotdown>0){
		q[3] = low_state.imu_state().rpy()[0]; q[4] = low_state.imu_state().rpy()[1]; q[5] = low_state.imu_state().rpy()[2];
		dq[3] = low_state.imu_state().gyroscope()[0]; dq[4] = low_state.imu_state().gyroscope()[1]; dq[5] = low_state.imu_state().gyroscope()[2];
		quat_to_R(low_state.imu_state().quaternion().data(),R);
	}else{	
		q[3] = LLData.att_euler[0]; q[4] = LLData.att_euler[1]; q[5] = LLData.att_euler[2];//+1.8;//+0.16);
		dq[3] = LLData.comp_angular_rate[0]; dq[4] = -LLData.comp_angular_rate[1]; dq[5] = -LLData.comp_angular_rate[2];
		// eigen_eul = {q[3],q[4],q[5]};
		R = IMURotation;//(eigen_eul,R);
	}

	// csvFile << motiontime << "," << state.imu.rpy[0] << "," << state.imu.rpy[1] << "," << state.imu.rpy[2] << ","
	// 		  	<< LLData.att_euler[0] << "," << LLData.att_euler[1] << "," << LLData.att_euler[2] << "\n";
	csvFile << motiontime << "," << tauEst[0] << "," << tauEst[1] << "," << tauEst[2] << ","
									<< tauEst[3] << "," << tauEst[4] << "," << tauEst[5] << ","
										<< tauEst[6] << "," << tauEst[7] << "," << tauEst[8] << ","
											<< tauEst[9] << "," << tauEst[10] << "," << tauEst[11] << ","
                                            << motorTemp[0] << "," << motorTemp[1] << "," << motorTemp[2] << ","
											<< motorTemp[3] << "," << motorTemp[4] << "," << motorTemp[5] << ","
											<< motorTemp[6] << "," << motorTemp[7] << "," << motorTemp[8] << ","
											<< motorTemp[9] << "," << motorTemp[10] << "," << motorTemp[11] << "\n";

	// quat_to_XYZ(state.imu.quaternion[0],state.imu.quaternion[1],state.imu.quaternion[2],state.imu.quaternion[3],
	//             q[3],q[4],q[5]);
	
	//Offset for Halo 
	// q[3]+= -3.0*MY_PI/180.0;
	// q[4] += 1.5*MY_PI/180.0;
	// q[5] += 0.0*MY_PI/180.0;
	/* //Offset for Willow
	q[3]+=  0*MY_PI/180.0;
	q[4] += 0*MY_PI/180.0;
	q[5] += 0*MY_PI/180.0;
	//dq[5] += 5.0*MY_PI/180.0; */

	rotMat = R.data();
	Eigen::Map< Eigen::Matrix<double, 3, 3> > rotE(rotMat, 3, 3);

	// Kinematic Estimation
	int footForce[4];
	footForce[0] = low_state.foot_force()[0]; footForce[1] = low_state.foot_force()[1];
	footForce[2] = low_state.foot_force()[2]; footForce[3] = low_state.foot_force()[3];
	static int contactIndex[4] = {1,1,1,1};
	// kinEst0(footForce,contactIndex,q,dq,R); // Defined in OtherUtils.hpp
	
	if(motiontime>2499){
        getStateEstimatefullll(q,dq,LLData.ind_LL,rotE,LLData.toePos,robotdown,motiontime,LLData.IMUacc);
        
    }else if(motiontime>0){
        kinestimatorrr(q,dq,LLData.ind_LL,rotE);
    }

    if(motiontime>=37000){
		q[2]=0.5;
		// dq[2]=0.0;
	}
	// std::cout << "[LL] Motiontime: " << motiontime << std::endl;
	// ===================================================== //
	// ============= Quad Initialization Time ============== //
	// ===================================================== //
	if ( (motiontime < settling) ){
		for(int i=0; i<12; ++i){
			low_cmd.motor_cmd()[i].dq() = 0.0;
			low_cmd.motor_cmd()[i].kp() = 0.0f;
			low_cmd.motor_cmd()[i].kd() = 0.0f;
			low_cmd.motor_cmd()[i].tau() = 0.0f;
		}
		
	}
	// ===================================================== //
	// ================== LL Controller ==================== //
	// ===================================================== //
	// std::cout << "Controller Executing \n";

    
    // std::cout<<"[Running LL]"<<std::endl;
    // Update the desired torques
    if(motiontime >= settling){ // Settle down
        // LLData.runMPC = 1;
		if (beginCommand && !softFall){
			if (setup){
				loco_obj->initStandVars(jointPosTotal.block(0,0,3,1),jointPosTotal(5),(int)duration);
				setup = false;
			}
			loco_obj->setIMUdata(LLData.att_euler);

			// double phaseVar;
			if(motiontime < loco_start){ // Start standing
				loco_obj->calcTau2(q,dq,rotMat,STAND,motiontime,LLData.solvetime);  
			}else{
				loco_obj->ExpWrapper(q,dq,rotMat,motiontime,LLData.solvetime,LLData.ind,LLData.comDes,LLData.fDes);
			}

			LLData.control_Tick = motiontime;
			memcpy(LLData.q,q,18*sizeof(double));
			memcpy(LLData.dq,dq,18*sizeof(double));

			LLData.toePos = loco_obj->getfootposition();
			LLData.QPforce = loco_obj->getpreviousQPforce();
			const int* ind_LL = loco_obj->getConDes();
			LLData.ind_LL[0] = ind_LL[0];
			LLData.ind_LL[1] = ind_LL[1];
			LLData.ind_LL[2] = ind_LL[2];
			LLData.ind_LL[3] = ind_LL[3];
	
			updateDataExp(SET_DATA, LL_DATA, &LLData);
		
			// discrete_butter_d(linearvelfilter,&LLData0.dq[0]);
			// moving_avg_filter->filter(&LLData.dq[0]);

			// Set the command
			for(int i = 0; i < 12; ++i){
				low_cmd.motor_cmd()[i].tau() = loco_obj->ll->tau[i + 6];
				low_cmd.motor_cmd()[i].q()  = loco_obj->ll->q(i + 6);
				low_cmd.motor_cmd()[i].dq() = loco_obj->ll->dq(i + 6);

				low_cmd.motor_cmd()[i].kp() = 10;
				low_cmd.motor_cmd()[i].kd() = 4;
				if(beginPose==true){
					for(int i=0; i<4; ++i){
						if(LLData.ind_LL[i]==0){
			//							cmd.motorCmd[i].tau = 0;
                            if(i<2){
							    low_cmd.motor_cmd()[3*i].kp() = 7;
							    low_cmd.motor_cmd()[3*i+1].kp() = 7;
							    low_cmd.motor_cmd()[3*i+2].kp() = 7;

							    low_cmd.motor_cmd()[3*i].kd() = 0.05;
							    low_cmd.motor_cmd()[3*i+1].kd() = 0.05;
							    low_cmd.motor_cmd()[3*i+2].kd() = 0.05;
                            }else{
                                low_cmd.motor_cmd()[3*i].kp() = 7;
							    low_cmd.motor_cmd()[3*i+1].kp() = 7;
							    low_cmd.motor_cmd()[3*i+2].kp() = 7;

							    low_cmd.motor_cmd()[3*i].kd() = 0.05;
							    low_cmd.motor_cmd()[3*i+1].kd() = 0.05;
							    low_cmd.motor_cmd()[3*i+2].kd() = 0.05;

                            }
                            // if(i>1){low_cmd.motor_cmd()[3*i].kp() = 20;}
						}
					}	
				}
			}
			
			// if(motiontime>switchtime*ctrlHz){
			// 	// cmd.motorCmd[6].q  = 0;
			// 	// cmd.motorCmd[9].q  = 0;
			// 	// cmd.motorCmd[6].dq = 0;
			// 	// cmd.motorCmd[9].dq = 0;
			// 	cmd.motorCmd[6].Kp = 20;
			// 	cmd.motorCmd[9].Kp = 20;
			// }

            // Saturate the command
			float hr_max = 30.0f, hr_min = -30.0f;
			float hp_max = 30.0f, hp_min = -30.0f;
			float kn_max = 33.0f, kn_min = -33.0f;
			for (int i = 0; i < 4; i++){
				low_cmd.motor_cmd()[3 * i + 0].tau() = (low_cmd.motor_cmd()[3 * i + 0].tau() > hr_max) ? hr_max : low_cmd.motor_cmd()[3 * i + 0].tau();
				low_cmd.motor_cmd()[3 * i + 0].tau() = (low_cmd.motor_cmd()[3 * i + 0].tau() < hr_min) ? hr_min : low_cmd.motor_cmd()[3 * i + 0].tau();
				low_cmd.motor_cmd()[3 * i + 1].tau() = (low_cmd.motor_cmd()[3 * i + 1].tau() > hp_max) ? hp_max : low_cmd.motor_cmd()[3 * i + 1].tau();
				low_cmd.motor_cmd()[3 * i + 1].tau() = (low_cmd.motor_cmd()[3 * i + 1].tau() < hp_min) ? hp_min : low_cmd.motor_cmd()[3 * i + 1].tau();
				low_cmd.motor_cmd()[3 * i + 2].tau() = (low_cmd.motor_cmd()[3 * i + 2].tau() > kn_max) ? kn_max : low_cmd.motor_cmd()[3 * i + 2].tau();
				low_cmd.motor_cmd()[3 * i + 2].tau() = (low_cmd.motor_cmd()[3 * i + 2].tau() < kn_min) ? kn_min : low_cmd.motor_cmd()[3 * i + 2].tau();
			}
        }
        else if(softFall){

			for (int i = 0; i < 12; ++i){
				low_cmd.motor_cmd()[i].tau() = 0.0;
				low_cmd.motor_cmd()[i].q()   = 0.0;
				low_cmd.motor_cmd()[i].dq()  = 0.0;

				low_cmd.motor_cmd()[i].kp() = 0;
				low_cmd.motor_cmd()[i].kd() = ( ( (i+1)%3 ) == 0 ) ? 6 : 3;
			}
		}
        else if ((!beginCommand) && (!softFall)){
			
			static bool printedA = false;
			if(!printedA){
				printf("\nPress A to continue\n");
				printedA = true;
			}
		}
    }

	motiontime += 1;
    // udpComp.SetSend(cmd);
    // udpComp.Send();
    low_cmd.crc() = crc32_core((uint32_t *)&low_cmd, (sizeof(unitree_go::msg::dds_::LowCmd_)>>2)-1);
    lowcmd_publisher->Write(low_cmd);

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

    std::array<mip::DescriptorRate, 6> sensor_descriptors = {{
        { mip::data_sensor::DATA_TIME_STAMP_GPS, sensor_decimation },
        { mip::data_sensor::DATA_ACCEL_SCALED,   sensor_decimation },
        { mip::data_sensor::DATA_GYRO_SCALED,    sensor_decimation },
		{mip::data_sensor::DATA_COMP_EULER_ANGLES, sensor_decimation},
		{mip::data_sensor::DATA_COMP_ORIENTATION_MATRIX, sensor_decimation},
		{mip::data_sensor::DATA_COMP_QUATERNION, sensor_decimation},
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
	device->registerExtractor(sensor_data_handlers[4], &sensor_comp_orientation_matrix);
	device->registerExtractor(sensor_data_handlers[5], &sensor_comp_quaternion);

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

        // IMUData.scaled_acc[0] = this->sensor_accel.scaled_accel[0];
		// IMUData.scaled_acc[1] = this->sensor_accel.scaled_accel[1];
		// IMUData.scaled_acc[2] = this->sensor_accel.scaled_accel[2];
		scaled_acc(0) = this->filter_comp_accel.accel[0];
		scaled_acc(1) = this->filter_comp_accel.accel[1];
		scaled_acc(2) = this->filter_comp_accel.accel[2];
		IMUData.comp_angular_rate[0] = this->filter_comp_angular_rate.gyro[0];
        IMUData.comp_angular_rate[1] = this->filter_comp_angular_rate.gyro[1];
        IMUData.comp_angular_rate[2] = this->filter_comp_angular_rate.gyro[2]; 
		for(size_t i = 0; i < 3; i++){
			for (size_t j = 0; j < 3; j++)
			{
				IMUR(i,j) = this->sensor_comp_orientation_matrix.m[3*i+j];
			}
		}

		// Eigen::Vector3d acc_world = IMUR*Eigen::Vector3d(IMUData.scaled_acc[0], IMUData.scaled_acc[1], IMUData.scaled_acc[2]);
		// Eigen::Vector3d acc_IMU0 = IMUR.transpose()*scaled_acc;//Eigen::Vector3d(scaled_acc[0], scaled_acc[1], scaled_acc[2]);
		IMUData.IMUacc = IMUframeoffset*IMUR.transpose()*scaled_acc;//acc_IMU0;
		IMUData.IMUacc(2) -= 9.81;

		// std::cout << motiontime << "\t" << IMUData.IMUacc(0) << "\t" << IMUData.IMUacc(1) << "\t" << IMUData.IMUacc(2) << "\t" 
		// 				<< scaled_acc(0) << "\t" << scaled_acc(1) << "\t" << scaled_acc(2) << std::endl; 

		IMURotation = IMUframeoffset*IMUR.transpose()*IMUframeoffset;
		// IMURotation = IMUR;//.transpose();
		// Eigen::Vector3d acc_world3 = IMURotation*Eigen::Vector3d(IMUData.scaled_acc[0], IMUData.scaled_acc[1], IMUData.scaled_acc[2]);
		// Eigen::Vector3d acc_world4 = IMURotation.transpose()*Eigen::Vector3d(IMUData.scaled_acc[0], IMUData.scaled_acc[1], IMUData.scaled_acc[2]);

		IMUData.att_euler[0] = -atan2(IMURotation(1,2),IMURotation(2,2));
		IMUData.att_euler[1] = asin(IMURotation(0,2));
		IMUData.att_euler[2] = -atan2(IMURotation(0,1),IMURotation(0,0));

		// omegabody[0] = IMUData.comp_angular_rate[0];
		// omegabody[1] = IMUData.comp_angular_rate[1];
		// omegabody[2] = IMUData.comp_angular_rate[2];
		// getthetadot(IMUData.att_euler,IMUData.comp_angular_rate);

		// std::cout << motiontime << "\t" << IMUData.att_euler[0] << "\t" << IMUData.att_euler[1] << "\t" << IMUData.att_euler[2] << "\t" 
		// 				<< IMUData.scaled_acc[0] << "\t" << IMUData.scaled_acc[1] << "\t" << IMUData.scaled_acc[2] << "\t" 
		// 				<< acc_world2(0) << "\t" << acc_world2(1) << "\t" << acc_world2(2) << "\t"  
		// 				<< omegabody[0] << "\t" << omegabody[1] << "\t" << omegabody[2] << "\t"
		// 				<< IMUData.comp_angular_rate[0] << "\t" << IMUData.comp_angular_rate[1] << "\t" << IMUData.comp_angular_rate[2] << "\t"
		// 				<< acc_world3(0) << "\t" << acc_world3(1) << "\t" << acc_world3(2) << "\t" 
		// 				<< acc_world4(0) << "\t" << acc_world4(1) << "\t" << acc_world4(2) << "\t"
		// 				<< acc_worldfinal(0) << "\t" << acc_worldfinal(1) << "\t" << acc_worldfinal(2) << std::endl;
		// IMUData.att_euler[0] = this->filter_euler_angles.roll; 
		// IMUData.att_euler[1] = this->filter_euler_angles.pitch; 
		// IMUData.att_euler[2] = this->filter_euler_angles.yaw;

		// std::cout << IMUData.att_euler[0] << "\t" << IMUData.att_euler[1] << "\t" << IMUData.att_euler[2] << std::endl;
		// csvFile << IMUData.att_euler[0] << "," << IMUData.att_euler[1] << "," << IMUData.att_euler[2] << ","
		// 			<< imurot_eul2(0) << "," << imurot_eul2(1) << "," << imurot_eul2(2) << ","
		// 				<< imurot_eul(0) << "," << imurot_eul(1) << "," << imurot_eul(2) << ","
		// 				//  << imuquat_eul3(0) << "," << imuquat_eul3(1) << "," << imuquat_eul3(2) << ","
		// 				 << imurot_eul3(0) << "," << imurot_eul3(1) << "," << imurot_eul3(2) << "\n";
 
        updateDataExp(SET_DATA, IMU_DATA, &IMUData);     
		// motiontime += 1;                     
    }
    
}



int main(int argc, char *argv[]) {

    // InitEnvironment();
    if (argc < 2){
        std::cout << "Usage: " << argv[0] << " networkInterface" << std::endl;
        exit(-1); 
    }

	ChannelFactory::Instance()->Init(0, argv[3]);
    ExternalComm extComm;

	extComm.loco_obj = std::unique_ptr<LocoWrapper>(new LocoWrapper(argc, argv));
    extComm.nmpc_obj  = std::unique_ptr<SRBNMPC>(new SRBNMPC(argc,argv,1,0));
    // extComm.moving_avg_filter = std::unique_ptr<MovingAverageFilter>(new MovingAverageFilter(1000,3));

    extComm.connectIMU();
    extComm.configureIMU();
    extComm.setupIMUfilter();

    extComm.Init();

    extComm.imuThreadPtr = CreateRecurrentThreadEx("imu_loop",   4, 1000,  &ExternalComm::getIMMUdata, &extComm);
    // sleep(1.0);
    extComm.mpcThreadPtr = CreateRecurrentThreadEx("mpc_loop",   2, 10000, &ExternalComm::HighLevel,   &extComm);
    // sleep(1.0);
    extComm.calcThreadPtr = CreateRecurrentThreadEx("calc_loop", 1, 1000,  &ExternalComm::Calc,        &extComm);
    // sleep(1.0);
	

    while ( (extComm.stop==0) ){
		sleep(0.1);
	};
	fclose(extComm.fid);
	return 0;
}

