// #include <mip/mip_all.hpp>
// #include <chrono>
// #include <iostream>
// #include <fstream>
// #include <thread>

// // Global to store last received CompAccel
// mip::data_filter::CompAccel g_comp_accel;
// bool g_got_accel = false;

// void handle_accel(void*, mip::PacketRef, mip::Timestamp, const mip::data_filter::CompAccel& data) {
//     g_comp_accel = data;
//     g_got_accel = true;
// }

// int main() {
//     const char* port = "/dev/ttyACM0";   // <- Update this to your device port
//     const uint32_t baudrate = 115200;

//     mip::DeviceInterface device;
//     mip::platform::SerialConnection connection(port, baudrate);

//     if (!device.connect(&connection)) {
//         std::cerr << "Failed to connect to device." << std::endl;
//         return 1;
//     }

//     std::cout << "Connected to IMU" << std::endl;

//     // Ping the device
//     if (!mip::commands_base::ping(&device)) {
//         std::cerr << "Ping failed!" << std::endl;
//         return 1;
//     }

//     // Set the filter message format to stream compensated acceleration
//     mip::commands_filter::MessageFormat::Entry accel_entry{
//         mip::data_filter::DATA_COMPENSATED_ACCELERATION, 100  // field descriptor, decimation = 100Hz
//     };

//     mip::commands_filter::writeMessageFormat(&device, {accel_entry});
//     mip::commands_filter::saveMessageFormat(&device);

//     // Register callback
//     mip::Dispatcher dispatcher;
//     dispatcher.registerExtractor<mip::data_filter::CompAccel>(&handle_accel, nullptr);
//     device.setPacketCallback(&dispatcher);

//     // Open log file
//     std::ofstream logfile("mip_accel_log.csv");
//     logfile << "Time(s),AccelX,AccelY,AccelZ\n";

//     auto t_start = std::chrono::steady_clock::now();

//     for (int i = 0; i < 500; ++i) {
//         g_got_accel = false;
//         device.update();

//         if (g_got_accel) {
//             auto t_now = std::chrono::steady_clock::now();
//             double t_sec = std::chrono::duration<double>(t_now - t_start).count();

//             logfile << t_sec << ","
//                     << g_comp_accel.accel[0] << ","
//                     << g_comp_accel.accel[1] << ","
//                     << g_comp_accel.accel[2] << "\n";

//             std::cout << "Accel [m/s^2]: " << g_comp_accel.accel[0] << ", "
//                       << g_comp_accel.accel[1] << ", "
//                       << g_comp_accel.accel[2] << std::endl;
//         }

//         std::this_thread::sleep_for(std::chrono::milliseconds(10));
//     }

//     logfile.close();
//     std::cout << "Logging complete." << std::endl;

//     return 0;
// }


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mip/mip_all.h>
#include <microstrain/connections/serial/serial_port.h>
#include <mip/definitions/data_filter.h>
#include <mip/definitions/commands_3dm.h>

int main() {
    mip_interface device;

    // Open serial port
    if (!mip_interface_init(&device, "/dev/ttyACM0", 115200)) {
        fprintf(stderr, "Failed to open port\n");
        return 1;
    }

    printf("Connected to device.\n");

    // Set data format for CompAccel (0x82, 0x04) at 100 Hz
    mip_3dm_message_format_descriptor descriptors[1];
    descriptors[0].descriptor = MIP_FILTER_DATA_COMPENSATED_ACCEL;
    descriptors[0].decimation = 1;  // 100 Hz if base rate is 100 Hz

    if (mip_3dm_write_message_format(&device, MIP_FILTER_DATA_DESC_SET, 1, descriptors) != MIP_STATUS_OK) {
        fprintf(stderr, "Failed to configure message format\n");
        return 1;
    }

    // Enable data stream (resume)
    if (mip_3dm_resume(&device) != MIP_STATUS_OK) {
        fprintf(stderr, "Failed to resume data stream\n");
        return 1;
    }

    printf("Streaming CompAccel...\n");

    // Loop to read and extract CompAccel
    for (int i = 0; i < 500; ++i) {
        mip_packet packet;
        mip_filter_comp_accel data;

        if (mip_interface_receive_packet(&device, &packet, 500) == MIP_STATUS_OK) {
            if (packet.descriptor_set == MIP_FILTER_DATA_DESC_SET &&
                packet.payload[0] == MIP_FILTER_DATA_COMPENSATED_ACCEL) {
                
                if (mip_filter_comp_accel_from_buffer(&data, packet.payload, packet.payload_length, 1) == MIP_STATUS_OK) {
                    printf("Accel [m/s^2]: %.3f, %.3f, %.3f\n", data.accel[0], data.accel[1], data.accel[2]);
                }
            }
        }
    }

    mip_interface_close(&device);
    return 0;
}

