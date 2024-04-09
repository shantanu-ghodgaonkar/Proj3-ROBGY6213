# Vision and IMU Fusion with Unscented Kalman Filter
## Description
This project is part of the ROB-6213 course and involves developing an Unscented Kalman Filter (UKF) to fuse inertial data (from an IMU) and vision-based pose and velocity estimation. The goal is to improve the overall estimation accuracy compared to using the IMU or vision data alone.

The project is divided into two parts:

Part 1: Using the visual pose estimation as the measurement in the UKF.
Part 2: Using only the velocity from the optical flow as the measurement in the UKF.
Dataset
The project uses the same sensor data format and type as the previous projects (Project 1 and Project 2). The datasets to be used are dataset1 and dataset4.

The sensor data is provided in the following formats:

IMU data: A struct with fields is_ready, t, omg, and acc.
Vision data: Matrices time, position, angle, linearVel, and angVel.
