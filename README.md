# Vision and IMU Fusion with Unscented Kalman Filter

## Project Overview

This project investigates the application of an **Unscented Kalman Filter (UKF)** for sensor fusion, combining data from an Inertial Measurement Unit (IMU) and a vision-based system for robot state estimation. The UKF is employed due to its ability to handle nonlinearities, potentially improving accuracy compared to traditional Kalman Filters.

The project is part of the course **Robot Localization and Navigation (ROB-GY 6213)** and aims to enhance robot state estimation by leveraging UKF's deterministic approach to propagate sigma points through the nonlinear process model.

## Key Features

- **Sensor Fusion**: Combines IMU and vision-based pose/velocity estimation for improved state estimation.
- **Nonlinearity Handling**: Utilizes UKF to handle nonlinearities that arise in the system's process model.
- **Two Scenarios Evaluated**:
  - Visual pose estimation as a measurement.
  - Optical flow-derived velocity as a measurement.
  
## Methodology

- The UKF framework is built by integrating IMU data and vision-based pose and velocity estimates.
- Two main components are evaluated:
  - **Part 1**: UKF with position and orientation measurements from visual data.
  - **Part 2**: UKF with velocity measurements from optical flow data.

### Prediction Step

- **Sigma Points**: Computation of sigma points to propagate through the nonlinear process.
- **Nonlinear Function**: Propagates sigma points through the non-linear process model.
- **Covariance Matrix**: Computes the predicted mean and covariance matrix.

### Update Step

- **Part 1**: Updates the state estimate using position and orientation from visual data.
- **Part 2**: Updates the state estimate using velocity from optical flow.

## Results

- **Part 1**: UKF effectively estimates position, orientation, and sensor biases using visual data.
- **Part 2**: UKF with velocity measurements shows a greater discrepancy in position estimation compared to Part 1 but maintains proficiency in velocity estimation.

## Conclusion

This project successfully implements the UKF for state estimation in a Micro Aerial Vehicle (MAV) simulation, demonstrating its capability in handling nonlinearities and sensor fusion. Future work could focus on adding additional sensor data or exploring alternative estimation techniques for enhanced robustness.

## Repository Structure

- **UKF Framework**: MATLAB code for prediction and update steps for the UKF.
- **Datasets**: Includes sample datasets for testing the UKF implementation.
- **Plots**: Visual results showcasing position, orientation, velocity, and bias estimation.

## References

This project leverages material from:
- Timothy D. Barfoot's *State Estimation for Robotics*.
- Previous course projects involving Extended Kalman Filter (EKF) and vision-based attitude estimation.

## License

This project does not currently have an associated license.
