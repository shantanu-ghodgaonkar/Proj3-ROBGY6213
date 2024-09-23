# Vision and IMU Fusion with Unscented Kalman Filter

## Project Overview

This project investigates the application of an **Unscented Kalman Filter (UKF)** for sensor fusion, combining data from an **Inertial Measurement Unit (IMU)** and a vision-based system for robot state estimation. The UKF is employed due to its ability to handle nonlinearities, which potentially improves accuracy compared to traditional Kalman Filters.

The project evaluates two scenarios:
1. **Part 1**: Using visual pose estimation (position and orientation) as measurements.
2. **Part 2**: Using optical flow-derived velocity as measurements.

The UKF framework is developed and tested using datasets, with performance measured by comparing estimated trajectories to actual sensor data.

## Key Features

- **Sensor Fusion**: Combines IMU and vision-based pose/velocity estimation.
- **Unscented Kalman Filter (UKF)**: Handles nonlinearities in the system model.
- **Pose and Velocity Estimation**: Tracks position, orientation, and velocity of the robot.
- **Nonlinear System Handling**: Uses sigma points to propagate state through a nonlinear process model.

## Methodology

### Part 1: UKF with Position and Orientation Measurements

In this part, the UKF estimates the robot’s position and orientation using data from the camera. The algorithm propagates the state estimate using IMU measurements and updates it with visual measurements for correction.

### Part 2: UKF with Velocity Measurements

Here, the UKF relies on velocity estimates derived from optical flow to update the state. Relying solely on velocity introduces trade-offs in position accuracy but demonstrates the UKF’s robustness with limited measurement data.

## Results

### Part 1: UKF with Position, Orientation, and Velocity Measurements

<p align="center">
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part1/dataset1/Position%20X.png" alt="Part 1: Position X Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part1/dataset1/Position%20Y.png" alt="Part 1: Position Y Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part1/dataset1/Position%20Z.png" alt="Part 1: Position Z Estimate" width="300"/>
</p>

<p align="center">
  <b>Fig 1:</b> Part 1 - Estimated Position X, Y, Z
</p>

<p align="center">
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part1/dataset1/Orientation%20X.png" alt="Part 1: Orientation X Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part1/dataset1/Orientation%20Y.png" alt="Part 1: Orientation Y Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part1/dataset1/Orientation%20Z.png" alt="Part 1: Orientation Z Estimate" width="300"/>
</p>

<p align="center">
  <b>Fig 2:</b> Part 1 - Estimated Orientation X, Y, Z
</p>

<p align="center">
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part1/dataset1/Velocity%20X.png" alt="Part 1: Velocity X Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part1/dataset1/Velocity%20Y.png" alt="Part 1: Velocity Y Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part1/dataset1/Velocity%20Z.png" alt="Part 1: Velocity Z Estimate" width="300"/>
</p>

<p align="center">
  <b>Fig 3:</b> Part 1 - Estimated Velocity X, Y, Z
</p>

### Part 2: UKF with Position, Orientation, and Velocity Measurements

<p align="center">
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part2/dataset1/Position%20X.png" alt="Part 2: Position X Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part2/dataset1/Position%20Y.png" alt="Part 2: Position Y Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part2/dataset1/Position%20Z.png" alt="Part 2: Position Z Estimate" width="300"/>
</p>

<p align="center">
  <b>Fig 4:</b> Part 2 - Estimated Position X, Y, Z
</p>

<p align="center">
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part2/dataset1/Orientation%20X.png" alt="Part 2: Orientation X Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part2/dataset1/Orientation%20Y.png" alt="Part 2: Orientation Y Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part2/dataset1/Orientation%20Z.png" alt="Part 2: Orientation Z Estimate" width="300"/>
</p>

<p align="center">
  <b>Fig 5:</b> Part 2 - Estimated Orientation X, Y, Z
</p>

<p align="center">
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part2/dataset1/Velocity%20X.png" alt="Part 2: Velocity X Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part2/dataset1/Velocity%20Y.png" alt="Part 2: Velocity Y Estimate" width="300"/>
  <img src="https://github.com/shantanu-ghodgaonkar/Proj3-ROBGY6213/blob/9ae138d5710731869b1e68067b401a620312cace/img/plots/part2/dataset1/Velocity%20Z.png" alt="Part 2: Velocity Z Estimate" width="300"/>
</p>

<p align="center">
  <b>Fig 6:</b> Part 2 - Estimated Velocity X, Y, Z
</p>

## Conclusion

This project successfully implements an **Unscented Kalman Filter (UKF)** for sensor fusion using IMU and vision-based systems. The UKF demonstrates its ability to handle nonlinearities, resulting in accurate state estimation in both configurations.

- **Part 1** effectively estimated the robot’s position and orientation using visual data.
- **Part 2** showed robust velocity estimation using optical flow, despite some trade-offs in position accuracy.

The UKF's performance highlights its potential for real-world applications where nonlinearities and noisy sensor data are common challenges.

## Future Enhancements

- **Extended Sensor Fusion**: Incorporate additional sensor modalities, such as GPS, to further improve state estimation.
- **Dynamic Environment Testing**: Perform testing in more complex environments to evaluate the robustness of the UKF.
- **Alternative Filters**: Explore alternative state estimation filters, such as the **Unscented Particle Filter (UPF)**, for improved performance.

## References

1. Timothy D. Barfoot. *State Estimation for Robotics*. Cambridge University Press, 2017. ISBN: 1107159393.
2. Ghodgaonkar, Shantanu. *Implementation of an Extended Kalman Filter (EKF) for State Estimation*. 2024. [Link](https://drive.google.com/file/d/1MTB0XjvszwhZJvR-bCVggfwa9PXGbHp4/view?usp=sharing).
3. Ghodgaonkar, Shantanu. *Vision-Based 3D Attitude Estimation Using AprilTags*. 2024. [Link](https://drive.google.com/file/d/1PdO5MHMplsqnmxrlCXx4kWT6H2ObG0Q3/view?usp=sharing).
4. MathWorks. Cholesky factorization. [Link](https://www.mathworks.com/help/matlab/ref/chol.html).
5. Bruno Siciliano et al. *Robotics: Modelling, Planning and Control*. Springer Publishing Company, 2010. ISBN: 1849966346.
6. Sebastian Thrun, Wolfram Burgard, and Dieter Fox. *Probabilistic Robotics*. The MIT Press, 2005. ISBN: 0262201623.
7. E.A. Wan and R. Van Der Merwe. “The unscented Kalman filter for nonlinear estimation”. In: *Proceedings of the IEEE 2000 Adaptive Systems for Signal Processing, Communications, and Control Symposium*. 2000. DOI: [10.1109/ASSPCC.2000.882463](https://doi.org/10.1109/ASSPCC.2000.882463).
