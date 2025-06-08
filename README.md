# Design-of-Whisker-Inspired-Sensors-for-Multi-Directional-Hydrodynamic-Sensing
## Overview
This repository is for the research paper "Design of Whisker-Inspired Sensors for Multi-Directional Hydrodynamic Sensing". The paper introduces a novel aquatic whisker-inspired flow sensor that accurately predicts both flow speed and direction. The sensor is modular, waterproof, and highly customizable where design changes have significant effects on specifications like sensitivity and range. Its capabilities are demonstrated through analytical models, experimental setups, and application tests. In this repository, we provide a MATLAB application based on our analytical model. The user interfce allows researchers to change the sensor design parameters and see how these changes affect the sensor's specifications. We hope this tool will enable future researchers to apply the work presented in this paper to new flow sensing tasks.

Paper Citation :
[1] Wang, T., Kent, T. A., & Bergbreiter, S. (2023). Design of whisker-inspired sensors for multi-directional hydrodynamic sensing. arXiv preprint arXiv:2307.09569.

### Authors
Teresa Kent Carnegie Mellon University Robotics Institute
Tuo Wang Carnegie Mellon Universtiy Department of Mechanical Engineering
Sarah Bergbreiter Carnegie Mellon University, Department of Mechanical Engineering

The authors would like to aknowledge Suhan Kim, a former member student in the Micro Robotics Lab at Carnegie Mellon University whose work and code developed for his papers [2] and [3] informed the development of this application.

[2] Kim, S., Kubicek, R., Paris, A., Tagliabue, A., How, J. P., & Bergbreiter, S. (2020, October). A whisker-inspired fin sensor for multi-directional airflow sensing. In 2020 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) (pp. 1330-1337). IEEE.
[3] Kim, S., Velez, C., Patel, D. K., & Bergbreiter, S. (2019, November). A magnetically transduced whisker for angular displacement and moment sensing. In 2019 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) (pp. 665-671). IEEE.

## Application for the design of Hydrodynamic Whisker Sensors

### Application Install
Download file folder in order to run the application. All files inside the folder contain functions needed to run the application.

### Sensor Design and Manufacturing
Steps to manufacture the sensor:
![Manufacturing (1)](https://github.com/user-attachments/assets/c4b360cf-5f03-4046-8816-6fb0b22d9949)
1. Make the sensor body using
  1. A water jet cutter on carbon fiber sheets [1]
  2. A CO2 laser cutter on
  3. A 3D printer as in [4]
2. 3D print the housing components
  1. a 2 mm x 2mm x 1 mm platform
  2. the housing for the Hall Effect Sensor
  3. the housing for the spring
3. Laser cut steel springs using a UV Laser 
4. Wire the Hall Effect sensor and cover the board with Silpoxy
5. Assembele the sensor
  1. Hall Effect Sensor in the housing
  2. Attach by glue and place in the housing: Whisker body onto the 2 mm platform onto the spring, magnet on the opposite side of the spring
  3. Attach the two housing components together and cover in SiylPoxy

Soon we will upload much of the sensor design CAD so others can alter the design to fit their application. Please reach out to the authors if you would like the files before that time.

