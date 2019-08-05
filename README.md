# Overview
This repository contains all the code needed to complete the final project for the Localization course in Udacity's Self-Driving Car Nanodegree.
# Particle filter implementation steps
all implemented codes can be found in `particle_filter.cpp`
1. initialization - using GPS input `ParticleFilter::init`
The most practical way to initialize our particles and generate real time output, is to make an initial estimate using GPS input. As with all sensor based operations, this step is impacted by noise.
2. Prediction - add control input (yaw rate & velocity) `ParticleFilter::prediction`
3. Update- update particle weights using map landmark positions and feature measurements `ParticleFilter::updateWeights`
	* Transform car coordinate to map coordinate
	* Data Associations -Associate  the nearest matched observed coordinate with map landmarks `ParticleFilter::dataAssociation`
	* Calculate the particle's final weight		
4.  Resample `ParticleFilter::resample()`

[//]: # (Image References)
[image1]: results_screenshot.PNG  "result screenshot"

## Running the Code
This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository includes two files that can be used to set up and install uWebSocketIO for either Linux or Mac systems. For windows you can use either Docker, VMware, or even Windows 10 Bash on Ubuntu to install uWebSocketIO.

Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./particle_filter

Alternatively some scripts have been included to streamline this process, these can be leveraged by executing the following in the top directory of the project:

1. ./clean.sh
2. ./build.sh
3. ./run.sh
## Results
The particle filter successfully passes the current grading code in the simulator.
1. **Accuracy**: error: x: .116, y: .108, yaw: .004
2. **Performance**:  69.28 sec
![alt text][image1]
