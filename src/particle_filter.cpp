/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 

#include "particle_filter.h"
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
  num_particles = 100;

  // This line creates a normal (Gaussian) distribution for x.
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_psi(theta, std[2]);

  particles.reserve(num_particles);
  for (int i = 0; i < num_particles; ++i) {
    Particle new_particle = {
      i, // id
      dist_x(gen), // x
      dist_y(gen), // y
      dist_psi(gen), // theta
      1 // weight
    };
    particles.push_back(new_particle);
  }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.

  normal_distribution<double> N_x_noise(0, std_pos[0]);
  normal_distribution<double> N_y_noise(0, std_pos[1]);
  normal_distribution<double> N_theta_noise(0, std_pos[2]);

  for (auto& particle : particles) {
    if (std::abs(yaw_rate) > negligible) {
      particle.x += velocity*(sin(particle.theta + yaw_rate*delta_t) - sin(particle.theta)) / yaw_rate;
      particle.y += velocity*(cos(particle.theta) - cos(particle.theta + yaw_rate*delta_t)) / yaw_rate;
      particle.theta += yaw_rate*delta_t;
    }
    else {
      particle.x += velocity*cos(particle.theta)*delta_t;
      particle.y += velocity*sin(particle.theta)*delta_t;
    }

    // add noise
    particle.x += N_x_noise(gen);
    particle.y += N_y_noise(gen);
    particle.theta += N_theta_noise(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
  for (auto& predicted_landmark : predicted) {
    double min_dist = dist(predicted_landmark.x, predicted_landmark.y, observations[0].x, observations[0].y);
    int min_dist_id = observations[0].id;

    for (auto& landmark: observations) {
      double distance = dist(predicted_landmark.x, predicted_landmark.y, landmark.x, landmark.y);
      if (distance < min_dist) {
        min_dist = distance;
        min_dist_id = landmark.id;
      }
    }

    predicted_landmark.id = min_dist_id;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.

  for (int i = 0; i < particles.size(); ++i) {
    Particle particle = particles[i];
    vector<LandmarkObs> predicted;

    // transplate observation into map's coordinate system
    for (auto& observation : observations) {
      LandmarkObs translated_observation;

      translated_observation.x = observation.x*cos(particle.theta) - observation.y*sin(particle.theta) + particle.x;
      translated_observation.y = observation.x*sin(particle.theta) + observation.y*cos(particle.theta) + particle.y;

      predicted.push_back(translated_observation);
    }

    // choose landmarks which are inside sensor range
    vector<LandmarkObs> landmarks_in_sensor_range;
    for (auto& landmark: map_landmarks.landmark_list) {
      if ((std::abs(landmark.x_f - particle.x) < sensor_range) &&
        (std::abs(landmark.y_f - particle.y) < sensor_range)) {

        landmarks_in_sensor_range.push_back(LandmarkObs{
          landmark.id_i,
          landmark.x_f,
          landmark.y_f
        });
      }
    }

    // find correspondence between observations and landmarks
    dataAssociation(predicted, landmarks_in_sensor_range);

    // update particle weight
    for (auto& observation : predicted) {
      auto corresponding_landmark = map_landmarks.landmark_list[observation.id - 1];

      weights[i] *= norm_probability(observation.x, observation.y,
        corresponding_landmark.x_f, corresponding_landmark.y_f,
        std_landmark[0], std_landmark[1]);
    }

  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
