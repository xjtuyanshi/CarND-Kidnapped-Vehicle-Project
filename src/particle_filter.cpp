/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::discrete_distribution;
std::default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	/**
	 * TODO: Set the number of particles. Initialize all particles to
	 *   first position (based on estimates of x, y, theta and their uncertainties
	 *   from GPS) and all weights to 1.
	 * TODO: Add random Gaussian noise Gaussian each particle.
	 * NOTE: Consult particle_filter.h for more information about this method
	 *   (and others in this file).
	 */
	 // creates a Gaussian distribution for x,y and theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	num_particles = 100;  // Set the number of particles(not sure how many need so set it to 100 for now)
	for (int i = 0; i < num_particles; i++) {
		Particle par;
		par.id = i;
		par.x = dist_x(gen);
		par.y = dist_y(gen);
		par.theta = dist_theta(gen);
		par.weight = 1;
		particles.push_back(par);
	}
	is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[],
	double velocity, double yaw_rate) {
	/**
	 * TODO: Add measurements to each particle and add random Gaussian noise.
	 * NOTE: When adding noise you may find std::normal_distribution
	 *   and std::default_random_engine useful.
	 *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	 *  http://www.cplusplus.com/reference/random/default_random_engine/
	 */
	for (int i = 0; i < num_particles; i++) {
		double theta_particle = particles[i].theta;
		//update theta
		double theta_t = theta_particle + yaw_rate * delta_t;
		normal_distribution<double> pred_dist_theta(theta_t, std_pos[2]);
		particles[i].theta = pred_dist_theta(gen);

		//get x
		double x_t = particles[i].x + velocity / yaw_rate * (sin(theta_t) - sin(theta_particle));
		normal_distribution<double> pred_dist_x(x_t, std_pos[0]);
		particles[i].x = pred_dist_x(gen);

		//get y
		double y_t = particles[i].x + velocity / yaw_rate * (cos(theta_particle) - cos(theta_t));
		normal_distribution<double> pred_dist_y(y_t, std_pos[1]);
		particles[i].y = pred_dist_y(gen);

	}


}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
	vector<LandmarkObs>& observations) {

	// Using a nearest - neighbors data association
	for (LandmarkObs& obs : observations) {
		int landmark_id = -1;
		//intial distance assume it's extremely large
		double min_dist = std::numeric_limits<double>::max();
		for (LandmarkObs& pre : predicted) {
			double distance = dist(obs.x, obs.y, pre.x, pre.y);
			//iterate the distance to get the min distance
			if (distance < min_dist) {
				landmark_id = pre.id;
				min_dist = distance;
			}
		}
		// get this obs' updated associated landmark id
		obs.id = landmark_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	const vector<LandmarkObs>& observations,
	const Map& map_landmarks) {
	/**
	 * TODO: Update the weights of each particle using a mult-variate Gaussian
	 *   distribution. You can read more about this distribution here:
	 *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	 * NOTE: The observations are given in the VEHICLE'S coordinate system.
	 *   Your particles are located according to the MAP'S coordinate system.
	 *   You will need to transform between the two systems. Keep in mind that
	 *   this transformation requires both rotation AND translation (but no scaling).
	 *   The following is a good resource for the theory:
	 *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	 *   and the following is a good resource for the actual equation to implement
	 *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
	 */

	 // my plan, observation in range (50m) are considered?
	 //for each particular - get the tranformed observations to landmark using fomular
	 //then for each particle,for each obs
		 // weight for each obs = multiv_prob(std_landmark[0],std_landmark[1],tran_x_obs,tran_y_obs,	
	 //map_landmarks.landmark_list[obs.id].x,map_landmarks.landmark_list[obs.id].y
	 // weight_obs *=weight_obs;
	 //p.weight = weight_obs;
	for (int i = 0; i < num_particles; i++) {
		double par_x = particles[i].x;
		double par_y = particles[i].y;
		double par_theta = particles[i].theta;

		//get the landmarks in a particular' sense range
		vector<LandmarkObs> landmark_in_range;
		for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			int landmark_id = map_landmarks.landmark_list[j].id_i;
			float landmark_x = map_landmarks.landmark_list[j].x_f;
			float landmark_y = map_landmarks.landmark_list[j].y_f;
			if (dist(par_x, par_y, landmark_x, landmark_y) < sensor_range) {
				LandmarkObs landmark = { landmark_id, landmark_x, landmark_y };
				landmark_in_range.push_back(landmark);
			}
		}
		// convert observations car coordinates to map coordinates
		vector<LandmarkObs> converted_obs;
		for (int k = 0; k < observations.size(); k++) {
			double x_m = par_x + cos(par_theta) * observations[k].x - sin(par_theta) * observations[k].y;
			double y_m = par_y + sin(par_theta) * observations[k].x + cos(par_theta) * observations[k].y;
			LandmarkObs obs = { observations[k].id, x_m, y_m };
			converted_obs.push_back(obs);
		}
		// data association 
		dataAssociation(landmark_in_range, converted_obs);

		//calculate weights
		for (LandmarkObs& obs : converted_obs) {
			// get the landmark which has same id with observation
			double weight = 1.0;
			for (LandmarkObs& landmark : landmark_in_range) {
				if (landmark.id == obs.id) {
					weight = multiv_prob(std_landmark[0], std_landmark[1], obs.x, obs.y, landmark.x, landmark.y);
					std::cout << obs.x << "," << obs.y << "<<"<<landmark.x<<","<<landmark.y<<std::endl;
					//FIX ME it seems the differrence between obs.y and landmark.y is quite large, consistently off by 4
					//obs.y is large than  landmark.y 
					break;
				}
			}
			if (weight == 0) {
				particles[i].weight *= 0.000001;
			}
			else {
				particles[i].weight *= weight;
			}
		}


	}
	std::cout << "Done" ;

}

void ParticleFilter::resample() {
	/**
	 * TODO: Resample particles with replacement with probability proportional
	 *   to their weigresampleht.
	 * NOTE: You may find std::discrete_distribution helpful here.
	 *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	 */
	 // get max of weight
	vector<double>weights;
	double max_weight = 0;
	for (Particle par : particles) {
		weights.push_back(par.weight);
		if (par.weight > max_weight) {
			max_weight = par.weight;
		}
	}
	//random distribution for index and beta
	discrete_distribution<int> index_dist(0, num_particles);
	std::uniform_real_distribution<double> beta_dist(0.0, 1.0);

	int index = index_dist(gen);
	vector<Particle>resampled_particles;

	// use wheel resample method
	double beta = 0.0;
	for (int i = 0; i < num_particles; i++) {
		beta += beta_dist(gen) * 2.0 * max_weight;
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		resampled_particles.push_back(particles[index]);
	}
	particles = resampled_particles;

}

void ParticleFilter::SetAssociations(Particle& particle,
	const vector<int>& associations,
	const vector<double>& sense_x,
	const vector<double>& sense_y) {
	// particle: the particle to which assign each listed association, 
	//   and association's (x,y) world coordinates mapping
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates
	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
	vector<int> v = best.associations;
	std::stringstream ss;
	copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
	vector<double> v;

	if (coord == "X") {
		v = best.sense_x;
	}
	else {
		v = best.sense_y;
	}

	std::stringstream ss;
	copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
