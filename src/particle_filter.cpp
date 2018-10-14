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
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define EPS 0.0001

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	std::default_random_engine gen;

	if (is_initialized) {
		return;
	}
	num_particles = 100;

	std::normal_distribution<double> noise_x(x, std[0]);
	std::normal_distribution<double> noise_y(y, std[1]);
	std::normal_distribution<double> noise_theta(theta, std[2]);

	for (int i = 0; i < num_particles; i++) {
		Particle particle;
		particle.id = i;
		particle.x = noise_x(gen);
		particle.y = noise_y(gen);
		particle.theta = noise_theta(gen);
		particle.weight = 1.0;

		particles.push_back(particle);
		weights.push_back(1);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	std::default_random_engine gen;
	double new_x, new_y, new_theta;
	

	for (int i = 0; i < num_particles; i++) {
		Particle p = particles[i];
		if (fabs(yaw_rate) < EPS) {
			new_x = particles[i].x + velocity * cos(particles[i].theta) * delta_t;
			new_y = particles[i].y + velocity * sin(particles[i].theta) * delta_t;
			new_theta = particles[i].theta;
		}
		else {
			new_x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			new_y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			new_theta = particles[i].theta +  yaw_rate * delta_t;
		}

		std::normal_distribution<double> dist_x(new_x, std_pos[0]);
		std::normal_distribution<double> dist_y(new_y, std_pos[1]);
		std::normal_distribution<double> dist_theta(new_theta, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = disty(gen);
		particles[i].theta = dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
		for(unsigned int i = 0; i < observations.size(); i++) {
		unsigned int nObs = observations.size();
		unsigned int nPred = predicted.size();
		for(unsigned int i = 0; i < nObs; i++) { // For each observation
			double minDist = numeric_limits<double>::max();
			int mapId = -1;
			for(unsigned j = 0; j < nPred; j++ ) { // For each predition.
				double xDist = observations[i].x - predicted[j].x;
				double yDist = observations[i].y - predicted[j].y;
				double distance = xDist * xDist + yDist * yDist;
				if(distance < minDist) {
					minDist = distance;
					mapId = predicted[j].id;
				}
				observations[i].id = mapId;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int k = 0; k < num_particles; k++) {
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;

		std::vector<LandmarkObs> trans_observations;
		LandmarkObs obs;
		LandmarkObs t_obs;
		Particle particle = particles[k];

		for (unsigned int j = 0; j < observations.size(); j++) {
			obs = observations[j];

			t_obs.x = particle.x + (obs.x*cos(particle.theta) - obs.y*sin(particle.theta));
			t_obs.y = particle.y + (obs.x*sin(particle.theta) - obs.y*cos(particle.theta));
			trans_observations.push_back(t_obs);
		}

		particles[k].weight = 1.0;
		double closest_value = sensor_range*sensor_range;

		for(unsigned int i = 0; i < trans_observations.size(); i++) {
			LandmarkObs obs = trans_observations[i];
			int minIndex = 0;

			for (unsigned int m = 0; m < map_landmarks.landmark_list.size(); m++) {
				double landmark_x = map_landmarks.landmark_list[m].x_f;
				double landmark_y = map_landmarks.landmark_list[m].y_f;

				double calc_value = pow(obs.x - landmark_x, 2.0) + pow(obs.y - landmark_y, 2.0);
				if (calc_value < closest_value) {
					closest_value = calc_value;
					minIndex = m;
				}
			}

			double x = obs.x;
			double y = obs.y;
			double mu_x = map_landmarks.landmark_list[minIndex].x_f;
			double mu_y = map_landmarks.landmark_list[minIndex].y_f;
			double std_x = std_landmark[0];
			double std_y = std_landmark[1];

			double mult_value = 1/(2*M_PI*std_x*std_y)*exp(-0.5*(pow((x-mu_x)/std_x, 2) + pow((y-mu_y)/std_y, 2)));
			if (mult_value > 0) {
				particles[k].weight *= mult_value;
			}
			associations.push_back(minIndex+1);
			sense_x.push_back(x);
			sense_y.push_back(y);
		}

		weights[k] = particles[k].weight;
		particles[k] = SetAssociations(particles[k],associations,sense_x, sense_y);
		
	}
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;

	vector<double> weights;
	for(int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);

	}

	// index for resampling wheel
	uniform_int_distribution<int> int_dist(0, num_particles-1);
	auto index = int_dist(gen);

	double max_weight = *max_element(weights.begin(), weights.end());
    uniform_real_distribution<double> real_dist(0.0, max_weight);

	vector<Particle> resampled_particles;

	double beta = 0.0;
	for(int i = 0; i < num_particles; i++) {
		beta += real_dist(gen)*2.0;
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		resampled_particles.push_back(particles[index]);
	}

	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
	const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;

	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
