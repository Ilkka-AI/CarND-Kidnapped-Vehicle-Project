/*
 * particle_filter.cpp
 *
 Edited by Ilkka
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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of 
	// x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	
	default_random_engine gen;
	normal_distribution<double> dist_x(x,std[0]);
	normal_distribution<double> dist_y(y,std[1]);
	num_particles=100;	
	normal_distribution<double> dist_theta(theta, std[2]);
	vector<double> weights;
	
	
	for(int i=0;i<num_particles;i++){
	particles.push_back(Particle());	
	particles[i].x = dist_x(gen);
	particles[i].y = dist_y(gen);
	particles[i].id = i;
	particles[i].theta = dist_theta(gen);	
	particles[i].weight=1;
	
	}
	is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	default_random_engine gen;	
	double minimum_yaw_rate=0.0001;
	normal_distribution<double> dist_x(0,std_pos[0]);
	normal_distribution<double> dist_y(0,std_pos[1]);
	normal_distribution<double> dist_theta(0,std_pos[2]);
	
	// Predict particle states
	for(int i=0;i<num_particles;i++){
	
	// If going straight, yaw rate is zero, divide by zero problem emerges, avoid it
	if (fabs(yaw_rate) < minimum_yaw_rate) {  
      	particles[i].x+=velocity*cos(particles[i].theta)*delta_t;
      	particles[i].y+=velocity*sin(particles[i].theta)*delta_t;
	particles[i].theta+=dist_theta(gen);	
   	 } 
	// In nonzero yaw rate, the motion equations are 
	else {
	particles[i].x+=velocity/yaw_rate*( sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta) )+dist_x(gen);
	particles[i].y+=velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t) )+dist_y(gen);
	particles[i].theta+=yaw_rate*delta_t+dist_theta(gen);	
}	
}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Find the predicted measurement that is closest to each observed measurement and assign the 
	// observed measurement to this particular landmark.
	
	// Loop through each observation and find nearest map landmark
	for(int i=0; i<observations.size();i++){
		std::vector<double> distances(predicted.size());
		int min_distance_element;	
		double min_distance;	
		//find nearest landmark for observation j
		for(int j=0;j<predicted.size();j++){
			distances[j]=pow(predicted[j].x-observations[i].x,2)+pow((predicted[j].y-observations[i].y),2);
			if(j==0){
				min_distance=distances[0];
				min_distance_element=0;
}			else{
				if(distances[j]<min_distance){
					min_distance_element=j;
					min_distance=distances[j];
}	
}	

		// observations[i].id=predicted[min_distance_element].id;
		// Make the code a bit faster by saving the index of closest match instead of the actual id
		observations[i].id=min_distance_element;

}}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	//  Update the weights of each particle using a mult-variate Gaussian distribution. 
	//  The observations are given in the VEHICLE'S coordinate system. Particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	
	std::vector<LandmarkObs> observations_map_coord( observations.size());
	double normalization_constant=1/(2*M_PI*std_landmark[0]*std_landmark[1]);
	
	for(int i=0;i<particles.size();i++){

	
	// Transform observation coordinates to map coordinates
		for(int j=0;j<observations.size();j++){
		observations_map_coord[j].x=particles[i].x+cos(particles[i].theta)*observations[j].x-sin(particles[i].theta)*observations[j].y;
		observations_map_coord[j].y=particles[i].y+sin(particles[i].theta)*observations[j].x+cos(particles[i].theta)*observations[j].y;
}	
		std::vector<LandmarkObs> withinrange_landmarks;
	
	for(int j=0;j<map_landmarks.landmark_list.size();j++){
	if(pow(map_landmarks.landmark_list[j].x_f-particles[i].x,2)+pow(map_landmarks.landmark_list[j].y_f-particles[i].y,2)<pow(sensor_range,2)){
//if((map_landmarks.landmark_list[j].x_f-particles[i].x)<sensor_range && (map_landmarks.landmark_list[j].y_f-particles[i].y)<sensor_range){
	
	withinrange_landmarks.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i,map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f});


}
}
	// Find closest map landmark to the sensor observations	
	dataAssociation(withinrange_landmarks, observations_map_coord);
	
	// Calculate probabilities of having observed the given landmarks, that is, the weights	
	double probs=1;
	for(int j=0;j<observations_map_coord.size();j++){

	int obs_id=observations_map_coord[j].id;
	int mu_x,mu_y;	
	
		
//	for(int k=0;k<withinrange_landmarks.size();k++){
//	if(obs_id==withinrange_landmarks[k].id){

	// With nearest neighbors, instead of actual landmark ids, use their index in the within_landmarks list
	int k=obs_id;	
	mu_x=withinrange_landmarks[k].x;
	mu_y=withinrange_landmarks[k].y;
//}
//}
 
    
	// Add likelihood of each single observation	
	probs = probs*exp(-(pow(observations_map_coord[j].x-mu_x,2)/2/std_landmark[0]+pow(observations_map_coord[j].y-mu_y,2)/2/std_landmark[1]))*normalization_constant;

}
	particles[i].weight=probs;
	}

}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight. 
	

	vector<double> weights;
	vector<Particle> particles2;
	double beta=0;

	//unload weight vector from particle list
	for(int i=0;i<particles.size();i++){
		weights.push_back(particles[i].weight);
	}

	// maximum weight
	double mw = *max_element(weights.begin(), weights.end());

	// Random start index
	int index=rand() % num_particles;
	default_random_engine gen;
	uniform_real_distribution<double> my_uniform_real_distribution(0.0, mw);

	// Sample with replacement by using the sampling wheel num_particles particles

	// Run the sampling wheel
	for(int i=0;i<particles.size();i++){

		beta+=my_uniform_real_distribution(gen)*2*mw;
		while(beta>weights[index]){
			beta-=weights[index];
			index=(index+1) % num_particles;
		}
		
		particles2.push_back(particles[index]);
	}
	particles=particles2;
	}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
