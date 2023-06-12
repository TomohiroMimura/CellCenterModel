#include <random>
#include "Random.hpp"

Random::Random(){
	std::random_device rnd;
	RandomSeed = rnd();
	mt.seed(RandomSeed);
}

Random::Random(int seed){
	RandomSeed = seed;
	mt.seed(RandomSeed);
}


int Random::GetRandomSeed() const{
	return RandomSeed;
}


double Random::RAND_P(){
	std::uniform_real_distribution<> rnd(0.0, 1.0);
	return rnd(mt);
}

double Random::RAND_PM(){
	std::uniform_real_distribution<> rnd(-1.0, 1.0);
	return rnd(mt);
}

double Random::RAND_R(double a, double b){
	std::uniform_real_distribution<> rnd(a, b);
	return rnd(mt);
}

int Random::RAND_N(int N){
	std::uniform_int_distribution<> rnd(0, N-1);
	return rnd(mt);
}

int Random::RAND_Z(int a, int b){
	std::uniform_int_distribution<> rnd(a, b);
	return rnd(mt);
}

double Random::RAND_Normal(double m,double n){
	std::normal_distribution<> rnd(m,n);
	return rnd(mt);
}