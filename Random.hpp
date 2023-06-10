#ifndef RANDOMCLASSHEADERDEF
#define RANDOMCLASSHEADERDEF

#include <random>

class Random{
private:
	int RandomSeed;
	std::mt19937 mt;

public:
	Random();
	Random(int seed);

	int GetRandomSeed() const;

	//一様乱数
	double RAND_P(); // [0.0, 1.0]の範囲の浮動小数点数の一様乱数
	double RAND_PM(); // [-1.0, 1.0]の範囲の浮動小数点数の一様乱数
	double RAND_R(double a, double b); // [a, b]の浮動小数点数の範囲の一様乱数
	int RAND_N(int N); // [0, N-1]の範囲の整数の一様乱数
	int RAND_Z(int a, int b); // [a, b]の範囲の整数の一様乱数

	//正規分布
	double RAND_Normal(double m,double s);//平均m、標準偏差nで分布させる
};

#endif
