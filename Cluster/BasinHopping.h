#pragma once
#include "Common.h"
#include"tool.h"
#include<ctime>




class BH_Individual
{
public:
	int N;
	double *cood;
	double energy;

	BH_Individual();
	BH_Individual(int N);
	void set4N(int N);
	~BH_Individual();
	static BH_Individual* Individual(int N);
	void Copy(BH_Individual *bh,int N);
	BH_Individual& operator=(BH_Individual&);
};

class BasinHopping
{
public:
	BasinHopping();
	BasinHopping(int FeN,int N,int popsized,double E);
	~BasinHopping();
	//BH_Individual *X;
	
	
	//单个种群初始扰动比较选择；
	void start();
    void initDate(int N);
	void disturbance();
	double Greedy();

	//多种群初始扰动，比较选择；
	void Initialization();
	void compare(int independent_times,int fate,double durtime);
	void remake();
	void SphereCut();
	void SphereCutCompare();
	double compare_end();
	void adjustment();
	void exchange();
	void Genetic(double *cood,int N);

	//void rechange();
	//bool EndingCondition(BasinHopping &);
	void check();
	double changeRate();

private:
	int popsize;
	int N;
	int FeN;
	double Ebest;
	int Localtimes;
	int Localsuccesstimes;
	int S_times;
	BH_Individual *X;
	BH_Individual *Y;
	BH_Individual *Z;
	BH_Individual current;
	BH_Individual change;
	BH_Individual best;
	BH_Individual prebest;
	BH_Individual readbest;
	BH_Individual GAbest;

	
	PEnergy P_E;
	PEnergy P_Enew;
	PForce P_F;
	PForce P_Fnew;
	float successRate;
	double k[10];
	int kIndex;
	

};	

