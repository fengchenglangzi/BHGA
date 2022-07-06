#pragma once
#include"Tool.h"
#include "PE_Tool.h"
#define RAND1 rand()/(RAND_MAX+1.0)


class GA_Individual
{
public:
	int N;
	//double *coor;
	double *cood;
	double energy;
	double value;

	GA_Individual();
	GA_Individual(int N);
	~GA_Individual();
	void setN(int N);
	static GA_Individual* Individual(int N);
	void Copy(GA_Individual *one,int N);
	void Copy1(GA_Individual *one,int N);
	GA_Individual& operator=(GA_Individual&);
	//double docking(double *cood,int N,int FeN);
	//void changeIndividual(GA_Individual one,GA_Individual two,int N);
};


class geneticTool
{
public:
	geneticTool(void);
	~geneticTool(void);
	geneticTool(int N,int FeN);

	//使用基因局部优化定义
	void Genetic(double *cood,int N,int FeN);
	//遗传算法的内容的定义
	void select_operator(double *coor_L,double BestEnergy,double WorstEnergy);
	void mutation_operator();
	void crossover_operator();
	void adjustment_operator(double *cood,double *R);
	double GAdocking(double *coorL,double *cood,double *coodL,int N,int FeN,double *R);
	void GetBestcood(double *coorL,double *cood,int N);

private:
	int N;
	int FeN;
	int GApopsize;
	GA_Individual *G;
	GA_Individual *A;
	GA_Individual *M;
	GA_Individual pbest;
	GA_Individual pop;
	GA_Individual best_pop;
	GA_Individual worst_pop;

	PEnergy P_E;
	PForce P_F;
};


