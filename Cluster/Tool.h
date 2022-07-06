#pragma once
#include "Common.h"

class Tool
{
public:
	Tool(void);
	~Tool(void);
	static void Distance(double *cood, double *R, int N);
	static int* RandPerm(int N, int K);
	static double Local(double *cood, int N, PEnergy PE_E, PForce PE_F);
	static double Localcoor(double *cood, int N, PEnergy PE_E, PForce PE_F);
	static double LocalWithR(double *cood,double *R, int N, PEnergy PE_E, PForce PE_F);
	static void SphereCutSplice(double *fr, double *mr, double *child1, double *child2, int natom);
	static double LJ_E(double *cood);
	static void printCood(double *cood, int N, char *root);
	static void readCood(double *cood, int N, char *root);
	static void printfDiamond(double *cood, int N, double E,char *path,int S,int H,int F,int FeN);
	static void printfDiamond1(double *cood, int N, double E,char *path);

	static void printfDistance(double *cood, double *R,int N,int FeN,int F);
	static void printfCoordinationNumber(double *cood,int *R,int N,int FeN,int F);
	static bool readDiamond(double *cood, int N,int differ = 0);
	static bool readDiamond1(double *cood, int N,int differ = 0);
	static bool readDiamond2(double *cood, int N,int differ);
	static bool readDiamond_best(double *cood, int N,int FeN,int differ);
	static bool readDiamond_atom_best(double *cood_atom,int N,int FeN,int differ);
	static void coodToDiamond(char *root, int N, char *path, PEnergy energy);
	static void printfEnergy(int keep,double E,int independent,char *path,int N,int FeN,int Localtimes,int Localsuccesstimes,double durtime);
	static void adjustment(double *cood,int N,int FeN);
	static void distool(double *cood,double *R,int N,int FeN,double temp,int F);
	static void distoolprintf(double *cood,double *R,int N,int FeN,double temp,int F);
	static void CoordinationNumber(double *cood,double *R,int N,int FeN,int F);
	//static void Genetic(double *cood1,double *cood2,int N);
	static void testR(int N);
	static bool testThree(double *R, int N);
	static void Quicksort(double *R, int N);
	static void similar_fun(double *cood1,double *cood2,int N,int FeN);

private:
	PEnergy P_E;
	PForce P_F;
};

