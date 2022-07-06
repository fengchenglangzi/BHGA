#include "BH_Manage.h"
#include "BasinHopping.h"
#include <direct.h>
#include "Common.h"
#include "Tool.h"
#include "time.h"
extern string EnergyName;

BH_Manage::BH_Manage(void)
{
	this->popNum = 0;
	//this->pop = NULL;
	
}

BH_Manage::~BH_Manage(void)
{
	cout<<"BH_Manage析构"<<endl;
}

void BH_Manage::process(){
	
	int m=0;
	clock_t start,finish;
	double duration;
	char energyFile[50];
	this->atomNum = 38;
	start = clock(); 
	this->popNum = 30;
	this->objbest = -1000;
	//for(this->atomNum = 24;this->atomNum <= 24;this->atomNum+=6){
	//for(this->FeNum = 22;this->FeNum <= 22;this->FeNum++){
	FeNum = 6;
	BasinHopping BH1(FeNum,atomNum,popNum,objbest);
	//BasinHopping BH1;
	//BH1.initDate(atomNum);

	/*
	
		数据持久化位置：Gupta_38_N

	*/

	sprintf(energyFile,"%s_%d_%d",EnergyName.c_str(),atomNum,FeNum);
	mkdir(energyFile);
	//BH1.check();
	//cout<<"第"<<Day<<"次独立实验"<<endl;
	
	
	// 表示独立20次实验 一次迭代200步
	for(int i=1;i< 20;i++)
	{  
		//disand(atomNum,FeNum,i);
		int k=200;
		BH1.Initialization();
		// BH1.initDate(atomNum);
		while(k--)

		{	
			finish = clock();
			duration = (double)(finish - start) / CLOCKS_PER_SEC;
			BH1.remake();//
			BH1.compare(i,k,duration);
			BH1.exchange();
			if (k%5 == 0)
			{
				BH1.SphereCut();

			}

			//if (!EndingCondition1(BH1,FeNum))//表示再某一次独立实验中找到了最优值，此时进入下一次循环
			//{
			//	//disand(atomNum,FeNum,i);
			//	cout<<"第"<<i<<"次独立实验结束"<<endl<<endl;
			//	break;
			//}
		}

	//	
	//}
	}

	/*while(1)
	{

	BH1.disturbance();
	BH1.Greedy();
	}

	while(EndingCondition(BH1))
	{
	BH1.disturbance();
	BH1.Greedy();
	}*/
}

bool BH_Manage::EndingCondition1(BasinHopping &BH,int FeNum){
	
	//double a[100]={-108.967,-109.767,-110.577,-111.390,-112.213,-113.031,-113.859,-114.691,-115.532,-116.213,-116.919,-117.603,-118.306,-118.963,-119};
	/*double *a;
	Tool::readDiamond2(a);
	int temp = 0;
	temp = this->atomNum - FeNum;*/
	//if(BH.compare_end() <= a[temp])
	/*if (Tool::readDiamond1(Z[0].cood,N,0))
		{
			cout<<"已经读到最优数据"<<endl;
			Z[0].energy = Tool::Local(Z[0].cood,N,P_E,P_F);
			cout<<Z[0].energy<<endl;
		}*/
	if(BH.compare_end() <= objbest)

		//cout<<"迭代到此结束"<<endl;
	//cout<<"最后的实验结果是："<<current.energy<<endl;
	return false;

}

void BH_Manage::checkin(){
	int Nin;
	//int success=0;
	//BasinHopping BH(10);
	double objEnergy;
	//double E,successRate = 0;
	cout<<"请输入原子数："<<endl;
	cin>>Nin;
	cout<<"请输入Fe原子数目："<<endl;
	cin>>FeNum;
	cout<<"请输入目标能量："<<endl;
	cin>>objEnergy;
	cout<<"请输入种群大小："<<endl;
	cin>>popNum;
	atomNum=Nin;
	this->FeNum=FeNum;
	objbest=objEnergy;
	this->popNum=popNum;
	/*current = 0;
	for (current = 0; current < repeat; current++ )
	{
		E = BH.Greedy();
		if (fabs(E - objEnergy) < 0.01)
		{
			success ++;
		}
	}
	successRate = success*100.0 / repeat;
	cout<<"一共进行"<<repeat<<"次，其中成功"<<success<<"次，成功率为"<<successRate<<"%"<<endl;*/
}

bool BH_Manage::EndingCondition(BasinHopping &BH){
	//BH_Individual current;
	//BH_Individual best;
	double m=0.01;
	//if(BH.changeRate()<m)//涉及到平均变化率、、有问题
	if(fabs((BH.Greedy()-objbest)/objbest)<=m)
			
		cout<<"迭代到此结束"<<endl;
		//cout<<"最后的实验结果是："<<current.energy<<endl;
		return false;
	
}


void BH_Manage::set_objbest(double setbest){

	objbest=setbest;
}

double BH_Manage::get_objbest(){
	return objbest;
}

//10.20两个问题：1.终止条件中没有输入两个能量值：2.扰动过程中生成的随机用时太短，不能随时间变化，
//第一次和第二次的时间间隔太短，没能随机到数组;今天解决



//原面原子偏距和配位数
void BH_Manage::disand(int N,int FeN,int F)
{
	int i;
	double temp=0;
	double *cood,*R,*Rcoor,*cood_atom;
	R = (double *)calloc(N*N,sizeof(double));
	Rcoor = (double *)calloc(N*N,sizeof(double));
	for( i = 0; i < N; i++)
	{	
		cood = (double *)malloc(4*N * sizeof(double));
		cood_atom = (double *)malloc(4*N * sizeof(double));
	}

	if(Tool::readDiamond_best(cood,N,FeN,0))
	{
		//cout<<"已经读到cood"<<endl;
		//Tool::CoordinationNumber(cood,Rcoor,N,FeN,F);
		//Tool::distoolprintf(cood,R,N,FeN,temp,F);
		//if(Tool::readDiamond_atom_best(cood_atom,N,0,0))//相似函数与Pt的比较
		if(Tool::readDiamond_atom_best(cood_atom,N,0,0))//相似函数与Fe的比较
		{
			//Tool::similar_fun(cood,cood_atom,N,FeN);
		}
	}
	
		
}