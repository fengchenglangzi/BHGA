#include "Tool.h"
#include "PE_Tool.h"
#include "geneticTool.h"
#include <iostream>
#include <stdlib.h>
extern string EnergyName;
using namespace std;

Tool::Tool(void)
{
	this->P_E = PE_Tool::FS_E;
	this->P_F = PE_Tool::FS_F;
}

Tool::~Tool(void)
{
}

void Tool:: printCood(double *cood, int N, char *root)
{
	FILE *fp = fopen(root,"w");
	for(int i=0;i<N;i++)
	{
		fprintf(fp,"%lf\n",cood[i]);
	}
	fclose(fp);
}

void Tool:: readCood(double *cood, int N, char *root)
{
	FILE *fp = fopen(root,"r");
	for(int i=0;i<N;i++)
	{
		fscanf(fp,"%lf",cood+i);
	}
	fclose(fp);
}

bool Tool::readDiamond(double *cood, int N,int differ /* = 0 */)
{
	char path[100];
	char line[100];
	int atomN;
	double E;
	int afterDiffer = N - differ;
	sprintf(path,"%s\\%s_%d.txt",EnergyName.c_str(),EnergyName.c_str(),afterDiffer);

	FILE *fp = fopen(path,"r");
	if (fp == NULL)
	{
		return false;
	}
	fscanf(fp,"%d",&atomN);
	fscanf(fp,"%s  %lf",line,&E);
	for (int i=0;i<afterDiffer;i++)
	{
		fscanf(fp,"%lf %lf %lf %lf",cood+i,cood+N+i,cood+2*N+i,cood+3*N+i);
	}
	fclose(fp);
	return true;
}
bool Tool::readDiamond1(double *cood, int N,int differ /* = 0 */)//读入Ir的最好数据
{
	char path[100];
	char line[100];
	int atomN;
	double E;
	int afterDiffer = N - differ;
	sprintf(path,"%s\\%s_%d\\1.txt",EnergyName.c_str(),EnergyName.c_str(),afterDiffer);

	FILE *fp = fopen(path,"r");
	if (fp == NULL)
	{
		return false;
	}
	fscanf(fp,"%d",&atomN);
	fscanf(fp,"%s %lf",line,&E);
	for (int i=0;i<afterDiffer;i++)
	{
		fscanf(fp,"%lf %lf %lf %lf",cood+i,cood+N+i,cood+2*N+i,cood+3*N+i);
	}
	fclose(fp);
	return true;
}

bool Tool::readDiamond2(double *cood, int N,int differ)//读入Pd最好数据；
{
char path[100];
	char line[100];
	int atomN;
	double E;
	int afterDiffer = N - differ;
	sprintf(path,"%s\\%s_%d\\2.txt",EnergyName.c_str(),EnergyName.c_str(),afterDiffer);

	FILE *fp = fopen(path,"r");
	if (fp == NULL)
	{
		return false;
	}
	fscanf(fp,"%d",&atomN);
	fscanf(fp,"%s %lf",line,&E);
	for (int i=0;i<afterDiffer;i++)
	{
		fscanf(fp,"%lf %lf %lf %lf",cood+i,cood+N+i,cood+2*N+i,cood+3*N+i);
	}
	fclose(fp);
	return true;
}

bool Tool::readDiamond_best(double *cood, int N,int FeN,int differ)//读入Pt最好数据；
{
	char path[100];
	char line[100];
	int atomN;
	double E;
	int afterDiffer = N - differ;
	sprintf(path,"%s\\%s\\%d-%d.txt",EnergyName.c_str(),EnergyName.c_str(),N,FeN);

	FILE *fp = fopen(path,"r");
	if (fp == NULL)
	{
		return false;
	}
	fscanf(fp,"%d",&atomN);
	fscanf(fp,"%s %lf",line,&E);
	for (int i=0;i<afterDiffer;i++)
	{
		fscanf(fp,"%lf %lf %lf %lf",cood+i,cood+N+i,cood+2*N+i,cood+3*N+i);
	}
	fclose(fp);
	return true;
}

bool Tool::readDiamond_atom_best(double *cood, int N,int FeN,int differ)//读到单金元素的最好结果用与做相似函数的比较
{
	char path[100];
	char line[100];
	int atomN;
	double E;
	int afterDiffer = N - differ;
	sprintf(path,"%s\\%s\\%d-%d.txt",EnergyName.c_str(),EnergyName.c_str(),N,FeN);

	FILE *fp = fopen(path,"r");
	if (fp == NULL)
	{
		return false;
	}
	fscanf(fp,"%d",&atomN);
	fscanf(fp,"%s %lf",line,&E);
	for (int i=0;i<afterDiffer;i++)
	{
		fscanf(fp,"%lf %lf %lf %lf",cood+i,cood+N+i,cood+2*N+i,cood+3*N+i);
	}
	fclose(fp);
	return true;
}

void Tool::printfDiamond(double *cood, int N, double E,char *path,int S,int H,int F,int FeN)
{
	FILE *fp;
	int i;
	double *coor,*x,*y,*z;
	coor = cood;x = cood + N; y = cood + 2*N; z = cood + 3*N;
	sprintf(path,"%s_%d_%d\\Pd %d.txt",EnergyName.c_str(),N,FeN,N);
	fp = fopen(path,"w+");
	//fprintf(fp,"第%d次独立重复试验结果：\n\n",F);
	////fprintf(fp,"1代表Fe，0代表Mo\n");
	//fprintf(fp,"%d\n",N);
	//fprintf(fp,"另一种元素占个数%d\n",FeN);
	//fprintf(fp,"run the best:  %lf\n",E);
	//fprintf(fp,"局部优化次数：%d.\t局部优化成功次数：%d.\t\n\n\n\n",S,H);
	for(i=0;i<3*N;i++)
	{
		fprintf(fp,"%lf\n",cood[i+N]);
	}

	fclose(fp);
}

void Tool::printfDiamond1(double *cood, int N, double E,char *path)
{
	FILE *fp0;
	int i;
	double *coor,*x,*y,*z;
	coor = cood;x = cood+N; y = cood + 2*N; z = cood + 3*N;
	fp0 = fopen("PdIr38AtomBestStructure.txt","w");
	fprintf(fp0,"%d\n",N);
	fprintf(fp0,"Pd\trun the best:  %lf\n",E);
	for(i=0;i<N;i++)
	{
		if(coor[i]<0.5){
			fprintf(fp0,"Pd\t%lf\t%lf\t%lf\n",x[i],y[i],z[i]);
		}
		else{
			fprintf(fp0,"Ir\t%lf\t%lf\t%lf\n",x[i],y[i],z[i]);
		}
	}
	fclose(fp0);
}

void Tool::printfDistance(double *cood, double *R,int N,int FeN,int F)
{
	FILE *fp0;
	int i;
	double *x;
	char path[200];
	double E0=0;
	double E1=0;
	x=cood;
	sprintf(path,"%s_%d_%d\\%d_%d_distance.txt",EnergyName.c_str(),N,FeN,N,FeN);
	fp0 = fopen(path,"w");
	fprintf(fp0,"%d\n",N);
//	fprintf(fp0,"Fe,run the best:  %lf\n",E);
	for(i=0;i<N;i++)
	{
		if(x[i]<0.1)
		{
			fprintf(fp0,"%1.0lf\t%lf\n",x[i],R[i]);
			E0 += R[i];
		}
	}

	for(i=0;i<N;i++)
	{
		if(x[i]>0.1)
		{
			fprintf(fp0,"%1.0lf\t%lf\n",x[i],R[i]);
			E1 += R[i];
		}
	}
	E0=E0/(N-FeN);
	E1=E1/(FeN);
	fprintf(fp0,"0的质心距离平均数是%lf\t1的质心距离平均数是%lf\n",E0,E1);
	fclose(fp0);
}

void Tool::printfEnergy(int keep,double E,int independent,char *path,int N,int FeN,int Localtimes,int Localsuccesstimes,double durtime)
{
	FILE *fp1;
	
	sprintf(path,"%s_%d_%d\\%s_%d_%dEnergy.txt",EnergyName.c_str(),N,FeN,EnergyName.c_str(),N,FeN);
	fp1 = fopen(path,"at+");

	//fprintf(fp1,"%d\t",S);
	fprintf(fp1,"%lf\t%d\t%d\t%d\t%lf\n",E,keep,Localtimes,Localsuccesstimes,durtime);
	fclose(fp1);
}


//void Tool::coodToDiamond(char *root, int N, char *path, PEnergy energy)
//{
//	double E;
//	int S=0,H=0,F=0;
//	double *cood,*R;
//	cood = (double *)malloc(3*N*sizeof(double));
//	R = (double *)malloc(N*N*sizeof(double));
//	readCood(cood,3*N,root);
//	Distance(cood,R,N);
//	E = energy(R,N);
//	printfDiamond(cood,N,E,path,S,H,F);
//	free(cood);
//	free(R);
//}

void Tool::Distance(double *cood, double *R, int N)
{
	int i,j;
	double *x,*y,*z;
    double Rmax=0;
	double Rmin = 10000;
	x = cood + N; y = cood + 2*N; z = cood + 3 * N; 
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			*(R+i*N+j)=(x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);
			*(R+i*N+j)= sqrt(*(R+i*N+j));
			*(R+j*N+i) = *(R+i*N+j);
			if(*(R+i*N+j)<Rmin)
				Rmin = *(R+i*N+j);
			if(*(R+i*N+j)>Rmax)
				Rmax = *(R+i*N+j);
		}
	}
}

//局部优化定义-最陡梯度下降法；
double Tool::Local(double *cood, int N,PEnergy PE_E, PForce PE_F)
{
	int i;
	int s = 1;
	double alpha = 0.01;
	double *R,*tempR,*tempCood;
	double *tempX,*tempY,*tempZ,*tempCoor;
	double E0,E1;
	double *FF;
	double *FX,*FY,*FZ,*coor,*x,*y,*z;
	double Fmax;

	/*
		
		最陡下降法迭代终止条件： 梯度是否接近于0，迭代次数是否达到预设值，
		
	*/

	double e = 1E-6;
	int smax = 3000;
	double alphaE = 1E-10;

	R = (double *)calloc(N*N,sizeof(double));
	tempR = (double *)calloc(N * N,sizeof(double));
	tempCood = (double *)calloc(4 * N,sizeof(double));
	
	FF = (double *)calloc(3 * N,sizeof(double));

	coor = cood;x = cood + N; y = cood + 2 * N; z = cood + 3 * N;
	FX = FF; FY = FF + N; FZ = FF + 2 * N; 
	tempCoor = tempCood;tempX = tempCood + N; tempY = tempCood + 2*N; tempZ = tempCood + 3*N;
//	Tool::printCood(cood,N);
//	Tool::readCood(cood,3* N,"s.txt");
	Distance(cood,R,N);
	///E0 = PE_Tool::FS_E(cood,R,N);
//	Fmax = PE_Tool::FS_F(cood,R,FF,N);
//	E0 = PE_Tool::Johnson_E(R,N);
//	Fmax = PE_Tool::Johnson_F(cood,R,FF,N);
//	Tool::printCood(FF,N,"FF.txt");
	E0 = PE_E(cood,R,N);

	// cout << "未进行梯度优化前的能量值 E0：" << E0 << endl;

//	cout<<"Without Localresearch energy is :"<<E0;
	Fmax = PE_F(cood,R,FF,N);	
	while(1)
	/*
		
		------------------------最陡下降法执行过程---------------------------------
		------------------------------------------------------------------------
		-------------迭代终止条件： 梯度是否接近于0，迭代次数是否达到预设值--------------
		------------------------------------------------------------------------
		------------------------------------------------------------------------
		
	*/
	{
		//cout<<"Fmax"<<Fmax<<endl;
		if(Fmax <= e || s > smax || alpha < e)
		{

			if(s>smax)
				//cout<<"Fmax"<<Fmax<<"  "<<s<<"  "<<alpha<<"  "<<E0<<endl;
				
				//getchar();
			//	if(E0>-55)
			//	{
				//	Tool::printfDiamond(cood,N,E0,"badCluster.txt");
				//	printf("badCluster");
				//	getchar();
			//	}
			break;
		}
	/*	for(i = 0; i < 3 * N; i++){
			FF[i] = (FF[i]>0.5)?0.5:((FF[i]<-0.5)?-0.5:0);
		}*/
		for(i=0; i<N ;i++)
		{
			tempX[i] = x[i] + alpha * FX[i]/Fmax;
			tempY[i] = y[i] + alpha * FY[i]/Fmax;
			tempZ[i] = z[i] + alpha * FZ[i]/Fmax;
		//	tempX[i] = x[i] + alpha * FX[i];
		//	tempY[i] = y[i] + alpha * FY[i];
		//	tempZ[i] = z[i] + alpha * FZ[i];
		}

		Distance(tempCood,tempR,N);
	//	E1 = PE_Tool::FS_E(tempR,N);
	//	E1 = PE_Tool::Johnson_E(tempR,N);
		E1 = PE_E(cood,tempR,N);
		if(E1 < E0)
		{
			alpha = alpha * 1.1;
			E0 = E1;
			for(i = 0 ;i < N;i ++)
			{
				x[i] = tempX[i];
				y[i] = tempY[i];
				z[i] = tempZ[i];
				cood[i + N] = x[i];
				cood[i + 2*N] = y[i];
				cood[i + 3*N] = z[i];
			}
			for(i=0;i<N*N;i++)
				R[i] = tempR[i];
		//	Fmax = PE_Tool::FS_F(cood,R,FF,N);
		//	Fmax = PE_Tool::Johnson_F(cood,R,FF,N);
			Fmax = PE_F(cood,R,FF,N);
		}
		else
			alpha = alpha * 0.6;
		s ++;
	}
	free(R);
	free(tempR);
	free(tempCood);
	free(FF);
	return E0;
}

//double Tool::Localcoor(double *cood, int N,int FeN,PEnergy PE_E, PForce PE_F)
//{
//	int i;
//	double *R;
//	double *coor,*x,*y,*z;
//
//	R = (double *)calloc(N*N,sizeof(double));
//	coor = cood;x = cood + N; y = cood + 2 * N; z = cood + 3 * N;
//	Distance(cood,R,N);
//	E0 = PE_E(cood,R,N);
//	geneticTool::Genetic(cood,N,FeN);
//	free(R);
//	return E0;
//}

double Tool::LocalWithR(double *cood,double *R,int N,PEnergy PE_E, PForce PE_F)
{
	int i,j;
	int s = 1;
	double alpha = 0.01;
	double *tempR;
	double E0,E1;
	double *FF;
	double Fmax;
	double e = 1E-6;
	int smax = 3000;
	double alphaE = 1E-10;

	tempR = (double *)calloc(N * N,sizeof(double));
	FF = (double *)calloc(N * N,sizeof(double));

	E0 = PE_E(cood,R,N);
	Fmax = PE_F(NULL,R,FF,N);
	while(1)
	{
		if(Fmax <= e || s > smax || alpha < e)
		{
			if(s>smax)
				cout<<"Fmax"<<Fmax<<"  "<<s<<"  "<<alpha<<"  "<<E0<<endl;
			break;
		}

		for(i=0; i<N-1 ;i++)
		{
			for(j=i+1;j<N;j++)
				*(tempR+i*N+j) = *(R+i*N+j) + alpha * (*(FF+i*N+j))/Fmax;
		}

		E1 = PE_E(cood,tempR,N);
		if(E1 < E0)
		{
			alpha = alpha * 1.1;
			E0 = E1;

			for(i=0; i<N-1 ;i++)
			{
				for(j=i+1;j<N;j++)
					*(R+i*N+j) = *(tempR+i*N+j);
			}

			Fmax = PE_F(NULL,R,FF,N);
		}
		else
			alpha = alpha * 0.6;
		s ++;
	}

	free(tempR);
	free(FF);
	return E0;
}

int *Tool::RandPerm(int N, int K)
{
	int *allP;
	int *p;
	allP = (int *)malloc(N * sizeof(int));
	p = (int *)malloc(K * sizeof(int));
	for(int i = 0; i < N; i++)
		allP[i] = i;
	for(int i=0;i < K;i++)
	{
		int point = rand()%(N-i);
		int temp = allP[i];
		allP[i] = allP[point+i];
		allP[point+i] = allP[i];
		p[i] = allP[i];
	}
	free(allP);
	return p;
}

void Tool::SphereCutSplice(double *fr, double *mr, double *child1, double *child2, int N)
{
	int j,k,kexist,kmid,ncross;
	double x,y,z,mid,frcentre[3],mrcentre[3];
	double *fr3_coor,*fr3_x,*fr3_y,*fr3_z,*mr3_coor,*mr3_x,*mr3_y,*mr3_z;
	int *existId;
	double *frd,*mrd,*frdrank,*mrdrank;

//	cout<<1;

	existId = (int *)malloc(N * sizeof(int));
	frd = (double *)malloc(N * sizeof(double));
	mrd = (double *)malloc(N * sizeof(double));
	frdrank = (double *)malloc(N * sizeof(double));
	mrdrank = (double *)malloc(N * sizeof(double));

	for(j=0;j < N;j++){//father
		fr3_coor = fr;
		fr3_x = fr + N;
		fr3_y = fr + 2*N;
		fr3_z = fr + 3 * N;
		mr3_coor = mr;
		mr3_x = mr + N;//mother
		mr3_y = mr + 2*N;
		mr3_z = mr + 3 * N;
	}
	x=0;y=0;z=0;
	for(j=0;j < N;j++){
		x = x + fr3_x[j];
		y = y + fr3_y[j];
		z = z + fr3_z[j];	
	}
	frcentre[0] = x/N;
	frcentre[1] = y/N;
	frcentre[2] = z/N;
	x=0;y=0;z=0;
	for(j = 0;j < N;j++){
		x = x + mr3_x[j];
		y = y + mr3_y[j];
		z = z + mr3_z[j];	
	}
	mrcentre[0] = x / N;
	mrcentre[1] = y / N;
	mrcentre[2] = z / N;
	//centre
	for(j = 0; j < N;j++){
		fr3_x[j] = fr3_x[j] - frcentre[0];
		fr3_y[j] = fr3_y[j] - frcentre[1];
		fr3_z[j] = fr3_z[j] - frcentre[2];
		mr3_x[j] = mr3_x[j] - mrcentre[0];
		mr3_y[j] = mr3_y[j] - mrcentre[1];
		mr3_z[j] = mr3_z[j] - mrcentre[2];
	}
	//distance
	for(j = 0;j < N;j++){
		frd[j] = sqrt(fr3_x[j] * fr3_x[j] + fr3_y[j] * fr3_y[j] + fr3_z[j] * fr3_z[j]);
		mrd[j] = sqrt(mr3_x[j] * mr3_x[j] + mr3_y[j] * mr3_y[j] + mr3_z[j] * mr3_z[j]);
		frdrank[j]=frd[j];
		mrdrank[j]=mrd[j];
	}
	//let not two ones are exactly the same
	for( k = 0;k < N-1;k++){
		for(j=k+1;j<N;j++){
			if(frd[k]==frd[j]){
				frd[k]=frd[k]+(2 + RAND)*(1e-16);
				frdrank[k]=frd[k];
			}
			if(mrd[k]==mrd[j]){
				mrd[k]=mrd[k]+(2 + RAND)*(1e-16);
				mrdrank[k]=mrd[k];
			}
		}
	}
	//rank
	for(k=N;k>1;k--){
		for(j=1;j<k;j++){
			if(frdrank[j-1]>frdrank[j]){
				mid=frdrank[j];
				frdrank[j]=frdrank[j-1];
				frdrank[j-1]=mid;
			}
		}
		for(j=1;j<k;j++){
			if(mrdrank[j-1]>mrdrank[j]){
				mid=mrdrank[j];
				mrdrank[j]=mrdrank[j-1];
				mrdrank[j-1]=mid;
			}
		}
	}
	//crossover by sphere
	//check if there exists such a sphere
	kexist=0;
	for(k=1;k<N-1;k++){
		if(!((frdrank[k-1]>mrdrank[k]) || (frdrank[k]<mrdrank[k-1]))){
			existId[kexist]=k;
			kexist=kexist+1;
		}
	}

	//cout<<"k="<<kexist<<";";
	if(kexist!=0){
		//random number of crossover atoms
		
		kmid=(int)ceil(kexist* RAND);
		if(kmid==0)
			kmid = 1;
	//	cout<<"kmid="<<kmid<<";";
		ncross=existId[kmid-1];
		k=0;
		for(j=0;j<N;j++){
			if(frd[j]<frdrank[ncross]-1e-16){
				//*(child1 + k) = fr3_coor[j];
				*(child1 + N + k) = fr3_x[j];
				*(child1 + 2*N + k) = fr3_y[j];
				*(child1 + 3*N + k) = fr3_z[j];
				k=k+1;
			}
			if(mrd[j]>mrdrank[ncross-1]+1e-16){
				//*(child1 + k) = mr3_coor[j];
				*(child1 +  N + k) = mr3_x[j];
				*(child1 + 2*N + k) = mr3_y[j];
				*(child1 + 3*N + k) = mr3_z[j];
				k=k+1;
			}
		}
		//cout<<"child1="<<k<<";";
		k=0;
		for(j=0;j<N;j++){
			if(frd[j]>frdrank[ncross-1]+1e-16){
				//*(child2 + k) = fr3_coor[j];
				*(child2 + N + k) = fr3_x[j];
				*(child2 + 2*N + k) = fr3_y[j];
				*(child2 + 3*N + k) = fr3_z[j];
				k=k+1;
			}
			if(mrd[j]<mrdrank[ncross]-1e-16){
				//*(child2 + k) = mr3_coor[j];
				*(child2 + N + k) = mr3_x[j];
				*(child2 + 2*N + k) = mr3_y[j];
				*(child2 + 3*N + k) = mr3_z[j];
				k=k+1;
			}
		}
		//cout<<"child2="<<k<<";";
	}
	else{
		for(j=0;j<4 * N;j++){
			*(child1 + j) = fr[j];
			*(child2 + j) = mr[j];
		}
	}

	for(j=0;j<N;j++)
	{
		*(child1 + j) = fr3_coor[j];//遗传交叉的函数调用
		*(child2 + j) = mr3_coor[j];//用来得到比较好的结果
	}
	//Genetic(child1.cood,child2.cood,N);
	//cout<<2;
	free(existId);
	free(frd);
	free(mrd);
	free(frdrank);
	free(mrdrank);
	//cout<<endl;
}

void Tool::testR(int N){
	double *cood = new double[3*N];
	double *R = new double[N*N];
	double E0,E1;

	//随机产生一个结构
	for (int i=0; i < 3*N; i++)
	{
		cood[i] = (RAND-0.5) * 2.75 * pow((double)N,1.0/3);
	}
	Distance(cood,R,N);
	E0 = PE_Tool::FS_E(cood,R,N);
	//对R进行最速下降法
	E1 = LocalWithR(cood,R,N,PE_Tool::FS_E,PE_Tool::FS_F);
	printf("原子总数为%d;\n",N);
	printf("初始能量:%lf\n",E0);
	printf("优化r后的能量:%lf\n",E1);

	//得到R下，三角形个数和无法围成合法三角形的个数
	testThree(R,N);
}

bool Tool::testThree(double *R,int N){
	int numOfAll = 0;
	int numOfBad = 0;
	bool rational = true;
	printf("无法构成三角形的节点和三条边的长度如下：\n");
	for (int i = 0; i < N-2; i++)
	{
		for(int j=i+1; j < N-1; j++){
			for (int m = j+1; m < N; m++)
			{
				numOfAll ++;
				double a = *(R + i * N + j);
				double b = *(R + i* N + m);
				double c = *(R + j * N + m);
				if (a + b<c || a + c < b || b + c <a)
				{
					rational = false;
					numOfBad ++;
					printf("%d,%d,%d\t%lf %lf %lf\n",i,j,m,a,b,c);
				} 
			}
		}
	}
	printf("一共%d三角形,其中坏的有%d个\n",numOfAll,numOfBad);
	return rational;
}


void Tool::adjustment(double *cood,int N,int FeN)
{
	int oneNum=0;
	int ren,g,o;
	int a[100];
	int temp;
	int zeroNum=0;
	cout<<"开始调整"<<endl;
	for (int i=0;i<50;i++)
	{
		ren = N;
		//产生随机的数组
		for (g=0;g<N;g++)
		{
			a[g] = g+1;
			//cout<<a[g]<<endl;
		}

		for (g=0;g<N-1;g++)
		{
			o = rand()%ren;
			temp=a[o];a[o]=a[ren-1];a[ren-1]=temp;
			ren--;		
		}

		for (g=0;g<N;g++)
		{
			if (a[g]<=FeN)
			{
				cood[g] = 1;
			}
			else 
			{
				cood[g] = 0;
			}
		}

	}

}

void Tool::distool(double *cood,double *R,int N,int FeN,double temp,int F){
	
	int i=0,j=0;
	double dcentre[3] = {0};
	double *cx,*cy,*cz;
	double r0;
	temp = 0;r0 = 2.75; //此时随机生成的范围应该比之前初始化的要小才合适;

			//if ((Ebest - Z[i].energy)/Ebest >= 0.07)
		
		//cout<<"第"<<i<<"个种群需要处理"<<"其能量是："<<Z[i].energy<<endl;
		cx = cood + N;
		cy = cood + 2*N;
		cz = cood + 3*N;
		for (j=0;j<N;j++)
			{
				dcentre[0] += cx[j];
				dcentre[1] += cy[j];
				dcentre[2] += cz[j];

				dcentre[0] /= N;
				dcentre[1] /= N;
				dcentre[2] /= N;

			}

		for (j=0;j<N;j++)
			{
				R[j] = sqrt((cx[j]-dcentre[0])*(cx[j]-dcentre[0]) + (cy[j]-dcentre[1])*(cy[j]-dcentre[1]) + (cz[j]-dcentre[2])*(cz[j]-dcentre[2]));
				temp = R[0];
				//cout<<"距离质心最远的距离是："<<distancecentre[j]<<endl;
			}
		//printfDistance(cood,R,N,FeN,F);
		
		for (j=1;j<N;j++)
		{
			if (temp<R[j])
			{
				temp = R[j];

			}


		}
		
	
}


void Tool::distoolprintf(double *cood,double *R,int N,int FeN,double temp,int F){
	
	int i=0,j=0;
	double dcentre[3] = {0};
	double *cx,*cy,*cz;
	double r0;
	temp = 0;r0 = 2.75; //此时随机生成的范围应该比之前初始化的要小才合适;

			//if ((Ebest - Z[i].energy)/Ebest >= 0.07)
		
		//cout<<"第"<<i<<"个种群需要处理"<<"其能量是："<<Z[i].energy<<endl;
		cx = cood + N;
		cy = cood + 2*N;
		cz = cood + 3*N;
		for (j=0;j<N;j++)
			{
				dcentre[0] += cx[j];
				dcentre[1] += cy[j];
				dcentre[2] += cz[j];

				dcentre[0] /= N;
				dcentre[1] /= N;
				dcentre[2] /= N;

			}

		for (j=0;j<N;j++)
			{
				R[j] = sqrt((cx[j]-dcentre[0])*(cx[j]-dcentre[0]) + (cy[j]-dcentre[1])*(cy[j]-dcentre[1]) + (cz[j]-dcentre[2])*(cz[j]-dcentre[2]));
				temp = R[0];
				//cout<<"距离质心最远的距离是："<<distancecentre[j]<<endl;
			}
		printfDistance(cood,R,N,FeN,F);
		
		for (j=1;j<N;j++)
			{
				if (temp<R[j])
				{
					temp = R[j];

				}


			}
		
	
}

void Tool::CoordinationNumber(double *cood,double *R,int N,int FeN,int F)
{
	int i,j;
	int *neighbour;
	//int mina0;
	neighbour = (int *)calloc(N,sizeof(int));
	//
	
	double *x,*y,*z,*co;
    double Rmax=0;
	double Rmin = 10000;
	co = cood;x = cood + N; y = cood + 2*N; z = cood + 3 * N; 
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			*(R+i*N+j)=(x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);
			*(R+i*N+j)= sqrt(*(R+i*N+j));
			*(R+j*N+i) = *(R+i*N+j);
			if(*(R+i*N+j)<Rmin)
				Rmin = *(R+i*N+j);
			if(*(R+i*N+j)>Rmax)
				Rmax = *(R+i*N+j);
		}
	}
	//
	for(i=0;i<N;i++){
		neighbour[i]=0;
	}
	for(i=0;i<N-1;i++){	
		for(j=i+1;j<N;j++){	
			if(*(R+i*N+j)< Rmin+0.4){
				neighbour[i]++;
				neighbour[j]++;
			}
			
		}//cout<<"第"<<i<<'\t'<<co[i]<<"个原子的配位是："<<neighbour[i]<<endl;
	}
	//cout<<"最短的原子距离是："<<Rmin<<endl;
	printfCoordinationNumber(cood,neighbour,N,FeN,F);
}

void Tool::printfCoordinationNumber(double *cood,int *R,int N,int FeN,int F)
{
	FILE *fp0;
	int i;
	double *x;
	char path[200];
	double E0=0;
	double E1=0;
	x=cood;
	sprintf(path,"%s_%d_%d\\%d_%d_CoordinationNum.txt",EnergyName.c_str(),N,FeN,N,FeN);
	fp0 = fopen(path,"w");
	fprintf(fp0,"%d\n",N);
//	fprintf(fp0,"Fe,run the best:  %lf\n",E);
	for(i=0;i<N;i++)
	{
		if(x[i]<0.1)
		{
			fprintf(fp0,"%1.0lf\t%d\n",x[i],R[i]);
			E0 += R[i];
		}
	}

	for(i=0;i<N;i++)
	{
		if(x[i]>0.1)
		{
			fprintf(fp0,"%1.0lf\t%d\n",x[i],R[i]);
			E1 += R[i];
		}
	}
	E0=E0/(N-FeN);
	E1=E1/(FeN);
	cout<<"0的配位数平均数是"<<E0<<endl;
	cout<<"1的配位数平均数是"<<E1<<endl;
	fprintf(fp0,"0的配位数平均数是%lf\t1的配位数平均数是%lf\n",E0,E1);
	fclose(fp0);
	
}


//描述相似函数部分

void Tool::similar_fun(double *cood1,double *cood2,int N,int FeN)
{
	int i,j;
	double s=0,q=0;
	double dcentre1[3] = {0};
	double dcentre2[3] = {0};
	double *cx1,*cy1,*cz1;
	double *cx2,*cy2,*cz2;
	double *R1,*R2;
	double r0;
	R1 = (double *)calloc(N*N,sizeof(double));
	R2 = (double *)calloc(N*N,sizeof(double));
	r0 = 2.75; //此时随机生成的范围应该比之前初始化的要小才合适;
		cx1 = cood1 + N;
		cy1 = cood1 + 2*N;
		cz1 = cood1 + 3*N;
		cx2 = cood2 + N;
		cy2 = cood2 + 2*N;
		cz2 = cood2 + 3*N;
		for (j=0;j<N;j++)
			{
				dcentre1[0] += cx1[j];
				dcentre1[1] += cy1[j];
				dcentre1[2] += cz1[j];

				dcentre1[0] /= N;
				dcentre1[1] /= N;
				dcentre1[2] /= N;

			}

		for (j=0;j<N;j++)
			{
				R1[j] = sqrt((cx1[j]-dcentre1[0])*(cx1[j]-dcentre1[0]) + (cy1[j]-dcentre1[1])*(cy1[j]-dcentre1[1]) + (cz1[j]-dcentre1[2])*(cz1[j]-dcentre1[2]));
			}

		for (j=0;j<N;j++)
			{
				dcentre2[0] += cx2[j];
				dcentre2[1] += cy2[j];
				dcentre2[2] += cz2[j];

				dcentre2[0] /= N;
				dcentre2[1] /= N;
				dcentre2[2] /= N;

			}

		for (j=0;j<N;j++)
			{
				R2[j] = sqrt((cx2[j]-dcentre2[0])*(cx2[j]-dcentre2[0]) + (cy2[j]-dcentre2[1])*(cy2[j]-dcentre2[1]) + (cz2[j]-dcentre2[2])*(cz2[j]-dcentre2[2]));
			}

		Quicksort(R1,N);
		Quicksort(R2,N);


		for(i=0;i<N;i++)
		{
			//cout<<"第11111得到的质心距离的排序"<<R1[i]<<endl;
			//cout<<"第22222得到的质心距离的排序"<<R2[i]<<endl;
			q += (R1[i]-R2[i])*(R1[i]-R2[i]);
		}
		q = sqrt(q/N);

		s = 1/(1+q);

		//cout<<"得到的相似函数的q是："<<q<<endl;
		cout<<"得到的相似函数的s是："<<s<<endl;
}


void Tool::Quicksort(double *R, int N)
{
	int i,j;
	double temp;
	for(i=0;i<N;i++)
	{
		for(j=i+1;j<N;j++)
		{
			if(R[i]<R[j])
			{
				temp = R[i];
				R[i] = R[j];
				R[j]  =temp;
			}
		}
	//cout<<"得到的质心距离的排序"<<R[i]<<endl;
	}
}