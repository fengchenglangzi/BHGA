#include "geneticTool.h"
#include"Tool.h"
#include "PE_Tool.h"




GA_Individual* GA_Individual::Individual(int N)
{
	GA_Individual *obj = (GA_Individual*)malloc(sizeof(GA_Individual));
	//obj->coor = (double *)malloc(N * sizeof(double));
	obj->cood = (double *)malloc(4*N * sizeof(double));
	return obj;
}

GA_Individual::GA_Individual(void){
	//this->coor = NULL;
	this->cood = NULL;
}

GA_Individual::GA_Individual(int N){
	this->setN(N);
}

GA_Individual::~GA_Individual(void){

	/*if(this->coor != NULL)
	{
		free(this->coor);
	}
*/
	if(this->cood != NULL)

	{free(this->cood);}
}

void GA_Individual::setN(int N)
{
	this->N = N;
	//this->coor = (double *)malloc(N * sizeof(double));
	this->cood = (double *)malloc(4*N * sizeof(double));
}
 
GA_Individual& GA_Individual::operator=(GA_Individual &ind)
{
	if (N != ind.N && N != 0)
	{
		throw "GA����ԭ�������Ե�";
		return *this;
	}
	if( cood == NULL || ind.cood == NULL){
		throw "δ������������";
		return *this;
	}

	//for (int i=0; i < N; i++)
	//{coor[i] = ind.coor[i];}
	//
	energy = ind.energy;
	for (int i=0; i < 4*N; i++)
	{cood[i] = ind.cood[i];}
	return *this;
}


void GA_Individual::Copy(GA_Individual *one,int N)
{
	int i;
	/*for(i=0;i<N;i++)
	{
		this->coor[i] =one->coor[i];
	}
*/
	for(i = 0; i < 4*N; i++)
	{
		this->cood[i] = one->cood[i];
	}
	this->energy = one->energy;
}

void GA_Individual::Copy1(GA_Individual *one,int N)
{
	int i;
	/*for(i=0;i<N;i++)
	{
		one->coor[i] = this->coor[i];
	}*/

	for(i = 0; i < 4*N; i++)
	{
		one->cood[i] = this->cood[i];
	}
	one->energy = this->energy;
}

geneticTool::geneticTool(void)
{

	cout<<"�Ŵ��ֲ��Ż���ʼ����ʼ"<<endl;
	//popsize = popsized;

	this->GApopsize = 100;

	//worst_pop.coor = calloc(N,sizeof(int));
	this->G = (GA_Individual *)malloc(this->GApopsize * sizeof(GA_Individual));
	this->A = (GA_Individual *)malloc(this->GApopsize * sizeof(GA_Individual));
	this->M = (GA_Individual *)malloc(this->GApopsize * sizeof(GA_Individual));
	for( int i = 0; i < this->GApopsize; i++)
	{	
		//this->G[i].coor = (double *)malloc(this->N * sizeof(double));
		this->G[i].cood = (double *)malloc(4*this->N * sizeof(double));
		///this->A[i].coor = (double *)malloc(this->N * sizeof(double));
		this->A[i].cood = (double *)malloc(4*this->N * sizeof(double));
		//this->M[i].coor = (double *)malloc(this->N * sizeof(double));
		this->M[i].cood = (double *)malloc(4*this->N * sizeof(double));
		//this->G[i].N = N;
		//this->A[i].N = N;
	}

	this->G->N=N;
	this->A->N=N;
	this->M->N=N;
}
	



geneticTool::~geneticTool(void)
{
	int i;
	cout<<"GA �ֲ��Ż�������"<<endl;
	for(i=0;i<this->GApopsize;i++){
	
		//if(G[i].coor != NULL)
		//{free(this->G[i].coor);}


		//if(A[i].coor != NULL)
		//{free(this->A[i].coor);}

		//if(M[i].coor != NULL)
		//{free(this->M[i].coor);}

		//if(G[i].cood != NULL)
		{free(this->G[i].cood);}
//
	//	if(A[i].cood != NULL)
		{free(this->A[i].cood);}

	//	if(M[i].cood != NULL)
		{free(this->M[i].cood);}

	}
	free(this->G);
	free(this->A);
	free(this->M);
	

}
//double GA_Individual::docking(double *cood,int N,int FeN)
//{
//	double E1=0;
//
//	return 0;
//}


geneticTool::geneticTool(int N,int FeN)
{
	cout<<"�Ŵ��ֲ��Ż���ʼ����ʼ"<<endl;
	//popsize = popsized;
	this->N = N;
	this->FeN = FeN;
	this->GApopsize = 100;
	this->best_pop.setN(N);
	this->pop.setN(N);
	this->pbest.setN(N);
	this->worst_pop.setN(N);
	this->P_E = PE_Tool::FS_E;
	this->P_F = PE_Tool::FS_F;
	//best_pop.coor = calloc(N,sizeof(int));
	//worst_pop.coor = calloc(N,sizeof(int));
	this->G = (GA_Individual *)malloc(this->GApopsize * sizeof(GA_Individual));
	this->A = (GA_Individual *)malloc(this->GApopsize * sizeof(GA_Individual));
	this->M = (GA_Individual *)malloc(this->GApopsize * sizeof(GA_Individual));
	for( int i = 0; i < this->GApopsize; i++)
	{	
		//this->G[i].coor = (double *)malloc(this->N * sizeof(double));
		this->G[i].cood = (double *)malloc(4*this->N * sizeof(double));
		//this->A[i].coor = (double *)malloc(this->N * sizeof(double));
		this->A[i].cood = (double *)malloc(4*this->N * sizeof(double));
		//this->M[i].coor = (double *)malloc(this->N * sizeof(double));
		this->M[i].cood = (double *)malloc(4*this->N * sizeof(double));
		//this->G[i].N = N;
		//this->A[i].N = N;
	}

	this->G->N=N;
	this->A->N=N;
	this->M->N=N;
}

void geneticTool::Genetic(double *cood,int N,int FeN)
{
	double E1 = 0,E2 = 0;
	int GArate = 10;//����������ʾ�Ŵ��ֲ��Ż����еĴ��� 
	this->N = N;
	this->FeN = FeN;
		
	this->P_E = PE_Tool::FS_E;
	this->P_F = PE_Tool::FS_F;
	int GApopsize = 10;
	int a[100]={0};
	int i,j,g,o=0,ren,temp;
	//double *R,*cood_new;
	//double *coor,*cood_L;
	//while(k--)
	
	//����ԭ�Ӽ�ľ��� ��Ϊ�̶�����ÿ��ֻ��Ҫ����һ��
	double *R;
	R = (double *)calloc(N*N,sizeof(double));
	Tool::Distance(cood,R,N);//�ŵ�ǰ������Լʱ��

	for (i=0;i<GApopsize;i++)
		{
			
	
			ren = N;
		///�������������
			for (g=0;g<N;g++)
			{
				a[g] = g+1;
				
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
					A[i].cood[g] = 1;
				}
				else 
				{
					A[i].cood[g] = 0;
				}
			}

		
		}
		for(i=0;i<N;i++)
		{
			A[0].cood[i] =cood[i];
		}

	//G��Ⱥ�ĳ�ʼ�����õ�G��Ⱥ���������ݣ������ҵ��������	
	for(i=0;i<GApopsize;i++){
		
			//A[i].cood = A[i].coor;
			for(int j=N;j<4*N;j++)
			{
				A[i].cood[j]=cood[j];
			}
			A[i].energy = Tool::Local(A[i].cood,N,P_E,P_F);
		//A[i].energy = GAdocking(A[i].coor,cood,A[i].cood,N,FeN,R);
		cout<<"��ʼ��"<<i<<"�������������"<<A[i].energy<<endl;
		}
	for(i=0;i<GApopsize;i++){
			
			this->M[i].Copy(&A[i],N);
		}

	this->best_pop.Copy(&M[0],N);
	this->worst_pop.Copy(&M[0],N);
	for(i=0;i<GApopsize;i++)
		{

			if(M[i].energy<best_pop.energy)
				this->best_pop.Copy(&M[i],N);

			if(M[i].energy>worst_pop.energy)
				this->worst_pop.Copy(&M[i],N);
		}
	cout<<"��õ�ֵ��"<<best_pop.energy<<endl;
	cout<<"����ֵ��"<<worst_pop.energy<<endl;
	
	while(GArate--)
	{	
		//cout<<"��0�������������"<<A[0].energy<<endl;
		//E0 = best_pop.energy;
		select_operator(cood,best_pop.energy,worst_pop.energy);//������õ�cood����������������Ӧ����Ӧ�ĸ���
		//cout<<"��N�������������"<<A[6].energy<<endl;
		////����
		crossover_operator();//���õ�cood�����ݵ����õ�ͬ����������
		//cout<<"��N�������������"<<A[6].energy<<endl;
		///����
		mutation_operator();
		//cout<<"��N�������������"<<A[6].energy<<endl;
		///�������Ÿ���
		
		//adjustment_operator(cood);
		//
		//this->best_pop.Copy(&A[0],N);
		//copy(&A[0],&best_pop,N);//�õ����ŵĸ���
		for(i=0;i<GApopsize;i++){
			//A[i].cood = A[i].coor;*/
		/*for(int j=N;j<4*N;j++)
			{
				A[i].cood[j]=cood[j];
			}
			*/
			A[i].energy = Tool::Local(A[i].cood,N,P_E,P_F);
			//A[i].energy = GAdocking(A[i].coor,cood,A[i].cood,N,FeN,R);//A[i].coodû�б仯
			//cout<<"�����仯֮��ĵ�"<<i<<"�������������"<<A[i].energy<<endl;
		}

		for(i=0;i<GApopsize;i++)
		{
			if(M[i].energy>A[i].energy)
				this->M[i].Copy(&A[i],N);
				//this->A[i].Copy(&M[i],N);
			//cout<<"�����������"<<i<<"����������"<<M[i].energy<<endl;
		}

		//this->best_pop.Copy(&M[0],N);
		//this->worst_pop.Copy(&M[0],N);
		//copy(&best_pop,&(A[0]),N);
		//copy(&worst_pop,&(A[0]),N);
		
		for(i=0;i<GApopsize;i++)
		{
			//this->G[i].Copy(A[i],N);

			if(M[i].energy<best_pop.energy)
				this->best_pop.Copy(&M[i],N);
				
			if(M[i].energy>worst_pop.energy)
				this->worst_pop.Copy(&M[i],N);
		}
		
		E1 = best_pop.energy;
		E2 = worst_pop.energy;
		//this->worst_pop.docking(cood,N,FeN);

		//E1 =  Tool::Local(best_pop.coor,N,P_E,P_F);//ʹ��cood����ͱ���Ҫ�ĵ�֮ǰ��*cood������
		//E2 =  Tool::Local(worst_pop.coor,N,P_E,P_F);
		cout<<endl;
		//cout<<"���ŵĽ���ǣ�"<<E1<<endl;
		//cout<<"���Ľ���ǣ�"<<E2<<endl;

		//cout<<"��"<<20-GArate<<"��GA�ֲ��Ż������ŵĸ������"<<best_pop.energy<<endl;
		//cout<<"��"<<20-GArate<<"��GA�ֲ��Ż������ĸ������"<<worst_pop.energy<<endl;
		//cout<<"��"<<20-GArate<<"��GA�ֲ��Ż���diliu�������"<<A[6].energy<<endl;

		/*for(i=0;i<0.2*GApopsize;i++){
			
			this->A[i].Copy(&M[i],N);
		}*/
	}	
	
	this->pbest.Copy(&best_pop,N);//pbest�ĳ�ʼ
	//this->pop.Copy(&best_pop,N);//pop�ĳ�ʼ��
	//E1 = best_pop.energy;
	pbest.energy = Tool::Local(pbest.cood,N,P_E,P_F);
	//pbest.energy = GAdocking(pbest.coor,cood,pbest.cood,N,FeN,R);
	cout<<endl;
	cout<<"shuchu �õ����Ÿ�����"<<pbest.energy<<endl;

	//pbest.cood[i] = 
	//�����Ÿ����ֵ����cood

	//�����Ÿ��崫��ȥ�ĺ���
	for(int h=0;h<4*N;h++)
	{cood[h] = pbest.cood[h];  }


				

	

}

/// ѡ������
void geneticTool::select_operator(double *cood,double BestEnergy,double WorstEnergy){
	int i,index;
	double p,sum=0;
	double *cvalue;

	cvalue = (double *)calloc(GApopsize,sizeof(double));
	

	for(i=0;i<GApopsize;i++)
		A[i].value=exp(-3*(A[i].energy-BestEnergy)/(WorstEnergy-BestEnergy+0.0001));;

	for(i=0;i<GApopsize;i++)
		cvalue[i] = A[i].value/sum;

	for(i=1;i<GApopsize;i++)
		cvalue[i] = cvalue[i-1]+cvalue[i];
	cout<<"��õ�ֵ��"<<BestEnergy<<endl;
	cout<<"����ֵ��"<<WorstEnergy<<endl;
	for(i=0;i<GApopsize;i++)
	{
		//cout<<"���Ե�ֵ"<<A[i].value<<endl;

		p = rand()%1000/1000.0;
		index = 0;
		while(p>cvalue[index])
			index++;
		this->G[i].Copy(&A[index],N);
		
	}

	for(i=0;i<GApopsize;i++)
		this->A[i].Copy(&G[i],N);
	

	free(cvalue);
	
}

///��������
void geneticTool::crossover_operator(){
	int i,j;
	int length;
	int *index;
	int dNum[2];
	int point,point1,point2,point_left,point_right,temp;
	double p;
	double pc = 0.4;//��Ҫ����Ч��ȥ�޸�
	double ch;
	
	index = (int *)calloc(GApopsize,sizeof(int));


	/////����õ�������Ӵ�
	for(i=0;i<GApopsize;i++)
		index[i] = i;
	for(i=0;i<GApopsize;i++)
	{
		point = rand()%(GApopsize-i);
		temp = index[i];
		index[i] = index[point+i];
		index[point+i] = temp;
		//cout<<"���ѡ����Ӵ����ǣ�"<<index[i]<<endl;

	}
	
	//ѡ�񽻲�ڵ�

	for(i=0;i<GApopsize-1;i+=2)
	{
		p = rand()%1000/1000.0;
		
		if(p<pc)
		{

			dNum[0] = 0; dNum[1] = 0;
			point_left = (int)(N*rand()/(RAND_MAX+1.0));
			point_right = (int)(N*rand()/(RAND_MAX+1.0));
			
			if(point_right<point_left)//�������ҳ���
			{
				point = point_left;
				point_left = point_right;
				point_right = point;
			}
			//�����жϸ���Ԫ�ص�ռ�ȸ���
			for(j=point_left;j<=point_right;j++)
			{
				(A[index[i]].cood[j] == 0)?(dNum[0]--):(dNum[0]);
				(A[index[i]].cood[j] == 1)?(dNum[1]--):(dNum[1]);
			
				(A[index[i+1]].cood[j] == 0)?(dNum[0]++):(dNum[0]);
				(A[index[i+1]].cood[j] == 1)?(dNum[1]++):(dNum[1]);//a�ж�A[index[i]]��A[index[i+1]]�Ķ�Ӧλ�õ�coor�Ƿ�һ��

		/////////////////////////////////////////////////////////////���һ�£���ôΪ0ʱdNum[0]Ϊ0��Ϊ1ʱdNum[1]Ϊ0�������һ�£���ôΪ0ʱdNum[0]Ϊ1��Ϊ1ʱdNum[1]Ϊ1����ʱ��Ϊ0������һ�����ж�
						///////////////////////////////////////////////�õ����  �������������01��������ƥ�����Ҫ���д���
				ch = A[index[i]].cood[j];
				A[index[i]].cood[j] = A[index[i+1]].cood[j];
				A[index[i+1]].cood[j] = ch;//////////////////////////����A[index[i]]��A[index[i+1]]��ѡ����ڵ��Ų�������õ���ͬ�Ӵ���
			}
			
			length = N - (point_right - point_left + 1);
			for(j=0;j<2;j++)
			{	
				while(dNum[j]>0)
				{
					point = rand()%length;////�����ȥ�����������ڵ�������������ȥ����  ����������
					if(point>=point_left)
						point1 = point_right+(point+1-point_left);
					else
						point1 = point;

					while(A[index[i]].cood[point1] != j)//
					{
						point = (point+1)%length;
						if(point>=point_left)
							point1 = point_right+(point+1-point_left);
						else
							point1 = point;
					}

					point = rand()%length;
					if(point>=point_left)
						point2 = point_right+(point+1-point_left);
					else
						point2 = point;

					temp = A[index[i+1]].cood[point2];
					while(dNum[temp]>=0)
					{
						point = (point+1)%length;
						if(point>=point_left)
							point2 = point_right+(point+1-point_left);
						else
							point2 = point;

						temp = A[index[i+1]].cood[point2]; 
					}

					A[index[i]].cood[point1] = temp;
					A[index[i+1]].cood[point2] = j;

					dNum[j]--;
					dNum[temp]++;
				}
			}
		}
	}

	free(index);
}

////��������
void geneticTool::mutation_operator(){
	int i,j;
	int point;
	double ch;
	double p;
	double pm=0.1;

	for(i=0;i<GApopsize;i++)
		for(j=0;j<N;j++)
		{
			p = rand()%1000/1000.0;
			if(p<pm)
			{
				ch = A[i].cood[j];
				point = (int)(N*rand()/(RAND_MAX+1.0));
				while(A[i].cood[point] == ch)
				{
					point++;
					if(point==N)
					point = 0;
				}	
				A[i].cood[j] = A[i].cood[point];
				A[i].cood[point] = ch;	
			}
		}
}


void geneticTool::adjustment_operator(double *cood,double *R)
{
	int i,randN1,randN2,tempcoor;
	double tempE,r,rate;
	
	for(i=0;i<this->GApopsize;i++)
	{
		r = RAND1;
		if(r < 1.0)
		{
			randN1 = (int)(N * RAND1);
			randN2 = (int)(N * RAND1);
			while(pbest.cood[randN1] == pbest.cood[randN2])
			{
				randN1 = (int)(N * RAND1);
				randN2 = (int)(N * RAND1);
			}
			tempcoor = pbest.cood[randN1];
			pbest.cood[randN1] = pbest.cood[randN2];
			pbest.cood[randN2] = tempcoor;

			//tempE = GAdocking(pbest.coor,cood,pbest.cood,N,FeN,R);
		//	tempE = Tool::Local(pbest.coor,N,PE_E,PE_F);
			if(tempE <= pbest.energy)
			{
				pbest.energy = tempE;
				continue;
			} 
			tempcoor = pbest.cood[randN1];
			pbest.cood[randN1] = pbest.cood[randN2];
			pbest.cood[randN2] = tempcoor;
		}
	}
}


//void geneticTool::changeIndividual(GA_Individual one,GA_Individual two,int N){
//	int i;
//	one.energy=two.energy;
//	for(i=0;i<N;i++){
//		one.coor[i] = two.coor[i];
//	}
//}

double geneticTool::GAdocking(double *coorL,double *cood,double *coodL,int N,int FeN,double *R)//�����Խ�coor��cood���صõ�����  ����Tool����Local��3.8������д����ԣ�
{
	double E0;
	//double E0,*coor_L,*cood_L,*temp,*R;
	//int	m;
	int i;
	E0 = 0;

	//coor_L = cood;
	//cood_L = cood + N;

	/*coor_L = temp;
	temp = coor;
	coor = coor_L;*/
	coodL = coorL;
	/*for(i=0;i<N;i++)
	{
		coodL[i] = coorL[i];
	}*/
	for(i=N;i<4*N;i++){
	coodL[i]=cood[i];
	}

	/*coodL = coorL;
	coodL+N = cood+N;*/
	E0= Tool::Local(coodL,N,P_E,P_F);
	//E0 = PE_Tool::FS_E(cood,R,N);
	cout<<"�����ǣ�"<<E0;
	return E0;
}

void geneticTool::GetBestcood(double *coorL,double *cood,int N)
{
	cood = coorL;
	/*int i;
	for(i=0;i<4*N;i++)
	{
		cood[i] = coorL[i];
	}*/
}
