#include "geneticTool.h"
#include"Tool.h"
#include "PE_Tool.h"




GA_Individual* GA_Individual::Individual(int N)
{
	GA_Individual *obj = (GA_Individual*)malloc(sizeof(GA_Individual));
	obj->coor = (double *)malloc(N * sizeof(double));
	obj->cood = (double *)malloc(4*N * sizeof(double));
	return obj;
}

GA_Individual::GA_Individual(void){
	cout<<"GA-ind"<<endl;
	this->coor = NULL;
	this->cood = NULL;
}

GA_Individual::GA_Individual(int N){
	cout<<"GA-ind"<<endl;
	this->setN(N);
}

GA_Individual::~GA_Individual(void){

	if(this->coor != NULL)
	{
		free(this->coor);
	}

	if(this->cood != NULL)

	{free(this->cood);}
}

void GA_Individual::setN(int N)
{
	this->N = N;
	this->coor = (double *)malloc(N * sizeof(double));
	this->cood = (double *)malloc(4*N * sizeof(double));
}
 
GA_Individual& GA_Individual::operator=(GA_Individual &ind)
{
	if (N != ind.N && N != 0)
	{
		throw "GA个体原子数不对等";
		return *this;
	}
	if( cood == NULL || ind.cood == NULL){
		throw "未开辟坐标数组";
		return *this;
	}

	for (int i=0; i < N; i++)
	{coor[i] = ind.coor[i];}
	
	energy = ind.energy;
	for (int i=0; i < 4*N; i++)
	{cood[i] = ind.cood[i];}
	return *this;
}


void GA_Individual::Copy(GA_Individual *one,int N)
{
	int i;
	for(i=0;i<N;i++)
	{
		this->coor[i] =one->coor[i];
	}

	for(i = 0; i < 4*N; i++)
	{
		this->cood[i] = one->cood[i];
	}
	this->energy = one->energy;
}

void GA_Individual::Copy1(GA_Individual *one,int N)
{
	int i;
	for(i=0;i<N;i++)
	{
		one->coor[i] = this->coor[i];
	}

	for(i = 0; i < 4*N; i++)
	{
		one->cood[i] = this->cood[i];
	}
	one->energy = this->energy;
}

geneticTool::geneticTool(void)
{

	cout<<"遗传局部优化初始化开始"<<endl;
	//popsize = popsized;

	this->GApopsize = 100;

	//worst_pop.coor = calloc(N,sizeof(int));
	this->G = (GA_Individual *)malloc(this->GApopsize * sizeof(GA_Individual));
	this->A = (GA_Individual *)malloc(this->GApopsize * sizeof(GA_Individual));
	this->M = (GA_Individual *)malloc(this->GApopsize * sizeof(GA_Individual));
	for( int i = 0; i < this->GApopsize; i++)
	{	
		this->G[i].coor = (double *)malloc(this->N * sizeof(double));
		this->G[i].cood = (double *)malloc(4*this->N * sizeof(double));
		this->A[i].coor = (double *)malloc(this->N * sizeof(double));
		this->A[i].cood = (double *)malloc(4*this->N * sizeof(double));
		this->M[i].coor = (double *)malloc(this->N * sizeof(double));
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
	cout<<"GA 局部优化的析构"<<endl;
	for(i=0;i<this->GApopsize;i++){
	
		//if(G[i].coor != NULL)
		{free(this->G[i].coor);}


		//if(A[i].coor != NULL)
		{free(this->A[i].coor);}

		//if(M[i].coor != NULL)
		{free(this->M[i].coor);}

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
	cout<<"遗传局部优化初始化开始"<<endl;
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
		this->G[i].coor = (double *)malloc(this->N * sizeof(double));
		this->G[i].cood = (double *)malloc(4*this->N * sizeof(double));
		this->A[i].coor = (double *)malloc(this->N * sizeof(double));
		this->A[i].cood = (double *)malloc(4*this->N * sizeof(double));
		this->M[i].coor = (double *)malloc(this->N * sizeof(double));
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
	int GArate = 10;//这里用来表示遗传局部优化进行的次数 
	this->N = N;
	this->FeN = FeN;
		
	this->P_E = PE_Tool::FS_E;
	this->P_F = PE_Tool::FS_F;
	int GApopsize = 50;
	int a[100]={0};
	int i,j,g,o=0,ren,temp;
	//double *R,*cood_new;
	//double *coor,*cood_L;
	//while(k--)
	
	//计算原子间的距离 此为固定数据每次只需要计算一次
	double *R;
	R = (double *)calloc(N*N,sizeof(double));
	Tool::Distance(cood,R,N);//放到前面会跟节约时间

	for (i=0;i<GApopsize;i++)
		{
			
	
			ren = N;
		///产生随机的数组
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
					A[i].coor[g] = 1;
				}
				else 
				{
					A[i].coor[g] = 0;
				}
			}

		
		}
		for(i=0;i<N;i++)
		{
			A[0].coor[i] =cood[i];
		}

	//G种群的初始化，得到G种群的所有数据，并且找到最优书记	
	for(i=0;i<GApopsize;i++){
		
			A[i].cood = A[i].coor;
			for(int j=N;j<4*N;j++)
			{
				A[i].cood[j]=cood[j];
			}
			A[i].energy = Tool::Local(A[i].cood,N,P_E,P_F);
		//A[i].energy = GAdocking(A[i].coor,cood,A[i].cood,N,FeN,R);
		cout<<"开始第"<<i<<"个个体的能量是"<<A[i].energy<<endl;
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
	cout<<"最好的值是"<<best_pop.energy<<endl;
	cout<<"最差的值是"<<worst_pop.energy<<endl;
	
	while(GArate--)
	{	
		//cout<<"第0个个体的能量是"<<A[0].energy<<endl;
		//E0 = best_pop.energy;
		select_operator(cood,best_pop.energy,worst_pop.energy);//这里会用到cood的内容用来计算相应的适应的概率
		//cout<<"第N个个体的能量是"<<A[6].energy<<endl;
		////交叉
		crossover_operator();//不用到cood的内容调整得到同比例的算子
		//cout<<"第N个个体的能量是"<<A[6].energy<<endl;
		///变异
		mutation_operator();
		//cout<<"第N个个体的能量是"<<A[6].energy<<endl;
		///保留最优个体
		
		//adjustment_operator(cood);
		//
		//this->best_pop.Copy(&A[0],N);
		//copy(&A[0],&best_pop,N);//得到最优的个体
		for(i=0;i<GApopsize;i++){
			A[i].cood = A[i].coor;
		for(int j=N;j<4*N;j++)
			{
				A[i].cood[j]=cood[j];
			}
			
			A[i].energy = Tool::Local(A[i].cood,N,P_E,P_F);
			//A[i].energy = GAdocking(A[i].coor,cood,A[i].cood,N,FeN,R);//A[i].cood没有变化
			cout<<"经过变化之后的第"<<i<<"个个体的能量是"<<A[i].energy<<endl;
		}

		for(i=0;i<GApopsize;i++)
		{
			if(M[i].energy>A[i].energy)
				this->M[i].Copy(&A[i],N);
				//this->A[i].Copy(&M[i],N);
			cout<<"保存的优良第"<<i<<"个的能量是"<<M[i].energy<<endl;
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

		//E1 =  Tool::Local(best_pop.coor,N,P_E,P_F);//使用cood计算就必须要的到之前的*cood的数据
		//E2 =  Tool::Local(worst_pop.coor,N,P_E,P_F);
		cout<<endl;
		//cout<<"最优的结果是："<<E1<<endl;
		//cout<<"最差的结果是："<<E2<<endl;

		cout<<"第"<<20-GArate<<"次GA局部优化的最优的个体输出"<<best_pop.energy<<endl;
		cout<<"第"<<20-GArate<<"次GA局部优化的最差的个体输出"<<worst_pop.energy<<endl;
		cout<<"第"<<20-GArate<<"次GA局部优化的diliu个体输出"<<M[6].energy<<endl;

		/*for(i=0;i<0.2*GApopsize;i++){
			
			this->A[i].Copy(&M[i],N);
		}*/
	}	
	
	this->pbest.Copy(&best_pop,N);//pbest的初始
	//this->pop.Copy(&best_pop,N);//pop的初始化
	//E1 = best_pop.energy;

	//pbest.energy = GAdocking(pbest.coor,cood,pbest.cood,N,FeN,R);
	cout<<endl;
	cout<<"得到最优个体是"<<pbest.energy<<endl;

	//pbest.cood[i] = 
	//将最优个体的值传给cood
	//将最优个体传出去的函数
	cood = pbest.cood;               ////???????有问题下午改？？？？？？/////
	//GetBestcood(pbest.cood,cood,N);
	

}

/// 选择算子
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
	cout<<"最好的值是"<<BestEnergy<<endl;
	cout<<"最差的值是"<<WorstEnergy<<endl;
	for(i=0;i<GApopsize;i++)
	{
		//cout<<"测试的值"<<A[i].value<<endl;

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

///交叉算子
void geneticTool::crossover_operator(){
	int i,j;
	int length;
	int *index;
	int dNum[2];
	int point,point1,point2,point_left,point_right,temp;
	double p;
	double pc = 0.4;//需要根据效率去修改
	double ch;
	
	index = (int *)calloc(GApopsize,sizeof(int));


	/////随机得到交叉的子代
	for(i=0;i<GApopsize;i++)
		index[i] = i;
	for(i=0;i<GApopsize;i++)
	{
		point = rand()%(GApopsize-i);
		temp = index[i];
		index[i] = index[point+i];
		index[point+i] = temp;
		//cout<<"随机选择的子代数是："<<index[i]<<endl;

	}
	
	//选择交叉节点

	for(i=0;i<GApopsize-1;i+=2)
	{
		p = rand()%1000/1000.0;
		
		if(p<pc)
		{

			dNum[0] = 0; dNum[1] = 0;
			point_left = (int)(N*rand()/(RAND_MAX+1.0));
			point_right = (int)(N*rand()/(RAND_MAX+1.0));
			
			if(point_right<point_left)//保持左右成立
			{
				point = point_left;
				point_left = point_right;
				point_right = point;
			}
			//用来判断各种元素的占比个数
			for(j=point_left;j<=point_right;j++)
			{
				(A[index[i]].coor[j] == 0)?(dNum[0]--):(dNum[0]);
				(A[index[i]].coor[j] == 1)?(dNum[1]--):(dNum[1]);
			
				(A[index[i+1]].coor[j] == 0)?(dNum[0]++):(dNum[0]);
				(A[index[i+1]].coor[j] == 1)?(dNum[1]++):(dNum[1]);//a判断A[index[i]]与A[index[i+1]]的对应位置的coor是否一致

		/////////////////////////////////////////////////////////////如果一致，那么为0时dNum[0]为0，为1时dNum[1]为0；如果不一致，那么为0时dNum[0]为1，为1时dNum[1]为1，此时不为0将做下一步的判断
						///////////////////////////////////////////////得到结果  如果在区间里面01个数不是匹配的着要进行处理
				ch = A[index[i]].coor[j];
				A[index[i]].coor[j] = A[index[i+1]].coor[j];
				A[index[i+1]].coor[j] = ch;//////////////////////////互换A[index[i]]与A[index[i+1]]在选择点内的排布情况，得到不同子代；
			}
			
			length = N - (point_right - point_left + 1);
			for(j=0;j<2;j++)
			{	
				while(dNum[j]>0)
				{
					point = rand()%length;////随机在去除交叉区间内的两边区域里面去交换  来调整比例
					if(point>=point_left)
						point1 = point_right+(point+1-point_left);
					else
						point1 = point;

					while(A[index[i]].coor[point1] != j)//
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

					temp = A[index[i+1]].coor[point2];
					while(dNum[temp]>=0)
					{
						point = (point+1)%length;
						if(point>=point_left)
							point2 = point_right+(point+1-point_left);
						else
							point2 = point;

						temp = A[index[i+1]].coor[point2]; 
					}

					A[index[i]].coor[point1] = temp;
					A[index[i+1]].coor[point2] = j;

					dNum[j]--;
					dNum[temp]++;
				}
			}
		}
	}

	free(index);
}

////变异算子
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
				ch = A[i].coor[j];
				point = (int)(N*rand()/(RAND_MAX+1.0));
				while(A[i].coor[point] == ch)
				{
					point++;
					if(point==N)
					point = 0;
				}	
				A[i].coor[j] = A[i].coor[point];
				A[i].coor[point] = ch;	
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
			while(pbest.coor[randN1] == pbest.coor[randN2])
			{
				randN1 = (int)(N * RAND1);
				randN2 = (int)(N * RAND1);
			}
			tempcoor = pbest.coor[randN1];
			pbest.coor[randN1] = pbest.coor[randN2];
			pbest.coor[randN2] = tempcoor;

			tempE = GAdocking(pbest.coor,cood,pbest.cood,N,FeN,R);
		//	tempE = Tool::Local(pbest.coor,N,PE_E,PE_F);
			if(tempE <= pbest.energy)
			{
				pbest.energy = tempE;
				continue;
			} 
			tempcoor = pbest.coor[randN1];
			pbest.coor[randN1] = pbest.coor[randN2];
			pbest.coor[randN2] = tempcoor;
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

double geneticTool::GAdocking(double *coorL,double *cood,double *coodL,int N,int FeN,double *R)//用来对接coor与cood返回得到能量  调用Tool：：Local，3.8晚，明天写完测试；
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
	cout<<"能量是："<<E0;
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
