//����дBH�㷨��������




#include "BasinHopping.h"
#include "BH_Manage.h"
#include "LocalTool.h"
#include "Tool.h"
#include "PE_Tool.h"
#include "LocalTool.h"
#include "geneticTool.h"
#include<ctime>


BH_Individual* BH_Individual::Individual(int N)
{
	BH_Individual *obj = (BH_Individual*)malloc(sizeof(BH_Individual));
	obj->cood = (double *)malloc(4 * N * sizeof(double));
	return obj;
}

BH_Individual::BH_Individual(void){
	this->cood = NULL;
}

BH_Individual::BH_Individual(int N){
	this->set4N(N);
}

BH_Individual::~BH_Individual(void){

	if(this->cood != NULL)
		free(this->cood);
}

void BH_Individual::set4N(int N)
{
	this->N = N;
	this->cood = (double *)malloc(4 * N * sizeof(double));
}

BH_Individual& BH_Individual::operator=(BH_Individual &ind)
{
	if (N != ind.N && N != 0)
	{
		throw "BH����ԭ�������Ե�";
		return *this;
	}
	if( cood == NULL || ind.cood == NULL){
		throw "δ������������";
		return *this;
	}

	for (int i=0; i < 4*N; i++)
	cood[i] = ind.cood[i];
	energy = ind.energy;

	return *this;
}

void BH_Individual::Copy(BH_Individual *bh, int N)
{
	int i;
	for(i = 0; i < 4 * N; i++)
	{
		this->cood[i] = bh->cood[i];
	}
	this->energy = bh->energy;
}


BasinHopping::BasinHopping(int FeN,int N,int popsized,double E)
	// Params:
		/*
			FeN: ��ԭ������
			N����ԭ����
			popsized: BH�㷨��Ⱥ��С�趨
			E����������Ԥ��ֵ
		*/
{
	cout<<"--------------------BH�㷨��ʼ������----------------------"<<endl;
	this->N = N;
	this->current.set4N(N);
	this->best.set4N(N);
	this->change.set4N(N);
	this->prebest.set4N(N);
	this->readbest.set4N(N);
	this->GAbest.set4N(N);
	//initDate(N);
	//�����ʼ������
	int i;
	cout<<"------------------����ʼ--------------------------------"<<endl;
	cout<<"Total Atom Number is��"<< N<<endl;
	cout<<"Ir Atom Number is:"<<FeN  << endl;
	cout << "Pt Atom Number is:" << N - FeN << endl;
	this->popsize = popsized;
	this->Ebest=E;
	this->FeN=FeN;
	
	/*
		
		�ƺ���ѡȡ��Gupta���ܺ���

	*/

	this->P_E = PE_Tool::FS_E;
	this->P_F = PE_Tool::FS_F;
	//this->P_Enew = PE_Tool::FS_Enew;
	//this->P_Fnew = PE_Tool::FS_Fnew;
	this->X = (BH_Individual *)malloc(this->popsize * sizeof(BH_Individual));
	this->Y = (BH_Individual *)malloc(this->popsize * sizeof(BH_Individual));
	this->Z = (BH_Individual *)malloc(this->popsize * sizeof(BH_Individual));
	this->best.cood = (double *)malloc(4*this->N * sizeof(double));
	for( i = 0; i < this->popsize; i++)
	{	
		this->X[i].cood = (double *)malloc(4*this->N * sizeof(double));
		this->Y[i].cood = (double *)malloc(4*this->N * sizeof(double));
		this->Z[i].cood = (double *)malloc(4*this->N * sizeof(double));
	}
	this->X->N=N;
	this->Y->N=N;
	this->Z->N=N;
	//
	//this->Initialization(N,popsize);
}

BasinHopping::BasinHopping(){
	
	int i;
	cout<<"N ="<<N<<endl;
	//this->popsize = 50;
	cout<<"popsize = "<<popsize<<endl;
	/*this->P_E = PE_Tool::FS_E;
	this->P_F = PE_Tool::FS_F;*/
	/*this->P_E = PE_Tool::FS_Enew;
	this->P_F = PE_Tool::FS_Fnew;*/
	this->X = (BH_Individual *)malloc(this->popsize * sizeof(BH_Individual));
	this->Y = (BH_Individual *)malloc(this->popsize * sizeof(BH_Individual));
	this->Z = (BH_Individual *)malloc(this->popsize * sizeof(BH_Individual));
	this->best.cood = (double *)malloc(4*this->N * sizeof(double));
	for( i = 0; i < this->popsize; i++)
	{	
		this->X[i].cood = (double *)malloc(4*this->N * sizeof(double));
		this->Y[i].cood = (double *)malloc(4*this->N * sizeof(double));
		this->Z[i].cood = (double *)malloc(4*this->N * sizeof(double));
	}
	//this->Initialization(N,popsize);
	
}

BasinHopping::~BasinHopping()
{
	int i;
	cout<<"BasinHopping����"<<endl;
	for(i=0;i<this->popsize;i++){
	
		free(this->X[i].cood); 
		free(this->Y[i].cood);
		free(this->Z[i].cood);

	}
	free(this->X);
	free(this->Y);
	free(this->Z);

}

void BasinHopping::initDate(int N){
	this->N=N;
	this->P_E = PE_Tool::FS_E;
	this->P_F = PE_Tool::FS_F;
	kIndex = 0;
	double r0= 2.75;

	for(int i=N; i<4 * N; i++)
		// ���һ����ʼ�ṹ
		current.cood[i] = (RAND-0.5) * r0 * pow((double)N,1.0/3);

	BH_Individual temp(N);
	// ��¼δ���оֲ���ǰ������������Ϣ
	temp = current;
	current.energy = Tool::Local(temp.cood,N,P_E,P_F);   //����Ѿ�����˳�ʼ״̬�µľֲ��Ż���Сֵ��
	best = current;
	prebest = best;

	// ͨ����½����仯�������
	cout<<"��ǰ�ֲ��Ż����������"<<current.energy<<endl;
}


//��Ⱥ��ʼ���Ĺ��̣�

void BasinHopping::Initialization()
{
	int i,j,g,o=0,ren,temp;
	int a[100];
	this->Localtimes = 0;
	this->Localsuccesstimes=0;
	this->S_times = 0;
	//static int g = 0;
	//static int o = 0;
	//double *R = new double[N*N];
	int rate1=10;
	int rate2=20;
	double r0= 2.75;
	BH_Individual *best1;
	best1 = this->Z;


	
	if (Tool::readDiamond(readbest.cood,N,0))
		{
			readbest.energy = Tool::Local(readbest.cood,N,P_E,P_F);
			cout<<readbest.energy<<endl;
		}

	
	for (i=0;i < popsize;i++)
	/*
		
		PdIr38�ṹ��ʼ��

		��Ⱥ30�����䣺

			1-10������Pt���Žṹ��������ã���������������
			10-20���������������
			20-30���� ��Fe���Žṹ������ã���������������
	
	*/
	{
		if(i >= rate2)
		{
			for (j=0;j<3*N;j++)
			{
			Z[i].cood[j+N] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
			}
		}
		else if(i < rate2 && i >= rate1)
		{
			
			if(Tool::readDiamond1(Z[i].cood,N,0))
			{
				cout<<"����Ir�����Žṹ�Ѵ��ڣ�ѡ�����Žṹ��Ϊ��ʼ���͡�"<<endl;
			}
			else
			{
				for (j=0;j<3*N;j++)
				{
			
					Z[i].cood[j+N] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
				}
			}
			
		
		}
		else
		{
			
			if(Tool::readDiamond2(Z[i].cood,N,0))
			{
				cout<<"����Pd�����Žṹ�Ѵ��ڣ�ѡ�����Žṹ��Ϊ��ʼ���͡�"<<endl;
			}
			else
			{
				for (j=0;j<3*N;j++)
				{
			
					Z[i].cood[j+N] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
				}
			}
			
		}

		ren = N;
		//�������������
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
				Z[i].cood[g] = 1;
			}
			else 
			{
				Z[i].cood[g] = 0;
			}
		}




		
		//X[i].energy = Tool::Local(X[i].cood,N,P_E,P_F);
		//Y[i].energy = Tool::Local(Y[i].cood,N,P_E,P_F);
		//Z[i].energy = PE_Tool::FS_E(Z[i].cood,R, N);

		// �Խṹ������½����ֲ��Ż�

		Z[i].energy = Tool::Local(Z[i].cood,N,P_E,P_F);
		//PE_Tool::FS_E(double *cood,double *R,int N)
		//Z[i].energy = Tool::Local(Z[i].cood,N,P_Enew,P_Fnew);
		
		this->X[i].Copy(&Z[i],this->N);

		if(Z[i].energy < best1->energy)
		best1 = this->Z+i;
		
		//best.energy=best1->energy;
		//cout<<"��ʼ��ȺX�е�"<<i+1<<"�����������"<<Z[i].energy<<endl;
		//this->Y[i]=this->X[i];
	}
	//cout<<"�õ�"<<popsize<<"����Ⱥ�е�����ֵ����"<<best1->energy<<endl;
	cout << endl;
	/*for(i=0;i<popsize;i++)
	{
	o++;
	cout<<"������ȺY�е�"<<o<<"�����������"<<Y[i].energy<<endl;

	}*/
}

/*void BasinHopping::rechange()
{
	int i=0,h=0,m=0;
	h=rand()%N + 1;
	m=rand()%popsize + 1;

	for (i=0;i<popsize;i++)
	{
		
	}

}*/


void BasinHopping::start(){
	successRate = 0;
	int success = 0;
	int all = 0;
	


}


void BasinHopping::disturbance()//�Ŷ���������
{
	//_sleep(0.1*1000);
	srand((unsigned)time(NULL));
	int i,a[3];
	static int k=0;
	double s=0.35;
	double *lx,*ly,*lz,*cx,*cy,*cz;
	change=current;
	cx = current.cood;
	cy = current.cood + N;
	cz = current.cood + 2*N;
	lx = change.cood;
	ly = change.cood + N;
	lz = change.cood + 2*N;

  

	    for(i = 0; i < N; i++){
		
		for(int j=0;j<3;j++)
		{
			
			a[j]=rand()%3-1;

			//cout<<a[j]<<'\t';
			
			
		}
			lx[i] = cx[i] + a[0]*s;
			ly[i] = cy[i] + a[1]*s;
			lz[i] = cz[i] + a[2]*s;
			change.cood[i] = lx[i];
			change.cood[i + N] = ly[i];
			change.cood[i + 2*N] = lz[i];
			//X[u].cood[i] = ly[i];
			//X[u].cood[i + N] = ly[i];
			//X[u].cood[i + 2*N] = ly[i];

 		
	//}
	}
	
	//cout<<a[0]<<'\t'<<a[1]<<'\t'<<a[2]<<endl;
	change.energy = Tool::Local(change.cood,N,P_E,P_F);//�Ŷ��ƺ�û��Ч�� ��һ���Ŷ��Ľ���͵ڶ����Ŷ����һ��
	++k;
	cout<<"��"<<k<<"���Ŷ���������"<<change.energy<<endl;
	//current=change;


}

void BasinHopping::remake()
{
	srand((unsigned)time(NULL));
	int i,j,l;
	double b[300];
	//BH_Individual *best2;
	const int M(100);
	double *a;
	double Max_dis=0;
	double s=0.45;
	//best2 = this->Z;
	static int Kate=0;
	a = (double *)calloc(N,sizeof(double));
	//static int Localtimes=0;
	//this->Localtimes = 0;
	for(i=0;i<popsize;i++)
	{	
		this->Y[i].Copy(&Z[i],this->N);
		//Tool:: distool(Z[i].cood,a,N,Max_dis);
		for(l=0;l<3*N;l++)
			{
				
				//b[l]=(rand()%21-10)/10.0;//���ǡ�-1,1��֮����0.1�ı��
				b[l]=rand()%2-1;//����ȡ��-1,1���������ı��
				
				if(b[l] == 0)
					b[l] = 1;
				//cout<<b[l]<<'\t';
			}
		
		for(j=0;j<N;j++)
		{
			
			Z[i].cood[j + N] += b[j]*s;
			Z[i].cood[j + 2*N] += b[j + N]*s;
			Z[i].cood[j + 3*N] += b[j + 2*N]*s;
			
			/*Z[i].cood[j + N] += b[j]*s*a[j];
			Z[i].cood[j + 2*N] += b[j + N]*s*a[j];
			Z[i].cood[j + 3*N] += b[j + 2*N]*s*a[j];*/
				//cout<<b[0]<<'\t'<<b[1]<<'\t'<<b[2]<<endl;
		}
	
		Z[i].energy = Tool::Local(Z[i].cood,N,P_E,P_F);
		//Z[i].energy = Tool::Local(Z[i].cood,N,P_Enew,P_Fnew);

		//���þֲ��Ż�����
		if(Y[i].energy != Z[i].energy)
		{
			this->Localtimes++;
		}

		//if(Z[i].energy < best2->energy)
		//{
		//	best2 = this->Z+i;
		//  
		//} 
		//cout<<"��Ⱥ�е�"<<i+1<<"�������Ŷ��������"<<Z[i].energy<<endl;
	}
	cout<<endl;
	++Kate;
	cout<<"�ֲ��Ż�����"<<Localtimes<<endl<<endl;
	//cout<<"��"<<Kate<<"��SphereCut�Ŷ������Ⱥ�е����Ž��"<<best2->energy<<endl;
	//cout<<endl;
	
	//delete [] lx1;
	//delete [] ly1;
	//delete [] lz1;
	//delete [] cx1;
	//delete [] cy1;
	//delete [] cz1;
}

void BasinHopping::compare(int independent_times,int fate,double durtime)//�������ʵ��Ĵ���
{
	BH_Individual *best3,*bestGA;
	int i,vote;
	static int j=0;
	double E0;
	int GArand=0;
	//int all = 0;
	char root[200];
	const int rat=1;
	static int keep=0;
	//static int Localsuccesstimes=0;
	int better = 0,bet=0;
	best3 = this->X;
	bestGA = this->X;
	for (i=0;i<popsize;i++)
	{
		if(X[i].energy>=Z[i].energy)
			{
				
				this->X[i].Copy(&Z[i],this->N);

				this->Localsuccesstimes++;
				better++;
			} 
		
		//this->X[i].cood=this->Z[i].cood;
			
		if(X[i].energy < best3->energy)
			{
				best3 = this->X+i;
				j=i;
			}
		
		//cout<<"��Ⱥ�е�"<<i+1<<"������ȽϺ������"<<X[i].energy<<endl;
	}

	//����z�ǲ��ϱ仯�ĸ���Ⱥ��������ʵ��ʱ��ɽ�z��Ϊ����  ����ÿ�����������ϱ仯�ģ���Z�ı仯����������ģ������Ƿ����������Եĸ��� �д�ʵ����֤
	for(i=0;i<popsize,i!=j;i++)
	{
		if(Z[i].energy < bestGA->energy);
		{
			//bestGAΪ���ϱ仯��Z��Ⱥ�е����Ÿ���ȥ
			bestGA = this->Z+i;
		}
	}

	keep++;
	//if(best.energy > best3->energy)
	//{	
		/*cout<<"------------------------------------*************_______________________________________"<<endl;

		cout<<"����GA�ֲ��Ż�qian���������������"<<'\t'<<best.energy<<endl;
		geneticTool GA(N,FeN);
		GA.Genetic(best.cood,N,FeN);
		best.energy = Tool::Local(best.cood,N,P_E,P_F);
		cout<<"����GA�ֲ��Ż�����������������"<<'\t'<<best.energy<<endl;

		cout<<"------------------------------------*************_______________________________________"<<endl;

		vote=rand()%49;
		this->Z[vote].Copy(&best,N)*/;
		//this->Z[j].Copy(&best,N);

	//	}
	//else
	//this->GAbest.Copy(bestGA,N);
	//cout<<"GAbest.energy="<<GAbest.energy<<endl;
	this->best.Copy(best3,N);
	this->GAbest.Copy(best3,N);
	//cout<<"best.energy="<<best.energy<<endl;
	//best.energy = best3->energy;
	//������ʹ���Ŵ��ֲ��Ż��õ�ÿ�εıȽ�֮��õ������Ž����Ԫ�����еľֲ����ţ�*****************************&&&&&&&&&&&&&^^^^^^^^^*********

	/*cout<<"׼������GA�ֲ��Ż����������������"<<j<<'\t'<<Z[j].energy<<endl;
	geneticTool GA(N,FeN);
	GA.Genetic(Z[j].cood,N,FeN);
	Z[j].energy = Tool::Local(Z[j].cood,N,P_E,P_F);
	cout<<"����GA�ֲ��Ż�����������������"<<j<<'\t'<<Z[j].energy<<endl;*/

	//������ʹ���Ŵ��ֲ��Ż��õ�ÿ�εıȽ�֮��õ������Ž����Ԫ�����еľֲ����ţ�****************************&&&&&&&&&&&&&^^^^^^^^^*********

	cout<<"��"<<keep<<"��SphereCut�ȽϺ����Ⱥ�е����Ž��"<<best3->energy<<endl<<endl;
	cout<<"�ֲ��Ż��ɹ�������"<<Localsuccesstimes<<endl;
	
	//���ÿ�αȽϺ���������Ժ����ݷ����Ƚ�ʵ���ʱ���ã�
	//Tool::printfEnergy(keep,best3->energy,independent_times,root,N,FeN);
	Tool::printfEnergy(keep,best3->energy,independent_times,root,N,FeN,Localtimes,Localsuccesstimes,durtime);

	/*geneticTool GA(N,FeN);
	GA.Genetic(bestcood,N,FeN);
	E0 = Tool::Local(best3->cood,this->N,this->P_E,this->P_F);
	cout<<"�Ŵ��ֲ��Ż������Ž����"<<E0<<endl;*/

	//����ҵ�������������Diamond��Ҫ��ԭ�ӵ����ꣻ
	if (best3->energy <= Ebest)
	{
		S_times ++;
		//��ʱ�� ����ģ�
		//Tool::printfDiamond(best3->cood,N,best3->energy,root,this->Localtimes,Localsuccesstimes,independent_times);
		//break;
	}

	if (S_times == 1||fate < 1)
	{
		Tool::printfDiamond(best3->cood,N,best3->energy,root,this->Localtimes,Localsuccesstimes,independent_times,FeN);
		Tool::printfDiamond1(best3->cood,N,best3->energy,root);
	}

	//���������ȥ��Ч����
	//��Ҫÿ�ζ�����𣿻���ֻҪ���һ�������

	/*if (better <= rat)
	{	
		cout<<"better:"<<better<<endl;
		cout<<"SphereCut����"<<endl;
		this->SphereCut();
		this->SphereCutCompare();
		better -= better;
	}

	else
	{
		better -= better;
	}*/
	

	/*if (better = 1)
	{
	bet++;	
	better--;
	}

	if (bet>=rat)
	{   
	cout<<"SphereCut����"<<endl;
	this->SphereCut();
	bet -= rat;

	}

	cout<<endl;*/
  
}

void BasinHopping::SphereCut()
{
	int j=0,h=0,i,k=0,ren,g,m=0,a[100]={0},temp;
	double E0=0;

	srand((unsigned)time(NULL));
	BH_Individual childBH1,childBH2;
	childBH1.cood = (double *)malloc(4 * this->N * sizeof(double));
	childBH2.cood = (double *)malloc(4 * this->N * sizeof(double));
	for (i=0;i<(this->popsize)/2;i++)
	{
		/*	ren=25;
		for(g=0;g<25;g++)
		{
		a[g]=g+25;
		}
		for(g=0;g<24;g++)
		{
		k=rand()%ren;
		temp=a[k];a[k]=a[ren-1];a[ren-1]=temp;
		ren--;
		}*/
		//m=rand()%25+26;//���ֻ��
		Tool::SphereCutSplice(this->Z[i].cood,this->Z[popsize - 1 - i].cood,childBH1.cood,childBH2.cood,this->N);
		//Tool::SphereCutSplice(this->Z[i].cood,this->Z[a[i]].cood,childBH1.cood,childBH2.cood,this->N);
		//Tool::adjustment(childBH1.cood,this->N,FeN);
		//Tool::adjustment(childBH2.cood,this->N,FeN);
		childBH1.energy = Tool::Local(childBH1.cood,this->N,this->P_E,this->P_F);
		childBH2.energy = Tool::Local(childBH2.cood,this->N,this->P_E,this->P_F);
		//Genetic(child2.cood,N);
		//geneticTool GA(N,FeN);
		//GA.Genetic(childBH1.cood,N,FeN);
		//Genetic(child1.cood,N);
		//cout<<"childBH1.energy"<<childBH1.energy<<endl;
		//cout<<"childBH2.energy"<<childBH2.energy<<endl;
		if(childBH1.energy < childBH2.energy)
		{
			this->Z[i].Copy(&childBH1,this->N);
			
		} 
		else
		{
			this->Z[i].Copy(&childBH2,this->N);
			
		}
		//�õ��ȽϺõĽ�����������Ŵ��ľֲ��Ż���
		

	}




	
	/*if(j>=49)
	{
		j=49-j;
	}
*/
	childBH1.cood = NULL;
	childBH2.cood = NULL;
	free(childBH1.cood);
	free(childBH2.cood);



	cout<<"------------------------------------*************_______________________________________"<<endl;
	//h=rand()%20;
	//cout<<"����GA�ֲ��Ż�qian���������������"<<'\t'<<GAbest.energy<<endl;
	geneticTool GA(N,FeN);
	GA.Genetic(GAbest.cood,N,FeN);
	GAbest.energy = Tool::Local(GAbest.cood,N,P_E,P_F);
	//cout<<"����GA�ֲ��Ż�����������������"<<'\t'<<GAbest.energy<<endl;

	cout<<"------------------------------------*************_______________________________________"<<endl;

	h=rand()%29;
	this->Z[h].Copy(&GAbest,this->N);


}

void BasinHopping::Genetic(double *cood,int N)
{
	
	double E1 = 0;
	this->N = N;
	int POPSIZE = 50;
	int a[100]={0};
	int i,j,g,o=0,ren,temp;
	//double *R,*cood_new;
	double *coor,*cood_L;

	for (i=0;i<popsize;i++)
	{
		//A[i].cood = cood
		ren = N;
		//�������������
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
				Z[i].cood[g] = 1;
			}
			else 
			{
				Z[i].cood[g] = 0;
			}
		}

		for(int m=N;m<4*N;m++){
			Z[i].cood[m] = cood[m];//������Լ�Ӹ�ֵ�ȽϷ��㡢

		}

	}

}

void BasinHopping::SphereCutCompare()
{
	static int SphereCuttime=0;
	for (int i=0;i<popsize;i += 2)
	{

		X[i].cood=Z[i].cood;
		X[i].energy=Z[i].energy;
		SphereCuttime++;
		cout<<endl;
		cout<<"�����Ŷ��Ĵ�����"<<SphereCuttime<<endl<<endl;
	}
}


double BasinHopping::Greedy()//�˴���Ҫ����������ѡ�񡣡�������
{	
	prebest = best;
	if(best.energy>=change.energy)
		{
			best = change;
			//current = change;
			cout<<"T�ȽϺ�����ֵ��"<<best.energy<<endl;
		}
	else{
		//current = change;
		cout<<"F�ȽϺ������ֵ��"<<best.energy<<endl;
	}
	
	return best.energy;
}

double BasinHopping::changeRate()
{
	double rate = 0;
	/*kIndex ++;
	k[kIndex] = fabs((prebest.energy - bestE) / bestE); ////����û��ȡ�����10������ֵ
	cout<<"�ֲ��Ż�����ֵ��"<<best.energy<<endl;
	for(int i= 0; i < 10; i++)
		rate += k[i];
	rate /= 10;*/
	rate = fabs((prebest.energy - Ebest) / Ebest);
	cout<<"��"<<rate<<endl;
	return rate;

}

/*void check(){
	double pNum;
	cout<<"��������Ⱥ����"<<endl;
	cin>>pNum;
    this->popsize=pNum;
}

/*bool BasinHopping::EndingCondition(BasinHopping &BH){
	//BH_Individual current;
	//BH_Manage BM;
	//this->objbest=objbest; 
	//BH_Individual best;
	double m=0.001;
	//if(BH.changeRate()<m)//rateû�и�ֵ����
	if((current.energy-objbest)/objbest<=m)
	{
		cout<<"�������˽���"<<endl;
		cout<<"����ʵ�����ǣ�"<<current.energy<<endl;
		return false;
	}
	else
	{
		cout<<"��ǰ�ֲ��Ż�����������ֵ��"<<current.energy<<endl;
		return true;
	}
}*/

double BasinHopping::compare_end()//�������ʵ��Ĵ���
{
	
	return best.energy;
  
}

void BasinHopping::adjustment()
{
	int oneNum=0;
	int ren,g,o;
	int a[100];
	int temp;
	int zeroNum=0;
	cout<<"��ʼ����"<<endl;
	for (int i=0;i<popsize;i++)
	{
		ren = this->N;
		//�������������
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
				Z[i].cood[g] = 1;
			}
			else 
			{
				Z[i].cood[g] = 0;
			}
		}

	}

}

void BasinHopping::exchange()
{
	int i=0,j=0;
	double dcentre[3] = {0};
	int K=0;
	//double *cx,*;cy,*cz;
	double *distancecentre,temp,r0;
	temp = 0;r0 = 2; //��ʱ������ɵķ�ΧӦ�ñ�֮ǰ��ʼ����ҪС�ź���;
	distancecentre = (double *)malloc(N * sizeof(double));
	/*for (i=0;i<popsize;i++)
	{*/
	//	temp = Z[0].energy;
	//	if(Z[i].energy > temp)
	//	{
	//		temp = Z[i].energy;
	//		//K = i;
	//	}
	//}
	//cout<<"���Ž��������"<<temp<<endl;
	for (i=0;i<popsize;i++)
	{
		if(abs(best.energy - Z[i].energy) >= 2.5)
		//if(i == K)	
		//if ((Ebest - Z[i].energy)/Ebest >= 0.07)
		{

 
			//cout<<"��"<<i<<"����Ⱥ��Ҫ����"<<"�������ǣ�"<<Z[i].energy<<endl;
			//cx = Z[i].cood + N;
			//cy = Z[i].cood + 2*N;
			//cz = Z[i].cood + 3*N;
			//for (j=0;j<N;j++)
			//{
			//	dcentre[0] += cx[j];
			//	dcentre[1] += cy[j];
			//	dcentre[2] += cz[j];

			//	dcentre[0] /= N;
			//	dcentre[1] /= N;
			//	dcentre[2] /= N;

			//}

			//for (j=0;j<N;j++)
			//{
			//	distancecentre[j] = (cx[j]-dcentre[0])*(cx[j]-dcentre[0]) + (cy[j]-dcentre[1])*(cy[j]-dcentre[1]) + (cz[j]-dcentre[2])*(cz[j]-dcentre[2]);
			//	temp = distancecentre[0];
			//	//cout<<"����������Զ�ľ����ǣ�"<<distancecentre[j]<<endl;
			//}

			//for (j=1;j<N;j++)
			//{
			//	if (temp<distancecentre[j])
			//	{
			//		temp = distancecentre[j];

			//	}
			//	

			//}
			Tool::distool(Z[i].cood,distancecentre,N,FeN,temp,K);
			//cout<<"��Զ�ľ��룺"<<temp<<endl;
			for (j=0;j<N;j++)
			{
				if (distancecentre[j]>=temp)
				{
					
					//cout<<"��"<<j<<"��ԭ����Ҫ������"<<endl;
					Z[i].cood[j+N] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
					Z[i].cood[j+2*N] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
					Z[i].cood[j+3*N] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
					
				}


			}
			
			cout<<"�仯֮ǰ������:"<<Z[i].energy<<endl;
			Z[i].energy = Tool::Local(Z[i].cood,N,P_E,P_F);
			cout<<"�仯֮�������:"<<Z[i].energy<<endl;

		}	
	}
}