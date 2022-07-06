#include "DE_Manage.h"
#include "Common.h"
 #include <direct.h>
#include "Tool.h"
extern string EnergyName;

DE_Manage::DE_Manage(void)
{
	this->popNum = 0;
	this->pop = NULL;
	this->numList = NULL;
}

DE_Manage::~DE_Manage(void)
{
	cout<<"DE_Manage����"<<endl;
}

void DE_Manage::AddPop(DE *p)
{
	int i;
	DE *temp;

	if(this->popNum == 0) 
		this->pop = p;
	else
	{
		temp = this->pop;
		for( i = 1; i < this->popNum; i++)
			temp = temp->next;
		temp->next = p;
	}
	p->next = this->pop;
	this->popNum ++;
}

void DE_Manage::removeAllPop()
{
	int i;
	DE *p,*q;
	p = q = this->pop;
	for(i = 0; i<this->popNum; i++)
	{
		q = p->next;
		delete p;
		p = q;
	}
	this->pop = NULL;
	this->popNum = 0;
}

double DE_Manage::MultiPop_Master(int N)
{
	int i,j,S;
	DE_Individual *best;
	DE *p;
	char root[20];
	double bestEnergy;
	FILE *fp;

	p = new DE(DERAND1,N);
	this->AddPop(p);
	p = new DE(DEBEST1,N);
	this->AddPop(p);
	//p = new DE(DEBEST2,N);
	//this->AddPop(p);
	//p = new DE(DERAND2,N);
	//this->AddPop(p);
	p = new DE(DERANDTOBEST1,N);
	this->AddPop(p);
	
	N = p->N;
	
	sprintf(root,"%s_%d.txt",EnergyName.c_str(),N);
	S = 10;
	best = DE_Individual::Individual(p->N);
	
	atomNum = N;
	readBestResult();

	char energyFile[50];
	sprintf(energyFile,"%s_%d",EnergyName.c_str(),N);
	mkdir(energyFile);
	sprintf(energyFile,"%s_%d\\%d.txt",EnergyName.c_str(),N,current);
	fp = fopen(energyFile,"w+");

	for(i = 0; i < iterations; i++)
	{
		p = this->pop;				//����һ��
		for(j = 0; j < this->popNum; j++)	
		{
			p->OneGeneration();
			p = p->next;
		}

		best->energy = 0;			//�õ�ÿ����Ⱥ���Ÿ���
		p = this->pop;
		for(j = 0; j < this->popNum; j++)	
		{
			if(p->bestPop.energy < best->energy)
				best->Copy(&p->bestPop,p->N);
			p = p->next;
		}
		cout<<i<<" "<<best->energy;
		fprintf(fp,"%d\t%f",i,best->energy);

		p = this->pop;
		for(j = 0; j < this->popNum; j++)	
		{
			cout<<" "<<p->bestPop.energy;
			fprintf(fp,"\t%f",p->bestPop.energy);
			p = p->next;
		}
		cout<<endl;
		fprintf(fp,"\n");

		if(i % S != 0)
			continue;
		cout<<"��Ⱥ�������Ÿ���"<<endl;
		p = this->pop;				//�����Ÿ��巢�͸�ÿ����Ⱥ
		for(j = 0; j < this->popNum; j++)	
		{
			p->RecIndividual(best);
			p = p->next;
		}

		if (best->energy < bestResult)
		{
			Tool::printfDiamond(best->cood,N,best->energy,root);
			bestResult = best->energy;
		}
		
		bestEnergy = best->energy;
	}

	fclose(fp);
	this->removeAllPop();
	free(best->cood);
	free(best);
	
	return bestEnergy;
}

void DE_Manage::inputPara(){
	int i,num;
	char c;
	char numStr[100];

	cout<<"������ԭ�Ӹ���(��13,14,15-18):"<<endl;
	do 
	{
		gets(numStr);
		i = 0;
		c = numStr[i];
		while (c!=0)
		{
			if(!(c==','||c=='-'||(c>='0'&&c<='9')))
				break;
			c = numStr[++i];
		}
		if (c == 0)
			break;
		cout<<"���ڷǷ��ַ�,����������:"<<endl;
	} while (1);

	i = 0;
	c = numStr[i];
	num = 0;
	while (c!='\0')
	{
		if (c >= '0' && c <= '9')
			num = num * 10 + (c - '0');

		if (c==','){
			AddToTail(&numList,num);
			num = 0;
		}

		if (c=='-')
		{
			int startNum = num;
			int endNum = 0;
			while (c!=',')
			{
				c = numStr[++i];
				if (c >= '0' && c <= '9')
					endNum = endNum * 10 + (c - '0');
			}
			for (int j = startNum; j <= endNum; j++)
				AddToTail(&numList,j);
			num = 0;
		}
		c = numStr[++i];
	}
	cout<<"����������������:"<<endl;
	cin>>iterations;
	cout<<"�������ظ�������"<<endl;
	cin>>repeat;
}

void DE_Manage::start(){
	if (numList==NULL)
	{
		cout<<"û������ԭ��������"<<endl;
		return;
	}
	ListNode *p;
	p = numList;
	while (p!=NULL)
	{	
		current = 0;
		for (current=0;current<repeat;current++)
		{
			MultiPop_Master(p->m_nValue);
		}
		
		p = p->m_pNext;
	}
}

void DE_Manage::Repeat(){
	int N,success =0;
	double E,objEnergy,successRate = 0;

	cout<<"������ԭ������"<<endl;
	cin>>N;
	cout<<"������������������"<<endl;
	cin>>iterations;
	cout<<"������Ŀ��������"<<endl;
	cin>>objEnergy;
	cout<<"�������ظ�������"<<endl;
	cin>>repeat;
	//�������ɹ��ʵĺ�����
	current = 0;
	for (current = 0; current < repeat; current++ )
	{
		E = MultiPop_Master(N);
		if (fabs(E - objEnergy) < 0.01)
		{
			success ++;
		}
	}
	successRate = success*100.0 / repeat;
	cout<<"һ������"<<repeat<<"�Σ����гɹ�"<<success<<"�Σ��ɹ���Ϊ"<<successRate<<"%"<<endl;
}

void DE_Manage::readBestResult(){
	char path[100];
	char line[100];
	int atomN;
	double E;
	sprintf(path,"%s_%d.txt",EnergyName.c_str(),atomNum);

	FILE *fp = fopen(path,"r");
	if (fp != NULL)
	{
		fscanf(fp,"%d",&atomN);
		fscanf(fp,"%s  %lf",line,&E);
		fclose(fp);
		bestResult = E;
	} 
	else
	{
		bestResult = 0;
	}
}