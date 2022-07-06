#include<stdio.h>
#include<math.h>
#define N 10
#define eps pow(10,-6)

/*double f(double x[],double g[],double t)
{
	double s;
	s = pow(x[0]-t * g[0],2)+4 * pow(x[1]-t * g[1],2);
	return s;
}

void st(double *a,double *b,double x[],double g[])
{
	double t0,t1,t,h,alpha,f0,f1;
	int k=0;
	t0=10;
	h=1;
	alpha=2;
	f0=f(x,g,t0);
	t1=t0+h;
	f1=f(x,g,t1);
	while(1)
	{	
		if(f1<f0)
		{
			h=alpha * h;
			t=t0;
			t0=t1;
			f0=f1;
			k++;
		}
		else
		{
			if(k==0)
			{h=-h;t=t1;}
			else
			{
				*a=t<t1?t:t1;
				*b=t<t1?t:t1;
				break;
			}
		}

		t1=t0+h;
		f1=f(x,g,t1);
	}
}

double ji(double x[],double g[])
{
	double beta,t1,t2,t;
	double f1,f2;
	double a=0,b=0;
	double *c,*d;
	c=&a;d=&b;
	st(c,d,x,g);
	printf("\n[a,b]=[%lf,%lf]",a,b);
	beta=(sqrt(5)-1.0)/2;
	t2=a+beta*(b-a);
	f2=f(x,g,t1);
	t1=a+b-t2;
	f1=f(x,g,t1);
	while(1)
	{
		if(fabs(t1-t2)<eps)
			break;
		else
		{
		if(f1<f2)
		{
			t=(t1+t2)/2;
			b=t2;
			t2=t1;
			f2=f1;
			t1=a+b-t2;
			f1=f(x,g,t1);
		}
		else
		{
			a=t1;
			t1=t2;
			f1=f2;
			t2=a+beta*(b-a);
			f2=f(x,g,t2);

		}
		}
	}
	t=(t1+t2)/2;
	return t;
}

void cool()
{
	double x[N],g[N],t=0,f0,mod;
	int i,n;
	printf("请输入n=");
	scanf("%d",&n);
	printf("\n请输入初始值：\n");
	for(i=0;i<n;i++)
		scanf("%lf",&x[i]);
	f0=f(x,g,t);
	g[0]=2*x[0];g[1]=8*x[1];
	t=ji(x,g);
	printf("t=%lf",mod);
	while(1)
	{
		mod=sqrt(pow(g[0],2)+pow(g[1],2));
		printf("mod=%lf\n",mod);
		if(mod<eps)
			break;
		else
	}

}*/