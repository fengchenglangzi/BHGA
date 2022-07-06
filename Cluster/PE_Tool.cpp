//定义两种势能函数（ＦＳ和Johnson势能），并且求得其导数，为后面的求局部优化做准备



#include "PE_Tool.h"


PE_Tool::PE_Tool(void)
{
}

PE_Tool::~PE_Tool(void)
{
}



//计算FS势能函数；
//double PE_Tool::FS_E(double *cood,double *R,int N)
//{
//	int i,j;
//	double d,A,beta,c,c0,c1,c2;
//	double d_1,A_1,beta_1,c_1,c0_1,c1_1,c2_1;
//	double d_2,A_2,beta_2,c_2,c0_2,c1_2,c2_2;
//	double r;
//	double tempV,tempP;
//	double *coor,*VEN,*PEN;
//	double E = 0;
//	double minR = 0,minR_2 = 0,tempPmin;
//
//	//Fe
//	d = 3.569745;A = 1.828905;beta = 1.8;c = 3.40;c0 = 1.2371147;c1 = -0.3592185; c2 = -0.0385607;
//	//d = 3.699579;A = 1.889846;beta = 1.8;c = 3.40;c0 = 1.2110601;c1 = -0.7510840; c2 = 0.1380773;
//	//Mo
//	d_1 = 4.114825;A_1 = 1.887117;beta_1 = 0;c_1 = 3.25;c0_1 = 43.4475218;c1_1 = -31.9332978; c2_1 = 6.084249;
//	//Fe-Mo
//	d_2 = 4.114825;A_2 = 1.887117;beta_2 = 0.9;c_2 = 3.25;c0_2 = 43.4475218;c1_2 = -31.9332978; c2_2 = 6.084249;
//	minR = 1.8;  minR_2 = d_2*(beta_2-1)/beta_2;//
//	coor = cood;
//	VEN = (double *)calloc(N,sizeof(double));
//	PEN = (double *)calloc(N,sizeof(double));
//	//tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
//	
////	tempVmin = (minR-c)*(minR-c)*(c0+c1*minR+c2*minR*minR);
//	
//	for(i = 0; i < N-1; i++)
//		for(j = i+1; j < N; j++)
//		{
//			r = *(R + i*N + j);
//			if (coor[i]==1&&coor[j]==1)
//			{
//				tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
//				tempV = (r>=c)?0 : (r-c)*(r-c)*(c0+c1*r+c2*r*r);
//				tempP = (r>=d)?0 : ((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);
//			}
//			else if (coor[i]==0&&coor[j]==0)
//			{
//				tempV = (r>=c_1)?0 : (r-c_1)*(r-c_1)*(c0_1+c1_1*r+c2_1*r*r);
//				tempP = (r>=d_1)?0 : (r-d_1)*(r-d_1);
//			}
//			else
//			{
//				tempPmin = (minR_2-d_2)*(minR_2-d_2)+beta_2*(minR_2-d_2)*(minR_2-d_2)*(minR_2-d_2)/d_2;
//				tempV = (r>=c_2)?0 : (r-c_2)*(r-c_2)*(c0_2+c1_2*r+c2_2*r*r);
//				tempP = (r>=d_2)?0 : ((r<=minR_2)?tempPmin:(r-d_2)*(r-d_2)+beta_2*(r-d_2)*(r-d_2)*(r-d_2)/d_2);
//			}
//
//			//tempV = (r>=c)?0 : (r-c)*(r-c)*(c0+c1*r+c2*r*r);
//			//tempP = (r>=d)?0 : ((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);
//
//			VEN[i] += tempV/2;
//			VEN[j] += tempV/2;
//			PEN[i] += tempP;
//			PEN[j] += tempP;
//		}
//
//	for(i=0 ;i < N; i++){
//		PEN[i] = sqrt(PEN[i]);
//		E += VEN[i] - A * PEN[i];
//	//	cout<<VEN[i]<<" "<<A * PEN[i]<<" "<<VEN[i]-A*PEN[i]<<endl;
//	}
////	Tool::printCood(PEN,N,"PEN.txt");
//	free(VEN);
//	free(PEN);
//	return E;
//}
//
////计算FS势能函数对应的梯度；
//double PE_Tool::FS_F(double *cood,double *R,double *FF, int N)
//{
//	int i,j;	
//	double d,A,beta,c,c0,c1,c2;
//	double d_1,A_1,beta_1,c_1,c0_1,c1_1,c2_1;
//	double d_2,A_2,beta_2,c_2,c0_2,c1_2,c2_2;
//	double r;
//	double tempdV,tempP,tempdP;
//	double *PEN;
//	double maxF = 0,FK;
//	double minR = 0,minR_2 = 0,tempPmin;
//	double *coor,*x,*y,*z,*FX,*FY,*FZ;
//
//	PEN = (double *)calloc(N,sizeof(double));
//
//	//Fe
//	d = 3.569745;A = 1.828905;beta = 1.8;c = 3.40;c0 = 1.2371147;c1 = -0.3592185; c2 = -0.0385607;
//	//d = 3.699579;A = 1.889846;beta = 1.8;c = 3.40;c0 = 1.2110601;c1 = -0.7510840; c2 = 0.1380773;
//	d_1 = 4.114825;A_1 = 1.887117;beta_1 = 0;c_1 = 3.25;c0_1 = 43.4475218;c1_1 = -31.9332978; c2_1 = 6.084249;
//	d_2 = 4.114825;A_2 = 1.887117;beta_2 = 0.9;c_2 = 3.25;c0_2 = 43.4475218;c1_2 = -31.9332978; c2_2 = 6.084249;
//	coor = cood;x = cood + N; y = cood + 2 * N; z = cood + 3 * N; FX = FF; FY = FF + N; FZ = FF + 2 * N;
//	minR = 1.8;  minR_2 = d_2*(beta_2-1)/beta_2;//此时得到的beta是不能为0；若beta_2为0则minR_2不要;
//
//	for(i = 0; i < 3 * N; i++)
//		FF[i] = 0;
//
//	
//	//tempdPmin = 2 * (minR -d) + 3 * beta * (minR - d) * (minR - d) / d;
//	//tempdVmin = (2*(minR - c)*(c0 + c1 * minR + c2 * minR * minR) + (c1 + 2 * c2 * minR)* (minR - c) * (minR - c)) * 2;
//	for(i=0;i<N-1; i++)
//		for(j=i+1;j<N;j++)
//		{
//			r = *(R + i*N + j);
//			if (coor[i]==1&&coor[j]==1)
//			{
//				tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
//				tempP = (r>=d)?0:((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);
//			}
//			else if (coor[i]==0&&coor[j]==0)
//			{
//				
//				tempP = (r>=d_1)?0:(r-d_1)*(r-d_1);
//			}
//			else
//			{
//				tempPmin = (minR-d_2)*(minR-d_2)+beta_2*(minR-d_2)*(minR-d_2)*(minR-d_2)/d_2;
//				tempP = (r>=d_2)?0:((r<=minR_2)?tempPmin:(r-d_2)*(r-d_2)+beta_2*(r-d_2)*(r-d_2)*(r-d_2)/d_2);
//			}
//			
//
//			PEN[i] += tempP;
//			PEN[j] += tempP;
//		}
//
//		for(i = 0; i < N ;i ++)
//			PEN[i] = (PEN[i] == 0)?0:1 / sqrt(PEN[i]) / 2;
//
//		for(i = 0; i < N-1;i ++)
//			for(j=i+1;j<N;j++)
//			{
//				r = *(R + i*N + j);
//				if (coor[i]==1&&coor[j]==1)
//				{
//					tempdV = (r>=c)? 0 : (2 * (r - c) * (c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r) * (r - c) * (r - c));
//					tempdP = (r>=d)? 0 : (2 * (r -d) + 3 * beta * (r - d) * (r - d) / d);
//				}
//				else if (coor[i]==0&&coor[j]==0)
//				{
//					tempdV = (r>=c_1)? 0 : (2 * (r - c_1) * (c0_1 + c1_1 * r + c2_1 * r * r) + (c1_1 + 2 * c2_1 * r) * (r - c_1) * (r - c_1));
//					tempdP = (r>=d_1)? 0 : (2 * (r-d_1));
//				}
//				else
//				{
//					tempdV = (r>=c_2)? 0 : (2*(r - c_2)*(c0_2 + c1_2 * r + c2_2 * r * r) + (c1_2 + 2 * c2_2 * r) * (r - c_2) * (r - c_2));
//					tempdP = (r>=d_2)? 0 : (2 * (r -d_2) + 3 * beta_2 * (r - d_2) * (r - d_2) / d_2);
//				}
//				//	tempdV = (r>=c)?0: (r<=minR)?tempdVmin : (2*(r - c)*(c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r)* (r - c) * (r - c)) * 2;
//				//	tempdP = (r>=d)?0 : (r<=minR)?tempdPmin : 2 * (r -d) + 3 * beta * (r - d) * (r - d) / d;
//				//tempdV = (r>=c)?0 : (2*(r - c)*(c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r)* (r - c) * (r - c)) * 2;
//				//tempdP = (r>=d)? 0 : (2 * (r -d) + 3 * beta * (r - d) * (r - d) / d);
//
//				tempdP = (PEN[i] + PEN[j]) * tempdP;
//
//				FK = -tempdV / 2 + A * tempdP / 2;
//
//				FX[i]=FX[i]+FK*(x[i]-x[j])/r;
//				//	cout<<j<<" "<<FX[i]<<" "<<FK<<" "<<r<<endl;
//				FX[j]=FX[j]-FK*(x[i]-x[j])/r;
//				FY[i]=FY[i]+FK*(y[i]-y[j])/r;
//				FY[j]=FY[j]-FK*(y[i]-y[j])/r;
//				FZ[i]=FZ[i]+FK*(z[i]-z[j])/r;
//				FZ[j]=FZ[j]-FK*(z[i]-z[j])/r;
//			}
//			for(i = 0;i < 3 * N; i++)
//				maxF = (FF[i]>maxF)?FF[i]:maxF;
//
//			free(PEN);
//			return maxF;
//}

//描述Gupta势能函数
double PE_Tool::FS_E(double *cood,double *R,int N){
	int i,j;
	//double e=2.7183;
	double tempV,tempP;
	double r;
	double E=0;
	double *coor,*VEN,*PEN;
	double r0,A,xi,P,q;
	double r0_1,A_1,xi_1,P_1,q_1;
	double r0_2,A_2,xi_2,P_2,q_2;

	VEN = (double *)calloc(N,sizeof(double));
	PEN = (double *)calloc(N,sizeof(double));

	//Ag
	//r0=2.8921;A=0.1031;xi=1.1895;P=-10.85;q=-3.18;
	//Cu
	//r0=2.5560;A=0.0855;xi=1.2240;P=-10.960;q=-2.2780;
	//pt
	//r0=2.7747;A=0.2795;xi=2.6950;P=-10.6210;q=-4.0040;
	//Au
	//r0_1=2.8840;A_1=0.2061;xi_1=1.7900;P_1=-10.2290;q_1=-4.0360;
	//Cu-Au
	//r0_2=2.5560;A_2=0.1539;xi_2=1.5605;P_2=-11.0500;q_2=-3.0475;
	//Au-Pt
	//r0_2=2.8294;A_2=0.2500;xi_2=2.2000;P_2=-10.4200;q_2=-4.0200;

	//Pt
	// r0=2.7746;A=0.2975;xi=2.695;P=-10.612;q=-4.004;
	//Co
	// r0_1=2.502;A_1=0.95;xi_1=1.488;P_1=-11.604;q_1=-2.286;
	//Pt _Co
	// r0_2 = 2.63;A_2 = 0.245; xi_2 = 2.386; P_2 = -9.97; q_2 = -3.32;
	
	/*
		
		----------------------Gupta势函数相关参数设定-------------------------
		-------------------------------------------------------------------
		-------------------------------------------------------------------
		-------------------------------------------------------------------
		-------------------------------------------------------------------

	
	*/

	//Pd
	r0 = 2.7486; A = 0.1746; xi = 1.718; P = -10.867; q = -3.742;
	//Ir
	r0_1 = 2.7146; A_1 = 0.1156; xi_1 = 2.289; P_1 = -16.980; q_1 = -2.7146;
	//Pd _Ir
	r0_2 = 2.7316; A_2 = 0.1421; xi_2 = 1.983; P_2 = -13.923; q_2 = -3.216;

	//Pt
	//r0=2.7746;A=0.2975;xi=2.695;P=-10.612;q=-4.004;
	////Fe
	//r0_1=2.553;A_1=0.13315;xi_1=1.6179;P_1=-10.50;q_1=-2.60;
	////Pt-Fe
	//r0_2=2.6638;A_2=0.19903;xi_2=2.0881;P_2=-10.556;q_2=-3.302;

		//Co
	//r0=2.502;A=0.95;xi=1.488;P=-11.604;q=-2.286;
	////Cu
//	r0_1=2.5560;A_1=0.0855;xi_1=1.224;P_1=-10.96;q_1=-2.278;
	////Co-Cu
///	r0_2=2.54;A_2=0.9;xi_2=1.3301;P_2=-11.282;q_2=-2.2820;
	coor = cood;

	for(i = 0; i < N-1; i++)
			for(j = i+1; j < N; j++)
				{
					r = *(R + i*N + j);
					if (coor[i]==1&&coor[j]==1)
					{
						tempV = A_1*exp(P_1*(r/r0_1 - 1));
						tempP = xi_1*xi_1*exp(2*q_1*(r/r0_1 - 1));
						
					}
					else if (coor[i]==0&&coor[j]==0)
					{
						tempV = A*exp(P*(r/r0 - 1));
						tempP = xi*xi*exp(2*q*(r/r0 - 1));
					}
					else
					{
						tempV = A_2*exp(P_2*(r/r0_2 - 1));
						tempP = xi_2*xi_2*exp(2*q_2*(r/r0_2 - 1));
					}
					
					VEN[i] += tempV;
					VEN[j] += tempV;
					PEN[i] += tempP;
					PEN[j] += tempP;
					
				}
				
				for(i=0 ;i < N; i++){
					PEN[i] = sqrt(PEN[i]);
					E += VEN[i] - PEN[i];
				}
				
				free(VEN);
				free(PEN);
				return E;
				


}

double PE_Tool::FS_F(double *cood,double *R,double *FF,int N){
	int i,j;
	double e=2.7183;
	double tempdV,tempP,tempdP;
	double r;
	double E=0,FK,Fmax = 0;
	double *coor,*VEN,*PEN;
	double *x,*y,*z,*FX,*FY,*FZ;
	double r0,A,xi,P,q;
	double r0_1,A_1,xi_1,P_1,q_1;
	double r0_2,A_2,xi_2,P_2,q_2;

	//Ag
	//r0=2.8921;A=0.1031;xi=1.1895;P=-10.85;q=-3.18;

	//Cu
	//r0=2.5560;A=0.0855;xi=1.2240;P=-10.960;q=-2.2780;
	//pt
	//r0=2.7747;A=0.2795;xi=2.6950;P=-10.6210;q=-4.0040;
	//Au
	//r0_1=2.8840;A_1=0.2061;xi_1=1.7900;P_1=-10.2290;q_1=-4.0360;
	//Cu-Au
	//r0_2=2.5560;A_2=0.1539;xi_2=1.5605;P_2=-11.0500;q_2=-3.0475;

	//Au-Pt
	//r0_2=2.8294;A_2=0.2500;xi_2=2.2000;P_2=-10.4200;q_2=-4.0200;

	/*
	*  -----------------Pt-Co合金晶格参数--------------------------
	//Pt
	r0=2.7746;A=0.2975;xi=2.695;P=-10.612;q=-4.004;
	//Co
	r0_1=2.502;A_1=0.95;xi_1=1.488;P_1=-11.604;q_1=-2.286;
	//Pt _Co
	r0_2 = 2.63;A_2 = 0.245; xi_2 = 2.386; P_2 = -9.97; q_2 = -3.32;
	*/
	// ----------------Pd-Ir合金晶格参数----------------------------
	//Pd
	r0 = 2.7486; A = 0.1746; xi = 1.718; P = -10.867; q = -3.742;
	//Ir
	r0_1 = 2.7146; A_1 = 0.1156; xi_1 = 2.289; P_1 = -16.980; q_1 = -2.7146;
	//Pd _Ir
	r0_2 = 2.7316; A_2 = 0.1421; xi_2 = 1.983; P_2 = -13.923; q_2 = -3.216;
	

	//VEN = (double *)calloc(N,sizeof(double));
	PEN = (double *)calloc(N,sizeof(double));
	//Pt
//	r0=2.7746;A=0.2975;xi=2.695;P=-10.612;q=-4.004;
	////Fe
//	r0_1=2.553;A_1=0.13315;xi_1=1.6179;P_1=-10.50;q_1=-2.60;
	////Pt-Fe
//	r0_2=2.6638;A_2=0.19903;xi_2=2.0881;P_2=-10.556;q_2=-3.302;

	//Co
	//r0=2.502;A=0.95;xi=1.488;P=-11.604;q=-2.286;
	////Cu
	//r0_1=2.5560;A_1=0.0855;xi_1=1.224;P_1=-10.96;q_1=-2.278;
	////Co-Cu
	//r0_2=2.54;A_2=0.9;xi_2=1.3301;P_2=-11.282;q_2=-2.2820;

	coor = cood;x = cood + N;y = cood + 2*N;z = cood + 3*N; FX = FF; FY = FF + N; FZ = FF + 2 * N;

	for(i = 0; i < 3 * N; i++)
		FF[i] = 0;

	for(i = 0; i < N-1; i++)
		for(j = i+1; j < N; j++)
		{
			r = *(R + i*N + j);
			if (coor[i]==1&&coor[j]==1)
			{
				tempP = xi_1*xi_1*exp(2*q_1*(r/r0_1 - 1));
			}
			else if (coor[i]==0&&coor[j]==0)
			{
				tempP = xi*xi*exp(2*q*(r/r0 - 1));
			}
			else
			{
				tempP = xi_2*xi_2*exp(2*q_2*(r/r0_2 - 1));
			}
			PEN[i] += tempP;
			PEN[j] += tempP;
		}

		for(i = 0; i < N ;i ++)
			PEN[i] = (PEN[i] == 0)?0:1 / sqrt(PEN[i]);

	for(i = 0; i < N-1; i++)
		for(j = i+1; j < N; j++)
		{
			r = *(R + i*N + j);
			if (coor[i]==1&&coor[j]==1)
			{
				tempdV = A_1*P_1*exp(P_1*(r/r0_1 - 1))/r0_1;
				tempdP = xi_1*xi_1*q_1*exp(2*q_1*(r/r0_1 - 1))/r0_1;

			}
			else if (coor[i]==0&&coor[j]==0)
			{
				tempdV = A*P*exp(P*(r/r0 - 1))/r0;
				tempdP = xi*xi*q*exp(2*q*(r/r0 - 1))/r0;
			}
			else
			{
				tempdV = A_2*P_2*exp(P_2*(r/r0_2 - 1))/r0_2;
				tempdP = xi_2*xi_2*q_2*exp(2*q_2*(r/r0_2 - 1))/r0_2;
			}

			tempdP = (PEN[i] + PEN[j]) * tempdP;

			FK = -tempdV + tempdP / 2;
			
			FX[i]=FX[i]+FK*(x[i]-x[j])/r;
			//	cout<<j<<" "<<FX[i]<<" "<<FK<<" "<<r<<endl;
			FX[j]=FX[j]-FK*(x[i]-x[j])/r;
			FY[i]=FY[i]+FK*(y[i]-y[j])/r;
			FY[j]=FY[j]-FK*(y[i]-y[j])/r;
			FZ[i]=FZ[i]+FK*(z[i]-z[j])/r;
			FZ[j]=FZ[j]-FK*(z[i]-z[j])/r;
		}
		for(i = 0;i < 3 * N; i++)
		Fmax = (FF[i]>Fmax)?FF[i]:Fmax;
		free(PEN);
		return Fmax;

}


//#####%%%%%*******计算新的FS势能，这个星期的重点。
//double PE_Tool::FS_Enew(double *cood,double *R,int N)
//{
//	int i,j;
//	double d,A,beta,c,c0,c1,c2,c3,c4;
//	double d_1,A_1,beta_1,c_1,c0_1,c1_1,c2_1,c3_1,c4_1;
//	double d_2,A_2,beta_2,c_2,c0_2,c1_2,c2_2,c3_2,c4_2;
//	double r;
//	double tempV,tempP;
//	double *coor,*VEN,*PEN;
//	double E = 0;
//	double minR = 0,minR_2 = 0,tempPmin;
//
//	//Ag
//	d = 4.41;A = 0.325514;beta = -1.293394;c = 4.76;c0 = 10.6812;c1 = -12.04517; c2 = 5.203072;c3 = -1.013304;c4 = 0.0742308;
//	
//	//Mo
//	d_1 = 4.1472;A_1 = 1.848648;beta_1 = 0;c_1 = 3.2572;c0_1 = 47.98066;c1_1 = -34.09924; c2_1 = 5.832293;c3_1 = 0.017494;c4_1 = 0.020393;
//	
//	//Ag-Mo
//	d_2 = 4.3;A_2 = 0.5452;beta_2 = 0;c_2 = 4.4;c0_2 = 29.0724;c1_2 = -29.2625; c2_2 = 9.8921;c3_2 = -1.1251;c4_2 = 0;
//	minR = 1.8; 
//	coor = cood;
//	VEN = (double *)calloc(N,sizeof(double));
//	PEN = (double *)calloc(N,sizeof(double));
//	//tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
//
//	//	tempVmin = (minR-c)*(minR-c)*(c0+c1*minR+c2*minR*minR);
//
//	for(i = 0; i < N-1; i++)
//		for(j = i+1; j < N; j++)
//		{
//			r = *(R + i*N + j);
//			if (coor[i]==1&&coor[j]==1)
//			{
//				//tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
//				tempV = (r>=c)?0 : (r-c)*(r-c)*(c0+c1*r+c2*r*r+c3*r*r*r+c4*r*r*r*r);
//				tempP = (r>=d)?0 : ((r-d)*(r-d)+beta*beta*(r-d)*(r-d)*(r-d)*(r-d));
//			}
//			else if (coor[i]==0&&coor[j]==0)
//			{
//				tempV = (r>=c_1)?0 : (r-c_1)*(r-c_1)*(c0_1+c1_1*r+c2_1*r*r+c3_1*r*r*r+c4_1*r*r*r*r);
//				tempP = (r>=d_1)?0 : (r-d_1)*(r-d_1);
//			}
//			else
//			{
//				//tempPmin = (minR_2-d_2)*(minR_2-d_2)+beta_2*(minR_2-d_2)*(minR_2-d_2)*(minR_2-d_2)/d_2;
//				tempV = (r>=c_2)?0 : (r-c_2)*(r-c_2)*(c0_2+c1_2*r+c2_2*r*r+c3_2*r*r*r+c4_2*r*r*r*r);
//				tempP = (r>=d_2)?0 : (r-d_2)*(r-d_2);
//			}
//
//			//tempV = (r>=c)?0 : (r-c)*(r-c)*(c0+c1*r+c2*r*r);
//			//tempP = (r>=d)?0 : ((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);
//
//			VEN[i] += tempV/2;
//			VEN[j] += tempV/2;
//			PEN[i] += tempP;
//			PEN[j] += tempP;
//		}
//
//		for(i=0 ;i < N; i++){
//			PEN[i] = sqrt(PEN[i]);
//			E += VEN[i] - A * PEN[i];
//			//	cout<<VEN[i]<<" "<<A * PEN[i]<<" "<<VEN[i]-A*PEN[i]<<endl;
//		}
//		//	Tool::printCood(PEN,N,"PEN.txt");
//		free(VEN);
//		free(PEN);
//		return E;
//}

//计算新的FSnew势能函数对应的梯度；
//double PE_Tool::FS_Fnew(double *cood,double *R,double *FF, int N)
//{
//	int i,j;	
//	double d,A,beta,c,c0,c1,c2,c3,c4;
//	double d_1,A_1,beta_1,c_1,c0_1,c1_1,c2_1,c3_1,c4_1;
//	double d_2,A_2,beta_2,c_2,c0_2,c1_2,c2_2,c3_2,c4_2;
//	double r;
//	double tempdV,tempP,tempdP;
//	double *PEN;
//	double maxF = 0,FK;
//	double minR = 0,minR_2 = 0,tempPmin;
//	double *coor,*x,*y,*z,*FX,*FY,*FZ;
//
//	PEN = (double *)calloc(N,sizeof(double));
//
//	//Ag
//	d = 4.41;A = 0.325514;beta = -1.293394;c = 4.76;c0 = 10.6812;c1 = -12.04517; c2 = 5.203072;c3 = -1.013304;c4 = 0.0742308;
//
//	//Mo
//	d_1 = 4.1472;A_1 = 1.848648;beta_1 = 0;c_1 = 3.2572;c0_1 = 47.98066;c1_1 = -34.09924; c2_1 = 5.832293;c3_1 = 0.017494;c4_1 = 0.020393;
//
//	//Ag-Mo
//	d_2 = 4.3;A_2 = 0.5452;beta_2 = 0;c_2 = 4.4;c0_2 = 29.0724;c1_2 = -29.2625; c2_2 = 9.8921;c3_2 = -1.1251;c4_2 = 0;
//
//
//	////Fe
//	//d = 3.569745;A = 1.828905;beta = 1.8;c = 3.40;c0 = 1.2371147;c1 = -0.3592185; c2 = -0.0385607;
//	////d = 3.699579;A = 1.889846;beta = 1.8;c = 3.40;c0 = 1.2110601;c1 = -0.7510840; c2 = -0.1380773;
//	//d_1 = 4.114825;A_1 = 1.887117;beta_1 = 0;c_1 = 3.25;c0_1 = 43.4475218;c1_1 = -31.9332978; c2_1 = 6.084249;
//	//d_2 = 4.114825;A_2 = 1.887117;beta_2 = 0.9;c_2 = 3.25;c0_2 = 43.4475218;c1_2 = -31.9332978; c2_2 = 6.084249;
//	coor = cood;x = cood + N; y = cood + 2 * N; z = cood + 3 * N; FX = FF; FY = FF + N; FZ = FF + 2 * N;
//	minR = 1.8;  //minR_2 = d_2*(beta_2-1)/beta_2;//此时得到的beta是不能为0；若beta_2为0则minR_2不要;
//
//	for(i = 0; i < 3 * N; i++)
//		FF[i] = 0;
//
//
//	//tempdPmin = 2 * (minR -d) + 3 * beta * (minR - d) * (minR - d) / d;
//	//tempdVmin = (2*(minR - c)*(c0 + c1 * minR + c2 * minR * minR) + (c1 + 2 * c2 * minR)* (minR - c) * (minR - c)) * 2;
//	for(i=0;i<N-1; i++)
//		for(j=i+1;j<N;j++)
//		{
//			r = *(R + i*N + j);
//			if (coor[i]==1&&coor[j]==1)
//			{
//				//tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
//				tempP = (r>=d)?0:((r-d)*(r-d)+beta*beta*(r-d)*(r-d)*(r-d)*(r-d));
//			}
//			else if (coor[i]==0&&coor[j]==0)
//			{
//
//				tempP = (r>=d_1)?0:(r-d_1)*(r-d_1);
//			}
//			else
//			{
//				tempP = (r>=d_2)?0:(r-d_2)*(r-d_2);
//			}
//
//
//			PEN[i] += tempP;
//			PEN[j] += tempP;
//		}
//
//		for(i = 0; i < N ;i ++)
//			PEN[i] = (PEN[i] == 0)?0:1 / sqrt(PEN[i]) / 2;
//
//		for(i = 0; i < N-1;i ++)
//			for(j=i+1;j<N;j++)
//			{
//				r = *(R + i*N + j);
//				if (coor[i]==1&&coor[j]==1)
//				{
//					tempdV = (r>=c)? 0 : (2*(r - c)*(c0 + c1*r + c2*r*r + c3*r*r*r + c4*r*r*r*r) + (c1 + 2*c2*r + 3*c3*r*r + 4*c4*r*r*r)*(r - c) * (r - c));
//					tempdP = (r>=d)? 0 : (2 * (r -d) + 4 * beta * beta * (r - d) * (r - d) * (r-d));
//				}
//				else if (coor[i]==0&&coor[j]==0)
//				{
//					tempdV = (r>=c_1)? 0 : (2*(r - c_1)*(c0_1 + c1_1*r + c2_1*r*r + c3_1*r*r*r + c4_1*r*r*r*r) + (c1_1 + 2 * c2_1*r + 3*c3_1*r*r + 4*c4_1*r*r*r) * (r - c_1) * (r - c_1));
//					tempdP = (r>=d_1)? 0 : 2 * (r-d_1);
//				}
//				else
//				{
//					tempdV = (r>=c_2)? 0 : (2*(r - c_2)*(c0_2 + c1_2*r + c2_2*r*r + c3_2*r*r*r + c4_2*r*r*r*r) + (c1_2 + 2 * c2_2*r + 3*c3_2*r*r + 4*c4_2*r*r*r) * (r - c_2) * (r - c_2));
//					tempdP = (r>=d_2)? 0 : 2 * (r -d_2);
//				}
//				//	tempdV = (r>=c)?0: (r<=minR)?tempdVmin : (2*(r - c)*(c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r)* (r - c) * (r - c)) * 2;
//				//	tempdP = (r>=d)?0 : (r<=minR)?tempdPmin : 2 * (r -d) + 3 * beta * (r - d) * (r - d) / d;
//				//tempdV = (r>=c)?0 : (2*(r - c)*(c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r)* (r - c) * (r - c)) * 2;
//				//tempdP = (r>=d)? 0 : (2 * (r -d) + 3 * beta * (r - d) * (r - d) / d);
//
//				tempdP = (PEN[i] + PEN[j]) * tempdP;
//
//				FK = -tempdV / 2 + A * tempdP;
//
//				FX[i]=FX[i]+FK*(x[i]-x[j])/r;
//				//	cout<<j<<" "<<FX[i]<<" "<<FK<<" "<<r<<endl;
//				FX[j]=FX[j]-FK*(x[i]-x[j])/r;
//				FY[i]=FY[i]+FK*(y[i]-y[j])/r;
//				FY[j]=FY[j]-FK*(y[i]-y[j])/r;
//				FZ[i]=FZ[i]+FK*(z[i]-z[j])/r;
//				FZ[j]=FZ[j]-FK*(z[i]-z[j])/r;
//			}
//			for(i = 0;i < 3 * N; i++)
//				maxF = (FF[i]>maxF)?FF[i]:maxF;
//
//			free(PEN);
//			return maxF;
//}
//
//
//
//double PE_Tool::FS_Fr(double *cood,double *R,double *FF, int N)
//{
//	int i,j;	
//	double d,A,beta,c,c0,c1,c2;
//	double r;
//	double tempdV,tempP,tempdP;
//	double *PEN;
//	double maxF = 0,FK;
//	double minR = 1.6,tempPmin;
//
//	PEN = (double *)calloc(N,sizeof(double));
//
//	//Fe
//	//d = 3.569745;A = 1.828905;beta = 1.8;c = 3.40;c0 = 1.2371147;c1 = -0.3592185; c2 = -0.0385607;
//	//Mo
//	d = 4.114825;A = 1.887117;beta = 0;c = 3.25;c0 = 43.4475218;c1 = -31.9332978; c2 = 6.084249;
//	for(i = 0; i < N * N; i++)
//		FF[i] = 0;
//
//	tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
//
//	for(i=0;i<N-1; i++)
//		for(j=i+1;j<N;j++)
//		{
//			r = *(R + i*N + j);
//			tempP = (r>=d)?0:((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);
//
//			PEN[i] += tempP;
//			PEN[j] += tempP;
//		}
//
//	for(i = 0; i < N ;i ++)
//		PEN[i] = (PEN[i] == 0)?0:1 / sqrt(PEN[i]) / 2;
//
//	for(i = 0; i < N-1;i ++)
//		for(j=i+1;j<N;j++)
//		{
//			r = *(R + i*N + j);
//
//			tempdV = (r>=c)?0 : (2*(r - c)*(c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r)* (r - c) * (r - c)) * 2;
//			tempdP = (r>=d)? 0 : (2 * (r -d) + 3 * beta * (r - d) * (r - d) / d);
//
//			tempdP = (PEN[i] + PEN[j]) * tempdP;
//
//			FK = -tempdV / 2 + A * tempdP;
//			*(FF + i*N + j) = FK;
//			maxF = (abs(FK)>maxF)?abs(FK):maxF;
//		}
//
//		free(PEN);
//		return maxF;
//}

double PE_Tool::Johnson_E(double *cood,double *R, int N)
{
	double rc = 4.5;

	double re = 2.481987;
	double fe = 1.885957;
	double roue = 20.041463;
	double alpha = 9.818270;
	double beta = 5.236411;
	double A = 0.392811; 
	double B = 0.646243;
	double kappa = 0.170306;
	double lambda = 0.340613;
	double Fn0 = -2.534992;
	double Fn1 = -0.059605;
	double Fn2 = 0.193065;
	double Fn3 = -2.282322;
	double F0 = -2.54;
	double F1 = 0;
	double F2 = 0.200269;
	double F3 = -0.148770;
	double eta = 0.391750;
	double Fe = -2.539945;

	double roun = 0.85 * roue;
	double rouo = 1.15 * roue;

	int i,j;
	double E,tempRou;
	double tempR,Phi1,Phi2;
	double *Phi,*rou;

	Phi = (double *)calloc(N,sizeof(double));
	rou = (double *)calloc(N,sizeof(double));

	for(i = 0; i < N-1; i++)
	{
		for(j = i+1; j < N; j++)
		{
			tempR = *(R + i*N + j);
			if(tempR > rc)
			{
				continue;
			} else {
			tempR= tempR / re;
			Phi1 = exp(-alpha * (tempR - 1)) / (1 + pow(tempR - kappa,20));
			Phi2 = exp(-beta * (tempR - 1)) / (1 + pow(tempR - lambda,20));
			Phi[i] += (A * Phi1 - B * Phi2) / 2;
			Phi[j] += (A * Phi1 - B * Phi2) / 2;
			rou[i] += fe * Phi2;
			rou[j] += fe * Phi2;
			}
		}
	}
	
	for(i = 0; i < N; i++)
	{
		tempRou = rou[i];
		if(tempRou < roun)
		{
			tempRou = tempRou / roun - 1;
			tempRou = Fn0 + Fn1 * tempRou + Fn2 * tempRou * tempRou + Fn3 * tempRou * tempRou * tempRou;	
		} else if(tempRou < rouo)
		{
			tempRou = tempRou / roue - 1;
			tempRou = F0 + F1 * tempRou + F2 * tempRou * tempRou + F3 * tempRou * tempRou * tempRou;
		} else
		{
			tempRou = pow(tempRou/roue,eta);
			tempRou = Fe * (1-log(tempRou)) * tempRou; 
		}
		rou[i] = tempRou;
	}

	E = 0;
	for(i = 0; i < N; i ++)
	{
		E += Phi[i] + rou[i];
	}

	free(Phi);
	free(rou);
	return E;
}

double PE_Tool::Johnson_F(double *cood, double *R, double *FF, int N)
{
	double rc = 4.5;

	double re = 2.481987;
	double fe = 1.885957;
	double roue = 20.041463;
	double alpha = 9.818270;
	double beta = 5.236411;
	double A = 0.392811; 
	double B = 0.646243;
	double kappa = 0.170306;
	double lambda = 0.340613;
	double Fn0 = -2.534992;
	double Fn1 = -0.059605;
	double Fn2 = 0.193065;
	double Fn3 = -2.282322;
	double F0 = -2.54;
	double F1 = 0;
	double F2 = 0.200269;
	double F3 = -0.148770;
	double eta = 0.391750;
	double Fe = -2.539945;

	double roun = 0.85 * roue;
	double rouo = 1.15 * roue;

	int i,j;
	double r,tempR,Phi2,dPhi1,dPhi2,dPhi,drou,tempRou,FK,maxF;
	double *FX,*FY,*FZ,*x,*y,*z;
	double *rou;
	
	rou = (double*)calloc(N,sizeof(double));
	x = cood; y = cood + N; z = cood + 2 * N; FX = FF; FY = FF + N; FZ = FF + 2 * N;
	for(i = 0; i < 3 * N; i++)
		FF[i] = 0;
	
	for(i = 0; i < N-1; i++)
	{
		for(j = i+1; j < N; j++)
		{
			tempR = *(R + i*N + j);
			tempR = tempR / re;
			Phi2 = exp(-beta * (tempR - 1)) / (1 + pow(tempR - lambda,20));
			
			rou[i] += fe * Phi2;
			rou[j] += fe * Phi2;
		}	
	}

	for(i = 0; i < N; i++)
	{
		tempRou = rou[i];
		if(tempRou < roun)
		{
			tempRou = tempRou / roun - 1;
			tempRou = (Fn1 + 2*Fn2 * tempRou + 3*Fn3 * tempRou * tempRou) / roun;
		} else if(tempRou < rouo)
		{
			tempRou = tempRou / roue - 1;
			tempRou = (F1 + 2*F2 * tempRou + 3*F2 * tempRou * tempRou) / roue;
		} else
		{
			tempRou = tempRou / roue;
			tempRou = -( Fe * eta * log(pow(tempRou,eta)) * pow(tempRou,(eta-1)) ) / roue; 
		}
		rou[i] = tempRou;
	}

	for( i = 0; i < N-1; i++)
	{
		for( j = i+1; j < N; j++)
		{
			r = *(R + i*N +j);
			tempR = r / re;
			dPhi2 = -(exp(-beta * (tempR-1)) * ( (beta + beta * pow((tempR-lambda),20) + 20 * pow((tempR-lambda),19)) / (1 + pow((tempR-lambda),40) + 2*pow((tempR-lambda),20)) ) ) / re;
			if(r<1.6){
				dPhi = -100 + 36.1698 * r;
			}
			else{
			dPhi1 = -(exp(-alpha * (tempR-1)) * ( (alpha + alpha * pow((tempR-kappa),20) + 20 * pow((tempR-kappa),19)) / (1 + pow((tempR-kappa),40) + 2*pow((tempR-kappa),20)) )) / re;
			dPhi = (r > rc)?0:A * dPhi1 - B * dPhi2;
			}
			drou = fe * dPhi2;
			drou = (r>rc)?0:(rou[i] + rou[j]) * drou;
			FK = -dPhi - drou;
			
			FX[i]=FX[i]+FK*(x[i]-x[j])/r;
			FX[j]=FX[j]-FK*(x[i]-x[j])/r;
			FY[i]=FY[i]+FK*(y[i]-y[j])/r;
			FY[j]=FY[j]-FK*(y[i]-y[j])/r;
			FZ[i]=FZ[i]+FK*(z[i]-z[j])/r;
			FZ[j]=FZ[j]-FK*(z[i]-z[j])/r;
		}
	}

	maxF = 0;
	for(i = 0;i < 3 * N; i++)
	{
	/*	if(FF[i]>10)
			FF[i] = 10;
		if(FF[i]<-10)
			FF[i] = -10;*/
		maxF = (FF[i]>maxF)?FF[i]:maxF;
	}
	
	free(rou);
	return maxF;
}

PEnergy PE_Tool::GetPEnergy(PotentialEnergyType type)
{
	switch(type)
	{
	case PotentialEnergyType_FS:
		return FS_E;
		break;
	case PotentialEnergyType_Johnson:
		return Johnson_E;
		break;
	}
	return NULL;
}

PForce PE_Tool::GetPForce(PotentialEnergyType type)
{
	switch(type)
	{
	case PotentialEnergyType_FS:
		return FS_F;
		break;
	case PotentialEnergyType_Johnson:
		return Johnson_F;
		break;
	}
	return NULL;
}