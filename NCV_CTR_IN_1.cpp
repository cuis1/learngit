#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctime>
#include <math.h>

#define MAX(a,b) ((a)>(b)?(a):(b))
#define SIZE  120			//��Ⱥ��ģ	
#define MAXGEN  50			//����������		
#define P_CORSS 0.75		//�������
#define P_MUTATION 0.05		//�������
#define P_CTR   0.33333333	//���ƶˣ�CTR��������ѡ�����ߡ��ϵĸ���
#define LIN 3
#define RAN 10				//LIN��RAN��Ⱦɫ��������к���
#define ORAN 8				// ORAN = pow(2.0,LIN) �ɹ��������������
#define NL 2
#define NR 3				//NCV�ſ����飨����
#define L NL-1
#define R NR-1					//L,R ��Ϊ��������������趨��

FILE *galog;

typedef struct node		//��Ⱥ�ṹ��
{
  //char *gen[LIN][RAN];	//��Ⱥ�еĸ��壨�Ի����ʾ��
	int gen[LIN][RAN];
	int n_gen[LIN][RAN];	//����Ⱦɫ��
	int output_1[LIN][ORAN];		//ԭʼ����ʱ������
	int output_2[LIN][ORAN];		//��ԭʼ���Ϊ����ʱ��������
	int mid[LIN];			//�м������
  //int	input[LIN],
	//	med[NR];
		//o_input[NR];
	int fitness,		//�����Ӧ�ȣ��������е�ĿǰΪֹ���������ֵ��
		fitsum;		//��ǰ��Ⱥ��Ӧ���ܺ�
}node;

node chr [SIZE],		//��ǰ��Ⱥ����
	 next [SIZE],		//��һ����Ⱥ����
	 max,				//��Ӧ�����ĸ���
	 min;				//��Ӧ����С�ĸ���

int input[LIN][ORAN];

//���ɲ���ӡ��������
void c_input()
{
	fprintf(galog,"\n��Ҫ���������\n");
	//����
	//int input[LIN][ORAN];
	for (int i=0;i<ORAN;i++)
	{
		int temp = i;
		for(int j=0;j<LIN;j++)
		{
			input[j][i] = temp%2;
			temp=temp>>1;
		}
	}
	//��ӡ
	for (int i=0;i<LIN;i++)
	{

		for(int j=0;j<ORAN;j++)
		{
			if (j%ORAN == 0)
				fprintf(galog,"\n");
			fprintf(galog,"%2d",input[i][j]);
		}
	}
	fprintf(galog,"\n*******************************\n");
}

//����0~1֮������С��
double randd()			
{
  return (double)rand()/RAND_MAX;
}

//����0~k֮����������
int randi(int k)		
{
	//+0.5��Ϊ�˱�֤��ȡ���ϱ߽��ϵ���
  return (int)(randd()*k+0.5);
}

//��ӡ
void print()
{
	for(int i=0;i<SIZE;i++)
	{
		fprintf(galog,"\n*******************************\n");
		fprintf(galog,"\n%d\n",i);
	//���ƻ���
	//if (i%10 == 0)			
	//	fprintf(galog,"\n");
		for(int k = 0; k < LIN;k++)
		{if (k%LIN == 0)
					fprintf(galog,"\n*******************************\n");
			for(int j = 0; j < RAN;j++)
			{ 
				if (j%RAN == 0)			   
					fprintf(galog,"\n");
				fprintf(galog,"%3d",chr[i].gen[k][j]);
				//switch (chr[i].gen[k][j])
				//	{
				//	case 0:
				//		fprintf(galog,"ZE  ");
				//		break;
				//	case 1:
				//		fprintf(galog,"CTR ");
				//		break;
				//	case 2:
				//		fprintf(galog,"CV+ ");
				//		break;
				//	case 3:
				//		fprintf(galog,"CV  ");
				//		break;
				//	case 4:
				//		fprintf(galog,"CN  ");
				//		break;
				//	default:
				//		fprintf(galog,"Null ");
				//		break;
				//}	
			 }
		}
	}
}
//void print_output()
//{
//	for(int k = 0; k < LIN;k++)
//	{
//		if (k%LIN == 0)
//				fprintf(galog,"\n*******************************\n");
//		for(int j = 0; j < ORAN;j++)
//		{ 
//			if (j%ORAN == 0)			   
//				fprintf(galog,"\n");
//			fprintf(galog,"%3d",chr[i].gen[k][j]);
//				//fprintf(galog,"%4s",chr[i].gen[k][j]);
//		 }
//	}
//}
void cal_fitness_2(int i)
{
	//for (int i = 0; i < SIZE; i++)
	//{
		chr[i].fitness = 0;
		chr[i].fitsum = 0;
		for (int a = 0; a < LIN; a++)
		{
			for (int b = 0; b < RAN; b++)
			{
				chr[i].n_gen[a][RAN-b-1] = chr[i].gen[a][b];
			}
		}
		for (int j = 0;j < ORAN; j++)
		{
			for (int m = 0; m < LIN; m++)
			{
				chr[i].mid[m] = chr[i].output_1[m][j];
			}
			for (int k = 0; k < RAN; k++)
			{	
				int temp[3] = {-1};
				int ctr = -1, ze = -1, ncv = -1;
				for (int l = 0; l < LIN; l++)				
				{
					if (chr[i].n_gen[l][k] == 1)
						ctr = l;
					else if (chr[i].n_gen[l][k] == 0)
						ze = l;
					else
						ncv = l;
				}
				chr[i].output_2[ze][j] = chr[i].mid[ze];
				if (chr[i].mid[ctr] == 0)
				{
					chr[i].output_2[ctr][j] = chr[i].mid[ctr];
					chr[i].output_2[ncv][j] = chr[i].mid[ncv];
				}
				//�������λ����Ϊ1ʱ�ܿ������
				else if(chr[i].mid[ctr] == 1)
				{
					//2��ʾCV+��3��ʾCV��4��ʾCN
					switch (chr[i].n_gen[ncv][k])
					{
						case 0:
							fprintf(galog,"chr[%d].n_gen[%d][%d] case 0 error\n",i,ncv,k);
							break;
						case 1:
							fprintf(galog,"chr[%d].n_gen[%d][%d] case 1 error\n",i,ncv,k);
							break;
						case 2:
							//2��ʾV0��3��ʾV1
							switch (chr[i].mid[ncv])
							{
								case 0:
									chr[i].output_2[ncv][j] = 3;
									chr[i].mid[ncv] = 3;
									break;
								case 1:
									chr[i].output_2[ncv][j] = 2;
									chr[i].mid[ncv] = 2;
									break;
								case 2:
									chr[i].output_2[ncv][j] = 0;
									chr[i].mid[ncv] = 0;
									break;
								case 3:
									chr[i].output_2[ncv][j] = 1;
									chr[i].mid[ncv] = 1;
									break;
							default:
								fprintf(galog,"chr[%d].mid[%d]=%d case 2 default error\n",i,ncv,chr[i].mid[ncv]);
								break;
							}break;
							case 3:
							//2��ʾV0��3��ʾV1
							switch (chr[i].mid[ncv])
							{
								
								case 0:
									chr[i].output_2[ncv][j] = 2;
									chr[i].mid[ncv] = 2;
									break;
								case 1:
									chr[i].output_2[ncv][j] = 3;
									chr[i].mid[ncv] = 3;
									break;
								case 2:
									chr[i].output_2[ncv][j] = 1;
									chr[i].mid[ncv] = 1;
									break;
								case 3:
									chr[i].output_2[ncv][j] = 0;
									chr[i].mid[ncv] = 0;
									break;
							default:
								fprintf(galog,"chr[%d].mid[%d]=%d case 3 default error\n",i,ncv,chr[i].mid[ncv]);
								break;
							}break;
							case 4:
							//2��ʾV0��3��ʾV1
							switch (chr[i].mid[ncv])
							{
								case 0:
									chr[i].output_2[ncv][j] = 1;
									chr[i].mid[ncv] = 1;
									break;
								case 1:
									chr[i].output_2[ncv][j] = 0;
									chr[i].mid[ncv] = 0;
									break;
								case 2:
									chr[i].output_2[ncv][j] = 3;
									chr[i].mid[ncv] = 3;
									break;
								case 3:
									chr[i].output_2[ncv][j] = 2;
									chr[i].mid[ncv] = 2;
									break;
							default:
								fprintf(galog,"chr[%d].mid[%d]=%d case 4 default error\n",i,ncv,chr[i].mid[ncv]);
								break;
							}break;
					default:
						fprintf(galog,"chr[%d].n_gen[%d][%d] default error\n",i,ncv,k);
						break;
					}
				}
				//����������������Ϊ����ֵ����Ϊ0��1��ʱ�����
				else 
				{   
					k = RAN - 1;
					for (int a = 0; a < LIN; a++)
					{
						chr[i].output_2[a][j] = -1;
					}
					////Ϊ�˵��Լӵ�ע��
					//fprintf(galog,"\nchr[%d].mid[%d]=%d  \n",i,ctr,chr[i].mid[ctr]);
					//fprintf(galog,"\nchr[%d].n_gen[%d][%d]=%d \n",i,ctr,k,chr[i].n_gen[ctr][k]);
					//int temp = -1;
					//for (int g = 2; g >= 0; g--)
					//{
					//	if (chr[i].mid[g] == 0||chr[i].mid[g] == 1)
					//	{
					//		//Ϊ�˵��Լӵ�ע��
					//		fprintf(galog,"\n����chr[%d].mid[%d]=%d  \n",i,g,chr[i].mid[g]);
					//		fprintf(galog,"\nCTR����chr[%d].n_gen[%d][%d]=%d \n",i,g,k,chr[i].n_gen[g][k]);
					//		temp = chr[i].gen[ctr][k];
					//		chr[i].gen[ctr][k] = chr[i].gen[g][k];
					//		chr[i].gen[g][k] = temp;
					//		k =-1;
					//		for (int m = 0; m < LIN; m++)
					//		{
					//			chr[i].mid[m] = input[m][j];
					//		}
					//		for (int l = 0; l < LIN; l++)				
					//		{
					//			if (chr[i].gen[l][k] == 1)
					//				ctr = l;
					//			else if (chr[i].gen[l][k] == 0)
					//				ze = l;
					//			else
					//				ncv = l;
					//		}
					//		break;
					//	}
					//}
				}
			}
		}
		//��ӡԭʼ����ʱ������
		//for(int k = 0; k < LIN;k++)
		//{
		//	for(int j = 0; j < ORAN;j++)
		//	{ 
		//		if (j%ORAN == 0)			   
		//			fprintf(galog,"\n");
		//		fprintf(galog,"%3d",chr[i].output_2[k][j]);
		//			//fprintf(galog,"%4s",chr[i].gen[k][j]);
		//	 }
		//	if (k%LIN == 2)
		//		fprintf(galog,"\n*******************************\n");
		//}
		for (int c = 0; c < ORAN; c++)
		{
			int temp = 0, fit = 0;
			for (int d = 0; d < LIN; d++)
			{
				if (chr[i].output_2[d][c] == input[d][c])
					temp++;
			}
			if (temp == 3)
				chr[i].fitness = chr[i].fitness + 2;   
			for (int e = 0; e < LIN; e++)
			{
				if (chr[i].output_1[e][c] == 0||chr[i].output_1[e][c] == 1)
					fit++;
			}
			chr[i].fitness = chr[i].fitness + fit-1;
		}
		chr[i].fitsum=i>0?(chr[i].fitness+chr[i-1].fitsum):(chr[0].fitness);
		fprintf(galog,"chr[%d].fitness=%d default error\n",i,chr[i].fitness);
		fprintf(galog,"chr[%d].fitsum=%d default error\n",i,chr[i].fitsum);
		fprintf(galog,"\n*******************************\n");
	//}
}

//����ԭʼ����ʱ������
void cal_fitness_1()
{
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0;j < ORAN; j++)
		{
			for (int m = 0; m < LIN; m++)
			{
				chr[i].mid[m] = input[m][j];
			}
			for (int k = 0; k < RAN; k++)
			{	
				int temp[3] = {-1};
				int ctr = -1, ze = -1, ncv = -1;
				for (int l = 0; l < LIN; l++)				
				{
					if (chr[i].gen[l][k] == 1)
						ctr = l;
					else if (chr[i].gen[l][k] == 0)
						ze = l;
					else
						ncv = l;
				}
				chr[i].output_1[ze][j] = chr[i].mid[ze];
				if (chr[i].mid[ctr] == 0)
				{
					chr[i].output_1[ctr][j] = chr[i].mid[ctr];
					chr[i].output_1[ncv][j] = chr[i].mid[ncv];
				}
				//�������λ����Ϊ1ʱ�ܿ������
				else if(chr[i].mid[ctr] == 1)
				{
					//2��ʾCV+��3��ʾCV��4��ʾCN
					switch (chr[i].gen[ncv][k])
					{
						case 0:
							fprintf(galog,"chr[%d].gen[%d][%d] case 0 error\n",i,ncv,k);
							break;
						case 1:
							fprintf(galog,"chr[%d].gen[%d][%d] case 1 error\n",i,ncv,k);
							break;
						case 2:
							//2��ʾV0��3��ʾV1
							switch (chr[i].mid[ncv])
							{
								case 0:
									chr[i].output_1[ncv][j] = 3;
									chr[i].mid[ncv] = 3;
									break;
								case 1:
									chr[i].output_1[ncv][j] = 2;
									chr[i].mid[ncv] = 2;
									break;
								case 2:
									chr[i].output_1[ncv][j] = 0;
									chr[i].mid[ncv] = 0;
									break;
								case 3:
									chr[i].output_1[ncv][j] = 1;
									chr[i].mid[ncv] = 1;
									break;
							default:
								fprintf(galog,"chr[%d].mid[%d]=%d case 2 default error\n",i,ncv,chr[i].mid[ncv]);
								break;
							}break;
							case 3:
							//2��ʾV0��3��ʾV1
							switch (chr[i].mid[ncv])
							{
								
								case 0:
									chr[i].output_1[ncv][j] = 2;
									chr[i].mid[ncv] = 2;
									break;
								case 1:
									chr[i].output_1[ncv][j] = 3;
									chr[i].mid[ncv] = 3;
									break;
								case 2:
									chr[i].output_1[ncv][j] = 1;
									chr[i].mid[ncv] = 1;
									break;
								case 3:
									chr[i].output_1[ncv][j] = 0;
									chr[i].mid[ncv] = 0;
									break;
							default:
								fprintf(galog,"chr[%d].mid[%d]=%d case 3 default error\n",i,ncv,chr[i].mid[ncv]);
								break;
							}break;
							case 4:
							//2��ʾV0��3��ʾV1
							switch (chr[i].mid[ncv])
							{
								case 0:
									chr[i].output_1[ncv][j] = 1;
									chr[i].mid[ncv] = 1;
									break;
								case 1:
									chr[i].output_1[ncv][j] = 0;
									chr[i].mid[ncv] = 0;
									break;
								case 2:
									chr[i].output_1[ncv][j] = 3;
									chr[i].mid[ncv] = 3;
									break;
								case 3:
									chr[i].output_1[ncv][j] = 2;
									chr[i].mid[ncv] = 2;
									break;
							default:
								fprintf(galog,"chr[%d].mid[%d]=%d case 4 default error\n",i,ncv,chr[i].mid[ncv]);
								break;
							}break;
					default:
						fprintf(galog,"chr[%d].gen[%d][%d] default error\n",i,ncv,k);
						break;
					}
				}
				//����������������Ϊ����ֵ����Ϊ0��1��ʱ�����
				else 
				{
					/*Ϊ�˵��Լӵ�ע��
					fprintf(galog,"\nchr[%d].mid[%d]=%d  \n",i,ctr,chr[i].mid[ctr]);
					fprintf(galog,"\nchr[%d].gen[%d][%d]=%d \n",i,ctr,k,chr[i].gen[ctr][k]);*/
					int temp = -1;
					for (int g = 2; g >= 0; g--)
					{
						if (chr[i].mid[g] == 0||chr[i].mid[g] == 1)
						{
							//Ϊ�˵��Լӵ�ע��
							//fprintf(galog,"\n����chr[%d].mid[%d]=%d  \n",i,g,chr[i].mid[g]);
							//fprintf(galog,"\nCTR����chr[%d].gen[%d][%d]=%d \n",i,g,k,chr[i].gen[g][k]);
							temp = chr[i].gen[ctr][k];
							chr[i].gen[ctr][k] = chr[i].gen[g][k];
							chr[i].gen[g][k] = temp;
							k =-1;
							for (int m = 0; m < LIN; m++)
							{
								chr[i].mid[m] = input[m][j];
							}
							//for (int l = 0; l < LIN; l++)				
							//{
							//	if (chr[i].gen[l][k] == 1)
							//		ctr = l;
							//	else if (chr[i].gen[l][k] == 0)
							//		ze = l;
							//	else
							//		ncv = l;
							//}
							break;
						}
					}
				}
			}
		}
		//��ӡԭʼ����ʱ������
		for(int k = 0; k < LIN;k++)
		{
			for(int j = 0; j < ORAN;j++)
			{ 
				if (j%ORAN == 0)			   
					fprintf(galog,"\n");
				fprintf(galog,"%3d",chr[i].output_1[k][j]);
					//fprintf(galog,"%4s",chr[i].gen[k][j]);
			 }
			if (k%LIN == 2)
				fprintf(galog,"\n*******************************\n");
		}
		cal_fitness_2(i);
	}
}

//���ɳ���Ⱦɫ�岢����
void init()
{
	//char *NCV[NL][NR] = {{"CTR","CTR","CTR"},{"CV+","CV","CN"}};
	char NCV[NL][NR] = {{1,1,1},{2,3,4}};
	//fprintf(galog,"������Ⱥ\n");
	int ran,lin,take;
	double p;
	for(int i=0;i<SIZE;i++)
	{
		for(int j = 0; j < RAN;j++)
		{   
			p = randd(); //���С��
			take = randi(L);//2��3��ȡֵ�ж�ֵ
			//lin = randi(L);//���ȡNCV����
			ran = randi(R);//���ȡNCV����
			int k = 0;
			if (p <= P_CTR)
			{
				chr[i].gen[k][j] = NCV[k][k];
				if (take == 1)
				{
					chr[i].gen[k+1][j] = NCV[take][ran];
					chr[i].gen[k+2][j] = 0;
				}
				else
				{
					chr[i].gen[k+1][j] = 0;
					chr[i].gen[k+2][j] = NCV[take+1][ran];
				}
			}
			else
			{
				chr[i].gen[k][j] = 0;
				if (take == 1)
				{
					chr[i].gen[k+1][j] = NCV[k][k];
					chr[i].gen[k+2][j] = NCV[k+1][ran];
				}
				else
				{
					chr[i].gen[k+1][j] = NCV[k+1][ran];
					chr[i].gen[k+2][j] = NCV[k][k];
				}
			}
		/*	
			k = randi(L)+1;
			chr[i].gen[k][j] = NCV[lin^1][ran];
			if (k == 1)
				chr[i].gen[k+1][j] = 0;
			 else
				chr[i].gen[k-1][j] = 0;*/
		}
	}
   // print();
  cal_fitness_1();
  //cal_fitness_2();
  print();
}
  int main()
{
  srand((unsigned)time(NULL));
  galog = fopen("galog.txt","w+");
  c_input();
  init();
  //GA();
  system("pause");
  return 0;
}