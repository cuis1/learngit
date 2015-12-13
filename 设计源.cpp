#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctime>
#include <math.h>


#define SIZE  120			//��Ⱥ��ģ	
#define MAXGEN  50			//����������		
#define P_CORSS 0.75		//�������
#define P_MUTATION 0.05		//�������
#define LIN 3
#define RAN 10				//LIN��RAN��Ⱦɫ��������к���
#define NL 2
#define NR 3				//NCV�ſ����飨����
#define L NL-1
#define R NR-1					//L,R ��Ϊ��������������趨��

FILE *galog;

char *NCV[NL][NR] = {{"CTR","CTR","CTR"},{"CV+","CV","CN"}};

typedef struct node		//��Ⱥ�ṹ��
{
  char *gen[LIN][RAN];	//��Ⱥ�еĸ��壨�Ի����ʾ��
  int	input[NR],
		output[NR],
		med[NR];
		//o_input[NR];
  double fitness,		//�����Ӧ�ȣ��������е�ĿǰΪֹ���������ֵ��
		 fitsum;		//��ǰ��Ⱥ��Ӧ���ܺ�
}node;

node chr [SIZE],		//��ǰ��Ⱥ����
	 next [SIZE],		//��һ����Ⱥ����
	 max,				//��Ӧ�����ĸ���
	 min;				//��Ӧ����С�ĸ���

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
	//���ƻ���
	//if (i%10 == 0)			
	//	fprintf(galog,"\n");
		for(int k = 0; k < LIN;k++)
		{if (k%3 == 0)
					fprintf(galog,"\n*******************************\n");
			for(int j = 0; j < RAN;j++)
			{ 
				if (j%10 == 0)			
					fprintf(galog,"\n");
		
				fprintf(galog,"%4s",chr[i].gen[k][j]);
			 }
		}
	}
}
//���ɳ���Ⱦɫ�岢����
void init()
{
	fprintf(galog,"������Ⱥ\n");
  int tmp,ran,lin;
  for(int i=0;i<SIZE;i++)
  {
    for(int j = 0; j < RAN;j++)
    {   
	  lin = randi(L);//���ȡ��
	  ran = randi(R);//���ȡ��
		int k = 0;
      chr[i].gen[k][j] = NCV[lin][ran];
	  k = randi(L)+1;
	  chr[i].gen[k][j] = NCV[lin^1][ran];
	  if (k == 1)
		  chr[i].gen[k+1][j] = "ZE";
	  else
		  chr[i].gen[k-1][j] = "ZE";
	}
  //cal_fitness();
	}
  print();
}
  int main()
{
  srand((unsigned)time(NULL));
  galog = fopen("galog.txt","w+");
  init();
  //GA();
  
  system("pause");
  return 0;
}