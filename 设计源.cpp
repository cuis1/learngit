#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctime>
#include <math.h>


#define SIZE  120			//种群规模	
#define MAXGEN  50			//最大迭代次数		
#define P_CORSS 0.75		//交叉概率
#define P_MUTATION 0.05		//变异概率
#define LIN 3
#define RAN 10				//LIN，RAN是染色体数组的行和列
#define ORAN 8				// ORAN = pow(2.0,LIN) 可供输入的数据组数
#define NL 2
#define NR 3				//NCV门库数组（基因）
#define L NL-1
#define R NR-1					//L,R 是为了生成随机数而设定的

FILE *galog;

typedef struct node		//种群结构体
{
  char *gen[LIN][RAN];	//种群中的个体（以基因表示）
  int output[LIN][RAN];		//输出数组
  int mid[LIN];			//中间计算结果
  //int	input[LIN],
	//	med[NR];
		//o_input[NR];
  double fitness,		//最大适应度（程序运行到目前为止函数的最大值）
		 fitsum;		//当前种群适应度总和
}node;

node chr [SIZE],		//当前种群数组
	 next [SIZE],		//下一代种群数组
	 max,				//适应度最大的个体
	 min;				//适应度最小的个体

//生成并打印输入数组
void c_input()
{
	fprintf(galog,"\n需要输入的数据\n");
	//生成
	int input[LIN][ORAN];
	for (int i=0;i<ORAN;i++)
	{
		int temp = i;
		for(int j=0;j<LIN;j++)
		{
			input[j][i] = temp%2;
			temp=temp>>1;
		}
	}
	//打印
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

//产生0~1之间的随机小数
double randd()			
{
  return (double)rand()/RAND_MAX;
}

//产生0~k之间的随机整数
int randi(int k)		
{
	//+0.5是为了保证能取到上边界上的数
  return (int)(randd()*k+0.5);
}

//打印
void print()
{
	for(int i=0;i<SIZE;i++)
	{
	//控制换行
	//if (i%10 == 0)			
	//	fprintf(galog,"\n");
		for(int k = 0; k < LIN;k++)
		{if (k%LIN == 0)
					fprintf(galog,"\n*******************************\n");
			for(int j = 0; j < RAN;j++)
			{ 
				if (j%RAN == 0)			   
					fprintf(galog,"\n");
		
				fprintf(galog,"%4s",chr[i].gen[k][j]);
			 }
		}
	}
}
//生成初代染色体并编码
void init()
{
	char *NCV[NL][NR] = {{"CTR","CTR","CTR"},{"CV+","CV","CN"}};
	fprintf(galog,"初代种群\n");
	int ran,lin;
	for(int i=0;i<SIZE;i++)
	{
		for(int j = 0; j < RAN;j++)
		{   
			lin = randi(L);//随机取NCV的行
			ran = randi(R);//随机取NCV的列
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
  c_input();
  init();
  //GA();
  
  system("pause");
  return 0;
}