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
#define NL 2
#define NR 3				//NCV门库数组（基因）
#define L NL-1
#define R NR-1					//L,R 是为了生成随机数而设定的

FILE *galog;

char *NCV[NL][NR] = {{"CTR","CTR","CTR"},{"CV+","CV","CN"}};

typedef struct node		//种群结构体
{
  char *gen[LIN][RAN];	//种群中的个体（以基因表示）
  int	input[NR],
		output[NR],
		med[NR];
		//o_input[NR];
  double fitness,		//最大适应度（程序运行到目前为止函数的最大值）
		 fitsum;		//当前种群适应度总和
}node;

node chr [SIZE],		//当前种群数组
	 next [SIZE],		//下一代种群数组
	 max,				//适应度最大的个体
	 min;				//适应度最小的个体

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
//生成初代染色体并编码
void init()
{
	fprintf(galog,"初代种群\n");
  int tmp,ran,lin;
  for(int i=0;i<SIZE;i++)
  {
    for(int j = 0; j < RAN;j++)
    {   
	  lin = randi(L);//随机取行
	  ran = randi(R);//随机取列
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