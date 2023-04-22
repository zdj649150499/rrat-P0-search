// also needed for this example...
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "cpgplot.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_statistics.h"

#define PI 3.141592653589793238462643383279

void help();
void Getdifference(double *in, double *out, double P0, int N);
void cutExtName(char *name,char *out);
int main(int argc, char* argv[])
{
  int i,j,k;
  char filename[1000];
  int MJDNUM=0;
  int M ;            // number of nonuniform points
  int N ;            // number of modes
  int def1mark=0;
  int def0mark=0;
  int bpmark=0;
  double bpl,bpr,bpstep;
  double minsig,maxsig;
  double bestP0;
  int DIFFNUM;

  FILE *fp;

  if(argc < 2) help();
  for(i=0;i<argc;i++)
  {
    if (strcmp(argv[i],"-h")==0) help();
    else if(strcmp(argv[i],"-bp")==0)
    {
      bpmark=1;
      sscanf(argv[++i],"%lf,%lf,%lf",&bpl,&bpr,&bpstep);
    }
    else if(strcmp(argv[i],"-df1")==0)
      def1mark=1;
    else if(strcmp(argv[i],"-df0")==0)
      def0mark=1;
  }

  strcpy(filename,argv[argc-1]);//获取文件名
  printf("Reading MJD from file: %s\n",filename);

  fp=fopen(filename,"r");
  while(!feof(fp))                  //获取list中的mjd数
	{
		if(fgetc(fp)=='\n')
		{
			MJDNUM++;
		}
	}
  printf("MJDs num: %d\n",MJDNUM);
  fclose(fp);
  //rewind(fp);

  M=MJDNUM*(MJDNUM-1)/2;                                            //两两相减后的数量
  N=(bpr-bpl)/bpstep + 1;                                           //周期数
  if(def0mark==0) DIFFNUM=M;                                        //使用两两相减得到sigma
  else DIFFNUM=MJDNUM-1;                                            //直接获取MJDs的减0号MJD得到sigma

  double* MJDs = (double *)malloc(sizeof(double)*MJDNUM);
  double* x = (double *)malloc(sizeof(double)*M);                       //两两相减
  double* x1 = (double *)malloc(sizeof(double)*(M-1)*M/2);              //两两相减，再两两相减
  double* period = (double *)malloc(sizeof(double)*N);
  double* Sigma = (double *)malloc(sizeof(double)*N);
  double* diff = (double *)malloc(sizeof(double)*DIFFNUM);          //每个周期中的点数=M-1
  float * xflot = (float *)malloc(sizeof(float)*M);
  double* diffMJDs=(double *)malloc(sizeof(double)*(MJDNUM-1));


  fp=fopen(filename,"r");
  for ( i = 0; i < MJDNUM; i++)
  {
    fscanf(fp,"%lf",&MJDs[i]);
    if(bpmark ==  0) printf("%.15f\n",MJDs[i]);

  }
  fclose(fp);
  gsl_sort(MJDs,1,MJDNUM);

  for(i=0;i<MJDNUM-1;i++)
  {
    diffMJDs[i] = (MJDs[i+1]-MJDs[0]);                          //获取MJDs的减0号MJD
  }

  k=0;
  int l=0;
  for ( i = l; i < MJDNUM-1; i++)
  {
    for (j = l+1; j< MJDNUM; j++)
    {
      x[k] = (MJDs[j] - MJDs[i])*24*3600;
      k++;
    }
    l++;
  }
  gsl_sort(x,1,M);
  for (k = 0; k < M; k++)
  {
    xflot[k] = (float)(x[k]);
  }
  
  l=0;
  k=0;
  for ( i = l; i < M-1; i++)
  {
    for (j = l+1; j< M; j++)
    {
      x1[k] = (x[j] - x[i]);
      k++;
    }
    l++;
  }
  gsl_sort(x1,1,(M-1)*M/2); 
  if(def1mark==1)
  for(i=0;i<(M-1)*M/2;i++)
  {
    printf("%f\n",x1[i]);
  }

  // for (k = 1; k < M; k++)
  // {
  //   x1[k-1] = x[k] - x[k-1];
  // }
  // gsl_sort(x1,1,M-1);

  if(bpmark ==  0)
  {
    for (k = 0; k < M; k++)
    {
      printf("%.15f\n",x[k]);
    }
  }
  else
  {
    printf("Caculate the sigma of each P0...\n");
    printf("Using Left-P0:%f,  Right-P0:%f,  P0-step:%f\n",bpl,bpr,bpstep);
    for(i=0;i<N;i++)
    {
      period[i] = bpl+bpstep*i;                                 //获取周期序列
    }

    printf("\n\n");

    if(def0mark==0)
      for(i=0;i<N;i++)
      {
        Getdifference(x, diff, period[i],DIFFNUM);                 //获取与对比周期的差值的百分比(两两相减的)
        Sigma[i] = sqrt(gsl_stats_variance(diff,1,DIFFNUM));       //获取与对比周期的的标准差
      }
    else
      for(i=0;i<N;i++)
      {
        Getdifference(diffMJDs, diff, period[i],DIFFNUM);                 //获取与对比周期的差值的百分比(获取MJDs的减0号MJD)
        Sigma[i] = sqrt(gsl_stats_variance(diff,1,DIFFNUM));       //获取与对比周期的的标准差
      }
    


    float* periodf = (float *)malloc(sizeof(float)*N);
    float* Sigmaf = (float *)malloc(sizeof(float)*N);
    bestP0 = 0;
    minsig = Sigma[0];
    maxsig = Sigma[0];
    for(i=0;i<N;i++)
    {
      if(minsig > Sigma[i]) 
      {
        bestP0 = period[i];
        minsig = Sigma[i];
        
      }
      if(maxsig < Sigma[i]) 
      {
        maxsig = Sigma[i];
      }
      periodf[i] = (float)period[i];
      Sigmaf[i] = (float)Sigma[i];
    }
    printf("Min Sigma: P0:%.10f  Sigma:%.10f\n",bestP0,minsig);

    Getdifference(x, diff, bestP0,DIFFNUM);                    //获取与对比周期的差值的百分比

    char output[1000];
    strcpy(output,filename);
    cutExtName(output,".");
    if(def0mark==0) 
    {
      printf("Output: %s-P0-sigma.ps\n",output);
      strcat(output,"-P0-sigma.ps/cps");
    }
    else
    {
      printf("Output: %s-df0-P0-sigma.ps\n",output);
      strcat(output,"-df0-P0-sigma.ps/cps");
    }
    

    cpgbeg(0,output,1,1);
    cpgsvp(0.1,0.9,0.1,0.9);
    cpgslw(2);
    cpgswin((float)bpl,(float)bpr,0,(float)(maxsig+(maxsig-minsig)*0.1));
    //cpgpt(N,periodf,Sigmaf,1);
    cpgline(N,periodf,Sigmaf);
    cpgbox("bcnst",0,0,"bcnst",0,0);
    cpglab("P0 (s)","Sigma","");
    cpgend();


    // char output[1000];
    // strcpy(output,filename);
    // printf("Output: %s-hist.ps\n",output);
    // strcat(output,"-hist.ps/cps");

    // cpgbeg(0,output,1,1);
    // cpgsvp(0.1,0.9,0.1,0.9);
    // //cpgswin(0,100,0,10);
    // cpghist(M,xflot,0.0,50.0,100,0);
    // cpgbox("bcnst",0,0,"bcnst",0,0);
    // cpglab("Distance (s)","No.","");
    // cpgend();

    //////////////////////////////////*画MJD-P0图*//////////////////////////////////////////

    // float* MJDsf = (float *)malloc(sizeof(float)*MJDNUM);
    // float* difff = (float *)malloc(sizeof(float)*MJDNUM);
    // double* diffP0 = (double *)malloc(sizeof(double)*MJDNUM);
    // double MJDmin,MJDmax;
    // double diffmin,diffmax;
    // float MJDminf,MJDmaxf;
    // float diffminf,diffmaxf;

    // for(i = 0; i < MJDNUM; i++)
    // {
    //   diffP0[i] = (MJDs[i]-MJDs[0])/bestP0;
    //   if(diffP0[i] > 0.5) diffP0[i] = 1.0 - diffP0[i];
      

    //   difff[i] = diffP0[i];
    //   MJDsf[i] = MJDs[i];
    //   printf("%f %f %f %f\n",MJDs[i],MJDsf[i],diffP0[i],difff[i]);
    // }
    // gsl_stats_minmax(&MJDmin,&MJDmax,MJDs,1,MJDNUM);
    // gsl_stats_minmax(&diffmin,&diffmax,diffP0,1,MJDNUM);
    // diffminf=diffmin;
    // diffmaxf=diffmax;
    // MJDminf=MJDmin;
    // MJDmaxf=MJDmax;
    // printf("min:%f max:%f\n",MJDminf,MJDmaxf);

    // strcpy(output,filename);
    // printf("Output: %s-MJD-P0.ps\n",output);
    // strcat(output,"-MJD-P0.ps/cps");
    // cpgbeg(0,output,1,1);
    // cpgsvp(0.1,0.9,0.1,0.9);
    // cpgslw(2);
    // cpgswin(MJDminf-(MJDmaxf-MJDminf)*0.1,MJDmaxf+(MJDmaxf-MJDminf)*0.1,diffminf-(diffmaxf-diffminf)*0.1,diffmaxf+(diffmaxf-diffminf)*0.1);
    // cpgpt(MJDNUM,MJDsf,difff,2);
    // cpgbox("bcnst",0,0,"bcnst",0,0);
    // cpglab("MJD","Residual (s)","");
    // cpgend();


    free(periodf);
    free(Sigmaf);
    // free(MJDsf);
    // free(difff);
    // free(diffP0);
  }
  


  free(MJDs);
  free(x);
  free(x1);
  free(period);
  free(Sigma);
  free(diff);
  free(xflot);
  free(diffMJDs);
  return 0;
}


void help()
{
  printf("# This is a profram for calculate Differences between each MJDs\n  and the Differences between adjacent two of them\n");
  printf("# Usage:\n");
  printf("#   rratbestP0 [options] MJDs-list-file\n");
  printf("#   -bp  a,b,c      Set left-P0 , right-P0 and step-P0 for caculate the sigma of each P0\n");
  printf("#   -df0            Use diff-MJDs0 series for search P0\n");
  // printf("#   -df1            Use diff-diff-MJDs0 series for search P0\n");
  printf("#   -h              Print this help\n");
  exit(0);
}

void Getdifference(double *in, double *out, double P0, int N)
{
  int i;
  for(i=0;i<N;i++)
  {
    out[i] = in[i]/P0;
    out[i] = out[i]- (int)out[i];
    if(out[i] > 0.5 ) out[i] = 1-out[i];
    //out[i] = (out[i] < 0.5) ? out[i]:(out[i]-0.5);
  }
}

void cutExtName(char *name,char *out)
{
  int i,len;
  len=strlen(name);
  for(i=len-1;i>=0;i--)
    if(name[i]==*out)
      break;
  if(i>0)
    name[i]='\0';
}