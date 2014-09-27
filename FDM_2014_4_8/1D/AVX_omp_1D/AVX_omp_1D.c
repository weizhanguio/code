/*
  purpuse: finite difference scheme for variable-order time fractional diffusion equation, fractional order is a function of space x and time t.
  reference paper:Sun Hongguang  "finite difference scheme for variable-order time fractional diffusion equation", page 1250085-12, example 1.
  date: 2014.3.5
  author: Zhang Wei
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<immintrin.h>
#include<time.h>  
#include<omp.h>
#define  K             0.000001                   /* diffusion efficient*/
#define  pi            3.1415926
#define  x_start       0.0                    /*starting point of space*/
#define  x_end         10.0                 /*end point of space*/
#define  t_start       0.0                    /*starting point of time*/
#define  t_end         1.0   	              /*end point of time*/
#define  L           (x_end-x_start)
#define  t(x)          (t_start+dt*(x))
#define  x(a)          (x_start+dx*(a))
#define  r(x)         r[x]
#define  s(x)         s[x] 
#define  alpha(x)       alpha[(x)]
#define  powfunction(x)        powfunction[(x)]
#define  gammafunction(x)      gammafunction[(x)]

#define  Jfunction(alpha,k)    pow((k),1.0-(alpha))   /*compute the coef*/
#define  real      double 





 
int 
main(int nargs, char** args)
{
	
	
  int i;
  int j;
  int k;
  int m;
  int n;
  
  	              	
  real begin;
  real finish;
  int  x_N, t_N;
  x_N=200;                   //number of space intervals
  t_N=200;                   //number of time intervals


  if(nargs>1){
    x_N=atof(args[1]);
    t_N=atof(args[2]);

  }
    
  int BlockWidth_x=4; 
  int BlockWidth_t=4; 
  int nBlock_x=(x_N-1)/BlockWidth_x;
  int nBlock_t=(t_N-1)/BlockWidth_t;
  int rem_x=(x_N-1)%BlockWidth_x;
  int rem_t=(t_N-1)%BlockWidth_t;
  int NextStart=x_N-1+4-(x_N-1)%4+1;    // padding, add one block 
  int x_Nnew=NextStart+3;       //next i start position: x_N-1+4-(x_N-1)%4+1
  
  real   dx;        dx= ((x_end-x_start)/x_N);          // space step
  real   dt;        dt= ((t_end-t_start)/t_N);          // time step
  
  printf("------omp  AVX version----------\n");
  printf("%d intervals in space, %d intervals in time \n\n",x_N, t_N);  
  /*x_N+1 points in space, t_N+1 points in time*/
	
	
	
  int threadNum;	
#pragma omp parallel
  {
#pragma omp master
    {
      threadNum=omp_get_num_threads();
      printf("Number of threads = %d\n",threadNum);
    }
  }
  	
	
	

#define  u(x,y)         u[1ul*(x)*(t_N+1)+(y)]
#define  coef(x,y)      coef[1ul*(x)*(t_N+1)+(y)]


  real *u=(real*)_mm_malloc(1ul* (x_Nnew+1)*(t_N+1)*sizeof(real), 64);   
  real *alpha=(real*)_mm_malloc((x_Nnew+1)*sizeof(real),64);
  real *coef=(real*)_mm_malloc(1ul* (t_N+1) *(x_Nnew+1)*sizeof(real),64);  /*coef is d in eq.6  */ 
  real *powfunction=(real*)_mm_malloc((x_Nnew+1)*sizeof(real),64);
  real *gammafunction=(real*)_mm_malloc((x_Nnew+1)*sizeof(real),64);
  real *r=(real*)malloc((x_Nnew+1)*sizeof(real));
  real *s=(real*)malloc((x_Nnew+1)*sizeof(real));   //store the value of 1-(  Jfunction(alpha(i),2)-Jfunction(alpha(i),1) ) 
   
 
   
  real* outputbuffer=(real*)_mm_malloc((x_Nnew+1)*BlockWidth_t*sizeof(real),64);
 

#pragma omp parallel for
  for(i=0;i<(x_Nnew+1)*BlockWidth_t;i++)
    outputbuffer[i]=0.0;



 
#pragma omp parallel for
  for(i=0;i<(x_Nnew+1)*(t_N+1);i++)
    u[i]=0.0;  

  
#pragma omp parallel for  
  for(i=1;i<=x_N-1;i++)    //coefficients alpha	   
    {	             
      alpha(i)= 0.8;
      powfunction(i)=pow(dt,alpha(i));
      gammafunction(i)=tgamma(2-alpha(i));
      r(i)=1.0*K*pow(dt,alpha(i))*tgamma(2-alpha(i))/(dx*dx);
      s(i)=1.0-(  Jfunction(alpha(i),2)-Jfunction(alpha(i),1) );
    }

  for(i=x_N;i<=x_Nnew;i++)    //coefficients alpha	   

    {	             
      alpha(i)= 0.0;
      powfunction(i)=0.0;
      gammafunction(i)=0.0;
      r(i)=0.0;
      s(i)=0.0;
    }	  
	  	  
	   
	 
  real Min=0.0;
#pragma omp parallel for	
  for(i=1;i<=x_N-1;i++)
	 
    {     
      real s1;
      s1= r(i)-0.25*(2- pow(2,1-alpha(i))) ;     
      Min=(s1<Min?s1:Min);
    }
  // printf("Min is %f\n",Min);
  if(Min>=0)
    {
      printf("error, do not satisfy the stability requirement.\n\n");
      
    }
   



#pragma omp parallel for
  for (i=1; i<x_N; i++) 
    for (k=1; k<=t_N; k++)
      {
	coef(i,k) = 2.0*Jfunction(alpha(i),k+1)
	  -Jfunction(alpha(i),k)-Jfunction(alpha(i),k+2);
      }


#pragma omp parallel for
  for(i=1;i<=x_N-1;i++)            /*intial value u(x,0)=0*/
    {
      u(i,0)=1.0;
    }
	
  u(0,0)=0.2;  		  /*left corner*/
  u(x_N,0)=0.2;                /*right corner*/
 
#pragma omp parallel for 
  for(i=1;i<=t_N;i++)
    {
      u(0,i)=0.2;                 /*left boundary*/ 
      u(x_N,i)=0.2;               /*right boundary*/	 	
    }


begin=omp_get_wtime();  

	
#pragma omp parallel for
  for(i=1;i<=x_N-1;i++)             /*compute  u(:,1)*/
    {        
      u(i,1)=(1-2*r(i))*u(i,0)+r(i)*u(i+1,0)+r(i)*u(i-1,0);
               
    }

#pragma omp barrier


  real sum;


  for (k=1; k<5; k++)              /*compute other value */
    for (i=1; i<x_N; i++) {    
      sum=0.0;
      for (j=1; j<k; j++){
	sum += u(i,j) * coef(i,k-j);
      }    
      u(i,k+1)= r(i)*(u(i+1,k)+u(i-1,k)) +(s(i)-2*r(i))*u(i,k) + sum + (Jfunction(alpha(i),(k+1))- Jfunction(alpha(i),k))*u(i,0) ;
                                                 
    }







#pragma omp parallel default(shared) private(k,i,j,m) 
  {    	   
    for (k=5; k<=nBlock_t*BlockWidth_t; k+=4){              /*compute other value */

        int id=omp_get_thread_num();
#pragma omp for 
      for (i=1; i<=NextStart-1; i+=4) {
    
	__m256d  sum1;
	__m256d  sum2;
	__m256d  sum3;
	__m256d  sum4;

	__m256d  load_u1;
	__m256d  load_u2;
	__m256d  load_u3;
	__m256d  load_u4;
 
	__m256d  mulp1;
	__m256d  mulp2;
	__m256d  mulp3;
	__m256d  mulp4;

	__m256d  coef1;
	__m256d  coef2;
	__m256d  coef3;
	__m256d  coef4;
          
	sum1=_mm256_setzero_pd();
	sum2=_mm256_setzero_pd();
	sum3=_mm256_setzero_pd();
	sum4=_mm256_setzero_pd();
	for (j=1; j<k; j++){
       
                
        load_u1=_mm256_broadcast_sd((real*)(u+i*(t_N+1)+j));        /*u(i,:) */
	load_u2=_mm256_broadcast_sd((real*)(u+(i+1)*(t_N+1)+j));    /*u(i+1,:) */
	load_u3=_mm256_broadcast_sd((real*)(u+(i+2)*(t_N+1)+j) );    /*u(i+2,:) */
	load_u4=_mm256_broadcast_sd((real*)(u+ (i+3)*(t_N+1)+j));    /*u(i+3,:) */    
																				
	coef1= _mm256_load_pd(coef+i*(t_N+1)+k-j);
	coef2= _mm256_load_pd(coef+(i+1)*(t_N+1)+k-j);
	coef3= _mm256_load_pd(coef+(i+2)*(t_N+1)+k-j);
	coef4= _mm256_load_pd(coef+(i+3)*(t_N+1)+k-j);       
										



	  mulp1=_mm256_mul_pd(load_u1 ,coef1);
	  mulp2=_mm256_mul_pd(load_u2 ,coef2);
	  mulp3=_mm256_mul_pd(load_u3 ,coef3);
	  mulp4=_mm256_mul_pd(load_u4 ,coef4);       
        
	  sum1=_mm256_add_pd(sum1,mulp1); 
	  sum2=_mm256_add_pd(sum2,mulp2);
	  sum3=_mm256_add_pd(sum3,mulp3); 
	  sum4=_mm256_add_pd(sum4,mulp4);       
 

    
	}    

       
	
        _mm256_store_pd(&outputbuffer[i*4],sum1);
	_mm256_store_pd(&outputbuffer[(i+1)*4],sum2);
	_mm256_store_pd(&outputbuffer[(i+2)*4],sum3);
	_mm256_store_pd(&outputbuffer[(i+3)*4],sum4);   
  
   
      



 
	for(m=0;m<4;m++)
	  {   
	   
	    u(i+m,k+1)= r(i+m)*(u(i+m+1,k)+u(i+m-1,k)) +(1-(  Jfunction(alpha(i+m),2)-Jfunction(alpha(i+m),1) )-2*r(i+m))*u(i+m,k) + (outputbuffer[(i+m)*4]) + (Jfunction(alpha(i+m),(k+1))- Jfunction(alpha(i+m),k))*u(i+m,0);
      
	                  
        
	  }       
                      
      }
 
#pragma omp single
      {
	u(x_N,k+1)=0.2;
      }
  
  
  
      /*compute u(i,k+2)*/

#pragma omp for 
      for (i=1; i<=x_N-1; i++) { 
           
	outputbuffer[i*4+1]+=u(i,k)*coef(i,1);

	u(i,k+2)=r(i)*(u(i+1,k+1)+u(i-1,k+1)) +(s(i)-2*r(i))*u(i,k+1) + outputbuffer[i*4+1] + (Jfunction(alpha(i),(k+2))- Jfunction(alpha(i),k+1))*u(i,0);
	      
   
      }
 
  
#pragma omp barrier 
 
 
    
      /*compute u(i,k+3)*/
#pragma omp for 
      for (i=1; i<=x_N-1; i++) { 
	
	outputbuffer[i*4+2]+=u(i,k)*coef(i,2)+u(i,k+1)*coef(i,1);
	u(i,k+3)=r(i)*(u(i+1,k+2)+u(i-1,k+2)) +(s(i)-2*r(i))*u(i,k+2) +outputbuffer[i*4+2] + (Jfunction(alpha(i),(k+3))- Jfunction(alpha(i),k+2))*u(i,0);
        
    
      }
    
#pragma omp barrier  
      /*compute u(i,k+4)*/
    
    
    
#pragma omp for 
      for (i=1; i<=x_N-1; i++) { 
      	
	outputbuffer[i*4+3]+=u(i,k)*coef(i,3)+u(i,k+1)*coef(i,2)+u(i,k+2)*coef(i,1);
	u(i,k+4)=r(i)*(u(i+1,k+3)+u(i-1,k+3)) +(s(i)-2*r(i))*u(i,k+3) +outputbuffer[i*4+3]+ (Jfunction(alpha(i),(k+4))- Jfunction(alpha(i),k+3))*u(i,0);
    
    
    
      }
#pragma omp barrier    
    
      /*compute remains of u(i,k+4)*/
    
    }
  }
  /*compute the remains in time*/


  real sum_t;

  for (k=nBlock_t*BlockWidth_t; k<t_N; k++)              /*compute other value */
    {
#pragma omp parallel for  
      for (i=1; i<x_N; i++) {    
	sum_t=0.0;
	for (j=1; j<k; j++){
	  sum_t += u(i,j) * coef(i,k-j);
	}    
	u(i,k+1)= r(i)*(u(i+1,k)+u(i-1,k)) +(s(i)-2*r(i))*u(i,k) + sum_t + (Jfunction(alpha(i),(k+1))- Jfunction(alpha(i),k))*u(i,0) ;
                                  
      }
#pragma omp barrier
    }
                  
  finish=omp_get_wtime();
  
  real time=1.0*(finish-begin);
     	
 	
  real numer=0.0;
                     
  for(i=1;i<=x_N-1;i++)    
    for(k=1;k<=t_N;k++)     
      {  
	numer +=u(i,k)*u(i,k);
	 
      }
     	
  printf("x_N: %d, t_N: %d ,Norm_n:  %e,threadNum: %d, timeusage is %6.2f seconds\n\n",x_N,t_N,sqrt(numer),threadNum, time);


  _mm_free(u);    
  _mm_free(alpha);      
  _mm_free(coef);
  _mm_free(powfunction); 
  _mm_free(gammafunction );
  _mm_free(outputbuffer);
  free(r);
  free(s);
  
 
  return EXIT_SUCCESS;     
}

