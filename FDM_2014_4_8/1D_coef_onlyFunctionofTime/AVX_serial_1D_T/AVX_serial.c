/*
   purpuse: finite difference scheme for variable-order time fractional diffusion equation, fractional order is a function of space x and time t.
  reference paper:Sun Hongguang  "finite difference scheme for variable-order time fractional diffusion equation", page 1250085-12, example 1.
  date: 2014.3.5
  author: Zhang Wei
*/
//#include<papi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<immintrin.h>
#include<time.h>  
#define  K             0.000001                   /* diffusion efficient*/
#define  pi            3.1415926
#define  x_start       0.0                    /*starting point of space*/
#define  x_end         10.0                 /*end point of space*/
#define  t_start       0.0                    /*starting point of time*/
#define  t_end         1.0     	              /*end point of time*/
#define  L           (x_end-x_start)
#define  t(x)          (t_start+dt*(x))
#define  x(a)          (x_start+dx*(a))

#define  alpha        0.8
#define  powfunction(x)        powfunction[(x)]
#define  gammafunction(x)      gammafunction[(x)]

#define  Jfunction(alpha,k)    pow((k),1.0-(alpha))   /*compute the coef*/
#define  real              double 




 
int 
main(int nargs, char** args)
{
	
	
  int i;
  int j;
  int k;
  int m;
  int n;
  
  	              	
  clock_t begin;
  clock_t finish;
  int  x_N, t_N;
  x_N=200;                   //number of space intervals
  t_N=200;                   //number of time intervals


  if(nargs>1){
    x_N=atof(args[1]);
    t_N=atof(args[2]);

  }
/*   
int Events[] = { PAPI_L2_DCR, PAPI_L3_DCR };
    int NUM_EVENTS = sizeof(Events)/sizeof(Events[0]);
    long long res_papi[NUM_EVENTS];
    char EventName[128];
    int num_hwcntrs = 0;
    int EventSet = PAPI_NULL;
    int retval;
 */
  int BlockWidth_x=4; 
  int BlockWidth_t=4; 
  int nBlock_x=(x_N-1)/BlockWidth_x;
  int nBlock_t=(t_N-1)/BlockWidth_t;
  int rem_x=(x_N-1)%BlockWidth_x;
  int rem_t=(t_N-1)%BlockWidth_t;
    
 
  real   dx;        dx= ((x_end-x_start)/x_N);          // space step
  real   dt;        dt= ((t_end-t_start)/t_N);          // time step
  
  printf("---------AVX version----------\n");
  printf("%d intervals in space, %d intervals in time \n\n",x_N, t_N);  
  /*x_N+1 points in space, t_N+1 points in time*/
	

#define  u(x,y)         u[1ul*(x)*(t_N+1)+(y)]
#define  coef(y)      coef[(y)]



  
  real *u=(real*)_mm_malloc(1ul*(x_N+1)*(t_N+1)*sizeof(real), 64);   
 
  real *coef=(real*)_mm_malloc((t_N+1)*sizeof(real),64);  /*coef is d in eq.6  */ 
  real *powfunction=(real*)_mm_malloc((x_N+1)*sizeof(real),64);
  real *gammafunction=(real*)_mm_malloc((x_N+1)*sizeof(real),64);
  
  real *sum1_rem=(real*)malloc((x_N+1)*sizeof(real));
  real *sum2_rem=(real*)malloc((x_N+1)*sizeof(real));
  real *sum3_rem=(real*)malloc((x_N+1)*sizeof(real)); 
  real *sum_rem=(real*)malloc((x_N+1)*sizeof(real));
   
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


  real* outputbuffer=(real*)_mm_malloc((x_N+1)*BlockWidth_t*sizeof(real),64);
  real r=1.0*K*pow(dt,alpha)*tgamma(2-alpha)/(dx*dx);
  real s=1-(pow(2,1.0-(alpha))-pow(1,1.0-(alpha)));   //store the value of 1-(  Jfunction(alpha,2)-Jfunction(alpha,1) )	
 
  /*
 retval = PAPI_library_init( PAPI_VER_CURRENT );
     retval = PAPI_create_eventset( &EventSet );

     if (PAPI_add_events( EventSet, Events, NUM_EVENTS) != PAPI_OK){
     printf("PAPI_add_events failed\n");
     }

     for (i=0; i<NUM_EVENTS; i++){
     res_papi[i] = 0;
     }


     if ( PAPI_start( EventSet ) != PAPI_OK){
     printf("PAPI_read_counters failed\n");
     }


 */


  for(i=0;i<(x_N+1)*BlockWidth_t;i++)
    outputbuffer[i]=0.0;



 
  for(i=0;i<(x_N+1)*(t_N+1);i++)
    u[i]=0.0;  

  
  
  for(i=1;i<=x_N-1;i++)    //coefficients alpha	   
    {	             

      powfunction(i)=pow(dt,alpha);
      gammafunction(i)=tgamma(2-alpha);
      
    }

	  
	   
	 
  
   




  for (i=1; i<x_N; i++) 
    for (k=1; k<=t_N; k++)
      {
	coef(k) = 2.0*Jfunction(alpha,k+1)
	  -Jfunction(alpha,k)-Jfunction(alpha,k+2);
      }



  for(i=1;i<=x_N-1;i++)            /*intial value u(x,0)=0*/
    {
      u(i,0)=1.0;
    }
	
  u(0,0)=0.2;  		  /*left corner*/
  u(x_N,0)=0.2;                /*right corner*/
 
  for(i=1;i<=t_N;i++)
    {
      u(0,i)=0.2;                 /*left boundary*/ 
      u(x_N,i)=0.2;               /*right boundary*/	 	
    }
	
 begin=clock();
  for(i=1;i<=x_N-1;i++)             /*compute  u(:,1)*/
    {        
      u(i,1)=(1-2*r)*u(i,0)+r*u(i+1,0)+r*u(i-1,0);//+pow(dt,alpha)*tgamma(2-alpha)*q(alpha,sinfunction(i), gammafunction2(i),t(1)) ;          
               
    }

 
  real sum;

  for (k=1; k<5; k++)              /*compute other value */
    for (i=1; i<x_N; i++) {    
      sum=0.0;
      for (j=1; j<k; j++){
	sum += u(i,j) * coef(k-j);
      }    
      u(i,k+1)= r*(u(i+1,k)+u(i-1,k)) +(s-2*r)*u(i,k) + sum + (Jfunction(alpha,(k+1))- Jfunction(alpha,k))*u(i,0) ;//+ powfunction(i)*gammafunction(i)*q(alpha,sinfunction(i), gammafunction2(i),t(k+1))  ;   
                                                 
    }






    	   
  for (k=5; k<=nBlock_t*BlockWidth_t; k+=4){              /*compute other value */
    for (i=1; i<=nBlock_x*BlockWidth_x; i+=4) {
          
      sum1=_mm256_setzero_pd();
      sum2=_mm256_setzero_pd();
      sum3=_mm256_setzero_pd();
      sum4=_mm256_setzero_pd();
      for (j=1; j<k; j++){
       /*
        load_u1=_mm256_load_pd((real*)(u+i*(t_N+1)+k-j)); 
        load_u2=_mm256_load_pd((real*)(u+(i+1)*(t_N+1)+k-j));
        load_u3=_mm256_load_pd((real*)(u+(i+2)*(t_N+1)+k-j)); 
        load_u4=_mm256_load_pd((real*)(u+ (i+3)*(t_N+1)+k-j)); 
              
       
        coef1= _mm256_broadcast_sd(coef+i*(t_N+1)+j);
        coef2= _mm256_broadcast_sd(coef+(i+1)*(t_N+1)+j);
        coef3= _mm256_broadcast_sd(coef+(i+2)*(t_N+1)+j);
        coef4= _mm256_broadcast_sd(coef+(i+3)*(t_N+1)+j);       
        */

        load_u1=_mm256_load_pd((real*)(coef+k-j)); 
        load_u2=_mm256_load_pd((real*)(coef+k-j));
        load_u3=_mm256_load_pd((real*)(coef+k-j)); 
        load_u4=_mm256_load_pd((real*)(coef+k-j)); 
              
       
        coef1= _mm256_broadcast_sd(u+i*(t_N+1)+j);
        coef2= _mm256_broadcast_sd(u+(i+1)*(t_N+1)+j);
        coef3= _mm256_broadcast_sd(u+(i+2)*(t_N+1)+j);
        coef4= _mm256_broadcast_sd(u+(i+3)*(t_N+1)+j);       

                


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
	   
	  u(i+m,k+1)= r*(u(i+m+1,k)+u(i+m-1,k)) +(1-(  Jfunction(alpha,2)-Jfunction(alpha,1) )-2*r)*u(i+m,k) + (outputbuffer[(i+m)*4]) + (Jfunction(alpha,(k+1))- Jfunction(alpha,k))*u(i+m,0) ;//+ powfunction(i+m)*gammafunction(i+m)*q(alpha,sinfunction(i+m), gammafunction2(i+m),t(k+1))  ;   
      
	                
        
	}       
    
                      
    }
 
    /*compute the remains of u(:,k+1)*/
 
 

 
    for (i=nBlock_x*BlockWidth_x+1; i<x_N; i++) {    
      sum_rem[i]=0.0;
      sum1_rem[i]=0.0;
      sum2_rem[i]=0.0;
      sum3_rem[i]=0.0;
      for (j=1; j<k; j++){
	sum_rem[i] +=u(i,j) *coef(k-j);
        sum1_rem[i]+=u(i,j)*coef(k-j+1);
        sum2_rem[i]+=u(i,j)*coef(k-j+2);    
        sum3_rem[i]+=u(i,j)*coef(k-j+3);     
      }    

   


      u(i,k+1)= r*(u(i+1,k)+u(i-1,k)) +(s-2*r)*u(i,k) + sum_rem[i] + (Jfunction(alpha,(k+1))- Jfunction(alpha,k))*u(i,0) ;//+ powfunction(i)*gammafunction(i)*q(alpha,sinfunction(i), gammafunction2(i),t(k+1))  ;   
                   
            
        
           
    }
 
  
  
    /*compute u(i,k+2)*/

 
    for (i=1; i<=nBlock_x*BlockWidth_x; i++) { 
           
      outputbuffer[i*4+1]+=u(i,k)*coef(1);

      u(i,k+2)=r*(u(i+1,k+1)+u(i-1,k+1)) +(s-2*r)*u(i,k+1) + outputbuffer[i*4+1] + (Jfunction(alpha,(k+2))- Jfunction(alpha,k+1))*u(i,0);// + powfunction(i)*gammafunction(i)*q(alpha,sinfunction(i), gammafunction2(i),t(k+2))  ; 
     
    }
 
  
 
 
    /*compute the remains of u(:,k+2)*/
 
    for (i=nBlock_x*BlockWidth_x+1; i<x_N; i++) {  
        
      sum1_rem[i]+=u(i,k)*coef(1);
      
      u(i,k+2)=r*(u(i+1,k+1)+u(i-1,k+1)) +(s-2*r)*u(i,k+1) + sum1_rem[i] + (Jfunction(alpha,(k+2))- Jfunction(alpha,k+1))*u(i,0) ;//+ powfunction(i)*gammafunction(i)*q(alpha,sinfunction(i), gammafunction2(i),t(k+2))  ; 
      
      
      

    }
 
    
    /*compute u(i,k+3)*/

    for (i=1; i<=nBlock_x*BlockWidth_x; i++) { 
	
      outputbuffer[i*4+2]+=u(i,k)*coef(2)+u(i,k+1)*coef(1);
      u(i,k+3)=r*(u(i+1,k+2)+u(i-1,k+2)) +(s-2*r)*u(i,k+2) +outputbuffer[i*4+2] + (Jfunction(alpha,(k+3))- Jfunction(alpha,k+2))*u(i,0) ;//+ powfunction(i)*gammafunction(i)*q(alpha,sinfunction(i), gammafunction2(i),t(k+3))  ;
        
    
    }
    
    /*compute remains of  u(i,k+3)*/   
 
    for (i=nBlock_x*BlockWidth_x+1; i<x_N; i++) { 
      sum2_rem[i]+=u(i,k)*coef(2)+u(i,k+1)*coef(1);
      u(i,k+3)=r*(u(i+1,k+2)+u(i-1,k+2)) +(s-2*r)*u(i,k+2) + sum2_rem[i] + (Jfunction(alpha,(k+3))- Jfunction(alpha,k+2))*u(i,0);// + powfunction(i)*gammafunction(i)*q(alpha,sinfunction(i), gammafunction2(i),t(k+3))  ;
         
       
    }
    
    
    /*compute u(i,k+4)*/
    
    
    
    for (i=1; i<=nBlock_x*BlockWidth_x; i++) { 
      	
      outputbuffer[i*4+3]+=u(i,k)*coef(3)+u(i,k+1)*coef(2)+u(i,k+2)*coef(1);
      u(i,k+4)=r*(u(i+1,k+3)+u(i-1,k+3)) +(s-2*r)*u(i,k+3) +outputbuffer[i*4+3]+ (Jfunction(alpha,(k+4))- Jfunction(alpha,k+3))*u(i,0);// + powfunction(i)*gammafunction(i)*q(alpha,sinfunction(i), gammafunction2(i),t(k+4))  ;
    
    
    
    }
    
    /*compute remains of u(i,k+4)*/
    
    for (i=nBlock_x*BlockWidth_x+1; i<x_N; i++) { 
      sum3_rem[i]+=u(i,k)*coef(3)+u(i,k+1)*coef(2)+u(i,k+2)*coef(1);
      u(i,k+4)=r*(u(i+1,k+3)+u(i-1,k+3)) +(s-2*r)*u(i,k+3) + sum3_rem[i] + (Jfunction(alpha,(k+4))- Jfunction(alpha,k+3))*u(i,0);// + powfunction(i)*gammafunction(i)*q(alpha,sinfunction(i), gammafunction2(i),t(k+4))  ;
    }

 
  }

  /*compute the remains in time*/


  real sum_t;

  for (k=nBlock_t*BlockWidth_t; k<t_N; k++)              /*compute other value */
    for (i=1; i<x_N; i++) {    
      sum_t=0.0;
      for (j=1; j<k; j++){
	sum_t += u(i,j) * coef(k-j);
      }    
      u(i,k+1)= r*(u(i+1,k)+u(i-1,k)) +(s-2*r)*u(i,k) + sum_t + (Jfunction(alpha,(k+1))- Jfunction(alpha,k))*u(i,0) ;//+ powfunction(i)*gammafunction(i)*q(alpha,sinfunction(i), gammafunction2(i),t(k+1))  ;   
                                         
    }

                  
  finish=clock();
/*
if ( PAPI_stop( EventSet, res_papi ) != PAPI_OK){
    printf("PAPI_accum_counters failed\n");
    }
  */
  
  real time=1.0*(finish-begin)/CLOCKS_PER_SEC;
     	
     		
  real numer=0.0;
     	
                    
  for(i=1;i<=x_N-1;i++)    
    for(k=1;k<=t_N;k++)     
      {  
	numer +=u(i,k)*u(i,k);
      }
     	
  printf("AVX:x_N:  %5d, t_N:  %5d   ,Norm_n:  %e, timeusage is %6.2f seconds\n\n",x_N,t_N,sqrt(numer), time);
 /*
 for (i = 0; i<NUM_EVENTS; i++){
    PAPI_event_code_to_name(Events[i], EventName);
    printf("PAPI Event name: %s, value: %lld\n", EventName, res_papi[i]);
    }

    PAPI_shutdown ();
  
*/


  _mm_free(u);    
  
  _mm_free(coef);
  _mm_free(powfunction); 
  _mm_free(gammafunction ); 
  _mm_free(outputbuffer); 
  free(sum_rem);  
  free(sum1_rem);  
  free(sum2_rem);
  free(sum3_rem);
 
 
  return EXIT_SUCCESS;     
}

