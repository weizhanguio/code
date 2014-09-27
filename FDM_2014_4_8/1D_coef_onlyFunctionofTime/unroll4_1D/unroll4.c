/*
  purpuse: finite difference scheme for variable-order time fractional diffusion equation, fractional order is a function of  time t.
  reference paper:Sun Hongguang  "finite difference scheme for variable-order time fractional diffusion equation", page 1250085-12, example 1.
  date: 2014.3.5
  author: Zhang Wei
*/
//#include<papi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
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
#define  powfunction(x)        powfunction[x]
#define  gammafunction(x)      gammafunction[x]
#define  Jfunction(alpha,k)    pow((k),1.0-(alpha))   /*compute the coef*/
#define  real           double 




 
int 
main(int nargs, char** args)
{
	
	
  int i;
  int j;
  int k;
	              	
  clock_t begin;
  clock_t finish;
  int  x_N, t_N;
  x_N=3000;                   //number of space intervals
  t_N=3000;                   //number of time intervals

  /*
    int Events[] = { PAPI_L2_DCR, PAPI_L3_DCR };
    int NUM_EVENTS = sizeof(Events)/sizeof(Events[0]);
    long long res_papi[NUM_EVENTS];
    char EventName[128];
    int num_hwcntrs = 0;
    int EventSet = PAPI_NULL;
    int retval;
  */

  if(nargs>1){
    x_N=atof(args[1]);
    t_N=atof(args[2]);

  }
    
 
  double   dx;        dx= ((x_end-x_start)/x_N);          // space step
  double   dt;        dt= ((t_end-t_start)/t_N);          // time step
  
  printf("---------unroll 4---coef(x,t)----------\n");
  printf("%d intervals in space, %d intervals in time \n\n",x_N, t_N);  
  /*x_N+1 points in space, t_N+1 points in time*/
	

#define  u(x,y)         u[1ul*(x)*(t_N+1)+(y)]
#define  coef(y)      coef[(y)]




  real *u=(real*)malloc(1ul*(x_N+1 )*(t_N+1)*sizeof(real));   
  
  real *coef=(real*)malloc((t_N+1)*sizeof(real));  /*coef is d in eq.6  */ 
  real *powfunction=(real*)malloc((x_N+1)*sizeof(real));
  real *gammafunction=(real*)malloc((x_N+1)*sizeof(real));
   
  real *sum1=(real*)malloc((x_N+1)*sizeof(real));
  real *sum2=(real*)malloc((x_N+1)*sizeof(real));
  real *sum3=(real*)malloc((x_N+1)*sizeof(real)); 
  real *sum=(real*)malloc((x_N+1)*sizeof(real));
  real *p_u; 
  real r=1.0*K*pow(dt,alpha)*tgamma(2-alpha)/(dx*dx);
  real s=1-(pow(2,1.0-(alpha))-pow(1,1.0-(alpha)));   //store the value of 1-(  Jfunction(alpha,2)-Jfunction(alpha,1) )	
 
  
  for(i=1;i<=x_N-1;i++)    
    {	             
       
      powfunction(i)=pow(dt,alpha);
      gammafunction(i)=tgamma(2-alpha);
      
    }
	  
	   
	 
  real Min=0.0;
	
  for(i=1;i<=x_N-1;i++)
	 
    {     
      real s1;
      s1= r-0.25*(2- pow(2,1-alpha)) ;     
      Min=(s1<Min?s1:Min);
    }

  if(Min>=0)
    {
      printf("error, do not satisfy the stability requirement.\n\n");
      
    }
   




 
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
      u(i,1)=(1-2*r)*u(i,0)+r*u(i+1,0)+r*u(i-1,0);
    }

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
  
  real  sum_t=0.0;

  for (k=1; k<5; k++)              
    for (i=1; i<x_N; i++) {    
      sum_t=0.0;
      p_u=&(u[i*(t_N+1)]);
      
      for (j=1; j<k; j++){
	
	//sum_t += p_u[j]*coef[k-j];
        sum_t+= u(i,j) * coef(k-j);
      }    
      u(i,k+1)= r*(u(i+1,k)+u(i-1,k)) +(s-2*r)*u(i,k) + sum_t + (Jfunction(alpha,(k+1))- Jfunction(alpha,k))*u(i,0);
                                                 
    }

    	   
  for (k=5; k<t_N; k+=4){              /*compute other value */
    for (i=1; i<x_N; i++) {    
      sum[i]=0.0;
      sum1[i]=0.0;
      sum2[i]=0.0;
      sum3[i]=0.0;
   //   p_u=&(u[i*(t_N+1)]);
      
      for (j=1; j<k; j++){
	/*
        sum[i] +=p_u[j]*coef[k-j];
        sum1[i]+=p_u[j]*coef[k-j+1];
        sum2[i]+=p_u[j]*coef[k-j+2];     
        sum3[i]+=p_u[j]*coef[k-j+3]; 
*/
    
        sum[i] +=u(i,j)*coef(k-j);
        sum1[i]+=u(i,j)*coef(k-j+1);
        sum2[i]+=u(i,j)*coef(k-j+2);     
        sum3[i]+=u(i,j)*coef(k-j+3);
  

      }    

    

      u(i,k+1)= r*(u(i+1,k)+u(i-1,k)) +(s-2*r)*u(i,k) + sum[i] + (Jfunction(alpha,(k+1))- Jfunction(alpha,k))*u(i,0) ;  
                   
                       
    }
  
    

    for (i=1; i<x_N; i++) {     
      sum1[i]+=u(i,k)*coef(1);
      u(i,k+2)=r*(u(i+1,k+1)+u(i-1,k+1)) +(s-2*r)*u(i,k+1) + sum1[i] + (Jfunction(alpha,(k+2))- Jfunction(alpha,k+1))*u(i,0); 
            
 
    }


    for (i=1; i<x_N; i++) { 
      sum2[i]+=u(i,k+1)*coef(1)+u(i,k)*coef(2);
      u(i,k+3)=r*(u(i+1,k+2)+u(i-1,k+2)) +(s-2*r)*u(i,k+2) + sum2[i] + (Jfunction(alpha,(k+3))- Jfunction(alpha,k+2))*u(i,0);
    }
   
   
    for (i=1; i<x_N; i++) { 
      sum3[i]+=u(i,k+2)*coef(1)+u(i,k+1)*coef(2)+u(i,k)*coef(3);
      u(i,k+4)=r*(u(i+1,k+3)+u(i-1,k+3)) +(s-2*r)*u(i,k+3) + sum3[i] + (Jfunction(alpha,(k+4))- Jfunction(alpha,k+3))*u(i,0);
    }




  }

                   
  finish=clock();
 
  
  
    double time=1.0*(finish-begin)/CLOCKS_PER_SEC;
  /*
    if ( PAPI_stop( EventSet, res_papi ) != PAPI_OK){
    printf("PAPI_accum_counters failed\n");
    }
  */	
     		
  real numer=0.0;
                    
  for(i=1;i<=x_N-1;i++)    
    for(k=1;k<=t_N;k++)     
      {  
	numer +=u(i,k)*u(i,k);
      }
     	
  printf("unroll 4: x_N:  %5d, t_N:  %5d   ,Norm_n:  %e, timeusage is %6.2f seconds\n\n",x_N,t_N,sqrt(numer), time);  

  /*
    for (i = 0; i<NUM_EVENTS; i++){
    PAPI_event_code_to_name(Events[i], EventName);
    printf("PAPI Event name: %s, value: %lld\n", EventName, res_papi[i]);
    }

    PAPI_shutdown ();
  */


  free(u);    
      
  free(coef);
  free(powfunction); 
  free(gammafunction );
 
  free(sum);  
  free(sum1);  
  free(sum2);
  free(sum3);
  

 
  return EXIT_SUCCESS;     
}

