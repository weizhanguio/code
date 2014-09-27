/*
  purpuse: finite difference scheme for  3D time fractional diffusion equation
  D_t^alpha(u(x,y,z,t))=K(d^2u/dx^2+d^2u/dy^2+d^2u/dz^2)
  boundary codition:
  u(0,:,:,:)=u(L,:,:,:)=u(:,0,:,:)=u(:,L,:,:)=u(:,:,0,:)=u(:,:,L,:)=0.2  
  intial condition:
  u(x,y,z,0)=1  
  date: 2014.3.5
  author: Zhang Wei
*/
//#include<papi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<time.h>  
#define  K             0.00001                   /* diffusion efficient*/
#define  x_start       0.0                    /*starting point of space*/
#define  x_end         10.0                      /*end point of space*/
#define  t_start       0.0                    /*starting point of time*/
#define  t_end         1.0     	       /*end point of time*/
#define  L           (x_end-x_start)
#define  alpha        0.9
#define  real         double





 
int 
main(int nargs, char** args)
{
	
	
  int i;
  int j;
  int k;
  int p;
  int n;	              	
  clock_t begin;
  clock_t finish;
  int  x_N,y_N,z_N, t_N;
  x_N=50;                   //number of space intervals
  y_N=50;                    
  z_N=50;                    
  /*we choose x_N=y_N=z_N*/
  t_N=200;                   //number of time intervals

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
    y_N=atof(args[1]);
    z_N=atof(args[1]);
    t_N=atof(args[2]);

  }
    
 
  real   dh;        dh= ((x_end-x_start)/x_N);          // space step
  /*we choose dx=dy=dz=dh*/
  real   dt;        dt= ((t_end-t_start)/t_N);          // time step
  
  
  printf("---------3D baseline version----------\n");
  printf("%d intervals in space(X,Y,Z), %d intervals in time \n\n",x_N, t_N);  
  /*x_N+1 points in space, t_N+1 points in time*/
	
  

  real A=pow(dt,-alpha)/tgamma(2-alpha);
  real B=K/A/dh/dh;

  real *u=(real*)calloc(1ul*(x_N+1 )*(y_N+1 )*(z_N+1 )*(t_N+1),sizeof(real));   
  real *coeff=(real*)calloc((t_N+1),sizeof(real));
   real *p_u;
   

  real sum;

#define  u(x,y,z,t)   u[1ul*(x)*(y_N+1)*(z_N+1)*(t_N+1)+1ul*(y)*(z_N+1)*(t_N+1)+1ul* (z)*(t_N+1) +(t) ]
#define  coeff(x)     coeff[x] 
  
   real s=(pow(2,1-alpha)- pow(1,1-alpha));  


  for(i=0;i<=t_N;i++)
    coeff(i)=pow(i+2,1-alpha)+pow(i,1-alpha)-2.0*pow(i+1,1-alpha);
			
   
  /*boundary value*/	  
  for(j=0;j<=y_N;j++)             
    for(k=0;k<=z_N;k++)
      for(p=0;p<=t_N;p++)
	{
	  u(0,j,k,p)=0.2;
	  u(x_N,j,k,p)=0.2;
	}
	
  for(i=0;i<=x_N;i++)            
    for(k=0;k<=z_N;k++)
      for(p=0;p<=t_N;p++)
	{
	  u(i,0,k,p)=0.2;
	  u(i,y_N,k,p)=0.2;
	}

  for(i=0;i<=x_N;i++)             
    for(j=0;j<=y_N;j++)
      for(p=0;p<=t_N;p++)
	{
	  u(i,j,0,p)=0.2;
	  u(i,j,z_N,p)=0.2;
	}


  /*initial value*/
  for(i=1;i<=x_N-1;i++)             
    for(j=1;j<=y_N-1;j++)
      for(k=1;k<=z_N-1;k++)
	{
	  u(i,j,k,0)=1.0;
	 
	}


  /*compute u(:,:,:,1)*/


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


  begin=clock();     
       
  for (i=1; i<=x_N-1; i++)     
    for (j=1; j<=y_N-1; j++)
      for(k=1;k<=z_N-1;k++)
	{
                
	  u(i,j,k,1)= (1-6.0*B)*u(i,j,k,0)+ B*(u(i+1,j,k,0)+u(i-1,j,k,0)+u(i,j+1,k,0)+u(i,j-1,k,0)+u(i,j,k+1,0)+u(i,j,k-1,0));
	}



  /*compute other u*/ 
    
//#define  u(x,y,z,t)   u[1ul*(x)*(y_N+1)*(z_N+1)*(t_N+1)+1ul*(y)*(z_N+1)*(t_N+1)+1ul* (z)*(t_N+1) +(t) ]	   
  for (n=1; n<=t_N-1; n++)              
    for (i=1; i<=x_N-1; i++)     
      for (j=1; j<=y_N-1; j++)
        for(k=1;k<=z_N-1;k++)
	  {
	    sum=0.0;
           // p_u=&(u[1ul*(i)*(y_N+1)*(z_N+1)*(t_N+1)+1ul*(j)*(z_N+1)*(t_N+1)+1ul* (k)*(t_N+1)]);
	    for(p=1; p<=n-1; p++)
	      {
            
              //    sum +=p_u[p]*coeff(n-p);
                 sum +=u(i,j,k,p)*coeff(n-p);
	      }
          
	    u(i,j,k,n+1)= (1-6.0*B)*u(i,j,k,n)-sum +u(i,j,k,0)*(pow(n+1,1-alpha)-pow(n,1-alpha))-u(i,j,k,n)*s+B*(u(i+1,j,k,n)+u(i-1,j,k,n)+u(i,j+1,k,n)+u(i,j-1,k,n)+u(i,j,k+1,n)+u(i,j,k-1,n));

	  }
  
           
  finish=clock();
 
  real time=1.0*(finish-begin)/CLOCKS_PER_SEC;

/*
  if ( PAPI_stop( EventSet, res_papi ) != PAPI_OK){
		printf("PAPI_accum_counters failed\n");
	}
 */    	
     
   		
  real numer=0.0;
     	
                    
  for(i=1;i<=x_N-1;i++) 
    for(j=1;j<=y_N-1;j++)
      for(k=1;k<=z_N-1;k++)     
	for(n=1;n<=t_N;n++) 
	  { 
	    numer +=u(i,j,k,n)*u(i,j,k,n);
	  }
     	
   printf("x_N:  %d, t_N:  %d ,time %6.2f seconds. Norm_num  %e \n\n",x_N, t_N,time, sqrt(numer));  
 
/*  
   for (i = 0; i<NUM_EVENTS; i++){
		PAPI_event_code_to_name(Events[i], EventName);
		printf("PAPI Event name: %s, value: %lld\n", EventName, res_papi[i]);
	}

	PAPI_shutdown ();  
*/ 

  free(u);    
  free(coeff);
  
  
 
  return EXIT_SUCCESS;     
}


