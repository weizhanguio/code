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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<time.h>  
#include<omp.h>
#include<immintrin.h>
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
  int m;	              	
  real begin;
  real finish;
  int  x_N,y_N,z_N, t_N;
  x_N=20;                   //number of space intervals
  y_N=20;                    
  z_N=20;                    
  /*we choose x_N=y_N=z_N*/
  t_N=20;                   //number of time intervals
  

  
  if(nargs>1){
    x_N=atof(args[1]);
    y_N=atof(args[1]);
    z_N=atof(args[1]);
    t_N=atof(args[2]);

  }
    
  int threadNum;	
#pragma omp parallel
  {
#pragma omp master
    {
      threadNum=omp_get_num_threads();
      printf("Number of threads = %d\n",threadNum);
    }
  }

  int BlockWidth_x=4; 

  int BlockWidth_t=4; 

  int nBlock_x=(x_N-1)/BlockWidth_x;    //x=0  x=x_N are known, x_N-1 points left

  int nBlock_t=(t_N-1)/BlockWidth_t;  // t=0  t=1 are computed separately, t_N-1 points left

  int rem_x=(x_N-1)%BlockWidth_x;

  int rem_t=(t_N-1)%BlockWidth_t;


  int NextStart=x_N-1+4-(x_N-1)%4+1;    // padding, add one block 

  int x_Nnew=NextStart+3;       //next i start position: x_N-1+4-(x_N-1)%4+1
  int NextStart_t=t_N-1+4-(t_N-1)%4+1;    // padding, add one block 

  int t_Nnew=NextStart_t+3;       //next i start position: x_N-1+4-(x_N-1)%4+1
  
   
 
  real   dh;        dh= ((x_end-x_start)/x_N);          // space step
  /*we choose dx=dy=dz=dh*/
  real   dt;        dt= ((t_end-t_start)/t_N);          // time step
  
  
  printf("---------3D omp version----------\n");
  printf("%d intervals in space(X,Y,Z), %d intervals in time \n\n",x_N, t_N);  
  /*x_N+1 points in space, t_N+1 points in time*/
	
  

  real A=pow(dt,-alpha)/tgamma(2-alpha);
  real B=K/A/dh/dh;

  real *u=(real*)_mm_malloc(1u*(x_Nnew+1 )*(y_N+1 )*(z_N+1 )*(t_Nnew)*sizeof(real),64);   
  real *coeff=(real*)_mm_malloc(t_Nnew*sizeof(real),64);
 
  real *sum=(real*)_mm_malloc(1ul* (x_Nnew+1)*(y_N+1 )*(z_N+1)*sizeof(real),64);
  real *sum1=(real*)_mm_malloc(1ul*(x_Nnew+1)*(y_N+1 )*(z_N+1)*sizeof(real),64);
  real *sum2=(real*)_mm_malloc(1ul*(x_Nnew+1 )*(y_N+1 )*(z_N+1)*sizeof(real),64);
  real *sum3=(real*)_mm_malloc(1ul*(x_Nnew+1)*(y_N+1 )*(z_N+1)*sizeof(real),64);

  real* outputbuffer=(real*)_mm_malloc(16*threadNum*64*sizeof(real),64);
  real s=(pow(2,1-alpha)- pow(1,1-alpha));  
 
#define  u(x,y,z,t)   u[1ul*(x)*(y_N+1)*(z_N+1)*(t_Nnew)+1ul* (y)*(z_N+1)*(t_Nnew)+(z)*(t_Nnew) +(t) ]
#define  coeff(t)     coeff[(t)] 

#define  sum(x,y,z)     sum[1ul*(x)*(y_N+1)*(z_N+1)+(y)*(z_N+1)+(z) ]
#define  sum1(x,y,z)    sum1[1ul*(x)*(y_N+1)*(z_N+1)+(y)*(z_N+1)+(z) ] 
#define  sum2(x,y,z)    sum2[1ul*(x)*(y_N+1)*(z_N+1)+(y)*(z_N+1)+(z) ]
#define  sum3(x,y,z)    sum3[1ul*(x)*(y_N+1)*(z_N+1)+(y)*(z_N+1)+(z) ]


  begin=omp_get_wtime();

#pragma omp parallel default(shared) private(i,n)
  {
#pragma omp for	
   // for(i=1;i<=x_Nnew;i++)
      for(n=0;n<=t_N;n++)
	coeff(n)=pow(n+2,1-alpha)+pow(n,1-alpha)-2.0*pow(n+1,1-alpha);
  } 

 // for(i=1;i<=x_Nnew;i++)
    for(n=t_N+1;n<=t_Nnew;n++)
      coeff(n)=0.0;			

#pragma omp parallel for
  for(i=0;i<(x_Nnew+1 )*(y_N+1 )*(z_N+1 )*t_Nnew;i++)
    u[i]=0.0;     



  /*boundary value*/

#pragma omp parallel default(shared) private(j,k,p)
  {
#pragma omp for	  
    for(j=0;j<=y_N;j++)             
      for(k=0;k<=z_N;k++)
	for(p=0;p<=t_N;p++)
	  {
	    u(0,j,k,p)=0.2;
	    u(x_N,j,k,p)=0.2;
	  }
  }



#pragma omp parallel default(shared) private(i,k,p)
  {
#pragma omp for		
    for(i=0;i<=x_N;i++)            
      for(k=0;k<=z_N;k++)
	for(p=0;p<=t_N;p++)
	  {
	    u(i,0,k,p)=0.2;
	    u(i,y_N,k,p)=0.2;
	  }
  }





#pragma omp parallel default(shared) private(i,j,p)
  {
#pragma omp for
    for(i=0;i<=x_N;i++)             
      for(j=0;j<=y_N;j++)
	for(p=0;p<=t_N;p++)
	  {
	    u(i,j,0,p)=0.2;
	    u(i,j,z_N,p)=0.2;
	  }
  }




  /*initial value*/

#pragma omp parallel default(shared) private(i,j,k)
  {
#pragma omp for
    for(i=1;i<=x_N-1;i++)             
      for(j=1;j<=y_N-1;j++)
	for(k=1;k<=z_N-1;k++)
	  {
	    u(i,j,k,0)=1.0;
	 
	  }

  }
  /*compute u(:,:,:,1)*/
#pragma omp parallel default(shared) private(i,j,k)
  {  
#pragma omp for          
    for (i=1; i<=x_N-1; i++)     
      for (j=1; j<=y_N-1; j++)
	for(k=1;k<=z_N-1;k++)
	  {
                
	    u(i,j,k,1)= (1-6.0*B)*u(i,j,k,0)+ B*(u(i+1,j,k,0)+u(i-1,j,k,0)+u(i,j+1,k,0)+u(i,j-1,k,0)+u(i,j,k+1,0)+u(i,j,k-1,0));//+q(i*dh,j*dh,k*dh,1.0*dt)/A; 
	  }
  } 
   
    
#pragma omp barrier
   
  real sum_t;
   
  for (n=1; n<5; n++)              
    for (i=1; i<=x_N-1; i++)     
      for (j=1; j<=y_N-1; j++)
        for(k=1;k<=z_N-1;k++)
	  {
	    sum_t=0.0;
	    for(p=1; p<=n-1; p++)
	      {
		sum_t +=coeff(p)*u(i,j,k,n-p);
	      }
          
	    u(i,j,k,n+1)= (1-6.0*B)*u(i,j,k,n)-sum_t +u(i,j,k,0)*(pow(n+1,1-alpha)-pow(n,1-alpha))-u(i,j,k,n)*s+B*(u(i+1,j,k,n)+u(i-1,j,k,n)+u(i,j+1,k,n)+u(i,j-1,k,n)+u(i,j,k+1,n)+u(i,j,k-1,n));//+q(i*dh,j*dh,k*dh,(n+1)*dt)/A;          

	  }
   
   
   
  #if 1 
   
  /*compute other u*/ 
    	 
#pragma omp parallel default(shared) private(n,p,i,j,k)
  {
    for (n=5; n<=NextStart_t-1; n+=4)  
      {        
        int id=omp_get_thread_num();
        
#pragma omp for
	for (i=1; i<=NextStart-1; i+=4)
	  for (j=1; j<=y_N-1; j++)
	    for(k=1;k<=z_N-1;k++)
	      {   
               
		__m256d  s1;
		__m256d  s2;
		__m256d  s3;
		__m256d  s4;

		__m256d  load_u1;
		__m256d  load_u2;
		__m256d  load_u3;
		__m256d  load_u4;

		__m256d  mulp1;
		__m256d  mulp2;
		__m256d  mulp3;
		__m256d  mulp4;

		__m256d  load_coeff1;
                __m256d  load_coeff2;
                __m256d  load_coeff3;
                __m256d  load_coeff4;



		s1=_mm256_setzero_pd();
		s2=_mm256_setzero_pd();
		s3=_mm256_setzero_pd();
		s4=_mm256_setzero_pd();
  
	      
		for(p=1; p<=n-1; p++)
		  {
		
		    load_coeff1=_mm256_load_pd((real*)(coeff+n-p)); 
                    load_coeff2=_mm256_load_pd((real*)(coeff+n-p));
		    load_coeff3=_mm256_load_pd((real*)(coeff+n-p)); 
                    load_coeff4=_mm256_load_pd((real*)(coeff+n-p));
 		 
		    load_u1= _mm256_broadcast_sd(&u(i,j,k,p));

		    load_u2= _mm256_broadcast_sd(&u(i+1,j,k,p));

		    load_u3= _mm256_broadcast_sd(&u(i+2,j,k,p));

		    load_u4= _mm256_broadcast_sd(&u(i+3,j,k,p));		 
			       

		    mulp1=_mm256_mul_pd(load_u1 ,load_coeff1);

		    mulp2=_mm256_mul_pd(load_u2 ,load_coeff2);

		    mulp3=_mm256_mul_pd(load_u3 ,load_coeff3);

		    mulp4=_mm256_mul_pd(load_u4 ,load_coeff4);	   
          
		    s1=_mm256_add_pd(s1,mulp1);  /*i*/

		    s2=_mm256_add_pd(s2,mulp2);  /*i+1*/

		    s3=_mm256_add_pd(s3,mulp3); /*i+2*/

		    s4=_mm256_add_pd(s4,mulp4);   /*i+3*/
  
     

		  }

		_mm256_store_pd(&outputbuffer[0+id*16*64],s1);
		_mm256_store_pd(&outputbuffer[4+id*16*64],s2);
		_mm256_store_pd(&outputbuffer[8+id*16*64],s3);
		_mm256_store_pd(&outputbuffer[12+id*16*64],s4);



		sum(i,j,k)=outputbuffer[0+id*16*64];
		sum1(i,j,k)=outputbuffer[1+id*16*64];
		sum2(i,j,k)=outputbuffer[2+id*16*64];
		sum3(i,j,k)=outputbuffer[3+id*16*64]; 
    
		sum(i+1,j,k)=outputbuffer[4+id*16*64];
		sum1(i+1,j,k)=outputbuffer[5+id*16*64];
		sum2(i+1,j,k)=outputbuffer[6+id*16*64];
		sum3(i+1,j,k)=outputbuffer[7+id*16*64]; 

		sum(i+2,j,k)=outputbuffer[8+id*16*64];
		sum1(i+2,j,k)=outputbuffer[9+id*16*64];
		sum2(i+2,j,k)=outputbuffer[10+id*16*64];
		sum3(i+2,j,k)=outputbuffer[11+id*16*64];
          
		sum(i+3,j,k)=outputbuffer[12+id*16*64];
		sum1(i+3,j,k)=outputbuffer[13+id*16*64];
		sum2(i+3,j,k)=outputbuffer[14+id*16*64];
		sum3(i+3,j,k)=outputbuffer[15+id*16*64];
          
  
		for(m=0;m<4;m++)
		  {
		    u(i+m,j,k,n+1)= (1-6.0*B)*u(i+m,j,k,n)-sum(i+m,j,k) +u(i+m,j,k,0)*(pow(n+1,1-alpha)-pow(n,1-alpha))-u(i+m,j,k,n)*(pow(2,1-alpha)- pow(1,1-alpha))+B*(u(i+m+1,j,k,n)+u(i+m-1,j,k,n)+u(i+m,j+1,k,n)+u(i+m,j-1,k,n)+u(i+m,j,k+1,n)+u(i+m,j,k-1,n));//+q((i+m)*dh,j*dh,k*dh,(n+1)*dt)/A;               
		    
            
		  }

	      }


#pragma omp barrier


#pragma omp for
	for(j=1; j<=y_N-1; j++)
	  for(k=1;k<=z_N-1;k++)
	    u(x_N,j,k,n+1)=0.2;
          



#pragma omp for
	for (i=1; i<=x_N-1; i++)     
	  for (j=1; j<=y_N-1; j++)
	    for(k=1;k<=z_N-1;k++)
	      {
		sum1(i,j,k)+=coeff(1)*u(i,j,k,n);    
		u(i,j,k,n+2)=(1-6.0*B)*u(i,j,k,n+1)-sum1(i,j,k) +u(i,j,k,0)*(pow(n+2,1-alpha)-pow(n+1,1-alpha))-u(i,j,k,n+1)*s+B*(u(i+1,j,k,n+1)+u(i-1,j,k,n+1)+u(i,j+1,k,n+1)+u(i,j-1,k,n+1)+u(i,j,k+1,n+1)+u(i,j,k-1,n+1));//+q(i*dh,j*dh,k*dh,(n+2)*dt)/A;             
	      }

#pragma omp barrier
#pragma omp for
	for (i=1; i<=x_N-1; i++)     
	  for (j=1; j<=y_N-1; j++)
	    for(k=1;k<=z_N-1;k++)
	      {
		sum2(i,j,k)+=coeff(2)*u(i,j,k,n)+coeff(1)*u(i,j,k,n+1);    
		u(i,j,k,n+3)=(1-6.0*B)*u(i,j,k,n+2)-sum2(i,j,k) +u(i,j,k,0)*(pow(n+3,1-alpha)-pow(n+2,1-alpha))-u(i,j,k,n+2)*s+B*(u(i+1,j,k,n+2)+u(i-1,j,k,n+2)+u(i,j+1,k,n+2)+u(i,j-1,k,n+2)+u(i,j,k+1,n+2)+u(i,j,k-1,n+2));//+q(i*dh,j*dh,k*dh,(n+3)*dt)/A;           
                  
	      }
 
#pragma omp barrier
#pragma omp for  
	for (i=1; i<=x_N-1; i++)     
	  for (j=1; j<=y_N-1; j++)
	    for(k=1;k<=z_N-1;k++)
	      {
		sum3(i,j,k)+=coeff(3)*u(i,j,k,n)+coeff(2)*u(i,j,k,n+1)+coeff(1)*u(i,j,k,n+2);    
		u(i,j,k,n+4)=(1-6.0*B)*u(i,j,k,n+3)-sum3(i,j,k) +u(i,j,k,0)*(pow(n+4,1-alpha)-pow(n+3,1-alpha))-u(i,j,k,n+3)*s+B*(u(i+1,j,k,n+3)+u(i-1,j,k,n+3)+u(i,j+1,k,n+3)+u(i,j-1,k,n+3)+u(i,j,k+1,n+3)+u(i,j,k-1,n+3));//+q(i*dh,j*dh,k*dh,(n+4)*dt)/A;             
                
	      }




      }

 
  }

#endif
  finish=omp_get_wtime();
  
  real time=1.0*(finish-begin);
     	
     
  		
  real numer=0.0;
  
                 
  for(i=1;i<=x_N-1;i++) 
    for(j=1;j<=y_N-1;j++)
      for(k=1;k<=z_N-1;k++)     
	for(n=1;n<=t_N;n++) 
	  { 
           
	    numer +=u(i,j,k,n)*u(i,j,k,n);
	  }
     	
  printf("x_N: %d, t_N: %d ,timeusage: %6.2f seconds.Threads:%d,  Norm_n:  %e\n\n",x_N, t_N,time,threadNum, sqrt(numer));

  _mm_free(u);    
  _mm_free(coeff);
  _mm_free(sum);
  _mm_free(sum1);
  _mm_free(sum2);
  _mm_free(sum3);
  _mm_free(outputbuffer);
 
 
 
  return EXIT_SUCCESS;     
}



