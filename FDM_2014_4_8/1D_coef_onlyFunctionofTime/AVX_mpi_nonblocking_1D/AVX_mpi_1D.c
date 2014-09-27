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
#include<mpi.h> 
#define  K             0.000001                   /* diffusion efficient*/
#define  pi            3.1415926
#define  x_start       0.0                    /*starting point of space*/
#define  x_end         10.0                 /*end point of space*/
#define  t_start       0.0                    /*starting point of time*/
#define  t_end         1.0   	              /*end point of time*/
#define  L           (x_end-x_start)
#define  t(x)          (t_start+dt*(x))
#define  x(a)          (x_start+dx*(a))

#define  alpha        0.8
#define  powfunction(x)        powfunction[(x)]
#define  gammafunction(x)      gammafunction[(x)]

#define  Jfunction(alpha,k)    pow((k),1.0-(alpha))   /*compute the coef*/
#define  MPI_type           MPI_DOUBLE
#define  real               double 


int taskid;                                                 
int numtasks;   
int left, right;
int npoints;
int first;
int RtoL=10;
int LtoR=20; 
int  x_N, t_N; 
real *u ;


 
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
  int NextStart=x_N-1+4-(x_N-1)%4+1;
  
  int num_p=(NextStart-1)/4+1;
  
  int x_Nnew=NextStart+3;       //next i start position: x_N-1+4-(x_N-1)%4+1
  
  
  MPI_Init(&nargs, &args);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  MPI_Status stats[4];
  MPI_Request reqs[4];  
  
    
 
  real   dx;        dx= ((x_end-x_start)/x_N);          // space step
  real   dt;        dt= ((t_end-t_start)/t_N);          // time step
  
   
  
  if(taskid==numtasks-1)
    right=0;
  else
    right=taskid+1;


  if(taskid==0)
    left=numtasks-1;
  else
    left=taskid-1;


  if(taskid==0)
    {
      printf("Number of MPI procs=%d\n",numtasks);
      printf("---------MPI non-blocking AVX Version----------\n");
      printf("%d intervals in space, %d intervals in time \n\n",x_N, t_N);
      /*x_N+1 points in space, t_N+1 points in time*/
    }


  int nmin,nleft;     /*nmin: min points for each proc, nleft: the remainder*/
  int npts;

  nmin=num_p/numtasks;                    /*x_N+1 points, 0 and x_N are boundary*/
  nleft=num_p%numtasks;

  for(i=0,k=0; i<numtasks;i++)
    {

      npts=(i<nleft)?nmin+1:nmin;
      if(taskid==i)
	{
	  first=k+1;   /*first index of u computed by  i  processor*/
	  npoints=npts;   /*npoints: number of points computed by the proc*/
	}
      else
	k+=npts;

    }

 
  first=(first-1)*4+1;
 
  
#define  u(x,y)         u[1ul*(x)*(t_N+1)+(y)]
#define  coef(y)      coef[(y)]



  
  u=(real*)_mm_malloc(1ul*(4*npoints+2)*(t_N+1)*sizeof(real), 64);   

  
  real *coef=(real*)_mm_malloc((t_N+1)*sizeof(real),64);  /*coef is d in eq.6  */ 
  real *powfunction=(real*)_mm_malloc((4*npoints+2)*sizeof(real),64);
  real *gammafunction=(real*)_mm_malloc((4*npoints+2)*sizeof(real),64);
  real r=1.0*K*pow(dt,alpha)*tgamma(2-alpha)/(dx*dx);
  real s=1-(pow(2,1.0-(alpha))-pow(1,1.0-(alpha)));   //store the value of 1-(  Jfunction(alpha,2)-Jfunction(alpha,1) ) 
  real* outputbuffer=(real*)_mm_malloc((x_N+1)*BlockWidth_t*sizeof(real),64);
   

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
          

  
  
  for(i=0;i<=(4*npoints+2)*(t_N+1)-1;i++)
    u[i]=0.0;  

 
 
  for(i=1;i<=npoints*4;i++)    //coefficients alpha	   
    {	             
      powfunction(i)=pow(dt,alpha);       
      gammafunction(i)=tgamma(2-alpha);    
    }

 

  int lastpoint;        
  if(taskid!=numtasks-1)
    {lastpoint=npoints*4 ;
    }
  else
    {
      lastpoint=x_N-1-first+1;
    }

 

    for (k=1; k<=t_N; k++)
      {
	coef(k) = 2.0*Jfunction(alpha,k+1)
	  -Jfunction(alpha,k)-Jfunction(alpha,k+2);
	  
      }


 

 
 
  for(i=1;i<= lastpoint;i++)            /*intial value u(x,0)=0*/
    {
      u(i,0)=1.0;
    }
	


  if(taskid==0)
    {
      u(0,0)=0.2;  		   /*left corner*/
      for(i=1;i<=t_N;i++) 
	u(0,i)=0.2;                /*left boundary*/
 
    }


  if(taskid==numtasks-1)
    {
      u(lastpoint+1,0)=0.2;                   /*right corner*/
      for(i=1;i<=t_N;i++)              
	u(lastpoint+1,i)=0.2;                /*right boundary*/
    }


  MPI_Status  status;
  /*exchange data with left neighbour*/
  if(taskid !=0)
    {
      MPI_Recv(&u(0,0),1,MPI_type,left,LtoR,MPI_COMM_WORLD,&status);
      MPI_Send(&u(1,0),1,MPI_type,left,RtoL,MPI_COMM_WORLD);
     
    }
  /*exchange data with right neighbour*/
  if(taskid !=(numtasks-1))
    {
      MPI_Send(&u(lastpoint,0),1,MPI_type,right,LtoR,MPI_COMM_WORLD);
      MPI_Recv(&u(lastpoint+1,0),1,MPI_type,right,RtoL,MPI_COMM_WORLD,&status);

    }


   
begin=MPI_Wtime();
  
  for(i=1;i<=lastpoint;i++)             /*compute  u(:,1)*/
    {
      u(i,1)=(1-2*r)*u(i,0)+r*u(i+1,0)+r*u(i-1,0);

    }
 
 
 
  if(taskid==(numtasks-1))  
    u(x_N-1-first+2,k+1)=0.2;

 
 
 
 
 
  real sum_t;

  for (k=1; k<5; k++)              /*compute other value */
    {
     
      if(taskid !=0)
        {
          MPI_Send(&u(1,k+0),1,MPI_type,left,RtoL,MPI_COMM_WORLD);
          MPI_Recv(&u(0,k+0),1,MPI_type,left,LtoR,MPI_COMM_WORLD,&status);
        }
      /*exchange data with right neighbour*/
      if(taskid !=(numtasks-1))
        {
          MPI_Send(&u(npoints*4,k+0),1,MPI_type,right,LtoR,MPI_COMM_WORLD);
          MPI_Recv(&u(npoints*4+1,k+0),1,MPI_type,right,RtoL,MPI_COMM_WORLD,&status);

        }

 
      for (i=1; i<=lastpoint; i++) {    
	sum_t=0.0;
	for (j=1; j<k; j++){
	  sum_t += u(i,j) * coef(k-j);
	}    
	u(i,k+1)=r*(u(i+1,k)+u(i-1,k)) +(s-2*r)*u(i,k) + sum_t + (Jfunction(alpha,(k+1))- Jfunction(alpha,k))*u(i,0);
     
                                                 
      } 

      if(taskid==(numtasks-1))  
	u(x_N-1-first+2,k+1)=0.2;
	

    }               
   

 
 	   
  for (k=5; k<=nBlock_t*BlockWidth_t; k+=4){              /*compute other value */
    if(taskid !=0)              
      {
        
        MPI_Irecv(&u(0,k+0),1,MPI_type,left,LtoR,MPI_COMM_WORLD,&reqs[1]);
        MPI_Isend(&u(1,k+0),1,MPI_type,left,RtoL,MPI_COMM_WORLD,&reqs[0]);      /*exchange k+4, last k-loop */
	
      }
     
    if(taskid !=(numtasks-1))
      {
	MPI_Isend(&u(lastpoint,k+0),1,MPI_type,right,LtoR,MPI_COMM_WORLD,&reqs[2]);
	MPI_Irecv(&u(lastpoint+1,k+0),1,MPI_type,right,RtoL,MPI_COMM_WORLD,&reqs[3]);

      }
    
    for (i=1; i<=lastpoint; i+=4) {
     
      sum1=_mm256_setzero_pd();
      sum2=_mm256_setzero_pd();
      sum3=_mm256_setzero_pd();
      sum4=_mm256_setzero_pd();
      for (j=1; j<k; j++){
        
                
        load_u1=_mm256_broadcast_sd((real*)(u+i*(t_N+1)+j));        /*u(i,:) */
	load_u2=_mm256_broadcast_sd((real*)(u+(i+1)*(t_N+1)+j));    /*u(i+1,:) */
	load_u3=_mm256_broadcast_sd((real*)(u+(i+2)*(t_N+1)+j) );    /*u(i+2,:) */
	load_u4=_mm256_broadcast_sd((real*)(u+ (i+3)*(t_N+1)+j));    /*u(i+3,:) */    
																				
																				
	coef1= _mm256_load_pd(coef+k-j);
	coef2= _mm256_load_pd(coef+k-j);
	coef3= _mm256_load_pd(coef+k-j);
	coef4= _mm256_load_pd(coef+k-j);     

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
     

	   
    }  
    
    
    if(taskid !=0)
      {
	MPI_Wait(&reqs[0], &stats[0]);
	MPI_Wait(&reqs[1], &stats[1]);
      }
    if(taskid !=(numtasks-1))
      {
	MPI_Wait(&reqs[2], &stats[2]);
	MPI_Wait(&reqs[3], &stats[3]);
      }
 
    
    
    for(i=1;i<=lastpoint;i++) 
      {
    
	u(i,k+1)=r*(u(i+1,k)+u(i-1,k)) +(s-2*r)*u(i,k) + (outputbuffer[(i)*4]) + (Jfunction(alpha,(k+1))- Jfunction(alpha,k))*u(i,0) ;      
        
    
      }
    
    
    
 
    if(taskid==(numtasks-1))
      { 
	u(x_N-1-first+2,k+1)=0.2; 
      }        
     
    //exchange(k+1,first);
    if(taskid !=0)
      {
        MPI_Recv(&u(0,k+1),1,MPI_type,left,LtoR,MPI_COMM_WORLD,&status);
	MPI_Send(&u(1,k+1),1,MPI_type,left,RtoL,MPI_COMM_WORLD);
	
      }
    /*exchange data with right neighbour*/
    if(taskid !=(numtasks-1))
      {
	MPI_Send(&u(npoints*4,k+1),1,MPI_type,right,LtoR,MPI_COMM_WORLD);
	MPI_Recv(&u(npoints*4+1,k+1),1,MPI_type,right,RtoL,MPI_COMM_WORLD,&status);

      }
 
 
    /*compute u(i,k+2)*/

 
    for (i=1; i<=lastpoint; i++) { 
           
      outputbuffer[i*4+1]+=u(i,k)*coef(1);

      u(i,k+2)=r*(u(i+1,k+1)+u(i-1,k+1)) +(s-2*r)*u(i,k+1) +  outputbuffer[i*4+1] + (Jfunction(alpha,(k+2))- Jfunction(alpha,k+1))*u(i,0); 
            
   
    }
 
    if(taskid==(numtasks-1))  
      u(x_N-1-first+2,k+2)=0.2;
    // exchange(k+2,first);
 
    if(taskid !=0)

      {
        MPI_Recv(&u(0,k+2),1,MPI_type,left,LtoR,MPI_COMM_WORLD,&status);
	MPI_Send(&u(1,k+2),1,MPI_type,left,RtoL,MPI_COMM_WORLD);
	
      }
    /*exchange data with right neighbour*/
    if(taskid !=(numtasks-1))
      {
	MPI_Send(&u(npoints*4,k+2),1,MPI_type,right,LtoR,MPI_COMM_WORLD);
	MPI_Recv(&u(npoints*4+1,k+2),1,MPI_type,right,RtoL,MPI_COMM_WORLD,&status);

      }
 
 
    
    /*compute u(i,k+3)*/
 
    for (i=1; i<=lastpoint; i++) { 
	
      outputbuffer[i*4+2]+=u(i,k)*coef(2)+u(i,k+1)*coef(1);
      u(i,k+3)=r*(u(i+1,k+2)+u(i-1,k+2)) +(s-2*r)*u(i,k+2) + outputbuffer[i*4+2] + (Jfunction(alpha,(k+3))- Jfunction(alpha,k+2))*u(i,0) ;
    
    }
    
    if(taskid==(numtasks-1))  
      u(x_N-1-first+2,k+3)=0.2;
    
    //   exchange(k+3,first);
 
    if(taskid !=0)
      {
        MPI_Recv(&u(0,k+3),1,MPI_type,left,LtoR,MPI_COMM_WORLD,&status);
	MPI_Send(&u(1,k+3),1,MPI_type,left,RtoL,MPI_COMM_WORLD);
	
      }
    /*exchange data with right neighbour*/
    if(taskid !=(numtasks-1))
      {
	MPI_Send(&u(npoints*4,k+3),1,MPI_type,right,LtoR,MPI_COMM_WORLD);
	MPI_Recv(&u(npoints*4+1,k+3),1,MPI_type,right,RtoL,MPI_COMM_WORLD,&status);

      }
    
    /*compute u(i,k+4)*/
    
    
    
    for (i=1; i<=lastpoint; i++) { 
      	
      outputbuffer[i*4+3]+=u(i,k)*coef(3)+u(i,k+1)*coef(2)+u(i,k+2)*coef(1);
      u(i,k+4)=r*(u(i+1,k+3)+u(i-1,k+3)) +(s-2*r)*u(i,k+3) +outputbuffer[i*4+3]+ (Jfunction(alpha,(k+4))- Jfunction(alpha,k+3))*u(i,0) ;
    
    }
    if(taskid==(numtasks-1))  
      u(x_N-1-first+2,k+4)=0.2;
    /*compute remains of u(i,k+4)*/
    
  }
  
  


  if(taskid !=0)
    {
      MPI_Recv(&u(0,k+4),1,MPI_type,left,LtoR,MPI_COMM_WORLD,&status);
      MPI_Send(&u(1,k+4),1,MPI_type,left,RtoL,MPI_COMM_WORLD);
      
    }
  /*exchange data with right neighbour*/
  if(taskid !=(numtasks-1))
    {
      MPI_Send(&u(npoints*4,k+4),1,MPI_type,right,LtoR,MPI_COMM_WORLD);
      MPI_Recv(&u(npoints*4+1,k+4),1,MPI_type,right,RtoL,MPI_COMM_WORLD,&status);

    }

  if(taskid==(numtasks-1))  
    u(x_N-1-first+2,k+4)=0.2;


  /*compute the remains in time*/


  for (k=nBlock_t*BlockWidth_t; k<t_N; k++)              /*compute other value */
    {
 
      for (i=1; i<=lastpoint; i++) {    
	sum_t=0.0;
	for (j=1; j<k; j++){
	  sum_t += u(i,j) * coef(k-j);
	}    
	u(i,k+1)=r*(u(i+1,k)+u(i-1,k)) +(s-2*r)*u(i,k) + sum_t + (Jfunction(alpha,(k+1))- Jfunction(alpha,k))*u(i,0);  
                                                 
      }
      if(taskid==(numtasks-1))  
	u(x_N-1-first+2,k+1)=0.2;
	
	
      if(taskid !=0)
	{
	  MPI_Send(&u(1,k+1),1,MPI_type,left,RtoL,MPI_COMM_WORLD);
	  MPI_Recv(&u(0,k+1),1,MPI_type,left,LtoR,MPI_COMM_WORLD,&status);
	}
      /*exchange data with right neighbour*/
      if(taskid !=(numtasks-1))
	{
	  MPI_Send(&u(npoints*4,k+1),1,MPI_type,right,LtoR,MPI_COMM_WORLD);
	  MPI_Recv(&u(npoints*4+1,k+1),1,MPI_type,right,RtoL,MPI_COMM_WORLD,&status);

	}

    }               
   

   
  finish=MPI_Wtime();
  real time=1.0*(finish-begin);
     	
  
  real numer=0.0;
   
  real numer_all=0.0;
   
  
  for(k=1;k<=t_N;k++)                    
    for(i=1;i<=lastpoint;i++)      
      {  
	numer +=u(i,k)*u(i,k);
	 
      }
   

  MPI_Reduce(&numer, &numer_all, 1, MPI_type ,  MPI_SUM,0,MPI_COMM_WORLD);
   
     	    
  if(taskid==0)
    {
      printf("x_N: %d, t_N: %d, numtasks %d  ,Norm_n: %e, timeusage is %6.2f seconds\n\n",x_N, t_N,numtasks,sqrt(numer_all) ,time);
    }
 


 
  _mm_free(u);       
  _mm_free(coef);

  _mm_free(powfunction); 
  _mm_free(gammafunction );
   
  _mm_free(outputbuffer);
 
 
 
  MPI_Finalize(); 
  
  return EXIT_SUCCESS;     
}

