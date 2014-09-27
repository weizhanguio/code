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
#include<mpi.h> 
#include<immintrin.h>
#define  K             0.00001                   /* diffusion efficient*/
#define  x_start       0.0                    /*starting point of space*/
#define  x_end         10.0                      /*end point of space*/
#define  t_start       0.0                    /*starting point of time*/
#define  t_end         1.0     	       /*end point of time*/
#define  L           (x_end-x_start)
#define  MPI_type       MPI_DOUBLE 
#define  alpha         0.9
#define   real         double 




int taskid;                                                 
int numtasks;   
int left, right;
int npoints;
int first;
int RtoL=10;
int LtoR=20; 
real *u; 
   



int 
main(int nargs, char** args)
{
	
	
  int i;
  int j;
  int k;
  int p;
  int n;
  int m;
  int  x_N,y_N,z_N, t_N;
  x_N=20;                   //number of space intervals
  y_N=20;                    
  z_N=20;                    
  /*we choose x_N=y_N=z_N*/
  t_N=20;  	
  
              	
  real begin;
  real finish;
 
 

  if(nargs>1){
    x_N=atof(args[1]);
    y_N=atof(args[1]);
    z_N=atof(args[1]);
    t_N=atof(args[2]);

  }
 

  MPI_Init(&nargs, &args);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  real   dh;        dh= ((x_end-x_start)/x_N);          // space step
  real   dt;        dt= ((t_end-t_start)/t_N);          // time step
  

  int BlockWidth_x=4; 

  int BlockWidth_t=4; 

  int nBlock_x=(x_N-1)/BlockWidth_x;    //x=0  x=x_N are known, x_N-1 points left

  int nBlock_t=(t_N-1)/BlockWidth_t;  // t=0  t=1 are computed separately, t_N-1 points left

  int rem_x=(x_N-1)%BlockWidth_x;

  int rem_t=(t_N-1)%BlockWidth_t;


  int NextStart=x_N-1+4-(x_N-1)%4+1;  

  int x_Nnew=NextStart+3;       //next i start position: x_N-1+4-(x_N-1)%4+1
  int NextStart_t=t_N-1+4-(t_N-1)%4+1;    // padding, add one block 

  int t_Nnew=NextStart_t+3;       //next i start position: x_N-1+4-(x_N-1)%4+1

  int num_x=(NextStart-1)/4+1;
  
  int num_t=(NextStart_t-1)/4+1;
  

  
  
  MPI_Status stats[4];
  MPI_Request reqs[4];  
  MPI_Datatype type;
  MPI_Type_vector((y_N+1)*(z_N+1), 1, t_Nnew, MPI_DOUBLE, &type);
  MPI_Type_commit(&type);
  MPI_Status  status;
 
  
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
      printf("------3D AVX MPI nonblocking  Version----------\n");
      printf("%d intervals in space, %d intervals in time \n\n",x_N, t_N);
      /*x_N+1 points in space, t_N+1 points in time*/
    }


  int nmin,nleft;     /*nmin: min points for each proc, nleft: the remainder*/
  int npts;

  nmin=num_x/numtasks;                    /*x_N+1 points, 0 and x_N are boundary*/
  nleft=num_x%numtasks;

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


  int lastpoint;        //last point to be dealt with
  if(taskid!=numtasks-1)
    {
      lastpoint=npoints*4 ;
    }
  else
    {
      lastpoint=x_N-1-first+1;
    }



  real A=pow(dt,-alpha)/tgamma(2-alpha);
  real B=K/A/dh/dh;

  u=(real*)_mm_malloc(1ul*(4*npoints+2)*(y_N+1 )*(z_N+1 )*(t_Nnew)*sizeof(real),64);   
  real *coeff=(real*)_mm_malloc(t_Nnew*sizeof(real),64);

 double number_total=1.0*(4*npoints+2)*(y_N+1 )*(z_N+1 )*(t_Nnew)+4.0*(4*npoints+2)*(y_N+1 )*(z_N+1);
// printf(" myrank %d, number_total %f\n",taskid, number_total); 
   
  real *sum=(real*)_mm_malloc(1ul*(4*npoints+2)*(y_N+1 )*(z_N+1)*sizeof(real),64);
  real *sum1=(real*)_mm_malloc(1ul*(4*npoints+2)*(y_N+1 )*(z_N+1)*sizeof(real),64);
  real *sum2=(real*)_mm_malloc(1ul*(4*npoints+2 )*(y_N+1 )*(z_N+1)*sizeof(real),64);
  real *sum3=(real*)_mm_malloc(1ul*(4*npoints+2)*(y_N+1 )*(z_N+1)*sizeof(real),64);


  real* outputbuffer=(real*)_mm_malloc(16*sizeof(real),64);
  real s=(pow(2,1-alpha)- pow(1,1-alpha));
  
#define  sum(x,y,z)     sum[1ul*(x)*(y_N+1)*(z_N+1)+(y)*(z_N+1)+(z) ]
#define  sum1(x,y,z)    sum1[1ul*(x)*(y_N+1)*(z_N+1)+(y)*(z_N+1)+(z) ] 
#define  sum2(x,y,z)    sum2[1ul*(x)*(y_N+1)*(z_N+1)+(y)*(z_N+1)+(z) ]
#define  sum3(x,y,z)    sum3[1ul*(x)*(y_N+1)*(z_N+1)+(y)*(z_N+1)+(z) ]
#define  u(x,y,z,t)     u[1ul*(x)*(y_N+1)*(z_N+1)*(t_Nnew)+1ul*(y)*(z_N+1)*(t_Nnew)+(z)*(t_Nnew) +(t) ]
#define  coeff(t)       coeff[(t)] 


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
	


  for(n=0;n<t_Nnew;n++)
    {
      coeff(n)=pow(n+2,1-alpha)+pow(n,1-alpha)-2.0*pow(n+1,1-alpha);
        
    }	


  /*initial value*/

  for(i=0;i<(4*npoints+2)*(y_N+1 )*(z_N+1 )*(t_Nnew);i++)
    {
      u[i]=1.0;
	 
    }


  for(i=1;i<=lastpoint;i++)                  
    for(k=1;k<=z_N-1;k++)
      for(p=0;p<=t_N;p++)
        {
	  u(i,0,k,p)=0.2;
	  u(i,y_N,k,p)=0.2;      
	}

  for(i=1;i<=lastpoint;i++)                  
    for(j=1;j<=y_N-1;j++)
      for(p=0;p<=t_N;p++)
        {
	  u(i,j,0,p)=0.2;
	  u(i,j,z_N,p)=0.2;      
	}


 

  /*boundary value*/
  if(taskid==0)
    {
      for(j=0;j<=y_N;j++)             
	for(k=0;k<=z_N;k++)
	  for(p=0;p<=t_N;p++)
	    u(0,j,k,p)=0.2;
	  
    }


  if(taskid==numtasks-1)
    {
      for(j=0;j<=y_N;j++)             
	for(k=0;k<=z_N;k++)
	  for(p=0;p<=t_N;p++)	  
	    u(lastpoint+1,j,k,p)=0.2;
	  
    }



  

  if(taskid !=0)
    {
      
      MPI_Recv(&u(0,0,0,0),1,type,left,LtoR,MPI_COMM_WORLD,&status);
      MPI_Send(&u(1,0,0,0),1,type,left,RtoL,MPI_COMM_WORLD);
   
    }
  /*exchange data with right neighbour*/
  if(taskid !=(numtasks-1))
    {
      MPI_Send(&u(lastpoint,0,0,0),1,type,right,LtoR,MPI_COMM_WORLD);
      MPI_Recv(&u(lastpoint+1,0,0,0),1,type,right,RtoL,MPI_COMM_WORLD,&status);

    }


begin=MPI_Wtime();

  for(i=1;i<=lastpoint;i++)             /*compute  u(:,1)*/     
    for (j=1; j<=y_N-1; j++)
      for(k=1;k<=z_N-1;k++)
	{
                
	  u(i,j,k,1)= (1-6.0*B)*u(i,j,k,0)+ B*(u(i+1,j,k,0)+u(i-1,j,k,0)+u(i,j+1,k,0)+u(i,j-1,k,0)+u(i,j,k+1,0)+u(i,j,k-1,0));
      

	}
 



  if(taskid !=0)
    {
      
      MPI_Recv(&u(0,0,0,1),1,type,left,LtoR,MPI_COMM_WORLD,&status);
      MPI_Send(&u(1,0,0,1),1,type,left,RtoL,MPI_COMM_WORLD);
   
    }
  /*exchange data with right neighbour*/
  if(taskid !=(numtasks-1))
    {
      MPI_Send(&u(lastpoint,0,0,1),1,type,right,LtoR,MPI_COMM_WORLD);
      MPI_Recv(&u(lastpoint+1,0,0,1),1,type,right,RtoL,MPI_COMM_WORLD,&status);

    }




  real sum_t;


  for (n=1; n<5; n++) {             
    for (i=1; i<=lastpoint; i++)     
      for (j=1; j<=y_N-1; j++)
        for(k=1;k<=z_N-1;k++)
	  {
	    sum_t=0.0;
	    for(p=1; p<=n-1; p++)
	      {
		sum_t +=u(i,j,k,p)*coeff(n-p);
	      }
          
	    u(i,j,k,n+1)= (1-6.0*B)*u(i,j,k,n)-sum_t +u(i,j,k,0)*(pow(n+1,1-alpha)-pow(n,1-alpha))-u(i,j,k,n)*s+B*(u(i+1,j,k,n)+u(i-1,j,k,n)+u(i,j+1,k,n)+u(i,j-1,k,n)+u(i,j,k+1,n)+u(i,j,k-1,n)); 
                       
	  }

    if(taskid !=0)
      {
      
	MPI_Recv(&u(0,0,0,n+1),1,type,left,LtoR,MPI_COMM_WORLD,&status);
	MPI_Send(&u(1,0,0,n+1),1,type,left,RtoL,MPI_COMM_WORLD);
   
      }
    /*exchange data with right neighbour*/
    if(taskid !=(numtasks-1))
      {
	MPI_Send(&u(lastpoint,0,0,n+1),1,type,right,LtoR,MPI_COMM_WORLD);
	MPI_Recv(&u(lastpoint+1,0,0,n+1),1,type,right,RtoL,MPI_COMM_WORLD,&status);

      }


  }

 

	   
  for (n=5; n<=NextStart_t-1; n+=4)  
    {   
       
      if(taskid !=0)
	{
	  
	  MPI_Irecv(&u(0,0,0,n+0),1,type,left,LtoR,MPI_COMM_WORLD,&reqs[1]);
          MPI_Isend(&u(1,0,0,n+0),1,type,left,RtoL,MPI_COMM_WORLD,&reqs[0]);
      

	}
      /*exchange data with right neighbour*/
      if(taskid !=(numtasks-1))
	{
	  MPI_Isend(&u(lastpoint,0,0,n+0),1,type,right,LtoR,MPI_COMM_WORLD,&reqs[2]);
	  MPI_Irecv(&u(lastpoint+1,0,0,n+0),1,type,right,RtoL,MPI_COMM_WORLD,&reqs[3]);
    
	}


      /*compute the internal points*/
      for(i=1; i<=npoints*4; i+=4)           
	for (j=1; j<=y_N-1; j++)
	  for(k=1;k<=z_N-1;k++)
	    {
                  
	      s1=_mm256_setzero_pd();
	      s2=_mm256_setzero_pd();
	      s3=_mm256_setzero_pd();
	      s4=_mm256_setzero_pd();
		
	      for(p=1; p<=n-1; p++)
		{
             
		  load_coeff1=_mm256_load_pd(&(coeff(n-p))); 
                  load_coeff2=_mm256_load_pd(&(coeff(n-p)));
                  load_coeff3=_mm256_load_pd(&(coeff(n-p))); 
                  load_coeff4=_mm256_load_pd(&(coeff(n-p)));


		  load_u1= _mm256_broadcast_sd(&u(i,j,k,p));
		  load_u2= _mm256_broadcast_sd(&u(i+1,j,k,p));
		  load_u3= _mm256_broadcast_sd(&u(i+2,j,k,p));
		  load_u4= _mm256_broadcast_sd(&u(i+3,j,k,p));

		  mulp1=_mm256_mul_pd(load_u1 ,load_coeff1);
		  mulp2=_mm256_mul_pd(load_u2 ,load_coeff2);
		  mulp3=_mm256_mul_pd(load_u3 ,load_coeff3);
		  mulp4=_mm256_mul_pd(load_u4 ,load_coeff4);	   
		  s1=_mm256_add_pd(s1,mulp1); 
		  s2=_mm256_add_pd(s2,mulp2);  
		  s3=_mm256_add_pd(s3,mulp3);  
		  s4=_mm256_add_pd(s4,mulp4);   

		}

       
	      _mm256_store_pd(&outputbuffer[0],s1);
	      _mm256_store_pd(&outputbuffer[4],s2);
	      _mm256_store_pd(&outputbuffer[8],s3);
	      _mm256_store_pd(&outputbuffer[12],s4);

	      sum(i,j,k)=outputbuffer[0];
	      sum1(i,j,k)=outputbuffer[1];
	      sum2(i,j,k)=outputbuffer[2];
	      sum3(i,j,k)=outputbuffer[3]; 
   
	      sum(i+1,j,k)=outputbuffer[4];
	      sum1(i+1,j,k)=outputbuffer[5];
	      sum2(i+1,j,k)=outputbuffer[6];
	      sum3(i+1,j,k)=outputbuffer[7]; 
  
	      sum(i+2,j,k)=outputbuffer[8];
	      sum1(i+2,j,k)=outputbuffer[9];
	      sum2(i+2,j,k)=outputbuffer[10];
	      sum3(i+2,j,k)=outputbuffer[11];
          
	      sum(i+3,j,k)=outputbuffer[12];
	      sum1(i+3,j,k)=outputbuffer[13];
	      sum2(i+3,j,k)=outputbuffer[14];
	      sum3(i+3,j,k)=outputbuffer[15];

         
        
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



      for(i=1; i<=lastpoint; i++)
	for (j=1; j<=y_N-1; j++)
	  for(k=1;k<=z_N-1;k++)
	    {

	      u(i,j,k,n+1)= (1-6.0*B)*u(i,j,k,n)-sum(i,j,k) +u(i,j,k,0)*(pow(n+1,1-alpha)-pow(n,1-alpha))-u(i,j,k,n)*s+B*(u(i+1,j,k,n)+u(i-1,j,k,n)+u(i,j+1,k,n)+u(i,j-1,k,n)+u(i,j,k+1,n)+u(i,j,k-1,n));
	    }
	

	 
	 
      if(taskid==numtasks-1)
	{
	  for (j=0; j<=y_N; j++)
	    for(k=0;k<=z_N;k++) 
	      u(lastpoint+1,j,k,n+1)=0.2;
	}
	 
	
 
      if(taskid !=0)
	{
	 
	  MPI_Recv(&u(0,0,0,n+1),1,type,left,LtoR,MPI_COMM_WORLD,&status);
          MPI_Send(&u(1,0,0,n+1),1,type,left,RtoL,MPI_COMM_WORLD);
	}
      /*exchange data with right neighbour*/
      if(taskid !=(numtasks-1))
	{
	  MPI_Send(&u(lastpoint,0,0,n+1),1,type,right,LtoR,MPI_COMM_WORLD);
	  MPI_Recv(&u(lastpoint+1,0,0,n+1),1,type,right,RtoL,MPI_COMM_WORLD,&status);

	}


      for (i=1; i<=lastpoint; i++)     
	for (j=1; j<=y_N-1; j++)
	  for(k=1;k<=z_N-1;k++)
	    {
	      sum1(i,j,k)+=coeff(1)*u(i,j,k,n);    
	      u(i,j,k,n+2)=(1-6.0*B)*u(i,j,k,n+1)-sum1(i,j,k) +u(i,j,k,0)*(pow(n+2,1-alpha)-pow(n+1,1-alpha))-u(i,j,k,n+1)*s+B*(u(i+1,j,k,n+1)+u(i-1,j,k,n+1)+u(i,j+1,k,n+1)+u(i,j-1,k,n+1)+u(i,j,k+1,n+1)+u(i,j,k-1,n+1)); 
	      
                           
	    }

 
      if(taskid !=0)
	{
	 
	  MPI_Recv(&u(0,0,0,n+2),1,type,left,LtoR,MPI_COMM_WORLD,&status);
          MPI_Send(&u(1,0,0,n+2),1,type,left,RtoL,MPI_COMM_WORLD);
	}
      /*exchange data with right neighbour*/
      if(taskid !=(numtasks-1))
	{
	  MPI_Send(&u(lastpoint,0,0,n+2),1,type,right,LtoR,MPI_COMM_WORLD);
	  MPI_Recv(&u(lastpoint+1,0,0,n+2),1,type,right,RtoL,MPI_COMM_WORLD,&status);

	}


      for (i=1; i<=lastpoint; i++)     
	for (j=1; j<=y_N-1; j++)
	  for(k=1;k<=z_N-1;k++)
	    {
	      sum2(i,j,k)+=coeff(2)*u(i,j,k,n)+coeff(1)*u(i,j,k,n+1);    
	      u(i,j,k,n+3)=(1-6.0*B)*u(i,j,k,n+2)-sum2(i,j,k) +u(i,j,k,0)*(pow(n+3,1-alpha)-pow(n+2,1-alpha))-u(i,j,k,n+2)*s+B*(u(i+1,j,k,n+2)+u(i-1,j,k,n+2)+u(i,j+1,k,n+2)+u(i,j-1,k,n+2)+u(i,j,k+1,n+2)+u(i,j,k-1,n+2)); 
	                
	    }


      if(taskid !=0)
	{
	 
	  MPI_Recv(&u(0,0,0,n+3),1,type,left,LtoR,MPI_COMM_WORLD,&status);
	  MPI_Send(&u(1,0,0,n+3),1,type,left,RtoL,MPI_COMM_WORLD);
	}
      /*exchange data with right neighbour*/
      if(taskid !=(numtasks-1))
	{
	  MPI_Send(&u(lastpoint,0,0,n+3),1,type,right,LtoR,MPI_COMM_WORLD);
	  MPI_Recv(&u(lastpoint+1,0,0,n+3),1,type,right,RtoL,MPI_COMM_WORLD,&status);

	}

   
      for (i=1; i<=lastpoint; i++)     
	for (j=1; j<=y_N-1; j++)
	  for(k=1;k<=z_N-1;k++)
	    {
	      sum3(i,j,k)+=coeff(3)*u(i,j,k,n)+coeff(2)*u(i,j,k,n+1)+coeff(1)*u(i,j,k,n+2);    
	      u(i,j,k,n+4)=(1-6.0*B)*u(i,j,k,n+3)-sum3(i,j,k) +u(i,j,k,0)*(pow(n+4,1-alpha)-pow(n+3,1-alpha))-u(i,j,k,n+3)*s+B*(u(i+1,j,k,n+3)+u(i-1,j,k,n+3)+u(i,j+1,k,n+3)+u(i,j-1,k,n+3)+u(i,j,k+1,n+3)+u(i,j,k-1,n+3)); 
            
                
	    }
 


    }  


                  
  finish=MPI_Wtime();


  real time=1.0*(finish-begin);
     	

       
     		
  real numer=0.0;
  real numer_all=0.0;
     	

  for(i=1;i<=lastpoint;i++) 
    for(j=1;j<=y_N-1;j++)
      for(k=1;k<=z_N-1;k++)     
	for(n=1;n<=t_N;n++) 
	  {  
	    numer +=u(i,j,k,n)*u(i,j,k,n);
	  }
     	

     	
  MPI_Reduce(&numer, &numer_all, 1, MPI_type ,  MPI_SUM,0,MPI_COMM_WORLD);
     	    
  if(taskid==0)
    {
       
      printf("x_N: %d, t_N: %d, numtasks %d  ,timeusage%6.2f seconds, Norm_n:  %e \n\n",x_N, t_N,numtasks, time,sqrt(numer_all));
    }
  
 


  _mm_free(u);    
  _mm_free(coeff);
  _mm_free(sum);
  _mm_free(sum1);
  _mm_free(sum2);
  _mm_free(sum3);
  _mm_free(outputbuffer);

 
  MPI_Finalize(); 

  return EXIT_SUCCESS;     
}


