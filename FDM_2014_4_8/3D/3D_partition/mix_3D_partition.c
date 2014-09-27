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
#include<omp.h> 
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
#pragma  message  "Make sure node_x*node_y*node_z=number of procs"



int taskid;                                                 
int numtasks;   
 
int npoints;
int first;
int RtoL=10;
int LtoR=20; 
int UptoDown=30;
int DowntoUp=40;
int FronttoBack=50;
int BacktoFront=60;
int left, right,up,down, front,back; 
int t_Nnew;

real *u; 
int x_length;
int y_length;
int z_length;
MPI_Comm comm;
MPI_Status  status;
int dim[3], period[3], reorder;
int coord[3], id;


#define sendbuffer_left(j,k)     sendbuffer_left[1ul*(z_length+1)*(j)+k]
#define recvbuffer_left(j,k)     recvbuffer_left[1ul*(z_length+1)*(j)+k]
#define sendbuffer_right(j,k)     sendbuffer_right[1ul*(z_length+1)*(j)+k]
#define recvbuffer_right(j,k)     recvbuffer_right[1ul*(z_length+1)*(j)+k]
#define sendbuffer_up(i,j)     sendbuffer_up[1ul*(y_length+2)*(i)+j]
#define recvbuffer_up(i,j)     recvbuffer_up[1ul*(y_length+2)*(i)+j]
#define sendbuffer_down(i,j)     sendbuffer_down[1ul*(y_length+2)*(i)+j]
#define recvbuffer_down(i,j)     recvbuffer_down[1ul*(y_length+2)*(i)+j]
#define sendbuffer_front(i,k)     sendbuffer_front[1ul*(z_length+1)*(i)+k]
#define recvbuffer_front(i,k)     recvbuffer_front[1ul*(z_length+1)*(i)+k]
#define sendbuffer_back(i,k)     sendbuffer_back[1ul*(z_length+1)*(i)+k]
#define recvbuffer_back(i,k)     recvbuffer_back[1ul*(z_length+1)*(i)+k]
#define  u(x,y,z,t)      u[1ul*(x)*(y_length+2)*(z_length+1)*(t_Nnew)+1ul*(y)*(z_length+1)*(t_Nnew)+1ul*(z)*(t_Nnew) +(t) ]
#define  coeff(t)       coeff[(t)]

real *sendbuffer_left;
real *recvbuffer_left;
real *sendbuffer_right;
real *recvbuffer_right;
real *sendbuffer_up;
real *recvbuffer_up;
real *sendbuffer_down;
real *recvbuffer_down;
real *sendbuffer_front;
real *recvbuffer_front;
real *sendbuffer_back;
real *recvbuffer_back;

int new_z;
real time_communication=0.0;
real time_pack=0.0;

void communication(int n)
{
    int i,j,k;
   real start;
   real finish;
   
   int id=omp_get_thread_num();

   start=omp_get_wtime();
#pragma omp for
 for(j=1;j<=y_length;j++)
     for(k=1;k<=new_z;k++)
        {
        sendbuffer_left(j,k)=u(1,j,k,n);
        sendbuffer_right(j,k)=u(x_length,j,k,n);
 
     }
   finish=omp_get_wtime();
#pragma omp master
{
time_pack+=(finish-start);
}

start=MPI_Wtime();//omp_get_wtime();

#pragma omp master
{
//x-axis -left
    MPI_Sendrecv(sendbuffer_left,(y_length+2)*(z_length+1),MPI_type,left,RtoL, recvbuffer_left, (y_length+2)*(z_length+1),MPI_type,left,LtoR,comm,&status);

 //x-axis -right

    MPI_Sendrecv(sendbuffer_right,(y_length+2)*(z_length+1),MPI_type,right,LtoR, recvbuffer_right, (y_length+2)*(z_length+1),MPI_type,right,RtoL,comm,&status);
}

  finish=MPI_Wtime();//omp_get_wtime();

#pragma omp master
{
time_communication+=(finish-start);
}


start=omp_get_wtime();

if(left !=MPI_PROC_NULL){
#pragma omp for
 for(j=1;j<=y_length;j++)
     for(k=1;k<=new_z;k++)
        u(0,j,k,n)=recvbuffer_left(j,k);
    }


if(right !=MPI_PROC_NULL){
#pragma omp for
 for(j=1;j<=y_length;j++)
     for(k=1;k<=new_z;k++)
        u(x_length+1,j,k,n)=recvbuffer_right(j,k);
    }

finish=omp_get_wtime();

#pragma omp master
{
time_pack+=(finish-start);
}


start=omp_get_wtime();

#pragma omp for
 for(i=1;i<=x_length;i++)
     for(k=1;k<=new_z;k++)
        {
        sendbuffer_front(i,k)=u(i,1,k,n);
        sendbuffer_back(i,k)=u(i,y_length,k,n);
     }

finish=omp_get_wtime();
#pragma omp master
{
time_pack+=(finish-start);
}


start=MPI_Wtime();//omp_get_wtime();
#pragma omp master
{
    MPI_Sendrecv(sendbuffer_front,(x_length+2)*(z_length+1),MPI_type,front,BacktoFront, recvbuffer_front,  (x_length+2)*(z_length+1),MPI_type,front,FronttoBack,comm,&status);

 //y-axis -back

    MPI_Sendrecv(sendbuffer_back,(x_length+2)*(z_length+1),MPI_type,back,FronttoBack, recvbuffer_back, (x_length+2)*(z_length+1),MPI_type,back,BacktoFront,comm,&status);

}


finish=MPI_Wtime();//omp_get_wtime();

#pragma omp master
{
time_communication+=(finish-start);
}
start=omp_get_wtime();
if(front !=MPI_PROC_NULL ){
#pragma omp for
for(i=1;i<=x_length;i++)
     for(k=1;k<=new_z;k++)
        u(i,0,k,n)=recvbuffer_front(i,k);
     }



if(back !=MPI_PROC_NULL ){
#pragma omp for
for(i=1;i<=x_length;i++)
     for(k=1;k<=new_z;k++)
        u(i,y_length+1,k,n)=recvbuffer_back(i,k);
     }

finish=omp_get_wtime();
#pragma omp master
{
time_pack+=(finish-start);
}

start=omp_get_wtime();

#pragma omp for
for(i=1;i<=x_length;i++)
     for(j=1;j<=y_length;j++)
        {
        sendbuffer_up(i,j)=u(i,j,new_z,n);
        sendbuffer_down(i,j)=u(i,j,1,n);
     }

finish=omp_get_wtime();
#pragma omp master
{
time_pack+=(finish-start);
}

start=MPI_Wtime();//omp_get_wtime();

#pragma omp master

{
    MPI_Sendrecv(sendbuffer_up,(x_length+2)*(y_length+2),MPI_type,up,DowntoUp, recvbuffer_up, (x_length+2)*(y_length+2),MPI_type,up,UptoDown,comm,&status);

 //z-axis -down

    MPI_Sendrecv(sendbuffer_down,(x_length+2)*(y_length+2),MPI_type,down,UptoDown, recvbuffer_down, (x_length+2)*(y_length+2),MPI_type,down,DowntoUp,comm,&status);

}


finish=MPI_Wtime();//omp_get_wtime();
#pragma omp master
{
time_communication+=(finish-start);
}
start=omp_get_wtime();
if(up !=MPI_PROC_NULL){
#pragma omp for
for(i=1;i<=x_length;i++)
     for(j=1;j<=y_length;j++)
        u(i,j,new_z+1,n)=recvbuffer_up(i,j);
}


if(down !=MPI_PROC_NULL){
#pragma omp for
for(i=1;i<=x_length;i++)
     for(j=1;j<=y_length;j++)
        u(i,j,0,n)=recvbuffer_down(i,j);
}

finish=omp_get_wtime();
#pragma omp master
{
time_pack+=(finish-start);
}
 
}
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
  x_N=200;                   //number of space intervals
  y_N=200;                    
  z_N=200;                    
  /*we choose x_N=y_N=z_N*/
  t_N=20;  	
  
              	
  real begin;
  real finish;
 
  //number of time intervals
  int node_x=3;
  int node_y=2;
  int node_z=3;

  if(nargs>1){
    x_N=atof(args[1]);
    y_N=atof(args[2]);
    z_N=atof(args[3]);
    t_N=atof(args[4]);
    node_x=atof(args[5]);
    node_y=atof(args[6]);
    node_z=atof(args[7]);
  }
 
 
  
 // MPI_Comm comm;
 // int dim[3], period[3], reorder;
 // int coord[3], id;
//  int left, right,up,down, front,back;

   int provided;
  //MPI_Init(&nargs, &args);
   MPI_Init_thread( &nargs, &args, MPI_THREAD_FUNNELED, &provided );
  dim[0]=node_x; dim[1]=node_y; dim[2]=node_z; 
  period[0]=0; period[1]=0; period[2]=0;
  reorder=1;
  MPI_Cart_create(MPI_COMM_WORLD, 3, dim, period, reorder, &comm);

  MPI_Comm_rank(comm, &taskid);
  MPI_Comm_size(comm,&numtasks);
  MPI_Cart_shift(comm, 0, 1, &left,&right);
  MPI_Cart_shift(comm, 1, 1,&front,&back);
  MPI_Cart_shift(comm, 2, 1, &down,&up );    


  
  real   dx;        dx= ((x_end-x_start)/x_N);          // space step x
  real   dy;         dy=  ((x_end-x_start)/y_N);         // space step y
  real   dz;         dz=  ((x_end-x_start)/z_N);          // space step z
  real   dt;        dt= ((t_end-t_start)/t_N);          // time step
  
  MPI_Cart_coords(comm, taskid, 3, coord);
#if 1
  //if (taskid==0)
    {
      //MPI_Cart_coords(comm, taskid, 3, coord);
      printf("\n");
   
      printf("node_x: %d, node_y: %d,  node_z: %d\n", node_x, node_y, node_z);
      printf("taskid: %d, left:%d, right:%d, down:%d, up:%d, front: %d, back:%d\n",taskid, left,right, down, up,front, back);
      printf("The processor at position (%d, %d, %d) has taskid %d\n", coord[0], coord[1], coord[2] , taskid);fflush(stdout);
     
      printf("Number of tasks %d\n", numtasks);fflush(stdout);
      printf("\n");
    }

#endif


 


  int threadNum;	
#pragma omp parallel
  {
#pragma omp master
    {
     threadNum=omp_get_num_threads();
      //printf("Number of threads = %d\n",threadNum);
    }
  }
  

  int nmin,nleft;     /*nmin: min points for each proc, nleft: the remainder*/
  int npts;

  int problem_size[3];
  problem_size[0]=x_N-1;
  problem_size[1]=y_N-1;
  problem_size[2]=z_N-1;

  int ii;
  int pp;

int my_start[3];
int my_finish[3]; 


  for(pp=0;pp<numtasks;pp++)
    {
      if(taskid==pp){

	for (ii=0;ii<3;ii++)
	  {
	    int num_x=problem_size[ii];
  
	    nmin=num_x/dim[ii];                    /*x_N+1 points, 0 and x_N are boundary*/
	    nleft=num_x%dim[ii];

	    for(i=0,k=0; i<dim[ii];i++)
	      {

		npts=(i<nleft)?nmin+1:nmin;
		if(coord[ii]==i)
		  {
		    first=k+1;   /*first index of u computed by  i  processor*/
		    npoints=npts;   /*npoints: number of points computed by the proc*/
                    my_start[ii]=first;
                    my_finish[ii]=first+npoints-1;

		   // printf("The processor at position (%d, %d, %d) has taskid %d, first: %d, npoints :%d\n", coord[0], coord[1], coord[2] , taskid,first, npoints);fflush(stdout);
                    
		  }
		else
		  k+=npts;

	      }
	  }
      }
    }


for (ii=0;ii<3;ii++)
 {

;//printf("The proc (%d, %d, %d), my_start[%d]: %d, my_finish[%d]: %d\n",coord[0], coord[1], coord[2],ii, my_start[ii],  ii, my_finish[ii]);

}



  new_z= my_finish[2]-my_start[2]+1;

 int BlockWidth_z=4; 

  int BlockWidth_t=4; 

  int nBlock_z=new_z/BlockWidth_z;    //x=0  x=x_N are known, x_N-1 points left

  int nBlock_t=(t_N-1)/BlockWidth_t;  // t=0  t=1 are computed separately, t_N-1 points left



  int NextStart=new_z+4-new_z%4+1;  

  int z_Nnew=NextStart+3;       //next i start position: x_N-1+4-(x_N-1)%4+1
  int NextStart_t=t_N-1+4-(t_N-1)%4+1;    // padding, add one block 

   t_Nnew=NextStart_t+3;       //next i start position: x_N-1+4-(x_N-1)%4+1

  int num_z=(NextStart-1)/4+1;      // z/4 consider padding
  
   //printf("num_z: %d\n", num_z);

  


  real A=pow(dt,-alpha)/tgamma(2-alpha);
  real B1=K/A/dx/dx;
  real B2=K/A/dy/dy;
  real B3=K/A/dz/dz;
  real B=1.0-2.0*(B1+B2+B3);

   x_length=my_finish[0]-my_start[0]+1;
   y_length=my_finish[1]-my_start[1]+1;
   z_length=num_z*4;
  
   printf("The processor at position (%d, %d, %d), x_length %d, y_length %d, z_length %d, rank%d\n", coord[0], coord[1], coord[2] , x_length, y_length, z_length,taskid);


  u=(real*)_mm_malloc(1ul*(x_length+2)*(y_length+2)*(z_length+1)*(t_Nnew)*sizeof(real),64); //  z_length+1: 1 more left point
  real *coeff=(real*)_mm_malloc(t_Nnew*sizeof(real),64);

 
   
  real *sum=(real*)_mm_malloc(1ul*(x_length+2)*(y_length+2)*(z_length+1)*sizeof(real),64);
  real *sum1=(real*)_mm_malloc(1ul*(x_length+2)*(y_length+2)*(z_length+1)*sizeof(real),64);
  real *sum2=(real*)_mm_malloc(1ul*(x_length+2)*(y_length+2)*(z_length+1)*sizeof(real),64);
  real *sum3=(real*)_mm_malloc(1ul*(x_length+2)*(y_length+2)*(z_length+1)*sizeof(real),64);
  real *outputbuffer=(real*)_mm_malloc(16*threadNum*64*sizeof(real),64);
//  real *outputbuffer=(real*)malloc(16*threadNum*64*sizeof(real));


 sendbuffer_left=(real*)malloc(1ul*(y_length+2)*(z_length+1)*sizeof(real));
recvbuffer_left=(real*)malloc(1ul*(y_length+2)*(z_length+1)*sizeof(real));
sendbuffer_right=(real*)malloc(1ul*(y_length+2)*(z_length+1)*sizeof(real));
recvbuffer_right=(real*)malloc(1ul*(y_length+2)*(z_length+1)*sizeof(real));
sendbuffer_up=(real*)malloc(1ul*(x_length+2)*(y_length+2)*sizeof(real));
recvbuffer_up=(real*)malloc(1ul*(x_length+2)*(y_length+2)*sizeof(real));
sendbuffer_down=(real*)malloc(1ul*(x_length+2)*(y_length+2)*sizeof(real));
recvbuffer_down=(real*)malloc(1ul*(x_length+2)*(y_length+2)*sizeof(real));
sendbuffer_front=(real*)malloc(1ul*(x_length+2)*(z_length+1)*sizeof(real));
recvbuffer_front=(real*)malloc(1ul*(x_length+2)*(z_length+1)*sizeof(real));
sendbuffer_back=(real*)malloc(1ul*(x_length+2)*(z_length+1)*sizeof(real));
recvbuffer_back=(real*)malloc(1ul*(x_length+2)*(z_length+1)*sizeof(real));


 
#define  sum(x,y,z)      sum[1ul*(x)*(y_length+2)*(z_length+1)+(y)*(z_length+1)+(z)]
#define  sum1(x,y,z)     sum1[1ul*(x)*(y_length+2)*(z_length+1)+(y)*(z_length+1)+(z)]
#define  sum2(x,y,z)     sum2[1ul*(x)*(y_length+2)*(z_length+1)+(y)*(z_length+1)+(z)]
#define  sum3(x,y,z)     sum3[1ul*(x)*(y_length+2)*(z_length+1)+(y)*(z_length+1)+(z)]



  real s=(pow(2,1-alpha)- pow(1,1-alpha));
	

#pragma omp parallel for default(shared) private(n)
  for(n=0;n<t_Nnew;n++)
    {
      coeff(n)=pow(n+2,1-alpha)+pow(n,1-alpha)-2.0*pow(n+1,1-alpha);
        
    }	

 

#pragma omp parallel for default(shared) private(i)
  for(i=0;i<(x_length+2)*(y_length+2)*(z_length+1)*(t_Nnew);i++)
    {
      u[i]=1.0;
	 
    }


if(left==MPI_PROC_NULL)
{ 
#pragma omp parallel for default(shared) private(j,k,p) collapse(3)
      for(j=0;j<=y_length+1;j++)             
	for(k=0;k<=new_z+1;k++)
	  for(p=0;p<=t_N;p++)
	    u(0,j,k,p)=0.2;
}	  

  
if(right==MPI_PROC_NULL)
{
#pragma omp parallel for default(shared) private(j,k,p) collapse(3)
      for(j=0;j<=y_length+1;j++)             
	for(k=0;k<=new_z+1;k++)
	  for(p=0;p<=t_N;p++)	  
	    u(x_length+1,j,k,p)=0.2;

}


if(front==MPI_PROC_NULL)

{
#pragma omp parallel for default(shared) private(i,k,p) collapse(3)
  for(i=0;i<=x_length+1;i++)                  
    for(k=0;k<=new_z+1;k++)
      for(p=0;p<=t_N;p++)
        {
	  u(i,0,k,p)=0.2;
	       
	}

}


if(back==MPI_PROC_NULL)
{
#pragma omp parallel for default(shared) private(i,k,p) collapse(3)
  for(i=0;i<=x_length+1;i++)                  
    for(k=0;k<=new_z+1;k++)
      for(p=0;p<=t_N;p++)
        {
	  u(i,y_length+1,k,p)=0.2;      
	}

}





if(up==MPI_PROC_NULL)
{
#pragma omp parallel for default(shared) private(i,j,p) collapse(3)
  for(i=0;i<=x_length+1;i++)                  
    for(j=0;j<=y_length+1;j++)
      for(p=0;p<=t_N;p++)
        {
	  
	  u(i,j,new_z+1,p)=0.2;      
	}
}



if(down==MPI_PROC_NULL)
{
#pragma omp parallel for default(shared) private(i,j,p) collapse(3)
  for(i=0;i<=x_length+1;i++)                  
    for(j=0;j<=y_length+1;j++)
      for(p=0;p<=t_N;p++)
        {
	  
	  u(i,j,0,p)=0.2;      
	}
}





 communication(0);
#pragma omp barrier
  begin=MPI_Wtime();

#pragma omp parallel for default(shared) private(i,j,k) collapse(3)
  for(i=1;i<=x_length;i++)                  
    for (j=1; j<=y_length; j++)
      for(k=1;k<=new_z;k++)
	{
                
	  u(i,j,k,1)= B*u(i,j,k,0)+ B1*(u(i+1,j,k,0)+u(i-1,j,k,0))+B2*(u(i,j+1,k,0)+u(i,j-1,k,0))+B3*(u(i,j,k+1,0)+u(i,j,k-1,0));
         }
  
real sum_t;

#pragma omp parallel default(shared) private(n,i,j,k,p,sum_t)
{
for(n=1;n<5;n++)
  {
   communication(n);
  #pragma omp barrier

             if(up==MPI_PROC_NULL){
#pragma omp for
 for(i=1;i<=x_length;i++)
  for(j=1;j<=y_length;j++)
        u(i,j,new_z+1,n)=0.2;
    }

   


 #pragma omp for
  for(i=1;i<=x_length;i++)            
    for (j=1; j<=y_length; j++)
      for(k=1;k<=new_z;k++)
        { 
           sum_t=0.0;
	  for(p=1; p<=n-1; p++)
	    {
	      sum_t +=coeff(n-p)*u(i,j,k,p);
                 
	    }  
              u(i,j,k,n+1)= B*u(i,j,k,n)-sum_t +u(i,j,k,0)*(pow(n+1,1-alpha)-pow(n,1-alpha))-u(i,j,k,n)*s+B1*(u(i+1,j,k,n)+u(i-1,j,k,n))+B2*(u(i,j+1,k,n)+u(i,j-1,k,n))+B3*(u(i,j,k+1,n)+u(i,j,k-1,n));   
                           
        }
}

}
#if 1
#pragma omp parallel default(shared) private(n,i,j,k,p)
  {	   
    for (n=5; n<=NextStart_t-1; n+=4)  
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
            

          int m;
          communication(n); 
        if(up==MPI_PROC_NULL){
#pragma omp for
 for(i=1;i<=x_length;i++)
  for(j=1;j<=y_length;j++)
        u(i,j,new_z+1,n)=0.2;
    }


     #pragma omp barrier     
     int id=omp_get_thread_num();
 #pragma omp for
         for(i=1;i<=x_length;i++)            
          for (j=1; j<=y_length; j++)
            for(k=1;k<=NextStart-1;k+=4)
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
      		load_u2= _mm256_broadcast_sd(&u(i,j,k+1,p));
		load_u3= _mm256_broadcast_sd(&u(i,j,k+2,p));
		load_u4= _mm256_broadcast_sd(&u(i,j,k+3,p));		 			       

		mulp1=_mm256_mul_pd(load_u1 ,load_coeff1);
		mulp2=_mm256_mul_pd(load_u2 ,load_coeff2);
		mulp3=_mm256_mul_pd(load_u3 ,load_coeff3);
		mulp4=_mm256_mul_pd(load_u4 ,load_coeff4);	   
          
		s1=_mm256_add_pd(s1,mulp1);   
		s2=_mm256_add_pd(s2,mulp2);   
		s3=_mm256_add_pd(s3,mulp3);  
		s4=_mm256_add_pd(s4,mulp4);    
  

	      }


	    _mm256_store_pd(&outputbuffer[0+id*16*64],s1);
	    _mm256_store_pd(&outputbuffer[4+id*16*64],s2);
	    _mm256_store_pd(&outputbuffer[8+id*16*64],s3);
	    _mm256_store_pd(&outputbuffer[12+id*16*64],s4);


	    sum(i,j,k)=outputbuffer[0+id*16*64];
	    sum1(i,j,k)=outputbuffer[1+id*16*64];
	    sum2(i,j,k)=outputbuffer[2+id*16*64];
	    sum3(i,j,k)=outputbuffer[3+id*16*64]; 
    
	    sum(i,j,k+1)=outputbuffer[4+id*16*64];
	    sum1(i,j,k+1)=outputbuffer[5+id*16*64];
	    sum2(i,j,k+1)=outputbuffer[6+id*16*64];
	    sum3(i,j,k+1)=outputbuffer[7+id*16*64]; 

	    sum(i,j,k+2)=outputbuffer[8+id*16*64];
	    sum1(i,j,k+2)=outputbuffer[9+id*16*64];
	    sum2(i,j,k+2)=outputbuffer[10+id*16*64];
	    sum3(i,j,k+2)=outputbuffer[11+id*16*64];
          
	    sum(i,j,k+3)=outputbuffer[12+id*16*64];
	    sum1(i,j,k+3)=outputbuffer[13+id*16*64];
	    sum2(i,j,k+3)=outputbuffer[14+id*16*64];
	    sum3(i,j,k+3)=outputbuffer[15+id*16*64];

    for(m=0;m<4;m++){
      u(i,j,k+m,n+1)= B*u(i,j,k+m,n)-sum(i,j,k+m)+u(i,j,k+m,0)*(pow(n+1,1-alpha)-pow(n,1-alpha))-u(i,j,k+m,n)*s+B1*(u(i+1,j,k+m,n)+u(i-1,j,k+m,n))+B2*(u(i,j+1,k+m,n)+u(i,j-1,k+m,n))+B3*(u(i,j,k+1+m,n)+u(i,j,k-1+m,n));
}
}


	
       communication(n+1); 
#pragma omp barrier
if(up==MPI_PROC_NULL){
#pragma omp for
 for(i=1;i<=x_length;i++)
  for(j=1;j<=y_length;j++)
        u(i,j,new_z+1,n+1)=0.2;
    }






#pragma omp for  collapse(3)
	for (i=1; i<=x_length; i++)     
	  for (j=1; j<=y_length; j++)
	    for(k=1;k<=new_z;k++)
	      {
		sum1(i,j,k)+=coeff(1)*u(i,j,k,n);    
		u(i,j,k,n+2)=B*u(i,j,k,n+1)-sum1(i,j,k) +u(i,j,k,0)*(pow(n+2,1-alpha)-pow(n+1,1-alpha))-u(i,j,k,n+1)*s+B1*(u(i+1,j,k,n+1)+u(i-1,j,k,n+1))+B2*(u(i,j+1,k,n+1)+u(i,j-1,k,n+1))+B3*(u(i,j,k+1,n+1)+u(i,j,k-1,n+1));
 
	      }


           communication(n+2);
#pragma omp barrier
#pragma omp for collapse(3)
	for (i=1; i<=x_length; i++)     
	  for (j=1; j<=y_length; j++)
	    for(k=1;k<=new_z;k++)
	      {
		sum2(i,j,k)+=coeff(2)*u(i,j,k,n)+coeff(1)*u(i,j,k,n+1);    
		u(i,j,k,n+3)=B*u(i,j,k,n+2)-sum2(i,j,k) +u(i,j,k,0)*(pow(n+3,1-alpha)-pow(n+2,1-alpha))-u(i,j,k,n+2)*s+B1*(u(i+1,j,k,n+2)+u(i-1,j,k,n+2))+B2*(u(i,j+1,k,n+2)+u(i,j-1,k,n+2))+B3*(u(i,j,k+1,n+2)+u(i,j,k-1,n+2)); 

	      }
           
       communication(n+3);
#pragma omp barrier
#pragma omp  for  collapse(3)   
	for (i=1; i<=x_length; i++)     
	  for (j=1; j<=y_length; j++)
	    for(k=1;k<=new_z;k++)
	      {
		sum3(i,j,k)+=coeff(3)*u(i,j,k,n)+coeff(2)*u(i,j,k,n+1)+coeff(1)*u(i,j,k,n+2);    
                u(i,j,k,n+4)=B*u(i,j,k,n+3)-sum3(i,j,k) +u(i,j,k,0)*(pow(n+4,1-alpha)-pow(n+3,1-alpha))-u(i,j,k,n+3)*s+B1*(u(i+1,j,k,n+3)+u(i-1,j,k,n+3))+B2*(u(i,j+1,k,n+3)+u(i,j-1,k,n+3))+B3*(u(i,j,k+1,n+3)+u(i,j,k-1,n+3)); 
                
	      }
 
  
}

}

#endif
  finish=MPI_Wtime();


  real time=1.0*(finish-begin);
     	
       
     		
  real numer=0.0;
  real numer_all=0.0;

  for(i=1;i<=x_length;i++) 
    for(j=1;j<=y_length;j++)
      for(k=1;k<=new_z;k++)     
	for(n=1;n<=t_N;n++) 
	  {  
	    numer +=u(i,j,k,n)*u(i,j,k,n);
	  }
//printf("Norm of numerical %e\n",sqrt(numer));     	

     	
  MPI_Reduce(&numer, &numer_all, 1, MPI_type ,  MPI_SUM,0,MPI_COMM_WORLD);
  

   	    
  if(taskid==0)
    {
 //     printf("%e\n",u(3,3,3,1851));   
      printf("x_N: %d,y_N: %d, z_N: %d, t_N: %d, numtasks %d,threads %d, total %d,time %5.2f secs, Norm_n: %e time_com %5.2f, time_pack %5.2f\n\n",x_N,y_N,z_N, t_N,numtasks,threadNum,numtasks*threadNum, time,sqrt(numer_all), time_communication,time_pack);
    }
  
  _mm_free(u);    
  _mm_free(coeff);
  _mm_free(sum);
  _mm_free(sum1);
  _mm_free(sum2);
  _mm_free(sum3);
  _mm_free(outputbuffer);
  free(sendbuffer_left);
  free(recvbuffer_left);

  free(sendbuffer_right);
  free(recvbuffer_right);

  free(sendbuffer_up);
  free(recvbuffer_up);

  free(sendbuffer_down);
  free(recvbuffer_down);

  free(sendbuffer_front);
  free(recvbuffer_front);

  free(sendbuffer_back);
  free(recvbuffer_back);
 


 
  MPI_Finalize(); 

  return EXIT_SUCCESS;     
}


