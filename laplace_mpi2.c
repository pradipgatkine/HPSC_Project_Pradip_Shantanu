#include <stdio.h>
#include <math.h>
#include "mpi.h"
//Taken from: http://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src/jacobi/C/main.html
/* This example handles a 12 x 12 mesh, on 4 processors only. */
//#define NC 202 
//#define NR 201 //keep NR greater than NC

int main( argc, argv )
int argc;
char **argv;
{
    int        taskid, value, numtasks, errcnt, toterr, i, j, k, offset;
    int        i_first, i_last, i_first0, i_last0;
    MPI_Status status;
    int NR=50;
    int NC=100;

    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
    MPI_Comm_size( MPI_COMM_WORLD, &numtasks );
	i_last = NR/numtasks;
	if(i_last<=2){
		if (taskid==0){
			printf("Too small problem\n");	
		}	
		MPI_Finalize( );
		return 0;
	}
		
    int extra = NR%numtasks;
    double xlocal[i_last+2+1][NC];
    double a[NR][NC];
    //if (numtasks != 4) MPI_Abort( MPI_COMM_WORLD, 1 );

    /* Note that top and bottom processes have one less row of interior points */
    i_first = 1;
    //i_last  = NR/numtasks;
    if (taskid == numtasks - 1) i_last=i_last--;    	
    if (taskid<extra) i_last=i_last++;
    if (taskid == 0) i_last=i_last--;
    
    if (taskid==0){
    	for (j=0; j<NC; j++){
    		xlocal[i_first-1][j]=100;//Top BC
    	}
    }
    if (taskid==numtasks-1){
    	for (j=0; j<NC; j++){
    		xlocal[i_last+1][j]=0;//Bottom BC
    	}
    }
    for (i=i_first; i<=i_last; i++){
    	xlocal[i][0]=0;//Left BC
    	xlocal[i][NC-1]=0;//Right BC
    }
    for (i=i_first; i<=i_last; i++) 
	for (j=1; j<NC-1; j++) 
	    xlocal[i][j] = 0.5*(xlocal[i][0]+xlocal[i][NC-1]);//Internal points initialization
    		
    for (k=0;k<5000;k++){
    MPI_Barrier(MPI_COMM_WORLD);
	/* Send up unless I'm at the top, then receive from below */
	/* Note the use of xlocal[i] for &xlocal[i][0] */
	if (taskid < numtasks - 1) 
	    MPI_Send( xlocal[i_last], NC, MPI_DOUBLE, taskid + 1, 0, MPI_COMM_WORLD );
	if (taskid > 0)
	    MPI_Recv( xlocal[i_first-1], NC, MPI_DOUBLE, taskid - 1, 0, MPI_COMM_WORLD, &status );
	/* Send down unless I'm at the bottom */
	if (taskid > 0) 
	    MPI_Send( xlocal[i_first], NC, MPI_DOUBLE, taskid - 1, 1, MPI_COMM_WORLD );
	if (taskid < numtasks - 1) 
	    MPI_Recv( xlocal[i_last+1], NC, MPI_DOUBLE, taskid + 1, 1, MPI_COMM_WORLD, &status );
	
	/* Compute new values (but not on boundary) */
	for (i=i_first; i<=i_last; i++){
	    for (j=1; j<NC-1; j++) {
			xlocal[i][j] = (xlocal[i][j+1] + xlocal[i][j-1] + xlocal[i+1][j] + xlocal[i-1][j]) / 4.0;
		}
	}
	if (taskid == 0) printf( "At iteration %d\n",k);
    } //Main iteration end
    
    if (taskid==numtasks-1){
    	printf("I_LAST: %f\n",xlocal[i_last-2][1]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (taskid!=0){
    	MPI_Send(&i_first, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    	MPI_Send(&i_last, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
    	MPI_Send(xlocal[i_first], (i_last-i_first+1)*NC, MPI_DOUBLE, 0, 0,MPI_COMM_WORLD);
    }
    if (taskid==0){
    	offset=0;
    	//printf("i_last: %i\n", i_last);
    	for(i=0;i<=i_last;i++){    	    
    		for (j=0;j<NC;j++){
    			a[offset][j]=xlocal[i][j];
    		}
    		offset++;
    		//printf("offset: %i\n", offset);
    	}
    	for (i=1;i<=numtasks-1;i++){
    		MPI_Recv(&i_first0, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
    		MPI_Recv(&i_last0, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
    		MPI_Recv(a[offset], (i_last0-i_first0+1)*NC, MPI_DOUBLE, i, 0,MPI_COMM_WORLD, &status);
    		if (i==numtasks-1){
    			printf("offset: %i\n", offset);
    			printf("i_first0: %i\n", i_first0);
    			printf("i_last0: %i\n", i_last0);
    		}
    		offset=offset+i_last0-i_first0+1;
    	}
    	printf("OK %f \n",a[NR-1][2]);
    	printf("Printing answer: \n");
		FILE *fp;
		fp=fopen("laplace_mpi2.dat","w"); //output will be stored in this file
		for(i=0;i<=NR-1;i++){
			for(j=0;j<=NC-1;j++){
				fprintf(fp,"%.3f\t",a[i][j]);
    		}
    		fprintf(fp,"\n");
		}
		fclose(fp);
		
    }	
    MPI_Finalize( );
    return 0;
}
