#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

double max(double a, double b, double c){
	double ans=a;
	if (ans<b){
		ans=b;
	}
	if (ans<c){
		ans=c;
	}
	return ans;
}

double absolute(double a){
	double ans;
	if (a>=0){
		ans=a;
	}
	else{
		ans=-a;
	}
	return ans;
}

int main( int argc, char** argv ){
	double Lx = 0.40; // length along X (m)
	double Ly = 0.40; // length along Y (m)
	double delta_x = 0.002; // length of CV along X(m)
	double delta_y = 0.002; // length of CV along Y(m)
	int nx = Lx/delta_x; // no of CVs along X
	int ny = Ly/delta_y; // no of CVs along Y
	double nu = 0.004; // kinematic viscosity (m^2/s)
	double mu = nu; // dynamic viscosity ()
	double rho = 1; // density (kg/m^3)
	double rel_u=0.5;
	double rel_v=0.5;
	double rel_p=0.5;
	int Nx=nx+2;	
	int Ny=ny+2;
	int i,j,k,u_no,v_no,b_count,p_no,i_u,i_v,i_p,i_ulocal,i_vlocal,i_plocal;
	double max_b, globalmax_b, delx_e, delx_w, dely_n, dely_s, De_u, Dw_u, Dn_u, Ds_u, De_v, Dw_v, Dn_v, Ds_v, temp;
	double Fe_u, Fw_u, Fn_u, Fs_u, Fe_v, Fw_v, Fn_v, Fs_v, ae_u, aw_u, an_u, as_u, ae_v, aw_v, an_v, as_v, ap_u, ap_v;
	double ae_p, aw_p, an_p, as_p, ap_p;
	double p[ny][nx], u[Ny][Nx-1], v[Ny-1][Nx];
	double converge=1e-6;//Convergance criterion
	int NR_u,NC_u;
	int taskid, value, numtasks, offset;
    int i_first_u, i_last_u, i_first_u0, i_last_u0;
    int i_first_v, i_last_v, i_first_v0, i_last_v0;
    int NR_v, NC_v;
    int i_first_p, i_last_p, i_first_p0, i_last_p0;
    int NR_p, NC_p;
    MPI_Status status;
	NR_u=Ny;
	NC_u=Nx-1;
	NR_v=Ny;
	NC_v=Nx-1;
	NR_p=ny;
	NC_p=nx;
	MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
    MPI_Comm_size( MPI_COMM_WORLD, &numtasks );
	i_last_u = NR_u/numtasks;
	if(i_last_u<=2){
		if (taskid==0){
			printf("Too small problem\n");	
		}	
		MPI_Finalize( );
		return 0;
	}
		
    int extra = NR_u%numtasks;
    double ulocal[i_last_u+2+1][NC_u],Ap_u[i_last_u+2+1][NC_u],uc[i_last_u+2+1][NC_u];
    
	//Initialiaztion of u starting
    i_first_u = 1;
    if (taskid == numtasks - 1) i_last_u=i_last_u--;    	
    if (taskid<extra) i_last_u=i_last_u++;
    if (taskid == 0) i_last_u=i_last_u--;
    
    if (taskid==0){
    	for (j=0; j<NC_u; j++){
    		ulocal[i_first_u-1][j]=4;//Top BC
    	}
    }
    if (taskid==numtasks-1){
    	for (j=0; j<NC_u; j++){
    		ulocal[i_last_u+1][j]=0;//Bottom BC
    	}
    }
    for (i=i_first_u; i<=i_last_u; i++){
    	ulocal[i][0]=0;//Left BC
    	ulocal[i][NC_u-1]=0;//Right BC
    }
    for (i=i_first_u; i<=i_last_u; i++){
		for (j=1; j<NR_u-1; j++){
	    	ulocal[i][j] = 0.5*(ulocal[i][0]+ulocal[i][NC_u-1]);//Internal points initialization
	    }
	}
	//Initialiaztion of u completed

	//Initialiaztion of v starting
	i_last_v = NR_v/numtasks;
	if(i_last_v<=2){
		if (taskid==0){
			printf("Too small problem\n");	
		}	
		MPI_Finalize( );
		return 0;
	}		
    extra = NR_v%numtasks;
    double vlocal[i_last_v+2+1][NC_v],Ap_v[i_last_v+2+1][NC_v],vc[i_last_v+2+1][NC_v];

    i_first_v = 1;
    if (taskid == numtasks - 1) i_last_v=i_last_v--;    	
    if (taskid<extra) i_last_v=i_last_v++;
    if (taskid == 0) i_last_v=i_last_v--;
    
    if (taskid==0){
    	for (j=0; j<NC_v; j++){
    		vlocal[i_first_v-1][j]=0;//Top BC
    	}
    }
    if (taskid==numtasks-1){
    	for (j=0; j<NC_v; j++){
    		vlocal[i_last_v+1][j]=0;//Bottom BC
    	}
    }
    for (i=i_first_v; i<=i_last_v; i++){
    	vlocal[i][0]=0;//Left BC
    	vlocal[i][NC_v-1]=0;//Right BC
    }
    for (i=i_first_v; i<=i_last_v; i++){ 
		for (j=1; j<NR_v-1; j++){
	    	vlocal[i][j] = 0.5*(vlocal[i][0]+vlocal[i][NC_v-1]);//Internal points initialization
	    }
	}
    //Initialiaztion of v completed

	//Initialiaztion of p starting   
	i_last_p = NR_p/numtasks;
	if(i_last_p<=2){
		if (taskid==0){
			printf("Too small problem\n");	
		}	
		MPI_Finalize( );
		return 0;
	}
	extra = NR_p%numtasks;
    double plocal[i_last_p+2+1][NC_p],pc[i_last_p+2+1][NC_p],b[i_last_p+2+1][NC_p];
    
    i_first_p = 1;
    //i_last_p  = NR_p/numtasks;
    if (taskid == numtasks - 1) i_last_p=i_last_p;    	
    if (taskid<extra) i_last_p=i_last_p++;
    if (taskid == 0) i_last_p=i_last_p--; i_first_p=0;
    
    for (i=0; i<=i_last_p+1; i++){ 
		for (j=0; j<NC_p; j++){
	    	plocal[i][j] = 0;
	    	pc[i][j] = 0;
	    	b[i][j] =0;
	    }
	}   
	//Initialiaztion of p completed

	int count=0;
	while (1){
		count=count+1;
	    // X mom GS
	    for (u_no=1;u_no<6;u_no++){	    	
	    	MPI_Barrier(MPI_COMM_WORLD);
	    	extra=NR_u%numtasks;
	    	if (taskid < numtasks - 1) 
	    		MPI_Send( ulocal[i_last_u], NC_u, MPI_DOUBLE, taskid + 1, 0, MPI_COMM_WORLD );//Send last row to lower process u
			if (taskid > 0)
	    		MPI_Recv( ulocal[i_first_u-1], NC_u, MPI_DOUBLE, taskid - 1, 0, MPI_COMM_WORLD, &status );//Receive first row from upper process u 			
			if (taskid > 0) 
	    		MPI_Send( ulocal[i_first_u], NC_u, MPI_DOUBLE, taskid - 1, 1, MPI_COMM_WORLD );//Send first row to upper process u
			if (taskid < numtasks - 1) 
	    		MPI_Recv( ulocal[i_last_u+1], NC_u, MPI_DOUBLE, taskid + 1, 1, MPI_COMM_WORLD, &status );//Receive last row from lower process u
	    	if (taskid < extra-1) 
	    		MPI_Send( vlocal[i_last_v], NC_v, MPI_DOUBLE, taskid + 1, 2, MPI_COMM_WORLD );//Send last row to lower process v
	    	if (taskid < extra && taskid > 0)
	    		MPI_Recv( vlocal[0], NC_v, MPI_DOUBLE, taskid - 1, 2, MPI_COMM_WORLD, &status );//Receive first row from upper process v
	    	if (taskid < extra-2) 
	    		MPI_Send( plocal[i_last_p], NC_p, MPI_DOUBLE, taskid + 1, 3, MPI_COMM_WORLD );//Send last row to lower process p
	    	if (taskid < extra-1 && taskid > 0)
	    		MPI_Recv( plocal[0], NC_p, MPI_DOUBLE, taskid - 1, 3, MPI_COMM_WORLD, &status );//Receive first row from upper process p
	    	if (taskid == extra) 
	    		MPI_Send( plocal[1], NC_p, MPI_DOUBLE, taskid - 1, 4, MPI_COMM_WORLD );//Send first row to upper process p
	    	if (taskid == extra-1)
	    		MPI_Recv( plocal[i_last_u], NC_p, MPI_DOUBLE, taskid + 1, 4, MPI_COMM_WORLD, &status );//Receive first row from upper process p
	    	if (taskid > extra) 
	    		MPI_Send( plocal[1], NC_p, MPI_DOUBLE, taskid - 1, 5, MPI_COMM_WORLD );//Send first row to upper process p
	    	if (taskid < numtasks-1 && taskid > extra-1)
	    		MPI_Recv( plocal[i_last_u+1], NC_p, MPI_DOUBLE, taskid + 1, 5, MPI_COMM_WORLD, &status );//Receive last row from lower process p
	    	for (i=i_first_u;i<=i_last_u;i++){ 
				for (j=1;j<NC_u-1;j++){
					delx_e=delta_x;
	                delx_w=delta_x;
	                dely_n=delta_y;
	                dely_s=delta_y;
	                if (taskid==0 && i==i_first_u){ 
	                    dely_n=delta_y/2;
	                }
	                else if (taskid==numtasks-1 && i==i_last_u){ 
	                    dely_s=delta_y/2;
	                }
	                if (taskid<extra){
	                	i_v=i-1;
	                }
	                if (taskid>=extra){
	                	i_v=i;	
	                }
	                if (taskid<extra-1){
	                	i_p=i-1; //take i_last row of upper process
	                }
	                if (taskid==extra-1){
	                	i_p=i; //One row less in plocal; take i_first row of taskid=extra	
	                }
	                if (taskid>extra-1){
	                	i_p=i+1; //	take i_first row of lower process
	                }
	                De_u=mu/delx_e;
	                Dw_u=mu/delx_w;
	                Dn_u=mu/dely_n;
	                Ds_u=mu/dely_s;
	                Fe_u=rho*(ulocal[i][j]+ulocal[i][j+1])/2;
	                Fw_u=rho*(ulocal[i][j]+ulocal[i][j-1])/2;
	                Fn_u=rho*(vlocal[i_v][j]+vlocal[i_v][j+1])/2;
	                Fs_u=rho*(vlocal[i][j]+vlocal[i][j+1])/2;
	                ae_u=max(-Fe_u,De_u-Fe_u/2,0.0);
	                aw_u=max(Fw_u,Dw_u+Fw_u/2,0);	
	                an_u=max(-Fn_u,Dn_u-Fn_u/2,0);
	                as_u=max(Fs_u,Ds_u+Fs_u/2,0);
	                ap_u=ae_u+aw_u+an_u+as_u;
	                Ap_u[i][j]=ap_u;
	                
	                if (ap_u==0){//This is just error condition, put to avoid division by zero
	                	printf("Error: ap_u is zero %i %i \n",i,j);
	                }
	                ulocal[i][j]=rel_u*(ae_u*ulocal[i][j+1]+aw_u*ulocal[i][j-1]+an_u*ulocal[i-1][j]+as_u*ulocal[i+1][j]+(plocal[i_p][j-1]-plocal[i_p][j])*delta_y)/ap_u+(1-rel_u)*ulocal[i][j];
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		

	    // Y mom GS
	    for (v_no=1;v_no<6;v_no++){
	    	MPI_Barrier(MPI_COMM_WORLD);
	    	extra=NR_v%numtasks;
	    	if (taskid < numtasks - 1) 
	    		MPI_Send( vlocal[i_last_v], NC_v, MPI_DOUBLE, taskid + 1, 6, MPI_COMM_WORLD );//Send last row to lower process v
			if (taskid > 0)
	    		MPI_Recv( vlocal[i_first_v-1], NC_v, MPI_DOUBLE, taskid - 1, 6, MPI_COMM_WORLD, &status );//Receive first row from upper process v 			
			if (taskid > 0) 
	    		MPI_Send( vlocal[i_first_v], NC_v, MPI_DOUBLE, taskid - 1, 7, MPI_COMM_WORLD );//Send first row to upper process v
			if (taskid < numtasks - 1) 
	    		MPI_Recv( vlocal[i_last_v+1], NC_v, MPI_DOUBLE, taskid + 1, 7, MPI_COMM_WORLD, &status );//Receive last row from lower process v
	    	if (taskid < extra+1 && taskid > 0 ) 
	    		MPI_Send( ulocal[1], NC_u, MPI_DOUBLE, taskid - 1, 8, MPI_COMM_WORLD );//Send first row to upper process u
	    	if (taskid < extra)
	    		MPI_Recv( ulocal[i_last_v+1], NC_u, MPI_DOUBLE, taskid + 1, 8, MPI_COMM_WORLD, &status );//Receive last row from lower process u
	    	if (taskid < extra-1) 
	    		MPI_Send( plocal[i_last_p], NC_p, MPI_DOUBLE, taskid + 1, 9, MPI_COMM_WORLD );//Send last row to lower process p
	    	if (taskid < extra && taskid > 0)
	    		MPI_Recv( plocal[0], NC_p, MPI_DOUBLE, taskid - 1, 9, MPI_COMM_WORLD, &status );//Receive first row from upper process p
	    	if (taskid >= extra+1) 
	    		MPI_Send( plocal[1], NC_p, MPI_DOUBLE, taskid - 1, 10, MPI_COMM_WORLD );//Send first row to upper process p
	    	if (taskid >= extra && taskid < numtasks-1)
	    		MPI_Recv( plocal[i_last_v+1], NC_p, MPI_DOUBLE, taskid + 1, 10, MPI_COMM_WORLD, &status );//Receive last row from lower process p
	    	for (i=i_first_v;i<=i_last_v;i++){
	    		for (j=1;j<NC_v-1;j++){
	    			extra=NR_v%numtasks;
					delx_e=delta_x;
	                delx_w=delta_x;
	                dely_n=delta_y;
	                dely_s=delta_y;
	                if (j==1){
	                    delx_w=delta_x/2;
	                }
	                else if (j==Nx-1){ 
	                    delx_e=delta_x/2;
	                }
	                if (taskid<=extra){
	                	i_u=i+1; // take first row of lower process till taskid<extra
	                }
	                if (taskid>extra){
	                	i_u=i;
	                }
	                if (taskid< extra){
	                	i_p=i-1; //take last row of upper process
	                }
	                if (taskid>=extra){
	                	i_p=i; //take first row of lower process for taskid>=extra
	                }
	                De_v=mu/delx_e;
	                Dw_v=mu/delx_w;
	                Dn_v=mu/dely_n;
	                Ds_v=mu/dely_s;
	                Fe_v=rho*(ulocal[i][j]+ulocal[i+1][j])/2;
	                Fw_v=rho*(ulocal[i][j-1]+ulocal[i+1][j-1])/2;
	                Fn_v=rho*(vlocal[i][j]+vlocal[i-1][j])/2;
	                Fs_v=rho*(vlocal[i][j]+vlocal[i+1][j])/2;
	                ae_v=max(-Fe_v,De_v-Fe_v/2,0);
	                aw_v=max(Fw_v,Dw_v+Fw_v/2,0);
	                an_v=max(-Fn_v,Dn_v-Fn_v/2,0);
	                as_v=max(Fs_v,Ds_v+Fs_v/2,0);
	                ap_v=ae_v+aw_v+an_v+as_v;
	                Ap_v[i][j]=ap_v;
	                if (ap_v==0){//This is just error condition, put to avoid division by zero
	                	printf("Error: ap_v is zero %i %i \n",i,j);	                
	                }
	                v[i][j]=rel_v*(ae_v*v[i][j+1]+aw_v*v[i][j-1]+an_v*v[i-1][j]+as_v*v[i+1][j]+(p[i_p+1][j-1]-p[i_p][j-1])*delta_x)/ap_v+(1-rel_v)*v[i][j];
				}
			}
		}	            
		MPI_Barrier(MPI_COMM_WORLD);  
		 
		max_b=0;
		extra=NR_p%numtasks;
		if (taskid < extra+1 && taskid > 0 ) 
	    	MPI_Send( ulocal[1], NC_u, MPI_DOUBLE, taskid - 1, 11, MPI_COMM_WORLD );//Send first row to upper process u
	    if (taskid < extra)
	    	MPI_Recv( ulocal[i_last_u+1], NC_u, MPI_DOUBLE, taskid + 1, 11, MPI_COMM_WORLD, &status );//Receive last row from lower process u
		if (taskid> extra && taskid < numtasks - 1) 
	    	MPI_Send( ulocal[i_last_u], NC_u, MPI_DOUBLE, taskid + 1, 12, MPI_COMM_WORLD );//Send last row to lower process u
		if (taskid > extra+1)
	    	MPI_Recv( ulocal[0], NC_u, MPI_DOUBLE, taskid - 1, 12, MPI_COMM_WORLD, &status );//Receive first row from upper process u
	    if (taskid> extra-1 && taskid < numtasks - 1) 
	    	MPI_Send( vlocal[i_last_v], NC_v, MPI_DOUBLE, taskid + 1, 13, MPI_COMM_WORLD );//Send last row to lower process v
		if (taskid > extra)
	    	MPI_Recv( vlocal[0], NC_v, MPI_DOUBLE, taskid - 1, 13, MPI_COMM_WORLD, &status );//Receive first row from upper process v 
	    if (taskid < extra+1 && taskid > 0 ) 
	    	MPI_Send( vlocal[1], NC_v, MPI_DOUBLE, taskid - 1, 14, MPI_COMM_WORLD );//Send first row to upper process v
	    if (taskid < extra)
	    	MPI_Recv( vlocal[i_last_v+1], NC_v, MPI_DOUBLE, taskid + 1, 14, MPI_COMM_WORLD, &status );//Receive last row from lower process v
		for (i=1;i<=i_last_p;i++){
			for (j=0;j<NC_p-1;j++){
				if (taskid<=extra){
					i_u=i+1;//append last row to taskid<extra
					i_v=i;//for i_v+1 append last row to taskid<extra
				}
				if (taskid==extra+1){
					i_u=i;
					i_v=i-1;//append first row to taskid==extra+1
				}
				if (taskid>extra+1){
					i_u=i-1;//append first row to taskid>extra+1
					i_v=i-1;//append first row to taskid>extra+1
				}
				b[i][j]=(ulocal[i_u][j]-ulocal[i_u][j+1])*delta_y+(vlocal[i_v+1][j+1]-vlocal[i_v][j+1])*delta_x;
				if (max_b<absolute(b[i][j])){
					max_b=absolute(b[i][j]);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		double globalmax[numtasks];
		if (taskid!=0){		
    		MPI_Send(&max_b, 1, MPI_DOUBLE, 0, 31, MPI_COMM_WORLD);
	    }  
		if (taskid==0){
			globalmax[0]=max_b;
			for (i=1;i<=numtasks-1;i++){
    			MPI_Recv(&globalmax[i], 1, MPI_DOUBLE, i, 31, MPI_COMM_WORLD, &status);
    		}
    		globalmax_b=globalmax[0];
    		for (i=1;i<numtasks;i++){
    			if (globalmax_b<globalmax[i]){
    				globalmax_b=globalmax[i];
    			}
    		}
    		for (i=1;i<=numtasks-1;i++){
    			MPI_Send(&globalmax_b, 1, MPI_DOUBLE, i, 32, MPI_COMM_WORLD);
    		}
    	}
    	if (taskid!=0){		
			MPI_Recv(&globalmax_b, 1, MPI_DOUBLE, 0, 32, MPI_COMM_WORLD, &status);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
		if (taskid==0){
			printf("Iteration no.: %i\n", count);
		}
		if (globalmax_b<converge){//Convergence Criterion
			break;
		}
		// Pressure correction GS  
		else{
			for (i=i_first_p-1;i<=i_last_p+1;i++){
				for (j=0;j<NC_p;j++){ //initialise as zeros in each iteration
					pc[i][j]=0;
				}
			} 
			extra=NR_p%numtasks;			
			for (p_no=1;p_no<2;p_no++){				
				MPI_Barrier(MPI_COMM_WORLD);
				if (taskid < extra+1 && taskid > 0){ 					
	    			MPI_Send( Ap_u[1], NC_u, MPI_DOUBLE, taskid - 1, 15, MPI_COMM_WORLD );}//Send first row to upper process Ap_u  		
				if (taskid < extra) 
	    			MPI_Recv( Ap_u[i_last_u+1], NC_u, MPI_DOUBLE, taskid + 1, 15, MPI_COMM_WORLD, &status );//Receive last row from lower process Ap_u	    		
	    		if (taskid > extra+2 && taskid>0) 
	    			MPI_Send( Ap_u[i_last_u], NC_u, MPI_DOUBLE, taskid - 1, 16, MPI_COMM_WORLD );//Send first row to upper process Ap_u
				if (taskid < numtasks-1 && taskid > extra+1) 
	    			MPI_Recv( Ap_u[0], NC_u, MPI_DOUBLE, taskid + 1, 16, MPI_COMM_WORLD, &status );//Receive last row from lower process Ap_u
	    		//printf("OK till here extra: %i, taskid: %i\n",extra, taskid);
	    		if (taskid <extra+1 && taskid > 0) 
	    			MPI_Send( Ap_v[1], NC_v, MPI_DOUBLE, taskid - 1, 17, MPI_COMM_WORLD );//Send first row to upper process Ap_v
				if (taskid < extra) 
	    			MPI_Recv( Ap_v[i_last_v+1], NC_v, MPI_DOUBLE, taskid + 1, 17, MPI_COMM_WORLD, &status );//Receive last row from lower process Ap_v	    		
	    		if (taskid < numtasks-1 && taskid > extra-1) 
	    			MPI_Send( Ap_v[i_last_v], NC_v, MPI_DOUBLE, taskid + 1, 18, MPI_COMM_WORLD );//Send first row to upper process Ap_v
				if (taskid > extra) 
	    			MPI_Recv( Ap_v[0], NC_v, MPI_DOUBLE, taskid - 1, 18, MPI_COMM_WORLD, &status );//Receive last row from lower process Ap_v
				if (taskid < numtasks - 1) 
	    			MPI_Send( pc[i_last_p], NC_p, MPI_DOUBLE, taskid + 1, 19, MPI_COMM_WORLD );//Send last row to lower process pc
				if (taskid > 0)
	    			MPI_Recv( pc[i_first_p-1], NC_p, MPI_DOUBLE, taskid - 1, 19, MPI_COMM_WORLD, &status );//Receive first row from upper process pc			
				if (taskid > 0) 
	    			MPI_Send( pc[i_first_p], NC_p, MPI_DOUBLE, taskid - 1, 20, MPI_COMM_WORLD );//Send first row to upper process pc
				if (taskid < numtasks - 1) 
	    			MPI_Recv( pc[i_last_p+1], NC_p, MPI_DOUBLE, taskid + 1, 20, MPI_COMM_WORLD, &status );//Receive last row from lower process pc
				for (i=2;i<=i_last_p;i++){
					for (j=1;j<NC_p-1;j++){ 	
						if (taskid<=extra){
							i_u=i+1;//append last row to taskid<extra
							i_v=i;//for i+1 append last row to taskid<extra
						}
						if (taskid==extra+1){
							i_u=i;
							i_v=i-1;//append first row to taskid==extra+1
						}		
						if (taskid>extra+1){
							i_u=i-1;//append first row to taskid>extra+1
							i_v=i-1;//append first row to taskid>extra+1
						}		
    					ae_p=delta_y*delta_y/Ap_u[i_u][j+1];
    					aw_p=delta_y*delta_y/Ap_u[i_u][j];
    					an_p=delta_x*delta_x/Ap_v[i_v][j+1];
    					as_p=delta_x*delta_x/Ap_v[i_v+1][j+1];
    					if (Ap_u[i_u][j]==0){
    						printf("Ap_u[%i][%i]=0 taskid: %i\n", i_u,j,taskid);
    					}
    					if (Ap_u[i_u][j+1]==0){
    						printf("Ap_u[%i][%i]=0 taskid: %i\n", i_u,j+1,taskid);
    					}
    					if (Ap_v[i_v][j]==0){
    						printf("Ap_v[%i][%i]=0 taskid: %i\n", i_v,j,taskid);
    					}
    					if (Ap_u[i_v][j+1]==0){
    						printf("Ap_v[%i][%i]=0 taskid: %i\n", i_v,j+1,taskid);
    					}
    					if (i==0 && taskid==0){
    						an_p=0;
    					}
    					else if (i==i_last_p+1 && taskid==numtasks-1){
    						as_p=0;
    					}
    					if (j==0){
    						aw_p=0;
    					}
    					else if (j==NC_p-1){
    						ae_p=0;
    					}
    					ap_p=ae_p+aw_p+an_p+as_p;
    					if (ap_p==0){
    						printf("ap_p=0 for i=%i and j=%i, taskid=%i\n", i,j,taskid);
    					}
    					if (i==0 && taskid==0){
    						if (j==0){
    							pc[i][j]=(ae_p*pc[i][j+1]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    						else if (j==NC_p-1){
    							pc[i][j]=(aw_p*pc[i][j-1]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    						else{
    							pc[i][j]=(ae_p*pc[i][j+1]+aw_p*pc[i][j-1]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    					}
    					else if (i==i_last_p+1 && taskid==numtasks-1){
    						if (j==0){
    							pc[i][j]=(ae_p*pc[i][j+1]+an_p*pc[i-1][j]+b[i][j])/ap_p;
    						}
    						else if (j==NC_p-1){
    							pc[i][j]=(aw_p*pc[i][j-1]+an_p*pc[i-1][j]+b[i][j])/ap_p;
    						}
    						else{
    							pc[i][j]=(ae_p*pc[i][j+1]+aw_p*pc[i][j-1]+an_p*pc[i-1][j]+b[i][j])/ap_p;
    						}
    					}
    					else{
    						if (j==0){
    							pc[i][j]=(ae_p*pc[i][j+1]+an_p*pc[i-1][j]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    						else if (j==NC_p-1){
    							pc[i][j]=(aw_p*pc[i][j-1]+an_p*pc[i-1][j]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    						else{
    							pc[i][j]=(ae_p*pc[i][j+1]+aw_p*pc[i][j-1]+an_p*pc[i-1][j]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    					}
    				}
    			}
    		}
    	}
        MPI_Barrier(MPI_COMM_WORLD);                        
        // corrections
        for (i=0;i<=i_last_p+1;i++){
        	for (j=0;j<NC_p;j++){
        		p[i][j]=p[i][j]+rel_p*pc[i][j];
        	}
        }
        MPI_Barrier(MPI_COMM_WORLD);
        extra=NR_u%numtasks;
        if (taskid < extra-2) 
	    	MPI_Send( pc[i_last_p], NC_p, MPI_DOUBLE, taskid + 1, 21, MPI_COMM_WORLD );//Send last row to lower process pc
		if (taskid>0 && taskid <extra-1)
   			MPI_Recv( pc[0], NC_p, MPI_DOUBLE, taskid - 1, 21, MPI_COMM_WORLD, &status );//Receive first row from upper process pc
   		if (taskid > extra-1) 
	    	MPI_Send( pc[0], NC_p, MPI_DOUBLE, taskid - 1, 22, MPI_COMM_WORLD );//Send last row to lower process pc
		if (taskid>=extra-1 && taskid <numtasks-1)
   			MPI_Recv( pc[i_last_p+1], NC_p, MPI_DOUBLE, taskid + 1, 22, MPI_COMM_WORLD, &status );//Receive first row from upper process pc		
        for (i=1;i<=i_last_u;i++){
        	if (taskid<extra-1){
        		i_p=i-1;//append first row for taskid<extra-1
        	}
        	if (taskid==extra-1){
        		i_p=i;//append last row for taskid==extra-1
        	}
        	if (taskid>extra-1){
        		i_p=i+1;//append last row for taskid>extra-1
        	}        		
        	for (j=1;j<NC_u-1;j++){
        		if (Ap_u[i][j]==0){
    				printf("Ap_u[%i][%i]=0 taskid: %i\n", i,j,taskid);
    			}
        		uc[i][j]=delta_y*(pc[i_p][j-1]-pc[i_p][j])/Ap_u[i][j];
                u[i][j]=u[i][j]+uc[i][j];
            }
        }
        extra=NR_v%numtasks;
        if (taskid <= extra-2) 
	    	MPI_Send( pc[i_last_p], NC_p, MPI_DOUBLE, taskid + 1, 23, MPI_COMM_WORLD );//Send last row to lower process pc
		if (taskid>0 && taskid <=extra-1)
   			MPI_Recv( pc[0], NC_p, MPI_DOUBLE, taskid - 1, 23, MPI_COMM_WORLD, &status );//Receive first row from upper process pc
   		if (taskid > extra-1) 
	    	MPI_Send( pc[0], NC_p, MPI_DOUBLE, taskid - 1, 24, MPI_COMM_WORLD );//Send last row to lower process pc
		if (taskid>=extra-1 && taskid <numtasks-1)
   			MPI_Recv( pc[i_last_p+1], NC_p, MPI_DOUBLE, taskid + 1, 24, MPI_COMM_WORLD, &status );//Receive first row from upper process pc		
        for (i=1;i<=i_last_v;i++){
        	if (taskid<extra-1){
        		i_p=i-1;//append first row for taskid<extra-1
        	}
        	if (taskid==extra-1){
        		i_p=i-1;//for i append last row, for i-1, append first row
        	}
        	if (taskid>extra-1){
        		i_p=i;//for i append last row
        	}
        	for (j=1;j<NC_v-1;j++){
        		vc[i][j]=delta_x*(pc[i][j-1]-pc[i-1][j-1])/Ap_v[i][j];
                v[i][j]=v[i][j]+vc[i][j];    
            }
        }
        if(count>476){
        	break;
        }
	}
	if (taskid!=0){		
    	MPI_Send(&i_first_u, 1, MPI_INT, 0, 25, MPI_COMM_WORLD);
    	MPI_Send(&i_last_u, 1, MPI_INT, 0, 26, MPI_COMM_WORLD);
    	MPI_Send(ulocal[i_first_u], (i_last_u-i_first_u+1)*NC_u, MPI_DOUBLE, 0, 27,MPI_COMM_WORLD);
    }        
    if (taskid==0){
    	printf("Total Number of Iterations: %i\n", count);
    	offset=0;
    	for(i=0;i<=i_last_u;i++){    	    
    		for (j=0;j<NC_u;j++){
    			u[offset][j]=ulocal[i][j];
    		}
    		offset++;
	   	}
    	for (i=1;i<=numtasks-1;i++){
    		MPI_Recv(&i_first_u0, 1, MPI_INT, i, 25, MPI_COMM_WORLD, &status);
    		MPI_Recv(&i_last_u0, 1, MPI_INT, i, 26, MPI_COMM_WORLD, &status);
    		MPI_Recv(u[offset], (i_last_u0-i_first_u0+1)*NC_u, MPI_DOUBLE, i, 27,MPI_COMM_WORLD, &status);
    		offset=offset+i_last_u0-i_first_u0+1;
    	}
    }
    if (taskid!=0){		
    	MPI_Send(&i_first_v, 1, MPI_INT, 0, 28, MPI_COMM_WORLD);
    	MPI_Send(&i_last_v, 1, MPI_INT, 0, 29, MPI_COMM_WORLD);
    	MPI_Send(vlocal[i_first_v], (i_last_v-i_first_v+1)*NC_v, MPI_DOUBLE, 0, 30,MPI_COMM_WORLD);
    }        
    if (taskid==0){
    	offset=0;
    	//printf("i_last_u: %i\n", i_last_u);
    	for(i=0;i<=i_last_v;i++){    	    
    		for (j=0;j<NC_v;j++){
    			v[offset][j]=ulocal[i][j];
    		}
    		offset++;
    	//printf("offset: %i\n", offset);
    	}
    	for (i=1;i<=numtasks-1;i++){
    		MPI_Recv(&i_first_v0, 1, MPI_INT, i, 28, MPI_COMM_WORLD, &status);
    		MPI_Recv(&i_last_v0, 1, MPI_INT, i, 29, MPI_COMM_WORLD, &status);
    		MPI_Recv(v[offset], (i_last_v0-i_first_v0+1)*NC_v, MPI_DOUBLE, i, 30,MPI_COMM_WORLD, &status);
    		offset=offset+i_last_v0-i_first_v0+1;
    	}
    }
    MPI_Barrier(MPI_COMM_WORLD);
	    MPI_Finalize( );
		return 0;
}     
