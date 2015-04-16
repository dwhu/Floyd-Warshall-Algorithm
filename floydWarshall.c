/*

Program kMeans.c
for K Means Lab

David Hughes
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef int bool;
enum { false, true };


/*
*   Helper Function to do a square root
*/
int sq_root(int n){
    int left, mid, right, mid2;

    right =1;
    while((right*right) <= n){
        right = (right << 1) + 1;
    }

    left = right;
    while((left * left) >= n){
        left >>=1;
    }

    while( (right-left) > 1){
        mid = (left + right) /2;
        mid2 = mid*mid;
        if(mid2 < n)
            left = mid;
        else if(mid2 == n)
            return mid;
        else
            right = mid;
    }

    return left;
}

/*
* Malloc Function that safely handles a malloc failure
* If Malloc fails, MPI is shutdown.
*/
void * safe_malloc(char * label, int bytes){

    void * p;
    p = malloc(bytes);
    if(p == NULL){
        fprintf(stderr, "malloc() failed on %s\n", label);
        MPI_Finalize();
        exit(3);
    }
    return p;
}

int main(int argc, char* argv[]){
    
    /* Local Info */
    int world_id;
    int nprocs;
    int n;
    int n_squared;
    double * local_section;
    int * my_row_group;
    int * my_column_group;
    MPI_Status status;


    /* MPI Group Info */
    MPI_Group world_group;
    MPI_Group row_group;
    MPI_Group column_group;
    MPI_Comm row_comm;
    MPI_Comm column_comm;

    /* Worker 0 */
    double * map;
    
    char * FNAME;

    if( argc < 1){
        fprintf(stderr, "Not Enough args\n");
        exit(1);
    }
    
    /***************************************************
                Setting Up Broadcast Groups
    ***************************************************/

    /* World */
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&world_id);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_group(MPI_COMM_WORLD,&world_group);

    /* Verify the correct number of workers to complete task
        ie. square_root(p) has no rounding error
    */
    int sqroot_nprocs = sq_root(nprocs);

    if(sqroot_nprocs*sqroot_nprocs != nprocs){
        if(world_id = 0){
            fprintf(stderr, "Number of procs '%d' us not a perfect square.\n", nprocs);
        }
        MPI_Finalize();
        exit(2);
    }

    int i;
    /* Row */
    my_row_group = safe_malloc("row group", sqroot_nprocs*sizeof(int));
    /* Find the workers in my row*/
    for(i =0; i < sqroot_nprocs; i++){
        my_row_group[i] = i + (world_id/sqroot_nprocs)*sqroot_nprocs;
    }
    /*Create the group and communication port*/
    MPI_Group_incl(world_group,sqroot_nprocs,my_row_group,&row_group);
    MPI_Comm_create(MPI_COMM_WORLD,row_group,&row_comm);

    /* Column */
    my_column_group = safe_malloc("column group", sqroot_nprocs*sizeof(int));
    // Find the workers in my column
    for(i =0; i < sqroot_nprocs; i++){
        my_column_group[i] = i*sqroot_nprocs + world_id %sqroot_nprocs;
    }
    /*Create the group and communication port*/
    MPI_Group_incl(world_group,sqroot_nprocs,my_column_group,&column_group);
    MPI_Comm_create(MPI_COMM_WORLD,column_group,&column_comm);

    MPI_Barrier(MPI_COMM_WORLD);

    // for(i = 0; i < nprocs; i++){
    //     if(i == world_id){
    //         printf("(%d) [",world_id);
    //         int x;
    //         for(x =0; x < sqroot_nprocs; x++){
    //             printf(" %d ",my_column_group[x]);
    //         }
    //         printf("]\n");
    //     }
    // }


    /***************************************************
                Read in and Distribute
    ***************************************************/
    if(world_id == 0 ){
        FNAME = argv[1];

        FILE * fp = fopen(FNAME,"r+");
        if(fp == NULL)
        {
            fprintf(stderr, "Unable to read file '%s'\n", FNAME);
            exit(2);
        }
        
        /* Read off the dimensions of the array */
        fscanf(fp,"%d",&n);

        /* Have to declare malloc.  Otherwise the OS might not give us enough space */
        map = (double*) safe_malloc("creating map",n*n*sizeof(double));

        int r = 1;
        n_squared = n*n;
        float tmp;
        int y;
        //Read in the map as a 1-D Array
        for(y = 0; y < n_squared && r != EOF; y++){
            r = fscanf(fp,"%f",&tmp);
            if ( r != EOF) {
                map[y] = tmp;  
            }
        }
        /* Close File */
        fclose(fp);
    }

    //Wait for Worker 0 to finish
    MPI_Barrier(MPI_COMM_WORLD);

    /***************************************************
                End Read in and Distribute
    ***************************************************/

    // /* Tell Everyone what n is */
    // MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // n_squared = n*n;

    // MPI_Barrier(MPI_COMM_WORLD);

    // double * buffer;
    // buffer = (double*) safe_malloc("creating buffer",n_squared/nprocs*sizeof(double));
    // MPI_Scatter((void *) map, n*n/nprocs, MPI_DOUBLE, (void *) buffer, n*n/nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // MPI_Barrier(MPI_COMM_WORLD);

    // for(i = 0; i < nprocs; i++){
    //     if(i == world_id){
    //         printf("(%d) n = %d\n",world_id,n);
    //         printf("(%d)",world_id);
    //         int x;
    //         for(x =0; x< n_squared/nprocs; x++){
    //             printf(" %f ",local_section[x]);
    //             if(x % n/sqroot_nprocs == 0){
    //                 printf("\n");
    //             }
    //         }
    //     }
    // }

    MPI_Finalize();
    
}

