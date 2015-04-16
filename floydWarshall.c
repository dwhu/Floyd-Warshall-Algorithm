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
    int p;
    int n;
    int nn;
    int nnp;
    int n_sqp;
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
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_group(MPI_COMM_WORLD,&world_group);

    /* Verify the correct number of workers to complete task
        ie. square_root(p) has no rounding error
    */
    int sqroot_p = sq_root(p);

    if(sqroot_p*sqroot_p != p){
        if(world_id = 0){
            fprintf(stderr, "Number of procs '%d' us not a perfect square.\n", p);
        }
        MPI_Finalize();
        exit(2);
    }

    int i;
    /* Row */
    my_row_group = safe_malloc("row group", sqroot_p*sizeof(int));
    /* Find the workers in my row*/
    for(i =0; i < sqroot_p; i++){
        my_row_group[i] = i + (world_id/sqroot_p)*sqroot_p;
    }
    /*Create the group and communication port*/
    MPI_Group_incl(world_group,sqroot_p,my_row_group,&row_group);
    MPI_Comm_create(MPI_COMM_WORLD,row_group,&row_comm);

    /* Column */
    my_column_group = safe_malloc("column group", sqroot_p*sizeof(int));
    // Find the workers in my column
    for(i =0; i < sqroot_p; i++){
        my_column_group[i] = i*sqroot_p + world_id %sqroot_p;
    }
    /*Create the group and communication port*/
    MPI_Group_incl(world_group,sqroot_p,my_column_group,&column_group);
    MPI_Comm_create(MPI_COMM_WORLD,column_group,&column_comm);

    MPI_Barrier(MPI_COMM_WORLD);


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
        nn = n*n;
        nnp = nn/p;
        n_sqp = n / sqroot_p;
        float tmp;
        int y;
        //Read in the map as a 1-D Array
        for(y = 0; y < nn && r != EOF; y++){
            r = fscanf(fp,"%f",&tmp);
            if ( r != EOF) {
                map[y] = tmp;  
                printf(" %f ", map[y]);
            }
        }
        printf("\n\n");
        /* Close File */
        fclose(fp);
    }

    //Wait for Worker 0 to finish
    MPI_Barrier(MPI_COMM_WORLD);

    /***************************************************
                End Read in and Distribute
    ***************************************************/

    /* Tell Everyone what n is */
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    nn = n*n;
    nnp = nn/p;
    n_sqp = n / sqroot_p;

    local_section = (double*) safe_malloc("creating buffer",nnp*sizeof(double));
    MPI_Scatter((void *) map, nnp, MPI_DOUBLE, (void *) local_section, nnp, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // for(i = 0; i < p; i++){
    //     if(i == world_id){
    //         printf("(%d):",world_id);
    //         int x;
    //         for(x =0; x< nnp; x++){
    //             printf(" %f ",local_section[x]);
    //             if(x % sqroot_p ==0){
    //                 printf("\n");
    //             }
    //         }
    //         printf("\n\n");
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    MPI_Finalize();
    
}

