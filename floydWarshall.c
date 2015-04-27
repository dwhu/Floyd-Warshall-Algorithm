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
    int i,j;
    int index;
    int p;
    int n;
    int nn;
    int nnp;
    int n_sqp;
    int starting_i;
    int ending_i;
    int starting_j;
    int ending_j;
    double * local_section;
    double * row_section;
    double * column_section;
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
            if( y % n == 0){
                printf("\n");
            }
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
                End Read In to Master Worker
    ***************************************************/

    /* Tell Everyone what n is */
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    nn = n*n;
    nnp = nn/p;
    n_sqp = n / sqroot_p;


    /***************************************************
                Distribute Blocks to Workers
    ***************************************************/
    int dest;
    for(dest = p-1; dest >= 0; dest--){

        //World id cuts up the map for the current processor
        if(world_id == 0){
            local_section = (double*) safe_malloc("creating local_section buffer",nnp*sizeof(double));

            starting_i = (dest % sqroot_p)*n_sqp;
            ending_i = starting_i+n_sqp;
            starting_j = (dest / sqroot_p)*n_sqp;
            ending_j = starting_j+n_sqp;

            index = 0;
            for(j = starting_j; j < ending_j;j++){
                int j_index = j*n;
                for(i = starting_i;i < ending_i; i++){
                    local_section[index] = map[i+j_index];
                    index++;
                }
            }

            if(dest != 0){
                MPI_Send(local_section,nnp,MPI_DOUBLE,dest,dest,MPI_COMM_WORLD);
            }else{
                printf("(%d)\n",world_id);
                for(j=0;j<n_sqp;j++){
                    for(i=0;i<n_sqp;i++){
                        printf(" %f ", local_section[i+j*n_sqp]);
                    }
                    printf("\n");
                }
                printf("\n\n");
            }

        }else if( dest == world_id){

            local_section = (double*) safe_malloc("creating local_section buffer",nnp*sizeof(double));
            MPI_Recv(local_section,nnp,MPI_DOUBLE,0,dest,MPI_COMM_WORLD,&status);

            printf("(%d)\n",world_id);
            for(j=0;j<n_sqp;j++){
                for(i=0;i<n_sqp;i++){
                    printf(" %f ", local_section[i+j*n_sqp]);
                }
                printf("\n");
            }
            printf("\n\n");
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    int row_block = world_id/sqroot_p;
    int column_block = world_id % sqroot_p;
    int kb,kc;

    //Iteration Block from 0 to root(p)
    for(kb=0;kb<sqroot_p;kb++){

        if(world_id == 0){
            printf("\nKB = %d\n\n",kb);
        }

        //Iteration over a column/row workers partial segment block
        //i to n/root(p)
        for(kc=0;kc<n_sqp;kc++){

            MPI_Barrier(MPI_COMM_WORLD);

            /****************************************
                  Distirubte the Row and Columns
            *****************************************/

            //Allocate Memory
            row_section = (double *) safe_malloc("creating buffer", n_sqp*sizeof(double));

            //if this worker holds the row section
            if(row_block == kb){
                
                //Copy the Data in
                int target_row = kc*n_sqp;
                for(i=0; i < n_sqp; i++){
                    row_section[i] = local_section[i+target_row];
                }
            }
            MPI_Bcast(row_section,n_sqp,MPI_DOUBLE,kb,column_comm);
            MPI_Barrier(MPI_COMM_WORLD);

            column_section = (double *) safe_malloc("creating buffer", n_sqp*sizeof(double));

            //if this worker holds this column section
            if(kb == column_block){

                //Copy the data in
                for(j=0; j < n_sqp; j++){
                    int row = j*n_sqp;
                    column_section[j] = local_section[kc + row];

                }
            }

            //Broadcast both row and column to the correct group
            MPI_Bcast(column_section,n_sqp,MPI_DOUBLE,kb,row_comm);
            MPI_Barrier(MPI_COMM_WORLD);

            //Update
            for(i =0; i < n_sqp; i++){
                for(j=0;j< n_sqp;j++){
                    if(row_section[i]+column_section[j] < local_section[i*n_sqp+j] && row_section[i]+column_section[j] > 0){
                        local_section[i*n_sqp+j] = row_section[i]+column_section[j];
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);


            if(world_id == 0){
                printf("KC = %d\n",kc);
             }

            //Free the row and Column after done
            free((void *) column_section);
            free((void *) row_section);

        }
    }

    if(world_id == 0){
        index = 0;
        for(j=0;j<n_sqp;j++){
            for(i=0;i<n_sqp;i++){
                map[index] = local_section[i+j*n_sqp];
                index++;
            }
        }
    }

    for(dest = 1; dest < p; dest++){

        MPI_Barrier(MPI_COMM_WORLD);

        /*****************************
            Send Data Back to Master
        ******************************/
        if(world_id == dest){

            MPI_Send(local_section,nnp,MPI_DOUBLE,0,0,MPI_COMM_WORLD);

        }else if(world_id == 0){

            /******************************************
                Master Recieves Data and Puts it back
             ******************************************/
            local_section = (double*) safe_malloc("creating local_section buffer",nnp*sizeof(double));
            MPI_Recv(local_section,nnp,MPI_DOUBLE,dest,0,MPI_COMM_WORLD,&status);

            starting_i = (dest % sqroot_p)*n_sqp;
            ending_i = starting_i+n_sqp;
            starting_j = (dest / sqroot_p)*n_sqp;
            ending_j = starting_j+n_sqp;

            index = 0;
            for(j = starting_j; j < ending_j;j++){
                int j_index = j*n;
                for(i = starting_i;i < ending_i; i++){
                    map[i+j_index] = local_section[index];
                    index++;
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(world_id==0){
        printf("\n\n");
        for(j=0;j<n;j++){
            for(i=0;i<n;i++){
                printf(" %f ",map[i+j*n_sqp]);
            }
            printf("\n");
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    
}

