/*
 * Author : Sai Kumar Manakan
 * UTA ID : 1001236131
 * Email  : saikumar.manakan@mavs.uta.edu
 *
 * Array shifting
 *
 * Working:-
 * 1. Processor 0 collects the 'Array size' input from user
 * 2. Processor 0 collects the data of each of the processor from user and then distributes to the
 *    respective processor
 * 3. Array data shift is then performed
 * 4. All Processors then send the resultant shifted array data to Processor 0
 * 5. Processor 0 then displays the resultant array
 *
 * Key MPI concepts used: MPI_Comm_size, MPI_Comm_rank, MPI_Bcast,
 *                        MPI_Send, MPI_Isend, MPI_Recv, MPI_Barrier
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void performArrayShift(int , int, int *, int, int);

void collectArrayData(int *, int, int);
void collectArrayDataFromProcessors(int , int *, int, int);

int main(int argc, char **argv)
{
    int rank, totalnodes, shiftFactor, i, p;
    int arraySize;

    int *myDataArray, *othersDataArray;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if(rank == 0)
    {
        printf("Array shift program");
        printf("\nEnter array size for each processor: ");
        fflush(stdout);

        scanf("%d",&arraySize);
        fflush(stdin);
    }

    /* Processor 0 broadcast 'arraySize' to all other processors */
    MPI_Bcast(&arraySize, 1, MPI_INT, 0, MPI_COMM_WORLD);


    // Quits in case of zero array size
    if(arraySize <= 0)
    {
        MPI_Finalize();
        return 0;
    }

    /* Processor 0 collects all the array data and sends to all processors */
    if(rank == 0)
    {
        myDataArray = (int *) malloc(arraySize * sizeof(int));
        collectArrayData(myDataArray, arraySize, 0);

        for(p = 1; p < totalnodes; ++p)
        {
            fflush(stdin);

            othersDataArray = (int *) malloc(arraySize * sizeof(int));
            collectArrayData(othersDataArray, arraySize, p);
            fflush(stdout);

            MPI_Send(othersDataArray, arraySize, MPI_INT, p, p, MPI_COMM_WORLD) ;
        }
    }
    else
    {
        MPI_Status status;

        /* Receive array data from Processor 0*/
        myDataArray = (int *) malloc(arraySize * sizeof(int));
        MPI_Recv(myDataArray, arraySize, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
        printf("\n######## Array ########");

    collectArrayDataFromProcessors(rank, myDataArray, arraySize, totalnodes);


    /* Array shifting */
    while(1)
    {
        if(rank == 0)
        {
            printf("\nEnter the shift factor (0 quits): ");
            fflush(stdout);
            scanf("%d", &shiftFactor);
        }

        /* Processor 0 broadcast 'shiftfactor' to all other processors */
        MPI_Bcast(&shiftFactor, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if(shiftFactor == 0)
            break;

        if(rank == 0)
            printf("\n######## Array before shift ########");

        collectArrayDataFromProcessors(rank, myDataArray, arraySize, totalnodes);

        performArrayShift(shiftFactor, rank, myDataArray, arraySize, totalnodes);

        MPI_Barrier(MPI_COMM_WORLD);

        if(rank == 0)
            printf("\n######## Array after shift ########");

        collectArrayDataFromProcessors(rank, myDataArray, arraySize, totalnodes);
    }

    MPI_Finalize();
    return 0;
}

void performArrayShift(int shiftFactor, int rank, int *myDataArray, int arraySize, int totalnodes)
{
    int data, i;
    MPI_Status status;
    MPI_Request request;

    //Positive shift
    if(shiftFactor > 0)
    {
        while(shiftFactor > 0)
        {
            if(rank == totalnodes - 1)
                MPI_Isend(myDataArray + arraySize - 1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
            else
                MPI_Isend(myDataArray + arraySize - 1, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &request);

            if(rank == 0)
            {
                MPI_Recv(&data, 1, MPI_INT, totalnodes - 1, 0, MPI_COMM_WORLD, &status);

                for(i = arraySize - 1; i > 0; --i)
                    myDataArray[i] = myDataArray[i - 1];

                myDataArray[0] = data;

            }
            else
            {
                MPI_Recv(&data, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);

                for(i = arraySize - 1; i > 0; --i)
                    myDataArray[i] = myDataArray[i - 1];

                myDataArray[0] = data;
            }

            MPI_Barrier(MPI_COMM_WORLD);
            --shiftFactor;
        }

    }
    else //Negative shift
    {
        while(shiftFactor < 0)
        {
            if(rank == 0)
                MPI_Isend(myDataArray, 1, MPI_INT, totalnodes - 1, 0, MPI_COMM_WORLD, &request);
            else
                MPI_Isend(myDataArray, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &request);

            if(rank == totalnodes - 1)
            {
                MPI_Recv(&data, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

                for(i = 0; i < arraySize - 1; ++i)
                    myDataArray[i] = myDataArray[i + 1];

                myDataArray[arraySize - 1] = data;

            }
            else
            {
                MPI_Recv(&data, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &status);

                for(i = 0; i < arraySize - 1; ++i)
                    myDataArray[i] = myDataArray[i + 1];

                myDataArray[arraySize - 1] = data;
            }

            MPI_Barrier(MPI_COMM_WORLD);
            ++shiftFactor;
        }
    }
}

void collectArrayDataFromProcessors(int rank, int *myDataArray, int arraySize, int totalnodes)
{
    /* Receive data from all other processors */
    if(rank == 0)
    {
        int index, k = 0, i = 0, p;
        int *dataArray = (int *) malloc(totalnodes * arraySize * sizeof(int));

        for(index = 0; index < (rank + 1) * arraySize; ++index)
        {
            dataArray[index] = myDataArray[k];
            ++k;
        }

        MPI_Status status;
        int *othersDataArray = (int *) malloc(arraySize * sizeof(int));


        for(p = 1; p < totalnodes; ++p)
        {
            MPI_Recv(othersDataArray, arraySize, MPI_INT, p, p, MPI_COMM_WORLD, &status) ;

            k = 0;

            for(; index < (p + 1) * arraySize; ++index)
            {
                dataArray[index] = othersDataArray[k];
                ++k;
            }
        }

        for(p = 0; p < totalnodes; ++p)
        {
            for(; i < (p + 1) * arraySize; ++i)
                printf("\n%d processor value : %d", p, dataArray[i]);
        }
    }
    else
    {
        MPI_Send(myDataArray, arraySize, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
}

void collectArrayData(int *dataArray, int arraySize, int rank)
{
    int i;
    printf("\nEnter values for processor %d", rank);

    for(i = 0; i < arraySize; ++i)
    {
        printf("\nEnter array value %d:", i);
        fflush(stdout);
        scanf("%d", &dataArray[i]);
    }
}
