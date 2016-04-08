/*
 * Author : Sai Kumar Manakan
 * UTA ID : 1001236131
 * Email  : saikumar.manakan@mavs.uta.edu
 *
 * MPI Program to compute Pi using Simpson's rule
 *
 * Working:-
 * 1. Processor 0 takes the 'number of intervals' input from user
 * 2. The summation part in the Simpson's rule is executed in parallel by all other processors
 * 3. The summation result is then sent to Processor 0 and the rest of the calculation is done by Processor 0
 * 4. Processor 0 then prints the final result
 *
 * Key MPI concepts used: MPI_Comm_size, MPI_Comm_rank, MPI_Bcast, MPI_Reduce
 */

#include <stdio.h>
#include <mpi.h>
#include <math.h>

/* Function to calculate Pi at an interva;l */
double f(double x) {return (4.0 / (1 + x * x));}

int main(int argc, char **argv)
{
    int rank, totalnodes, intervals, i;
    double sum, calc, partial_sum, pi, subIntervalSize, odd_pos, even_pos;
    double PI25DT = 3.141592653589793238462643;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0)
    {
        printf("Program to compute Pi using Simpson's rule");
        printf("\nEnter the number of intervals : ");
        fflush(stdout);
        scanf("%d", &intervals);

        if(intervals <= 0)
            printf("\nInterval cannot be 0 or negative. Exiting.\n");
    }

    /* Processor 0 broadcast 'intervals' to all other processors */
    MPI_Bcast(&intervals, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(intervals <= 0)
    {
        MPI_Finalize();
        return 0;
    }

    subIntervalSize = 1.0 / intervals;
    sum = 0.0;

    /* Simpson's rule
     * v = (1 / 3n) (f(0) + f(1) + summation from 1 to n/2 (4f(i-1) + 2f(i)))
     * summation from 1 to n/2 (4f(i-1) + 2f(i)) is calculated in parallel
     */
    for(i = rank + 1; i <= (intervals / 2); i += totalnodes)
    {
        even_pos = (2.0 * i);
        odd_pos = (2.0 * i) - 1;
        sum += (4.0 * f(odd_pos * subIntervalSize)) + (2.0 * f(even_pos * subIntervalSize));
    }

    MPI_Reduce(&sum, &partial_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
        /* Simpson's rule
         * v = (1 / 3n) (f(0) + f(1) + summation from 1 to n/2 (4f(i-1) + 2f(i)))
         * summation from 1 to n/2 (4f(i-1) + 2f(i)) is calculated in parallel
         */
        pi = (1.0 / (3 * intervals)) * (f(0.0) + partial_sum - f(1.0));

        printf("\nCalculated pi : %.16f, Error of %.16f\n", pi, fabs(pi - PI25DT));
    }

    MPI_Finalize();
    return 0;
}
