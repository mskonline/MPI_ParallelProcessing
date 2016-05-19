/*
 *
 * Serial and Parallel 1D BFS
 *	
 * Authors : 
 *		Gaurav Prakash,  <gaurav.prakash@mavs.uta.edu>	
 *		UTA ID : 1001236154
 *		
 *		Sai Kumar Manakan, <saikumar.manakan@mavs.uta.edu>
 *		UTA ID : 1001236131
 *
 *		Akash Gupta, <akash.gupta@mavs.uta.edu>
 *		UTA ID : 1001122031
 *
 *		Pavan Akshay Abhange, <pavanakshay.abhange@mavs.uta.edu>
 *		UTA ID : 1001103771
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <mpi.h>

// list
typedef struct node
{
    int data;
    struct node *next;
    struct node *prev;
}node;

typedef struct list
{
    int size;
    struct node *head;
    struct node *rear;
}list;

int exists_in_list(list *lst, int d){
    struct node *cur = lst->head;
    while(cur != NULL){
        if(cur->data == d)
            return 1;
        cur = cur->next;
    }
    return 0;
}

void print_list(list *lst){
    struct node *cur = lst->head;
    while(cur != NULL){
        printf("%d ", cur->data);
        cur = cur->next;
    }
    printf("\n");
}

int *array_from_list(list *lst){
    if(lst->size == 0){
        free(lst);
        return NULL;
    }
    int i = 0;
    int *arr = (int *)malloc(sizeof(int) * lst->size);
    struct node *cur = lst->head;
    struct node *garbage;
    while(cur != NULL){
        arr[i] = cur->data;
        garbage = cur;
        cur = cur->next;
        free(garbage);
    }
    free(lst);
    return arr;
}

void destroy_list(list *lst){
    struct node *garbage = lst->head;
    while(lst->head != NULL){
        lst->head = lst->head->next;
        free(garbage);
        garbage = lst->head;
    }
    free(lst);
}

list *create_list(){
    struct list *lst = (list *)malloc(sizeof(list));
    lst->head = NULL;
    lst->rear = NULL;
    lst->size = 0;
    
    return lst;
}

void push(list *lst, int d)
{
    struct node *new_node = (struct node *)malloc(sizeof(struct node));
    new_node->data = d;
    new_node->next = NULL;
    new_node->prev = lst->rear;

    if(!lst->rear){
        lst->head = new_node;
        lst->rear = new_node;
    }
    else
        lst->rear->next = new_node;

    lst->rear = new_node;
    lst->size += 1;
}

int pop(list *lst)
{
    if(!lst->rear)
        return -1;

    int d = lst->rear->data;

    struct node *tmp = lst->rear->prev;
    free(lst->rear);

    lst->rear = tmp;

    lst->size--;
    
    if(lst->size == 0)
        lst->head = NULL;

    return d;
}

int len(list *lst)
{
    if(lst)
        return lst->size;
    else
        return 0;
}

// queue

struct qnode
{
    int data;
    struct qnode *next;
};

struct qnode *HEAD = NULL;
struct qnode *REAR = NULL;

int qSize = 0;

void Qpush(int d)
{
    struct qnode *newNode = (struct qnode *) malloc(sizeof(struct qnode));
    newNode->data = d;
    newNode->next = NULL;

    if(!HEAD)
        HEAD = newNode;
    else
        REAR->next = newNode;

    REAR = newNode;
    ++qSize;
}

int Qpop()
{
    if(!HEAD)
        return -1;

    int d = HEAD->data;

    struct qnode *tmp = HEAD->next;
    free(HEAD);

    HEAD = tmp;

    --qSize;
    return d;
}

int Qsize()
{
    return qSize;
}

// Global variables
int world_size; // total number of processors
int num_of_vertices_per_proc; // number of vertices per processor
int num_of_vertices; // total number of vertices of the graph
double t_start, t_end;


// Prints the submatrix where num_of_vertices_per_proc is num of rows
// and num_of_vertices is num of columns
void print_sub_adj_matrix(unsigned char *sub_adj_matrix, int world_rank){
    int i, j, rank;
    MPI_Barrier(MPI_COMM_WORLD);
    for(rank = 0; rank < world_size; rank++){
        if(rank == world_rank){
            for(i = 0; i < num_of_vertices_per_proc; i++){
                for(j = 0; j < num_of_vertices ; j++){
                    printf("%d  ", *(sub_adj_matrix + i*num_of_vertices + j));
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// generates the frontier for the current level of BFS from
// level[] and current level l
// frontier is stored as a list
void generate_frontier(int *level, int l, int first_vertex, list *F){
    int i; // for iteration
    for(i = 0; i < num_of_vertices_per_proc; i++){
        if(*(level + i) == l)
            push(F, first_vertex + i);
    }
}

// finds the neighbors of a vertex
void get_neighbors(int vertex, int first_vertex, unsigned char *sub_adj_matrix, list *N){
    int row;
    int j;
    row = vertex - first_vertex;
    for(j = 0; j < num_of_vertices; j++){
        if(*(sub_adj_matrix + row * num_of_vertices + j) == 1){
            if(!exists_in_list(N, j))
                push(N, j);
        }
    }
}

// prints the distance of each vertex from source
void print_level_info(int *level, int world_rank){
    int i, rank;
    MPI_Barrier(MPI_COMM_WORLD);
    for(rank = 0; rank < world_size; rank++){
        if(rank == world_rank){
            for(i = 0; i < num_of_vertices_per_proc; i++){
                printf("%d  ", level[i]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

unsigned char* generateAdjMatrix(int);
void bfs(unsigned char*, int, int);

int main(int argc, char** argv){
    

    // declarations
    int i, j, k; // for iterations
    int world_rank; // rank of the processor
    int source_vertex; // source vertex for BFS
    int first_vertex; // index of first vertex that belong to the processor
    unsigned char *sub_adj_matrix; // store sub-adjacency matrix local to processor
    int *level; // store the length of path from start_vertex
    int l; // to keep track of level information for BFS
    struct list* F; // to store the Frontier in a list
    struct list* N; // to store the Neighbors in a list
    int frontier_has_vertices, continue_search; // for terminating the BFS
    int vertex; // to store the popped vertex from Frontier
    int neighbor, owner; // owner represents the owner of the vertex neighbor
    struct list **map; // stores neighbors that need to be sent to each proc
    int *send_counts; // number of neighbors that will be sent to each proc
    int *send_buf; // contains the data that need to be sent to each proc
    int *send_disps; // start index of the data that has to be sent to each proc
    int *recv_counts; // number of neighbors to be received from each proc
    int *recv_buf; // to store the data to be received from each proc
    int *recv_disps; // start index of storage for data received from each proc
    int send_buf_size = 0; // size of the send_buf
    int recv_buf_size = 0; // size of the recv_buf
    list *lst; // temporary storeage
    int rel_vertex; // to store the relative number of proc vertex in a proc


    if(argc != 3){
        fprintf(stderr, "Usage: a.out number_of_vertices source_vertex\n");
        exit(1);
    }

    srand(2);

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    unsigned char *adjMatrix;
    int p;

    if(world_rank == 0)
    {
        if(argc < 2)
        {
            printf("Please enter number of vertices per processor and start vertex\n");
            return 0;
        }

        num_of_vertices_per_proc = atoi(argv[1]);
        printf("Total vertices in the graph : %d\n", num_of_vertices_per_proc * world_size);
    }

    MPI_Bcast(&num_of_vertices_per_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(num_of_vertices_per_proc < 0)
    {
        MPI_Finalize();
        return 0;
    }

    source_vertex = atoi(argv[2]);

    //num_of_vertices_per_proc = num_of_vertices / world_size;
    num_of_vertices = num_of_vertices_per_proc * world_size;
    unsigned int dataLength = num_of_vertices * num_of_vertices_per_proc;


    first_vertex = world_rank * num_of_vertices_per_proc;

    if(world_rank == 0)
    {
        adjMatrix = generateAdjMatrix(num_of_vertices);

        // Perform serial BFS
        bfs(adjMatrix, num_of_vertices, source_vertex);

        sub_adj_matrix = (unsigned char *) malloc(num_of_vertices_per_proc * num_of_vertices * sizeof(unsigned char));
        memcpy(sub_adj_matrix, adjMatrix, dataLength);

        unsigned char* other_proc_adj_data = (unsigned char *) malloc(num_of_vertices_per_proc * num_of_vertices * sizeof(unsigned char));

        for(p = 1; p < world_size; ++p)
        {
            memcpy(other_proc_adj_data, (adjMatrix + p * dataLength), dataLength);
            MPI_Send(other_proc_adj_data, dataLength, MPI_UNSIGNED_CHAR, p, p, MPI_COMM_WORLD);
        }

        free(other_proc_adj_data);
        free(adjMatrix);
    }
    else
    {
        MPI_Status status;

        /* Receive sub adjacency matrix from Processor 0 */
        sub_adj_matrix = (unsigned char *) malloc(num_of_vertices_per_proc * num_of_vertices * sizeof(unsigned char));
        MPI_Recv(sub_adj_matrix, dataLength, MPI_UNSIGNED_CHAR, 0, world_rank, MPI_COMM_WORLD, &status);
    }

    MPI_Barrier(MPI_COMM_WORLD);
//    print_sub_adj_matrix(sub_adj_matrix, world_rank);

    // BFS

    // level[] stores the level information for each vertex
    // which represents the length of the path from start_vertex
    level = (int *)malloc(sizeof(int) * num_of_vertices_per_proc);

    // initialize level[]
    for(i = 0; i < num_of_vertices_per_proc; i++){
        // comparing source_vertex number with the vertex that the proc owns
        if(first_vertex + i == source_vertex){
            level[i] = 0;
        }
        else
            level[i] = -1;
    }

    // level-synchronous iteration for BFS traversal
    // iterating at the max for num_of_vertices - 1 times in the complete
    // graphs because the max length of BFS can be num_of_vertices - 1
    F = create_list();
    N = create_list();
    // for distributing neighbors from N to respective map[owner]
    map = (struct list**)malloc(sizeof(struct list *) * world_size);

    // count the total number of vertices that need to be sent to each proc
    // store the number of vertices that need to be sent to each proc
    send_counts = (int *)malloc(sizeof(int) * world_size);
    send_disps = (int *)malloc(sizeof(int) * world_size);

    // inform each porcessor, how many vertices to expect in the
    // next communication
    recv_counts = (int *)malloc(sizeof(int) * world_size);
    recv_disps = (int *)malloc(sizeof(int) * world_size);

    t_start = MPI_Wtime();
    for(l = 0; l < num_of_vertices - 1; l++){
        MPI_Barrier(MPI_COMM_WORLD);
        // frontier
        generate_frontier(level, l, first_vertex, F);
        if(len(F) > 0)
            frontier_has_vertices = 1;
        else
            frontier_has_vertices = 0;

        MPI_Allreduce(&frontier_has_vertices, &continue_search, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        // terminate if all frontiers are empty
        if(!continue_search){
            break;
        }


        // find neighbors of all the vertex in the frontier
        vertex = pop(F);
        while(vertex != -1){
            get_neighbors(vertex, first_vertex, sub_adj_matrix, N);
            vertex = pop(F);
        }


        // initialize the list for each index of the map
        for(i = 0; i < world_size; i++){
            map[i] = create_list();
        }

        // distribute
        neighbor = pop(N);
        while(neighbor != -1){
            owner = neighbor / num_of_vertices_per_proc;
            push(map[owner], neighbor);
            neighbor = pop(N);
        }

        // generate send_disps, send_counts and send_buf_size
        for(i = 0; i < world_size; i++){
            send_disps[i] = send_buf_size;
            send_buf_size += len(map[i]);
            send_counts[i] = len(map[i]);
        }

        // creating send_buf
        k = 0;
        send_buf = (int *)malloc(sizeof(int) * send_buf_size);
        for(i = 0; i < world_size; i++){
            lst = map[i];
            for(j = 0; j < send_counts[i] ; j++){
                send_buf[k] = pop(lst);
                k++;
            }
            free(lst);
            map[i] = NULL;
        }


        MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);
        // generate recv_disps and recv_buf_size
        for(i = 0; i < world_size; i++){
            recv_disps[i] = recv_buf_size;
            recv_buf_size += recv_counts[i];
        }

        // recv_buf stores the vertices that the proc owns, but were
        // in the Neighbor list of various processors and its own,
        // for the current level
        recv_buf = (int *)malloc(sizeof(int) * recv_buf_size);
        MPI_Alltoallv(send_buf, send_counts, send_disps, MPI_INT,
                recv_buf, recv_counts, recv_disps, MPI_INT, MPI_COMM_WORLD);


        for(i = 0; i < recv_buf_size; i++){
            // not the actual vertex number, but relative number on the proc
            rel_vertex = recv_buf[i] % num_of_vertices_per_proc;
            if(level[rel_vertex] == -1)
                level[rel_vertex] = l + 1; // increasing the level by 1
        }


        free(send_buf);
        free(recv_buf);
        // reinitializing the size counter
        send_buf_size = 0;
        recv_buf_size = 0;
    } // end for
    t_end = MPI_Wtime();

    if(world_rank == 0)
    {
        printf("Parallel BFS complete. Time taken = %lf\n", t_end - t_start);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // print_level_info(level, world_rank);
    // Clean up
    destroy_list(F);
    destroy_list(N);
    free(level);
    free(sub_adj_matrix);
    free(map);
    free(send_counts);
    free(send_disps);
    free(recv_counts);
    free(recv_disps);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void bfs(unsigned char* adjMatrix, int num_of_vertices, int source_vertex)
{
    int currentVertex, i;
    unsigned char *label = (unsigned char *) malloc(num_of_vertices * sizeof(unsigned char));
    unsigned char *visited = (unsigned char *) malloc(num_of_vertices * sizeof(unsigned char));

    for(i = 0; i < num_of_vertices; ++i)
        label[i] = -1;

    for(i = 0; i < num_of_vertices; ++i)
        visited[i] = 0;

    label[source_vertex] = 0;
    visited[source_vertex] = 1;
    Qpush(source_vertex);

    t_start = MPI_Wtime();

    while(Qsize() != 0)
    {
        currentVertex = Qpop();

        int j;
        for(j = 0; j < num_of_vertices; ++j)
        {
            if(*(adjMatrix + currentVertex * num_of_vertices + j) == 1)
            {
                if(visited[j] != 1)
                {
                    label[j] = label[currentVertex] + 1;
                    Qpush(j);
                    visited[j] = 1;
                }
            }
        }
    }

    t_end = MPI_Wtime();

    printf("Serial BFS complete. Time taken = %lf\n", t_end - t_start);

    free(label);
    free(visited);
}

unsigned char* generateAdjMatrix(int numVertices)
{
    int i, j;
    unsigned char *adjMatrix = (unsigned char *) malloc(numVertices * numVertices *  sizeof(unsigned char*));

    for (i = 0;i < numVertices; ++i)
    {
        for(j = 0; j < numVertices; ++j)
        {
            if(j >= i)
            {
                if(i == j)
                    *(adjMatrix + i * numVertices + j) = 0;
                else
                    *(adjMatrix + i * numVertices + j) = (rand() % 2);
            }
            else
                *(adjMatrix + i * numVertices + j) = *(adjMatrix + j * numVertices + i);
        }
    }

    return adjMatrix;
}
