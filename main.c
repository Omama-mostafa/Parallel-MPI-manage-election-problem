#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
struct Voters
{
    int *candidates;

} voters[0];
void printRandoms(FILE* fileobj,int lower, int upper, int voterNum)
{
    //FILE *fileobj ;
    int arr[upper];
    fileobj = fopen("Elections.txt", "w");
    fprintf(fileobj, "%d ", upper);
    fprintf(fileobj, "%d ", voterNum);
    fprintf(fileobj,"\n");

    int i,j,k;
    int counter;
    int flag;

    for(j = 0; j<voterNum ; j++)
    {
        for(k = 0 ; k<upper ; k++)
        {
            arr[k] = 0;
        }
        counter =0;
        for (i = 0; i < upper; i++)
        {
            int num = (rand() %(upper - lower + 1)) + lower;

            if(arr[num-1] == 1 && counter < upper) i--;

            else
            {
                fprintf(fileobj, "%d ", num);
                arr[num-1] = 1;
                counter++;
            }

        }
        fprintf(fileobj,"\n");
    }
    fclose(fileobj);
}

void Sort_arr(double arr[], int size, int candidateid[])
{
    int i, j, id;
    double temp;
    for(i=0; i<size; i++)
    {
        for(j=0; j<size; j++)
        {
            if(arr[i] > arr[j])
            {
                temp = arr[i];
                id = candidateid[i];
                candidateid[i]=candidateid[j];
                candidateid[j]=id;
                arr[i] = arr[j];
                arr[j] = temp;
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int rank, noPro, tag = 0;
    int i, j, k;
    FILE *myFile;

    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &noPro);


    int voterSize, candidateSize;
    int portion, rem;
    double maxi = 0.0;
    if(rank == 0)
    {   srand(time(NULL));
        int lower = 1, upper = (rand()% 7) + 2, voterNum = (rand()% 1300) + 100;
       // printf("voters %d \n",voterNum);

        printRandoms(myFile,lower, upper, voterNum);

        myFile = fopen("Elections.txt", "r");
        if (myFile == NULL)
        {
            printf("Error Reading File\n");
            exit (0);
        }
        fscanf(myFile, "%d", &candidateSize);
        fscanf(myFile, "%d", &voterSize);
        fclose(myFile);

        portion = voterSize / (noPro-1) ;
        rem = voterSize % (noPro-1);

        printf("\nvoter size %d ",voterSize);
        printf("\ncandidate size %d ",candidateSize);
        printf("\nportion %d ", portion);
        printf("\nRem %d\n", rem);
    }

    MPI_Bcast(&portion, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&candidateSize,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&voterSize,1,MPI_INT,0,MPI_COMM_WORLD);
    voters[portion];

    int candidateID[candidateSize];
    int candidateScore[candidateSize];
    for(i=0; i<candidateSize; i++)
    {
        candidateScore[i]=0;
        candidateID[i]=i+1;
    }
    if(rank !=0)
    {
        //ReadInSlaves(rank, i, j, k, myFile, voterSize, candidateSize, portion,voters);
        myFile = fopen("Elections.txt", "r");
        i = rank;
        fscanf(myFile, "%*s");
        int l=0;
        char line[256];
        for(k=0; k<portion; k++)
        {
            voters[k].candidates = (int*)malloc(candidateSize * sizeof(int));
        }
        int counter = 0;
        while (fgets(line, sizeof(line), myFile))
        {
            if(counter==voterSize)break;
            l++;
            if(l == (portion*(rank-1)+1))
            {
                l=0;
                for (i = 0; i < portion; i++)
                {
                    for(j=0; j<candidateSize; j++)
                    {
                        fscanf(myFile, "%d",&voters[i].candidates[j]);
                    }
                    if(i == portion)
                    {
                        l=0;
                        break;
                    }
                }
                break;
            }
            counter++;
        }
        fclose(myFile);

        for(k=0; k<candidateSize; k++)
        {
            candidateID[k] = k+1;
            candidateScore[k] =  0;
        }
        for(i=0; i<portion; i++)
        {
            int index = voters[i].candidates[0];

            candidateScore[index-1]++;
        }
    }

    int globalCandidateScore[candidateSize];
    for(i=0; i<candidateSize; i++)
    {
        globalCandidateScore[i]=0;
    }

    MPI_Reduce(candidateScore, globalCandidateScore, candidateSize, MPI_INT, MPI_SUM, 0,MPI_COMM_WORLD);
    if(rank == 0)
    {
        if(rem != 0)
        {
            //ReadInMaster(rank, i, j, k, myFile, voterSize, candidateSize, rem);
            myFile = fopen("Elections.txt", "r");
            fscanf(myFile, "%*s");
            for(k=0; k<rem; k++)
            {
                voters[k].candidates = (int*)malloc(candidateSize * sizeof(int));
                for(i=0; i<candidateSize; i++)
                {
                    voters[k].candidates[i]=0;
                }
            }

            int l=0;
            char line[256];

            while (fgets(line, sizeof(line), myFile))
            {
                l++;
                if(l == (voterSize-rem + 1))
                {

                    l=0;
                    for (i = 0; i < rem; i++)
                    {
                        for(j=0; j<candidateSize; j++)
                        {
                            fscanf(myFile, "%d",&voters[i].candidates[j]);

                        }
                        if(i == rem)break;
                    }
                    break;
                }
                if(l==(voterSize-rem+2))break;
            }
            fclose(myFile);
        }
        for(k=0; k<candidateSize; k++)
        {
            candidateID[k] = k+1;
            candidateScore[k] =  0;
        }
        for(i=0; i<rem; i++)
        {
            int index = voters[i].candidates[0];
            candidateScore[index-1]++;
        }

        for(i=0; i<candidateSize; i++)
        {
            globalCandidateScore[i] += candidateScore[i];
        }

        for(i=0; i<candidateSize; i++)
        {
            printf("\nFinalScore in first round = %d ,, ID : %d ", globalCandidateScore[i],i+1);
        }

        printf("\n");
        double probCandidate[candidateSize];
        for(i=0; i<candidateSize; i++)
            probCandidate[i] = 0.0;

        for(i=0; i<candidateSize; i++)
        {
            probCandidate[i] = ((double)globalCandidateScore[i]/ voterSize)*100.0;
        }

        Sort_arr(probCandidate, candidateSize,candidateID);
        for(i=0; i<candidateSize; i++)
        {
            printf("\nProbability[%d] = %16f , ID : %d \n", i, probCandidate[i],candidateID[i]);
        }
        maxi = probCandidate[0];

        if(maxi > 50.0)
        {
            printf("\n\nWinner with Score in first round with ID = %d and probablity = %16f \n",  candidateID[0], maxi);
        }
        else
        {
            printf("no winner at first round \n");
        }

    }
    MPI_Bcast(candidateID,candidateSize,MPI_INT,0,MPI_COMM_WORLD);

    if(maxi <= 50.0)
    {
        for(k=0; k<candidateSize; k++)
        {
            //candidateID[k] = k+1;
            candidateScore[k] =  0;
            globalCandidateScore[k] = 0;
        }
        if(rank != 0)
        {
            myFile = fopen("Elections.txt", "r");
            i = rank;
            fscanf(myFile, "%*s");
            int l=0;
            char line[256];
            for(k=0; k<portion; k++)
            {
                voters[k].candidates = (int*)malloc(candidateSize * sizeof(int));
            }
            int counter = 0;
            while (fgets(line, sizeof(line), myFile))
            {
                if(counter==voterSize)break;
                l++;
                if(l == (portion*(rank-1)+1))
                {
                    l=0;
                    for (i = 0; i < portion; i++)
                    {
                        for(j=0; j<candidateSize; j++)
                        {
                            fscanf(myFile, "%d",&voters[i].candidates[j]);
                        }
                        if(i == portion)
                        {
                            l=0;
                            break;
                        }
                    }
                    break;
                }
                counter++;
            }
            fclose(myFile);

            for(k=0; k<candidateSize; k++)
            {
                //candidateID[k] = k+1;
                candidateScore[k] =  0;
            }

            int**secondroundvote=malloc(portion*sizeof(int));
            for(i=0; i<portion; i++)
            {
                secondroundvote[i]=malloc(2*sizeof(int));
            }
            int flag =0;
            for(i=0; i<portion; i++)
            {
                flag =0;
                for(j=0; j<candidateSize; j++)
                {
                    if(voters[i].candidates[j]==candidateID[0]&& !flag)
                    {
                        int index = voters[i].candidates[j];
                        candidateScore[index-1]++;
                        flag =1;
                    }
                    else if(voters[i].candidates[j]==candidateID[1] && !flag)
                    {
                        int index = voters[i].candidates[j];
                        candidateScore[index-1]++;
                        flag = 1;
                    }
                }

            }
        }

        MPI_Reduce(candidateScore, globalCandidateScore, candidateSize, MPI_INT, MPI_SUM, 0,MPI_COMM_WORLD);
        if(rank == 0)
        {
            if(rem != 0)
            {
                //ReadInMaster(rank, i, j, k, myFile, voterSize, candidateSize, rem);
                myFile = fopen("Elections.txt", "r");
                fscanf(myFile, "%*s");
                for(k=0; k<rem; k++)
                {
                    voters[k].candidates = (int*)malloc(candidateSize * sizeof(int));
                    for(i=0; i<candidateSize; i++)
                    {
                        voters[k].candidates[i]=0;
                    }
                }

                int l=0;
                char line[256];
                while (fgets(line, sizeof(line), myFile))
                {
                    l++;
                    if(l == (voterSize-rem + 1))
                    {
                  //      printf("\nL = %d", l);
                        l=0;
                        for (i = 0; i < rem; i++)
                        {
                            for(j=0; j<candidateSize; j++)
                            {
                                fscanf(myFile, "%d",&voters[i].candidates[j]);

                            }
                            if(i == rem)break;
                        }
                        break;
                    }
                    if(l==(voterSize-rem+2))break;
                }
                fclose(myFile);

                for(k=0; k<candidateSize; k++)
                {
                    candidateScore[k] =  0;
                }

                int**secondroundvote=malloc(portion*sizeof(int));
                for(i=0; i<portion; i++)
                {
                    secondroundvote[i]=malloc(2*sizeof(int));
                }
                int flag =0;
                //printf("candidate id in rem: %d \n", candidateID[0]);
                for(i=0; i<rem; i++)
                {
                    flag =0;
                    for(j=0; j<candidateSize; j++)
                    {
                        if(voters[i].candidates[j]==candidateID[0]&& !flag)
                        {
                            int index = voters[i].candidates[j];
                            candidateScore[index-1]++;
                            flag =1;
                        }

                        else if(voters[i].candidates[j]==candidateID[1] && !flag)
                        {
                            int index = voters[i].candidates[j];
                            candidateScore[index-1]++;
                            flag = 1;
                        }
                    }


                }

                for(i=0; i<candidateSize; i++)
                {
                    globalCandidateScore[i] += candidateScore[i];
                }
            }
            for(i=0; i<candidateSize; i++)
            {if(globalCandidateScore[i]!=0)
                printf("\nFinalScore in second round = %d , ID = %d", globalCandidateScore[i],i+1);
            }

            printf("\n");
            double probCandidate[candidateSize];
            for(i=0; i<candidateSize; i++)
                probCandidate[i] = 0.0;

            for(i=0; i<candidateSize; i++)
            {
                probCandidate[i] = ((double)globalCandidateScore[i]/ voterSize)*100.0;
            }


            for(i=0; i<candidateSize; i++)
            {if(probCandidate[i] != 0.0)
                printf("\nProbability[%d] = %16f , ID : %d \n", i, probCandidate[i],i+1);
            }
            double max=-1.00;
            int id=0;
            for(i=0;i<candidateSize;i++)
            {
             if(probCandidate[i]>max)
             {max=probCandidate[i];
              id=i+1;
             }
            }
            printf("Winner in second round = %d ,, with Probablity = %16f \n", id, max);

        }
    }
    MPI_Finalize();
    return 0;
}
