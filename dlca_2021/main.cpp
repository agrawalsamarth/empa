#include "abdusalamov.h"

int main()
{
    dlca test;
    test.runSimulation();

    int N = dlca::N;
    int D = dlca::D;

    FILE *f;
    f = fopen("cluster.csv", "w");

    if (D==3){
        for (int i = 0; i < N; i++)
            fprintf(f, "%lf,%lf,%lf\n", test.pos(i,0),test.pos(i,1),test.pos(i,2));
    }

    else {
        for (int i = 0; i < N; i++)
            fprintf(f, "%lf,%lf\n", test.pos(i,0),test.pos(i,1));
    }

    fclose(f);

    return 0;
}
