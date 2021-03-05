#include <cmath>
#include <stdlib.h>
#include <iostream>

class dlca
{
  public:

    static const int N  = 10000;
    static const int D  = 2;
    const double rho    = 0.1;
    const double s_fraction = 0.1;
    const int N_s = (int)(N * s_fraction);
    const double R = 1.;


    dlca()
    {
        pos_           = (double*)malloc(sizeof(double) * N * D);
        clusters_S_    = (double*)malloc(sizeof(double) * N_s * D);
        seed_          = (int*)malloc(sizeof(int) * N);
        clusterNumber_ = (int*)malloc(sizeof(int) * N);
        srand48(1);

        L = std::pow((1. * N)/(rho), 1./(1. * D));
        initArrays();
        //std::cout<<"L ="<<L<<std::endl;

        totalClusters_ = N_s;
        totalWalkers_  = N - N_s;
        halfL_         = 0.5 * L;

        s[0] = 0.1 * R;
        s[1] = 0.05 * R;

        //std::cout<<"done"<<std::endl;


    }

    ~dlca()
    {
        free(pos_);
        free(seed_);
        free(clusterNumber_);
    }

    double
    pos(const int i, const int j) const
    { return pos_[i * D + j]; }
    double&
    pos(const int i, const int j)
    { return pos_[i * D + j]; }

    double
    clusters_S(const int i, const int j) const
    { return clusters_S_[i * D + j]; }
    double&
    clusters_S(const int i, const int j)
    { return clusters_S_[i * D + j]; }

    double
    periodic(const double x) const
    { return x + L * (x < 0) - L * (x > L); }

    void
    initArrays()
    {
        for (int i = 0; i < (N * D); i++)
            pos_[i] = drand48() * 1. * L;

        for (int i = 0; i < (N-N_s); i++) {
            seed_[i] = 0;
            clusterNumber_[i] = N;
        }

        for (int i = (N-N_s); i < N; i++) {
            seed_[i] = 1;
            clusterNumber_[i] = (i-(N-N_s));
        }

    }

    void
    buildClusterDisplacement2D() {

        for (int i = 0; i < N_s; i++) {

            theta_ = 2.0 * M_PI * drand48();
            //phi_   = 2.0 * M_PI * drand48();
            clusters_S(i, 0) = s[1] * std::cos(theta_);
            clusters_S(i, 1) = s[1] * std::sin(theta_);

        }

    }

    void
    buildClusterDisplacement3D() {

        for (int i = 0; i < N_s; i++) {

            theta_ = 2.0 * M_PI * drand48();
            phi_   = 2.0 * M_PI * drand48();
            clusters_S(i, 0) = s[1] * std::cos(theta_);
            clusters_S(i, 1) = s[1] * std::sin(theta_) * std::cos(phi_);
            clusters_S(i, 2) = s[1] * std::sin(theta_) * std::sin(phi_);

        }

    }


    void
    updatePositions2D () {

        for (int i = 0; i < N; i++) {

            if (seed_[i] == 0) {

                theta_    = 2.0 * M_PI * drand48();
                pos(i,0)  = periodic(pos(i,0) + (s[0] * std::cos(theta_)));
                pos(i,1)  = periodic(pos(i,1) + (s[0] * std::sin(theta_)));

            }

            else {

                clusterIndex_ = clusterNumber_[i];

                for (int j = 0; j < D; j++)
                    pos(i, j) = periodic(pos(i,j) + clusters_S(clusterIndex_, j));

            }

        }

    }

    void
    updatePositions3D () {

        for (int i = 0; i < N; i++) {

            if (seed_[i] == 0) {

                theta_    = 2.0 * M_PI * drand48();
                phi_      = 2.0 * M_PI * drand48();
                pos(i,0)  = periodic(pos(i,0) + (s[0] * std::cos(theta_)));
                pos(i,1)  = periodic(pos(i,1) + (s[0] * std::sin(theta_) * std::cos(phi_)));
                pos(i,2)  = periodic(pos(i,1) + (s[0] * std::sin(theta_) * std::sin(phi_)));

            }

            else {

                clusterIndex_ = clusterNumber_[i];

                for (int j = 0; j < D; j++)
                    pos(i, j) = periodic(pos(i,j) + clusters_S(clusterIndex_, j));

            }

        }

    }


    void
    bindWalker (const int walker, const int cluster)
    {
        seed_[walker] = 1;
        clusterNumber_[walker] = clusterNumber_[cluster];
    }

    void
    bindClusters (const int cluster_1, const int cluster_2)
    {
        for (int i = 0; i < N; i++) {
            if (clusterNumber_[i] == cluster_1)
                clusterNumber_[i] = cluster_2;
        }
    }

    void
    checkForAggregations() {

        for (int i = 0; i < N; i++) {
            for (int j = i+1; j < N; j++) {

                r2 = 0.;

                for (int axis = 0; axis < D; axis++) {
                    posDiff[axis] = periodic(pos(i,axis) - pos(j, axis));
                    r2           += posDiff[axis] * posDiff[axis];
                }

                if (r2 <= R) {

                    if ((seed_[i] == 0) && (seed_[j] == 1)) {
                        bindWalker(i,j);
                        totalWalkers_--;
                    }

                    if ((seed_[i] == 1) && (seed_[j] == 0)) {
                        bindWalker(j,i);
                        totalWalkers_--;
                    }

                    if ((seed_[i] == 1) && (seed_[j] == 1) && (clusterNumber_[i] != clusterNumber_[j])) {
                        clusterIndex_1 = clusterNumber_[i];
                        clusterIndex_2 = clusterNumber_[j];
                        bindClusters(clusterIndex_1, clusterIndex_2);
                        totalClusters_--;
                    }

                }

            }
        }

    }

    void
    iterationStep() {

	switch(D){
        	case 2: buildClusterDisplacement2D(); updatePositions2D(); break;
            case 3: buildClusterDisplacement3D(); updatePositions3D(); break;
	}

        checkForAggregations();

    }

    void
    runSimulation() {

        while ((totalWalkers_ != 0) || (totalClusters_ != 1)) {
            iterationStep();
            //std::cout<<"clusters = "<<totalClusters_<<std::endl;
        }

    }

private:

    double L;
    double *pos_;
    int    *seed_;
    int    *clusterNumber_;
    int    totalClusters_;
    int    totalWalkers_;
    double halfL_;
    double s[2];
    double posDiff[D];
    double *clusters_S_;
    int    clusterIndex_;
    int    clusterIndex_1;
    int    clusterIndex_2;
    double r2;
    double theta_;
    double phi_;

};
