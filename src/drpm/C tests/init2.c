#include <stdio.h>

void initialize(int *Si_iter, int *gamma_iter, int *nh, int *nclus_iter, int *nsubject, int *ntime) {
    int ntime1 = *ntime + 1;  // Additional time period for scratch memory
    int j, t;

    // Initialize Si_iter, gamma_iter, and nh
    for (j = 0; j < *nsubject; j++) {
        for (t = 0; t < ntime1; t++) {
            Si_iter[j * ntime1 + t] = 1;
            gamma_iter[j * ntime1 + t] = 0;
            nh[j * ntime1 + t] = 0;
            if (t == 1) Si_iter[j * ntime1 + t] = 1; 
            if (t == *ntime) Si_iter[j * ntime1 + t] = 0;
        }
    }

    // Initial enumeration of number of subjects per cluster
    for (j = 0; j < *nsubject; j++) {
        for (t = 0; t < *ntime; t++) {
            nh[(Si_iter[j * ntime1 + t] - 1) * ntime1 + t]++;
        }
    }

    // Initialize the number of clusters
    for (t = 0; t < *ntime; t++) {
        nclus_iter[t] = 0;
        for (j = 0; j < *nsubject; j++) {
            if (nh[j * ntime1 + t] > 0) nclus_iter[t]++;
        }
    }
    nclus_iter[*ntime] = 0;
}

int main() {
    int nsubject = 3;
    int ntime = 4;
    int ntime1 = ntime + 1;

    int Si_iter[nsubject * ntime1];
    int gamma_iter[nsubject * ntime1];
    int nh[nsubject * ntime1];
    int nclus_iter[ntime1];

    // Call the initialization function
    initialize(Si_iter, gamma_iter, nh, nclus_iter, &nsubject, &ntime);

    // Print Si_iter
    printf("Si_iter:\n");
    for (int j = 0; j < nsubject; j++) {
        for (int t = 0; t < ntime1; t++) {
            printf("%d ", Si_iter[j * ntime1 + t]);
        }
        printf("\n");
    }

    // Print gamma_iter
    printf("gamma_iter:\n");
    for (int j = 0; j < nsubject; j++) {
        for (int t = 0; t < ntime1; t++) {
            printf("%d ", gamma_iter[j * ntime1 + t]);
        }
        printf("\n");
    }

    // Print nh
    printf("nh:\n");
    for (int j = 0; j < nsubject; j++) {
        for (int t = 0; t < ntime1; t++) {
            printf("%d ", nh[j * ntime1 + t]);
        }
        printf("\n");
    }

    // Print nclus_iter
    printf("nclus_iter:\n");
    for (int t = 0; t <= ntime; t++) {
        printf("%d ", nclus_iter[t]);
    }
    printf("\n");

    return 0;
}
