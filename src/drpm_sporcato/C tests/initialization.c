#include <stdio.h>
#include <stdlib.h>

void initialize_and_print(int *nsubject, int *ntime, int ntime1) {
    int *Si_iter, *gamma_iter, *nh, *nclus_iter;
    int j, t;

    // Allocate memory for the arrays
    Si_iter = (int *)malloc((*nsubject) * ntime1 * sizeof(int));
    gamma_iter = (int *)malloc((*nsubject) * ntime1 * sizeof(int));
    nh = (int *)malloc((*nsubject) * ntime1 * sizeof(int));
    nclus_iter = (int *)malloc((*ntime + 1) * sizeof(int));

    // Initialize arrays
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
            nh[(Si_iter[j * ntime1 + t] - 1) * ntime1 + t] = nh[(Si_iter[j * ntime1 + t] - 1) * ntime1 + t] + 1;
        }
    }

    // Initialize the number of clusters
    for (t = 0; t < *ntime; t++) {
        nclus_iter[t] = 0;
        for (j = 0; j < *nsubject; j++) {
            if (nh[j * ntime1 + t] > 0) nclus_iter[t] = nclus_iter[t] + 1;
        }
    }
    nclus_iter[*ntime] = 0;

    // Print the results
    printf("Si_iter:\n");
    for (j = 0; j < *nsubject; j++) {
        for (t = 0; t < ntime1; t++) {
            printf("%d ", Si_iter[j * ntime1 + t]);
        }
        printf("\n");
    }

    printf("gamma_iter:\n");
    for (j = 0; j < *nsubject; j++) {
        for (t = 0; t < ntime1; t++) {
            printf("%d ", gamma_iter[j * ntime1 + t]);
        }
        printf("\n");
    }

    printf("nh:\n");
    for (j = 0; j < *nsubject; j++) {
        for (t = 0; t < ntime1; t++) {
            printf("%d ", nh[j * ntime1 + t]);
        }
        printf("\n");
    }

    printf("nclus_iter:\n");
    for (t = 0; t <= *ntime; t++) {
        printf("%d ", nclus_iter[t]);
    }
    printf("\n");

    // Free allocated memory
    free(Si_iter);
    free(gamma_iter);
    free(nh);
    free(nclus_iter);
}

int main() {
    int nsubject = 3;
    int ntime = 4;
    int ntime1 = 5;

    initialize_and_print(&nsubject, &ntime, ntime1);

    return 0;
}
