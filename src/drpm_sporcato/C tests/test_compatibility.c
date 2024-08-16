#include <stdio.h>
#include <string.h>

// Function prototypes
int compatibility(int *rho1, int *rho2, int n);
void relabel(int *Si, int n, int *Sirelab, int *nhrelab, int *oldLab);

void print_vec(int vec[], int n){
	for (int i=0; i<n; ++i){
		printf("%2d",vec[i]);
	}
	printf("\n");
}
// Function definitions
int compatibility(int *rho1, int *rho2, int n) {
    /*******************************************************************************
     * PURPOSE:
     * Compare rhot | gammat = 1 and rhot_1 | gammat = 1 to determine compatibility
     * I will relabel before comparison so that the first unit is always assigned to
     * cluster 1, the first unit not in cluster 1 will be assigned cluster 2, etc.
     *
     * INPUTS:
     * rho1 - n dimensional integer array holding cluster labels from reduced rho_t
     * rho2 - n dimensional integer array holding cluster labels from reduced rho_t-1
     * n - integer indicating length of rho1 and rho2
     *
     * OUTPUTS:
     * comp - binary indicating if two partitions are comatible
     *      0 - not compatible
     *      1 - compatible
     *******************************************************************************/
    int i, comp = 1, scr1[n], scr2[n], scr3[n], scr4[n];

    relabel(rho1, n, scr1, scr3, scr4);
    relabel(rho2, n, scr2, scr3, scr4);


	printf("original rho1 :"); print_vec(rho1,n);
	printf("relabeled rho1:"); print_vec(scr1,n);
	printf("original rho2 :"); print_vec(rho2,n);
	printf("relabeled rho2:"); print_vec(scr2,n);
    for (i = 0; i < n; i++) {
        if (scr1[i] != scr2[i]) {
            comp = 0;
            break;
        }
    }
    return comp;
}

void relabel(int *Si, int n, int *Sirelab, int *nhrelab, int *oldLab) {
    /*******************************************************************************
     * PURPOSE:
     * Relabel groups so that group with first unit is always 1, then group with
     * closest unit to the first that is not a part of the first group is 2, etc.
     *
     *
     * INPUTS:
     * Si - n dimensional integer array holding cluster labels, this will be d
     * n - integer indicating length of rho1 and rho2
     *
     * OUTPUTS:
     * Sirelab - n dimensional integer that will hold the relabeled group labels
     * nhrelab - n dimensional integer that will contain reorder cluster sizes
     *
     *******************************************************************************/
    int j;
    int shuffle = n, loc = 0, lab = 1;

    for (j = 0; j < n; j++) {
        nhrelab[j] = 0;
        Sirelab[j] = 0;
    }

    while (shuffle > 0) {
        for (j = 0; j < n; j++) {
            if (Si[j] == Si[loc]) {
                Sirelab[j] = lab;
                shuffle = shuffle - 1;
                nhrelab[lab - 1] = nhrelab[lab - 1] + 1;
            }
        }

        oldLab[lab - 1] = Si[loc];

        lab = lab + 1;
        for (j = 0; j < n; j++) {
            if (Sirelab[j] == 0) {
                loc = j;
                break;
            }
        }
    }
}

int main() {
    int rho1[] = {1, 2, 1, 3, 2};
    int rho2[] = {1, 1, 2, 3, 2};
    int rho3[] = {3, 1, 3, 2, 1}; // relabel-equivalent to rho1
    int n = 5;

    int comp1 = compatibility(rho1, rho2, n);
	printf("Compatibility between rho1 and rho2: %d\n", comp1);

    int comp2 = compatibility(rho1, rho3, n);
    printf("Compatibility between rho1 and rho3: %d\n", comp2);

	int rho4[] = {1,3,6};
    int rho5[] = {1,2,4};
    n = 3;
    comp1 = compatibility(rho4, rho5, n);
	printf("Compatibility between rho1 and rho2: %d\n", comp1);

    return 0;
}
