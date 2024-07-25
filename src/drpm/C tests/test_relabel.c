#include <stdio.h>
#include <stdlib.h>

void relabel(int *Si, int n, int *Sirelab, int *reorder, int *oldLab){
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
 * reorder - n dimensional integer that will contain reorder cluster sizes
 *
 *******************************************************************************/
	int j;
	int shuffle = n, loc = 0, lab = 1;

	for(j = 0; j < n; j++){
		reorder[j] = 0, Sirelab[j] = 0;
	}

	while(shuffle > 0){
		for(j = 0; j < n; j++){
			if(Si[j] == Si[loc]){
				Sirelab[j] = lab;
				shuffle = shuffle - 1;
				reorder[lab-1] = reorder[lab-1] + 1;
			}
		}

		oldLab[lab-1] = Si[loc];

		lab = lab+1;
		for(j = 0; j < n; j++){
			if(Sirelab[j] == 0){
				loc = j;
				break;
			}
		}
	}
}

int main() {
    int n = 10;
    // int Si[10] = {1,3,3,5,5,4,4,1,2,1};
    // int Si[10] = {1,1,1,1,1,3,1,1,2,3};
    int Si[10] = {4,2,1,3,1,4,5,0,0,0};
  	int Sirelab[10] = {0};
    int reorder[10] = {0};
    int oldLab[10]  = {0};

    relabel(Si, n, Sirelab, reorder, oldLab);

    printf("             Original labels (Si): ");
    for (int i = 0; i < n; i++) {
        printf("%d ", Si[i]);
    }
    printf("\n");

    printf("       Relabeled groups (Sirelab): ");
    for (int i = 0; i < n; i++) {
        printf("%d ", Sirelab[i]);
    }
    printf("\n");

    printf("Reordered cluster sizes (reorder): ");
    for (int i = 0; i < n; i++) {
        printf("%d ", reorder[i]);
    }
    printf("\n");

    printf("              Old labels (oldLab): ");
    for (int i = 0; i < n; i++) {
        printf("%d ", oldLab[i]);
    }
    printf("\n");

    return 0;
}