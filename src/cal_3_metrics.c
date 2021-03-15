#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/***************
 Calculate indices to measure the agreement between two partitions.

 This code is originially from the paper below and written by the authors. I adopt the
 calculation for ARI and Jaccard. Jaccard is usually computed for two set of lists, however,
 in this implementation, it can be generalized in different sets of lists.

@article{Milligan:1986,
  author="{Milligan, G. W.} and {Cooper, M. C.}",
  title={A study of the comparability of external criteria for hierarchical
cluster analysis},
  journal={Multivariate Behavioral Research},
  volume={21},
  pages={441-458},
  year={1986}
}
***************/

/***************
  cl1 --- partition 1 of the data set. 'cl1' is a 'n' by 1 vector.
  cl1u --- unique values of elements in 'cl1'. 'cl1u' is a 'm1' by 1 vector.
  cl2 --- partition 2 of the data set. 'cl2' is a 'n' by 1 vector.
  cl2u --- unique values of elements in 'cl2'. 'cl2u' is a 'm2' by 1 vector.
  m1 --- number of clusters in partition 1
  m2 --- number of clusters in partition 2 (m2 can be not equal to m1)
  n --- number of data points
  flag = 1 --- Rand index
  flag = 2 --- Fowlkes and Mallows's index
  flag = 3 --- Jaccard index
  r12 --- the output index
***************/
void cal_3_metrics(int *cl1, int *cl1u, int *cl2, int *cl2u, int *m1, int *m2,
		  int *n, int *flag, double *r12)
{
    int i, j, t, r, *nmatrix;
    int mm1, mm2, nn, fflag;
    double a, b, c, d, denom;
    double *nc, *nr, ni_2, n_j2, nt, nij_2;

    mm1 = *m1; mm2 = *m2; nn = *n; fflag = *flag;

    nmatrix = (int *)malloc((size_t)(mm1 * mm2 * sizeof(int)));
    nc = (double *)malloc((size_t)(mm2 * sizeof(double)));
    nr = (double *)malloc((size_t)(mm1 * sizeof(double)));

    a = 0.0; b = 0.0; c = 0.0; d = 0.0;
    denom = 0.0;
    for(t = 0; t < nn ; t ++){
        for(r = t + 1; r < nn; r ++){
            if((cl1[t] == cl1[r]) && (cl2[t] == cl2[r])){
                a = a + 1.0;
            } else if((cl1[t] == cl1[r]) && (cl2[t] != cl2[r])){
                b = b + 1.0;
            } else if((cl1[t] != cl1[r]) && (cl2[t] == cl2[r])){
                c = c + 1.0;
            } else{
                d = d + 1.0;
            }
        }
    }
    // get nij
    for(t = 0; t < mm1; t ++){
        for(r = 0; r < mm2; r ++){
            nmatrix[t * mm2 + r] = 0;
            for(i = 0; i < nn; i ++){
                if((cl1[i] == cl1u[t]) && (cl2[i] == cl2u[r])){
                    nmatrix[t * mm2 + r] += 1;
                }
            }
        }
    }

    nij_2 = 0.0;
    for(i = 0; i < mm1; i ++){
        nr[i] = 0;
        for(j = 0; j < mm2; j ++){
            nr[i] += nmatrix[i * mm2 + j];
            nij_2 += pow(nmatrix[i * mm2 + j],2);
        }
    }

    for(i = 0; i < mm2; i ++){
        nc[i] = 0;
        for(j = 0; j < mm1; j ++){
            nc[i] += nmatrix[j * mm2 + i];
        }
    }

    ni_2 = 0.0; n_j2 = 0.0; nt = 0.0;
    for(i = 0; i < mm1; i ++){
        nt += nr[i];
        ni_2 += nr[i] * nr[i];
    }
    for(j = 0; j < mm2; j ++){
        n_j2 += nc[j] * nc[j];
    }

    if(fflag == 1){// Rand
        *r12 = (a + d) / (a + b + c + d);
    } else if(fflag == 2) {  //Fowlkes and Mallows
        *r12 = a / sqrt((a + b) * (a + c));
    } else if(fflag == 3) { //Jaccard
        *r12 = a / (a + b + c);
    }

    return;
}






