#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

struct answer {
    int u;
    int v;
};

void print_arr(int a[], int n, char s[]) {
    int i;
    printf("%s is\n", s);
    for (i=0; i<n; ++i) {
        printf("%d\n", a[i]);
    }
    printf("------------------------\n");
}

int int_pow(int a, int x) {
    int ret_val = 1;
    int i;
    for (i=0; i<x; ++i) {
        ret_val *= a;
    }
    return ret_val;
}

int rand_range(int min, int max) {
    return min + (int)(rand() * (max - min + 1.0) / (1.0 + RAND_MAX));
}

int min(int a, int b) {
    if (a < b) {
        return a;
    } else {
        return b;
    }
}

struct answer maximum_sub_sequence(int x[], int n) {
    int global_max = x[0];
    int u = 0;
    int v = 0;
    int current_max = x[0];
    int q = 0;
    int i;

    for (int i=0; i<n; ++i) {
        if (current_max >= 0) {
            current_max += x[i];
        } else {
            current_max = x[i];
            q = i;
        }

        if (current_max > global_max) {
            global_max = current_max;
            u = q;
            v = i;
        }
    }

    struct answer ret_val;
    ret_val.u = u;
    ret_val.v = v;

    return ret_val;
}

void parallel_prefix_sum(int x[], int p[], int n) {
    int i;
    int j;
    int k;
    int exp_2_k;
    int start;

    #pragma omp parallel for
    for (i=0; i<n; ++i) {
        p[i] = x[i];
    }

    for (k=0; k<=(int)log2((double)n); k++) {
        exp_2_k = int_pow(2, k);
        #pragma omp parallel for private(start)
        for (j=0; j < (n - exp_2_k) / (2 * exp_2_k) + 1; ++j) {
            start = exp_2_k + 2 * j * exp_2_k;
            #pragma omp parallel for
            for (i=0; i<min(exp_2_k, n - start); ++i) {
                p[start + i] += p[start - 1];
            }
        }
    }
}

void parallel_prefix_sum2(int x[], int p[], int n) {
    int log2n = (int)log2((double)n) + 1;
    int i,j;
    int *x_sub, *p_sub;
    x_sub = (int *)malloc((n/log2n+1) * sizeof(int));
    p_sub = (int *)malloc((n/log2n+1) * sizeof(int));

    #pragma omp parallel for
    for (i=0; i<n; ++i) {
        p[i] = x[i];
    }

    #pragma omp parallel for
    for (i=0; i<n/log2n+1; ++i) {
        #pragma omp parallel for ordered
        for (j=i*log2n; j<min((i+1)*log2n-1, n-1); ++j) {
            #pragma omp ordered
            p[j+1] += p[j];
        }
    }

    #pragma omp parallel for
    for (i=0; i<n/log2n+1; ++i) {
        x_sub[i] = p[min((i+1)*log2n-1, n-1)];
    }

    parallel_prefix_sum(x_sub, p_sub, n/log2n+1);

    #pragma omp parallel for
    for (i=0; i<n/log2n+1; ++i) {
        p[min((i+1)*log2n-1, n-1)] = p_sub[i];
    }

    #pragma omp parallel for
    for (i=1; i<n/log2n+1; ++i) {
        #pragma omp parallel for
        for (j=i*log2n; j<min((i+1)*log2n-1, n-1); ++j) {
            p[j] += p[i*log2n-1];
        }
    }
}

void parallel_postfix_maximum(int p[], int m[], int a[], int n) {
    int i;
    int j;
    int k;
    int tmp, tmp_;
    int exp_2_k, start;

    #pragma omp parallel for
    for (i=0; i<n; ++i) {
        m[i] = p[n - i - 1];
        a[i] = n - i - 1;
    }

    for (k=0; k<=(int)log2((double)n); k++) {
        exp_2_k = int_pow(2, k);
        #pragma omp parallel for private(start)
        for (j=0; j < (n - exp_2_k) / (2 * exp_2_k) + 1; ++j) {
            start = exp_2_k + 2 * j * exp_2_k;
            #pragma omp parallel for
            for (i=0; i<min(exp_2_k, n - start); ++i) {
                if (m[start + i] < m[start - 1]) {
                    m[start + i] = m[start - 1];
                    a[start + i] = a[start - 1];
                }
            }
        }
    }

    #pragma omp parallel for
    for (i=0; i<n/2 ; ++i) {
        {
            tmp = m[i];
            m[i] = m[n - i - 1];
            m[n - i - 1] = tmp;
        }

        {
            tmp_ = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = tmp_;
        }
    }
}

void parallel_prefix_maximum(int p[], int q[], int m[], int a[], int n) {
    int i;
    int j;
    int k;
    int exp_2_k, start;

    #pragma omp parallel for
    for (i=0; i<n; ++i) {
        m[i] = p[i];
        a[i] = q[i];
    }

    for (k=0; k<=(int)log2((double)n); k++) {
        exp_2_k = int_pow(2, k);
        #pragma omp parallel for private(start)
        for (j=0; j < (n - exp_2_k) / (2 * exp_2_k) + 1; ++j) {
            start = exp_2_k + 2 * j * exp_2_k;
            #pragma omp parallel for
            for (i=0; i<min(exp_2_k, n - start); ++i) {
                if (m[start + i] < m[start - 1]) {
                    m[start + i] = m[start - 1];
                    a[start + i] = a[start - 1];
                }
            }
        }
    }
}

void parallel_postfix_maximum2(int p[], int m[], int a[], int n) {
    int log2n = (int)log2((double)n) + 1;
    int i, j, tmp, tmp_;
    int *p_sub, *q_sub, *m_sub, *a_sub;
    p_sub = (int *)malloc((n/log2n+1) * sizeof(int));
    m_sub = (int *)malloc((n/log2n+1) * sizeof(int));
    a_sub = (int *)malloc((n/log2n+1) * sizeof(int));
    q_sub = (int *)malloc((n/log2n+1) * sizeof(int));

    #pragma omp parallel for
    for (i=0; i<n; ++i) {
        m[i] = p[n - i - 1];
        a[i] = n - i - 1;
    }

    #pragma omp parallel for
    for (i=0; i<n/log2n+1; ++i) {
        #pragma omp parallel for ordered
        for (j=i*log2n; j<min((i+1)*log2n-1, n-1); ++j) {
            #pragma omp ordered
            if (m[j+1] < m[j]) {
                m[j+1] = m[j];
                a[j+1] = a[j];
            }
        }
    }

    #pragma omp parallel for
    for (i=0; i<n/log2n+1; ++i) {
        p_sub[i] = m[min((i+1)*log2n-1, n-1)];
        q_sub[i] = a[min((i+1)*log2n-1, n-1)];
    }

    parallel_prefix_maximum(p_sub, q_sub, m_sub, a_sub, n/log2n+1);

    #pragma omp parallel for
    for (i=0; i<n/log2n+1; ++i) {
        m[min((i+1)*log2n-1, n-1)] = m_sub[i];
        a[min((i+1)*log2n-1, n-1)] = a_sub[i];
    }

    #pragma omp parallel for
    for (i=1; i<n/log2n+1; ++i) {
        #pragma omp parallel for
        for (j=i*log2n; j<min((i+1)*log2n-1, n-1); ++j) {
            if (m[j] < m[i*log2n-1]) {
                m[j] = m[i*log2n-1];
                a[j] = a[i*log2n-1];
            }
        }
    }

    #pragma omp parallel for
    for (i=0; i<n/2 ; ++i) {
        {
            tmp = m[i];
            m[i] = m[n - i - 1];
            m[n - i - 1] = tmp;
        }

        {
            tmp_ = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = tmp_;
        }
    }
}

struct answer maximum_sub_sequence_parallel(int x[], int n) {
    int *p;
    p = (int *)malloc(n * sizeof(int));
    int *m;
    m = (int *)malloc(n * sizeof(int));
    int *a;
    a = (int *)malloc(n * sizeof(int));
    int u;
    int *mss;
    mss = (int *)malloc(n * sizeof(int));
    int u_star;
    int max_val;
    struct answer ret_val;

    parallel_prefix_sum2(x, p, n);
    parallel_postfix_maximum2(p, m, a, n);

    #pragma omp parallel for
    for (u=0; u<n; ++u) {
        mss[u] = m[u] - p[u] + x[u];
    }

    max_val = mss[0];
    u_star = 0;

    #pragma omp parallel for
    for (u=1; u<n; ++u) {
        #pragma omp flush(max_val)
        if (mss[u] > max_val) {
            #pragma omp critical
            {
                if (mss[u] > max_val) {
                    max_val = mss[u];
                    u_star = u;
                }
            }
        }
    }

    ret_val.u = u_star;
    ret_val.v = a[u_star];
    free(p);
    free(m);
    free(a);
    free(mss);
   return ret_val;
}

int main(int argc, char *argv[]) {

    if (argc != 2) {
        fprintf(stderr, "Usage: program <maximum exponent of # of integers to sort>\n");
        exit(1);
    }

    int nthread;
    int mthread = omp_get_num_procs();
    int max_exponent = atoi(argv[1]);
    //int seed = atoi(argv[3]);
    int i;
    int min = -5;
    int max = 5;
    FILE *fp;

    srand(10);
    fp = fopen("results.csv", "w");
    if( fp == NULL ){
        printf( "%s cannot be opened.¥n", "sequential.csv");
        return -1;
    }
    fprintf(fp, "nthread,size,seq_time,par_time");

    for (nthread=1; nthread<=mthread; ++nthread){
    omp_set_num_threads(nthread);
    printf("------------------------------\nnthread = %d\n", nthread);
    for (i=1; i<=max_exponent; ++i) {
        int na = int_pow(10, i);
        int *a;
        a = (int *)malloc(na * sizeof(int));
        if (a == 0) {
            fprintf(stderr, "malloc error\n");
            exit(1);
        }
        int j;
        struct answer ret_seq, ret_par;
        double start_seq, end_seq, start_par, end_par;

        for (j=0; j<na; ++j) {
            a[j] = rand_range(min, max);
        }

        if (nthread == 1 && i == 1) {
            print_arr(a, na, "a");
        }

        start_seq = omp_get_wtime();
        ret_seq = maximum_sub_sequence(a, na);
        end_seq = omp_get_wtime();

        start_par = omp_get_wtime();
        ret_par = maximum_sub_sequence_parallel(a, na);
        end_par = omp_get_wtime();

        fprintf(fp, "%d,%d,%f,%f\n", nthread, na, (end_seq-start_seq)*1000.0, (end_par-start_par)*1000.0);
        fprintf(stdout, "(u, v) = (%d, %d)\n", ret_seq.u, ret_seq.v);
        fprintf(stdout, "(u, v) = (%d, %d)\n", ret_par.u, ret_par.v);

        int sum_seq=0, sum_par=0, hoge;
        for (hoge=ret_seq.u; hoge<ret_seq.v+1; ++hoge) {
            sum_seq += a[hoge];
        }
        for (hoge=ret_par.u; hoge<ret_par.v+1; ++hoge) {
            sum_par += a[hoge];
        }
        printf("(seq, par) = (%d, %d)\n", sum_seq, sum_par);

        printf("------------------------\n");
        free(a);
    }
    }
    fclose(fp);

    return 0;
}
