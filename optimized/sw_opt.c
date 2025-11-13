#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <immintrin.h>

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MATCH     2
#define MISMATCH -1
#define GAP      -2

void generate_sequence(char *seq, int n) {
    const char alphabet[] = "ACGT";
    for (int i = 0; i < n; i++)
        seq[i] = alphabet[rand() % 4];
    seq[n] = '\0';
}

static inline int sw_cell(char a, char b, int diag, int up, int left) {
    int match = diag + (a == b ? MATCH : MISMATCH);
    int del = up + GAP;
    int ins = left + GAP;
    return MAX(0, MAX(match, MAX(del, ins)));
}

int smith_waterman_optimized(const char *seq1, const char *seq2, int len1, int len2) {
    // Use two rows only (previous and current) to reduce memory and improve cache behavior.
    int *prev = calloc(len2 + 1, sizeof(int));
    int *curr = calloc(len2 + 1, sizeof(int));
    if (!prev || !curr) {
        fprintf(stderr, "Allocation failed\n");
        free(prev); free(curr);
        return 0;
    }

    int max_score = 0;

    for (int i = 1; i <= len1; ++i) {
        const char a = seq1[i-1];
        // curr[0] is always 0 for Smith-Waterman
        curr[0] = 0;
        for (int j = 1; j <= len2; ++j) {
            int diag = prev[j-1] + (a == seq2[j-1] ? MATCH : MISMATCH);
            int up   = prev[j]   + GAP;
            int left = curr[j-1] + GAP;

            int t = diag;
            if (up > t) t = up;
            if (left > t) t = left;
            if (t < 0) t = 0;

            curr[j] = t;
            if (t > max_score) max_score = t;
        }
        // swap prev and curr pointers (avoid copying)
        int *tmp = prev;
        prev = curr;
        curr = tmp;
    }

    free(prev);
    free(curr);
    return max_score;
}


int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <sequence_length>\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    srand(42);

    char *seq1 = malloc((N + 1) * sizeof(char));
    char *seq2 = malloc((N + 1) * sizeof(char));

    generate_sequence(seq1, N);
    generate_sequence(seq2, N);

    clock_t start = clock();
    int score = smith_waterman_optimized(seq1, seq2, N, N);
    clock_t end = clock();

    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Sequence length: %d\n", N);
    printf("Smith-Waterman optimized score: %d\n", score);
    printf("Execution time: %.6f seconds\n", elapsed);

    free(seq1);
    free(seq2);
    return 0;
}
