/*
Generate partially permuted data triples [A1_perm, B1_perm, C1_perm]
 Input:
       M in 4 x 4 x 1
       N in 4 x 4 x n
       P in 4 x 4 x n
       r : scrambling rate
 Output:
       M_perm in 4 x 4 x n (n copies of M)
       N_perm in 4 x 4 x n (permuated N)
      P_perm in 4 x 4 x n (samme as P)
*/

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

void permFixABC(double M[4][4], double N[4][4][100], double P[4][4][100], double r, double M_perm[4][4][100], double N_perm[4][4][100], double P_perm[4][4][100])
{
    int n = 100; // update with the actual size of N
    for (int i = 0; i < n; ++i)
    {
        // Copy M into M_perm
        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < 4; ++k)
            {
                M_perm[j][k][i] = M[j][k];
            }
        }
    }

    // Scramble data in N
    for (int i = 0; i < n; ++i)
    {
        std::vector<int> perm(16);
        for (int j = 0; j < 16; ++j)
        {
            perm[j] = j;
        }
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(perm.begin(), perm.end(), g);
        
        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < 4; ++k)
            {
                int idx = j * 4 + k;
                int new_idx = perm[idx];
                int new_j = new_idx / 4;
                int new_k = new_idx % 4;
                N_perm[j][k][i] = N[new_j][new_k][i];
            }
        }
    }

    // Copy P into P_perm
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < n; ++k)
            {
                P_perm[i][j][k] = P[i][j][k];
            }
        }
    }
}
