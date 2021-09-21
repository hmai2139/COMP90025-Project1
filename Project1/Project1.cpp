
// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ 
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include "sha512.hh"
#include <omp.h>
#include <numa.h>
#include <cmath>

using namespace std;

std::string getMinimumPenalties(std::string* genes, int k, int pxy, int pgap, int* penalties);
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int* xans, int* yans);

/*
Examples of sha512 which returns a std::string
sw::sha512::calculate("SHA512 of std::string") // hash of a string, or
sw::sha512::file(path) // hash of a file specified by its path, or
sw::sha512::calculate(&data, sizeof(data)) // hash of any block of data
*/

// Return current wallclock time, for performance measurement
uint64_t GetTimeStamp() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
}

int main(int argc, char** argv) {
	int misMatchPenalty;
	int gapPenalty;
	int k;
	std::cin >> misMatchPenalty;
	std::cin >> gapPenalty;
	std::cin >> k;
	std::string genes[k];

	for (int i = 0; i < k; i++)
	{
		std::cin >> genes[i];
	}
	int numPairs = k * (k - 1) / 2;

	int penalties[numPairs];

	uint64_t start = GetTimeStamp();

	// return all the penalties and the hash of all allignments
	std::string alignmentHash = getMinimumPenalties(genes,
		k, misMatchPenalty, gapPenalty,
		penalties);

	// print the time taken to do the computation
	printf("Time: %ld us\n", (uint64_t)(GetTimeStamp() - start));

	// print the alginment hash
	std::cout << alignmentHash << std::endl;

	for (int i = 0; i < numPairs; i++) {
		std::cout << penalties[i] << " ";
	}
	std::cout << std::endl;
	return 0;
}

int min3(int a, int b, int c) {
	if (a <= b && a <= c) {
		return a;
	}
	else if (b <= a && b <= c) {
		return b;
	}
	else {
		return c;
	}
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
int** new2d(int width, int height) {
	int** dp = new int* [width];
	size_t size = width;
	size *= height;
	int* dp0 = new int[size];
	if (!dp || !dp0) {
		std::cerr << "getMinimumPenalty: new failed" << std::endl;
		exit(1);
	}
	dp[0] = dp0;
	for (int i = 1; i < width; i++)
		dp[i] = dp[i - 1] + height;

	return dp;
}

std::string getMinimumPenalties(std::string* genes, int k, int pxy, int pgap,
	int* penalties) {

	std::string alignmentHash = "";
	int probNum = 0;

	for (int i = 1; i < k; ++i) {
		for (int j = 0; j < i; ++j) {
			std::string gene1 = genes[i];
			std::string gene2 = genes[j];
			int m = gene1.length();
			int n = gene2.length();
			int l = m + n;
			int xans[l + 1], yans[l + 1];

			penalties[probNum] = getMinimumPenalty(gene1, gene2, pxy, pgap, xans, yans);

			// Since we have assumed the answer to be n+m long,
			// we need to remove the extra gaps in the starting
			// id represents the index from which the arrays
			// xans, yans are useful
			int id = 1;
			int a;
			for (a = l; a >= 1; --a)
			{
				if ((char)yans[a] == '_' && (char)xans[a] == '_')
				{
					id = a + 1;
					break;
				}
			}
			std::string align1 = "";
			std::string align2 = "";
			for (a = id; a <= l; a++)
			{
				align1.append(1, (char)xans[a]);
			}
			for (a = id; a <= l; ++a)
			{
				align2.append(1, (char)yans[a]);
			}

			std::string align1hash = sw::sha512::calculate(align1);
			std::string align2hash = sw::sha512::calculate(align2);
			std::string problemhash = sw::sha512::calculate(align1hash.append(align2hash));

			alignmentHash = sw::sha512::calculate(alignmentHash.append(problemhash));
			probNum++;
		}
	}
	return alignmentHash;
}
bool isValid(int i, int j, int m, int n)
{
	if (i < 0 || i >= m
		|| j >= n || j < 0)
		return false;
	return true;
}
// Function to find out the minimum penalty
// Returns the minimum penalty and put the aligned sequences in xans and yans
int getMinimumPenalty(std::string gene1, std::string gene2, int pGene1Gene2, int pgap, int* xans, int* yans)
{
	int i, j;

	int m = gene1.length();
	int n = gene2.length();

	int rows = m + 1;
	int cols = n + 1;

	// Table for storing optimal substructure answers
	int** dp = new2d(m + 1, n + 1);
	size_t size = m + 1;
	size *= n + 1;
	memset(dp[0], 0, size);

	// Intialising the table
	for (i = 0; i <= m; ++i) {
		dp[i][0] = i * pgap;
	}
	for (i = 0; i <= n; ++i) {
		dp[0][i] = i * pgap;
	}
	
	// Reference: https://www.tutorialspoint.com/zigzag-or-diagonal-traversal-of-matrix-in-cplusplus

	int threads = omp_get_max_threads();

	// Dimensions of a given block of cells for a given thread.
	int blockWidth = (int) ceil((1.0 * m) / threads);
	int blockLength = (int) ceil((1.0 * n) / threads);

	for (int traversalNum = 1; traversalNum <= (2 * threads - 1); traversalNum++)
	{
		// Column index of the starting cell of the current diagonal traversal.
		int startCol = max(1, traversalNum - threads + 1);

		// Number of cells on the current diagonal traversal.
		int cellNum = min(traversalNum, threads);

		#pragma omp parallel for
		for (int currentCol = startCol; currentCol <= cellNum; currentCol++) {

			// Given the current column index, the block assigned to the current thread are defined by the following indices.
			int rowStart = (currentCol - 1) * blockWidth + 1;
            int rowEnd = min(rowStart + blockWidth, rows); // Prevents out of bound.
            int colStart = (traversalNum - currentCol) * blockLength + 1;
            int colEnd = min(colStart + blockLength, cols); // Prevents out of bound.
			for (int i = rowStart; i < rowEnd; ++i) {
				for (int j = colStart; j < colEnd; ++j) {
					if (gene1[i - 1] == gene2[j - 1]) {
						dp[i][j] = dp[i - 1][j - 1];
					}
					else {
						dp[i][j] = min(min(dp[i - 1][j - 1] + pGene1Gene2,
							dp[i - 1][j] + pgap),
							dp[i][j - 1] + pgap);
					}
				}
			}
		}
	}

	// Reconstructing the solution.
	int l = n + m; // maximum possible length

	i = m; 
	j = n;

	int xpos = l;
	int ypos = l;

	while (!(i == 0 || j == 0))
	{
		if (gene1[i - 1] == gene2[j - 1])
		{
			xans[xpos--] = (int)gene1[i - 1];
			yans[ypos--] = (int)gene2[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j - 1] + pGene1Gene2 == dp[i][j])
		{
			xans[xpos--] = (int)gene1[i - 1];
			yans[ypos--] = (int)gene2[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)gene1[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if (dp[i][j - 1] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)'_';
			yans[ypos--] = (int)gene2[j - 1];
			j--;
		}
	}
	while (xpos > 0)
	{
		if (i > 0) xans[xpos--] = (int)gene1[--i];
		else xans[xpos--] = (int)'_';
	}
	while (ypos > 0)
	{
		if (j > 0) yans[ypos--] = (int)gene2[--j];
		else yans[ypos--] = (int)'_';
	}

	int ret = dp[m][n];

	delete[] dp[0];
	delete[] dp;

	return ret;
}

