
// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ 
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include "sha512.hh"
#include <omp.h>
#include <cmath>

using namespace std;

std::string getMinimumPenalties(std::string* genes, int sequenceNum, int mismatchPenalty, int gapPenalty, int* penalties);
int getMinimumPenalty(std::string gene1, std::string gene2, int mismatchPenalty, int gapPenalty, int* gene1Ans, int* gene2Ans);

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
	int sequenceNum;
	std::cin >> misMatchPenalty;
	std::cin >> gapPenalty;
	std::cin >> sequenceNum;
	std::string genes[sequenceNum];

	for (int i = 0; i < sequenceNum; i++)
	{
		std::cin >> genes[i];
	}
	int numPairs = sequenceNum * (sequenceNum - 1) / 2;

	int penalties[numPairs];

	uint64_t startTime = GetTimeStamp();

	// Return all the penalties and the hash of all alignments.
	std::string alignmentHash = getMinimumPenalties(genes,
		sequenceNum, misMatchPenalty, gapPenalty,
		penalties);

	// Print the time taken to do the computation.
	printf("Time: %ld us\n", (uint64_t)(GetTimeStamp() - startTime));

	// Print the alignment hash.
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

// Equivalent of  int *dp[width] = new int[height][width],
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

std::string getMinimumPenalties(std::string* genes, int sequenceNum, int misMatchPenalty, int gapPenalty,
	int* penalties) {
	std::string alignmentHash = "";
	int probNum = 0;

	for (int i = 1; i < sequenceNum; ++i) {
		for (int j = 0; j < i; ++j) {
			std::string gene1 = genes[i];
			std::string gene2 = genes[j];

			int length1 = gene1.length();
			int length2 = gene2.length();
			int totalLength = length1 + length2;

			int gene1Ans[totalLength + 1], gene2Ans[totalLength + 1];

			penalties[probNum] = getMinimumPenalty(gene1, gene2, misMatchPenalty, gapPenalty, gene1Ans, gene2Ans);

			// Since we have assumed the answer to be n+m long,
			// we need to remove the extra gaps in the starting
			// id represents the index from which the arrays
			// xans, yans are useful.
			int id = 1;
			int a;
			for (a = totalLength; a >= 1; --a)
			{
				if ((char) gene1Ans[a] == '_' && (char) gene2Ans[a] == '_')
				{
					id = a + 1;
					break;
				}
			}
			std::string align1 = "";
			std::string align2 = "";
			for (a = id; a <= totalLength; a++)
			{
				align1.append(1, (char) gene1Ans[a]);
			}
			for (a = id; a <= totalLength; ++a)
			{
				align2.append(1, (char) gene2Ans[a]);
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

// Function to find out the minimum penalty.
// Returns the minimum penalty and put the aligned sequences in xans and yans.
int getMinimumPenalty(std::string gene1, std::string gene2, int mismatchPenalty, int gapPenalty, int* gene1Ans, int* gene2Ans)
{
	int i, j;

	int length1 = gene1.length();
	int length2 = gene2.length();

	int rows = length1 + 1;
	int cols = length2 + 1;

	// Table for storing optimal substructure answers.
	int** dp = new2d(length1 + 1, length2 + 1);
	size_t size = length1 + 1;
	size *= length2 + 1;
	memset(dp[0], 0, size);

	// Intialising the table.
	for (i = 0; i <= length1; ++i) {
		dp[i][0] = i * gapPenalty;
	}
	for (i = 0; i <= length2; ++i) {
		dp[0][i] = i * gapPenalty;
	}
	
	/* References for my parallelisation approach: 
	   1. Sequential algorithm, adapted from: 
	        - https://www.tutorialspoint.com/zigzag-or-diagonal-traversal-of-matrix-in-cplusplus.
	   2. Parallel algorithm, idea and approach from: 
			- https://etd.ohiolink.edu/apexprod/rws_etd/send_file/send?accession=kent1429528937&disposition=inline.
			- https://cse.buffalo.edu/~vipin/book_Chapters/2006/2006_2.pdf.
	*/

	int threads = omp_get_max_threads();

	// Dimensions in (matrix) cells of a given block of cells for a given thread, min = 1.
	int blockWidth = (int) ceil((1.0 * length1) / threads);
	int blockLength = (int) ceil((1.0 * length2) / threads);

	// Instead of (rows + cols - 1) like the source material, the number of traversal necessary is 
	// (threads + threads - 1), because threads are moving block by block, not cell by cell.
	for (int traversalNum = 1; traversalNum <= (2 * threads - 1); traversalNum++) {

		// Column index of the starting cell of the current diagonal traversal.
		int startCol = max(1, traversalNum - threads + 1);

		// Number of cells on the current diagonal traversal.
		int cells = min(traversalNum, threads);

		omp_set_dynamic(0);
		omp_set_num_threads(threads);

		#pragma omp parallel for
		for (int currentCol = startCol; currentCol <= cells; currentCol++) {

			// Given the current column index, the block assigned to the current thread
			// consists of cells from (rowStart, colStart) to (colStart, colEnd).
			int rowStart = (currentCol - 1) * blockWidth + 1;
            int colStart = (traversalNum - currentCol) * blockLength + 1;

			// Prevents out-of-bound traversal.
			int rowEnd = min(rowStart + blockWidth, rows);
            int colEnd = min(colStart + blockLength, cols); 

			// Start from cell at (currentRow, currentRol), visit its northwest, north, west neighbours.
			for (int currentRow = rowStart; currentRow < rowEnd; ++currentRow) {

				/*std::cout << "traversalNum: " << traversalNum << endl;
				std::cout << "currentCol: " << currentCol << endl;
				std::cout << "rowStart: " << rowStart << endl;
				std::cout << "rowEnd: " << rowEnd << endl;
				std::cout << "colStart: " << colStart << endl;
				std::cout << "colEnd: " << colEnd << endl;*/

				for (int currentCol = colStart; currentCol < colEnd; ++currentCol) {
					if (gene1[currentRow - 1] == gene2[currentCol - 1]) {
						dp[currentRow][currentCol] = dp[currentRow - 1][currentCol - 1];
					}
					else {
						dp[currentRow][currentCol] = min3(
							dp[currentRow - 1][currentCol - 1] + mismatchPenalty,
							dp[currentRow - 1][currentCol] + gapPenalty,
							dp[currentRow][currentCol - 1] + gapPenalty);
					}
				}
			}
		}
	}

	// Reconstructing the solution.
	int l = length2 + length1; // maximum possible length

	i = length1; 
	j = length2;

	int xPos = l;
	int yPos = l;

	while (!(i == 0 || j == 0)) {
		if (gene1[i - 1] == gene2[j - 1]) {
			gene1Ans[xPos--] = (int) gene1[i - 1];
			gene2Ans[yPos--] = (int) gene2[j - 1];
			i--; j--;
		}

		else if (dp[i - 1][j - 1] + mismatchPenalty == dp[i][j]) {
			gene1Ans[xPos--] = (int) gene1[i - 1];
			gene2Ans[yPos--] = (int) gene2[j - 1];
			i--; j--;
		}

		else if (dp[i - 1][j] + gapPenalty == dp[i][j]) {
			gene1Ans[xPos--] = (int) gene1[i - 1];
			gene2Ans[yPos--] = (int) '_';
			i--;
		}

		else if (dp[i][j - 1] + gapPenalty == dp[i][j]) {
			gene1Ans[xPos--] = (int) '_';
			gene2Ans[yPos--] = (int) gene2[j - 1];
			j--;
		}
	}

	while (xPos > 0) {
		if (i > 0) gene1Ans[xPos--] = (int) gene1[--i];
		else gene1Ans[xPos--] = (int) '_';
	}

	while (yPos > 0) {
		if (j > 0) gene2Ans[yPos--] = (int) gene2[--j];
		else gene2Ans[yPos--] = (int) '_';
	}

	int ret = dp[length1][length2];

	delete[] dp[0];
	delete[] dp;

	return ret;
}

