#include "stdafx.h"
#include<limits.h>
#include<stdio.h>
#include<list>
#include <iostream>
#include<stdlib.h>
#include <time.h>
#include <vector>

using namespace std;

#define LINUX

#define DEBUG 0

typedef unsigned long int  bitvector_t;
typedef bitvector_t        *cs_t;

typedef enum { identity_b, complement_b }  cmode_t;

typedef struct row_s {
	long  originalRowNumber;
	cs_t  columnSet;
} row_t;

int          bitsPerBV;
int          noBVs;
bitvector_t  bitMaskLastBV;
int          biclusterCount;
long   noRows;
long   noColumns;
long   minNoRows;
long   minNoColumns;
long   maxLevels;
row_t  *rows;
cs_t   *consideredColumns;
cs_t   *mandatoryColumns;
cs_t   columnIntersection;
std::vector<std::vector<int> > store;
int maxedged = 0;
int maxnotfound = 0;

// definitions
struct Coord {
	int x, y;
};


int  newrn, newln, newm, newn, maxln, maxrn, maxn, maxm;
int n, m, ln, rn;
int       maxobjval = 0;
list<Coord> clist; // coordinates list 


int  isSet(cs_t  columnSet, long  column)
{
	bitvector_t  bv;

	if (column >= 0L && column < noColumns) {
		bv = 1U << (column % bitsPerBV);
		return ((columnSet[column / bitsPerBV] & bv) != 0);
	}
	return 0;
} /* setColumn */

void  setColumn(cs_t  columnSet, long  column)
{
	bitvector_t  bv;

	if (column >= 0L && column < noColumns) {
		bv = 1U << (column % bitsPerBV);
		columnSet[column / bitsPerBV] |= bv;
	}
} /* setColumn */

void  unsetColumn(cs_t  columnSet, long  column)
{
	bitvector_t  bv;

	if (column >= 0L && column < noColumns) {
		bv = ~(~(columnSet[column / bitsPerBV]) | (1U << (column % bitsPerBV)));
		columnSet[column / bitsPerBV] &= bv;
	}
} /* unsetColumn */

long  columnCount(cs_t  columnSet)
{
	long         i, j;
	long         counter;
	bitvector_t  bv;

	columnSet[noBVs - 1] &= bitMaskLastBV;
	counter = 0L;
	for (i = noBVs - 1; i >= 0; i--) {
		bv = columnSet[i];
		if (bv != 0U) {
			for (j = 0L; j < bitsPerBV; j++) {
				if (bv & 1U)  counter++;
				bv >>= 1;
			}
		}
	}
	return counter;
} /* columnCount */

int  compareColumns(cs_t  columnSetA, cs_t  columnSetB, cs_t  mask)
{
	int          i;
	int          contained, disjoint;
	bitvector_t  bitMask, sharedColumns;

	contained = disjoint = 1;
	bitMask = bitMaskLastBV;
	for (i = noBVs - 1; i >= 0; i--) {
		sharedColumns = ((columnSetA[i] & columnSetB[i]) & mask[i]) & bitMask;
		if ((sharedColumns | columnSetB[i]) != sharedColumns)
			contained = 0;
		if (sharedColumns != 0)
			disjoint = 0;
		bitMask = ~0U;
	}
	if (contained && disjoint)
		return -2; /* either set is empty */
	if (contained)
		return -1; /* set A contains set B */
	if (disjoint)
		return 1; /* set A and set B are disjoint */
	return 0; /* set B is larger than set A and the intersection is not empty */
} /* compareColumns */

void  copyColumnSet(cs_t  columnSet, cs_t  columnSetDest, cmode_t  copyMode)
{
	int  i;

	for (i = noBVs - 1; i >= 0; i--)
		if (copyMode == complement_b)
			columnSetDest[i] = ~columnSet[i];
		else
			columnSetDest[i] = columnSet[i];
} /* copyColumnSet */

void  intersectColumnSets(cs_t  columnSetA, cs_t  columnSetB, cs_t  columnSetDest)
{
	int  i;

	for (i = noBVs - 1; i >= 0; i--)
		columnSetDest[i] = (columnSetA[i] & columnSetB[i]);
} /* intersectColumnSets */

void  determineColumnsInCommon(long  firstRow, long  lastRow, cs_t  sharedColumnSet)
{
	int   i;
	long  j;

	if (firstRow >= 0L && lastRow >= firstRow && lastRow < noRows) {
		for (i = noBVs - 1; i >= 0; i--) {
			sharedColumnSet[i] = ~0U;
			for (j = firstRow; j <= lastRow; j++)
				sharedColumnSet[i] &= rows[j].columnSet[i];
		}
	}
} /* determineColumnsInCommon */

int  containsMandatoryColumns(cs_t  columnSet, int  noSets)
{
	int   contains, j;
	long  i;

	contains = 1;
	for (i = 0; i < noSets; i++) {
		if ((mandatoryColumns[i][noBVs - 1] & columnSet[noBVs - 1] & bitMaskLastBV) == 0U) {
			j = noBVs - 2;
			while (j >= 0 && (mandatoryColumns[i][j] & columnSet[j]) == 0U)
				j--;
			if (j < 0) {
				contains = 0;
				i = noSets;
			}
		}
	}
	return contains;
} /* containsMandatoryColumns */

void  swapRows(long  a, long  b)
{
	long   tempOriginalRowNumber;
	cs_t  tempColumnSet;

	if (a != b && a >= 0L && a < noRows && b >= 0L && b < noRows) {
		tempOriginalRowNumber = rows[a].originalRowNumber;
		tempColumnSet = rows[a].columnSet;
		rows[a].originalRowNumber = rows[b].originalRowNumber;
		rows[a].columnSet = rows[b].columnSet;
		rows[b].originalRowNumber = tempOriginalRowNumber;
		rows[b].columnSet = tempColumnSet;
	}
} /* swapRows */

long  chooseSplitRow(long  firstRow, long  lastRow, int  level)
{
	long  i;

	for (i = firstRow; i <= lastRow &&
		compareColumns(rows[i].columnSet, consideredColumns[level],
			consideredColumns[0]) < 0; i++);
	if (i <= lastRow)
		return i;
	return firstRow;
} /* chooseSplitRow */

long  selectRows(long  firstRow, long  lastRow, long  level, int  *overlapping)
{
	long  selected;

	selected = 0L;
	*overlapping = 0;
	while (firstRow <= lastRow) {
		switch (compareColumns(consideredColumns[level], rows[firstRow].columnSet,
			consideredColumns[level - 1L])) {
		case -2:
		case 1:
			swapRows(lastRow, firstRow);
			lastRow--;
			break;
		case 0:
			*overlapping = 1;
		default:
			selected++;
			firstRow++;
			break;
		}
	}
	return selected;
} /* selectRows */

void  writeBicluster(long  firstRow, long  lastRow, cs_t  columnSet)
{
	long  biclusterCounter = 0;
	long  i;
	FILE *fptr;
	bool flag = false;
	long  j;
	/* to check bicluster content */
	for (i = firstRow; i <= lastRow; i++)
	{
		for (j = 0L; j < noColumns; j++) {
			if (isSet(columnSet, j)) {
				if (store[rows[i].originalRowNumber][j]) {
					flag = true;
					break;
				}
			}
		}
		if (flag == true)
			break;
	}
	flag = false;

	if (flag == false) {

		biclusterCounter++;
		//cout << "Bicluster " << biclusterCount << " is calculated for BIMAX\n";
		//   printf("\n%ld\n", biclusterCounter);
		int count = 0;
		for (i = 0; i < noColumns; i++)
			if (isSet(columnSet, i))
				count++;

		if (maxedged < count * (lastRow - firstRow + 1)) {
			cout << " Found " << lastRow - firstRow + 1 << " x " << count << endl;
			maxedged = count * (lastRow - firstRow + 1);
			fopen_s(&fptr, "BIMAXResult.txt", "w");
			//cout << "K:" << lastRow - firstRow + 1 << " - " << "L:" << count << endl;
			fprintf(fptr, "%d\t%d\n", lastRow - firstRow + 1, count);
			for (i = firstRow; i <= lastRow; i++) {
				fprintf(fptr, "%d\t", rows[i].originalRowNumber + 1L);
				//     printf("%ld\t", rows[i].originalRowNumber + 1L);
			}
			//   printf("\n");
			fprintf(fptr, "\n");
			for (i = 0; i < noColumns; i++)
				if (isSet(columnSet, i)) {
					//       printf("%ld\t", i + 1L);
					fprintf(fptr, "%d\t", i);
				}
			fprintf(fptr, "\n");
			fclose(fptr);


			// Remove following block not to show bicluster
			//#ifndef DEBUG
			long  j;
			for (i = firstRow; i <= lastRow; i++)
			{
				cout << i << ":";
				for (j = 0L; j < noColumns; j++) {
					if (isSet(columnSet, j))
						cout << store[rows[i].originalRowNumber][j];
				}
				cout << endl;
			}
			//#endif
		}

		// rows[i].columnSet


		biclusterCount--;
	}
} /* writeBicluster */

bool conquer(long  firstRow, long  lastRow, long  level, long noMandatorySets)
{
	int   overlapping = 0;
	long  splitRow, noSelectedRows;
	if (biclusterCount != 0) {
		determineColumnsInCommon(firstRow, lastRow, columnIntersection);
		if (compareColumns(columnIntersection, consideredColumns[level],
			consideredColumns[level]) == -1) {
			writeBicluster(firstRow, lastRow, columnIntersection);
		}
		else {
			splitRow = chooseSplitRow(firstRow, lastRow, level);
			intersectColumnSets(consideredColumns[level], rows[splitRow].columnSet,
				consideredColumns[level + 1L]);
			if (columnCount(consideredColumns[level + 1L]) >= minNoColumns &&
				containsMandatoryColumns(consideredColumns[level + 1L], noMandatorySets)) {
				noSelectedRows = selectRows(firstRow, lastRow, level + 1L, &overlapping);
				if (noSelectedRows >= minNoRows)
					conquer(firstRow, firstRow + noSelectedRows - 1L, level + 1L, noMandatorySets);
			}
			copyColumnSet(consideredColumns[level + 1L], consideredColumns[level + 1L],
				complement_b);
			intersectColumnSets(consideredColumns[level], consideredColumns[level + 1L],
				consideredColumns[level + 1L]);
			if (overlapping) {
				copyColumnSet(consideredColumns[level + 1L], mandatoryColumns[noMandatorySets],
					identity_b);
				noMandatorySets++;
			}
			noSelectedRows = selectRows(firstRow, lastRow, level + 1L, &overlapping);
			copyColumnSet(consideredColumns[level], consideredColumns[level + 1L], identity_b);
			if (noSelectedRows >= minNoRows)
				conquer(firstRow, firstRow + noSelectedRows - 1L, level + 1L, noMandatorySets);
		}
	}

	if (minNoRows * minNoColumns <= maxedged) {
		return true;
	}
	else {
		return false;
	}
} /* conquer */

int  initialize()
{
	bitvector_t  dummy;
	int          failed;
	long         i;

	/* initilization for handling bit vectors */
	dummy = 1;
	bitsPerBV = 0;
	while (dummy != 0) {
		dummy <<= 1;
		bitsPerBV++;
	}
	bitMaskLastBV = (~0U >> (bitsPerBV - (noColumns % bitsPerBV)));
	noBVs = (noColumns / bitsPerBV) + ((noColumns % bitsPerBV) == 0 ? 0 : 1);

	/* memory allocation */
	failed = 0;
	rows = (row_t*)malloc(sizeof(row_t) * noRows);
	if (rows == NULL)  failed = 1;
	for (i = 0L; i < noRows; i++) {
		rows[i].originalRowNumber = i;
		rows[i].columnSet = (bitvector_t*)calloc(sizeof(bitvector_t), noBVs);
		if (rows[i].columnSet == NULL)
			failed = 1;
	}
	maxLevels = (noRows + 2L);
	consideredColumns = (cs_t*)calloc(sizeof(cs_t), maxLevels);
	if (consideredColumns == NULL)  failed = 1;
	else {
		for (i = 0L; i < maxLevels; i++) {
			consideredColumns[i] = (bitvector_t*)calloc(sizeof(bitvector_t), noBVs);
			if (consideredColumns[i] == NULL)  failed = 1;
		}
		if (!failed)
			for (i = 0L; i < noColumns; i++)
				setColumn(consideredColumns[0], i);
	}
	mandatoryColumns = (cs_t*)calloc(sizeof(cs_t), maxLevels);
	if (mandatoryColumns == NULL)  failed = 1;
	else {
		for (i = 0L; i < maxLevels; i++) {
			mandatoryColumns[i] = (bitvector_t*)calloc(sizeof(bitvector_t), noBVs);
			if (mandatoryColumns[i] == NULL)  failed = 1;
		}
	}
	columnIntersection = (bitvector_t*)calloc(sizeof(bitvector_t), noBVs);
	if (columnIntersection == NULL)  failed = 1;

	return !failed;
} /* initializeMemory */

void  readInDataMatrix(FILE  *fp)
{
	long  i, j, cell;
	// 	cout << endl;
	for (i = 0L; i < noRows; i++) {
		std::vector<int> tmp_vStoreRow;
		for (j = 0L; j < noColumns; j++) {
			fscanf_s(fp, "%d", &cell);
			tmp_vStoreRow.push_back(cell);
#ifndef DEBUG
			if (cell != 0)
				cout << i << " - " << j << endl;
#endif
			// 	cout << cell << "\t";
			if (cell == 0)
				unsetColumn(rows[i].columnSet, j);
			else
				setColumn(rows[i].columnSet, j);
		}
		store.push_back(tmp_vStoreRow);
		// 	cout << endl;
	}
} /* readInDataMatrix */


double diffclock(clock_t clock1, clock_t clock2)
{
	double diffticks = clock1 - clock2;
	double diffms = (diffticks * 10) / CLOCKS_PER_SEC;
	return diffms;
}


long timediff(clock_t t1, clock_t t2) {
	long elapsed;
	elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000;
	return elapsed;
}

int quadsearch(
	int IL,
	int JL,
	int IR,
	int JR)
{
	int maxedges;
	int m1 = 0, m2 = 0, m3 = 0;
	int I, J;
	int exists;

	printf("quadsearch(%d,%d,%d,%d)\n", IL, JL, IR, JR);
	if (IL > IR) return(0);
	if (JL > JR) return(0);
	I = (IL + IR) / 2;
	J = (JL + JR) / 2;

	exists = false;
	if (maxobjval < I*J) {
		minNoRows = I;
		minNoColumns = J;

#ifndef DEBUG_ENABLE
		cout << "Conquer " << I << " x " << J << endl;
#endif
		exists = conquer(0L, noRows - 1L, 0L, 0L);
	}
	maxedges = I*J;
	if (exists) {
		if (maxedged < maxedges && maxnotfound > maxedges ) m1 = quadsearch(I + 1, J + 1, IR, JR);
		if (maxedged < maxedges && maxnotfound > maxedges) m2 = quadsearch(IL, J + 1, I, JR);
		if (maxedged < maxedges && maxnotfound > maxedges) m3 = quadsearch(I + 1, JL, IR, J);
	}
	else {
		if (maxedged < maxedges && maxnotfound > maxedges) m1 = quadsearch(IL, JL, I - 1, J - 1);
		if (maxedged < maxedges && maxnotfound > maxedges) m2 = quadsearch(IL, J, I - 1, JR);
		if (maxedged < maxedges && maxnotfound > maxedges) m3 = quadsearch(I, JL, IR, J - 1);

		if (maxnotfound > maxedges)
			maxnotfound = maxedges;
		maxedges = 0;
	}

	if (maxedges < m1) {
		maxedges = m1;
	}
	if (maxedges < m2) {
		maxedges = m2;
	}
	if (maxedges < m3) {
		maxedges = m3;
	}
	return(maxedges);
}

int shallwesolve(
	int I,
	int J)
{
	list<Coord>::iterator li;

	if (newln*newrn < I*J) return(0);
	for (li = clist.begin(); li != clist.end(); li++) {
		if (((*li).x <= I) && ((*li).y <= J)) return(0);
	}
	return(1);
}

int main(int argc, char *argv[])
{
	double solnTime = 0;
	int numberOfBiclusters = 0;
	int maxDim1 = 0, maxDim2 = 0;
	long int timeBefore = 0;
	long int timeAfter = 0;

	biclusterCount = 10000000;
	FILE  *fp;
	FILE *fptr, *fptr2;
	fopen_s(&fptr, "BIMAXResult.txt", "w");

	cout << "/**************************************************/" << endl;
	cout << "	       BIMAX'S RUN " << endl;
	cout << "/**************************************************/" << endl;
	fclose(fptr);
	if (argv[1] != NULL)
		fopen_s(&fp, argv[1], "r");
	else
		fopen_s(&fp, "example.txt", "r");

	if (fp == NULL)  exit(0);
	fscanf_s(fp, "%ld", &noRows);
	fscanf_s(fp, "%ld", &noColumns);
	fscanf_s(fp, "%ld", &minNoRows);
	fscanf_s(fp, "%ld", &minNoColumns);
	maxnotfound = noRows * noColumns;
	//   cout << noRows << " - " << noColumns << " - " << minNoRows << " - " << minNoColumns << endl;
	if (minNoRows < 1L)
		minNoRows = 1L;
	if (minNoColumns < 1L)
		minNoColumns = 1L;
	if (noColumns > 0L && noRows > 0L && initialize()) {

		readInDataMatrix(fp);
		int rc;
		rc = quadsearch(1, 1, noRows, noColumns);

		if (fopen_s(&fptr2, "BIMAXResult.txt", "r"))
		{
			char strToRead[16] = "";
			long maxCase = 0;
			int dim1 = 0, dim2 = 0;
			cout << "Check File until the end" << endl;

			while (!feof(fptr2))
			{
				fscanf_s(fptr2, "%d%d", &dim1, &dim2);
				//cout << dim1 << " - " << dim2 << endl;
				for (int count = 0; count < dim1; count++)
				{
					fscanf_s(fptr2, "%s", strToRead);
				}

				for (int count = 0; count < dim2; count++)
				{
					fscanf_s(fptr2, "%s", strToRead);
				}

				if (dim1 * dim2 > maxCase)
				{
					maxDim1 = dim1;
					maxDim2 = dim2;
					maxCase = dim1 * dim2;
				}
				numberOfBiclusters++;

			}
			numberOfBiclusters--;
			fclose(fptr2);
		}
		else
		{
			cout << "File Could not Be opened" << endl;
		}


		fopen_s(&fptr2, "statusBimax.txt", "w");

		fprintf_s(fptr2, "%d\t%d\t%d\t%d\t%lf\t%lf", numberOfBiclusters, maxDim1, maxDim2);
		fclose(fptr2);

	}

	return 0;
} /* main */