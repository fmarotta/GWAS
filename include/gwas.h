/* Interface header for an ADT structure to perform a genome-wide association
 * study (GWAS)
 */

// TODO pass arguments to assay and qc-filter to specify what to assay or what
// to filter on. Good opportunity to use variable-length stuff.

// TODO give the chance to choose the penetrance model.

// TODO implement covariates, at least sex and perhaps population.

#ifndef GWAS_H_
#define GWAS_H_

#include <stdbool.h>
#include <stdio.h>

#define FIELDLEN 30

/* data structures
 *******************/

typedef struct Gwas_freq {
	struct Alleles {
		char seq;
		unsigned int count;
	} alleles[2];
	struct Genotypes {
		char seq[3]; // i.e. 0/1, 1/1 or 0|0 (if phased, like in VCFs)
		unsigned int count;
	} genotypes[3];
	unsigned int n_invalid_data;
} GWAS_FREQ;

/* node of a linked list */
typedef struct Gwas_marker {
	char id[30];
	unsigned int chr;
	unsigned long pos;
	GWAS_FREQ cases;
	GWAS_FREQ controls;
	float OR;	// Odds Ratio
	float E;	// Expected OR
	float RR;	// Relative Risk
	float Pvalue;
	struct Gwas_marker *next;
} GWAS_MARKER;

/* node of another linked list */
typedef struct Gwas_sample {
	char family_id[30];
	char individual_id[30];
	char father_id[30];
	char mother_id[30];
	int sex;
	int condition;
	unsigned int n_invalid_markers;
	struct Gwas_sample *next;
} GWAS_SAMPLE;

typedef struct Gwas_cohort {
	GWAS_MARKER *markers;
	GWAS_SAMPLE *samples;
	unsigned int n_markers;
	unsigned int n_samples;
	unsigned int n_cases;
	unsigned int n_controls;
	unsigned int n_missing_phenotype;
} GWAS_COHORT;

/* operations
 **************/

/* operation:		Reads the .ped and .map files into a COHORT structure.
 * precondition:	The relevant files are fopen'd; a pointer to COHORT is
 *					defined.
 * postcondition:	Stores all the data in the structure. */
void Initialize_cohort(FILE *ped_file, FILE *map_file, GWAS_COHORT *pcohort);

/* operation:		Perform preliminar analysis of data in the cohort.
 * precondition:	pcohort points to an initialized cohort.
 * postcondition:	Basic statistical information is produced about the data,
 * 					which can give hints about the further steps of analysis. */
void Assay_cohort(GWAS_COHORT *pcohort);

/* operation:		Apply some quality control filters to the cohort.
 * precondition:	pcohort points to an initialized cohort; users know what
 *					they are doing because they have run Assay_cohort().
 * postcondition:	Removes from the cohort markers or individuals that fail
 *					some quality control. */
void QC_filter_cohort(GWAS_COHORT *pcohort);

/* operation:		Perform a chi-square test of independence in order to find
 *					SNPs associated to the phenotype, using allelic model.
 * precondition:	pcohort points to an initialized cohort.
 * postcondition:	A chi-square test is performed and the P-value (already
 *					adjusted for multiple testing) is obtained. */
void Test_association(GWAS_COHORT *pcohort);

/* operation:		Prints a table with the results of the associations.
 * precondition:	pcohort points to an initialized cohort.
 * postcondition:	A txt table is provided for interpretation of results. */
void Print_results(GWAS_COHORT *pcohort);

#endif
