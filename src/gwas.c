/* Interface implementation for an ADT structure to perform a
 * genome-wide association study (GWAS)
 */

// TODO support counting of genotypes.

#include "../include/gwas.h"
#include "../utils/io_utils.h"
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef LOGFILE
#define LOGFILE stderr
#endif

static GWAS_MARKER *create_marker(unsigned int, char *, float, unsigned long);
static GWAS_SAMPLE *create_sample(char *, char *, char *, char*, int, int);
static void append_marker(GWAS_MARKER *, GWAS_COHORT *);
static void append_sample(GWAS_SAMPLE *, GWAS_COHORT *);
static void get_alleles(FILE *, char *);
static int count_alleles(char *, int, GWAS_MARKER *);
static void count_genotypes(char *, int, GWAS_MARKER *);

// Initialize_cohort {{{
void Initialize_cohort(FILE *ped_file, FILE *map_file, GWAS_COHORT *pcohort)
{
	GWAS_MARKER *pmarker;
	GWAS_SAMPLE *psample;

	unsigned int m_chr;
	char m_id[30];
	float m_dist;
	unsigned long m_pos;

	char family_id[30];
	char individual_id[30];
	char father_id[30];
	char mother_id[30];
	int sex;
	int condition;
	char alleles[2];
	unsigned int errors;

	/* Initialize cohort structure members */
	pcohort->n_markers = 0;
	pcohort->n_samples = 0;
	pcohort->n_cases = 0;
	pcohort->n_controls = 0;
	pcohort->n_missing_phenotype = 0;
	pcohort->markers = NULL;
	pcohort->samples = NULL;

	/* While reading the map file, allocate space for the loci. */
	while (fscanf(map_file, "%d%s%f%d", &m_chr, m_id, &m_dist, &m_pos) != EOF)
	{
		if ((pmarker = create_marker(m_chr, m_id, m_dist, m_pos)) == NULL)
		{
			fprintf(stderr, "Error: run out of memory at marker %s\n",
					m_id);
			exit(EXIT_FAILURE);
		}
		append_marker(pmarker, pcohort);
	}

	/* Read the constant fields of the ped file, then the alleles. */
	while (fscanf(ped_file, "%s%s%s%s%d%d", family_id, individual_id,
				father_id, mother_id, &sex, &condition) != EOF)
	{
		/* Read a sample */
		if ((psample = create_sample(family_id, individual_id,
						father_id, mother_id, sex, condition)) == NULL)
		{
			fprintf(stderr, "Error: run out of memory at sample %s\n",
					individual_id);
			exit(EXIT_FAILURE);
		}
		append_sample(psample, pcohort);

		/* Read a pair of alleles for each marker */
		errors = 0;
		pmarker = pcohort->markers;
		while (pmarker != NULL)
		{
			get_alleles(ped_file, alleles);
			errors += count_alleles(alleles, condition, pmarker);
			//count_genotypes(alleles, condition, pmarker);
			pmarker = pmarker->next;
		}
		psample->n_invalid_markers = errors;
	}

	/* Compute allele frequencies. */
	pmarker = pcohort->markers;
	while (pmarker != NULL)
	{
		pmarker->cases.alleles[0].freq =
			(float) pmarker->cases.alleles[0].count / pcohort->n_cases;
		pmarker->controls.alleles[1].freq =
			(float) pmarker->controls.alleles[1].count / pcohort->n_controls;
		pmarker = pmarker->next;
	}

	/* Logging: print cohort stats. */
	fprintf(LOGFILE, "Cohort initialized with %d individuals, of which %d cases "
			"and %d controls (%d missing phenotypes), and %d markers.\n",
			pcohort->n_samples, pcohort->n_cases, pcohort->n_controls,
			pcohort->n_missing_phenotype, pcohort->n_markers);
}

static GWAS_MARKER *create_marker(unsigned int m_chr, char *m_id,
		float m_dist, unsigned long m_pos)
{
	/* NOTE: m_dist (position of the marker in centimorgans) is ignored
	 * here. */
	GWAS_MARKER *pmarker;

	pmarker = (GWAS_MARKER *) malloc(sizeof(GWAS_MARKER));
	if (pmarker == NULL)
		return NULL;

	strncpy(pmarker->id, m_id, 30);
	pmarker->chr = m_chr;
	pmarker->pos = m_pos;
	pmarker->cases.alleles[0].seq = '\0';
	pmarker->cases.alleles[0].count = 0;
	pmarker->cases.alleles[0].freq = 0;
	pmarker->cases.alleles[1].seq = '\0';
	pmarker->cases.alleles[1].count = 0;
	pmarker->cases.alleles[1].freq = 0;
	/*pmarker->OR = 0; // I leave them undefined
	pmarker->E = 0;
	pmarker->Pvalue = 0;*/
	pmarker->next = NULL;

	return pmarker;
}

static GWAS_SAMPLE *create_sample(char *family_id, char *individual_id,
		char *father_id, char *mother_id, int sex, int condition)
{
	GWAS_SAMPLE *psample;

	psample = (GWAS_SAMPLE *) malloc(sizeof(GWAS_SAMPLE));
	if (psample == NULL)
		return NULL;

	strncpy(psample->family_id, family_id, FIELDLEN);
	strncpy(psample->individual_id, individual_id, FIELDLEN);
	strncpy(psample->father_id, father_id, FIELDLEN);
	strncpy(psample->mother_id, mother_id, FIELDLEN);
	psample->sex = sex;
	psample->condition = condition;
	psample->next = NULL;

	return psample;
}

static void append_marker(GWAS_MARKER *pmarker, GWAS_COHORT *pcohort)
{
	GWAS_MARKER *pnode;

	if (pcohort->markers == NULL)
		pcohort->markers = pmarker;
	else
	{
		pnode = pcohort->markers;
		while (pnode->next != NULL)
			pnode = pnode->next;
		pnode->next = pmarker;
	}

	pcohort->n_markers++;
}

static void append_sample(GWAS_SAMPLE *psample, GWAS_COHORT *pcohort)
{
	GWAS_SAMPLE *pnode;

	if (pcohort->samples == NULL)
		pcohort->samples = psample;
	else
	{
		pnode = pcohort->samples;
		while (pnode->next != NULL)
			pnode = pnode->next;
		pnode->next = psample;
	}

	/* Update cohort stats */
	if (psample->condition == 1)
		pcohort->n_controls++;
	else if (psample->condition == 2)
		pcohort->n_cases++;
	else
		pcohort->n_missing_phenotype++;
	pcohort->n_samples++;
}

static void get_alleles(FILE *ped_file, char *alleles)
{
	char ch;

	/* Read the next two non-blank characters */
	while (isspace(ch = fgetc(ped_file)))
		continue;
	alleles[0] = ch;
	while (isspace(ch = fgetc(ped_file)))
		continue;
	alleles[1] = ch;
}

/* Returns 0 on success, 1 if alleles are invalid */
static int count_alleles(char *alleles, int condition, GWAS_MARKER *pmarker)
{
	GWAS_FREQ *group;

	if (condition == 2)
		group = &(pmarker->cases);
	else
		group = &(pmarker->controls);

	/* Handle missing alleles. */
	if (strchr("ACGTacgt", alleles[0]) == NULL ||
			strchr("ACTGacgt", alleles[1]) == NULL)
	{
		group->n_invalid_data++;
		return 1;
	}

	/* Count the valid alleles. */
	if (group->alleles[0].seq == '\0')
	{
		if (alleles[0] == alleles[1])
		{
			group->alleles[0].seq = alleles[0];
			group->alleles[0].count += 2;
		}
		else
		{
			group->alleles[0].seq = alleles[0];
			group->alleles[0].count++;
			group->alleles[1].seq = alleles[1];
			group->alleles[1].count++;
		}
	}
	else
	{
		if (alleles[0] == group->alleles[0].seq)
			group->alleles[0].count++;
		else if (alleles[0] == group->alleles[1].seq)
			group->alleles[1].count++;
		else if (group->alleles[1].seq == '\0')
		{
			group->alleles[1].seq = alleles[0];
			group->alleles[1].count++;
		}
		else
			fprintf(stderr, "Warning: apparently marker %s is multiallelic, "
					"but we only support biallelic loci.\n", pmarker->id);

		if (alleles[1] == group->alleles[0].seq)
			group->alleles[0].count++;
		else if (alleles[1] == group->alleles[1].seq)
			group->alleles[1].count++;
		else if (group->alleles[1].seq == '\0')
		{
			group->alleles[1].seq = alleles[1];
			group->alleles[1].count++;
		}
		else
			fprintf(stderr, "Warning: apparently marker %s is multiallelic, "
					"but we only support biallelic loci.\n", pmarker->id);
	}

	return 0;
}

static void count_genotypes(char *alleles, int condition, GWAS_MARKER *pmarker)
{
	return;
}
// }}}
