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

static GWAS_MARKER *create_marker(unsigned int, char *, float, unsigned long);
static void append_marker(GWAS_MARKER *, GWAS_COHORT *);
static void get_alleles(FILE *, char *);
static void count_alleles(char *, int, GWAS_MARKER *);
static void count_genotypes(char *, int, GWAS_MARKER *);

void Initialize_cohort(FILE *ped_file, FILE *map_file, GWAS_COHORT *pcohort)
{
	GWAS_MARKER *pmarker;

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

	/* Initialize cohort */
	pcohort->n_markers = 0;
	pcohort->n_samples = 0;
	pcohort->n_cases = 0;
	pcohort->n_controls = 0;
	pcohort->markers = NULL;

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
		//printf("condition at beginning: %d ", condition);
		if (condition == 1)
			pcohort->n_controls++;
		else if (condition == 2)
			pcohort->n_cases++;
		else
		{
			EATLINE(ped_file);
			continue;
		}

		/* Read a pair of alleles for each marker */
		int i = 0;
		pmarker = pcohort->markers;
		while (pmarker != NULL)
		{
			//fscanf(ped_file, "%s%s", &alleles[0], &alleles[1]);
			get_alleles(ped_file, alleles);

			if (strchr("ACGTacgt", alleles[0]) == NULL)
			{
				// FIXME keep track of missing data.
				pmarker = pmarker->next;
				i++;
				continue;
			}

			count_alleles(alleles, condition, pmarker);
			//count_genotypes(alleles, condition, pmarker);

			pmarker = pmarker->next;
			i++;
		}
	}

	/* Debug: print all the markers */
	/*
	pmarker = pcohort->markers;
	while (pmarker != NULL)
	{
		printf("marker %s: %c %c\n", pmarker->id,
				pmarker->cases.alleles[0].seq, pmarker->cases.alleles[1].seq);
		pmarker = pmarker->next;
	}
	*/
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


static void count_alleles(char *alleles, int condition, GWAS_MARKER *pmarker)
{
	GWAS_FREQ *group;

	if (condition == 2)
		group = &(pmarker->cases);
	else
		group = &(pmarker->controls);

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
		if (alleles[0] == alleles[1])
		{
			if (alleles[0] == group->alleles[0].seq)
			{
				group->alleles[0].seq = alleles[0];
				group->alleles[0].count += 2;
			}
			else
			{
				group->alleles[1].seq = alleles[0];
				group->alleles[1].count += 2;
			}
		}
		else
		{
			if (alleles[0] = group->alleles[0].seq)
			{
				group->alleles[0].seq = alleles[0];
				group->alleles[0].count++;
				group->alleles[1].seq = alleles[1];
				group->alleles[1].count++;
			}
			else if (alleles[1] = group->alleles[0].seq)
			{
				group->alleles[0].seq = alleles[1];
				group->alleles[0].count++;
				group->alleles[1].seq = alleles[0];
				group->alleles[1].count++;
			}
		}
	}
}

static void count_genotypes(char *alleles, int condition, GWAS_MARKER *pmarker)
{
	return;
}
