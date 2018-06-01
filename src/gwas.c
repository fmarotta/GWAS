/* Interface implementation for an ADT structure to perform a
 * genome-wide association study (GWAS)
 */

#include "../include/gwas.h"
#include "../utils/io_utils.h"
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static GWAS_MARKER *create_marker(unsigned int, char *, float, unsigned long);
static void append_marker(GWAS_MARKER *, GWAS_COHORT *);

void Initialize_cohort(FILE *ped_file, FILE *map_file, GWAS_COHORT *pcohort)
{
	GWAS_MARKER *pmarker;

	unsigned int m_chr;
	char m_id[30];
	float m_dist;
	unsigned long m_pos;

	int condition;
	char field[10];
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

	/* Debug: print all the markers
	pmarker = pcohort->markers;
	while (pmarker != NULL)
	{
		printf("marker %s\n", pmarker->id);
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
	/*pmarker->OR = 0;
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
