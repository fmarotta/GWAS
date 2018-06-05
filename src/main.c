/*Exploits an ADT interface to perform a genome-wide association study*/

// TODO support output files?

#include "../include/gwas.h"
#include <argp.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>

const char * argp_program_version = "GWAS 1.0";
const char * argp_program_bug_address = "<federicomarotta@mail.com>";

/* Command-line options parsing {{{
 ********************************/

/* Program documentation. */
static char doc[] = "GWAS -- an educational attempt to perform a GWAS"
"\v"
"This program uses an ADT interface to perform:"
"	1) Reading of .ped and .map files into a suitable data structure;"
"	2) Preliminar analysis of the data;"
"	3) Quality control of samples and markers;"
"	4) Testing for association and displaying of the results.";

/* A description of the arguments we accept. */
static char args_doc[] = "PED_FILE MAP_FILE";

/* The options we understand. */
static struct argp_option options[] = {
	{0, 0, 0, 0, "General options" }, /* options header */
	{"output", 'o', "FILE", 0, "Output to FILE instead of standard output" },
	{ 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
	char *args[2]; /* ped_file & map_file */
	char *output_file;
};

/* Parse each option. */
static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	/* Get the input argument from argp_parse, which we know is a
	 * pointer to our arguments structure. */
	struct arguments *arguments = state->input;

	switch (key)
    {
		case 'o':
			arguments->output_file = arg;
			break;

		case ARGP_KEY_ARG:
			if (state->arg_num >= 2)
			/* Too many arguments. */
			argp_usage (state);
			arguments->args[state->arg_num] = arg;
			break;

		case ARGP_KEY_END:
			if (state->arg_num < 2)
			/* Not enough arguments. */
			argp_usage (state);
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

// }}}

int main(int argc, char *argv[])
{
	struct arguments arguments;
	FILE *ped_file, *map_file;
	GWAS_COHORT cohort;

	/* Default arguments values. */
	arguments.args[0] = "";
	arguments.args[1] = "";
	arguments.output_file = "-";

	/* Parse our arguments; every option seen by parse_opt will be
	 * reflected in arguments. */
	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	/* Open ped and map files */
	if ((ped_file = fopen(arguments.args[0], "r")) == NULL)
	{
		fprintf(stderr, "Error: could not open ped file %s.\n",
				arguments.args[0]);
		exit(EXIT_FAILURE);
	}
	if ((map_file = fopen(arguments.args[1], "r")) == NULL)
	{
		fprintf(stderr, "Error: could not open map file %s.\n",
				arguments.args[0]);
		exit(EXIT_FAILURE);
	}

	Initialize_cohort(ped_file, map_file, &cohort);

	Test_association(&cohort, "allelic");

	return 0;
}

