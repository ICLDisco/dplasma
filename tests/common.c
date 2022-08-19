/*
 * Copyright (c) 2009-2021 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 */
#include "parsec/runtime.h"
#include "parsec/execution_stream.h"
#include "parsec/utils/mca_param.h"
#include "parsec/utils/show_help.h"
#include "dplasma.h"

#include "common.h"

#include <stdlib.h>
#include <stdio.h>
#ifdef PARSEC_HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef PARSEC_HAVE_LIMITS_H
#include <limits.h>
#endif
#if defined(PARSEC_HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(PARSEC_HAVE_GETOPT_H) */
#ifdef PARSEC_HAVE_MPI
#include <mpi.h>
#endif
#if defined(DPLASMA_HAVE_CUDA)
#include "dplasmaaux.h"
#include <cublas.h>
#include <cusolverDn.h>
#endif

char *PARSEC_SCHED_NAME[] = {
    "", /* default */
    "lfq",
    "ltq",
    "ap",
    "lhq",
    "gd",
    "pbq",
    "ip",
    "rnd",
};

/*******************************
 * globals and argv set values *
 *******************************/
#if defined(PARSEC_HAVE_MPI)
MPI_Datatype SYNCHRO = MPI_BYTE;
#endif  /* PARSEC_HAVE_MPI */

const int   side[2]  = { dplasmaLeft,    dplasmaRight };
const int   uplo[2]  = { dplasmaUpper,   dplasmaLower };
const int   diag[2]  = { dplasmaNonUnit, dplasmaUnit  };
const int   trans[3] = { dplasmaNoTrans, dplasmaTrans, dplasmaConjTrans };
const int   norms[4] = { dplasmaMaxNorm, dplasmaOneNorm, dplasmaInfNorm, dplasmaFrobeniusNorm };

const char *sidestr[2]  = { "Left ", "Right" };
const char *uplostr[2]  = { "Upper", "Lower" };
const char *diagstr[2]  = { "NonUnit", "Unit   " };
const char *transstr[3] = { "N", "T", "H" };
const char *normsstr[4] = { "Max", "One", "Inf", "Fro" };

double time_elapsed = 0.0;
double sync_time_elapsed = 0.0;
double alpha = 1.;

/**********************************
 * Command line arguments
 **********************************/
void print_usage(void)
{
    fprintf(stderr,
            "Mandatory argument:\n"
            " number            : dimension (N) of the matrices (required)\n"
            "Optional arguments:\n"
            " -p -P --grid-rows : rows (P) in the PxQ process grid   (default: NP)\n"
            " -q -Q --grid-cols : columns (Q) in the PxQ process grid (default: NP/P)\n"
            "\n"
            " -N                : dimension (N) of the matrices (required)\n"
            " -M                : dimension (M) of the matrices (default: N)\n"
            " -K --NRHS C<-A*B+C: dimension (K) of the matrices (default: N)\n"
            "           AX=B    : columns in the right hand side (default: 1)\n"
            " -A --LDA          : leading dimension of the matrix A (default: full)\n"
            " -B --LDB          : leading dimension of the matrix B (default: full)\n"
            " -C --LDC          : leading dimension of the matrix C (default: full)\n"
            " -i --IB           : inner blocking     (default: autotuned)\n"
            " -t --MB           : rows in a tile     (default: autotuned)\n"
            " -T --NB           : columns in a tile  (default: autotuned)\n"
            " -s --SMB --kp     : repetitions in the process grid row cyclicity (formerly supertiles) (default: 1)\n"
            " -S --SNB --kq     : repetitions in the process grid column cyclicity (formerly supertiles) (default: 1)\n"
            " -z --HNB --HMB    : Inner NB/MB used for recursive algorithms (default: MB)\n"
            " -x --check        : verify the results\n"
            " -X --check_inv    : verify the results against the inverse\n"
            " -b --sync         : call the step by step version of the algorithm if exists\n"
            "\n"
            "    --qr_a         : Size of TS domain. (specific to xgeqrf_param)\n"
            "    --qr_p         : Size of the high level tree. (specific to xgeqrf_param)\n"
            " -d --domino       : Enable/Disable the domino between upper and lower trees. (specific to xgeqrf_param) (default: 1)\n"
            " -r --tsrr         : Enable/Disable the round-robin on TS domain. (specific to xgeqrf_param) (default: Disabled)\n"
            "    --treel        : Tree used for low level reduction inside nodes. (specific to xgeqrf_param)\n"
            "    --treeh        : Tree used for high level reduction between nodes, only if qr_p > 1. (specific to xgeqrf_param)\n"
            "                      (0: Flat, 1: Greedy, 2: Fibonacci, 3: Binary)\n"
            "\n"
            "    --criteria     : Choice of the criteria to switch between LU and QR\n"
            "                      (0: Alternate, 1: Higham, 2: MUMPS (specific to xgetrf_qrf)\n"
            " -a --alpha        : Threshold to swith back to QR. (specific to xgetrf_qrf)\n"
            "    --seed         : Set the seed for pseudo-random generator\n"
            "    --mtx          : Set the matrix generator (Default: 0, random)\n"
            "\n"
            " -y --butlvl       : Level of the Butterfly (starting from 0).\n"
            "\n"
            " --nruns            : Number of times to run the kernel.\n"
            "\n"
            " -v --verbose      : extra verbose output\n"
            " -h --help         : this message\n"
            "\n"
            );
    fprintf(stderr,
            "\n"
            " -c --cores        : number of concurent threads (default: number of physical hyper-threads)\n"
            " -g --gpus         : number of GPU (default: 0)\n"
            " -m --thread_multi : initialize MPI_THREAD_MULTIPLE (default: no)\n"
            " -o --scheduler    : select the scheduler (default: LFQ)\n"
            "                     Accepted values:\n"
            "                       LFQ -- Local Flat Queues\n"
            "                       LTQ -- Local Tree Queues\n"
            "                       AP  -- Absolute Priorities\n"
            "                       LHQ -- Local Hierarchical Queues\n"
            "                       GD  -- Global Dequeue\n"
            "                       PBQ -- Priority Based Local Flat Queues\n"
            "                       IP  -- Inverse Priorities\n"
            "                       RND -- Random\n"
            "\n"
            "    --dot          : create a dot output file (default: don't)\n"
            "\n"
            "    --ht nbth      : enable a SMT/HyperThreadind binding using nbth hyper-thread per core.\n"
            "                     This parameter must be declared before the virtual process distribution parameter\n"
            " -V --vpmap        : select the virtual process map (default: flat map)\n"
            "                     Accepted values:\n"
            "                       flat  -- Flat Map: all cores defined with -c are under the same virtual process\n"
            "                       hwloc -- Hardware Locality based: threads up to -c are created and threads\n"
            "                                bound on cores that are under the same socket are also under the same\n"
            "                                virtual process\n"
            "                       rr:n:p:c -- create n virtual processes per real process, each virtual process with p threads\n"
            "                                   bound in a round-robin fashion on the number of cores c (overloads the -c flag)\n"
            "                       file:filename -- uses filename to load the virtual process map. Each entry details a virtual\n"
            "                                        process mapping using the semantic  [mpi_rank]:nb_thread:binding  with:\n"
            "                                        - mpi_rank : the mpi process rank (empty if not relevant)\n"
            "                                        - nb_thread : the number of threads under the virtual process\n"
            "                                                      (overloads the -c flag)\n"
            "                                        - binding : a set of cores for the thread binding. Accepted values are:\n"
            "                                          -- a core list          (exp: 1,3,5-6)\n"
            "                                          -- a hexadecimal mask   (exp: 0xff012)\n"
            "                                          -- a binding range expression: [start];[end];[step] \n"
            "                                             wich defines a round-robin one thread per core distribution from start\n"
            "                                             (default 0) to end (default physical core number) by step (default 1)\n"
            "\n"
            "\n"
            "ENVIRONMENT\n"
            "  [SDCZ]<FUNCTION> : defines the priority limit of a given function for a given precision\n"
            "\n");
            parsec_usage();
}

#define GETOPT_STRING "bc:mo:g::p:P:q:Q:N:M:K:A:B:C:i:t:T:s:S:xXv::hd:ry:V:a:R:G:"

#if defined(PARSEC_HAVE_GETOPT_LONG)
static struct option long_options[] =
{
    /* PaRSEC specific options */
    {"cores",       required_argument,  0, 'c'},
    {"c",           required_argument,  0, 'c'},
    {"thread_multi",no_argument,        0, 'm'},
    {"m",           no_argument,        0, 'm'},
    {"o",           required_argument,  0, 'o'},
    {"scheduler",   required_argument,  0, 'o'},
    {"gpus",        required_argument,  0, 'g'},
    {"g",           required_argument,  0, 'g'},
    {"V",           required_argument,  0, 'V'},
    {"vpmap",       required_argument,  0, 'V'},
    {"ht",          required_argument,  0, 'H'},
    {"dot",         required_argument,  0, '.'},

    /* Generic Options */
    {"grid-rows",   required_argument,  0, 'p'},
    {"p",           required_argument,  0, 'p'},
    {"P",           required_argument,  0, 'p'},
    {"grid-cols",   required_argument,  0, 'q'},
    {"q",           required_argument,  0, 'q'},
    {"Q",           required_argument,  0, 'q'},

    {"N",           required_argument,  0, 'N'},
    {"M",           required_argument,  0, 'M'},
    {"K",           required_argument,  0, 'K'},
    {"NRHS",        required_argument,  0, 'K'},
    {"LDA",         required_argument,  0, 'A'},
    {"A",           required_argument,  0, 'A'},
    {"LDB",         required_argument,  0, 'B'},
    {"B",           required_argument,  0, 'B'},
    {"LDC",         required_argument,  0, 'C'},
    {"C",           required_argument,  0, 'C'},
    {"IB",          required_argument,  0, 'i'},
    {"i",           required_argument,  0, 'i'},
    {"MB",          required_argument,  0, 't'},
    {"t",           required_argument,  0, 't'},
    {"NB",          required_argument,  0, 'T'},
    {"T",           required_argument,  0, 'T'},
    {"kp",          required_argument,  0, 's'},
    {"SMB",         required_argument,  0, 's'},
    {"s",           required_argument,  0, 's'},
    {"kq",          required_argument,  0, 'S'},
    {"SNB",         required_argument,  0, 'S'},
    {"S",           required_argument,  0, 'S'},
    {"check",       no_argument,        0, 'x'},
    {"x",           no_argument,        0, 'x'},
    {"check_inv",   no_argument,        0, 'X'},
    {"X",           no_argument,        0, 'X'},

    {"sync",        no_argument,        0, 'b'},
    {"b",           no_argument,        0, 'b'},

    /* HQR options */
    {"qr_a",        required_argument,  0, '0'},
    {"qr_p",        required_argument,  0, '1'},
    {"d",           required_argument,  0, 'd'},
    {"domino",      required_argument,  0, 'd'},
    {"r",           no_argument,        0, 'r'},
    {"tsrr",        no_argument,        0, 'r'},
    {"treel",       required_argument,  0, 'l'},
    {"treeh",       required_argument,  0, 'L'},

    /* LU-QR options */
    {"criteria",    required_argument,  0, '1'},
    {"alpha",       required_argument,  0, 'a'},
    {"seed",        required_argument,  0, 'R'},
    {"mtx",         required_argument,  0, 'G'},

    /* Recursive options */
    {"z",           required_argument,  0, 'z'},
    {"HNB",         required_argument,  0, 'z'},
    {"HMB",         required_argument,  0, 'z'},

    /* HERBT options */
    {"butlvl",      required_argument,  0, 'y'},
    {"y",           required_argument,  0, 'y'},

    /* Auxiliary options */
    {"verbose",     optional_argument,  0, 'v'},
    {"v",           optional_argument,  0, 'v'},
    {"help",        no_argument,        0, 'h'},
    {"h",           no_argument,        0, 'h'},

    {"nruns",       required_argument,  0, '3'},

    {0, 0, 0, 0}
};
#endif  /* defined(PARSEC_HAVE_GETOPT_LONG) */

extern char **environ;

static void read_arguments(int *_argc, char*** _argv, int* iparam)
{
    int opt = 0;
    int rc, c;
    int argc = *_argc;
    char **argv = *_argv;
    char *add_dot = NULL;
    char *value;

    /* Default seed */
    iparam[IPARAM_RANDOM_SEED] = 3872;
    iparam[IPARAM_MATRIX_INIT] = dplasmaMatrixRandom;
    iparam[IPARAM_NRUNS] = 3;

    do {
#if defined(PARSEC_HAVE_GETOPT_LONG)
        c = getopt_long_only(argc, argv, "",
                        long_options, &opt);
#else
        c = getopt(argc, argv, GETOPT_STRING);
        (void) opt;
#endif  /* defined(PARSEC_HAVE_GETOPT_LONG) */

        // printf("%c: %s = %s\n", c, long_options[opt].name, optarg);
        switch(c)
        {
            case 'c': iparam[IPARAM_NCORES] = atoi(optarg); break;
            case '3': iparam[IPARAM_NRUNS] = atoi(optarg); break;
            case 'm': iparam[IPARAM_THREAD_MT] = 1; break;
            case 'o':
                if( !strcmp(optarg, "LFQ") )
                    iparam[IPARAM_SCHEDULER] = PARSEC_SCHEDULER_LFQ;
                else if( !strcmp(optarg, "LTQ") )
                    iparam[IPARAM_SCHEDULER] = PARSEC_SCHEDULER_LTQ;
                else if( !strcmp(optarg, "AP") )
                    iparam[IPARAM_SCHEDULER] = PARSEC_SCHEDULER_AP;
                else if( !strcmp(optarg, "LHQ") )
                    iparam[IPARAM_SCHEDULER] = PARSEC_SCHEDULER_LHQ;
                else if( !strcmp(optarg, "GD") )
                    iparam[IPARAM_SCHEDULER] = PARSEC_SCHEDULER_GD;
                else if( !strcmp(optarg, "PBQ") )
                    iparam[IPARAM_SCHEDULER] = PARSEC_SCHEDULER_PBQ;
                else if( !strcmp(optarg, "IP") )
                    iparam[IPARAM_SCHEDULER] = PARSEC_SCHEDULER_IP;
                else if( !strcmp(optarg, "RND") )
                    iparam[IPARAM_SCHEDULER] = PARSEC_SCHEDULER_RND;
                else {
                    fprintf(stderr, "#!!!!! malformed scheduler value %s (accepted: LFQ LTQ PBQ AP GD RND LHQ IP). Reverting to default LFQ\n",
                            optarg);
                    iparam[IPARAM_SCHEDULER] = PARSEC_SCHEDULER_LFQ;
                }
                parsec_setenv_mca_param( "mca_sched", PARSEC_SCHED_NAME[iparam[IPARAM_SCHEDULER]], &environ );
                break;

            case 'g':
#if !defined(DPLASMA_HAVE_CUDA)
                iparam[IPARAM_NGPUS] = DPLASMA_ERR_NOT_SUPPORTED; /* force an error message */
#endif
                if(iparam[IPARAM_NGPUS] == DPLASMA_ERR_NOT_SUPPORTED) {
                    fprintf(stderr, "#!!!!! This test does not have GPU support. GPU disabled.\n");
                    break;
                }
                if(optarg)  iparam[IPARAM_NGPUS] = atoi(optarg);
                else        iparam[IPARAM_NGPUS] = INT_MAX;

                rc = asprintf(&value, "%d", iparam[IPARAM_NGPUS]);
                parsec_setenv_mca_param( "device_cuda_enabled", value, &environ );
                free(value);
                break;

            case 'p': case 'P': iparam[IPARAM_P] = atoi(optarg); break;
            case 'q': case 'Q': iparam[IPARAM_Q] = atoi(optarg); break;
            case 'N': iparam[IPARAM_N] = atoi(optarg); break;
            case 'M': iparam[IPARAM_M] = atoi(optarg); break;
            case 'K': iparam[IPARAM_K] = atoi(optarg); break;
            case 'A': iparam[IPARAM_LDA] = atoi(optarg); break;
            case 'B': iparam[IPARAM_LDB] = atoi(optarg); break;
            case 'C': iparam[IPARAM_LDC] = atoi(optarg); break;

            case 'i': iparam[IPARAM_IB] = atoi(optarg); break;
            case 't': iparam[IPARAM_MB] = atoi(optarg); break;
            case 'T': iparam[IPARAM_NB] = atoi(optarg); break;
            case 's': iparam[IPARAM_KP] = atoi(optarg); break;
            case 'S': iparam[IPARAM_KQ] = atoi(optarg); break;

            case 'X': iparam[IPARAM_CHECKINV] = 1;
                      /* Fall through */
            case 'x': iparam[IPARAM_CHECK] = 1; iparam[IPARAM_VERBOSE] = max(2, iparam[IPARAM_VERBOSE]); break;

            case 'b': iparam[IPARAM_ASYNC] = 0; break;

                /* HQR parameters */
            case '0': iparam[IPARAM_QR_TS_SZE]    = atoi(optarg); break;
            case '1': iparam[IPARAM_QR_HLVL_SZE]  = atoi(optarg); break;

            case 'R': iparam[IPARAM_RANDOM_SEED]  = atoi(optarg); break;
            case 'G': iparam[IPARAM_MATRIX_INIT]  = atoi(optarg); break;

            case 'd': iparam[IPARAM_QR_DOMINO]    = atoi(optarg) ? 1 : 0; break;
            case 'r': iparam[IPARAM_QR_TSRR]      = 1; break;

            case 'l': iparam[IPARAM_LOWLVL_TREE]  = atoi(optarg); break;
            case 'L': iparam[IPARAM_HIGHLVL_TREE] = atoi(optarg); break;

                /* GETRF/QRF parameters */
            case 'a': alpha = atof(optarg); break;

                /* Butterfly parameters */
            case 'y': iparam[IPARAM_BUT_LEVEL] = atoi(optarg); break;

                /* Recursive parameters */
            case 'z': iparam[IPARAM_HNB] = iparam[IPARAM_HMB] = atoi(optarg); break;

            case 'v':
                if(optarg)  iparam[IPARAM_VERBOSE] = atoi(optarg);
                else        iparam[IPARAM_VERBOSE] = 2;
                break;

            case 'H':
                if( 0 == iparam[IPARAM_RANK] ) fprintf(stderr, "#!!!!! option '%s' deprecated in testing programs, it should be passed to PaRSEC instead in the exact same format after --\n", long_options[opt].name);
                exit(-10);  /* No kidding! */
                break;  /* because some compilers are just too annoying */

            case 'V':
                if( 0 == iparam[IPARAM_RANK] ) fprintf(stderr, "#!!!!! option '%s' deprecated in testing programs, it should be passed to PaRSEC instead in the exact same format after --\n", long_options[opt].name);
                exit(-10);  /* No kidding! */
                break;  /* because some compilers are just too annoying */

            case '.':
                add_dot = optarg;
                break;

            case 'h': print_usage(); exit(0);
                break;

            case '?': /* getopt_long already printed an error message. */
                exit(1);
                break;  /* because some compilers are just too annoying */

            default:
                break; /* Assume anything else is parsec/mpi stuff */
        }
    } while(-1 != c);

    if( NULL != add_dot ) {
        int i, has_dashdash = 0, has_parsecdot = 0;
        for(i = 1; i < argc; i++) {
            if( !strcmp( argv[i], "--") ) {
                has_dashdash = 1;
            }
            if( has_dashdash && !strncmp( argv[i], "--parsec_dot", 11 ) ) {
                has_parsecdot = 1;
                break;
            }
        }
        if( !has_parsecdot ) {
            char **tmp;
            int  tmpc;
            if( !has_dashdash ) {
                tmpc = *(_argc) + 2;
                tmp = (char **)malloc((tmpc+1) * sizeof(char*));
                tmp[ tmpc - 2 ] = strdup("--");
            } else {
                tmpc = *(_argc) + 1;
                tmp = (char **)malloc((tmpc+1) * sizeof(char*));
            }
            for(i = 0; i < (*_argc);i++)
                tmp[i] = (*_argv)[i];

            rc = asprintf( &tmp[ tmpc - 1 ], "--parsec_dot=%s", add_dot );
            tmp[ tmpc     ] = NULL;

            *_argc = tmpc;
            *_argv = tmp;
        }
    }
    
    /* Set matrices dimensions to default values if not provided */
    /* Search for N as a bare number if not provided by -N */
    if(0 == iparam[IPARAM_N] && optind < argc) {
        iparam[IPARAM_N] = atoi(argv[optind++]);
    }

    /* If we ask for check (-x), then we do '0' run, i.e. only the 'warmup'
     * run, which is sometimes used to initialize the data and do the initial
     * computation */
    if(iparam[IPARAM_CHECK]) {
        iparam[IPARAM_NRUNS] = 0;
    }
    (void)rc;
}

static void parse_arguments(int *iparam) {
    int verbose = iparam[IPARAM_RANK] ? 0 : iparam[IPARAM_VERBOSE];

    if(iparam[IPARAM_NGPUS] < 0) iparam[IPARAM_NGPUS] = 0;
    if(iparam[IPARAM_NGPUS] > 0) {
        if (iparam[IPARAM_VERBOSE] > 3) {
            parsec_setenv_mca_param( "device_show_capabilities", "1", &environ );
        }
        if (iparam[IPARAM_VERBOSE] > 2) {
            parsec_setenv_mca_param( "device_show_statistics", "1", &environ );
        }
    }

    /* Check the process grid */
    if(0 == iparam[IPARAM_P])
        iparam[IPARAM_P] = iparam[IPARAM_NNODES];
    else if(iparam[IPARAM_P] > iparam[IPARAM_NNODES])
    {
        fprintf(stderr, "#XXXXX There are only %d nodes in the world, and you requested P=%d\n",
                iparam[IPARAM_NNODES], iparam[IPARAM_P]);
        exit(2);
    }
    if(0 == iparam[IPARAM_Q])
        iparam[IPARAM_Q] = iparam[IPARAM_NNODES] / iparam[IPARAM_P];
    int pqnp = iparam[IPARAM_Q] * iparam[IPARAM_P];
    if(pqnp > iparam[IPARAM_NNODES])
    {
        fprintf(stderr, "#XXXXX the process grid PxQ (%dx%d) is larger than the number of nodes (%d)!\n", iparam[IPARAM_P], iparam[IPARAM_Q], iparam[IPARAM_NNODES]);
        exit(2);
    }
    if(verbose && (pqnp < iparam[IPARAM_NNODES]))
    {
        fprintf(stderr, "#!!!!! the process grid PxQ (%dx%d) is smaller than the number of nodes (%d). Some nodes are idling!\n", iparam[IPARAM_P], iparam[IPARAM_Q], iparam[IPARAM_NNODES]);
    }

    if(0 == iparam[IPARAM_N])
    {
        fprintf(stderr, "#XXXXX the matrix size (N) is not set!\n");
        exit(2);
    }
    if(0 == iparam[IPARAM_M]) iparam[IPARAM_M] = iparam[IPARAM_N];
    if(0 == iparam[IPARAM_K]) iparam[IPARAM_K] = iparam[IPARAM_N];

    /* Set some sensible defaults for the leading dimensions */
    if(-'m' == iparam[IPARAM_LDA]) iparam[IPARAM_LDA] = iparam[IPARAM_M];
    if(-'n' == iparam[IPARAM_LDA]) iparam[IPARAM_LDA] = iparam[IPARAM_N];
    if(-'k' == iparam[IPARAM_LDA]) iparam[IPARAM_LDA] = iparam[IPARAM_K];
    if(-'m' == iparam[IPARAM_LDB]) iparam[IPARAM_LDB] = iparam[IPARAM_M];
    if(-'n' == iparam[IPARAM_LDB]) iparam[IPARAM_LDB] = iparam[IPARAM_N];
    if(-'k' == iparam[IPARAM_LDB]) iparam[IPARAM_LDB] = iparam[IPARAM_K];
    if(-'m' == iparam[IPARAM_LDC]) iparam[IPARAM_LDC] = iparam[IPARAM_M];
    if(-'n' == iparam[IPARAM_LDC]) iparam[IPARAM_LDC] = iparam[IPARAM_N];
    if(-'k' == iparam[IPARAM_LDC]) iparam[IPARAM_LDC] = iparam[IPARAM_K];

    /* Set no defaults for IB, NB, MB, the algorithm have to do it */
    assert(iparam[IPARAM_IB]); /* check that defaults have been set */
    if(iparam[IPARAM_NB] <= 0 && iparam[IPARAM_MB] > 0) iparam[IPARAM_NB] = iparam[IPARAM_MB];
    if(iparam[IPARAM_MB] <= 0 && iparam[IPARAM_NB] > 0) iparam[IPARAM_MB] = iparam[IPARAM_NB];
    if(iparam[IPARAM_NGPUS] && iparam[IPARAM_MB] < 0) iparam[IPARAM_MB] = -384;
    if(iparam[IPARAM_NGPUS] && iparam[IPARAM_NB] < 0) iparam[IPARAM_NB] = -384;
    if(iparam[IPARAM_MB] < 0) iparam[IPARAM_MB] = -iparam[IPARAM_MB];
    if(iparam[IPARAM_NB] < 0) iparam[IPARAM_NB] = -iparam[IPARAM_NB];

    /* No supertiling by default */
    if(-'p' == iparam[IPARAM_KP]) iparam[IPARAM_KP] = (iparam[IPARAM_M]/iparam[IPARAM_MB])/iparam[IPARAM_P];
    if(-'q' == iparam[IPARAM_KQ]) iparam[IPARAM_KQ] = (iparam[IPARAM_N]/iparam[IPARAM_NB])/iparam[IPARAM_Q];
    if(0 == iparam[IPARAM_KP]) iparam[IPARAM_KP] = 1;
    if(0 == iparam[IPARAM_KQ]) iparam[IPARAM_KQ] = 1;
    if(0 == iparam[IPARAM_HMB]) iparam[IPARAM_HMB] = iparam[IPARAM_MB];
    if(0 == iparam[IPARAM_HNB]) iparam[IPARAM_HNB] = iparam[IPARAM_NB];

    /* HQR */
    if(-'P' == iparam[IPARAM_QR_HLVL_SZE]) iparam[IPARAM_QR_HLVL_SZE] = iparam[IPARAM_P];
    if(-'Q' == iparam[IPARAM_QR_HLVL_SZE]) iparam[IPARAM_QR_HLVL_SZE] = iparam[IPARAM_Q];
}

static void print_arguments(int* iparam)
{
    int verbose = iparam[IPARAM_RANK] ? 0 : iparam[IPARAM_VERBOSE];

    if(verbose)
        fprintf(stderr, "#+++++ cores detected       : %d\n", iparam[IPARAM_NCORES]);

    if(verbose > 1) fprintf(stderr, "#+++++ nodes x cores + gpu  : %d x %d + %d (%d+%d)\n"
                                    "#+++++ thread mode          : %s\n"
                                    "#+++++ P x Q                : %d x %d (%d/%d)\n",
                            iparam[IPARAM_NNODES],
                            iparam[IPARAM_NCORES],
                            iparam[IPARAM_NGPUS],
                            iparam[IPARAM_NNODES] * iparam[IPARAM_NCORES],
                            iparam[IPARAM_NNODES] * iparam[IPARAM_NGPUS],
                            iparam[IPARAM_THREAD_MT]? "THREAD_MULTIPLE": "THREAD_SERIALIZED",
                            iparam[IPARAM_P], iparam[IPARAM_Q],
                            iparam[IPARAM_Q] * iparam[IPARAM_P], iparam[IPARAM_NNODES]);

    if(verbose)
    {
        fprintf(stderr, "#+++++ M x N x K|NRHS       : %d x %d x %d\n",
                iparam[IPARAM_M], iparam[IPARAM_N], iparam[IPARAM_K]);
    }

    if(verbose > 2)
    {
        if(iparam[IPARAM_LDB] && iparam[IPARAM_LDC])
            fprintf(stderr, "#+++++ LDA , LDB , LDC      : %d , %d , %d\n", iparam[IPARAM_LDA], iparam[IPARAM_LDB], iparam[IPARAM_LDC]);
        else if(iparam[IPARAM_LDB])
            fprintf(stderr, "#+++++ LDA , LDB            : %d , %d\n", iparam[IPARAM_LDA], iparam[IPARAM_LDB]);
        else
            fprintf(stderr, "#+++++ LDA                  : %d\n", iparam[IPARAM_LDA]);
    }

    if(verbose)
    {
        if(iparam[IPARAM_IB] > 0)
            fprintf(stderr, "#+++++ MB x NB , IB         : %d x %d , %d\n",
                            iparam[IPARAM_MB], iparam[IPARAM_NB], iparam[IPARAM_IB]);
        else
            fprintf(stderr, "#+++++ MB x NB              : %d x %d\n",
                    iparam[IPARAM_MB], iparam[IPARAM_NB]);
        if(iparam[IPARAM_KQ] * iparam[IPARAM_KP] != 1)
            fprintf(stderr, "#+++++ KP x KQ              : %d x %d\n", iparam[IPARAM_KP], iparam[IPARAM_KQ]);
        if((iparam[IPARAM_HNB] != iparam[IPARAM_NB]) || (iparam[IPARAM_HMB] != iparam[IPARAM_MB]))
            fprintf(stderr, "#+++++ HMB x HNB            : %d x %d\n", iparam[IPARAM_HMB], iparam[IPARAM_HNB]);
    }
}




static void iparam_default(int* iparam)
{
    /* Just in case someone forget to add the initialization :) */
    memset(iparam, 0, IPARAM_SIZEOF * sizeof(int));
    iparam[IPARAM_NNODES] = 1;
    iparam[IPARAM_ASYNC]  = 1;
    iparam[IPARAM_QR_DOMINO]    = -1;
    iparam[IPARAM_LOWLVL_TREE]  = DPLASMA_GREEDY_TREE;
    iparam[IPARAM_HIGHLVL_TREE] = -1;
    iparam[IPARAM_QR_TS_SZE]    = -1;
    iparam[IPARAM_QR_HLVL_SZE]  = -'P';
}

void iparam_default_ibnbmb(int* iparam, int ib, int nb, int mb)
{
    iparam[IPARAM_IB] = ib ? ib : -1;
    iparam[IPARAM_NB] = -nb;
    iparam[IPARAM_MB] = -mb;
}


void iparam_default_facto(int* iparam)
{
    iparam_default(iparam);
    iparam[IPARAM_K] = 1;
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = 0;
    iparam[IPARAM_LDC] = 0;
}

void iparam_default_solve(int* iparam)
{
    iparam_default(iparam);
    iparam[IPARAM_K] = 1;
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'n';
    iparam[IPARAM_LDC] = 0;
    iparam[IPARAM_M] = -'n';
}

void iparam_default_gemm(int* iparam)
{
    iparam_default(iparam);
    iparam[IPARAM_K] = 0;
    /* no support for transpose yet */
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'k';
    iparam[IPARAM_LDC] = -'m';
}

#ifdef PARSEC_PROF_TRACE
static char* argvzero;
char cwd[1024];
int unix_timestamp;
#endif

#if defined(DPLASMA_HAVE_CUDA)
static void destroy_cuda_handles(void *_h, void *_n)
{
    dplasma_cuda_handles_t *handles = (dplasma_cuda_handles_t*)_h;
    (void)_n;
    cublasDestroy_v2(handles->cublas_handle);
    cusolverDnDestroy(handles->cusolverDn_handle);
    free(handles);
}
#endif

parsec_context_t* setup_parsec(int argc, char **argv, int *iparam)
{
#ifdef PARSEC_PROF_TRACE
    argvzero = argv[0];
    unix_timestamp = time(NULL);
    getcwd(cwd, sizeof(cwd));
#endif
    read_arguments(&argc, &argv, iparam);
#ifdef PARSEC_HAVE_MPI
    {
        int requested = iparam[IPARAM_THREAD_MT]? MPI_THREAD_MULTIPLE: MPI_THREAD_SERIALIZED;
        int provided;
        MPI_Init_thread(&argc, &argv, requested, &provided);
        if( requested > provided ) {
            fprintf(stderr, "#XXXXX User requested %s but the implementation returned a lower thread\n", requested==MPI_THREAD_MULTIPLE? "MPI_THREAD_MULTIPLE": "MPI_THREAD_SERIALIZED");
            exit(2);
        }
    }
    MPI_Comm_size(MPI_COMM_WORLD, &iparam[IPARAM_NNODES]);
    MPI_Comm_rank(MPI_COMM_WORLD, &iparam[IPARAM_RANK]);
#else
    iparam[IPARAM_NNODES] = 1;
    iparam[IPARAM_RANK] = 0;
#endif
    parse_arguments(iparam);
    int verbose = iparam[IPARAM_VERBOSE];
    if(iparam[IPARAM_RANK] > 0 && verbose < 4) verbose = 0;

    TIME_START();

    /* Once we got out arguments, we should pass whatever is left down */
    int parsec_argc, idx;
    char** parsec_argv = (char**)calloc(argc, sizeof(char*));
    parsec_argv[0] = argv[0];  /* the app name */
    for( idx = parsec_argc = 1;
         (idx < argc) && (0 != strcmp(argv[idx], "--")); idx++);
    if( idx != argc ) {
        for( parsec_argc = 1, idx++; idx < argc;
             parsec_argv[parsec_argc] = argv[idx], parsec_argc++, idx++);
    }
    parsec_context_t* ctx = parsec_init(iparam[IPARAM_NCORES],
                                        &parsec_argc, &parsec_argv);
    free(parsec_argv);
    if( NULL == ctx ) {
        /* Failed to correctly initialize. In a correct scenario report
         * upstream, but in this particular case bail out.
         */
        exit(-1);
    }

    /* In order to find the help-dplasma.txt file both in case
     * this is run within the dplasma source trunk or after
     * installation, we add both paths to the show_help directories.
     * We add the install one first so that in normal operations (installed),
     * the number of file access is minimized. */
    parsec_show_help_add_dir(DPLASMA_SHARE_DIRECTORY_INSTALL);
    parsec_show_help_add_dir(DPLASMA_SHARE_DIRECTORY_SOURCE);

    /* If the number of cores has not been defined as a parameter earlier
     update it with the default parameter computed in parsec_init. */
    if(iparam[IPARAM_NCORES] <= 0)
    {
        int p, nb_total_comp_threads = 0;
        for(p = 0; p < ctx->nb_vp; p++) {
            nb_total_comp_threads += ctx->virtual_processes[p]->nb_cores;
        }
        iparam[IPARAM_NCORES] = nb_total_comp_threads;
    }
    print_arguments(iparam);

#if defined(DPLASMA_HAVE_CUDA)
    int dev, nbgpu = 0;
    for(dev = 0; dev < (int)parsec_nb_devices; dev++) {
        parsec_device_module_t *device = parsec_mca_device_get(dev);
        if( PARSEC_DEV_CUDA == device->type ) {
            nbgpu++;
        }
    }
    if( nbgpu > 0 ) {
        cublasStatus_t status = cublasInit();
        assert(CUBLAS_STATUS_SUCCESS == status);
        parsec_info_register(&parsec_per_stream_infos, "DPLASMA::CUDA::HANDLES",
                             destroy_cuda_handles, NULL,
                             dplasma_create_cuda_handles, NULL,
                             NULL);
    }
#endif

    if(verbose > 2) TIME_PRINT(iparam[IPARAM_RANK], ("PaRSEC initialized\n"));
    return ctx;
}

void cleanup_parsec(parsec_context_t* parsec, int *iparam)
{
#if defined(DPLASMA_HAVE_CUDA)
    parsec_info_id_t CuHI = parsec_info_lookup(&parsec_per_stream_infos, "DPLASMA::CUDA::HANDLES", NULL);
    parsec_info_unregister(&parsec_per_stream_infos, CuHI, NULL);
    cublasShutdown();
#endif

    parsec_fini(&parsec);

#ifdef PARSEC_HAVE_MPI
    MPI_Finalize();
#endif
    (void)iparam;
}

void dplasma_warmup(parsec_context_t *parsec)
{
    int Aseed = 3872;
    int Bseed = 4674;
    int Cseed = 2873;
    int tA = dplasmaNoTrans;
    int tB = dplasmaNoTrans;

    // DPLASMA might have been compiled with only one precision, or any subset of the 4 possible
    // precisions. The following logic tries to find /a/ kernel we can use. We check in the
    // arbitrary order d, s, c, z, but really we just want any of them.
#if defined(DPLASMA_DGEMM_NN)
    double alpha =  0.51;
    double beta  = -0.42;
    parsec_matrix_type_t mtype = PARSEC_MATRIX_DOUBLE;
#define KERNEL_NEW      dplasma_dgemm_New
#define KERNEL_DESTRUCT dplasma_dgemm_Destruct
#define INIT_MATRIX     dplasma_dplrnt
#elif defined(DPLASMA_SGEMM_NN)
    float alpha =  0.51;
    float beta  = -0.42;
    parsec_matrix_type_t mtype = PARSEC_MATRIX_FLOAT;
#define KERNEL_NEW      dplasma_sgemm_New
#define KERNEL_DESTRUCT dplasma_sgemm_Destruct
#define INIT_MATRIX     dplasma_splrnt
#elif defined(DPLASMA_CGEMM_NN)
    dplasma_complex32_t alpha =  0.51 + I * 0.32;
    dplasma_complex32_t beta  = -0.42 + I * 0.21;
    parsec_matrix_type_t mtype = PARSEC_MATRIX_COMPLEX_FLOAT;
#define KERNEL_NEW      dplasma_cgemm_New
#define KERNEL_DESTRUCT dplasma_cgemm_Destruct
#define INIT_MATRIX     dplasma_cplrnt
#elif defined(DPLASMA_ZGEMM_NN)
    dplasma_complex64_t alpha =  0.51 + I * 0.32;
    dplasma_complex64_t beta  = -0.42 + I * 0.21;
    parsec_matrix_type_t mtype = PARSEC_MATRIX_COMPLEX_DOUBLE;
#define KERNEL_NEW      dplasma_zgemm_New
#define KERNEL_DESTRUCT dplasma_zgemm_Destruct
#define INIT_MATRIX     dplasma_zplrnt
#else
#warning "DPLASMA is configured without any of the sdcz precisions... Warmup will be no-op"
#endif

#if defined(KERNEL_NEW)
    int M, N, K;
    int MB;
    int LDA, LDB, LDC;
    int P, Q;

    int rank = parsec->my_rank;
    int nodes = parsec->nb_nodes;

    int gpus = 0;

#if defined(DPLASMA_HAVE_CUDA)
    int devid;
    for(devid = 0; devid < (int)parsec_nb_devices; devid++) {
        parsec_device_module_t *device = parsec_mca_device_get(devid);
        if( PARSEC_DEV_CUDA == device->type ) {
            gpus++;
        }
    }
#endif

    if(0 == gpus) {
        MB = 512;
        N = nodes * parsec->virtual_processes[0]->nb_cores * 3 * MB;
        P = nodes;
        Q = 1;
    } else {
        MB = 64;
        N = nodes * gpus * 3 * MB;
        P = nodes;
        Q = 1;
    }
    M = MB;
    K = MB;

    LDA = MB;
    LDB = MB;
    LDC = MB;

    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        parsec_matrix_block_cyclic, (&dcC, mtype, PARSEC_MATRIX_TILE,
                               rank, MB, MB, LDC, N, 0, 0,
                               M, N, nodes, P, Q, 1, 0, 0));

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, mtype, PARSEC_MATRIX_TILE,
                                rank, MB, MB, LDA, K, 0, 0,
                                M, K, nodes, P, Q, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcB, 1,
        parsec_matrix_block_cyclic, (&dcB, mtype, PARSEC_MATRIX_TILE,
                                rank, MB, MB, LDB, N, 0, 0,
                                K, N, nodes, P, Q, 1, 0, 0));

    /* matrix generation */
    INIT_MATRIX( parsec, 0, (parsec_tiled_matrix_t *)&dcA, Aseed);
    INIT_MATRIX( parsec, 0, (parsec_tiled_matrix_t *)&dcB, Bseed);
    INIT_MATRIX( parsec, 0, (parsec_tiled_matrix_t *)&dcC, Cseed);

    parsec_devices_release_memory();
    parsec_taskpool_t* PARSEC_gemm = KERNEL_NEW (tA, tB, alpha, (parsec_tiled_matrix_t *)&dcA, (parsec_tiled_matrix_t *)&dcB, beta, (parsec_tiled_matrix_t *)&dcC); 
    PARSEC_CHECK_ERROR(parsec_context_add_taskpool(parsec, PARSEC_gemm), "parsec_context_add_taskpool"); 
    PARSEC_CHECK_ERROR(parsec_context_start(parsec), "parsec_context_start"); 
    PARSEC_CHECK_ERROR(parsec_context_wait(parsec), "parsec_context_wait"); 
    KERNEL_DESTRUCT ( PARSEC_gemm );
    parsec_devices_reset_load(parsec);

    parsec_data_free(dcA.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    parsec_data_free(dcB.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);
    parsec_data_free(dcC.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC);

#endif
}
