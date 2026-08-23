Attribute VB_Name = "Module4"

    


'branch.h
'#ifndef __BRANCH
'#define __BRANCH

'/* Structure and functions concerned with setup, prior, and proposal on
' * branches.  Currently only handles integrated branch lengths with a single
' * parameter, the expected branch length.
' */

#include <stdarg.h>
#include "settings.h"
#include "sampler.h"    ' Basic probability calculations

struct _branch {
Type branch
    mu_mean As Double                     ' Hierarchical mean of average branch length
    mu_variance As Double                 ' Hierarchical variance of average branch length
    double (*Log_Prior) (const branch *, double);       ' Function to compute log prior of current average branch length (arg 2)
    double (*Propose) (settings *, ...);            ' Function to propose a new average branch length
    double (*Propose_Merge) (double *, double, double, settings *, ...);    ' Function to merge two segments with two separate average branch lengths
    double (*Propose_Split) (double *, double *, double, settings *, ...);  ' Function to propose split of one segment into two with separate average branch lengths
    double (*Initialize) (settings *set);                   ' Function to return an initial value for average branch length
End Type 'branch
};

' Global function declarations: Create and destroy
void BranchMake(branch **, settings *);
void BranchDelete(branch *);

''#End If

'branch.c
#include "branch.h"

' Following functions are the first and only fulfillment of the function descriptions made in header file.
' Alternative functions could be used an linked in instead.
' These functions are made externally available via the function pointers.
static double Initialize(settings *);                       ' Assume initial distribution is standard uniform (seems a little under-dispersed)!
static double Propose(settings *, ...);                     ' Proposal distribution is c*exp(U - 0.5( where U is standard uniform (bigger moves for bigger values)
static double HierarchicalLogPrior(const branch *, double);         ' Current prior is iid against hierarchical log normal
static double Propose_Merge(double *, double, double, settings *, ... );    ' Propose merge when deleting parameter change point (weighted average on log-scale)
static double Propose_Split(double *, double *, double, settings *, ... );  ' Propose split when adding parameter change point (inverse Green rjMCMC move of delete)

void BranchMake(branch **br, settings *set) {
    *br = (branch *)malloc(sizeof(branch));
    (*br)->mu_mean = set->mu_hyper_mean;
    (*br)->mu_variance = set->mu_hyper_variance;
    (*br)->Log_Prior = &HierarchicalLogPrior;
    (*br)->Initialize = &Initialize;
    (*br)->Propose = &Propose;
    (*br)->Propose_Merge = &Propose_Merge;
    (*br)->Propose_Split = &Propose_Split;
}' BranchMake

'static double Initialize(settings *set) {
Public Function Initialize(setx As settings) As Double
    'return set->rng->nextStandardUniform(set->rng);
    Initialize = nextStandardUniform(setx.rng)

End Function '}' Initialize
   
static double Propose(settings *set, ...) {
    va_list vargs;
    Dim current_value As Double, lambda As Double, ret As Double

    va_start(vargs, set);
    current_value = va_arg(vargs, double);
    lambda = va_arg(vargs, double);
    va_end(vargs);

    ret = current_value*exp(lambda*(set->rng->nextStandardUniform(set->rng) - 0.5));
    return ret;
}' Propose

static double Propose_Merge(double *pMu, double cMu1, double cMu2, settings *set, ... ) {
    Dim left As Long, middle As Long, right As Long
    Dim weight1 As Double, weight2 As Double, zMu As Double
    Dim vargs As va_list

    va_start(vargs, set)
    left = va_arg(vargs, int);
    middle = va_arg(vargs, int);
    right = va_arg(vargs, int);
    va_end(vargs);

    weight1 = (middle - left) / (right - left)
    weight2 = (right - middle) / (right - left)

    *pMu = exp( weight1 * log(cMu1) + weight2 * log(cMu2) )
    zMu = log( cMu1 / *pMu ) / weight2 / set->sigmaMu
    return logStandardNormalDensity(zMu);
}' Propose_Merge

static double Propose_Split(double *pMu1, double *pMu2, double cMu, settings *set, ... ) {
    Dim left As Long, middle As Long, right As Long
    Dim weight1 As Double, weight2 As Double, zMu As Double
    Dim vargs As va_list

    va_start(vargs, set);
    left = va_arg(vargs, int);
    middle = va_arg(vargs, int);
    right = va_arg(vargs, int);
    va_end(vargs);

    weight1 = (double) (middle - left) / (right - left);
    weight2 = (double) (right - middle) / (right - left);

    zMu = set->rng->nextStandardNormal(set->rng);
    *pMu1 = cMu * exp( weight2 * set->sigmaMu * zMu );
    *pMu2 = cMu * exp( -weight1 * set->sigmaMu * zMu );
    return logStandardNormalDensity(zMu);
}' Propose_Split

static double HierarchicalLogPrior(const branch *br, double mu) {
    Dim log_mu As Double
    log_mu = Log(mu)
    return -1.0*(log_mu - br->mu_mean)*(log_mu - br->mu_mean)/2.0/br->mu_variance - log_mu - 0.5*log(br->mu_variance) - normal_const;
}' HierarchicalLogPrior

void BranchDelete(branch *br) {
    if(br) free(br);
}' BranchDelete

'cbrother.c

#include "constants.h"
#include "settings.h"
#include "seqdata.h"
#include "cpsampler.h"
#include "dcpsampler.h"
'#include "multistddata.h"

' GLOBAL VARIABLE DEFINITIONS
int global_debug = 0;
boolean compute_likelihood = true;

static const char *prog_name = "cbrother";
static int debug = 0;

void PrintCorrectUsage() {

    fprintf(stderr, "ERROR: Invalid usage.  Correct usage below.\n");
    fprintf(stderr, "specify model:<recomb|diverge|mcp_recomb> in cmdfile\n\n");
    fprintf(stderr, "recomb:\n\t%s <seed> <cmdfile> <phylipfile> <postfile>\n", prog_name);
    fprintf(stderr, "diverge: \n\t%s <seed> <cmdfile> <phylipfile1> <phylipfile2> ... <postfile>\n", prog_name);
    exit(EXIT_FAILURE);
}

'int main(int argc, char * argv[]) {
Public Sub DualBrother(argc As Long, argv() As Byte)


    Dim setx As settings
    Dim sqd As seqdata    ' sequences
    Dim seed As Byte, cmdfile As String, postfile As Byte
    ' not implemented multistddata msd;

    'if( argc < 3 ) PrintCorrectUsage();
    
    seed = argv(1)
    cmdfile = argv(2)
    postfile = argv(argc - 1)

    'if( debug || global_debug ) fprintf(stderr, "Posterior file: %s\n", postfile);

    Call ReadCmdfile(setx, cmdfile)
    
    'Set_Seed(setx ->rng, setx ->cmdfile_seed ? setx ->seed : (unsigned) atol(seed));  ' WARNING: Ignore command-line seed if there is a seed setx  in cmdfile
    If setx.cmdfile_seed = 1 Then
        Call Set_Seed(setx.rng, setx.seed)   ' WARNING: Ignore command-line seed if there is a seed setx  in cmdfile
    Else
        Call Set_Seed(setx.rng, atol(seed))
    End If
    
    
    'compute_likelihood = setx ->compute_likelihood || compute_likelihood;
    
    If compute_likelihood = 1 Or setx.compute_likelihood = 1 Then compute_likelihood = 1
    
    'global_debug = setx ->debug ? setx ->debug : global_debug;
    If setx.debug = 1 Then
        global_debug = 1
    Else
        global_debug = global_debug
    End If
    
    If setx.simulate_data = 0 Then      ' Read Phylip File(s)

        'if( argc != 5 ) PrintCorrectUsage();

        'ReadPhylip(&sqd, argv[3], false, setx ->model);
        Call ReadPhylip(sqd(), argv(3), 0, setx.model)
    Else
        'if( argc != 4 ) PrintCorrectUsage();
        'SimulateAlignment(&sqd, setx );
        Call SimulateAlignment(sqd, setx)
        'PrintSequences(sqd);
        'Print_Distances(sqd);
    End If

       ' Make data objects
       
       ' not implemented if (nfamilies > 1) Make_MultiStdData(&msd, sqd, nfamilies, ntaxa);
              
    If setx.model = SCP_RECOMB Then ' Run MCMC
        'cpsampler *scp;
        Dim scp As cpsampler
        'CPSamplerSetup(&scp, sqd, setx , postfile);
        Call CPSamplerSetup(scp, sqd, setx, postfile)
        'scp.smp.run(scp.smp);
        CPSamplerRun (scp.smp)
        'CPSamplerDelete(&scp);
        Call CPSamplerDelete(scp)
    End If
    ' Not implemented yet
    'else if (setx ->model == DIVERGE) { ' Run Gu
    '  Setup_Gu(&msd, setx , postfile);
    '  RunGu();
    '}
    else if( setx ->model == DCP_RECOMB ) {
        dcpsampler *dcp = NULL;
        DCPSamplerSetup(&dcp, sqd, setx , postfile);
        dcp->smp->run(dcp->smp);
        'DCPSamplerDelete(&dcp);
    }
    
    
    ' Cleanup
    SeqDataDelete(sqd);
    ' Not implemented yet
    'else if(setx ->model == DIVERGE) {
    '  MultiStdDataDelete(&msd);
    '  GuDelete();
    '}
    
    Settings_Cleanup(setx );
    return 1;
    
End Sub '} 'main



    ' to make stddata one file at a time
       /******************************
    {
        int i, nfamilies = 0, nchars = 0;
        if( setx .model == DIVERGE ) {
            nfamilies = argc - 3;
            nchars = 20;
            m_std = malloc (sizeof(stddata) * nfamilies);
        }
        else if( set.model == RECOMB ) {
            nfamilies = 1;
            nchars = 4;
        }
        
        for( i = 0; i < nfamilies; i++ ) {
            
            ReadPhylip(&sqd, argv[i+3], false, set.model);
            Make_StdData(&m_std[i], &sqd);
        }

    }
    *******************/
    
/*****************
' test main
int main() {
    
    seqdata sqd;
    stddata std;
    settings set;
    char *str1 = "ATAAAA";
    char *str2 = "AATATA";
    char *str3 = "AAAAAA";
    int *data1 = malloc(sizeof(int) * strlen(str1));
    int *data2 = malloc(sizeof(int) * strlen(str2));
    int *data3 = malloc(sizeof(int) * strlen(str3));
    sequence seqs[3];
    int length = strlen(str1);
    
    set.model = RECOMB;
    Setup_Sequence(&seqs[0], str1, "seq1", data1, RECOMB);
    Setup_Sequence(&seqs[1], str2, "seq2", data2, RECOMB);
    Setup_Sequence(&seqs[2], str3, "seq3", data3, RECOMB);
    
    Make_SeqData(&sqd, seqs, 3,length, 4);
    
'  PrintSequenceInfo(&seqs[0], false);
'  PrintSequenceInfo(&seqs[1], false);
'  PrintSequenceInfo(&seqs[2], false);
    
    Make_StdData(&std, &sqd);
    PrintCounts(sqd);
    PrintMap(sqd);
    PrintRMap(sqd);
    printf("sequences:\n");
    PrintSequences(sqd);
    printf("\nsorted seq:\n");
    PrintSortedSequences(sqd);
    printf("\ncompressed:\n");
    PrintCompressedSequences(sqd);
    'PrintSeqData(&sqd);
        

return 1;
}
************/



' cpsampler.c

#include "cpsampler.h"

'#define SWAP_SIZE 8

static int debug = 1;

' TEMPORARY
static int total_adds = 0;
static int total_deletes = 0;
static double poisson[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
static double avg_nchpts = 0.0;
static double std_nchpts = 0.0;
static double avg_mu = 0.0;
static double std_mu = 0.0;
static double avg_alpha = 0.0;
static double std_alpha = 0.0;
static int sample_nchpts = 0;
static int sample_poisson = 0;
static int sample_param = 0;
static void CPOutputLine(const sampler *);
static void CPSamplerRun(sampler *);
static void CPRecordStatistics(cpsampler *);
static void CPPrintTopologies(sampler *, const char *);

' Setup functions
static void CPInitialState(cpsampler *);
static void CPSetParameters(cpsampler *);
static void CPSetPostTreeMultipleQueries(cpsampler *);

' Proposal mechanisms
static void ProposeNewTopologies(cpsampler *, partition **, const int, const double, const boolean);
static void ProposeNewQandMu(cpsampler *, partition **, const int, const double, const boolean);
static void CPProposeNewXi(cpsampler *, partition **, const int, const double, const boolean);
static void ProposeAddSegment(cpsampler *);
static void ProposeDeleteSegment(cpsampler *);

' Auxiliary functions
static int determineSameT(const int, partition *[]);
static int updateSameT(const cpsampler *, const int, const int);
static double CPLogJacobian(double, double, double, double, double, double, double, double);
static void CPSRemovePartition(cpsampler *, int);
static void CPSAddPartition(cpsampler *, partition *, int);
static double DrawnAlpha1(double, double, double, double);
static double DrawnAlpha2(double, double, double, double);
static double DrawnMu1(double, double, double, double);
static double DrawnMu2(double, double, double, double);
static double CondensedAlpha(double, double, double, double);
static double CondensedMu(double, double, double, double);
static double InverseZMu(double, double, double, double);
static double InverseZAlpha( double, double, double, double);
static double Logit(const double);
static double InvLogit(const double);

' I/O
static void PrintInitialValues(const cpsampler *);


' Extra space allocated to avoid constant reallocation (calculation space)
' TODO: how to use this memory pool with the new setup?
/*
static double *memory_PartialLogLikelihood;
static double *memory_PartialLogHyperParameterPrior;
static double *memory_HyperParameter;
static double *memory_Alpha;
static int **memory_counts;
static int *memory_Tree;
static int *memory_BreakPoints;
static int *memory_testcounts;
static boolean *memory_doUpdate;
*/

static void CPParseLastLine(void *smp) {
    char *tmp;
    cpsampler *cps = (cpsampler *) smp;
    char *str = cps->set->init_string;
    int count = 0, current_region = 0, var_num = 0, i;
    partition *part;
    
    tmp = strtok(str, " ");

    while( tmp != NULL ) {
        ' Count tokens
        count++;
        ' First word is the sample number
        if( count == 1) cps->smp->JumpNumber = atoi(tmp);
        ' Second word is the number of segments
        else if( count == 2 ) {
            cps->npartitions = atoi(tmp);
            ' Allocate memory for all other parameters
            cps->part_list = (partition **) malloc(sizeof(partition *)*cps->npartitions);
            part = cps->part_list[0];
            PartitionMake(&part, cps->lenunique, 0, 0, false, false);
            ' TODO: make independent of matrix type
            iHKYNoBoundFixPiMatrixMakeDefault(&part->cmatrix, ALPHA);
        }
        ' The number of evolutionary parameter change-points
        else if( count == 3 ) {
            cps->cSameT = atoi(tmp);
        }
        ' Segment-specific parameters
        else if( count <= 3 + cps->npartitions * 8 ) {
            ' topology
            if( var_num % 8 == 0 ) {
                tree *ttree;
                Make_Tree(&ttree, tmp, cps->sqd->num_chars);
                Balance_Tree(ttree);
                ' Find the matching tree in the list of possibilities
                for(i=0; i<cps->numTrees; i++) {
                    if( !SameTrees(&cps->PostTree[i], ttree, false) ) {
                         part->cTree = i;
                         break;
                    }
                }
                if(ttree) free(ttree);
            ' log likelihood
            } else if( var_num % 8 == 1 ) {
                part->cPartialLogLikelihood = atof(tmp);
                cps->cLogLikelihood += part->cPartialLogLikelihood;
            ' pi_A
            } else if ( var_num % 8 == 2 ) {
                part->cmatrix->pi[0] = atof(tmp);
            ' pi_C
            } else if ( var_num % 8 == 3 ) {
                part->cmatrix->pi[1] = atof(tmp);
            ' pi_G
            } else if ( var_num % 8 == 4 ) {
                part->cmatrix->pi[2] = atof(tmp);
            ' pi_T
            } else if ( var_num % 8 == 5 ) {
                part->cmatrix->pi[3] = atof(tmp);
            ' alpha
            } else if ( var_num % 8 == 6 ) {
                Set_Alpha(part->cmatrix, atof(tmp));
            ' mu
            } else if ( var_num % 8 == 7 ) {
                part->cHyperParameter = atof(tmp);
                part->cPartialLogHyperParameterPrior = - part->cHyperParameter;
                current_region++;
                if( current_region < cps->npartitions )
                    part = cps->part_list[current_region];
                    PartitionMake(&part, cps->lenunique, 0, 0, false, false);
                    ' TODO: make independent of matrix type
                    iHKYNoBoundFixPiMatrixMakeDefault(&part->cmatrix, ALPHA);
            }
            var_num++;
        }
        ' The last words are the change-point locations
        else {
            if(current_region >= cps->npartitions) {
                 current_region = 0;
                 part = cps->part_list[current_region];
            }
            part->left = atoi(tmp);
            if(current_region) {
                 cps->part_list[current_region-1]->right = part->left - 1;
                 PartitionCopySegmentCounts(cps->part_list[current_region-1], cps->sqd, cps->part_list[current_region-1]->left, cps->part_list[current_region-1]->right+1);
            }
            current_region++;
            part = cps->part_list[current_region];
        }
        ' Get next token
        tmp = strtok(NULL, " ");
    }
    ' Set the last segment
    part->right = cps->lenSeq - 1;
    PartitionCopySegmentCounts(part, cps->sqd, part->left, cps->lenSeq);
    cps->cSameT = determineSameT(cps->npartitions, cps->part_list);
} ' CPParseLastLine

static void CPInitialState(cpsampler *cps) {
    partition *part = NULL;
    ' BUG: assume only one partition

    if( cps->set->debug>2 || debug>2 ) fprintf(stderr, "Enter Initial State\n");
    
    /*memory_PartialLogLikelihood = (double *) malloc(sizeof(double)*SWAP_SIZE);
    memory_PartialLogHyperParameterPrior = (double *) malloc(sizeof(double)*SWAP_SIZE);
    memory_HyperParameter = (double *) malloc(sizeof(double)*SWAP_SIZE);
    memory_Alpha = (double *) malloc(sizeof(double)*SWAP_SIZE);
    memory_Tree = (int *) malloc(sizeof(int)*SWAP_SIZE);
    memory_BreakPoints = (int *) malloc(sizeof(int)*SWAP_SIZE);
    memory_counts = (int **) malloc(sizeof(int *)*SWAP_SIZE);
    memory_testcounts = (int *) malloc(sizeof(int)*lenunique);
    memory_doUpdate = (boolean *) malloc(sizeof(boolean)*SWAP_SIZE);
    for(i=0; i<SWAP_SIZE; i++) memory_counts[i] = (int *) malloc(sizeof(int)*lenunique);
    */

    PartitionMake(&(cps->part_list[0]), cps->lenunique, 0, cps->lenSeq - 1, true, true);
    part = cps->part_list[0];
    PartitionCopyCounts(part, cps->sqd);
    ' Select random starting parameters
    part->cTree = cps->set->rng->nextStandardUniform(cps->set->rng)*cps->numTrees;'ChooseStartingTree(cps, 0);
    part->cHyperParameter = - log( cps->set->rng->nextStandardUniform(cps->set->rng) );
    part->cPartialLogHyperParameterPrior = - part->cHyperParameter;
    cps->set->ctmc_model = HKY;
    cps->set->ctmc_parameterization = ALPHA;
    ' TODO: make independent of matrix type
    iHKYNoBoundFixPiMatrixMakeInitial(&(part->cmatrix), ALPHA, cps->sqd, cps->set);
    
    part->cPartialLogLikelihood = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[part->cTree], cps->smp, part->cmatrix, part->counts, part->cHyperParameter, false) : 0.0;
    cps->cLogLikelihood = part->cPartialLogLikelihood;
    cps->cSameT = determineSameT(cps->npartitions, cps->part_list);

    if( cps->set->debug>2 || debug>2 ) {
        PrintInitialValues(cps);
        fprintf(stderr, "Leaving Initial State\n");
    }
}' CPInitialState

static void CPSetParameters(cpsampler *cps) {
    settings *set = cps->set;

    cps->jumpClasses = set->jump_classes;
    cps->lambda = set->scp->lambda;
    cps->mix1 = set->mix1;

    cps->logLambda = log(cps->lambda);
    cps->logTwo = log(2.0);
    cps->sigmaA = set->sigmaAlpha;
    cps->sigmaM = set->sigmaMu;
    cps->C = set->C;

    cps->logWSameT = set->weight;   ' log(w)
    cps->logWNotSameT = 0;
    cps->wSameT  = exp(cps->logWSameT);
    cps->wNotSameT = (1.0 - cps->wSameT);
    if( cps->wNotSameT == 0 ) cps->wNotSameT = -cps->logWSameT;
    cps->wNotSameT /=  (double) cps->numTrees - 1.0 ;
    if( cps->logWNotSameT == 0 ) cps->logWNotSameT = log(cps->wNotSameT);
}' CPSetParameters


' for normal set up, use NULL cfileName, otherwise setup will continue from last line in cfileName
'void CPSamplerSetup(cpsampler **in_cps, seqdata *in_sqd, settings *in_set, char *ofilename) {
Public Sub CPSamplerSetup(in_cps As cpsampler, in_sqd As seqdata, in_set As settings, ofilename As String)


    'char tmp[MAX_TREE_STRING];  ' Potential BUG
    'strcpy(tmp, "");
    Dim tmp As String
    
    'sampler *smp;
    Dim smp As sampler
    
    'cpsampler *cps;
    Dim cps As cpsampler
    
    
    '(*in_cps) = (cpsampler *) malloc(sizeof(cpsampler));
    'cps = *in_cps;
    cps = in_cps
    
    
    'SamplerMake(&cps->smp);
    Call SamplerMake(cps.smp)
    'smp = cps->smp;
    smp = cps.smp
    'smp->derived_smp = (void *)cps;
    smp.derived_smp = cps
    smp->set = in_set;
    smp->sqd = in_sqd;

    smp->max_move_name_length = strlen("evolutionary parameter");
    SamplerSetNumberMoves(smp, 7);
    SamplerAddMoveName(smp, 0, "topology");
    SamplerAddMoveName(smp, 1, "evolutionary parameter");
    SamplerAddMoveName(smp, 2, "move changepoint");
    SamplerAddMoveName(smp, 3, "add changepoint");
    SamplerAddMoveName(smp, 4, "add breakpoint");
    SamplerAddMoveName(smp, 5, "delete changepoint");
    SamplerAddMoveName(smp, 6, "delete breakpoint");

    ' Sampler output
    smp->fout = fopen(ofilename, "w");
    if(!smp->fout) {
        fprintf(stderr, "ERROR:  Could not open output file %s\n", ofilename);
        exit(EXIT_FAILURE);
    }
    smp->OutputLine = CPOutputLine;
    smp->run = CPSamplerRun;

    cps->sqd = in_sqd;
    cps->set = smp->set;
    
    cps->lenSeq = cps->sqd->lenseq;
    cps->indexSeq = cps->sqd->map;
    cps->lenunique = cps->sqd->lenunique;

    Make_Tree(&cps->start_tree, cps->set->pTree[0], cps->sqd->num_chars);
    CPSetPostTreeMultipleQueries(cps);
    
    CPSetParameters(cps);
    smp->JumpNumber = 1;

    if( !cps->set->init_string ) {
        cps->part_list = (partition **) malloc(sizeof(partition *));
        PartitionMake(&cps->part_list[0], cps->lenunique, 0, cps->lenSeq - 1, true, true);
        PartitionCopyCounts(cps->part_list[0], cps->sqd);
        cps->npartitions = 1;

        CPInitialState(cps);
    }
    else {
        CPParseLastLine(cps);
    }

    if( smp->JumpNumber > cps->set->burnin )
        SamplerSaveEstimates(cps->smp, cps->set->length);
    
End Sub ' CPSamplerSetup

static void CPSetPostTreeMultipleQueries(cpsampler *cps) {
    int num_queries = cps->sqd->ntaxa - cps->start_tree->nleaves;
    int queries_added = 0;
    int numTreesPerEnumeration = 0;
    int numCurrentTrees = 0;
    int numNextTrees = 0;
    int numCurrentLeaves = 0;
    tree *NewPostTree = NULL;
    tree *tree_ptr = NULL;
    tree *current_tree = NULL;
    int i =0;
    char tmp[MAX_TREE_STRING];  ' Potential BUG
    
    if (num_queries < 1) {
        fprintf(stderr, "No query sequences, cannot run\n");
        exit(EXIT_FAILURE);
    }
    
    if( cps->set->debug>2 || debug>2 ) {
        fprintf(stderr, "\nnum seq: %d\n", cps->sqd->ntaxa);
        fprintf(stderr, "num total trees: %d\n", Number_All_Trees(cps->start_tree, num_queries));
    }

    ' Enumerate Single Initial Tree for first query
    numTreesPerEnumeration = Number_Parental_Trees(cps->start_tree); ' 2* nleaves - 3
    numCurrentTrees = numTreesPerEnumeration;
    numCurrentLeaves = cps->start_tree->nleaves;
    cps->numTrees = 2 * numCurrentLeaves - 3;   ' Number of unrooted trees
    cps->PostTree = NULL;'(tree *) malloc(sizeof(tree)*numCurrentTrees);

    EnumerateLastTaxon(&(cps->PostTree), cps->start_tree);
    queries_added++;
    
    while (queries_added < num_queries) {
        
        numCurrentLeaves++;
        numNextTrees = numCurrentTrees * (2 * numCurrentLeaves - 3);
        numTreesPerEnumeration = 2 * (numCurrentLeaves + 1) - 3; 'Number_Parental_Trees(current_tree);  (complete enumeration, ok to hard-code the number of trees)
        
        NewPostTree = (tree *) malloc(sizeof(tree)*numNextTrees);
        
        ' for each tree in PostTree EnumerateLastTaxon()
        for (i = 0; i < numCurrentTrees; i++) {
            current_tree = &cps->PostTree[i];
            tree_ptr = &NewPostTree[i*numTreesPerEnumeration];
            EnumerateLastTaxon(&tree_ptr, current_tree);
            
        }' for each current post tree

        free(cps->PostTree);
        
        ' set the new PostTree
        numCurrentTrees = numNextTrees;
        
        cps->PostTree = NewPostTree;
        cps->numTrees = numNextTrees;

        queries_added++;
        if( cps->set->debug>2 || debug>2 ) fprintf(stderr, "added query: %d\n", queries_added);
    }

    if( cps->set->debug>2 || debug>2 ) {
        fprintf(stderr, "\n\ntotal # of post trees: %d\n", numCurrentTrees);
        for( i = 0; i < numCurrentTrees; i++ ) {
            tmp[0] = '\0';
            toString(tmp, cps->PostTree[i].root, false);
            fprintf(stderr, "tree #%d: %s\n", i + 1, tmp);
        }
        fprintf(stderr, "\n");
    }

}' CPSetPostTreeMultipleQueries

void CPSamplerDelete(cpsampler **cps) {
    int i;
    if(!*cps) return;
    if((*cps)->part_list) {
        for(i=0; i<(*cps)->npartitions; i++)
            PartitionDelete((*cps)->part_list[i]);
        free((*cps)->part_list);
    }
    if((*cps)->PostTree) free((*cps)->PostTree);
    'if(memory_counts) free(memory_counts);
    'if(memory_PartialLogLikelihood) free(memory_PartialLogLikelihood);
    'if(memory_PartialLogHyperParameterPrior) free(memory_PartialLogHyperParameterPrior);
    'if(memory_HyperParameter) free(memory_HyperParameter);
    'if(memory_Alpha) free(memory_Alpha);
    'if(memory_Tree) free(memory_Tree);
    'if(memory_BreakPoints) free(memory_BreakPoints);
    'for(i=0;i<SWAP_SIZE;i++) if(memory_counts[i]) free(memory_counts[i]);
    'if(memory_counts) free(memory_counts);
    'if(memory_testcounts) free(memory_testcounts);
    'if(memory_doUpdate) free(memory_doUpdate);
    free(*cps);
}' CPSamplerDelete


static void CPOutputLine(const sampler *smp) {
    char tmp[MAX_TREE_STRING];  ' Potential BUG
    int i;
    cpsampler *cps = (cpsampler *)smp->derived_smp;

    strcpy(tmp, "");

    fprintf(cps->smp->fout, "%-6d%3d%3d", cps->smp->JumpNumber, cps->npartitions, cps->cSameT);
    for(i=0; i<cps->npartitions; i++) {
        const partition *part = cps->part_list[i];
        tmp[0] = '\0';
        toString(tmp, cps->PostTree[part->cTree].root, false);
        fprintf(cps->smp->fout," %s %10.2f %5.4f %5.4f %5.4f %5.4f %5.4f %7.4f",
            tmp,
            part->cPartialLogLikelihood,
            part->cmatrix->v[0],
            part->cmatrix->pi[0],part->cmatrix->pi[1],part->cmatrix->pi[2],part->cmatrix->pi[3],
            part->cHyperParameter);
    }
    if( cps->jumpClasses ) for(i=0; i<cps->npartitions; i++) fprintf(cps->smp->fout," %d", cps->part_list[i]->left);
    fprintf(cps->smp->fout, "\n");
}'method CPOutputLine

'static void CPSamplerRun(sampler *smp) {
Public Sub CPSamplerRun(smp As sampler)

    
    Dim sincePrint As Long, i As Long
    sincePrint = 0
    
    Dim BirthOrDeathOrMove As Double
    BirthOrDeathOrMove = 0#
    
    'cpsampler *cps = (cpsampler *)smp->derived_smp;
    Dim cps As cpsampler
    cps = smp.derived_smp
    
    'while( true ) {
    Do While (True)
        ' These were checked on Oct 30, 2001 */
        'cps->bk = (double)cps->lambda / (double)(cps->npartitions);
        cps.bk = cps.lambda / cps.npartitions
        'if( cps->bk > 1.0 ) cps->bk = 1.0;
        If cps.bk > 1# Then cps.bk = 1#
        'cps->dk = (double)(cps->npartitions-1) / (double)cps->lambda;
        cps.dk = (cps.npartitions - 1) / cps.lambda
        'if( cps->dk > 1.0 ) cps->dk = 1.0;
        If cps.dk > 1# Then cps.dk = 1#
        'cps->bkm1 = (double)cps->lambda / (double)(cps->npartitions-1);
        cps.bkm1 = cps.lambda / (cps.npartitions - 1)
        'if( cps->bkm1 > 1.0 ) cps->bkm1 = 1.0;
        If cps.bkm1 > 1# Then cps.bkm1 = 1#
        'cps->dkp1 = (double)(cps->npartitions) / (double)cps->lambda;
        cps.dkp1 = (cps.npartitions) / cps.lambda
        'if( cps->dkp1 > 1.0 ) cps->dkp1 = 1.0;
        If cps.dkp1 > 1# Then cps.dkp1 = 1#
        
        'cps->bk *= cps->C;
        cps.bk = cps.bk * cps.C
        'cps->dk *= cps->C;
        cps.dk = cps.dk * cps.C
        'cps->bkm1 *= cps->C;
        cps.bkm1 = cps.bkm1 * cps.C
        'cps->dkp1 *= cps->C;
        cps.dkp1 = cps.dkp1 * cps.C

        BirthOrDeathOrMove = 0#
        
        'if( cps->jumpClasses ) BirthOrDeathOrMove = cps->set->rng->nextStandardUniform(cps->set->rng);
        If cps.jumpClasses <> 0 Then BirthOrDeathOrMove = nextStandardUniform(cps.setx.rng)
        'if( cps->set->debug>1 || debug>1 ) fprintf(stderr, "Move (%d): %.4f ", cps->set->rng->useRnList?cps->set->rng->current_rn:cps->smp->JumpNumber, BirthOrDeathOrMove);
        
        'if( BirthOrDeathOrMove < (1.0 - cps->bk - cps->dk) ) {
        If BirthOrDeathOrMove < 1# - cps.bk - cps.dk Then
            'if( cps->set->debug>1 || debug>1 ) fprintf(stderr, " fixed dimension sampler\n");
            
            'ProposeNewTopologies(cps, cps->part_list, cps->npartitions, 1.0, true);
            Call ProposeNewTopologies(cps, cps.part_list, cps.npartitions, 1#, 1)
            'ProposeNewQandMu(cps, cps->part_list, cps->npartitions, 1.0, true);
            Call ProposeNewQandMu(cps, cps.part_list, cps.npartitions, 1#, 1)
            if( (cps->npartitions > 1) && cps->jumpClasses  ) {
                CPProposeNewXi(cps, cps->part_list, cps->npartitions, 1.0, true);
            }
        '} else if( BirthOrDeathOrMove < (1.0 - cps->dk) ) {
         ElseIf BirthOrDeathOrMove < 1# - cps.dk Then
            if( cps->set->debug>1 || debug>1 ) fprintf(stderr, " add segment\n");
            if( cps->setx->debug>1 || debug>1 ) fprintf(stderr, " add segment\n");
            ProposeAddSegment(cps);
        '} else if( cps->npartitions > 1 ) {
        ElseIf cps.npartitions > 1 Then
            'if( cps->set->debug>1 || debug>1 ) fprintf(stderr, " delete segment\n");
            if( cps->setx->debug>1 || debug>1 ) fprintf(stderr, " delete segment\n");
            ProposeDeleteSegment(cps);
        End If
        cps->cLogLikelihood = 0.0;
        'for(i=0; i<cps->npartitions; i++) if( cps->set->compute_likelihood ) cps->cLogLikelihood += cps->part_list[i]->cPartialLogLikelihood;
        for(i=0; i<cps->npartitions; i++) if( cps->setx->compute_likelihood ) cps->cLogLikelihood += cps->part_list[i]->cPartialLogLikelihood;
        cps->smp->JumpNumber++;
        sincePrint++;
        'if( cps->smp->JumpNumber > cps->set->burnin ) CPRecordStatistics(cps);
        if( cps->smp->JumpNumber > cps->setx->burnin ) CPRecordStatistics(cps);
        'if( (sincePrint >= cps->set->subsample) ) {
        if( (sincePrint >= cps->setx->subsample) ) {
            sincePrint = 0;
            'fprintf(stderr, "%d ", cps->smp->JumpNumber);
            'if( cps->smp->JumpNumber > cps->set->burnin ) {
            if( cps->smp->JumpNumber > cps->setx->burnin ) {
                CPPrintTopologies(cps->smp, "CPSamplerRun");
                'SamplerSaveEstimates(cps->smp, cps->set->length);
                SamplerSaveEstimates(cps->smp, cps->setx->length);
            }
        }
    Loop
End Sub 'method CPSamplerRun

static void CPPrintTopologies(sampler *smp, const char *where) {
    int i;
    double ll = 0.0;
    cpsampler *cps = (cpsampler *) smp->derived_smp;
    partition **part = cps->part_list;

    if( cps->set->debug || debug ) {
        fprintf(stderr, "%15s (%8d; %.2f): ", where, smp->JumpNumber, (double)avg_nchpts/sample_nchpts);
        for( i=0; i<cps->npartitions; i++ ) {
            ll += part[i]->cPartialLogLikelihood;
            if( i ) fprintf(stderr, " [%4d]", part[i]->left);
            fprintf(stderr, " %.2f %.2f %2d", part[i]->cHyperParameter, part[i]->cmatrix->v[0], part[i]->cTree);
        }
        fprintf(stderr, " :|: %.2f\n", ll);
    }
    if( smp->JumpNumber >= cps->set->length ) {
        double la = cps->lambda;
        int fact = 1;
        fprintf(stderr, "\nChange point statistics : avg(%.4f == %.4f) std(%.4f == %.4f)", (double)avg_nchpts/sample_nchpts, la,
                (std_nchpts - avg_nchpts*avg_nchpts/sample_nchpts)/(double) (sample_nchpts-1), la);
        fprintf(stderr, "\nChange point Poisson:");
        for( i=0; i<10; i++ ) fprintf(stderr, " %.4f", (double)poisson[i]/sample_poisson);
        fprintf(stderr, "\nChange point Poisson:");
        for( i=0; i<10; i++ ) {
            fprintf(stderr, " %.4f", exp(-la)*pow(la,i)/fact);
            fact *= (i+1);
        }
        fprintf(stderr, "\n   Mu: %.4f +/- %.4f", avg_mu/sample_param, (std_mu - avg_mu*avg_mu/sample_param)/(double) (sample_param-1));
        fprintf(stderr, "\nAlpha: %.4f +/- %.4f", avg_alpha/sample_param, (std_alpha - avg_alpha*avg_alpha/sample_param)/(double) (sample_param-1));
        fprintf(stderr, "\n");
    }
}' CPPrintTopologies

Public Sub ProposeNewTopologies(cps As cpsampler, plist() As partition, np As Long, like_factor As Double, record_try As Byte)


    Dim whichTreeJump As Long, t As Long, pTree As Long, i As Long, pUpdateSameT As Long
    'double logRatio = 0.0, pPartialLogLikelihood;
    Dim logRatio As Double, pPartialLogLikelihood As Double, mixing As Double
    Dim part As partition
    logRatio = 0
    'boolean goodTree = false;
    Dim goodTree As Byte
    goodTree = 0
    i = 0
    For t = 0 To np - 1
        'partition *part = plist[t];
        part = plist(t)
        'if( !part->doUpdate ) continue;
        If (part.doUpdate = 0) Then
            ' Propose a new topology for this region
            goodTree = False
            pTree = 0
            whichTreeJump = 0
            Do While (goodTree = 0)
                'double mixing = cps->set->rng->nextStandardUniform(cps->set->rng);
                mixing = nextStandardUniform(cps.setx.rng)
                If mixing < cps.mix1 Then   ' Why do this?
                    'pTree = (int) (cps->set->rng->nextStandardUniform(cps->set->rng)*cps->numTrees);
                    pTree = nextStandardUniform(cps.setx.rng) * cps.numTrees
                    whichTreeJump = 0
                End If
                'if( pTree != part->cTree ) goodTree = true;
                If pTree <> part.cTree Then goodTree = 1
            Loop
            'if( true ) {
            If True Then
    
                'pPartialLogLikelihood = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[pTree], cps->smp, part->cmatrix, part->counts, part->cHyperParameter, false) : 0.0;
                If cps.setx.compute_likelihood = 1 Then
                    pPartialLogLikelihood = TreeLogLikelihood(cps.PostTree(pTree), cps.smp, part.cmatrix, part.counts, part.cHyperParameter, 0)
                Else
                    pPartialLogLikelihood = 0#
                End If
                pUpdateSameT = updateSameT(cps, pTree, t);
    
                if(record_try) cps->smp->tries[whichTreeJump]++;
    
                logRatio = like_factor*pPartialLogLikelihood - like_factor*part->cPartialLogLikelihood  ' likelihood ratio
                    + pUpdateSameT*cps->logWSameT - pUpdateSameT*cps->logWNotSameT;         ' q(\tau^*) / q(\tau)
    
                if( cps->set->debug>1 || debug>1 ) fprintf(stderr, "ProposeNewTopologies (%d): propose tree %d for segment %d -> %.4f", cps->set->rng->useRnList?cps->set->rng->current_rn:cps->smp->JumpNumber, pTree, t, logRatio);
    
                if( LogMHAccept( logRatio, cps->set->rng->nextStandardUniform(cps->set->rng) ) ) {
    
                    if( cps->set->debug>1 || debug>1 ) fprintf(stderr, " ACCEPTED");
    
                    if(record_try) cps->smp->acceptancerate[whichTreeJump]++;
    
                    part->cTree = pTree;
                    part->cPartialLogLikelihood = pPartialLogLikelihood;
                    cps->cSameT += pUpdateSameT;
                }
            Else  ' Update a whole segment (not used)
                'fprintf(stderr, "SHOULD NEVER GET HERE!!!!");
                'exit(EXIT_FAILURE);
            End If
            if( cps->set->debug>1 || debug>1 ) fprintf(stderr, "\n");
        End If
    Next t
End Sub 'method ProposeNewTopologies

static void ProposeNewQandMu(cpsampler *cps, partition **plist, const int np, const double like_factor, const boolean record_try) {
    double pPartialLogLikelihood, pLogHyperParameterPrior, pHyperParameter, logRatio;
    int i;
    qmatrix *pmatrix = NULL;

    for(i=0; i<np; i++) {
        partition *part = plist[i];
        if( !part->doUpdate ) continue;
        part->cmatrix->Matrix_Proposer(&pmatrix, part->cmatrix, cps->set);          ' First time allocates memory
        pHyperParameter = fabs(part->cHyperParameter + cps->set->rng->nextNormal(cps->set->rng, cps->set->sdMu));   ' Propose new average branchlength
        pLogHyperParameterPrior = - pHyperParameter;
        pPartialLogLikelihood = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[part->cTree], cps->smp, pmatrix, part->counts, pHyperParameter, false) : 0.0;
        if(record_try) cps->smp->tries[1]++;
        logRatio = like_factor*pPartialLogLikelihood - like_factor*part->cPartialLogLikelihood      ' likelihood ratio
             + pLogHyperParameterPrior - part->cPartialLogHyperParameterPrior;          ' prior ratio q(\mu*)/q(\mu); q(\alpha*)/q(\alpha) = 1
        if( cps->set->debug>1 || debug>1 ) fprintf(stderr, "ProposeNewQandMu (%d): propose %.4f and %.4f for region %d -> %.4f", cps->set->rng->useRnList?cps->set->rng->current_rn:cps->smp->JumpNumber, pHyperParameter, pmatrix->v[0], i, logRatio);
        if( LogMHAccept( logRatio, cps->set->rng->nextStandardUniform(cps->set->rng) ) ) {
            if( cps->set->debug>1 || debug>1 ) fprintf(stderr, " ACCEPTED");
            part->cPartialLogLikelihood = pPartialLogLikelihood;
            part->cHyperParameter = pHyperParameter;
            part->cPartialLogHyperParameterPrior = pLogHyperParameterPrior;
            part->cmatrix->Matrix_Copy(part->cmatrix, pmatrix);
            if(record_try) cps->smp->acceptancerate[1]++;
        }
        if( cps->set->debug>1 || debug>1 ) fprintf(stderr, "\n");
    }
    if(pmatrix) QMatrixDelete(pmatrix);
    pmatrix = NULL;
}'method ProposeNewQandMu

static void CPProposeNewXi(cpsampler *cps, partition **plist, const int np, const double like_factor, boolean record_try) {
    int rnd, upperBound, lowerBound, proposed, current;
    double leftLikelihood, rightLikelihood, logRatio;
    partition *middle_part;

    for(rnd = 1; rnd < np; rnd++) {
        partition *prev_part = plist[rnd-1];
        partition *part = plist[rnd];
        if( !part->doUpdate ) continue;
        lowerBound = prev_part->left+1;
        upperBound = part->right;
        current = part->left;
        proposed = current;
        if( (upperBound - lowerBound) > 1 ) {
            proposed = ProposeNewChangePointPosition(cps->smp, current, lowerBound, upperBound);
            ' Now proposed is symmetric and reflected and does not equal the original value
            'if( cps->set->debug || debug ) fprintf(stderr, "CPProposeNewXi: %d -> %d on region %d\n", (*_cBreakPoints)[rnd], proposed, rnd);
            if( proposed > current ) {
                PartitionMake(&middle_part, cps->lenunique, current, proposed-1, false, true);
                PartitionCopySegmentCounts(middle_part, cps->sqd, current, proposed);
                ' Calculate likelihood of sites under both models
                leftLikelihood  = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[prev_part->cTree], cps->smp, prev_part->cmatrix, middle_part->counts, prev_part->cHyperParameter, false) : 0.0;
                rightLikelihood = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[part->cTree], cps->smp, part->cmatrix, middle_part->counts, part->cHyperParameter, false) : 0.0;
                logRatio = like_factor*(leftLikelihood - rightLikelihood);
                if( cps->set->debug>1 || debug>1 ) fprintf(stderr, "CPProposeNewXi (%d): move %d to %d (%d,%d) -> %.4f", cps->set->rng->useRnList?cps->set->rng->current_rn:cps->smp->JumpNumber, current, proposed, lowerBound, upperBound, logRatio);
                if( record_try ) cps->smp->tries[2]++;
                if( LogMHAccept( logRatio, cps->set->rng->nextStandardUniform(cps->set->rng) ) ) {
                    if( cps->set->debug>1 || debug>1 ) fprintf(stderr, " ACCEPTED1");
                    PartitionSubtractPartition(part, middle_part);
                    PartitionAddPartition(prev_part, middle_part);
                    prev_part->cPartialLogLikelihood += leftLikelihood;
                    part->cPartialLogLikelihood -= rightLikelihood;
                    part->left = proposed;
                    prev_part->right = proposed - 1;
                    if( record_try ) cps->smp->acceptancerate[2]++;
                }
                if( cps->set->debug>1 || debug>1 ) fprintf(stderr, "\n");
            } else { ' proposed < current
                PartitionMake(&middle_part, cps->lenunique, proposed, current-1, false, true);
                PartitionCopySegmentCounts(middle_part, cps->sqd, proposed, current);
                leftLikelihood  = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[prev_part->cTree], cps->smp, prev_part->cmatrix, middle_part->counts, prev_part->cHyperParameter, false) : 0.0;
                rightLikelihood = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[part->cTree], cps->smp, part->cmatrix, middle_part->counts, part->cHyperParameter, false) : 0.0;
                logRatio = like_factor*(rightLikelihood - leftLikelihood);
                'if( cps->set->debug || debug ) fprintf(stderr, "CPProposeNewXi: likelihood %f - %f (%d, %f, %f)", leftLikelihood, rightLikelihood, _cTree[rnd], _cHyperParameter[rnd], _cAlpha[rnd]);
                if( cps->set->debug>1 || debug>1 ) fprintf(stderr, "CPProposeNewXi (%d): move %d to %d (%d,%d) -> %.4f", cps->set->rng->useRnList?cps->set->rng->current_rn:cps->smp->JumpNumber, current, proposed, lowerBound, upperBound, logRatio);
                if( record_try) cps->smp->tries[2]++;
                if( LogMHAccept( logRatio, cps->set->rng->nextStandardUniform(cps->set->rng) ) ) {
                    if( cps->set->debug>1 || debug>1 ) fprintf(stderr, " ACCEPTED2");
                    PartitionSubtractPartition(prev_part, middle_part);
                    PartitionAddPartition(part, middle_part);
                    prev_part->cPartialLogLikelihood -= leftLikelihood;
                    part->cPartialLogLikelihood += rightLikelihood;
                    part->left = proposed;
                    prev_part->right = proposed - 1;
                    if( record_try ) cps->smp->acceptancerate[2]++;
                }
                if( cps->set->debug>1 || debug>1 ) fprintf(stderr, "\n");
            }
            PartitionDelete(middle_part);
        }
    }
}'method CPProposeNewXi

' TODO: need to write Al-Awadhi stuff to match new setup
static void ProposeAddSegment(cpsampler *cps) {
    ' Add a new segment.
    ' Choose location U[1,lenSeq) and add new EH to right portion of intersected segment
    int i, k, proposed_cp, rnd, local_left_region, local_right_region, start_region, end_region, right_cp, left_cp, end_cp;
    int num_affected_regions, added_region, sameTree;
    double l1, l2, weight1, weight2, zAlpha, zMu, logRatio = 0.0, log_tau_proposal_prob;
    int p1UpdateSameT, p1UpdateNotSameT;
    partition **p1 = NULL, *land_part = NULL, **p2 = NULL;
    boolean abut_left, abut_right, priorDraw, newOnRight, local_debug = false;
    const char *function_name = "ProposeAddSegment";

    ' Propose the new changepoint location wp 1 / (cps->lenSeq - cps->npartitions)
    proposed_cp = 0;       ' New changepoint location
    while( proposed_cp == 0 ) {
        proposed_cp = (int) ( cps->set->rng->nextStandardUniform(cps->set->rng)*cps->lenSeq );
        ' Make sure it isn't overlapping with an existing change point
        for(i=0; i<cps->npartitions; i++) if( proposed_cp == cps->part_list[i]->left ) proposed_cp = 0;
    }
    ' Determine the region containing the new changepoint location
    rnd = 0;            ' Current segment to contain new changepoint location
    while( (rnd < (cps->npartitions-1)) && (proposed_cp > cps->part_list[rnd+1]->left) ) rnd++;
    ' Is the region constrained by alignment ends?
    abut_left = false;
    abut_right = false;
    if( rnd == 0 ) abut_left = true;
    if( rnd == (cps->npartitions-1) ) abut_right = true;
    land_part = cps->part_list[rnd];

    ' Determine number of regions that will be altered by proposal
    ' Must not exceed 2*SWAP_SPACE
    num_affected_regions = 2 + (cps->set->alawadhi ? ((abut_left?0:1) + (abut_right?0:1)) : 0);

    ' Set aside some space to store the proposed state
    /*p1BreakPoints = memory_BreakPoints;
    p1Tree = memory_Tree;
    p1Alpha = memory_Alpha;
    p1HyperParameter = memory_HyperParameter;
    p1PartialLogHyperParameterPrior = memory_PartialLogHyperParameterPrior;
    p1PartialLogLikelihood = memory_PartialLogLikelihood;
    p1testcounts = memory_counts;*/

    if(cps->set->alawadhi) {
        local_left_region = abut_left ? 0 : 1;                  ' region left of new change point
        local_right_region = local_left_region + 1;             ' region right of new change point
        start_region = rnd - (abut_left ? 0 : 1);               ' first region in locally affected region
        end_region = rnd + (abut_right ? 0 : 1);                ' last region in locally affected region
        right_cp = abut_right ? cps->lenSeq : cps->part_list[end_region]->left; ' change point right of proposed one
        end_cp = cps->part_list[end_region]->right+1;               ' right boundary of local region
    } else {
        local_left_region = 0;
        local_right_region = 1;
        start_region = rnd;
        end_region = rnd;
        right_cp = land_part->right + 1;
        end_cp = right_cp;
    }
    left_cp = cps->part_list[start_region]->left;

    ' Set local parameter values equal to the current values, except in two new regions (these params must be generated)
    'p1BreakPoints[0] = cBreakPoints[start_region];
    p1 = (partition **) malloc(sizeof(partition *)*num_affected_regions);
    PartitionMake(&p1[local_right_region], cps->lenunique, proposed_cp, right_cp-1, true, true);
    PartitionCopySegmentCounts(p1[local_right_region], cps->sqd, proposed_cp, right_cp);
    cps->part_list[start_region]->cmatrix->Matrix_Make_Copy(&(p1[local_right_region]->cmatrix), cps->part_list[start_region]->cmatrix);
    added_region = 0;
    for(i=0; i<num_affected_regions; i++) {
        partition *part = NULL;
        if(i==local_right_region) { ' skip one of added region (arbitrarily right)
            added_region = 1;
            continue;
        }
        part = cps->part_list[start_region+i-added_region];
        PartitionMakeCopy(&p1[i], part);
    }
    PartitionSubtractPartition(p1[local_left_region], p1[local_right_region]);
    p1[local_left_region]->right = proposed_cp - 1;

    ' Calculate the weights associated with each new region
    l1 = (double) (proposed_cp - left_cp);
    l2 = (double) (right_cp - proposed_cp);
    weight1 = l1 / (l1 + l2);
    weight2 = l2 / (l1 + l2);

    /*
    logEdgeEffect = 0;
    if( abut_left || abut_right ) logEdgeEffect = cps->logTwo;
    */

    ' Propose new values for the evolutionary parameters in the two new regions
    zAlpha = cps->set->rng->nextStandardNormal(cps->set->rng);
    zMu    = cps->set->rng->nextStandardNormal(cps->set->rng);
    p1[local_left_region]->cHyperParameter = DrawnMu1(land_part->cHyperParameter, weight2, cps->sigmaM, zMu);
    p1[local_left_region]->cPartialLogHyperParameterPrior = - p1[local_left_region]->cHyperParameter;
    p1[local_right_region]->cHyperParameter = DrawnMu2(land_part->cHyperParameter, weight1, cps->sigmaM, zMu);
    p1[local_right_region]->cPartialLogHyperParameterPrior = - p1[local_right_region]->cHyperParameter;
    Set_Alpha(p1[local_left_region]->cmatrix, DrawnAlpha1(land_part->cmatrix->v[0], weight2, cps->sigmaA, zAlpha));
    Set_Alpha(p1[local_right_region]->cmatrix, DrawnAlpha2(land_part->cmatrix->v[0], weight1, cps->sigmaA, zAlpha));

    ' Propose new topologies for the two new regions
    priorDraw = true;
    p1[local_left_region]->cTree = land_part->cTree;
    sameTree = 0;
    if( cps->set->rng->nextStandardUniform(cps->set->rng) < cps->wSameT ) {
        priorDraw = false;
        sameTree = 1;
        log_tau_proposal_prob = cps->logWSameT;
    } else {
        while( p1[local_left_region]->cTree == land_part->cTree ) p1[local_left_region]->cTree = cps->set->rng->nextStandardUniform(cps->set->rng)*cps->numTrees;
        log_tau_proposal_prob = cps->logWNotSameT - cps->logTwo;    ' logTwo is for newOnRight probability
    }
    newOnRight = false;
    if( cps->set->rng->nextStandardUniform(cps->set->rng) < 0.5 ) newOnRight = true;
    if(newOnRight) {
        p1[local_right_region]->cTree = p1[local_left_region]->cTree;
        p1[local_left_region]->cTree = land_part->cTree;
    } else {
        p1[local_right_region]->cTree = land_part->cTree;
    }

    ' Calculate the relative change to count of same neighbor trees (for prior on tau)
    ' |---left neighbor---|---left - right---|---right neighbor---|
    p1UpdateSameT = 0;
    p1UpdateNotSameT = 0;
    if( p1[local_left_region]->cTree == p1[local_right_region]->cTree ) p1UpdateSameT++;
    if( p1[local_left_region]->cTree != p1[local_right_region]->cTree ) p1UpdateNotSameT++;
    if( newOnRight ) {
        if( !abut_right ) {
            if( land_part->cTree == cps->part_list[rnd+1]->cTree )          p1UpdateSameT--;    ' middle was == right neighbor
            if( land_part->cTree != cps->part_list[rnd+1]->cTree )          p1UpdateNotSameT--; ' middle was == right neighbor
            if( p1[local_right_region]->cTree == cps->part_list[rnd+1]->cTree ) p1UpdateSameT++;    ' right == right neighbor
            if( p1[local_right_region]->cTree != cps->part_list[rnd+1]->cTree ) p1UpdateNotSameT++; ' right == right neighbor
        }
    } else {
        if( !abut_left ) {
            if( land_part->cTree == cps->part_list[rnd-1]->cTree )          p1UpdateSameT--;    ' middle was == left neighbor
            if( land_part->cTree != cps->part_list[rnd-1]->cTree )          p1UpdateNotSameT--; ' middle was == left neighbor
            if( p1[local_left_region]->cTree == cps->part_list[rnd-1]->cTree )  p1UpdateSameT++;    ' left == left neighbor
            if( p1[local_left_region]->cTree != cps->part_list[rnd-1]->cTree )  p1UpdateNotSameT++; ' left == left neighbor
        }
    }
if( p1UpdateSameT + p1UpdateNotSameT != 1 ) {
fprintf(stderr, "Add %d: %d %d %d %d -> %d, %d = %d\n", proposed_cp, abut_left?-1:cps->part_list[rnd-1]->cTree, p1[local_left_region]->cTree, p1[local_right_region]->cTree, abut_right?-1:cps->part_list[rnd+1]->cTree, p1UpdateSameT, p1UpdateNotSameT, p1UpdateSameT + p1UpdateNotSameT);
exit(0);
}
    'if( cps->set->debug || debug ) printf("%s proposal 1 location %d-%d-%d-%d-%d\n", function_name, p1BreakPoints[0], left_cp, proposed_cp, right_cp, end_cp);
    for(i=0;i<num_affected_regions;i++)
        p1[i]->cPartialLogLikelihood = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[p1[i]->cTree], cps->smp, p1[i]->cmatrix, p1[i]->counts, p1[i]->cHyperParameter, false) : 0.0;
    'if( cps->set->debug || debug ) fprintf(stderr, "\tmu: %f -> (%f %f)\n", cHyperParameter[rnd], p1HyperParameter[local_left_region], p1HyperParameter[local_right_region]);
    'if( cps->set->debug || debug ) fprintf(stderr, "\talpha: %f -> (%f %f)\n", cAlpha[rnd], p1Alpha[local_left_region], p1Alpha[local_right_region]);
    'if( cps->set->debug || debug ) fprintf(stderr, "\ttrees: %d -> (%d %d) (%d)\n", cTree[rnd], p1Tree[local_left_region], p1Tree[local_right_region], p1UpdateSameT);
    ' Finished collecting information on proposal 1; Start collecting info on proposal 2

    if(cps->set->alawadhi) {
        int local_cSameT;

        if( cps->set->debug>2 || debug>2 || local_debug ) fprintf(stderr, "...using Al-Awadhi add...\n");
        'p2BreakPoints = &memory_BreakPoints[num_affected_regions];
        'p2Tree = &memory_Tree[num_affected_regions];
        'p2Alpha = &memory_Alpha[num_affected_regions];
        'p2HyperParameter = &memory_HyperParameter[num_affected_regions];
        'p2PartialLogHyperParameterPrior = &memory_PartialLogHyperParameterPrior[num_affected_regions];
        'p2PartialLogLikelihood = &memory_PartialLogLikelihood[num_affected_regions];
        'p2testcounts = memory_counts + num_affected_regions;

        p2 = (partition **) malloc(sizeof(partition *)*num_affected_regions);

        ' Except for surrounding changepoints; allow the neighboring (left and right) and proposed changepoints to move too
        ' |---left neighbor---|---left---|---right---|---right neighbor--|
        '                     ^          ^           ^
        for(i=0; i<num_affected_regions; i++) {
             PartitionMakeCopy(&p2[i], p1[i]);
             p2[i]->doUpdate = false;
        }
        ' Only update two regions around new changepoint
        p2[local_left_region]->doUpdate = true;
        p2[local_right_region]->doUpdate = true;
        p2[0]->doXiUpdate = false;

        local_cSameT = determineSameT(num_affected_regions, p2);

        ' Go through alawadhi_k cycles of updates within dimension
        for(k=0;k<cps->set->alawadhi_k;k++) {
            ProposeNewTopologies(cps, p2, num_affected_regions, cps->set->alawadhi_factor, false);
            ProposeNewQandMu(cps, p2, num_affected_regions, cps->set->alawadhi_factor, false);
            CPProposeNewXi(cps, p2, num_affected_regions, cps->set->alawadhi_factor, false);
        }
        if( cps->set->debug>2 || debug>2 || local_debug ) {
            for(i=0;i<num_affected_regions;i++) fprintf(stderr, "%s proposal 2 Loglikelihood[%d]: %f\n", function_name, i, p2[i]->cPartialLogLikelihood);
            fprintf(stderr, "%s: proposal 2 location %d-%d-%d\n", function_name, p2[local_left_region]->left, p2[local_right_region]->left, end_cp);
            fprintf(stderr, "\tmu1: %f -> (%f %f)\n", land_part->cHyperParameter, p2[local_left_region]->cHyperParameter, p2[local_right_region]->cHyperParameter);
            fprintf(stderr, "\talpha: %f -> (%f %f)\n", land_part->cmatrix->v[0], p2[local_left_region]->cmatrix->v[0], p2[local_right_region]->cmatrix->v[0]);
            fprintf(stderr, "\ttrees: %d -> (%d %d)\n", land_part->cTree, p2[local_left_region]->cTree, p2[local_right_region]->cTree);
        }
        ' Calculate likelihood ratio for move theta^{(k^*)} to theta^{(k')} (2/12/04 notes)
        for(i=0;i<num_affected_regions;i++) {
            logRatio += cps->set->alawadhi_factor*( p1[i]->cPartialLogLikelihood - p2[i]->cPartialLogLikelihood );
            if( cps->set->debug>2 || debug>2 || local_debug )
                fprintf(stderr, "%s Al-Awadhi LR[%d]: %f - %f\n", function_name, i, p1[i]->cPartialLogLikelihood, p2[i]->cPartialLogLikelihood);
        }
        for(i=0;i<num_affected_regions;i++) {
            logRatio += p2[i]->cPartialLogLikelihood;
            if( cps->set->debug>2 || debug>2 || local_debug ) fprintf(stderr, "%s LR-top[%d]: %f\n", function_name, i, p2[i]->cPartialLogLikelihood);
        }
        for(i=0;i<num_affected_regions-1;i++) {
            logRatio -= cps->part_list[start_region+i]->cPartialLogLikelihood;
            if( cps->set->debug>2 || debug>2 || local_debug )
                fprintf(stderr, "%s LR-bottom[%d]: %f\n", function_name, start_region+i, cps->part_list[start_region+i]->cPartialLogLikelihood);
        }
    } else {
        for(i=0;i<num_affected_regions;i++) {
            logRatio += p1[i]->cPartialLogLikelihood;
            'if( cps->set->debug || debug ) fprintf(stderr, "%s proposal 1 LogLikelihood[%d]: %f\n", function_name, i, p1[i]->cPartialLogLikelihood);
        }
        logRatio -= land_part->cPartialLogLikelihood;
        'if( cps->set->debug || debug ) fprintf(stderr, "%s proposal 1 cLogLikelihood[%d]: %f\n", function_name, rnd, land_part->cPartialLogLikelihood);
    }
    ' Add prior ratio on mu (alpha has uniform prior so its prior cancels)
    for(i=0;i<num_affected_regions;i++) logRatio += - p1[i]->cHyperParameter;
    for(i=0;i<num_affected_regions-1;i++) logRatio -= cps->part_list[start_region+i]->cPartialLogHyperParameterPrior;

    ' q(J) = exp(-lambda)(lambda)^J/J! ==> q(J+1)/q(J) = lambda/(J+1)
    ' q(\xi|J) = \frac{J!}{(L-1)(L-2)...(L-J)} ==> q(\xi|J+1)/q(\xi|J) = (J+1)/(L-J-1)
    ' q(\tau|J=K+I) = (1/ntrees)*((1-w)/(ntrees-1))^Kw^I
    ' p(\tau*|\tau) = w or (1-w)/(2*(ntrees-1))
    /* ORIGINAL */
    logRatio += cps->logLambda          ' PRIOR: q(J+1,\xi)/q(J,\xi)
'      - log((double) cps->npartitions)    ' CORRECTION KSD
'      - log((double) (cps->npartitions + 1))          ' Prior on I (Model Mi) CORRECTION KSD
'      + logEdgeEffect             ' CORRECTION KSD
    '  - cps->logTwo
        - log( cps->npartitions ) - log_tau_proposal_prob   ' KSD added
'      + ((double) (p1UpdateSameT-sameTree))*cps->logWSameT + ((double) (sameTree - p1UpdateSameT))*cps->logWNotSameT  ' Prior on tau
        + p1UpdateSameT*cps->logWSameT + p1UpdateNotSameT*cps->logWNotSameT ' Prior on tau KSD CORRECTION
        + log( cps->dkp1 ) - log( cps->bk )                                     ' Dimensional move proposal
        - logStandardNormalDensity(zAlpha) - logStandardNormalDensity(zMu)                      ' Parameter proposal
        + CPLogJacobian(land_part->cmatrix->v[0], p1[local_left_region]->cmatrix->v[0], p1[local_right_region]->cmatrix->v[0], cps->sigmaA,
            land_part->cHyperParameter, p1[local_left_region]->cHyperParameter, p1[local_right_region]->cHyperParameter, cps->sigmaM);
    /**/
    /* KSD
    logRatio += cps->logLambda - log(cps->lenSeq - cps->npartitions)                            ' Prior on \xi,J
        + p1UpdateSameT*cps->logWSameT + p1UpdateNotSameT*cps->logWNotSameT                     ' Prior on tau
        + Clog( cps->dkp1 ) - cps->logTwo - Clog( cps->npartitions )                            ' Proposal: death
        - Clog( cps->bk ) - log_tau_proposal_prob - ClogStandardNormalDensity(zAlpha) - ClogStandardNormalDensity(zMu)  ' Proposal: birth -> move, \tau, \alpha, \mu
        + Clog( cps->lenSeq - cps->npartitions )                                    ' Proposal: birth -> \xi
        + CPLogJacobian(land_part->cmatrix->v[0], p1[local_left_region]->cmatrix->v[0], p1[local_right_region]->cmatrix->v[0], cps->sigmaA,
            land_part->cHyperParameter, p1[local_left_region]->cHyperParameter, p1[local_right_region]->cHyperParameter, cps->sigmaM);
    */

    if( cps->set->debug>1 || debug>1 || local_debug )
        fprintf(stderr, "tagProposeAddSegment (%d): propose %d with trees (%d,%d), alphas (%.4f,%.4f) and mus (%.4f,%.4f) -> %.4f %f %f\n", cps->set->rng->useRnList?cps->set->rng->current_rn:cps->smp->JumpNumber,proposed_cp, p1[local_left_region]->cTree, p1[local_right_region]->cTree, p1[local_left_region]->cmatrix->v[0], p1[local_right_region]->cmatrix->v[0], p1[local_left_region]->cHyperParameter, p1[local_right_region]->cHyperParameter, exp(logRatio), log(cps->npartitions+1), - log( cps->npartitions+1 )- log_tau_proposal_prob - logStandardNormalDensity(zAlpha) - logStandardNormalDensity(zMu)+ log( cps->lenSeq - cps->npartitions ) );

    'if( cps->set->debug || debug ) fprintf(stderr, "%s LogRatio: %f %f %f %f %f %d %f %f %f %f %f %f %f %f %f %f\n", function_name, logRatio, -p1HyperParameter[0], -p1HyperParameter[1], cPartialLogHyperParameterPrior[rnd], logEdgeEffect, (p1UpdateSameT-sameTree), logLambda, log(numDatasets), log(numDatasets+1), logWSameT, logWNotSameT, log(dkp1), log(bk), logStandardNormalDensity(zAlpha), logStandardNormalDensity(zMu), zMu);
    if( cps->set->debug>1 || debug>1 || local_debug ) fprintf(stderr, "%s (%d): %f", function_name, cps->smp->JumpNumber, logRatio);
    if( priorDraw ) cps->smp->tries[3]++;
    else cps->smp->tries[4]++;

    if( LogMHAccept( logRatio, cps->set->rng->nextStandardUniform(cps->set->rng) ) ) {
        partition *new_part;

        ' Treat right region of landing partition as new partition (and
        ' modify the left region for new right change point)
        PartitionMakeCopy(&new_part, cps->set->alawadhi ? p2[local_right_region] : p1[local_right_region]);
        CPSAddPartition(cps, new_part, rnd+1);

        ' Copy other updated regions
        added_region = 0;
        for(i=0;i<num_affected_regions;i++) {
            if(i == start_region + local_right_region) {
                added_region = 1;
                continue;
            }
            if(cps->set->alawadhi)
                PartitionCopy(cps->part_list[start_region+i+added_region], p2[i+added_region]);
            Else
                PartitionCopy(cps->part_list[start_region+i+added_region], p1[i+added_region]);
        }

        cps->cSameT = determineSameT(cps->npartitions, cps->part_list);
        if( priorDraw ) cps->smp->acceptancerate[3]++;
        else cps->smp->acceptancerate[4]++;
        total_adds++;
        if( cps->set->debug>1 || debug>1 || local_debug ) fprintf(stderr, "Added region.\n");
    } else {
        if( cps->set->debug>1 || debug>1 || local_debug ) fprintf(stderr, "\n");
    }

    ' Clear temporary partitions
    for(i=0; i<num_affected_regions; i++) {
         PartitionDelete(p1[i]);
         if(cps->set->alawadhi) PartitionDelete(p2[i]);
    }
    free(p1);
    if( cps->set->alawadhi ) free(p2);
}'method ProposeAddSegment

static void ProposeDeleteSegment(cpsampler *cps) {
'  const char *function_name = "ProposeDeleteSegment";
    ' Remove a segment
    ' Choose segment U[1,numDatasets-1) and collapse to left EH
    int i, rnd, num_affected_regions, left_region, right_region, start_region, tree_region, start, end, midpt;
    int local_left_region, local_right_region, local_tree_region, sameTree;
    double weight1, weight2, zAlpha, zMu, logRatio, log_tau_proposal_prob;
    partition **p1 = NULL, **p2 = NULL, *new_part = NULL;
    int pUpdateSameT = 0, pUpdateNotSameT = 0;
    boolean abut_left, abut_right, keepLeft, local_debug = false;

    ' Select regions to merge wp 1/cps->npartitions/2 (divide by 2 iff select internal region)
    'rnd = (int) ( cps->set->rng->nextStandardUniform(cps->set->rng)*cps->npartitions );
    rnd = (int) ( cps->set->rng->nextStandardUniform(cps->set->rng)*(cps->npartitions - 1));    ' Select randomly one change point to delete
    if( cps->set->rng->nextStandardUniform(cps->set->rng) < 0.5 ) keepLeft = true;
    else keepLeft = false;

    abut_left = false;
    abut_right = false;
    /*
    logEdgeEffect = 0;
    if( rnd == 0 ) {
        keepLeft = false;
        abut_left = true;
        logEdgeEffect = cps->logTwo;
    } else if( rnd == (cps->npartitions - 1) ) {
        keepLeft = true;
        abut_right = true;
        logEdgeEffect = cps->logTwo;
    }
    */
    if(rnd==0) abut_left = true;
    if(rnd==cps->npartitions-2) abut_right = true;
    num_affected_regions = 2 + (cps->set->alawadhi ? ((abut_left?0:1) + (abut_right?0:1)) : 0);
    start_region = rnd - ((cps->set->alawadhi && !abut_left)?1:0);
    left_region = rnd;
    right_region = rnd + 1;
    tree_region = keepLeft ? left_region : right_region;
    local_left_region = (abut_left?0:1);
    local_right_region = (abut_left?1:2);
    local_tree_region = keepLeft ? local_left_region : local_right_region;

    if(cps->set->alawadhi) {
        int local_cSameT;
        if( cps->set->debug>2 || debug>2 || local_debug ) fprintf(stderr, "...using Al-Awadhi delete...\n");

        p1 = (partition **) malloc(sizeof(partition *)*num_affected_regions);

        ' First we need to take several steps in this higher dimension
        ' Record current state in p1

        for(i=0; i<num_affected_regions; i++)
            PartitionMakeCopy(&p1[i], cps->part_list[start_region+i]);
        end = (start_region+num_affected_regions==cps->npartitions) ? cps->lenSeq : cps->part_list[start_region+num_affected_regions]->left;
        ' We have recorded the original state of the system

        ' Now copy over into p2, which we'll update several steps

        p2 = (partition **) malloc(sizeof(partition *)*num_affected_regions);

        for(i=0; i<num_affected_regions; i++) {
            PartitionMakeCopy(&p2[i], p1[i]);
            p2[i]->doUpdate = false;
        }
        p2[local_left_region]->doUpdate = true;
        p2[local_right_region]->doUpdate = true;
        p2[0]->doXiUpdate = false;

        local_cSameT = determineSameT(num_affected_regions, p2);

        ' Go through alawadhi_k cycles of updates within dimension
        for(i=0;i<cps->set->alawadhi_k;i++) {
            ProposeNewTopologies(cps, p2, num_affected_regions, cps->set->alawadhi_factor, false);
            ProposeNewQandMu(cps, p2, num_affected_regions, cps->set->alawadhi_factor, false);
            CPProposeNewXi(cps, p2, num_affected_regions, cps->set->alawadhi_factor, false);
        }

        if( cps->set->debug>2 || debug>2 || local_debug ) {
            fprintf(stderr, "ProposeDeleteSegment proposal:\n");
            for(i=0;i<num_affected_regions;i++) {
                fprintf(stderr, "\tBreakpoint: %d -> %d\n", p1[i]->left, p2[i]->left);
                fprintf(stderr, "\tTree: %d -> %d\n", p1[i]->cTree, p2[i]->cTree);
                fprintf(stderr, "\tMu: %f -> %f\n", p1[i]->cHyperParameter, p2[i]->cHyperParameter);
                fprintf(stderr, "\tAlpha: %f -> %f\n", p1[i]->cmatrix->v[0], p2[i]->cmatrix->v[0]);
                fprintf(stderr, "\tLL: %f -> %f\n", p1[i]->cPartialLogLikelihood, p2[i]->cPartialLogLikelihood);
            }
        }

        ' We now have the second state of the system

        ' Calculate likelihood ratio for move theta^{(k^*)} to theta^{(k')} (2/12/04 notes)
        logRatio = 0;
        for(i=0;i<num_affected_regions;i++) {
            if( cps->set->debug>2 || debug>2 || local_debug )
                fprintf(stderr, "ProposeDeleteSegment Al-Awadhi LR[%d]: %f %f\n", i, p1[i]->cPartialLogLikelihood, p2[i]->cPartialLogLikelihood);
            logRatio += cps->set->alawadhi_factor*( p1[i]->cPartialLogLikelihood - p2[i]->cPartialLogLikelihood );
        }

        ' Make the dimension-changing proposal

        ' Calculate the change on number of neighboring same trees
        pUpdateSameT = 0;
        '     |---left neighbor---|---left---|---right---|---right neighbor---|
        if( p2[left_region]->cTree == p2[right_region]->cTree ) pUpdateSameT--;         ' left == right
        if(!abut_left && p2[left_region-1]->cTree == p2[left_region]->cTree) pUpdateSameT--;    ' left neighbor == left
        if(!abut_left && p2[left_region-1]->cTree == p2[right_region]->cTree) pUpdateSameT++;   ' left neighbor == right
        if(!abut_right && p2[right_region]->cTree == p2[right_region+1]->cTree) pUpdateSameT--; ' right == right neighbor
        if(!abut_right && p2[left_region]->cTree == p2[right_region+1]->cTree) pUpdateSameT++;  ' left == right neighbor

        sameTree = 0;
        if( cps->part_list[left_region]->cTree == cps->part_list[right_region]->cTree ) {
            sameTree++;
            log_tau_proposal_prob = cps->logWSameT;
        } else {
            log_tau_proposal_prob = cps->logWNotSameT;' - cps->logTwo;
        }

        start = p2[local_left_region]->left;
        midpt = p2[local_right_region]->left;
        end = (right_region+1==cps->npartitions) ? cps->lenSeq : p2[local_right_region+1]->left;

        PartitionMake(&new_part, cps->lenunique, p2[local_left_region]->left, p2[local_right_region]->right, p2[local_left_region]->topchange, p2[local_left_region]->parchange);
        PartitionCopyPartitionSum(new_part, p2[local_left_region], p2[local_right_region]);

        weight1 = (double) (midpt - start) / (end - start);
        weight2 = (double) (end - midpt) / (end - start);

        Set_Alpha(new_part->cmatrix, CondensedAlpha( p2[local_left_region]->cmatrix->v[0], weight1, p2[local_right_region]->cmatrix->v[0], weight2 ));
        new_part->cHyperParameter = CondensedMu( p2[local_left_region]->cHyperParameter, weight1, p2[local_right_region]->cHyperParameter, weight2 );
        new_part->cPartialLogHyperParameterPrior = - new_part->cHyperParameter;
        zAlpha = InverseZAlpha( p2[local_left_region]->cmatrix->v[0], new_part->cmatrix->v[0], weight2, cps->sigmaA );
        zMu = InverseZMu( p2[local_left_region]->cHyperParameter, new_part->cHyperParameter, weight2,  cps->sigmaM );
        'if(cps->set->debug) fprintf(stderr, "ProposeDeleteSegment dimension change proposal\n");
        'if(cps->set->debug) fprintf(stderr, "\tAlpha: %f %f -> %f\n", p2Alpha[local_left_region], p2Alpha[local_right_region], new_part->cmatrix->v[0]);
        'if(cps->set->debug) fprintf(stderr, "\tMu: %f %f -> %f\n", p2HyperParameter[local_left_region], p2HyperParameter[local_right_region], new_part->cHyperParameter);
        'if(cps->set->debug) fprintf(stderr, "\tSuppl: %f %f\n", zAlpha, zMu);
        new_part->cPartialLogLikelihood = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[p2[local_tree_region]->cTree], cps->smp, new_part->cmatrix, new_part->counts, new_part->cHyperParameter, false) : 0.0;
        logRatio += new_part->cPartialLogLikelihood + new_part->cPartialLogHyperParameterPrior;
        'if(cps->set->debug) fprintf(stderr, "ProposeDeleteSegment LR-top[middle]: %f\n", new_part->cPartialLogLikelihood);
        if(!abut_left) logRatio += p2[0]->cPartialLogLikelihood + p2[0]->cPartialLogHyperParameterPrior;
        'if(!abut_left && set->debug) fprintf(stderr, "ProposeDeleteSegment LR-top[left]: %f\n", p2PartialLogLikelihood[0]);
        if(!abut_right) logRatio += p2[local_right_region+1]->cPartialLogLikelihood + p2[local_right_region+1]->cPartialLogHyperParameterPrior;
        'if(!abut_right && set->debug) fprintf(stderr, "ProposeDeleteSegment LR-top[right]: %f\n", p2PartialLogLikelihood[local_right_region+1]);

        for(i=0;i<num_affected_regions;i++) {
            logRatio -= cps->part_list[start_region+i]->cPartialLogLikelihood - cps->part_list[start_region+i]->cPartialLogHyperParameterPrior;
            'if(cps->set->debug) fprintf(stderr, "ProposeDeleteSegment LR-bottom[%d]: %f\n", start_region+i, cPartialLogLikelihood[start_region+i]);
        }
        logRatio -= CPLogJacobian(new_part->cmatrix->v[0], p2[local_left_region]->cmatrix->v[0], p2[local_right_region]->cmatrix->v[0], cps->sigmaA, new_part->cHyperParameter, p2[local_left_region]->cHyperParameter, p2[local_right_region]->cHyperParameter, cps->sigmaM);

    } else {
        ' Make the dimension-changing proposal

        ' Calculate the change on number of neighboring same trees
        pUpdateSameT = 0;
        pUpdateNotSameT = 0;
        '     |---left neighbor---|---left---|---right---|---right neighbor---|
        if( cps->part_list[left_region]->cTree == cps->part_list[right_region]->cTree ) pUpdateSameT--;                 ' left == right
        if( cps->part_list[left_region]->cTree != cps->part_list[right_region]->cTree ) pUpdateNotSameT--;              ' left != right
        if(!abut_left && !keepLeft && cps->part_list[left_region-1]->cTree == cps->part_list[left_region]->cTree) pUpdateSameT--;   ' left neighbor == left
        if(!abut_left && !keepLeft && cps->part_list[left_region-1]->cTree == cps->part_list[right_region]->cTree) pUpdateSameT++;  ' left neighbor == right
        if(!abut_left && !keepLeft && cps->part_list[left_region-1]->cTree != cps->part_list[left_region]->cTree) pUpdateNotSameT--;    ' left neighbor == left
        if(!abut_left && !keepLeft && cps->part_list[left_region-1]->cTree != cps->part_list[right_region]->cTree) pUpdateNotSameT++;   ' left neighbor == right
        if(!abut_right && keepLeft && cps->part_list[right_region]->cTree == cps->part_list[right_region+1]->cTree) pUpdateSameT--; ' right == right neighbor
        if(!abut_right && keepLeft && cps->part_list[left_region]->cTree == cps->part_list[right_region+1]->cTree) pUpdateSameT++;  ' left == right neighbor
        if(!abut_right && keepLeft && cps->part_list[right_region]->cTree != cps->part_list[right_region+1]->cTree) pUpdateNotSameT--;  ' right == right neighbor
        if(!abut_right && keepLeft && cps->part_list[left_region]->cTree != cps->part_list[right_region+1]->cTree) pUpdateNotSameT++;   ' left == right neighbor
if( pUpdateSameT + pUpdateNotSameT != -1 ) {
fprintf(stderr, "Delete %d keeping %5s: %d %d %d %d -> %d, %d = %d\n", cps->part_list[right_region]->left, keepLeft?"left":"right", abut_left?-1:cps->part_list[left_region-1]->cTree, cps->part_list[left_region]->cTree, cps->part_list[right_region]->cTree, abut_right?-1:cps->part_list[right_region+1]->cTree, pUpdateSameT, pUpdateNotSameT, pUpdateSameT + pUpdateNotSameT);
exit(0);
}

        sameTree = 0;
        if( cps->part_list[left_region]->cTree == cps->part_list[right_region]->cTree ) {
            sameTree++;
            log_tau_proposal_prob = cps->logWSameT;
        } else {
            log_tau_proposal_prob = cps->logWNotSameT;' - cps->logTwo;
        }

        PartitionMake(&new_part, cps->lenunique, cps->part_list[left_region]->left, cps->part_list[right_region]->right, cps->part_list[left_region]->topchange, cps->part_list[left_region]->parchange);
        PartitionCopyPartitionSum(new_part, cps->part_list[left_region], cps->part_list[right_region]);
        cps->part_list[right_region]->cmatrix->Matrix_Make_Copy(&(new_part->cmatrix), cps->part_list[right_region]->cmatrix);

        start = new_part->left;
        midpt = cps->part_list[right_region]->left;
        end = new_part->right + 1;

        weight1 = (double) (midpt - start) / (end - start);
        weight2 = (double) (end - midpt) / (end - start);

        Set_Alpha(new_part->cmatrix, CondensedAlpha( cps->part_list[left_region]->cmatrix->v[0], weight1, cps->part_list[right_region]->cmatrix->v[0], weight2 ));
        new_part->cHyperParameter = CondensedMu( cps->part_list[left_region]->cHyperParameter, weight1, cps->part_list[right_region]->cHyperParameter, weight2 );
        new_part->cPartialLogHyperParameterPrior = - new_part->cHyperParameter;
        zAlpha = InverseZAlpha( cps->part_list[left_region]->cmatrix->v[0], new_part->cmatrix->v[0], weight2, cps->sigmaA );
        zMu = InverseZMu( cps->part_list[left_region]->cHyperParameter, new_part->cHyperParameter, weight2,  cps->sigmaM );
        new_part->cTree = cps->part_list[tree_region]->cTree;

        'if(cps->set->debug) fprintf(stderr, "ProposeDeleteSegment dimension change proposal (%d, %d, %d, %d)\n", left_region, right_region, tree_region, numDatasets);
        'if(cps->set->debug) fprintf(stderr, "\tAlpha: %f %f -> %f\n", cAlpha[left_region], cAlpha[right_region], new_part->cmatrix->v[0]);
        'if(cps->set->debug) fprintf(stderr, "\tMu: %f %f -> %f\n", cHyperParameter[left_region], cHyperParameter[right_region], new_part->cHyperParameter);
        'if(cps->set->debug) fprintf(stderr, "\tSuppl: %f %f\n", zAlpha, zMu);
        'if(cps->set->debug) fprintf(stderr, "\tTree: %d %d -> %d\n", cTree[left_region], cTree[right_region], cTree[tree_region]);
        new_part->cPartialLogLikelihood = cps->set->compute_likelihood ? TreeLogLikelihood(&cps->PostTree[new_part->cTree], cps->smp, new_part->cmatrix, new_part->counts, new_part->cHyperParameter, false) : 0.0;
        'if(cps->set->debug) fprintf("ProposeDeleteSegment LR-top[middle]: %f\n", pPartialLogLikelihood);
        logRatio = new_part->cPartialLogLikelihood - cps->part_list[left_region]->cPartialLogLikelihood - cps->part_list[right_region]->cPartialLogLikelihood;
        'if(cps->set->debug) fprintf(stderr, "ProposeDeleteSegment LR-bottom: [%d] %f [%d] %f\n", left_region, cPartialLogLikelihood[left_region], right_region, cPartialLogLikelihood[right_region]);

        logRatio -= CPLogJacobian(new_part->cmatrix->v[0], cps->part_list[left_region]->cmatrix->v[0], cps->part_list[right_region]->cmatrix->v[0], cps->sigmaA, new_part->cHyperParameter, cps->part_list[left_region]->cHyperParameter, cps->part_list[right_region]->cHyperParameter, cps->sigmaM);
        logRatio += new_part->cPartialLogHyperParameterPrior - cps->part_list[left_region]->cPartialLogHyperParameterPrior - cps->part_list[right_region]->cPartialLogHyperParameterPrior;
    }

    /* ORIGINAL */
    logRatio += - log( cps->dk ) + log( cps->bkm1 )                                         ' Birth/death proposal ratio
        - cps->logLambda
        + log( cps->npartitions - 1)
'      + log( cps->npartitions )   ' CORRECTION KSD
'      - logEdgeEffect                 ' Prior ratio on K CORRECTION KSD
'      + cps->logTwo
        + log_tau_proposal_prob ' KSD add
'      + ((double) (pUpdateSameT+sameTree))*cps->logWSameT - ((double) (pUpdateSameT+sameTree))*cps->logWNotSameT      ' prior ratio on tau
        + pUpdateSameT*cps->logWSameT + pUpdateNotSameT*cps->logWNotSameT       ' prior ratio on tau KSD CORRECTION
        + logStandardNormalDensity(zAlpha) + logStandardNormalDensity(zMu);                         ' proposal ratio
    /* */
    /* KSD
    logRatio += - Clog( cps->dk ) + Clog( cps->npartitions - 1 ) + cps->logTwo                          ' Proposal: death -> move, partition, merge direction
        + Clog( cps->bkm1 ) + ClogStandardNormalDensity(zAlpha) + ClogStandardNormalDensity(zMu)                ' Proposal: birth -> move, \mu, \alpha
        + log_tau_proposal_prob - Clog(cps->lenSeq - cps->npartitions + 1)                          ' Proposal: birth -> \tau, \xi
        - Ccps->logLambda + Clog( cps->lenSeq - cps->npartitions + 1)                               ' Prior ratio on K
        + C((double) (pUpdateSameT+sameTree))*cps->logWSameT - ((double) (pUpdateSameT+sameTree))*cps->logWNotSameT     ' Prior ratio on tau
        ;
    */

    if( cps->set->debug>1 || debug>1 || local_debug )
        fprintf(stderr, "tagProposeDeleteSegment (%d): propose to delete %d with alpha %.4f, mu %.4f, and tree %d -> %.4f %f %f", cps->set->rng->useRnList?cps->set->rng->current_rn:cps->smp->JumpNumber, midpt, new_part->cmatrix->v[0], new_part->cHyperParameter, new_part->cTree, exp(logRatio), log( cps->npartitions) ,  + log( cps->npartitions ) + logStandardNormalDensity(zAlpha) + logStandardNormalDensity(zMu)+ log_tau_proposal_prob - log(cps->lenSeq - cps->npartitions + 1));
    if(sameTree) cps->smp->tries[5]++;
    else cps->smp->tries[6]++;
    if( LogMHAccept( logRatio, cps->set->rng->nextStandardUniform(cps->set->rng) ) ) {

        'if(cps->set->debug) fprintf(stderr, "%d %d %d %d\n", start_region, left_region, right_region, num_affected_regions);
        CPSRemovePartition(cps, left_region);
        PartitionCopy(cps->part_list[left_region], new_part);
        if( cps->set->alawadhi && !abut_left ) {
            PartitionCopy(cps->part_list[start_region], p2[0]);
            if( cps->set->debug>2 || debug>2 || local_debug ) fprintf(stderr, "%d: %f\n", start_region, p2[0]->cPartialLogLikelihood);
        }
        'fprintf(stderr, "Setting %d\n", left_region);
        if( cps->set->alawadhi && !abut_right ) {
            PartitionCopy(cps->part_list[right_region], p2[local_right_region+1]);
            'if( cps->set->debug || debug ) fprintf(stderr, "%d: %f\n", right_region, p2[local_right_region+1]->cPartialLogLikelihood);
        }

        cps->cSameT = determineSameT(cps->npartitions, cps->part_list);
        if(sameTree) cps->smp->acceptancerate[5]++;
        else cps->smp->acceptancerate[6]++;
        total_deletes++;
        if( cps->set->debug>2 || debug>2 || local_debug ) fprintf(stderr, "Removed region.\n");
    } else {
        if( cps->set->debug>1 || debug>1 || local_debug ) fprintf(stderr, "\n");
    }
    ' Clear temporary partitions
    if( cps->set->alawadhi) {
        for(i=0; i<num_affected_regions; i++) {
             PartitionDelete(p1[i]);
             PartitionDelete(p2[i]);
        }
        free(p1);
        free(p2);
    }
    PartitionDelete(new_part);
}' ProposeDeleteSegment

static void CPSAddPartition(cpsampler *cps, partition *new_part, int insert_after) {
     partition **new_part_list = (partition **) malloc(sizeof(partition *)*(cps->npartitions+1));
     int i;
     for(i=0; i<cps->npartitions; i++) {
          if( i<insert_after )
               new_part_list[i] = cps->part_list[i];
          Else
               new_part_list[i+1] = cps->part_list[i];
     }
      new_part_list[insert_after] = new_part;
      cps->npartitions++;
      free(cps->part_list);
      cps->part_list = new_part_list;
}' CPSAddPartition

' WORKING: to update with new partition setup
static void CPSRemovePartition(cpsampler *cps, int del_segment) {
     partition **new_part_list = (partition **) malloc(sizeof(partition *)*(cps->npartitions-1));
     int i;
     for(i=0; i<cps->npartitions-1; i++) {
        if( i < del_segment )
             new_part_list[i] = cps->part_list[i];
        Else
             new_part_list[i] = cps->part_list[i+1];
     }
    PartitionDelete(cps->part_list[del_segment]);
    cps->npartitions--;
    free(cps->part_list);
    cps->part_list = new_part_list;
}' CPSRemovePartition

static int determineSameT(const int np, partition *p[]) {
    int a = 0, i, lastTree = p[0]->cTree;
    for(i=1; i<np; i++) {
        int thisTree = p[i]->cTree;
        if( thisTree == lastTree ) a++;
        else lastTree = thisTree;
    }
    return a;
}' determineSameT

static int updateSameT(const cpsampler *cps, const int pTree, const int i) {
    int update = 0;

    if( pTree == cps->part_list[i]->cTree ) return update;
    if( i > 0 ) {
        if( (cps->part_list[i-1]->cTree == cps->part_list[i]->cTree) && (cps->part_list[i-1]->cTree != pTree) ) update--;
        else if( (cps->part_list[i-1]->cTree != cps->part_list[i]->cTree) && (cps->part_list[i-1]->cTree == pTree) ) update++;
    }
    if( i < (cps->npartitions - 1) ) {
        if( (cps->part_list[i+1]->cTree == cps->part_list[i]->cTree) && (cps->part_list[i+1]->cTree != pTree) ) update--;
        else if( (cps->part_list[i+1]->cTree != cps->part_list[i]->cTree) && (cps->part_list[i+1]->cTree == pTree) ) update++;
    }
    return update;
}' updateSameT

static double CPLogJacobian(double alpha, double alpha1, double alpha2, double sigmaa, double mu, double mu1, double mu2, double sigmam) {
    return log( alpha1 ) + log( 1.0 - alpha1 ) + log( alpha2 ) + log( 1.0 - alpha2 ) + log( mu1 ) + log( mu2 ) + log( sigmaa ) + log( sigmam )
        - log( alpha ) - log( 1.0 - alpha ) - log( mu );
}' CPLogJacobian

static void PrintInitialValues(const cpsampler *cps) {
    int i, j;
    char tmp[MAX_TREE_STRING];  ' Potential BUG
    printf("\nTotal Starting Regions: %d", cps->npartitions);
    for(i=0; i<cps->npartitions; i++) {
        printf("\n\nREGION #%d\n", i+1);
        for(j=0;j<4;j++) printf("%f ", cps->part_list[i]->cmatrix->pi[j]);
        printf("\n");
        printf("Starting Mu & Alpha: %lf %lf\n", cps->part_list[i]->cHyperParameter, cps->part_list[i]->cmatrix->v[0]);
        tmp[0] = '\0';
        toString(tmp, cps->PostTree[cps->part_list[i]->cTree].root, false);
        printf("Starting Tree: %s\n", tmp);
        printf("Starting likelihood: %lf\n", cps->part_list[i]->cPartialLogLikelihood);
    }
    printf("\n\n");
}' PrintInitialValues

static double DrawnAlpha1(double alpha, double weight, double sigma, double z) {
    double pert = exp(sigma * weight * z);
    return alpha * pert / ( 1.0 - alpha + alpha * pert);
}' DrawnAlpha1

static double DrawnAlpha2(double alpha, double weight, double sigma, double z) {
    double pert = exp(-sigma * weight * z);
    return alpha * pert / ( 1.0 - alpha + alpha * pert);
}' DrawnAlpha2

static double DrawnMu1(double mu, double weight, double sigma, double z) {
    return mu * exp(weight * sigma * z);
}' DrawnMu1

static double DrawnMu2(double mu, double weight, double sigma, double z) {
    return mu * exp(- weight * sigma * z);
}' DrawnMu2

static double CondensedAlpha(double alpha1, double weight1, double alpha2, double weight2 ) {
    return InvLogit( weight1 * Logit(alpha1) + weight2 * Logit(alpha2) );
}' CondensedAlpha

static double CondensedMu(double mu1, double weight1, double mu2, double weight2) {
    return exp( weight1 * log(mu1) + weight2 * log(mu2) );
}' CondensedMu

static double InverseZMu(double mu1, double mu0, double weight1, double sigma) {
    return log( mu1 / mu0 ) / weight1 / sigma;
}' InverseZMu

static double InverseZAlpha( double alpha1, double alpha0, double weight1, double sigma) {
    return log( alpha1 * (1.0 - alpha0) / (1.0 - alpha1) / alpha0 ) / weight1 / sigma;
}'InverseZAlpha

/**
* Computes the logit transform of x.
*/
static double Logit(const double x) {
    return log( x / (1.0 - x) );
}' Logit

/**
* Computes the inverse of the logit transform of x.
*/
static double InvLogit(const double x) {
    double tmp = exp(x);
    return tmp / (1.0 + tmp);
}' InvLogit

static void CPRecordStatistics(cpsampler *cps) {
    int i;
    avg_nchpts += cps->npartitions-1;
    std_nchpts += (cps->npartitions-1)*(cps->npartitions-1);
    sample_nchpts++;
    for( i=0;i<cps->npartitions;i++ ) {
        avg_mu += cps->part_list[i]->cHyperParameter;
        std_mu += cps->part_list[i]->cHyperParameter * cps->part_list[i]->cHyperParameter;
        avg_alpha += cps->part_list[i]->cmatrix->v[0];
        std_alpha += cps->part_list[i]->cmatrix->v[0] * cps->part_list[i]->cmatrix->v[0];
        sample_param++;
    }
    if( cps->npartitions < 11) {
            poisson[cps->npartitions-1]++;
        sample_poisson++;
    }
}' CPRecordStatistics

' dcpsampler.c
#include "dcpsampler.h"

static int debug = 1;               ' Set to positive integer for local debugging output
'static const char *file_name = "dcpsampler.c";

' debugging variables
static double avg_ntopol = 0;
static double std_ntopol = 0;
static double avg_nparam = 0;
static double std_nparam = 0;
static int cnt = 0;
static double top_poisson[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
static double par_poisson[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
static double first_posn_tree[3] = {0.0,0.0,0.0};
static int mono_prob = 0;
static double mean_log_mu = 0.0;
static double mean_log_kappa = 0.0;
static double sd_log_mu = 0.0;
static double sd_log_kappa = 0.0;
static double *rolling_ntopol_avg = NULL;
static double *rolling_nparam_avg = NULL;
static int rolling_size = 10000;
static int rolling_index = 999; ' Set to rolling_size - 1
static double iact_J = 0;
static int *record_J = NULL;
static int mean_J = 0;
static int cnt_J = 0;
static double iact_kappa = 0;
static double *record_kappa = NULL;
static double mean_kappa = 0.0;
static int cnt_kappa = 0;

' Local function predeclarations:
static void RecordStatistics(dcpsampler *);
static void DCPPartitionListMakeFrom(partition_list **, char *, sampler *);
static void DCPPartitionListInitialize(dcpsampler *);
static void DCPReportState(sampler *, const char *, double, double);
static void DCPReportProposedState(sampler *, const char *, tree **, int);
static void DCPReportStatistics(sampler *, const char *, double, double, double, double, double, int, int, int);
static double DCPLogJacobian(const sampler *, ...);
static void DCPRun(sampler *);
static void DCPOutputLine(const sampler *);
static void Fixed_Dimension_Sampler(sampler *, int, int);
static void Alawadhi_Copy_State(sampler *, ...);
static void Alawadhi_Accept(sampler *);
static void Alawadhi_Reject(sampler *);
static double Log_Prior(const sampler *, boolean);
static void Exit_Condition(const dcpsampler *);
static void Compute_Auto_Correlation(const dcpsampler *);

void DCPSamplerSetup(dcpsampler **in_dcp, seqdata *sqd, settings *set, char *ofilename) {
    const char *fxn_name = "DCPSamplerSetup";
    int i;
    sampler *smp;
    dcpsampler *dcp;

    *in_dcp = (dcpsampler *) malloc(sizeof(dcpsampler));
    if( *in_dcp == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    dcp = *in_dcp;

    dcp->set = set; ' For ease of access
    dcp->mc = 0.0;
    dcp->alawadhi_part_list = NULL;
    dcp->part_list = NULL;

    ' Setup base sampler object (things common to all samplers ... not much!)
    SamplerMake(&dcp->smp);

    ' Set up object pointers (reach right into the guts of the struct!!)
    smp = dcp->smp;
    smp->derived_smp = (void *)dcp;
    smp->set = set;
    smp->rng = set->rng;
    smp->sqd = sqd;

    ' Set up the move types and names
    smp->max_move_name_length = strlen("parameter update");
    SamplerSetNumberMoves(smp, 11);
    SamplerAddMoveName(smp, UPDATE_TAU, "topology update");
    SamplerAddMoveName(smp, UPDATE_KAPPA_AND_MU, "parameter update");   ' longest name
    SamplerAddMoveName(smp, UPDATE_XI, "xi update");
    SamplerAddMoveName(smp, UPDATE_RHO, "rho update");
    SamplerAddMoveName(smp, ADD_XI, "add one xi");
    SamplerAddMoveName(smp, ADD_TWO_XI, "add two xi");
    SamplerAddMoveName(smp, ADD_RHO, "add rho");
    SamplerAddMoveName(smp, DELETE_XI, "remove one xi");
    SamplerAddMoveName(smp, DELETE_TWO_XI, "remove two xi");
    SamplerAddMoveName(smp, DELETE_RHO, "remove rho");
    SamplerAddMoveName(smp, FIXED_DIMENSION, "fixed dimension");

    ' Setup function pointers
    smp->OutputLine = DCPOutputLine;
    smp->logJacobian = DCPLogJacobian;
    smp->run = DCPRun;
    smp->Fixed_Dimension_Sampler = &Fixed_Dimension_Sampler;
    smp->Alawadhi_Copy_State = &Alawadhi_Copy_State;
    smp->Alawadhi_Accept = &Alawadhi_Accept;
    smp->Alawadhi_Reject = &Alawadhi_Reject;
    smp->Log_Prior = &Log_Prior;
    smp->Report_State = &DCPReportState;
    smp->Report_Proposed_State = &DCPReportProposedState;
    smp->Report_Proposal_Statistics = &DCPReportStatistics;

    if( set->ctmc_model == HKY ) {
        smp->Matrix_Make_Initial = &iHKYNoBoundFixPiMatrixMakeInitial;
        smp->Matrix_Make_Default = &iHKYNoBoundFixPiMatrixMakeDefault;
        smp->Matrix_Make_and_Set = &iHKYNoBoundFixPiMatrixMakeAndSet;
        set->ctmc_parameterization = KAPPA; ' Only current choice for dcpsampler
    }

    ' Open output file
    smp->fout = fopen(ofilename, "w");
    if( !dcp->smp->fout ) {
        fprintf(stderr, "ERROR:  Could not open output file %s\n", ofilename);
        exit(EXIT_FAILURE);
    }

    ' DEBUG output
    rolling_ntopol_avg = (double *) malloc(sizeof(double)*rolling_size);
    rolling_nparam_avg = (double *) malloc(sizeof(double)*rolling_size);
    for( i=0; i<rolling_size; i++ ) {
        rolling_ntopol_avg[i] = 0.0;
        rolling_nparam_avg[i] = 0.0;
    }

    record_J = (int *) malloc(sizeof(int)*set->length);
    record_kappa = (double *) malloc(sizeof(double)*set->length);

    ' Make tree vector object that will handle update the full vector \tau (will hold current state information and, so will be initialized below)
    TreeVectorMake(&dcp->tree_vec, smp);

    ' Make branch object to handle hierarchical dependence in topology branch lengths across segments (contains no current state information)....
    ' Because there is no current state info stored within, the prior object information is initalized once here.
    ' If branch structure becomes more complex, model after CTMC matrix setup (below) or tree vector setup (above).
    BranchMake(&dcp->smp->br, smp->set);

    ' Make an HKY-specific global qmatrix_prior object to handle hierarchical dependence across segments.
    ' Will hold current state information and, so will be initialized later (as part of Partition_List initialization; OK there because of conditional independence)
    ' Per-segment information is stored in individual ihkynoboundfixpimatrix objects.
    ' Global hierarhical information is stored in a single qmatrix_prior object which the individual objects access and manipulate smartly (i.e. they know its hidden structure)
    ' to calculate their own priors.
    if( set->ctmc_model == HKY ) iHKYNoBoundFixPiMatrixGlobalInitialize(dcp->smp->set);

    if( !set->init_string ) {
        DCPPartitionListInitialize(dcp);
    } else {
        DCPPartitionListMakeFrom(&dcp->part_list, set->init_string, smp);'ParseLastLine(dcp);  ' Prepares partitions
        TreeVectorInitialize(dcp->tree_vec, dcp->part_list);
    }

    ' Override default partition_list choices to retain required structure on topologies across partitions
    dcp->part_list->Propose_Topology_Change_Point_To_Delete = &Propose_Topology_Change_Point_To_Delete;
    dcp->part_list->Propose_Two_Topology_Change_Points_To_Delete = &Propose_Two_Topology_Change_Points_To_Delete;

    if( smp->JumpNumber > set->burnin ) SamplerSaveEstimates(smp, set->length);

}' DCPSamplerSetup

static void DCPOutputLine(const sampler *smp) {
    char tmp[MAX_TREE_STRING];' Potential BUG - ksd
    const dcpsampler *dcp = (dcpsampler *) smp->derived_smp;
    int i;

    fprintf(dcp->smp->fout, "%-6d%3d%3d", smp->JumpNumber, dcp->part_list->npartitions, dcp->part_list->npartitions - 1 - dcp->part_list->topology_changes);
    for(i=0; i<dcp->part_list->npartitions; i++) {
        const partition *part = dcp->part_list->part[i];
        tmp[0] = '\0';
        toString(tmp, part->ctree->root, false);
        fprintf(dcp->smp->fout," %s %10.2f %5.4f %5.4f %5.4f %5.4f %5.4f %7.4f",
            tmp,
            part->cPartialLogLikelihood,
            part->cmatrix->v[0],
            part->cmatrix->pi[0],part->cmatrix->pi[1],part->cmatrix->pi[2],part->cmatrix->pi[3],
            part->cHyperParameter);
    }
    if( dcp->set->jump_classes ) for(i=0; i< dcp->part_list->npartitions; i++) fprintf(dcp->smp->fout," %d", dcp->part_list->part[i]->left);
    fprintf(dcp->smp->fout, "\n");
}'method DCPOutputLine

static void DCPRun(sampler *smp) {
    const char *fxn_name = "DCPRun";
    dcpsampler *dcp = (dcpsampler *)smp->derived_smp;
    int sincePrint = 0;
    double BirthOrDeathOrMove = 0.0;
    settings *set = dcp->set;

    while(true) {
        dcp->part_list->Update_Move_Probabilities(dcp->part_list, set);

        if( dcp->set->jump_classes )
            BirthOrDeathOrMove = set->rng->nextStandardUniform(set->rng);
        Else    ' Only fixed-dimension moves
            BirthOrDeathOrMove = 1.0;

        if( debug>1 || global_debug>1 )
            fprintf(stderr, "Next move (%d): %.4f -> ", dcp->set->rng->useRnList?dcp->set->rng->current_rn:dcp->smp->JumpNumber, BirthOrDeathOrMove);

        dcp->smp->tries[FIXED_DIMENSION]++;

        if( BirthOrDeathOrMove < dcp->part_list->top_one_bk ) {
            if( debug>1 || global_debug>1 ) fprintf(stderr, "TopologyAddOne\n");
            if( dcp->part_list->topology_changes < dcp->smp->sqd->lenseq - 1 ) {
                dcp->tree_vec->Add_One(dcp->tree_vec, dcp->part_list, dcp->smp);
            }
        } else if( BirthOrDeathOrMove < dcp->part_list->top_one_bk + dcp->part_list->top_two_bk ) {
            if( debug>1 || global_debug>1 ) fprintf(stderr, "TopologyAddTwo\n");
            if( dcp->part_list->topology_changes < dcp->smp->sqd->lenseq - 2 ) {
                    dcp->tree_vec->Add_Two(dcp->tree_vec, dcp->part_list, dcp->smp);
            }
        } else if( BirthOrDeathOrMove < dcp->part_list->top_one_bk + dcp->part_list->top_two_bk + dcp->part_list->top_one_dk ) {
            if( debug>1 || global_debug>1 ) fprintf(stderr, "TopologyDeleteOne\n");
            if( dcp->part_list->topology_changes > 0 ) {
                    dcp->tree_vec->Delete_One(dcp->tree_vec, dcp->part_list, dcp->smp);
            }
        } else if( BirthOrDeathOrMove < dcp->part_list->top_one_bk + dcp->part_list->top_two_bk + dcp->part_list->top_one_dk + dcp->part_list->top_two_dk ) {
            if( debug>1 || global_debug>1 ) fprintf(stderr, "TopologyDeleteTwo\n");
            if( dcp->part_list->topology_changes > 1 ) {
                    dcp->tree_vec->Delete_Two(dcp->tree_vec, dcp->part_list, dcp->smp);
            }
        } else if( BirthOrDeathOrMove < dcp->part_list->top_one_bk + dcp->part_list->top_two_bk + dcp->part_list->top_one_dk + dcp->part_list->top_two_dk + dcp->part_list->par_bk ) {
            if( debug>1 || global_debug>1 ) fprintf(stderr, "ParameterAddOne\n");
            if( dcp->part_list->parameter_changes < dcp->smp->sqd->lenseq - 1 ) {
                    ParChgPtAdd(dcp->part_list, dcp->smp);
            }
        } else if( BirthOrDeathOrMove < dcp->part_list->top_one_bk + dcp->part_list->top_two_bk + dcp->part_list->top_one_dk + dcp->part_list->top_two_dk + dcp->part_list->par_bk + dcp->part_list->par_dk ) {
            if( debug>1 || global_debug>1 ) fprintf(stderr, "ParameterDeleteOne\n");
            if( dcp->part_list->parameter_changes > 0 ) {
                    ParChgPtDelete(dcp->part_list, dcp->smp);
            }
        } else {
            if( debug>1 || global_debug>1 ) fprintf(stderr, "FixedDimensionSampler\n");
            dcp->smp->acceptancerate[FIXED_DIMENSION]++;

            smp->Fixed_Dimension_Sampler(dcp->smp, 1, 0);
        }

        sincePrint++;
        dcp->smp->JumpNumber++;
        if( dcp->smp->JumpNumber > dcp->set->burnin ) RecordStatistics(dcp);
        if( sincePrint >= dcp->set->subsample ) {
            sincePrint = 0;
            if( dcp->smp->JumpNumber > dcp->set->burnin ) {
                PrintTopologies(dcp->smp, fxn_name, false, false);
                SamplerSaveEstimates(dcp->smp, dcp->set->length);
            }
        }
        if( dcp->set->exit_condition ) Exit_Condition(dcp);
    }

}' DCPRun

static void Fixed_Dimension_Sampler(sampler *smp, int ntimes, int alawadhi) {
    dcpsampler *dcp = (dcpsampler *)smp->derived_smp;
    partition_list *pl = alawadhi ? dcp->alawadhi_part_list : dcp->part_list;
    topology_vector *tv = dcp->tree_vec;
    int i;

    for( i=0; i<ntimes; i++ ) {
        if( !alawadhi || alawadhi&KAPPA_MU ) UpdateParameters(pl, dcp->smp, alawadhi);
        if( tv->numTrees > 1 && ( !alawadhi || alawadhi&TAU ) ) tv->Update_Topologies(tv, pl, dcp->smp, alawadhi);
            
'      if( set->update_hyperparameters ) __INSERT_SOMETHING__;

        if( pl->parameter_changes > 0 && ( !alawadhi || alawadhi&RHO ) ) UpdateChangePointLocations(pl, dcp->smp, false, alawadhi);
        if( pl->topology_changes > 0 && ( !alawadhi || alawadhi&XI ) ) UpdateChangePointLocations(pl, dcp->smp, true, alawadhi);
    }
}' Fixed_Dimension_Sampler

static void Alawadhi_Copy_State(sampler *smp, ...) {
    dcpsampler *dcp = (dcpsampler *) smp->derived_smp;
    partition_list **tpl = NULL;
    va_list vargs;

    va_start(vargs, smp);
    tpl = va_arg(vargs, partition_list **);
    va_end(vargs);

    PartitionListMakeCopy(&dcp->alawadhi_part_list, dcp->part_list);
    *tpl = dcp->alawadhi_part_list;
}' Alawadhi_Copy_State

static void Alawadhi_Accept(sampler *smp) {
    dcpsampler *dcp = (dcpsampler *) smp->derived_smp;

    PartitionListDelete(dcp->part_list, true);
    dcp->part_list = dcp->alawadhi_part_list;
    dcp->alawadhi_part_list = NULL;
}' Alawadhi_Accept

static void Alawadhi_Reject(sampler *smp) {
    dcpsampler *dcp = (dcpsampler *) smp->derived_smp;

    PartitionListDelete(dcp->alawadhi_part_list, true);
    dcp->alawadhi_part_list = NULL;
    TreeVectorInitialize(dcp->tree_vec, dcp->part_list);
}' Alawadhi_Reject

Public Function Log_Prior(smp As sampler, alawadhi As Byte) As Double


    
    'dcpsampler *dcp = (dcpsampler *) smp->derived_smp;
    dcp As dcpsampler
    dcp = smp.derived_smp
    
    'partition_list *pl = alawadhi ? dcp->alawadhi_part_list : dcp->part_list;
    pl As partition_list
    
    If alawadhi = 1 Then
        pl = dcp.alawadhi_part_list
    Else
        dcp.part_list
    End If
    
    'topology_vector *tv = dcp->tree_vec;
    tv As topology_vector
    tv = dcp.tree_vec
    
    Dim i As Long
    Dim mup As Double, kap As Double, lp As Double
    mup = 0: kap = 0
    lp = tv.log_prior_prob + pl->K_Log_Prior(pl) + pl->J_Log_Prior(pl); ' tau, K, J
'fprintf(stderr, "Log_Prior(\\tau): %.4f\n", tv->log_prior_prob); fprintf(stderr, "Log_Prior(   K): %.4f\n", pl->K_Log_Prior(pl)); fprintf(stderr, "Log_Prior(   J): %.4f\n", pl->J_Log_Prior(pl));
    lp += pl->Xi_Log_Prior(pl) + pl->Rho_Log_Prior(pl);             ' xi, rho
'fprintf(stderr, "Log_Prior( \\xi): %.4f\n", pl->Xi_Log_Prior(pl)); fprintf(stderr, "Log_Prior(\\rho): %.4f\n", pl->Rho_Log_Prior(pl));
    for( i=0; i<pl->npartitions; i++ ) {
        if( !pl->part[i]->parchange ) continue;
        lp += pl->part[i]->cmatrix->log_prior;                  ' kappa
        lp += pl->part[i]->cPartialLogHyperParameterPrior;          ' mu
        kap += pl->part[i]->cmatrix->log_prior;
        mup += pl->part[i]->cPartialLogHyperParameterPrior;
'fprintf(stderr, "Log_Prior(\\kappa[%d]): %.4f\n", i, pl->part[i]->cmatrix->log_prior); fprintf(stderr, "Log_Prior(\\mu[%d]): %.4f\n", i, pl->part[i]->cPartialLogHyperParameterPrior);
    }
'fprintf(stderr, "Log_Prior(\\kappa): %.4f\n", kap); fprintf(stderr, "Log_Prior(\\mu): %.4f\n", mup);
    return lp;
End Function ' Log_Prior

static double DCPLogJacobian(const sampler *smp, ...) {
    double alpha1, alpha2, alpha;
    qmatrix *cmatrix, *pmatrix1, *pmatrix2;
    double sigmaa, sigmam;
    double mu, mu1, mu2;
    va_list vargs;

    va_start(vargs, smp);
    cmatrix = va_arg(vargs, qmatrix *);
    pmatrix1 = va_arg(vargs, qmatrix *);
    pmatrix2 = va_arg(vargs, qmatrix *);
    sigmaa = va_arg(vargs, double);
    mu = va_arg(vargs, double);
    mu1 = va_arg(vargs, double);
    mu2 = va_arg(vargs, double);
    sigmam = va_arg(vargs, double);
    va_end(vargs);

    alpha1 = pmatrix1->v[0];
    alpha2 = pmatrix2->v[0];
    alpha = cmatrix->v[0];
    return log(alpha1) + log(alpha2) + log(mu1) + log(mu2) + log(sigmaa) + log(sigmam) - log(alpha) - log(mu);
}' DCPLogJacobian

' DCP-specific functions for creating partition_list and initializing partition_list
static void DCPPartitionListMakeFrom(partition_list **pl, char *str, sampler *smp) {
    const char *fxn_name = "DCPPartitionListMakeFrom";
    char *tmp;
    int count = 0, current_region = 0, var_num = 0, i;
    dcpsampler *dcp = (dcpsampler *) smp->derived_smp;
    partition *part = NULL;

    tmp = strtok(str, " ");

    while( tmp != NULL ) {
        ' Count tokens
        count++;
        ' First word is the sample number
        if( count == 1) smp->JumpNumber = atoi(tmp);
        ' Second word is the number of segments
        else if( count == 2 ) {
            int nparts = atoi(tmp);
            PartitionListMake(pl, smp, nparts);
            if( *pl == NULL ) {
                fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                exit(EXIT_FAILURE);
            }
            ' Allocate memory for all other parameters
            PartitionMake(&((*pl)->part[0]), smp->sqd->lenunique, 0, 0, false, false);
            part = (*pl)->part[0];
            smp->Matrix_Make_Default(&part->cmatrix, KAPPA);
        }
        ' The number evolutionary parameter change-points
        else if( count == 3 ) {
            (*pl)->topology_changes = (*pl)->npartitions - 1 - atoi(tmp);
            (*pl)->parameter_changes = 0;   ' Computed later
        }
        ' Segment-specific parameters
        else if( count <= 3 + (*pl)->npartitions * 8 ) {
            ' topology
            if( var_num % 8 == 0 ) {
                tree *ttree = NULL;
                Make_Tree(&ttree, tmp, smp->sqd->num_chars);
                Balance_Tree(ttree);
                ' Find the matching tree in the list of possibilities
                for(i=0; i<dcp->tree_vec->numTrees; i++) {
                    if( SameTrees(dcp->tree_vec->tree_list[i], ttree, false) ) {    ' WAS BUG !SameTrees
                        part->ctree = dcp->tree_vec->tree_list[i];
                        break;
                    }
                }
                if(ttree) free(ttree);
            ' log likelihood
            } else if( var_num % 8 == 1 ) {
                part->cPartialLogLikelihood = atof(tmp);
            ' kappa
            } else if ( var_num % 8 == 2 ) {    ' WAS BUG: out of order
                Set_Kappa(part->cmatrix, atof(tmp));
            ' pi_A
            } else if ( var_num % 8 == 3 ) {
                part->cmatrix->pi[0] = atof(tmp);
            ' pi_C
            } else if ( var_num % 8 == 4 ) {
                part->cmatrix->pi[1] = atof(tmp);
            ' pi_G
            } else if ( var_num % 8 == 5 ) {
                part->cmatrix->pi[2] = atof(tmp);
            ' pi_T
            } else if ( var_num % 8 == 6 ) {
                part->cmatrix->pi[3] = atof(tmp);
            ' mu
            } else if ( var_num % 8 == 7 ) {
                part->cHyperParameter = atof(tmp);
                part->cPartialLogHyperParameterPrior = dcp->smp->br->Log_Prior(dcp->smp->br, part->cHyperParameter);
                current_region++;
                if( current_region < (*pl)->npartitions ) { ' WAS BUG
                    PartitionMake(&((*pl)->part[current_region]), smp->sqd->lenunique, 0, 0, false, false);
                    part = (*pl)->part[current_region];
                    smp->Matrix_Make_Default(&part->cmatrix, KAPPA);
                }
            }
            var_num++;
        }
        ' The last words are the change-point locations
        else {
            if(current_region >= (*pl)->npartitions) current_region = 0;
            (*pl)->part[current_region]->left = atoi(tmp);
            if(current_region) (*pl)->part[current_region-1]->right = (*pl)->part[current_region]->left - 1;
            current_region++;
        }
        ' Get next token
        tmp = strtok(NULL, " ");
    }
    ' Set the right end of the last segment
    (*pl)->part[(*pl)->npartitions-1]->right = smp->sqd->lenseq - 1;

    ' Set up the partitions boundaries, count data, types of change points, and likelihoods
    for(i=0; i<(*pl)->npartitions; i++) {
        partition *part = (*pl)->part[i];
        boolean btc = i ? (part->ctree != (*pl)->part[i-1]->ctree) : true;  ' WAS BUG ==
        boolean bpc = i
            ? (fabs(part->cHyperParameter - (*pl)->part[i-1]->cHyperParameter) > tolerance || fabs(part->cmatrix->v[0] - (*pl)->part[i-1]->cmatrix->v[0]) > tolerance)  ' WAS BUG?
            : true;
        if(bpc && i) (*pl)->parameter_changes++;    ' WAS BUG && i added

        PartitionReset(part, part->left, part->right, btc, bpc);
        PartitionCopySegmentCounts(part, smp->sqd, part->left, part->right+1);
        part->cmatrix->Matrix_Sync(part->cmatrix);

        part->cPartialLogHyperParameterPrior = dcp->smp->br->Log_Prior(dcp->smp->br, part->cHyperParameter);
    }
    PrintTopologies(smp, "DCPPartitionListMakeFrom", false, false);
'  VerifyCounts(smp, "DCPPartitionListMakeFrom", false);
'  VerifyLikelihood(smp, false);
}' DCPPartitionListMakeFrom

void DCPPartitionListInitialize(dcpsampler *dcp) {
    int i;
    sampler *smp = (sampler *)dcp->smp;
    char tmp[MAX_TREE_STRING];  ' BUGGY!
    partition_list *pl = NULL;
    'int bp[3] = {0, 1375, 2221};

    ' Make partition list (our initial distribution puts all weight on 1 segment)
    'PartitionListMake(&pl, smp, 3);   ' DEBUG
    PartitionListMake(&pl, smp, 1);
    dcp->part_list = pl;

    'dcp->part_list->topology_changes = 2; ' DEBUG

    ' Make and initialize each partition (set tree, \mu, qmatrix, likelihood)
    for( i=0; i<pl->npartitions; i++ ) {
        partition *part = NULL;
        PartitionMake(&(pl->part[i]), smp->sqd->lenunique, 0, smp->sqd->lenseq - 1, true, true);    ' DEBUG
        'if(i<2) PartitionMake(&(pl->part[i]), smp->sqd->lenunique, bp[i], bp[i+1] - 1, true, false);  ' DEBUG
        'else PartitionMake(&(pl->part[i]), smp->sqd->lenunique, bp[i], smp->sqd->lenseq - 1, true, false);    ' DEBUG
        part = pl->part[i];
        'if(!i) part->parchange = true;    ' DEBUG
        PartitionCopyCounts(part, smp->sqd);    ' DEBUG
        'if(i<2) PartitionCopySegmentCounts(part, smp->sqd, bp[i], bp[i+1]);
        'else PartitionCopySegmentCounts(part, smp->sqd, bp[i], smp->sqd->lenseq - 1);
        part->ctree = dcp->tree_vec->Draw_Initial_Tree(dcp->tree_vec, smp->set);' Only works if \tau_i depends at most on \tau_j, j<i  ' DEBUG
        'if(!i || i==2) part->ctree = dcp->tree_vec->tree_list[4]; ' DEBUG
        'if(i==1) part->ctree = dcp->tree_vec->tree_list[2];   ' DEBUG
        tmp[0] = '\0';
        toString(tmp, part->ctree->root, false);
fprintf(stderr, "Starting tree: %s\n", tmp);
        part->cHyperParameter = smp->br->Initialize(smp->set);          ' Only works if \mu_i depends at most on \mu_j, j<i
fprintf(stderr, "Starting hyperparameter: %f\n", part->cHyperParameter);
        part->cPartialLogHyperParameterPrior = smp->br->Log_Prior(smp->br, part->cHyperParameter);
fprintf(stderr, "Starting log hyperparameter: %f\n", part->cPartialLogHyperParameterPrior);
        smp->Matrix_Make_Initial(&(part->cmatrix), KAPPA, smp->sqd, smp->set);          ' Only works if CTMC matrices are independent (for some reason hierarchical struct handled diff)
fprintf(stderr, "Starting ep: %f (%.4f)\n", part->cmatrix->v[0], part->cmatrix->log_prior);

        ' By now, everything should be set to compute likelihood
        part->cPartialLogLikelihood = TreeLogLikelihood(part->ctree, smp, part->cmatrix, part->counts, part->cHyperParameter, false);
fprintf(stderr, "Starting likelihood: %f\n", part->cPartialLogLikelihood);
    }
    VerifyLikelihood(dcp->smp, false);      ' DEBUG
    VerifyCounts(dcp->smp, "DCPPartitionListInitialize", false);        ' DEBUG

    ' Initialize vector \tau now that the trees have been selected for each region
    ' In fact, this would handle tree selection if there were multiple regions and a complex dependence structure
    TreeVectorInitialize(dcp->tree_vec, dcp->part_list);
fprintf(stderr, "Starting prior on tree: %.4f\n", dcp->tree_vec->log_prior_prob);

}' DCPPartitionListInitialize

void DCPSamplerDelete(dcpsampler *dcp) {
    int i;
    if(!dcp) return;
    if(dcp->part_list) {
        for( i=0; i<dcp->part_list->npartitions; i++ ) {
            QMatrixDelete(dcp->part_list->part[i]->cmatrix);
            PartitionDelete(dcp->part_list->part[i]);
        }
        PartitionListDelete(dcp->part_list, true);
    }
    if(dcp->tree_vec) {
        TopologyVectorDelete(dcp->tree_vec, true);
    }
    if(dcp->smp->br) BranchDelete(dcp->smp->br);
    free(dcp);
}' DCPSamplerDelete

' Debugging functions
static void RecordStatistics(dcpsampler *dcp) {
    int i;
    boolean all_mono = true;
    topology_gmodel_prior *tgp = NULL;

    if( dcp->set->gmodel ) tgp = (topology_gmodel_prior *) dcp->tree_vec->top_prior;

    for( i=0; i<(dcp->part_list->topology_changes+1); i++ ) if( dcp->set->gmodel && !tgp->monophyletic[dcp->tree_vec->current_trees[i]->tree_index] ) all_mono = false;
    for( i=0; i<dcp->part_list->npartitions; i++ ) {
        if( dcp->part_list->part[i]->cPartialLogLikelihood - 0.0 > tolerance ) {
            PrintTopologies(dcp->smp, "RecordStatistics", false, false);
            fprintf(stderr, "ERROR: positive log likelihood %f in partition %d.\n", dcp->part_list->part[i]->cPartialLogLikelihood, i);
            exit(EXIT_FAILURE);
        }
    }
    if( all_mono ) mono_prob ++;
    mean_log_mu += log(dcp->part_list->part[0]->cHyperParameter);
    mean_log_kappa += log(dcp->part_list->part[0]->cmatrix->v[0]);
    sd_log_mu += log(dcp->part_list->part[0]->cHyperParameter) * log(dcp->part_list->part[0]->cHyperParameter);
    sd_log_kappa += log(dcp->part_list->part[0]->cmatrix->v[0]) * log(dcp->part_list->part[0]->cmatrix->v[0]);
    avg_ntopol += dcp->part_list->topology_changes;
    std_ntopol += dcp->part_list->topology_changes*dcp->part_list->topology_changes;
    avg_nparam += dcp->part_list->parameter_changes;
    std_nparam += dcp->part_list->parameter_changes*dcp->part_list->parameter_changes;
    rolling_ntopol_avg[rolling_index] = 0.0;
    rolling_nparam_avg[rolling_index] = 0.0;
    if( ++rolling_index == rolling_size ) rolling_index = 0;
    for( i=0; i<rolling_size; i++ ) {
        rolling_ntopol_avg[i] += dcp->part_list->topology_changes;
        rolling_nparam_avg[i] += dcp->part_list->parameter_changes;
    }
    if( dcp->part_list->topology_changes < 10 ) {
        top_poisson[dcp->part_list->topology_changes]++;
    }
    if( dcp->part_list->parameter_changes < 10 ) {
        par_poisson[dcp->part_list->parameter_changes]++;
    }
    if(dcp->part_list->part[0]->ctree->tree_index < 3) first_posn_tree[dcp->part_list->part[0]->ctree->tree_index]++;
    cnt++;
    dcp->mc = (double)mono_prob/cnt;
    record_J[cnt_J] = dcp->part_list->parameter_changes;
    mean_J += record_J[cnt_J++];
    record_kappa[cnt_kappa] = dcp->part_list->part[0]->cmatrix->v[0];
    mean_kappa += record_kappa[cnt_kappa++];
}' RecordStatistics

static void Exit_Condition(const dcpsampler *dcp) {
    int i, j;
    partition_list *pl = dcp->part_list;

    ' EDIT: to stop on trigger
    const int num_cops = 3;
    const double true_cops[3] = {500,1000,1500};
    const double cop_interval[3] = {30,30,30};
    const boolean topology_change_points = false;
    boolean got_cop[3] = {false,false,false};
    ' EDIT: end
    boolean got_all_cops = true;
    for( i=1; i<pl->npartitions; i++ ) {
        if( (topology_change_points && pl->part[i]->topchange) || (!topology_change_points && pl->part[i]->parchange) ) {
            for( j=0; j<num_cops; j++ ) {
                if( abs( pl->part[i]->left - true_cops[j] ) < cop_interval[j] ) got_cop[j] = true;
            }
        }
    }
    for( j=0; j<num_cops; j++ ) if( !got_cop[j] ) got_all_cops = false;
    if( got_all_cops ) PrintTopologies(dcp->smp, "Exit_Condition", false, true);
}' Exit_Condition

static void Compute_Auto_Correlation(const dcpsampler *dcp) {
    double mean = (double) mean_J / cnt_J;
    double mean_k = mean_kappa / cnt_kappa;
    int lag_size = 100000;
    double var = 0, lag = 0, lag_kappa, var_kappa=0;
    int current_len = dcp->smp->JumpNumber;
    int i, k;

    if( lag_size > current_len ) lag_size = current_len;
    iact_J = 1.0;
    iact_kappa = 1.0;
    for( k=0; k<lag_size; k++ ) {
        lag = 0;
        lag_kappa = 0;
        for( i=0; i<current_len - k; i++ ) {
            lag += (record_J[i] - mean) * (record_J[i+k] - mean);
            lag_kappa += (record_kappa[i] - mean_k) * (record_kappa[i+k] - mean_k);
        }
        if( !k ) {
            var = lag;
            var_kappa = lag_kappa;
        }
        iact_J += 2*lag/var;
        iact_kappa += 2*lag_kappa/var_kappa;
'      if( k+2 > lag_size ) fprintf(stderr, "%d: %.2f %.2f %.2f %.2f\n", k, iact_J, iact_kappa, lag_kappa, mean_k);
    }
}' Compute_Auto_Correlation

void PrintTopologies(const sampler *smp, const char *where, boolean alawadhi, boolean exit_now) {
    int i;
    topology_gmodel_prior *tgp = NULL;
    boolean all_mono = true;
    dcpsampler *dcp = (dcpsampler *) smp->derived_smp;
    double ll = 0.0;
    partition_list *pl = alawadhi ? dcp->alawadhi_part_list : dcp->part_list;

    if( dcp->set->gmodel ) tgp = (topology_gmodel_prior *) dcp->tree_vec->top_prior;
    if( dcp->set->report_iact ) Compute_Auto_Correlation(dcp);

    if( debug || global_debug ) {
        fprintf(stderr, "%15s (%7d; %.2f[%.2f]; %.2f[%.2f]", where, dcp->smp->JumpNumber, (double) avg_ntopol/cnt, (double) rolling_ntopol_avg[rolling_index]/rolling_size, (double) avg_nparam/cnt, (double) rolling_nparam_avg[rolling_index]/rolling_size);
        if( dcp->set->gmodel ) fprintf(stderr, "; %.2f", dcp->mc);
        fprintf(stderr, "): ");
        for( i=0; i<pl->npartitions; i++ ) {
            ll += pl->part[i]->cPartialLogLikelihood;
            if(i) {
                fprintf(stderr, "[%4d:%c]", pl->part[i]->left, (pl->part[i]->topchange && pl->part[i]->parchange) ? 'B' : (pl->part[i]->topchange ? 'T' : 'P'));
                if( !pl->part[i]->topchange && pl->part[i-1]->ctree != pl->part[i]->ctree ) {
                    fprintf(stderr, "inconsistency %d NOT topology change point\n", pl->part[i]->ctree->tree_index);
                    exit(EXIT_FAILURE);
                }
                if( pl->part[i]->topchange && pl->part[i-1]->ctree == pl->part[i]->ctree ) {
                    fprintf(stderr, "inconsistency %d IS topology change point\n", pl->part[i]->ctree->tree_index);
                    exit(EXIT_FAILURE);
                }
            }
            fprintf(stderr, " %.2f %.2f %2d", pl->part[i]->cHyperParameter, pl->part[i]->cmatrix->v[0], pl->part[i]->ctree->tree_index);', pl->part[i]->cPartialLogLikelihood);
            if( dcp->set->gmodel ) {
                fprintf(stderr, "(%d)", tgp->monophyletic[pl->part[i]->ctree->tree_index]);
                if( ! tgp->monophyletic[pl->part[i]->ctree->tree_index] ) all_mono = false;
            }
            if(i) fprintf(stderr, "\t");
            else fprintf(stderr, " ");
        }
        if( dcp->set->gmodel ) fprintf(stderr, "%s", all_mono?"***":"");
        fprintf(stderr, " :|: %.2f %.2f", ll, (double)smp->acceptancerate[ADD_TWO_XI]/smp->tries[ADD_TWO_XI]);
        if( dcp->set->report_iact ) fprintf(stderr, " %.2f %.2f", iact_J, iact_kappa);
        fprintf(stderr, "\n");
    }
    if( exit_now ) {
        double tl = pl->top_lambda;
        double pla = pl->par_lambda;
        fprintf(stderr, "\nTopology change statistics : avg(%.4f == %.4f) std(%.4f == %.4f)", ((double) avg_ntopol/cnt), tl,
                (std_ntopol - avg_ntopol*avg_ntopol/cnt)/(double) (cnt-1), tl);
        fprintf(stderr, "\nParameter change statistics: avg(%.4f == %.4f) std(%.4f == %.4f)", ((double) avg_nparam/cnt), pla,
                (std_nparam - avg_nparam*avg_nparam/cnt)/(double) (cnt-1), pla);
        fprintf(stderr, "\nFraction monophyletic: %f\n", (double)mono_prob/cnt);
        fprintf(stderr, "\nTop Poisson: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f", ((double) top_poisson[0]/cnt),((double) top_poisson[1]/cnt),
                ((double) top_poisson[2]/cnt),((double) top_poisson[3]/cnt),((double) top_poisson[4]/cnt),((double) top_poisson[5]/cnt),
                ((double) top_poisson[6]/cnt),((double) top_poisson[7]/cnt),((double) top_poisson[8]/cnt),((double) top_poisson[9]/cnt));
        fprintf(stderr, "\nTop Poisson: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f", exp(-tl), exp(-tl)*tl, exp(-tl)*pow(tl,2)/2, exp(-tl)*pow(tl,3)/6, exp(-tl)*pow(tl,4)/24,
                exp(-tl)*pow(tl,5)/24/5, exp(-tl)*pow(tl,6)/24/30, exp(-tl)*pow(tl,7)/24/30/7, exp(-tl)*pow(tl,8)/24/30/56, exp(-tl)*pow(tl,9)/24/30/56/9);
        fprintf(stderr, "\nPar Poisson: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f", ((double) par_poisson[0]/cnt),((double) par_poisson[1]/cnt),
                ((double) par_poisson[2]/cnt),((double) par_poisson[3]/cnt),((double) par_poisson[4]/cnt),((double) par_poisson[5]/cnt),
                ((double) par_poisson[6]/cnt),((double) par_poisson[7]/cnt),((double) par_poisson[8]/cnt),((double) par_poisson[9]/cnt));
        fprintf(stderr, "\nPar Poisson: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f", exp(-pla), exp(-pla)*pla, exp(-pla)*pow(pla,2)/2, exp(-pla)*pow(pla,3)/6, exp(-pla)*pow(pla,4)/24,
                exp(-pla)*pow(pla,5)/24/5, exp(-pla)*pow(pla,6)/24/30, exp(-pla)*pow(pla,7)/24/30/7, exp(-pla)*pow(pla,8)/24/30/56, exp(-pla)*pow(pla,9)/24/30/56/9);
        fprintf(stderr, "\nTrees: %.4f %.4f %.4f", ((double) first_posn_tree[0]/cnt),((double) first_posn_tree[1]/cnt),((double) first_posn_tree[2]/cnt));
        fprintf(stderr, "\nMean log    mu: %.4f %.4f", mean_log_mu / cnt, (sd_log_mu - mean_log_mu*mean_log_mu/cnt)/(double) (cnt-1));
        fprintf(stderr, "\nMean log kappa: %.4f %.4f", mean_log_kappa / cnt, (sd_log_kappa - mean_log_kappa*mean_log_kappa/cnt)/(double) (cnt-1));
        fprintf(stderr, "\n");

        CloseSampler(smp);
    }

}' PrintTopologies

static void DCPReportState(sampler *smp, const char *location, double logRatio, double log_tau_ratio) {
'  const char *fxn_name = "DCPReportState";
    dcpsampler *dcp = (dcpsampler *) smp->derived_smp;
    
    PrintTopologies(smp, location, dcp->alawadhi_part_list != NULL, false);
    fprintf(stderr, "%.4f (%.4f)\n", logRatio, log_tau_ratio);
}' DCPReportState

static void DCPReportProposedState(sampler *smp, const char *where, tree **trees, int dim) {
    const char *fxn_name = "DCPReportProposedState";
    int i;
    dcpsampler *dcp = (dcpsampler *) smp->derived_smp;
    topology_gmodel_prior *tgp;

    if( smp->set->gmodel ) {
        tgp = (topology_gmodel_prior *) dcp->tree_vec->top_prior;
    }

    fprintf(stderr, "%s (%s): ", fxn_name, where);
    for( i=0; i<dim; i++ ) {
        fprintf(stderr, "[x] %d ", trees[i]->tree_index);
    }
    if( smp->set->gmodel ) fprintf(stderr, " (%.4f)", dcp->mc);
    fprintf(stderr, "\n");
}' DCPReportProposedState

static void DCPReportStatistics(sampler *smp, const char *transition_name, double priorRatio, double birthProb, double deathProb, double proposed_llike, double current_llike, int stat1, int stat2, int stat3) {
    const char *fxn_name = "ReportStatistics";
    dcpsampler *dcp = (dcpsampler *) smp->derived_smp;
    partition_list *pl = dcp->alawadhi_part_list ? dcp->alawadhi_part_list : dcp->part_list;
    topology_gmodel_prior *tgp = NULL;

    if( smp->set->gmodel ) tgp = (topology_gmodel_prior *) dcp->tree_vec->top_prior;
    fprintf(stderr, "%s(%d): ", fxn_name, (dcp->set->rng->useRnList?dcp->set->rng->current_rn:smp->JumpNumber));
    if( !strcmp(transition_name, "AddOne") ) {
        fprintf(stderr, " %s adding %d with state %d on %s: ap=%6.4e = (lpl=%6.4e - lcl=%6.4e) * pr=%6.4e * dp=%6.4e / bp=%6.4e (iap=%.4f)\n", transition_name, stat1, stat2,
                stat3?"right":"left",
                exp(priorRatio+deathProb-birthProb+proposed_llike-current_llike),
                proposed_llike, current_llike,
                exp(priorRatio), exp(deathProb), exp(birthProb), exp(-priorRatio-deathProb+birthProb));
    } else if( !strcmp(transition_name, "AddTwo") ) {
        fprintf(stderr, " %s adding (%d,%d) with state %d: ap=%6.4e = (lpl=%6.4e - lcl=%6.4e) * pr=%5.4e * dp=%6.4e / bp=%6.4e (iap=%.4f)\n", transition_name,
            stat1, stat2, stat3,
            exp(priorRatio+deathProb-birthProb),
            proposed_llike, current_llike,
            exp(priorRatio), exp(deathProb), exp(birthProb), exp(-priorRatio-deathProb+birthProb));
    } else if( !strcmp(transition_name, "DeleteOne") ) {
        fprintf(stderr, " %s deleting %d keeping %s: ap=%6.4e = (lpl=%6.4e - lcl=%6.4e) * pr=%6.4e * bp=%6.4e / dp=%6.4e (iap=%.4f)\n", transition_name,
                pl->part[stat1]->left, stat2?"left":"right",
                exp(priorRatio+deathProb-birthProb+proposed_llike-current_llike),
                proposed_llike, current_llike, exp(priorRatio), exp(birthProb), exp(deathProb),
                exp(-priorRatio-birthProb+deathProb));
    } else if( !strcmp(transition_name, "UpdateTopology") ) {
        fprintf(stderr, " %s updating %dth tree from %d", transition_name, stat1, stat2);
        if( smp->set->gmodel ) fprintf(stderr, "(%d)", tgp->monophyletic[stat2]);
        fprintf(stderr, " to %d", stat3);
        if( smp->set->gmodel ) fprintf(stderr, "(%d)", tgp->monophyletic[stat3]);
        fprintf(stderr, ": ap=%6.4e = (lpl=%6.4e - lcl=%6.4e) * pr=%6.4e * bp=%6.4e / dp=%6.4e (iap=%.4f)\n",
                exp(priorRatio+deathProb-birthProb+proposed_llike-current_llike),
                proposed_llike, current_llike,
                exp(priorRatio), exp(birthProb), exp(deathProb), exp(-priorRatio-birthProb+deathProb));
    } else if( !strcmp(transition_name, "UpdateXi") || !strcmp(transition_name, "UpdateRho") ) {
        fprintf(stderr, " %s moving %d to %d: ap=%6.4e = (lpl=%6.4e - lcl=%6.4e) * pr=%6.4e * bp=%6.4e / dp=%6.4e (iap=%.4f)\n", transition_name, stat1, stat2,
                exp(priorRatio+deathProb-birthProb+proposed_llike-current_llike),
                proposed_llike, current_llike,
                exp(priorRatio), exp(birthProb), exp(deathProb), exp(-priorRatio-birthProb+deathProb));
    } else {
        fprintf(stderr, " %s deleting (%d,%d): ap=%6.4e = (lpl=%6.4e - lcl=%6.4e) * pr=%6.4e * bp=%6.4e / dp=%6.4e (iap=%.4f)\n", transition_name,
                pl->part[stat1]->left, pl->part[stat2]->left,
                exp(priorRatio-deathProb+birthProb+proposed_llike-current_llike),
                    proposed_llike, current_llike,
                exp(priorRatio), exp(birthProb), exp(deathProb),
                exp(-priorRatio-birthProb+deathProb));
    }
}' DCPReportStatistics

boolean VerifyLikelihood(sampler *smp, boolean alawadhi) {
    dcpsampler *dcp = (dcpsampler *)smp->derived_smp;
    partition_list *pl = alawadhi ? dcp->alawadhi_part_list : dcp->part_list;
    int i;
    if( ! compute_likelihood ) return true;
    for( i=0; i<pl->npartitions; i++ ) {
        double ll1 = pl->part[i]->cPartialLogLikelihood;
        double ll2 = TreeLogLikelihood(pl->part[i]->ctree, smp, pl->part[i]->cmatrix, pl->part[i]->counts, pl->part[i]->cHyperParameter, false);
        if( fabs(ll1 - ll2) > tolerance ) {
            fprintf(stderr, "PROBLEM in VerifyLikelihood region %d (%d, %d) %f != %f (should be)\n", i, pl->part[i]->left, pl->part[i]->right, ll1, ll2);
            exit(EXIT_FAILURE);
        }
        if( ll1 > 0.0 ) {
            fprintf(stderr, "PROBLEM IN VerifyLikelihood region %d (%d, %d) has positive likelihood %f\n", i, pl->part[i]->left, pl->part[i]->right, ll1);
            exit(EXIT_FAILURE);
        }
    }
    return true;
}' VerifyLikelihood

boolean VerifyCounts(sampler *smp, const char *caller, boolean alawadhi) {
    dcpsampler *dcp = (dcpsampler *) smp->derived_smp;
    partition_list *pl = alawadhi ? dcp->alawadhi_part_list : dcp->part_list;
    int i, j, sum = 0, lsum;
    if( !alawadhi && dcp->smp->acceptancerate[ADD_RHO] - dcp->smp->acceptancerate[DELETE_RHO] != pl->parameter_changes ) {
        fprintf(stderr, "%s lead to acceptance rate error %d != %d\n", caller, dcp->smp->acceptancerate[ADD_RHO] - dcp->smp->acceptancerate[DELETE_RHO], pl->parameter_changes);
        exit(EXIT_FAILURE);
    }
    if( !alawadhi && dcp->smp->acceptancerate[ADD_XI] - dcp->smp->acceptancerate[DELETE_XI] != pl->topology_changes ) {
        fprintf(stderr, "%s lead to acceptance rate error %d != %d\n", caller, dcp->smp->acceptancerate[ADD_XI] - dcp->smp->acceptancerate[DELETE_XI], pl->topology_changes);
        exit(EXIT_FAILURE);
    }
    for( i=0; i<pl->npartitions; i++ ) {
        partition *part = pl->part[i];
        lsum = 0;
        if( fabs(part->cPartialLogHyperParameterPrior - dcp->smp->br->Log_Prior(dcp->smp->br, part->cHyperParameter)) > tolerance ) {
            fprintf(stderr, "%s lead to mismatch in prior on \\mu %f != %f\n", caller, part->cPartialLogHyperParameterPrior, dcp->smp->br->Log_Prior(dcp->smp->br, part->cHyperParameter));
            exit(EXIT_FAILURE);
        }
        for( j=0; j<part->lenunique; j++ ) {
            if( part->counts[j] < 0 ) {
                fprintf(stderr, "%s produced negative counts %d in region %d\n", caller, part->counts[j], i);
                exit(EXIT_FAILURE);
            }
            lsum += part->counts[j];
        }
        sum += lsum;
        if( lsum != part->right - part->left + 1 ) {
            fprintf(stderr, "%s produced inaccurate counts in region %d (%d != %d in (%d,%d))\n", caller, i, lsum, part->right - part->left + 1, part->left, part->right);
            exit(EXIT_FAILURE);
        }
    }
    if( sum != smp->sqd->lenseq ) {
        fprintf(stderr, "%s lost some counts (%d != %d)\n", caller, sum, smp->sqd->lenseq);
        exit(EXIT_FAILURE);
    }
    return true;
}' VerifyCounts

'discrete_gamma.c
#include "discrete_gamma.h"

    /* The following routines were copied from Yang (1994) */

/*-------------------------------------------------------------------------------
|                                                                               |
|  Discretization of gamma distribution with equal proportions in each          |
|  category.                                                                    |
|                                                                               |
-------------------------------------------------------------------------------*/
 int DiscreteGamma(double *rK, double alfa, double beta, int K, boolean median)
{
'System.err.println("DiscreteGamma("+alfa+", "+alfa+", "+K+")");
    int i;
    double  gap05 = gap05 = 1.0/(2.0*(double)K), t, factor = alfa/beta*K, lnga1;
    if (median)
    {
        for (i=0; i<K; i++)
            rK[i] = PointChi2(((double)i*2.0+1.0)*gap05, 2*alfa)/(2*beta);
        for (i=0,t=0; i<K; i++)
            t += rK[i];
        for (i=0; i<K; i++)
            rK[i] *= factor/t;
    }
    Else
    {
        double *freqK = (double *) malloc (sizeof(double) * K);
        lnga1 = LnGamma(alfa+1);
        for (i=0; i<K-1; i++)
        {
            freqK[i] = PointChi2(((double)i+1.0)/K, 2*alfa)/(2*beta);
        }
        for (i=0; i<K-1; i++)
            freqK[i] = IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
        rK[0] = freqK[0]*factor;
        rK[K-1] = (1-freqK[K-2])*factor;
        for (i=1; i<K-1; i++)
            rK[i] = (freqK[i]-freqK[i-1])*factor;
        free(freqK);
    }
    
    return (0);
}' DiscreteGamma

/*-------------------------------------------------------------------------------
|                                                                               |
|  Returns z so That Prob{x<z} = prob where x is Chi2 distributed with df=v.    |
|  Returns -1 if in error.   0.000002 < prob < 0.999998.                        |
|                                                                               |
|  RATNEST FORTRAN by                                                           |
|  Best, D. J. and D. E. Roberts.  1975.  The percentage points of the          |
|     Chi2 distribution.  Applied Statistics 24:385-388.  (AS91)                |
|                                                                               |
|  Converted into C by Ziheng Yang, Oct. 1993.                                  |
|                                                                               |
-------------------------------------------------------------------------------*/
double PointChi2 (double prob, double v)
{
    double  e = 0.5e-6, aa = 0.6931471805, p = prob, g,
                xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0,
                x = 0.0;

    if (p < 0.000002 || p > 0.999998 || v <= 0.0)
        return (-1.0);
    g = LnGamma (v/2.0);
    xx = v/2.0;
    c = xx - 1.0;
    if (v >= -1.24*log(p))
    {
        if (v > 0.32)
        {
            x = PointNormal (p);
            p1 = 0.222222/v;
            ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
            if (ch > 2.2*v+6.0)
                ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
            return( Getch(ch,xx,g,p,aa,c) );
        }
        Else
        {
            ch = 0.4;
            a = log(1.0-p);
            Do
            {
                q = ch;
                p1 = 1.0+ch*(4.67+ch);
                p2 = ch*(6.73+ch*(6.66+ch));
                t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
                ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
                if (fabs(q/ch-1.0)-0.01 <= 0.0)
                    return( Getch(ch,xx,g,p,aa,c) );
            }
            while(true);
        }
    }
    Else
    {
        ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
        if (ch-e<0)
            return (ch);
        return( Getch(ch,xx,g,p,aa,c) );
    }
}' PointChi2

double Getch(double ch, double xx, double g, double p, double aa, double c)
{
    double  q, p1, t, p2, b, a, s1, s2, s3, s4, s5, s6;
    double  e = 0.5e-6;
    Do
    {
        q = ch;
        p1 = 0.5*ch;
        if ((t = IncompleteGamma (p1, xx, g)) < 0.0)
        {
            fprintf(stderr, "\nerr IncompleteGamma\n");
            return (-1.0);
        }
        p2 = p-t;
        t = p2*exp(xx*aa+g+p1-c*log(ch));
        b = t/ch;
        a = 0.5*t-b*c;
        s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
        s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
        s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
        s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
        s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
        s6 = (120.0+c*(346.0+127.0*c))/5040.0;
        ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
    }
    while(fabs(q/ch-1.0) > e);
    return ch;
}' Getch

/*-------------------------------------------------------------------------------
|                                                                               |
|  Returns the incomplete gamma ratio I(x,alpha) where x is the upper           |
|  limit of the integration and alpha is the shape parameter.  Returns (-1)     |
|  if in error.                                                                 |
|  LnGamma_alpha = ln(Gamma(alpha)), is almost redundant.                      |
|  (1) series expansion     if (alpha>x || x<=1)                                |
|  (2) continued fraction   otherwise                                           |
|                                                                               |
|  RATNEST FORTRAN by                                                           |
|  Bhattacharjee, G. P.  1970.  The incomplete gamma integral.  Applied         |
|     Statistics, 19:285-287 (AS32)                                             |
|                                                                               |
-------------------------------------------------------------------------------*/
double IncompleteGamma (double x, double alpha, double LnGamma_alpha)
{
    int         i;
    double  p = alpha, g = LnGamma_alpha,
                accurate = 1e-8, overflow = 1e30,
                factor, gin = 0.0, rn = 0.0, a = 0.0, b = 0.0, an = 0.0,
                dif = 0.0, term = 0.0;
    double  pn[6];

    if (x == 0.0)
        return (0.0);
    if (x < 0 || p <= 0)
        return (-1.0);

    factor = exp(p*log(x)-x-g);
    if (x>1 && x>=p)
    {
        a = 1.0-p;
        b = a+x+1.0;
        term = 0.0;
        pn[0] = 1.0;
        pn[1] = x;
        pn[2] = x+1;
        pn[3] = x*b;
        gin = pn[2]/pn[3];
        Do
        {
            a++;
            b += 2.0;
            term++;
            an = a*term;
            for (i=0; i<2; i++)
                pn[i+4] = b*pn[i+2]-an*pn[i];
            if (pn[5] == 0)
            {
                for (i=0; i<4; i++)
                    pn[i] = pn[i+2];
                if (fabs(pn[4]) < overflow)
                    continue;
                for (i=0; i<4; i++)
                    pn[i] /= overflow;
                continue;
            }
            rn = pn[4]/pn[5];
            dif = fabs(gin-rn);
            if (dif>accurate)
            {
                gin = rn;
                for (i=0; i<4; i++)
                    pn[i] = pn[i+2];
                if (fabs(pn[4]) < overflow)
                    continue;
                for (i=0; i<4; i++)
                    pn[i] /= overflow;
                continue;
            }
            Else
            {
                gin = 1.0-factor*gin;
                return (gin);
            }
        }
        while(true);
    }
    Else
    {
        gin = 1.0;
        term = 1.0;
        rn = p;
        Do
        {
            rn++;
            term *= x/rn;
            gin += term;
        }
        while(term > accurate);
        gin *= factor/p;
        return (gin);
    }

}' IncompleteGamma

/*-------------------------------------------------------------------------------
|                                                                               |
|  Returns ln(gamma(alpha)) for alpha > 0, accurate to 10 decimal places.       |
|  Stirling's formula is used for the central polynomial part of the procedure. |
|                                                                               |
|  Pike, M. C. and I. D. Hill.  1966.  Algorithm 291: Logarithm of the gamma    |
|     function.  Communications of the Association for Computing                |
|     Machinery, 9:684.                                                         |
|                                                                               |
-------------------------------------------------------------------------------*/
double LnGamma (double alpha)
{
    double  x = alpha, f = 0.0, z;

    if (x < 7)
    {
        f = 1.0;
        z = x-1.0;
        While (Z < 7#)
            f *= z;
        x = z;
        f = -log(f);
    }
    z = 1.0/(x*x);
    'total =  f + (x-0.5)*log(x) - x + 0.918938533204673 + (((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z + 0.083333333333333)/x;
    
    'total = 0.789;
    'return (double) total;
    return  f + (x-0.5)*log(x) - x + 0.918938533204673 + (((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z + 0.083333333333333)/x;

}' LnGamma


/* ------------------------------------------------------------------------------
|                                                                               |
|  Returns z so That Prob{x<z} = prob where x ~ N(0,1) and                      |
|  (1e-12) < prob < 1-(1e-12).  Returns (-9999) if in error.                    |
|                                                                               |
|  Odeh, R. E. and J. O. Evans.  1974.  The percentage points of the normal     |
|     distribution.  Applied Statistics, 22:96-97 (AS70)                        |
|                                                                               |
|  Newer methods:                                                               |
|                                                                               |
|  Wichura, M. J.  1988.  Algorithm AS 241: The percentage points of the        |
|     normal distribution.  37:477-484.                                         |
|  Beasley, JD & S. G. Springer.  1977.  Algorithm AS 111: The percentage       |
|     points of the normal distribution.  26:118-121.                           |
|                                                                               |
-------------------------------------------------------------------------------*/
double PointNormal (double prob)
{

    double      a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
            a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
            b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
            y, z = 0.0, p = prob, p1;

    p1 = (p<0.5 ? p : 1-p);
    if (p1<1e-20)
       return (-9999);
    y = sqrt (log(1/(p1*p1)));
    z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
    return (p<0.5 ? -z : z);

}' PointNormal


' evol_param.c
#include "evol_param.h"

' FILE-WIDE VARIABLE DEFINITIONS
static int debug = 0;               ' Set to positive integer to turn on local debug output.
'static const char *file_name = "evol_param.c";
static const boolean fang_correction = false;   ' Use Fang's 0-value correction to acceptance probability calculation.

static void ParChgPtAddAccept(sampler *, partition_list *, partition *, partition *, int, int, int, int, qmatrix *, qmatrix *, double, double, double, double, double *, double, double);
/**
 * Sequentially updates evolutionary parameters at all partitions of the alignment
 */
    
void UpdateParameters(partition_list *pl, sampler *smp, boolean alawadhi) {
    const char *fxn_name = "UpdateParameters";
    double pLogHyperParameterPrior = 0;
    double pHyperParameter = 0; ' Proposed average branch length
    double cEP, pEP;        ' Current and updated evolutionary parameter values
    double pLogLikelihood = 0.0, cLogLikelihood = 0.0, logRatio;
        double *pPartialLogLikelihood = (double *) malloc(sizeof(double)*pl->npartitions);
    int i;
    int curr_index = 0;     ' Index of current region
    int prev_index = 0;     ' Index of previous region
    int next_index;         ' Index of next region
    qmatrix *pMatrix = NULL;    ' Proposed matrix
    partition *curr_part = NULL;
    boolean local_debug = false;

    ' Start with the first partition
    
    if( pPartialLogLikelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    do {
        curr_part = pl->part[curr_index];
    
        next_index = curr_index + 1;
        while( next_index < pl->npartitions && !pl->part[next_index]->parchange ) next_index++;

        ' Propose new evolutionary parameters for this region
        curr_part->cmatrix->Matrix_Proposer(&pMatrix, curr_part->cmatrix, smp->set);    ' Allocates memory
        pHyperParameter = smp->br->Propose(smp->set, curr_part->cHyperParameter, 1.0);
        pLogHyperParameterPrior = smp->br->Log_Prior(smp->br, pHyperParameter);
        cEP = curr_part->cmatrix->v[0];
        pEP = pMatrix->v[0];

        for( i=curr_index; i<next_index; i++ ) {    ' WAS BUG: was from prev_index !!
            partition *cpart = pl->part[i];

            if(compute_likelihood) pPartialLogLikelihood[i] = TreeLogLikelihood(cpart->ctree, smp, pMatrix, cpart->counts, pHyperParameter, false);
            else pPartialLogLikelihood[i] = 0.0;
            pLogLikelihood += pPartialLogLikelihood[i];
            if(compute_likelihood) cLogLikelihood += cpart->cPartialLogLikelihood;  ' WAS BUG: was cLogLikelihood =
        }

        logRatio = (alawadhi ? smp->set->alawadhi_factor : 1.0) *
            (
            pLogLikelihood -cLogLikelihood                      ' Likelihood ratio
            + pMatrix->log_prior - curr_part->cmatrix->log_prior            ' Prior on CTMC matrix
            + pLogHyperParameterPrior - curr_part->cPartialLogHyperParameterPrior   ' Prior on mu
            )
            + log(pHyperParameter) - log(curr_part->cHyperParameter)        ' Proposal on mu
            + log(pEP) - log(cEP)                           ' Proposal on kappa
            ;
        if( debug>1 || global_debug>1 || local_debug ) {
            fprintf(stderr, "U(%d): mu: %.4f -> %.4f; kappa: %.4f -> %.4f: ap=%6.4e = LR=%6.4e * pk=%6.4e / ck=%6.4e * pm=%6.4e / cm=%6.4e *%.2f/%.2f*%.2f/%.2f\n",
                (smp->set->rng->useRnList?smp->set->rng->current_rn:smp->JumpNumber),
                curr_part->cHyperParameter, pHyperParameter, cEP, pEP,
                exp(logRatio),
                exp(pLogLikelihood - cLogLikelihood),
                exp(pMatrix->log_prior), exp(curr_part->cmatrix->log_prior),
                exp(pLogHyperParameterPrior), exp(curr_part->cPartialLogHyperParameterPrior),
                pHyperParameter, curr_part->cHyperParameter, pEP, cEP);
        }
        if( !alawadhi ) smp->tries[UPDATE_KAPPA_AND_MU]++;
  
        ' Update all previous partitions since  last parameter change point

        if( logRatio > 0 || smp->set->rng->nextStandardUniform(smp->set->rng) < exp(logRatio) ) {

            if( !alawadhi ) smp->acceptancerate[UPDATE_KAPPA_AND_MU]++;
        
            ' Update all previous partitions since last parameter change point
            for( i=curr_index; i<next_index; i++ ) {
                partition *update_part = pl->part[i];
            
                update_part->cmatrix->Matrix_Copy(update_part->cmatrix, pMatrix);
                update_part->cHyperParameter = pHyperParameter;
                update_part->cPartialLogHyperParameterPrior = pLogHyperParameterPrior;
                update_part->cPartialLogLikelihood = pPartialLogLikelihood[i];
            }
            if( debug>3 || global_debug>3 || global_debug==-1 ) {
                VerifyLikelihood(smp, false);
                VerifyCounts(smp, "UpdateParameters", false);
            }
        }

        prev_index = curr_index;
        curr_index++;
        while( curr_index < pl->npartitions && !pl->part[curr_index]->parchange ) curr_index++;
        if( debug>0 || global_debug>0 ) smp->Report_State(smp, "UpdateParameters", logRatio, 0.0);
    } while( curr_index < pl->npartitions );

    
    QMatrixDelete(pMatrix);'pMatrix->Matrix_Delete(pMatrix);
    if( pPartialLogLikelihood ) free(pPartialLogLikelihood);
}' UpdateParameters

/**
 * Adds a new evolutionary change-point to the parameter space
 */
    
void ParChgPtAdd(partition_list *pl, sampler *smp) {
    const char *fxn_name = "ParChgPtAdd";
    int proposed;           ' Candidate new parameter change point
    int land_index;         ' Index of segment in which proposed falls
    int left_change_index;      ' Next partition over to left that has parameter changepoint on its left
    int right_change_index;     ' Next partition over to right that has parameter changepoint on its left
    int left_change_pos;        ' Left boundary of left_par_change partition
    int right_change_pos;       ' Right boundary of right_par_change partition
    int nth_part = 0;       ' The parameter partition to change is the nth
    partition *land_part;       ' Partition where proposed point lands
    partition *new_left;        ' Partition to hold new left region
    partition *new_right;       ' Partition to hold new right region
    partition *cpart;
    double pHyperParameter1, pHyperParameter2;  ' Proposed hyperparameters
    double pPartialLogHyperParameterPrior1, pPartialLogHyperParameterPrior2;
    qmatrix *pMatrix1, *pMatrix2;   ' New qmatrices
    double *pPartialLogLikelihood = (double *) malloc(sizeof(double)*pl->npartitions);
    double pLogLikelihood = 0.0, cLogLikelihood = 0.0, left_likelihood = 0.0, right_likelihood = 0.0;
    double logRatio = 0.0, birthProb, deathProb, priorRatio, logJacob;
    double matrix_log_proposal_prob, branch_log_proposal_prob;
    settings *set = smp->set;
    int i;
    partition_list *tpl = NULL; ' Just a pointer
    boolean local_debug = false;

    if( pPartialLogLikelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    ' Move x n-->m x' (from n- to m-dimensional space)
    proposed = pl->Propose_Change_Point(pl, set, false);

    ' Find partition where new point landed
    land_index = PartitionContaining(pl, proposed, 0, pl->npartitions - 1);
    land_part = pl->part[land_index];

    if( smp->set->alawadhi_parameter || smp->set->alawadhi ) {
        i = 1;
        while( i<=land_index ) if( pl->part[i++]->parchange ) nth_part++;
        nth_part++;
    }

    ' Find left and right parameter change point
    left_change_index = land_index;
    while( !pl->part[left_change_index]->parchange ) left_change_index--;
    right_change_index = land_index + 1;
    while( right_change_index < pl->npartitions && !pl->part[right_change_index]->parchange ) right_change_index++;

    ' Compute start and end of the constant parameter interval
    left_change_pos = pl->part[left_change_index]->left;
    right_change_pos = pl->part[right_change_index-1]->right + 1;

    ' Compute counts on both sides of proposed new change point
    if( land_part->left != proposed ) {
        PartitionMake(&new_left, smp->sqd->lenunique, land_part->left, proposed-1, land_part->topchange, land_part->parchange); ' Allocates memory
        PartitionCopySegmentCounts(new_left, smp->sqd, land_part->left, proposed);
        PartitionMake(&new_right, smp->sqd->lenunique, proposed, land_part->right, false, true);    ' Allocates memory
        PartitionCopyPartitionCountDifferences(new_right, land_part, new_left);
    } else {
        new_right = land_part;
    }
    
    ' Propose new evolutionary parameters (qmatrix should handle this by itself!)
    branch_log_proposal_prob = smp->br->Propose_Split(&pHyperParameter1, &pHyperParameter2, land_part->cHyperParameter, set, left_change_pos, proposed, right_change_pos);
    pPartialLogHyperParameterPrior1 = smp->br->Log_Prior(smp->br, pHyperParameter1);
    pPartialLogHyperParameterPrior2 = smp->br->Log_Prior(smp->br, pHyperParameter2);
    matrix_log_proposal_prob = land_part->cmatrix->Matrix_Propose_Split(&pMatrix1, &pMatrix2, land_part->cmatrix, set, left_change_pos, proposed, right_change_pos);
        
    ' Compute current and proposal likelihood on the left
    for( i = left_change_index; i < land_index; i++ ) {
        cpart = pl->part[i];
        if( compute_likelihood ) pPartialLogLikelihood[i] = TreeLogLikelihood(cpart->ctree, smp, pMatrix1, cpart->counts, pHyperParameter1, false);
        else pPartialLogLikelihood[i] = 0.0;

        pLogLikelihood += pPartialLogLikelihood[i];
        if( compute_likelihood ) cLogLikelihood += cpart->cPartialLogLikelihood;
    }

    ' and on the right
    for( i = land_index + 1; i < right_change_index; i++ ) {
        cpart = pl->part[i];
        if( compute_likelihood ) pPartialLogLikelihood[i] = TreeLogLikelihood(cpart->ctree, smp, pMatrix2, cpart->counts, pHyperParameter2, false);
        else pPartialLogLikelihood[i] = 0.0;

        pLogLikelihood += pPartialLogLikelihood[i];
        if( compute_likelihood ) cLogLikelihood += cpart->cPartialLogLikelihood;
    }

    ' Compute likelihood in landing region
    if( land_part->left != proposed && compute_likelihood )
        left_likelihood = TreeLogLikelihood(land_part->ctree, smp, pMatrix1, new_left->counts, pHyperParameter1, false);

    if( compute_likelihood ) right_likelihood = TreeLogLikelihood(land_part->ctree, smp, pMatrix2, new_right->counts, pHyperParameter2, false);

    pLogLikelihood += left_likelihood + right_likelihood;
    if( compute_likelihood ) cLogLikelihood += land_part->cPartialLogLikelihood;

    ' pi(x') / pi(x) =
    ' ln[ l(x') / l(x) ] = pLogLikelihood - cLogLikelihood
    ' + ln[ q(x') / q(x) ] = priorRatio
    priorRatio = pPartialLogHyperParameterPrior1 + pPartialLogHyperParameterPrior2 - land_part->cPartialLogHyperParameterPrior  ' \mu
        + pMatrix1->log_prior + pMatrix2->log_prior - land_part->cmatrix->log_prior                     ' \kappa
        ;

    ' ln[ q_{nm}(x',x) / q_{mn}(x,x') ] = deathProb - birthProb
    birthProb =
        + branch_log_proposal_prob                  ' \mu
        + matrix_log_proposal_prob                  ' \kappa
            ;
    deathProb = 0.0;                            ' cancels with birth/death probs

    logJacob = smp->logJacobian(smp, land_part->cmatrix, pMatrix1, pMatrix2, set->sigmaAlpha, land_part->cHyperParameter, pHyperParameter1, pHyperParameter2, set->sigmaMu);

    ' Compute acceptance ratio
    if( smp->set->alawadhi_parameter || smp->set->alawadhi ) {
        double ologRatio = pLogLikelihood - cLogLikelihood + priorRatio - birthProb + deathProb + logJacob;
        ' x m-->n x' -> ... -> x*: acceptance ratio = min{ 1, pi_n(x*) pi_n^*(x') q_{nm}(x',x) / pi_m(x) / pi_n^*(x*) / q_{mn}(x,x') }
        ' Convenient ratios to compute: pi_n(x*) / pi_n^*(x*) and pi_n^*(x') / pi_m(x)
        double cll = 0.0, pll = 0.0, ill = 0.0, clp, plp, ilp;

        ' Put the terms back in that cancel with the prior ratio (q(x')/q(x)) in the regular calculation
        deathProb += - log(pl->parameter_changes+1) + log(pl->par_dkp1);
        birthProb += log(pl->par_bk) - log(pl->alignment_length - 1 - pl->parameter_changes);
        logRatio = deathProb - birthProb + logJacob;

        ' q_{nm}(x',x) / q_{mn}(x,x') : The proposal ratio is unchanged, but was incomplete since probability of choice for new \rho^* canceled with prior

        ' pi(x)
        if( compute_likelihood ) for( i=0; i<pl->npartitions; i++ ) cll += pl->part[i]->cPartialLogLikelihood;
        clp = smp->Log_Prior(smp, false);

        ' x' -> ... -> x* (in tpl)
        smp->Alawadhi_Copy_State(smp, &tpl);        ' Sets up copies of pl on which we will do fixed dimension sampling in dimension n
        ParChgPtAddAccept(smp, tpl, new_left, new_right, proposed, left_change_index, land_index, right_change_index, pMatrix1, pMatrix2, pHyperParameter1,
                pHyperParameter2, pPartialLogHyperParameterPrior1, pPartialLogHyperParameterPrior2, pPartialLogLikelihood, left_likelihood, right_likelihood);

        ' pi(x')
        if( compute_likelihood ) for( i=0; i<tpl->npartitions; i++ ) ill += tpl->part[i]->cPartialLogLikelihood;
        ilp = smp->Log_Prior(smp, true);

        ' Fixed dimension sampler: x' -> x*
        smp->Fixed_Dimension_Sampler(smp, smp->set->alawadhi_k, KAPPA_MU|RHO);

        ' pi(x*)
        if( compute_likelihood ) for( i=0; i<tpl->npartitions; i++ ) pll += tpl->part[i]->cPartialLogLikelihood;
        plp = smp->Log_Prior(smp, true);

        logRatio += smp->set->alawadhi_factor*(ill + ilp - pll - plp) + pll + plp - cll - clp;
        if( smp->set->alawadhi_debug ) {
            fprintf(stderr, "ALAWADHI_DEBUG (%20s): logRatio: %e -> %e\n", smp->move_names[ADD_RHO], ologRatio, logRatio);
        }
    } else {
        logRatio = pLogLikelihood - cLogLikelihood + priorRatio - birthProb + deathProb + logJacob;
        if( fang_correction ) logRatio += log(pl->par_dkp1) - log(pl->par_bk) + pl->log_par_lambda - log(pl->parameter_changes+1);

        if( debug>1 || global_debug>1 || local_debug ) {
            fprintf(stderr, "A(%d): %d, %d, %d; mu: %.4f -> %.4f, %.4f; kappa: %.4f -> %.4f, %.4f : ap=%6.4e",
                smp->set->rng->useRnList?smp->set->rng->current_rn:smp->JumpNumber,
                left_change_pos, proposed, right_change_pos,
                land_part->cHyperParameter, pHyperParameter1, pHyperParameter2, 'sqrt(-2*(branch_log_proposal_prob - logOneOverSqrtTwoPi)),
                land_part->cmatrix->v[0], pMatrix1->v[0], pMatrix2->v[0], exp(logRatio)); 'sqrt(-2*(matrix_log_proposal_prob - logOneOverSqrtTwoPi)), exp(logRatio));
            if( debug>2 || global_debug>2 )
                fprintf(stderr, " = LR = (%f - %f) * pm1=%6.4e * pm2=%6.4e / cm=%6.4e * pa1=%6.4e * pa2=%6.4e / ca=%6.4e / zM=%6.4e / zA=%6.4e * j=%6.4e\n",
                    pLogLikelihood, cLogLikelihood,
                    exp(pPartialLogHyperParameterPrior1), exp(pPartialLogHyperParameterPrior2),
                    exp(land_part->cPartialLogHyperParameterPrior),
                    exp(pMatrix1->log_prior), exp(pMatrix2->log_prior),
                    exp(land_part->cmatrix->log_prior),
                    exp(branch_log_proposal_prob),
                    exp(matrix_log_proposal_prob),
                    exp(smp->logJacobian(smp, land_part->cmatrix, pMatrix1, pMatrix2, set->sigmaAlpha, land_part->cHyperParameter, pHyperParameter1, pHyperParameter2, set->sigmaMu)));
            else fprintf(stderr, "\n");
        }
    }

    smp->tries[ADD_RHO]++;
    if( logRatio > 0 || set->rng->nextStandardUniform(set->rng) < exp(logRatio) ) {

        smp->acceptancerate[ADD_RHO]++;

        if( smp->set->alawadhi_parameter || smp->set->alawadhi ) {
            smp->Alawadhi_Accept(smp);
        } else {
            ParChgPtAddAccept(smp, pl, new_left, new_right, proposed, left_change_index, land_index, right_change_index, pMatrix1, pMatrix2, pHyperParameter1,
                pHyperParameter2, pPartialLogHyperParameterPrior1, pPartialLogHyperParameterPrior2, pPartialLogLikelihood, left_likelihood, right_likelihood);
        }

        if( debug>0 || global_debug>0 || local_debug>0 ) smp->Report_State(smp, "AddParChange", logRatio, 0.0);
        if( debug>3 || global_debug>3 || global_debug==-1 ) {
            VerifyLikelihood(smp, false);
            VerifyCounts(smp, "AddParChange", false);
        }
    } else {
        if( smp->set->alawadhi_parameter || smp->set->alawadhi ) smp->Alawadhi_Reject(smp);
        else {
            if( land_part->left != proposed ) {
                PartitionDelete(new_left);
                PartitionDelete(new_right);
            }
            QMatrixDelete(pMatrix1);'pMatrix1->Matrix_Delete(pMatrix1);
            QMatrixDelete(pMatrix2);'pMatrix2->Matrix_Delete(pMatrix2);
        }
    }
    if(pPartialLogLikelihood) free(pPartialLogLikelihood);
}' ParChgPtAdd

/**
 * Processes the addition of a parameter change point
 */

static void ParChgPtAddAccept(sampler *smp, partition_list *pl, partition *new_left, partition *new_right, int proposed, int left_change_index, int land_index, int right_change_index, qmatrix *pMatrix1, qmatrix *pMatrix2, double pHyperParameter1, double pHyperParameter2, double pPartialLogHyperParameterPrior1, double pPartialLogHyperParameterPrior2, double *pPartialLogLikelihood, double left_likelihood, double right_likelihood) {
    int i;
    partition *land_part = pl->part[land_index];

    ' Update intermediate likelihoods and evolutionary parameters
    for( i = left_change_index; i < land_index; i++ ) {
        partition *cpart = pl->part[i];

        cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
        cpart->cmatrix->Matrix_Copy(cpart->cmatrix, pMatrix1);
        cpart->cHyperParameter = pHyperParameter1;
        cpart->cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior1;
    }
    
    for( i = land_index + 1; i < right_change_index; i++ ) {
        partition *cpart = pl->part[i];

        cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
        cpart->cmatrix->Matrix_Copy(cpart->cmatrix, pMatrix2);
        cpart->cHyperParameter = pHyperParameter2;
        cpart->cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior2;
    }
        
    ' Insert new parameter change point

    ' This segment already exists as a topology change point
    if( land_part->left == proposed ) {

        land_part->parchange = true;
        land_part->cPartialLogLikelihood = right_likelihood;
        ' Delete the old matrix associated with this partition
        QMatrixDelete(land_part->cmatrix);'PartitionDeleteMatrix(land_part);
        ' And point to the new matrix instead
        land_part->cmatrix = pMatrix2;
        land_part->cHyperParameter = pHyperParameter2;
        land_part->cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior2;

        QMatrixDelete(pMatrix1);
    }
    else{
        ' Make landing partition into the new left partition
        memcpy(land_part->counts, new_left->counts, sizeof(int)*smp->sqd->lenunique);
        land_part->cPartialLogLikelihood = left_likelihood;

        ' Delete the old matrix associated with this partition
        QMatrixDelete(land_part->cmatrix);

        ' And point to the new left matrix instead
        land_part->cmatrix = pMatrix1;
        land_part->cHyperParameter = pHyperParameter1;
        land_part->cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior1;
        land_part->right = new_left->right;

        ' We don't need new_left anymore.  Clear the memory
        PartitionDelete(new_left);

        ' Update right side
        new_right->cPartialLogLikelihood = right_likelihood;
        ' Give the new partition a matrix
        new_right->cmatrix = pMatrix2;
        new_right->cHyperParameter = pHyperParameter2;
        new_right->cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior2;
        new_right->ctree = land_part->ctree;    ' tree unchanged

        PartitionListAddPartition(pl, new_right, land_index + 1);
    }
    pl->parameter_changes++;
}' ParChgPtAddAccept

/**
 * Removes an evolutionary change-point from the parameter space
 */

void ParChgPtDelete(partition_list *pl, sampler *smp) {
    const char *fxn_name = "ParChgPtDelete";
    int propose_to_delete_index;        ' Index of partition with partition change point that is proposed to delete
    int left_change_index;          ' Index of leftward partition with parameter change point
    int right_change_index;         ' Index of rightward partition with parameter change point
    int left_change_pos;            ' Beginning of constant parameter interval
    int right_change_pos;           ' End of constant parameter interval
    int nth_parameter_segment;      ' Which parameter interval (right half) is proposed for deletion
    partition *propose_to_delete_part;  ' Partition proposed to delete
    partition *left_part;           ' Partition just left of the change point proposed to delete
    partition *cpart;
    double cHyperParameter1, cHyperParameter2;
    double pHyperParameter;
    double cPartialLogHyperParameterPrior1, cPartialLogHyperParameterPrior2, pPartialLogHyperParameterPrior;
    qmatrix *pMatrix = NULL;        ' Proposed matrix for new merged segments
    double *pPartialLogLikelihood = NULL;
    double pLogLikelihood = 0.0, cLogLikelihood = 0.0, logRatio, birthProb, deathProb, priorRatio;
    double matrix_log_proposal_prob, branch_log_proposal_prob;
    partition_list *tpl = NULL;
    double cll=0.0, clp=0.0, ill=0.0, ilp=0.0;
    boolean local_debug = false;
    int i;

    ' Before we propose to move to another dimension: do some fix dimension moves
    if( smp->set->alawadhi_parameter || smp->set->alawadhi ) {
        ' Move: x* -> ... -> x' (in tpl)
        for( i=0; i<pl->npartitions; i++ ) cll += pl->part[i]->cPartialLogLikelihood;   ' log l(x*)
        clp = smp->Log_Prior(smp, false);   ' log q(x*)
        smp->Alawadhi_Copy_State(smp, &tpl);
        smp->Fixed_Dimension_Sampler(smp, smp->set->alawadhi_k, KAPPA_MU|RHO);
        for( i=0; i<tpl->npartitions; i++ ) ill += tpl->part[i]->cPartialLogLikelihood; ' log l(x')
        ilp = smp->Log_Prior(smp, true);    ' log q(x')
    } else {
        tpl = pl;
    }

    pPartialLogLikelihood = (double *) malloc(sizeof(double)*tpl->npartitions);
    if( pPartialLogLikelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    
    ' Move: x' n-->m x (in tpl)
    propose_to_delete_index = tpl->Propose_Parameter_Change_Point_To_Delete(tpl, smp->set);
    nth_parameter_segment = 0;
    i = 1;
    while( i<propose_to_delete_index ) if( tpl->part[i++]->parchange ) nth_parameter_segment++;

    ' Find left and right parameter change point
    left_change_index = propose_to_delete_index - 1;
    while( !tpl->part[left_change_index]->parchange ) left_change_index--;
    right_change_index = propose_to_delete_index + 1;
    while( right_change_index < tpl->npartitions && !tpl->part[right_change_index]->parchange ) right_change_index++;

    ' Compute start and end of the constant parameter interval
    left_change_pos = tpl->part[left_change_index]->left;
    right_change_pos = smp->sqd->lenseq;
    if( right_change_index < tpl->npartitions ) right_change_pos = tpl->part[right_change_index]->left;

    propose_to_delete_part = tpl->part[propose_to_delete_index];
    left_part = tpl->part[propose_to_delete_index - 1];

    ' Propose new evolutionary parameters
    cHyperParameter1 = left_part->cHyperParameter;
    cHyperParameter2 = propose_to_delete_part->cHyperParameter;
    cPartialLogHyperParameterPrior1 = left_part->cPartialLogHyperParameterPrior;
    cPartialLogHyperParameterPrior2 = propose_to_delete_part->cPartialLogHyperParameterPrior;

    branch_log_proposal_prob = smp->br->Propose_Merge(&pHyperParameter, cHyperParameter1, cHyperParameter2, smp->set, left_change_pos, propose_to_delete_part->left, right_change_pos);
    pPartialLogHyperParameterPrior = smp->br->Log_Prior(smp->br, pHyperParameter);
    
    matrix_log_proposal_prob = left_part->cmatrix->Matrix_Propose_Merge(&pMatrix, left_part->cmatrix, propose_to_delete_part->cmatrix, smp->set, left_change_pos, propose_to_delete_part->left, right_change_pos);

    ' Compute current and proposal likelihood on the left and on the right
    for( i = left_change_index; i < right_change_index; i++ ) {
        cpart = tpl->part[i];
        if( compute_likelihood) pPartialLogLikelihood[i] = TreeLogLikelihood(cpart->ctree, smp, pMatrix, cpart->counts, pHyperParameter, false);
        else pPartialLogLikelihood[i] = 0.0;
        
        pLogLikelihood += pPartialLogLikelihood[i];
        if( compute_likelihood ) cLogLikelihood += cpart->cPartialLogLikelihood;
    }

    ' q_{mn}(x,x') / q_{nm}(x',x)
    birthProb = matrix_log_proposal_prob + branch_log_proposal_prob;
    deathProb = 0.0;

    ' Prior ratio: q(x) / q(x')
    priorRatio = pMatrix->log_prior
        - left_part->cmatrix->log_prior - propose_to_delete_part->cmatrix->log_prior
        + pPartialLogHyperParameterPrior
        - cPartialLogHyperParameterPrior1 - cPartialLogHyperParameterPrior2
        ;

    logRatio = birthProb - deathProb + priorRatio + pLogLikelihood - cLogLikelihood
        - smp->logJacobian(smp, pMatrix, left_part->cmatrix, propose_to_delete_part->cmatrix, smp->set->sigmaAlpha, pHyperParameter, cHyperParameter1, cHyperParameter2, smp->set->sigmaMu);
    if( fang_correction) logRatio += log(tpl->par_bkm1) - log(tpl->par_dk) - tpl->log_par_lambda + log(tpl->parameter_changes);

    ' Compute acceptance ratio
    if( smp->set->alawadhi_parameter || smp->set->alawadhi ) {
        ' pl          tpl
        ' x* -> ... -> x' -> x
        ' l*(x*) / l*(x')
        double ologRatio = logRatio;

        ' Already includes factor:  pll + plp - ill - ilp
        ' pi(x') * pi*(x*) / pi(x*) / pi*(x')
        logRatio += smp->set->alawadhi_factor*(cll + clp - ill - ilp) + ill + ilp - cll - clp;
        
        if( smp->set->alawadhi_debug ) {
            fprintf(stderr, "ALAWADHI_DEBUG (%20s): logRatio: %e -> %e\n", smp->move_names[DELETE_RHO], ologRatio, logRatio);
        }
    } else {

        if( debug>1 || global_debug>1 || local_debug ) {
            fprintf(stderr, "D(%d): %d, %d, %d; mu: %.4f %.4f -> %.4f; kappa: %.4f %.4f -> %.4f: ap=%.4f",
                smp->set->rng->useRnList?smp->set->rng->current_rn:smp->JumpNumber,
                left_change_pos, propose_to_delete_part->left, right_change_pos,
                left_part->cHyperParameter, propose_to_delete_part->cHyperParameter, pHyperParameter,
                left_part->cmatrix->v[0], propose_to_delete_part->cmatrix->v[0], pMatrix->v[0],
                exp(logRatio));
            if( debug>2 || global_debug>2 )
                fprintf(stderr, " = LR = %f * pk=%6.4e / ck1=%6.4e / ck2=%6.4e * pm=%6.4e / cm1=%6.4e / cm2=%6.4e * zA=%6.4e * zM=%6.4e / j=%6.4e %f %f\n",
                    exp(pLogLikelihood - cLogLikelihood),
                    exp(pMatrix->log_prior),
                    exp(left_part->cmatrix->log_prior), exp(propose_to_delete_part->cmatrix->log_prior),
                    exp(pPartialLogHyperParameterPrior),
                    exp(cPartialLogHyperParameterPrior1), exp(cPartialLogHyperParameterPrior2),
                    exp(matrix_log_proposal_prob), exp(branch_log_proposal_prob),
                    exp(smp->logJacobian(smp, pMatrix, left_part->cmatrix, propose_to_delete_part->cmatrix, smp->set->sigmaAlpha, pHyperParameter, cHyperParameter1, cHyperParameter2, smp->set->sigmaMu)),
                    pLogLikelihood, cLogLikelihood);
            else fprintf(stderr, "\n");
        }
    }

    smp->tries[DELETE_RHO]++;

    if( logRatio > 0 || smp->set->rng->nextStandardUniform(smp->set->rng) < exp(logRatio) ) {
        
        smp->acceptancerate[DELETE_RHO]++;

        ' Update intermediate likelihoods and evolutionary parameters
        for( i = left_change_index; i < right_change_index; i++ ) {
            if( i == propose_to_delete_index ) continue;
            cpart = tpl->part[i];
        
            cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
            cpart->cmatrix->Matrix_Copy(cpart->cmatrix, pMatrix);
            cpart->cHyperParameter = pHyperParameter;
            cpart->cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior;
        }
        
        ' Can't remove the change point because it is also a topology change point
        if( propose_to_delete_part->topchange ) {
            propose_to_delete_part->parchange = false;
            propose_to_delete_part->cPartialLogLikelihood = pPartialLogLikelihood[propose_to_delete_index];
            propose_to_delete_part->cmatrix->Matrix_Copy(propose_to_delete_part->cmatrix, pMatrix);
            propose_to_delete_part->cHyperParameter = pHyperParameter;
            propose_to_delete_part->cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior;
        }
        ' Remove parameter change point since it is not also topology change point
        else{
            PartitionAddPartition(left_part, propose_to_delete_part);
            left_part->cPartialLogLikelihood += pPartialLogLikelihood[propose_to_delete_index];
            left_part->right = propose_to_delete_part->right;
        
            PartitionListRemovePartition(tpl, propose_to_delete_index);
        }

        tpl->parameter_changes--;
        if( smp->set->alawadhi_parameter || smp->set->alawadhi ) smp->Alawadhi_Accept(smp);
        if( debug>0 || global_debug>0 || local_debug>0 ) smp->Report_State(smp, "DeleteParChange", logRatio, 0.0);
        if( debug>3 || global_debug>3 || global_debug==-1 ) {
            VerifyLikelihood(smp, false);
            VerifyCounts(smp, "DeleteParChange", false);
        }
    } else {
        if( smp->set->alawadhi_parameter || smp->set->alawadhi ) smp->Alawadhi_Reject(smp);
    }
    QMatrixDelete(pMatrix);
    if(pPartialLogLikelihood) free(pPartialLogLikelihood);
}' ParChgPtDelete

' ihkynoboundfixpimatrix.c

#include "ihkynoboundfixpimatrix.h"

' TODO: review for appropriate positioning of calls to Matrix_Update_Cache

' Global variables:
extern void *qmatrix_prior; ' Shared prior across instances of qmatrices

' File local variables:
'static int debug = 0;                                 ' Set to positive integer to turn on local debugging output
'static const char *file_name = "ihkynoboundfixpimatrix.c";

' File local function predeclarations:
static void MakePrior(settings *);                          ' Setup the global prior structure
static void Default_iHKYNoBoundFixPiMatrix(qmatrix **);                 ' Makes default HKY matrix of 1st arg qmatrix
static void iHKYNoBoundFixPiMatrixReset(qmatrix *, const double *, const double *); ' Resets 1st arg HKY matrix using 2nd arg vector of parameter values, and 3rd arg stationary distn
static void UpdateVariables(qmatrix *);                         ' Internal function for updating internal state of matrix object upon change
static void iHKYNoBoundFixPiMatrixMakeCopy(qmatrix **, const qmatrix *);        ' Make 1st arg qmatrix by copying 2nd arg qmatrix
static void iHKYNoBoundFixPiMatrixCopy(qmatrix *, const qmatrix *);         ' Copy 2nd arg qmatrix into 2nd arg qmatrix (assumes both args pre-exist)
static void iHKYNoBoundFixPiMatrixUpdateCache(qmatrix *, ...);              ' Update 1st arg qmatrix cache
static double iHKYNoBoundFixPiMatrixProposer(qmatrix **, const qmatrix *, settings *);  ' See qmatrix.h Matrix_Proposer function pointer
static double iHKYNoBoundFixPiMatrixLogPrior(qmatrix *, ...);

' Split alignment segment into two contiguous, adjacent segments (See qmatrix.h Matrix_Propose_Split function pointer)
static double iHKYNoBoundFixPiMatrixProposeSplit(qmatrix **, qmatrix **, const qmatrix *, settings *, ...);

' Merge two physically adjacent alignments segments (See qmatrix.h Matrix_Propose_Merge function pointer)
static double iHKYNoBoundFixPiMatrixProposeMerge(qmatrix **, const qmatrix *, const qmatrix *, settings *, ...);

static void Default_iHKYNoBoundFixPiMatrix(qmatrix **qmt) {
    const char *fxn_name = "Default_iHKYNoBoundFixPiMatrix";
    ihkynoboundfixpimatrix *ihky = NULL;

    QMatrixMake(qmt, 4, 1);
    (*qmt)->v[0] = 1.0/3.0;
    (*qmt)->Matrix_Update_Cache = &iHKYNoBoundFixPiMatrixUpdateCache;
    (*qmt)->Matrix_Copy = &iHKYNoBoundFixPiMatrixCopy;
    (*qmt)->Matrix_Make_Copy = &iHKYNoBoundFixPiMatrixMakeCopy;
    (*qmt)->Matrix_Delete = NULL;   ' No allocations in ihky except the derived_mt itself so no special function needed
    (*qmt)->Matrix_Proposer = &iHKYNoBoundFixPiMatrixProposer;
    (*qmt)->Matrix_Propose_Split = &iHKYNoBoundFixPiMatrixProposeSplit;
    (*qmt)->Matrix_Propose_Merge = &iHKYNoBoundFixPiMatrixProposeMerge;
    (*qmt)->Matrix_Sync = &UpdateVariables;

    ' Setup prior
    (*qmt)->Matrix_Log_Prior = &iHKYNoBoundFixPiMatrixLogPrior;

    ihky = (ihkynoboundfixpimatrix *) malloc(sizeof(ihkynoboundfixpimatrix));   ' BUG: this change creates bugs elsewhere
    if( ihky == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    (*qmt)->derived_mt = (void *)ihky;
    ihky->cached_avg_brlen = 0;
    ihky->cached_ep = -1;
    ihky->model = ALPHA;
    UpdateVariables(*qmt);
}' Default_iHKYNoBoundFixPiMatrix

void iHKYNoBoundFixPiMatrixGlobalInitialize(settings *set) {
    MakePrior(set);
}' iHKYNoBoundFixPiMatrixGlobalInitialize

void iHKYNoBoundFixPiMatrixMakeDefault(qmatrix **qmt, int model) {
'  const char *fxn_name = "iHKYNoBoundFixPiMatrixMake1";
    ihkynoboundfixpimatrix *ihky = NULL;

    Default_iHKYNoBoundFixPiMatrix(qmt);
    ihky = (ihkynoboundfixpimatrix *) (*qmt)->derived_mt;
    ihky->model = model;
    (*qmt)->v[0] = 0.5;
    (*qmt)->pi[0] = (*qmt)->pi[1] = (*qmt)->pi[2] = (*qmt)->pi[3] = 0.25;
    UpdateVariables(*qmt);
}' iHKYNoBoundFixPiMatrixMakeDefault

void iHKYNoBoundFixPiMatrixMakeInitial(qmatrix **qmt, int model, const seqdata *sqd, settings *set) {
'  const char *fxn_name = "iHKYNoBoundFixPiMatrixMake2";
    ihkynoboundfixpimatrix *ihky = NULL;

    Default_iHKYNoBoundFixPiMatrix(qmt);
    ihky = (ihkynoboundfixpimatrix *) (*qmt)->derived_mt;
    ihky->model = model;
    if( model == ALPHA) (*qmt)->v[0] = set->rng->nextStandardUniform(set->rng);
    else if( model == KAPPA ) (*qmt)->v[0] = set->rng->nextDispersedUniform(set->rng);
    Alignment_Composition(sqd, &((*qmt)->pi));
    UpdateVariables(*qmt);
}' iHKYNoBoundFixPiMatrixMakeInitial

void iHKYNoBoundFixPiMatrixMakeAndSet(qmatrix **qmt, int model, const double *inV, const double *inPi) {
'  const char *fxn_name = "iHKYNoBoundFixPiMatrixMakeAndSet";
    ihkynoboundfixpimatrix *ihky = NULL;

    Default_iHKYNoBoundFixPiMatrix(qmt);
    ihky = (ihkynoboundfixpimatrix *) (*qmt)->derived_mt;
    ihky->model = model;
    iHKYNoBoundFixPiMatrixReset(*qmt, inV, inPi);
}' iHKYNoBoundFixPiMatrixMakeAndSet

static void iHKYNoBoundFixPiMatrixMakeCopy(qmatrix **new_qmt, const qmatrix *old_qmt) {
    ihkynoboundfixpimatrix *old_ihky = (ihkynoboundfixpimatrix *) old_qmt->derived_mt;
    iHKYNoBoundFixPiMatrixMakeAndSet(new_qmt, old_ihky->model, old_qmt->v, old_qmt->pi);
}' iHKYNoBoundFixPiMatrixMakeCopy

static void iHKYNoBoundFixPiMatrixCopy(qmatrix *update_qmt, const qmatrix *copied_qmt) {
    ihkynoboundfixpimatrix *update_hky = (ihkynoboundfixpimatrix *) update_qmt->derived_mt;
    ihkynoboundfixpimatrix *copied_hky = (ihkynoboundfixpimatrix *) copied_qmt->derived_mt;
    QMatrixCopy(update_qmt, copied_qmt);
    update_hky->cached_ep = copied_hky->cached_ep;
    update_hky->cached_avg_brlen = copied_hky->cached_avg_brlen;
    update_hky->alpha = copied_hky->alpha;
    update_hky->kappa = copied_hky->kappa;
    update_hky->model = copied_hky->model;
    update_hky->pR = copied_hky->pR;
    update_hky->pY = copied_hky->pY;
    update_hky->b = copied_hky->b;
    update_hky->c = copied_hky->c;
    update_hky->d = copied_hky->d;
}' iHKYNoBoundFixPiMatrixCopy

static void iHKYNoBoundFixPiMatrixReset(qmatrix *qmt, const double *inV, const double *inPi) {
    ihkynoboundfixpimatrix *ihky = (ihkynoboundfixpimatrix *) qmt->derived_mt;
    QMatrixUpdateParameters(qmt, inV, inPi);
    UpdateVariables(qmt);
    ihky->cached_avg_brlen = 0;
    ihky->cached_ep = -1;
}' iHKYNoBoundFixPiMatrixReset

static void UpdateVariables(qmatrix *qmt) {
'  const char *fxn_name = "UpdateVariables";
    ihkynoboundfixpimatrix *ihky = (ihkynoboundfixpimatrix *)qmt->derived_mt;

    if(ihky->model == KAPPA) {
        ihky->kappa = qmt->v[0];
        qmt->log_prior = qmt->Matrix_Log_Prior(qmt, ihky->kappa);
    } else if (ihky->model == ALPHA) {
        ihky->alpha = qmt->v[0];
'      qmt->log_prior = qmt->Matrix_Log_Prior(qmt, ihky->alpha);
    }
    ihky->pR = qmt->pi[0] + qmt->pi[1];
    ihky->pY = qmt->pi[2] + qmt->pi[3];
    if( ihky->model == ALPHA ) {
        ihky->b = (1-ihky->alpha)/2;
        ihky->c = ihky->alpha*ihky->pR + ihky->b*ihky->pY;
        ihky->d = ihky->b*ihky->pR + ihky->alpha*ihky->pY;
    } else if( ihky->model == KAPPA ) {
        ihky->b = 0.5 / (ihky->kappa*qmt->pi[0]*qmt->pi[1] + ihky->kappa*qmt->pi[2]*qmt->pi[3] + ihky->pR*ihky->pY);
        ihky->c = (ihky->pY + ihky->kappa*ihky->pR)*ihky->b;
        ihky->d = (ihky->pR + ihky->kappa*ihky->pY)*ihky->b;
    }
}' UpdateVariables

void Set_Alpha(qmatrix *qmt, double alpha) {
    ihkynoboundfixpimatrix *ihky = (ihkynoboundfixpimatrix *)qmt->derived_mt;
    ihky->model = ALPHA;
    qmt->v[0] = alpha;
    UpdateVariables(qmt);
    qmt->Matrix_Update_Cache(qmt, ihky->cached_avg_brlen);
}' Set_Alpha

void Set_Kappa(qmatrix *qmt, double kappa) {
    ihkynoboundfixpimatrix *ihky = (ihkynoboundfixpimatrix *) qmt->derived_mt;
    ihky->model = KAPPA;
    qmt->v[0] = kappa;
    UpdateVariables(qmt);
    qmt->Matrix_Update_Cache(qmt, ihky->cached_avg_brlen);
}' Set_Kappa

' Assumption: pi is not changing
static void iHKYNoBoundFixPiMatrixUpdateCache(qmatrix *qmt, ... ) {' double t) {
'  const char *fxn_name = "iHKYNoBoundFixPiMatrixUpdateCache";
    ihkynoboundfixpimatrix *ihky = (ihkynoboundfixpimatrix *) qmt->derived_mt;
    double *pi = qmt->pi;
    va_list vargs;
    double t;
    int i, j;

    ' Retrieve arguments
    va_start(vargs, qmt);
    t = va_arg(vargs, double);
    va_end(vargs);

    if( t != ihky->cached_avg_brlen || qmt->v[0] != ihky->cached_ep ) {
        double pY, pR, bt, ct, dt;
        ihky->cached_avg_brlen = t;
        pY = ihky->pY;
        pR = ihky->pR;
        bt = 1 + ihky->b*t;
        ct = 1 + ihky->c*t;
        dt = 1 + ihky->d*t;

' Have transposed cached_qmatrix to speed up likelihood calculation
        if( ihky->model == KAPPA ) {
            ihky->cached_ep = ihky->kappa;
            qmt->cached_qmatrix[0][0] = pi[0] + pi[0]*pY/pR/bt + pi[1]/pR/ct;
            qmt->cached_qmatrix[1][0] = pi[1] + pi[1]*pY/pR/bt - pi[1]/pR/ct;
            qmt->cached_qmatrix[2][0] = pi[2] - pi[2]/bt;
            qmt->cached_qmatrix[3][0] = pi[3] - pi[3]/bt;
            qmt->cached_qmatrix[0][1] = pi[0] + pi[0]*pY/pR/bt - pi[0]/pR/ct;
            qmt->cached_qmatrix[1][1] = pi[1] + pi[1]*pY/pR/bt + pi[0]/pR/ct;
            qmt->cached_qmatrix[2][1] = pi[2] - pi[2]/bt;
            qmt->cached_qmatrix[3][1] = pi[3] - pi[3]/bt;
            qmt->cached_qmatrix[0][2] = pi[0] - pi[0]/bt;
            qmt->cached_qmatrix[1][2] = pi[1] - pi[1]/bt;
            qmt->cached_qmatrix[2][2] = pi[2] + pi[2]*pR/pY/bt + pi[3]/pY/dt;
            qmt->cached_qmatrix[3][2] = pi[3] + pi[3]*pR/pY/bt - pi[3]/pY/dt;
            qmt->cached_qmatrix[0][3] = pi[0] - pi[0]/bt;
            qmt->cached_qmatrix[1][3] = pi[1] - pi[1]/bt;
            qmt->cached_qmatrix[2][3] = pi[2] + pi[2]*pR/pY/bt - pi[2]/pY/dt;
            qmt->cached_qmatrix[3][3] = pi[3] + pi[3]*pR/pY/bt + pi[2]/pY/dt;
        } else if( ihky->model == ALPHA ) {
            ihky->cached_ep = ihky->alpha;
                    qmt->cached_qmatrix[0][0] = pi[0] + pi[0]*pY/pR/bt + pi[1]/pR/ct;
                    qmt->cached_qmatrix[1][0] = pi[1] + pi[1]*pY/pR/bt - pi[1]/pR/ct;
                    qmt->cached_qmatrix[2][0] = pi[2] - pi[2]/bt;
                    qmt->cached_qmatrix[3][0] = pi[3] - pi[3]/bt;
                    qmt->cached_qmatrix[0][1] = pi[0] + pi[0]*pY/pR/bt - pi[0]/pR/ct;
                    qmt->cached_qmatrix[1][1] = pi[1] + pi[1]*pY/pR/bt + pi[0]/pR/ct;
                    qmt->cached_qmatrix[2][1] = pi[2] - pi[2]/bt;
                    qmt->cached_qmatrix[3][1] = pi[3] - pi[3]/bt;
                    qmt->cached_qmatrix[0][2] = pi[0] - pi[0]/bt;
                    qmt->cached_qmatrix[1][2] = pi[1] - pi[1]/bt;
                    qmt->cached_qmatrix[2][2] = pi[2] + pi[2]*pR/pY/bt + pi[3]/pY/dt;
                    qmt->cached_qmatrix[3][2] = pi[3] + pi[3]*pR/pY/bt - pi[3]/pY/dt;
                    qmt->cached_qmatrix[0][3] = pi[0] - pi[0]/bt;
                    qmt->cached_qmatrix[1][3] = pi[1] - pi[1]/bt;
                    qmt->cached_qmatrix[2][3] = pi[2] + pi[2]*pR/pY/bt - pi[2]/pY/dt;
                    qmt->cached_qmatrix[3][3] = pi[3] + pi[3]*pR/pY/bt + pi[2]/pY/dt;
        }
    }
    for( i=0; i<4; i++ )
        for( j=0; j<4; j++ )
            if( qmt->cached_qmatrix[j][i]>1.0) fprintf(stderr, "qmatrix (%d, %d; %f): %f\n", i, j, ihky->cached_ep, qmt->cached_qmatrix[j][i]);
}'iHKYNoBoundFixPiMatrixUpdateCache

'static void iHKYNoBoundFixPiMatrixUpdateCache(qmatrix *qmt, ... ) {' double t) {
'static void Matrix_Update_Cache(qmatrix *qmt, ... ) {' double t) {
Public Sub Matrix_Update_Cache(qmt As qmatrix, xxx)


'  const char *fxn_name = "iHKYNoBoundFixPiMatrixUpdateCache";
    'ihkynoboundfixpimatrix *ihky = (ihkynoboundfixpimatrix *) qmt->derived_mt;
    Dim ihky As ihkynoboundfixpimatrix
    ihky = qmt.derived_mt
    
    'double *pi = qmt->pi;
    Dim pi As Double, t As Double, i As Long, j As Long
    pi = qmt.pi
    vargs As va_list
    
    ' Retrieve arguments
    'va_start(vargs, qmt);
    Call va_start(vargs, qmt)
    t = va_arg(vargs, double);
    
    va_end(vargs);

    if( t != ihky->cached_avg_brlen || qmt->v[0] != ihky->cached_ep ) {
        double pY, pR, bt, ct, dt;
        ihky->cached_avg_brlen = t;
        pY = ihky->pY;
        pR = ihky->pR;
        bt = 1 + ihky->b*t;
        ct = 1 + ihky->c*t;
        dt = 1 + ihky->d*t;

' Have transposed cached_qmatrix to speed up likelihood calculation
        if( ihky->model == KAPPA ) {
            ihky->cached_ep = ihky->kappa;
            qmt->cached_qmatrix[0][0] = pi[0] + pi[0]*pY/pR/bt + pi[1]/pR/ct;
            qmt->cached_qmatrix[1][0] = pi[1] + pi[1]*pY/pR/bt - pi[1]/pR/ct;
            qmt->cached_qmatrix[2][0] = pi[2] - pi[2]/bt;
            qmt->cached_qmatrix[3][0] = pi[3] - pi[3]/bt;
            qmt->cached_qmatrix[0][1] = pi[0] + pi[0]*pY/pR/bt - pi[0]/pR/ct;
            qmt->cached_qmatrix[1][1] = pi[1] + pi[1]*pY/pR/bt + pi[0]/pR/ct;
            qmt->cached_qmatrix[2][1] = pi[2] - pi[2]/bt;
            qmt->cached_qmatrix[3][1] = pi[3] - pi[3]/bt;
            qmt->cached_qmatrix[0][2] = pi[0] - pi[0]/bt;
            qmt->cached_qmatrix[1][2] = pi[1] - pi[1]/bt;
            qmt->cached_qmatrix[2][2] = pi[2] + pi[2]*pR/pY/bt + pi[3]/pY/dt;
            qmt->cached_qmatrix[3][2] = pi[3] + pi[3]*pR/pY/bt - pi[3]/pY/dt;
            qmt->cached_qmatrix[0][3] = pi[0] - pi[0]/bt;
            qmt->cached_qmatrix[1][3] = pi[1] - pi[1]/bt;
            qmt->cached_qmatrix[2][3] = pi[2] + pi[2]*pR/pY/bt - pi[2]/pY/dt;
            qmt->cached_qmatrix[3][3] = pi[3] + pi[3]*pR/pY/bt + pi[2]/pY/dt;
        } else if( ihky->model == ALPHA ) {
            ihky->cached_ep = ihky->alpha;
                    qmt->cached_qmatrix[0][0] = pi[0] + pi[0]*pY/pR/bt + pi[1]/pR/ct;
                    qmt->cached_qmatrix[1][0] = pi[1] + pi[1]*pY/pR/bt - pi[1]/pR/ct;
                    qmt->cached_qmatrix[2][0] = pi[2] - pi[2]/bt;
                    qmt->cached_qmatrix[3][0] = pi[3] - pi[3]/bt;
                    qmt->cached_qmatrix[0][1] = pi[0] + pi[0]*pY/pR/bt - pi[0]/pR/ct;
                    qmt->cached_qmatrix[1][1] = pi[1] + pi[1]*pY/pR/bt + pi[0]/pR/ct;
                    qmt->cached_qmatrix[2][1] = pi[2] - pi[2]/bt;
                    qmt->cached_qmatrix[3][1] = pi[3] - pi[3]/bt;
                    qmt->cached_qmatrix[0][2] = pi[0] - pi[0]/bt;
                    qmt->cached_qmatrix[1][2] = pi[1] - pi[1]/bt;
                    qmt->cached_qmatrix[2][2] = pi[2] + pi[2]*pR/pY/bt + pi[3]/pY/dt;
                    qmt->cached_qmatrix[3][2] = pi[3] + pi[3]*pR/pY/bt - pi[3]/pY/dt;
                    qmt->cached_qmatrix[0][3] = pi[0] - pi[0]/bt;
                    qmt->cached_qmatrix[1][3] = pi[1] - pi[1]/bt;
                    qmt->cached_qmatrix[2][3] = pi[2] + pi[2]*pR/pY/bt - pi[2]/pY/dt;
                    qmt->cached_qmatrix[3][3] = pi[3] + pi[3]*pR/pY/bt + pi[2]/pY/dt;
        }
    }
    'for( i=0; i<4; i++ )
    '    for( j=0; j<4; j++ )
    '        if( qmt->cached_qmatrix[j][i]>1.0) fprintf(stderr, "qmatrix (%d, %d; %f): %f\n", i, j, ihky->cached_ep, qmt->cached_qmatrix[j][i]);
End Sub 'Matrix_Update_Cache'iHKYNoBoundFixPiMatrixUpdateCache


static double iHKYNoBoundFixPiMatrixProposer(qmatrix **new_qmt, const qmatrix *qmt, settings *set) {
    ihkynoboundfixpimatrix *ihky = (ihkynoboundfixpimatrix *) qmt->derived_mt;
    double lambda = 1.0;
    double new_param, log_proposal_prob = 0.0;
    if( ihky->model == KAPPA ) {
        new_param = ihky->kappa*exp(lambda*(set->rng->nextStandardUniform(set->rng) - 0.5));
'fprintf(stderr, "iHKY Propose: %f from %f %f\n", new_param, ihky->kappa, lambda);
        log_proposal_prob = 0.0;
    } else if(ihky->model == ALPHA ) {
        new_param = ihky->alpha + set->rng->nextNormal(set->rng, set->sdEP);
        log_proposal_prob = logNormalDensity(new_param, ihky->alpha, set->sdEP);
        while( new_param > 1.0 || new_param < 0.00 ) {
            if( new_param > 1.0 ) new_param = 2.0 - new_param;
            else new_param = -new_param;
        }
    }
    if(*new_qmt) {
        iHKYNoBoundFixPiMatrixReset(*new_qmt, &new_param, qmt->pi);
    } else {
        iHKYNoBoundFixPiMatrixMakeAndSet(new_qmt, ihky->model, &new_param, qmt->pi);
    }
    return log_proposal_prob;
}' iHKYNoBoundFixPiMatrixProposer

' Assumes KAPPA model
static double iHKYNoBoundFixPiMatrixProposeSplit(qmatrix **qmt1, qmatrix **qmt2, const qmatrix *qmt, settings *set, ... ) {'int left, int middle, int right) {
    ihkynoboundfixpimatrix *ihky = (ihkynoboundfixpimatrix *) qmt->derived_mt;
    double zKappa = set->rng->nextStandardNormal(set->rng); ' Random variable used to supplement anemic dimension in rjMCMC's one-to-one mapping between dimensions
    double weight1 = 0, weight2 = 0;    ' Weights used in determining new kappa's from existing kappa
    double pKappa1 = 0, pKappa2 = 0;    ' Proposed kappas for resulting split matrices
    va_list vargs;              ' Used to process variable number of terminal arguments
    int left, middle, right;        ' The variables that will hold the terminal arguments (ends and split point of alignment segment)

    ' Read variable argument list
    va_start(vargs, set);
    left = va_arg(vargs, int);
    middle = va_arg(vargs, int);
    right = va_arg(vargs, int);
    va_end(vargs);

    weight1 = (double)(middle - left) / (right - left);
    weight2 = (double)(right - middle) / (right - left);
    pKappa1 = ihky->kappa * exp( weight2 * set->sigmaAlpha * zKappa );
    pKappa2 = ihky->kappa * exp( - weight1 * set->sigmaAlpha * zKappa );

'fprintf(stderr, "pEPs: %.4f %.4f based on %.4f %.4f %.4f %.4f %.4f and with proposal prob %.4e\n", pKappa1, pKappa2, ihky->kappa, weight1, weight2, set->sigmaAlpha, zKappa, exp(logStandardNormalDensity(zKappa)));
    iHKYNoBoundFixPiMatrixMakeAndSet(qmt1, KAPPA, &pKappa1, qmt->pi);
    iHKYNoBoundFixPiMatrixMakeAndSet(qmt2, KAPPA, &pKappa2, qmt->pi);
    return logStandardNormalDensity(zKappa);
}' iHKYNoBoundFixPiMatrixProposeSplit

static double iHKYNoBoundFixPiMatrixProposeMerge(qmatrix **merged_qmt, const qmatrix *qmt1, const qmatrix *qmt2, settings *set, ... ) {'int left, int middle, int right) {
    ihkynoboundfixpimatrix *ihky1 = (ihkynoboundfixpimatrix *) qmt1->derived_mt;
    ihkynoboundfixpimatrix *ihky2 = (ihkynoboundfixpimatrix *) qmt2->derived_mt;
    double weight1 = 0, weight2 = 0;    ' Weights used in invertible function that maps dimensions in rjMCMC
    double pKappa = 0;          ' Proposed kappa of merged region, invertible function of qmt1->kappa and qmt2->kappa
    va_list vargs;              ' Used to process variable number of terminal arguments
    int left, middle, right;        ' The variables that will hold the terminal arguments (ends and split point of alignment segment)

    ' Read variable argument list
    va_start(vargs, set);
    left = va_arg(vargs, int);  ' Left-most end of split region
    middle = va_arg(vargs, int);    ' Point where region is to be split
    right = va_arg(vargs, int); ' Right-most end of split region
    va_end(vargs);

    weight1 = (double)(middle - left) / (right - left);
    weight2 = (double)(right - middle) / (right - left);
    pKappa = exp(  weight1 * log(ihky1->kappa) + weight2 * log(ihky2->kappa) );

    iHKYNoBoundFixPiMatrixMakeAndSet(merged_qmt, KAPPA, &pKappa, qmt1->pi);
    return logStandardNormalDensity( log( ihky1->kappa / pKappa ) / weight2 / set->sigmaAlpha );
}' iHKYNoBoundFixPiMatrixProposeMerge

static double iHKYNoBoundFixPiMatrixLogPrior(qmatrix *qmt, ... ) { 'double kappa) {
    hky_kappa_hierarchical_prior *khier = (hky_kappa_hierarchical_prior *) qmatrix_prior;
    double kappa, log_kappa;
    va_list vargs;

    if( !qmatrix_prior ) return 0.0;

    va_start(vargs, qmt);
    kappa = va_arg(vargs, double);
    va_end(vargs);
    log_kappa = log(kappa);

    return -1.0*(log_kappa - khier->kappa_mean)*(log_kappa - khier->kappa_mean)/(2.0*khier->kappa_variance) - log_kappa - 0.5*log(khier->kappa_variance) - normal_const;
}' iHKYNoBoundFixPiMatrixLogPrior

static void MakePrior(settings *set) {
    if( !qmatrix_prior ) {
        hky_kappa_hierarchical_prior *pr = (hky_kappa_hierarchical_prior *) malloc(sizeof(hky_kappa_hierarchical_prior));
        if( pr == NULL ) {
            fprintf(stderr, "MakePrior: memory allocation error\n");
            exit(EXIT_FAILURE);
        }
        pr->kappa_mean = set->titv_hyper_mean;
        pr->kappa_variance = set->titv_hyper_variance;
        qmatrix_prior = (void *) pr;
    }
}' MakePrior

/*
void Update_Matrix_SimpleAA(double t, const double lambda, const double *pi) {
    int x,y;

    float Pii = (1.0 / 20.0) + (19.0 / 20.0)* exp(-lambda * t);
    float Pij = (1.0 / 20.0) - (1.0 / 20.0) * exp(-lambda * t);

    if (t <= 0.0000000000000001) {
        Pii = 1.0;
        Pij = 0.0;
    }
    
    
    for (x = 0; x < qmatrix_nchars; x++) {
        for (y = 0; y < qmatrix_nchars; y++) {
            cached_qmatrix[x][y] = Pij;
        }
    }

    for (x = 0; x < qmatrix_nchars; x++) {
        cached_qmatrix[x][x] = Pii;
    }

}
*/

'multiseqdata.c
#include "multiseqdata.h"

' only assigns data array to multiseqdata; does not assign counts, map, or rmap. Currently these are not used' in diverge
void Make_MultiSeqData(multiseqdata *msd, seqdata *all_seq, int in_num_data_sets, int *in_ntaxa) {
        int ***split_data;
    int i,j,n,m,l=0;
        sequence **alignment;
    
    msd->num_datasets = in_num_data_sets;
        msd->total_taxa = 0;

    for (i = 0; i < in_num_data_sets; i++) msd->total_taxa += in_ntaxa[i];

        split_data = (int ***) malloc(sizeof(int **) * msd->num_datasets);
    alignment = (sequence **) malloc(sizeof(sequence *) * msd->num_datasets);
    msd->seq_sets = malloc(sizeof(seqdata) * msd->num_datasets);
    
    for (n=0; n<msd->num_datasets; n++) {
        split_data[n] = (int **) malloc (sizeof(int *) * all_seq->lenunique);
        for (m = 0; m < all_seq->lenunique; m++) {
            split_data[n][m] = (int *) malloc (sizeof(int) * in_ntaxa[n]);
        }

        for (i = 0; i<all_seq->lenunique; i++) {    ' for each site i:
            for (j = 0; j<in_ntaxa[n]; j++) {  ' for each taxa in align n:
                    split_data[n][i][j] = all_seq->data[i][j+l];
            }
        }
    
        alignment[n] = (sequence *) malloc(sizeof(sequence) * in_ntaxa[n]);
        for (i = 0; i < in_ntaxa[n]; i++) {
            alignment[n][i] = all_seq->alignment[i+l];
        }
        l += in_ntaxa[n];

        Set_SeqDataFromSeq(&msd->seq_sets[n], in_ntaxa[n], all_seq, alignment[n]);
        
    }' chop data into num_datasets int[][] holding compressed data for each align
                    
'      msd->data_sets = malloc(sizeof(seqdata) * msd->num_datasets);
'  msd->seq_sets = malloc(sizeof(seqdata) * msd->num_datasets);
'  for (i=0; i<msd->num_datasets; i++) {
'      Set_SeqDataFromSeq(&msd->seq_sets[i], in_ntaxa[i], all_seq, alignment[i]);
'  }

    
} ' Make_MultiSeqData

void Delete_MultiSeqData(multiseqdata *msd) {
    int i;

    for (i = 0; i < msd->num_datasets; i++) {
        SeqDataDelete(&msd->seq_sets[i]);
    }
    free(msd->seq_sets);
} ' Delete_MultiSeqData

' pi should be allocated beforehand; use nchars + 1 to include -9's
' NOTE: this may lead to values of 0 for some elements of pi
void GetPiFromData(multiseqdata *msd, double *pi, int nchars) {
    int i, j, k, total = 0;
    
    for (k = 0; k < nchars; k++) pi[k] = 0;
    
    for (i = 0; i < msd->num_datasets; i++) {
        for (j = 0; j < msd->seq_sets->ntaxa; j++) {
            total = 0; ' final total will represent pi[k] after all sequences
            for (k = 0; k < nchars; k++) {
                pi[k] += msd->seq_sets[i].alignment[j].count[k];
                total += pi[k];
            }
        }
    }
    
    for (k = 0; k < nchars; k++) pi[k] /= total;
}

' for debug
void PrintMultiSeqDataInfo (multiseqdata * msd) {
    int i, j;
    printf("datasets: %d, total taxa: %d\n", msd->num_datasets, msd->total_taxa);
    for (i = 0; i < msd->num_datasets; i++) {
        printf("for set %d\n", i+1);
        for (j = 0; j < msd->seq_sets[i].ntaxa; j++) {
            printf("%s\n", msd->seq_sets[i].alignment[j].name);
        }
    }

}

'node.c
#include "node.h"

' Code herein only works for trees input as ((0,1),(2,3)), where the number corresponds to the 0-based sequence index in the PHYLIP file.
' Tree arrangement not important, e.g. (0,(1,(2,3))) also works
' I think unrooted trees also work, e.g. (0,1,(2,3)) also works?

static int debug = 0;               ' Set to positive integer to turn on local debugging output
static const char this_filename[7] = "node.c";

void printstuff(char *i, char *j) {
    char *ptr1 = i, *ptr2 = j;
    do {
        fputc(*ptr1, stderr);
    } while(ptr1++ != ptr2);
        fprintf(stderr, "\n");
}

node *Make_Subtree_From_Tree(const tree *tr, node *nodes, int *nnodes, double **like, boolean **islike, node *parent) {
    node *cnode = &nodes[*nnodes];
    node *onode = &tr->node_list[*nnodes];
    if( onode->is_branch ) {
        Make_Node(cnode, parent, NULL, NULL, (*nnodes)++, -1, onode->branch_length, *like, (*islike)++);
        (*like) += tr->nchars;
        cnode->left = Make_Subtree_From_Tree(tr, nodes, nnodes, like, islike, cnode);
        cnode->right = Make_Subtree_From_Tree(tr, nodes, nnodes, like, islike, cnode);
        cnode->left->brother = cnode->right;
        cnode->right->brother = cnode->left;
    } else {
        Make_Node(cnode, parent, NULL, NULL, (*nnodes)++, onode->id, onode->branch_length, *like, (*islike)++);
        (*like) += tr->nchars;
    }
    return cnode;
}' Make_Subtree_From_Tree

node *Make_Subtree(node *nodes, char *lsptr, char *rsptr, int *nnodes, double **like, boolean **islike, node *parent, const boolean has_branches, const int nchars) {
    double branch_length=0;
    int comma_position, id;
    const char this_fxnname[13] = "Make_Subtree";
    char *tstr = lsptr;'(char *) malloc(sizeof(char)*strlen(lsptr));
    node *cnode = &nodes[*nnodes];
'  strcpy(tstr, lsptr);
'  tstr[rsptr - lsptr] = '\0';


    if( debug>5 || global_debug>5 ) {
        fprintf(stderr, "%s::%s: Entering with string ", this_filename, this_fxnname);
        printstuff(tstr, rsptr);
        do {
            fputc(*tstr, stderr);
        } while(tstr++ != rsptr);
        fprintf(stderr, "\n");
    }

    tstr = lsptr;

    if( *lsptr == '(' ) {   ' This is a branch
        char *newrsptr;

        branch_length = has_branches ? Read_Branch_Length(rsptr) : -9;
        if( debug>5 || global_debug>5 ) fprintf(stderr, "%s::%s: read branch length %f\n", this_filename, this_fxnname, branch_length);
        comma_position = Find_Branch_Split(lsptr);
        Make_Node(cnode, parent, NULL, NULL, *nnodes, -1, branch_length, *like, (*islike)++);
        (*nnodes)++;
        (*like) += nchars;
        cnode->left = Make_Subtree(nodes, lsptr+1, lsptr + comma_position - 1, nnodes, like, islike, cnode, has_branches, nchars);

        newrsptr = rsptr;
        while(newrsptr != lsptr && newrsptr != NULL && *newrsptr != ')') {
            newrsptr--;
        }
        if (newrsptr == lsptr) newrsptr = rsptr;
        
        if( debug>5 || global_debug>5 ) printstuff(lsptr + comma_position + 1, newrsptr -1);
        cnode->right = Make_Subtree(nodes, lsptr + comma_position + 1, newrsptr-1, nnodes, like, islike, cnode, has_branches, nchars);
        cnode->left->brother = cnode->right;
        cnode->right->brother = cnode->left;
    } else {
        if( debug>5 || global_debug>5 ) {
            fprintf(stderr, "Make node from else\n");
            fprintf(stderr, "current: \n");
            do {
                fputc(*tstr, stderr);
            } while(tstr++ != rsptr);
            fprintf(stderr, "\n");
        
            while(rsptr != NULL && *rsptr != ')') {
                fputc(*rsptr, stderr);
                rsptr--;
            }
        }
        
        branch_length = has_branches ? Read_Branch_Length(rsptr) : -9;
        if( debug>5 || global_debug>5 ) fprintf(stderr, "%s::%s: read branch length %f\n", this_filename, this_fxnname, branch_length);
        id = atoi(lsptr);
        Make_Node(cnode, parent, NULL, NULL, (*nnodes), id, branch_length, *like, (*islike)++);
        (*nnodes)++;
        (*like) += nchars;
    }
    return cnode;
}' Make_Subtree


void Make_Node(node *this, node *up, node *left, node *right, const int uid, const int id, const double branch, double *like, boolean *blike) {
    if( debug>5 || global_debug>5 ) {
        fprintf(stderr, "Make_Node: making %d (%d) with ", uid, id);
        if(up) fprintf(stderr, " parent %d and ", up->uid);
        if(left) fprintf(stderr, " left child %d ", left->uid);
        if(right) fprintf(stderr, " right child %d", right->uid);
        fprintf(stderr, "\n");
    }
    this->is_root = false;
    this->up = up;
    this->left = left;
    this->right = right;
    this->id = id;
    this->uid = uid;
    if(id == -1) this->is_branch = true;
    else this->is_branch = false;
    this->branch_length = branch;
    this->is_likelihood_done = blike;
    this->clikelihood = like;
    this->descending_clades = NULL;
    this->ascending_clades = NULL;
    this->descending_clade = -1;
    this->ascending_clade = -1;
    if( debug>5 || global_debug>5 ) fprintf(stderr, "node %d\tuid %d\tstock %f\n", this->id, this->uid, this->branch_length);
     
}' Make_Node

void Delete_Node(node *this) {
          if(this->clikelihood) free(this->clikelihood);
}' Delete_Node

double Read_Branch_Length(char *rtreeptr) {
    char *treeptr = rtreeptr;
    double branch;
    int a;
    const char this_fxnname[19] = "Read_Branch_Length";

    if (*treeptr == ')') return -9;
    
    while( *treeptr != ':' ) treeptr--;
    if((a = sscanf(treeptr, ": %lf", &branch)) != 1) {
        fprintf(stderr, "%s::%s: Invalid tree string format (%s)", this_filename, this_fxnname, treeptr);
        exit(EXIT_FAILURE);
    }
    return branch;
}' Read_Branch_Length

int Find_Branch_Split(char *ltreeptr) { ' Find dividing ',' for this node
    int nLevel = 0;
    int comma_position = 0;
    char *treeptr = ltreeptr;
    const char this_fxnname[18] = "Find_Branch_Split";

    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s::%s: entering...\n", this_filename, this_fxnname);
    while( (*treeptr != ',') || (nLevel != 1) ) {
      if( *treeptr == '(' ) nLevel++;
      if( *treeptr == ')' ) nLevel--;
      treeptr++;
      comma_position++;
    }

    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s::%s: leaving with comma position at %d\n", this_filename, this_fxnname, comma_position);
    return comma_position;
}' Find_Branch_Split

void Clear_Likelihood(node *n) {
    (*n->is_likelihood_done) = false;
    if( n->left ) {
        Clear_Likelihood(n->left);
        Clear_Likelihood(n->right);
    }
}' Clear_Likelihood

' Very dangerous (assumes str is sufficiently allocated)
void toString(char *str, const node *n, boolean branch) {
    if( n->is_branch ) {
        str = strcat(str, "(");
        toString(str, n->left, branch);
        str = strcat(str, ",");
        toString(str, n->right, branch);
        str = strcat(str, ")");
        if (branch && n->branch_length != -9) sprintf(str,"%s:%f", str,n->branch_length);
    } else {
        char tstr[100];
        if (branch && n->branch_length != -9) sprintf(tstr, "%d:%f", n->id, n->branch_length);
        else sprintf(tstr, "%d", n->id);
        strcat(str, tstr);
    }
}' toString

node *Brother(const node n) {
    if(n.up == NULL) {
        fprintf(stderr, "node.c::Brother: cannot call this function on root\n");
        exit(EXIT_FAILURE);
    }
    if (n.uid == n.up->right->uid) return (n.up->left);
    return (n.up->right);
}' Brother

boolean isLeftChild(const node *n) {
    if (n->up == NULL) {
        fprintf(stderr, "node.c::isLeftChild: cannot call this function on root\n");
        exit(EXIT_FAILURE);
    }
    if (n->uid == n->up->left->uid) return true;
    return false;
}' isLeftChild

int Number_Nodes(const node n) {
    int rtn = 0;
    const node *c = &n;
    while( c->up != NULL ) {
        rtn++;
        c = c->up;
    }
    return rtn;
}' Number_Nodes

int Number_Children(const node n) {
    if( n.left != NULL ) return(1 + Number_Children((*(n.left))) + Number_Children((*(n.right))) );
    return 0;
}' Number_Children

int Number_Leaves(const node *n) {
    if(n->left != NULL) return(Number_Leaves(n->left) + Number_Leaves(n->right));
    return 1;
}' Number_Leaves

int CountDown(const node n) {
    if( n.left != NULL ) return(1 + CountDown(*(n.right)) + CountDown((*(n.left))));
    return 1;
}' CountDown

int Balance(node *n) {
    if(n->is_branch) {
        int i = Balance(n->left);
        int j = Balance(n->right);
        if(i < j) {
            return i;
        } else {
            node *tmp = n->left;
            n->left = n->right;
            n->right = tmp;
            return j;
        }
    } else {
        return n->id;
    }
}' Balance

void Cladify(node *n, const int nleaves, const int nclades, const int *clades) {
    int i, num_clades;
    if( n->descending_clades ) return;
    n->descending_clades = (boolean *)malloc(sizeof(boolean)*nclades);
    for( i=0; i<nclades; i++ )
        n->descending_clades[i] = false;
    if( n->is_branch ) {
        Cladify(n->left, nleaves, nclades, clades);
        Cladify(n->right, nleaves, nclades, clades);
        for( i=0, num_clades=0; i<nclades; i++ ) {
            if( n->left->descending_clades[i] || n->right->descending_clades[i] ) {
                n->descending_clades[i] = true;
                num_clades++;
                n->descending_clade = i;
            }
        }
        if( num_clades > 1 ) n->descending_clade = -1;
    } else {
        n->descending_clade = clades[n->id];
        n->descending_clades[n->descending_clade] = true;
    }
}' Cladify

void Cladify_Up(node *n, const int nleaves, const int nclades, const int *clades) {
    int i, num_clades;
    if( n->up ) {
        boolean *brother, *ascend;
        if(!n->ascending_clades) {
            n->ascending_clades = (boolean *)malloc(sizeof(boolean)*nclades);
            for( i=0; i<nclades; i++ ) n->ascending_clades[i] = false;
        }
        if( isLeftChild(n) ) {
            Cladify(n->up->right, nleaves, nclades, clades);
            brother = n->up->right->descending_clades;
        } else {
            Cladify(n->up->left, nleaves, nclades, clades);
            brother = n->up->left->descending_clades;
        }
        if( n->up->up ) {
            ascend = n->up->ascending_clades;
            for( i=0, num_clades=0; i<nclades; i++ ) {
                if( brother[i] || ascend[i] ) {
                    n->ascending_clades[i] = true;
                    num_clades++;
                    n->ascending_clade = i;
                }
            }
            if( num_clades > 1 ) n->ascending_clade = -1;
        } else {
            memcpy(n->ascending_clades, brother, sizeof(boolean));
            n->ascending_clade = n->brother->descending_clade;
        }
    }
    if( n->is_branch ) {
        Cladify_Up(n->left, nleaves, nclades, clades);
        Cladify_Up(n->right, nleaves, nclades, clades);
    }
}' Cladify_Up

'      0
'   1     5
'       2   6
'      3 4
boolean Clade_Ancestor_Exists(node *n, int nleaves, int clade_id, int nclade_members) {
    if( Clade_Ancestor(n, nleaves, clade_id, nclade_members) )
        return true;
    else if( n->is_branch ) {
        return Clade_Ancestor_Exists(n->left, nleaves, clade_id, nclade_members) || Clade_Ancestor_Exists(n->right, nleaves, clade_id, nclade_members);
    }
    return false;
}' Clade_Ancestor_Exists

boolean Clade_Ancestor(node *n, int nleaves, int clade_id, int nclade_members) {
    return ( ( n->descending_clade == clade_id && Number_Leaves(n) == nclade_members ) ||
        ( n->up && n->ascending_clade == clade_id && nleaves - Number_Leaves(n) == nclade_members ) );
}' Clade_Ancestor

void Clear_Clades(node *n) {
    if(n->descending_clades) free(n->descending_clades);
    if(n->ascending_clades) free(n->ascending_clades);
    n->descending_clade = -1;
    n->ascending_clade = -1;
    if( n->is_branch ) {
        Clear_Clades(n->left);
        Clear_Clades(n->right);
    }
}' Clear_Clades

' To simulate on a tree with branches, first check if branches are defined.
' If so, use the branches in the tree.
' If not, read extra arguments (...) from the argument list and generate a branch according to assumptions in settings arg
void Simulate_Down_Branch(node *n, node *parent, double mu, qmatrix *qmt, rngen *rng) {
    int i;
    double brlen = rng->nextExponential(rng, parent->is_root&&parent->right==n?0.0:mu);
    double rn, cum = 0.0;

    qmt->Matrix_Update_Cache(qmt, brlen);
    rn = rng->nextStandardUniform(rng);
    for( i=0; i<qmt->nchars; i++ ) {
        if( rn < cum + qmt->cached_qmatrix[i][parent->state] ) {
            n->state = i;
            break;
        }
        cum += qmt->cached_qmatrix[i][parent->state];
    }
    if( n->is_branch ) {
        Simulate_Down_Branch(n->left, n, mu, qmt, rng);
        Simulate_Down_Branch(n->right, n, mu, qmt, rng);
    }
}' Simulate_Down_Branch


' DEBUG FUNCTIONS *****************
'
' for debug; dfs traversal on tree
void PrintNodesDFS(const node n) {
    if (!n.left && !n.right) {
        printf("leaf node: %d (leaf %d) with branch length: %f\n\t", n.uid, n.id, n.branch_length);
                'for (i = 0; i < 20; i++) printf("%f ", n.clikelihood);
        printf("\n");
    }
    if (n.left) {
        printf("go left from %d with branch length: %f\n\t", n.uid, n.branch_length);
                'for (i = 0; i < 20; i++) printf("%f ", n.clikelihood);
        printf("\n");
        PrintNodesDFS((*(n.left)));
    }
    if (n.right) {
        printf("go right from %d with branch length: %f\n\t", n.uid, n.branch_length);
                'for (i = 0; i < 20; i++) printf("%f ", n.clikelihood);
        printf("\n");
        PrintNodesDFS((*(n.right)));
    }
}' PrintNodeDFS

partition.C
#include "partition.h"

'static int debug = 0;             ' Set to positive integer for local debugging output
'static const char *file_name = "partition.c";

' File-local function declarations
static void PartitionDefaults(partition *);                     ' Set structure elements to default values

static void PartitionDefaults(partition *prt) {
'  const char *fxn_name = "PartitionDefaults";
    prt->left = 0;
    prt->right = 0;                 ' Right end of a partition/segment
    prt->cPartialLogLikelihood = 0.0;
    prt->cHyperParameter = 0.0;
    prt->cPartialLogHyperParameterPrior = 0.0;
    prt->ctree = NULL;
    prt->cmatrix = NULL;
    prt->counts = NULL;
    prt->topchange = true;
    prt->parchange = true;
    prt->doUpdate = true;
    prt->doXiUpdate = true;
}' PartitionDefaults

void PartitionMake(partition **prt, int lenuniq, int lt, int rt, boolean bpt, boolean bpc) {
'  const char *fxn_name = "PartitionMake";
    *prt = (partition *) malloc(sizeof(partition)); ' Creates BUGS elsewhere
    PartitionDefaults(*prt);
    (*prt)->lenunique = lenuniq;
    (*prt)->counts = (int *) malloc(sizeof(int)*lenuniq);
    PartitionReset(*prt, lt, rt, bpt, bpc);
}' PartitionMake

void PartitionMakeCopy(partition **new_prt, const partition *old_prt) {
'  const char *fxn_name = "PartitionMakeCopy";
    *new_prt = (partition *) malloc(sizeof(partition));
    (*new_prt)->left = old_prt->left;
    (*new_prt)->right = old_prt->right;
    (*new_prt)->cPartialLogLikelihood = old_prt->cPartialLogLikelihood;
    (*new_prt)->cHyperParameter = old_prt->cHyperParameter;
    (*new_prt)->cPartialLogHyperParameterPrior = old_prt->cPartialLogHyperParameterPrior;
    (*new_prt)->cTree = old_prt->cTree;
    (*new_prt)->ctree = old_prt->ctree;
    (*new_prt)->topchange = old_prt->topchange;
    (*new_prt)->parchange = old_prt->parchange;
    (*new_prt)->doUpdate = old_prt->doUpdate;
    (*new_prt)->doXiUpdate = old_prt->doXiUpdate;
    (*new_prt)->lenunique = old_prt->lenunique;
    (*new_prt)->cmatrix = NULL;
    old_prt->cmatrix->Matrix_Make_Copy(&(*new_prt)->cmatrix, old_prt->cmatrix);
    (*new_prt)->counts = (int *) malloc(sizeof(int)*(*new_prt)->lenunique);
    memcpy((*new_prt)->counts, old_prt->counts, (*new_prt)->lenunique*sizeof(int));
}' PartitionMakeCopy

void PartitionCopy(partition *new_prt, partition *old_prt) {
'  const char *fxn_name = "PartitionCopy";

    new_prt->left = old_prt->left;
    new_prt->right = old_prt->right;
    new_prt->cPartialLogLikelihood = old_prt->cPartialLogLikelihood;
    new_prt->cHyperParameter = old_prt->cHyperParameter;
    new_prt->cPartialLogHyperParameterPrior = old_prt->cPartialLogHyperParameterPrior;
    new_prt->cTree = old_prt->cTree;
    new_prt->ctree = old_prt->ctree;
    new_prt->topchange = old_prt->topchange;
    new_prt->parchange = old_prt->parchange;
    new_prt->doUpdate = old_prt->doUpdate;
    new_prt->lenunique = old_prt->lenunique;
    new_prt->cmatrix = old_prt->cmatrix;    ' BUGGY: copying a pointer, but memory managment of matrices elsewhere
    memcpy(new_prt->counts, old_prt->counts, new_prt->lenunique*sizeof(int));
}' PartitionCopy

void PartitionReset(partition *prt, int lt, int rt, boolean bpt, boolean bpc) {
    prt->left = lt;
    prt->right = rt;
    prt->topchange = bpt;
    prt->parchange = bpc;
}' PartitionReset

void PartitionCopyCounts(partition *prt, const seqdata *sqd) {
    memcpy(prt->counts, sqd->counts, sizeof(int)*sqd->lenunique);
}' PartitionCopyCounts

void PartitionCopySegmentCounts(partition *prt, const seqdata *sqd, int left, int right) {
    int *indexSeq = sqd->map;
    int i;
    for(i=0; i<sqd->lenunique; i++) prt->counts[i] = 0; ' Reinitialize
    for(i=left; i<right; i++) {
        prt->counts[indexSeq[i]]++;
    }
}' PartitionCopySegmentCounts

void PartitionCopyPartitionCountDifferences(partition *prt, const partition *fprt, const partition *pprt) {
    int i;
    for(i=0; i<prt->lenunique; i++) prt->counts[i] = fprt->counts[i] - pprt->counts[i];
}' PartitionCopyInverseSegmentCounts

void PartitionCopyPartitionCountDifferences2(partition *prt, const partition *fprt, const partition *pprt1, const partition *pprt2) {
    int i;
    for(i=0; i<prt->lenunique; i++) prt->counts[i] = fprt->counts[i] - pprt1->counts[i] - pprt2->counts[i];
}' PartitionCopyInverseSegmentCounts2

void PartitionCopyPartitionSum(partition *prt, const partition *prt1, const partition *prt2) {
    int i;
    for(i=0; i<prt->lenunique; i++) prt->counts[i] = prt1->counts[i] + prt2->counts[i];
}' PartitionCopyPartitionSum

void PartitionAddPartition(partition *prt, const partition *aprt) {
    int i;
    for(i=0; i<prt->lenunique; i++) prt->counts[i] += aprt->counts[i];
}' PartitionAddPartition

void PartitionSubtractPartition(partition *prt, const partition *aprt) {
    int i;
    for(i=0; i<prt->lenunique; i++) prt->counts[i] -= aprt->counts[i];
}' PartitionSubtractPartition

void PartitionPrintCounts(partition *prt) {
    int i;
    for( i=0; i<prt->lenunique; i++ ) fprintf(stderr, "%d\n", prt->counts[i]);
}' PartitionPrintCounts

void PartitionDelete(partition *prt) {
    if(!prt) return;
    if(prt->counts) free(prt->counts);
    free(prt);
    prt = NULL;
}' PartitionDelete

'partition_list.c
#include "partition_list.h"

#include "dcpsampler.h"

' BEGIN FILE-WIDE VARIABLE DECLARATIONS/DEFINITIONS

static int debug = 0;       ' Set to positive integer for file-specific debugging output
static double current_llike;    ' For DEBUGGING only
static double proposed_llike;   ' For DEBUGGING only
static sampler *gsmp = NULL;

' BEGIN FUNCTION DECLARATIONS

/**
 ** The following functions are used externally (via function pointers) to add/delete change points
 ** in various ways.
 **/

/* function: Propose new change point of specified type (4th argument) that does not land
 * on an existing change point of the same type.
 */
static int Propose_Change_Point(const partition_list *, settings *, boolean);

/* function: Propose new change point of specified type (6th argument) within partitions
 * indexed by 3rd and 4th argument that is NOT the same as the specified location
 * (5th argument) or any other existing change point in the region.
 */
static int Propose_Second_Change_Point(const partition_list *, settings *, int, int, int, boolean);

/* function: Select one of the existing parameter change points to delete.
 * Corresponding functions for topology change points are defined elsewhere (tree_vec)
 * since there are topology dependencies between segments.
 */
static int Propose_Parameter_Change_Point_To_Delete(const partition_list *, settings *);

/* function: Set move probabilities for various add/delete moves and fixed dimension moves.
 */
static void Update_Move_Probabilities(partition_list *, const settings *);

/**
 ** The following functions are used locally to update change point locations.
 **/

/* function: Compute log ratio (actually just likelihood ratio) for the proposed change point
 * location move.  Creates a new partition (4th argument) if the landing partition is split.
 * 5th argument will hold landing index
 * 6th argument will hold vector of partial log likelihoods for all segments that will be wholesale swapped to new parameter values
 * 7th argument will hold log likelihood of partial land part after update (with new parameter values)
 * 8th argument will hold log likelihood of partial land part before update (with old parameter values)
 * 9-10th arguments are the proposed and current location of change point, respectively
 * 11-13th arguments are the leftmost (or rightmost for next) index of the previous, current, and next partition with different parameter values
 * 14th argument indicates whether we are moving topology or parameter change point.
 */
static double LogRatioForMovingChgPt(const sampler *, const partition_list *, partition **, int *, double *, double *, double *, int, int, int, int, int, boolean);

/* function: Reset partitions, counts, and data after a proposed change point location move is accepted.
 * 2nd argument is the new partition created in LogRatioForMovingChgPt function if necessary.
 * 3th argument is current index of moved change point
 * 4th argument is land index of moved change point
 * 5th argument is vector of proposed partial log likelihoods for all regions affected in their entirety
 * 6th argument is the proposed log likelihood of the partial landing segment (if the landing partition will be split)
 * 7th argument is the current log likelihood of the partial landing segment (if the landing partition will be split)
 * 8th argument indicates if we are moving topology or parameter change point
 */
static int MoveCPToLeft(partition_list *, partition *, int, int, int, double *, double, double, boolean);

/* function: reset partitions, counts, and data after a proposed change point location move is accepted.
 * Reverse function to above, but arguments are the same, though the logic is by asymmetry of partitions, slightly different.
 * Extra int argument (before current_index) is index of previous region with different parameter values.
 */
static int MoveCPToRight(partition_list *, partition *, int, int, int, int, double *, double, double, boolean);

static double log_factorial(int);
static double log_factorial_stop(int, int);
static double K_Log_Prior(const partition_list *);
static double J_Log_Prior(const partition_list *);
static double Xi_Log_Prior(const partition_list *);
static double Rho_Log_Prior(const partition_list *);


' BEGIN FUNCTION DEFINITIONS

void PartitionListMake(partition_list **pl, const sampler *smp, int nparts) {
    const char *fxn_name = "PartitionListMake";
    int i;

    *pl = (partition_list *)malloc(sizeof(partition_list));
    if( pl == NULL ) {
        fprintf(stderr, "%s: Memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    (*pl)->part = (partition **)malloc(sizeof(partition *)*nparts);
    if( (*pl)->part == NULL ) {
        fprintf(stderr, "%s: Memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    (*pl)->alignment_length = smp->sqd->lenseq;
    (*pl)->npartitions = nparts;
    (*pl)->parameter_changes = 0;
    (*pl)->topology_changes = 0;
    (*pl)->top_lambda = smp->set->dcp->top_lambda;
    (*pl)->top_lambda_squared = smp->set->dcp->top_lambda * smp->set->dcp->top_lambda;
    (*pl)->log_top_lambda = log(smp->set->dcp->top_lambda);
    (*pl)->par_lambda = smp->set->dcp->par_lambda;
    (*pl)->log_par_lambda = log(smp->set->dcp->par_lambda);

    (*pl)->Update_Move_Probabilities = &Update_Move_Probabilities;
    (*pl)->Propose_Change_Point = &Propose_Change_Point;
    (*pl)->Propose_Second_Change_Point = &Propose_Second_Change_Point;
    (*pl)->Propose_Parameter_Change_Point_To_Delete = &Propose_Parameter_Change_Point_To_Delete;
    (*pl)->K_Log_Prior = &K_Log_Prior;
    (*pl)->J_Log_Prior = &J_Log_Prior;
    (*pl)->Xi_Log_Prior = &Xi_Log_Prior;
    (*pl)->Rho_Log_Prior = &Rho_Log_Prior;

    for( i=0; i<nparts; i++ ) (*pl)->part[i] = NULL;

    ' Use this function pointer to swap in different priors on K (or J too); currently hard-coded in acceptance probability calculations
'  (*pl)->Topology_Change_Point_Prior_Ratio = &KPrior;
    
    ' Must rely on other assistance to make topology partitions due to dependency structure across partitions
'handled in topology_vector    (*pl)->Propose_Topology_Change_Point_To_Delete = &Propose_Topology_Change_Point_To_Delete;
'handled in topology_vector    (*pl)->Propose_Second_Topology_Change_Point_To_Delete = &Propose_Second_Topology_Change_Point_To_Delete;
}' PartitionListMake

void PartitionListMakeCopy(partition_list **npl, const partition_list *opl) {
    int i;
    *npl = (partition_list *) malloc(sizeof(partition_list));
    (*npl)->part = (partition **) malloc(sizeof(partition *)*opl->npartitions);
    for( i=0; i<opl->npartitions; i++ ) {
        (*npl)->part[i] = NULL;
        PartitionMakeCopy(&((*npl)->part[i]), opl->part[i]);
    }
    (*npl)->npartitions = opl->npartitions;
    (*npl)->topology_changes = opl->topology_changes;
    (*npl)->parameter_changes = opl->parameter_changes;
    (*npl)->top_lambda = opl->top_lambda;
    (*npl)->log_top_lambda = opl->log_top_lambda;
    (*npl)->top_lambda_squared = opl->top_lambda_squared;
    (*npl)->par_lambda = opl->par_lambda;
    (*npl)->log_par_lambda = opl->log_par_lambda;
    (*npl)->alignment_length = opl->alignment_length;
    (*npl)->top_one_bk = opl->top_one_bk;
    (*npl)->top_two_bk = opl->top_two_bk;
    (*npl)->top_one_dk = opl->top_one_dk;
    (*npl)->top_two_dk = opl->top_two_dk;
    (*npl)->par_bk = opl->par_bk;
    (*npl)->par_dk = opl->par_dk;
    (*npl)->top_one_bkm1 = opl->top_one_bkm1;
    (*npl)->top_two_bkm2 = opl->top_two_bkm2;
    (*npl)->top_one_dkp1 = opl->top_one_dkp1;
    (*npl)->top_two_dkp2 = opl->top_two_dkp2;
    (*npl)->Propose_Parameter_Change_Point_To_Delete = opl->Propose_Parameter_Change_Point_To_Delete;
    (*npl)->Propose_Topology_Change_Point_To_Delete = opl->Propose_Topology_Change_Point_To_Delete;
    (*npl)->Propose_Two_Topology_Change_Points_To_Delete = opl->Propose_Two_Topology_Change_Points_To_Delete;
    (*npl)->Propose_Change_Point = opl->Propose_Change_Point;
    (*npl)->Propose_Second_Change_Point = opl->Propose_Second_Change_Point;
    (*npl)->Update_Move_Probabilities = opl->Update_Move_Probabilities;
    (*npl)->K_Log_Prior = opl->K_Log_Prior;
    (*npl)->J_Log_Prior = opl->J_Log_Prior;
    (*npl)->Xi_Log_Prior = opl->Xi_Log_Prior;
    (*npl)->Rho_Log_Prior = opl->Rho_Log_Prior;
}' PartitionListMakeCopy

static void Update_Move_Probabilities(partition_list *pl, const settings *set) {
    const char *fxn_name = "Update_Move_Probabilities";
    double factor = set->C;

    ' Dimension change probabilities
    pl->top_one_bk = pl->top_lambda / (pl->topology_changes + 1);
    pl->top_two_bk = pl->top_lambda_squared / (pl->topology_changes + 1) / (pl->topology_changes + 2);
    pl->top_one_dk = (double) pl->topology_changes / pl->top_lambda;
    pl->top_two_dk = (double) pl->topology_changes * (pl->topology_changes - 1) / pl->top_lambda_squared;
    pl->par_bk = pl->par_lambda / (pl->parameter_changes + 1);
    pl->par_dk = (double) pl->parameter_changes / pl->par_lambda;

    ' Reverse move probabilities
    pl->top_one_bkm1 = pl->topology_changes ? pl->top_lambda / (double) pl->topology_changes : 0.0;
    pl->top_two_bkm2 = pl->topology_changes>1 ? pl->top_lambda_squared / (pl->topology_changes - 1) / pl->topology_changes : 0.0;
    pl->top_one_dkp1 = (double) (pl->topology_changes + 1) / pl->top_lambda;
    pl->top_two_dkp2 = (double) (pl->topology_changes + 2) * (pl->topology_changes + 1) / pl->top_lambda_squared;
    pl->par_bkm1 = pl->parameter_changes ? pl->par_lambda / (double) pl->parameter_changes : 0.0;
    pl->par_dkp1 = (double) (pl->parameter_changes + 1) / pl->par_lambda;

    ' Fix 'em up so they are probabilities
    if( pl->top_one_bk > 1.0 ) pl->top_one_bk = 1.0;
    if( pl->top_two_bk > 1.0 ) pl->top_two_bk = 1.0;
    if( pl->top_one_dk > 1.0 ) pl->top_one_dk = 1.0;
    if( pl->top_two_dk > 1.0 ) pl->top_two_dk = 1.0;
    if( pl->par_bk > 1.0 ) pl->par_bk = 1.0;
    if( pl->par_dk > 1.0 ) pl->par_dk = 1.0;
    if( pl->top_one_bkm1 > 1.0 ) pl->top_one_bkm1 = 1.0;
    if( pl->top_two_bkm2 > 1.0 ) pl->top_two_bkm2 = 1.0;
    if( pl->top_one_dkp1 > 1.0 ) pl->top_one_dkp1 = 1.0;
    if( pl->top_two_dkp2 > 1.0 ) pl->top_two_dkp2 = 1.0;
    if( pl->par_bkm1 > 1.0 ) pl->par_bkm1 = 1.0;
    if( pl->par_dkp1 > 1.0 ) pl->par_dkp1 = 1.0;

    ' Multiply by a mixing parameter that should insure all move probabilities do not exceed 1
    pl->top_one_bk *= factor;
    pl->top_two_bk *= factor;
    pl->top_one_dk *= factor;
    pl->top_two_dk *= factor;
    pl->par_bk *= factor;
    pl->par_dk *= factor;
    pl->top_one_bkm1 *= factor;
    pl->top_two_bkm2 *= factor;
    pl->top_one_dkp1 *= factor;
    pl->top_two_dkp2 *= factor;
    pl->par_bkm1 *= factor;
    pl->par_dkp1 *= factor;
    
    ' DEBUGGING:
    if( !set->add_xi ) {
        pl->top_two_bk = 0.0;
        pl->top_two_dk = 0.0;
        pl->top_one_bk = 0.0;
        pl->top_one_dk = 0.0;
    }
    if( !set->add_rho ) {
        pl->par_bk = 0.0;
        pl->par_dk = 0.0;
    }

    ' Double-check that the user has made a smart choice on the tuning parameter.
    if( pl->top_one_bk + pl->top_two_bk + pl->top_one_dk + pl->top_two_dk + pl->par_bk + pl->par_dk > 0.9 ) {
        fprintf(stderr, "%s: You probably don't want the sum of all move probabilities to exceed 0.9 (%f).\n", fxn_name, pl->top_one_bk + pl->top_two_bk + pl->top_one_dk + pl->top_two_dk + pl->par_bk + pl->par_dk);
        exit(EXIT_FAILURE);
    }
}' Update_Move_Probabilities

void PartitionListAddPartition(partition_list *pl, partition *new_part, int insert_index) {
    const char *fxn_name = "PartitionListAddPartition";
    partition **new_partitions = (partition **) malloc(sizeof(partition *)*(pl->npartitions + 1));
    int i;

    if( new_partitions == NULL ) {
        fprintf(stderr, "%s: memory allocation error", fxn_name);
        exit(EXIT_FAILURE);
    }

    if( insert_index > pl->npartitions ) {
        fprintf(stderr, "PartitionListAddPartition: exceed array bounds\n");
        exit(EXIT_FAILURE);
    }

    for( i=0; i<pl->npartitions; i++ ) {
        if( i<insert_index )    new_partitions[i] = pl->part[i];
        else            new_partitions[i+1] = pl->part[i];
    }
    new_partitions[insert_index] = new_part;
    pl->npartitions++;
    free(pl->part);
    pl->part = new_partitions;
}' PartitionListAddPartition

void PartitionListRemovePartition(partition_list *pl, int delete_index) {
    const char *fxn_name = "PartitionListRemovePartition";
    partition **new_partitions = (partition **) malloc(sizeof(partition *)*(pl->npartitions - 1));
    int i;

    if( new_partitions == NULL ) {
        fprintf(stderr, "%s: memory allocation error", fxn_name);
        exit(EXIT_FAILURE);
    }
    if( delete_index >= pl->npartitions || delete_index < 0 ) {
        fprintf(stderr, "PartitionListRemovePartition: exceed array bounds\n");
        exit(EXIT_FAILURE);
    }

    for( i=0;i<pl->npartitions-1; i++ ) {
        if( i < delete_index )  new_partitions[i] = pl->part[i];
        else            new_partitions[i] = pl->part[i+1];
    }
    if( pl->part[delete_index]->cmatrix ) QMatrixDelete(pl->part[delete_index]->cmatrix);
    if( pl->part[delete_index] ) PartitionDelete(pl->part[delete_index]);
    pl->npartitions--;
    free(pl->part);
    pl->part = new_partitions;
}' PartitionListRemovePartition

int PartitionContaining(const partition_list *pl, int new_point, int start, int end) {
    const char *fxn_name = "PartitionContaining";
    int cpart_index = start;
    partition *cpart = NULL;

    if( start < 0 || start >= pl->npartitions || start > end || end >= pl->npartitions ) {
        fprintf(stderr, "%s: exceeded array bounds\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    do {
        cpart = pl->part[cpart_index];
        if( new_point >= cpart->left && new_point <= cpart->right ) return cpart_index;
    } while ( ++cpart_index <= end );

    fprintf(stderr, "%s: could not find partition containing %d in specified range (%d, %d)\n", fxn_name, new_point, start, end);
    exit(EXIT_FAILURE);

    return -1;
}' PartitionContaining

static double K_Log_Prior(const partition_list *pl) {
    int k = pl->topology_changes;
    return - pl->top_lambda + k*pl->log_top_lambda - log_factorial(k);
}' K_Log_Prior

static double J_Log_Prior(const partition_list *pl) {
    int j = pl->parameter_changes;
    return - pl->par_lambda + j*pl->log_par_lambda - log_factorial(j);
}' J_Log_Prior

static double Xi_Log_Prior(const partition_list *pl) {
    int k = pl->topology_changes;
    return log_factorial(k) - log_factorial_stop(pl->alignment_length-1, pl->alignment_length-k);
}' Xi_Log_Prior

static double Rho_Log_Prior(const partition_list *pl) {
    int j = pl->parameter_changes;
    return log_factorial(j) - log_factorial_stop(pl->alignment_length-1, pl->alignment_length-j);
}' Rho_Log_Prior

static double log_factorial(int n) {
    if( n<=1 ) return 0.0;
    else return log((double)n) + log_factorial(n-1);
}' log_factorial

static double log_factorial_stop(int n, int t) {
    if( n<t ) return 0.0;
    else if( n==t ) return log((double)n);
    else return log((double)n) + log_factorial_stop(n-1, t);
}' log_factorial_stop

static int Propose_Change_Point(const partition_list *pl, settings *set, boolean topchange) {
'  const char *fxn_name = "Propose_Change_Point";
    boolean good_point = false;
    int new_point = 0, i;

    while( !good_point ) {
        new_point = (int) (set->rng->nextStandardUniform(set->rng)*pl->alignment_length);
        good_point = true;
        for( i=0; i<pl->npartitions; i++ ) {
            partition *cpart = pl->part[i];
            if( ((topchange && cpart->topchange) || (!topchange && cpart->parchange)) && new_point == cpart->left ) good_point = false;
        }
    }

    return new_point;
}' Propose_Change_Point

static int Propose_Second_Change_Point(const partition_list *pl, settings *set, int left_index, int right_index, int first_xi, boolean topchange) {
    boolean good_point = false;
    int new_point = 0;
    int i;

    do {
        good_point = true;
        new_point = (int) ( set->rng->nextStandardUniform(set->rng) *
                ( pl->part[right_index-1]->right - pl->part[left_index]->left ) + pl->part[left_index]->left + 1);  ' WAS BUG: could proposed lenSeq
        if( new_point == first_xi ) {
            good_point = false;
            continue;
        }
        for( i=left_index; i<right_index; i++ ) {
            partition *cpart = pl->part[i];
            if( ((topchange && cpart->topchange) || (!topchange && cpart->parchange)) && new_point == cpart->left ) {
                good_point = false;
                break;
            }
        }
    } while( ! good_point );

    return new_point;
}' Propose_Second_Change_Point

/**
 * Updates the locations of change points
 */

void UpdateChangePointLocations(partition_list *pl, sampler *smp, boolean topchange, boolean alawadhi) {
    const char *fxn_name = "UpdateChangePointLocations";
    int lowerBound;             ' Left parameter change point
    int upperBound;             ' Right parameter change point
    int current;                ' Current parameter change point to move (within above bounds)
    int proposed;               ' New parameter change point location (within bounds)
    int land_index;             ' Index of partition proposed to receive newly located change point
    double proposed_small_like = 0.0;   ' Proposed likelihood of portion of landing partition over which change point moves
    double current_small_like = 0.0;    ' Current likelihood of portion of landing parititon over which change point moves
    double logRatio;
    int curr_index;             ' Index of region just right of change point we're moving
    int prev_index = 0;         ' Index of region just right of change point previous to one we are moving
    int next_index;             ' Index of region just right of change point next to one we are moving
    partition *prev_part;           ' Indicates the next parameter region to the left of the change point we're moving
    partition *curr_part;           ' Indicates the region just right of the change point we're moving
    partition *new_part = NULL;     ' If needed, this will hold the new partition counts
    boolean local_debug = false;
    int old_nchanges, new_nchanges, i;

    curr_index = 1;
    if( topchange ) {
        while( curr_index < pl->npartitions && !pl->part[curr_index]->topchange ) curr_index++;
    } else {
        while( curr_index < pl->npartitions && !pl->part[curr_index]->parchange ) curr_index++;
    }
    ' If are no change points of the right type (nothing to do and this shouldn't happen because of earlier checks)
    if( curr_index == pl->npartitions ) return;

    do {
        double *pPartialLogLikelihood = (double *) malloc(sizeof(double)*pl->npartitions);  ' VERIFIED 2/11/05
        if( pPartialLogLikelihood == NULL ) {
            fprintf(stderr, "%s: Memory allocation error\n", fxn_name);
            exit(EXIT_FAILURE);
        }
        old_nchanges = pl->parameter_changes;

        prev_part = pl->part[prev_index];
        curr_part = pl->part[curr_index];

        next_index = curr_index + 1;
        if( topchange ) {
            while( next_index < pl->npartitions && !pl->part[next_index]->topchange ) next_index++;
        } else {
            while( next_index < pl->npartitions && !pl->part[next_index]->parchange ) next_index++;
        }

        ' Proposed likelihoods for segments between old and new change point location
        current = curr_part->left;
   
        ' Determine the bounds for the proposed value
        lowerBound = prev_part->left + 1;
        upperBound = pl->part[next_index - 1]->right;   ' Do this because next_index may be pl->npartitions

        ' There is room to propose a new change point
        if( upperBound - lowerBound > 0 ) {

            '  Propose a new change point location that is symmetric and reflected and does not equal the original value
            proposed = ProposeNewChangePointPosition(smp, current, lowerBound, upperBound);

            logRatio = (alawadhi ? smp->set->alawadhi_factor : 1.0) * LogRatioForMovingChgPt(smp, pl, &new_part, &land_index, pPartialLogLikelihood, &proposed_small_like, &current_small_like, proposed, current,
                    prev_index, curr_index, next_index - 1, topchange);

            if( debug>1 || global_debug>1 || local_debug )
                smp->Report_Proposal_Statistics(smp, (topchange?"UpdateXi":"UpdateRho"), 0.0, 0.0, 0.0, proposed_llike, current_llike, current, proposed, 0);

            if( !alawadhi ) smp->tries[(topchange?UPDATE_XI:UPDATE_RHO)]++;

            if( logRatio > 0 || smp->set->rng->nextStandardUniform(smp->set->rng) < exp(logRatio) ) {
                if( !alawadhi ) smp->acceptancerate[(topchange?UPDATE_XI:UPDATE_RHO)]++;
                if( proposed > current ) {  ' Change point moves to right
                    curr_index = MoveCPToRight(pl, new_part, prev_index, curr_index, land_index, proposed, pPartialLogLikelihood, proposed_small_like, current_small_like, topchange);
                } else {            ' Change point moves to left
                    gsmp = smp;
                    curr_index = MoveCPToLeft(pl, new_part, curr_index, land_index, proposed, pPartialLogLikelihood, proposed_small_like, current_small_like, topchange);
                }
                if( debug>0 || global_debug>0 ) smp->Report_State(smp, (topchange?"UpdateXi":"UpdateRho"), logRatio, 0.0);
                if( debug>3 || global_debug>3 || global_debug==-1 ) {
                    VerifyLikelihood(smp, false);
                    VerifyCounts(smp, "UpdateXi|Rho", false);
                }
            } else if( pl->part[land_index]->left != proposed ) {   ' Clean up the new part that was created in LogRatioForMovingChgPt
                PartitionDelete(new_part);
            }
            new_nchanges = 0;
            for( i=1; i<pl->npartitions; i++ ) if(pl->part[i]->parchange) {
                new_nchanges++;
                }
            if( new_nchanges > old_nchanges && local_debug ) {
                fprintf(stderr, "***%d***%d****\n", old_nchanges, new_nchanges);
                smp->Report_State(smp, topchange?"UpdateXi":"UpdateRho", logRatio, 0.0);
                exit(0);
            }
        }

        ' Reset topology change points
        prev_index = curr_index;
        curr_index++;
        if( topchange ) {
            while( curr_index < pl->npartitions && !pl->part[curr_index]->topchange ) curr_index++;
        } else {
            while( curr_index < pl->npartitions && !pl->part[curr_index]->parchange ) curr_index++;
        }
        if(pPartialLogLikelihood) free(pPartialLogLikelihood);
        pPartialLogLikelihood = NULL;
    } while( curr_index < pl->npartitions );
}' UpdateChangePointLocations

int ProposeNewChangePointPosition(sampler *smp, int current, int left, int right) {
    int proposed = current;

    while( proposed == current ) {
        int prnd = (int) (smp->set->rng->nextStandardUniform(smp->set->rng)*smp->set->recomb->lenWindow + 1);
        if( smp->set->rng->nextStandardUniform(smp->set->rng) < 0.5 ) prnd *= -1;
        proposed = current + prnd;

        while( proposed < left || proposed > right ) {
            if( proposed < left )
                proposed = 2*left - proposed;
            else    proposed = 2*right - proposed;
        }
    }
    return proposed;
}' ProposeNewChangePointPosition

static double LogRatioForMovingChgPt(const sampler *smp, const partition_list *pl, partition **npart, int *land_index, double *pPartialLogLikelihood, double *proposed_small_like, double *current_small_like, int proposed, int current, int prev_part_index, int curr_part_index, int last_part_index, boolean top_change_point) {
    const char *fxn_name = "LogRatioForMovingTopChgPt";
    partition *land_part;
    partition *prev_part = pl->part[prev_part_index];
    partition *curr_part = pl->part[curr_part_index];
    partition *new_part, *old_part;
    int start_index, end_index, new_region_start, new_region_end, i;
    boolean parch, topch;
    double proposed_likelihood=0, current_likelihood=0, logRatio;

    *land_index = PartitionContaining(pl, proposed, prev_part_index, last_part_index);
    land_part = pl->part[*land_index];

    '             *                                     *>
    ' prev_part---|---curr_part---|---...---|------land_part------|
    '                                       |new segment|
    if( proposed > current ) {
        start_index = curr_part_index;
        end_index = *land_index;
        new_part = prev_part;
        old_part = curr_part;
        new_region_start = land_part->left;
        new_region_end = proposed - 1;
        parch = land_part->parchange;
        topch = land_part->topchange;
    '       <*                                     *
    ' |-----land_part----|---...---|---prev_part---|---curr_part---|
    '        |new segment|
    } else {
        start_index = ( land_part->left == proposed ) ? *land_index : *land_index + 1;
        end_index = curr_part_index;
        new_part = curr_part;
        old_part = prev_part;
        new_region_start = proposed;
        new_region_end = land_part->right;
        topch = top_change_point?true:false;
        parch = top_change_point?false:true;
    }

    if( start_index < 0 || end_index > pl->npartitions ) {
        fprintf(stderr, "%s: debugging write access violation\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    ' Sum intermediate likelihoods of all full-length parameter partitions under both models between old and new topology change point location
    for( i=start_index; i < end_index; i++ ) {
        partition *tpart = pl->part[i];
        ' Likelihood computed with new tree, old evolutionary parameter
        if( compute_likelihood ) {
            if( top_change_point ) pPartialLogLikelihood[i] = TreeLogLikelihood(new_part->ctree, smp, tpart->cmatrix, tpart->counts, tpart->cHyperParameter, false);
            else pPartialLogLikelihood[i] = TreeLogLikelihood(tpart->ctree, smp, new_part->cmatrix, tpart->counts, new_part->cHyperParameter, false);   ' WAS BUG: new_part->counts
        } else pPartialLogLikelihood[i] = 0.0;
        proposed_likelihood += pPartialLogLikelihood[i];
        if( compute_likelihood ) current_likelihood += tpart->cPartialLogLikelihood;
    }

    ' Unless proposed move is onto another change point, a new partition will need to be created
    if( land_part->left != proposed ) {
        PartitionMake(npart, smp->sqd->lenunique, new_region_start, new_region_end, topch, parch);
        PartitionCopySegmentCounts(*npart, smp->sqd, new_region_start, new_region_end + 1);
        if( compute_likelihood ) {
            if( top_change_point ) *proposed_small_like = TreeLogLikelihood(new_part->ctree, smp, land_part->cmatrix, (*npart)->counts, land_part->cHyperParameter, false);
            else *proposed_small_like = TreeLogLikelihood(land_part->ctree, smp, new_part->cmatrix, (*npart)->counts, new_part->cHyperParameter, false);
        } else *proposed_small_like = 0.0;
        proposed_likelihood += *proposed_small_like;
        if( compute_likelihood ) {
            if( top_change_point ) *current_small_like = TreeLogLikelihood(old_part->ctree, smp, land_part->cmatrix, (*npart)->counts, land_part->cHyperParameter, false);
            else *current_small_like = TreeLogLikelihood(land_part->ctree, smp, old_part->cmatrix, (*npart)->counts, old_part->cHyperParameter, false);
        } else *current_small_like = 0.0;
        if( compute_likelihood ) current_likelihood += *current_small_like;
    }
    proposed_llike = proposed_likelihood;   ' for DEBUGGING only
    current_llike = current_likelihood; ' for DEBUGGING only

    logRatio = proposed_likelihood - current_likelihood;
    return(logRatio);
}' LogRatioForMovingChgPt


'             *                                     *>
' prev_part---|---curr_part---|---...---|---land_part---|
'                                       |new segment|
' For nice graphics of each case see MoveParXiToRight and MoveParXiToLeft in parchpt.c
static int MoveCPToRight(partition_list *pl, partition *new_part, int prev_index, int curr_index, int land_index, int prop, double *pPartialLogLikelihood, double psmalllike, double csmalllike, boolean top_change_point) {
    const char *fxn_name = "MoveCPToRight";
    partition *curr_part = pl->part[curr_index];
    partition *prev_part = pl->part[prev_index];
    partition *land_part = pl->part[land_index];
    int new_curr_index, i;

    if( curr_index < 0 || land_index > pl->npartitions ) {
        fprintf(stderr, "%s: debugging memory access violation\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    
    ' Update trees and likelihoods on intermediate parameter partitions
    for( i = curr_index; i < land_index; i++ ) {
        partition *cpart = pl->part[i];
        if( top_change_point ) cpart->ctree = prev_part->ctree;
        else {
            cpart->cmatrix->Matrix_Copy(cpart->cmatrix, prev_part->cmatrix);
            cpart->cHyperParameter = prev_part->cHyperParameter;
            cpart->cPartialLogHyperParameterPrior = prev_part->cPartialLogHyperParameterPrior;
        }
        cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
    }

    ' Change point we're moving was a both parameter and topology change point (so it's got to stay)
    if( (top_change_point && curr_part->parchange) || (!top_change_point && curr_part->topchange) ) {

        ' Change point landed on old parameter change point
        if( land_part->left == prop ) {
' before:    P&T                        T
' prev_part---|---curr_part---|---...---|---land_part---|  top_change_point = false
' after :     T                        P&T
            if( top_change_point ) land_part->topchange = true; ' It is now a topology change point
            else land_part->parchange = true;
            new_curr_index = land_index;
            if( top_change_point ) curr_part->topchange = false;    ' Remove old topology change point
            else curr_part->parchange = false;
        }

        ' Insert a new topology change point (left segment of land_part)
        else {
            if( land_index == curr_index ) {
'            P&T             *>
' prev_part---|---curr_part==land_part---| ' new_part is setup with properties of land_part at left boundary
'             T---new_part---P
                if( top_change_point ) {
                    new_part->topchange = false;
                    land_part->parchange = false;
                } else {
                    new_part->parchange = false;
                    land_part->topchange = false;
                }
            } else {
' before:    P&T                                       *>
' prev_part---|---curr_part---|---...---|----land_part----|
' after :     T                         |---new_part---P
                if( top_change_point ) {
                    land_part->topchange = true;
                    land_part->parchange = false;
                    curr_part->topchange = false;
                } else {
                        land_part->parchange = true;
                        land_part->topchange = false;
                    curr_part->parchange = false;   ' WAS BUG: parchange = true
                }
            }
            PartitionSubtractPartition(land_part, new_part);
            land_part->cPartialLogLikelihood -= csmalllike;
            land_part->left = prop;

            ' Update likelihood
            new_part->cPartialLogLikelihood = psmalllike;

            if( top_change_point ) {
                ' Tree updated
                new_part->ctree = prev_part->ctree;

                ' Other parameters are unchanged
                new_part->cHyperParameter = land_part->cHyperParameter;
                new_part->cPartialLogHyperParameterPrior = land_part->cPartialLogHyperParameterPrior;
                land_part->cmatrix->Matrix_Make_Copy(&new_part->cmatrix, land_part->cmatrix);
            } else {
                ' Tree unchanged
                new_part->ctree = land_part->ctree;

                ' Other parameters updated
                new_part->cHyperParameter = prev_part->cHyperParameter;
                new_part->cPartialLogHyperParameterPrior = prev_part->cPartialLogHyperParameterPrior;
                prev_part->cmatrix->Matrix_Make_Copy(&new_part->cmatrix, prev_part->cmatrix);
            }

            PartitionListAddPartition(pl, new_part, land_index);

            new_curr_index = land_index + 1;
        }
    }

'             *                                     *>
' prev_part---|---curr_part---|---...---|---land_part---|
'                                       |new segment|
    ' Moved point is topology change point only
    else {
        partition *left_part = pl->part[curr_index - 1];

        ' Topology change point moves within the same segment
        if( curr_index == land_index ) {
'             P              *>
' left_part---|---curr_part==land_part---|
' left_part------------------P
            ' Update change points
            curr_part->left = prop;
            left_part->right = prop - 1;

            ' Update counts
            PartitionSubtractPartition(curr_part, new_part);
            PartitionAddPartition(left_part, new_part);

            ' Update likelihoods
            curr_part->cPartialLogLikelihood -= csmalllike;
            left_part->cPartialLogLikelihood += psmalllike;

            new_curr_index = curr_index;
            PartitionDelete(new_part);
        }

        ' Topology change point moves to another segment
        else {

            if( top_change_point ) land_part->topchange = true;
            else land_part->parchange = true;

            ' Topology change point lands on an existing parameter change point
            if( land_part->left == prop ) {
' left_part---P                         T>
' left_part---|---curr_part---|---...---|---land_part---|
' left_part-------------------|---...--P&T
                new_curr_index = land_index - 1;    ' An upstream segment will be removed
            }

            ' Insert new topology change point
            else {
' left_part---P                                        *>
' left_part---|---curr_part---|---...---|----land_part----|
' left_part-------------------|---...---|---new_part---P
                ' Update likelihoods
                land_part->cPartialLogLikelihood -= csmalllike;
                land_part->left = prop;
                PartitionSubtractPartition(land_part, new_part);    ' WAS BUG: forgot to update counts (this line) on land_part
                new_part->cPartialLogLikelihood = psmalllike;

                if( top_change_point ) {
                    land_part->parchange = false;

                    ' Update tree
                    new_part->ctree = prev_part->ctree;

                    ' Other parameters are unchanged
                    new_part->cHyperParameter = land_part->cHyperParameter;
                    new_part->cPartialLogHyperParameterPrior = land_part->cPartialLogHyperParameterPrior;
                    land_part->cmatrix->Matrix_Make_Copy(&new_part->cmatrix, land_part->cmatrix);
                } else {
                        land_part->topchange = false;

                    ' Tree unchanged
                    new_part->ctree = land_part->ctree;

                    ' Update evolutionary parameters
                    new_part->cHyperParameter = prev_part->cHyperParameter;
                    new_part->cPartialLogHyperParameterPrior = prev_part->cPartialLogHyperParameterPrior;
                    prev_part->cmatrix->Matrix_Make_Copy(&new_part->cmatrix, prev_part->cmatrix);
                }

                PartitionListAddPartition(pl, new_part, land_index);

                new_curr_index = land_index;        ' An upstream segment will be removed
            }

            ' Remove old topology change point (physically remove segment)
            left_part->right = curr_part->right;
            PartitionAddPartition(left_part, curr_part);
            left_part->cPartialLogLikelihood += curr_part->cPartialLogLikelihood;' aka pPartialLogLikelihood[curr_index] (already updated)

            PartitionListRemovePartition(pl, curr_index);
        }
    }
    return new_curr_index;
}' MoveCPToRight


'      <*                                     *
' |----land_part----|---...---|---prev_part---|---curr_part---|
'       |new segment|
static int MoveCPToLeft(partition_list *pl, partition *new_part, int curr_index, int land_index, int prop, double *pPartialLogLikelihood, double psmalllike, double csmalllike, boolean top_change_point) {
    const char *fxn_name = "MoveCPToLeft";
    partition *land_part = pl->part[land_index];
    partition *curr_part = pl->part[curr_index];
    int new_curr_index, i;

    if( land_index < -1 || curr_index > pl->npartitions ) {
        fprintf(stderr, "%s: debugging memory access violation\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    ' Update trees and likelihoods on intermediate parameter partitions
    for( i = land_index + 1; i < curr_index; i++ ) {
        partition *cpart = pl->part[i];
        if( top_change_point ) {
            cpart->ctree = curr_part->ctree;
        } else {
            cpart->cmatrix->Matrix_Copy(cpart->cmatrix, curr_part->cmatrix);
            cpart->cHyperParameter = curr_part->cHyperParameter;
            cpart->cPartialLogHyperParameterPrior = curr_part->cPartialLogHyperParameterPrior;
        }
        cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
    }
       
    ' Moved point was topology and parameter change point (it must stay and a new topology change point added)
    if( (top_change_point && curr_part->parchange) || (!top_change_point && curr_part->topchange) ) {
        ' New topology change point lands on existing change point
        if( land_part->left == prop ) {
'<*                                        T&P
' |---land_part---|---...---|---prev_part---|---curr_part---|
' P                                         T
            if( top_change_point ) {
                land_part->topchange = true;
                land_part->ctree = curr_part->ctree;
            } else {
                land_part->parchange = true;
                land_part->cmatrix->Matrix_Copy(land_part->cmatrix, curr_part->cmatrix);
                land_part->cHyperParameter = curr_part->cHyperParameter;
                land_part->cPartialLogHyperParameterPrior = curr_part->cPartialLogHyperParameterPrior;
            }
            land_part->cPartialLogLikelihood = pPartialLogLikelihood[land_index];   ' WAS BUG: used to set to psmalllike, but psmalllike was not set for land_part->left == prop
            new_curr_index = land_index;
            }
        ' Insert a new topology change point
        else {
'     <*                                   T&P
' |---land_part---|---...---|---prev_part---|---curr_part---|
'      P-new_part-|                         T
            ' Update likelihoods
            land_part->cPartialLogLikelihood -= csmalllike;
            PartitionSubtractPartition(land_part, new_part);    ' WAS BUG: forgot this line
            land_part->right = prop - 1;
            new_part->cPartialLogLikelihood = psmalllike;

            if( top_change_point ) {
                ' Update tree
                new_part->ctree = curr_part->ctree;

                ' Other parameters unchanged
                new_part->cHyperParameter = land_part->cHyperParameter;
                new_part->cPartialLogHyperParameterPrior = land_part->cPartialLogHyperParameterPrior;
                land_part->cmatrix->Matrix_Make_Copy(&new_part->cmatrix, land_part->cmatrix);
            } else {
                ' Tree unchanged
                new_part->ctree = land_part->ctree;

                ' Update other parameters
                new_part->cHyperParameter = curr_part->cHyperParameter;
                new_part->cPartialLogHyperParameterPrior = curr_part->cPartialLogHyperParameterPrior;
                curr_part->cmatrix->Matrix_Make_Copy(&new_part->cmatrix, curr_part->cmatrix);
            }

            PartitionListAddPartition(pl, new_part, land_index + 1);

            new_curr_index = land_index + 1;
        }

        ' Remove old topology point
        if( top_change_point ) curr_part->topchange = false;
        else curr_part->parchange = false;
    }
    ' Moved point was topology change point only (it will be removed and a new one added at landing spot)
'    <*                                     *
' |---land_part---|---...---|---prev_part---|---curr_part---|
'     |new segment|          ...left_part---|
    else {
        partition *left_part = pl->part[curr_index - 1];

        ' Moves left into immediate neighbor region
        if( curr_index - 1 == land_index ) {
        
            ' Lands on existing change point
            if( left_part->left == prop ) {
'<*    left_part--P
' |---land_part---|---curr_part---|         ' land_part == left_part
' P-------------------------------|
                ' Merge regions
                left_part->right = curr_part->right;

                ' Update counts
                PartitionAddPartition(left_part, curr_part);

                ' Update likelihood
                left_part->cPartialLogLikelihood = pPartialLogLikelihood[land_index] + curr_part->cPartialLogLikelihood;    ' WAS BUG: used unset psmalllike

                if( top_change_point ) {
                    left_part->topchange = true;

                    ' Update tree
                    left_part->ctree = curr_part->ctree;

                    ' Other parameters unchanged
                } else {
                    left_part->parchange = true;
                    left_part->cHyperParameter = curr_part->cHyperParameter;
                    left_part->cPartialLogHyperParameterPrior = curr_part->cPartialLogHyperParameterPrior;
                    left_part->cmatrix->Matrix_Copy(left_part->cmatrix, curr_part->cmatrix);
                }

                ' Remove old topology change point
                PartitionListRemovePartition(pl, curr_index);
            
                new_curr_index = curr_index - 1;
            }
            ' Does not land on existing change point
            else {
'    <*left_part--P
' |---land_part---|---curr_part---|
' |---P---------------------------|
                ' Update change points
                curr_part->left = prop;
                left_part->right = prop - 1;    ' aka land_part

                ' Update counts
                PartitionAddPartition(curr_part, new_part);
                PartitionSubtractPartition(left_part, new_part);

                ' Update likelihoods
                curr_part->cPartialLogLikelihood += psmalllike;
                left_part->cPartialLogLikelihood -= csmalllike; ' aka land_part

                    new_curr_index = curr_index;
                PartitionDelete(new_part);
            }
        }
        ' There is at least one intermediate change point between start and end position
'      <*                                     *
' |----land_part----|---...---|---prev_part---|---curr_part---|
'       |new segment|
        else {
            int removal_index;

            ' Lands on a parameter change point
            if( land_part->left == prop ) {
'<*                                         P
' |---land_part---|---...---|---prev_part---|---curr_part---|
'                           |---left_part-------------------|
                if( top_change_point ) {
                    land_part->topchange = true;

                    ' Update tree
                    land_part->ctree = curr_part->ctree;
                } else {
                    land_part->parchange = true;
                    ' Update evolutionary parameter
                    land_part->cHyperParameter = curr_part->cHyperParameter;
                    land_part->cPartialLogHyperParameterPrior = curr_part->cPartialLogHyperParameterPrior;
                    land_part->cmatrix->Matrix_Copy(land_part->cmatrix, curr_part->cmatrix);
                }

                ' Update likelihood
                land_part->cPartialLogLikelihood = pPartialLogLikelihood[land_index];   ' WAS BUG: used unset psmalllike
                removal_index = curr_index;

                new_curr_index = land_index;
            }
            ' Lands in a new spot
            else {
'    <*                                     P
' |---land_part---|---...---|---prev_part---|---curr_part---|
'     P--new_part-|         |---left_part-------------------|
                ' Update likelihoods
                new_part->cPartialLogLikelihood = psmalllike;
                land_part->right = prop - 1;
                land_part->cPartialLogLikelihood -= csmalllike;
                PartitionSubtractPartition(land_part, new_part);    ' WAS BUG: forgot this line

                if( top_change_point ) {
                    ' Update tree
                    new_part->ctree = curr_part->ctree;

                    ' Evolutionary parameters are unchanged from land_part
                    new_part->cHyperParameter = land_part->cHyperParameter;
                    new_part->cPartialLogHyperParameterPrior = land_part->cPartialLogHyperParameterPrior;
                    land_part->cmatrix->Matrix_Make_Copy(&new_part->cmatrix, land_part->cmatrix);
                } else {
                    ' Tree unchanged
                    new_part->ctree = land_part->ctree;

                    ' Evolutionary parameters updated
                    new_part->cHyperParameter = curr_part->cHyperParameter;
                    new_part->cPartialLogHyperParameterPrior = curr_part->cPartialLogHyperParameterPrior;
                    curr_part->cmatrix->Matrix_Make_Copy(&new_part->cmatrix, curr_part->cmatrix);
                }

                PartitionListAddPartition(pl, new_part, land_index + 1);

                removal_index = curr_index + 1;
                new_curr_index = land_index + 1;
            }

            ' Merge left_part and curr_part to remove topology change point
            left_part->right = curr_part->right;
            PartitionAddPartition(left_part, curr_part);
            left_part->cPartialLogLikelihood += curr_part->cPartialLogLikelihood;   ' left_part likelihood already updated at top of function

            ' Remove old topology point
            PartitionListRemovePartition(pl, removal_index);
        }
    }
    return new_curr_index;
}' MoveCPToLeft

double KPriorRatio(partition_list *pl, int dim_change) {
    double log_prior_ratio;
    if( !dim_change) return 0;
    log_prior_ratio = dim_change==1 ? pl->log_top_lambda - log(pl->topology_changes + dim_change) :
        2*pl->log_top_lambda  - log(pl->topology_changes + dim_change) - log(pl->topology_changes + dim_change - 1);
    return log_prior_ratio;
}' KPriorRatio

static int Propose_Parameter_Change_Point_To_Delete(const partition_list *pl, settings *set) {
    int i, j=0;
    int rnd = (int) (set->rng->nextStandardUniform(set->rng)*pl->parameter_changes);

    for( i=1; i<pl->npartitions; i++ ) {
        if( pl->part[i]->parchange && j++ == rnd ) return i;
    }
    return -1;
}' Propose_Parameter_Change_Point_To_Delete

void PartitionListDelete(partition_list *pl, boolean delete_qmatrix) {
'  const char *fxn_name = "PartitionListDelete";
    int i;
    if( !pl ) return;
    for( i=0; i<pl->npartitions; i++ ) {
        if(delete_qmatrix) QMatrixDelete(pl->part[i]->cmatrix);
        PartitionDelete(pl->part[i]);
    }
    free(pl->part);
    free(pl);
    pl = NULL;
}' PartitionListDelete

'qmatrix.c
#include "qmatrix.h"

' Global variables
void *qmatrix_prior = NULL;

' File local variables:
'static int debug = 0;         ' Set to positive integer to turn on local debugging output

' File local function declarations:
static void QMatrixSetParameters(qmatrix *, double *, double *);

void QMatrixMake(qmatrix **qmt, int nchars, int nvariables) {
    int i;
    const char *fxn_name = "QMatrixMake";

    (*qmt) = (qmatrix *) malloc(sizeof(qmatrix));
    if( *qmt == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    (*qmt)->nchars = nchars;
    (*qmt)->nvariables = nvariables;
    if( !(*qmt)->nchars ) return;
    (*qmt)->pi = (double *) malloc(sizeof(double)*(*qmt)->nchars);

    if( (*qmt)->pi == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    for( i=0; i<(*qmt)->nchars; i++ )
        (*qmt)->pi[i] = 1.0/(*qmt)->nchars;

    (*qmt)->cached_qmatrix = (double **) malloc(sizeof(double *)*(*qmt)->nchars);
    if( (*qmt)->cached_qmatrix == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    for(i=0; i<(*qmt)->nchars; i++) {
        (*qmt)->cached_qmatrix[i] = (double *) malloc(sizeof(double)*(*qmt)->nchars);
        if( (*qmt)->cached_qmatrix[i] == NULL ) {
            fprintf(stderr, "%s: memory allocation error\n", fxn_name);
            exit(EXIT_FAILURE);
        }
    }
    if( !(*qmt)->nvariables ) return;
    (*qmt)->v = (double *) malloc(sizeof(double)*(*qmt)->nvariables);
    if( (*qmt)->v == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    for(i=0; i<(*qmt)->nvariables; i++)
        (*qmt)->pi[i] = 1.0/(*qmt)->nvariables;
}' QMatrixMake

void QMatrixCopy(qmatrix *update_qmt, const qmatrix *copied_qmt) {
    int i;
    update_qmt->Matrix_Update_Cache = copied_qmt->Matrix_Update_Cache;
    update_qmt->Matrix_Proposer = copied_qmt->Matrix_Proposer;
    update_qmt->Matrix_Copy = copied_qmt->Matrix_Copy;
    update_qmt->Matrix_Make_Copy = copied_qmt->Matrix_Make_Copy;
    update_qmt->Matrix_Delete = copied_qmt->Matrix_Delete;
    update_qmt->nchars = copied_qmt->nchars;
    update_qmt->nvariables = copied_qmt->nvariables;
    update_qmt->log_prior = copied_qmt->log_prior;
    QMatrixUpdateParameters(update_qmt, copied_qmt->v, copied_qmt->pi);
    for( i=0; i<update_qmt->nchars; i++ ) {
        memcpy(update_qmt->cached_qmatrix[i], copied_qmt->cached_qmatrix[i], sizeof(double)*update_qmt->nchars);
    }
}' QMatrixCopy

void QMatrixUpdateParameters(qmatrix *qmt, const double *inV, const double *inPi) {
    memcpy(qmt->v, inV, sizeof(double)*qmt->nvariables);
    memcpy(qmt->pi, inPi, sizeof(double)*qmt->nchars);
}' QMatrixUpdateParameters

void QMatrixSetParameters(qmatrix *qmt, double *inV, double *inPi) {
    qmt->v = inV;
    qmt->pi = inPi;
}' QMatrixSetParameters

void QMatrixDelete(qmatrix *qmt) {
    int i;
    if(!qmt) return;
    if(qmt->v) free(qmt->v);
    if(qmt->pi) free(qmt->pi);
    if(qmt->cached_qmatrix) {
        for(i=0; i<qmt->nchars; i++) {
            free(qmt->cached_qmatrix[i]);
        }
        free(qmt->cached_qmatrix);
    }
    ' new
    if( qmt->derived_mt ) {
        if( qmt->Matrix_Delete ) qmt->Matrix_Delete(qmt->derived_mt);
        else free(qmt->derived_mt);
    }
    ' end new
    free(qmt);
    qmt = NULL;
}' QMatrixDelete

'sampler.c
#include "sampler.h"

'static int debug = 0; ' Set to positive integer for local debugging output

'void SamplerMake(sampler **smp ) {
Public Sub SamplerMake(smp As sampler)
    
    'const char *fxn_name = "SamplerMake";
    Dim fxn_name As String
    fxn_name = "SamplerMake"

    '*smp = (sampler *) malloc(sizeof(sampler));
    
    'if( *smp == NULL ) {
    '
    '    fprintf(stderr, "%s: memory allocation error\n", fxn_name);
    '    exit(EXIT_FAILURE);
    '}
    
    (*smp)->JumpNumber = 0;
    (*smp)->sincePrint = 0;
    (*smp)->nmoves = 0;
    (*smp)->max_move_name_length = 50;  ' Don't depend on this!! --> SEGFAULT
    (*smp)->tries = NULL;
    (*smp)->acceptancerate = NULL;
    (*smp)->move_names = NULL;
    (*smp)->set = NULL;
    (*smp)->sqd = NULL;
    (*smp)->br = NULL;
    (*smp)->derived_smp = NULL;
    
    smp.JumpNumber = 0
    smp.sincePrint = 0
    smp.nmoves = 0
    smp.max_move_name_length = 50  ' Don't depend on this!! --> SEGFAULT
    smp.tries = 0
    smp.acceptancerate = 0
    smp.move_names = 0
    smp.set = 0
    smp.sqd = 0
    smp.br = 0
    smp.derived_smp = 0
End Sub ' SamplerMake

/*
void ChooseStartingPi(const seqdata *sqd, double *pi) {
    Alignment_Composition(sqd, &pi);
}' ChooseStartingPi
*/

' Only used by cpsampler (hard-coded in dcpsampler for speed)
boolean LogMHAccept(const double logRatio, const double random) {
    if( (logRatio > 0) || (random < exp(logRatio)) ) return true;
    return false;
}' LogMHAccept

void SamplerSaveEstimates(const sampler *smp, int length) {
    smp->OutputLine(smp);
    fflush(smp->fout);
    if( smp->JumpNumber >= length ) {
        CloseSampler(smp);
    }
}' SamplerSaveEstimates

void CloseSampler(const sampler *smp) {
    int i;

    for( i=0; i<smp->nmoves; i++ ) {
        if( smp->tries[i] > 0 )
            fprintf(stderr, "# Acceptance rate #%2d (%20s) %f %6d => %6d\n", i, smp->move_names[i], (double) smp->acceptancerate[i]/(double) smp->tries[i], smp->tries[i], smp->acceptancerate[i]);
        else  fprintf(stderr, "# Acceptance rate #%2d (%20s) not used\n", i, smp->move_names[i]);
    }
    fclose(smp->fout);
    exit(EXIT_SUCCESS);
}'method CloseSampler

void SamplerSetNumberMoves(sampler *smp, int n) {
    const char *fxn_name = "SamplerSetNumberMoves";
    int i;
    smp->nmoves = n;
    smp->move_names = (char **) malloc(sizeof(char *)*smp->nmoves);
    if( smp->move_names == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    for( i=0; i<smp->nmoves; i++ ) {
        smp->move_names[i] = (char *)malloc(sizeof(char)*(smp->max_move_name_length+1));
        if( smp->move_names[i] == NULL ) {
            fprintf(stderr, "%s: memory allocation error\n", fxn_name);
            exit(EXIT_FAILURE);
        }
    }
    smp->acceptancerate = (int *) malloc(sizeof(int)*smp->nmoves);
    if( smp->acceptancerate == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    for( i=0; i<smp->nmoves; i++ ) smp->acceptancerate[i] = 0;
    smp->tries = (int *) malloc(sizeof(int)*smp->nmoves);
    if( smp->tries == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    for( i=0; i<smp->nmoves; i++ ) smp->tries[i] = 0;
}' SamplerSetNumberMoves

void SamplerAddMoveName(sampler *smp, int i, const char *name) {
    strcpy(smp->move_names[i], name);
}' SamplerAddMoveName

/**
 * Calculates the natural log density at x under Standard Normal distribution.
 */
double logStandardNormalDensity(const double x) {
    return logOneOverSqrtTwoPi - 0.5 * x * x;
}' logStandardNormalDensity

/**
 * Calculates the natural log density at x under Normal distribution.
 */
double logNormalDensity(const double x, const double mean, const double std_dev) {
    return logOneOverSqrtTwoPi - 0.5 * ( x - mean ) * ( x - mean ) / std_dev / std_dev - log(std_dev);
}' logNormalDensity
    
void SamplerDelete(sampler *smp) {
    if(!smp) return;
    if(smp->tries) free(smp->tries);
    if(smp->move_names[0]) free(smp->move_names[0]);
    if(smp->move_names) free(smp->move_names);
    if(smp->acceptancerate) free(smp->acceptancerate);
    free(smp);
}' SamplerDelete

'seqdata.c
#include "seqdata.h"

static int debug = 0;               ' Set to positive integer for local debugging output
static const char *file_name = "seqdata.c";

void Compress_Data(seqdata *);

void Make_SeqData(seqdata **sd, sequence *inseq, const int ntaxa, const int lenseq, const int num_chars) {
    const char *fxn_name = "Make_SeqData";
    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): entering and setting %d taxa...\n", fxn_name, file_name, ntaxa);
    *sd = (seqdata *)malloc(sizeof(seqdata));
    (*sd)->ntaxa = ntaxa;
    (*sd)->lenseq = lenseq;
    (*sd)->alignment = inseq;
    (*sd)->num_chars = num_chars;
    (*sd)->lenunique = 0;
'  (*sd)->Compress_Data = &Compress_Data;
        Compress_Data(*sd);
}' Make_SeqData

' Sets the member data having to do with unique patterns
' note: we may not have to set counts, map, and rmap if we set them in derived class
void Set_SeqData(seqdata *sd, const int inlenunique, int incounts[], int inmap[], int inrmap[]) {
    sd->lenunique = inlenunique;
    if( debug>3 || global_debug>3 ) fprintf(stderr, "Number of unique sites: %d\n", sd->lenunique);
    sd->counts = incounts;
    sd->map = inmap;
    sd->rmap = inrmap;
}' Set_SeqData

' used to set a cluster from an alignment of all sequences
void Set_SeqDataFromSeq(seqdata *sd, int ntaxa, seqdata *all_seq, sequence *in_alignment) {
    sd->ntaxa = ntaxa;
    sd->alignment = in_alignment;
    sd->lenunique = all_seq->lenunique;
    sd->counts = all_seq->counts;
    sd->map = all_seq->map;
    sd->rmap = all_seq->rmap;
    sd->num_chars = all_seq->num_chars;
}' Set_SeqDataFromSeq

void Compress_Data(seqdata *sqd) {
    const char *fxn_name = "Make_StdData";
    int lenseq = sqd->lenseq;
    int ntaxa = sqd->ntaxa;
    int *a = NULL, *oc = NULL, *c = NULL, *rm = NULL, *m = NULL;
    int i, j, k, l, temp, lenuniq;
    int *y = NULL, **data, *ydata;
    boolean flip, tied;

    if( debug || global_debug ) fprintf(stderr, "%s(%s): alignment has length: %d and %d sequences\n", fxn_name, file_name, lenseq, ntaxa);

    a = (int *) malloc(sizeof(int)*lenseq);     /* ALLOCATE_MEMORY */
    oc = (int *) malloc(sizeof(int)*lenseq);    /* ALLOCATE_MEMORY */

    if( a == NULL || oc == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    for( i=0; i<lenseq; i++ ) {
        a[i] = i;
        oc[i] = 1;
    }

    ' note: previous version assumed indata->alignment[l].data are continuous in memory, which may no longer be true
    
    ' Shellsort
    for( i=1; i<=lenseq/9; i = (3*i + 1) );
    for( ; i>0; i/=3) {
        for( j=i+1; j<=lenseq; j+=1 ) {
            temp = a[j-1];
            k = j;
            do {
                flip = false;
                tied = true;
                l = 0;
                while(l<ntaxa && tied) {
            
                    y = sqd->alignment[l].data;
                    
                    flip = (y[a[k-i-1]] > y[temp]);
                    tied = tied && (y[a[k-i-1]] == y[temp]);
                    'flip = (y[l*lenseq + a[k-i-1]] > y[l*lenseq + temp]);
                    'tied = tied && (y[l*lenseq + a[k-i-1]] == y[l*lenseq + temp]);
                    l++;
                }
                if(flip) {
                    'printf("flip a[%d] and a[%d]\n", k-1, k-i-1);
                    a[k-1] = a[k-i-1];
                    k -= i;
                }
            } while( (k>i) && flip);
            a[k-1] = temp;
        }
    }


    ' Compress sites
    i = 0;
    lenuniq = 0;
    c = (int *) malloc(sizeof(int)*lenseq);     /* ALLOCATE_MEMORY */
    rm = (int *) malloc(sizeof(int)*lenseq);    /* ALLOCATE_MEMORY */
    m = (int *) malloc(sizeof(int)*lenseq);

    if( c == NULL || rm == NULL || m == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    while(i<lenseq) {
        tied = true;
        j = i+1;
        m[a[i]] = lenuniq;
        while(j < lenseq && tied) {
            k = 0;
            while(k<ntaxa && tied) {
                y = sqd->alignment[k].data;
                tied = (y[a[i]] == y[a[j]]);
                'tied = (y[k*lenseq+a[i]] == y[k*lenseq+a[j]]);
                k++;
            }
            if(!tied) break;
            oc[i] = oc[i] + oc[j];
            m[a[j++]] = lenuniq;
        }
        c[lenuniq] = oc[i];
        rm[lenuniq++] = a[i];
        i = j;
    }
    c = (int *) realloc(c, sizeof(int)*lenuniq);
    rm = (int *) realloc(rm, sizeof(int)*lenuniq);
    if(oc) free(oc);
    if(a) free(a);

    ydata = (int *) malloc(sizeof(int)*lenuniq*ntaxa);  /* ALLOCATE_MEMORY */
    data = (int **) malloc(sizeof(int *)*lenuniq);      /* ALLOCATE_MEMORY */

    if( ydata == NULL || data == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    for(i=0; i<lenuniq; i++) {
        data[i] = &ydata[i*ntaxa];
        for(j=0; j<ntaxa; j++) {
            y = sqd->alignment[j].data;
            data[i][j] = y[rm[i]];
        }
    }

    
    Set_SeqData(sqd, lenuniq, c, m, rm);
    sqd->data = data;
}' Compress_Data

' Calcualte empirical pi (maybe less efficient than Composition)
void Uncompressed_Alignment_Composition(const seqdata *sd, double **freq) {
    const char *fxn_name = "Uncompressed_Alignment_Composition";
    double sum=0;
    int i,j,k;

    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): entering...\n", fxn_name, file_name);
    for( i=0;i<sd->num_chars;i++ ) (*freq)[i] = 0;
    for( i=0; i<sd->lenseq; i++ ) {         ' For each sequence
        for( j=0; j<sd->ntaxa; j++ ) {      ' For each taxon
            for( k=0; k<sd->num_chars; k++ ) {' For each character type (e.g. nucleotide)
                if( sd->alignment[j].data[i] == k ) {
                    (*freq)[k]++;
                    sum++;
                }
            }
        }
    }
    for (k=0; k<sd->num_chars; k++) (*freq)[k] /= sum;
}' Uncompressed_Alignment_Composition

void Alignment_Composition(const seqdata *sd, double **freq) {
    const char *fxn_name = "Alignment_Composition";
    int i, j, k, sum=0;
    
    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): entering with %d taxa...\n", fxn_name, file_name, sd->ntaxa);
    if(!sd->lenunique) return Uncompressed_Alignment_Composition(sd, freq);
    for(i=0;i<sd->num_chars;i++) (*freq)[i] = 0;
    for(i=0;i<sd->lenunique;i++) {
        for(j=0;j<sd->ntaxa;j++) {
            for(k=0;k<sd->num_chars;k++) {
                if(sd->alignment[j].data[sd->rmap[i]] == k) {
                    (*freq)[k] += sd->counts[i];
                    sum += sd->counts[i];
                }
            }
        }
    }
    for(k=0; k<sd->num_chars; k++) (*freq)[k] /= sum;
    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): exiting with %lf, %lf, %lf, %lf\n", fxn_name, file_name, (*freq)[0], (*freq)[1], (*freq)[2], (*freq)[3]);
}' Alignment_Composition

FILE *ReadPhylipHeader(const char *fname, int *nseqs, int *nsites) {
    const char *fxn_name = "ReadPhylipHeader";
    FILE *fin;
    char *line = NULL;
    int line_size = 10, line_read, nm;

    line = (char *) malloc (sizeof(char) * line_size);
    if( line == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): Reading header of file: %s\n", fxn_name, file_name, fname);
    
    fin = fopen(fname, "r");
        if (!fin) {
        fprintf(stderr, "Could not open file %s\n", fname);
        exit(EXIT_FAILURE);
    }
    if((line_read = getline(&line, &line_size, fin)) == -1) {   /* MEMORY_ALLOCATION */
        fprintf(stderr, "Error reading file %s\n", fname);
        exit(EXIT_FAILURE);
    }
    nm = sscanf(line, " %d %d \n", nseqs, nsites);
    if(nm != 2) {
        fprintf(stderr, "Failed to read first line of phylip file %s\n", fname);
        exit(EXIT_FAILURE);
    }
    
    'fprintf(stderr, "line not freed in seqdata.c : ReadPhylip");
    if(line) free(line);
    
    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): Number of sequences: %d\n", fxn_name, file_name, *nseqs);
    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): Number of sites: %d\n", fxn_name, file_name, *nsites);
    return fin;
}


/**************/
Public Sub ReadPhylip(sd() As seqdata, fname As String, noGaps As Byte, model As Long)

'    const char *fxn_name = "ReadPhylip";
'    char *line, *strSeq, *allNames, *ptr;
'    int line_size = 0, line_read, i, num_chars = -1;
'    int cSeq, *data, numseqs, numsites;
'    sequence *s;
'    FILE *fin = ReadPhylipHeader(fname, &numseqs, &numsites);
'fprintf(stderr, "Read %d %d\n", numseqs, numsites);
'    if( model == SCP_RECOMB || model == DCP_RECOMB ) num_chars = 4;
'    else if (model == DIVERGE) num_chars = 20;
'    else {
'        fprintf(stderr, "%s(%s) - Model is unknown\n", file_name, fxn_name);
'        exit(EXIT_FAILURE);
'    }
'
'    fprintf(stderr, "In ReadPhylip...\n");'
'
'    s = (sequence *) malloc(sizeof(sequence)*numseqs);          ' MEMORY_ALLOCATED
'    strSeq = (char *) malloc(sizeof(char)*numseqs*(numsites+1));    ' MEMORY_ALLOCATED
'    allNames = (char *) malloc(sizeof(char)*numseqs*(NAME_LENGTH+1));   ' MEMORY_ALLOCATED
'    data = (int *) malloc(sizeof(int)*numseqs*numsites);        ' MEMORY_ALLOCATED
'
'    cSeq = 0;
'
'    for(i=0; i<numseqs; i++) {
'        line = NULL;
'        line_read = getdelim(&line, &line_size, ' ', fin);
'        if(line_read == -1 || line[strlen(line)-1] != ' ') {
'            fprintf(stderr, "Failed to read in names from phylip file\n");
'            if(allNames) free(allNames);
'            if(strSeq) free(strSeq);
'            exit(EXIT_FAILURE);
'        }
'        line[strlen(line)-1] = '\0';
'        strcpy(&allNames[(NAME_LENGTH+1)*cSeq], line);
'        if(debug) fprintf(stdout, "Read in sequence name %s\n", &allNames[(NAME_LENGTH+1)*cSeq]);
'        ' Read in the sequence
'        line_read = getline(&line, &line_size, fin);
'        if(line_read == -1) {
'            fprintf(stderr, "Failed to read in sequences from phylip file.\n");
'            if(allNames) free(allNames);
'            if(strSeq) free(strSeq);
'            exit(EXIT_FAILURE);
'        }
'        ' Chomp newline
'        if(line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
'        ' Remove spaces
'        ptr = line;
'        while(*ptr != '\0') {
'            if(*ptr == ' ') memmove(ptr, ptr+1, strlen(ptr+1)+1);
'            else ptr++;
'        }
'        ' Move the buffer into final sequence resting place
'        strcpy(&strSeq[(numsites+1)*cSeq], line);
'        cSeq++;
'        if( line ) free(line);
'        line = NULL;
'    }
'    while((line_read = getline(&line, &line_size, fin)) != -1) {
'        char *ptr = line;
'        ' Chomp newline
'        if(line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
'        ' Remove spaces
'        ptr = line;
'        while( *ptr != '\0' ) {
'            if(*ptr == ' ') memmove(ptr, ptr+1, strlen(ptr+1)+1);
'            else ptr++;
'        }
'        if(strlen(line) == 0) continue;
'        if(cSeq == numseqs) cSeq = 0;
'        strcat(&strSeq[(numsites+1)*cSeq], line);
'        cSeq++;
'        if( line ) free(line);
'        line = NULL;
'    }
'    if( line ) free(line);
'    line = NULL;
'    if(noGaps) {
'        fprintf(stderr, "ReadPhylip does not yet implement noGaps = 1\n");
'        if(allNames) free(allNames);
'        if(strSeq) free(strSeq);
'        exit(EXIT_FAILURE);
'    } else if(debug) {
'        printf("\tLoaded phylip file with gaps included\n");
'        fflush(stdout);
'    }
'
'    for(i=0; i<numseqs; i++) {
'        Setup_Sequence(&s[i], strSeq, allNames, data, model);
'        strSeq += numsites+1;
'        allNames += NAME_LENGTH+1;
'        data += numsites;
'    }
'    if(debug) printf("Done with ReadPhylip\n");
'    Make_SeqData(sd, s, numseqs, numsites, num_chars);'

End Sub ' ReadPhylip
/***************/

/******** to treat data from multiple files as 1; e.g. 1 alignment ********
void ReadPhylip(seqdata **sd, const char **fname, boolean noGaps, int in_num_files, int model) {
    const char *fxn_name = "ReadPhylip";
    char *line = NULL, *strSeq, *allNames, *ptr, num_chars = 0;
    int line_size = 0, line_read, i;
    int cSeq, *data, file_numseqs, file_numsites = 0, last_numsites = 0;
    int total_numseqs = 0;
    int index = 0;
    sequence *s = NULL;
    int last_seq_num =0;
    int *ntaxa = (int *)malloc(sizeof(int)*in_num_files);

    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): entering with %d files to read\n", fxn_name, file_name, in_num_files);

    if( model == SCP_RECOMB || model == DCP_RECOMB ) num_chars = 4;
    else if( model == DIVERGE ) num_chars = 20;
    else {
        fprintf(stderr, "In ReadPhylip - Model is unknown\n");
        exit(EXIT_FAILURE);
    }

    ' do this first so we can allocate memory for all sequences; had problems getting realloc to work
    for (index = 0; index < in_num_files; index++) {
        FILE *fin = ReadPhylipHeader(fname[index], &file_numseqs, &file_numsites);
        if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): opening %s\n", fxn_name, file_name, fname[index]);
        ntaxa[index] = file_numseqs;
        total_numseqs += file_numseqs;
        fclose(fin);
    }
    
    s = (sequence *) malloc(sizeof(sequence)*total_numseqs);    ' MEMORY_ALLOCATED
    if( s == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    if( debug>5 || global_debug>5 ) fprintf(stderr, "%s(%s): Allocated space for %d sequences\n", fxn_name, file_name, total_numseqs);
    
    for( index = 0; index < in_num_files; index++ ) {
        FILE *fin = ReadPhylipHeader(fname[index], &file_numseqs, &file_numsites);
        if( index != 0 && file_numsites != last_numsites ) {
            fprintf(stderr, "%s: Files do not have the same number of sites. Cannot be processed.\n", fxn_name);
            exit(EXIT_FAILURE);
        }
        last_numsites = file_numsites;
        strSeq = (char *) malloc(sizeof(char)*file_numseqs*(file_numsites+1));  ' MEMORY_ALLOCATED
        if( strSeq == NULL ) {
            fprintf(stderr, "%s: memory allocation error\n", fxn_name);
            exit(EXIT_FAILURE);
        }

        allNames = (char *) malloc(sizeof(char)*file_numseqs*(NAME_LENGTH+1));  ' MEMORY_ALLOCATED
        if( allNames == NULL ) {
            fprintf(stderr, "%s: memory allocation error\n", fxn_name);
            exit(EXIT_FAILURE);
        }

            data = (int *) malloc(sizeof(int)*file_numseqs*file_numsites);          ' MEMORY_ALLOCATED
        if( data == NULL ) {
            fprintf(stderr, "%s: memory allocation error\n", fxn_name);
            exit(EXIT_FAILURE);
        }

      
        cSeq = 0;

        ' Read in the first interleave set with names
        for( i=0; i<file_numseqs; i++ ) {
            line_read = getdelim(&line, &line_size, ' ', fin);
            if(line_read == -1 || line[strlen(line)-1] != ' ') {
                fprintf(stderr, "Failed to read in names from phylip file\n");
                if(allNames) free(allNames);
                if(strSeq) free(strSeq);
                if(data) free(data);
                exit(EXIT_FAILURE);
            }
            line[strlen(line)-1] = '\0';
            strcpy(&allNames[(NAME_LENGTH+1)*cSeq], line);
            if( debug>5 || global_debug>5 ) fprintf(stdout, "%s(%s): Read in sequence name %s\n", fxn_name, file_name, &allNames[(NAME_LENGTH+1)*cSeq]);
            ' Read in the sequence
            line_read = getline(&line, &line_size, fin);
            if(line_read == -1) {
                fprintf(stderr, "Failed to read in sequences from phylip file.\n");
                if(allNames) free(allNames);
                if(strSeq) free(strSeq);
                if(data) free(data);
                exit(EXIT_FAILURE);
            }
            ' Chomp newline
            if(line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
            ' Remove spaces
            ptr = line;
            while( *ptr != '\0' ) {
                if(*ptr == ' ') memmove(ptr, ptr+1, strlen(ptr+1)+1);
                else ptr++;
            }
            ' Move the buffer into final sequence resting place
            strcpy(&strSeq[(file_numsites+1)*cSeq], line);
            cSeq++;
            if(line) free(line);
            line = NULL;
        }
        ' Read in rest of aligned data
        while( (line_read = getline(&line, &line_size, fin)) != -1 ) {
            char *ptr = line;
            char *nptr = ptr;
            ' Chomp newline
            if(line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
            ' Remove spaces
            while( *ptr ) {
                if( !isspace(*ptr) ) *nptr++ = *ptr;
                ptr++;
            }
            *nptr = 0;
            if(strlen(line) == 0) {
                if(line) free(line);
                line = NULL;
                continue;
            }
            if(cSeq == file_numseqs) cSeq = 0;
            strcat(&strSeq[(file_numsites+1)*cSeq], line);
            cSeq++;
            if(line) {
                free(line);
                line = NULL;
            }
        }
        if(line) free(line);
        if(noGaps) {
            fprintf(stderr, "ERROR - %s: does not yet implement noGaps = 1\n", fxn_name);
            if(allNames) free(allNames);
            if(strSeq) free(strSeq);
            if(data) free(data);
            exit(EXIT_FAILURE);
        } else if( debug>5 || global_debug>5 ) {
            fprintf(stdout, "%s(%s): Loaded phylip file with gaps included\n", fxn_name, file_name);
            fflush(stdout);
        }

        ' set up sequences from current file
        for( i = last_seq_num; i <last_seq_num + file_numseqs; i++ ) {
            Setup_Sequence(&s[i], strSeq, allNames, data, model);
            strSeq += file_numsites+1;
            allNames += NAME_LENGTH+1;
            data += file_numsites;
        }
    
        last_seq_num = i;
        fclose(fin);
    }

    if( debug>5 || global_debug>5 ) {
        fprintf(stderr, "%s(%s): %d sequences: \n\n", fxn_name, file_name, total_numseqs);
        for( i = 0; i < total_numseqs; i++ ) {
            PrintSequenceInfo(&s[i], false);
        }
    }
        
    if( debug>5 || global_debug>5 ) fprintf(stdout, "%s(%s): calling Make_SeqData with %d seqs, %d sites, %d chars\n", fxn_name, file_name, total_numseqs, file_numsites, num_chars);
    Make_SeqData(sd, s, total_numseqs, file_numsites, num_chars);
    
}' ReadPhylip
******************/

' TODO: assumes nucleotide sequences (search nuc assumption)
' Assumes same pi throughout (model does too)

'void SimulateAlignment(seqdata **sqd, settings *set) {
Public Sub SimulateAlignment(sqd() As seqdata, setx As settings)
    'const char *fxn_name = "SimulateAlignment";
    Dim fxn_name As String
    fxn_name = "SimulateAlignment"
    
    'static const int SEQ_NAME_LEN = 8;
    Dim SEQ_NAME_LEN As Long
    SEQ_NAME_LEN = 8
    
    Dim i As Long, j As Long, k As Long, l As Long, nleaves
    l = 0
    'tree *tr = NULL;
    Dim tr As tree
    
    'int *sim_data = NULL;
    Dim sim_data As Long
    
    'qmatrix *qmt = NULL;
    Dim qmt As qmatrix
    
    'sequence *s = NULL;
    Dim s As sequence
    
    'char *n = NULL;
    Dim n As Byte
    
    nleaves = 0

    ' Generate the sequences
'    for( k=0; k<set->sim->segments; k++ ) {
'        iHKYNoBoundFixPiMatrixMakeAndSet(&qmt, KAPPA, &set->sim->kappa[k], set->sim->pi);
'        Make_Tree(&tr, set->sim->tree[k], 4);   ' nuc assumption
'        if( !nleaves ) nleaves = tr->nleaves;
'        else if( nleaves != tr->nleaves ) {
'            fprintf(stderr, "%s(%s): invalid trees in sim_tree %d != %d\n", fxn_name, file_name, nleaves, tr->nleaves);
'            exit(EXIT_FAILURE);
'        }
'        if( !s ) s = (sequence *) malloc(sizeof(sequence)*nleaves);             ' MEMORY_ALLOCATED
'        if( !sim_data ) sim_data = (int *) malloc(sizeof(int)*nleaves*set->sim->total_length);  ' MEMORY_ALLOCATED
'        for( i=l; i<l+set->sim->length[k]; i++ ) {
'            Simulate_Position(tr, qmt, set->sim->pi, set->sim->mu[k], set->rng);
'            for( j=0; j<nleaves; j++ ) {
'                sim_data[set->sim->total_length*j + i] = tr->leaf_list[j]->state;
'            }
'        }
'        QMatrixDelete(qmt);
'        TreeDelete(tr);
'        if(tr) free(tr);
'        tr = NULL;
'        l += set->sim->length[k];
'    }

    ' Name the sequences
'    if( !n ) n = (char *) malloc(sizeof(char)*SEQ_NAME_LEN*nleaves);
'    for( i=0; i<nleaves; i++ ) {
'        sprintf(&n[i*SEQ_NAME_LEN], "Sim. %2d", i);
'        Setup_Sequence_From_Data(&s[i], &sim_data[i*set->sim->total_length], set->sim->total_length, &n[i*SEQ_NAME_LEN]);' This function internally assumes nucleotide sequence
'    }

'    Make_SeqData(sqd, s, nleaves, set->sim->total_length, 4);   ' nuc assumption
End Sub ' SimulateAlignment

void SeqDataDelete(seqdata *sd) {
    if( debug>5 || global_debug>5 ) fprintf(stderr, "Entering SeqDataDelete...\n");
    if(sd->counts) free(sd->counts);
    if(sd->map) free(sd->map);
    if(sd->rmap) free(sd->rmap);
    SequenceDelete(sd->alignment);
    if(sd->alignment) free(sd->alignment);
    if(sd->data) {
        if(sd->data[0]) free(sd->data[0]);
        if(sd->data) free(sd->data);
    }
    if( debug>5 || global_debug>5 ) fprintf(stderr, "Leaving SeqDataDelete...\n");
}' SeqDataDelete

void PrintSequences(const seqdata *sd) {
    int i;
    for( i = 0; i < sd->ntaxa; i++ ) {
        printf("%d: %-10s%s\n", i, sd->alignment[i].name, sd->alignment[i].strand);
    }
}' PrintSequences

void PrintSortedSequences(const seqdata *sd) {
    int i, j;
    for (i=0; i<sd->ntaxa; i++) {
        printf("%-10s", sd->alignment[i].name);
        for (j=0; j<sd->lenunique; j++) {
            putchar(sd->alignment[i].strand[sd->rmap[j]]);
        }
        printf("\n");
    }
}' PrintSortedSequences

void PrintCompressedSequences(const seqdata *sd) {
    int i, j, k;
    for (i=0; i<sd->ntaxa; i++) {
        printf("%-10s", sd->alignment[i].name);
        for (j=0; j<sd->lenunique; j++) {
            for (k=0; k<sd->counts[j]; k++) {
                putchar(sd->alignment[i].strand[sd->rmap[j]]);
            }
        }
        printf("\n");
    }
}' PrintCompressedSequences

void PrintCounts(const seqdata *sd) {
    int i;
    printf("counts array: ");
    for (i = 0; i < sd->lenunique; i++) {
        printf("%d ", sd->counts[i]);
    }
    
    printf("\n");
}' PrintCounts

void PrintMap(const seqdata *sd) {
    int i;
    printf("map array: ");
    for (i = 0; i < sd->lenseq; i++) {
        printf("%d ", sd->map[i]);
    }
    printf("\n");
}' PrintMap

void PrintRMap(const seqdata *sd) {
    int i;
    printf("rmap array: ");
    for( i = 0; i < sd->lenunique; i++ ) {
        printf("%d ", sd->rmap[i]);
    }
    printf("\n");
}' PrintMap

' For debugging or information output
void Print_Summary(const seqdata *sd) {
    int max = 0, i;
    for( i=0; i<sd->lenunique; i++ ) {
        if(max < sd->counts[i]) max = sd->counts[i];
    }
    printf("Number of iid data sites     : %d\n", sd->lenseq);
    printf("Number of unique data sites  : %d\n", sd->lenunique);
    printf("Maximum count per unique site: %d\n", max);
}'method Print_Summary

void Print_Data(const seqdata *sqd) {
    int i, j;
    printf("lenunique: %d\n", sqd->lenunique);
    printf("taxa: %d\n", sqd->ntaxa);
    printf("data:\n");
    for(i=0; i<sqd->lenunique; i++) {
        printf("%d: ", i);
        for(j=0; j<sqd->ntaxa; j++) {
            printf("%d", sqd->data[i][j]);
        }
        printf("\n");
    }
}' Print_Data

void Print_Distances(const seqdata *sqd) {
    int i, j, k;
    double **dis = (double **) malloc(sizeof(double *)*sqd->ntaxa);
    printf("Print_Distances:\n");
    for( i=0; i<sqd->ntaxa; i++ ) {
        dis[i] = (double *) malloc(sizeof(double)*sqd->ntaxa);
        for( j=0; j<sqd->ntaxa; j++ ) {
            dis[i][j] = 0.0;
            for( k=0; k<sqd->lenunique; k++ ) {
                if( sqd->alignment[i].strand[sqd->rmap[k]] != sqd->alignment[j].strand[sqd->rmap[k]] ) dis[i][j] += (double) sqd->counts[k];
            }
            dis[i][j] /= (double) sqd->lenseq;
            printf(" %6.4f", dis[i][j]);
        }
        printf("\n");
        free(dis[i]);
    }
    free(dis);
}' Print_Distances
'sequence.c
#include "sequence.h"

static boolean debug = false;
static int char_code_DNA[] = {'A','G','C','T'}; ' Translation taken from DNASequence.java
static int char_code_AA[] =                     ' Translation taken from AASequence.java
        {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};


void PrintSequenceInfo(sequence *s, boolean print_data) {
    int i;
    fprintf(stderr, "%s: %s", s->name, s->strand);

    if( print_data ) {
        printf("\n%s: ", s->name);
        for( i = 0; i < s->length; i++ ) printf("%d", s->data[i]);
    }
    printf("\n");
}

' TODO: assumes DNA sequences; assumes no gaps
void Setup_Sequence_From_Data(sequence *s, int *data, int len, char *name) {
    int i;

    s->length = len;
    s->name = name;
    s->data = data;
    s->strand = (char *) malloc(sizeof(char)*(len+1));
    s->count = (int *) malloc(sizeof(int)*5);
    for( i=0; i<5; i++ ) {
        s->count[i] = 0;
    }
    for( i=0; i<len; i++ ) {
        s->strand[i] = char_code_DNA[data[i]];
        s->count[data[i]]++;
    }
    s->strand[i] = '\0';
}' Make_Sequence_From_Data

void Setup_Sequence(sequence *s, char *inStr, char *name, int *d, int model) {
' Initializer to allocate a sequence in memory from a string
    int i;
    s->length = strlen(inStr);
    if( debug ) printf("set length: %d\n", s->length);
'  s->strand = (char *) malloc(sizeof(char)*(s->length+1));    /* MEMORY_ALLOCATED */
    s->strand = inStr;
    if( debug ) printf("set strand: %s\n", s->strand);
'  s->name = (char *) malloc(sizeof(char)*(NAME_LENGTH+1));
    s->name = name;
    if( debug ) printf("set name: %s\n", s->name);
    s->data = d;
'  if( debug ) printf("set data: %d\n", s->data);
'  s->data = (int *) malloc(sizeof(int)*s->length);        /* MEMORY_ALLOCATED */
    
    if( model == SCP_RECOMB || model == DCP_RECOMB ) {
        s->count = malloc(sizeof(int) * (4+1));
        for (i = 0; i < 5; i++) s->count[i] = 0;
        for(i=0; i<s->length; i++ ) {
            s->data[i] = To_Int_DNA(inStr[i]);
            if (s->data[i] != -9) s->count[s->data[i]]++;
            else s->count[4]++;
        }
    }
    else if( model == DIVERGE ) {
        s->count = malloc(sizeof(int) * (20 + 1));
        for (i = 0; i < 21; i++) s->count[i] = 0;
        for(i=0; i<s->length; i++ ) {
            s->data[i] = To_Int_AA(inStr[i]);
            if (s->data[i] != -9) s->count[s->data[i]]++;
            else s->count[20]++;
        }
    }
    if( debug ) printf("done with setup sequence\n");
}' Setup_Sequence

int To_Int_DNA(char c) {
    int i;
    for(i=0; i<4; i++) {
        if(toupper(char_code_DNA[i]) == c) return i;
    }
    return -9;
}' To_Int_DNA

int To_Int_AA(char c) {
        int i;
        for(i=0; i<20; i++) {
            if(toupper(char_code_AA[i]) == c) return i;
    }
    return -9;
}' To_Int_AA

void toAAString(sequence *s, char *str) {
    int i, length = s->length;
    if (!str) {
        printf("allocate memory of str before calling toAAString\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < length; i++) {
        'printf("s->data[i]: %d\n", s->data[i]);
        if (s->data[i] == -9) str[i] = '-';
            else str[i] = char_code_AA[s->data[i]];
        'printf("set str[i] to %c\n", str[i]);
    }
    'printf("in file str: %s\n", str);
}

void Copy_Sequence(sequence *s, const sequence d) {
    int i;
    s->length = d.length;
    s->name = d.name;
    s->strand = (char *) malloc(sizeof(char)*(s->length+1));    /* MEMORY_ALLOCATED */
    strcpy(s->strand, d.strand);
    s->data = (int *) malloc(sizeof(int)*s->length);        /* MEMORY_ALLOCATED */
    for (i = 0; i < s->length; i++) {
        s->data[i] = d.data[i];
    }
}' Copy_Sequence

void SequenceDelete(sequence *s) {
    if(!s) return;
    if(debug) printf("Clearing sequences\n");
    if(s->data) free(s->data);
    if(s->strand) free(s->strand);
    if(s->name) free(s->name);
    s->length = 0;
}' SequenceDelete

int getBase(sequence sd, const int inBase) {
    return sd.data[inBase];
}' getBase
  
' Composition:  Proportion of integer (code for character) in sequence, excluding non-[AGCT] characters
double Sequence_Composition(sequence sd, const int in) {
    int i;
    double rtn = 0.0;
    double len = 0.0;
    for(i=0; i<sd.length; i++) {
        if( sd.data[i] == in ) rtn++;
        if( sd.data[i] != -9 ) len++;
    }
    return (rtn/len);
}' Sequence_Composition


' AllComposition:  Frequency of character inChar in sequence
double Sequence_Absolute_Composition(sequence sd, const int in) {
    int rtn = 0, i;
    for (i=0; i<sd.length; i++) {
        if (sd.data[i] == in) rtn++;
    }
    return rtn;
}' Sequence_Absolute_Composition

settings.C
#include "settings.h"

static int debug = 1;       ' Set to positive integer for local debugging output
static const char *model_types = "<recomb|diverge|dcp_recomb>";
static const char *file_name = "settings.c";

static double nextStandardUniform(rngen *);
static double nextDispersedUniform(rngen *);
static double nextStandardNormal(rngen *);
static double nextNormal(rngen *, double);
static double nextGamma(rngen *, double, double);
static double nextExponential(rngen *, double);
static void Read_List(void *, char *, int, const char *, int);
static void Process_Alawadhi_Option(settings *, const char *);
static void Process_Change_Point_Option(settings *, const char *);
static char *strip_white(char *str);

void UnknownOption(char *cmd, char *model) {
    fprintf(stderr, "Option %s unknown (for %s) .\n", cmd, model);
    exit(EXIT_FAILURE);
}

/**
* Reads in configuration file and initializes pseudo-random number generator.
*/

void CheckBoundsGr0(const double x, const char *vname) {
    if( x <= 0 ) {
        fprintf(stderr, "Paremeter %s must be > 0.\n", vname);
        exit(-1);
    }
}' CheckBoundsGr0

void CheckBoundsGrEqual0(const double x, const char *vname) {
    if( x < 0 ) {
        fprintf(stderr, "Paremeter %s must be >=0.\n", vname);
        exit(-1);
    }
}' CheckBoundsGrEqual0

void CheckBounds01(const double x, const char *vname) {
    if( x < 0 || x > 1) {
        fprintf(stderr, "Parameter %s must be >=0 and <=1.\n", vname);
        exit(-1);
    }
}' CheckBounds01

void Set_Defaults(settings *set) {
    set->model = -1;
    set->ctmc_model = HKY;
    set->ctmc_parameterization = KAPPA;
    set->recomb = NULL;
    set->diverge = NULL;
    set->dcp = NULL;
    set->scp = NULL;
    set->gmodel = false;
    set->simulate_data = false;
    set->sim = NULL;
    set->rng = (rngen *) malloc(sizeof(rngen));
    set->cmdfile_seed = false;
    
    set->rng->useRnList = false;
    set->rng->nextStandardUniform = nextStandardUniform;
    set->rng->nextDispersedUniform = nextDispersedUniform;
    set->rng->nextStandardNormal = nextStandardNormal;
    set->rng->nextNormal = nextNormal;
    set->rng->nextGamma = nextGamma;
    set->rng->nextExponential = nextExponential;

    (*set).num_pTrees = 0;
    (*set).length = 10000000;
    
    (*set).burnin = 100;
    (*set).subsample = 50;

    set->alawadhi = false;
    set->alawadhi_topology = false;
    set->alawadhi_topology_one = false;
    set->alawadhi_topology_two = false;
    set->alawadhi_parameter = false;
    set->alawadhi_k = 100;
    set->alawadhi_factor = 0.5;

    set->debug = 0;
    set->exit_condition = false;
    set->report_iact = false;
    set->alawadhi_debug = 0;
    set->compute_likelihood = true;
    set->pTree = NULL;

    set->sdEP = 0.1;
    set->weight = log(0.5);
    set->sigmaAlpha = 0.75;
    set->sigmaMu = 0.75;
    set->C = 0.45;
    set->jump_classes = false;
    set->add_rho = false;
    set->add_xi = false;

    set->update_hyperparameters = true;

    ' What are these parameters?
    set->sdP = 0.02;
    set->sdNP = 1;
    set->sdT1 = 0.1;
    set->sdT2 = 0.01;
    set->sdHyperEP = 1;
    set->sdUV = 1;
    set->mix1 = 1.0;

    set->init_string = NULL;

}' Set_Defaults

void Set_Seed(rngen *rng, const unsigned long int seed) {
    if(debug) fprintf(stdout, "Setting random number seed: %ld\n", seed);
    rng->rengine = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng->rengine, seed);
}' Set_Seed

'static double nextStandardUniform(rngen *rng) {
Public Function nextStandardUniform(rng As rngen) As Double

    If rng.useRnList <> 0 Then
        If rng.current_rn = rng.total_rn Then
            'fprintf(stderr, "ERROR: nextStandardUniform\n");
            'exit(EXIT_FAILURE);
        End If
        'return rng->rnList[rng->current_rn++];
        nextStandardUniform = rng.rnList(rng.current_rn + 1)
    Else
        'return gsl_rng_uniform(rng->rengine);
        nextStandardUniform = gsl_rng_uniform(rng.rengine)
    End If
End Function ' nextStandardUniform

static double nextDispersedUniform(rngen *rng) {
    if( rng->useRnList ) {
        if( rng->current_rn == rng->total_rn ) {
            fprintf(stderr, "ERROR: nextDispersedUniform\n");
            exit(EXIT_FAILURE);
        }
        return 100*rng->rnList[rng->current_rn++];
    } else
        return 100*gsl_rng_uniform(rng->rengine);
}' nextDispersedUniform

static double nextStandardNormal(rngen *rng) {
    if( rng->useRnList ) {
        if( rng->current_rn == rng->total_rn ) {
            fprintf(stderr, "ERROR: nextStandardNormal\n");
            exit(EXIT_FAILURE);
        }
        return rng->rnList[rng->current_rn++];  ' BUG - in some sense
    } else
        return gsl_ran_gaussian(rng->rengine, 1);
}' nextStandardNormal

static double nextNormal(rngen *rng, double sd) {
    if( rng->useRnList ) {
        if(rng->current_rn == rng->total_rn) {
            fprintf(stderr, "ERROR: nextNormal\n");
            exit(EXIT_FAILURE);
        }
        return rng->rnList[rng->current_rn++];  ' BUG - in some sense
    } else
        return gsl_ran_gaussian(rng->rengine, sd);
}' nextNormal

static double nextGamma(rngen *rng, double alpha, double beta) {
    if( rng->useRnList ) {
        if( rng->current_rn == rng->total_rn ) {
            fprintf(stderr, "ERROR: nextGamma\n");
            exit(EXIT_FAILURE);
        }
        return rng->rnList[rng->current_rn++];  ' BUG - in some sense
    } else
        return gsl_ran_gamma(rng->rengine, alpha, beta);
}' nextGamma

static double nextExponential(rngen *rng, double mean) {
    return -mean * log( rng->nextStandardUniform(rng) );
}' nextExponential

void Set_Recomb_Defaults(settings *set) {
    
    
    '(*set).recomb->mixBD = 0.90;
    set->recomb->lenWindow = 5;

'  set->recomb->init_alpha_string = NULL;
'  set->recomb->init_mu_string = NULL;
'  set->recomb->init_tree_string = NULL;
'  set->recomb->init_changepoint_string = NULL;
} ' Set_Recomb_Defaults

void Set_SCP_Defaults(settings *set) {
    set->scp->lambda = 1;
    (*set).sdMu = 0.05;
}' Set_SCP_Defaults

void Set_Diverge_Defaults(settings *set) {
    (*set).sdMu = 0.01;
    (*set).diverge->M = 1;
    (*set).diverge->sdTheta = 0.001;
    (*set).diverge->sdAlpha = 0.01;
    (*set).diverge->numMu = 0;
    (*set).diverge->Mu_Default = 0.5;
} ' Set_Diverge_Defaults

void Set_DCP_Defaults(settings *set) {
    set->dcp->par_lambda = 1;
    set->dcp->top_lambda = 1;
'  set->dcp->P = 0.95;
}' Set_DCP_Defaults

void ReadCmdfile(settings **set_ptr, const char *fname) {
    const char *fxn_name = "ReadCmdfile";
    FILE *fin;
    char *cmd = NULL, *line = NULL;
    char *value = NULL, *value_ptr = NULL, *cmd_ptr = NULL;
    int line_num = 0, i;
    boolean set_titv_mean=false, set_titv_variance=false, set_mu_mean=false, set_mu_variance=false;
    settings *set;

    line = (char *) malloc( sizeof(char)* (MAX_LINE_LENGTH+1));
    cmd = (char *) malloc( sizeof(char)* (MAX_LINE_LENGTH+1));
    value = (char *) malloc( sizeof(char)* (MAX_LINE_LENGTH+1));
    if( line == NULL || cmd == NULL || value == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    *set_ptr = (settings *) malloc(sizeof(settings));
    set = *set_ptr;
    Set_Defaults(set);

    fin = fopen(fname, "r");
    if(fin == NULL) {
        fprintf(stderr, "Error opening cmdfile: %s\n", fname);
        exit(EXIT_FAILURE);
    }
    ' Read in cmd (up to and including :)
    while( (line = fgets(line, MAX_LINE_LENGTH, fin)) != NULL ) {

        ' Ignore empty reads (shouldn't happen)
        if( !strlen(line) ) continue;

        ' Remove trailing/leading blank space (including newlines)
        line = strip_white(line);

        ' Ignore comments
        if( strlen(line) && line[0] == '#' )  continue;

        ' Read option name
        strcpy(cmd, line);
        if( (cmd_ptr = strchr(cmd, ':')) != NULL ) {
            cmd_ptr++;
            *cmd_ptr = '\0';
        } else {
            fprintf(stderr, "%s: error reading the following line of the cmdfile:\n%s\n", fxn_name, line);
            exit(EXIT_FAILURE);
        }

        ' Read option value
        strcpy(value, ++cmd_ptr);
        
        ' Remove trailing comments on value (PROBLEM: prevents use of character '#' within option value; is that OK?)
        value_ptr = strchr(value, '#');
        if( value_ptr != NULL ) *value_ptr = '\0';

        ' Remove trailing/leading blank space (including newlines)
        value = strip_white(value);
        cmd = strip_white(cmd);

        line_num++;
        ' read model
        if( line_num == 1 && !strcmp(cmd, "model:") ) {
            if (!strcmp(value, "recomb") ) {
                set->model = SCP_RECOMB;
                set->scp = (set_scp *)malloc(sizeof(set_scp));
                set->recomb = (set_recomb *)malloc(sizeof(set_recomb));
                if( set->scp == NULL || set->recomb == NULL ) {
                    fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                    exit(EXIT_FAILURE);
                }
                Set_Recomb_Defaults(set);
                Set_SCP_Defaults(set);
            }
            else if (!strcmp(value, "diverge") ) {
                set->model = DIVERGE;
                set->diverge = (set_diverge *)malloc(sizeof(set_diverge));
                if( set->diverge == NULL) {
                    fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                    exit(EXIT_FAILURE);
                }
                Set_Diverge_Defaults(set);
            }
            else if (!strcmp(value, "dcp_recomb") ) {
                set->model = DCP_RECOMB;
                set->dcp = (set_dcp *)malloc(sizeof(set_dcp));
                set->recomb = (set_recomb *)malloc(sizeof(set_recomb));
                if( set->dcp == NULL || set->recomb == NULL ) {
                    fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                    exit(EXIT_FAILURE);
                }
                Set_Recomb_Defaults(set);
                Set_DCP_Defaults(set);
            }
            else {
                fprintf(stderr,"ERROR: must specify model: %s in cmdfile\n", model_types);
                exit(EXIT_FAILURE);
            }
            
        } else if (line_num == 1) {
            fprintf(stderr,"ERROR: first line in cmd file must specify model: %s\n", model_types);
            exit(EXIT_FAILURE);
        }
        
        ' read common parameters for all models
        else if( !strcmp(cmd, "length:") ) {
            ' note: the max length will be 2,147,483,647, even if value exceeds it
            set->length = atol(value);
            CheckBoundsGr0(set->length, "length");
            if(debug) fprintf(stdout, "length: %d\n", set->length);
        } else if( !strcmp(cmd, "ctmc_model:") ) {
            if( !strcmp(value, "HKY") )
                set->ctmc_model = HKY;
            else {
                fprintf(stderr, "ERROR: unrecognized option for ctmc_model (Valid options: HKY)\n");
                exit(EXIT_FAILURE);
            }
            if(debug) fprintf(stdout, "ctmc_model: %s\n", value);
        } else if( !strcmp(cmd, "ctmc_parameterization:") ) {
            if( !strcmp(value, "ALPHA") )
                set->ctmc_parameterization = ALPHA;
            else if( !strcmp(value, "KAPPA") )
                set->ctmc_parameterization = KAPPA;
            else {
                fprintf(stderr, "ERROR: unrecognized option for ctmc_model (Valid options: kappa|alpha; kappa default)\n");
                exit(EXIT_FAILURE);
            }
            if(debug) fprintf(stdout, "ctmc_parameterization: %s\n", value);
        } else if( !strcmp(cmd, "gmodel:") ) {
            double mono_prob = atof(value);
            set->gmodel = true;
            CheckBounds01(mono_prob, "gmodel");
            set->log_mono_prob = log(mono_prob);
            if(debug) fprintf(stdout, "gmodel: %f\n", set->log_mono_prob);
        } else if( !strcmp(cmd, "simulate_data:") && !set->simulate_data ) {
            char *ptr = value;
            ' Remove possible comment
            while( !isspace(*ptr) && *ptr != '\0' ) ptr++;
            *ptr = '\0';
            set->simulate_data = !strcmp(value, "false") ? false : true;
            if( !set->simulate_data ) continue;
            if( !set->sim ) set->sim = (set_sim *)malloc(sizeof(set_sim));
            set->sim->segments = 1;
            if( strcmp(value, "true") )
                set->sim->segments = atoi(value);
            if(debug) fprintf(stdout, "simulate_data: true\n");
        } else if( !strcmp(cmd, "sim_mu:") && set->simulate_data ) {
            if( !set->sim ) {
                set->sim = (set_sim *)malloc(sizeof(set_sim));
                set->sim->segments = 1;
            }
            set->sim->mu = (double *) malloc(sizeof(double)*set->sim->segments);
            if(debug) fprintf(stdout, "sim_mu:");
            Read_List(set->sim->mu, value, set->sim->segments, "sim_mu", DOUBLE);
            for( i=0; i<set->sim->segments; i++ ) CheckBoundsGr0(set->sim->mu[i], "sim_mu");
        } else if( !strcmp(cmd, "sim_kappa:") && set->simulate_data ) {
            if( !set->sim ) {
                set->sim = (set_sim *)malloc(sizeof(set_sim));
                set->sim->segments = 1;
            }
            set->sim->kappa = (double *) malloc(sizeof(double)*set->sim->segments);
            if(debug) fprintf(stdout, "sim_kappa:");
            Read_List(set->sim->kappa, value, set->sim->segments, "sim_kappa", DOUBLE);
            for( i=0; i<set->sim->segments; i++ ) CheckBoundsGr0(set->sim->kappa[i], "sim_kappa");
        } else if( !strcmp(cmd, "sim_pi:") && set->simulate_data ) {    ' TODO: assumes nucleotide and trusts the user
            double sum = 0.0;
            if( !set->sim ) {
                set->sim = (set_sim *)malloc(sizeof(set_sim));
                set->sim->segments = 1;
            }
            set->sim->pi = (double *) malloc(sizeof(double)*4); ' assume nucleotide
            if(debug) fprintf(stdout, "sim_pi:");
            Read_List(set->sim->pi, value, 4, "sim_pi", DOUBLE);' assume nucleotide
            for( i=0; i<4; i++ ) {                  ' assume nucleotide
                CheckBounds01(set->sim->pi[i], "sim_pi");
                sum += set->sim->pi[i];
            }
            if( fabs(sum - 1.0) > tolerance ) { ' tolerance is set in constants.h
                fprintf(stderr, "ERROR: stationary distribution (sim_pi) does not sum to 1\n");
                exit(EXIT_FAILURE);
            }
        } else if( !strcmp(cmd, "sim_length:") && set->simulate_data ) {
            if( !set->sim ) {
                set->sim = (set_sim *)malloc(sizeof(set_sim));
                set->sim->segments = 1;
            }
            set->sim->total_length = 0;
            set->sim->length = (int *) malloc(sizeof(int)*set->sim->segments);
            if(debug) fprintf(stdout, "sim_length:");
            Read_List(set->sim->length, value, set->sim->segments, "sim_length", INT);
            for( i=0; i<set->sim->segments; i++ ) {
                CheckBoundsGr0(set->sim->length[i], "sim_length");
                set->sim->total_length += set->sim->length[i];
            }
        } else if( !strcmp(cmd, "sim_tree:") && set->simulate_data ) {
            char *ptr, *rptr;
            if( !set->sim ) {
                set->sim = (set_sim *)malloc(sizeof(set_sim));
                set->sim->segments = 1;
            }
            set->sim->tree = (char **) malloc(sizeof(char *)*set->sim->segments);
            if( set->sim->tree == NULL ) {
                fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                exit(EXIT_FAILURE);
            }
            if(debug) fprintf(stdout, "sim_tree:");
            ptr = value;
            for( i=0; i<set->sim->segments; i++ ) {
                rptr = ptr;
                while( *ptr != ' ' && *ptr != '\0' ) ptr++;
                *ptr = '\0';
                set->sim->tree[i] = (char *) malloc(sizeof(char)*(strlen(rptr)+1));
                strcpy(set->sim->tree[i], rptr);
                ptr++;
            }
        } else if( !strcmp(cmd, "burnin:") ) {
            set->burnin = atoi(value);
            CheckBoundsGrEqual0(set->burnin, "burnin");
            if(debug) fprintf(stdout, "burnin: %d\n", set->burnin);
        } else if( !strcmp(cmd, "subsample:") ) {
            set->subsample = atoi(value);
            CheckBoundsGr0(set->subsample, "subsample");
            if(debug) fprintf(stdout, "subsample: %d\n", set->subsample);
        } else if( !strcmp(cmd, "start_tree:") || !strcmp(cmd, "parent_tree:") ) {
            char **newpTree, *ptr = value;
            ' Remove possible comment
            while( !isspace(*ptr) && *ptr != '\0' ) ptr++;
            *ptr = '\0';
            ' couldn't get realloc to work
            newpTree = (char **) malloc(sizeof(char*) * (set->num_pTrees+1));
            if( newpTree == NULL ) {
                fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                exit(EXIT_FAILURE);
            }
            for (i = 0; i < set->num_pTrees; i++) {
                newpTree[i] = set->pTree[i];
            }
            newpTree[set->num_pTrees] = (char *) malloc(sizeof(char)*(strlen(value) + 1));  ' BUGGY, BUGGY THIS WAS!!
            if( newpTree[set->num_pTrees] == NULL ) {
                fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                exit(EXIT_FAILURE);
            }
            strcpy(newpTree[set->num_pTrees], value);
            if( set->pTree ) free( set->pTree );                        ' BUGGY, BUGGY (used to free unmalloced space)
            set->pTree = newpTree;
            set->num_pTrees++;
            if(debug) fprintf(stdout, "start_tree %d: %s\n", i+1, value);

        } else if( !strcmp(cmd, "compute_likelihood:") ) {
            if( !strcmp(value, "false") ) set->compute_likelihood = false;
            if(debug && set->compute_likelihood) fprintf(stdout, "compute_likelihood: true\n");
            else if(debug) fprintf(stdout, "compute_likelihood: false\n");
            compute_likelihood = set->compute_likelihood;
        } else if( !strcmp(cmd, "mu_hyper_mean:") ) {
            set->mu_hyper_mean = atof(value);
'          CheckBoundsGrEqual0(set->mu_hyper_mean, "mu_hyper_mean");
            if(debug) fprintf(stdout, "mu_hyper_mean: %f\n", set->mu_hyper_mean);
            set_mu_mean = true;
        } else if( !strcmp(cmd, "mu_hyper_variance:") ) {
            set->mu_hyper_variance = atof(value);
            CheckBoundsGr0(set->mu_hyper_variance, "mu_hyper_variance");
            if(debug) fprintf(stdout, "mu_hyper_variance: %f\n", set->mu_hyper_variance);
            set_mu_variance = true;
        } else if( !strcmp(cmd, "titv_hyper_mean:") ) {
            set->titv_hyper_mean = atof(value);
'          CheckBoundsGrEqual0(set->titv_hyper_mean, "titv_hyper_mean");
            if(debug) fprintf(stdout, "titv_hyper_mean: %f\n", set->titv_hyper_mean);
            set_titv_mean = true;
        } else if( !strcmp(cmd, "titv_hyper_variance:") ) {
            set->titv_hyper_variance = atof(value);
            CheckBoundsGr0(set->titv_hyper_variance, "titv_hyper_variance");
            if(debug) fprintf(stdout, "titv_hyper_variance: %f\n", set->titv_hyper_variance);
            set_titv_variance = true;
        } else if( !strcmp(cmd, "sd_mu:") ) {
            set->sdMu = atof(value);
            CheckBoundsGr0(set->sdMu, "sd_mu");
            if(debug) fprintf(stdout, "sd_mu: %f\n", set->sdMu);
        } else if( !strcmp(cmd, "sd_p:") ) {
            set->sdP = atof(value);
            CheckBoundsGr0(set->sdP, "sd_p");
        } else if( !strcmp(cmd, "sd_np:") ) {
            set->sdNP = atof(value);
            CheckBoundsGr0(set->sdNP, "sd_np");
        } else if( !strcmp(cmd, "sd_t1:") ) {
            set->sdT1 = atof(value);
            CheckBoundsGr0(set->sdT1, "sd_t1");
        } else if( !strcmp(cmd, "sd_t2:") ) {
            set->sdT2 = atof(value);
            CheckBoundsGr0(set->sdT2, "sd_t2");
        } else if( !strcmp(cmd, "sd_ep:") ) {
            set->sdEP = atof(value);
            CheckBoundsGr0(set->sdEP, "sd_ep");
            if(debug) fprintf(stdout, "sd_ep: %f\n", set->sdEP);
        } else if( !strcmp(cmd, "sd_hyep:") ) {
            set->sdHyperEP = atof(value);
            CheckBoundsGr0(set->sdHyperEP, "sd_hyep");
            if(debug) fprintf(stdout, "sd_hyep: %f\n", set->sdHyperEP);
        } else if( !strcmp(cmd, "sd_uv:") ) {
            set->sdUV = atof(value);
            CheckBoundsGr0(set->sdUV, "sd_uv");
        } else if( !strcmp(cmd, "exp_weight:") ) {
            set->weight = atof(value);
            CheckBoundsGr0(-set->weight, "exp_weight");
            if(debug) fprintf(stdout, "exp_weight: %f\n", set->weight);
        } else if( !strcmp(cmd, "sigma_alpha:") ) {
            set->sigmaAlpha = atof(value);
            CheckBoundsGr0(set->sigmaAlpha, "sigma_alpha");
            if(debug) fprintf(stdout, "sigma_alpha: %f\n", set->sigmaAlpha);
        } else if( !strcmp(cmd, "sigma_mu:") ) {
            set->sigmaMu = atof(value);
            CheckBoundsGr0(set->sigmaMu, "sigma_mu");
            if(debug) fprintf(stdout, "sigma_mu: %f\n", set->sigmaMu);
        } else if( !strcmp(cmd, "mix1:") ) {
            set->mix1 = atof(value);
            CheckBounds01(set->mix1, "mix1");
            if(debug) fprintf(stdout, "mix1: %f\n", set->mix1);
        } else if( !strcmp(cmd, "C:") ) {
            set->C = atof(value);
            CheckBounds01(set->C, "C");
            if(debug) fprintf(stdout, "C: %f\n", set->C);
        } else if( !strcmp(cmd, "jump_classes:") || !strcmp(cmd, "change_points:") ) {
            char *str_ptr = value;
            if( !strcmp(value, "true") ) {  ' Legacy, backwards compatibility
                set->jump_classes = true;
                set->add_xi = true;
                set->add_rho = true;
                if(debug && set->jump_classes) printf("change_points: true\n");
                continue;
            } else if( !strcmp(value, "false") ) {
                if( debug ) printf("change_points: false\n");
                continue;
            }
            Process_Change_Point_Option(set, str_ptr);
            while( (str_ptr = strchr(str_ptr, '|')) != NULL ) {
                str_ptr++;
                Process_Change_Point_Option(set, str_ptr);
            }
        } else if( !strcmp(cmd, "seed:") ) {
            set->cmdfile_seed = true;
            set->seed = atol(value);
        } else if( !strcmp(cmd, "report_iact:") ) {
            if( !strcmp(value, "true") ) set->report_iact = true;
            else set->report_iact = false;
            if( debug && set->report_iact ) printf("report_iact: true\n");
        } else if( !strcmp(cmd, "exit_condition:") ) {
            if( !strcmp(value, "true") ) set->exit_condition = true;
            else set->exit_condition = false;
            if( debug && set->exit_condition ) printf("exit_condition: true\n");
        } else if( !strcmp(cmd, "random:") ) {
            FILE *finr = fopen(value, "r");
            if(finr == NULL) {
                fprintf(stderr, "Error opening random file: %s\n", value);
                exit(EXIT_FAILURE);
            }
            set->rng->useRnList = true;
            set->rng->total_rn = 1000000;   ' TODO: get this information from the file (yucky to do this portably in C)
            set->rng->rnList = (double *) malloc(sizeof(double)*set->rng->total_rn);
            if( set->rng->rnList == NULL ) {
                fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                exit(EXIT_FAILURE);
            }
            set->rng->current_rn = 0;
            i=0;
            while( (line = fgets(line, MAX_LINE_LENGTH, fin)) != NULL ) {
                set->rng->rnList[i++] = atof(line);
            }
            fclose(finr);
        } else if( !strcmp(cmd, "initial_values:") ) {
            set->init_string = (char *) malloc(sizeof(char)* (strlen(value) + 1) );
            if( set->init_string == NULL ) {
                fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                exit(EXIT_FAILURE);
            }
            strcpy(set->init_string, value);
            if(debug) fprintf(stdout, "initial_values: %s\n", set->init_string);
        } else if( !strcmp(cmd, "alawadhi_debug:") ) {
            set->alawadhi_debug = atoi(value);
            if(debug) fprintf(stdout, "alawadhi_debug: %d\n", set->alawadhi_debug);
        } else if( !strcmp(cmd, "debug:") ) {
            set->debug = atoi(value);
            if(debug) fprintf(stdout, "debug: %d\n", set->debug);
            debug = set->debug ? set->debug : debug;
        } else if( !strcmp(cmd, "window_length:") ) {
            set->recomb->lenWindow = atoi(value);
            CheckBoundsGr0(set->recomb->lenWindow, "window_length");
            if(debug) fprintf(stdout, "window_length: %d\n", set->recomb->lenWindow);
        } else if( !strcmp(cmd, "alawadhi:") ) {
            char *str_ptr = value;
            if( !strcmp(value, "true") ) {  ' Legacy, backwards compatibility
                set->alawadhi = true;
                if(debug && set->alawadhi) printf("alawadhi: true\n");
                continue;
            } else if( !strcmp(value, "false") ) {
                if( debug ) printf("alawadhi: false\n");
                continue;
            }
            Process_Alawadhi_Option(set, str_ptr);
            while( (str_ptr = strchr(str_ptr, '|')) != NULL ) {
                str_ptr++;
                Process_Alawadhi_Option(set, str_ptr);
            }
        } else if( !strcmp(cmd, "alawadhi_k:") ) {
            set->alawadhi_k = atoi(value);
            CheckBoundsGrEqual0((double) set->alawadhi_k, "alawadhi_k");
            if(debug) printf("alawadhi_k: %d\n", set->alawadhi_k);
        } else if( !strcmp(cmd, "alawadhi_factor:") ) {
            set->alawadhi_factor = atof(value);
            CheckBounds01((double) set->alawadhi_factor, "alawadhi_factor");
            if(debug) printf("alawadhi_factor: %f\n", set->alawadhi_factor);
        }

        ' read single changepoint recomb parameters
        if( set->model == SCP_RECOMB ) {
        
            if( !strcmp(cmd, "lambda:") ) {
                set->scp->lambda = atoi(value);
                CheckBoundsGr0(set->scp->lambda, "lambda");
                if(debug) fprintf(stdout, "lambda: %f\n", set->scp->lambda);
            }
            ' } else if( !strcmp(cmd, "mix_bd:") ) {
              '    set->mixBD = atof(value);
              '    CheckBounds01(set->mixBD, "mix_bd");
              '    if(debug) fprintf(stdout, "mix_bd: %f\n", set->mixBD);
            '} else if( !strcmp(cmd, "mix2:") ) {
            '  set->mix2 = atof(value);
            '  CheckBounds01(set->mix2, "mix2");
            '  if(debug) fprintf(stdout, "mix2: %f\n", set->mix2);
            /*
            } else if( !strcmp(cmd, "initial_alphas:")) {
                set->recomb->init_alpha_string = (char *) malloc(sizeof(char)* (strlen(value) + 1) );
                strcpy(set->recomb->init_alpha_string, value);
                if(debug) fprintf(stdout, "initial_values: %s\n", set->recomb->init_alpha_string);
            } else if( !strcmp(cmd, "initial_mus:")) {
                set->recomb->init_mu_string = (char *) malloc(sizeof(char)* (strlen(value) + 1));
                strcpy(set->recomb->init_mu_string, value);
                if(debug) fprintf(stdout, "initial_values: %s\n", set->recomb->init_mu_string);
            } else if( !strcmp(cmd, "initial_trees:")) {
                set->recomb->init_tree_string = (char *) malloc(sizeof(char)*(strlen(value)+1));
                strcpy(set->recomb->init_tree_string, value);
                if(debug) fprintf(stdout, "initial_values: %s\n", set->recomb->init_tree_string);
            } else if( !strcmp(cmd, "initial_changepoints:")) {
                set->recomb->init_changepoint_string = (char *) malloc(sizeof(char)*(strlen(value)+1));
                strcpy(set->recomb->init_changepoint_string, value);
                if(debug) fprintf(stdout, "initial_values: %s\n", set->recomb->init_changepoint_string);
            */
        } ' end recomb parameters

        ' read diverge parameters
        else if (set->model == DIVERGE) {
                if( !strcmp(cmd, "mix3:") ) {
                set->diverge->mix3 = atof(value);
                CheckBounds01(set->diverge->mix3, "mix3");
                if(debug) fprintf(stdout, "mix3: %f\n", set->diverge->mix3);
            '} else if (!strcmp(cmd, "Mu:") ) {
            '  set->diverge->Mu = atof(value);
            '  CheckBounds01(set->diverge->Mu, "Mu");
            '  if(debug) printf("sdT2: %f\n", set->diverge->sdT2);
            } else if (!strcmp(cmd, "Alpha:") ) {
                set->diverge->Alpha = atof(value);
                CheckBounds01(set->diverge->Alpha, "Alpha");
                if(debug) printf("Alpha: %f\n", set->diverge->Alpha);
            } else if (!strcmp(cmd, "Theta:") ) {
                set->diverge->Theta = atof(value);
                CheckBounds01(set->diverge->Theta, "Theta");
                if (debug) printf("Theta: %f\n", set->diverge->Theta);
            } else if (!strcmp(cmd, "sdTheta:") ) {
                set->diverge->sdTheta = atof(value);
                CheckBounds01(set->diverge->sdTheta, "sdTheta");
                if(debug) printf("sdTheta: %f\n", set->diverge->sdTheta);
            } else if (!strcmp(cmd, "sdAlpha:") ) {
                set->diverge->sdAlpha = atof(value);
                CheckBounds01(set->diverge->sdAlpha, "sdAlpha");
                if(debug) printf("sdAlpha: %f\n", set->diverge->sdAlpha);
            } else if( !strcmp(cmd, "Mu:")) {
                double *newMu;
                int i = 0;
                set->diverge->numMu++;
                ' couldn't get realloc to work
                newMu = (double *) malloc(sizeof(double) * set->diverge->numMu);
                if( newMu == NULL ) {
                    fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                    exit(EXIT_FAILURE);
                }
                for (i = 0; i < set->diverge->numMu-1; i++) {
                    newMu[i] = set->diverge->Mu[i];
                }
                'set->pTree = (char **) realloc(**(set->pTree), sizeof(char*) * set->num_pTrees);
                newMu[i] = atof(value);
                CheckBounds01(newMu[i], "Mu");
                free(set->diverge->Mu);
                set->diverge->Mu = newMu;
                if(debug) fprintf(stdout, "Mu %d: %f\n", i+1, set->diverge->Mu[set->diverge->numMu-1]);
            } else UnknownOption(cmd, "diverge");
        } ' end diverge parameters

        ' read DCPM parameters
        else if (set->model == DCP_RECOMB) {
            if( !strcmp(cmd, "par_lambda:") ) {
                set->dcp->par_lambda = atof(value);
                CheckBoundsGr0(set->dcp->par_lambda, "par_lambda");
            } else if( !strcmp(cmd, "top_lambda:") ) {
                set->dcp->top_lambda = atof(value);
                CheckBoundsGr0(set->dcp->top_lambda, "top_lambda");
            '} else if( !strcmp(cmd, "prob:") ) {
            '  set->dcp->P = atof(value);
            '  CheckBounds01(set->dcp->P, "prob");
            }
        }
    }

    if( line ) free(line);
    if( cmd ) free(cmd);
    if( value ) free(value);
    fclose(fin);

    if (set->model == -1) {
        fprintf(stderr, "\tERROR: must specify model: %s in cmdfile\n\n", model_types);
        exit(EXIT_FAILURE);
    }

    if( set_mu_mean && set_mu_variance && set_titv_mean && set_titv_variance ) {
        set->update_hyperparameters = false;
    } else if (set_mu_mean || set_mu_variance || set_titv_mean || set_titv_variance ) {
        fprintf(stderr, "ERROR: you must specify _all_ options: mu_hyper_mean, mu_hyper_variance, titv_hyper_mean, and titv_hyper_variance (%d, %d, %d, %d)\n", set_mu_mean, set_mu_variance, set_titv_mean, set_titv_variance);
        exit(EXIT_FAILURE);
    } else {
        set->update_hyperparameters = true;
    }

    ' diverge error checks
    if (set->model == DIVERGE) {
        if (set->num_pTrees != 2) {
            fprintf(stderr, "\tERROR: must specify 2 starting trees\n");
            exit(EXIT_FAILURE);
        }
        if (set->diverge->numMu == 0) {
            set->diverge->numMu = set->num_pTrees;
            set->diverge->Mu = malloc (sizeof (double) * set->diverge->numMu);
            if( set->diverge->Mu == NULL ) {
                fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                exit(EXIT_FAILURE);
            }
            for (i = 0; i < set->diverge->numMu; i++) {
                set->diverge->Mu[i] = set->diverge->Mu_Default;
            }
            if (debug) printf("setting Mu: %f for each cluster\n", set->diverge->Mu_Default);
        }
        if (set->diverge->numMu == 1) {
            float m = set->diverge->Mu[0];
            free (set->diverge->Mu);
            set->diverge->numMu = set->num_pTrees;
            set->diverge->Mu = malloc (sizeof (double) * set->diverge->numMu);
            if( set->diverge->Mu == NULL ) {
                fprintf(stderr, "%s: memory allocation error\n", fxn_name);
                exit(EXIT_FAILURE);
            }
            for (i = 0; i < set->diverge->numMu; i++) {
                set->diverge->Mu[i] = m;
            }
            if (debug) printf("setting Mu: %f for each cluster\n", m);
        }
        else if (set->diverge->numMu != set->num_pTrees) {
            fprintf(stderr, "\tERROR: must specify a starting Mu for each cluster\n");
            exit(EXIT_FAILURE);
        }
    }
    if(debug) fprintf(stdout, "\tRead command file.\n");
}' ReadCmdfile

static void Read_List(void *list, char *str, int len, const char *option, int type) {
    const char *fxn_name = "Read_List";
    int i;
    char *ptr = str;
    char *rptr = ptr;
    int *ilist = NULL;
    double *dlist = NULL;

    if( type == INT )
        ilist = (int *) list;
    else if( type == DOUBLE )
        dlist = (double *) list;
    
    for( i=0; i<len; i++ ) {
        rptr = ptr;
        while( *ptr != ' ' && *ptr != '\0' ) ptr++;
        if( *ptr == '\0' && i<len-1 ) {
            fprintf(stderr, "%s(%s): invalid %s entry\n", fxn_name, file_name, option);
            exit(EXIT_FAILURE);
        }
        *ptr = '\0';
        if( type == DOUBLE ) {
            dlist[i] = atof(rptr);
            if(debug) fprintf(stdout, " %f", dlist[i]);
        } else if( type == INT ) {
            ilist[i] = atoi(rptr);
            if(debug) fprintf(stdout, " %d", ilist[i]);
        }
        ptr++;
    }
    if(debug) fprintf(stdout, "\n");
}' Read_List

static void Process_Change_Point_Option(settings *set, const char *str_ptr) {
    if( !strcmp(str_ptr, "parameter") ) {
        set->jump_classes = true;
        set->add_rho = true;
        if( debug ) printf("change_point: parameter\n");
    } else if( !strcmp(str_ptr, "topology") ) {
        set->jump_classes = true;
        set->add_xi = true;
        if( debug ) printf("change_point: topology\n");
    } else {
        fprintf(stderr, "ERROR settings.c::Process_Change_Point_Option: invalid value for cmdfile option change_points or jump_classes\n");
        exit(EXIT_FAILURE);
    }
}' Process_Jump_Classes_Option

static void Process_Alawadhi_Option(settings *set, const char *str_ptr) {
    if( ! strncmp( str_ptr, "parameter", 9) ) {
        set->alawadhi_parameter = true;
        if(debug && set->alawadhi_parameter) printf("alawadhi: parameter dimension J update with alawadhi\n");
    } else if( ! strncmp( str_ptr, "topology_one", 12) ) {
        set->alawadhi_topology_one = true;
        if(debug && set->alawadhi_topology_one) printf("alawadhi: topology dimension K update via AddOne/DeleteOne with alawadhi\n");
    } else if( ! strncmp( str_ptr, "topology_two", 12) ) {
        set->alawadhi_topology_two = true;
        if(debug && set->alawadhi_topology_two) printf("alawadhi: topology dimension K update via AddTwo/DeleteTwo with alawadhi\n");
    } else if( ! strncmp( str_ptr, "topology", 8) ) {
        set->alawadhi_topology = true;
        if(debug && set->alawadhi_topology) printf("alawadhi: topology dimension K update with alawadhi\n");
    }
}' Process_Alawadhi_Option

static char *strip_white(char *str) {
    char *p = NULL;
    char *result = str;
    int i = 0;
    while(isspace(str[i]) && i<(int)strlen(str)) i++;
    p = &str[i];
    i = strlen(p)-1;
    while(i>=0 && isspace(p[i])) i--;
    p[i+1] = '\0';
    memmove(result, p, strlen(p) + 1);
    return result;
}' strip_white

void Settings_Cleanup(settings *set) {
    if(set->pTree) free(set->pTree);
    if(set->init_string) free(set->init_string);
    
    if (set->recomb != NULL) {
        'if(set->recomb->rnList) free(set->recomb->rnList);
        /*
        if(set->recomb->init_alpha_string) free(set->recomb->init_alpha_string);
        if(set->recomb->init_mu_string) free(set->recomb->init_mu_string);
        if(set->recomb->init_tree_string) free(set->recomb->init_tree_string);
        if(set->recomb->init_changepoint_string) free(set->recomb->init_changepoint_string);
        */
        free(set->recomb);
    }
    else if (set->diverge != NULL) {
        if (set->diverge->Mu) free(set->diverge->Mu);
        free(set->diverge);
    }
}' Settings_Cleanup

'tree.c
#include "tree.h"

static int debug = 0;               ' set to positive integer to turn on local debugging output
static const char *file_name = "tree.c";

static double *tCL;

static void Initialize_Tree(tree *);
static void Make_Leaf_List(tree *);
static void Make_Bigger_Tree(tree *, char *, int, int, int);
static void Setup_Counts(tree *, char *);

static void Initialize_Tree(tree *tr) {
    tr->has_branches = false;
    tr->nnodes = 0;
    tr->likelihood = NULL;
    tr->node_list = NULL;
    tr->nleaves = 0;
    'tr->nbranches = 0;
    tr->root = NULL;
    tr->nchars = 0;
}' Initialize_Tree

void TreeMakeCopy(tree **ntr, const tree *otr) {
    const char *fxn_name = "TreeMakeCopy";
    node *nptr;
    double *dptr;
    boolean *bptr;
    int cnodes = 0;

    *ntr = (tree *) malloc(sizeof(tree));
    (*ntr)->nleaves = otr->nleaves;
    (*ntr)->nnodes = otr->nnodes;
    (*ntr)->tree_index = otr->tree_index;
    (*ntr)->nchars = otr->nchars;
    (*ntr)->has_branches = otr->has_branches;
    (*ntr)->is_likelihood_done = (boolean *) malloc(sizeof(boolean)*otr->nnodes);
    memcpy((*ntr)->is_likelihood_done, otr->is_likelihood_done, sizeof(boolean)*otr->nnodes);
    (*ntr)->is_likelihood_done_blank = (boolean *) malloc(sizeof(boolean)*otr->nnodes);
    memcpy((*ntr)->is_likelihood_done_blank, otr->is_likelihood_done_blank, sizeof(boolean)*otr->nnodes);
    (*ntr)->likelihood = (double *) malloc(sizeof(double)*otr->nnodes*otr->nchars);
    memcpy((*ntr)->likelihood, otr->likelihood, sizeof(double)*otr->nnodes*otr->nchars);
    (*ntr)->node_list = (node *) malloc(sizeof(node)*otr->nnodes);  ' MEMORY_ALLOCATED
    (*ntr)->root = (*ntr)->node_list;

    nptr = (*ntr)->node_list;
    dptr = (*ntr)->likelihood;
    bptr = (*ntr)->is_likelihood_done;

    Make_Subtree_From_Tree(otr, nptr, &cnodes, &dptr, &bptr, NULL);

    (*ntr)->leaf_list = (node **) malloc(sizeof(node *)*(*ntr)->nleaves);
    if( (*ntr)->leaf_list == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    Make_Leaf_List((*ntr));
    (*ntr)->root->is_root = true;
}' TreeMakeCopy

void Make_Tree(tree **tr, char *treestr, int nchars) {
    const char *fxn_name = "make_tree";
    int cnodes=0;
    char *lptr=treestr, *rptr=treestr + strlen(treestr) - 1;
    double *dptr;
    node *nptr;
    boolean *bptr;
    int i;
fprintf(stderr, "%s(%s): making tree %s of length %d\n", fxn_name, file_name, treestr, strlen(treestr));
    if( debug>5 || global_debug>5 ) printf("%s(%s): entering...\n", fxn_name, file_name);

    *tr = (tree *) malloc(sizeof(tree));
    if( *tr == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    Initialize_Tree(*tr);
    (*tr)->nchars = nchars;

    Setup_Counts(*tr, treestr);
    
    (*tr)->node_list = (node *) malloc(sizeof(node)*(*tr)->nnodes); /* memory_allocated */
    if( (*tr)->node_list == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    nptr = (*tr)->node_list;
    (*tr)->root = (*tr)->node_list;
    (*tr)->likelihood = (double *) malloc(sizeof(double)*(*tr)->nnodes*nchars); /* memory_allocated */
    if( (*tr)->likelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    
    /******************/
    
    dptr = (*tr)->likelihood;
    (*tr)->is_likelihood_done = (boolean *) malloc(sizeof(boolean)*(*tr)->nnodes);
    if( (*tr)->is_likelihood_done == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    (*tr)->is_likelihood_done_blank = (boolean *) malloc(sizeof(boolean)*(*tr)->nnodes);
    if( (*tr)->is_likelihood_done_blank == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    for( i=0; i<(*tr)->nnodes; i++ ) (*tr)->is_likelihood_done_blank[i] = false;

    bptr = (*tr)->is_likelihood_done;

    Make_Subtree(nptr, lptr, rptr, &cnodes, &dptr, &bptr, NULL, (*tr)->has_branches, nchars);
    
    if ((*tr)->has_branches) {
        for (i = 1; i < (*tr)->nnodes; i++) {
            if ((*tr)->node_list[i].branch_length == -9) {
                fprintf(stderr, "error: not all branch lengths are specified\n");
                exit(EXIT_FAILURE);
            }
        }
    
        /* so tree becomes  \_____./
         *                  /      \
         */
        
        (*tr)->root->left->branch_length += (*tr)->root->right->branch_length;
        (*tr)->root->right->branch_length = 0.00;
    }
    
    (*tr)->leaf_list = (node **) malloc(sizeof(node *)*(*tr)->nleaves);
    if( (*tr)->leaf_list == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    Make_Leaf_List((*tr));
    (*tr)->root->is_root = true;

    ' allocate memory for cached q matrix
    if(tCL == NULL) tCL = (double *) malloc(sizeof(double)*nchars);
    if(tCL == NULL) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    if( debug>5 || global_debug>5 ) printf("%s(%s): leaving after making %d==%d nodes\n", fxn_name, file_name, cnodes, (*tr)->nnodes);
/*********/
    
}' Make_Tree

static void Make_Leaf_List(tree *tr) {
    int i;
    for(i=0; i<tr->nnodes; i++) {
        if(tr->node_list[i].is_branch) continue;
        tr->leaf_list[tr->node_list[i].id] = &(tr->node_list[i]);
    }
}' Make_Leaf_List

static void Make_Bigger_Tree(tree *tr, char *treestr, int nchars, int more_nbranch, int more_nleaves) {
    int i;
    int cnodes=0, bigger = more_nbranch+more_nleaves;
    char *lptr=treestr, *rptr=treestr + strlen(treestr) - 1;
    const char *fxn_name = "make_bigger_tree";
    double *dptr;
    node *nptr;
    boolean *bptr;

    if( debug>5 || global_debug>5 ) printf("%s(%s): entering...\n", fxn_name, file_name);

    Initialize_Tree(tr);
    tr->nchars = nchars;

    Setup_Counts(tr, treestr);
    tr->node_list = (node *) malloc(sizeof(node)*(tr->nnodes+bigger));  /* memory_allocated */
    if( tr->node_list == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    nptr = tr->node_list;
    tr->root = tr->node_list;
    tr->likelihood = (double *) malloc(sizeof(double)*(tr->nnodes+bigger)*nchars);  /* memory_allocated */
    if( tr->likelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    dptr = tr->likelihood;
    tr->is_likelihood_done = (boolean *) malloc(sizeof(boolean)*(tr->nnodes+bigger));
    if( tr->is_likelihood_done == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    tr->is_likelihood_done_blank = (boolean *) malloc(sizeof(boolean)*(tr->nnodes+bigger));
    if( tr->is_likelihood_done_blank == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    for( i=0; i<(tr->nnodes+bigger); i++ ) tr->is_likelihood_done_blank[i] = false;


    bptr = tr->is_likelihood_done;

    Make_Subtree(nptr, lptr, rptr, &cnodes, &dptr, &bptr, NULL, tr->has_branches, nchars);

    if (tr->has_branches) {
        for (i = 1; i < tr->nnodes; i++) {
            if (tr->node_list[i].branch_length == -9) {
                fprintf(stderr, "error: not all branch lengths are specified\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    
    tr->leaf_list = (node **) malloc(sizeof(node *)*(tr->nleaves+more_nleaves));
    if( tr->leaf_list == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    Make_Leaf_List(tr);
    tr->root->is_root = true;

    if( debug>5 || global_debug>5 ) printf("%s(%s): leaving after making %d==%d nodes\n", fxn_name, file_name, cnodes, tr->nnodes);
}' Make_Bigger_Tree

static void Setup_Counts(tree *tr, char *treestr) {
    int nopen = 0, nclose = 0, ncommas = 0, i, nsplit = 0, len = strlen(treestr);
    const char *fxn_name = "Setup_Counts";

    if( debug>5 || global_debug>5 ) printf("%s(%s): entering...\n", fxn_name, file_name);

    ' traverse string looking for commas and parentheses
    for(i=0;i<len;i++) {
        char current = treestr[i];
        if( current == ',' && (nopen-nclose)==1) {
            if(!nsplit) nsplit = i;
            ncommas++;
        } else if( current == '(' ) nopen++;
        else if( current == ')' ) nclose++;
        else if( current == ':' ) tr->has_branches = true;
    }
    if(nopen!=nclose || ncommas>=2 || ncommas<1) {
        fprintf(stderr, "%s(%s): invalid tree structure (%s)\n", fxn_name, file_name, treestr);
        exit(EXIT_FAILURE);
    }

    ' this seg faults during realloc and malloc - not sure why - gmd
    ' root what was an unrooted tree string
    /**************************************
    if(ncommas == 2) {
        char *ttreestr = NULL;
        int add = 3;
        if (tr->has_branches) add = 5;
    
        printf("realloc treestr\n");
        treestr = (char *) realloc(treestr, sizeof(char) * (strlen(treestr) + add));
        printf("done\n");
        printf("malloc ttreestr\n");
        ttreestr = (char *) malloc(sizeof(char) * (strlen(treestr) + add));
        printf("done\n");
                
        printf("start copying\n");
        strcpy(ttreestr, "");
        ttreestr = strncat(ttreestr, treestr, nsplit);
        ttreestr = strcat(ttreestr, "(");
        ttreestr = strncat(ttreestr, treestr+nsplit+1, strlen(treestr) - nsplit);
        if(tr->has_branches) ttreestr = strcat(ttreestr, ":0)");
        else ttreestr = strcat(ttreestr, ")");
        strcpy(treestr, ttreestr);
        free(ttreestr);
        nopen++;
    }
    *******************/
    tr->nleaves = nopen+1;  ' (0,1,(2,3)) has 3 leaves; ((0,1),(2,3)) has 4 leaves
    tr->nnodes = 2*tr->nleaves - 1;

    if( debug>5 || global_debug>5 ) printf("%s(%s): found tree with branches? (%d).\n", fxn_name, file_name, tr->has_branches);
}' Setup_Counts

boolean SameTrees(tree *tr1, tree *tr2, boolean branch) {
    char tmp1[MAX_TREE_STRING], tmp2[MAX_TREE_STRING]; ' dangerous
    tmp1[0] = '\0'; tmp2[0] = '\0';
    toString(tmp1, tr1->root, branch);
    toString(tmp2, tr2->root, branch);
    if (strcmp(tmp1, tmp2) == 0) return true;
    return false;
}' SameTrees


void TreeDelete(tree *tr) {
    if(tr->node_list) free(tr->node_list);
    if(tr->likelihood) free(tr->likelihood);
    if(tr->is_likelihood_done) free(tr->is_likelihood_done);
    if(tr->is_likelihood_done_blank) free(tr->is_likelihood_done_blank);
    if(tr->leaf_list) free(tr->leaf_list);
}' TreeDelete


/** original cbrother loglikelihood function ****************

' idea: to allow different models, use a struct with a function pointer and a void * argument to a struct of needed variables; update is either alpha (recomb) or lambda (diverge)
' you can pass in a structure, but no easy way to pass in structure without copying it (unless you set a pointer to the struct each time you want to use it)
'update = alpha or lambda
double loglikelihood(tree *tr, seqdata *sqd, stddata *std, const int *counts, const double avg_brlen, const double update, const double *pi, const int model){
    int numsites = sqd->lenunique;
    int numleaves = sqd->ntaxa;
    int nchars = sqd->num_chars;
    int i,j,k,l;
    double ll = 0.0, l = 0;
    node *cnode = NULL;
    double stock;
    int **data = std->data;
    double *lvector, *rvector;
    const char *fxn_name = "loglikelihood";

    update_matrix_ihky85(avg_brlen, update, pi);

'  if(debug) printf("%s(%s): entering...\n", fxn_name, file_name);
    for(i=0; i<numsites; i++) {
'  if(debug) printf("%s(%s): processing site %d\n", fxn_name, file_name, i);
        for(j=0;j<tr->nnodes;j++) tr->is_likelihood_done[j] = false;
        if( counts[i] != 0 ) {
'printf("%s(%s): %d %d %d %d %d %d %d\n", fxn_name, file_name, data[i][0], data[i][1], data[i][2], data[i][3], data[i][4], data[i][5], data[i][6]);
            for(j=0; j<numleaves; j++) {    ' here, let's calculate the likelihood of a single word
                cnode = tr->leaf_list[j];
'  if(debug) printf("%s(%s): processing leaf node %d with data %d\n", fxn_name, file_name, cnode->uid, data[i][j]);
                        
                if( (data[i][j] >= 0) ) {
                    for(k=0; k<nchars; k++)
                        cnode->clikelihood[k] = cached_qmatrix[k][data[i][j]];
                } else {
                    for(k=0; k<nchars; k++) cNode->clikelihood[k] = 1;
                }
'if(debug) printf("%s(%s): node %d clikelihood %18.16lf %18.16lf %18.16lf %18.16lf\n", fxn_name, file_name, cNode->uid, cNode->clikelihood[0], cNode->clikelihood[1], cNode->clikelihood[2], cNode->clikelihood[3]);
                *(cNode->is_likelihood_done) = true;
                while( (cNode != tr->root) && *(cNode->brother->is_likelihood_done) ) { ' Propogate up tree until we hit the root or an unfinished section
                    cNode = cNode->up;
                    lvector = cNode->left->clikelihood;
                    rvector = cNode->right->clikelihood;
'  if(debug) printf("Moving up to branch %d\n", cNode->uid);
'if(debug) printf("%s(%s): node %d clikelihood %18.16lf %18.16lf %18.16lf %18.16lf\n", fxn_name, file_name, lChild->uid, lChild->clikelihood[0], lChild->clikelihood[1], lChild->clikelihood[2], lChild->clikelihood[3]);
                    if (cNode != tr->root) {
'  if(debug) printf("%s(%s): processing branch %d\n", fxn_name, file_name, cNode->uid);
                        for(k=0; k<nchars; k++) tCL[k] = lvector[k] * rvector[k];
                        if(cNode->up == tr->root && *(cNode->brother->is_likelihood_done)) {
                            for(k=0; k<nchars; k++) cNode->clikelihood[k] = tCL[k];
                        } else {
                            for(k=0; k<nchars; k++) {
                                cNode->clikelihood[k] = 0;
                                for(l=0; l<nchars; l++) cNode->clikelihood[k] += cached_qmatrix[k][l] * tCL[l];
                            }
                        }
                    } else for(k=0; k<nchars; k++) cNode->clikelihood[k] = lvector[k] * rvector[k];
'  if(debug) printf("%s(%s): processing root %d\n", fxn_name, file_name, cNode->uid);
                    *(cNode->is_likelihood_done) = true;
'if(debug) printf("%s(%s): node %d clikelihood %18.16lf %18.16lf %18.16lf %18.16lf\n", fxn_name, file_name, cNode->uid, cNode->clikelihood[0], cNode->clikelihood[1], cNode->clikelihood[2], cNode->clikelihood[3]);
                }
            }
        ' By this point all conditional likelihoods are calculated and cNode = Root. The conditionals
        ' for the Root are still in tCL
        ' Posterior likelihood = SUM( conditional * PI_i ) of all i in span.
            L = 0.0;
            for(k=0; k<nchars; k++) L += pi[k] * cNode->clikelihood[k];
'  if(debug) printf("%s(%s): counts for this site %d\n", fxn_name, file_name, counts[i]);
            LL += (double) log(L) * counts[i];
'if(debug) printf("%s(%s): loglike thru position %d so far %lf = %lf (log(%lf)) * %d\n", fxn_name, file_name, i, LL, log(L), L, counts[i]);
        }
    }
    return LL;
}' LogLikelihood

****************/


' IDEA: to allow different models, use a struct with a function pointer and a VOID * argument to a struct of needed variables; update is either alpha (recomb) or lambda (DIVERGE)
' you can pass in a structure, but no easy way to pass in structure without copying it (unless you set a pointer to the struct each time you want to use it)
'update = alpha or lambda

' wrapper functions
'double TreeLogLikelihood(tree *tr, const sampler *smp, qmatrix *qmat, const int *counts, const double avg_brlen, boolean update_each_branch) {
Public Function TreeLogLikelihood(tr As tree, smp As sampler, qmat As qmatrix, counts As Long, avg_brlen As Double, update_each_branch As Byte) As Double
    'return CalcLikelihood(tr, smp->sqd, qmat, counts, avg_brlen, update_each_branch, -1, true);
    TreeLogLikelihood = CalcLikelihood(tr, smp.sqd, qmat, counts, avg_brlen, update_each_branch, -1, 1)
End Function '}

double SiteLikelihood(tree *tr, const seqdata *sqd, qmatrix *qmat, boolean update_each_branch, int site) {
    return CalcLikelihood(tr, sqd, qmat, NULL, -1, update_each_branch, site, false);
}


' SET PARAMETER                    TO
' Log = true           return loglikelihood
' Log = false              return likelihood
' site >= 0            compute for given site
' site < 0             compute for all unique sites
' update_each_branch = true    update matrix for each new branch length (DIVERGE)
' update_each_branch = false   update matrix once (for avg. branch length, RECOMB)

double CalcLikelihood(tree *tr, const seqdata *sqd, qmatrix *qmat, const int *counts, const double avg_brlen, boolean update_each_branch, int site, boolean Log) {
Public Function CalcLikelihood(tr As tree, sqd As seqdata, qmat As qmatrix, counts() As Long, avg_brlen As Double, update_each_branch As Byte, site As Long, Log As Byte) As Double


    'const char *fxn_name = "CalcLikelihood";
    Dim fxn_name As String, numSites As Long, numLeaves As Long, nchars As Long, LL As Double, l As Double
    fxn_name = "CalcLikelihood"
    
    'char tree_vis[MAX_TREE_STRING];
    Dim tree_vis As String
        
    Dim cNode As node
        
    'int numSites = site>=0 ? 1 : sqd->lenunique;
    If site >= 0 Then
        numSites = 1
    Else
        numSites = sqd.lenunique
    End If
    
    'int numLeaves = tr->nleaves; 'sqd->ntaxa; ntaxa only works if enumerated trees are passed in
    numLeaves = tr.nleaves 'sqd->ntaxa; ntaxa only works if enumerated trees are passed in
    
    'nchars = sqd->num_chars;
    nchars = sqd.num_chars
    'double LL = 0.0, L = 0.0;
    LL = 0#: l = 0#
    cNode = Null
    
    'int **data = sqd->data;
    Dim data As Long
    data = sqd.data
    
    'double *lvector, *rvector;
    Dim lvector As Double, rvector As Double
    'int i = 0, j,k,l;
    Dim i As Long, j As Long, k As Long, l As Long
    i = 0
    
    'double *pi = qmat->pi;
    Dim pi As Double
    pi = qmat.pi
    
    'boolean display = false;
    Dim display As Byte
    display = 0
    
    'tree_vis[0] = '\0';
    tree_vis = "\0"
    
    'toString(tree_vis, tr->root, false);
    Call toString(tree_vis, tr.root, 0)
'fprintf(stderr, "CalcLikelihood: %f %f %s", qmat->v[0], avg_brlen, tree_vis);

    'if (display) PrintTreeInfo(tr);
    
    'if (site >= 0) {  ' only calculate for one site
    If site >= 0 Then  ' only calculate for one site
        i = site
        numSites = site + 1
    End If

    ' For this to work for updating each branch, a vector of matrices needs to be stored.
    'if (!update_each_branch) {
    If update_each_branch = 0 Then
        'qmat->Matrix_Update_Cache(qmat, avg_brlen);
        Call Matrix_Update_Cache(qmat, avg_brlen)
    End If
    
    For i = 0 To numSites - 1
'fprintf(stderr, "%d ", counts[i]);
        If site >= 0 Then i = site
        'memcpy(tr->is_likelihood_done, tr->is_likelihood_done_blank, sizeof(boolean)*tr->nnodes);
        For X = 0 To tr.nnodes
            tr.is_likelihood_done_blank(X) = tr.is_likelihood_done(X)
        Next X
        'if( counts == NULL || counts[i] != 0 ) {  ' do we need this?
        If counts(0) = 0 Or counts(i) <> 0 Then ' do we need this?
            
            For j = 0 To numLeaves - 1 ' Here, let's calculate the likelihood of a single word
                cNode = tr.leaf_list(j)
                    'if (display) printf("%s(%s): processing leaf %d with data %d\n", fxn_name, file_name, cNode->id, data[i][j]);
                'if (update_each_branch) qmat->Matrix_Update_Cache(qmat, cNode->branch_length);
                If update_each_branch = 1 Then Call Matrix_Update_Cache(qmat, cNode.branch_length)
                if( (data[i][j] >= 0) ) {
                    memcpy(cNode->clikelihood, qmat->cached_qmatrix[data[i][j]], sizeof(double)*nchars);
'                  for(k=0; k<nchars; k++) cNode->clikelihood[k] = qmat->cached_qmatrix[k][data[i][j]];
                    for(k=0; k<nchars; k++) if( cNode->clikelihood[k] > 1.0 ) fprintf(stderr, "%f!!!!!\n", cNode->clikelihood[k]);
                } else {
                    for(k=0; k<nchars; k++) cNode->clikelihood[k] = 1.0;
                }
                *(cNode->is_likelihood_done) = true;
                while( (cNode != tr->root) && *(cNode->brother->is_likelihood_done) ) { ' Propogate up tree until we hit the root or an unfinished section
                    
                    cNode = cNode->up;

                    lvector = cNode->left->clikelihood;
                    rvector = cNode->right->clikelihood;

                    if (display) {
                      /**/
                      printf("lvector for node %d: ", cNode->left->uid);
                      for (k = 0; k < nchars; k++) printf("%f ", lvector[k]);
                      printf("\n");
                        
                      printf("rvector for node %d: ", cNode->right->uid);
                      for (k = 0; k < nchars; k++) printf("%f ", rvector[k]);
                      printf("\n");
                      /**/
                    }
                    
                    if (cNode != tr->root) {
                        
                        for (k = 0; k < nchars; k++) {
                            tCL[k] = lvector[k]*rvector[k];
'fprintf(stderr, "tCL[%d]: %f\n", k, tCL[k]);
                        }
                        if (display) printf("\nclike for node %d: ", cNode->uid);
                        ' handle avg brlen of 0 from root to right child
                        ' How come this doesn't apply if update_each_branch is true?  update_each_branch will have the branch length set to 0 explicitly
                        if( !update_each_branch && cNode->up == tr->root && *(cNode->brother->is_likelihood_done) ) {
                            for( k=0; k<nchars; k++) cNode->clikelihood[k] = tCL[k];
                        } else {

                            if (update_each_branch) qmat->Matrix_Update_Cache(qmat, cNode->branch_length);
                            for (k = 0; k < nchars; k++) {
                                cNode->clikelihood[k] = 0;
                                for (l=0; l<nchars; l++) cNode->clikelihood[k] += qmat->cached_qmatrix[l][k]*tCL[l];
                                if (display) printf("%f ", cNode->clikelihood[k]);
                            }
                            if (display) printf("\n\n");
                        }

                    } else {
                        for(k=0; k<nchars; k++) cNode->clikelihood[k] = lvector[k] * rvector[k];
                    }
                    *(cNode->is_likelihood_done) = true;
                }
            Next j
        ' By this point all conditional likelihoods are calculated and cNode = Root. The conditionals
        ' for the Root are still in tCL
        ' Posterior likelihood = SUM( conditional * PI_i ) of all i in span.
            L = 0.0;
            for(k=0; k<nchars; k++) {
                if (display) printf("%f x %f\n", pi[k], cNode->clikelihood[k]);
                L += pi[k] * cNode->clikelihood[k];
'fprintf(stderr, "L = %f += %f*%f\n", L, pi[k], cNode->clikelihood[k]);
            }
            ' Oh no, but counts might be NULL!
'fprintf(stderr, "LL = %f = log(%f)*%d\n", LL+log(L)*counts[i], L, counts[i]);
            if (Log) LL += (double) log(L) * counts[i];
            'else Likelihood *= L;
        End If
    Next i

' Note: You should probably output a warning message in the event of L > 1
    'fprintf(stderr, " %f\n", L);
    if( L - 1.0 > tolerance ) {
        fprintf(stderr, "WARNING %f!!\n", L);
        exit(EXIT_FAILURE);
        L = 1;
        LL = 0;
    }

    if (Log) return LL;

    if (display) printf("L: %f\n", L);
' But this doesn't make sense.  If you calculate over multiple sites L will only contain data on the last site.  There is no internal mechanism to catch this outcome.
    return L;
End Function ' CalcLikelihood

void Simulate_Position(tree *tr, qmatrix *qmt, const double *pi, double mu, rngen *rng) {
    int i;
    node *cn = NULL;
    double rn, cum=0.0;

    cn = tr->root;
    rn = rng->nextStandardUniform(rng);
    for( i=0; i<tr->nchars; i++ ) {
        if( rn < cum + pi[i] ) {
            cn->state = i;
            break;
        }
        cum += pi[i];
    }
    if( cn->is_branch ) {
        Simulate_Down_Branch(cn->left, cn, mu, qmt, rng);
        Simulate_Down_Branch(cn->right, cn, mu, qmt, rng);
    }
}' Simulate_Position


/*
' I really don't know what this does, or if it does whatever that is correctly (I think not)
void YankHardOnTree(const node *parent, node *child){
    node *t1;
    if( child->up == parent ) return;
    t1 = child->up;
    child->up = parent;
    if( child->right == parent ) child->right = t1;
    else child->left = t1;
    if( child->right != NULL ) {
        if( child->right == t1 ) YankHardOnTree(child, child->right);
        else YankHardOnTree(child, child->right);
    }
    if( child->left != NULL ) {
        if( child->left == t1 ) YankHardOnTree(child, child->Left);
        else YankHardOnTree(child, child->Left);
    }
    return;
}' YankHardOnTree
*/


/*********** OLD VERSION FOR ONLY 1 QUERY *******************
void EnumerateLastTaxon(tree *trlist, const tree *tr) {
    int ntrees = 2 * tr->nleaves - 3;     ' The number of trees when parental tree fixed and one query sequence
    int i, j;
    char temp[200]="";
    const char *fxn_name = "EnumerateLastTaxon";

'  if(debug) printf("%s(%s): entering...\n", fxn_name, file_name);

    toString(temp, tr->root);
    for(i=0; i<ntrees; i++){  ' Fill each tree with the user-specified start tree (relating parental sequences)
        Make_Bigger_Tree(&trlist[i], temp, tr->nchars, 1, 1);
    }
    j = 0;
    for(i=0; i<tr->nnodes; i++) {   ' Make all possible trees formed by growing query out of all possible branches
        if(!tr->node_list[i].is_branch) {
            continue;
        } else if(tr->node_list[i].uid == tr->root->uid) {
            GrowLeaf(&trlist[j++], i, false);
            Balance_Tree(trlist[j-1]);
            temp[0] = '\0';
            toString(temp, trlist[j-1].root);
            if(debug) printf("Tree %s\n", temp);
        } else {
            GrowLeaf(&trlist[j++], i, false);   ' Grow query out of right descendent branch
            Balance_Tree(trlist[j-1]);
            temp[0] = '\0';
            toString(temp, trlist[j-1].root);
            if(debug) printf("Tree %s\n", temp);
            GrowLeaf(&trlist[j++], i, true);    ' Grow query out of left descendent branch
            Balance_Tree(trlist[j-1]);
            temp[0] = '\0';
            toString(temp, trlist[j-1].root);
            if(debug) printf("Tree %s\n", temp);
        }
    }
}' EnumerateLastTaxon
********************************/


/******************/
void EnumerateLastTaxon(tree **trlist, const tree *tr) {
'  const char *fxn_name = "EnumerateLastTaxon";
    int ntrees = 2 * tr->nleaves - 3;     ' The number of trees when parental tree fixed and one query sequence
    int i, j;
    char temp[200]="";  ' BUG: hard-coded length cannot handle very large tree
    int trees_added = 0;

    if( (*trlist) == NULL ) {
        *trlist = (tree *) malloc(sizeof(tree)*ntrees);
    }

    toString(temp, tr->root, false);
    
    for(i=0; i<ntrees; i++){  ' Fill each tree with the user-specified start tree (relating parental sequences)
        Make_Bigger_Tree(&((*trlist)[i]), temp, tr->nchars, 1, 1);
    }
    
    j = 0;
    for(i=0; i<tr->nnodes; i++) {   ' Make all possible trees formed by growing query out of all possible branches

        if (trees_added >= ntrees) break;

        if ( (j >= ntrees) || ((*trlist)[j].node_list == NULL) ) continue; ' added by gmd
    
        ' note: tr and trees in trlist may have differing internal nodes even though same structure
        
        'if(!tr->node_list[i].is_branch) continue;
        if (!(*trlist)[j].node_list[i].is_branch) continue;  ' gmd: use trlist instead why?  Is this a BUG? - ksd
        
        'else if(tr->node_list[i].uid == tr->root->uid) {
        else if((*trlist)[j].node_list[i].uid == (*trlist)[j].root->uid) {  ' gmd: use trlist instead why?  Is this a BUG? - ksd
            GrowLeaf(&((*trlist)[j++]), i, false);
            Balance_Tree(&((*trlist)[j-1]));
'          temp[0] = '\0';
'          toString(temp, trlist[j-1].root, false);
        }
        else {
            GrowLeaf(&((*trlist)[j++]), i, true);       ' Grow query out of left descendent branch
            Balance_Tree(&((*trlist)[j-1]));
'          temp[0] = '\0';
'          toString(temp, trlist[j-1].root, false);
            GrowLeaf(&((*trlist)[j++]), i, false);   ' Grow query out of right descendent branch
            Balance_Tree(&((*trlist)[j-1]));
'          temp[0] = '\0';
'          toString(temp, trlist[j-1].root, false);
        }
    }

}' EnumerateLastTaxon
/*************/


void GrowLeaf(tree *tr, const int n, const boolean left) {
    node *B, *Q, *O=&(tr->node_list[n]);
    double *dptr;
    boolean *bptr;
    int nnodes = tr->nnodes;
    const char *fxn_name = "GrowLeaf";

    if( debug>5 || global_debug>5 ) printf("%s(%s): entering...\n", fxn_name, file_name);

    dptr = tr->likelihood;
    dptr += nnodes*tr->nchars;

    bptr = tr->is_likelihood_done;
    bptr += nnodes;

    tr->nnodes += 2;
    tr->nleaves++;
'  tr->nbranches++;

    B = &(tr->node_list[nnodes]);
    Q = &(tr->node_list[nnodes+1]);

    if(left) {
        Make_Node(B, O, Q, O->left, nnodes++, -1, -9, dptr++, bptr++);
        O->left = B;
        B->brother = O->right;
        O->right->brother = B;
    } else {
        Make_Node(B, O, Q, O->right, nnodes++, -1, -9, dptr++, bptr++);
        O->right = B;
        B->brother = O->left;
        O->left->brother = B;
    }
    B->right->up = B;
    Make_Node(Q, B, NULL, NULL, nnodes, tr->nleaves-1, -9, dptr++, bptr++);
    tr->leaf_list[tr->nleaves-1] = Q;
    B->left->brother = B->right;
    B->right->brother = B->left;
}' GrowLeaf

int Number_Parental_Trees(tree *tr) {
    return(2*tr->nleaves - 3);
}' Number_Parental_Trees

int Number_All_Trees(tree *tr, const int numQueries) {
        int total = 1, numLeaves = tr->nleaves, i = 1;
        for (i = 0; i < numQueries; i++) {
            total = total * (2 * (numLeaves + i) - 3);
    }
    return total;
}' Number_All_Trees

boolean Monophyletic(tree *tr, int *clades) {
    int i, group_size;
    Cladify(tr->root, tr->nleaves, 2, clades);
    Cladify_Up(tr->root, tr->nleaves, 2, clades);
    for( i=0, group_size=0; i<tr->nleaves; i++ ) {
        if( clades[i] ) group_size++;
    }
    return Clade_Ancestor_Exists(tr->root, tr->nleaves, 1, group_size);
}' Monophyletic

void Balance_Tree(tree *tr) {
    Balance(tr->root);
}' Balance_Tree


void PrintTreeInfo(tree *tr) {
        int i = 0;
        char tmp[500];
        tmp[0] = '\0';
        toString(tmp, tr->root, true);
        printf("tree: %s\n", tmp);
        printf("num leaves: %d, num nodes: %d, num nchars: %d\n", tr->nleaves, tr->nnodes, tr->nchars);
        printf("node list:\n");
        for (i  = 0; i < tr->nnodes; i++) {
            printf("  %d: %d, %d,   %f\n", i, tr->node_list[i].uid, tr->node_list[i].id, tr->node_list[i].branch_length);
        }
        printf("\n");
        PrintNodesDFS(*(tr->root));
        printf("\n\n");
}' PrintTreeInfo

void PrintLeafList(tree *tr) {
    int i;
        for (i  = 0; i < tr->nnodes; i++) {
        if (tr->node_list[i].is_branch) continue;
        printf("%d ",tr->node_list[i].uid);
    }
    printf("\n");

}' PrintLeafList
    
    
    '                               normal(0, sdT2),             normal(0, sdT1),    unif
       'public Tree JointBranchAndTopology(AbstractDistribution nglobal, AbstractDistribution nlocal, AbstractDistribution unif, double pmix) {

' treestr is string for current tree; new tree will be contained in tr
void JointBranchAndTopology(char *treestr, tree *tr, settings *set, int nchars, double pmix) {
        
    'bTree pTree = new bTree(this);
    int i;
    node *selected = NULL;
    node *parent = NULL;
    node *moving = NULL;
    node *aunt = NULL;
    boolean SelIsLeftKid;

    Make_Tree(&tr, treestr, nchars);
    
    if( gsl_rng_uniform(set->rng->rengine) < pmix ) {   ' unif.nextDouble() < pmix
        for(i=1; i<tr->nleaves; i++)
            tr->leaf_list[i]->branch_length = fabs(tr->leaf_list[i]->branch_length + gsl_ran_gaussian(set->rng->rengine, set->sdT2)); 'nglobal.nextDouble());
            tr->node_list[1].branch_length = fabs(tr->node_list[1].branch_length + gsl_ran_gaussian(set->rng->rengine, set->sdT2));  ' nglobal.nextDouble();
        '  for(i=2; i<tr->nbranches; i++)
            for(i=2; i<tr->nnodes; i++) {
                if (!tr->node_list[i].is_branch) { ' added by gmd
                    continue;
                }
                selected = &tr->node_list[i];
                selected->branch_length += gsl_ran_gaussian(set->rng->rengine, set->sdT2); 'nglobal.nextDouble
                
                if( selected->branch_length < 0 ) {
                    selected->branch_length = - selected->branch_length;
                    parent = selected->up;
                    moving = NULL;
                    SelIsLeftKid = false;
                    if( selected == parent->left )
                        SelIsLeftKid = true;
                    if( gsl_rng_uniform(set->rng->rengine) < 0.5 ) { ' unif.nextDouble < 0.5
                        moving = selected->left;
                        if( SelIsLeftKid ) {
                            ' reattach as right
                            aunt = parent->right;
                            aunt->up = selected;
                            selected->left = aunt;
                            parent->right = moving;
                            moving->up = parent;
                        } else {
                            ' reattach as left;
                            aunt = parent->left;
                            aunt->up = selected;
                            selected->left = aunt;
                            parent->left = moving;
                            moving->up = parent;
                        }
                    } else {
                        moving = selected->right;
                        if( SelIsLeftKid ) {
                            aunt = parent->right;
                            aunt->up = selected;
                            selected->right = aunt;
                            parent->right = moving;
                            moving->up = parent;
                        } else {
                            aunt = parent->left;
                            aunt->up = selected;
                            selected->right = aunt;
                            parent->left = moving;
                            moving->up = parent;
                        }
                    }
                    Balance_Tree(tr);
                }
            }
    } else {
        int r = 0;
        r = (int) (gsl_rng_uniform(set->rng->rengine)*(2*tr->nleaves - 3)); ' unif.nextDouble * 2*numLeaves - 3
        if( r == 0 )
            tr->root->right->branch_length = fabs(tr->root->right->branch_length +  gsl_ran_gaussian(set->rng->rengine, set->sdT1)/*nlocal.nextDouble()*/ );
        else if( r < tr->nleaves )
            tr->leaf_list[r]->branch_length = fabs(tr->leaf_list[r]->branch_length + gsl_ran_gaussian(set->rng->rengine, set->sdT1)/*nlocal.nextDouble()*/ );
        else { ' update an internal branch
            r -= tr->nleaves;
            r += 2;
            selected = &tr->node_list[r];
            selected->branch_length = selected->branch_length + gsl_ran_gaussian(set->rng->rengine, set->sdT1)/*nlocal.nextDouble()*/;
            if( selected->branch_length < 0 ) {
                selected->branch_length = - selected->branch_length;
                parent = selected->up;
                moving = NULL;
                SelIsLeftKid = false;
                if( selected == parent->left )
                    SelIsLeftKid = true;
                if( gsl_rng_uniform(set->rng->rengine) < 0.5 ) {  ' unif.nextDouble() < 0.5
                    moving = selected->left;
                    if( SelIsLeftKid ) {
                        ' reattach as right
                        aunt = parent->right;
                        aunt->up = selected;
                        selected->left = aunt;
                        parent->right = moving;
                        moving->up = parent;
                    } else {
                        ' reattach as left;
                        aunt = parent->left;
                        aunt->up = selected;
                        selected->left = aunt;
                        parent->left = moving;
                        moving->up = parent;
                    }
                } else {
                    moving = selected->right;
                    if( SelIsLeftKid ) {
                        aunt = parent->right;
                        aunt->up = selected;
                        selected->right = aunt;
                        parent->right = moving;
                        moving->up = parent;
                    } else {
                        aunt = parent->left;
                        aunt->up = selected;
                        selected->right = aunt;
                        parent->left = moving;
                        moving->up = parent;
                    }
                }
                Balance_Tree(tr);
            }
        }
    }
        
}' JointBranchAndTopology


' Can's version is commented out; did not include 1st leaf branch_length
double SumOfBranchLengths(tree *tr) {
    int i;
    double rtn = 0;
/*
    for(i=1; i<tr->nleaves; i++)
        rtn += tr->leaf_list[i]->branch_length;
    for(i=1; i<tr->nbranches; i++)
            rtn += tr->node_list[i].branch_length;
*/

    for (i = 1; i < tr->nnodes; i++) rtn += tr->node_list[i].branch_length;
    return rtn;
    
}' SumOfBranchLengths


void PrintBrothers(tree *tr) {
    node cNode;
    int i;
    for (i = 0; i < tr->nnodes; i++) {
        cNode = tr->node_list[i];
        if (cNode.brother != NULL) {
            printf("brother of %d is %d\n", cNode.uid, cNode.brother->uid);
        }
        else printf("brother of %d is NULL\n", cNode.uid);
    }
}
'tree_vector.c
#include "tree_vector.h"

' FILE-WIDE VARIABLE DEFINITIONS
static int debug = 0;                                                   ' Set to positive integer to turn debugging on in this file only
static const char *file_name = "tree_vector.c";
static void *top_prior;                                                 ' Object for computing prior on vector \tau (handles elt depencies)

' LOCAL FUNCTION PRE-DECLARATIONS

' External linkage via function pointer: i.e. there can be multiple versions
' of these functions for different MCMC transition kernels.
static tree *Draw_Initial_Tree(const topology_vector *, settings *);                ' Initialize MCMC: select tree from initial distribution
static void TopologyAddOne(topology_vector *, partition_list *, sampler *);         ' Move: Add one topology change point
static void TopologyAddTwo(topology_vector *, partition_list *, sampler *);         ' Move: Add two topology change points
static void TopologyDeleteOne(topology_vector *, partition_list *, sampler *);          ' Move: Delete one topology change point
static void TopologyDeleteTwo(topology_vector *, partition_list *, sampler *);          ' Move: Delete two topology change points
static void Update_Topologies(topology_vector *, partition_list *, sampler *, boolean);     ' Move: change topology

' Setup functions:
static void Make_Tree_List(topology_vector *, const sampler *);                 ' Setup list of possible trees based on fixed parental tree
static void TopologyPriorMake(topology_vector *, const sampler *);              ' Setup prior; makes top_prior if DNE and necessary/relevant
static void TopologyGModelPriorMake(const topology_vector *, const sampler *);          ' Make prior object for Gmodel (case where top_prior necessary)
static double UniformMarkovPriorRatio(const topology_vector *, int, ...);           ' Compute prior ratio on \tau (modified Vladimir style)
static double MonophyleticMarkovPriorRatio(const topology_vector *, int, ...);          ' Compute prior ratio on \tau (Fang style)

' Local utility (help move functions): actually first three are externally linked
' via function pointer, but there is no need for this at the moment.
static tree *Propose_New_Tree(const topology_vector *, settings *, tree *, tree *, tree *);             ' Propose a new tree given left, current, and right tree
static double Reverse_Add_One_Log_Proposal_Probability(topology_vector *, const partition_list *, int, boolean, tree *);' Move auxiliary: Probability of reverse move for add one
static double Reverse_Add_Two_Log_Proposal_Probability(topology_vector *, const partition_list *, int, tree *);     ' Move auxiliary: Probability of reverse move for add two
static int ProposedTreeStructureWithDeleteOne(topology_vector *, const partition_list *, int, boolean);         ' Move auxiliary: Prepare proposed tree structure during delete one
static int ProposedTreeStructureWithDeleteTwo(topology_vector *, const partition_list *, int, int);             ' Move auxiliary: Prepare proposed tree structure during delete two

static void AddOneAccept(sampler *, partition_list *, partition *, partition *, topology_vector *, int, int, int, int, boolean, tree *, tree *, double *, double, double, double);
static void AddTwoAccept(sampler *, partition_list *, topology_vector *, int, int, partition *, int, int, tree *, double, double, double, double *, double, int);
static void DeleteOneAccept(sampler *, partition_list *, topology_vector *, int, int, int, boolean, tree *, double *, double);

void TreeVectorMake(topology_vector **tv, const sampler *smp) {
    const char *fxn_name = "TreeVectorMake";

    *tv = (topology_vector *) malloc(sizeof(topology_vector));

    if( *tv == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    (*tv)->gmodel = smp->set->gmodel;
    (*tv)->proposed_trees = NULL;
    (*tv)->current_trees = NULL;

    ' Create the start_tree structure (fixed parental tree)
    Make_Tree(&((*tv)->start_tree), smp->set->pTree[0], smp->sqd->num_chars);

    ' Create the tree list
    Make_Tree_List(*tv, smp);

    ' Setup prior
    TopologyPriorMake(*tv, smp);

    (*tv)->Draw_Initial_Tree = &Draw_Initial_Tree;
    (*tv)->Propose_New_Tree = &Propose_New_Tree;
    (*tv)->Add_One = &TopologyAddOne;
    (*tv)->Reverse_Add_One_Log_Proposal_Probability = &Reverse_Add_One_Log_Proposal_Probability;
    (*tv)->Add_Two = &TopologyAddTwo;
    (*tv)->Reverse_Add_Two_Log_Proposal_Probability = &Reverse_Add_Two_Log_Proposal_Probability;
    (*tv)->Delete_One = &TopologyDeleteOne;
    (*tv)->Delete_Two = &TopologyDeleteTwo;
    (*tv)->Update_Topologies = &Update_Topologies;
}' TreeVectorMake

void TreeVectorMakeCopy(topology_vector **ntv, const topology_vector *otv, int dim) {
    *ntv = (topology_vector *) malloc(sizeof(topology_vector));
    (*ntv)->numTrees = otv->numTrees;
    (*ntv)->log_numTrees = otv->log_numTrees;
    (*ntv)->log_numTrees_minus1 = otv->log_numTrees_minus1;
    (*ntv)->log_numTrees_minus2 = otv->log_numTrees_minus2;
    (*ntv)->log_prior_prob = otv->log_prior_prob;
    (*ntv)->top_prior = otv->top_prior;
    (*ntv)->gmodel = otv->gmodel;
    (*ntv)->Draw_Initial_Tree = otv->Draw_Initial_Tree;
    (*ntv)->Propose_New_Tree = otv->Propose_New_Tree;
    (*ntv)->Add_One = otv->Add_One;
    (*ntv)->Reverse_Add_One_Log_Proposal_Probability = otv->Reverse_Add_One_Log_Proposal_Probability;
    (*ntv)->Add_Two = otv->Add_Two;
    (*ntv)->Reverse_Add_Two_Log_Proposal_Probability = otv->Reverse_Add_Two_Log_Proposal_Probability;
    (*ntv)->Delete_One = otv->Delete_One;
    (*ntv)->Delete_Two = otv->Delete_Two;
    (*ntv)->Update_Topologies = otv->Update_Topologies;
    (*ntv)->Log_Prior_Ratio = otv->Log_Prior_Ratio;
    (*ntv)->start_tree = NULL;  ' Doesn't need start tree, since obtaining tree_list from otv
    (*ntv)->tree_list = otv->tree_list; ' Just copy pointer
    (*ntv)->current_trees = (tree **) malloc(sizeof(tree *)*dim);
    memcpy((*ntv)->current_trees, otv->current_trees, sizeof(tree *)*dim);
    (*ntv)->proposed_trees = (tree **) malloc(sizeof(tree *)*dim);
    memcpy((*ntv)->proposed_trees, otv->proposed_trees, sizeof(tree *)*dim);
}' TreeVectorMakeCopy

'static tree *Draw_Initial_Tree(const topology_vector *tv, settings *set) {
Public Function Draw_Initial_Tree(tv As topology_vector, setx As settings) As tree
    'return tv->tree_list[ (int) (setx->rng->nextStandardUniform(set->rng) * tv->numTrees) ];
    Draw_Initial_Tree = tv.tree_list(CLng(setx.rng.nextStandardUniform(setx.rng) * tv.numTrees))
}' Draw_Initial_Tree

void TreeVectorInitialize(topology_vector *tv, const partition_list *pl) {
    const char *fxn_name = "TreeVectorInitialize";
    int i=0, j=0;

    if( tv->gmodel || debug>3 || global_debug>3 ) {
        if( tv->current_trees ) free(tv->current_trees);
        tv->current_trees = (tree **) malloc(sizeof(tree *)*(pl->topology_changes + 1));    ' VERIFIED 2/11/05
        if( tv->current_trees == NULL ) {
            fprintf(stderr, "%s: memory allocation error\n", fxn_name);
            exit(EXIT_FAILURE);
        }

        while( i < pl->npartitions ) {
            tv->current_trees[j++] = pl->part[i++]->ctree;
            while( i<pl->npartitions && !pl->part[i]->topchange ) i++;
        }
    }
    tv->log_prior_prob = 0.0;
    tv->log_prior_prob = tv->Log_Prior_Ratio(tv, pl->topology_changes + 1, false);

}' TreeVectorInitialize

static void Make_Tree_List(topology_vector *tv, const sampler *smp) {
    const char *fxn_name = "Make_Tree_List";
    int i;
    int queries_added = 0;                      ' How many queries have been added to fixed parental tree
    int numTreesPerEnumeration = 0;                 ' How many possible trees given the current tree?
    int numCurrentLeaves = 0;                   ' The current number of leaves in the tree
    tree *tree_list = NULL;                     ' Will contain our list of possible topologies
    tree *new_tree_list = NULL;                 ' Memory management
    tree *tree_ptr = NULL;                      ' Just a pointer to a tree in one of the above lists
    int num_queries = smp->sqd->ntaxa - tv->start_tree->nleaves;    ' The number of queries (sequences not in the fixed parental tree)
    char tmp[MAX_TREE_STRING];                  ' BUGGY: hard code length of string (to hold tree). Change constant MAX_TREE_STRING in constants.h.

    if ( num_queries < 1 && smp->set->add_xi ) {
        fprintf(stderr, "%s(%s): No query sequences, cannot run\n", fxn_name, file_name);
        exit(EXIT_FAILURE);
    }

    if( num_queries == 0 ) {
        tv->numTrees = 1;
        tree_list = tv->start_tree;
    } else {

    numCurrentLeaves = tv->start_tree->nleaves;
    tv->numTrees = 2 * numCurrentLeaves - 3;    ' Number of unrooted trees
    
        EnumerateLastTaxon(&tree_list, tv->start_tree); ' Now handles allocation of memory if necessary
    queries_added++;
    
    while (queries_added < num_queries) {
        
        numCurrentLeaves++;
        numTreesPerEnumeration = 2 * numCurrentLeaves - 3;
        
        new_tree_list = (tree *) malloc(sizeof(tree)*tv->numTrees*numTreesPerEnumeration);

        if( new_tree_list == NULL ) {
            fprintf(stderr, "%s: 1memory allocation error\n", fxn_name);
            exit(EXIT_FAILURE);
        }
        
        ' for each tree in PostTree EnumerateLastTaxon()
        for( i = 0; i < tv->numTrees; i++ ) {
            tree_ptr = &new_tree_list[i*numTreesPerEnumeration];
            EnumerateLastTaxon(&tree_ptr, &tree_list[i]);
        }

        tv->numTrees *= numTreesPerEnumeration;
        free(tree_list);
        tree_list = new_tree_list;
        queries_added++;
    }
    }

    tv->tree_list = (tree **) malloc(sizeof(tree *)*tv->numTrees);

    if( tv->tree_list == NULL ) {
        fprintf(stderr, "%s: 2memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    /* For DEBUGGING
    tv->tree_list[0] = &tree_list[2];
    tv->tree_list[0]->tree_index = 0;
    tv->tree_list[1] = &tree_list[6];
    tv->tree_list[1]->tree_index = 1;
    tv->tree_list[2] = &tree_list[4];
    tv->tree_list[2]->tree_index = 2;
    tv->tree_list[3] = &tree_list[5];
    tv->tree_list[3]->tree_index = 3;
    tv->tree_list[4] = &tree_list[7];
    tv->tree_list[4]->tree_index = 4;
    tv->tree_list[5] = &tree_list[8];
    tv->tree_list[5]->tree_index = 5;
    tv->tree_list[6] = &tree_list[1];
    tv->tree_list[6]->tree_index = 6;
    tv->tree_list[7] = &tree_list[3];
    tv->tree_list[7]->tree_index = 7;
    tv->tree_list[8] = &tree_list[0];
    tv->tree_list[8]->tree_index = 8;
    for( i=0; i<tv->numTrees; i++ ) {
        tmp[0] = '\0';
        toString(tmp, tv->tree_list[i]->root, false);
        fprintf(stderr, "Tree %d: %s\n", tv->tree_list[i]->tree_index, tmp);
    }
    */

    /**/
    for( i=0; i<tv->numTrees; i++ ) {
        tmp[0] = '\0';
        tree_list[i].tree_index = i;
        toString(tmp, tree_list[i].root, false);
        fprintf(stderr, "Tree %d: %s\n", i, tmp);
        tv->tree_list[i] = &tree_list[i];
    }
    /**/

    tv->log_numTrees = log(tv->numTrees);
    tv->log_numTrees_minus1 = log(tv->numTrees - 1);
    tv->log_numTrees_minus2 = log(tv->numTrees - 2);
}' Make_Tree_List

' Only handles gmodel prior
static void TopologyPriorMake(topology_vector *tv, const sampler *smp) {
    if( smp->set->gmodel )  {
        if( !top_prior ) TopologyGModelPriorMake(tv, smp);
        tv->top_prior = (void *) top_prior;
        tv->Log_Prior_Ratio = &MonophyleticMarkovPriorRatio;
    } else {
        tv->top_prior = (void *) top_prior;
        tv->Log_Prior_Ratio = &UniformMarkovPriorRatio;
    }
}' TopologyPriorMake

void TopologyGModelPriorMake(const topology_vector *tv, const sampler *smp) {
    const char *fxn_name = "TopologyGModelPriorMake";
    int i;
    int *query_list = (int *)malloc(sizeof(int)*smp->sqd->ntaxa);   ' VERIFIED 2/11/05
    topology_gmodel_prior *tgp = (topology_gmodel_prior *) malloc(sizeof(topology_gmodel_prior));

    if( query_list == NULL || tgp == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    for( i=0; i<smp->sqd->ntaxa; i++ ) {
        if( i < tv->start_tree->nleaves ) query_list[i] = 0;
        else query_list[i] = 1;
    }
    
    tgp->monophyletic = (boolean *) malloc(sizeof(boolean)*tv->numTrees);   ' VERIFIED 2/11/05

    if( tgp->monophyletic == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    tgp->num_mono_trees = 0;
    for( i=0; i<tv->numTrees; i++ ) {
'      char output_string[500] = "";
'      toString(output_string, tv->tree_list[i]->root, false);
        tgp->monophyletic[i] = Monophyletic(tv->tree_list[i], query_list);
        if( tgp->monophyletic[i] ) tgp->num_mono_trees++;
'      fprintf(stderr, "%dth tree %s is %s\n", i, output_string, tgp->monophyletic[i]?"monophyletic":"not monophyletic");
    }
    tgp->num_nonmono_trees = tv->numTrees - tgp->num_mono_trees;
    tgp->log_num_nmtrees[0] = log(tgp->num_nonmono_trees);
    tgp->log_num_nmtrees[1] = log(tgp->num_mono_trees);
    tgp->log_num_nmtrees_minus1[0] = log(tgp->num_nonmono_trees - 1);
    tgp->log_num_nmtrees_minus1[1] = log(tgp->num_mono_trees - 1);
    tgp->log_mono_prob = smp->set->log_mono_prob;
    top_prior = (void *) tgp;
    if(query_list) free(query_list);
}' TopologyGModelPriorMake

static double UniformMarkovPriorRatio(const topology_vector *tv, int dim, ...) {
'fprintf(stderr, "UniformMarkovPriorRatio: %f - %f\n", exp(- (dim - 1)*tv->log_numTrees_minus1 - tv->log_numTrees), exp(tv->log_prior_prob));
    return - (dim - 1)*tv->log_numTrees_minus1 - tv->log_numTrees - tv->log_prior_prob;
}' UniformMarkovPriorRatio

static double MonophyleticMarkovPriorRatio(const topology_vector *tv, int dim, ...) {
    const char *fxn_name = "MonophyleticMarkovPrior";
    boolean local_debug = 0;
    int i;
    double log_mono_prob[2] = {0,0};
    double log_prior;
    tree **trees = NULL;
    topology_gmodel_prior *tgp = (topology_gmodel_prior *) tv->top_prior;
    boolean proposed = false;
    va_list vargs;

    va_start(vargs, dim);
    proposed = va_arg(vargs, int);
    va_end(vargs);
    trees = proposed ? tv->proposed_trees : tv->current_trees;

    log_mono_prob[1] = tgp->log_mono_prob / dim;
    log_mono_prob[0] = log( 1 - exp(log_mono_prob[1]) );

    ' Initial state distribution
    log_prior = log_mono_prob[tgp->monophyletic[trees[0]->tree_index]] - tgp->log_num_nmtrees[tgp->monophyletic[trees[0]->tree_index]];
    if(local_debug) fprintf(stderr, "%s: region 0 with tree %d which is %d: %.4f\n", fxn_name, trees[0]->tree_index, tgp->monophyletic[trees[0]->tree_index], exp(log_mono_prob[tgp->monophyletic[trees[0]->tree_index]] - tgp->log_num_nmtrees[tgp->monophyletic[trees[0]->tree_index]]));
    
    ' Transition probabilities
    for( i=1; i<dim; i++ ) {
        if( tgp->monophyletic[trees[i-1]->tree_index] == tgp->monophyletic[trees[i]->tree_index] ) {
            if(local_debug) fprintf(stderr, "%s: region %d with tree %d which is %d: %.4f\n", fxn_name, i, trees[i]->tree_index, tgp->monophyletic[trees[i]->tree_index], exp(log_mono_prob[tgp->monophyletic[trees[i]->tree_index]] - tgp->log_num_nmtrees_minus1[tgp->monophyletic[trees[i]->tree_index]]));
            log_prior += log_mono_prob[tgp->monophyletic[trees[i]->tree_index]] - tgp->log_num_nmtrees_minus1[tgp->monophyletic[trees[i]->tree_index]];
        } else {
            if(local_debug) fprintf(stderr, "%s: region %d with tree %d which is %d: %.4f\n", fxn_name, i, trees[i]->tree_index, tgp->monophyletic[trees[i]->tree_index], exp(log_mono_prob[tgp->monophyletic[trees[i]->tree_index]] - tgp->log_num_nmtrees[tgp->monophyletic[trees[i]->tree_index]]));
            log_prior += log_mono_prob[tgp->monophyletic[trees[i]->tree_index]] - tgp->log_num_nmtrees[tgp->monophyletic[trees[i]->tree_index]];
        }
    }
    return log_prior - tv->log_prior_prob;
}' MonophyleticMarkovPriorRatio

static tree *Propose_New_Tree(const topology_vector *tv, settings *set, tree *left_tree, tree *curr_tree, tree *right_tree) {
    tree *new_tree = NULL;

    if( left_tree && curr_tree && right_tree && left_tree != curr_tree && left_tree != right_tree && curr_tree != right_tree && tv->numTrees == 3 ) return NULL;

    while( 1 ) {
        new_tree = tv->tree_list[ (int) (set->rng->nextStandardUniform(set->rng)*tv->numTrees) ];
        if( (!left_tree || new_tree != left_tree) && (!curr_tree || new_tree != curr_tree) && (!right_tree || new_tree != right_tree) ) return new_tree;
    }
}' Propose_New_Tree

static void TopologyAddOne(topology_vector *tv, partition_list *pl, sampler *smp) {
    const char *fxn_name = "MarkovTopologyProposeAddOne";
    int proposed;               ' Candidate for location of new topology changepoint
    int land_index;             ' Index of partition where new location lands
    int left_top_change;            ' Next partition to left with different topology
    int right_top_change;           ' Next partition to right with different topology
    int nth_topology_part = 0;      ' Which topology segment are we proposing to split (0-index)?
    partition *land_part;           ' Partition where new position lands
    partition *left_part, *right_part;  ' Partition to left and right of the target segment with different topologies
    partition *new_left_part;       ' New partition for holding left half of split partition
    partition *new_right_part;      ' New partition for holding right half of split partition
    partition *cpart;           ' Temporary holder of a partition
    tree *left_tree;            ' Left tree (index)
    tree *right_tree;           ' Right tree (index)
    double *pPartialLogLikelihood = (double *) malloc(sizeof(double)*pl->npartitions);      ' List of new likelihoods for all affected regions
    double pLogLikelihood = 0.0;        ' Sum of above
    double cLogLikelihood = 0.0;        ' Current log likelihood of the same region
    double proposed_left_likelihood=0.0, proposed_right_likelihood=0.0;
    double current_left_likelihood=0.0, current_right_likelihood=0.0;
    double priorRatio, deathProb, birthProb, logRatio;
    boolean keep_on_left;           ' Keep the tree on the left
    double log_number_tree_choices;
    double selection_prob;
    double log_tau_prior_ratio;
    partition_list *tpl;
    int i;

    ' Catch and handle possible memory allocation error
    if( pPartialLogLikelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    ' x \in \Re^n --> x' \in \Re^m
    proposed = pl->Propose_Change_Point(pl, smp->set, true);            ' Select a new topology change point               WP: 1/(lenSeq - topology_changes - 1)

    ' Find partition where new point landed (index arguments are inclusive)
    land_index = PartitionContaining(pl, proposed, 0, pl->npartitions - 1);
    land_part = pl->part[land_index];
    
    ' Find relative position of topology part with proposed as its left boundary
    if( smp->set->alawadhi_debug ) {
        nth_topology_part = 0;
        i = 1;
        while( i<=land_index ) if( pl->part[i++]->topchange ) nth_topology_part++;
        nth_topology_part++;
    }

    ' Find left and right parameter change point
    left_top_change = land_index;
    while( !pl->part[left_top_change]->topchange ) left_top_change--;
    right_top_change = land_index + 1;
    while( right_top_change < pl->npartitions && !pl->part[right_top_change]->topchange ) right_top_change++;

    ' Obtain a pointer to the boundary partitions
    left_part = land_part;
    if( left_top_change > 0 ) left_part = pl->part[left_top_change - 1];
    right_part = land_part;
    if( right_top_change < pl->npartitions ) right_part = pl->part[right_top_change];

    ' Propose new topologies
    if( smp->set->rng->nextStandardUniform(smp->set->rng) < 0.5 ) {                     ' Keep current tree on the left, new tree on right WP: 0.5
        keep_on_left = true;
        left_tree = land_part->ctree;

        ' Propose new right tree (and proposal probability)
        if( right_top_change == pl->npartitions ) {                         ' Adding to right end                  WP: handled
            right_tree = tv->Propose_New_Tree(tv, smp->set, left_tree, NULL, NULL);         ' Propose any tree not equivalent to left_tree     WP: 1/(numTrees-1)
            log_number_tree_choices = tv->log_numTrees_minus1;
        } else {                                            '                          WP: handled
            right_tree = tv->Propose_New_Tree(tv, smp->set, left_tree, NULL, right_part->ctree);    ' Propose any tree not like left or right      WP: 1/(numTrees-2)
            log_number_tree_choices = tv->log_numTrees_minus2;
        }

        ' Compute current and proposed likelihood on the right
        for( i=land_index+1; i<right_top_change; i++ ) {
            cpart = pl->part[i];
            if(compute_likelihood) pPartialLogLikelihood[i] = TreeLogLikelihood(right_tree, smp, cpart->cmatrix, cpart->counts, cpart->cHyperParameter, false);
            else pPartialLogLikelihood[i] = 0.0;
            pLogLikelihood += pPartialLogLikelihood[i];
            if( compute_likelihood ) cLogLikelihood += cpart->cPartialLogLikelihood;
        }
    }
    else {                                                  ' Keep current tree on right               WP: 0.5
        keep_on_left = false;
        right_tree = land_part->ctree;

        ' Propose new tree on left (and compute proposal probability)
        if( !left_top_change ) {                                    ' Adding to left end                   WP: handled
            left_tree = tv->Propose_New_Tree(tv, smp->set, NULL, NULL, right_tree);         ' Propose new tree on the end              WP: 1/(numTrees-1)
            log_number_tree_choices = tv->log_numTrees_minus1;
        } else {                                            ' Adding in middle                 WP: handled
            left_tree = tv->Propose_New_Tree(tv, smp->set, left_part->ctree, NULL, right_tree); ' Propose new tree in middle               WP: 1/(numTrees-2)
            log_number_tree_choices = tv->log_numTrees_minus2;
        }

        if( left_top_change < 0 || land_index > pl->npartitions ) {
            fprintf(stderr, "%s: debugging memory write error\n", fxn_name);
            exit(EXIT_FAILURE);
        }
        ' Compute current and proposal likelihood on left
        for( i=left_top_change; i<land_index; i++) {
            cpart = pl->part[i];
            if(compute_likelihood) pPartialLogLikelihood[i] = TreeLogLikelihood(left_tree, smp, cpart->cmatrix, cpart->counts, cpart->cHyperParameter, false);
            else pPartialLogLikelihood[i] = 0.0;
            pLogLikelihood += pPartialLogLikelihood[i];
            if( compute_likelihood ) cLogLikelihood += cpart->cPartialLogLikelihood;
        }
    }

    ' Construct new partitions if proposal splits an existing one
    if( land_part->left != proposed ) {             ' It can if !land_part->topchange
        PartitionMake(&new_left_part, smp->sqd->lenunique, land_part->left, proposed-1, land_part->topchange, land_part->parchange);
        PartitionCopySegmentCounts(new_left_part, smp->sqd, land_part->left, proposed);
        PartitionMake(&new_right_part, smp->sqd->lenunique, proposed, land_part->right, true, false); ' BUG
        PartitionCopyPartitionCountDifferences(new_right_part, land_part, new_left_part);
    }

    ' Compute current and proposal likelihood on the landing partition
    if(land_part->left != proposed) {
        ' TODO: there is a slightly wasted computation on the event of no acceptance, but trying to eliminate it has generated bugs in the past!
        if( compute_likelihood )  { ' WAS BUG: computed one of these on mismatch of left tree with new_right
            current_right_likelihood = TreeLogLikelihood(land_part->ctree, smp, land_part->cmatrix, new_right_part->counts, land_part->cHyperParameter, false);
            proposed_right_likelihood = keep_on_left ? TreeLogLikelihood(right_tree, smp, land_part->cmatrix, new_right_part->counts, land_part->cHyperParameter, false)
                : current_right_likelihood;
        }
        pLogLikelihood += proposed_right_likelihood;
        cLogLikelihood += current_right_likelihood;
        if( compute_likelihood )  { ' WAS BUG: computed one of these on mismatch of left tree with new_right
            current_left_likelihood = TreeLogLikelihood(land_part->ctree, smp, land_part->cmatrix, new_left_part->counts, land_part->cHyperParameter, false);
            proposed_left_likelihood = keep_on_left ? current_left_likelihood
                : TreeLogLikelihood(left_tree, smp, land_part->cmatrix, new_left_part->counts, land_part->cHyperParameter, false);
        }
        pLogLikelihood += proposed_left_likelihood;
        cLogLikelihood += current_left_likelihood;
    }
    else if(keep_on_left) {
        if( compute_likelihood ) proposed_right_likelihood = TreeLogLikelihood(right_tree, smp, land_part->cmatrix, land_part->counts, land_part->cHyperParameter, false);
        pLogLikelihood += proposed_right_likelihood;
        if( compute_likelihood ) {
            current_right_likelihood = land_part->cPartialLogLikelihood;
            cLogLikelihood += land_part->cPartialLogLikelihood;
        }
    }

    ' Compute proposal probability for reverse move (DeleteOne)
    selection_prob = tv->Reverse_Add_One_Log_Proposal_Probability(tv, pl, land_index, keep_on_left, (keep_on_left?right_tree:left_tree));

    if( debug>3 || global_debug>3 ) smp->Report_Proposed_State(smp, "AddOne", tv->proposed_trees, pl->topology_changes + 1 + 1);

    ' Compute prior ratio on topology \tau (based on proposed_trees)
    log_tau_prior_ratio = tv->Log_Prior_Ratio(tv, pl->topology_changes + 1 + 1, true);
    
    ' x n-->m x' -> ... -> x*
    ' ln[ pi(x') / pi(x) ] =
    ' ln[ l(x') / l(x) ] = pLogLikelihood - cLogLikelihood
    ' + ln[ q(x') / q(x) ] = priorRatio
    priorRatio = pl->log_top_lambda                 ' K
        ' - log(dcp->topology_changes + 1)         ' K (cancels with new prior on change points)
        - log(pl->alignment_length - pl->topology_changes - 1)  ' \xi
        + log_tau_prior_ratio                   ' \tau
        ;
    ' ln[ q_{nm}(x',x) / q_{mn}(x,x') ] = deathProb - birthProb
    birthProb =
        + log(pl->top_one_bk)                   ' BIRTH: probability of AddOne
        - log(pl->alignment_length - pl->topology_changes - 1)  ' BIRTH: number of choices for new \xi
        - log_number_tree_choices               ' BIRTH: number of choices for new \tau
        - log_two                       ' BIRTH: KSD: choice of keep_on_left vs. !keep_on_left during birth (no longer cancels when use prob below)
        ;
    deathProb = log(pl->top_one_dkp1)               ' DEATH: death prob
        + log(selection_prob)                   ' DEATH: KSD: probability of selecting this change point for deletion
        ;

    if( smp->set->alawadhi_topology_one || smp->set->alawadhi_topology || smp->set->alawadhi ) {
        ' x m-->n x' -> ... -> x* : min{ 1, pi_n(x*) pi_n^*(x') q_{nm}(x',x) / pi_m(x) / pi_n^*(x*) / q_{mn}(x,x') }
        double ologRatio = pLogLikelihood - cLogLikelihood + priorRatio + deathProb - birthProb;
        double cll = 0.0, pll = 0.0, ill=0.0, clp, ilp, plp;
        ' q_{nm}(x',x) / q_{mn}(x,x') : proposal ratio unchanged
        logRatio = deathProb - birthProb;

        ' l(x): current log likelihood
        if(compute_likelihood) for( i=0; i<pl->npartitions; i++ ) cll += pl->part[i]->cPartialLogLikelihood;    ' l(x)
        ' q(x): current log prior
        clp = smp->Log_Prior(smp, false);   ' false means don't use temporary state vector

        ' x' -> ... -> x* (in tpl)
        smp->Alawadhi_Copy_State(smp, &tpl);    ' Copy current state; tpl points to copy
        AddOneAccept(smp, tpl, new_left_part, new_right_part, tv, proposed, left_top_change, land_index, right_top_change, keep_on_left, left_tree, right_tree, pPartialLogLikelihood, proposed_left_likelihood, proposed_right_likelihood, log_tau_prior_ratio);

        ' l(x'): intermediate log likelihood
        if(compute_likelihood) for( i=0; i<tpl->npartitions; i++ ) ill += tpl->part[i]->cPartialLogLikelihood;  ' l(x)
        ' q(x'): intermediate log prior
        ilp = smp->Log_Prior(smp, true);    ' true means use temporary state vector

        ' Fixed dimension sampler: x' -> x*
        smp->Fixed_Dimension_Sampler(smp, smp->set->alawadhi_k, TAU|XI);

        ' l(x*): Total proposed log likelihood
        if(compute_likelihood) for( i=0; i<tpl->npartitions; i++ ) pll += tpl->part[i]->cPartialLogLikelihood;
        ' q(x*): intermediate log prior
        plp = smp->Log_Prior(smp, true);    ' true means use temporary state vector (now updated by fixed dimension sampling)

        ' pi_n(x*) / pi_n^*(x*) = (1-AF)*pi_n(x*) = (1-AF)*l(x*)*q(x*)
        logRatio += smp->set->alawadhi_factor*(ill + ilp - pll - plp) + pll + plp - cll - clp;

        if( smp->set->alawadhi_debug && logRatio > log(0.2) ) {
            int seg_index = 1; i = 0; while( i<nth_topology_part ) if( tpl->part[seg_index++]->topchange ) i++; seg_index--;
            fprintf(stderr, "ALAWADHI_DEBUG (%20s): %d:%d -> %d:%d move prob: %e; prior: %e -> %e -> %e; logRatio: %e -> %e\n", smp->move_names[ADD_XI], proposed, keep_on_left?right_tree->tree_index:left_tree->tree_index, tpl->part[seg_index]->left, keep_on_left?tpl->part[seg_index]->ctree->tree_index:tpl->part[seg_index-1]->ctree->tree_index, (deathProb - birthProb), (clp+cll), (ilp+ill), (plp+pll), (ologRatio), (logRatio));
            if( pll+plp-ill-ilp > 0 ) fprintf(stderr, "+++++++ ");
            else fprintf(stderr, "------- ");
            if( logRatio > ologRatio ) fprintf(stderr, "+++++++ %f\n", (double)smp->acceptancerate[ADD_XI]/smp->tries[ADD_XI]);
            else fprintf(stderr, "------- %f\n", (double)smp->acceptancerate[ADD_XI]/smp->tries[ADD_XI]);
        }
    } else {
        logRatio = pLogLikelihood - cLogLikelihood + priorRatio + deathProb - birthProb;

    }
    smp->tries[ADD_XI]++;

    ' TODO: this report is not correct if Al-Awadhi is turned on
    if( debug>1 || global_debug>1 )
        smp->Report_Proposal_Statistics(smp, "AddOne", priorRatio, birthProb, deathProb, pLogLikelihood, cLogLikelihood, proposed, keep_on_left?right_tree->tree_index:left_tree->tree_index, keep_on_left);

    if( logRatio > 0 || smp->set->rng->nextStandardUniform(smp->set->rng) < exp(logRatio) ) {
        smp->acceptancerate[ADD_XI]++;

        if( smp->set->alawadhi_topology_one || smp->set->alawadhi_topology || smp->set->alawadhi ) {
            smp->Alawadhi_Accept(smp);
        } else {
            AddOneAccept(smp, pl, new_left_part, new_right_part, tv, proposed, left_top_change, land_index, right_top_change, keep_on_left, left_tree, right_tree, pPartialLogLikelihood, proposed_left_likelihood, proposed_right_likelihood, log_tau_prior_ratio);
        }
        if( debug>0 || global_debug>0 ) smp->Report_State(smp, "AddOne", logRatio, log_tau_prior_ratio);
        if( debug>3 || global_debug>3 || global_debug==-1 ) {
            VerifyLikelihood(smp, false);
            VerifyCounts(smp, "TopologyAddOne", false);
        }
    } else {
        if( smp->set->alawadhi_topology_one || smp->set->alawadhi_topology || smp->set->alawadhi ) {
            smp->Alawadhi_Reject(smp);  ' Resets tree_vec
        } else {
            if( smp->set->gmodel || debug>3 || global_debug>3 ) {
                if( tv->proposed_trees ) free(tv->proposed_trees);
                tv->proposed_trees = NULL;
            }
                if( land_part->left != proposed ) {
                PartitionDelete(new_right_part);
                PartitionDelete(new_left_part);
            }
        }
    }
    if(pPartialLogLikelihood) free(pPartialLogLikelihood);
}' TopologyAddOne

static void AddOneAccept(sampler *smp, partition_list *pl, partition *new_left_part, partition *new_right_part, topology_vector *tv, int proposed, int left_top_change, int land_index, int right_top_change, boolean keep_on_left, tree *left_tree, tree *right_tree, double *pPartialLogLikelihood, double proposed_left_likelihood, double proposed_right_likelihood, double log_tau_prior_ratio) {
    int i;
    partition *land_part = pl->part[land_index];

    ' Handle acceptance for \tau
    tv->log_prior_prob += log_tau_prior_ratio;
    if( smp->set->gmodel || debug>3 || global_debug>3 ) {
        if( tv->current_trees ) free(tv->current_trees);
        tv->current_trees = tv->proposed_trees;
        tv->proposed_trees = NULL;
    }

    ' Update intermediate likelihoods and topologies
    if( keep_on_left ) {
        for( i=land_index+1; i<right_top_change; i++) {
            partition *cpart = pl->part[i];
            cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
            cpart->ctree = right_tree;
        }
    }
    else {
        for( i=left_top_change; i<land_index; i++) {
            partition *cpart = pl->part[i];
            cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
            cpart->ctree = left_tree;
        }
    }

    ' Insert new topology change point : on an existing parameter change point
    if( land_part->left == proposed ) {
        land_part->topchange = true;
        if( keep_on_left ) {
            land_part->cPartialLogLikelihood = proposed_right_likelihood;
            land_part->ctree = right_tree;
        }
    }
    ' Insert new topology change point
    else {
        ' Update left side
        memcpy(land_part->counts, new_left_part->counts, sizeof(int)*land_part->lenunique);
        land_part->right = proposed - 1;
        land_part->cPartialLogLikelihood = proposed_left_likelihood;
            'keep_on_left ? TreeLogLikelihood(land_part->ctree, smp, land_part->cmatrix, land_part->counts, land_part->cHyperParameter, false) : proposed_left_likelihood;
        land_part->ctree = left_tree;

        ' Update right side
        new_right_part->cPartialLogLikelihood = proposed_right_likelihood;
            'keep_on_left ? proposed_right_likelihood : TreeLogLikelihood(right_tree, smp, land_part->cmatrix, new_right_part->counts, land_part->cHyperParameter, false);
        land_part->cmatrix->Matrix_Make_Copy(&new_right_part->cmatrix, land_part->cmatrix);
        new_right_part->cHyperParameter = land_part->cHyperParameter;
        new_right_part->cPartialLogHyperParameterPrior = land_part->cPartialLogHyperParameterPrior;
        new_right_part->ctree = right_tree;

        PartitionListAddPartition(pl, new_right_part, land_index+1);
        PartitionDelete(new_left_part);
    }

    pl->topology_changes++;
}' AddOneAccept

static void TopologyAddTwo(topology_vector *tv, partition_list *pl, sampler *smp) {
    const char *fxn_name = "TopologyAddTwo";
    int left_proposed = 0;          ' Left end point
    int right_proposed = 0;         ' Right end point
    int left_land_index = 0;        ' Index of region on left boundary
    int right_land_index = 0;       ' Index of region on right boundary
    int left_boundary_index = 0;        ' Index of left-most region in this topology segment
    int right_boundary_index = 0;       ' Index of left-most region in next topology segment (or npartitions)
    int number_xi_choices = 0;      ' Length of region in which to propose second topology changepoint
    tree *middle_tree;          ' New tree for inserted region
    partition *land_part_left=NULL;     ' Landing partition of left_proposed
    partition *land_part_right=NULL;    ' Landing partition of right_proposed
    partition *middle_left_part=NULL;   ' New partition to hold data in right portion of land_part_left (if applicable)
    partition *middle_right_part=NULL;  ' New partition to hold data in left portion of land_part_right (if applicable)
    double *pPartialLogLikelihood = (double *) malloc(sizeof(double)*pl->npartitions);
    double pLeftLogLikelihood = 0.0, cLeftLogLikelihood = 0.0, pRightLogLikelihood = 0.0, cRightLogLikelihood = 0.0;
    double pLogLikelihood = 0.0;        ' Total proposed log likelihood
    double cLogLikelihood = 0.0;        ' Total current log likelihood
    settings *set = smp->set;
    double priorRatio, deathProb, birthProb, logRatio, selection_prob, log_tau_prior_ratio;
    int i;
    boolean local_debug = 0;

    if( pPartialLogLikelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    ' TODO: The following should compute its own proposal probability, rather than relying on a match with assumptions in this code (harder to debug though)
    left_proposed = pl->Propose_Change_Point(pl, set, true);    ' Uniform over remaining choices               WP: 1/(lenSeq - topology_changes - 1)
    'first_proposed = pl->Propose_Change_Point(pl, set, true);
    left_land_index = PartitionContaining(pl, left_proposed, 0, pl->npartitions - 1);
    'first_land_index = PartitionContaining(pl, first_proposed, 0, pl->npartitions - 1);
    left_boundary_index = left_land_index;
    'left_boundary_index = first_land_index;
    while( left_boundary_index > 0 && !pl->part[left_boundary_index]->topchange ) left_boundary_index--;
    right_boundary_index = left_land_index+1;
    'right_boundary_index = first_land_index+1;
    while( right_boundary_index < pl->npartitions && !pl->part[right_boundary_index]->topchange ) right_boundary_index++;
    number_xi_choices = pl->part[right_boundary_index - 1]->right - pl->part[left_boundary_index]->left - 1;    ' include right, not include left, not include proposed
    /*
    if( set->rng->nextStandardUniform(set->rng) < 0.5 ) {
        number_xi_choices = first_proposed - pl->part[left_boundary_index]->left - 1;   ' exclusive both ends
        if( number_xi_choices < 1 ) return;
        right_proposed = first_proposed;
        right_land_index = first_land_index;
        left_proposed = pl->part[left_boundary_index]->left + 1 + (int) ((double)set->rng->nextStandardUniform(set->rng)*number_xi_choices);
        left_land_index = PartitionContaining(pl, left_proposed, left_boundary_index, right_boundary_index-1);
    } else {
        number_xi_choices = pl->part[right_boundary_index - 1]->right - first_proposed; ' exclusive left end
        if( number_xi_choices < 1 ) return;
        left_proposed = first_proposed;
        left_land_index = first_land_index;
        right_proposed = left_proposed + 1 + (int) ((double)set->rng->nextStandardUniform(set->rng)*number_xi_choices);
        right_land_index = PartitionContaining(pl, right_proposed, left_boundary_index, right_boundary_index-1);
    }

    number_xi_choices = (double) 1 / (right_proposed - pl->part[left_boundary_index]->left - 1);    ' exclusive both ends
    number_xi_choices += (double) 1 / (pl->part[right_boundary_index-1]->right - first_proposed);   ' exclusive both ends
    */
    if( number_xi_choices < 1 ) return; ' No room for two change points
    
    ' Should add to proposal probability
    right_proposed = pl->Propose_Second_Change_Point(pl, set, left_boundary_index, right_boundary_index, left_proposed, true);
    right_land_index = PartitionContaining(pl, right_proposed, left_boundary_index, right_boundary_index - 1);  ' Inclusive indices
    ' Should add to proposal probability
    
    middle_tree = tv->Propose_New_Tree(tv, set, pl->part[left_land_index]->ctree, NULL, NULL);  ' Propose tree for new region      WP: 1/(numTrees-1)
    selection_prob = tv->Reverse_Add_Two_Log_Proposal_Probability(tv, pl, left_land_index, middle_tree);

    if( left_proposed > right_proposed ) {
        int tmp_proposed = right_proposed;
        int tmp_land_index = right_land_index;
        right_proposed = left_proposed;
        left_proposed = tmp_proposed;
        right_land_index = left_land_index;
        left_land_index = tmp_land_index;
    }

    log_tau_prior_ratio = tv->Log_Prior_Ratio(tv, pl->topology_changes + 1 + 2, true);

    priorRatio = 2*pl->log_top_lambda               ' K: going K to K+2
'      - log(pl->topology_changes + 2)         ' K
'      - log(pl->topology_changes + 1)         ' K
        - log(pl->alignment_length - pl->topology_changes - 1)  ' \xi
        - log(pl->alignment_length - pl->topology_changes - 2)  ' \xi
        + log_tau_prior_ratio
'      - 2*dcp->log_numTrees_minus1                ' \tau
        ;
    deathProb = log(pl->top_two_dkp2)           ' DEATH: death probability
        + log(selection_prob)               ' DEATH: KSD: probably of selecting this pair of changepoints for deletion
        ;
    birthProb =
        + log(pl->top_two_bk)               ' BIRTH: birth probability
        - log(pl->alignment_length - pl->topology_changes - 1)' BIRTH: number of choices for first topology change point
        - log(number_xi_choices)            ' BIRTH: number of choices for second topology change point
'      + log(number_xi_choices)            ' BIRTH: number of choices for second topology change point
        - tv->log_numTrees_minus1           ' BIRTH: number of choices for new topology during birth
'      - log_two                   ' BIRTH: whether second change point is left or right of first
        + log_two                   ' BIRTH: there are two ways to select these two topology change points: left first or right first
        ;

    land_part_left = pl->part[left_land_index];
    land_part_right = pl->part[right_land_index];

    if( land_part_left->left != left_proposed ) {
        PartitionMake(&middle_left_part, smp->sqd->lenunique, left_proposed, (left_land_index==right_land_index)?right_proposed-1:land_part_left->right, true, false);
        PartitionCopySegmentCounts(middle_left_part, smp->sqd, left_proposed, (left_land_index==right_land_index)?right_proposed:land_part_left->right+1);
    } else middle_left_part = land_part_left;
    if( land_part_right->left != right_proposed && left_land_index != right_land_index ) {
        PartitionMake(&middle_right_part, smp->sqd->lenunique, land_part_right->left, right_proposed-1, land_part_right->topchange, land_part_right->parchange);
        PartitionCopySegmentCounts(middle_right_part, smp->sqd, land_part_right->left, right_proposed);
    }

    for(i=left_land_index+1; i<right_land_index; i++) {
        partition *cpart = pl->part[i];
        if(compute_likelihood) pPartialLogLikelihood[i] = TreeLogLikelihood(middle_tree, smp, cpart->cmatrix, cpart->counts, cpart->cHyperParameter, false);
        else pPartialLogLikelihood[i] = 0.0;
        pLogLikelihood += pPartialLogLikelihood[i];
        if( compute_likelihood ) cLogLikelihood += cpart->cPartialLogLikelihood;
    }
    if( compute_likelihood ) {
        if( compute_likelihood)
            pLeftLogLikelihood = TreeLogLikelihood(middle_tree, smp, land_part_left->cmatrix, middle_left_part->counts, land_part_left->cHyperParameter, false);
        pLogLikelihood += pLeftLogLikelihood;
        if( compute_likelihood ) cLeftLogLikelihood = TreeLogLikelihood(land_part_left->ctree, smp, land_part_left->cmatrix, middle_left_part->counts, land_part_left->cHyperParameter, false);
        cLogLikelihood += cLeftLogLikelihood;
        if( land_part_right->left != right_proposed && left_land_index != right_land_index ) {
            if( compute_likelihood )
                pRightLogLikelihood = TreeLogLikelihood(middle_tree, smp, land_part_right->cmatrix, middle_right_part->counts, land_part_right->cHyperParameter, false);
            pLogLikelihood += pRightLogLikelihood;
            if( compute_likelihood )
                cRightLogLikelihood = TreeLogLikelihood(land_part_right->ctree, smp, land_part_right->cmatrix, middle_right_part->counts, land_part_right->cHyperParameter, false);
            cLogLikelihood += cRightLogLikelihood;
        }
    }

    if( smp->set->alawadhi_topology_two || smp->set->alawadhi_topology || smp->set->alawadhi ) {
        ' x m-->n x' -> ... -> x* : min{ 1, pi_n(x*) pi_n^*(x') q_{nm}(x',x) / pi_m(x) / pi_n^*(x*) / q_{mn}(x,x') }
        double ologRatio = pLogLikelihood - cLogLikelihood + priorRatio + deathProb - birthProb;
        double cll = 0.0, pll = 0.0, ill=0.0, clp, ilp, plp;
        partition_list *tpl = NULL;
        int nthone = 0;

        ' q_{nm}(x',x) / q_{mn}(x,x') : proposal ratio unchanged
        logRatio = deathProb - birthProb;

        ' l(x): current log likelihood
        if(compute_likelihood) for( i=0; i<pl->npartitions; i++ ) cll += pl->part[i]->cPartialLogLikelihood;    ' l(x)
        ' q(x): current log prior
        clp = smp->Log_Prior(smp, false);   ' false means use pl

        if( smp->set->alawadhi_debug ) {
            i = 1;
            while( i<=left_land_index ) if( pl->part[i++]->topchange ) nthone++;
        }

        ' x' -> ... -> x* (in tpl)
        smp->Alawadhi_Copy_State(smp, &tpl);
        AddTwoAccept(smp, tpl, tv, left_land_index, right_land_index, middle_left_part, left_proposed, right_proposed, middle_tree, pLeftLogLikelihood, cLeftLogLikelihood, pRightLogLikelihood, pPartialLogLikelihood, log_tau_prior_ratio, local_debug);

        ' l(x'): intermediate log likelihood
        if(compute_likelihood) for( i=0; i<tpl->npartitions; i++ ) ill += tpl->part[i]->cPartialLogLikelihood;  ' l(x)
        ' q(x'): intermediate log prior
        ilp = smp->Log_Prior(smp, true);    ' true means use tpl

        ' Fixed dimension sampler: x' -> x*
        smp->Fixed_Dimension_Sampler(smp, smp->set->alawadhi_k, TAU|XI);

        ' l(x*): Total proposed log likelihood
        if(compute_likelihood) for( i=0; i<tpl->npartitions; i++ ) pll += tpl->part[i]->cPartialLogLikelihood;
        ' q(x*): intermediate log prior
        plp = smp->Log_Prior(smp, true);    ' true means use tpl (now updated by fixed dimension sampling)

        ' pi_n(x*) / pi_n^*(x*) = (1-AF)*pi_n(x*) = (1-AF)*l(x*)*q(x*)
        logRatio += smp->set->alawadhi_factor*(ill + ilp - pll - plp) + pll + plp - cll - clp;

        if( smp->set->alawadhi_debug && logRatio > log(0.2) ) {
            int leftone, rightone;
            leftone = 1; i = 0; while( i <= nthone ) if( tpl->part[leftone++]->topchange ) i++; leftone--;
            rightone = 1; i= 0; while( i <= nthone+1) if( tpl->part[rightone++]->topchange ) i++; rightone--;
            fprintf(stderr, "ALAWADHI_DEBUG (%20s): (%d,%d:%d) -> (%d,%d:%d) move prob: %e; %e -> %e -> %e; logRatio: %e -> %e\n", smp->move_names[ADD_TWO_XI], left_proposed, right_proposed, middle_tree->tree_index, tpl->part[leftone]->left, tpl->part[rightone]->left, tpl->part[leftone]->ctree->tree_index, deathProb - birthProb, clp+cll, ilp+ill, plp+pll, (ologRatio), (logRatio));
            if( pll+plp-ill-ilp > 0 ) fprintf(stderr, "+++++++ ");
            else fprintf(stderr, "------- ");
            if( logRatio > ologRatio ) fprintf(stderr, "+++++++ %f\n", (double)smp->acceptancerate[ADD_TWO_XI]/smp->tries[ADD_TWO_XI]);
            else fprintf(stderr, "------- %f\n", (double)smp->acceptancerate[ADD_TWO_XI]/smp->tries[ADD_TWO_XI]);
        }
    } else {
        logRatio = pLogLikelihood - cLogLikelihood + priorRatio + deathProb - birthProb;
    }

    if( debug>1 || global_debug>1 || local_debug>1 )
        smp->Report_Proposal_Statistics(smp, "AddTwo", priorRatio, birthProb, deathProb, pLogLikelihood, cLogLikelihood, left_proposed, right_proposed, middle_tree->tree_index);

        smp->tries[ADD_TWO_XI]++;

    if( logRatio > 0.0 || set->rng->nextStandardUniform(set->rng) < exp(logRatio) ) {
        smp->acceptancerate[ADD_TWO_XI]++;

        if( smp->set->alawadhi_topology_two || smp->set->alawadhi_topology || smp->set->alawadhi )
            smp->Alawadhi_Accept(smp);
        Else
            AddTwoAccept(smp, pl, tv, left_land_index, right_land_index, middle_left_part, left_proposed, right_proposed, middle_tree, pLeftLogLikelihood, cLeftLogLikelihood, pRightLogLikelihood, pPartialLogLikelihood, log_tau_prior_ratio, local_debug);

        if( debug>0 || global_debug>0 || local_debug>0 ) smp->Report_State(smp, "AddTwo", logRatio, log_tau_prior_ratio);
        if( debug>3 || global_debug>3 || global_debug==-1 ) {
            VerifyLikelihood(smp, false);
            VerifyCounts(smp, "TopologyAddTwo", false);
        }
    } else {
        if( smp->set->alawadhi_topology_two || smp->set->alawadhi_topology || smp->set->alawadhi ) {
            smp->Alawadhi_Reject(smp);
        } else {
            if( set->gmodel || debug>3 || global_debug>3 || local_debug>3 ) {
                free(tv->proposed_trees);
                tv->proposed_trees = NULL;
            }
            if( land_part_left->left != left_proposed ) PartitionDelete(middle_left_part);
            if( land_part_right->left != right_proposed && left_land_index != right_land_index ) PartitionDelete(middle_right_part);
        }
    }
    if( pPartialLogLikelihood ) free(pPartialLogLikelihood);
}' TopologyAddTwo

static void AddTwoAccept(sampler *smp, partition_list *pl, topology_vector *tv, int left_land_index, int right_land_index, partition *middle_left_part, int left_proposed, int right_proposed, tree *middle_tree, double pLeftLogLikelihood, double cLeftLogLikelihood, double pRightLogLikelihood, double *pPartialLogLikelihood, double log_tau_prior_ratio, int local_debug) {
    int i;
    partition *land_part_left = pl->part[left_land_index];
    partition *land_part_right = pl->part[right_land_index];

    tv->log_prior_prob += log_tau_prior_ratio;
    if( smp->set->gmodel || debug>3 || global_debug>3 || local_debug>3 ) {
        if( tv->current_trees ) free(tv->current_trees);
        tv->current_trees = tv->proposed_trees;
        tv->proposed_trees = NULL;
    }

    ' Update intermediate likelihoods and topologies
    for(i = left_land_index+1; i < right_land_index; i++) {
        partition *cpart = pl->part[i];
        cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
        cpart->ctree = middle_tree;
    }

    ' Insert first new toplogy change point
    ' Turn an existing parameter change point into a topology change point also
    if( land_part_right->left != right_proposed ) {
        partition *far_right_part = NULL;
        PartitionMake(&far_right_part, smp->sqd->lenunique, right_proposed, land_part_right->right, true, false);
        PartitionCopySegmentCounts(far_right_part, smp->sqd, right_proposed, land_part_right->right+1);
        land_part_right->cmatrix->Matrix_Make_Copy(&far_right_part->cmatrix, land_part_right->cmatrix);
        far_right_part->cHyperParameter = land_part_right->cHyperParameter;
        far_right_part->cPartialLogHyperParameterPrior = land_part_right->cPartialLogHyperParameterPrior;
        far_right_part->ctree = land_part_right->ctree;
        if( compute_likelihood )
            far_right_part->cPartialLogLikelihood = TreeLogLikelihood(far_right_part->ctree, smp, far_right_part->cmatrix, far_right_part->counts, far_right_part->cHyperParameter, false);
        PartitionListAddPartition(pl, far_right_part, right_land_index+1);

        if( left_land_index == right_land_index ) {
            land_part_left->cPartialLogLikelihood -= far_right_part->cPartialLogLikelihood;
            land_part_left->right = right_proposed-1;
            PartitionSubtractPartition(land_part_left, far_right_part);
        } else {
            land_part_right->cPartialLogLikelihood = pRightLogLikelihood;
            land_part_right->ctree = middle_tree;
            land_part_right->right = right_proposed - 1;
            PartitionSubtractPartition(land_part_right, far_right_part);
        }
    } else {
        land_part_right->topchange = true;
    }
    if( land_part_left->left != left_proposed ) {
        ' middle_left_part is either affected region (if left_land_index == right_land_index) else the leftmost partition of middle segment
        middle_left_part->cPartialLogLikelihood = pLeftLogLikelihood;
        land_part_left->cmatrix->Matrix_Make_Copy(&middle_left_part->cmatrix, land_part_left->cmatrix);
        middle_left_part->cHyperParameter = land_part_left->cHyperParameter;
        middle_left_part->cPartialLogHyperParameterPrior = land_part_left->cPartialLogHyperParameterPrior;
        middle_left_part->ctree = middle_tree;
        PartitionListAddPartition(pl, middle_left_part, left_land_index+1);

        land_part_left->cPartialLogLikelihood -= cLeftLogLikelihood;
        PartitionSubtractPartition(land_part_left, middle_left_part);
        land_part_left->right = left_proposed-1;
    } else {
        land_part_left->topchange = true;
        land_part_left->cPartialLogLikelihood = pLeftLogLikelihood;
        land_part_left->ctree = middle_tree;
    }
    pl->topology_changes += 2;
    if( debug>3 || global_debug>3 || local_debug>3 || global_debug==-1 ) {
        VerifyCounts(smp, "AddTwoAccept", smp->set->alawadhi_topology_two || smp->set->alawadhi_topology || smp->set->alawadhi);
        VerifyLikelihood(smp, smp->set->alawadhi_topology_two || smp->set->alawadhi_topology || smp->set->alawadhi);
    }
}' AddTwoAccept

static void TopologyDeleteOne(topology_vector *tv, partition_list *pl, sampler *smp) {
    const char *fxn_name = "TopologyDeleteOne";
    int i;
    int propose_to_delete_index;        ' region to delete
    int left_top_change_index;      ' first region on left with topology change
    int right_top_change_index;     ' first region on right with topology change
    partition *propose_to_delete_part;  ' Pointer to partition to delete
    partition *left_part;           ' Region just left of the one proposed to delete
    tree *after_collapse_tree = NULL;   ' The topology that will spread in to replace old topology
    boolean keep_on_left = false;       ' Determine which tree (on left or right) that will spread to fill removed topology's segments
    double log_number_tree_choices;
    double *pPartialLogLikelihood = NULL;
    double pLogLikelihood = 0.0, cLogLikelihood = 0.0, logRatio, priorRatio, birthProb, deathProb;
    double selection_prob = 0.0, log_tau_prior_ratio;
    double cll=0.0, clp=0.0, ill=0.0, ilp=0.0;
    partition_list *tpl = NULL;


    ' Before we propose to move to another dimension, do some fixed dimension moves
    if( smp->set->alawadhi_topology_one || smp->set->alawadhi_topology || smp->set->alawadhi ) {
        for( i=0; i<pl->npartitions; i++ ) cll += pl->part[i]->cPartialLogLikelihood;   ' log l(x*)
        clp = smp->Log_Prior(smp, false);   ' log q(x*)
        smp->Alawadhi_Copy_State(smp, &tpl);
        smp->Fixed_Dimension_Sampler(smp, smp->set->alawadhi_k, TAU|XI);            ' Produces proposed_trees during UpdateTopologies, but then deleted
        for( i=0; i<tpl->npartitions; i++ ) ill += tpl->part[i]->cPartialLogLikelihood; ' log l(x')
        ilp = smp->Log_Prior(smp, true);    ' log q(x')

    } else {
        tpl = pl;
    }

    ' Check for and handle possible memory allocation error
    pPartialLogLikelihood = (double *) malloc(sizeof(double)*tpl->npartitions);
    if( pPartialLogLikelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    
    ' Return index (>0) of partition with left end a topology changepoint
    propose_to_delete_index = tpl->Propose_Topology_Change_Point_To_Delete(tpl, smp->set, &selection_prob, &keep_on_left);
    if( tv->gmodel || debug>3 || global_debug>3 ) ProposedTreeStructureWithDeleteOne(tv, tpl, propose_to_delete_index, keep_on_left);   ' Allocates and sets proposed_trees

    ' |---left_part(A)---||---proposed_to_delete_part(B)---
    propose_to_delete_part = tpl->part[propose_to_delete_index];
    left_part = tpl->part[propose_to_delete_index - 1];

    ' Find left and right neighboring topology partitions
    left_top_change_index = propose_to_delete_index - 1;
    while( ! tpl->part[left_top_change_index]->topchange ) left_top_change_index--;
    right_top_change_index = propose_to_delete_index + 1;
    while( right_top_change_index < tpl->npartitions && ! tpl->part[right_top_change_index]->topchange ) right_top_change_index++;

    ' Choose tree and compute proposal probability
    log_number_tree_choices = tv->log_numTrees_minus2;
    if( keep_on_left ) {                                ' Keep tree on left                            WP: handled
        after_collapse_tree = left_part->ctree;
        if( right_top_change_index == tpl->npartitions )            ' Choose to delete right terminal segment              WP: handled
            log_number_tree_choices = tv->log_numTrees_minus1;      ' Choose new tree                          WP: 1/(numTrees - 1)
    } else {                                    ' Keep tree on right                           WP: handled
        after_collapse_tree = propose_to_delete_part->ctree;
        if( left_top_change_index == 0 )                    ' Choose to delete left terminal segment               WP: handled
            log_number_tree_choices = tv->log_numTrees_minus1;      ' Choose new tree                          WP: 1/(numTrees - 1)
    }

    ' Compute current and proposal likelihoods
    if( keep_on_left ) {
    '  <--left_top_change_index                                     propose_to_delete_index                           <----right_top_change_index
    ' ||---left_top_change_part--|--...--|--left_part(A:after_collapse_tree)--||--propose_to_delete_part(B)--|--...--||--right_top_change_part(C)--|--...--||--(D)

        ' Compute current and proposal likelihood on the right
        for(i = propose_to_delete_index; i < right_top_change_index; i++) {
            partition *cpart = tpl->part[i];
            if(compute_likelihood) pPartialLogLikelihood[i] = TreeLogLikelihood(after_collapse_tree, smp, cpart->cmatrix, cpart->counts, cpart->cHyperParameter, false);
            else pPartialLogLikelihood[i] = 0.0;
            
            pLogLikelihood += pPartialLogLikelihood[i];
            if( compute_likelihood ) cLogLikelihood += cpart->cPartialLogLikelihood;
        }
    }
    else {
    '  <--left_top_change_index                                       propose_to_delete_index                            <----right_top_change_index
    ' ||---left_top_change_part(C)--|--...--|--left_part(A:after_collapse_tree)--||--propose_to_delete_part(B)--|--...--||--right_top_change_part--|--...--||--(D)

        ' Compute current and proposal likelihood on the left
        for(i = left_top_change_index; i < propose_to_delete_index; i++) {
            partition *cpart = tpl->part[i];
            if(compute_likelihood) pPartialLogLikelihood[i] = TreeLogLikelihood(after_collapse_tree, smp, cpart->cmatrix, cpart->counts, cpart->cHyperParameter, false);
            else pPartialLogLikelihood[i] = 0.0;

            pLogLikelihood += pPartialLogLikelihood[i];
            if( compute_likelihood ) cLogLikelihood += cpart->cPartialLogLikelihood;
        }
    }

    if( debug>3 || global_debug>3 ) smp->Report_Proposed_State(smp, "DeleteOne", tv->proposed_trees, tpl->topology_changes);

    ' Compute prior ratio for \tau (based on proposed_trees)
    log_tau_prior_ratio = tv->Log_Prior_Ratio(tv, tpl->topology_changes + 1 - 1, true);

    ' Compute acceptance ratio
    priorRatio = - tpl->log_top_lambda              ' K
'      + log(tpl->topology_changes)                ' K (now cancels with \xi prior)
        + log(tpl->alignment_length - tpl->topology_changes)    ' \xi
        + log_tau_prior_ratio                   ' replacement \tau
'      + log(dcp->numTrees - 1)                ' old \tau
        ;
    birthProb = log(tpl->top_one_bkm1)              ' BIRTH: birth probability
        - log(tpl->alignment_length - tpl->topology_changes)    ' BIRTH: places to put new \xi during birth
        - log_number_tree_choices               ' BIRTH: choices for new \tau during birth
        - log_two                       ' BIRTH: KSD: choice if new \tau on left or right during birth (no longer cancels with deathProb elt)
        ;

    deathProb =
        + log(tpl->top_one_dk)                  ' DEATH: death probability
        + log(selection_prob)                   ' DEATH: choice of \xi to delete (new version)
'      + log(tpl->topology_changes) - log_two          ' DEATH: choice of \xi to delete and keep left or right (old version)
        ;

    ' q_{mn}(x,x') / q_{nm}(x',x) * pi(x) / pi(x')
    logRatio = birthProb - deathProb + priorRatio + pLogLikelihood - cLogLikelihood;

    if( smp->set->alawadhi_topology_one || smp->set->alawadhi_topology || smp->set->alawadhi ) {
        ' x* -> ... -> x' -> x
        ' Already included in ratio: pi(x) / pi(x') == pll*plp / ill / ilp
        ' pi(x') / pi(x*) * pi*(x*) / pi*(x')
        'double ologRatio = logRatio;
        logRatio += smp->set->alawadhi_factor*(cll + clp - ill - ilp) + ill + ilp - cll - clp;
        if( smp->set->alawadhi_debug && logRatio > log(0.2) ) {
'          fprintf(stderr, "ALAWADHI_DEBUG (%20s): %d move prob: %e; %e -> %e -> g%e; logRatio: %e -> %e\n", smp->move_names[DELETE_XI], pl->part[propose_to_delete_index]->left, (birthProb - deathProb), (clp+cll), (ilp+ill), (priorRatio+pLogLikelihood-cLogLikelihood)+ill+ilp, (ologRatio), (logRatio));
            /*if( logRatio > ologRatio ) fprintf(stderr, "+++++++\n");
            else fprintf(stderr, "-------\n");*/
        }
    }

    if( debug>1 || global_debug>1 )
        smp->Report_Proposal_Statistics(smp, "DeleteOne", priorRatio, birthProb, deathProb, pLogLikelihood, cLogLikelihood, propose_to_delete_index, (int) keep_on_left, 0);

    smp->tries[DELETE_XI]++;

    if( logRatio > 0 || smp->set->rng->nextStandardUniform(smp->set->rng) < exp(logRatio) ) {
    '  <--left_top_change_index                                     propose_to_delete_index                           <----right_top_change_index
    ' ||---left_top_change_part--|--...--|--left_part(A:after_collapse_tree)--||--propose_to_delete_part(B)--|--...--||--right_top_change_part(C)--|--...--||--(D)
        
        smp->acceptancerate[DELETE_XI]++;

        DeleteOneAccept(smp, tpl, tv, left_top_change_index, propose_to_delete_index, right_top_change_index, keep_on_left, after_collapse_tree, pPartialLogLikelihood, log_tau_prior_ratio);
        if( smp->set->alawadhi_topology_one || smp->set->alawadhi_topology || smp->set->alawadhi ) smp->Alawadhi_Accept(smp);

        if( debug>0 || global_debug>0 ) smp->Report_State(smp, "DeleteOne", logRatio, log_tau_prior_ratio);
        if( debug>3 || global_debug>3 || global_debug == -1 ) {
            VerifyLikelihood(smp, false);
            VerifyCounts(smp, "TopologyDeleteOne", false);
        }
    } else {
        if( smp->set->alawadhi_topology_one || smp->set->alawadhi_topology || smp->set->alawadhi ) {
            smp->Alawadhi_Reject(smp);
        } else if( smp->set->gmodel || debug>3 || global_debug>3 ) {
            if( tv->proposed_trees ) free(tv->proposed_trees);
            tv->proposed_trees = NULL;
        }
    }
    if(pPartialLogLikelihood) free(pPartialLogLikelihood);
}' TopologyDeleteOne

static void DeleteOneAccept(sampler *smp, partition_list *pl, topology_vector *tv, int left_top_change_index, int propose_to_delete_index, int right_top_change_index, boolean keep_on_left, tree *after_collapse_tree, double *pPartialLogLikelihood, double log_tau_prior_ratio) {
    int i;
    partition *propose_to_delete_part = pl->part[propose_to_delete_index];
    partition *left_part = pl->part[propose_to_delete_index-1];

    tv->log_prior_prob += log_tau_prior_ratio;
    if( smp->set->gmodel || debug>3 || global_debug>3 ) {
        if( tv->current_trees ) free(tv->current_trees);
        tv->current_trees = tv->proposed_trees;
        tv->proposed_trees = NULL;
    }

    if( keep_on_left ) {

        ' Update likelihoods and topology in the affected segments to right
        for(i = propose_to_delete_index + 1; i < right_top_change_index; i++) {
            partition *cpart = pl->part[i];
            cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
            cpart->ctree = after_collapse_tree;
        }

        ' Delete topology change point
        ' Change point must stay because it is also a parameter change point
        if( propose_to_delete_part->parchange ) {
            propose_to_delete_part->topchange = false;
            propose_to_delete_part->cPartialLogLikelihood = pPartialLogLikelihood[propose_to_delete_index];
            propose_to_delete_part->ctree = after_collapse_tree;
        }
        else {
            ' Extend left partition and remove proposed_to_delete_part entirely
            PartitionAddPartition(left_part, propose_to_delete_part);
            left_part->cPartialLogLikelihood += pPartialLogLikelihood[propose_to_delete_index];
            left_part->right = propose_to_delete_part->right;

            PartitionListRemovePartition(pl, propose_to_delete_index);
        }
    }
    else {
    
        ' Update intermediate likelihoods and topology on the left
        for(i = left_top_change_index; i < propose_to_delete_index; i++) {
            partition *cpart = pl->part[i];
            cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
            cpart->ctree = after_collapse_tree;
        }
    
        ' Delete topology change point
        ' Change point must stay because it was also a parameter change point
        if( propose_to_delete_part->parchange ) {
            propose_to_delete_part->topchange = false;
        }
        else {
            ' Extend left partition
            PartitionAddPartition(left_part, propose_to_delete_part);
            left_part->cPartialLogLikelihood += propose_to_delete_part->cPartialLogLikelihood;
            left_part->right = propose_to_delete_part->right;
    
            PartitionListRemovePartition(pl, propose_to_delete_index);
        }
    }
    pl->topology_changes--;
}' DeleteOneAccept

static void TopologyDeleteTwo(topology_vector *tv, partition_list *pl, sampler *smp) {
    const char *fxn_name = "TopologyDeleteTwo";
    int left_top_change_index;      ' first region on left with topology change
    int right_top_change_index;     ' first region on right with topology change
    tree *after_collapse_tree = NULL;   ' The topology that will spread in to replace old topology
    int left_to_delete_index = 0;       ' Region left of region to delete
    int right_to_delete_index = 0;      ' Region right of region to delete
    int number_xi_choices;          ' Number of places to put topology change points
    int i;
    double *pPartialLogLikelihood = NULL;
    double pLogLikelihood = 0.0, cLogLikelihood = 0.0, logRatio, priorRatio, deathProb, birthProb;
    partition *cpart;
    double selection_prob, log_tau_prior_ratio;
    boolean local_debug = 0;
    double cll=0.0, clp=0.0, ill=0.0, ilp=0.0;
    partition_list *tpl = NULL;


    ' Before we propose to move to another dimension, do some fixed dimension moves
    if( smp->set->alawadhi_topology_two || smp->set->alawadhi_topology || smp->set->alawadhi ) {
        ' x* -> ... -> x' -> x
        for( i=0; i<pl->npartitions; i++ ) cll += pl->part[i]->cPartialLogLikelihood;   ' log l(x*)
        clp = smp->Log_Prior(smp, false);   ' log q(x*)
        smp->Alawadhi_Copy_State(smp, &tpl);
        smp->Fixed_Dimension_Sampler(smp, smp->set->alawadhi_k, TAU|XI);
        for( i=0; i<tpl->npartitions; i++ ) ill += tpl->part[i]->cPartialLogLikelihood; ' log l(x')
        ilp = smp->Log_Prior(smp, true);    ' log q(x')
    } else {
        tpl = pl;
    }

    pPartialLogLikelihood = (double *) malloc(sizeof(double)*tpl->npartitions);

    if( pPartialLogLikelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    if( !tpl->Propose_Two_Topology_Change_Points_To_Delete(tpl, smp->set, &left_to_delete_index, &right_to_delete_index, &selection_prob) ) return;

    if( tv->gmodel || debug>3 || global_debug>3 || local_debug>3 ) ProposedTreeStructureWithDeleteTwo(tv, tpl, left_to_delete_index, right_to_delete_index);

    ' Find left and right topology change point
    left_top_change_index = left_to_delete_index - 1;
    while( left_top_change_index && ! tpl->part[left_top_change_index]->topchange ) left_top_change_index--;
    right_top_change_index = right_to_delete_index + 1;
    while( right_top_change_index < tpl->npartitions && ! tpl->part[right_top_change_index]->topchange ) right_top_change_index++;

    after_collapse_tree = tpl->part[left_top_change_index]->ctree;

    number_xi_choices = tpl->part[right_top_change_index-1]->right - tpl->part[left_top_change_index]->left - 1;    ' include right, not include left, not include other
    'number_xi_choices = (double)1/ (tpl->part[right_top_change_index-1]->right - tpl->part[left_to_delete_index]->left);  ' non-exclusive on one end
    'number_xi_choices += (double)1/ (tpl->part[right_to_delete_index]->left - tpl->part[left_top_change_index]->left - 1);    ' non-exclusive on both ends

    '  <--left_top_change_index                                     propose_to_delete_index                           <----right_top_change_index
    ' ||---left_top_change_part--|--...--|--left_part(A:after_collapse_tree)--||--propose_to_delete_part(B)--|--...--||--right_top_change_part(C)--|--...--||--(D)

    ' Compute current and proposal likelihood
    for(i = left_to_delete_index; i < right_to_delete_index; i++) {
        cpart = tpl->part[i];
        if(compute_likelihood) pPartialLogLikelihood[i] = TreeLogLikelihood(after_collapse_tree, smp, cpart->cmatrix, cpart->counts, cpart->cHyperParameter, false);
        else pPartialLogLikelihood[i] = 0.0;

        pLogLikelihood += pPartialLogLikelihood[i];
        cLogLikelihood += cpart->cPartialLogLikelihood;
    }

    log_tau_prior_ratio = tv->Log_Prior_Ratio(tv, tpl->topology_changes + 1 - 2, true);

    ' Compute acceptance ratio
    priorRatio = 'log(tpl->topology_changes)
        '+ log(tpl->topology_changes - 1)
        + log(tpl->alignment_length - tpl->topology_changes)
        + log(tpl->alignment_length - tpl->topology_changes + 1)
        + log_tau_prior_ratio
'      + 2*dcp->log_numTrees_minus1
        - 2*tpl->log_top_lambda             ' K: K -> K-2
        ;
    birthProb = log(tpl->top_two_bkm2)          ' BIRTH: birth probability
        - log(tpl->alignment_length - tpl->topology_changes + 1)    ' BIRTH: number of choices for first topology change point during birth L - 1 - (K-2) = L - K + 1
        - log(number_xi_choices)            ' BIRTH (uniform): number of choices for second topology change point during birth
'      + log(number_xi_choices)            ' BIRTH (equal left/right): _probability_ of selecting second topology change point (depends on who's first)
        - tv->log_numTrees_minus1           ' choices for tree in new region during birth
        + log_two
'      - log_two                   ' BIRTH: probability left is selected first (or right)
        ;
    deathProb =
        + log(tpl->top_two_dk)              ' death probability
        + log(selection_prob)               ' selection of these change points to delete during death
        ;

    ' q_{mn}(x,x') / q_{nm}(x',x) * pi(x) / pi(x')
    logRatio = birthProb - deathProb + priorRatio + pLogLikelihood - cLogLikelihood;

    if( smp->set->alawadhi_topology_two || smp->set->alawadhi_topology || smp->set->alawadhi ) {
        ' x* -> ... -> x' -> x
'      double ologRatio = logRatio;

        ' pi(x') / pi(x*) * pi*(x*) / pi*(x')
        logRatio += smp->set->alawadhi_factor*(cll + clp - ill - ilp) + ill + ilp - cll - clp;
        if( smp->set->alawadhi_debug && logRatio > log(0.2) ) {
'          fprintf(stderr, "ALAWADHI_DEBUG (%20s): (%d,%d), move prob: %e; %e -> %e -> g%e; logRatio: %e -> %e\n", smp->move_names[DELETE_TWO_XI], tpl->part[left_to_delete_index]->left, tpl->part[right_to_delete_index]->left, (birthProb - deathProb), (clp+cll), (ilp+ill), (priorRatio+pLogLikelihood-cLogLikelihood)+ill+ilp, (ologRatio), (logRatio));
            /*if( logRatio > ologRatio ) fprintf(stderr, "+++++++\n");
            else fprintf(stderr, "-------\n");*/
        }
    }

    if( debug>1 || global_debug>1 || local_debug>1 )
        smp->Report_Proposal_Statistics(smp, "DeleteTwo", priorRatio, birthProb, deathProb, pLogLikelihood, cLogLikelihood, left_to_delete_index, right_to_delete_index, 0);

    smp->tries[DELETE_TWO_XI]++;

    if( logRatio > 0.0 || smp->set->rng->nextStandardUniform(smp->set->rng) < exp(logRatio) ) {
        partition *first_propose_to_delete_part = tpl->part[left_to_delete_index];
        partition *left_of_propose_to_delete_part = tpl->part[left_to_delete_index - 1];
        partition *right_of_propose_to_delete_part = tpl->part[right_to_delete_index];
        partition *last_propose_to_delete_part = tpl->part[right_to_delete_index - 1];

        smp->acceptancerate[DELETE_TWO_XI]++;

        tv->log_prior_prob += log_tau_prior_ratio;
        if( smp->set->gmodel || debug>3 || global_debug>3 || local_debug>3 ) {
            if( tv->current_trees ) free(tv->current_trees);
            tv->current_trees = tv->proposed_trees;
            tv->proposed_trees = NULL;
        }

        ' Update intermediate likelihoods and topology
        for( i = left_to_delete_index; i < right_to_delete_index; i++ ) {
            cpart = tpl->part[i];
            cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
            cpart->ctree = after_collapse_tree;
        }

        ' Delete RIGHT topology change point
        ' The next region on the right is a parameter change point so segment must stay
        if( right_of_propose_to_delete_part->parchange ) {
            right_of_propose_to_delete_part->topchange = false;
            ' Next region is not parameter change point, so will go
        } else {
            ' Extend left partition
            PartitionAddPartition(last_propose_to_delete_part, right_of_propose_to_delete_part);
            last_propose_to_delete_part->cPartialLogLikelihood += right_of_propose_to_delete_part->cPartialLogLikelihood;
            last_propose_to_delete_part->right = right_of_propose_to_delete_part->right;

            PartitionListRemovePartition(tpl, right_to_delete_index);
        }

        ' Delete LEFT topology change point
        ' The first partition in deleted region is a parameter change point and must stay
        if( first_propose_to_delete_part->parchange ) {
            first_propose_to_delete_part->topchange = false;
        }
        ' The first partition (of more than one) in deleted region was only topology change point and will go
        else {
            ' Extend left partition
            PartitionAddPartition(left_of_propose_to_delete_part, first_propose_to_delete_part);
            left_of_propose_to_delete_part->cPartialLogLikelihood += first_propose_to_delete_part->cPartialLogLikelihood;
            left_of_propose_to_delete_part->right = first_propose_to_delete_part->right;

            PartitionListRemovePartition(tpl, left_to_delete_index);
        }
        tpl->topology_changes -= 2;

        if( smp->set->alawadhi_topology_two || smp->set->alawadhi_topology || smp->set->alawadhi ) smp->Alawadhi_Accept(smp);

        if( debug>0 || global_debug>0 || local_debug>0 ) smp->Report_State(smp, "DeleteTwo", logRatio, log_tau_prior_ratio);
        if( debug>3 || global_debug>3 || global_debug==-1 ) {
            VerifyLikelihood(smp, false);
            VerifyCounts(smp, "TopologyDeleteTwo", false);
        }
    } else {
        if( smp->set->alawadhi_topology_two || smp->set->alawadhi_topology || smp->set->alawadhi ) {
            smp->Alawadhi_Reject(smp);
        } else {
            if( smp->set->gmodel || debug>3 || global_debug>3 || local_debug>3 ) {
                if( tv->proposed_trees ) free(tv->proposed_trees);
                tv->proposed_trees = NULL;
            }
        }
    }
    if(pPartialLogLikelihood) free(pPartialLogLikelihood);
}' TopologyDeleteTwo


static void Update_Topologies(topology_vector *tv, partition_list *pl, sampler *smp, boolean alawadhi) {
    const char *fxn_name = "Update_Topologies";
    boolean local_debug = false;
    int i, j;
    partition *curr_part;       ' Current part
    tree *pTree;            ' Proposed tree
    tree *prev_tree;        ' Tree of previous region, if any
    tree *next_tree;        ' Tree of next region, if any
    int curr_index=0;       ' Index of current topology change point
    int next_index = 0;     ' Index of next topology change point
    int prev_index = 0;     ' Index of previous topology change point
    double *pPartialLogLikelihood = (double *) malloc(sizeof(double)*pl->npartitions);
    double pLogLikelihood = 0.0;    ' Likelihood of region affected by proposed change
    double cLogLikelihood = 0.0;    ' Current likelihood of region affected by proposed change
    double log_tau_prior_ratio = 0.0;
    double logRatio = 0.0;
       
    if( pPartialLogLikelihood == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    if( smp->set->gmodel || debug>3 || global_debug>3 ) {
        tv->proposed_trees = (tree **)malloc(sizeof(tree *)*(pl->topology_changes+1));
        if( tv->proposed_trees == NULL ) {
            fprintf(stderr, "%s: memory allocation error\n", fxn_name);
            exit(-1);
        }
        for( i=0; i<(pl->topology_changes+1); i++ ) tv->proposed_trees[i] = tv->current_trees[i];
    }

    curr_index = 0;
    j = 0;
    do {
        next_index = curr_index + 1;
        while( next_index < pl->npartitions && !pl->part[next_index]->topchange ) next_index++;
        ' WAS INCONSEQUENTIAL BUG: recompute prev_index

        prev_tree = pl->part[prev_index]->ctree;
        curr_part = pl->part[curr_index];
        next_tree = next_index<pl->npartitions ? pl->part[next_index]->ctree : NULL;

        pTree = tv->Propose_New_Tree(tv, smp->set, prev_tree, curr_part->ctree, next_tree);

        if( pTree != NULL ) {   ' There is a new tree we could propose for this region

            if( smp->set->gmodel || debug>3 || global_debug>3 ) tv->proposed_trees[j] = pTree;

            if( local_debug || debug>3 || global_debug>3 ) smp->Report_Proposed_State(smp, "UpdateTopology", tv->proposed_trees, pl->topology_changes + 1);

            log_tau_prior_ratio = tv->Log_Prior_Ratio(tv, pl->topology_changes + 1, true);

            pLogLikelihood = 0.0;
            cLogLikelihood = 0.0;

            for( i=curr_index; i<next_index; i++ ) {
                partition *cpart = pl->part[i];

                if(compute_likelihood) pPartialLogLikelihood[i] = TreeLogLikelihood(pTree, smp, cpart->cmatrix, cpart->counts, cpart->cHyperParameter, false);
                else pPartialLogLikelihood[i] = 0.0;

                pLogLikelihood += pPartialLogLikelihood[i];
                if(compute_likelihood) cLogLikelihood += cpart->cPartialLogLikelihood;
            }
            logRatio = (alawadhi ? smp->set->alawadhi_factor : 1.0) * (pLogLikelihood - cLogLikelihood + log_tau_prior_ratio);
            if( !alawadhi ) smp->tries[UPDATE_TAU]++;
            if( local_debug || debug>1 || global_debug>1 )
                smp->Report_Proposal_Statistics(smp, "UpdateTopology", log_tau_prior_ratio, 0.0, 0.0, pLogLikelihood, cLogLikelihood, j, curr_part->ctree->tree_index, pTree->tree_index);

            if( logRatio > 0 || smp->set->rng->nextStandardUniform(smp->set->rng) < exp(logRatio) ) {
                if( !alawadhi ) smp->acceptancerate[UPDATE_TAU]++;

                tv->log_prior_prob += log_tau_prior_ratio;
                if( smp->set->gmodel || debug>3 || global_debug>3 ) tv->current_trees[j] = pTree;
        
                ' Update all partitions
                for( i = curr_index; i < next_index; i++ ) {
                    partition *cpart = pl->part[i];
                    cpart->ctree = pTree;
                    cpart->cPartialLogLikelihood = pPartialLogLikelihood[i];
                }
                if( local_debug>0 || debug>0 || global_debug>0 ) smp->Report_State(smp, "UpdateTopology", logRatio, log_tau_prior_ratio);
                if( debug>3 || global_debug>3 || global_debug==-1 ) {
                    VerifyLikelihood(smp, false);
                    VerifyCounts(smp, "Update_Topologies", false);
                }
            } else {
                if( smp->set->gmodel || debug>3 || global_debug>3 ) tv->proposed_trees[j] = tv->current_trees[j];
            }
        }

        j++;    ' Next topology
        prev_index = curr_index;
        curr_index++;
        while( curr_index < pl->npartitions && !pl->part[curr_index]->topchange ) curr_index++;
    } while( curr_index < pl->npartitions );

    if(pPartialLogLikelihood) free(pPartialLogLikelihood);
    if( smp->set->gmodel || debug>3 || global_debug>3 ) {
            if( tv->proposed_trees ) free(tv->proposed_trees);
        tv->proposed_trees = NULL;
    }
}' Update_Topologies

static double Reverse_Add_One_Log_Proposal_Probability(topology_vector *tv, const partition_list *pl, int insert_index, boolean keep_left, tree *nTree) {
    const char *fxn_name = "Reverse_Add_One_Log_Proposal_Probability";
    int new_topology_changes = pl->topology_changes + 1;
    int curr_index, next_index, i, nchoices;
    boolean local_debug = false;

    tv->proposed_trees = (tree **) malloc(sizeof(tree *)*(new_topology_changes+1));

    if( tv->proposed_trees == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    i = 0;
    curr_index =  1;
    while( curr_index < pl->npartitions && !pl->part[curr_index]->topchange ) curr_index++;
    if( insert_index < curr_index ) {
        if( keep_left ) {
            tv->proposed_trees[i++] = pl->part[0]->ctree;
            if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
            tv->proposed_trees[i++] = nTree;
            if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
        } else {
            tv->proposed_trees[i++] = nTree;
            if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
            tv->proposed_trees[i++] = pl->part[0]->ctree;
            if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
        }
    } else {
        tv->proposed_trees[i++] = pl->part[0]->ctree;
        if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
    }

    ' Addition:                             *
    ' |---prev_prev---|---prev---|---curr---|---new----|---next---|
    ' |---prev_prev---|---prev---|---new----|---curr---|---next---|
    while( curr_index < pl->npartitions ) {
        next_index = curr_index + 1;
        while( next_index < pl->npartitions && !pl->part[next_index]->topchange ) next_index++;
        if( insert_index < curr_index || insert_index >= next_index ) {
            tv->proposed_trees[i++] = pl->part[curr_index]->ctree;
            if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
        } else if( keep_left ) {
            tv->proposed_trees[i++] = pl->part[curr_index]->ctree;
            if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
            tv->proposed_trees[i++] = nTree;
            if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
        } else if( !keep_left ) {
            tv->proposed_trees[i++] = nTree;
            if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
            tv->proposed_trees[i++] = pl->part[curr_index]->ctree;
            if(local_debug) fprintf(stderr, "%d", tv->proposed_trees[i-1]->tree_index);
        }
        curr_index = next_index;
    }

    if( i!= new_topology_changes+1 ) {
        fprintf(stderr, "%s: debugging write access violation\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    nchoices = 2;                       ' end top change points can always go
    for( i=1; i<new_topology_changes; i++ ) {
        if( tv->proposed_trees[i-1] != tv->proposed_trees[i+1] ) nchoices += 2; ' see below
    }
    ' Delete 1, keep_left, delete 2, !keep_left are possible given the same condition A != C
    '         1      2
    ' |---A---|------|---C---|------|

    if( nchoices == 0 ) exit(EXIT_FAILURE);

    if( !tv->gmodel && debug<=3 && global_debug<=3 && tv->proposed_trees ) {
        free(tv->proposed_trees);
        tv->proposed_trees = NULL;
    }

    return (double) 1 / nchoices;
}' Reverse_Add_One_Log_Proposal_Probability

static double Reverse_Add_Two_Log_Proposal_Probability(topology_vector *tv, const partition_list *pl, int insert_left_index, tree *nTree) {
    const char *fxn_name = "ProposeTwoTopChangePointDeleteProbability";
    int i, nchoices, curr_index, next_index;

    tv->proposed_trees = (tree **) malloc(sizeof(tree *)*(pl->topology_changes+3));

    if( tv->proposed_trees == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    i = 0;
    tv->proposed_trees[i++] = pl->part[0]->ctree;
    curr_index = 1;
    while( curr_index < pl->npartitions && !pl->part[curr_index]->topchange ) curr_index++;
    if( insert_left_index < curr_index ) {
        tv->proposed_trees[i++] = nTree;
        tv->proposed_trees[i++] = pl->part[0]->ctree;
    }

    while( curr_index < pl->npartitions ) {
        next_index = curr_index + 1;
        while( next_index < pl->npartitions && !pl->part[next_index]->topchange) next_index++;
        tv->proposed_trees[i++] = pl->part[curr_index]->ctree;
        if( insert_left_index >= curr_index && insert_left_index < next_index ) {
            tv->proposed_trees[i++] = nTree;
            tv->proposed_trees[i++] = pl->part[curr_index]->ctree;
        }
        curr_index = next_index;
    }

    nchoices = 0;
    for( i=1; i<pl->topology_changes+2; i++ ) if( tv->proposed_trees[i-1] == tv->proposed_trees[i+1] ) nchoices++;

    if( !tv->gmodel && debug<=3 && global_debug<=3 && tv->proposed_trees ) {
        free(tv->proposed_trees);
        tv->proposed_trees = NULL;
    }

    if( nchoices == 0 ) exit(EXIT_FAILURE);

    return (double) 1 /nchoices;
}' Reverse_Add_Two_Log_Proposal_Probability

int Propose_Topology_Change_Point_To_Delete(const partition_list *pl, settings *set, ... ) {
    const char *fxn_name = "Propose_Topology_Change_Point_To_Delete";
    int *choices = (int *) malloc(sizeof(int)*2*(pl->topology_changes+1));
    int curr_index, prev_index, prev_prev_index, next_index;
    int nchoices=0, returned_choice, last_one;
    double *prob;
    boolean *keep_left;
    va_list vargs;

    if( choices == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    va_start(vargs, set);
    prob = va_arg(vargs, double *);
    keep_left = va_arg(vargs, boolean *);
    va_end(vargs);

    curr_index = 1;
    while( curr_index < pl->npartitions && !pl->part[curr_index]->topchange ) curr_index++;

    ' Deletion:
    ' |---prev_prev---|---prev---|---curr---|---next---|
    ' |---prev_prev---|<-------------curr---|---next---|   keep right
    ' |---prev_prev---|---prev------------->|---next---|   keep left
    while( curr_index < pl->npartitions ) {
        prev_index = curr_index - 1;
        while( prev_index > 0 && ! pl->part[prev_index]->topchange ) prev_index--;
        prev_prev_index = prev_index ? prev_index - 1 : prev_index;
        while( prev_prev_index > 0 && ! pl->part[prev_prev_index]->topchange ) prev_prev_index--;
        next_index = curr_index + 1;
        while( next_index < pl->npartitions && ! pl->part[next_index]->topchange ) next_index++;

        if( prev_prev_index < prev_index && (pl->part[prev_prev_index]->ctree != pl->part[curr_index]->ctree) )
            choices[nchoices++] = curr_index;
        else if( prev_prev_index == prev_index )
            choices[nchoices++] = curr_index;
        if( next_index < pl->npartitions && (pl->part[prev_index]->ctree != pl->part[next_index]->ctree) )
            choices[nchoices++] = -curr_index;
        else if( next_index == pl->npartitions )
            choices[nchoices++] = -curr_index;
        last_one = curr_index;
        curr_index = next_index;
    }
    
    *prob = (double) 1 / nchoices;
    
    returned_choice = choices[ (int) (set->rng->nextStandardUniform(set->rng)*nchoices) ];
    free(choices);

    *keep_left = false;
    if( returned_choice < 0 ) *keep_left = true;
    returned_choice = abs(returned_choice);

    return(returned_choice);

}' Propose_Topology_Change_Point_To_Delete

boolean Propose_Two_Topology_Change_Points_To_Delete(const partition_list *pl, settings *set, ... ) {'int *left, int *right, double *prob) {
    const char *fxn_name = "Propose_Two_Topology_Change_Points_To_Delete";
    int *left_pair = (int *) malloc(sizeof(int)*pl->topology_changes);
    int *right_pair = (int *) malloc(sizeof(int)*pl->topology_changes);
    int curr_index, prev_index, next_index, nchoices=0, choice;
    int *left, *right;
    double *prob;
    va_list vargs;

    if( left_pair == NULL || right_pair == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    va_start(vargs, set);
    left = va_arg(vargs, int *);
    right = va_arg(vargs, int *);
    prob = va_arg(vargs, double *);
    va_end(vargs);

    curr_index = 1;
    while( curr_index < pl->npartitions && !pl->part[curr_index]->topchange ) curr_index++;

    while( curr_index < pl->npartitions ) {
            prev_index = curr_index - 1;
            while( prev_index && !pl->part[prev_index]->topchange ) prev_index--;
            next_index = curr_index + 1;
            while( next_index<pl->npartitions && !pl->part[next_index]->topchange ) next_index++;
            if( next_index == pl->npartitions) break;
            if( pl->part[prev_index]->ctree == pl->part[next_index]->ctree ) {
                left_pair[nchoices] = curr_index;
                right_pair[nchoices++] = next_index;
            }
            curr_index = next_index;
    }

    if( ! nchoices ) {
        if(left_pair) free(left_pair);
        if(right_pair) free(right_pair);
        return false;
    }

    choice = (int) (set->rng->nextStandardUniform(set->rng)*nchoices);
    *left = left_pair[choice];
    *right = right_pair[choice];
    *prob = (double) 1 / nchoices;
    if( left_pair ) free(left_pair);
    if( right_pair ) free(right_pair);
    return true;
}' Propose_Two_Topology_Change_Point_To_Delete

static int ProposedTreeStructureWithDeleteOne(topology_vector *tr, const partition_list *pl, int delete_index, boolean keep_left) {
    const char *fxn_name = "ProposedTreeStructureWithDeleteOne";
    int i, j=0, next_index;

    tr->proposed_trees = (tree **) malloc(sizeof(tree *)*pl->topology_changes);

    if( tr->proposed_trees == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    for( i=0; i<pl->npartitions; i++ ) {
        partition *cpart = pl->part[i];
        if( !cpart->topchange || ( keep_left && i == delete_index ) ) continue;
        next_index = i+1;
        while( next_index < pl->npartitions && !pl->part[next_index]->topchange ) next_index++;
            if( !keep_left && next_index == delete_index ) continue;
        tr->proposed_trees[j++] = cpart->ctree;
    }
    if( j != pl->topology_changes ) {
        fprintf(stderr, "%s: debugging memory access violation\n", fxn_name);
        exit(EXIT_FAILURE);
    }
    return -1;
}' ProposedTreeStructureWithDeleteOne

static int ProposedTreeStructureWithDeleteTwo(topology_vector *tr, const partition_list *pl, int left, int right) {
    const char *fxn_name = "ProposedTreeStructureWithDeleteTwo";
    int i = 0, curr_index;

    tr->proposed_trees = (tree **) malloc(sizeof(tree *)*(pl->topology_changes-1));

    if( tr->proposed_trees == NULL ) {
        fprintf(stderr, "%s: memory allocation error\n", fxn_name);
        exit(EXIT_FAILURE);
    }

    tr->proposed_trees[i++] = pl->part[0]->ctree;
    curr_index = 1;
    while( curr_index < pl->npartitions && ! pl->part[curr_index]->topchange ) curr_index++;

    while( curr_index < pl->npartitions ) {
        if( curr_index < left ) tr->proposed_trees[i++] = pl->part[curr_index]->ctree;
        else if( curr_index > right ) tr->proposed_trees[i++] = pl->part[curr_index]->ctree;
        curr_index++;
        while( curr_index < pl->npartitions && ! pl->part[curr_index]->topchange ) curr_index++;
    }
    return -1;
}' ProposedTreeStructureWithDeleteTwo

void TopologyVectorDelete(topology_vector *tv, boolean delete_tree_list) {
    if( !tv ) return;
    if( tv->current_trees ) free(tv->current_trees);
    if( tv->top_prior ) free(tv->top_prior);
    if( tv->tree_list && delete_tree_list ) {
        free(tv->tree_list[0]);
        free(tv->tree_list);
    }
    free(tv);
    tv = NULL;
}' TopologyVectorDelete


'Constants.H
'#ifndef __CONSTANTS
'#define __CONSTANTS

' Define boolean data type
typedef unsigned char boolean;
'#define true 1
'#define false 0

' Constants
enum {SCP_RECOMB, DIVERGE, DCP_RECOMB}; ' possible models
enum {HKY};             ' possible ctmc models
enum {ALPHA, KAPPA};            ' possible parameterizations of ctmc models
enum {ADD_RHO, DELETE_RHO, ADD_XI,  ' 11 move types
    DELETE_XI, ADD_TWO_XI, DELETE_TWO_XI, FIXED_DIMENSION, UPDATE_TAU, UPDATE_KAPPA_AND_MU, UPDATE_XI, UPDATE_RHO};

static const int KAPPA_MU = 1;
static const int TAU = 2;
static const int XI = 4;
static const int RHO = 8;

static const int MAX_TREE_STRING = 100;             ' maximum length of string representation of phylogenetic trees
static const int MAX_LINE_LENGTH = 200;             ' maximum length of line in cmdfile
static const double log_two = 0.693147180559945309417232121;    ' log(2)
static const double normal_const = 0.9189385;           ' ln(2*pi)/2
static const double tolerance = 1e-5;

' Predefine "classes"
typedef struct _sampler sampler;
typedef struct _dcpsampler dcpsampler;
typedef struct _qmatrix qmatrix;
typedef struct _tree tree;
typedef struct _partition partition;
typedef struct _partition_list partition_list;
typedef struct _branch branch;
typedef struct _topology_vector topology_vector;
typedef struct _seqdata seqdata;
typedef struct _node node;
typedef struct _rngen rngen;

' Global variables
extern int global_debug;        ' Override or change without compiling using cmdfile option debug:
extern boolean compute_likelihood;  ' Override or change without compiling using cmdfile option compute_likelihood:

' Global debugging functions (in dcpsampler.c)
boolean VerifyLikelihood(sampler *, boolean);
boolean VerifyCounts(sampler *, const char *, boolean);


''#End If
cpsampler.H
'#ifndef __CPSAMPLER
'#define __CPSAMPLER

#include <math.h>
#include "constants.h"
#include "tree.h"
#include "seqdata.h"
#include "settings.h"
#include "sampler.h"
#include "partition.h"
#include "partition_list.h"
#include "ihkynoboundfixpimatrix.h"

'typedef struct {
Type cpsampler
    'seqdata *sqd;       '
    'settings *set;      ' Object with all setting information
    'sampler *smp;
    
     sqd As seqdata       '
     setx As settings      ' Object with all setting information
     smp As sampler
    
     cLogLikelihood As Double   ' Current loglikelihood

    'tree *start_tree;   ' Fixed parental tree
    'tree *PostTree;     ' List of possible trees
    'int numTrees;       ' Number of possible trees
    'int cSameT;     ' Number of change points that are parameter-only change points
    'int npartitions;    ' Number of segments
    'partition **part_list;
    
    start_tree As tree   ' Fixed parental tree
    PostTree As tree     ' List of possible trees
    numTrees As Long       ' Number of possible trees
    cSameT As Long     ' Number of change points that are parameter-only change points
    npartitions As Long    ' Number of segments
    part_list() As partition
    
    ' MCMC parameters
    lenWindow As Long
    lambda As Double
    mix1 As Double
    logLambda As Double
    logTwo As Double
    sigmaA As Double
    sigmaM As Double
    logWSameT As Double
    logWNotSameT As Double
    wSameT As Double
    wNotSameT As Double
    C As Double
    bk As Double
    dk As Double
    bkm1 As Double
    dkp1 As Double

    jumpClasses As Byte

    'int *indexSeq;      ' Map from alignment position to unique pattern index
    indexSeq() As Long      ' Map from alignment position to unique pattern index
    lenSeq As Long     ' Length of alignment
    lenunique As Long      ' Number of unique patterns in data
End Type ' cpsampler;

void CPSamplerSetup(cpsampler **, seqdata *, settings *, char *);
void CPSamplerDelete(cpsampler **);



'#End If
'dcpsampler.h
'#ifndef __DCPSAMPLER
'#define __DCPSAMPLER

#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "constants.h"
#include "tree.h"
#include "tree_vector.h"
#include "evol_param.h"
#include "settings.h"
#include "partition_list.h"
#include "sampler.h"
#include "ihkynoboundfixpimatrix.h"


'struct _dcpsampler {
Public Type dcpsampler
    'sampler *smp;           ' Common parts of all samplers
    'settings *set;          ' Settings (pointer for easier access)
    smp As sampler      ' Common parts of all samplers
    set As settings        ' Settings (pointer for easier access)
    
    mc As Double           ' Global probability of monophyletic tree through alignment (see gmodel)

    'partition_list *part_list;  ' K, J, \xi, \rho
    'topology_vector *tree_vec;  ' \tau
    'partition_list *alawadhi_part_list;
    'topology_vector *alawadhi_tree_vec;
    
    part_list As partition_list  ' K, J, \xi, \rho
    tree_vec As topology_vector  ' \tau
    alawadhi_part_list As partition_list
    alawadhi_tree_vec As topology_vector
End Type 'dcpsampler

' Setup
void DCPSamplerSetup(dcpsampler **, seqdata *, settings *, char *);

' Output
void PrintTopologies(const sampler *, const char *, boolean, boolean);

' Cleanup
void DCPSamplerDelete(dcpsampler *);

'#End If

'discrete_gamma.h
' discrete_gamma.h
'
'#ifndef __discrete_gamma
'#define __discrete_gamma

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"

double LnGamma (double alpha);
double PointNormal (double prob);
double Getch(double ch, double xx, double g, double p, double aa, double c);
double PointChi2 (double prob, double v);
int DiscreteGamma(double *rK, double alfa, double beta, int K, boolean median);
double IncompleteGamma (double x, double alpha, double LnGamma_alpha);

'#End If

'evol_param.h
'#ifndef __EVOL_PARAM
'#define __EVOL_PARAM

/*
 * These functions are concerned with updating evolutionary parameters and adding/deleting
 * evolutionary change points.  They heavily on partition_list for proposing change points
 * to add or delete, qmatrix and derivatives for proposing new CTMC evolutionary parameters,
 * and branch for proposing new \mu parameters.
 *
 * If we want to consider more complicated relationships among parameters, structures would
 * be required to hold information.  See tree_vector.c for inspiration.
 */

#include "partition.h"
#include "partition_list.h"
#include "sampler.h"

' Function declarations
void UpdateParameters(partition_list *, sampler *, boolean);    ' Update parameters in each parameter segment, one segment at a time
void ParChgPtAdd(partition_list *, sampler *);          ' Propose and accept via MH step a new parameter change point
void ParChgPtDelete(partition_list *, sampler *);       ' Propose and accept via MH step to delete an existing parameter change point

'#End If

'ihkynoboundfixpimatrix.h
'#ifndef __IHKYNOBOUNDFIXPIMATRIX
'#define __IHKYNOBOUNDFIXPIMATRIX

/* These structures and functions are concerned with handling the HKY CTMC model
 * of evolution.
 *
 * There is also a structure to hold information about the hierarchical prior.
 * One file-wide instance of this is all that is needed since HKY parameters are
 * conditionally iid in the hierarchical structure.
 */

#include <stdarg.h>
#include "qmatrix.h"
#include "seqdata.h"
#include "settings.h"
#include "sampler.h"    ' Basic probability calculations

'typedef struct {
Public Type ihkynoboundfixpimatrix
    cached_ep As Double      ' Cached evolutionary parameter
    cached_avg_brlen As Double     ' Cached average branchlength
    alpha As Double            ' Easy place to store alpha in ALPHA version
    kappa As Double            ' Easy place to store kappin in KAPPA version
    model As Long          ' choice of ALPHA or KAPPA parameterized model
    pR As Double           ' Prob. of purine
    pY As Double           ' Prob. of pyrimidine
    b As Double          ' Auxiliary variable
    C As Double
    d As Double
End Type ' ihkynoboundfixpimatrix;

'typedef struct {
Public Type hky_kappa_hierarchical_prior
    kappa_mean As Double       ' Hierarchical mean on HKY parameters \kappa
    kappa_variance As Double      ' Hierarchical variance on HKY parameters \kappa
End Type '} hky_kappa_hierarchical_prior;

' Global function declarations:  these are the bare minimum to pull ourselves up by our bootstraps

' Construction: various constructors are needed for different situations (see dcpsampler.c)

' Construct default matrix with reasonable, but fixed parameter values: make default HKY
' matrix object in 1st arg qmatrix ptr, assuming 2nd arg model parameterization.
' The matrix can be further specified with calls to Set_Alpha or Set_Kappa.
void iHKYNoBoundFixPiMatrixMakeDefault(qmatrix **, int);
' Construct initially-distributed matrix: make HKY matrix object in 1st arg qmatrix ptr,
' asuming 2nd arg model parameterization, then initialize using data from 3rd arg seqdata,
' and default initial distribution on parameters.
void iHKYNoBoundFixPiMatrixMakeInitial(qmatrix **, int, const seqdata *, settings *);
' Construct fully-specified matrix: make HKY matrix object in 1st arg qmatrix ptr, assuming
' 2nd arg model parameterization, 3rd arg parameter vector, and 4th arg equilibrium distribution
void iHKYNoBoundFixPiMatrixMakeAndSet(qmatrix **, int, const double *, const double *);

' Setup: main class should call these along with constructors to get everything started
' Most subsequent calls to qmatrix functions are via the function pointers (see qmatrix.h)

' Set 1st arg HKY qmatrix to use ALPHA parameterization and set parameter alpha to 2nd arg
void Set_Alpha(qmatrix *, double);
' Set 1st arg HKY qmatrix to use KAPPA parameterization and set parameter kappa to 2nd arg
void Set_Kappa(qmatrix *, double);
' Do any global initializations required (such as setting up dependencies, like hierarchical prior on multiple qmatrix objects)
void iHKYNoBoundFixPiMatrixGlobalInitialize(settings *);

'#End If
multiseqdata.H
'#ifndef __MULTISEQDATA
'#define __MULTISEQDATA

#include "seqdata.h"


' used for DIVERGE; data_sets contain only data array;
' seq_sets contains counts, map, rmap but is not used;
' data will be compressed in terms of alignment for all sequences, and counts, map, etc are the same for each' set

'typedef struct {
Public Type multiseqdata
    num_datasets As Long     ' number of alignments (clusters)
    total_taxa As Long       ' number of sequences of all the alignments

    seq_sets As seqdata ' counts, map, rmap for each alignment

End Type '} multiseqdata;

void Delete_MultiStdData(multiseqdata *msd);

int Get_num_datasets(multiseqdata *);
int **Get_Data_set (multiseqdata *, int);
void Make_MultiStdData(multiseqdata *, seqdata *, int, int *);

void PrintMultiStdDataInfo (multiseqdata *);

'#End If

'node.h
'#ifndef __NODE
'#define __NODE

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "constants.h"
#include "settings.h"
#include "qmatrix.h"
#include "tree.h"


struct _node {
    node *up;
    node *right;
    node *left;
    node *brother;
    int id;
    int uid;
    int state;
    boolean is_branch;
    boolean *is_likelihood_done;
    boolean is_root;
    double *clikelihood;
    double branch_length;
    int descending_clade;
    int ascending_clade;
    boolean *descending_clades;
    boolean *ascending_clades;
};

' Node construction
node *Make_Subtree(node *, char *, char *, int *, double **, boolean **, node *, const boolean, const int);
node *Make_Subtree_From_Tree(const tree *, node *, int *, double **, boolean **, node *);
void Make_Node(node *, node *, node *, node *, const int, const int, const double, double *, boolean *);
double Read_Branch_Length(char *);
int Find_Branch_Split(char *);

' Simulation
void Simulate_Down_Branch(node *, node *, double, qmatrix *, rngen *);

' Auxiliary
void toString(char *, const node *, boolean);
boolean isLeftChild(const node *);
node *Brother(const node);
int Number_Nodes(const node);
int Number_Children(const node);
int Number_Leaves(const node *);
int CountDown(const node);
int Balance (node *);
void Cladify(node *, const int, const int, const int *);
void Cladify_Up(node *, const int, const int, const int *);
void Clear_Clades(node *);
boolean Clade_Ancestor(node *, int, int, int);
boolean Clade_Ancestor_Exists(node *, int, int, int);

' Output
void PrintNodesDFS(const node);

' Node cleanup
void Clear_Likelihood(node *);
void Delete_Node(node *);

'#End If

'partition.h
'#ifndef __PARTITION
'#define __PARTITION

#include "constants.h"
#include "qmatrix.h"
#include "seqdata.h"

struct _partition {
    int left;               ' Left end (inclusive) of a partition/segment
    int right;              ' Right end (inclusive) of a partition/segment
    boolean topchange;          ' Is left end point a topology change point?
    boolean parchange;          ' Is left end point a parameter change point?
    double cPartialLogLikelihood;       ' Current log likelihood of this partition
    double cHyperParameter;         ' Current expected divergence $\mu$ in this partition
    double cPartialLogHyperParameterPrior;  ' Current log prior density of $\mu$ for this partition
    qmatrix *cmatrix;           ' Current evolutionary matrix for this partition
    tree *ctree;                ' Pointer to current tree in this region
    int *counts;                ' Sufficient statistics
    int lenunique;              ' Number of unique patterns (length of counts)

    ' Only used by legacy cpsample code
    boolean doUpdate;           ' Whether to update this segment (legacy: cpsampler)
    boolean doXiUpdate;         ' Update its left boundary (legacy: cpsampler)
    int cTree;              ' Current topology for this partition (legacy: cpsampler)
};

' Partition creation::
void PartitionMake(partition **, int, int, int, boolean, boolean);  ' Create a partition with specified values of lenunique, left, right, topchange, parchange
void PartitionMakeCopy(partition **, const partition *);        ' Create a partition by copying another partition contents
void PartitionCopy(partition *, partition *);               ' Copy one partition onto another (memory assumed allocated)
void PartitionReset(partition *, int, int, boolean, boolean);       ' Reset partition left, right, topchange, and parchange

' Partition count creation/setting::
void PartitionCopyCounts(partition *, const seqdata *);         ' Set partition data to that of seqdata object
void PartitionCopySegmentCounts(partition *, const seqdata *, int, int);' Set partition data to that of args left, right-defined segment of seqdata alignment
void PartitionAddPartition(partition *, const partition *);     ' Supplement partition (arg 1) data with partition (arg 2) data
void PartitionSubtractPartition(partition *, const partition *);    ' Subtract data of partition (arg 2) from partition (arg 1) data

' Set partition data to be that of argument 2 data minus argument 3 data (the first partition should contain the second partition or invalid data can result; no checks)
void PartitionCopyPartitionCountDifferences(partition *, const partition *, const partition *);
' Set partition data to be that of argument 2 data minus argument 3 and 4 data (the first partition should contain the next two partitions; no checks)
void PartitionCopyPartitionCountDifferences2(partition *, const partition *, const partition *, const partition *);
' Set partition data to be that of argument 2 data plus argument 3 data (the first partition should be union of the other two partitions)
void PartitionCopyPartitionSum(partition *, const partition *, const partition *);

' Partition elimination::
void PartitionDelete(partition *);                  ' Delete the partition (does NOT delete matrix since there's no guarantee it's allocated)
void PartitionDeleteMatrix(partition *);                ' Delete the matrix belonging to this (arg 1) partition

'#End If
 'partition_list.h
'#ifndef __PARTITION_LIST
'#define __PARTITION_LIST

/* Information and functions for the number of topology change points K, the number of parameter change points J,
 * the location of topology change points \xi, and the location of parameter change points \rho.
 *
 * Since there is only one prior (poisson) so far used for both K and J, we do not establish
 * the machinery (functions and function pointers) to compute prior ratios via function calls.
 * The reason is it is faster and easier to verify log ratios if the computations are not
 * hidden in function calls.  However, some clues as to how such function pointers might be
 * established are included in comments.
 *
 * The same argument applies to the prior on topology and parameter change point locations.
 */

#include "constants.h"
#include "partition.h"
#include "sampler.h"
#include "tree_vector.h"

struct _partition_list {
    partition **part;                   ' List of partitions
    int npartitions;                    ' Current number of partitions
    int topology_changes;                   ' Current number of topology change points
    int parameter_changes;                  ' Current number of parameter change points

    ' Prior constants
    double top_lambda;                  ' Prior expected number of topology change points (set in settings)
    double log_top_lambda;                  ' Pre-computed log
    double top_lambda_squared;              ' Pre-computed square
    double par_lambda;                  ' Prior expected number of parameter change points (set in settings)
    double log_par_lambda;                  ' Pre-computed log
    int alignment_length;                   ' L

    ' Move probabilities:
    double top_one_bk, top_two_bk, top_one_dk, top_two_dk;                  ' Topology change point move probabilities
    double par_bk, par_dk;                                  ' Parameter change point move probabilities
    double top_one_bkm1, top_two_bkm2, top_one_dkp1, top_two_dkp2, par_dkp1, par_bkm1;  ' Reverse move probabilities

    ' Function pointers:  currently not necessary to use function pointers since there is only one way to propose change points to add/delete.
    ' However, there already exist multiple mechanisms for these proposals, though they are not coded here (see Vladimir's java code).
    ' TODO: to be really useful these functions must also compute proposal probabilities (topology change point delete functions do).
    int (*Propose_Parameter_Change_Point_To_Delete) (const partition_list *, settings *);           ' Propose an existing parameter change point to delete
    int (*Propose_Topology_Change_Point_To_Delete) (const partition_list *, settings *, ...);       ' Propose an existing topology change point to delete
    boolean (*Propose_Two_Topology_Change_Points_To_Delete) (const partition_list *, settings *, ...);  ' Propose an existing contiguous topolog segment to delete
    int (*Propose_Change_Point) (const partition_list *, settings *, boolean);              ' Propose location of new change point to add
    int (*Propose_Second_Change_Point) (const partition_list *, settings *, int, int, int, boolean);    ' Propose the location of a second topology change point (in AddTwo)
    void (*Update_Move_Probabilities) (partition_list *, const settings *);                 ' Set move probabilities based on current state

    ' TODO: ideas for functions that will allow easy coding of alternative priors
    double (*K_Log_Prior) ();       ' Computes prior for the number K of topology change points (don't use ratio because no log_prior member)
    double (*J_Log_Prior) ();       ' Computes prior for the number J of parameter change points (don't use ratio because no log_prior member)
    double (*Xi_Log_Prior) ();
    double (*Rho_Log_Prior) ();
};

' Global Functions

' Construction:
void PartitionListMake(partition_list **, const sampler *, int);    ' Make a partition with length given by 3rd argument
void PartitionListMakeCopy(partition_list **, const partition_list *);  ' Makes a complete copy of a partition_list

' Manipulation:
void PartitionListAddPartition(partition_list *, partition *, int); ' Insert new partition (2nd argument) into list at index location of 3rd argument
void PartitionListRemovePartition(partition_list *, int);       ' Remove partition at index position 2nd argument

' Information:
int PartitionContaining(const partition_list *, int, int, int);     ' Identify the partition in inclusive range (3rd and 4th arguments) that contains location (2nd argument)

' MCMC Functions: these may be better as function pointers for flexibility in changing MCMC or coding multiple MCMC in same program

/* function: propose and possibility accept (via MH step) new change point locations (4th argument indicates
 * what kind of change point in multiple change point model)
 */
void UpdateChangePointLocations(partition_list *, sampler *, boolean, boolean);

/* function: Propose new change point different from the current location (2nd argument) within a window
 * surrounding the current change point location, but not beyond the boundary positions (3rd and 4th arguments).
 * This declaration is here because cpsampler uses it.
 */
int ProposeNewChangePointPosition(sampler *, int, int, int);

void PartitionListDelete(partition_list *, boolean);            ' Boolean is to delete (or not) the associated qmatrix objects

'#End If

'qmatrix.h
'#ifndef __QMATRIX
'#define __QMATRIX

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "settings.h"

struct _qmatrix {
    ' Data:
    int nchars;             ' Number of states in this CTMC model
    int nvariables;             ' Number of variables in v
    double log_prior;           ' Log prior probability of matrix
    double *v;              ' Evolutionary parameters
    double *pi;             ' Stationary distribution
    double **cached_qmatrix;        ' Cached matrix to speed up calculations
    double **cached_qmatrix_transpose;  ' Cached matrix to speed up calculations
    void *derived_mt;           ' Pointer to derived matrix (model-specific)
    ' Functions:
    void (*Matrix_Update_Cache) (qmatrix *, ...);                       ' Function updates transition probabilities of 1st arg qmatrix
    double (*Matrix_Proposer) (qmatrix **, const qmatrix *, settings *);            ' Function proposes new 1st arg qmatrix based on 2nd arg qmatrix; returns proposal probability
    ' Function proposes new 1st & 2nd arg qmatrices from 3rd qmatrix arg; returns proposal prob
    double (*Matrix_Propose_Split) (qmatrix **, qmatrix **, const qmatrix *, settings *, ...);
    ' Function proposes a new 1st arg qmatrix from 2nd and 3rd qmatrix args; returns proposal prob
    double (*Matrix_Propose_Merge) (qmatrix **, const qmatrix *, const qmatrix *, settings *, ...);
    void (*Matrix_Copy) (qmatrix *, const qmatrix *);                   ' Function copies 2nd qmatrix arg to 1st qmatrix arg
    void (*Matrix_Make_Copy) (qmatrix **, const qmatrix *);                 ' Function makes its 1st arg qmatrix as a copy of its 2nd argument qmatrix
    void (*Matrix_Sync) (qmatrix *);                            ' Function to sync internal state of qmatrix
    void (*Matrix_Delete) (qmatrix *);                          ' Function for deleting its qmatrix argument
    double (*Matrix_Log_Prior) (qmatrix *, ...);                        ' Function computes log prior of 1st arg qmatrix
};

void QMatrixMake(qmatrix **, int, int);     ' Function: makes 1st arg qmatrix on 2nd arg state space size (e.g. 4 for nucs) with 3rd arg free parameters (e.g. 1 HKY)
void QMatrixCopy(qmatrix *, const qmatrix *);   ' Function: copy one matrix onto another (assumes both are allocated already)
void QMatrixDelete(qmatrix *);          ' Function: deletes qmatrix passed in at 1st arg
' Used by derived matrix objects to update the base object
void QMatrixUpdateParameters(qmatrix *, const double *, const double *);

'#End If
sampler.H
'#ifndef __SAMPLER
'#define __SAMPLER

#include <math.h>
#include "constants.h"
#include "settings.h"
#include "seqdata.h"
#include "branch.h"

'struct _sampler {

Public Type sampler
   ' int JumpNumber;         ' Current iteration number
   ' settings *set;          ' Settings
   ' seqdata *sqd;           ' Sequence data
   ' rngen *rng;         ' Random number generator (or wrapper to one)
   ' branch *br;         ' Structure for branches
   ' void *derived_smp;      ' Pointer to derived sampler
    
    JumpNumber As Long         ' Current iteration number
    set As settings          ' Settings
    sqd As seqdata           ' Sequence data
    rng As rngen          ' Random number generator (or wrapper to one)
    br As branch         ' Structure for branches
    derived_smp As sampler      ' Pointer to derived sampler
    
    ' Information about moves
    'int nmoves;         ' Number of different moves
    'int max_move_name_length;   ' Max length of a move name
    'int *tries, *acceptancerate;    ' Information about move attempts and acceptance rates
    'char **move_names;      ' User-friendly names for each move type for reporting results

    'int sincePrint;         ' Number of iterations since last print
    
    nmoves As Long         ' Number of different moves
    max_move_name_length As Long    ' Max length of a move name
    tries() As Long
    acceptancerate() As Long    ' Information about move attempts and acceptance rates
    move_names() As String      ' User-friendly names for each move type for reporting results

    sincePrint As Long         ' Number of iterations since last print
    
    
    ' Output
    fout As String         ' Output file

    'void (*OutputLine) (const sampler *);       ' Function to prepare one line of output
    'double (*logJacobian) (const sampler *, ...);   ' Overridable function pointer for computation of jacobian in change point add/delete
    'void (*run)(sampler *);             ' Run MCMC
    'void (*Matrix_Make_Default) (qmatrix **, int);  ' Make default matrix of given type
    'void (*Matrix_Make_Initial) (qmatrix **, int, const seqdata *, settings *); ' Make matrix of given type using available data
    'void (*Matrix_Make_and_Set) (qmatrix **, int, const double *, const double *);  ' Make matrix of given type with parameters provided
    'void (*Fixed_Dimension_Sampler) (sampler *, int, int);  ' As name indicates
    'void (*Alawadhi_Copy_State) (sampler *, ...);   ' Al-Awadhi: copy state to temporary state vector to allow updates without permanently accepting
    'void (*Alawadhi_Accept) (sampler *);        ' Al-Awadhi: accept temporary state vector
    'void (*Alawadhi_Reject) (sampler *);        ' Al-Awadhi: reject temporary state vector and restore initial state
    'double (*Log_Prior) (const sampler *, boolean);
    
    ' DEBUG:
    'void (*Report_State) (sampler *, const char *, double, double);                             ' Report current state (after an accept)
    'void (*Report_Proposed_State) (sampler *, const char *, tree **, int);                          ' Report proposed state (before accept or reject)
    'void (*Report_Proposal_Statistics) (sampler *, const char *, double, double, double, double, double, int, int, int);    ' Report information about a proposal (before accept/reject)
End Type 'sampler

static const double logOneOverSqrtTwoPi = -0.9189;

' Auxiliary
double logStandardNormalDensity(const double);
double logNormalDensity(const double, const double, const double);

'void ChooseStartingPi(const seqdata *, double *);

boolean LogMHAccept(const double, const double);

' Information about move types:
void SamplerSetNumberMoves(sampler *, int);
void SamplerAddMoveName(sampler *, int, const char *);
void SamplerSaveEstimates(const sampler *, int);
void CloseSampler(const sampler *);
void SamplerMake(sampler **);

' Clean up:
void SamplerDelete(sampler *);

'#End If
'seqdata.h
'#ifndef __SeqData
'#define __SeqData

'#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "constants.h"
#include "sequence.h"
#include "settings.h"
#include "tree.h"
#include "ihkynoboundfixpimatrix.h"

/**
* This is an abstract super class for the alignment data object where positions
* are iid.
*
* It declares all generic variables that are present regardless of the exact
* data actually used.  For example, the length of the alignment, the number of
* unique patterns in the alignment, and the sufficient statistics (counts of
* patterns) are stored here.
*/

' ntaxa by cluster was added so seqdata can be used like MultiSeqData in Gu.java

struct _seqdata {
    int ntaxa;          ' The total number of taxa in the alignment
    int lenseq;     ' The length of the alignment
    int lenunique;      ' The number of unique patterns in the alignment
    int num_chars;      ' The number of distinct characters in alignment
    int *counts;        ' The number of times each unique pattern appears in alignment

    sequence *alignment;    ' The minimally processed sequence alignment
    int *map;       ' Maps the original site location to compressed site location
    int *rmap;      ' Maps the compressed site location to AN original site location
'  void (*Compress_Data) (seqdata *);
    int **data;     ' Data for this class is a matrix of integers
};

' Creation
void Make_SeqData(seqdata **, sequence *, const int, const int, const int);
void Set_SeqData(seqdata *, const int, int[], int[], int[]);
void Set_SeqDataFromSeq(seqdata *, int, seqdata *, sequence *);
void ReadPhylip(seqdata **, const char *, boolean, int);    ' Creates a seqdata object, allocating memory via Make_SeqData
void SimulateAlignment(seqdata **, settings *);                     ' Creates a seqdata object by simulation, allocating memory via Make_SeqData

' Deletion
void SeqDataDelete(seqdata *);

' Utility
void Uncompressed_Alignment_Composition(const seqdata *, double **);    ' Composition function before data compressed
void Alignment_Composition(const seqdata *, double **);
FILE *ReadPhylipHeader(const char *, int *, int *);

' Screen I/O (debug)
void PrintSequences(const seqdata *);
void PrintSortedSequences(const seqdata *);
void PrintCompressedSequences(const seqdata *);
void PrintCounts(const seqdata *);
void PrintMap(const seqdata *);
void PrintRMap(const seqdata *);
void Print_Summary(const seqdata *);
void Print_Data(const seqdata *);
void Print_Distances(const seqdata *);

'#End If
'sequence.h
'#ifndef __Sequence
'#define __Sequence

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "constants.h"

'typedef struct {
'    int length;
'    char *strand;
'    char *name;
'    int *data;
'    int *count;  ' array for the number of of each character in the sequence; last element in array is count of -9
'} sequence;

Type sequence
    length As Long
    strand() As Byte
    namex() As String
    data() As Long
    count() As Long  ' array for the number of of each character in the sequence; last element in array is count of -9
End Type ' sequence;


'#define NAME_LENGTH 10



'void Setup_Sequence(sequence *, int instrand, const int); ' Needed?
'void Setup_Sequence(sequence *);              ' Needed?

void Setup_Sequence_From_Data(sequence *, int *, int, char *);  ' Setup an existing sequence obj given integer representatio of data (arg 2), length (arg 3), name (arg 4)
void Setup_Sequence(sequence *, char *, char *, int *, int);    ' Setup an existing sequence object given sequence in char * (arg 2), name (arg 3), integer representation (arg 4), model (arg 4)

void PrintSequenceInfo(sequence *, boolean);
void Copy_Sequence(sequence *, const sequence);
void SequenceDelete(sequence *);
int getBase(sequence, const int);
double Sequence_Composition(sequence, const int);
double Sequence_Absolute_Composition(sequence, const int);
int To_Int_DNA(char);
int To_Int_AA(char);

void toAAString(sequence *, char *);
    
'#End If
settings.H
'#ifndef __Settings
'#define __Settings

'#define _GNU_SOURCE

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "constants.h"

enum {INT, DOUBLE, STRING}; ' possible data types

' setting variables unique to RECOMB
typedef struct {

    /* Mixing Parameters */
    'double sdPi;      ' Standard deviation for nucleotide frequency updater
    double mix2;        ' not used
    int lenWindow;      ' Mixing parameter: Changepoint moves are uniformly proposed within a window of this size
    
} set_recomb;

typedef struct {
    double lambda;      ' Total expected number of changepoints
} set_scp;

' settings variables unique to dcp model
typedef struct {
    double par_lambda;  ' Prior expected number of parameter change points
    double top_lambda;  ' Prior expected number of topology change points
} set_dcp;

struct _rngen {
    gsl_rng *rengine;   ' Random number generator

    ' Random number generation functions
    double (*nextNormal) (rngen *, double);     ' Function that produces the next N(0,sd), where sd is arg 1
    double (*nextStandardUniform) (rngen *);        ' Function that produces the next standard uniform deviate
    double (*nextDispersedUniform) (rngen *);       ' Function that produces the next dispersed uniform deviate
    double (*nextStandardNormal) (rngen *);     ' Function that produces the next standard normal deviate
    double (*nextGamma) (rngen *, double, double);
    double (*nextExponential) (rngen *, double);

    ' For debug purposes:
    boolean useRnList;  ' Use a pre-fabricated list of random numbers
    double *rnList;     ' The list of random numbers
    int current_rn;     ' The current position in the list of random numbers
    int total_rn;       ' The total number of numbers in the random list
};


' settings variables unique to DIVERGE
'typedef struct {
Type set_diverge
    ' for Gu only
    mix3 As Double
    sdAlpha As Double
    sdTheta As Double
    alpha As Double
    Theta As Double
    M As Double
    numMu As Long
    'double *Mu;
    mu() As Double
    Mu_Default As Double
    
End Type ' set_diverge;

'typedef struct {
Type set_sim
    segments As Long   ' Number of segments
    'double *mu; ' Average branch length
    mu() As Double  ' Average branch length
    'double *kappa;  ' Ti/Tv ratio
    kappa() As Double  ' Ti/Tv ratio
    'double *pi; ' Stationary distribution
    pi() As Double ' Stationary distribution
    'char **tree;    ' Simulation tree
    tree() As Byte    ' Simulation tree
    'int *length;    ' Length of alignment
    length()   As Long    ' Length of alignment
    total_length As Long   ' Total length of alignment
End Type ' set_sim;


'typedef struct {
Type settings
    model As Long              ' current model (see constants.h for list)
    ctmc_model As Long          ' HKY is only option
    ctmc_parameterization As Long   ' Parameterization, when multiple possible, of CTMC model (e.g. ALPHA vs. KAPPA)

    '/* Run Parameters */
    length As Long      ' The length of the MCMC run
    burnin As Long      ' The length of the MCMC burnin period
    subsample As Long       ' The MCMC subsampling rate (sample very ____)

    'char *init_string;  ' String of output with which to initialize chain
    init_string() As Byte   ' String of output with which to initialize chain
    
    num_pTrees As Long     ' the number of parental trees (clusters)
    pTree()   As Byte       ' Fixed parental tree(s) stored as string

    debug As Long
    alawadhi_debug As Long
    compute_likelihood As Byte
    exit_condition As Byte
    report_iact As Byte

    C As Double       ' rjMCMC mixing parameter: The multiplier of the birth/death probabilities
    sdMu As Double        ' Standard deviation for average branch length $\mu$ updater (recomb)
    sdEP As Double        ' Standard deviation for evolutionary parameter updater
    weight As Double      ' log probability that two adjacent segments have the same topology
    jump_classes As Byte   ' If set, the number of changepoints will be estimated
    add_xi As Byte     ' If set, add/delete topology change points
    add_rho As Byte    ' If set, add/delete parameter change points
    gmodel As Byte
    simulate_data As Byte  ' If set, generate data by simulation according to other arguments
    log_mono_prob As Byte

    'set_recomb *recomb; ' Recomb settings
    'set_diverge *diverge;   ' Gu settings
    'set_dcp *dcp;       ' Dual change point settings
    'set_scp *scp;       ' Single change point settings
    'set_sim *sim;   ' Simulation settings
    
    recomb As set_recomb ' Recomb settings
    diverge As set_diverge    ' Gu settings
    dcp As set_dcp        ' Dual change point settings
    scp As set_scp       ' Single change point settings
    sim As set_sim   ' Simulation settings
    
    cmdfile_seed As Byte
    'unsigned long int seed;
    seed As Long
    rng As rngen

    ' Hyperparameters
    mu_hyper_mean As Double
    mu_hyper_variance As Double
    titv_hyper_mean As Double
    titv_hyper_variance As Double
    update_hyperparameters As Byte

    ' Al-Awadhi
    alawadhi As Byte       ' If set, use Al-Awadhi mechanism for increasing dimension jumps (turn alawadhi on across the board)
    alawadhi_topology As Byte  ' If set, use Al-Awadhi mechanism for increasing dimension jumps for K (number of topology change points)
    alawadhi_topology_one As Byte  ' If set, use Al-Awadhi mechanism for increasing dimension jumps for K (number of topology change points)
    alawadhi_topology_two As Byte  ' If set, use Al-Awadhi mechanism for increasing dimension jumps for K (number of topology change points)
    alawadhi_parameter As Byte ' If set, use Al-Awadhi mechanism for increasing dimension jumps for J (number of parameter change points)
    alawadhi_k As Long         ' TODO: add different k and factor for parameter vs. topology (first get parameter working!)
    alawadhi_factor As Double

    ' What are these?
    sdP As Double     ' Transition kernel standard deviation for updating population-level stationary distribution
    dNP As Double        ' Transition kernel standard deviation to update pseudo-counts on population-level stationary distribution
    sdT1 As Double        ' Transition kernel standard deviation to update one branch length at a time
    sdT2 As Double        ' Transition kernel standard deviation to update all branch lengths at a time
    sdHyperEP As Double   ' Transition kernel standard deviation to update population-level evolutionary parameter
    sdUV As Double        ' Transition kernel standard deviation to update population-level evolutionary parameters
    mix1 As Double        ' Mixing parameter
    sigmaAlpha As Double  ' rjMCMC mixing parameter: The standard deviation for generating a new $\alpha$ during a birth step
    sigmaMu As Double     ' rjMCMC mixing parameter: The standard deviation for generating a new $\mu$ during a birth step

End Type '} settings;

' Main interface function for reading options
void ReadCmdfile(settings **, const char *);

' Random number generation
void Set_Seed(rngen *, const unsigned long int);

' Setup
void Set_Defaults(settings *);
void Set_DCP_Defaults(settings *);
void Set_Diverge_Defaults(settings *);
void Set_Recomb_Defaults(settings *);
void UnknownOption(char *, char *);

' Verification of cmdfile options
void CheckBoundsGr0(const double, const char *);
void CheckBoundsGrEqual0(const double, const char *);
void CheckBounds01(const double, const char *);

' Delete
void Settings_Cleanup(settings *);

' Not used...
/*
boolean SettingsParseInitialParameters(const settings *, int *, double **, double **, char ***, int **);
void Settings_ParseLastLine(settings *, int *, int *, int *, char ***, double **, double **, double **, double **, double **, double **, double **, int **);
*/

'#End If
'tree.h
'#ifndef __TREE
'#define __TREE

#include <math.h>
#include <stdlib.h>
#include "node.h"
#include "seqdata.h"
#include "qmatrix.h"
#include "settings.h"
#include "constants.h"
#include "sampler.h"


struct _tree {
Public Type tree
    nleaves As Long
    nnodes As Long
    tree_index As Long         ' Unique identifier for each tree (set externally; used externally to id trees)
    'int nbranches;
    nchars As Long
    is_likelihood_done() As Byte
    is_likelihood_done_blank() As Byte
    has_branches As Byte
    'node *root;
    'node *node_list;
    'node **leaf_list;
    'double *likelihood;
    root As node
    node_list As node
    leaf_list()   As node
    likelihood() As Double
End Type 'tree

' Tree construction
void Make_Tree(tree **, char *, int);
void EnumerateLastTaxon(tree **, const tree *);
void GrowLeaf(tree *, const int, const boolean);

' Likelihood calculation
double TreeLogLikelihood(tree *, const sampler *, qmatrix *, const int *, const double, boolean);
double SiteLikelihood(tree *, const seqdata *, qmatrix *, boolean, int);
double CalcLikelihood(tree *, const seqdata *, qmatrix *, const int *, const double, boolean, int, boolean);

' Simulation
void Simulate_Position(tree *, qmatrix *, const double *, double, rngen *);

' Auxiliary
boolean SameTrees(tree *, tree *, boolean);
void Balance_Tree(tree *);
double SumOfBranchLengths(tree *);
int Number_Parental_Trees(tree *);
int Number_All_Trees(tree *, const int);
boolean Monophyletic(tree *, int *);

' Proposers
void JointBranchAndTopology(char *, tree *, settings *, int, double);

' Output
void PrintBrothers(tree *);
void PrintTreeInfo(tree *);
void PrintLeafList(tree *);

' Cleanup
void TreeDelete(tree *);

'#End If
tree_vector.H
'#ifndef __TREE_VECTOR
'#define __TREE_VECTOR

/* Structures and functions holding information and computing things related to the $\tau$ vector, i.e. vector of topologies.
 * topology_vector: basic structure storing information about allowable topologies, number of topologies, and pre-computed logs
 * topology_gmodel_prior: place to store monophyletic/nonmonophyletic prior auxiliary data
 */

#include <stdarg.h>
#include "constants.h"
#include "partition_list.h"
#include "sampler.h"
#include "seqdata.h"
#include "tree.h"

' See Fang's dissertation for information about this prior
typedef struct {
    boolean *monophyletic;          ' Vector of indicator variables corresponding to tree_list and indicating whether tree is monophyletic
    int num_nonmono_trees;          ' The total number of nonmonophyletic trees
    int num_mono_trees;         ' The total number of monophyletic trees (with above should sum to numTrees)
    double log_mono_prob;           ' Pre-computed log
    double log_num_nmtrees[2];      ' Pre-computed log
    double log_num_nmtrees_minus1[2];   ' Pre-computed log
} topology_gmodel_prior;

Type topology_vector
'struct _topology_vector {
    ' Information about possible trees:
    'tree **tree_list;           ' List of possible trees
    tree_list() As tree
    numTrees As Long               ' The total number of possible trees
    start_tree As tree           ' The fixed parental tree (queries excluded)
    'tree **current_trees;           ' The list of current topologies along the alignment (i.e. \tau)
    current_trees() As tree
    'tree **proposed_trees;          ' A list of proposed topologies along the alignment
    proposed_trees() As tree
    ' Pre-computed numbers for faster computing
    log_numTrees As Double            ' Pre-computed log
    log_numTrees_minus1 As Double     ' Pre-comptued log
    log_numTrees_minus2 As Double     ' Pre-computed log

    ' Prior information:
    log_prior_prob As Double          ' Log prior probability of current state
    void *top_prior;            ' (optional) Pointer to a structure with information about the prior
    gmodel As Byte             ' Whether the gmodel prior is being used

    ' Function pointers:
    ' Draw a starting tree from the assume initial distribution (allows
    ' multiple possible initial distributions or even dependencies among
    ' \tau elements in initial distribution)
    
    'tree * (*Draw_Initial_Tree) (const topology_vector *, settings *);
    (*Draw_Initial_Tree) (const topology_vector *, settings *) as tree
    ' Propose a new tree where the three tree objects are supposed to be
    ' tree to left, current tree, and tree to right, but NULL allowed.
    ' Returns NULL if there are no allowable proposals.
    tree * (*Propose_New_Tree) (const topology_vector *, settings *, tree *, tree *, tree *);
    ' Propose and accept via MH step a new topology change point
    void (*Add_One) (topology_vector *, partition_list *, sampler *);
    ' Used to compute probability of reverse move to Add_One given arg 3
    ' is new topology change point, arg 4 indicates whether new tree (arg
    ' 4) is on left or right.  May wish to use variable argument list to
    ' handle more complex situations.
    double (*Reverse_Add_One_Log_Proposal_Probability) (topology_vector *, const partition_list *, int, boolean, tree *);
    ' Propose and accept via MH step a new topology segment (two
    ' neighboring topology change points).
    void (*Add_Two) (topology_vector *, partition_list *, sampler *);
    ' Used to compute probability of reverse move to Add_Two given arg 3
    ' index of segment to be spliced by new topology segment and the new
    ' tree (arg 4).
    double (*Reverse_Add_Two_Log_Proposal_Probability) (topology_vector *, const partition_list *, int, tree *);
    ' Propose and accept via MH step to delete an existing topology
    ' change point.
    void (*Delete_One) (topology_vector *, partition_list *, sampler *);
    'double (*Reverse_Delete_One_Log_Proposal_Probability) ();
    ' Propose and accept via MH step to delete an existing topology
    ' segment.
    void (*Delete_Two) (topology_vector *, partition_list *, sampler *);
    'double (*Reverse_Delete_Two_Log_Proposal_Probability) ();
    ' Propose and accept via MH step new topologies for all topology partitions,
    ' one at a time, in existing partition list.
    void (*Update_Topologies) (topology_vector *, partition_list *, sampler *, boolean);
    ' Compute prior ratio assuming (3rd argument) topology partitions in
    ' proposed structure (it knows the number of topology partitions in
    ' current structure) for some priors.  May rely on proposed_trees and
    ' current_trees for some priors with more complex dependencies among
    ' \tau elements.
    double (*Log_Prior_Ratio) (const topology_vector *, int, ...);
End Type 'struct _topology_vector

' GLOBAL FUNCTIONS PRE-DECLARATIONS:

' Construction:
void TreeVectorMake(topology_vector **, const sampler *);
void TreeVectorMakeCopy(topology_vector **, const topology_vector *, int);
void TreeVectorInitialize(topology_vector *, const partition_list *);

' MCMC: the following are used in updates that partition_list handles. Some
' aspects of the proposal require knowledge about the whole \tau vector and
' constraints thereon (i.e. no identical neighbors) so they are handled here.
' TODO: partition_list should have naive versions of these functions that are
' overwritten by specialized functions, such as these, in particular applications
int Propose_Topology_Change_Point_To_Delete(const partition_list *, settings *, ...);       ' Propose a topology change point to delete
boolean Propose_Two_Topology_Change_Points_To_Delete(const partition_list *, settings *, ...);  ' Propose a topology segment to delete

void TopologyVectorDelete(topology_vector *, boolean);                      ' Boolean is to delete (or not) the tree_list

'#End If

'Makefile
# Compiler/Linker variables
CC = gcc
VERSION = 1#
# -DDBG
ifdef DBG
    CFLAGS = -G - fno - inline
Else
    CFLAGS = -O3
End If

# Suggested gcc flags
CFLAGS=-Wall -W -O3 -ffast-math -fexpensive-optimizations -malign-double -march=i486 -funroll-all-loops -finline-functions
# For debugging
#CFLAGS=-g -fno-inline
# For debugging with gprof
#CFLAGS=-g -fno-inline -pg
LDFLAGS = -lm - lgsl - lgslcblas

# Local variables
srcs = $(wildcard *.c)
objs = $(srcs:.c=.o)
deps = $(srcs:.c=.d)

cbrother : $(objs)
    $(CC) $(CFLAGS) -o cbrother $(objs) $(LDFLAGS)

include $(deps)

%.d : %.c
    -@$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) $< \
        | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@'

.PHONY : backup restore clean public

backup:
    tar czvf cbrother.tar.gz *.c *.h Makefile TODO README.local INSTALL cmdfile

restore:
    tar xzvf cbrother.tar.gz

clean:
    rm *.o *.d cbrother

public:
    tar czvf cbrother$(VERSION).tar.gz *.c *.h Makefile INSTALL README examples
INSTALL
1. Install GSL (GNU Scientific Library) available at www.gnu.org/software/gsl which is required for random number generation.
2. Then, starting with the cbrother$VERSION.tar.gz file, issue these commands in order
    tar xzvf cbrother$VERSION.tar.gz    # Uncompress and dump archive contents into the current directory (no subdirectories)
    make                    # Compile the code
3. You may wish to move (install) the cbrother executable to a central directory in your path.

'README
Overview:
cbrother is software for inferring recombination when recombination is rare.  This is a c version of the code originally written in Java and available elsewhere (http:'www.biomath.medsch.ucla.edu/msuchard/).  Unlike the Java version, this version does not estimate the hierarchical parameters.

Features:
* Estimate recombination using the single multiple change point (SMCP) model (see reference 4).
* Estimate recombination using the dual multiple change point (DMCP) model (see reference 2).
* Simulate data under the model of evolution assumed by the software.
* In development: Test recombinants with similar mosaic structure (e.g. circulating recombinant forms) for evidence that they result from multiple recombination events, rather than a single recombination event (see reference 1).

Installation:
See the INSTALL file for information on installation.  Contact K. S. Dorman at http:'www.biomath.org/dormanks if you have trouble or questions.

Manual: At this time there is no manual.  Please read the rest of this document for some hints.  Also, see the example command files provided, as they contain important comments.

Usage Summary:
* Prepare an alignment in interleaved PHYLIP format containing reference sequences (representing putative parental sequences involved in the recombination event) and the putative recombinant sequences(s).
* Prepare a command file (use cmdfile as a template) for your data analysis.
* Determine the phylogenetic relationship of the reference sequences without branch lengths.  Enter this tree as the value of the cmdfile option start_tree.  The name of each reference sequence in this tree is its 0-index position in the PHYLIP file.
* Issue the command:
    cbrother $RANDOM <cmdfile> <phylip_file> <post_file>

Examples: This distribution comes with two examples.
1. AB073841 is a recombinant HBV sequence.  To run it type:
    cbrother $RANDOM examples/cmdfile.AB073841 examples/AB073841.phy examples/AB073841.post
2. A simulation example where three parameter change points and no topology change points are simulated.  To run it type:
    cbrother $RANDOM examples/cmdfile.simulation examples/simulation.post

References:
1. F. Fang, M. A. Suchard, K. S. Dorman (2005) Distinguishing multiple vs. single recombination events in a collection of similar mosaic structures. In preparation.
2. V. N. Minin, K. S. Dorman, M. A. Suchard. (2005) Dual multiple change-point model leads to more accurate recombination detection. Bioinformatics. Accepted.
3. M. A. Suchard, R. E. Weiss, K. S. Dorman, J. S. Sinsheimer. (2003) Inferring spatial phylogenetic variation along nucleotide sequences. JASA. 98:427-437.
4. M. A. Suchard, R. E. Weiss, K. S. Dorman, J. S. Sinsheimer. (2002) Oh brother where Art thou? A bayes factor test for recombination with uncertain heritage.  Systematic Biology. 51(5): 1-14.
examples/
examples/AB073841.phy
 7 3221
consD_fnl    CACCTCTGCC TAATCATCTC TTGTTCATGT CCTACTGTTC AAGCCTCCAA
consA_fnl    CACCTCTGCC TAATCATCTC TTGTACATGT CCCACTGTTC AAGCCTCCAA
consC_fnl    CACCTCTGCC TAATCATCTC ATGTTCATGT CCTACTGTTC AAGCCTCCAA
consB_fnl    CACCTCTGCC TAGTCATCTC TTGTTCATGT CCTACTGTTC AAGCCTCCAA
consF_fnl    CACCTCTGCC TAATCATCTT TTGTTCATGT CC-ACTGTTC AAGCCTCCAA
consH_fnl    CACCTCTGCC TAATCATCTT TTGTTCATGT CCCACTGTTC AAGCCTCCAA
AB073841     CACCTCTGCC TAATCATCTC ATGTTCATGT CCTACTGTTC AAGCCTCCAA

             GCTGTGCCTT GGGTGGCTTT GGGGCATGGA CATTGACC-T TATAAAGAAT
             GCTGTGCCTT GGGTGGCTTT GGGGCATGGA CATTGACCCT TATAAAGAAT
             GCTGTGCCTT GGGTGGCTTT GGGGCATGGA CATTGACCCG TATAAAGAAT
             GCTGTGCCTT GGGTGGCTTT AGG-CATGGA CATTGACCCT TATAAAGAAT
             GCTGTGCCTT GGGTGGCTTT GGGGCATGGA CATTGACCCT TATAAAGAAT
             GCTGTGCCTT GGGTGGCTTT GGGGCATGGA CATTGACCCT TATAAAGAAT
             GCTGTGCCTT GGGTGGCTTT AGGGCATGGA CATTGACCCG TATAAAGAAT

             TTGGAGCTAC TGTGGAGTTA CTCTCGTTTT TGCCTTCTGA CTTCTTTCCT
             TTGGAGCTAC TGTGGAGTTA CTCTCGTTTT TGCCTTCTGA CTTCTTTCCT
             TTGGAGCTTC TGTGGAGTTA CTCTCTTTTT TGCCTTCTGA CTTCTTTCCT
             TTGGAGCTAC TGTGGAGTTA CTCTCTTTTT TGCCTTCTGA CTTCTTTCCG
             TTGGAGCTTC TGTGGAATT- CTCTCTTTTT TGCCTTCTGA TTTCTTCCCG
             TTGGAGCTTC TGTGGAGTTA CTCTCATTTT TGCCTTCTGA CTTCTTCCCG
             TTGGAGCTTC TACGGAGTTA CTCTCTTTTT TGCCTTCTGA CTTCTTTCCT

             TCAGTACGAG ATCTTCTAGA TACCGCCTCA GCTCTGTATC GGGAAGCCTT
             TCCGTCAGAG ATCTC-TAGA CACCGCCTCA GCTCTGTATC GG-AAGCCTT
             TCTATTCGAG ATCTCCTCGA CACCGCCTCT GCTCTGTATC GGGAGGCCTT
             TCGGTGCGAG ACCTCCTAGA TACCGCTGCT GCTCTGTATC GGGAAGCCTT
             TC-GTTCGGG ACCTACTCGA CACCGCTTCA GC-CT-TACC GGGATGC-TT
             TCTGTCCGGG ACCTACTCGA CACCGCTTCA GCCCTCTACC GAGATGCCTT
             TCTATTCGAG ATCTCCTCGA CACCGCCACT GCTCTGTATC GGGAGGCCTT

             AGAGTCTCCT GAGCATTGTT CACCTCACCA TACTGCACTC AGGCAAGCAA
             AGAGTCTCCT GAGCATTGCT CACCTCACCA TACTGCACTC AGGCAAGCCA
             AGAGTCTCCG GAACATTGTT CACCTCACCA TACAGCACTC AGGCAAGCTA
             AGAATCTCCT GAACATTGCT CACCTCACCA CACAGCACTC AGGCAAGCTA
             AG-ATCACC- GAACATTGCA CCC---AACA TACCGCTCTC AGGCAAGCTA
             AGAATCACCC GAACATTGCA CCCCCAACCA TACTGCTCTC AGGCAAGCTA
             AGAGTCTCCG GAACATTGTA CACCTCACCA TACGGCACTC AGGCAAGCTA

             TTCTTTGCTG GGGGGAACTA ATGACTCTAG CTACCTGGGT GGGT-GTAAT
             TTCTCTGCTG GGGGGAATTG ATGACTCTAG CTACCTGGGT GGGTAATAAT
             TTCTGTGTTG GGGTGAGTTG ATGAATCTGG CCACCTGGGT GGGAAGTAAT
             TTCTGTGCTG GGGGGAATTA ATGACTCTAG CTACCTGGGT GGGTAATAAT
             TTT-GTGCTG GGGTGAGTTA ATGACTTTGG CTTCCTGGGT GGG-AATAAT
             TTTTGTGCTG GGGTGAGTTG ATGACCTTGG CTTCCTGGGT GGGCAATAAT
             TTCTGTGTTG GGGTGAGTTA ATGAATCTAG CCACCTGGGT GGGAAGTAAT

             TTGGAAGATC CA-CATC-AG GGACCTAGTA GTCAGTTATG TCAACACTAA
             TTGGAAGATC CAGCATCCAG GGATCTAGTA GTCAATTATG TTAATACTAA
             TTGGAAGACC CAGCATCCAG GGAATTAGTA GTCAGCTATG TCAATGTTAA
             TTACAAGATC CAGCGTCCAG GGATCTAGTA GTCAATTATG TTAACACTAA
             TTGGAAGA-C CTGCAGC-TA GGGA-TAGT- GTTAACTATG TAAA--CTAA
             TTAGAGGATC CTGCAGC-AG AGATCTAGTA GTTAATTATG TCAATACTAA
             TTGGAAGATC CAGCATCCAG GGAATTAGTA GTCGGCTATG TCAACGTTAA

             TATGGGCCTA AAGTTCAGGC AACTATTGTG GTTTCACATT TCTTGTCTCA
             CATGGGT-TA AAGATCAGGC AACTATTGTG GTTTCATATA TCTTGCCTTA
             TATGGGCCTA AAAATCAGAC AACTATTGTG GTTTCACATT TCCTGTCTTA
             CATGGGCCTA AAGATCAGGC AATTATTGTG GTTTCACATT TCCTGTCTTA
             -ATGGGCCTA AAAATTAGAC AATT---GTG GTTTCACATT TCCTGCCTTA
             TATGGGCCTA AAAATTAGAC AATTATTATG GTTTCATATT TCCTGCCTTA
             TATGGGAATA AAACTAAGAC AATTATTGTG GTTTCACCTT TCCTGTCTTA

             CTTTTGGAAG AGAAACG-TA -TAGAGTATT TGGTGTCTTT CGGAGTGTGG
             CTTTTGGAAG AGAA-CTGTA CTTGAATATT TGGTCTCTTT CGGAGTGTGG
             CTTTTGGAAG AGAAACTGTT CTTGAGTATT TGGTGTCTTT TGGAGTGTGG
             CTTTTGGAAG AGAAACTGTT CTTGAATATT TGGTGTCTTT TGGAGTGTGG
             CTTTTGGAAG AGAAACAGTT CTTGAGTATT TGGTGTCTTT -GGAGTGTGG
             CATTTGGAAG AGATACTGTT CTTGAGTATT TGGTGTCTTT TGGAGTGTGG
             TGTTTGGAAG AGACACTGTT CTTGAATATT TGGTGTCTTT TGGAGTGTGG

             ATTCGCACTC CTCCAGCTTA TAGACCACCA AATGCCCCTA TCTTATCAAC
             ATTCGCACTC CTCCAGCCTA TAGACCACCA AATGCCCCTA TCTTATCAAC
             ATTCGCACTC CTCC-GCTTA CAGACCACCA AATGCCCCTA TCTTATCAAC
             ATTCGCACTC CTCCTGCCTA CAGACCACCA AATGCCCCTA TCTTATCAAC
             ATTCGCACTC CTCCTGCTTA TAGACCACCA AATGCCCCTA TC-TATCCAC
             ATTCGCACTC CACCTGCTTA TAGACCACCA AATGCCCCTA TCCTATCAAC
             ATTCGCACTC CTCCTGCATA TAGACCACCA AATGCCCCTA TCTTATCAAC

             ACTTCCGGAG ACTACTGTTG TTAGACGACG ------AGGC AGGTCCCCTA
             ACTTCCGGAA ACTACTGTTG TTAGACGACG GGACCGAGGC AGGTCCCCTA
             ACTTCCGGAA ACTACTGTTG TTAGACGACG ------AGGC AGGTCCCCTA
             ACTTCCGGAA ACTACTGTTG TTAGACGACG ------AGGC AGGTCCCCTA
             ACTTCCGGAA ACTACTGTTG TTAGACGACG ------AGGC AGGTCCCCT-
             ACTTCCGGAG ACTACTGTTG TTAGACAACG ------AGGC AGGGCCCCTA
             ACTTCCGGAA ACTACTGTTG TTAGACGACG ------AGGC AGGTCCCCTA

             GAAGAAGAAC TCCCTCGCCT CGCAGACGAA GGTCTCAATC GCCGCGTCGC
             GAAGAAGAAC TCCCTCGCCT CGCAGACGCA GATCTCAATC GCCGCGTCGC
             GAAGAAGAAC TCCCTCGCCT CGCAGACGAA GGTCTCAATC GCCGCGTCGC
             GAAGAAGAAC TCCCTCGCCT CGCAGACGAA GGTCTCAATC ACCGCGTCGC
             GAAGAAGAAC TCCCTCGCCT CGCAGACGAA GGTCTCAATC GCCGCGTCGC
             GAAGAAGAAC TCCCTCGCCT CGCAGACGAA GATCTCAATC ACCGCGTCGC
             GAAGAAGAAC TCCCTCGCCT CGCAGACGAA GGTCTCAATC GCCGCGTCGC

             AGAAGATCTC AATCTCGGGA ATCTCAATGT TAGTATTCCT TGGACTCATA
             AGAAGATCTC AATCTCGGGA ATCTCAATGT TAGTATTCCT TGGACTCATA
             AGAAGATCTC AATCTCGGGA ATCTCAATGT TAGTATCCCT TGGACTCATA
             AGAAGATCTC AATCTCGGGA ACC-CAATGT TAGTATCCCT TGGACTCATA
             AGAAGATCTC AATCTCCAGC TTCCCAATGT TAGTATTCCT TGGACTCATA
             AGAAGATCTC AATCTCCAGC TTCCCAATGT TAGTATTCCT TGGACTCATA
             AGAAGATCTC AATCTCGGGA ACCTCAATGT TAGTATTCCT TGGACACATA

             AGGTGGGAAA CTTTACGGGG CTTTATTCTT CTACTGTACC TGTCTTTAAC
             AGGTGGGAAA CTTTACTGGG CTTTATTCCT CTACAGTACC TATCTTTAAT
             AGGTGGGAAA CTTTACTGGG CTTTATTCTT CTACTGTACC TGTCTTTAAT
             AGGTGGGAAA CTTTACGGGG CTT-ATTCTT CTACAGTACC TGTCTTTAAT
             AGGTGGGAAA -TTTACGGG - CTCTAT - ct - CTACTGT - CC - -TCTTTAAT
             AGGTGGGAAA CTTTACCGGT CTTTACTCCT CTACTGTACC TGTTTTCAAT
             AGGTGGGAAA CTTTACGGGG CTTTATTCTT CTACGGTACC TTGCTTTAAT

             CCTCATTGGA AAACACCCTC TTTTCCTAAT ATACATTTAC ACCAAGACAT
             CCTGAATGGC AAACTCCTTC CTTTCCTAAG ATTCATTTAC AAGAGGACAT
             CCTGAGTGGC AAACTCCCTC CTTTCCT-AC ATTCATTTAC AGGAGGACAT
             CCTGAATGGC AAACTCCTTC TTTTCCAGAC ATTCATTTGC AGGAGGACAT
             CCT-ACTGGT TAACTCCTTC TTTTCCTGAT ATTCATTTAC ATCAAGA--T
             CCTGACTGGT TAACTCCTTC TTTTCCTGAC ATTCACTTGC ATCAAGATTT
             CCTAAATGGC AAACTCCTTC TTTTCCTAAC ATTCATTTGC AGGAGGACAT

             TATCAAAAAA TGTGAACAAT TTGTAGGCCC ACTCACAGTC AATGAGAAAA
             TATTAATAGG TGTCAACAAT TTGTGGGCCC TCTCACTGTA AATGAAAAGA
             TATTAATAGA TGTCAACAAT ATGTGGGCCC TCTTACAGTT AATGAAAAAA
             TGTTGATAGA TGTAAGCAAT TTGTGGGACC CCTTACAGTA AATGAAAACA
             TGAT-A-AAA TGTGAACAAT TTGTAGGCCC CTACA---AA AATGAATTGA
             GATACAAAAA TGTGAACAAT TTGTAGGCCC ACTCACTAAA AATGAAAGGA
             TGTTGATAGA TGTAAGCAAT TTGTGGGACC CCTTACAGTA AATGAAAACA

             GAAGACTGCA ATTGATTATG CCTGCTAGGT TTTATCCAAA TGTTACCAAA
             GAAGATTGAA ATTAATTATG CCTGCTAGAT TT-ATCCTAC CCACACTAAA
             GGAGATTAAA ATTAATTATG CCTGCTAGGT TCTATCCTAA CCTTACCAAA
             GGAGACTAAA ATTAATAATG CCTGCTAGAT TTTATCCTAA TGTTACCAAA
             GAAGATTAAA ATTGGTTATG CCA-C-AGAT TTT-TCCTAA GGTTACCAAA
             GACGATTGAA ACTAATTATG CCAGCTAGGT TTTATCCCAA AGTTACTAAA
             GGAGACTAAA ATTAATTATG CCTGCTAGGT TTTATCCCAA AGTTACTAAA

             TATTTGCCAT TGGATAAGGG TATTAAACCT TATTATCCAG AACATCTAGT
             TATTTGCCCT TAGACAAAGG AATTAAACCT TATTATCCAG ATCAGGTAGT
             TATTTGCCCT TGGA-AAAGG CATTAAACCT -ATTATCCTG AACATGCAGT
             TATTTGCCCT TAGATAAAGG GATCAAACCT TATTATCCAG AGCATGTAGT
             TATTTTCCTA TGGA--AAGG -ATTAAACCC TATTATCCTG AA-----GGT
             TACTTCCCTT TGGATAAAGG TATTAAACCT TACTATCCAG AGAATGTGGT
             TATTTGCCCT TAGACAAAGG GATCAAACCG TATTATCCAG AGTATGTGGT

             TAATCATTAC TTCCAAACCA GACATTATTT ACACACTCTA TGGAAGGCGG
             TAATCATTAC TTCCAAACCA GACATTATTT ACATACTCTT TGGAAGGCTG
             TAATCATTAC TTCAAAACTA GGCATTATTT ACATACTCTG TGGAAGGCTG
             TAATCATTAC TTCCAGAC-A GACATTATTT GCATACTCTT TGGAAGGCGG
             TAATCATTAT TTTAAA-C-A GACA-TATTT GCATACTTTA TGGAAGGCGG
             TAATCATTAC TTCAAAACTA GACACTATTT ACATACTTTG TGGAAGGCAG
             TAATCATTAC TTCCAGACGC GACATTATTT ACACACTCTT TGGAAGGCGG

             GTATATTATA TAAGAGAGAA ACAACACATA GCGCCTCATT TTGTGGGTCA
             GTATTCTATA TAAGAGGGAA ACCACACGTA GCGCATCATT TTGCGGGTCA
             GCATTCTATA TAAGAGAGAA ACTACACGCA GCGCCTCATT TTGTGGGTCA
             GTATCTTATA TAAAAGAGAG TCAACACATA GCGCCTCATT TTGCGGGTCA
             A---TTTATA TAAGAGAGAA TCCACACGTA GCGCCTCATT TTGTGGGTCA
             GAATTCTATA TAAGAGAGAA TCCACACATA GCGCCTCATT TTGTGGGTCA
             GGATCTTATA TAAAAGAGAA TCCACACGTA GCGCTTCATT TTGCGGGTCA

             CCATATTCTT GGGAACAAGA GCTACAGCAT GGG------- ----------
             CCATATTCTT GGGAACAAGA GCTACAGCAT GGGAGGTTGG TCATCAAAAC
             CCATATTCTT GGGAACAAGA GCTACAGCAT GGGAGGTTGG TCTTCCAAAC
             CCATATTCTT GGGAACAAGA TCTACAGCAT GGGAGGTTGG TCTTCCAAAC
             CCATATTC-T GGGAACAAGA GCTACAGCAT GGGAGCACCT CT-TCAACGA
             CCATATTCCT GGGAACAAGA GCTACAGCAT GGGAGCACCT CTCTCAACGG
             CCATATTCTT GGGAACAAGA TCTACAGCAT GGGAGGTGTG TCTTCCAAAC

             ---------- ------GCAG AATCTTTCCA CCAGCAATCC TCTGGGATTC
             CTCGCAAAGG CATGGGGACG AATCTTTCTG TTCCCAACCC TCTGGGATTC
             CTCGACAAGG CATGGGGACG AATCTTTCTG TTCCCAATCC TCTGGGATTC
             CTCGAAAAGG CATGGGGACA AATCTTTCTG TCCCCAATCC CCTGGGATTC
             C-AGAAGGGG CATGGGACAG AATCTCTCTG T-CCCAATCC TCTGGGATTC
             CGAGAAGGGG CATGGGACAG AATCTTTCTG TGCCCAATCC TCTGGGATTC
             CTCGAAAAGG CATGGGGACA AATCTTGCTG TCCCCAATCC CCTGGGATTC

             TTTCCCGACC ACCAGTTGGA TCCAGCCTTC AGAGCAAACA CCGCAAATCC
             TTTCCCGATC ATCAGTTGGA CCCTGCATTC GGAGCCAACT CAAACAATCC
             TTTCCCGATC ACCAGTTGGA CCCTGCGTTC GGAGCCAACT CAAACAATCC
             TTCCCCGATC ATCAGTTGGA CCCTGCATTC AAAGCCAACT CAGAAAATCC
             TT-CCAGACC ATCAGCTGGA TCC-CTATTC AG-GCAAATT CCAGCAGTCC
             TTTCCAGACC ACCAGTTGGA TCCACTATTC AGAGCAAATT CCAGCAGTCC
             TTCCCCGATC ATCAGTTGGA CCCTCTATTC AAAGCCAACT CAGAAAATCC

             AGATTGGGAC TTCAATCCCA ACAAGGACAC CTGGCCAGAC GCCAACAAGG
             AGATTGGGAC TTCAACCCCA TCAAGGACCA CTGGCCAGCA GCCAACCAGG
             AGATTGGGAC TTCAACCCCA ACAAGGATCA CTGGCCAGAG GCAAATCAGG
             AGATTGGGAC CTCAACCCAC ACAAGGACAA CTGGCCGGAC GCCCACAAGG
             CGACTGGGAC TTCAACA-AA ACAAGGACA- TTGGCCAATG GCAAACAAGG
             CGATTGGGAC TTCAACACAA ACAAGGACAA TTGGCCAATG GCAAACAAGG
             AGATTGGGAC CTCAACCCGC ACAAGGACAA CTGGCCGGAC GCCAACAAGG

             TAGGAGCTGG AGCATTCGGG CTGGGATTCA CCCCACCGCA CGGAGGCCTT
             TAGGAGTGGG AGCATTCGGG CCAGGGCTCA CCCCTCCACA CGGCGG-TTT
             TAGGAGCGGG AGCATTCGGG CCAGGGTTCA CCCCACCACA CGGCGGTCTT
             TGGGAGTGGG AGCATTCGGG CCAGGGTTCA CCCCTCCCCA TGGGGGACTG
             TAGGAGTGGG AGG-TACGGT CC-GGGTTCA CACCCCCACA CGGTGGCCTG
             TAGGAGTGGG AGGCTTCGGT CCAGGGTTCA CACCCCCACA CGGTGGCCTT
             TGGGAGTGGG AGCATTCGGG CCAGGGTTCA CCCCTCCCCA TGGGGGACTG

             TTGGGGTGGA GCCCTCAGGC TCAGGGCATA CTACAAACCT TGCCAGCAAA
             TTGGGGTGGA GCCCTCAGGC TCAGGGCATA TTGACCACAG TGTCAACAAT
             TTGGGGTGGA GCCCTCAGGC TCAGGGCATA TTGACAACAG TGCCAGCAGC
             TTGGGGTGGA GCCCTCAGGC TCAGGGCATA CTCACATCTG TGCCAGCAGC
             CTGGGGTGGA GCCC-CAGGC ACAGGGTGTT ACAAC----T TGCCAGCAGA
             CTGGGGTGGA GCCCTCAGGC ACAGGGCATT CTGACAACCT CGCCACCAGA
             TTGGGGTGGA GCCCTCAGGC TCAGGGCCTA CTCACAACTG TGCCAGCAGC

             TCCGCCTCCT GCCTCTACCA ATCGCCAGTC AGGAAGGCAG CCTACCCCGC
             TCCTCCTCCT GCCTCCACCA ATCGGCAGTC AGGAAGGCAG CCTACTCCCA
             -CCTCCTCCT GCCTCCACCA ATCGGCAGTC AGGAAGACAG CCTACTCCCA
             TCCTCCTCCT GCCTCCACCA ATCGGCAGTC AGGAAGGCAG CCTACTCCCT
             TCCGCCTCCT GCTTCCACCA ATCGGCGGTC CGGGAGA-AA CCAACCCCAG
             TCCACCTCCT GCTTCCACCA ATCGGAGGTC -AGGAAG-AA CCAACCCCAG
             TCCTCCTCCT GCCTCCACCA ATCGGCAGTC AGGAAGGCAG CCTACTCCCT

             TGTCTCCACC TTTGAGAAAC ACTCATCCTC AGGCCATGCA GTGGAACTCC
             TCTCTCCACC TCTAAGAGAC AGTCATCCTC AGGCCATGCA GTGGAATTCC
             TCTCTCCACC TCTAAGAGAC AGTCATCCTC AGGCCATGCA GTGGAACTCC
             TATCTCCACC TCTAAGGGAC ACTCATCCTC AGGCCATGCA GTGGAACTCC
             TCTCTCCACC TCTAAGAGAC AC-CATCCAC AGGCCATGCA GTGGAACTCA
             TCTCTCCACC TCTAAGGGAC ACACATCCAC AGGCCATGCA GTGGAACTCA
             TATCTCCACC TCTAAGGGAC ACTCATCCTC AGGCCATGCA GTGGAACTCC

             ACAACCTTCC ACCAAACTCT GCAAGATCCC AGAGTGAGAG GCCTGTATTT
             AC--GCTTCC ACCAAGCTCT GCA-GATCCC AGAGTCAGGG GTCTGTA-TT
             ACAACATTCC ACCAAGCTCT GCTAGA-CCC AGAGTGAGGG GCCTATA-TT
             ACCACTTTCC ACCAAACTCT TCAAGATCCC AGAGTCAGGG CTCTGTACTT
             AC-CAGTTCC ACCAGGCTCT GTTG-ATCCG AGGGTAAGGG CTCTGTATTT
             ACACAGTTCC ACCAAGCACT GTTGGATCCG AGAGTAAGGG GTCTGTATTT
             ACCACTTTCC ACCAAACTCT TCAAGATCCC AAAGTCAGGG CCCTGTACTT

             CCCTGCTGGT GGCTCCAGTT CAGGAACAGT AAACCCTGTT CCGACTACTG
             TCCTGCTGGT GGCTCCAGTT CAGGAACAGT AAACCCTGCT CCGAATATTG
             TCCTGCTGGT GGCTCCAGTT CCGGAACAGT AAACCCTGTT CCGACTACTG
             TCCTGCTGGT GGCTCCAGTT CAGGAACAGT AAGCCCTGCT CAGAATACTG
             -CCTGCTGGT GGCTCCAGTT CAGGA-CACA GAACCCTGCT CCGACTATTG
             TCCTGCTGGT GGCTCCAGTT CAGAAACACA GAACCCTGCT CCGACTATTG
             TCCTGCTGGT GGCTCCAGTT CAGGAACAGT GAGCCCTGCT CAGAATACTG

             -CTCTCCCAT ATCGTCAATC TTCTCGAGGA TTGGGGACCC TGCGCTGAAC
             CCTCTCACAT CTCGTCAATC TCCGCGAGGA CTGGGGACCC TGTGACGAAC
             CCTCACCCAT ATCGTCAATC TTCTCGAGGA CTGGGGACCC TGCACCGAAC
             TCTCTGCCAT ATCGTCAATC TTATCGAAGA CTGGGGACCC TGTGCCGAAC
             CCTCTCTCAC ATCATCAATC TTCT-GAAGA CTGGGGGCCC TGCTATGAAC
             CCTCTCTCAC ATCATCAATC TTCTCGAAGA CTGGGGACCC TGCTATGAAC
             TCTCTCCCAT ATCGTCAATC TTATCGAAGA CTGGGGACCC TGTACCGAAC

             ATGGAGAACA TCACATCAGG ATTCCTAGGA CCCCTGCTCG TGTTACAGGC
             ATGGAGAACA TCACATCAGG ATTCCTAGGA CCCCTGCTCG TGTTACAGGC
             ATGGAGA-CA CAACATCAGG ATTCCTAGGA CCCCTGCTCG TGTTACAGGC
             ATGGAGAACA TCGCATCAGG ACTCCTAGGA CCCCTGCTCG TGTTACAGGC
             ATGGACAACA TCACATCAGG ACTCCTAGGA CCCCTGCTCG TGTTACAGGC
             ATGGAGAACA TCACATCAGG ACTCCTAGGA CCCCTTCTCG TGTTACAGGC
             ATGGAGAACA TCGCATCAGG ACTCCTAGGA CCCCTGCTCG TGTTACAGGC

             GGGGTTTTTC TTGTTGACAA GAATCCTCAC AATACCGCAG AGTCTAGACT
             GGGGTTTTTC TTGTTGACAA GAATCCTCAC AATACCGCAG AGTCTAGACT
             GGGGTTTTTC TTGTTGACAA GAATCCTCAC AATACCACAG AGTCTAGACT
             GGGGTTTTTC T-GTTGACAA AAATCCTCAC AATACCACAG AGTCTAGACT
             GGTGTGTTTC TTGTTGACAA AAATCCTCAC AATACCACAG AGTCTAGACT
             GGTGTGTTTC TTGTTGACAA AAATCCTCAC AATACCACAG AGTCTAGACT
             GGGGTTTTTC TTGTTGACAA AAATCCTCAC AATACCACAG AGTCTAGACT

             CGTGGTGGAC TTCTCTCAAT TTTCTAGGGG GAACTACCGT GTGTCTTGGC
             CGTGGTGGAC TTCTCTCAAT TTTCTAGGGG GATCACCCGT GTGTCTTGGC
             CGTGGTGGAC TTCTCTCAAT TTTCTAGGGG GAGCACCCAC GTGTCCTGGC
             CGTGGTGGAC TTCTCTCAAT TTTCTAGGGG GAACACCCGT GTGTCTTGGC
             CGTGGTGGAC TTCTCTCAAT TTTCTAGGGG GA--ACCCGG GTGTCCTGGC
             CGTGGTGGAC TTCTCTCAAT TTTCTAGGGG TACCACCCGG GTGTCCTGGC
             CGTGGTGGAC TTCTCTCAAT TTTCTAGGGG AAACACCCGT GTGTCTTGGC

             CAAAATTCGC AGTCCCCAAC CTCCAATCAC TCACCAACCT CCTGTCCTCC
             CAAAATTCGC AGTCCCCAAC CTCCAATCAC TCACCAACCT CCTGTCCTCC
             CAAAATTCGC AGTCCCCAAC CTCCAATCAC TCACCAACCT CTTGTCCTCC
             CAAAATTCGC AGTCCCAAAT CTCCAGTCAC TCACCAACCT GTTGTCCTCC
             CAAAATTCGC AGTCCCCAAC CTCCAATCAC TTACCAACCT CCTGTCCTCC
             CAAAATTCGC AGTCCCCAAT CTCCAATCAC TTACCAACCT CCTGTCCTCC
             CAAAATTCGC AGTCCCAAAT CTCCAGTCAC TCACTAACCT GTTGTCCTCC

             AACTTGTCCT GGTTATCGCT GGATGTGTCT GCGGCGTTTT ATCATCTTCC
             AATTTGTCCT GGTTATCGCT GGATGTGTCT GCGGCGTTTT ATCATATTCC
             AATTTGTCCT GGCTATCGCT GGATGTGTCT GCGGCGTTTT ATCATATTCC
             AATTTGTCCT GGTTATCGCT GGATGTGTCT GCGGCGTTTT ATCATCTTCC
             AACTTGTCCT GGCTATCG-T GGATGTGTCT GCGGCGTTTT ATCATCTTCC
             AACTTGTCCT GGCTATCGTT GGATGTGTCT GCGGCGTTTT ATCATCTTCC
             AATTTGTCCT GGTTATCGCT GGATGTGTCT GCGGCGTTTT ATCATCTTCC

             TCTTCATCCT GCTGCTATGC CTCATCTTCT TGTTGGTTCT TCTGGACTAT
             TCTTCATCCT GCTGCTATGC CTCATCTTCT TATTGGTTCT TCTGGATTAT
             TCTTCATCCT GCTGCTATGC CTCATCTTCT TGTTGGTTCT TCTGGACTAC
             TCTGCATCCT GCTGCTATGC CTCATCTTCT TGTTGGTTCT TCTGGACTAT
             TCTTCATCCT GCTGCTATGC CTCATCTTCT TGTTGGTTCT TCTGGACTAT
             TCTTCATCCT GCTGCTATGC CTCATCTTCT TGTTGGTTCT TCTGGACTAT
             TCTGCATCCT GCTGCTATGC CTCATCTTCT TGTTGGTTCT TCTGGACTAT

             CAAGGTATGT TGCCCGTTTG TCCTCTAATT CCAGGATCTT CAACCACCAG
             CAAGGTATGT TGCCCGTTTG TCCTCTAATT CCAGGATCAA CAACAACCAG
             CAAGGTATGT TGCCCGTTTG TCCTCTACTT CCAGGAACAT CAACTACCAG
             CAAGGTATGT TGCCCGTTTG TCCTCTAATT CCAGGATCAT CAACCACCAG
             CAAGGTATGT TGCCCGTTTG TCCTCTACTT CCAGGATCCA CGACCACCAG
             CAAGGTATGT TGCCCGTGTG TCCTCTACTT CCAGGATCTA CAACCACCAG
             CAAGGTATGT TGCCCGTTTG TCCTCTAATT CCAGGATCAT CAACAACCAG

             CACGGGACCA TGCAGAACCT GCACGACTCC TGCTCAAGGA ACCTCTATGT
             TACGGGACCA TGCAAAACCT GCACGACTCC TGCTCAAGGC AACTCTATGT
             CACGGGACCA TGCAAGACCT GCACGATTCC TGCTCAAGGA ACCTCTATGT
             CACGGGACCA TGCAAGACCT GCACAACTCC TGCTCAAGGA ACCTCTATGT
             CACGGGACCA TGCAAAACCT GCACAACTCT TGCTCAAGGA ACCTCTATGT
             CACGGGACCC TGCAAAACCT GCACCACTCT TGCTCAAGGA ACCTCTATGT
             CACCGGACCA TGCAAGACCT GCACAACTCC TGCTCAAGGA ACCTCTATGT

             ATCCCTCCTG TTGCTGTACC AAACCTTCGG ACGGAAATTG CACCTGTATT
             TTCCCTCATG TTGCTGTACA AAACCTACGG ATGGAAATTG CACCTGTATT
             TTCCCTCTTG TTGCTGTACA AAACCTTCGG ACGGAAACTG CACTTGTATT
             TTCCCTCATG TTGCTGTACA AAACCTACGG ACGGAAACTG CACCTGTATT
             TTCCCTCCTG -TGCTGTTCC AAACCCTCGG ACGGAAACTG CAC-TGTATT
             TTCCCTCCTG CTGCTGTACC AAACCTTCGG ACGGAAATTG CACCTGTATT
             TTCCCTCATG TTGCTGTACA AAACCTACGG ACGGAAACTG CACCTGTATT

             CCCATCCCAT CATCCTGGGC TTTCGGAAAA TTCCTATGGG AGTGGGCCTC
             CCCATCCCAT CGTCT-GGGC TTTCGCAAAA TACCTATGGG AGTGGGCCTC
             CCCATCCCAT CATCCTGGGC TTTCGCAAGA TTCCTATGGG AGTGGGCCTC
             CCCATCCCAT CATCTTGGGC TTTCGCAAAA TACCTATGGG AGTGGGCCTC
             CCCATCCCAT CATCT-GGGC TTTAGGAAAA TACCTATGGG AGTGGGCCTC
             CCCATCCCAT CATCTTGGGC TTTCGGAAAA TACCTATGGG AGTGGGCCTC
             CCCATCCCAT CATCTTGGGC TTTCGCAAAA TACCTATGGG AGTGGGCCTC

             AGCCCGTTTC TCCTGGCTCA GTTTACTAGT GCCATTTGTT CAGTGGTTCG
             AGTCCGTTTC TCTTGGCTCA GTTTACTAGT GCCATTTGTT CAGTGGTTCG
             AGTCCGTTTC TCCTGGCTCA GTTTACTAGT GCCATTTGTT CAGTGGTTCG
             AGTCCGTTTC TCTTGGCTCA GTTTACTAGT GCCATTTGTT CAGTGGTTCG
             AGCCCGTTTC TCCTGGCTCA GTTTACTAGT GCAATTTGTT CAGTGGTGCG
             AGCCCGTTTC TCTTGGCTCA GTTTACTAGT GCAATTTGTT CAGTGGTGCG
             AGTCCGTTTC TCTTGGCTCA GTTTACTAGT GCCATTTGTT CAGTGGTTCG

             TAGGGCTTTC CCCCACTGTT TGGCTTTCAG TTATATGGAT GATGTGGTAT
             TAGGGCTTTC CCCCACTGTT TGGCTTTCAG CTATATGGAT GATGTGGTAT
             TAGGGCTTTC CCCCACTGTT TGGCTTTCAG TTATATGGAT GATGTGGTAT
             TAGGGCTTTC CCCCACTGTC TGGCTTTCAG TTATATGGAT GATGTGGTAT
             TAGGGCTTTC CCCCACTGT- TGGCTTTTAG TTATATGGAT GATCTGGTAT
             TAGGGCTTTC CCCCACTGTC TGGCTTTTAG TTATATGGAT GATTTGGTAT
             TAGGGCTTTC CCCCACTGTC TGGCTTTCAG TTACATGGAT GATGTGGTTT

             TGGGGGCCAA GTCTGTACAG CATCTTGAGT CCCTTTTTAC CGCTGTTACC
             TGGGGGCCAA GTCTGTACAG CATCGTGAGC CC-TTTATAC CGCTGTTACC
             TGGGGGCCAA GTCTGTACAA CATCTTGAGT CCCTTTTTAC CTCTATTACC
             TGGGGGCCAA GTCTGTACAA CATCTTGAGT CCCTTTATGC CGCTGTTACC
             TGGGGGCCAA ATCTGTGCAG CATCTTGAGT CCCTTTATAC CGCTGTTACC
             TGGGGGCCAA ATCTGTGCAG CATCTTGAGT CCCTTTATAC CGCTGTTACC
             TGGGGGCCAA GTCTGTACAA CATCTTGAGT CCCTTTATGC CGCTGTTACC

             AATTTTCTTT TGTCTTTGGG TATACATTTA AACCCTAACA AAACAAAAAG
             AATTTTCTTT TGTCTCTGGG TATACATTTA AACCCTAACA AAACAAAAAG
             AATTTTCTTT TGTCTTTGGG TATACATTTG AACCCTAATA AAACCAAACG
             AATTTTCTTT TGTCTTTGGG TATACATTTA AACCCTCACA AAACAAAAAG
             AATTTTCTGT TATCTGTGGG TATCCATTTA AATACCTCTA AAACAAAAAG
             AATTTTTTGT TATCTGTGGG CATCCATTTG AACACAGCTA AAACAAAATG
             AATTTTCTTT TGTCTTTGGG TATACATTTA AACCCTCACA AAACAAAACG

             ATGGGGTTAC TCTTTACATT TCATGGGCTA TGTCATTGGA TGTTATGGGT
             ATGGGGTTAT TCCCTAAACT TCATGGGTTA C-TAATTGGA AGTTGGGGAA
             TTGGGGCTAC TCCCTTAACT TCATGGGATA TGTAATTGGA AGTTGGGGTA
             ATGGGGATAT TCCCTTAACT TCATGGGATA TGTAATTGGG AGTTGGGGCA
             ATGGGGTTAA CTA----ATT TCATGGGTTA TGTTATTGG- AGTTGGGG--
             GTGGGGTTAT TCCTTACACT TTATGGGTTA TATCATTGGG AGTTGGGGGA
             ATGGGGATAT TCCCTTAACT TCATGGGATA TGTAATTGGG AGTTGGGGCA

             CATTGCCACA AGATCACATC ATACAGAAAA TCAAAGAATG TTTTAGAAAA
             CATTGCCACA GGATCATATT GTACAAAAGA TCAAACACTG TTTTAGAAAA
             CTTTA-CCCA G-AACATATT GTAC-TAAAA TCAAGCAATG TTTTCG-AAA
             CATTGCCACA GGAACATATT GTACAAAAAA TCAAACTATG TTTTAGGAAA
             ACTTACC-CA AGATCA-ATT GTAC-AAAAA TCAAAGA-TG TTTTCG-AAA
             CATTGCCTCA GGAACATATT GTGCATAAAA TCAAAGATTG CTTTCGCAAA
             CATTGCCACA GAACCATATT GTACAAAAAA TCAAAATGTG TTTTAGGAAA

             CTTCCTGTTA ACAGGCCTAT TGATTGGAAA GTCTGTCAAC GTATTGTGGG
             CTTCCTGTTA ACAG-CCTAT TGATTGGAAA GTATGTCAAA GAATTGTGGG
             CTGCCTGTAA ATAGACCTAT TGATTGGAAA GTATGTCAA- GAATTGTGGG
             CTTCCTGTAA ACAGGCCTAT TGATTGGAAA GT-TGTCAAC GAATTGTGGG
             CTTCCTGTAA ATCG-CC-AT TGATTGGAAA GTTTGTCAAC GCATTGTGGG
             CTTCCCGTGA ATAGACCCAT TGATTGGAAG GTTTGTCAAC GAATTGTGGG
             CTTCCTGTCA ACAGGCCTAT TGATTGGAAA GTATGTCAAC GAATTGTGGG

             TCTTTTGGGT TTTGCTGCCC CTTTTACACA ATGTGGTTAT CCTGCTTTAA
             TCTTTTGGGC TTTGCTGCTC CATTTACACA ATGTGGATAT CCTGCCTTAA
             TCTTTTGGGC TTTGCTGCCC CTTTTACACA ATGTGGCTAT CCTGCCTTA-
             TCTTTTGGGG TTTGCTGCCC CTTTTACGCA ATGTGGATAT CCTGCTTTAA
             TCTTTTGGGC TTTGC-GCCC C-TT-AC-CA ATGTGGTTAT CCTGCTCTCA
             TCTTTTGGGC TTTGCAGCCC CTTTTACTCA ATGTGGTTAT CCTGCTCTCA
             TCTTTTGGGG TTTGCCGCCC CTTTCACGCA ATGTGGATAT CCTGCTTTAA

             TGCCCTTGTA TGCATGTATT CAATCTAAGC AGGCTTTCAC TTTCTCGCCA
             TGCCTTTGTA TGCATGTATA CAAGCTAAAC AGGCTTTCAC TTTCTCGCCA
             TGCCTTTATA TGCATGTATA CAATCTAAGC AGGCTTTCAC TTTCTCGCCA
             TGCCTTTATA TGCATGTATA CAAGCAAAAC AGGCTTTTAC TTTCTCGCCA
             TGCCTCTGTA T-CCTGTATA CTGCTAAA-C AGGCTTTTGT CTTTTCGCCA
             TGCCCTTGTA TGCCTGTATT ACCGCTAA-C AGGCTTTTGT TTTCTCGCCA
             TGCCTTTATA TGCATGTATA CAAGCAAAAC AGGCTTTTAC TTTCTCGCCA

             ACTTACAAGG CCTTTCTGTG TAAACAATAC CTGAACCTTT ACCCCGTTGC
             ACTTACAAGG CCTTTCTAAG TAAACAGTAC ATGAACCTTT ACCCCGTTGC
             ACTTACAAGG CCTTTCTGTG TAAACAATAT CTGAACCTTT ACCCCGTTGC
             ACTTACAAGG CCTTTCTAAG TAAACAGTAT CTGACCCTTT ACCCCGTTGC
             ACTTACAAGG CCTTTCT-TG TAA-CAATAC ATGAACCTTT ACCCCGTTGC
             ACTTACAAGG CCTTTCTCTG TAAACAATAC ATGAACCTTT ACCCCGTTGC
             ACTTACAAGG CCTTTCTAAG TCAACAGTAT CTGAACCTTT ACCCCGTTGC

             CCGGCAACGG CCAGGTCTGT GCCAAGTGTT TGCTGACGCA ACCCCCACTG
             TCGGCAACGG CCTGGTCTGT GCCAAGTGTT TGCTGACGCA ACCCCCACTG
             CCGGCAACGG TCAGGTCTCT GCCAAGTGTT TGCTGACGCA ACCCCCACTG
             TCGGCAACGG CCTGGTCTGT GCCAAGTGTT TGCTGACGCA ACCCCCACTG
             TCGGCAACGG CCAGGCCTGT GCCAAGTGTT TGCTGACGCA ACCCCCACTG
             TCGGCAACGG CCAGGCCTTT GCCAAGTGTT TGCTGACGCA ACCCCCACTG
             TCGGCAACGG CCTGGTCTGT GCCAAGTGTT TGCTGACGCA ACCCCCACTG

             GCTGGGGCTT GGTCATGGGC CATCAGCGCA TGCGTGGAAC CTTTCTGGCT
             GCTGGGGCTT GGCCATAGGC CATCAGCGCA TGCGTGGAAC CTTTGTGGCT
             GATGGGGCTT GGC-AT-GGC CATCGGCGCA TGCGTGGAAC CTTTGTGGCT
             GTTGGGGCTT GGCCATAGGC CATCAGCGCA TGCGTGGAAC CTTTGTGTCT
             GTTGGGGCTT GGCCATTGGC CATCAGCGCA TGCGTGGAAC CTTTGTGGCT
             GCTGGGGCTT GGCGATTGGC CATCAGCGCA TGCGCGGAAC CTTTGTGGCT
             GTTGGGGCTT GGCCATAGGC CATCGGCGCA TGCGTGGGAC CTTTGTGTCT

             CCTCTGCCGA TCCATACTGC GGAACTCCTA GCCGCTTGTT TTGCTCGCAG
             CCTCTGCCGA TCCATACTGC GGAACTCCTA GCCGCTTGTT TTGCTCGCAG
             CCTCTGCCGA TCCATACTGC GGAACTCCTA GCAGCTTGTT TTGCTCGCAG
             CCTCTGCCGA TCCATACTGC GGAACTCCTA GCCGCTTGTT TTGCTCGCAG
             CCTCTGCCGA TCCATACTGC GGAACTCCTT GCAGCTTGTT TCGCTCGCAG
             CCTCTGCCGA TCCATACTGC GGAACTCCTA GCAGCTTGTT TCGCTCGCAG
             CCTCTGCCGA TCCATACTGC GGAACTCCTA GCCGCTTGTT TTGCTCGCAG

             CAGGTCTGGA GCAAACATTC TCGGGACGGA TAACTCTGTT GTTCTCTCCC
             CCGGTCTGGA GC-AAACTCA TCGGAACTGA CAATTCTGTC GTCCTCTCGC
             CCGGTCTGGA GCGAAACTTA TCGGA-C-GA CAACTCTGTT GTCCTCTCTC
             CAGGTCTGGA GCGAAACTCA TCGGGACTGA CAATTCTGTC GTGCTCTCCC
             CCGGTCTGGA GCGAA--TCA TCGGCACAGA CAACTCTGTT GTCCTCTCTA
             CCGGTCTGGA GCGGACATTA TCGGCACTGA CAACTCTGTT GTCCTTTCTC
             CAGGTCTGGG GCAAAACTCA TCGGGACTGA CAATTCTGTC GTGCTCTCCC

             GCAAATATAC ATCGTTTCCA TGGCTGCTAG GCTGTGCTGC CAACTGGATC
             GGAAATATAC ATCGTTTCCA TGGCTGCTAG GCTGTACTGC CAACTGGATC
             GGAAATACAC CTCC-TTCCA TGGCTGCTGG G-TGTGCTGC CAACTGGATC
             GCAAGTATAC ATC-TTTCCA TGGCTGCTAG GCTGTGCTGC CAACTGGATC
             GGAAGTACAC CTCCTTCC-A TGGCTGCTCG GTTGTGCTGC CAACTGGATC
             GGAAGTACAC CTCCTTCCCA TGGCTGCTAG GCTGTGCTGC CAACTGGATC
             GCAAGTATAC ATCATTTCCA TGGCTGCTAG GCTGTGCTGC CAACTGGATC

             CTGCGCGGGA CGTCCTTTGT TTACGTCCCG TCGGCGCTGA ATCCCGCGGA
             CTTCGCGGGA CGTCCTTTGT TTACGTCCCG TCGGCGCTGA ATCCCGCGGA
             CTGCGCGGGA CGTCCTTTGT CTACGTCCCG TCGGCGCTGA ATCCCGCGGA
             CTGCGCGGGA CGTCCTTTGT TTACGTCCCG TCGGCGCTGA ATCCCGCGGA
             CTGCGCGGGA CGTCCTTTGT TTACGTCCCG TCGGCGCTGA ATCCCGCGGA
             CTGCGCGGGA CGTCCTTTGT CTACGTCCCG TCGGCGCTGA ATCCTGCGGA
             CTGCGCGGGA CGTCCTTTGT TTACGTCCCG TCGGCGCTGA ATCCCGCGGA

             CGACCCTTCT CGGGGCCGCT TGGGACTCTC TCGTCCCCTT CTCCGTCTGC
             CGACCCCTCT CGGGGCCGCT TGGGACTCTA TCGTCCCCTT CTCCGTCTGC
             CGACCCGTCT CGGGGCCGTT TGGG-CTCTA CCGTCCCCTT CTTC-TCTGC
             CGACCCCTCC CGGGGCCGCT TGGGGCTCTA CCGCCCGCTT CTCCGTCTGC
             CGACCC-TCC -GGGGTCGCT TGGGGCTGTA CCGCCCCCTT CTC-GTCTGC
             CGACCCCTCT CGTGGTCGCT TGGGGCTCTG CCGCCCTCTT CTCCGCCTGC
             CGACCCCTCC CGGGGCCGCT TGGGGCTCTA CCGCCCGCTT CTCCGCCTGT

             CGTTTCGACC GACCACGGGG CGCACCTCTC TTTACGCGGA CTCCCCGTCT
             CGTTCCAGCC GACCACGGGG CGCACCTCTC TTTACGCGGT CTCCCCGTCT
             CGTTCCGGCC GACCACGGGG CGCACCTCTC TTTACGCGGT CTCCCCGTCT
             CGTACCGACC GACCACGGGG CGCACCTCTC TTTACGCGGA CTCCCCGTCT
             CGTTCCAGCC GACGACGGGT CGCACCTCTC TTTACGCGG- CTCCCCGTCT
             CGTTCCGGCC GACGACGGGT CGCACCTCTC TTTACGCGGA CTCCCCGCCT
             TGTACCAACC GACCACGGGG CGCACCTCTC TTTACGCGGA CTCCCCGTCT

             GTGCCTTCTC ATCTGCCGGA CCGTGTGCAC TTCGCTTCAC CTCTGCACGT
             GTGCCTTCTC ATCTGCCGGT CCGTGTGCAC TTCGCTTCAC CTCTGCACGT
             GTGCCTTCTC ATCTGCCGGA CCGTGTGCAC TTCGCTTCAC CTCTGCACGT
             GTGCCTTCTC -TCTGCCGGA CCGTGTGCAC TTCGCTTCAC CTCTGCACGT
             GTTCCTTCTC ATCTGCCGGA CCGTGTGCAC TTCGCTTCAC CTCTGCACGT
             GTGCCTTCTC ATCTGCCGGC CCGTGTGCAC TTCGCTTCAC CTCTGCACGT
             GTGCCTTCTC ATCTGCCGGA CCGTGTGCAC TTCGCTTCAC CTCTGCACGT

             CGCATGGAGA CCACCGTGAA CGCCCACCAA TTCTTGCCCA AGGTCTTACA
             TGCATGGAGA CCACCGTGAA CGCCCATCAG ATCCTGCCCA AGGTCTTACA
             CGCATGGAGA CCACCGTGAA CGCCCACCAG GTCTTGCCCA AGGTCTTACA
             CGCATGGAGA CCACCGTGAA CGCCCACCGG AACCTGCCCA AGGTCTTGCA
             CGCATGGAGA CCACCGTGAA CGCCCC---G AG-TTGCCAA CAGTCTTACA
             CGCATGGAGA CCACCGTGAA CGCCCCTCGG AGCTTGCCAA CAACCTTACA
             CGCATGGAGA CCACCGTGAA CGCCCACGGG AACCTGCCCA AGGTCTTGCA

             TAAGAGGACT CTTGGACTCT CTGTAATGTC AACGACCGAC CTTGAGGCAT
             TAAGAGGACT CTTGGACTCC CAGCAATGTC AACGACCGAC CTTGAGGCCT
             TAAGAGGACT CTTGGACTCT CAGCAATGTC AACGACCGAC CTTGAGGCAT
             TAAGAGGACT CTTGGACTTT CAGCAATGTC AACGACCGAC CTTGAGGCAT
             TAAG-GGACT CTTGGACTTT CAGGACGGTC AATGACCTGG ATCGAAGA-T
             TAAGAGGACT CTTGGACTTT CGCCCCGGTC AACGACCTGG ATTGAGGAAT
             TAAGAGGACT CTTGGACTCT CAGCAATGTC AACGACCGAC CTTGAGGCAT

             ACTTCAAAGA CTGTTTGTTT AAAGACTGGG AGGAGTTGGG GGAGGAGATT
             ACTTCAAAGA CTGTGTGTTT AAGGACTGGG AGGAGCTGGG GGAGGAGATT
             ACTTCAAAGA CTGTTTGTTT AA-GACTGGG AGGAGTTGGG GGAGGAGATT
             ACTTCAAAGA CTGTGTGTTT ACTGAGTGGG AGGAGCTGGG GGAGGAGATT
             ACATCAAAGA CTGTGTATTT AAGGACTGGG AGGAGCTGGG GGAGGAGATC
             ACATCAAAGA CTGTGTATTT AAGGACTGGG AGGAGTCGGG GGAGGAGTTG
             ACTTCAAAGA CTGTGTGTTT ACTGAGTGGG AGGAGTTGGG GGAGGAGGTT

             AGATTAAAGG TCTTTGTA-T AGGAGGCTGT AGGCATAAAT TGGTCTGCGC
             AGGTTAAG-- TCTTTGTATT AGGAGGCTGT AGGCATAAAT TGGTCTGCGC
             AGGTTAATGA TCTTTGTACT AGGAGGCTGT AGGCATAAAT TGGTCTGTTC
             AGGTTAAAGG TCTTTGTACT AGGAGGCTGT AGGCATAAAT TGGTCTGTTC
             AGGTTAAAGG TCTTTGTA-T AGGAGGCTGT AGGCATAAAT TGGTCTGTTC
             AGGTTAAAGG TCTTTGTATT AGGAGGCTGT AGGCATAAAT TGGTCTGTTC
             AGGTTAAAGG TCTTTGTACT AGGAGGCTGT AGGCATAAAT TGGTGTGTTC

             ACCAGCACCA TGCAACTTTT T
             ACCAGCACCA TGCAACTTTT T
             ACCAGCACCA TGCAACTTTT T
             ACCAGCACCA TGCAACTTTT T
             ACCAGCACCA TGCAACTTTT T
             ACCAGCACCA TGCAACTTTT T
             ACCAGCACCA TGCAACTTTT T











examples/cmdfile.AB073841
model: dcp_recomb       # Select the model to run:
    # dcp_recomb: dual change point (DMCP) model for recombination
    # recomb: single change point (SMCP) model for recombination
length: 10000000        # RUN: Length of MCMC run
burnin: 0           # RUN: Number of initial samples to discard
subsample: 1000         # RUN: Report state to screen and posterior file every nth sample (specify n here)
change_points: true     # RUN (DMCP): Add/Delete what kind of change points
    # true: normal behavior
    # parameter: add/delete parameter change points only
    # topology: add/delete topology change points only
par_lambda: 6           # PRIOR (DMCP): expected number of parameter change points (truncated Poisson prior distribution)
top_lambda: 0.693       # PRIOR (DMCP): expected number of topology change points (truncated Poisson prior distribution)
#lambda: 3          # PRIOR (SMCP): expected number of change points (Poisson assumed)
titv_hyper_mean: 2      # PRIOR (DMCP): mean of log kappa
titv_hyper_variance: 1      # PRIOR (DMCP): variance of log kappa
mu_hyper_mean: -2       # PRIOR (DMCP): mean of log mu
mu_hyper_variance: 2        # PRIOR (DMCP): variance of log mu
#gmodel: 0.5            # PRIOR (MREM): Fang's prior on \tau
window_length: 10       # PROPOSAL: Size of window to randomly move \xi or \rho within during fixed dimension update
sigma_alpha: 1.75       # PROPOSAL: Random perturbation size for fixed dimension update of kappa (the transition/transversion ratio)
sigma_mu: 1.75          # PROPOSAL: Random perturbation size for fixed dimension update for mu (the average branch length)
C: 0.10             # PROPOSAL: Multiplier of move probabilities (add one, delete one, etc.)
    # Reset this option if you receive the error message: "You probably don't want the sum of all move probabilities to exceed 0.9"
start_tree: ((0,(2,(1,3))),(5,4))   # PROPOSAL: Fixed reference tree on which recombinant queries move
alawadhi: false         # PROPOSAL (Al-Awadhi): use alawadhi updates
    # false: turn off Al-Awadhi updates
    # true or topology|parameter: on all dimension change updates
    # topology: when changing the number of topology change points only
    # parameter: when changing the number of parameter change pointsn only
alawadhi_k: 10          # PROPOSAL (Al-Awadhi): Number of fixed dimension sampler iterations during alawadhi update x' -> x*
alawadhi_factor: 0.7        # PROPOSAL (Al-Awadhi): Multiplicative factor of posterior distribution during alawadhi update x' -> x*
alawadhi_debug: false       # DEBUG (Al-Awadhi): Turn on debugging output (limited) of alawadhi updates
debug: false            # DEBUG: Turn on verbose debugging output; the higher the number, the more verbose
compute_likelihood: true    # DEBUG: Set to false to test recovery of the prior
report_iact: true       # DEBUG: Turn on output of expected sample size
exit_condition: false       # DEBUG: set true/false to trigger exit as defined in dcpsampler.c (requires compile to set exit condition)
#random: ../rn          # DEBUG: Specify a file with random standard uniform numbers rather than use internal rng (for programmers)
simulate_data: false        # SIMULATION: set to false if you don't want to simulate data, otherwise it is the number of simulated segments
sim_length: 500 500 200 500 # SIMULATION: length of each simulated alignment segment
sim_mu: 0.03 0.05 0.02 0.06 # SIMULATION: mu for each segment
sim_kappa: 2.0 2.0 2.0 2.0  # SIMULATION: kappa for each segment; next line is tree for each segment
sim_tree: ((3,(2,1)),((4,5),0)) ((3,(2,1)),((4,5),0)) ((3,(2,1)),((4,5),0)) ((3,(2,1)),((4,5),0))
sim_pi: 0.22 0.23 0.28 0.27 # SIMULATION: stationary distribution for simulating nucleotide sequence (frequency of A,C,G,T)
# Any full-line comment that starts with # in first position is OK
examples/cmdfile.simulation
model: dcp_recomb       # Select the model to run:
    # dcp_recomb: dual change point (DMCP) model for recombination
    # recomb: single change point (SMCP) model for recombination
length: 10000000        # RUN: Length of MCMC run
burnin: 0           # RUN: Number of initial samples to discard
subsample: 1000         # RUN: Report state to screen and posterior file every nth sample (specify n here)
change_points: parameter    # RUN (DMCP): Add/Delete what kind of change points
    # true: normal behavior
    # parameter: add/delete parameter change points only
    # topology: add/delete topology change points only
par_lambda: 3           # PRIOR (DMCP): expected number of parameter change points (truncated Poisson prior distribution)
top_lambda: 3           # PRIOR (DMCP): expected number of topology change points (truncated Poisson prior distribution)
#lambda: 3          # PRIOR (SMCP): expected number of change points (Poisson assumed)
titv_hyper_mean: 2      # PRIOR (DMCP): mean of log kappa
titv_hyper_variance: 1      # PRIOR (DMCP): variance of log kappa
mu_hyper_mean: -2       # PRIOR (DMCP): mean of log mu
mu_hyper_variance: 2        # PRIOR (DMCP): variance of log mu
#gmodel: 0.5            # PRIOR (MREM): Fang's prior on \tau
window_length: 10       # PROPOSAL: Size of window to randomly move \xi or \rho within during fixed dimension update
sigma_alpha: 1.75       # PROPOSAL: Random perturbation size for fixed dimension update of kappa (the transition/transversion ratio)
sigma_mu: 1.75          # PROPOSAL: Random perturbation size for fixed dimension update for mu (the average branch length)
C: 0.10             # PROPOSAL: Multiplier of move probabilities (add one, delete one, etc.)
    # Reset this option if you receive the error message: "You probably don't want the sum of all move probabilities to exceed 0.9"
start_tree: ((0,(2,(1,3))),4)   # PROPOSAL: Fixed reference tree on which recombinant queries move
alawadhi: false         # PROPOSAL (Al-Awadhi): use alawadhi updates
    # false: turn off Al-Awadhi updates
    # true or topology|parameter: on all dimension change updates
    # topology: when changing the number of topology change points only
    # parameter: when changing the number of parameter change pointsn only
alawadhi_k: 10          # PROPOSAL (Al-Awadhi): Number of fixed dimension sampler iterations during alawadhi update x' -> x*
alawadhi_factor: 0.7        # PROPOSAL (Al-Awadhi): Multiplicative factor of posterior distribution during alawadhi update x' -> x*
alawadhi_debug: false       # DEBUG (Al-Awadhi): Turn on debugging output (limited) of alawadhi updates
debug: false            # DEBUG: Turn on verbose debugging output; the higher the number, the more verbose
compute_likelihood: true    # DEBUG: Set to false to test recovery of the prior
report_iact: true       # DEBUG: Turn on output of expected sample size
exit_condition: false       # DEBUG: set true/false to trigger exit as defined in dcpsampler.c (requires compile to set exit condition)
#random: ../rn          # DEBUG: Specify a file with random standard uniform numbers rather than use internal rng (for programmers)
simulate_data: 4        # SIMULATION: set to false if you don't want to simulate data, otherwise it is the number of simulated segments
sim_length: 500 500 200 500 # SIMULATION: length of each simulated alignment segment
sim_mu: 0.03 0.05 0.02 0.06 # SIMULATION: mu for each segment
sim_kappa: 2.0 2.0 2.0 2.0  # SIMULATION: kappa for each segment; next line is tree for each segment
sim_tree: ((3,(2,1)),((4,5),0)) ((3,(2,1)),((4,5),0)) ((3,(2,1)),((4,5),0)) ((3,(2,1)),((4,5),0))
sim_pi: 0.22 0.23 0.28 0.27 # SIMULATION: stationary distribution for simulating nucleotide sequence (frequency of A,C,G,T)
# Any full-line comment that starts with # in first position is OK


