// MathFuncsDll.h

#ifdef MATHFUNCSDLL_EXPORTS
#define MATHFUNCSDLL_API __declspec(dllexport) 
#else
#define MATHFUNCSDLL_API __declspec(dllimport) 
#endif


#include <float.h>
#include <math.h>
#include <string>
#include <comutil.h>
//#include "dna.h"

#define NJ_NAME_STATE  100
#define NJ_FLOAT_STATE 101
#define NJ_WS_STATE    102
#define NJ_EOF_STATE   103

#define NJ_PARSE_SYMMETRIC 100
#define NJ_PARSE_LOWER     101
#define NJ_PARSE_UPPER     102
#define NJ_PARSE_UNKNOWN   103
#define NJ_NCELLS(a)       ( ((a)*(a+1))/2 )

#define NJ_INTERNAL_NODE -1
#define NJ_LAST 101

#define NJ_INPUT_MODE_UNKNOWN             0
#define NJ_INPUT_MODE_DISTANCE            100
#define NJ_INPUT_MODE_UNALIGNED_SEQUENCES 101
#define NJ_INPUT_MODE_ALIGNED_SEQUENCES   102

#define NJ_MODEL_NONE    100
#define NJ_MODEL_JUKES   101
#define NJ_MODEL_KIMURA  102

/* some data structures */
typedef struct _NJ_DIST_TOKEN_STRUCT {

	char *buf;
	long int bufsize;
	int type;

} NJ_DIST_TOKEN;


//struct nodex_struct {
//
//	int neighbourindex;
//	struct nodex_struct *next;
//};
//
//typedef struct nodex_struct nodex;

typedef struct XOVERDEFINE {
	unsigned char OutsideFlag;
	unsigned char  MissIdentifyFlag;
	unsigned char  ProgramFlag;
	unsigned char  SBPFlag;
	unsigned char  Accept;
	short int MajorP;
	short int MinorP;
	short int Daughter;
	int Beginning;
	int Ending;
	int LHolder;
	int Eventnumber;
	float PermPVal;
	int BeginP;
	int EndP;
	double Probability;
	double DHolder;

}XOVERDEFINE;

typedef struct POINTAPI {
	int x;
	int y; 
} POINTAPI;

typedef struct node {
	struct node *next, *back;
	bool tip;
	long number;
	short int nayme;
	float v;
	long xcoord, ycoord, ymin, ymax;
} node;

typedef struct treeX {
	node **nodep;
	node *start;
} treeX;

typedef struct DMAT {

	long int ntaxa;   /* the original size of the distance matrix */
	long int size;    /* the current/effective size of the distance matrix */

	char **taxaname;  /* a pointer to an array of taxa name strings */

	float *val;       /* the distances */
	float *valhandle; /* to track the orig. pointer to free memory */

	float *r, *r2;    /* r and r2 vectors (used to compute transformed dists) */
	float *rhandle, *r2handle;  /* track orig. pointers to free memory */

} DMAT;

typedef struct NJ_TREE {

	struct NJ_TREE *left;  /* left child  */
	struct NJ_TREE *right; /* right child */

	float dist;  /* branch length.  i.e. dist from node to parent */

	long int taxa_index; /* for terminal nodes, track the taxon index */

} NJ_TREE;

typedef struct _STRUCT_NJ_VERTEX {

	NJ_TREE **nodes;
	NJ_TREE **nodes_handle;  /* original memory handle for freeing */
	long int nactive;  /* number of active nodes in the list */
	long int size;     /* the total size of the vertex */

} NJ_VERTEX;


/*
* NJ_MAP() -
*
* Thus function maps i, j coordinates to the correct offset into
* the distance matrix
*
*/
static inline
long int
NJ_MAP(long int i,
	long int j,
	long int ntaxa) {

	return((i*(2 * ntaxa - i - 1)) / 2 + j);
}


static inline
int
NJ_FLT_EQ(float x,
	float y) {

	if (fabs(x - y)<FLT_EPSILON) {
		return(1);
	}
	else {
		return(0);
	}
}



static inline
int
NJ_FLT_LT(float x,
	float y) {

	if (NJ_FLT_EQ(x, y)) {
		return(0);
	}
	else {
		if (x < y) {
			return(1);
		}
		else {
			return(0);
		}
	}
}


static inline
int
NJ_FLT_GT(float x,
	float y) {

	if (NJ_FLT_EQ(x, y)) {
		return(0);
	}
	else {
		if (x > y) {
			return(1);
		}
		else {
			return(0);
		}
	}
}



namespace MathFuncs
{
	// This class is exported from the MathFuncsDll.dll
	class MyMathFuncs
	{
	public:
		// Returns a + b
		static MATHFUNCSDLL_API double _stdcall Add(double a, double b);

		// Returns a - b
		static MATHFUNCSDLL_API double _stdcall Subtract(double a, double b);

		// Returns a * b
		static MATHFUNCSDLL_API double _stdcall Multiply(double a, double b);

		// Returns a / b
		// Throws const std::invalid_argument& if b is 0
		static MATHFUNCSDLL_API double _stdcall Divide(double a, double b);





		static MATHFUNCSDLL_API double _stdcall SuperDistP(int X, int Nextno, int UB14, int UB04, int UB13, int UB03, int UB12, int UB02, int UB11, double *avdst, float *pd, float *pv, float *dist, short int *redodist, int *SeqCatCount, short int *ISeq14, short int *ISeq04, short int *ISeq13, short int *ISeq03, short int *ISeq12, short int *ISeq02, short int *ISeq11, char *CompressValid14, char *CompressDiffs14, char *CompressValid13, char *CompressDiffs13, char *CompressValid12, char *CompressDiffs12, char *CompressValid11, char *CompressDiffs11, char *CompressDiffs04, char *CompressDiffs03, char *CompressDiffs02);
		static MATHFUNCSDLL_API double _stdcall SuperDistP2(int XX, int Nextno, int UB14, int UB04, int UB13, int UB03, int UB12, int UB02, int UB11, double *avdst, float *pd, float *pv, float *dist, short int *redodist, int *SeqCatCount, short int *ISeq14, short int *ISeq04, short int *ISeq13, short int *ISeq03, short int *ISeq12, short int *ISeq02, short int *ISeq11, char *CompressValid14, char *CompressDiffs14, char *CompressValid13, char *CompressDiffs13, char *CompressValid12, char *CompressDiffs12, char *CompressValid11, char *CompressDiffs11, char *CompressDiffs04, char *CompressDiffs03, char *CompressDiffs02);

		static MATHFUNCSDLL_API int _stdcall CMaxD2P(int incnum, int Seq1, int Seq2, int Seq3, int SBP, int EBP, int Nextno, int SLen, short int *SeqNum, short int *SeqnumX, int *IdenticalR, int *IdenticalF, unsigned char *NucMat, int *IncSeq3, unsigned char *IncSeq2, unsigned char *IncSeq, float *E, float *d, float *VScoreMat, float *AvDist, int *TotCount);
		static MATHFUNCSDLL_API int _stdcall CMaxD2P3(int incnum, int Seq1, int Seq2, int Seq3, int SBP, int EBP, int Nextno, int SLen, short int *SeqNum, short int *SeqnumX, int *IdenticalR, int *IdenticalF, unsigned char *NucMat, int *IncSeq3, unsigned char *IncSeq2, unsigned char *IncSeq, float *E, float *d, float *VScoreMat, float *AvDist, int *TotCount);
		static MATHFUNCSDLL_API int _stdcall GetFragsP(short int CircularFlag, int LenXoverSeq, int LSeq, int maxcount, char *SubSeq, int *FragSt, int *FragEn, int *FragScore, int *FragCount);
		static MATHFUNCSDLL_API int _stdcall GetFragsP3(short int CircularFlag, int LenXoverSeq, int LSeq, int maxcount, char *SubSeq, int *FragSt, int *FragEn, int *FragScore, int *FragCount);

		static MATHFUNCSDLL_API int _stdcall GetMaxFragScoreP(int LenXoverSeq, int lseq, short int CircularFlag, short int GCMissmatchPen, double *MissPen, int *MaxScorePos, int *FragMaxScore, int *FragScore, int *FragCount, int *hiscore);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqGCAP(char gcindelflag, int LSeq, int seq1, int seq2, int seq3, short int *SeqNum, char *SubSeq, int *XPosDiff, int *XDiffPos, int *NDiff);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqGCAP5(char gcindelflag, int LSeq, int seq1, int seq2, int seq3, short int *SeqNum, char *SubSeq, int *NDiff);
		static MATHFUNCSDLL_API double _stdcall GCGetHiPValP(int lseq, int LenXoverSeq, int *FragCount, double *PVals, int *MaxY, int *MaxX, int *highenough);
		static MATHFUNCSDLL_API double _stdcall CalcChiValsP(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals);
		static MATHFUNCSDLL_API int _stdcall GrowMChiWinP(int LO, int RO, int LenXoverSeq, int HWindowWidth, int TWin, int MaxY, int LS, int A, int C, int MaxFailCount, double *MPV, int *WinWin, double *MChi, int *TopL, int *TopR, int *TopLO, int *TopRO, unsigned char *Scores);
		static MATHFUNCSDLL_API int _stdcall GrowMChiWin2P(int LO, int RO, int LenXoverSeq, int HWindowWidth, int TWin, int MaxY, int LS, int a, int c, int MaxFailCount, double *MPV, int *WinWin, double *MChi, int *TopL, int *TopR, int *TopLO, int *TopRO, unsigned char *Scores, unsigned char *mdmap);
		static MATHFUNCSDLL_API int _stdcall GrowMChiWinP2(int MaxABWin, int LO, int RO, int LenXoverSeq, int HWindowWidth, int TWin, int MaxY, int LS, int A, int C, int MaxFailCount, double *MPV, int *WinWin, double *MChi, int *TopL, int *TopR, int *TopLO, int *TopRO, unsigned char *Scores, float *chitable, int *chimap);
		static MATHFUNCSDLL_API int _stdcall GrowMChiWin2P2(int MaxABWin, int LO, int RO, int LenXoverSeq, int HWindowWidth, int TWin, int MaxY, int LS, int a, int c, int MaxFailCount, double *MPV, int *WinWin, double *MChi, int *TopL, int *TopR, int *TopLO, int *TopRO, unsigned char *Scores, unsigned char *mdmap, float *chitable, int *chimap);
		static MATHFUNCSDLL_API int _stdcall FindMChiP(int LenSeq, int LenXoverSeq, int *MaxX, short int *MaxY, double *MChi, double *ChiVals);
		static MATHFUNCSDLL_API int _stdcall WinScoreCalcP(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int Seq1, int Seq2, int Seq3, unsigned char *Scores, int *XDiffPos, short int *SeqNum, int *WinScores);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqCP(int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, short int *seqnum, int *xdiffpos, int *xposdiff);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqDP(int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, short int *seqnum, int *xdiffpos, int *xposdiff);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqDP2(int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, short int *seqnum, int *xdiffpos);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqDP3(int UBFSS, int ubcs1, unsigned char *FSSRDP, unsigned char *CS, int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, int *LXOS, int UBXDP, int *XDP);

		static MATHFUNCSDLL_API int _stdcall WinScoreCalc4P(int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int Seq1, int Seq2, int Seq3, unsigned char *Scores, int *XDiffPos, short int *SeqNum, int *WinScores);
		static MATHFUNCSDLL_API int _stdcall WinScoreCalc4P2(int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int Seq1, int Seq2, int Seq3, unsigned char *Scores, int *XDiffPos, short int *SeqNum, int *WinScores);

		static MATHFUNCSDLL_API int _stdcall SmoothChiValsP(int LenXoverSeq, int LenSeq, double *ChiVals, double *SmoothChi);
		static MATHFUNCSDLL_API double _stdcall CalcChiVals4P(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, int *BanWins);
		static MATHFUNCSDLL_API double _stdcall FastBootDistP(int df, int reps, int nextno, int lenseq, float *diffsx, float *validx, int *weightmod, short int *seqnum, float *distance);
		static MATHFUNCSDLL_API double _stdcall CalcChiVals3P(int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals);
		static MATHFUNCSDLL_API double _stdcall FastBootDistIP(int df, int reps, int nextno, int lenseq, float *diffsx, float *validx, int *weightmod, short int *seqnum, float *distance);
		static MATHFUNCSDLL_API double _stdcall FastBootDistIP7(int df, int reps, int nextno, int lenseq, float *diffsx, float *validx, int *weightmod, short int *seqnum, float *distance);

		static MATHFUNCSDLL_API double _stdcall FastBootDistIP6(int dfx, int repsx, int nextnox, int lenseqx, float *dx, float *vx, int UBWM1, int UBWM2, int *wm, int UBSN1, int UBSN2, short int *sn, int UBD1, int UBD2, float *dist);
		static MATHFUNCSDLL_API double _stdcall FastBootDistIP5(int df, int reps, int nextno, int lenseq, int *weightmod, short int *seqnum, float *distance);
		static MATHFUNCSDLL_API double _stdcall ProbCalcP(double *fact, int xoverlength, int numincommon, double indprob, int lenxoverseq);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqP(int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int seq3, short int spacerno, short int *seqnum, short int *xoverseqnum, char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer);
		static MATHFUNCSDLL_API int _stdcall XOHomologyP(short int inlyer, int lenstrainseq, int lenxoverseq, short int xoverwindow, char *xoverseqnumw, int *xoverhomologynum);
		static MATHFUNCSDLL_API int _stdcall XOHomologyP2(short int inlyer, int lenstrainseq, int lenxoverseq, short int xoverwindow, char *xoverseqnumw, int *xoverhomologynum);

		static MATHFUNCSDLL_API int _stdcall FindFirstCOP(int x, int MedHomol, int  HighHomol, int LenXOverSeq, int UBXOHN, int *XOverHomologyNum);
		static MATHFUNCSDLL_API double _stdcall ViterbiCP(int SLen, int NumberAB, int NumberXY, double *OptXY, unsigned char *RecodeB, double *LaticeXY, double *TransitionM2, double *EmissionM2, double *LaticeAB);
		static MATHFUNCSDLL_API float _stdcall ViterbiCPF(int SLen, int NumberAB, int NumberXY, float *OptXY, unsigned char *RecodeB, float *LaticeXY, float *TransitionM2, float *EmissionM2, float *LaticeAB);
		static MATHFUNCSDLL_API int _stdcall FindNextP(int UBXOHN, int StartPosX, int HighHomol, int MedHomol, int LowHomol, int LenXoverSeq, int xoverwindow, int *XOverHomologyNum);
		static MATHFUNCSDLL_API int _stdcall FindNextPB(int UBXOHN, int StartPosX, int HighHomol, int MedHomol, int LowHomol, int LenXoverSeq, int xoverwindow, int *XOverHomologyNum);
		static MATHFUNCSDLL_API double _stdcall DefineEventP(int ShortOutFlag, int LongWindedFlag, int MedHomol, int HighHomol, int LowHomol, int TargetX, int CircularFlag, int XX, int  XOverWindow, int  lenseq, int  LenXoverSeq, int  SeqDaughter, int  SeqMinorP, int *EndFlag, int  *Be, int  *En, int  *NCommon, int  *XOverLength, char *XOverSeqNum, int *XDiffPos, int *XOverHomologyNum);
		static MATHFUNCSDLL_API int _stdcall CalcKMaxP(short int GCMissmatchPen, int XOLen, short int MCFlag, int MCCorrection, double LowestProb, double *pco, int *HiFragScore, double *critval, double *MissPen, double *LL, double *KMax, int *NDiff, int *highenough);
		static MATHFUNCSDLL_API int _stdcall CollapseNodesXP(int NextNo, float CutOff, int *T, float *DLen, float *CMat, float *AMat, float *TraceBak);
		static MATHFUNCSDLL_API int _stdcall CollapseNodesXP2(int NextNo, float CutOff, int *T, float *DLen, float *CMat, float *PAMat, float *TraceBak);
		static MATHFUNCSDLL_API int _stdcall CollapseNodesXP3(int NextNo, float CutOff, int *T, float *DLen, float *CMat, float *PAMat, float *TraceBak);

		static MATHFUNCSDLL_API double _stdcall ProbCalcP2(double *fact3x3, int ub3x3, int xoverlength, int numincommon, double indprob, int lenxoverseq);
		static MATHFUNCSDLL_API double _stdcall FastSimilarityBP(int df, int reps, int ISDim, int Nextno, int UBX, float *Valid, float *Diffs, short int *XCVal, short int *IntegerSeq, unsigned char *CompressValid, unsigned char *CompressDiffs, float *DistCheckB, int *weightmod);
		static MATHFUNCSDLL_API int _stdcall DoPerms3P(int LS, int SSWinLen, int SSNumPerms, int SSNumPerms2, int *PScores, char *VRandTemplate, char *VRandConv, int *PermPScores);
		static MATHFUNCSDLL_API int _stdcall MakeBanWinP(int UBBW, int Seq1, int Seq2, int Seq3, int HWindowWidth, int LS, int LenXoverSeq, int *BanWin, unsigned char *MDMap, unsigned char *MissingData, int *XPosDiff, int *XDiffPos);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqP2(int UBXSN, int UBVO, int UBCS, int XoverWindow, int lenseq, int A, int B, int C, int *AH, unsigned char *CompressedSeqs3, unsigned char *XoverSeqNumW, int *XDP, int *XPD, unsigned char *SkipTrip, int *FindSS0);
		static MATHFUNCSDLL_API int _stdcall ClearDeleteArray(int ls, int *da);
		static MATHFUNCSDLL_API int _stdcall ClearDeleteArrayB(int ls, unsigned char *da);
		static MATHFUNCSDLL_API int _stdcall DelPValsP(short int GCMaxOverlapFrags, int Y, int X, int LS, double *PVals, int *FragCount, int *FragSt, int *FragEn, int *MaxScorePos, int *DeleteArray);
		static MATHFUNCSDLL_API double _stdcall MakeSubProbP(int X, int LS, int LenXoverSeq, int BTarget, int ETarget, char *SubSeq, double *LL, double *KMax, double *MissPen, double *critval);
		static MATHFUNCSDLL_API int _stdcall MakeDeleteArrayP(int FragSt, int FragEn, int FragCount, int *DeleteArray);
		static MATHFUNCSDLL_API int _stdcall FindMissingP(int LS, int Seq1, int Seq2, int Seq3, int Z, int En, unsigned char *MissingData);
		static MATHFUNCSDLL_API int _stdcall CheckSplitP(int step, int LS, int Be, int En, int Seq1, int Seq2, int Seq3, int *Split, unsigned char *MissingData);
		static MATHFUNCSDLL_API double _stdcall GCCalcPValP(int lseq, long LXover, long *FragMaxScore, double *PVals, long *FragCount, double *KMax, double *LL, int *highenough, double *critval);
		static MATHFUNCSDLL_API double _stdcall ChiPVal2P(double X);
		static MATHFUNCSDLL_API int _stdcall GetACP(int LenXoverSeq, int LS, int MaxY, int MaxX, int TWin, int *A, int *C, unsigned char *Scores);
		static MATHFUNCSDLL_API double _stdcall CalcChiVals4P2(int UBCT, int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, int *BanWins, float *chitable);
		static MATHFUNCSDLL_API double _stdcall CalcChiVals4P3(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, int *BanWins, float *chitable);
		static MATHFUNCSDLL_API double _stdcall CalcChiValsP2(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, float *chitable);
		static MATHFUNCSDLL_API double _stdcall GetLaticePathP(int SLen, int NumberXY, double *LaticeXY, double *LaticeAB, int *LaticePath);
		static MATHFUNCSDLL_API float _stdcall GetLaticePathPF(int SLen, int NumberXY, float *LaticeXY, float *LaticeAB, int *LaticePath);

		static MATHFUNCSDLL_API double _stdcall UpdateCountsP(int SLen, int NumberABC, int NumberXY, int *LaticePath, unsigned char *RecodeB, double *TransitionCount, double *StateCount);
		static MATHFUNCSDLL_API float _stdcall UpdateCountsPF(int SLen, int NumberABC, int NumberXY, int *LaticePath, unsigned char *RecodeB, float *TransitionCount, float *StateCount);

		static MATHFUNCSDLL_API double _stdcall ForwardCP(int SLen, int NumberAB, int NumberXY, double *ValXY, double *OptXY, unsigned char *RecodeB, double *LaticeXY, double *TransitionM2, double *EmissionM2);
		static MATHFUNCSDLL_API float _stdcall ForwardCPF(int SLen, int NumberAB, int NumberXY, float *ValXY, float *OptXY, unsigned char *RecodeB, float *LaticeXY, float *TransitionM2, float *EmissionM2);

		static MATHFUNCSDLL_API double _stdcall ReverseCP(int SLen, int NumberAB, int NumberXY, double *ValXY, double *OptXY, unsigned char *RecodeB, double *LaticeXY, double *TransitionM2, double *EmissionM2);
		static MATHFUNCSDLL_API float _stdcall ReverseCPF(int SLen, int NumberAB, int NumberXY, float *ValXY, float *OptXY, unsigned char *RecodeB, float *LaticeXY, float *TransitionM2, float *EmissionM2);

		static MATHFUNCSDLL_API int _stdcall MakeBinArrayP(int Seq1, int LSeq, int Nextno, short int *Maskseq, short int *SeqNum, unsigned char *BinArray);
		static MATHFUNCSDLL_API int _stdcall MakeBinArray2P(int UBPV1, float *permvalid, int UBDP1, unsigned char *dopairs, int Seq1, int LSeq, int Nextno, short int *Maskseq, short int *SeqNum, unsigned char *BinArray, int *slookupnum, int *slookup, unsigned char *isin, int *tracesub, int *actualsize, int MinSeqSize);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqP4(int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int seq3, short int spacerno, short int *seqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray);
		//static MATHFUNCSDLL_API int _stdcall MakeBinArrayP(int ubba, int Seq1, int LSeq, int Nextno, short int *Maskseq, short int *SeqNum, unsigned char *BinArray);
		static MATHFUNCSDLL_API int _stdcall DefineEventP2(int UBXOHN, int ShortOutFlag, int LongWindedFlag, int MedHomol, int HighHomol, int LowHomol, int TargetX, int CircularFlag, int XX, int  XOverWindow, int  lenseq, int  LenXoverSeq, int  SeqDaughter, int  SeqMinorP, int *EndFlag, int  *Be, int  *En, int  *NCommon, int  *XOverLength, char *XOverSeqNum, int *XOverHomologyNum);
		static MATHFUNCSDLL_API int _stdcall MakeXPD2(int lenxoseq, int *xdiffpos, int *xposdiff);
		static MATHFUNCSDLL_API int _stdcall CopyCharArray(int ub1, int ub2, unsigned char *From, unsigned char *To);
		//static MATHFUNCSDLL_API int _stdcall FindSubSeqP5(int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int seq3, short int spacerno, short int *seqnum, short int *xoverseqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray, int *SSC);
		static MATHFUNCSDLL_API int _stdcall MakeISeq4P(int Nextno, int UBNS, int UBIS4, short int *SeqCompressor4, short int *ISeq4, char *NumSeq);
		static MATHFUNCSDLL_API int _stdcall MakeNumSeqP(int Nextno, int SLen, int UBNS, int StartPosInAlign, int EndPosInAlign, unsigned char *ConvNumSeq, short int *SeqNum, unsigned char *NumSeq);
		static MATHFUNCSDLL_API int _stdcall MakeBinArrayP4(int Seq1, int  Nextno, int UBIS4, int UBBC, int UBBA, short int *Maskseq, short int *ISeq4, unsigned char *BinArray, unsigned char *BinConverter4);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqP5(int UBXSN, int UBBA, int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int seq3, short int spacerno, short int *seqnum, short int *xoverseqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray3, unsigned char *binarray4, unsigned char *PairwiseTripArray, unsigned char *SS255RDP);
		static MATHFUNCSDLL_API int _stdcall MakeISeq3P(int Nextno, int UBNS, int UBIS4, short int *SeqCompressor4, short int *ISeq4, char *NumSeq);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqP6(int UBXSN, int UBBA, int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int seq3, short int spacerno, short int *seqnum, short int *xoverseqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray3, unsigned char *binarray4, unsigned char *binarray5, unsigned char *PairwiseTripArray, unsigned char *SS255RDP);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqP7(int UBXO1, int UBXO2, int UBXO3, int *lenxoverseq, int en, unsigned char *goong, int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int *elementseq, short int spacerno, short int *seqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray);
		static MATHFUNCSDLL_API int _stdcall GoRightP(int Seq1, int Seq2, int Seq3, int CircularFlag, int startpos, int LS, int UBMD, unsigned char *MissingData);
		static MATHFUNCSDLL_API int _stdcall GoLeftP(int Seq1, int Seq2, int Seq3, int CircularFlag, int startpos, int LS, int UBMD, unsigned char *MissingData);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqP8(int UBXO1, int UBXO2, int UBXO3, int *lenxoverseq, int en, unsigned char *goong, int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int *elementseq2, int *elementseq, short int spacerno, short int *seqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray);
		static MATHFUNCSDLL_API int _stdcall MakePairsP(int Nextno, int Da, int Ma, int Mi, int WinPP, int *RNum, int *Rlist, unsigned char *DoPairs);
		static MATHFUNCSDLL_API int _stdcall MarkRemovalsP(int Nextno, int WinPP, int Redolistsize, int *RedoList, int *RNum, int *Rlist, unsigned char *DoPairs);
		static MATHFUNCSDLL_API int _stdcall MakeBinArray3P(int SNextNo, int UBDP1, unsigned char *dopairs, int Seq1, int LSeq, int Nextno, short int *Maskseq, short int *SeqNum, unsigned char *BinArray, int *slookupnum, int *slookup, int *tracesub, int *actualsize, int MinSeqSize);
		//static MATHFUNCSDLL_API int _stdcall MakeBinArray2P(int Seq1, int LSeq, int Nextno, short int *Maskseq, short int *SeqNum, unsigned char *BinArray, int *slookupnum, int *slookup, unsigned char *isin, int *tracesub, int *actualsize, int MinSeqSize);
		static MATHFUNCSDLL_API int _stdcall SignalCountC(int Nextno, int UBXO1, int UBXO2, int AddNum, double LowestProb, short int *CurrentXOver, XOVERDEFINE *XOverlist, int *oRecombNo);
		static MATHFUNCSDLL_API int _stdcall DoSetsAP(int Nextno, int UBCX, int UBXO1, int SZ1, int lseq, int *RI, char *OLSeq, char *Sets, short int *CurrentXOver, XOVERDEFINE *XOverlist, unsigned char *DoIt, int *ISeqs);
		static MATHFUNCSDLL_API int _stdcall FillSetsP(int SZ1, int lseq, int Nextno, int UBCX, int UBXO1, int UBRL1, int UBS, int *RNum, int *RList, short int *CurrentXOver, XOVERDEFINE *XOverlist, int *RI, char *OLSeq, unsigned char *Sets);
		static MATHFUNCSDLL_API int _stdcall FillSetsP2(int UBXO2, int SZ1, int lseq, int Nextno, int UBCX, int UBXO1, int UBRL1, int UBS, int *RNum, int *RList, short int *CurrentXOver, XOVERDEFINE *XOverlist, int *RI, char *OLSeq, unsigned char *Sets);
		static MATHFUNCSDLL_API int _stdcall FillSetsP3(int BE, int EN, int UBXO2, int SZ1, int lseq, int Nextno, int UBCX, int UBXO1, int UBRL1, int UBS, int *RNum, int *RList, short int *CurrentXOver, XOVERDEFINE *XOverlist, int *RI, char *OLSeq, unsigned char *Sets);

		static MATHFUNCSDLL_API int _stdcall MakeNodeDepthC(int Nextno, int PermNextno, int UBND1, int UBDM1, int UBDD, unsigned char *DoneDist, short int *NodeDepth, float *DMat, float *TraceBak);
		static MATHFUNCSDLL_API int _stdcall MaketFSMat(int Nextno, int UBFM, int UBtFM, float *FMat, float *tFMat);
		static MATHFUNCSDLL_API int _stdcall FillRmat(int Y, int Nextno, int UBRM1, int UBRM2, int UBDM1, int UBDM2, int UBDM3, double *RMat, double *DistMat, unsigned char *ZP);
		static MATHFUNCSDLL_API int _stdcall FastBootDistIP4(int df, int reps, int nextno, int lenseq, float *diffsx, float *validx, int *weightmod, short int *seqnum, float *distance, unsigned char *fd, unsigned char *fv);
		static MATHFUNCSDLL_API unsigned char _stdcall FixOverlapsP(unsigned char DoneThisOne, int CurBegin, int CurEnd, int CurProg, int X, int Y, int MSX, float LSAdjust, int UBPD, int UBXONC1, int UBXONC2, unsigned char *ProgDo, short int *XOverNoComponent, short int *MaxXONo);
		static MATHFUNCSDLL_API int _stdcall CMaxD2P2(int incnum, int Seq1, int Seq2, int Seq3, int SBP, int EBP, int Nextno, int SLen, short int *SeqNum, short int *SeqnumX, const int *IdenticalR, int *IdenticalF, unsigned char *NucMat, int *IncSeq3, unsigned char *IncSeq2, unsigned char *IncSeq, float *E, float *d, float *VScoreMat, float *AvDist, int *TotCount);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqGCAP2(char gcindelflag, int LSeq, int seq1, int seq2, int seq3, short int *SeqNum, char *SubSeq, int *XPosDiff, int *XDiffPos, int *NDiff, unsigned char *binarray);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqGCAP3(int UBND, int UBXPD1, int UBSS1, int UBSS2, int elementnum, int *lxos, char gcindelflag, int LSeq, int seq1, int seq2, int *elementseq, short int *SeqNum, char *SubSeq, int *XPosDiff, int *XDiffPos, int *NDiff, unsigned char *binarray);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqGCAP4(int UBND, int UBXPD1, int UBSS1, int UBSS2, int elementnum, int *lxos, char gcindelflag, int LSeq, int seq1, int *elementseq2, int *elementseq, short int *SeqNum, char *SubSeq, int *XPosDiff, int *XDiffPos, int *NDiff, unsigned char *binarray);
		static MATHFUNCSDLL_API int _stdcall GetFragsP2(char *goon, int elementnum, int UBFC, int UBFS1, int UBFS2, int UBSS1, int UBSS2, short int CircularFlag, int *LenXoverSeq, int LSeq, int maxcount, char *SubSeq, int *FragSt, int *FragEn, int *FragScore, int *FragCount);
		static MATHFUNCSDLL_API int _stdcall MakeImageDataP(int bkr, int bkg, int bkb, int SX, int SY, int PosE1, int PosE0, int PosS1, int PosS0, int StS, int StSX, int CurScale, float XAD, float Min, float MR, int UBID1, int UBID2, int UBID3, int UBRM1, int UBHM1, int *HeatMap, float *RegionMat, unsigned char *ImageData);
		static MATHFUNCSDLL_API int _stdcall GetMaxFragScoreP2(int elementnum, int *LenXoverSeq, int lseq, short int CircularFlag, short int GCMissmatchPen, double *MissPen, int *MaxScorePos, int *FragMaxScore, int *FragScore, int *FragCount, int *hiscore, int *NDiffG);
		static MATHFUNCSDLL_API int _stdcall MakeSeqCatCount2P(int Nextno, int LSeq, int UBSN1, int UBSCC1, int StartPosInAlign, int EndPosInAlign, int *SeqCatCount, int *AA, short int *SeqNum, unsigned char *NucMat, unsigned char *SeqSpace, unsigned char  *NucMatB, unsigned char  *NucMat2, unsigned char *flp, unsigned char *ml, unsigned char *nl);
		static MATHFUNCSDLL_API int _stdcall ConvSimToDistP(int SLen, int Nextno, int UBDistance, int UBPermvalid, int UBFubvalid, int UBSubvalid, short int *RedoDist, float *Distance, float *FMat, float *SMat, float *PermValid, float *PermDiffs, float *Fubvalid, float *Fubdiffs, float *SubValid, float *SubDiffs);
		static MATHFUNCSDLL_API int _stdcall EraseEmptiesP(int Nextno, int UB, int UBFM, int SCO, int *ISeqs, float *FMat, float *FubValid, float *SMat, float *SubValid);
		static MATHFUNCSDLL_API int _stdcall DoAABlocksP(int xRes, int yRes, int UBPC, int UBIX, int UBID22, int UBID23, int UBID12, int UBID13, int *PCount, int *ImageX, unsigned char *ImageData, unsigned char *ImageData2);
		static MATHFUNCSDLL_API int _stdcall MakeBigMap(int IStart, int XRes, int YRes, int MBN, int TType, int TNum, int TSH, float XSize, float TSingle, int UBMB1, int UBMB2, int UBMB3, int UBMB4, HDC pict, float *MapBlocks);
		static MATHFUNCSDLL_API int _stdcall MakeBigMapB(int IStart, int XRes, int YRes, int MBN, int TType, int TNum, int TSH, float XSize, float TSingle, int UBMB1, int UBMB2, int UBMB3, int UBMB4, HDC pict, float *MapBlocks);
		static MATHFUNCSDLL_API int _stdcall DrawTreeLines(HDC Pict, int IStart, int TSHx, int TargetA, int TNum, int TType, int TDL1, int UBTD1, int UBTD2, int UBTD3, int UBTD4, float PRat, float TSingle, int *OS, float *TreeDrawB);
		static MATHFUNCSDLL_API float _stdcall GetMaxXPos(int CharLen, int TNum, int TType, int TDL0, int UBON, int UBTD1, int UBTD2, int UBTD3, int UBTD4, float PRat, short int *ONameLen, float *TreeDrawB);
		static MATHFUNCSDLL_API float _stdcall GetMaxXPosB(int UBTTS1, int UBTTS2, int UBTT, int CharLen, int TNum, int TType, int TDL0, int UBON, int UBTD1, int UBTD2, int UBTD3, int UBTD4, float PRat, short int *ONameLen, float *TreeDrawB, int *TreeTraceSeqs, int *TreeTrace);
		static MATHFUNCSDLL_API int _stdcall MakeTreeDrawB2(int UBA, int UBB, int UBC, int UBD, int UBE, float *TreeDraw, float *TreeDrawB);
		static MATHFUNCSDLL_API int _stdcall TSeqPermsP(int Seq1, int Seq2, int Seq3, int lseq, int *THold, unsigned char *tMissingData, short int *SeqNum, short int *SeqRnd);
		static MATHFUNCSDLL_API int _stdcall DoPermsXP(int LS, int SSWinLen, int SSNumPerms, char *PScores, char *VRandTemplate, char *VRandConv, int *PermPScores);
		static MATHFUNCSDLL_API int _stdcall FillIntTD(int UB, float *mindist, float *mintdist, float *adjustd, float *adjusttd, short int *IntTD, float *Distance, float *TreeDistance);
		static MATHFUNCSDLL_API int _stdcall ReadIntTD(int UB, float MinDist, float MinTDist, float AdjustD, float AdjustTD, short int *IntTD, float *Distance, float *TreeDistance);
		static MATHFUNCSDLL_API int _stdcall SeqColBlocks(HDC Pict, int UBSL, double tTYF, float TCA, double XConA, int NumSeqLines, int Targ, int VSV, int *SeqLines);
		static MATHFUNCSDLL_API int _stdcall SeqColBlocksP(HDC Pict, int UBSL, double tTYF, float TCA, double XConA, int NumSeqLines, int Targ, int VSV, int *SeqLines);
		static MATHFUNCSDLL_API int _stdcall PrintSeqs(int X1, HDC Pict, int UBST, int LOS, int Targ, int UBSL, int NumSeqLines, int StartX, int VSV, int SLFS, int SeqSpaceIncrement, int FirstSeq, int *SeqLines, char* SeqText);
		static MATHFUNCSDLL_API int _stdcall PrintSeqsP(int X1, HDC Pict, int UBST, int LOS, int Targ, int UBSL, int NumSeqLines, int StartX, int VSV, int SLFS, int SeqSpaceIncrement, int FirstSeq, int *SeqLines, char* SeqText);
		static MATHFUNCSDLL_API int _stdcall ExtraRemovalsP(int Nextno, int UBF, int UBS, int *ISeqs, int *ExtraRemove, float *FMat, float *SMat);
		static MATHFUNCSDLL_API int _stdcall FindNewX(int WinPPY, int WinPP, int Seq3, int Nextno, int *RNum);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqPB(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, int *xdiffpos, int *xposdiff, unsigned char *FSSRDP);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqPB2(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, unsigned char *FSSRDP);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqMCPB(int UBFSS, int ubcs1, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int *xdiffpos, int *xposdiff, unsigned char *FSSRDP);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqMCPB2(int UBFSS, int ubcs1, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, unsigned char *FSSRDP, int *XDiffPos);
		static MATHFUNCSDLL_API int _stdcall MakeWindowSizeP(int BEP, int ENP, int *CriticalDiff, int LenXoverSeq, double CWinFract, int CWinSize, int *HWindowWidth, int lHWindowWidth, short int CProportionFlag);
		static MATHFUNCSDLL_API int _stdcall AddPVal(int Prog, double *mtP, int HWindowWidth, int LenXoverSeq, int MCCorrection, short int MCFlag, double MChi, double LowestProb);
		static MATHFUNCSDLL_API int _stdcall MakeTWinP(unsigned char FindallFlag, int HWindowWidth, int *TWin, int LenXoverSeq);
		static MATHFUNCSDLL_API int _stdcall DestroyPeakP(int MaxY, int LS, int RO, int LO, int LenXoverSeq, double *LOT, double *SmoothChi, double *ChiVals);
		static MATHFUNCSDLL_API int _stdcall AddToMapP(int A, int S, double Win, int LS, int *APos, short int *Map);
		static MATHFUNCSDLL_API int _stdcall FastRecCheckMC(int SEN, int LongWindedFlag,double *BQPV, int EarlyBale, double UCTHresh, short int MCFlag, short int ShortOutFlag, int MCCorrection, double LowestProb, short int CircularFlag, int NextNo, int MaxABWin, int HWindowWidthX, int lHWindowWidth, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainSeq0, int CriticalDiffX, unsigned char FindallFlag, int UBFSSMC, int UBCS, int Seq1, int Seq2, int Seq3, unsigned char *CS, unsigned char *FSSMC, short int *SeqNum, int UBWS, unsigned char *Scores, int *Winscores, int *XDiffPos, double *Chivals,  int *BanWin,  unsigned char *MDMap, float *ChiTable2, int *Chimap, double *mtP,  double *SmoothChi);
		static MATHFUNCSDLL_API int _stdcall FastRecCheckMC2(int SEN, int LongWindedFlag, double *BQPV, int EarlyBale, double UCTHresh, short int MCFlag, short int ShortOutFlag, int MCCorrection, double LowestProb, short int CircularFlag, int NextNo, int MaxABWin, int HWindowWidthX, int lHWindowWidth, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainSeq0, int CriticalDiffX, unsigned char FindallFlag, int UBFSSMC, int UBCS, int Seq1, int Seq2, int Seq3, unsigned char *CS, unsigned char *FSSMC, short int *SeqNum, unsigned char *MissingData, int UBWS, unsigned char *Scores, int *Winscores, int *XDiffPos, int *XPosDiff, double *Chivals, int *BanWin, unsigned char *MDMap, float *ChiTable2, int *Chimap, double *mtP, double *SmoothChi);

		static MATHFUNCSDLL_API int _stdcall FastRecCheckChim(unsigned char *MissingData, int *XPosDiff, int *LXOS,int YP, int SEN, int LongWindedFlag, double *BQPV, int EarlyBale, double UCTHresh, short int MCFlag, short int ShortOutFlag, int MCCorrection, double LowestProb, short int CircularFlag, int NextNo, int MaxABWin, int HWindowWidthX, int lHWindowWidth, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainSeq0, int CriticalDiffX, unsigned char FindallFlag, int UBFSSRDP, int UBCS, int Seq1, int Seq2, int Seq3, unsigned char *CS, unsigned char *FSSRDP, short int *SeqNum, int UBWS, unsigned char *Scores, int *Winscores, int UBXDP, int *XDP, double *Chivals, int *BanWin, unsigned char *MDMap, float *ChiTable2, int *Chimap, double *mtP, double *SmoothChi);
		static MATHFUNCSDLL_API int _stdcall FindMChi3P(int LenSeq, int LenXoverSeq, int *MaxX, short int *MaxY, double *MChi, double *ChiVals);
		static MATHFUNCSDLL_API int _stdcall CleanChiVals(int LenXoverSeq, int LenSeq, double *ChiVals);
		static MATHFUNCSDLL_API int _stdcall CleanChiVals2(int LenXoverSeq, int LenSeq, double *ChiVals);
		static MATHFUNCSDLL_API int _stdcall SmoothChiVals3P(int LenXoverSeq, int LenSeq, double *ChiVals, double *SmoothChi);
		static MATHFUNCSDLL_API int _stdcall CleanXOSNW(int lenxoseq, int xoverwindow, int UBXO1, char *xoverseqnumw);
		static MATHFUNCSDLL_API int _stdcall GetBestMatch(int Nextno, int NumSeeds, int UBD, float *Dist, int *BestMatch);
		static MATHFUNCSDLL_API int _stdcall GetBestMatch2(int Nextno, int NumSeeds, int UBD, float *Dist, int *BestMatch, int *NIY);
		static MATHFUNCSDLL_API int _stdcall GetClosestTo(int A, int Nextno, int UBD, int *Done, float *ClosestTo, float *Dist);
		static MATHFUNCSDLL_API int _stdcall SuperDist14(int X, int Y, int UB14, int *tvd, short int *ISeq14A, short int *ISeq14B, char *CompressValid14, char *CompressDiffs14);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqPB3(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, unsigned char *FSSRDP);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqPB5(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, unsigned char *FSSRDP);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqPB6( int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, unsigned char *FSSRDP);

		static MATHFUNCSDLL_API float _stdcall Clearcut(int outlyer, int NextNo, int treetype, int nlen, int nseed, int RJ, int UBD, float *dists, char *outtree);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqPB4(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, int *xdiffpos, int *xposdiff, unsigned char *FSSRDP);
		static MATHFUNCSDLL_API int _stdcall UpdatePlotsCP(int UBAD,float ff, HDC Pict, int LSeq, short int P1, short int P2, short int P3, short int P4, int StepSize, float XFactor, float oDMax, float oPMax, int MaxHits, int *Decompress, float *PDistPlt, float *ProbPlt, int *HitPlt);
		static MATHFUNCSDLL_API int _stdcall UpdatePlotsCP2(int UBAD,float ff, HDC Pict, int LSeq, short int P1, short int P2, short int P3, short int P4, int StepSize, float XFactor, float oDMax, float oPMax, int MaxHits, int *Decompress, float *PDistPlt, float *ProbPlt, int *HitPlt, float *ll1, float *ll2, float *ll3);
		static MATHFUNCSDLL_API double _stdcall UpdateDonePVCO(double NPVal, double LPV, int Prg, int s1, int SIP, int UBXOL1, int UBDPV, short int *CurrentXOver, XOVERDEFINE *XoverList, double *DonePVCO);
		static MATHFUNCSDLL_API int _stdcall MarkDones(int Nextno, int lseq, int STA, int ENA, int A1, int A2, int A3, int UBDS1, int UBPXO, unsigned char *DoneSeq, short int *PCurrentXOver, XOVERDEFINE *PXOList);
		static MATHFUNCSDLL_API int _stdcall XOHomologyPB5(short int inlyer, int lenstrainseq, int lenxoverseq, short int xoverwindow, char *xoverseqnumw, int *xoverhomologynum);
		static MATHFUNCSDLL_API double _stdcall MakeRCompatP(int *ISeqs, int *CompMat, int WinPP, int Nextno, int *RCompat, int *RCompatB, int *InPen, int *RCats, int *RNum, int *NRNum, int *GoodC, int *DoneX, int *Rlist, int *NRList, float *FAMat, double *LDist);
		static MATHFUNCSDLL_API int _stdcall CompressTE(int lseq, unsigned char *DecompressSeq, unsigned char *TEString, int *Decompress);
		static MATHFUNCSDLL_API int _stdcall MakeVarSites(int lseq, int BPos, int EPos, int sa, int SB, int SX, int SY, int UBSN1, short int *SeqNum, int *VXPos, short int *VarSiteMap, int *VSBE);
		static MATHFUNCSDLL_API int _stdcall FakeMissing(int seq1, int lseq, int UBSN1, int UBTSN1, int *SeqnumBak, short int *Seqnum, short int *tSN);
		static MATHFUNCSDLL_API int _stdcall CountVSites(int x, int lseq, int SA, int SB, int SX, int Epos, int UBSN, short int *Seqnum);
		static MATHFUNCSDLL_API double _stdcall FastGC(int lseq, int Nextno, double PCO, int UBSN, int UBTP, int UBDP, short int *SeqNum, int *VarSites, unsigned char *Mask, unsigned char *TestPairs, unsigned char *DP);
		static MATHFUNCSDLL_API int _stdcall FastRecCheckP(int CircularFlag, int DoAllFlag, int MCCorrection, int MCFlag, int EarlyBale, double UCThresh, double LowestProb, int NextNo, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, int UBXSNW, unsigned char *CS, short int *Seqnum, int Seq1, int Seq2, int Seq3, int LenStrainSeq, int XoverWindow, short int XOverWindowX, short int *XoverSeqNum, char *XoverSeqNumW, int UBXOHN, int *XOverHomologyNum, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact, double *BQPV);
		static MATHFUNCSDLL_API int _stdcall FastRecCheckPB(int CircularFlag, int DoAllFlag, int MCCorrection, int MCFlag, int EarlyBale, double UCThresh, double LowestProb, int NextNo, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, int UBXSNW, unsigned char *CS, short int *Seqnum, int Seq1, int Seq2, int Seq3, int LenStrainSeq, int XoverWindow, short int XOverWindowX, short int *XoverSeqNum, char *XoverSeqNumW, int UBXOHN, int *XOverHomologyNum, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact, double *BQPV);

		static MATHFUNCSDLL_API int _stdcall AEFirstRDP(int Seq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DP, unsigned char *DP2,int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseq, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact);
		static MATHFUNCSDLL_API int _stdcall AEFirstRDP2(int Seq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DPO, unsigned char *DP2, int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseqO, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CSO, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDPO, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact);
		static MATHFUNCSDLL_API int _stdcall AEFirstRDP3(int oSeq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DP, unsigned char *DP2, int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseq, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact);
		static MATHFUNCSDLL_API float _stdcall DoHMMCycles(int nseed, int SLen, int HMMCycles, int LenStrainSeq0, int NumberXY, int NumberABC, unsigned char *RecodeB, double *TransitionBak, double *EmissionBak, double *InitPBak, int *LaticePathBak);
		static MATHFUNCSDLL_API float _stdcall DoHMMCyclesDetermin(int nseed, int SLen, int HMMCycles, int LenStrainSeq0, int NumberXY, int NumberABC, unsigned char *RecodeB, double *TransitionBak, double *EmissionBak, double *InitPBak, int *LaticePathBak);
		static MATHFUNCSDLL_API float _stdcall DoHMMCyclesSerial(int nseed, int SLen, int HMMCycles, int LenStrainSeq0, int NumberXY, int NumberABC, unsigned char *RecodeB, double *TransitionBak, double *EmissionBak, double *InitPBak, int *LaticePathBak);
		static MATHFUNCSDLL_API float _stdcall DoHMMCyclesPara(int nseed, int SLen, int HMMCycles, int LenStrainSeq0, int NumberXY, int NumberABC, unsigned char *RecodeB, double *TransitionBak, double *EmissionBak, double *InitPBak, int *LaticePathBak);

		static MATHFUNCSDLL_API int _stdcall AESecondRDP3(int oSeq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DP, unsigned char *DP2, int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseq, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact);
		static MATHFUNCSDLL_API int _stdcall PrimaryRDP3(int oSeq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DP, unsigned char *DP2, int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseq, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact);
		static MATHFUNCSDLL_API int _stdcall AlistRDP3(short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact);
		static MATHFUNCSDLL_API int _stdcall AlistRDP4(int ubslpv, double *StoreLPV, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact);

		static MATHFUNCSDLL_API int _stdcall AEFirstAlistRDP3(int UBDP, unsigned char *DP2, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact);
		static MATHFUNCSDLL_API int _stdcall MakeCompressSeqP(int NextNo, int UBR, unsigned char *Recoded, int UBCS, unsigned char *CompressSeq, int UBCR1, int UBCR2, unsigned char *CompressorRDP);
		static MATHFUNCSDLL_API int _stdcall MakeAListISE(int *rs, int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs, int UBPV, float *PermValid);
		static MATHFUNCSDLL_API int _stdcall MakeAListISE2(int BAL, float probdo,  int *rs, int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs, int UBPV, float *PermValid);

		static MATHFUNCSDLL_API int _stdcall MakeBALISE(int PermNextno, int NextNo, int MinSeqSize,  int UBDP, unsigned char *DoPairs, int UBPV, float *PermValid);

		static MATHFUNCSDLL_API int _stdcall MakeAListIS(int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs);
		static MATHFUNCSDLL_API int _stdcall MakeAListISP(int prg, int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs);
		static MATHFUNCSDLL_API int _stdcall MakeAListISP2(int prg, int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize, int UBAL, int UBAL2, short int *AList, int UBDP, unsigned char *DoPairs);
		static MATHFUNCSDLL_API int _stdcall MakeAListISP3(float probdo, int prg, int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize, int UBAL, int UBAL2, short int *AList, int UBDP, unsigned char *DoPairs);

		static MATHFUNCSDLL_API int _stdcall MakeBALISP3(int prg, int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize,  int UBDP, unsigned char *DoPairs);

		static MATHFUNCSDLL_API int _stdcall MakeAListOSP(int UBAL2, int BusyWithExcludes, int UBSV, float *SubValid, int sNextno, int UBTS1, int prg, int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs);
		static MATHFUNCSDLL_API int _stdcall MakeSubSeqPermP(int LS, int D, int P1, int P2, short int *PermSeqNum, int *XDiffPos, int *XPosDiff);
		static MATHFUNCSDLL_API int _stdcall MakeAListASEF(int BAL, int *rs, int oNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DP, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs);
		static MATHFUNCSDLL_API int _stdcall MakeAListASEF2(float propdo, int BAL, int *rs, int oNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DP, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs);
		static MATHFUNCSDLL_API int _stdcall GetBALF(int oNextno, int NextNo, int MinSeqSize,int UBDP, unsigned char *DP, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs);
		static MATHFUNCSDLL_API int _stdcall GetBALS(int oNextno, int NextNo, int MinSeqSize, int UBDP, unsigned char *DP, unsigned char *DP2, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs);
		                                            
		static MATHFUNCSDLL_API int _stdcall MakeAListASES(int BAL, int *rs, int oNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DP, unsigned char *DP2, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs);
		static MATHFUNCSDLL_API int _stdcall MakeAListASES2(float propdo, int BAL, int *rs, int oNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DP, unsigned char *DP2, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs);

		static MATHFUNCSDLL_API int _stdcall CheckMatrixP(int * MinS, int *ISeqs, int NextNo, int SCO, int MinSeqSize, int UBMP, unsigned char *MissPair, int UBPV, float *PermValid, int UBSV, float *SubValid, int UBF, float *FMat, float *SMat, int *ValtotF, int *ValtotS);
		static MATHFUNCSDLL_API int _stdcall DoRecode(int NextNo, int LenStrainSeq0, int UBSN, short int *SeqNum, int UBRec, unsigned char *Recoded, unsigned char *NucMat, int UBRep, unsigned char *Replace);
		static MATHFUNCSDLL_API int _stdcall DoRecodeP(int NextNo, int LenStrainSeq0, int UBSN, short int *SeqNum, int UBRec, unsigned char *Recoded, unsigned char *NucMat, int UBRep, unsigned char *Replace);
		static MATHFUNCSDLL_API int _stdcall MakeLowCI(int Y, int TargetNum, int oPermNum, int PNA, int UBMS, float *MapS, int UBPVM, float *PValMap);
		static MATHFUNCSDLL_API int _stdcall MakeHighCI(int PermNum, int Y, int TargetNum, int oPermNum, int PNA, int UBMS, float *MapS, int UBPVM, float *PValMap);

		static MATHFUNCSDLL_API int _stdcall CountNucs(int NextNo, int LenStrainSeq0, int UBSN, short int *SeqNum, unsigned char *NucMat, int UBNC, int *NucCount);
		static MATHFUNCSDLL_API int _stdcall FtoFA(int NSeqs, int LenStrainSeq0, int UBTS, int *TraceSeqs, int UBTFA, float *tFAMat, int UBFA, float *FAMat);
		static MATHFUNCSDLL_API int _stdcall MakeVarSiteMap(int SWin, int LenVarSeq, short int *VarSiteMap, float *VarSiteSmooth);
		static MATHFUNCSDLL_API int _stdcall FindOverlapP(int lenseq, int BPos2, int EPos2, int *RSize, int *OLSeq);
		static MATHFUNCSDLL_API float _stdcall FindOverlapP2(int lenseq, int BPos2, int EPos2, int *RSize, int *OLSeq);

		static MATHFUNCSDLL_API int _stdcall MakeCollecteventsC(int NextNo, int lenstrainseq0, int WinPP, int *RSize, int *OLSeq, int UBCM, int *CompMat, int UBRL, int *RList, int *RNum, int Addnum, int UBSM, float *SMatSmall, int *ISeqs, int *Trace, short int *PCurrentXOver, int UBPXO, XOVERDEFINE *PXOList, int UBCE, XOVERDEFINE *collectevents);
		static MATHFUNCSDLL_API int _stdcall TreeGroupsXP(int NextNo, char *THolder, int TLen, int NLen, char *TMatch, float *DLen);
		static MATHFUNCSDLL_API int _stdcall TransferDistP(int NSeqs, int cr, int Reps, float *tFMat, float *DstMat);
		static MATHFUNCSDLL_API int _stdcall TreeReps(int NSeqs, int Reps, int BSRndNumSeed, int NameLen, float *DstMat, int LTI, int *LTree, char *tMatch, float *DL);
		static MATHFUNCSDLL_API int _stdcall TreeRepsP(int NSeqs, int Reps, int BSRndNumSeed, int NameLen, float *DstMat, int LTI, int *LTree, char *tMatch, float *DL);
		//static MATHFUNCSDLL_API int _stdcall NEIGHBOURP(short int njoin, short int jumble, int nseed, int outgrno, int numsp, float *x, char *ot, float *coltotals);
		static MATHFUNCSDLL_API int _stdcall MaketFSMatL(int Nextno, int UBFM, int UBtFM, float *FMat, float *tFMat, int *LR);
		static MATHFUNCSDLL_API int _stdcall MakeTreeArrayXP(int nextno, float *tmat2, float *tmat2bak);
		static MATHFUNCSDLL_API int _stdcall MakeTreeArrayXP2(int nextno, float *tmat2, float *tmat2bak);
		static MATHFUNCSDLL_API int _stdcall UltraTreeDistP(double MD, int MaxCurPos, int NumberOfSeqs, float *TMat2, double *NumDone, int *DoneThis, int *AbBe, int *NodeOrder, double *MidNode, double *NodeLen);
		static MATHFUNCSDLL_API int _stdcall UltraTreeDistP2(int rr, int MaxCurPos, int NumberOfSeqs, float *TMat2, double *NumDone, int *DoneThis, int *AbBe, int *NodeOrder, double *MidNode, double *NodeLen, float *TMat);
		static MATHFUNCSDLL_API int _stdcall UltraTreeDistP3(int rr, int MaxCurPos, int NumberOfSeqs, float *TMat2, double *NumDone, int *DoneThis, int *AbBe, int *NodeOrder, double *MidNode, double *NodeLen, float *TMat);

		static MATHFUNCSDLL_API double _stdcall TreeMidP(int MaxCurPos, int NumberOfSeqs, double *NumDone, float *TMat2, int *TB, int *NodeOrder, double *MidNode, double *NodeLen);
		static MATHFUNCSDLL_API double _stdcall TreeToArrayP(short int nlen2, int nextno, int treelen, char *sholder, float *tmat, int  *nodeorder, int *donenode, int *tempnodeorder, unsigned char *rootnode, double *nodelen, double *numdone);
		static MATHFUNCSDLL_API double _stdcall TreeToArrayP2(short int nlen2, int nextno, int treelen, char *sholder, float *tmat, int  *nodeorder, int *donenode, int *tempnodeorder, double *nodelen, double *numdone);
		static MATHFUNCSDLL_API int _stdcall Tree2ArrayP(unsigned char EarlyExitFlag, int NameLen, int NumberOfSeqs, int LTree, char *T2Holder, int UBTM2, float *TMat2);
		static MATHFUNCSDLL_API int _stdcall Tree2ArrayP2(int rr, int NameLen, int NumberOfSeqs, int LTree, char *T2Holder, int UBTM2, float *TMat2);

		static MATHFUNCSDLL_API int _stdcall CleanFCMat2P(int Nextno, int UBFM, int UBFC, float *FCMat, int *invol);
		static MATHFUNCSDLL_API int _stdcall MakeNJTrees(int NSeqs, int NextNo, int *ISeqs, unsigned char *MinPair, unsigned char *SeqPair, int BSRndNumSeed, int NameLen, int LenStrainSeq0, int UBTS, int *Outlyer, int *TraceSeqs, int UBFM, float *FMat, int UBSM, float *SMat, int UBFAM, float *FAMat, int UBSAM, float *SAMat, int *LR, char *FHolder, char *SHolder, float *tFAMat, float *tSAMat);
		static MATHFUNCSDLL_API int _stdcall MakeNJTreesP(int NSeqs, int NextNo, int *ISeqs, unsigned char *MinPair, unsigned char *SeqPair, int BSRndNumSeed, int NameLen, int LenStrainSeq0, int UBTS, int *Outlyer, int *TraceSeqs, int UBFM, float *FMat, int UBSM, float *SMat, int UBFAM, float *FAMat, int UBSAM, float *SAMat, int *LR, char *FHolder, char *SHolder, float *tFAMat, float *tSAMat);
		static MATHFUNCSDLL_API int _stdcall MakeNJTreesP2(int RR, int NSeqs, int NextNo, int *ISeqs, unsigned char *MinPair, unsigned char *SeqPair, int BSRndNumSeed, int NameLen, int LenStrainSeq0, int UBTS, int *Outlyer, int *TraceSeqs, int UBFM, float *FMat, int UBSM, float *SMat, int UBFAM, float *FAMat, int UBSAM, float *SAMat, int *LR, char *FHolder, char *SHolder, float *tFAMat, float *tSAMat);
		static MATHFUNCSDLL_API int _stdcall MakeNJTreesP3(int RR, int NSeqs, int NextNo, int *ISeqs, unsigned char *MinPair, unsigned char *SeqPair, int BSRndNumSeed, int NameLen, int LenStrainSeq0, int UBTS, int *Outlyer, int *TraceSeqs, int UBFM, float *FMat, int UBSM, float *SMat, int UBFAM, float *FAMat, int UBSAM, float *SAMat, int *LR, char *FHolder, char *SHolder, float *tFAMat, float *tSAMat);
		static MATHFUNCSDLL_API double _stdcall GCXoverDP(int MCFlag, int UBPV, double *PVals, double LowestProb,  int MCCorrection, int ShortOutFlag, int CircularFlag, int GCDimSize, int lenstrainseq0, short int GCMissmatchPen, char GCIndelFlag, int Seq1, int Seq2, int Seq3, int UBFST, int *FragSt, int *FragEn, int UBFS, int *FragScore, short int *SeqNum, int UBSS, char *SubSeq, int UBMSP, int *MaxScorePos, int UBFMS, int *FragMaxScore, int *HighEnough);
		static MATHFUNCSDLL_API double _stdcall GCXoverDP2(double *BQPV, int ubcs, unsigned char *cs, int ubfss, unsigned char *fssgc, int MCFlag, int UBPV, double *PVals, double LowestProb, int MCCorrection, int ShortOutFlag, int CircularFlag, int GCDimSize, int lenstrainseq0, short int GCMissmatchPen, char GCIndelFlag, int Seq1, int Seq2, int Seq3, int UBFST, int *FragSt, int *FragEn, int UBFS, int *FragScore,  int UBSS, char *SubSeq, int UBMSP, int *MaxScorePos, int UBFMS, int *FragMaxScore, int *HighEnough);

		static MATHFUNCSDLL_API double _stdcall GCCalcPValP2(int lseq, long LXover, int *FragMaxScore, double *PVals, int *FragCount, double *KMax, double *LL, int *highenough, double *critval);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqGCAP6(int ubcs, unsigned char *cs, int ubfss, unsigned char *fssgc, char gcindelflag, int LSeq, int seq1, int seq2, int seq3, char *SubSeq, int *NDiff);
		static MATHFUNCSDLL_API int _stdcall AlistGC(char GCIndelFlag, short int GCMissmatchPen, int GCDimSize, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC);
		static MATHFUNCSDLL_API int _stdcall AlistGC2(int ubslpv, double *StoreLPV, char GCIndelFlag, short int GCMissmatchPen, int GCDimSize, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqGCAP7(int ubcs, unsigned char *cs, int ubfss, unsigned char *fssgc, char gcindelflag, int LSeq, int seq1, int seq2, int seq3, char *SubSeq, int *NDiff, int *XDiffPos, int *XPosDiff);
		static MATHFUNCSDLL_API int _stdcall cleanss(int y, int UBSS, char *SubSeq);
		static MATHFUNCSDLL_API int _stdcall AEFirstAlistGC(char GCIndelFlag, short int GCMissmatchPen, int GCDimSize, int UBDP, unsigned char *DP2, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC);
		static MATHFUNCSDLL_API int _stdcall AEFirstAlistMC(int SEN, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth, int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int UBDP, unsigned char *DP2, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainseq0, int UBCS, unsigned char *CS, int UBFSSMC, unsigned char *FSSMC, short int *SeqNum, unsigned char *MissingData, int *Chimap, float *ChiTable2);
		static MATHFUNCSDLL_API int _stdcall AEFirstAlistChi(int SEN, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth, int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int UBDP, unsigned char *DP2, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainseq0, int UBCS, unsigned char *CS, int UBFSSMC, unsigned char *FSSMC, short int *SeqNum, unsigned char *MissingData, int *Chimap, float *ChiTable2);

		
		static MATHFUNCSDLL_API int _stdcall FindBestRecSignalP(char DoneTarget, int NextNo, int UB, int UB2, double *LowP, char *DoneSeq, int *Trace, short int *PCurrentXOver, XOVERDEFINE *PXOList);
		static MATHFUNCSDLL_API int _stdcall FindBestRecSignalP2(char DoneTarget, int NextNo, int UB, int UB2, double *LowP, char *DoneSeq, int *Trace, short int *PCurrentXOver, double *TestPVs);
		static MATHFUNCSDLL_API int _stdcall MakeTestPVs(int UBDS, unsigned char *DoneSeq,int NextNo, int UB, int UB2, short int *PCurrentXOver, XOVERDEFINE *PXOList, double *TestPVs);
		static MATHFUNCSDLL_API int _stdcall UFDist(int LenStrainSeq0, int BPos3, int EPos3, int UBPV, float *PermValid, float *PermDIffs, float *BT, float *RT, int *ISeqs, int UBSN, short int *SeqNum);
		static MATHFUNCSDLL_API int _stdcall MarkOutsides(int UBDS, unsigned char *DoneSeq, int NextNo, int UB, short int *PCurrentXOver, XOVERDEFINE *PXOList);
		static MATHFUNCSDLL_API int _stdcall CheckYannP(int SEN, int NextNo, int LenStrainSeq0, int BPos, int Epos, int *ISeqs, int UBSN, short int *SeqNum, unsigned char *IsPresent, int *TraceSub, int UBXH, int UBXHMi, int UBXHMa, unsigned char *ExtraHits, unsigned char *ExtraHitsMi, unsigned char *ExtraHitsMa, int *A, int *b);
		static MATHFUNCSDLL_API int _stdcall AddjustCXO(int NextNo, int WinPP, double LowestProb, int UBDS1, int UBDS2, unsigned char *DoneSeq, int UBTD1, int UBTD2, unsigned char *TempDone, int *oRecombNo, int *RNum, int*RList, unsigned char *DoPairs, int UBTS, int *TraceSub, short int *tCurrentxover, int UBTXOL1, int UBTXOL2, XOVERDEFINE *TempXOList, short int *PCurrentXOver, int UBPXO1, int UBPXO2, XOVERDEFINE *PXOList);
		static MATHFUNCSDLL_API int _stdcall MakeLenFrag(int LenStrainSeq0, int NextNo, int ABPos, int AEPos, int *BCycle, int *BoundX, int UBSN21, short int *SeqNum2, int UBSN, short int *SeqNum);
		static MATHFUNCSDLL_API int _stdcall CleanRedo(int RWNN, int x, int NextNo, int *RedoLS, int UBRL, int *RedoList);
		static MATHFUNCSDLL_API int _stdcall MakeMoveDist(int NextNo, float *MoveDistF, float *MoveDistS, int UBFM, float *FMat, int UBSM, float *SMat);
		static MATHFUNCSDLL_API int _stdcall MakeSDMP(int NextNo, int SLen, int *SP, int *EP, int *ISeqs, int *CompMat, unsigned char *MissingData, short int *SeqNum, double *SDM, double *DistMat);
		static MATHFUNCSDLL_API int _stdcall MakeSDMP2(int NextNo, int SLen, int *SP, int *EP, int *ISeqs, int *CompMat, unsigned char *MissingData, short int *SeqNum, double *SDM, double *DistMat);
		static MATHFUNCSDLL_API int _stdcall MakeMatchMatX2P(int LSeq, int NextNo, int X, char *ContainSite, float *SMat, float *MatchMat, float *BMatch, int *BPMatch, short int *SeqNum, int *iseqs);
		static MATHFUNCSDLL_API int _stdcall MakeMatchMatX2P2(int LSeq, int NextNo, int X, char *ContainSite, float *SMat, float *MatchMat, float *BMatch, int *BPMatch, short int *SeqNum, int *iseqs);
		static MATHFUNCSDLL_API int _stdcall AlistMC2(unsigned char *wws, int StartP, int EndP, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth, int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int NextNo, int ubslpv, double *StoreLPV, short int *AList, int AListLen, unsigned char *RL, short int CircularFlag, int MCCorrection, short int MCFlag, double UCTHresh, double LowestProb, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainseq0, int UBCS, unsigned char *CS, int UBFSSMC, unsigned char *FSSMC, short int *SeqNum, int *Chimap, float *ChiTable2);
		static MATHFUNCSDLL_API int _stdcall AlistMC3(int SEN, unsigned char *wws, int StartP, int EndP, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth, int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int NextNo, int ubslpv, double *StoreLPV, short int *AList, int AListLen, unsigned char *RL, short int CircularFlag, int MCCorrection, short int MCFlag, double UCTHresh, double LowestProb, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainseq0, int UBCS, unsigned char *CS, int UBFSSMC, unsigned char *FSSMC, short int *SeqNum, unsigned char *MissingData, int *Chimap, float *ChiTable2);
		static MATHFUNCSDLL_API int _stdcall AlistChi(int SEN, unsigned char *MissingData, unsigned char *wws, int StartP, int EndP, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth, int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int NextNo, int ubslpv, double *StoreLPV, short int *AList, int AListLen, unsigned char *RL, short int CircularFlag, int MCCorrection, short int MCFlag, double UCTHresh, double LowestProb, double CWinFract, int CWinSize, short int CProportionFlag, int LenStrainseq0, int UBCS, unsigned char *CS, int UBFSSRDP, unsigned char *FSSRDP, short int *SeqNum, int *Chimap, float *ChiTable2);
		static MATHFUNCSDLL_API double _stdcall CalcChiVals5P(int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, int *BanWins);
		static MATHFUNCSDLL_API int _stdcall RecodeNucsLong(int Y,int LSeq, int UBRecoded, int UBReplace, short int *tSeqnum, unsigned char *NucMat, unsigned char *Replace, unsigned char *Recoded);

		static MATHFUNCSDLL_API int _stdcall MakeAListP(int PropTrips, int NextNo, short int *MaskSeq, int UBAL1, short int *Analysislist);
		static MATHFUNCSDLL_API int _stdcall MakeAListP2(float PropTrips, int NextNo, short int *MaskSeq, int UBAL1, short int *Analysislist);
		static MATHFUNCSDLL_API int _stdcall FindSubSeqDP6(int UBFSS, int ubcs1, unsigned char *FSSRDP, unsigned char *CS, int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, int *LXOS, int UBXDP, int *XDP, int *XPD);
		static MATHFUNCSDLL_API int _stdcall RecodeNucs(int NextNo, int LS, int UBNC, int *NucCount, int UBR, unsigned char *Replace);
		static MATHFUNCSDLL_API int _stdcall countSetBits(unsigned int n);
		static MATHFUNCSDLL_API int _stdcall countSetBits2(unsigned int n);
		static MATHFUNCSDLL_API int _stdcall MakeSSInCount(int NextNo, int UBSSIC, int *SSIncount, unsigned int *BinarySeqsL, unsigned int *BinarySeqsM, unsigned int *BinarySeqsH);
		static MATHFUNCSDLL_API int _stdcall MakeSSInCount2(int NextNo, unsigned int *v1, unsigned int *v2, unsigned int *v3, unsigned int *v4, int UBSSIC, int *SSIncount, unsigned int *BinarySeqsL, unsigned int *BinarySeqsM, unsigned int *BinarySeqsH);
		static MATHFUNCSDLL_API int _stdcall MakeSplits(int NextNo, int BinsLeft, int MaxSplitSize, int NumSplits, int x, int Y, int *BinCounts, int *BinPos, unsigned char *CloseSeq, int UBSSIC, int *SSIncount, unsigned int *BinarySeqsL, unsigned int *BinarySeqsM, unsigned int *BinarySeqsH);
		static MATHFUNCSDLL_API void _stdcall StripUnwanted(BSTR strainSeq);
		static MATHFUNCSDLL_API int _stdcall MakeListOfEdits(int UB, unsigned char *StrainSeq, int *EditList);
		static MATHFUNCSDLL_API int _stdcall CompressIt(int Element, int LenTString, int UBC, short int *Compressed, short int *SeqCompressor14, unsigned char *CS04, unsigned char *Charconvert, unsigned char *CurChar, unsigned char *TS2);
		static MATHFUNCSDLL_API int _stdcall MakeOLSeqP(int LS, int BPos, int EPos, int *RSize, int *OLSeq);
		static MATHFUNCSDLL_API int _stdcall MakeOLSeqP2(int LS, int BPos, int EPos, int *RSize, int *OLSeq, unsigned char *MDSeq);

	};

	// Isolated literal Clearcut/NJ path used by MakeNJTreesP2.  This keeps
	// the legacy tree routine's ABI-sensitive evaluation separate from the
	// surrounding matrix and event code.
	MATHFUNCSDLL_API float rdp_nj_clearcut_compat(int outlyer, int NextNo,
		int treetype, int nlen, int nseed, int RJ, int UBD, float *dists,
		char *outtree);
}
