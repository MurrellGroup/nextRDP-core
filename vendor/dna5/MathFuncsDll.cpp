// MathFuncsDll.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "MathFuncsDll.h"
//#include "dna.h"
#include "omp.h"

//#include "NEIGHBOUR.h"


//#include "ppl.h"
#include <stdexcept>
#include <malloc.h>
#include <windows.h>
#include <random>
#include <bits.h>
#include <intrin.h>
#include <string>
#include <comutil.h>
//#include <math.h>

//typedef float *vector;
//typedef long *intvector;
//
//static treeX curtree;
//
//
//struct nodex_struct {
//
//	int neighbourindex;
//	struct nodex_struct *next;
//};
//
//typedef struct nodex_struct nodex;

using namespace std;
int num_to_bits[16] = { 0, 1, 1, 2, 1, 2, 2, 3,
1, 2, 2, 3, 2, 3, 3, 4 };

namespace MathFuncs
{
	double MyMathFuncs::Add(double a, double b)
	{
		return a + b;
	}

	double MyMathFuncs::Subtract(double a, double b)
	{
		return a - b;
	}

	double MyMathFuncs::Multiply(double a, double b)
	{
		return a * b;
	}

	double MyMathFuncs::Divide(double a, double b)
	{
		if (b == 0)
		{
			throw invalid_argument("b cannot be zero!");
		}

		return a / b;
	}






	int MyMathFuncs::SuperDist14(int X, int Y, int UB14, int *tvd, short int *ISeq14A, short int *ISeq14B, char *CompressValid14, char *CompressDiffs14) {
		int Z;
		int TValid, TDiffs;
		int cv14;

		cv14 = 626;
		
		TValid = 0;
		TDiffs = 0;
		
			
	
		for (Z = 0; Z <= UB14; Z++) 
			TValid = TValid + CompressValid14[ISeq14A[Z] + ISeq14B[Z] * cv14];

		for (Z = 0; Z <= UB14; Z++) 
			TDiffs = TDiffs + CompressDiffs14[ISeq14A[Z] + ISeq14B[Z] * cv14];
		

		tvd[0] = TValid;
		tvd[1] = TDiffs;
			
		return(1);
	}
	
	
	

	

	////Phylip Neigbor
	//int MyMathFuncs::NEIGHBOURP(short int njoin, short int jumble, int nseed, int outgrno, int numsp, float *x, char *ot, float *coltotals)
	//{
	//	int tpos;
	//	float xx;
	//	xx = 0;
	//	//tpos = NEIGHBOURP2(njoin, jumble, nseed, outgrno, numsp, x, ot, coltotals);
	//	tpos = NEIGHBOUR(njoin, jumble, nseed, outgrno, numsp, x, ot, coltotals, &xx);
	//	
	//	return(tpos + 1);
	//}  /* setuptree */


	int MyMathFuncs::MakeNJTrees(int NSeqs, int NextNo, int *ISeqs, unsigned char *MinPair, unsigned char *SeqPair, int BSRndNumSeed, int NameLen, int LenStrainSeq0, int UBTS, int *Outlyer, int *TraceSeqs, int UBFM, float *FMat, int UBSM, float *SMat, int UBFAM, float *FAMat, int UBSAM, float *SAMat, int *LR, char *FHolder, char *SHolder, float *tFAMat, float *tSAMat) {

		float *tFMat, *tSMat, *TotMat, *ColTotalsF, *ColTotalsS;

		float MaxSc;
		int UBTFM, UBTSM, x, Y, Outie, Dummy, Z, LTreeF, LTreeS;
		float MinDistF, MinDistS;

		UBTFM = NSeqs;
		UBTSM = NSeqs;
		tFMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));
		tSMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));
		TotMat = (float*)calloc((NSeqs + 1), sizeof(float));


		Dummy = MaketFSMatL(NextNo, UBFM, UBTFM, FMat, tFMat, LR);


		for (x = 0; x <= NSeqs; x++) {
			for (Y = x + 1; Y <= NSeqs; Y++) {
				if (tFMat[x + x*(NSeqs + 1)] < 3) {
					if (tFMat[Y + Y*(NSeqs + 1)] < 3) {
						TotMat[x] = TotMat[x] + tFMat[x + Y*(NSeqs + 1)];
						TotMat[Y] = TotMat[Y] + tFMat[x + Y*(NSeqs + 1)];
					}
				}
			}
		}

		MaxSc = 0;
		for (x = 0; x <= NSeqs; x++) {
			if (TotMat[x] > MaxSc) {
				Outie = x;
				MaxSc = TotMat[x];
			}
		}
		free(TotMat);


		ColTotalsF = (float*)calloc((NSeqs + 1), sizeof(float));

		//LTreeF = NEIGHBOURP(1, 0, BSRndNumSeed, Outie + 1, NSeqs + 1, tFMat, FHolder, ColTotalsF);
		LTreeF =  rdp_nj_clearcut_compat(0, NSeqs, 1, 100, BSRndNumSeed, 1, NSeqs, tFMat, FHolder);
		free(ColTotalsF);
		free(tFMat);




		//tFAMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));


		Dummy = Tree2ArrayP(1, NameLen, NSeqs, LTreeF, FHolder, NSeqs, tFAMat);



		Dummy = FtoFA(NSeqs, LenStrainSeq0, UBTS, TraceSeqs, NSeqs, tFAMat, UBFAM, FAMat);

		//free(tFAMat);


		//this fills in any blanks in FAmat (in many cases only some of the sequences will be used to make famat)
		Dummy = CleanFCMat2P(NextNo, NextNo, UBFAM, FAMat, LR);


		//Use famat to make minpair 1
		MinDistF = 1000000;
		Outlyer[0] = 2;
		Outlyer[1] = 1;
		Outlyer[2] = 0;
		Z = 0;
		for (x = 0; x <= 1; x++) {
			for (Y = x + 1; Y <= 2; Y++) {
				if (FAMat[ISeqs[x] + ISeqs[Y] * (UBFAM + 1)] < MinDistF) {
					MinDistF = FAMat[ISeqs[x] + ISeqs[Y] * (UBFAM + 1)];
					MinPair[0] = Z;
					SeqPair[0] = x;
					SeqPair[1] = Y;
					SeqPair[2] = Outlyer[Z];
				}
				Z++;
			}
		}

		for (x = 0; x <= NextNo; x++)
			FAMat[x + x*(UBFAM + 1)] = 0;




		ColTotalsS = (float*)calloc((NSeqs + 1), sizeof(float));

		Dummy = MaketFSMatL(NextNo, UBSM, UBTSM, SMat, tSMat, LR);



		//LTreeS = NEIGHBOURP(1, 0, BSRndNumSeed, Outie + 1, NSeqs + 1, tSMat, SHolder, ColTotalsS);
		LTreeS = rdp_nj_clearcut_compat(0, NSeqs, 1, 100, BSRndNumSeed, 1, NSeqs, tSMat, SHolder);

		free(ColTotalsS);
		free(tSMat);

		//tSAMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));


		Dummy = Tree2ArrayP(1, NameLen, NSeqs, LTreeS, SHolder, NSeqs, tSAMat);

		Dummy = FtoFA(NSeqs, LenStrainSeq0, UBTS, TraceSeqs, NSeqs, tSAMat, UBSAM, SAMat);

		//free(tSAMat);

		Dummy = CleanFCMat2P(NextNo, NextNo, UBSAM, SAMat, LR);

		MinDistS = 1000000;
		Z = 0;
		for (x = 0; x <= 1; x++) {
			for (Y = x + 1; Y <= 2; Y++) {
				if (SAMat[ISeqs[x] + ISeqs[Y] * (UBSAM + 1)] < MinDistS) {
					MinDistS = SAMat[ISeqs[x] + ISeqs[Y] * (UBSAM + 1)];
					MinPair[1] = Z;
				}
				Z++;
			}
		}

		for (x = 0; x <= NextNo; x++)
			SAMat[x + x*(UBSAM + 1)] = 0;



		return(1);
	}
	int MyMathFuncs::MakeNJTreesP(int NSeqs, int NextNo, int *ISeqs, unsigned char *MinPair, unsigned char *SeqPair, int BSRndNumSeed, int NameLen, int LenStrainSeq0, int UBTS, int *Outlyer, int *TraceSeqs, int UBFM, float *FMat, int UBSM, float *SMat, int UBFAM, float *FAMat, int UBSAM, float *SAMat, int *LR, char *FHolder, char *SHolder, float *tFAMat, float *tSAMat) {

		float *tFMat, *tSMat;// , *TotMat;

		float MaxSc;
		int UBTFM, UBTSM, x, Y, Outie, Dummy, Z, LTreeF, LTreeS;
		float MinDistF, MinDistS;
		int SP;

		UBTFM = NSeqs;
		UBTSM = NSeqs;
		tFMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));
		tSMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));
		//TotMat = (float*)calloc((NSeqs + 1), sizeof(float));





		/*for (x = 0; x <= NSeqs; x++){
		for (Y = x + 1; Y <= NSeqs; Y++){
		if (tFMat[x + x*(NSeqs +1)] < 3){
		if (tFMat[Y + Y*(NSeqs + 1)] < 3){
		TotMat[x] = TotMat[x] + tFMat[x + Y*(NSeqs + 1)];
		TotMat[Y] = TotMat[Y] + tFMat[x + Y*(NSeqs + 1)];
		}
		}
		}
		}

		MaxSc = 0;
		for (x = 0; x <= NSeqs; x++){
		if (TotMat[x] > MaxSc){
		Outie = x;
		MaxSc = TotMat[x];
		}
		}
		free(TotMat);*/
		omp_set_num_threads(2);

#pragma omp parallel 
		{
#pragma omp sections private (x, Y, Z, Dummy)
			{
#pragma omp section
				{


					Dummy = MaketFSMatL(NextNo, UBFM, UBTFM, FMat, tFMat, LR);
					float *ColTotalsF;
					ColTotalsF = (float*)calloc((NSeqs + 1), sizeof(float));

					//LTreeF = NEIGHBOURP(1, 0, BSRndNumSeed, Outie + 1, NSeqs + 1, tFMat, FHolder, ColTotalsF);
					LTreeF = rdp_nj_clearcut_compat(0, NSeqs, 1, 100, BSRndNumSeed, 1, NSeqs, tFMat, FHolder);
					free(ColTotalsF);

					//tFAMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));


					Dummy = Tree2ArrayP(1, NameLen, NSeqs, LTreeF, FHolder, NSeqs, tFAMat);



					Dummy = FtoFA(NSeqs, LenStrainSeq0, UBTS, TraceSeqs, NSeqs, tFAMat, UBFAM, FAMat);

					//free(tFAMat);


					//this fills in any blanks in FAmat (in many cases only some of the sequences will be used to make famat)
					Dummy = CleanFCMat2P(NextNo, NextNo, UBFAM, FAMat, LR);


					//Use famat to make minpair 1
					MinDistF = 1000000;
					Outlyer[0] = 2;
					Outlyer[1] = 1;
					Outlyer[2] = 0;
					Z = 0;
					for (x = 0; x <= 1; x++) {
						for (Y = x + 1; Y <= 2; Y++) {
							if (FAMat[ISeqs[x] + ISeqs[Y] * (UBFAM + 1)] < MinDistF) {
								MinDistF = FAMat[ISeqs[x] + ISeqs[Y] * (UBFAM + 1)];
								MinPair[0] = Z;
								SeqPair[0] = x;
								SeqPair[1] = Y;
								SeqPair[2] = Outlyer[Z];
							}
							Z++;
						}
					}

					for (x = 0; x <= NextNo; x++)
						FAMat[x + x*(UBFAM + 1)] = 0;

				}
#pragma omp section
				{
					float *ColTotalsS;
					ColTotalsS = (float*)calloc((NSeqs + 1), sizeof(float));

					Dummy = MaketFSMatL(NextNo, UBSM, UBTSM, SMat, tSMat, LR);



					//LTreeS = NEIGHBOURP(1, 0, BSRndNumSeed, Outie + 1, NSeqs + 1, tSMat, SHolder, ColTotalsS);
					LTreeS = rdp_nj_clearcut_compat(0, NSeqs, 1, 100, BSRndNumSeed, 1, NSeqs, tSMat, SHolder);

					free(ColTotalsS);


					//tSAMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));


					Dummy = Tree2ArrayP(1, NameLen, NSeqs, LTreeS, SHolder, NSeqs, tSAMat);

					Dummy = FtoFA(NSeqs, LenStrainSeq0, UBTS, TraceSeqs, NSeqs, tSAMat, UBSAM, SAMat);

					//free(tSAMat);

					Dummy = CleanFCMat2P(NextNo, NextNo, UBSAM, SAMat, LR);

					MinDistS = 1000000;
					Z = 0;

					for (x = 0; x <= 1; x++) {
						for (Y = x + 1; Y <= 2; Y++) {
							if (SAMat[ISeqs[x] + ISeqs[Y] * (UBSAM + 1)] < MinDistS) {
								MinDistS = SAMat[ISeqs[x] + ISeqs[Y] * (UBSAM + 1)];
								SP = Z;
							}
							Z++;
						}
					}

					for (x = 0; x <= NextNo; x++)
						SAMat[x + x*(UBSAM + 1)] = 0;
				}
			}
		}
		MinPair[1] = SP;
		free(tFMat);
		free(tSMat);

		return(1);
	}
	int MyMathFuncs::MakeNJTreesP2(int RR, int NSeqs, int NextNo, int *ISeqs, unsigned char *MinPair, unsigned char *SeqPair, int BSRndNumSeed, int NameLen, int LenStrainSeq0, int UBTS, int *Outlyer, int *TraceSeqs,int UBFM, float *FMat, int UBSM, float *SMat, int UBFAM, float *FAMat, int UBSAM, float *SAMat, int *LR, char *FHolder, char *SHolder, float *tFAMat, float *tSAMat) {
		
		float *tFMat, *tSMat;// , *TotMat;
		
		float MaxSc;
		int UBTFM,UBTSM, x, Y, Outie, Dummy, Z, LTreeF, LTreeS;
		float MinDistF, MinDistS;
		int SP;

		UBTFM = NSeqs;
		UBTSM = NSeqs;
		tFMat = (float*)calloc((NSeqs+1)*(NSeqs+1), sizeof(float));
		tSMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));
		//TotMat = (float*)calloc((NSeqs + 1), sizeof(float));
		

		
		
    
		/*for (x = 0; x <= NSeqs; x++){
			for (Y = x + 1; Y <= NSeqs; Y++){
				if (tFMat[x + x*(NSeqs +1)] < 3){
					if (tFMat[Y + Y*(NSeqs + 1)] < 3){
						TotMat[x] = TotMat[x] + tFMat[x + Y*(NSeqs + 1)];
						TotMat[Y] = TotMat[Y] + tFMat[x + Y*(NSeqs + 1)];
					}
				}
			}
		}

		MaxSc = 0;
		for (x = 0; x <= NSeqs; x++){
			if (TotMat[x] > MaxSc){
				Outie = x;
				MaxSc = TotMat[x];
			}
		}
		free(TotMat);*/
		omp_set_num_threads(2);

#pragma omp parallel 
		{
#pragma omp sections private (x, Y, Z, Dummy)
			{
#pragma omp section
				{
					

					Dummy = MaketFSMatL(NextNo, UBFM, UBTFM, FMat, tFMat, LR);
					float *ColTotalsF;
					ColTotalsF = (float*)calloc((NSeqs + 1), sizeof(float));

					//LTreeF = NEIGHBOURP(1, 0, BSRndNumSeed, Outie + 1, NSeqs + 1, tFMat, FHolder, ColTotalsF);
					LTreeF = rdp_nj_clearcut_compat(0, NSeqs, 1, 100, BSRndNumSeed, 1, NSeqs, tFMat, FHolder);
					free(ColTotalsF);

					//tFAMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));


					Dummy = Tree2ArrayP2(RR, NameLen, NSeqs, LTreeF, FHolder, NSeqs, tFAMat);



					Dummy = FtoFA(NSeqs, LenStrainSeq0, UBTS, TraceSeqs, NSeqs, tFAMat, UBFAM, FAMat);

					//free(tFAMat);


					//this fills in any blanks in FAmat (in many cases only some of the sequences will be used to make famat)
					Dummy = CleanFCMat2P(NextNo, NextNo, UBFAM, FAMat, LR);


					//Use famat to make minpair 1
					MinDistF = 1000000;
					Outlyer[0] = 2;
					Outlyer[1] = 1;
					Outlyer[2] = 0;
					Z = 0;
					for (x = 0; x <= 1; x++) {
						for (Y = x + 1; Y <= 2; Y++) {
							if (FAMat[ISeqs[x] + ISeqs[Y] * (UBFAM + 1)] < MinDistF) {
								MinDistF = FAMat[ISeqs[x] + ISeqs[Y] * (UBFAM + 1)];
								MinPair[0] = Z;
								SeqPair[0] = x;
								SeqPair[1] = Y;
								SeqPair[2] = Outlyer[Z];
							}
							Z++;
						}
					}

					for (x = 0; x <= NextNo; x++)
						FAMat[x + x*(UBFAM + 1)] = 0;

				}
#pragma omp section
				{
					float *ColTotalsS;
					ColTotalsS = (float*)calloc((NSeqs + 1), sizeof(float));

					Dummy = MaketFSMatL(NextNo, UBSM, UBTSM, SMat, tSMat, LR);



					//LTreeS = NEIGHBOURP(1, 0, BSRndNumSeed, Outie + 1, NSeqs + 1, tSMat, SHolder, ColTotalsS);
					LTreeS = rdp_nj_clearcut_compat(0, NSeqs, 1, 100, BSRndNumSeed, 1, NSeqs, tSMat, SHolder);

					free(ColTotalsS);
					

					//tSAMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));


					Dummy = Tree2ArrayP2(RR, NameLen, NSeqs, LTreeS, SHolder, NSeqs, tSAMat);

					Dummy = FtoFA(NSeqs, LenStrainSeq0, UBTS, TraceSeqs, NSeqs, tSAMat, UBSAM, SAMat);

					//free(tSAMat);

					Dummy = CleanFCMat2P(NextNo, NextNo, UBSAM, SAMat, LR);

					MinDistS = 1000000;
					Z = 0;
					
					for (x = 0; x <= 1; x++) {
						for (Y = x + 1; Y <= 2; Y++) {
							if (SAMat[ISeqs[x] + ISeqs[Y] * (UBSAM + 1)] < MinDistS) {
								MinDistS = SAMat[ISeqs[x] + ISeqs[Y] * (UBSAM + 1)];
								SP = Z;
							}
							Z++;
						}
					}

					for (x = 0; x <= NextNo; x++)
						SAMat[x + x*(UBSAM + 1)] = 0;
				}
			}
		}
		MinPair[1] = SP;
		free(tFMat);
		free(tSMat);

		return(1);
	}


	int MyMathFuncs::CleanFCMat2P(int Nextno, int UBFM, int UBFC, float *FCMat, int *invol) {
		int X, Y, os, os2, os3, os4;
		float rv, nn;
		os = UBFC + 1;
		os2 = UBFM + 1;
		nn = (float)(Nextno);
		rv = (float)(((nn * 3) - 1) / 1000);

		for (X = 0; X <= Nextno; X++) {
			if (invol[X] == 1) {
				os3 = X*os;
				for (Y = 0; Y <= Nextno; Y++) {
					os4 = Y*os;
					FCMat[X + os4] = rv;
					FCMat[Y + os3] = rv;

				}
			}

		}
		/*
		for (X = 0; X <= Nextno; X++){
		For Y = X + 1 To Nextno
		FCMat(X, Y) = (CLng(FCMat(X, Y) * 10000000)) / 10000000
		FCMat(Y, X) = FCMat(X, Y)
		SCMat(X, Y) = (CLng(SCMat(X, Y) * 10000000)) / 10000000
		SCMat(Y, X) = SCMat(X, Y)
		Next Y
		}*/


		return(1);

	}
	int MyMathFuncs::MakeTreeArrayXP(int nextno, float *tmat2, float *tmat2bak) {
		double lowd = 0;
		double lld;
		int y, x;
		double cs = 1.0;
		double holder1, holder2;
		lld = 1000;
		while (lowd < 1000) {
			lowd = 1000;
			for (y = 0; y <= nextno; y++) {
				//offset = y*(nextno+1);
				for (x = y + 1; x <= nextno; x++) {
					if (tmat2[x + y*(nextno + 1)] < lowd)
						lowd = tmat2[x + y*(nextno + 1)];
				}
			}
			if (lowd < 1000) {
				holder1 = lowd*lowd;
				for (y = 0; y <= nextno; y++) {
					//offset = y*(nextno+1);

					for (x = y + 1; x <= nextno; x++) {
						holder2 = tmat2[x + y*(nextno + 1)] * tmat2[x + y*(nextno + 1)];
						if (holder2 / 0.99999 >= holder1 && holder2*0.99999 <= holder1) {
							tmat2bak[x + y*(nextno + 1)] = (float)(cs);
							tmat2[x + y*(nextno + 1)] = 10000;
							tmat2bak[y + x*(nextno + 1)] = (float)(cs);
							tmat2[y + x*(nextno + 1)] = 10000;
						}
					}
				}
			}
			else {

				for (y = 0; y <= nextno; y++) {

					for (x = y + 1; x <= nextno; x++) {
						tmat2[x + y*(nextno + 1)] = tmat2bak[x + y*(nextno + 1)] / 1000;
						tmat2[y + x*(nextno + 1)] = tmat2[x + y*(nextno + 1)];
					}
				}

			}
			if (lld != lowd)
				cs++;
			lld = lowd;
		}
		return(1);

	}
	
	int MyMathFuncs::MakeTreeArrayXP2(int nextno, float *tmat2, float *tmat2bak) {
		float lowd = 0;
		float lld;
		int y, x, tx, ty, winx;
		float cs = 1.0;
		float holder1, holder2;
		float *cm;
		//int *cmx;

		cm = (float*)calloc((nextno+1), sizeof(float));
		//cmx = (int*)calloc((nextno + 1), sizeof(int));
		for (x = 0; x <= nextno; x++) 
			cm[x] = 1000;
		
		for (x = 0; x <= nextno; x++) {
			for (y =0; y <= nextno; y++) {
				tmat2[x + y*(nextno + 1)] = (float)(round(tmat2[x + y*(nextno + 1)] * 10000)) / 10000;
				if (y != x) {
					if (cm[x] > tmat2[x + y*(nextno + 1)]) {
						cm[x] = tmat2[x + y*(nextno + 1)];
						//cmx[x] = y;
					}
				}
			}
		}

		lld = 1000;
		int winy;
		winy = 0;
		while (lowd < 1000) {
			lowd = 1000;
			for (y = 0; y <= nextno; y++) {
				if (cm[y] < lowd) {
					lowd = cm[y];
					//winx = cmx[y];
					winy = y;
				}
				//offset = y*(nextno+1);
				/*for (x = y + 1; x <= nextno; x++) {
					if (tmat2[x + y*(nextno + 1)] < lowd)
						lowd = tmat2[x + y*(nextno + 1)];
				}*/
			}
			if (lowd < 1000) {
				holder1 = lowd;// *lowd;
				for (y = 0; y <= nextno; y++) {
					//offset = y*(nextno+1);

					for (x = y + 1; x <= nextno; x++) {
						holder2 = tmat2[x + y*(nextno + 1)];// *tmat2[x + y*(nextno + 1)];
						if (holder2 == holder1){// && (holder2 == tmat2[x + winy*(nextno + 1)] || holder2 == tmat2[winy + y*(nextno + 1)]) ) {
							tmat2bak[x + y*(nextno + 1)] = (float)(cs);
							tmat2[x + y*(nextno + 1)] = 10000;
							tmat2bak[y + x*(nextno + 1)] = (float)(cs);
							tmat2[y + x*(nextno + 1)] = 10000;
						}
					}
				}
			}
			else {

				for (y = 0; y <= nextno; y++) {

					for (x = y + 1; x <= nextno; x++) {
						tmat2[x + y*(nextno + 1)] = tmat2bak[x + y*(nextno + 1)] / 1000;
						tmat2[y + x*(nextno + 1)] = tmat2[x + y*(nextno + 1)];
					}
				}

			}
			if (lld != lowd)
				cs++;
			lld = lowd;
			
			//y = winy;
			cm[winy] = 1000;
			for (x = 0; x <= nextno; x++) {
				/*for (y = 0; y <= nextno; y++) {
					tmat2[x + y*(nextno + 1)] = (float)(round(tmat2[x + y*(nextno + 1)] * 10000)) / 10000;*/
					if (winy != x) {
						if (cm[winy] > tmat2[x + winy*(nextno + 1)]) {
							cm[winy] = tmat2[x + winy*(nextno + 1)];
							//cmx[winy] = x;
						}

					}
				//}

			}
		
		}
		free(cm);
		//free(cmx);
		return(1);

	}

	double FAR pascal TreeMid(int MaxCurPos, int NumberOfSeqs,double *NumDone, float *TMat2, int *TB, int *NodeOrder, double *MidNode, double *NodeLen){
	
	int Seq1,Seq2,nextno, Y,  Z, C, Pos,Inside;
	double MD, tdist;
	tdist = 0.0;
	//Find greatest distance between seqs in tree
	nextno = NumberOfSeqs+1;
	MD = 0.0;
	for (Seq1 = 0; Seq1 < NumberOfSeqs; Seq1++){
		for (Seq2 = Seq1+1; Seq2 <= NumberOfSeqs; Seq2++){
			if (MD < TMat2[Seq1 + Seq2*nextno]) {
				MD = TMat2[Seq1 + Seq2*nextno];
				TB[0] = Seq1;
				TB[1] = Seq2;
			}
		}
	}

	MD = MD/2;
	for (Y = 0; Y<= MaxCurPos; Y++)
		NumDone[Y] = 1;
	

	//Find Midpoint of the Tree

	for (Y = 0; Y<= MaxCurPos; Y++){
		if (NodeOrder[Y] == TB[0]){ 
			Pos = Y + 1;
			do{
				//mark the route
				if (NodeOrder[Pos] > NumberOfSeqs)
					NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
				else{
					if (NodeOrder[Pos] == TB[1])
						break;
				}
				Pos++;
			} while (NodeOrder[Pos] != TB[1]);
			
			
			
			tdist = tdist + *(NodeLen + *(NodeOrder +Y));
			if (tdist < MD){
				for (Z = Y + 1; Z <= Pos; Z++){
					if (NodeOrder[Z] > NumberOfSeqs || NodeOrder[Z] == NodeOrder[Pos]){ 
						if (NumDone[NodeOrder[Z]] == -1.0 || NodeOrder[Z] == NodeOrder[Pos]){
							NumDone[NodeOrder[Z]] = -NumDone[NodeOrder[Z]]; 
							if (tdist + NodeLen[NodeOrder[Z]] < MD)
								tdist = tdist + NodeLen[NodeOrder[Z]];
							else{
								Inside = 0;
								C = Z;
								while (C > 0){
									C--;
									if(NodeOrder[C] == NodeOrder[Z]){
										Inside = 1;
										break;
									}
								}
                            
								//midpoint found
								if (Inside == 0){
									MidNode[0] = (double)(Z);
									MidNode[1] = MD - tdist;
									MidNode[2] = NodeLen[NodeOrder[Z]] - MidNode[1];
									break;
								}
								else if (Inside == 1){
									MidNode[0] = (double)(Z);
									MidNode[2] = MD - tdist;
									MidNode[1] = NodeLen[NodeOrder[Z]] - MidNode[2];
									break;
								}
							}
						}
					}
				}
			}
			else{
				MidNode[0] = (double)(Y);
				MidNode[2] = MD;
				MidNode[1] = NodeLen[NodeOrder[Y]] - MidNode[2];
			}
			break;
		}
		else if (NodeOrder[Y] == TB[1]){
			Pos = Y + 1;
			//mark the route between the most distant sequences
			do{
				if (NodeOrder[Pos] > NumberOfSeqs)
					NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
				else if (NodeOrder[Pos] == TB[0])
					break;
				Pos++;
			}while (NodeOrder[Pos] != TB[0]);
                             
			tdist = *(NodeLen+NodeOrder[Y]);
			if (tdist < MD){
				for (Z = Y + 1; Z<=Pos; Z++){
					if (NodeOrder[Z] > NumberOfSeqs || NodeOrder[Z] == NodeOrder[Pos]){
						if (NumDone[NodeOrder[Z]] == -1 || NodeOrder[Z] == NodeOrder[Pos]){
							NumDone[NodeOrder[Z]] = -NumDone[NodeOrder[Z]];
							if (tdist + NodeLen[NodeOrder[Z]] < MD)
								tdist = tdist + NodeLen[NodeOrder[Z]];
							else{
								Inside = 0;
								C = Z;
								while (C > 0){
									C --;
									if (NodeOrder[C] == NodeOrder[Z]) {
										Inside = 1;
										break;
									}
								}
								//midpoint found
								if (Inside == 0){
									MidNode[0] = (double)(Z);
									MidNode[1] = MD - tdist;
									MidNode[2] = NodeLen[NodeOrder[Z]] - MidNode[1];
									break;
									
								}
								else if (Inside == 1){
									MidNode[0] = (double)(Z);
									MidNode[2] = MD - tdist;
									MidNode[1] = NodeLen[NodeOrder[Z]] - MidNode[2];
									break;
								}
							}
						}
					}
				}
			}
			else{
				MidNode[0] = Y;
				MidNode[2] = MD;
				MidNode[1] = NodeLen[NodeOrder[Y]] - MidNode[2];
			}
			break;
		}
	}
	return(MD);
}






int MyMathFuncs::countSetBits(unsigned int n)
{
	int nibble = 0;
	if (0 == n)
		return num_to_bits[0];

	// Find last nibble
	nibble = n & 0xf;

	// Use pre-stored values to find count
	// in last nibble plus recursively add
	// remaining nibbles.
	return num_to_bits[nibble] + countSetBits(n >> 4);
}

int MyMathFuncs::countSetBits2(unsigned int n)
{
	
	return  __popcnt(n);
}

int MyMathFuncs::MakeSSInCount(int NextNo,int UBSSIC, int *SSIncount, unsigned int *BinarySeqsL,  unsigned int *BinarySeqsM, unsigned int *BinarySeqsH)
{
	int X, Y,  DiffCount, off;
	unsigned int  L, M, H, A, B, C, D, E, F;
	off = UBSSIC + 1;
	for (X = 0; X < NextNo; X++) {
		L = BinarySeqsL[X];
		M = BinarySeqsM[X];
		H = BinarySeqsH[X];
		for (Y = X + 1; Y <= NextNo; Y++) {

			A = BinarySeqsL[Y] ^ L;
			B = BinarySeqsM[Y] ^ M;
			C = BinarySeqsH[Y] ^ H;
			
			D = __popcnt(A);	
			E = __popcnt(B);
			F = __popcnt(C);

			DiffCount = D + E + F;

			if (DiffCount > 50)
				DiffCount = 50;

			SSIncount[DiffCount + X*off]++;
			SSIncount[DiffCount + Y*off]++;


		}
			
	}
	
	
	return (1);
}





int MyMathFuncs::MakeSSInCount2(int NextNo, unsigned int *v1, unsigned int *v2, unsigned int *v3, unsigned int *v4,int UBSSIC, int *SSIncount, unsigned int *BinarySeqsL, unsigned int *BinarySeqsM, unsigned int *BinarySeqsH)
{
	int X, Y, DiffCount, off;
	unsigned int Difflong, x1, x2, x3;
	off = UBSSIC + 1;
	
	for (X = 0; X < NextNo; X++) {
		v2[0] = BinarySeqsL[X];
		v2[1] = BinarySeqsM[X];
		v2[2] = BinarySeqsH[X];
		for (Y = X + 1; Y <= NextNo; Y++) {
			v3[0] = BinarySeqsL[Y];
			v3[1] = BinarySeqsM[Y];
			v3[2] = BinarySeqsH[Y];

			v1[0] = v3[0] ^ v2[0];
			v1[1] = v3[1] ^ v2[1];
			v1[2] = v3[2] ^ v2[2];



			v4[0] = __popcnt(v1[0]);
			v4[1] = __popcnt(v1[1]);
			v4[2] = __popcnt(v1[2]);

			DiffCount = v4[0] + v4[1] + v4[2];

			if (DiffCount > 50)
				DiffCount = 50;

			SSIncount[DiffCount + X*off]++;
			SSIncount[DiffCount + Y*off]++;


		}

	}


	return (1);
}


int MyMathFuncs::MakeSplits(int NextNo, int BinsLeft, int MaxSplitSize, int NumSplits,int x, int Y, int *BinCounts, int *BinPos, unsigned char *CloseSeq, int UBSSIC, int *SSIncount, unsigned int *BinarySeqsL, unsigned int *BinarySeqsM, unsigned int *BinarySeqsH)
{
	int off;
	unsigned int Difflong, DiffCount, L, M, H, A, B, C, D, E, F, LZ, MZ, HZ;
	off = UBSSIC + 1;
	
	//for (int Y = 0; Y <= NextNo; Y++) {
	//	if (CloseSeq[Y] == 0) {
	//		if (SSIncount[x + Y*off] > 0) {
				L = BinarySeqsL[Y];
				M = BinarySeqsM[Y];
				H = BinarySeqsH[Y];
				for (int Z = 0; Z <= NextNo; Z++) {
					if (Z != Y) {
						if (CloseSeq[Z] == 0) {
							if (SSIncount[x + Z*off] > 0) {

								A = BinarySeqsL[Z] ^ L;
								B = BinarySeqsM[Z] ^ M;
								C = BinarySeqsH[Z] ^ H;
								D = __popcnt(A);
								E = __popcnt(B);
								F = __popcnt(C);
								DiffCount = D + E + F;

								if (DiffCount == x) {

									SSIncount[x + Y*off] -= 1;
									BinsLeft -= 1;

									// Remove distances to Z from SSIncount
									LZ = BinarySeqsL[Z];
									MZ = BinarySeqsM[Z];
									HZ = BinarySeqsH[Z];
									for (int b = 0; b <= NextNo; b++) {
										if (BinPos[b] == Z + 1) {
											BinPos[b] = Y + 1;
										}
										if (b != Y && b != Z) {
											A = BinarySeqsL[b] ^ LZ;
											B = BinarySeqsM[b] ^ MZ;
											C = BinarySeqsH[b] ^ HZ;

											D = __popcnt(A);
											E = __popcnt(B);
											F = __popcnt(C);
											DiffCount = D + E + F;

											if (DiffCount > 50) DiffCount = 50;
											SSIncount[DiffCount + b*off] -= 1;

										}
									}
									CloseSeq[Z] = 1;

									BinCounts[Y] += BinCounts[Z];
									BinCounts[Z] = 0;

									if (BinsLeft == NumSplits-1) {
										break;
									}
									if (BinCounts[Y] > MaxSplitSize) {
										// Remove distances to Y from all SSIncount
										BinsLeft -= 1;
										for (int b = 0; b <= NextNo; b++) {
											if (b != Y && b != Z) {
												A = L ^ BinarySeqsL[b];
												B = M ^ BinarySeqsM[b];
												C = H ^ BinarySeqsH[b];

												D = __popcnt(A);
												E = __popcnt(B);
												F = __popcnt(C);
												DiffCount = D + E + F;

												if (DiffCount > 50) DiffCount = 50;
												SSIncount[DiffCount + b*off] -= 1;
											}
										}
										CloseSeq[Y] = 1;
										break;
									}
									if (SSIncount[x + Y*off] == 0) break;
									if (BinsLeft <= NumSplits-1) {
										break;
									}
								}
								if (BinsLeft <= NumSplits-1) {
									break;
								}
							}
						}
					}
				}
		//	}
		//}
		//if (BinsLeft <= NumSplits) {
		//	break;
		//}
	//}
	return(BinsLeft);
}

int MyMathFuncs::CompressIt(int Element, int LenTString, int UBC, short int *Compressed, short int *SeqCompressor14, unsigned char *CS04, unsigned char *Charconvert, unsigned char *CurChar, unsigned char *TS2) {
	unsigned char tChar;
	int X, Y, off;
	off = UBC + 1;
	Y = -1;
	for (X = 0; X <= LenTString-1; X = X + 4) {
		Y = Y + 1;

		tChar = TS2[X];
		CurChar[0] = Charconvert[tChar];
		CS04[X + 1] = CurChar[0];
		tChar = TS2[X+1];
		CurChar[1] = Charconvert[tChar];
		CS04[X + 2] = CurChar[1];
		tChar = TS2[X + 2];
		CurChar[2] = Charconvert[tChar];
		CS04[X + 3] = CurChar[2];
		tChar = TS2[X + 3];
		CurChar[3] = Charconvert[tChar];
		CS04[X + 4] = CurChar[3];
		Compressed[Y + Element*off] = SeqCompressor14[CurChar[0] + CurChar[1]*5 + CurChar[2]*25 + CurChar[3]*125];
		//break;
	}
	return 1;
}

int MyMathFuncs::MakeListOfEdits(int UB, unsigned char *StrainSeq, int *EditList) {
	int EL = 0;  // Initialize the edit count
	unsigned char TestChar;
				 // Iterate through the string from position 1 to UB (0-based in C++)
	for (int Z = 0; Z < UB; Z++) {
		TestChar = StrainSeq[Z];  // Get the character at position Z

									   // Check if the character is '-', 'A', 'C', 'T', or 'G'
		if (TestChar == 45 || TestChar == 65 || TestChar == 67 || TestChar == 84 || TestChar == 71) {
			continue;  // Skip these characters
		}
		else {
			EL++;  // Increase the edit count
			EditList[1 + EL * 2] = Z;
			if (TestChar == 85) {
				EditList[EL*2] = 2;  // 'U' is represented by 2
			}
			else {
				EditList[EL * 2] = 1; ;  // Other invalid characters are represented by 1
			}
		}
	}

	// Add the end marker to the list
	EL++;
	EditList[EL * 2] = -1;
	return EL;  // Return the number of edits + 1 (for the end marker)
}



void MyMathFuncs::StripUnwanted(BSTR strainSeq) {
	if (strainSeq == nullptr) return;

	// Convert BSTR to a wstring for easier manipulation
	std::wstring seq(strainSeq, SysStringLen(strainSeq));

	// Perform the character replacement as required
	for (int Y = 1; Y <= 255; Y++) {
		if (Y != 45 && Y != 65 && Y != 71 && Y != 84 && Y != 67) {
			size_t Pos = 0;
			while ((Pos = seq.find(static_cast<wchar_t>(Y), Pos)) != std::wstring::npos) {
				seq[Pos] = L'-';  // Replace the character with '-'
				Pos++;
			}
		}
	}

	// Update the original BSTR with the modified string
	SysReAllocString(&strainSeq, seq.c_str());
}




double MyMathFuncs::TreeToArrayP(short int nlen2, int nextno, int treelen, char *sholder, float *tmat, int  *nodeorder, int *donenode, int *tempnodeorder, unsigned char *rootnode, double *nodelen, double *numdone) {

	int  zz, done0, dh, tcpos, maxcurpos, tnode, x, z, y, lpos, currentpos, tpos, totcount, currentnode;
	double tallydist, th1, th2, th3;
	//Set the various position counters


	maxcurpos = nextno * 3;


	lpos = 0;
	currentpos = 0;
	totcount = 0;

	currentnode = nextno;


	//Move through the treefile one character at a time and work out tree
	//distances between sequences

	while (lpos < treelen && currentpos <= maxcurpos) {
		lpos++;
		if (sholder[lpos] == 83 && currentpos <= maxcurpos) { //If character is "S" - indicates a sequence name

															  //Read in the current sequence number
			for (x = 1; x <= nlen2; x++) {
				nodeorder[currentpos] += (int)(0.1 + (sholder[lpos + x] - 48)*pow(10, nlen2 - x));//*10^(1);//nlen2 - x);// Val(Mid$(TreeOut, lpos + 1, NameLen))

			}
			//return((double)(nodeorder[currentpos]));
			//return nodeorder[currentpos]; 
			//Update position to that of the "P" character +2
			tpos = lpos + 2;

			//Find the decimal that indicaties the position of the branch length
			while (tpos < treelen) {
				if (sholder[tpos] == 46) { //  'If caracter is "."
					if (sholder[tpos - 2] != 45) { // ' ie if number is not negative
												   //NodeLen(NodeOrder(currentpos)) = Val(Mid$(TreeOut, tpos - 1, 6))
						dh = 2;
						for (x = 0; x <= 6; x++) {
							if (sholder[tpos - 2 + x] > 47 && sholder[tpos - 2 + x] < 58) {
								th1 = 6 - x - dh;
								th2 = pow(10, th1);
								th3 = sholder[tpos - 2 + x] - 48;
								nodelen[nodeorder[currentpos]] += th3*th2;//(sholder[tpos - 2 + x] - 48)*pow(10,6-x-dh);
							}
							else
								dh--;
						}
					}
					else {
						dh = 2;
						for (x = 0; x <= 6; x++) {
							if (sholder[tpos - 2 + x] > 47 && sholder[tpos - 2 + x] < 58) {
								th1 = 6 - x - dh;
								th2 = pow(10, th1);
								th3 = sholder[tpos - 2 + x] - 48;
								nodelen[nodeorder[currentpos]] += th3*th2;//(sholder[tpos - 2 + x] - 48)*pow(10,6-x-dh);
							}

							else
								dh--;
						}
						//nodelen[nodeorder[currentpos]] *= -1.0;//uncomment if neg branches alowed
						nodelen[nodeorder[currentpos]] = 0.0;

					}
					nodelen[nodeorder[currentpos]] /= 10000;

					donenode[nodeorder[currentpos]] = 1;
					break;

				}
				tpos++;

			}

			//Find next available Node slot
			while (currentpos <= maxcurpos) {
				currentpos++;
				if (nodeorder[currentpos] == 0)
					break;
			}

		}
		else if (sholder[lpos] == 40 && currentpos <= maxcurpos) { //"(" - ie a new node

			currentnode++;
			nodeorder[currentpos] = currentnode;
			//Find next available Node slot
			while (nodeorder[currentpos] != 0) {
				currentpos++;
				if (currentpos > maxcurpos)
					break;
			}

			if (currentnode != nextno + 1) {  //If its not the first internal node
											  //The idea here is to find the matching ")" to get the branch length for this node
				tnode = currentnode;
				tcpos = currentpos;
				tpos = lpos;
				do {
					tpos++;
					if (sholder[tpos] == 40) {//'"("
						tnode++; //'Increace node count
						tcpos++;
						if (tcpos > maxcurpos)
							break;
					}
					else if (sholder[tpos] == 83) {//   '"S"
						tpos += nlen2 + 7;
						tcpos++;
						if (tcpos > maxcurpos)
							break;
					}
					else if (sholder[tpos] == 41) {// '")"
						tcpos++;
						if (tcpos > maxcurpos)
							break;
						tnode--; // ' Decrease node count until the count matches currentnode-1
						if (tnode == currentnode - 1) {
							//Get the brnach length
							//tpos += 3;
							while (tpos < treelen) {
								tpos++;
								if (sholder[tpos] == 46)
									break;
							}

							if (sholder[tpos - 2] != 45) { // ' ie if number is not negative
														   //NodeLen(NodeOrder(currentpos)) = Val(Mid$(TreeOut, tpos - 1, 6))
								dh = 2;
								for (x = 0; x <= 6; x++) {
									if (sholder[tpos - 2 + x] > 47 && sholder[tpos - 2 + x] < 58) {
										th1 = 6 - x - dh;
										th2 = pow(10, th1);
										th3 = sholder[tpos - 2 + x] - 48;
										nodelen[currentnode] += th3*th2;//(sholder[tpos - 2 + x] - 48)*pow(10,6-x-dh);
									}
									else
										dh--;
								}
							}
							else {
								dh = 2;
								for (x = 0; x <= 6; x++) {
									if (sholder[tpos - 2 + x] > 47 && sholder[tpos - 2 + x] < 58) {
										th1 = 6 - x - dh;
										th2 = pow(10, th1);
										th3 = sholder[tpos - 2 + x] - 48;
										nodelen[currentnode] += th3*th2;//(sholder[tpos - 2 + x] - 48)*pow(10,6-x-dh);
									}
									else
										dh--;
								}
								//nodelen[currentnode] *= -1.0;//uncomment if negative branch lengths allowed
								nodelen[currentnode] = 0.0;
							}
							//return nodeorder[currentpos];
							nodelen[currentnode] /= 10000;
							donenode[currentnode] = 1;
							nodeorder[tcpos - 1] = currentnode;
						}
					}
					else if (sholder[tpos] == 59 || sholder[tpos] == 0)// '";" ' If end of the tree is reached
						break;

				} while (donenode[currentnode] == 0);
				//return 4;
			}
		}
		else if (sholder[lpos] == 59 || sholder[lpos] == 0) {  //'";" ' End of the tree reached so calculate the matrix

			currentnode = nextno;
			currentpos = 0;
			break;
		}

	}//do

	//finished reading in vals in nh formatted tree string

	for (z = 0; z <= maxcurpos; z++)
		tempnodeorder[z] = nodeorder[z];
	//make rootnode (a binary nextno*2 x nextno*2 encoding of the tree topology
	for (y = 0; y <= nextno * 2; y++) {
		if (tempnodeorder[y] > nextno + 1) {
			z = y + 1;
			while (tempnodeorder[z] != tempnodeorder[y]) {
				if (tempnodeorder[z] > -1)
					*(rootnode + tempnodeorder[y] + tempnodeorder[z] * (maxcurpos + 1)) = 1;
				z++;
				if (z > maxcurpos) {
					z--;
					break;
				}
			}
			//return 5;
			tempnodeorder[z] = -1;
			tempnodeorder[y] = -1;
		}
	}

	//Work out pairwise distances and write them to a matrix

	for (x = 0; x <= maxcurpos; x++) {
		if (nodelen[x] < 0.0)
			nodelen[x] = 0.0;
		else if (nodelen[x] > 1.0)
			nodelen[x] = 1.0;
	}

	done0 = 0;

	for (zz = 0; zz <= maxcurpos; zz++) {
		if (nodeorder[zz] == 0) {
			if (done0 == 0)
				done0 = 1;
			else if (done0 == 1) {
				nodeorder[zz] = nextno + 1;
				break;
			}
		}
	}



	for (y = 0; y <= maxcurpos; y++) {
		if (nodeorder[y] == nextno + 1) {
			if (y > 2)
				break;
		}
		for (z = 0; z <= maxcurpos; z++)
			numdone[z] = 1;

		if (nodeorder[y] <= nextno) {
			tallydist = nodelen[nodeorder[y]];
			for (z = y + 1; z <= maxcurpos; z++) {
				if (nodeorder[z] == nextno + 1)
					break;
				if (nodeorder[z] > nextno) {
					tallydist += nodelen[nodeorder[z]] * numdone[nodeorder[z]];
					numdone[nodeorder[z]] = numdone[nodeorder[z]] * -1.0;
				}

				else {
					tmat[nodeorder[y] + nodeorder[z] * (nextno + 1)] = (float)(tallydist + nodelen[nodeorder[z]]);
					tmat[nodeorder[z] + nodeorder[y] * (nextno + 1)] = tmat[nodeorder[y] + nodeorder[z] * (nextno + 1)];
				}
			}
		}

	}


	//return 1;

	//for (x = 0; x < nextno; x++){
	//	for(z = x + 1; z <= nextno; z++){
	//		for (y = nextno + 1; y <= nextno * 2;y++){
	//			if ((*(rootnode + y + x*(maxcurpos*2+1) ) == 1 && *(rootnode + y + z*(maxcurpos*2+1)) == 0) || (*(rootnode + y + x*(maxcurpos*2+1)) == 0 && *(rootnode + y + z*(maxcurpos*2+1)) == 1))
	//				*(tmat + x + (nextno+1) * z) +=  nodelen[y];
	//		}
	//		*(tmat + x + (nextno+1)*z) += (nodelen[x] + nodelen[z]);
	//		*(tmat + z + (nextno+1)*x) = *(tmat + x + (nextno+1)*z);
	//	}
	//}

	//free (nodelen);
	//free (tempnodeorder);
	//free (donenode);
	//free (nodeorder);

	return 6;
}

double MyMathFuncs::TreeToArrayP2(short int nlen2, int nextno, int treelen, char *sholder, float *tmat, int  *nodeorder, int *donenode, int *tempnodeorder, double *nodelen, double *numdone) {

	int  zz, done0, dh, tcpos, maxcurpos, tnode, x, z, y, lpos, currentpos, tpos, totcount, currentnode;
	double tallydist, th1, th2, th3;
	//Set the various position counters


	maxcurpos = nextno * 3+100;


	lpos = 0;
	currentpos = 0;
	totcount = 0;

	currentnode = nextno;


	//Move through the treefile one character at a time and work out tree
	//distances between sequences

	while (lpos < treelen && currentpos <= maxcurpos) {
		lpos++;
		if (sholder[lpos] == 83 && currentpos <= maxcurpos) { //If character is "S" - indicates a sequence name

															  //Read in the current sequence number
			for (x = 1; x <= nlen2; x++) {
				nodeorder[currentpos] += (int)(0.1 + (sholder[lpos + x] - 48)*pow(10, nlen2 - x));//*10^(1);//nlen2 - x);// Val(Mid$(TreeOut, lpos + 1, NameLen))

			}
			//return((double)(nodeorder[currentpos]));
			//return nodeorder[currentpos]; 
			//Update position to that of the "P" character +2
			tpos = lpos + 2;

			//Find the decimal that indicaties the position of the branch length
			while (tpos < treelen) {
				if (sholder[tpos] == 46) { //  'If caracter is "."
					if (sholder[tpos - 2] != 45) { // ' ie if number is not negative
												   //NodeLen(NodeOrder(currentpos)) = Val(Mid$(TreeOut, tpos - 1, 6))
						dh = 2;
						for (x = 0; x <= 6; x++) {
							if (sholder[tpos - 2 + x] > 47 && sholder[tpos - 2 + x] < 58) {
								th1 = 6 - x - dh;
								th2 = pow(10, th1);
								th3 = sholder[tpos - 2 + x] - 48;
								nodelen[nodeorder[currentpos]] += th3*th2;//(sholder[tpos - 2 + x] - 48)*pow(10,6-x-dh);
							}
							else
								dh--;
						}
					}
					else {
						dh = 2;
						for (x = 0; x <= 6; x++) {
							if (sholder[tpos - 2 + x] > 47 && sholder[tpos - 2 + x] < 58) {
								th1 = 6 - x - dh;
								th2 = pow(10, th1);
								th3 = sholder[tpos - 2 + x] - 48;
								nodelen[nodeorder[currentpos]] += th3*th2;//(sholder[tpos - 2 + x] - 48)*pow(10,6-x-dh);
							}

							else
								dh--;
						}
						//nodelen[nodeorder[currentpos]] *= -1.0;//uncomment if neg branches alowed
						nodelen[nodeorder[currentpos]] = 0.0;

					}
					nodelen[nodeorder[currentpos]] /= 10000;

					donenode[nodeorder[currentpos]] = 1;
					break;

				}
				tpos++;

			}

			//Find next available Node slot
			while (currentpos <= maxcurpos) {
				currentpos++;
				if (nodeorder[currentpos] == 0)
					break;
			}

		}
		else if (sholder[lpos] == 40 && currentpos <= maxcurpos) { //"(" - ie a new node

			currentnode++;
			nodeorder[currentpos] = currentnode;
			//Find next available Node slot
			while (nodeorder[currentpos] != 0) {
				currentpos++;
				if (currentpos > maxcurpos)
					break;
			}

			if (currentnode != nextno + 1) {  //If its not the first internal node
											  //The idea here is to find the matching ")" to get the branch length for this node
				tnode = currentnode;
				tcpos = currentpos;
				tpos = lpos;
				do {
					tpos++;
					if (sholder[tpos] == 40) {//'"("
						tnode++; //'Increace node count
						tcpos++;
						if (tcpos > maxcurpos)
							break;
					}
					else if (sholder[tpos] == 83) {//   '"S"
						tpos += nlen2 + 7;
						tcpos++;
						if (tcpos > maxcurpos)
							break;
					}
					else if (sholder[tpos] == 41) {// '")"
						tcpos++;
						if (tcpos > maxcurpos)
							break;
						tnode--; // ' Decrease node count until the count matches currentnode-1
						if (tnode == currentnode - 1) {
							//Get the brnach length
							//tpos += 3;
							while (tpos < treelen) {
								tpos++;
								if (sholder[tpos] == 46)
									break;
							}

							if (sholder[tpos - 2] != 45) { // ' ie if number is not negative
														   //NodeLen(NodeOrder(currentpos)) = Val(Mid$(TreeOut, tpos - 1, 6))
								dh = 2;
								for (x = 0; x <= 6; x++) {
									if (sholder[tpos - 2 + x] > 47 && sholder[tpos - 2 + x] < 58) {
										th1 = 6 - x - dh;
										th2 = pow(10, th1);
										th3 = sholder[tpos - 2 + x] - 48;
										nodelen[currentnode] += th3*th2;//(sholder[tpos - 2 + x] - 48)*pow(10,6-x-dh);
									}
									else
										dh--;
								}
							}
							else {
								dh = 2;
								for (x = 0; x <= 6; x++) {
									if (sholder[tpos - 2 + x] > 47 && sholder[tpos - 2 + x] < 58) {
										th1 = 6 - x - dh;
										th2 = pow(10, th1);
										th3 = sholder[tpos - 2 + x] - 48;
										nodelen[currentnode] += th3*th2;//(sholder[tpos - 2 + x] - 48)*pow(10,6-x-dh);
									}
									else
										dh--;
								}
								//nodelen[currentnode] *= -1.0;//uncomment if negative branch lengths allowed
								nodelen[currentnode] = 0.0;
							}
							//return nodeorder[currentpos];
							nodelen[currentnode] /= 10000;
							donenode[currentnode] = 1;
							nodeorder[tcpos - 1] = currentnode;
						}
					}
					else if (sholder[tpos] == 59 || sholder[tpos] == 0)// '";" ' If end of the tree is reached
						break;

				} while (donenode[currentnode] == 0);
				//return 4;
			}
		}
		else if (sholder[lpos] == 59 || sholder[lpos] == 0) {  //'";" ' End of the tree reached so calculate the matrix

			currentnode = nextno;
			currentpos = 0;
			break;
		}

	}//do
	
	 //finished reading in vals in nh formatted tree string

		
	//Work out pairwise distances and write them to an outfile
	//make sure all branch lengths are between 0 and 1
	for (x = 0; x <= maxcurpos; x++) {
		if (nodelen[x] < 0.0)
			nodelen[x] = 0.0;
		else if (nodelen[x] > 1.0)
			nodelen[x] = 1.0;
	}

	done0 = 0;
	//label the first empty node after node 0,1 etc
	for (zz = 0; zz <= maxcurpos; zz++) {
		if (nodeorder[zz] == 0) {
			if (done0 == 0)
				done0 = 1;
			else if (done0 == 1) {
				nodeorder[zz] = nextno+1;
				break;
			}
		}
	}



	for (y = 0; y <= maxcurpos; y++) {
		if (nodeorder[y] == maxcurpos) {
			if (y > 2)
				break;
		}
		for (z = 0; z <= maxcurpos; z++)
			numdone[z] = 1;

		if (nodeorder[y] <= nextno) {
			tallydist = nodelen[nodeorder[y]];
			for (z = y + 1; z <= maxcurpos; z++) {
				if (nodeorder[z] == nextno + 1)
					break;
				if (nodeorder[z] > nextno) {
					tallydist += nodelen[nodeorder[z]] * numdone[nodeorder[z]];
					numdone[nodeorder[z]] = numdone[nodeorder[z]] * -1.0;
				}

				else {
					tmat[nodeorder[y] + nodeorder[z] * (nextno + 1)] = (float)(tallydist + nodelen[nodeorder[z]]);
					tmat[nodeorder[z] + nodeorder[y] * (nextno + 1)] = tmat[nodeorder[y] + nodeorder[z] * (nextno + 1)];
				}
			}
		}

	}
	

	
	


	//return 1;

	//for (x = 0; x < nextno; x++){
	//	for(z = x + 1; z <= nextno; z++){
	//		for (y = nextno + 1; y <= nextno * 2;y++){
	//			if ((*(rootnode + y + x*(maxcurpos*2+1) ) == 1 && *(rootnode + y + z*(maxcurpos*2+1)) == 0) || (*(rootnode + y + x*(maxcurpos*2+1)) == 0 && *(rootnode + y + z*(maxcurpos*2+1)) == 1))
	//				*(tmat + x + (nextno+1) * z) +=  nodelen[y];
	//		}
	//		*(tmat + x + (nextno+1)*z) += (nodelen[x] + nodelen[z]);
	//		*(tmat + z + (nextno+1)*x) = *(tmat + x + (nextno+1)*z);
	//	}
	//}

	//free (nodelen);
	//free (tempnodeorder);
	//free (donenode);
	//free (nodeorder);

	return 6;
}


double MyMathFuncs::TreeMidP(int MaxCurPos, int NumberOfSeqs, double *NumDone, float *TMat2, int *TB, int *NodeOrder, double *MidNode, double *NodeLen) {

		int Seq1, Seq2, nextno, Y, Z, C, Pos, Inside;
		double MD, tdist;
		tdist = 0.0;
		//Find greatest distance between seqs in tree
		nextno = NumberOfSeqs + 1;
		MD = 0.0;
		for (Seq1 = 0; Seq1 < NumberOfSeqs; Seq1++) {
			for (Seq2 = Seq1 + 1; Seq2 <= NumberOfSeqs; Seq2++) {
				if (MD < TMat2[Seq1 + Seq2*nextno]) {
					MD = TMat2[Seq1 + Seq2*nextno];
					TB[0] = Seq1;
					TB[1] = Seq2;
				}
			}
		}
		MidNode[3] = MD;
		MD = MD / 2;
		for (Y = 0; Y <= MaxCurPos; Y++)
			NumDone[Y] = 1;


		//Find Midpoint of the Tree

		for (Y = 0; Y <= MaxCurPos; Y++) {
			if (NodeOrder[Y] == TB[0]) {
				Pos = Y + 1;
				do {
					//mark the route
					if (NodeOrder[Pos] > NumberOfSeqs)
						NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
					else {
						if (NodeOrder[Pos] == TB[1])
							break;
					}
					Pos++;
				} while (NodeOrder[Pos] != TB[1]);



				tdist = tdist + *(NodeLen + *(NodeOrder + Y));
				if (tdist < MD) {
					for (Z = Y + 1; Z <= Pos; Z++) {
						if (NodeOrder[Z] > NumberOfSeqs || NodeOrder[Z] == NodeOrder[Pos]) {
							if (NumDone[NodeOrder[Z]] == -1.0 || NodeOrder[Z] == NodeOrder[Pos]) {
								NumDone[NodeOrder[Z]] = -NumDone[NodeOrder[Z]];
								if (tdist + NodeLen[NodeOrder[Z]] < MD)
									tdist = tdist + NodeLen[NodeOrder[Z]];
								else {
									Inside = 0;
									C = Z;
									while (C > 0) {
										C--;
										if (NodeOrder[C] == NodeOrder[Z]) {
											Inside = 1;
											break;
										}
									}

									//midpoint found
									if (Inside == 0) {
										MidNode[0] = (double)(Z);
										MidNode[1] = MD - tdist;
										MidNode[2] = NodeLen[NodeOrder[Z]] - MidNode[1];
										break;
									}
									else if (Inside == 1) {
										MidNode[0] = (double)(Z);
										MidNode[2] = MD - tdist;
										MidNode[1] = NodeLen[NodeOrder[Z]] - MidNode[2];
										break;
									}
								}
							}
						}
					}
				}
				else {
					MidNode[0] = (double)(Y);
					MidNode[2] = MD;
					MidNode[1] = NodeLen[NodeOrder[Y]] - MidNode[2];
				}
				break;
			}
			else if (NodeOrder[Y] == TB[1]) {
				Pos = Y + 1;
				//mark the route between the most distant sequences
				do {
					if (NodeOrder[Pos] > NumberOfSeqs)
						NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
					else if (NodeOrder[Pos] == TB[0])
						break;
					Pos++;
				} while (NodeOrder[Pos] != TB[0]);

				tdist = *(NodeLen + NodeOrder[Y]);
				if (tdist < MD) {
					for (Z = Y + 1; Z <= Pos; Z++) {
						if (NodeOrder[Z] > NumberOfSeqs || NodeOrder[Z] == NodeOrder[Pos]) {
							if (NumDone[NodeOrder[Z]] == -1 || NodeOrder[Z] == NodeOrder[Pos]) {
								NumDone[NodeOrder[Z]] = -NumDone[NodeOrder[Z]];
								if (tdist + NodeLen[NodeOrder[Z]] < MD)
									tdist = tdist + NodeLen[NodeOrder[Z]];
								else {
									Inside = 0;
									C = Z;
									while (C > 0) {
										C--;
										if (NodeOrder[C] == NodeOrder[Z]) {
											Inside = 1;
											break;
										}
									}
									//midpoint found
									if (Inside == 0) {
										MidNode[0] = (double)(Z);
										MidNode[1] = MD - tdist;
										MidNode[2] = NodeLen[NodeOrder[Z]] - MidNode[1];
										break;

									}
									else if (Inside == 1) {
										MidNode[0] = (double)(Z);
										MidNode[2] = MD - tdist;
										MidNode[1] = NodeLen[NodeOrder[Z]] - MidNode[2];
										break;
									}
								}
							}
						}
					}
				}
				else {
					MidNode[0] = Y;
					MidNode[2] = MD;
					MidNode[1] = NodeLen[NodeOrder[Y]] - MidNode[2];
				}
				break;
			}
		}
		return(MD);
	}

	int MyMathFuncs::UltraTreeDistP(double MD, int MaxCurPos, int NumberOfSeqs, float *TMat2, double *NumDone, int *DoneThis, int *AbBe, int *NodeOrder, double *MidNode, double *NodeLen) {
		int Y, Z, C, A;
		double TallyDist, Modi;

		for (Y = 0; Y <= MaxCurPos; Y++) {
			if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])])
				break;
			AbBe[Y] = 1;
		}

		for (Y = MaxCurPos; Y >= 0; Y--) {
			if (Y == MaxCurPos) {
				while (NodeOrder[Y] == 0)
					Y--;

			}
			if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])])
				break;
			AbBe[Y] = 3;
		}

		//find the last node that was added to the tree
		for (Y = MaxCurPos; Y >= 0; Y--) {
			if (Y == MaxCurPos) {
				while (NodeOrder[Y] == 0)
					Y--;
			}
			if (AbBe[Y] == 0)
				AbBe[Y] = 2;

		}

		//now modify the tree distance matrix to reflect equal distances from the root
		//In effect its lengthening the terminal branches so that the distance between
		//sequences in a distance matrix will reflect their reletive positions in the tree
		//in relation to the midpoint "root"
		//first do the "left" part of the tree

		for (Y = 0; Y <= MaxCurPos; Y++) {
			if (AbBe[Y] > 0) {
				if (NodeOrder[Y] <= NumberOfSeqs) {
					if (DoneThis[NodeOrder[Y]] == 0) {
						DoneThis[NodeOrder[Y]] = 1;
						TallyDist = 0.0;
						for (Z = 0; Z <= MaxCurPos; Z++)
							NumDone[Z] = 1.0;

						if (NodeOrder[Y] != NodeOrder[(int)(MidNode[0])]) {
							TallyDist = TallyDist + NodeLen[NodeOrder[Y]];
							NumDone[NodeOrder[Y]] = -NumDone[NodeOrder[Y]];
							if (AbBe[Y] == 2 || AbBe[Y] == 1) {
								for (Z = Y + 1; Z <= MaxCurPos; Z++) {
									if (NodeOrder[Z] == NodeOrder[(int)(MidNode[0])]) {
										TallyDist = TallyDist + MidNode[AbBe[Y]];
										Modi = MD - TallyDist;

										C = NodeOrder[Y];
										for (A = 0; A <= NumberOfSeqs; A++) {
											if (A != C) {
												TMat2[A + C*(NumberOfSeqs + 1)] = TMat2[A + C*(NumberOfSeqs + 1)] + (float)(Modi);
												TMat2[C + A*(NumberOfSeqs + 1)] = TMat2[A + C*(NumberOfSeqs + 1)];
											}
										}
										break;
									}
									else if (NodeOrder[Z] > NumberOfSeqs) {
										TallyDist = TallyDist + NodeLen[NodeOrder[Z]] * NumDone[NodeOrder[Z]];
										NumDone[NodeOrder[Z]] = -NumDone[NodeOrder[Z]];
									}
								}
							}
							else {
								for (Z = Y - 1; Z >= 0; Z--) {

									if (NodeOrder[Z] == NodeOrder[(int)(MidNode[0])]) {
										TallyDist = TallyDist + MidNode[1];
										Modi = MD - TallyDist;
										C = NodeOrder[Y];
										for (A = 0; A <= NumberOfSeqs; A++) {
											if (A != C) {
												TMat2[A + C*(NumberOfSeqs + 1)] = TMat2[A + C*(NumberOfSeqs + 1)] + (float)(Modi);
												TMat2[C + A*(NumberOfSeqs + 1)] = TMat2[A + C*(NumberOfSeqs + 1)];
											}
										}
										break;
									}
									else if (NodeOrder[Z] > NumberOfSeqs) {
										TallyDist = TallyDist + NodeLen[NodeOrder[Z]] * NumDone[NodeOrder[Z]];
										NumDone[NodeOrder[Z]] = -NumDone[NodeOrder[Z]];
									}
								}
							}
						}
					}
				}
			}
			else if (AbBe[Y] == 0)
				break;

		}
		return(1);
	}


	int MyMathFuncs::UltraTreeDistP2(int rr, int MaxCurPos, int NumberOfSeqs, float *TMat2, double *NumDone, int *DoneThis, int *AbBe, int *NodeOrder, double *MidNode, double *NodeLen, float *TMat) {
		int X, Y, Z, C, A, Pos, lrp, rrp, dim, sf, B;
		//float TallyDist, Modi;
		//mark everything from beginning to midnode as 1
		lrp = -1;
		rrp = -1;
		int *list1, *list2, *list3, *done, *no;//, *tl, *mask
		int ln1, ln2, ln3, mcp;
		list1 = (int*)calloc(MaxCurPos + 1, sizeof(int));
		list2 = (int*)calloc(MaxCurPos + 1, sizeof(int));
		list3 = (int*)calloc(MaxCurPos + 1, sizeof(int));
		no = (int*)calloc(MaxCurPos + 1, sizeof(int));
		/*sf = 100000;
		dim = (int)((MidNode[3]) * sf);
		if (dim < 10) {
		dim = dim * 10000;
		sf = sf * 10000;
		}
		else if (dim < 100) {
		dim = dim * 1000;
		sf = sf * 1000;
		}
		else if (dim < 1000) {
		dim = dim * 100;
		sf = sf * 100;
		}
		else if (dim < 1000) {
		dim = dim * 10;
		sf = sf * 10;
		}
		dim = dim + 1;*/
		//tl = (int*)calloc(dim + 1, sizeof(int));
		done = (int*)calloc(MaxCurPos, sizeof(int));
		//mask = (int*)calloc(MaxCurPos, sizeof(int));


		ln1 = 0;
		ln2 = 0;
		ln3 = 0;
		mcp = MaxCurPos;


		for (Y = MaxCurPos; Y >= 0; Y--) {
			if (NodeOrder[Y] != 0) {
				mcp = Y;
				break;
			}
		}

		for (Y = 0; Y <= mcp; Y++) {
			if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])]) {
				lrp = Y;
				break;

			}
			//else if (NodeOrder[Y] <= NumberOfSeqs) {
			//	for (Z = 0; Z <= mcp; Z++)
			//		NumDone[NodeOrder[Z]] = 1.0;

			//	//Pos = Y + 1;
			//	for (Pos = Y + 1; Pos < mcp; Pos++) {
			//		//mark the route
			//		if (NodeOrder[Pos] == NodeOrder[(int)(MidNode[0])]) {
			//			TallyDist = 0;
			//			//add up the branch lengths to the midnode
			//			for (Z = Y + 1; Z < Pos; Z++) {
			//				if (NumDone[NodeOrder[Z]] < 0) {
			//					TallyDist = TallyDist + NodeLen[NodeOrder[Z]];
			//				}
			//			}
			//			NodeLen[NodeOrder[Y]] = 10 - TallyDist;
			//			break;
			//		}
			//		else if (NodeOrder[Pos] > NumberOfSeqs)
			//			NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
			//	}
			//}

		}
		if (NodeOrder[(int)(MidNode[0])] <= NumberOfSeqs) {
			/*ln3 = 1;
			list3[1] = NodeOrder[(int)(MidNode[0])];*/
			rrp = lrp;
			/*for (Y = mcp; Y >= 0; Y--) {
			if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])]) {

			rrp = Y;
			lrp = Y;
			break;
			}
			}*/
			/*rrp = (int)MidNode[0];
			lrp = (int)MidNode[0];
			ln1 = 0;
			for (Y = 0; Y <= NumberOfSeqs; Y++) {
			if (list3[1] != Y) {
			ln1++;
			list1[ln1] = Y;
			}

			}*/
		}
		else {
			for (Y = mcp; Y >= 0; Y--) {
				if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])]) {

					rrp = Y;
					break;
				}
				//else if (NodeOrder[Y] <= NumberOfSeqs) {
				//	for (Z = 0; Z <= mcp; Z++)
				//		NumDone[NodeOrder[Z]] = 1.0;

				//	//Pos = Y + 1;
				//	for (Pos = Y - 1; Pos > 0; Pos--) {
				//		//mark the route
				//		if (NodeOrder[Pos] == NodeOrder[(int)(MidNode[0])]) {
				//			TallyDist = 0;
				//			//add up the branch lengths to the midnode
				//			for (Z = Y - 1; Z > Pos; Z--) {
				//				if (NumDone[NodeOrder[Z]] < 0) {
				//					TallyDist = TallyDist + NodeLen[NodeOrder[Z]];
				//				}
				//			}
				//			NodeLen[NodeOrder[Y]] = 10 - TallyDist;
				//			break;
				//		}
				//		else if (NodeOrder[Pos] > NumberOfSeqs)
				//			NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
				//	}
				//}

			}
			if (lrp == -1 || rrp == -1)
				lrp = 0;
			/*for (Y = 1; Y < lrp; Y++) {
			if (NodeOrder[Y] <= NumberOfSeqs) {
			ln1++;
			list1[ln1] = NodeOrder[Y];
			}
			}
			for (Y = mcp; Y >= rrp; Y--) {
			if (NodeOrder[Y] <= NumberOfSeqs) {
			ln2++;
			list2[ln2] = NodeOrder[Y];
			}
			}


			for (Y = lrp + 1; Y < rrp; Y++) {
			if (NodeOrder[Y] <= NumberOfSeqs) {
			ln3++;
			list3[ln3] = NodeOrder[Y];
			}
			}*/
		}
		//float maxd, maxdh;
		//maxd = (MidNode[3])*sf; //NumberOfSeqs * 2;



		//if (maxd > dim)
		//	maxd = dim;
		//if (maxd > sf / 2)
		//	maxd = sf / 2;
		//if (maxd > dim)
		//	maxd = dim;
		//tl[(int)(maxd)] = 1;
		//label sequence clusters branching off the root with the max distance
		float md;
		md = NumberOfSeqs + 1;

		//for (Y = 1; Y <= ln3; Y++) {
		//	for (Z = 1; Z <= ln1; Z++) {
		//		TMat2[list3[Y] + list1[Z] * (NumberOfSeqs + 1)] = md;// maxd;
		//		TMat2[list1[Z] + list3[Y] * (NumberOfSeqs + 1)] = md;//maxd;
		//	}
		//}
		//for (Y = 1; Y <= ln3; Y++) {
		//	for (Z = 1; Z <= ln2; Z++) {
		//		TMat2[list3[Y] + list2[Z] * (NumberOfSeqs + 1)] = md;//maxd;
		//		TMat2[list2[Z] + list3[Y] * (NumberOfSeqs + 1)] = md;//maxd;
		//	}
		//}




		/*for (Y = lrp; Y <= rrp; Y++)
		mask[NodeOrder[Y]] = 1;*/

		//

		////label sequences branching off the next branch with the next lowest distance
		//maxd = maxd - 1;
		//for (Y = 1; Y <= ln1; Y++) {
		//	for (Z = 1; Z <= ln2; Z++) {
		//		TMat2[list1[Y] + list2[Z] * (NumberOfSeqs + 1)] = maxd;
		//		TMat2[list2[Z] + list1[Y] * (NumberOfSeqs + 1)] = maxd;
		//	}
		//}

		//Now rearrange nodeorder to reflect the rerooting
		//this requires finding the node for the biggest subtree with the brackets on either side of the root node brackets (excluding the old rootnode brackets)
		//int *lnl, *rnl;
		int tp, nl1, nl2, st1, st2, en1, en2, mt, cpos, wnr, wnl, winn;
		//lnl = (int*)calloc(mcp, sizeof(int));
		//rnl = (int*)calloc(mcp, sizeof(int));
		//for (X = lrp - 1; X > 0; X--)
		//	lnl[NodeOrder[X]] = 1;

		//wnr = -1;
		//for (X = rrp + 1; X < mcp; X++) {
		//	rnl[NodeOrder[X]] = 1;
		//	if (lnl[NodeOrder[X]] == 1)
		//		wnr = X;
		//}
		//if (wnr > -1) {
		//	wnl = -1;
		//	for (X =  1; X < lrp; X++) {
		//		if (rnl[NodeOrder[X]] == 1) {
		//			wnl = X;
		//			break;
		//		}
		//	}
		//}
		//free(lnl);
		//free(rnl);
		//rearrange nodeorder
		cpos = -1;
		cpos++;
		no[cpos] = NumberOfSeqs * 2 + 3;
		for (A = lrp; A <= rrp; A++) {
			cpos++;
			no[cpos] = NodeOrder[A];
			////mask the "root" clade
			//mask[cpos] = 1;
		}
		//if (wnr > -1) {//if the root branch/node is bounded by brackets need to "invert" the bracket contents
		//	for (A = rrp + 1; A <=wnr; A++) {
		//		cpos++;
		//		no[cpos] = NodeOrder[A];
		//	}
		//	for (A = lrp - 1; A >= wnl; A--) {
		//		cpos++;
		//		no[cpos] = NodeOrder[A];
		//	}
		//	for (A = wnl-1; A > 0; A--) {
		//		cpos++;
		//		no[cpos] = NodeOrder[A];
		//	}
		//	for (A = wnr+1; A < mcp; A++) {
		//		cpos++;
		//		no[cpos] = NodeOrder[A];
		//	}

		//}
		//else {
		//invert the tree around the root node
		cpos++;
		no[cpos] = NumberOfSeqs * 2 + 2;
		for (A = lrp - 1; A > 0; A--) {
			cpos++;
			no[cpos] = NodeOrder[A];
		}
		for (A = mcp - 1; A >= rrp + 1; A--) {
			cpos++;
			no[cpos] = NodeOrder[A];
		}
		cpos++;
		no[cpos] = NumberOfSeqs * 2 + 2;
		cpos++;
		no[cpos] = NumberOfSeqs * 2 + 3;
		//for (A = rrp + 1; A < mcp; A++) {
		//	cpos++;
		//	no[cpos] = NodeOrder[A];
		//}
		//for (A = 1; A <= lrp - 1; A++) {
		//	cpos++;
		//	no[cpos] = NodeOrder[A];
		//}


		//}

		mcp = cpos;
		//for (A = 0; A <= mcp; A++) {
		//	
		//	NodeOrder[A] = no[A];
		//}
		tp = lrp;

		mt = 0;

		int fp;

		if (rr == 0)
			fp = 0;
		else
			fp = 1;

		//for (mt = 0; mt <= 1; mt++) {
		for (Y = mcp; Y > 0; Y--) {
			//if (mask[Y] == mt) {
			st1 = -1;
			st2 = -1;
			en1 = -1;
			en2 = -1;
			if (no[Y] > NumberOfSeqs) {//on a node
				if (done[no[Y]] == 0) {
					done[no[Y]] = 1;
					nl1 = no[Y];
					st1 = Y;
					for (X = st1 - 1; X > 0; X--) {
						if (no[X] > NumberOfSeqs && done[no[X]] == 0) {// && mask[X]==mt) {// on next most nested node
							nl2 = no[X];
							//done[no[X]] = 1;
							st2 = X;
							break;
						}
					}
					//scan till en2
					for (X = st2 - 1; X > 0; X--) {
						if (no[X] == nl2) {// at the end of the nested node
							en2 = X;
							break;
						}
					}
					if (en2 > -1) {//if en2 is -1 it means that an inner tip pair has been reached - these get sorted at the end.
						for (X = en2 - 1; X >= 0; X--) {
							if (no[X] == nl1) {// on end next outer node
								en1 = X;
								break;
							}
						}
						if (en1 > -1) {
							//add everthing from en1 to en2 and from st2 to st1 to list 1 and everything from en2 to st2 to list 2
							ln1 = 0;
							ln2 = 0;
							for (X = en1; X < en2; X++) {
								if (no[X] <= NumberOfSeqs) {// && mask[X]==mt) {// on a sequence
									ln1++;
									list1[ln1] = no[X];
								}
							}
							for (X = st2 + 1; X < st1; X++) {
								if (no[X] <= NumberOfSeqs) {// && mask[X] == mt) {// on a sequence
									ln1++;
									list1[ln1] = no[X];
								}
							}
							for (X = en2 + 1; X < st2; X++) {
								if (no[X] <= NumberOfSeqs) {// && mask[X] == mt) {// on a sequence
									ln2++;
									list2[ln2] = no[X];
								}
							}
							if (ln1 > 0 && ln2 > 0) {
								//make average distance between list 1 and list 2 sequences
								//maxd = 0;
								//for (X = 1; X <= ln1; X++) {
								//	for (Z = 1; Z <= ln2; Z++) {
								//		maxd = maxd + TMat[list1[X] + list2[Z] * (NumberOfSeqs + 1)];
								//	}
								//}

								//maxd = maxd / (float)(ln1*ln2);
								//if (maxd >= 0.5)
								//	maxd = 0.5;
								//
								//maxd = maxd * sf;

								//if ((int)(maxd) > dim)
								//	maxd = dim;

								//maxdh = maxd;

								//if ((int)(maxd) >= 1) {

								//	if (tl[(int)maxd] != 0) {
								//		while (tl[(int)maxd] == 1 && (int)(maxd) >= 1) {
								//			maxd--;
								//			if ((int)(maxd) < 1)
								//				break;
								//		}

								//	}
								//}
								//if ((int)(maxd) <= 1) {

								//	maxd = maxdh;
								//	if ((int)(maxd) < 1)
								//		maxd = 1;
								//	//if (tl[(int)maxd] != 0) {
								//	while (tl[(int)(maxd)] == 1 && (int)(maxd) < dim) {
								//		maxd++;
								//		if ((int)(maxd) > dim)
								//			break;
								//	}

								//	//}


								//}
								//if ((int)(maxd) > dim || (int)(maxd) <1)
								//	maxd = dim;

								//tl[(int)maxd] = 1;
								if (fp == 1)
									md--;
								for (X = 1; X <= ln1; X++) {
									for (Z = 1; Z <= ln2; Z++) {
										//if(list1[X] == 8 || list1[X])
										TMat2[list1[X] + list2[Z] * (NumberOfSeqs + 1)] = md;//maxd;
										TMat2[list2[Z] + list1[X] * (NumberOfSeqs + 1)] = md;// maxd;
									}
								}

								fp = 1;
							}
						}
					}
				}
			}
			//}

		}
		//}

		//now finish off the terminal branch pairs
		//for (Y = 0; Y <= mcp; Y++) {
		//	if (no[Y] == NumberOfSeqs)
		//		no[Y] = NumberOfSeqs;
		//}

		for (Y = 0; Y < mcp; Y++) {
			/*if (Y == mcp - 10)
			Y == mcp - 10;*/
			if (no[Y] <= NumberOfSeqs && no[Y + 1] <= NumberOfSeqs) {

				//maxd = TMat[no[Y] + no[Y + 1] * (NumberOfSeqs + 1)];
				//if (maxd >= 0.5)
				//	maxd = 0.5;
				//maxd = maxd * sf;
				//if ((int)(maxd) > dim)
				//	maxd = dim;
				//maxdh = maxd;
				//if ((int)(maxd) > 1) {

				//	if (tl[(int)maxd] != 0) {
				//		while (tl[(int)maxd] == 1 && (int)(maxd) >= 1) {
				//			maxd--;
				//			if ((int)(maxd) < 1)
				//				break;
				//		}

				//	}
				//}
				//if ((int)(maxd) <= 1) {

				//	maxd = maxdh;
				//	if ((int)(maxd) < 1)
				//		maxd = 1;
				//	//if (tl[(int)maxd] != 0) {
				//	while (tl[(int)maxd] == 1 && (int)(maxd) < dim) {
				//		maxd++;
				//		if ((int)(maxd)> dim)
				//			break;
				//	}

				//	//}


				//}
				//if ((int)(maxd)> dim)
				//	maxd = dim;

				//if ((int)(maxd)<0)
				//	maxd = 1;

				//tl[(int)maxd] = 1;
				md--;
				TMat2[no[Y] + no[Y + 1] * (NumberOfSeqs + 1)] = md; //maxd;
				TMat2[no[Y + 1] + no[Y] * (NumberOfSeqs + 1)] = md;// maxd;

			}

		}
		//for (X = mcp; X >= mcp - 10; X--) {
		//	if (no[X] == 2000)
		//		no[X] = 1000;
		//}
		//for (X = 0; X <= NumberOfSeqs; X++) {
		//	for (Y =X+1; Y <= NumberOfSeqs; Y++) {
		//		if (TMat2[Y + X * (NumberOfSeqs + 1)] < 1) {
		//			TMat2[Y + X * (NumberOfSeqs + 1)] =0;
		//		}

		//	}

		//}
		//Z = 0;
		//if (sf < 1)
		//	sf = 1;

		//int cnt;
		//cnt = 0;
		for (Y = 0; Y <= NumberOfSeqs; Y++) {
			for (Z = Y + 1; Z <= NumberOfSeqs; Z++) {
				//if ((int)(TMat2[Z + Y * (NumberOfSeqs + 1)]) == X) {

				TMat[Y + Z * (NumberOfSeqs + 1)] = TMat2[Y + Z * (NumberOfSeqs + 1)] / 1000;
				TMat[Z + Y * (NumberOfSeqs + 1)] = TMat2[Y + Z * (NumberOfSeqs + 1)] / 1000;
				//}
			}

		}
		/*for (X = 0; X <= dim; X++) {
		if (tl[X] == 1) {
		cnt++;
		for (Y = 0; Y <= NumberOfSeqs; Y++) {
		for (Z = Y + 1; Z <= NumberOfSeqs; Z++) {
		if ((int)(TMat2[Z + Y * (NumberOfSeqs + 1)]) == X) {

		TMat[Y + Z * (NumberOfSeqs + 1)] = (float)cnt / 1000;
		TMat[Z + Y * (NumberOfSeqs + 1)] = (float)cnt / 1000;
		}
		}

		}

		}

		}*/


		//free(tl);
		free(done);
		//free(mask);
		free(list1);
		free(list2);
		free(list3);
		free(no);



		return(1);
	}
	int MyMathFuncs::UFDist(int LenStrainSeq0,int BPos3, int EPos3, int UBPV, float *PermValid, float *PermDIffs, float *BT, float *RT, int *ISeqs, int UBSN, short int *SeqNum) {
		int T,x,Y, I0,I1, Z, os1, os2, os3, os4;
		float D, V;
		T = -1;
		os1 = UBSN + 1;
		os4 = UBPV + 1;
		for (x = 0; x <= 1; x++){
			if (ISeqs[x] <= UBPV) {
				I0 = ISeqs[x];
				os2 = os1*I0;
				for (Y = x + 1; Y <= 2; Y++) {
					T++;
					if (ISeqs[Y] <= UBPV) {
						I1 = ISeqs[Y];
						os3 = I1*os1;
						V = 0;
						D = 0;
						if (BPos3 < EPos3) {
							for (Z = BPos3; Z <= EPos3; Z++) {
								if (SeqNum[Z + os2] != 46) {
									if (SeqNum[Z + os3] != 46) {
										if (SeqNum[Z + os2] != SeqNum[Z + os3])
											D++;

										V++;
									}
								}
							}
						}
						else {
							for (Z = BPos3; Z <= LenStrainSeq0; Z++) {
								if (SeqNum[Z + os2] != 46) {
									if (SeqNum[Z + os3] != 46) {
										if (SeqNum[Z + os2] != SeqNum[Z + os3])
											D++;

										V++;
									}
								}
							}
							for (Z = 1; Z <= EPos3; Z++) {
								if (SeqNum[Z + os2] != 46) {
									if (SeqNum[Z + os3] != 46) {
										if (SeqNum[Z + os2] != SeqNum[Z + os3])
											D++;

										V++;
									}
								}
							}
						}
						if (V > 0)
							BT[T] = D / V;
						else
							BT[T] = 10;

						if (PermValid[I0 + I1*os4] - V > 0)
							RT[T] = (PermDIffs[I0 + I1*os4] - D) / (PermValid[I0 + I1*os4] - V);
						else
							RT[T] = 10;
					}
				}
			}
		}
    


		return(1);
	}
	int MyMathFuncs::UltraTreeDistP3(int rr, int MaxCurPos, int NumberOfSeqs, float *TMat2, double *NumDone, int *DoneThis, int *AbBe, int *NodeOrder, double *MidNode, double *NodeLen, float *TMat) {
		int X, Y, Z, C, A, Pos, lrp, rrp, dim, sf, B;
		//float TallyDist, Modi;
		//mark everything from beginning to midnode as 1
		lrp = -1;
		rrp = -1;
		int *list1, *list2, *list3, *done, *no;//, *tl, *mask
		int ln1, ln2, ln3, mcp;
		list1 = (int*)calloc(MaxCurPos + 1, sizeof(int));
		list2 = (int*)calloc(MaxCurPos + 1, sizeof(int));
		list3 = (int*)calloc(MaxCurPos + 1, sizeof(int));
		no = (int*)calloc(MaxCurPos + 1, sizeof(int));
		/*sf = 100000;
		dim = (int)((MidNode[3]) * sf);
		if (dim < 10) {
			dim = dim * 10000;
			sf = sf * 10000;
		}
		else if (dim < 100) {
			dim = dim * 1000;
			sf = sf * 1000;
		}
		else if (dim < 1000) {
			dim = dim * 100;
			sf = sf * 100;
		}
		else if (dim < 1000) {
			dim = dim * 10;
			sf = sf * 10;
		}
		dim = dim + 1;*/
		//tl = (int*)calloc(dim + 1, sizeof(int));
		done = (int*)calloc(MaxCurPos, sizeof(int));
		//mask = (int*)calloc(MaxCurPos, sizeof(int));


		ln1 = 0;
		ln2 = 0;
		ln3 = 0;
		mcp = MaxCurPos;


		for (Y = MaxCurPos; Y >= 0; Y--) {
			if (NodeOrder[Y] != 0) {
				mcp = Y;
				break;
			}
		}

		for (Y = 0; Y <= mcp; Y++) {
			if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])]) {
				lrp = Y;
				break;

			}
			//else if (NodeOrder[Y] <= NumberOfSeqs) {
			//	for (Z = 0; Z <= mcp; Z++)
			//		NumDone[NodeOrder[Z]] = 1.0;

			//	//Pos = Y + 1;
			//	for (Pos = Y + 1; Pos < mcp; Pos++) {
			//		//mark the route
			//		if (NodeOrder[Pos] == NodeOrder[(int)(MidNode[0])]) {
			//			TallyDist = 0;
			//			//add up the branch lengths to the midnode
			//			for (Z = Y + 1; Z < Pos; Z++) {
			//				if (NumDone[NodeOrder[Z]] < 0) {
			//					TallyDist = TallyDist + NodeLen[NodeOrder[Z]];
			//				}
			//			}
			//			NodeLen[NodeOrder[Y]] = 10 - TallyDist;
			//			break;
			//		}
			//		else if (NodeOrder[Pos] > NumberOfSeqs)
			//			NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
			//	}
			//}

		}
		if (NodeOrder[(int)(MidNode[0])] <= NumberOfSeqs) {
			/*ln3 = 1;
			list3[1] = NodeOrder[(int)(MidNode[0])];*/
			rrp = lrp;
			/*for (Y = mcp; Y >= 0; Y--) {
				if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])]) {

					rrp = Y;
					lrp = Y;
					break;
				}
			}*/
			/*rrp = (int)MidNode[0];
			lrp = (int)MidNode[0];
			ln1 = 0;
			for (Y = 0; Y <= NumberOfSeqs; Y++) {
				if (list3[1] != Y) {
					ln1++;
					list1[ln1] = Y;
				}

			}*/
		}
		else {
			for (Y = mcp; Y >= 0; Y--) {
				if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])]) {

					rrp = Y;
					break;
				}
				//else if (NodeOrder[Y] <= NumberOfSeqs) {
				//	for (Z = 0; Z <= mcp; Z++)
				//		NumDone[NodeOrder[Z]] = 1.0;

				//	//Pos = Y + 1;
				//	for (Pos = Y - 1; Pos > 0; Pos--) {
				//		//mark the route
				//		if (NodeOrder[Pos] == NodeOrder[(int)(MidNode[0])]) {
				//			TallyDist = 0;
				//			//add up the branch lengths to the midnode
				//			for (Z = Y - 1; Z > Pos; Z--) {
				//				if (NumDone[NodeOrder[Z]] < 0) {
				//					TallyDist = TallyDist + NodeLen[NodeOrder[Z]];
				//				}
				//			}
				//			NodeLen[NodeOrder[Y]] = 10 - TallyDist;
				//			break;
				//		}
				//		else if (NodeOrder[Pos] > NumberOfSeqs)
				//			NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
				//	}
				//}

			}
			if (lrp == -1 || rrp == -1)
				lrp = 0;
			/*for (Y = 1; Y < lrp; Y++) {
				if (NodeOrder[Y] <= NumberOfSeqs) {
					ln1++;
					list1[ln1] = NodeOrder[Y];
				}
			}
			for (Y = mcp; Y >= rrp; Y--) {
				if (NodeOrder[Y] <= NumberOfSeqs) {
					ln2++;
					list2[ln2] = NodeOrder[Y];
				}
			}


			for (Y = lrp + 1; Y < rrp; Y++) {
				if (NodeOrder[Y] <= NumberOfSeqs) {
					ln3++;
					list3[ln3] = NodeOrder[Y];
				}
			}*/
		}
		//float maxd, maxdh;
		//maxd = (MidNode[3])*sf; //NumberOfSeqs * 2;



		//if (maxd > dim)
		//	maxd = dim;
		//if (maxd > sf / 2)
		//	maxd = sf / 2;
		//if (maxd > dim)
		//	maxd = dim;
		//tl[(int)(maxd)] = 1;
		//label sequence clusters branching off the root with the max distance
		
		int tp, nl1, nl2, st1, st2, en1, en2, mt, cpos, wnr, wnl, winn;
		
		cpos = -1;
		cpos++;
		no[cpos] = NumberOfSeqs * 2+3;
		for (A = lrp; A <= rrp; A++) {
			cpos++;
			no[cpos] = NodeOrder[A];
			////mask the "root" clade
			//mask[cpos] = 1;
		}
		
		cpos++;
		no[cpos] = NumberOfSeqs*2+2;
		for (A = lrp - 1; A > 0; A--) {
			cpos++;
			no[cpos] = NodeOrder[A];
		}
		for (A = mcp - 1; A >= rrp + 1; A--) {
			cpos++;
			no[cpos] = NodeOrder[A];
		}
		cpos++;
		no[cpos] = NumberOfSeqs*2+2;
		cpos++;
		no[cpos] = NumberOfSeqs * 2 +3;
		

		mcp = cpos;
		
		tp = lrp;
		
		mt = 0;
		
		int fp;
		float *mdist;
		float md;
		md = NumberOfSeqs + 1;
		mdist = (float*)calloc(NumberOfSeqs + 2, sizeof(float));
		if (rr == 0)
			fp = 0;
		else
			fp = 1;

		//for (mt = 0; mt <= 1; mt++) {
			for (Y = mcp; Y > 0; Y--) {
				//if (mask[Y] == mt) {
					st1 = -1;
					st2 = -1;
					en1 = -1;
					en2 = -1;
					if (no[Y] > NumberOfSeqs) {//on a node
						if (done[no[Y]] == 0) {
							done[no[Y]] = 1;
							nl1 = no[Y];
							st1 = Y;
							for (X = st1 - 1; X > 0; X--) {
								if (no[X] > NumberOfSeqs && done[no[X]] == 0){// && mask[X]==mt) {// on next most nested node
									nl2 = no[X];
									//done[no[X]] = 1;
									st2 = X;
									break;
								}
							}
							//scan till en2
							for (X = st2 - 1; X > 0; X--) {
								if (no[X] == nl2) {// at the end of the nested node
									en2 = X;
									break;
								}
							}
							if (en2 > -1) {//if en2 is -1 it means that an inner tip pair has been reached - these get sorted at the end.
								for (X = en2 - 1; X >= 0; X--) {
									if (no[X] == nl1) {// on end next outer node
										en1 = X;
										break;
									}
								}
								if (en1 > -1) {
									//add everthing from en1 to en2 and from st2 to st1 to list 1 and everything from en2 to st2 to list 2
									ln1 = 0;
									ln2 = 0;
									for (X = en1; X < en2; X++) {
										if (no[X] <= NumberOfSeqs){// && mask[X]==mt) {// on a sequence
											ln1++;
											list1[ln1] = no[X];
										}
									}
									for (X = st2 + 1; X < st1; X++) {
										if (no[X] <= NumberOfSeqs){// && mask[X] == mt) {// on a sequence
											ln1++;
											list1[ln1] = no[X];
										}
									}
									for (X = en2 + 1; X < st2; X++) {
										if (no[X] <= NumberOfSeqs){// && mask[X] == mt) {// on a sequence
											ln2++;
											list2[ln2] = no[X];
										}
									}
									if (ln1 > 0 && ln2 > 0) {
										//make average distance between list 1 and list 2 sequences
										//maxd = 0;
										//for (X = 1; X <= ln1; X++) {
										//	for (Z = 1; Z <= ln2; Z++) {
										//		maxd = maxd + TMat[list1[X] + list2[Z] * (NumberOfSeqs + 1)];
										//	}
										//}

										//maxd = maxd / (float)(ln1*ln2);
										//if (maxd >= 0.5)
										//	maxd = 0.5;
										//
										//maxd = maxd * sf;

										//if ((int)(maxd) > dim)
										//	maxd = dim;

										//maxdh = maxd;

										//if ((int)(maxd) >= 1) {

										//	if (tl[(int)maxd] != 0) {
										//		while (tl[(int)maxd] == 1 && (int)(maxd) >= 1) {
										//			maxd--;
										//			if ((int)(maxd) < 1)
										//				break;
										//		}

										//	}
										//}
										//if ((int)(maxd) <= 1) {

										//	maxd = maxdh;
										//	if ((int)(maxd) < 1)
										//		maxd = 1;
										//	//if (tl[(int)maxd] != 0) {
										//	while (tl[(int)(maxd)] == 1 && (int)(maxd) < dim) {
										//		maxd++;
										//		if ((int)(maxd) > dim)
										//			break;
										//	}

										//	//}


										//}
										//if ((int)(maxd) > dim || (int)(maxd) <1)
										//	maxd = dim;

										//tl[(int)maxd] = 1;
										if (fp==1) 
											md--;
										for (X = 1; X <= ln1; X++) {
											for (Z = 1; Z <= ln2; Z++) {
												//if(list1[X] == 8 || list1[X])
												TMat2[list1[X] + list2[Z] * (NumberOfSeqs + 1)] = md;//maxd;
												TMat2[list2[Z] + list1[X] * (NumberOfSeqs + 1)] = md;// maxd;
											}
										}
										
										fp = 1;
									}
								}
							}
						}
					}
				//}

			}
		//}

		//now finish off the terminal branch pairs
		//for (Y = 0; Y <= mcp; Y++) {
		//	if (no[Y] == NumberOfSeqs)
		//		no[Y] = NumberOfSeqs;
		//}
		
		for (Y = 0; Y < mcp; Y++) {
			/*if (Y == mcp - 10)
				Y == mcp - 10;*/
			if (no[Y] <= NumberOfSeqs && no[Y + 1] <= NumberOfSeqs) {
				
				//maxd = TMat[no[Y] + no[Y + 1] * (NumberOfSeqs + 1)];
				//if (maxd >= 0.5)
				//	maxd = 0.5;
				//maxd = maxd * sf;
				//if ((int)(maxd) > dim)
				//	maxd = dim;
				//maxdh = maxd;
				//if ((int)(maxd) > 1) {

				//	if (tl[(int)maxd] != 0) {
				//		while (tl[(int)maxd] == 1 && (int)(maxd) >= 1) {
				//			maxd--;
				//			if ((int)(maxd) < 1)
				//				break;
				//		}

				//	}
				//}
				//if ((int)(maxd) <= 1) {

				//	maxd = maxdh;
				//	if ((int)(maxd) < 1)
				//		maxd = 1;
				//	//if (tl[(int)maxd] != 0) {
				//	while (tl[(int)maxd] == 1 && (int)(maxd) < dim) {
				//		maxd++;
				//		if ((int)(maxd)> dim)
				//			break;
				//	}

				//	//}


				//}
				//if ((int)(maxd)> dim)
				//	maxd = dim;

				//if ((int)(maxd)<0)
				//	maxd = 1;

				//tl[(int)maxd] = 1;
				md--;
				TMat2[no[Y] + no[Y + 1] * (NumberOfSeqs + 1)] = md; //maxd;
				TMat2[no[Y + 1] + no[Y] * (NumberOfSeqs + 1)] = md;// maxd;

			}

		}
		//for (X = mcp; X >= mcp - 10; X--) {
		//	if (no[X] == 2000)
		//		no[X] = 1000;
		//}
		//for (X = 0; X <= NumberOfSeqs; X++) {
		//	for (Y =X+1; Y <= NumberOfSeqs; Y++) {
		//		if (TMat2[Y + X * (NumberOfSeqs + 1)] < 1) {
		//			TMat2[Y + X * (NumberOfSeqs + 1)] =0;
		//		}

		//	}

		//}
		//Z = 0;
		//if (sf < 1)
		//	sf = 1;

		//int cnt;
		//cnt = 0;
		for (Y = 0; Y <= NumberOfSeqs; Y++) {
			for (Z = Y + 1; Z <= NumberOfSeqs; Z++) {
				//if ((int)(TMat2[Z + Y * (NumberOfSeqs + 1)]) == X) {

					TMat[Y + Z * (NumberOfSeqs + 1)] = TMat2[Y + Z * (NumberOfSeqs + 1)] / 1000;
					TMat[Z + Y * (NumberOfSeqs + 1)] = TMat2[Y + Z * (NumberOfSeqs + 1)] / 1000;
				//}
			}

		}
		/*for (X = 0; X <= dim; X++) {
			if (tl[X] == 1) {
				cnt++;
				for (Y = 0; Y <= NumberOfSeqs; Y++) {
					for (Z = Y + 1; Z <= NumberOfSeqs; Z++) {
						if ((int)(TMat2[Z + Y * (NumberOfSeqs + 1)]) == X) {
							
							TMat[Y + Z * (NumberOfSeqs + 1)] = (float)cnt / 1000;
							TMat[Z + Y * (NumberOfSeqs + 1)] = (float)cnt / 1000;
						}
					}

				}

			}

		}*/

		
		//free(tl);
		free(done); 
		//free(mask);
		free(list1);
		free(list2);
		free(list3);
		free(no);
		
		

		return(1);
	}

	//int MyMathFuncs::UltraTreeDistP3(double MD, int MaxCurPos, int NumberOfSeqs, float *TMat2, double *NumDone, int *DoneThis, int *AbBe, int *NodeOrder, double *MidNode, double *NodeLen, float *TMat) {
	//	int X, Y, Z, C, A, Pos, lrp, rrp, dim, sf;
	//	//float TallyDist, Modi;
	//	//mark everything from beginning to midnode as 1
	//	lrp = -1;
	//	rrp = -1;
	//	int *list1, *list2, *list3, *done, *mask, *tl;
	//	int ln1, ln2, ln3, mcp;
	//	list1 = (int*)calloc(MaxCurPos + 1, sizeof(int));
	//	list2 = (int*)calloc(MaxCurPos + 1, sizeof(int));
	//	list3 = (int*)calloc(MaxCurPos + 1, sizeof(int));
	//	sf = 100000;
	//	dim = (int)((MidNode[3]) * sf);
	//	if (dim < 10) {
	//		dim = dim * 10000;
	//		sf = sf * 10000;
	//	}
	//	else if (dim < 100) {
	//		dim = dim * 1000;
	//		sf = sf * 1000;
	//	}
	//	else if (dim < 1000) {
	//		dim = dim * 100;
	//		sf = sf * 100;
	//	}
	//	else if (dim < 1000) {
	//		dim = dim * 10;
	//		sf = sf * 10;
	//	}
	//	dim = dim + 1;
	//	tl = (int*)calloc(dim + 1, sizeof(int));
	//	done = (int*)calloc(MaxCurPos, sizeof(int));
	//	mask = (int*)calloc(MaxCurPos, sizeof(int));


	//	ln1 = 0;
	//	ln2 = 0;
	//	ln3 = 0;
	//	mcp = MaxCurPos;


	//	for (Y = MaxCurPos; Y >= 0; Y--) {
	//		if (NodeOrder[Y] != 0) {
	//			mcp = Y;
	//			break;
	//		}
	//	}

	//	for (Y = 0; Y <= mcp; Y++) {
	//		if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])]) {
	//			lrp = Y;
	//			break;

	//		}
	//		//else if (NodeOrder[Y] <= NumberOfSeqs) {
	//		//	for (Z = 0; Z <= mcp; Z++)
	//		//		NumDone[NodeOrder[Z]] = 1.0;

	//		//	//Pos = Y + 1;
	//		//	for (Pos = Y + 1; Pos < mcp; Pos++) {
	//		//		//mark the route
	//		//		if (NodeOrder[Pos] == NodeOrder[(int)(MidNode[0])]) {
	//		//			TallyDist = 0;
	//		//			//add up the branch lengths to the midnode
	//		//			for (Z = Y + 1; Z < Pos; Z++) {
	//		//				if (NumDone[NodeOrder[Z]] < 0) {
	//		//					TallyDist = TallyDist + NodeLen[NodeOrder[Z]];
	//		//				}
	//		//			}
	//		//			NodeLen[NodeOrder[Y]] = 10 - TallyDist;
	//		//			break;
	//		//		}
	//		//		else if (NodeOrder[Pos] > NumberOfSeqs)
	//		//			NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
	//		//	}
	//		//}

	//	}
	//	if (NodeOrder[(int)(MidNode[0])] <= NumberOfSeqs) {
	//		ln3 = 1;
	//		list3[1] = NodeOrder[(int)(MidNode[0])];
	//		rrp = (int)MidNode[0];
	//		lrp = (int)MidNode[0];
	//		ln1 = 0;
	//		for (Y = 0; Y <= NumberOfSeqs; Y++) {
	//			if (list3[1] != Y) {
	//				ln1++;
	//				list1[ln1] = Y;
	//			}

	//		}
	//	}
	//	else {
	//		for (Y = mcp; Y >= 0; Y--) {
	//			if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])]) {

	//				rrp = Y;
	//				break;
	//			}
	//			//else if (NodeOrder[Y] <= NumberOfSeqs) {
	//			//	for (Z = 0; Z <= mcp; Z++)
	//			//		NumDone[NodeOrder[Z]] = 1.0;

	//			//	//Pos = Y + 1;
	//			//	for (Pos = Y - 1; Pos > 0; Pos--) {
	//			//		//mark the route
	//			//		if (NodeOrder[Pos] == NodeOrder[(int)(MidNode[0])]) {
	//			//			TallyDist = 0;
	//			//			//add up the branch lengths to the midnode
	//			//			for (Z = Y - 1; Z > Pos; Z--) {
	//			//				if (NumDone[NodeOrder[Z]] < 0) {
	//			//					TallyDist = TallyDist + NodeLen[NodeOrder[Z]];
	//			//				}
	//			//			}
	//			//			NodeLen[NodeOrder[Y]] = 10 - TallyDist;
	//			//			break;
	//			//		}
	//			//		else if (NodeOrder[Pos] > NumberOfSeqs)
	//			//			NumDone[NodeOrder[Pos]] = -NumDone[NodeOrder[Pos]];
	//			//	}
	//			//}

	//		}
	//		if (lrp == -1 || rrp == -1)
	//			lrp = 0;
	//		for (Y = 1; Y < lrp; Y++) {
	//			if (NodeOrder[Y] <= NumberOfSeqs) {
	//				ln1++;
	//				list1[ln1] = NodeOrder[Y];
	//			}
	//		}
	//		for (Y = mcp; Y >= rrp; Y--) {
	//			if (NodeOrder[Y] <= NumberOfSeqs) {
	//				ln2++;
	//				list2[ln2] = NodeOrder[Y];
	//			}
	//		}


	//		for (Y = lrp + 1; Y < rrp; Y++) {
	//			if (NodeOrder[Y] <= NumberOfSeqs) {
	//				ln3++;
	//				list3[ln3] = NodeOrder[Y];
	//			}
	//		}
	//	}
	//	float maxd, maxdh;
	//	maxd = (MidNode[3])*sf; //NumberOfSeqs * 2;


	//	if (maxd > dim)
	//		maxd = dim;
	//	if (maxd > sf / 2)
	//		maxd = sf / 2;
	//	if (maxd > dim)
	//		maxd = dim;
	//	tl[(int)(maxd)] = 1;
	//	//label sequence clusters branching off the root with the max distance
	//	for (Y = 1; Y <= ln3; Y++) {
	//		for (Z = 1; Z <= ln1; Z++) {
	//			TMat2[list3[Y] + list1[Z] * (NumberOfSeqs + 1)] = maxd;
	//			TMat2[list1[Z] + list3[Y] * (NumberOfSeqs + 1)] = maxd;
	//		}
	//	}
	//	for (Y = 1; Y <= ln3; Y++) {
	//		for (Z = 1; Z <= ln2; Z++) {
	//			TMat2[list3[Y] + list2[Z] * (NumberOfSeqs + 1)] = maxd;
	//			TMat2[list2[Z] + list3[Y] * (NumberOfSeqs + 1)] = maxd;
	//		}
	//	}


	//	//mask the "root" clade

	//	for (Y = lrp; Y <= rrp; Y++)
	//		mask[NodeOrder[Y]] = 1;

	//	////label sequences branching off the next branch with the next lowest distance
	//	//maxd = maxd - 1;
	//	//for (Y = 1; Y <= ln1; Y++) {
	//	//	for (Z = 1; Z <= ln2; Z++) {
	//	//		TMat2[list1[Y] + list2[Z] * (NumberOfSeqs + 1)] = maxd;
	//	//		TMat2[list2[Z] + list1[Y] * (NumberOfSeqs + 1)] = maxd;
	//	//	}
	//	//}

	//	//Now do relationships within each of the three subtrees separately
	//	int tp, nl1, nl2, st1, st2, en1, en2, mt;
	//	tp = lrp;
	//	for (mt = 0; mt <= 1; mt++) {
	//		for (Y = mcp; Y > 0; Y--) {
	//			if (mask[NodeOrder[Y]] == mt) {
	//				st1 = -1;
	//				st2 = -1;
	//				en1 = -1;
	//				en2 = -1;
	//				if (NodeOrder[Y] > NumberOfSeqs) {//on a node
	//					if (done[NodeOrder[Y]] == 0) {
	//						done[NodeOrder[Y]] == 1;
	//						nl1 = NodeOrder[Y];
	//						st1 = Y;
	//						for (X = st1 - 1; X > 0; X--) {
	//							if (NodeOrder[X] > NumberOfSeqs && done[NodeOrder[X]] == 0 && mask[NodeOrder[X]] == mt) {// on next most nested node
	//								nl2 = NodeOrder[X];
	//								//done[NodeOrder[X]] = 1;
	//								st2 = X;
	//								break;
	//							}
	//						}
	//						//scan till en2
	//						for (X = st2 - 1; X > 0; X--) {
	//							if (NodeOrder[X] == nl2) {// at the end of the nested node
	//								en2 = X;
	//								break;
	//							}
	//						}
	//						if (en2 > -1) {//if en2 is -1 it means that an inner tip pair has been reached - these get sorted at the end.
	//							for (X = en2 - 1; X >= 0; X--) {
	//								if (NodeOrder[X] == nl1) {// on end next outer node
	//									en1 = X;
	//									break;
	//								}
	//							}
	//							if (en1 > -1) {
	//								//add everthing from en1 to en2 and from st2 to st1 to list 1 and everything from en2 to st2 to list 2
	//								ln1 = 0;
	//								ln2 = 0;
	//								for (X = en1; X < en2; X++) {
	//									if (NodeOrder[X] <= NumberOfSeqs && mask[NodeOrder[X]] == mt) {// on a sequence
	//										ln1++;
	//										list1[ln1] = NodeOrder[X];
	//									}
	//								}
	//								for (X = st2 + 1; X < st1; X++) {
	//									if (NodeOrder[X] <= NumberOfSeqs && mask[NodeOrder[X]] == mt) {// on a sequence
	//										ln1++;
	//										list1[ln1] = NodeOrder[X];
	//									}
	//								}
	//								for (X = en2 + 1; X < st2; X++) {
	//									if (NodeOrder[X] <= NumberOfSeqs && mask[NodeOrder[X]] == mt) {// on a sequence
	//										ln2++;
	//										list2[ln2] = NodeOrder[X];
	//									}
	//								}
	//								if (ln1 > 0 && ln2 > 0) {
	//									//make average distance between list 1 and list 2 sequences
	//									maxd = 0;
	//									for (X = 1; X <= ln1; X++) {
	//										for (Z = 1; Z <= ln2; Z++) {
	//											maxd = maxd + TMat[list1[X] + list2[Z] * (NumberOfSeqs + 1)];
	//										}
	//									}

	//									maxd = maxd / (float)(ln1*ln2);
	//									if (maxd >= 0.5)
	//										maxd = 0.5;

	//									maxd = maxd * sf;

	//									if ((int)(maxd) > dim)
	//										maxd = dim;

	//									maxdh = maxd;

	//									if ((int)(maxd) >= 1) {

	//										if (tl[(int)maxd] != 0) {
	//											while (tl[(int)maxd] == 1 && (int)(maxd) >= 1) {
	//												maxd--;
	//												if ((int)(maxd) < 1)
	//													break;
	//											}

	//										}
	//									}
	//									if ((int)(maxd) <= 1) {

	//										maxd = maxdh;
	//										if ((int)(maxd) < 1)
	//											maxd = 1;
	//										//if (tl[(int)maxd] != 0) {
	//										while (tl[(int)(maxd)] == 1 && (int)(maxd) < dim) {
	//											maxd++;
	//											if ((int)(maxd) > dim)
	//												break;
	//										}

	//										//}


	//									}
	//									if ((int)(maxd) > dim || (int)(maxd) <1)
	//										maxd = dim;

	//									tl[(int)maxd] = 1;

	//									for (X = 1; X <= ln1; X++) {
	//										for (Z = 1; Z <= ln2; Z++) {
	//											//if(list1[X] == 8 || list1[X])
	//											TMat2[list1[X] + list2[Z] * (NumberOfSeqs + 1)] = maxd;
	//											TMat2[list2[Z] + list1[X] * (NumberOfSeqs + 1)] = maxd;
	//										}
	//									}
	//								}
	//							}
	//						}
	//					}
	//				}
	//			}

	//		}
	//	}

	//	//now finish off the terminal branch pairs
	//	for (Y = 0; Y <= mcp; Y++) {
	//		if (NodeOrder[Y] == NumberOfSeqs)
	//			NodeOrder[Y] = NumberOfSeqs;
	//	}
	//	for (Y = 0; Y < mcp; Y++) {
	//		/*if (Y == mcp - 10)
	//		Y == mcp - 10;*/
	//		if (NodeOrder[Y] <= NumberOfSeqs && NodeOrder[Y + 1] <= NumberOfSeqs) {

	//			maxd = TMat[NodeOrder[Y] + NodeOrder[Y + 1] * (NumberOfSeqs + 1)];
	//			if (maxd >= 0.5)
	//				maxd = 0.5;
	//			maxd = maxd * sf;
	//			if ((int)(maxd) > dim)
	//				maxd = dim;
	//			maxdh = maxd;
	//			if ((int)(maxd) > 1) {

	//				if (tl[(int)maxd] != 0) {
	//					while (tl[(int)maxd] == 1 && (int)(maxd) >= 1) {
	//						maxd--;
	//						if ((int)(maxd) < 1)
	//							break;
	//					}

	//				}
	//			}
	//			if ((int)(maxd) <= 1) {

	//				maxd = maxdh;
	//				if ((int)(maxd) < 1)
	//					maxd = 1;
	//				//if (tl[(int)maxd] != 0) {
	//				while (tl[(int)maxd] == 1 && (int)(maxd) < dim) {
	//					maxd++;
	//					if ((int)(maxd)> dim)
	//						break;
	//				}

	//				//}


	//			}
	//			if ((int)(maxd)> dim)
	//				maxd = dim;

	//			if ((int)(maxd)<0)
	//				maxd = 1;

	//			tl[(int)maxd] = 1;

	//			TMat2[NodeOrder[Y] + NodeOrder[Y + 1] * (NumberOfSeqs + 1)] = maxd;
	//			TMat2[NodeOrder[Y + 1] + NodeOrder[Y] * (NumberOfSeqs + 1)] = maxd;

	//		}

	//	}
	//	//for (X = mcp; X >= mcp - 10; X--) {
	//	//	if (NodeOrder[X] == 2000)
	//	//		NodeOrder[X] = 1000;
	//	//}
	//	//for (X = 0; X <= NumberOfSeqs; X++) {
	//	//	for (Y =X+1; Y <= NumberOfSeqs; Y++) {
	//	//		if (TMat2[Y + X * (NumberOfSeqs + 1)] < 1) {
	//	//			TMat2[Y + X * (NumberOfSeqs + 1)] =0;
	//	//		}

	//	//	}

	//	//}
	//	//Z = 0;
	//	if (sf < 1)
	//		sf = 1;

	//	int cnt;
	//	cnt = 0;
	//	for (X = 0; X <= dim; X++) {
	//		if (tl[X] == 1) {
	//			cnt++;
	//			for (Y = 0; Y <= NumberOfSeqs; Y++) {
	//				for (Z = Y + 1; Z <= NumberOfSeqs; Z++) {
	//					if ((int)(TMat2[Z + Y * (NumberOfSeqs + 1)]) == X) {

	//						TMat[Y + Z * (NumberOfSeqs + 1)] = (float)cnt / 1000;
	//						TMat[Z + Y * (NumberOfSeqs + 1)] = (float)cnt / 1000;
	//					}
	//				}

	//			}

	//		}

	//	}

	//	//for (X = 0; X <= NumberOfSeqs; X++) {
	//	//	for (Y = X+1; Y <= NumberOfSeqs; Y++) {
	//	//		if (TMat2[Y + X * (NumberOfSeqs + 1)] < sf) {
	//	//			TMat[Y + X * (NumberOfSeqs + 1)] = TMat2[Y + X * (NumberOfSeqs + 1)] / (float)(sf);
	//	//			TMat[X + Y * (NumberOfSeqs + 1)] = TMat[Y + X * (NumberOfSeqs + 1)];
	//	//		}
	//	//		else {
	//	//			TMat[Y + X * (NumberOfSeqs + 1)] = 0.99;
	//	//			TMat[X + Y * (NumberOfSeqs + 1)] = 0.99;
	//	//		}
	//	//		//Z++;
	//	//	}

	//	//}
	//	free(tl);
	//	free(done);
	//	free(mask);
	//	free(list1);
	//	free(list2);
	//	free(list3);
	//	return(1);
	//	//for (Y = 0; Y <= MaxCurPos; Y++) {
	//	//	if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])])
	//	//		break;
	//	//	AbBe[Y] = 1;
	//	//}

	//	////mark everything from end to midnode as 1
	//	//for (Y = MaxCurPos; Y >= 0; Y--) {
	//	//	if (Y == MaxCurPos) {
	//	//		while (NodeOrder[Y] == 0)
	//	//			Y--;

	//	//	}
	//	//	if (NodeOrder[Y] == NodeOrder[(int)(MidNode[0])])
	//	//		break;
	//	//	AbBe[Y] = 3;
	//	//}

	//	////find the last node that was added to the tree
	//	//for (Y = MaxCurPos; Y >= 0; Y--) {
	//	//	if (Y == MaxCurPos) {
	//	//		while (NodeOrder[Y] == 0)
	//	//			Y--;
	//	//	}
	//	//	if (AbBe[Y] == 0)
	//	//		AbBe[Y] = 2;

	//	//}

	//	////now modify the tree distance matrix to reflect equal distances from the root
	//	////In effect its lengthening the terminal branches so that the distance between
	//	////sequences in a distance matrix will reflect their reletive positions in the tree
	//	////in relation to the midpoint "root"
	//	////first do the "left" part of the tree
	//	//unsigned char *rn;
	//	//int tc, bc;
	//	//tc = 0;
	//	//bc = 0;
	//	//rn = (unsigned char*)calloc((NumberOfSeqs + 1)*(NumberOfSeqs + 1), sizeof(unsigned char));
	//	//for (Y = 0; Y <= MaxCurPos; Y++) {
	//	//	if (AbBe[Y] > 0) {
	//	//		if (NodeOrder[Y] <= NumberOfSeqs) {
	//	//			if (DoneThis[NodeOrder[Y]] == 0) {
	//	//				DoneThis[NodeOrder[Y]] = 1;
	//	//				TallyDist = 0.0;
	//	//				for (Z = 0; Z <= MaxCurPos; Z++)
	//	//					NumDone[Z] = 1.0;

	//	//				if (NodeOrder[Y] != NodeOrder[(int)(MidNode[0])]) {
	//	//					TallyDist = TallyDist + NodeLen[NodeOrder[Y]];
	//	//					NumDone[NodeOrder[Y]] = -NumDone[NodeOrder[Y]];
	//	//					if (AbBe[Y] == 2 || AbBe[Y] == 1) {
	//	//						for (Z = Y + 1; Z <= MaxCurPos; Z++) {
	//	//							if (NodeOrder[Z] == NodeOrder[(int)(MidNode[0])]) {
	//	//								TallyDist = TallyDist + MidNode[AbBe[Y]];
	//	//								Modi = MD - TallyDist;

	//	//								C = NodeOrder[Y];
	//	//								for (A = 0; A <= NumberOfSeqs; A++) {
	//	//									if (A != C) {
	//	//										//if (TMat2[A + C*(NumberOfSeqs + 1)] < 1) {
	//	//											//TMat2[A + C*(NumberOfSeqs + 1)] = TMat2[A + C*(NumberOfSeqs + 1)] / 100000;
	//	//										rn[A + C*(NumberOfSeqs + 1)] = rn[A + C*(NumberOfSeqs + 1)] + 1;
	//	//										rn[C + A*(NumberOfSeqs + 1)] = rn[A + C*(NumberOfSeqs + 1)];
	//	//											TMat2[A + C*(NumberOfSeqs + 1)] = TMat2[A + C*(NumberOfSeqs + 1)] + (float)(Modi);
	//	//											//TMat2[A + C*(NumberOfSeqs + 1)] = (int)(TMat2[A + C*(NumberOfSeqs + 1)] * 100000);
	//	//											TMat2[C + A*(NumberOfSeqs + 1)] = TMat2[A + C*(NumberOfSeqs + 1)];
	//	//											tc++;
	//	//										//}
	//	//									}
	//	//								}
	//	//								break;
	//	//							}
	//	//							else if (NodeOrder[Z] > NumberOfSeqs) {
	//	//								TallyDist = TallyDist + NodeLen[NodeOrder[Z]] * NumDone[NodeOrder[Z]];
	//	//								NumDone[NodeOrder[Z]] = -NumDone[NodeOrder[Z]];
	//	//							}
	//	//						}
	//	//					}
	//	//					else {
	//	//						for (Z = Y - 1; Z >= 0; Z--) {

	//	//							if (NodeOrder[Z] == NodeOrder[(int)(MidNode[0])]) {
	//	//								TallyDist = TallyDist + MidNode[1];
	//	//								Modi = MD - TallyDist;
	//	//								C = NodeOrder[Y];
	//	//								for (A = 0; A <= NumberOfSeqs; A++) {
	//	//									if (A != C) {
	//	//										//if (TMat2[A + C*(NumberOfSeqs + 1)] < 1) {
	//	//											TMat2[A + C*(NumberOfSeqs + 1)] = TMat2[A + C*(NumberOfSeqs + 1)] + (float)(Modi);
	//	//											//TMat2[A + C*(NumberOfSeqs + 1)] = (int)(TMat2[A + C*(NumberOfSeqs + 1)] * 100000);
	//	//											TMat2[C + A*(NumberOfSeqs + 1)] = TMat2[A + C*(NumberOfSeqs + 1)];
	//	//											rn[A + C*(NumberOfSeqs + 1)] = rn[A + C*(NumberOfSeqs + 1)] +4;
	//	//											rn[C + A*(NumberOfSeqs + 1)] = rn[A + C*(NumberOfSeqs + 1)];
	//	//											bc++;
	//	//										//}
	//	//									}
	//	//								}
	//	//								break;
	//	//							}
	//	//							else if (NodeOrder[Z] > NumberOfSeqs) {
	//	//								TallyDist = TallyDist + NodeLen[NodeOrder[Z]] * NumDone[NodeOrder[Z]];
	//	//								NumDone[NodeOrder[Z]] = -NumDone[NodeOrder[Z]];
	//	//							}
	//	//						}
	//	//					}
	//	//				}
	//	//			}
	//	//		}
	//	//	}
	//	//	else if (AbBe[Y] == 0)
	//	//		break;

	//	//}

	//	//
	//	////Need to resolve bad TMat2 "ties" here
	//	//short int *dp;
	//	//int cs, yy, x, z,y,zn;
	//	//dp = (short int*)calloc((NumberOfSeqs + 1)*(NumberOfSeqs + 1), sizeof(short int));
	//	//
	//	//
	//	//
	//	//

	//	///*zn = 0;
	//	//for (Z = 0; Z <= MaxCurPos; Z++) {
	//	//	if (NodeLen[Z] > 0)
	//	//		zn++;

	//	//}*/
	//	////detect pairs on opposite sides of the midpoint
	//	//bc = 0;
	//	//for (x = 0; x < NumberOfSeqs; x++) {
	//	//	for (z = x+1; z <= NumberOfSeqs; z++) {
	//	//		TMat2[x + z*(NumberOfSeqs + 1)] = (int)(TMat2[x + z*(NumberOfSeqs + 1)] * 10000);
	//	//		TMat2[x + z*(NumberOfSeqs + 1)] = TMat2[x + z*(NumberOfSeqs + 1)] / 10000;
	//	//		if (rn[x + z*(NumberOfSeqs + 1)] != 2 && rn[x + z*(NumberOfSeqs + 1)] != 8 && rn[x + z*(NumberOfSeqs + 1)] != 0) {
	//	//			bc++;
	//	//			dp[x + z*(NumberOfSeqs + 1)] = 1;
	//	//			dp[z + x*(NumberOfSeqs + 1)] = 1;
	//	//			TMat2[x + z*(NumberOfSeqs + 1)] = 10;
	//	//			TMat2[z + x*(NumberOfSeqs + 1)] = 10;

	//	//		}

	//	//		

	//	//	}
	//	//}
	//	//cs = 0;
	//	//for (x = 0; x < NumberOfSeqs; x++) {
	//	//	for (z = x + 1; z <= NumberOfSeqs; z++) {
	//	//		if (dp[z + x*(NumberOfSeqs + 1)] == 0 && TMat2[z + x*(NumberOfSeqs + 1)] < 5) {
	//	//			cs++;
	//	//			dp[x + z*(NumberOfSeqs + 1)] = cs;
	//	//			dp[z + x*(NumberOfSeqs + 1)] = cs;
	//	//			for (y =0; y <= NumberOfSeqs; y++) {
	//	//				if (y != x && y != z) {
	//	//					if ((TMat2[x + z*(NumberOfSeqs + 1)]) == (TMat2[x + y*(NumberOfSeqs + 1)])) {
	//	//						dp[x + y*(NumberOfSeqs + 1)] = cs;
	//	//						dp[y + x*(NumberOfSeqs + 1)] = cs;
	//	//						for (yy = y+1; yy <= NumberOfSeqs; yy++) {
	//	//							if (yy != z && yy != x && yy != y) {
	//	//								if ((TMat2[x + y*(NumberOfSeqs + 1)]) == (TMat2[yy + y*(NumberOfSeqs + 1)])) {
	//	//									dp[y + yy*(NumberOfSeqs + 1)] = cs;
	//	//									dp[yy + y*(NumberOfSeqs + 1)] = cs;
	//	//									dp[x+ yy*(NumberOfSeqs + 1)] = cs;
	//	//									dp[yy + x*(NumberOfSeqs + 1)] = cs;
	//	//								}
	//	//							}
	//	//						}
	//	//					}
	//	//					else if ((TMat2[x + z*(NumberOfSeqs + 1)]) == (TMat2[z + y*(NumberOfSeqs + 1)])) {
	//	//						dp[y + z*(NumberOfSeqs + 1)] = cs;
	//	//						dp[z + y*(NumberOfSeqs + 1)] = cs;
	//	//						for (yy = y+1; yy <= NumberOfSeqs; yy++) {
	//	//							if (yy != z && yy != x && yy != y) {
	//	//								if ((TMat2[z + y*(NumberOfSeqs + 1)]) == (TMat2[yy + y*(NumberOfSeqs + 1)])) {
	//	//									dp[y + yy*(NumberOfSeqs + 1)] = cs;
	//	//									dp[yy + y*(NumberOfSeqs + 1)] = cs;
	//	//									dp[z + yy*(NumberOfSeqs + 1)] = cs;
	//	//									dp[yy + z*(NumberOfSeqs + 1)] = cs;
	//	//								}
	//	//							}
	//	//						}
	//	//					}
	//	//				}

	//	//			}
	//	//		}

	//	//	}
	//	//}



	//	//free(dp);
	//	//free(rn);

	//	return(1);
	//}
	int MyMathFuncs::Tree2ArrayP(unsigned char EarlyExitFlag, int NameLen, int NumberOfSeqs, int LTree, char *T2Holder, int UBTM2, float *TMat2) {
		/*Dim CS As Long, XYP(1) As Long, LowD As Double, TMat2Bak() As Single
		Dim MaxCurPos As Long, tMat() As Single, TB(1) As Long, MidNode(2) As Double
		Dim DoneNode() As Long, NodeOrder() As Long, RootNode() As Byte, AbBe() As Long, DoneThis() As Long, TempNodeOrder() As Long, NodeLen() As Double, NumDone() As Double
		Dim Z As Long, C As Long, Inside As Byte, MD As Double, TallyDist As Double, Done0 As Byte
		Dim Y As Long, x As Long, Dummy As Long, NoS As Long*/
		int MaxCurPos, NoS, Dummy, x, Y;
		double MD;
		double *NodeLen, *NumDone, *MidNode;
		float *tMat, *TMat2Bak;
		int *TempNodeOrder, *NodeOrder, *DoneNode, *DoneThis, *TB, *AbBe;
		unsigned char *RootNode;
		MaxCurPos = NumberOfSeqs * 3;


		NodeLen = (double*)calloc(MaxCurPos + 1, sizeof(double));
		NumDone = (double*)calloc(MaxCurPos + 1, sizeof(double));
		MidNode = (double*)calloc(4, sizeof(double));
		TempNodeOrder = (int*)calloc(MaxCurPos + 1, sizeof(int));
		AbBe = (int*)calloc(MaxCurPos + 1, sizeof(int));
		TB = (int*)calloc(3, sizeof(int));
		NodeOrder = (int*)calloc(MaxCurPos + 1, sizeof(int));
		DoneNode = (int*)calloc(MaxCurPos + 1, sizeof(int));
		RootNode = (unsigned char*)calloc((MaxCurPos + 1)*(MaxCurPos + 1), sizeof(unsigned char));
		tMat = (float*)calloc((MaxCurPos + 1)*(MaxCurPos + 1), sizeof(float));
		DoneThis = (int*)calloc(MaxCurPos + 1, sizeof(int));
		TMat2Bak = (float*)calloc((NumberOfSeqs + 1)*(NumberOfSeqs + 1), sizeof(float));

		NoS = NumberOfSeqs;

		//Make a matrix of distances within the tree

		Dummy = TreeToArrayP(NameLen, NumberOfSeqs, LTree, T2Holder, TMat2, NodeOrder, DoneNode, TempNodeOrder, RootNode, NodeLen, NumDone);

		for (x = 0; x <= NumberOfSeqs; x++) {
			if (TMat2[x + x*(UBTM2 + 1)] != 0)
				TMat2[x + x*(UBTM2 + 1)] = 0;

		}
		for (x = 0; x <= MaxCurPos; x++)
			NumDone[x] = 0;


		for (x = 0; x <= NumberOfSeqs; x++) {
			for (Y = x + 1; Y <= NumberOfSeqs; Y++) {
				TMat2[x + Y*(UBTM2 + 1)] = round(TMat2[x + Y*(UBTM2 + 1)] * 100000) / 100000;
				TMat2[Y + x*(UBTM2 + 1)] = TMat2[x + Y*(UBTM2 + 1)];
			}
		}

		for (x = 0; x <= MaxCurPos; x++)
			NodeLen[x] = round(NodeLen[x] * 100000) / 100000;



		MD = TreeMidP(MaxCurPos, NumberOfSeqs, NumDone, TMat2, TB, NodeOrder, MidNode, NodeLen);


		for (x = 0; x <= MaxCurPos; x++)
			NodeLen[x] = round(NodeLen[x] * 100000) / 100000;

		MidNode[1] = round(MidNode[1] * 100000) / 100000;
		MidNode[2] = round(MidNode[2] * 100000) / 100000;


		Dummy = UltraTreeDistP(MD, MaxCurPos, NoS, TMat2, NumDone, DoneThis, AbBe, NodeOrder, MidNode, NodeLen);


		for (x = 0; x <= NumberOfSeqs; x++) {
			for (Y = x + 1; Y <= NumberOfSeqs; Y++) {
				TMat2[x + Y*(UBTM2 + 1)] = round(TMat2[x + Y*(UBTM2 + 1)] * 100000) / 100000;
				TMat2[Y + x*(UBTM2 + 1)] = TMat2[x + Y*(UBTM2 + 1)];
			}
		}


		Dummy = MakeTreeArrayXP2(NumberOfSeqs, TMat2, TMat2Bak);



		free(NodeLen);
		free(NumDone);
		free(MidNode);
		free(TempNodeOrder);
		free(AbBe);
		free(TB);
		free(NodeOrder);
		free(DoneNode);
		free(RootNode);
		free(tMat);
		free(DoneThis);
		free(TMat2Bak);

		return(1);
	}

	int MyMathFuncs::Tree2ArrayP2(int rr, int NameLen, int NumberOfSeqs, int LTree, char *T2Holder, int UBTM2, float *TMat2) {
		/*Dim CS As Long, XYP(1) As Long, LowD As Double, TMat2Bak() As Single
		Dim MaxCurPos As Long, tMat() As Single, TB(1) As Long, MidNode(2) As Double
		Dim DoneNode() As Long, NodeOrder() As Long,  AbBe() As Long, DoneThis() As Long, TempNodeOrder() As Long, NodeLen() As Double, NumDone() As Double
		Dim Z As Long, C As Long, Inside As Byte, MD As Double, TallyDist As Double, Done0 As Byte
		Dim Y As Long, x As Long, Dummy As Long, NoS As Long*/
		int MaxCurPos, NoS, Dummy, x, Y;
		double MD;
		double *NodeLen, *NumDone, *MidNode;
		//float *tMat, 
		float	*TMat2Bak;
		int *TempNodeOrder, *NodeOrder, *DoneNode, *DoneThis, *TB, *AbBe;
		
		MaxCurPos = NumberOfSeqs * 3 + 100;
		

		NodeLen = (double*)calloc(MaxCurPos+1, sizeof(double));
		NumDone = (double*)calloc(MaxCurPos +1, sizeof(double));
		MidNode = (double*)calloc(4, sizeof(double));
		TempNodeOrder = (int*)calloc(MaxCurPos + 1, sizeof(int));
		AbBe = (int*)calloc(MaxCurPos + 1, sizeof(int));
		TB = (int*)calloc(3, sizeof(int));
		NodeOrder = (int*)calloc(MaxCurPos + 1, sizeof(int));
		DoneNode = (int*)calloc(MaxCurPos + 1, sizeof(int));
		
		//tMat = (float*)calloc((MaxCurPos + 1)*(MaxCurPos + 1), sizeof(float));
		DoneThis = (int*)calloc(MaxCurPos + 1, sizeof(int));
		TMat2Bak = (float*)calloc((NumberOfSeqs + 1)*(NumberOfSeqs + 1), sizeof(float));

		NoS = NumberOfSeqs;

		//Make a matrix of distances within the tree

		Dummy = TreeToArrayP2(NameLen, NumberOfSeqs, LTree, T2Holder, TMat2, NodeOrder, DoneNode, TempNodeOrder, NodeLen, NumDone);

		for (x = 0; x<= NumberOfSeqs; x++){
			if (TMat2[x + x*(UBTM2 + 1)] != 0)
				TMat2[x + x*(UBTM2 + 1)] = 0;
		   
		}
		for (x = 0; x <= MaxCurPos; x++)
			NumDone[x] = 0;
		

		for (x = 0; x <= NumberOfSeqs; x++) {
			for (Y = x + 1; Y <= NumberOfSeqs; Y++){
				TMat2[x + Y*(UBTM2 + 1)] = round(TMat2[x + Y*(UBTM2 + 1)] * 100000) / 100000;
				TMat2[Y + x*(UBTM2 + 1)] = TMat2[x + Y*(UBTM2 + 1)];
			}
		}

		for (x = 0; x <= MaxCurPos; x++)
			NodeLen[x] = round(NodeLen[x] * 100000) / 100000;
		

		
		MD = TreeMidP(MaxCurPos, NumberOfSeqs, NumDone, TMat2, TB, NodeOrder, MidNode, NodeLen);

		
		for (x = 0; x <= MaxCurPos; x++)
			NodeLen[x] = round(NodeLen[x] * 100000) / 100000;

		MidNode[1] = round(MidNode[1] * 100000) / 100000;
		MidNode[2] = round(MidNode[2] * 100000) / 100000;


		Dummy = UltraTreeDistP2(rr, MaxCurPos, NoS, TMat2Bak, NumDone, DoneThis, AbBe, NodeOrder, MidNode, NodeLen, TMat2);
		//Dummy = UltraTreeDistP(MD, MaxCurPos, NoS, TMat2, NumDone, DoneThis, AbBe, NodeOrder, MidNode, NodeLen );

		/*for (x = 0; x <= NumberOfSeqs; x++) {
			for (Y = x + 1; Y <= NumberOfSeqs; Y++) {
				TMat2[x + Y*(UBTM2 + 1)] = round(TMat2[x + Y*(UBTM2 + 1)] * 100000) / 100000;
				TMat2[Y + x*(UBTM2 + 1)] = TMat2[x + Y*(UBTM2 + 1)];
			}
		}
*/
		
		//Dummy = MakeTreeArrayXP2(NumberOfSeqs, TMat2, TMat2Bak);

		

		free(NodeLen);
		free(NumDone);
		free(MidNode);
		free(TempNodeOrder);
		free(AbBe);
		free(TB);
		free(NodeOrder);
		free(DoneNode);
		
		//free(tMat);
		free(DoneThis);
		free(TMat2Bak);

		return(1);
	}
	
	double MyMathFuncs::SuperDistP(int X, int Nextno, int UB14, int UB04, int UB13,int UB03, int UB12, int UB02, int UB11, double *avdst, float *pd, float *pv, float *dist, short int *redodist, int *SeqCatCount, short int *ISeq14, short int *ISeq04, short int *ISeq13, short int *ISeq03, short int *ISeq12, short int *ISeq02, short int *ISeq11, char *CompressValid14, char *CompressDiffs14, char *CompressValid13, char *CompressDiffs13, char *CompressValid12, char *CompressDiffs12, char *CompressValid11, char *CompressDiffs11, char *CompressDiffs04, char *CompressDiffs03, char *CompressDiffs02) {
		int Z, Y;
		int TValid, TDiffs, os14, os04, os13, os03, os12, os02, os11, ps14, ps04, ps13, ps03, ps12, ps02, ps11, qs14, qs04, qs13, qs03, qs12, qs02, qs11;
		int cv14, cv04, cv13, cv03, cv12, cv02, cv11;


		double v1, d1, dX, ad, upper;
		int o3, o4;
		ad = *avdst;
		upper = 0.0;

		cv14 = 626;
		cv04 = 1025;
		cv13 = 1025;
		cv03 = 730;
		cv12 = 730;
		cv02 = 1025;
		cv11 = 1025;

		os14 = UB14 + 1;
		os04 = UB04 + 1;
		os13 = UB13 + 1;
		os03 = UB03 + 1;
		os12 = UB12 + 1;
		os02 = UB02 + 1;
		os11 = UB11 + 1;

		qs14 = os14*X;
		qs04 = os04*X;
		qs13 = os13*X;
		qs03 = os03*X;
		qs12 = os12*X;
		qs02 = os02*X;
		qs11 = os11*X;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 -1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel for private (Z, TValid, TDiffs, ps14,ps04,ps13,ps03,ps12,ps02,ps11,v1,d1, o3, o4, dX)
		for (Y = X + 1; Y <= Nextno; Y++) {
			if (redodist[X] + redodist[Y] > 0) {

				TValid = 0;
				TDiffs = 0;

				ps14 = os14*Y;
				ps04 = os04*Y;
				ps13 = os13*Y;
				ps03 = os03*Y;
				ps12 = os12*Y;
				ps02 = os02*Y;
				ps11 = os11*Y;
				v1 = 0;
				d1 = 0;
				//#pragma omp parallel{
								//14
				//#pragma omp parallel sections
				//				{
				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))		
				for (Z = 1; Z <= UB14; Z++) {
					TValid = TValid + CompressValid14[ISeq14[Z + qs14] + ISeq14[Z + ps14] * cv14];
					TDiffs = TDiffs + CompressDiffs14[ISeq14[Z + qs14] + ISeq14[Z + ps14] * cv14];
				}
				//					}
								//13

				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))
				for (Z = 1; Z <= UB13; Z++) {
					TValid = TValid + CompressValid13[ISeq13[Z + qs13] + ISeq13[Z + ps13] * cv13];
					TDiffs = TDiffs + CompressDiffs13[ISeq13[Z + qs13] + ISeq13[Z + ps13] * cv13];
				}
				//					}
								//04
				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))
				for (Z = 1; Z <= UB04; Z++) {
					TDiffs = TDiffs + CompressDiffs04[ISeq04[Z + qs04] + ISeq04[Z + ps04] * cv04];
				}

				TValid = TValid + SeqCatCount[8];
				//					}
								//12
				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))
				for (Z = 1; Z <= UB12; Z++) {

					TValid = TValid + CompressValid12[ISeq12[Z + qs12] + ISeq12[Z + ps12] * cv12];
					TDiffs = TDiffs + CompressDiffs12[ISeq12[Z + qs12] + ISeq12[Z + ps12] * cv12];
				}
				//					}
								//03
				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))
				for (Z = 1; Z <= UB03; Z++) {
					TDiffs = TDiffs + CompressDiffs03[ISeq03[Z + qs03] + ISeq03[Z + ps03] * cv03];
				}
				TValid = TValid + SeqCatCount[6];
				//					}
								//11
				//#pragma omp section
				//#pragma loop(hint_parallel(8))
				//					{
				for (Z = 1; Z <= UB11; Z++) {
					TValid = TValid + CompressValid11[ISeq11[Z + qs11] + ISeq11[Z + ps11] * cv11];
					TDiffs = TDiffs + CompressDiffs11[ISeq11[Z + qs11] + ISeq11[Z + ps11] * cv11];
				}
				//					}
								//02
				//#pragma omp section
				//#pragma loop(hint_parallel(8))
				//					{
				for (Z = 1; Z <= UB02; Z++)
					TDiffs = TDiffs + CompressDiffs02[ISeq02[Z + qs02] + ISeq02[Z + ps02] * cv02];

				TValid = TValid + SeqCatCount[4];
				TValid = TValid + SeqCatCount[2];
				//					}
				//				}
				//				}


				if (TValid >= 1) {
					v1 = (double)(TValid);
					d1 = (double)(TDiffs);
					if (v1 >= 1.0)
						dX = (double)((v1 - d1) / v1);
					else
						dX = 0.0;
				}
				else
					dX = 0.0;

				o3 = X + Y*(Nextno + 1);
				o4 = Y + X*(Nextno + 1);
				dist[o3] = (float)(dX);
				dist[o4] = (float)(dX);
				pv[o3] = (float)(v1);
				pv[o4] = (float)(v1);
				pd[o3] = (float)(d1);
				pd[o4] = (float)(d1);
#pragma omp critical
				{
					ad += (1.0 - dX);
					if (dX > upper)
						upper = dX;
					
				}
			}

		}
		*avdst = ad;
		omp_set_num_threads(2);
		return(upper);
	}

	double MyMathFuncs::SuperDistP2(int XX, int Nextno, int UB14, int UB04, int UB13, int UB03, int UB12, int UB02, int UB11, double *avdst, float *pd, float *pv, float *dist, short int *redodist, int *SeqCatCount, short int *ISeq14, short int *ISeq04, short int *ISeq13, short int *ISeq03, short int *ISeq12, short int *ISeq02, short int *ISeq11, char *CompressValid14, char *CompressDiffs14, char *CompressValid13, char *CompressDiffs13, char *CompressValid12, char *CompressDiffs12, char *CompressValid11, char *CompressDiffs11, char *CompressDiffs04, char *CompressDiffs03, char *CompressDiffs02) {
		int X,Z, Y;

		int procs;
		procs = omp_get_num_procs();
		procs = procs / 2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

		int *xy;
		int c;
		xy = (int *)calloc((Nextno + 1)*Nextno, sizeof(int));
		c = 0;
		for (X = 0; X < Nextno; X++) {
			for (Y = X + 1; Y <= Nextno; Y++) {
				if (redodist[X] + redodist[Y] > 0) {
					xy[c] = X;
					xy[c + 1] = Y;
					c = c + 2;
				}
			}
		}
		c = c - 2;
		c = c / 2;
		double ad, upper;
		ad = *avdst;
		upper = 0.0;
#pragma omp parallel
		{
			int TValid, TDiffs, os14, os04, os13, os03, os12, os02, os11, ps14, ps04, ps13, ps03, ps12, ps02, ps11, qs14, qs04, qs13, qs03, qs12, qs02, qs11;
			int cv14, cv04, cv13, cv03, cv12, cv02, cv11;


			double v1, d1, dX;
			int o3, o4;
			

			cv14 = 626;
			cv04 = 1025;
			cv13 = 1025;
			cv03 = 730;
			cv12 = 730;
			cv02 = 1025;
			cv11 = 1025;

			os14 = UB14 + 1;
			os04 = UB04 + 1;
			os13 = UB13 + 1;
			os03 = UB03 + 1;
			os12 = UB12 + 1;
			os02 = UB02 + 1;
			os11 = UB11 + 1;

			

			int d;
#pragma omp for private (X,Y,d,Z, TValid, TDiffs, ps14,ps04,ps13,ps03,ps12,ps02,ps11,v1,d1, o3, o4, dX)
			for (d = 0; d <= c; d++) {
				X = xy[d * 2];
				Y = xy[1 + d * 2];

				qs14 = os14*X;
				qs04 = os04*X;
				qs13 = os13*X;
				qs03 = os03*X;
				qs12 = os12*X;
				qs02 = os02*X;
				qs11 = os11*X;

				//for (Y = X + 1; Y <= Nextno; Y++) {
			//	if (redodist[X] + redodist[Y] > 0) {

				TValid = 0;
				TDiffs = 0;

				ps14 = os14*Y;
				ps04 = os04*Y;
				ps13 = os13*Y;
				ps03 = os03*Y;
				ps12 = os12*Y;
				ps02 = os02*Y;
				ps11 = os11*Y;
				v1 = 0;
				d1 = 0;
				//#pragma omp parallel{
				//14
				//#pragma omp parallel sections
				//				{
				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))		
				for (Z = 1; Z <= UB14; Z++) {
					TValid = TValid + CompressValid14[ISeq14[Z + qs14] + ISeq14[Z + ps14] * cv14];
					TDiffs = TDiffs + CompressDiffs14[ISeq14[Z + qs14] + ISeq14[Z + ps14] * cv14];
				}
				//					}
				//13

				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))
				for (Z = 1; Z <= UB13; Z++) {
					TValid = TValid + CompressValid13[ISeq13[Z + qs13] + ISeq13[Z + ps13] * cv13];
					TDiffs = TDiffs + CompressDiffs13[ISeq13[Z + qs13] + ISeq13[Z + ps13] * cv13];
				}
				//					}
				//04
				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))
				for (Z = 1; Z <= UB04; Z++) {
					TDiffs = TDiffs + CompressDiffs04[ISeq04[Z + qs04] + ISeq04[Z + ps04] * cv04];
				}

				TValid = TValid + SeqCatCount[8];
				//					}
				//12
				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))
				for (Z = 1; Z <= UB12; Z++) {

					TValid = TValid + CompressValid12[ISeq12[Z + qs12] + ISeq12[Z + ps12] * cv12];
					TDiffs = TDiffs + CompressDiffs12[ISeq12[Z + qs12] + ISeq12[Z + ps12] * cv12];
				}
				//					}
				//03
				//#pragma omp section
				//					{
				//#pragma loop(hint_parallel(8))
				for (Z = 1; Z <= UB03; Z++) {
					TDiffs = TDiffs + CompressDiffs03[ISeq03[Z + qs03] + ISeq03[Z + ps03] * cv03];
				}
				TValid = TValid + SeqCatCount[6];
				//					}
				//11
				//#pragma omp section
				//#pragma loop(hint_parallel(8))
				//					{
				for (Z = 1; Z <= UB11; Z++) {
					TValid = TValid + CompressValid11[ISeq11[Z + qs11] + ISeq11[Z + ps11] * cv11];
					TDiffs = TDiffs + CompressDiffs11[ISeq11[Z + qs11] + ISeq11[Z + ps11] * cv11];
				}
				//					}
				//02
				//#pragma omp section
				//#pragma loop(hint_parallel(8))
				//					{
				for (Z = 1; Z <= UB02; Z++)
					TDiffs = TDiffs + CompressDiffs02[ISeq02[Z + qs02] + ISeq02[Z + ps02] * cv02];

				TValid = TValid + SeqCatCount[4];
				TValid = TValid + SeqCatCount[2];
				//					}
				//				}
				//				}


				if (TValid >= 1) {
					v1 = (double)(TValid);
					d1 = (double)(TDiffs);
					if (v1 >= 1.0)
						dX = (double)((v1 - d1) / v1);
					else
						dX = 0.0;
				}
				else
					dX = 0.0;

				o3 = X + Y*(Nextno + 1);
				o4 = Y + X*(Nextno + 1);
				dist[o3] = (float)(dX);
				dist[o4] = (float)(dX);
				pv[o3] = (float)(v1);
				pv[o4] = (float)(v1);
				pd[o3] = (float)(d1);
				pd[o4] = (float)(d1);
//#pragma omp critical
//				{
//					ad += (1.0 - dX);
//					if (dX > upper)
//						upper = dX;
//
//				}
				//}

			}
		}
		*avdst = ad;
		free(xy);
		omp_set_num_threads(2);
		return(upper);
	}

	double NormalZ(double Z) {

		long double Y, X, w, temp, Temp2;
		double Z_MAX, NormalZx, WinP;
		NormalZx = 0.0;
		Z_MAX = 6;

		if (fabs(Z) < 5.9999999) {
			if (Z == 0.0)
				X = 0.0;
			else {
				Y = 0.5 * fabs(Z);
				if (Y >= (Z_MAX * 0.5))
					X = 1.0;
				else if (Y < 1.0) {

					w = Y * Y;
					X = ((((((((0.000124818987 * w - 0.001075204047) * w + 0.005198775019) * w - 0.019198292004) * w + 0.059054035642) * w - 0.151968751364) * w + 0.319152932694) * w - 0.5319230073) * w + 0.797884560593) * Y * 2.0;

				}
				else {

					Y = Y - 2.0;
					X = (((((((((((((-0.000045255659 * Y
						+ 0.00015252929) * Y - 0.000019538132) * Y
						- 0.000676904986) * Y + 0.001390604284) * Y
						- 0.00079462082) * Y - 0.002034254874) * Y
						+ 0.006549791214) * Y - 0.010557625006) * Y
						+ 0.011630447319) * Y - 0.009279453341) * Y
						+ 0.005353579108) * Y - 0.002141268741) * Y
						+ 0.000535310849) * Y + 0.999936657524;
				}

				if ((X + 1.0) < (1.0 - X))
					NormalZx = ((double)(X)+1.0);
				else
					NormalZx = (1.0 - (double)(X));

			}
		}
		else {
			temp = ((fabs(Z) - 5.999999) * 10);
			Temp2 = pow(1.6, temp);
			WinP = pow(10, -9);
			WinP = WinP / Temp2;
			NormalZx = WinP;

		}

		return (NormalZx);

	}

	double ChiPVal(double X) {
		long double PValHolder;
		double ChiPValx;
		if (X == 0)
			ChiPValx = 1;
		else {

			PValHolder = (NormalZ(-sqrt(X)));

			if (PValHolder == 0) {//< 0.0000000001){
								  //if (X > 35)
								  //	Y=X;
								  //else
								  //	X=35.1;
				PValHolder = 0.0;
				PValHolder = pow(10, -9);
				PValHolder = PValHolder / (X - 34);
			}

			ChiPValx = (double)(PValHolder);
		}
		return(ChiPValx);

	}




	int MyMathFuncs::GrowMChiWinP(int LO, int RO, int LenXoverSeq, int HWindowWidth, int TWin, int MaxY, int LS, int A, int C, int MaxFailCount, double *MPV, int *WinWin, double *MChi, int *TopL, int *TopR, int *TopLO, int *TopRO, unsigned char *Scores) {
		//scores LS,2 
		int off1, FailCount;
		double hMChi, tpv, lxdh, lx, h, t, E, htpv, dA, dB, dC, dD, dTWin;

		off1 = MaxY*(LS + 1);
		FailCount = 0;
		lx = LenXoverSeq;
		h = HWindowWidth;
		lxdh = lx / h;
		htpv = *MChi;
		//TW2 = TWin*2;
//#pragma loop(no_vector)
		dA = (double)(A);
		dC = (double)(C);
		dTWin = (double)(TWin);
		while (FailCount <= MaxFailCount) {
			dA = dA + Scores[LO + off1];
			dC = dC + Scores[RO + off1];
			dB = TWin - A;
			dD = TWin - C;
			if (dB + dD > 0) {
				E = (double)(dA*dD - dB*dC);
				E = E*E;
				E = E * 2;
				hMChi = E / (double)((TWin*(dA + dC)*(dB + dD)));
				
				if (hMChi >= htpv) {
					htpv = hMChi;
					*TopLO = LO;
					*TopRO = RO;
					*MChi = hMChi;
					*WinWin = TWin;
					//*MPV = tpv;
					FailCount = 0;
					*TopL = (int)dA;
					*TopR = (int)dC;
				}
				else {
					if (TWin == HWindowWidth && *TopL == 0 && *TopR == 0) {
						*TopL = A;
						*TopR = C;
					}
					FailCount++;
					if (FailCount > MaxFailCount)
						return(1);
				}
			}
			else {
				FailCount++;
				if (FailCount > MaxFailCount)
					return(1);
				hMChi = 0;
			}
			RO++;
			LO--;
			if (LO < 1)
				LO = LenXoverSeq;
			if (RO > LenXoverSeq*2)
				RO = RO - LenXoverSeq*2;
			if (RO > LenXoverSeq)
				RO = RO - LenXoverSeq;
			TWin++;
		}
		return(1);
	}


	int MyMathFuncs::GrowMChiWinP2(int MaxABWin, int LO, int RO, int LenXoverSeq, int HWindowWidth, int TWin, int MaxY, int LS, int A, int C, int MaxFailCount, double *MPV, int *WinWin, double *MChi, int *TopL, int *TopR, int *TopLO, int *TopRO, unsigned char *Scores, float *chitable, int *chimap) {
		//scores LS,2 
		int off1, FailCount, ctoff, os1, os2, B, D;
		double hMChi, tpv, lxdh, lx, h, t, E, dB, dD, htpv, dA, dC, dTWin;
		
		off1 = MaxY*(LS + 1);
		FailCount = 0;
		lx = LenXoverSeq;
		h = HWindowWidth;
		lxdh = lx / h;
		htpv = *MChi;
		//TW2 = TWin*2;
		//TW2 = TWin*2;
		//#pragma loop(no_vector)
		while (FailCount <= MaxFailCount) {
			A = A +  Scores[LO + off1];
			C = C + Scores[RO + off1];
			B = TWin - A;
			D = TWin - C;
			if (B + D > 0) {
				if (TWin <= MaxABWin) {

					hMChi = chitable[chimap[TWin] + A + C*(TWin + 1)];
					//tpv = ChiPVal(hMChi) * (h / t);
					if (hMChi >= htpv) {
						htpv = hMChi;
						*TopLO = LO;
						*TopRO = RO;
						*MChi = hMChi;
						*WinWin = TWin;
						//*MPV = tpv;
						FailCount = 0;
						*TopL = A;
						*TopR = C;
					}
					else {
						if (TWin == HWindowWidth && *TopL == 0 && *TopR == 0) {
							*TopL = A;
							*TopR = C;
						}
						FailCount++;
						if (FailCount > MaxFailCount)
							return(1);
					}
				}

				else {
					dA = (double)A;
					dC = (double)C;
					dB = (double)B;
					dD = (double)D;
					dTWin = (double)TWin;
					//if (B + D > 0) {
						E = (double)(dA*dD - dB*dC);
						E = E*E;
						E = E * 2;
						hMChi = E / (double)((dTWin*(dA + dC)*(dB + dD)));
						if (hMChi < 0)
							hMChi = fabs(hMChi);
						
						if (hMChi >= htpv) {
							htpv = hMChi;
							*TopLO = LO;
							*TopRO = RO;
							*MChi = hMChi;
							*WinWin = TWin;
							//*MPV = tpv;
							FailCount = 0;
							*TopL = (int)(dA);
							*TopR = (int)(dC);
						}
						else {
							if (TWin == HWindowWidth && *TopL == 0 && *TopR == 0) {
								*TopL = (int)(dA);
								*TopR = (int)(dC);
							}
							FailCount++;
							if (FailCount > MaxFailCount)
								return(1);
						}
					/*}
					else {
						FailCount++;
						if (FailCount > MaxFailCount)
							return(1);
						hMChi = 0;
					}*/
				}
			}
			else {
				FailCount++;
				if (FailCount > MaxFailCount)
					return(1);
				hMChi = 0;
			}
			RO++;
			LO--;
			if (LO < 1)
				LO = LenXoverSeq;
			if (RO > LenXoverSeq * 2)
				RO = RO - LenXoverSeq * 2;
			if (RO > LenXoverSeq)
				RO = RO - LenXoverSeq;
			TWin++;
		}
		return(1);
	}

	int MyMathFuncs::FindMChiP(int LenSeq, int LenXoverSeq, int *MaxX, short int *MaxY, double *MChi, double *ChiVals)
	{
		int X, SO, tMaxX;
		short int Y, tMaxY;
		double ChiV, tMChi;
		tMaxX = -1;
		tMaxY = -1;
		tMChi = 0;
		ChiV = 0;
		SO = LenSeq + 1;
		for (X = 0; X < LenXoverSeq; X++) {
			for (Y = 0; Y <= 2; Y++) {
				ChiV = *(ChiVals + X + Y*SO);
				if (ChiV > tMChi) {
					tMChi = ChiV;
					tMaxX = X;
					tMaxY = Y;
				}
			}
		}
		*MChi = tMChi;
		*MaxX = tMaxX;
		*MaxY = tMaxY;
		return 1;
	}



	int MyMathFuncs::TSeqPermsP(int Seq1, int Seq2, int Seq3, int lseq, int *THold, unsigned char *tMissingData, short int *SeqNum, short int *SeqRnd) {
		int Z, off1, off2, os3, os1, os2, NewPos, o1, o2, no1, no2;
		double rn, rm, rf, ld, dh;
		ld = (double)(lseq);

		rn = rand();
		rm = RAND_MAX;
		dh = ld / rm;


		off1 = lseq + 1;
		off2 = off1 * 2;
		os1 = Seq1*off1;
		os2 = Seq2*off1;
		os3 = Seq3*off1;
		for (Z = 1; Z <= lseq; Z++) {

			SeqRnd[Z] = SeqNum[Z + os1];
			SeqRnd[Z + off1] = SeqNum[Z + os2];
			SeqRnd[Z + off2] = SeqNum[Z + os3];
		}

		for (Z = 1; Z <= lseq; Z++) {

			o1 = Z + off1;
			o2 = Z + off2;
			if (tMissingData[Z] == 0) {
				if (tMissingData[o1] == 0) {
					if (tMissingData[o2] == 0) {
						rn = rand();
						rf = (rn*dh) + 0.49;
						NewPos = (int)(rf);
						no1 = NewPos + off1;
						no2 = NewPos + off2;
						if (tMissingData[NewPos] == 0) {
							if (tMissingData[no1] == 0) {
								if (tMissingData[no2] == 0) {
									THold[0] = SeqRnd[Z];
									THold[1] = SeqRnd[o1];
									THold[2] = SeqRnd[o2];
									SeqRnd[Z] = SeqRnd[NewPos];
									SeqRnd[o1] = SeqRnd[no1];
									SeqRnd[o2] = SeqRnd[no2];
									SeqRnd[NewPos] = THold[0];
									SeqRnd[no1] = THold[1];
									SeqRnd[no2] = THold[2];
								}
							}
						}
					}
				}
			}

		}

		return(1);
	}


	int  MyMathFuncs::MakeBanWinP(int UBBW, int Seq1, int Seq2, int Seq3, int HWindowWidth, int LS, int LenXoverSeq, int *BanWin, unsigned char *MDMap, unsigned char *MissingData, int *XPosDiff, int *XDiffPos) {

		int X, Y, s1o, s2o, s3o, lx, xpd;
		//empty banwin and mdmap
		s1o = Seq1*(LS + 1);
		s2o = Seq2*(LS + 1);
		s3o = Seq3*(LS + 1);

		for (X = 0; X < UBBW; X++)
			BanWin[X] = 0;

		for (X = 0; X < LenXoverSeq + 2; X++)
			MDMap[X] = 0;
		lx = 0;
		for (X = 1; X <= LS; X++) {
			if (lx > X)
				break;
			lx = X;
			if (MissingData[X + s1o] == 1 || MissingData[X + s2o] == 1 || MissingData[X + s3o] == 1) {
				xpd = XPosDiff[X];
				//half windows are allowed to end on this position but not start on it or traverse it
				MDMap[xpd] = 1;
				if (xpd + HWindowWidth - 1 <= LenXoverSeq) {
					for (Y = xpd; Y < xpd + HWindowWidth; Y++)
						BanWin[Y] = 1;
				}
				else {
					for (Y = xpd; Y <= LenXoverSeq; Y++)
						BanWin[Y] = 1;

					for (Y = 0; Y < xpd + HWindowWidth - 1 - LenXoverSeq; Y++)
						BanWin[Y] = 1;
				}
				if (xpd < LenXoverSeq) {
					//half windows are allowed to start on this position but not end on it or traverse it
					MDMap[xpd + 1] = 1;

					if (xpd + 2 - HWindowWidth > 0) {
						for (Y = xpd + 2 - HWindowWidth; Y < xpd + 2; Y++)
							BanWin[Y] = 1;
					}
					else {
						for (Y = 0; Y < xpd + 2; Y++)
							BanWin[Y] = 1;
						for (Y = xpd + 2 - HWindowWidth + LenXoverSeq; Y <= LenXoverSeq; Y++)
							BanWin[Y] = 1;

					}
					if (XDiffPos[xpd + 1] > XDiffPos[xpd]) {
						if (X < XDiffPos[xpd + 1]) {
							X = XDiffPos[xpd + 1];
							if (X >= LS)
								break;
						}
					}
					else {

						MDMap[1] = 1;
						break;
					}
				}

				else {
					MDMap[1] = 1;
					break;

				}
			}
		}
		if (MDMap[LenXoverSeq] == 1 || MDMap[1] == 1) {
			for (X = LenXoverSeq - HWindowWidth + 2; X <= LenXoverSeq; X++)
				BanWin[X] = 1;

		}
		return(1);
	}

	int MyMathFuncs::WinScoreCalcP(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int Seq1, int Seq2, int Seq3, unsigned char *Scores, int *XDiffPos, short int *SeqNum, int *WinScores)
	{

		int goon, X, LO, RO, FO, SO, TO, TO2, target, ls1, ls2, ls3, xdp, s1, s2, s3;
		TO = LenSeq;
		TO2 = LenSeq * 2;
		FO = UBWS+1;  //LenSeq + HWindowWidth * 2;
		SO = FO * 2;
		goon = 0;

		target = (FO + 1) * 2;

		for (X = 0; X <= target; X++)
			WinScores[X] = 0;
		////Calculate scores per position
		//X = 0;
		//while (X <= LenXoverSeq) {

		//	*(WinScores + X) = 0;
		//	X++;
		//}

		//X = 0;
		//while (X <= LenXoverSeq) {
		//	
		//	*(WinScores + X + FO) = 0;
		//	X++;
		//	
		//}


		//X = 0;
		//while (X <= LenXoverSeq) {
		//	*(WinScores + X + SO) = 0;
		//	X++;
		//}
		ls1 = LenSeq*Seq1;
		ls2 = LenSeq*Seq2;
		ls3 = LenSeq*Seq3;
//#pragma omp parallel
//		{
//			
//#pragma omp  for private (xdp, s1, s2, s3)
			for (X = 1; X <= LenXoverSeq; X++) {
				xdp = XDiffPos[X];
				s1 = SeqNum[xdp + ls1];
				s2 = SeqNum[xdp + ls2];
				s3 = SeqNum[xdp + ls3];
				Scores[X] = (unsigned char)(s1 == s2);
				//for (X = 1; X <= LenXoverSeq; X++)
				Scores[X + TO] = (unsigned char)(s1 == s3);
				//for (X = 1; X <= LenXoverSeq; X++)
				Scores[X + TO2] = (unsigned char)(s2 == s3);

			}
		//}
		//calculate score for 0 window (actually the last window)
		for (X = (LenXoverSeq - HWindowWidth + 1); X <= LenXoverSeq; X++) {
			*(WinScores) = *(WinScores)+*(Scores + X);
			*(WinScores + FO) = *(WinScores + FO) + *(Scores + X + TO);
			*(WinScores + SO) = *(WinScores + SO) + *(Scores + X + TO2);
		}

		//Calculate scores for windows traversing the left end
		for (X = 1; X <= HWindowWidth; X++) {
			LO = ((LenXoverSeq - HWindowWidth) + X);
			*(WinScores + X) = *(WinScores + X - 1) - *(Scores + LO) + *(Scores + X);
			*(WinScores + X + FO) = *(WinScores + X + FO - 1) - *(Scores + LO + TO) + *(Scores + X + TO);
			*(WinScores + X + SO) = *(WinScores + X + SO - 1) - *(Scores + LO + TO2) + *(Scores + X + TO2);

		}

		//Calculate scores for internal windows
		for (X = HWindowWidth + 1; X <= LenXoverSeq; X++) {
			LO = X - HWindowWidth;
			*(WinScores + X) = *(WinScores + X - 1) - *(Scores + LO) + *(Scores + X);
			*(WinScores + X + FO) = *(WinScores + X + FO - 1) - *(Scores + LO + TO) + *(Scores + X + TO);
			*(WinScores + X + SO) = *(WinScores + X + SO - 1) - *(Scores + LO + TO2) + *(Scores + X + TO2);
		}



		//Calculate scores for windows traversing the right end
		for (X = LenXoverSeq + 1; X < LenXoverSeq + HWindowWidth; X++) {
			LO = (X - HWindowWidth);
			RO = X - LenXoverSeq;
			*(WinScores + X) = *(WinScores + X - 1) - *(Scores + LO) + *(Scores + RO);
			*(WinScores + X + FO) = *(WinScores + X + FO - 1) - *(Scores + LO + TO) + *(Scores + RO + TO);
			*(WinScores + X + SO) = *(WinScores + X + SO - 1) - *(Scores + LO + TO2) + *(Scores + RO + TO2);
		}



		return (1);
	}

	unsigned char MyMathFuncs::FixOverlapsP(unsigned char DoneThisOne, int CurBegin, int CurEnd, int CurProg, int X, int Y, int MSX, float LSAdjust, int UBPD, int UBXONC1, int UBXONC2, unsigned char *ProgDo, short int *XOverNoComponent, short int *MaxXONo) {
		int Z, RN, os1, os2, os3, os4, os5, holder;
		unsigned char dto;
		double J, K, L;
		os1 = UBXONC1 + 1;
		os2 = os1*(UBXONC2 + 1);
		os3 = CurProg + MSX*(UBPD + 1);
		os4 = X*os1;
		dto = DoneThisOne;
		//os5=
		//srand((int)(CurProb));
		L = (double)(RAND_MAX);
		os5 = CurProg + os4;
		for (Z = CurBegin; Z <= CurEnd; Z++) {

			//if (CurProb > 0){

			K = (double)(rand());
			J = (K / L) * 3;

			RN = (int)(J)+1;
			//RN = (int)((rand() / RAND_MAX)*100 +1);//Int((3 * Rnd) + 1)
			//return(RN);
			if (RN == 2) {
				holder = (int)(Z * LSAdjust);
				if (XOverNoComponent[os5 + holder*os2] == 0)
					XOverNoComponent[os5 + holder*os2] = Y;

				ProgDo[os3] = 1;

				if (dto == 0) {
					MaxXONo[MSX] = MaxXONo[MSX] + 1;
					dto = 1;
				}
			}

			// }
		}

		return(dto);
	}

	int MyMathFuncs::FindSubSeqDP(int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, short int *seqnum, int *xdiffpos, int *xposdiff)
	{
		int s1, s2, s3, y = 0;
		int os1, os2, os3;
		int x;
		os1 = seq1*lenseq;
		os2 = seq2*lenseq;
		os3 = seq3*lenseq;
		for (x = 1; x < lenseq; x++) {
			xposdiff[x - 1] = y;
			s1 = *(seqnum + x + os1);


			if (s1 != 46) {
				s2 = *(seqnum + x + os2);

				if (s2 != 46) {
					s3 = *(seqnum + x + os3);
					if (s3 != 46) {

						if (s1 != s2 || s1 != s3) {
							if (s1 == s2 || s1 == s3) {
								y++;
								xdiffpos[y] = x;
							}
						}
					}
				}
			}
		}
		return(y);
	}

	int MyMathFuncs::FindSubSeqDP2(int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, short int *seqnum, int *xdiffpos)
	{
		int s1, s2, s3, y = 0;
		int os1, os2, os3;
		int x;
		os1 = seq1*lenseq;
		os2 = seq2*lenseq;
		os3 = seq3*lenseq;
		for (x = 1; x < lenseq; x++) {
			
			s1 = *(seqnum + x + os1);


			if (s1 != 46) {
				s2 = *(seqnum + x + os2);

				if (s2 != 46) {
					s3 = *(seqnum + x + os3);
					if (s3 != 46) {

						if (s1 != s2 || s1 != s3) {
							if (s1 == s2 || s1 == s3) {
								y++;
								xdiffpos[y] = x;
							}
						}
					}
				}
			}
		}
		return(y);
	}

	int MyMathFuncs::FindSubSeqDP6(int UBFSS, int ubcs1, unsigned char *FSSRDP, unsigned char *CS, int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, int *LXOS, int UBXDP, int *XDP, int *XPD)
	{
		int s1, s2, s3, y = 0, s1o, s2o, s3o, se1, se2, se3, osf, h, uo1, uo2, y0, y1, y2, xh;
		int os1, os2, os3, os4, os5;
		int x, z;
		os1 = seq1*lenseq;
		os2 = seq2*lenseq;
		os3 = seq3*lenseq;

		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);

		uo1 = UBXDP + 1;
		uo2 = (UBXDP + 1) * 2;

		os4 = (UBFSS + 1) * 4;
		os5 = (UBFSS + 1)*(UBFSS + 1) * 4;
		y0 = 0;
		y1 = 0;
		y2 = 0;
		for (x = 1; x < ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];
			osf = se1 * 4 + os4 * se2 + os5 * se3;
			if (FSSRDP[3 + osf] > 0) {
				for (int z = 0; z <= 2; z++) {
					h = FSSRDP[z + osf];
					xh = (x - 1) * 3 + z + 1;
					if (h == 1) {
						y0++;
						y1++;
						XDP[y0] = xh;
						XDP[y1 + uo1] = xh;
					}
					else if (h == 2) {
						y0++;
						y2++;
						XDP[y0] = xh;
						XDP[y2 + uo2] = xh;
					}
					else if (h == 3) {
						y1++;
						y2++;
						XDP[y1 + uo1] = xh;
						XDP[y2 + uo2] = xh;
					}
					XPD[xh] = y0;
					XPD[xh + uo1] = y1;
					XPD[xh + uo2] = y2;
				}
			}
			else {
				xh = (x - 1) * 3 + 1;
				for (int z = 0; z <= 2; z++) {
					XPD[xh] = y0;
					XPD[xh + uo1] = y1;
					XPD[xh + uo2] = y2;
					XPD[xh + 1] = y0;
					XPD[xh + 1 + uo1] = y1;
					XPD[xh + 1 + uo2] = y2;
					XPD[xh + 2] = y0;
					XPD[xh + 2 + uo1] = y1;
					XPD[xh + 2 + uo2] = y2;
				}

			}

		}
		//for (x = 0; x <= lenseq; x++)

		LXOS[0] = y0;
		LXOS[1] = y1;
		LXOS[2] = y2;

		return(y);
	}
	int MyMathFuncs::FindSubSeqDP3(int UBFSS, int ubcs1, unsigned char *FSSRDP, unsigned char *CS, int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, int *LXOS, int UBXDP, int *XDP)
	{
		int s1, s2, s3, y = 0, s1o, s2o,s3o, se1,se2,se3, osf, h, uo1, uo2, y0, y1,y2, xh;
		int os1, os2, os3, os4, os5;
		int x, z;
		os1 = seq1*lenseq;
		os2 = seq2*lenseq;
		os3 = seq3*lenseq;

		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);

		uo1 = UBXDP + 1;
		uo2 = (UBXDP + 1) * 2;

		os4 = (UBFSS + 1) * 4;
		os5 = (UBFSS + 1)*(UBFSS + 1) * 4;
		y0 = 0;
		y1 = 0;
		y2 = 0;
		for (x = 1; x < ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];
			osf = se1 * 4 + os4 * se2 + os5 * se3;
			if (FSSRDP[3 + osf] > 0) {
				for (int z = 0; z <= 2; z++) {
					h = FSSRDP[z + osf];
					xh = (x-1)*3 + z + 1;
					
					if (h == 1) {
						y0++;
						y1++;
						XDP[y0] = xh;
						XDP[y1 + uo1] = xh;
					}
					else if (h == 2) {
						y0++;
						y2++;
						XDP[y0] = xh;
						XDP[y2 + uo2] = xh;
					}
					else if (h == 3) {
						y1++;
						y2++;
						XDP[y1 + uo1] = xh;
						XDP[y2 + uo2] = xh;
					
					}
				}
			}
		}
		LXOS[0] = y0;
		LXOS[1] = y1;
		LXOS[2] = y2;

		return(y);
	}

	int MyMathFuncs::FindSubSeqCP(int lenseq, short int nextno, short int seq1, short int seq2, short int seq3, short int *seqnum, int *xdiffpos, int *xposdiff)
	{
		short int s1, s2, s3;
		int x, y, o1, o2, o3;
		y = 0;
		o1 = seq1*lenseq;
		o2 = seq2*lenseq;
		o3 = seq3*lenseq;
		for (x = 1; x < lenseq; x++) {
			*(xposdiff + x - 1) = y;
			s1 = *(seqnum + x + o1);
			if (s1 != 46) {
				s2 = *(seqnum + x + o2);
				if (s2 != 46) {
					if (s1 != s2) {
						if (*(seqnum + x + o3) != 46) {
							y++;
							*(xdiffpos + y) = x;
						}
					}
					else {
						s3 = *(seqnum + x + o3);
						if (s3 != 46) {
							if (s1 != s3) {
								y++;
								*(xdiffpos + y) = x;
							}
						}
					}
				}
			}
		}
		return(y);
	}

	int MyMathFuncs::SmoothChiValsP(int LenXoverSeq, int LenSeq, double *ChiVals, double *SmoothChi)
	{
		short int Y;
		int OS, RO, X;
		double RunCount;
		int qWindowSize = 5;
		OS = LenSeq + 1;
		RO = qWindowSize * 2 + 1;
		for (Y = 0; Y <= 2; Y++) {
			RunCount = 0;
			for (X = 0 - qWindowSize; X <= 1 + qWindowSize; X++) {
				if (X < 1)
					RunCount += *(ChiVals + LenXoverSeq + X + Y*OS);
				else
					RunCount += *(ChiVals + X + Y*OS);
			}
			*(SmoothChi + Y*OS) = RunCount / RO;
			for (X = 1 - qWindowSize; X < LenXoverSeq - qWindowSize; X++) {
				if (X > 0) {
					if (X + RO <= LenXoverSeq)
						RunCount = RunCount - *(ChiVals + X + Y*OS) + *(ChiVals + X + RO + Y*OS);
					else
						RunCount = RunCount - *(ChiVals + X + Y*OS) + *(ChiVals + X + RO - LenXoverSeq + Y*OS);
				}
				else
					RunCount = RunCount - *(ChiVals + LenXoverSeq + X + Y*OS) + *(ChiVals + X + RO + Y*OS);

				*(SmoothChi + X + qWindowSize + Y*OS) = RunCount / RO;
			}

		}

		return 1;
	}

	double MyMathFuncs::ProbCalcP(double *fact, int xoverlength, int numincommon, double indprob, int lenxoverseq)
	{
		int z;// , target;
		long double nfactorial;
		long double mfactorial, nmfactorial;
		long double hold, hold2;

		double probability = 0.0;
		double td1, td2, prob = 0.0;
		td1 = (double)(xoverlength);
		td2 = (double)(lenxoverseq);
		//addjustfactor = 1;

		nfactorial = (long double)(fact[xoverlength]);
		//for(y=1; y<=xoverlength;y++){
		//	nfactorial = nfactorial * y;
		//}
//#pragma omp parallel for private ( target, mfactorial, nmfactorial, hold, hold2)
//#pragma omp parallel
//		{
			for (z = numincommon; z <= xoverlength; z++) {
				//target = ;

				mfactorial = (long double)(fact[z]);
				//for (y = 2; y <= z;y++){
				//    mfactorial = mfactorial * y;
				//}
				nmfactorial = (long double)(fact[xoverlength - z]);
				//for (a = 2; a <= target;a++){
				//	nmfactorial = nmfactorial * a;
				//}

				hold = mfactorial* nmfactorial;
				hold2 = nfactorial / hold;
				probability = probability + pow(indprob, z) * pow((1 - indprob), (xoverlength - z)) * hold2;

			}
//		}
		//td2 = td2/td1;
		probability = probability * td2 / td1;

		/*if (1 - probability == 1) {
		if (probability < 1)
		probability = probability * (lenxoverseq / xoverlength);

		}
		else{
		if (probability < 0)
		probability = -1;
		else if (probability < 1){
		td1 = lenxoverseq / xoverlength;
		td2 = 1 - probability;
		probability = 1 - pow(td1,td2);//((1 - probability) ^ ((lenxoverseq / xoverlength)));

		}
		}*/


		//if (probability < 1){
		//	if 1-probability >
		//	probability = 1 - pow((1 - probability),(lenxoverseq / xoverlength));
		//                
		//}
		//else{
		//	probability=1;
		//}
		//if (probability < 0.000000000000001){
		//	probability = 0;
		//}

		return (probability);
	}


	double MyMathFuncs::ProbCalcP2(double *fact3x3, int ub3x3, int xoverlength, int numincommon, double indprob, int lenxoverseq)
	{
		int os1, os2;// , target;
		int z;
		double dz;
		
		double hold2, hold3;

		double probability = 0.0;
		double td1, td2, prob = 0.0;
		td1 = (double)(xoverlength);
		td2 = (double)(lenxoverseq);
		os1 = ub3x3 + 1;
		os2 = os1*os1;
		
//#pragma omp parallel for private (z, hold2, h3) 
		for (z = numincommon; z <= xoverlength; z++) {
			hold2 = fact3x3[xoverlength + z * os1 + (xoverlength - z) * os2]; //nfactorial / (fact[z] * fact[xoverlength - z]);
			dz = (double)(z);
			hold3 = pow(indprob, dz);
			dz = (double)(xoverlength - z);
			hold3 = hold3*pow((1 - indprob), dz);
			probability += hold3 * hold2;


		}

		probability = probability * td2 / td1;



		return (probability);
	}



	double  MyMathFuncs::FastSimilarityBP(int df, int reps, int ISDim, int Nextno, int UBX, float *Valid, float *Diffs, short int *XCVal, short int *IntegerSeq, unsigned char *CompressValid, unsigned char *CompressDiffs, float *DistCheckB, int *weightmod)
	{

		int  X, A, Y, B, C, E, off5, off2, off1, off3, o3, o4, off6, off7;
		double dX, v1, d1, th2, th3;

		off6 = reps + 1;
		off7 = (Nextno + 1)*(reps + 1);

		for (X = 0; X <= Nextno; X++) {

			off1 = (ISDim + 1)*X;
			for (A = 1; A <= UBX; A++)
				XCVal[A] = IntegerSeq[A + off1];


			for (Y = X + 1; Y <= Nextno; Y++) {
				for (E = 0; E <= reps; E++) {
					Valid[E] = 0;
					Diffs[E] = 0;
				}
				off2 = (ISDim + 1)*Y;
				for (A = 1; A <= UBX; A++) {


					C = IntegerSeq[A + off2];
					B = XCVal[A];

					off3 = B + 626 * C;
					off5 = (reps + 1)*(A - 1);



					for (E = df; E <= reps; E++)
						Valid[E] = Valid[E] + weightmod[E + off5] * CompressValid[off3];

					if (B != C) {
						for (E = df; E <= reps; E++)
							Diffs[E] = Diffs[E] + weightmod[E + off5] * CompressDiffs[off3];

					}
				}
				for (E = df; E <= reps; E++) {


					v1 = (double)(Valid[E]);
					d1 = (double)(Diffs[E]);
					if (v1 > 0) {
						dX = (double)((v1 - d1) / v1);
						if (dX > 0.25) {
							th2 = (4.0 * dX - 1.0) / 3.0;
							th3 = log(th2);
							dX = -0.75*th3;
						}
						else
							dX = 10.0;
					}
					else
						dX = 10.0;

					o3 = E + X*off6 + Y*off7;
					o4 = E + Y*off6 + X*off7;
					DistCheckB[o3] = (float)(dX);
					DistCheckB[o4] = (float)(dX);
				}


			}


		}


		return(1);
	}


int MyMathFuncs::CheckMatrixP(int * MinS,int *ISeqs, int NextNo, int SCO, int MinSeqSize, int UBMP, unsigned char *MissPair, int UBPV,  float *PermValid, int UBSV, float *SubValid,int UBF, float *FMat, float *SMat, int *ValtotF, int *ValtotS){
	int x, Y, Z, GoOn, MinNum, LoopNo, MinSize;
	GoOn = 0;
	for (x = 0; x<= NextNo; x++){
        
        
        for (Y = x + 1; Y <= NextNo; Y++){
            
            if (FMat[Y + x*(UBF+1)] < 2.99 || FMat[Y + x*(UBF + 1)] > 3.01){
                if (PermValid[Y + x*(UBPV+1)] - SubValid[Y + x*(UBSV + 1)] < MinSeqSize){
					MissPair[x + Y*(UBMP+1)] = 1;
					MissPair[Y + x*(UBMP + 1)] = 1;
					for (Z = 0; Z <= 2; Z++){
						ValtotF[x] = ValtotF[x] + PermValid[ISeqs[Z] + x*(UBPV + 1)] - SubValid[ISeqs[Z] + x*(UBSV + 1)];
						ValtotF[Y] = ValtotF[Y] + PermValid[ISeqs[Z] + Y*(UBPV + 1)] - SubValid[ISeqs[Z] + Y*(UBSV + 1)];
                    }
					GoOn = 1;
                }
                
				if (SubValid[Y + x*(UBSV + 1)]  < SCO && (SMat[Y + x*(UBF + 1)] < 2.99 || SMat[Y + x*(UBF + 1)] > 3.01)){
					MissPair[x + Y*(UBMP + 1)] = 1;
					MissPair[Y + x*(UBMP + 1)] = 1;
					for (Z = 0; Z <= 2; Z++) {
						ValtotS[x] = ValtotS[x] + SubValid[ISeqs[Z] + x*(UBSV + 1)];
						ValtotS[Y] = ValtotS[Y] + SubValid[ISeqs[Z] + Y*(UBSV + 1)];
					}
					GoOn = 1;
                    
                }
            }
        }
    }
	
	MinNum = 0;
	LoopNo = 0;
    
	while (GoOn == 1) {

		MinSize = 1000000000;

		for (x = 0; x <= NextNo; x++) {
			if (ValtotS[x] > 0) {
				if (ValtotS[x] < MinSize) {
					MinNum = 0;
					MinSize = ValtotS[x];
					MinS[0] = x;
				}
				else if (ValtotS[x] == MinSize) {
					MinNum++;
					MinS[MinNum] = x;
				}
			}
			if (ValtotF[x] > 0) {
				if (ValtotF[x] < MinSize) {
					MinNum = 0;
					MinSize = ValtotF[x];
					MinS[0] = x;
				}
				else if (ValtotF[x] == MinSize) {
					MinNum++;
					MinS[MinNum] = x;
				}
			}
		}
		for (x = 0; x <= MinNum; x++) {
			ValtotS[MinS[x]] = 0;
			ValtotF[MinS[x]] = 0;
            for (Y = 0; Y <= NextNo; Y++){
				FMat[MinS[x] + Y*(UBF + 1)] = 3.0;
                FMat[Y + MinS[x]*(UBF+1)] = 3.0;
                SMat[MinS[x] + Y*(UBF + 1)] = 3.0;
                SMat[Y + MinS[x] * (UBF + 1)] = 3.0;
				MissPair[Y + MinS[x] * (UBMP + 1)] = 0;
				MissPair[MinS[x] + Y * (UBMP + 1)] = 0;
               
            }
        }
        
        GoOn = 0;
        for (x = 0; x <= NextNo; x++){
            for (Y = 0; Y <= NextNo; Y++){
            
                if (MissPair[x + Y*(UBMP+1)] == 1){
					GoOn = 1;
					break;
                }
           }
            
			if (GoOn == 1)
				break;
            
           
        }
        
        LoopNo++;
		if (LoopNo > NextNo * 10)
			break;
    }
	
	for (x = 0; x <= NextNo; x++){
		GoOn = 0;
		if (FMat[x + x*(UBF + 1)] < 2.999 || FMat[x + x*(UBF + 1)] > 3.001) {
			for (Y = 0; Y <= NextNo; Y++) {
				if (FMat[x + Y*(UBF + 1)] < 2.999 && x != Y) {
					GoOn = 1;
					break;
				}
					
			}

			if (GoOn==0)
				FMat[x + x*(UBF + 1)] = 3.000;

				
		}
    }
	
	return(1);
}


int MyMathFuncs::DoRecode(int NextNo, int LenStrainSeq0, int UBSN, short int *SeqNum, int UBRec, unsigned char *Recoded, unsigned char *NucMat, int UBRep, unsigned char *Replace){
int Y,x, NN;
	for (Y = 0; Y <= NextNo; Y++){
		for (x = 1; x <= LenStrainSeq0; x++){
			NN = SeqNum[x + Y*(UBSN + 1)];
			NN = NucMat[NN];
			Recoded[x + Y*(UBRec + 1)] = Replace[x + NN*(UBRep + 1)];
            
        }
   }
	return(1);
}


int MyMathFuncs::DoRecodeP(int NextNo, int LenStrainSeq0, int UBSN, short int *SeqNum, int UBRec, unsigned char *Recoded, unsigned char *NucMat, int UBRep, unsigned char *Replace) {
	int Y, x, NN;


	int procs;
	procs = omp_get_num_procs();
	procs = procs / 2 - 1;
	if (procs < 3)
		procs = 3;
	omp_set_num_threads(procs);

#pragma omp parallel for private (Y, x,NN)
	for (Y = 0; Y <= NextNo; Y++) {
		for (x = 1; x <= LenStrainSeq0; x++) {
			NN = SeqNum[x + Y*(UBSN + 1)];
			NN = NucMat[NN];
			Recoded[x + Y*(UBRec + 1)] = Replace[x + NN*(UBRep + 1)];

		}
	}
	return(1);
}


int MyMathFuncs::CountNucs(int NextNo, int LenStrainSeq0, int UBSN, short int *SeqNum, unsigned char *NucMat, int UBNC, int *NucCount) {
	int Y, x, NN;
	for (Y = 0; Y <= NextNo; Y++) {
		for (x = 1; x <= LenStrainSeq0; x++) {
			NN = SeqNum[x + Y*(UBSN + 1)];
			NN = NucMat[NN];
			NucCount[x + NN*(UBNC + 1)] = NucCount[x + NN*(UBNC + 1)] + 1;
			//Recoded[x + Y*(UBRec + 1)] = Replace[x + NN*(UBRep + 1)];

		}
	}
	return(1);
}

int MyMathFuncs::MakeVarSiteMap(int SWin, int LenVarSeq, short int *VarSiteMap, float *VarSiteSmooth) {
	int x, Z;
	float Tot;
	
	Tot = 0.0;
	for (x = 1 - SWin; x <= 1 + SWin; x++){

		if (x < 1)
			Z = LenVarSeq + x;
		else if (x > LenVarSeq)
			Z = x - LenVarSeq;
		else
			Z = x;
		
		Tot = Tot + (float)(VarSiteMap[Z]);
		/*if (VarSiteMap[Z] != 2)
			ZZ = VXPos[Z];*/
		
	}

	VarSiteSmooth[1] = Tot / (float)((SWin * 2 + 1) * 2);
	for (x = 2; x <= LenVarSeq; x++){
		Z = x - SWin - 1;
		if (Z < 1)
			Tot = Tot - (float)(VarSiteMap[LenVarSeq + Z]);
		else if (Z > LenVarSeq)
			Tot = Tot - (float)(VarSiteMap[Z - LenVarSeq]);
		else
			Tot = Tot - (float)(VarSiteMap[Z]);
		
		Z = x + SWin;
		if (Z > LenVarSeq)
			Tot = Tot + (float)(VarSiteMap[Z - LenVarSeq]);
		else
			Tot = Tot + (float)(VarSiteMap[Z]);
		
		VarSiteSmooth[x] = Tot / (float)((SWin * 2 + 1) * 2);
		
	}
	return(1);

}

int MyMathFuncs::FtoFA(int NSeqs, int LenStrainSeq0, int UBTS,int *TraceSeqs,int UBTFA, float *tFAMat, int UBFA, float *FAMat) {
	//this is also good for smat version
	int x, Y;
     for (x = 0; x <= NSeqs; x++){
        for (Y = x + 1; Y <= NSeqs; Y++){
			//tFAMat[x + Y*(UBTFA+1)] = round(tFAMat[x + Y*(UBTFA + 1)] * 10000) / 10000;
			//tFAMat[Y + x*(UBTFA + 1)] = tFAMat[x + Y*(UBTFA + 1)];
			FAMat[TraceSeqs[1 + x*(UBTS+1)] + TraceSeqs[1 + Y*(UBTS+1)]*(UBFA+1)] = tFAMat[x + Y*(UBTFA + 1)];
			FAMat[TraceSeqs[1 + Y*(UBTS + 1)] + TraceSeqs[1 + x*(UBTS + 1)] * (UBFA + 1)] = tFAMat[x + Y*(UBTFA + 1)];
       }
    }       
        
	return(1);
}

	double MyMathFuncs::FastBootDistIP(int df, int reps, int nextno, int lenseq, float *dx, float *vx, int *weightmod, short int *seqnum, float *distance)
	{
		//weightmod - reps,len
		//distance -reps,nextno,nextno

		float th2, th3;
		float th1, dst;
		
		int *validx, *diffsx;

		int h, r, a, s1, x, y, z, off1, off0, off3;
		int xoff, yoff;
		
		r = reps + 1;
		
		off0 = reps + 1;
		off1 = (nextno + 1)*(reps + 1);
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
		
#pragma omp parallel for private (xoff, yoff, z, s1, a, off3, th1, th2, th3, dst, validx, diffsx, y, x)
		for (x = 0; x < nextno; x++) {
			xoff = x*(lenseq + 1);
			validx = (int*)calloc(r, sizeof(int));
			diffsx = (int*)calloc(r, sizeof(int));
			for (y = x + 1; y <= nextno; y++) {
				yoff = y*(lenseq + 1);
				
				for (a = df; a <= reps; a++)
					diffsx[a] = 0.0;

				for (a = df; a <= reps; a++)
					validx[a] = 0.0;

				for (z = 1; z <= lenseq; z++) {
					s1 = seqnum[z + xoff];

					if (s1 != 46) {

						if (s1 == seqnum[z + yoff]) {
							off3 = (off0)*(z - 1);
							for (a = df; a <= reps; a++)
								validx[a] = validx[a] + weightmod[a + off3];
						}
						else if (seqnum[z + yoff] != 46) {
							off3 = (off0)*(z - 1);
							for (a = df; a <= reps; a++) {
								validx[a] = validx[a] + weightmod[a + off3];
								diffsx[a] = diffsx[a] + weightmod[a + off3];
							}
						}
					}

				}
				for (a = df; a <= reps; a++) {
					//h=diffsx[a];
					if (validx[a] > 0) {
						th1 = (float)(diffsx[a]);
						th2 = (float)(validx[a]);
						dst = (th2 - th1) / th2;
						//return(dst);
						if (dst > 0.25) {
							th2 = (float)((4.0 * dst - 1.0) / 3.0);
							th3 = (float)(log(th2));
							dst = (float)(-0.75*th3);
						}
						else
							dst = 10.0;

					}

					else
						dst = 10.0;


					distance[a + x*off0 + y*off1] = (float)(dst);
					distance[a + y*off0 + x*off1] = (float)(dst);

				}
				
			}
			free(validx);
			free(diffsx);
		}
		omp_set_num_threads(2);
		return(1);
	}

	double MyMathFuncs::FastBootDistIP7(int df, int reps, int nextno, int lenseq, float *dx, float *vx, int *weightmod, short int *seqnum, float *distance)
	{
		//weightmod - reps,len
		//distance -reps,nextno,nextno

		float th2, th3;
		float th1, dst;

		

		int h, r, a, s1, x, y, z, off1, off0, off3;
		int xoff, yoff;

		r = reps + 1;

		off0 = reps + 1;
		off1 = (nextno + 1)*(reps + 1);


		int *xy;
		int c, xx, yy;
		xy = (int*)calloc((nextno + 1)*(nextno)+10, sizeof(int));
		c = 0;
		for (xx = 0; xx < nextno; xx++) {
			for (yy = xx + 1; yy <= nextno; yy++) {
				xy[c] = xx;
				xy[c + 1] = yy;
				c = c + 2;
			}
		}
		c = c - 2;
		
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel
		{
			int *validx, *diffsx;
			int d;
			validx = (int*)calloc(r, sizeof(int));
			diffsx = (int*)calloc(r, sizeof(int));
	#pragma omp for private (xoff, yoff, z, s1, a, off3, th1, th2, th3, dst, y, x)
			for (d = 0; d <= c/2; d++) {
				x = xy[d*2];
				y = xy[d*2 + 1];
				xoff = x*(lenseq + 1);
			

				yoff = y*(lenseq + 1);

				for (a = df; a <= reps; a++)
					diffsx[a] = 0.0;

				for (a = df; a <= reps; a++)
					validx[a] = 0.0;

				for (z = 1; z <= lenseq; z++) {
					s1 = seqnum[z + xoff];

					if (s1 != 46) {

						if (s1 == seqnum[z + yoff]) {
							off3 = (off0)*(z - 1);
							for (a = df; a <= reps; a++)
								validx[a] = validx[a] + weightmod[a + off3];
						}
						else if (seqnum[z + yoff] != 46) {
							off3 = (off0)*(z - 1);
							for (a = df; a <= reps; a++) {
								validx[a] = validx[a] + weightmod[a + off3];
								diffsx[a] = diffsx[a] + weightmod[a + off3];
							}
						}
					}

				}
				for (a = df; a <= reps; a++) {
					//h=diffsx[a];
					if (validx[a] > 0) {
						th1 = (float)(diffsx[a]);
						th2 = (float)(validx[a]);
						dst = (th2 - th1) / th2;
						//return(dst);
						if (dst > 0.25) {
							th2 = (float)((4.0 * dst - 1.0) / 3.0);
							th3 = (float)(log(th2));
							dst = (float)(-0.75*th3);
						}
						else
							dst = 10.0;

					}

					else
						dst = 10.0;


					distance[a + x*off0 + y*off1] = (float)(dst);
					distance[a + y*off0 + x*off1] = (float)(dst);

				}

			}
			free(validx);
			free(diffsx);
		}
		free (xy);
		omp_set_num_threads(2);
		return(1);
	}
	double MyMathFuncs::FastBootDistIP6(int dfx, int repsx, int nextnox, int lenseqx, float *dx, float *vx, int UBWM1, int UBWM2,int *wm, int UBSN1, int UBSN2,short int *sn, int UBD1, int UBD2,float *dist)
	{
		//weightmod - reps,len
		//distance -reps,nextno,nextno
		int *xy;
		int c, xx, yy;
		xy = (int*)calloc((nextnox+1)*(nextnox)+nextnox, sizeof(int));
		c = 0;
		for (xx = 0; xx < nextnox; xx++) {
			for (yy = xx + 1; yy <= nextnox; yy++) {
				xy[c] = xx;
				xy[c + 1] = yy;
				c = c + 2;
			}
		}
		c = c - 2;
		c = c / 2;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2;
		if (procs < 2)
			procs = 2;
		omp_set_num_threads(procs);
#pragma omp parallel 
		{
			float th2, th3;
			float th1, dst;
			int h, r, a, s1, x, y, z, off1, off0, off3, df, reps, nextno, lenseq;
			int xoff, yoff, dummy;
			short int *seqnum;
			int *validx, *diffsx;
			int *weightmod;
			float *distance;
			reps = repsx;
			df = dfx;
			nextno = nextnox;
			lenseq = lenseqx;
			r = reps + 1;
			off0 = reps + 1;
			off1 = (nextno + 1)*(reps + 1);
			
			validx = (int*)calloc(r+1, sizeof(int));
			diffsx = (int*)calloc(r+1, sizeof(int));
			//distance = (float*)calloc((UBD1+1)*(UBD2+1)*(UBD2 + 1), sizeof(float));
			seqnum = (short int*)calloc((UBSN1+1)*(UBSN2 + 1), sizeof(short int));
			weightmod = (int*)calloc((UBWM1+1)*(UBWM2 + 1), sizeof(int));
			memcpy (seqnum, sn, (UBSN1 + 1)*(UBSN2 + 1)*sizeof(short int));
			memcpy(weightmod, wm, (UBWM1+1)*(UBWM2 + 1) * sizeof(int));
			int d;
#pragma omp for private(d)
			for (d = 0; d <= c; d++) {
				x = xy[d*2];
				y = xy[d*2 + 1];
				xoff = x*(UBSN1 + 1);
				yoff = y*(UBSN1 + 1);

					for (a = df; a <= reps; a++)
						diffsx[a] = 0.0;

					for (a = df; a <= reps; a++)
						validx[a] = 0.0;

					for (z = 1; z <= lenseq; z++) {
						s1 = seqnum[z + xoff];

						if (s1 != 46) {

							if (s1 == seqnum[z + yoff]) {
								off3 = (off0)*(z - 1);
								for (a = df; a <= reps; a++)
									validx[a] = validx[a] + weightmod[a + off3];
							}
							else if (seqnum[z + yoff] != 46) {
								off3 = (off0)*(z - 1);
								for (a = df; a <= reps; a++) {
									validx[a] = validx[a] + weightmod[a + off3];
									diffsx[a] = diffsx[a] + weightmod[a + off3];
								}
							}
						}

					}
					for (a = df; a <= reps; a++) {
						//h=diffsx[a];
						if (validx[a] > 0) {
							th1 = (float)(diffsx[a]);
							th2 = (float)(validx[a]);
							dst = (th2 - th1) / th2;
							//return(dst);
							if (dst > 0.25) {
								th2 = (float)((4.0 * dst - 1.0) / 3.0);
								th3 = (float)(log(th2));
								dst = (float)(-0.75*th3);
							}
							else
								dst = 10.0;

						}

						else
							dst = 10.0;

#pragma omp critical
						{
							//distance[a + x*off0 + y*off1] = (float)(dst);
							//distance[a + y*off0 + x*off1] = (float)(dst);
							dist[a + x*off0 + y*off1] = (float)dst;
							dist[a + y*off0 + x*off1] = (float)dst;
						}

					}

				
				
			}
			free(validx);
			free(diffsx);
			free(seqnum);
			free(weightmod);
			/*for (x = 0; x < nextno; x++) {
				xoff = x*(lenseq + 1);
				for (y = x + 1; y <= nextno; y++) {
					yoff = y*(lenseq + 1);
					for (a = df; a <= reps; a++) {
						dst = distance[a + x*off0 + y*off1];
						if (dst != 0) {
							dist[a + x*off0 + y*off1] = dst;
							dist[a + y*off0 + x*off1] = dst;
						}

					}
				}
			}*/
			//free(distance);
		}
		
		omp_set_num_threads(2);
		
		free(xy);
		return(1);
	}


	double MyMathFuncs::FastBootDistIP5(int df, int reps, int nextno, int lenseq, int *weightmod, short int *seqnum, float *distance)
	{
		//weightmod - reps,len
		//distance -reps,nextno,nextno

		float th2, th3;
		float th1, dst;

		

		int h, r, a, s1, x, y, z, off1, off0, off3,n,ls,lsp1,apo;
		int xoff, yoff;
		int *validx, *diffsx;
		r = reps + 1;

		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel private (xoff, yoff, z, s1, a, off3, th1, th2, th3, dst, validx, diffsx, y, x, off0, off1,n,ls,apo)
		{
		
		n = nextno;
		ls = lenseq;
		lsp1 = ls + 1;
		off0 = reps + 1;
		off1 = (n + 1)*(reps + 1);
		
		validx = (int*)calloc(r, sizeof(int));
		diffsx = (int*)calloc(r, sizeof(int));
		
#pragma omp for
			for (x = 0; x < n; x++) {
				xoff = x*lsp1;
				
				for (y = x + 1; y <= n; y++) {
					yoff = y*lsp1;

					for (a = df; a <= reps; a++)
						diffsx[a] = 0.0;

					for (a = df; a <= reps; a++)
						validx[a] = 0.0;

					for (z = 1; z <= ls; z++) {
						//s1 = seqnum[z + xoff];

						if (seqnum[z + xoff] != 46) {

							if (seqnum[z + xoff] == seqnum[z + yoff]) {
								off3 = (off0)*(z - 1);
								for (a = df; a <= reps; a++)
									validx[a] = validx[a] + weightmod[a + off3];
							}
							else if (seqnum[z + yoff] != 46) {
								off3 = (off0)*(z - 1);
								for (a = df; a <= reps; a++) {
									apo = a + off3;
									validx[a] = validx[a] + weightmod[apo];
									diffsx[a] = diffsx[a] + weightmod[apo];
								}
							}
						}

					}
					for (a = df; a <= reps; a++) {
						//h=diffsx[a];
						if (validx[a] > 0) {
							th1 = (float)(diffsx[a]);
							th2 = (float)(validx[a]);
							dst = (th2 - th1) / th2;
							//return(dst);
							if (dst > 0.25) {
								th2 = (float)((4.0 * dst - 1.0) / 3.0);
								th3 = (float)(log(th2));
								dst = (float)(-0.75*th3);
							}
							else
								dst = 10.0;

						}

						else
							dst = 10.0;


						distance[a + x*off0 + y*off1] = (float)(dst);
						distance[a + y*off0 + x*off1] = (float)(dst);

					}

				}
				
			}
			free(validx);
			free(diffsx);
		}
		omp_set_num_threads(2);
		return(1);
	}


	int MyMathFuncs::FindNewX(int WinPPY, int WinPP, int Seq3, int Nextno, int *RNum) {
		int X;

		for (X = 0; X <= WinPPY; X++) {
			if (Seq3 == Nextno - RNum[WinPP] + X)
				break;
		}

		return(X);
	}

	/*double MyMathFuncs::FastBootDistIP(int df, int reps, int nextno, int lenseq, float *diffsx, float *validx, int *weightmod, short int *seqnum, float *distance)
	{
		//weightmod - reps,len
		//distance -reps,nextno,nextno

		float th2, th3;

		float th1, dst;
		int nn, r, sl, sublen, df2;
		int a, s1, s2, x, y, z, off1, off0, off3, off4, b, c;
		int* dx;
		int* vx;
		nn = nextno;
		r = reps+1;
		df2 = df;
		int xoff, yoff;
		//sublen = IdenticalF[lenseq];
		s2 = 66;
		off0 = reps + 1;
		off1 = (nextno + 1)*(reps + 1);
		
//#pragma omp parallel for private (xoff, yoff, b, c, z, s1, a, off3, th1, th2, dst, dx, vx, y, x)
		for (x = 0; x < nextno; x++) {
			xoff = x*(lenseq + 1);
				for (y = x + 1; y <= nextno; y++) {
					
					yoff = y*(lenseq + 1);
					
					dx = (int*)calloc(r, sizeof(int));
					vx = (int*)calloc(r, sizeof(int));
					
					for (b = df2; b <= r; b++)
						dx[b] = 0;

					for (c = df2; c <= r; c++)
						vx[c] = 0;

					for (z = 1; z <= lenseq; z++) {
					//for (z = 1; z <= sublen; z++) {
						//off4 = IdenticalR[z];
						
						s1 = seqnum[z + xoff];

						if (s1 != 46) {

							if (s1 == seqnum[z + yoff]) {
								off3 = off0*(z-1);
								for (a = df2; a <= r; a++)
									vx[a] = vx[a] + weightmod[a + off3];
							}
							else if (seqnum[z + yoff] != 46) {
								off3 = off0*(z-1);
								for (a = df2; a <= r; a++) {
									vx[a] = vx[a] + weightmod[a + off3];


									//for (a = df; a <= reps; a++)	
									dx[a] = dx[a] + weightmod[a + off3];
								}
							}
						}

					}
					for (a = df2; a <= r; a++) {
						if (vx[a] > 0) {
							th1 = (float)(dx[a]);
							th2 = (float)(vx[a]);
							dst = (th2 - th1) / th2;
							//return(dst);
							if (dst > 0.25) {
								th2 = (float)((4.0 * dst - 1.0) / 3.0);
								th3 = (float)(log(th2));
								dst = (float)(-0.75*th3);
							}
							else
								dst = 10.0;

						}

						else
							dst = 10.0;


						distance[a + x*off0 + y*off1] = (float)(dst);
						distance[a + y*off0 + x*off1] = (float)(dst);
						
					}
					free(vx);
					free(dx);
					
				}
				
			}
		
		

		return(1);
	}*/

	int MyMathFuncs::FastBootDistIP4(int df, int reps, int nextno, int lenseq, float *diffsx, float *validx, int *weightmod, short int *seqnum, float *distance, unsigned char *fd, unsigned char *fv)
	{
		//weightmod - reps,len
		//distance -reps,nextno,nextno

		float th2, th3;

		float th1, dst;
		int nn, r, sl, sublen, df2;
		int a, s1, s2, x, y, z, off1, off0, off3, off4, b, c, d, v, os,wm;
		int* dx;
		int* vx;
		nn = nextno;
		r = reps + 1;
		df2 = df;
		int xoff, yoff;
		//sublen = IdenticalF[lenseq];
		s2 = 66;
		off0 = reps + 1;
		off1 = (nextno + 1)*(reps + 1);
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel private (xoff, yoff, b, c, z,  a, off3, th1, th2, dst, dx, vx, y, d, v,os,wm, x)
		{
			dx = (int*)calloc(r, sizeof(int));
			vx = (int*)calloc(r, sizeof(int));
#pragma omp for 
			for (x = 0; x < nextno; x++) {
				xoff = x*(lenseq + 1);



				for (y = x + 1; y <= nextno; y++) {

					yoff = y*(lenseq + 1);



					for (b = df2; b <= r; b++)
						dx[b] = 0;

					for (c = df2; c <= r; c++)
						vx[c] = 0;

					for (z = 1; z <= lenseq; z++) {
						off3 = off0*(z - 1);
						os = seqnum[z + xoff] + seqnum[z + yoff] * 86;
						d = (int)fd[os];
						v = (int)fv[os];
						for (a = df2; a <= r; a++){
							wm = weightmod[a + off3];
							vx[a] += d*wm;
							dx[a] += v*wm; 
						
						}

					}
					for (a = df2; a <= r; a++) {
						if (vx[a] > 0) {
							th1 = (float)(dx[a]);
							th2 = (float)(vx[a]);
							dst = (th2 - th1) / th2;
							//return(dst);
							if (dst > 0.25) {
								th2 = (float)((4.0 * dst - 1.0) / 3.0);
								th3 = (float)(log(th2));
								dst = (float)(-0.75*th3);
							}
							else
								dst = 10.0;

						}

						else
							dst = 10.0;


						distance[a + x*off0 + y*off1] = (float)(dst);
						distance[a + y*off0 + x*off1] = (float)(dst);

					}

				}
			}
			free(vx);
			free(dx);
		}
		omp_set_num_threads(2);
		return(1);
	}

	double MyMathFuncs::FastBootDistP(int df, int reps, int nextno, int lenseq, float *diffsx, float *validx, int *weightmod, short int *seqnum, float *distance)
	{
		//weightmod - reps,len
		//distance -reps,nextno,nextno

		float th2, th3;
		float th1, dst;

		int a, s1, s2, x, y, z, off1, off0, off3;
		int xoff, yoff;

		s2 = 66;
		off0 = reps + 1;
		off1 = (nextno + 1)*(reps + 1);
		for (x = 0; x < nextno; x++) {
			xoff = x*(lenseq + 1);
			for (y = x + 1; y <= nextno; y++) {
				yoff = y*(lenseq + 1);
				for (a = df; a <= reps; a++)
					diffsx[a] = 0.0;

				for (a = df; a <= reps; a++)
					validx[a] = 0.0;

				for (z = 1; z <= lenseq; z++) {
					s1 = seqnum[z + xoff];

					if (s1 != 46) {

						if (s1 == seqnum[z + yoff]) {
							off3 = (off0)*(z - 1);
							for (a = df; a <= reps; a++)
								validx[a] = validx[a] + weightmod[a + off3];
						}
						else if (seqnum[z + yoff] != 46) {
							off3 = (off0)*(z - 1);
							for (a = df; a <= reps; a++) {
								validx[a] = validx[a] + weightmod[a + off3];


								//for (a = df; a <= reps; a++)	
								diffsx[a] = diffsx[a] + weightmod[a + off3];
							}
						}
					}

				}
				for (a = df; a <= reps; a++) {
					if (validx[a] > 0) {
						th1 = (float)(diffsx[a]);
						th2 = (float)(validx[a]);
						dst = (th2 - th1) / th2;
						//return(dst);
						if (dst > 0.25) {
							th2 = (float)((4.0 * dst - 1.0) / 3.0);
							th3 = (float)(log(th2));
							dst = (float)(-0.75*th3);
						}
						else
							dst = 10.0;

					}

					else
						dst = 10.0;


					distance[a + x*off0 + y*off1] = (float)(dst);
					distance[a + y*off0 + x*off1] = (float)(dst);

				}
			}

		}

		return(1);
	}
	


	double MyMathFuncs::CalcChiVals4P(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, int *BanWins)
	{

		double ChiH, MChi;
		float E, A, B, C, D;
		int  X, LO, FO, SO, bp1;
		//float ;
		FO = UBWS + 1;// LenSeq + (HWindowWidth)* 2 + 1;
		SO = LenSeq + 1;
		LO = HWindowWidth;
		ChiH = 0;
		MChi = 0;
		for (X = 0; X < LenXoverSeq; X++) {
			bp1 = X - HWindowWidth;
			if (bp1 < 1)
				bp1 = bp1 + LenXoverSeq;
			if (BanWins[X] == 0 && BanWins[bp1] == 0) {
				/*for (Y = 0; Y <= 2;Y++){
				A = *(WinScores + X + FO*Y);
				C = *(WinScores + X + LO + FO*Y);
				if (A - C > criticaldiff || A-C < -criticaldiff){
				B = HWindowWidth - A;
				D = HWindowWidth - C;

				if (A + C > 0 && B + D > 0){
				E = A*D - B*C;
				ChiH = E*E*2/(HWindowWidth*(A + C)*(B + D));
				*(ChiVals + X + SO*Y) = ChiH;
				if (ChiH > MChi)
				MChi = ChiH;
				}
				else
				*(ChiVals + X + SO*Y) = 0;
				}
				else
				*(ChiVals + X + SO*Y) = 0;//abs(A-C)/HWindowWidth;
				}*/

				//0

				A = (float)(*(WinScores + X));
				C = (float)(*(WinScores + X + LO));
				if (A - C > criticaldiff || A - C < -criticaldiff) {
					B = HWindowWidth - A;
					D = HWindowWidth - C;

					if (A + C > 0 && B + D > 0) {
						E = A*D - B*C;
						ChiH = (double)(E*E * 2) / (double)(HWindowWidth*(A + C)*(B + D));
						ChiVals [X] = ChiH;
						if (ChiH > MChi)
							MChi = ChiH;
					}
					else
						ChiVals [X] = 0;
				}
				else
					ChiVals [X] = 0;//abs(A-C)/HWindowWidth;

									   //1
				A = (float)(*(WinScores + X + FO));
				C = (float)(*(WinScores + X + LO + FO));
				if (A - C > criticaldiff || A - C < -criticaldiff) {
					B = HWindowWidth - A;
					D = HWindowWidth - C;

					if (A + C > 0 && B + D > 0) {
						E = A*D - B*C;

						ChiH = (double)(E*E * 2) / (double)(HWindowWidth*(A + C)*(B + D));
						*(ChiVals + X + SO) = ChiH;
						if (ChiH > MChi)
							MChi = ChiH;
					}
					else
						*(ChiVals + X + SO) = 0;
				}
				else
					*(ChiVals + X + SO) = 0;//abs(A-C)/HWindowWidth;

											//2
				A = (float)(*(WinScores + X + FO * 2));
				C = (float)(*(WinScores + X + LO + FO * 2));
				if (A - C > criticaldiff || A - C < -criticaldiff) {
					B = HWindowWidth - A;
					D = HWindowWidth - C;

					if (A + C > 0 && B + D > 0) {
						E = A*D - B*C;
						ChiH = (double)(E*E * 2) / (double)(HWindowWidth*(A + C)*(B + D));
						*(ChiVals + X + SO * 2) = ChiH;
						if (ChiH > MChi)
							MChi = ChiH;
					}
					else
						*(ChiVals + X + SO * 2) = 0;
				}
				else
					*(ChiVals + X + SO * 2) = 0;//abs(A-C)/HWindowWidth;
			}
			else {

				*(ChiVals + X) = 0;
				*(ChiVals + X + SO) = 0;
				*(ChiVals + X + SO * 2) = 0;
			}
		}
		return (MChi);
	}
	double MyMathFuncs::CalcChiVals4P2(int UBCT, int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, int *BanWins, float *chitable)
	{

		double ChiH, MChi;
		int A, B, C, D;
		int  X, LO, FO, SO, bp1, os1, os2, S1, F1;
		//float ;
		FO = UBWS + 1;// LenSeq + (HWindowWidth)* 2 + 1;
		SO = LenSeq + 1;
		S1 = SO * 2;
		F1 = FO * 2;
		LO = HWindowWidth;
		ChiH = 0;
		MChi = 0;

		os1 = UBCT + 1;
		os2 = os1*os1*LO;
		for (X = 0; X < LenXoverSeq; X++) {
			bp1 = X - HWindowWidth;
			if (bp1 < 1)
				bp1 = bp1 + LenXoverSeq;
			if (BanWins[X] == 0 && BanWins[bp1] == 0) {

				
				ChiVals[X] = chitable[WinScores[X] + WinScores[X + LO]*os1 + os2];
				ChiVals[X + SO] = chitable[WinScores[X + FO] + WinScores[X + LO + FO] * os1 + os2];
				ChiVals[X + S1] = chitable[WinScores[X + F1] + WinScores[X + LO + F1] * os1 + os2];
			}
			else {

				*(ChiVals + X) = 0;
				*(ChiVals + X + SO) = 0;
				*(ChiVals + X + S1) = 0;
			}
		}

		for (X = 0; X < LenXoverSeq; X++) {
			if (ChiVals[X] > MChi)
				MChi = ChiVals[X];
		}
		for (X = 0; X < LenXoverSeq; X++) {
			if (ChiVals[X+SO] > MChi)
				MChi = ChiVals[X+SO];
		}
		for (X = 0; X < LenXoverSeq; X++) {
			if (ChiVals[X + S1] > MChi)
				MChi = ChiVals[X + S1];
		}
		return (MChi);
	}

	double MyMathFuncs::CalcChiVals4P3(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, int *BanWins, float *chitable)
	{

		double ChiH, MChi;
		int A, B, C, D;
		int  X, LO, FO, SO, bp1, os1, os2, S1, F1;
		//float ;
		FO = UBWS + 1;// LenSeq + (HWindowWidth)* 2 + 1;
		SO = LenSeq + 1;
		S1 = SO * 2;
		F1 = FO * 2;
		LO = HWindowWidth;
		ChiH = 0;
		MChi = 0;

		os1 = HWindowWidth + 1;
		
		for (X = 0; X < LenXoverSeq; X++) {
			bp1 = X - HWindowWidth;
			if (bp1 < 1)
				bp1 = bp1 + LenXoverSeq;
			if (BanWins[X] == 0 && BanWins[bp1] == 0) {


				ChiVals[X] = chitable[WinScores[X] + WinScores[X + LO] * os1];
				ChiVals[X + SO] = chitable[WinScores[X + FO] + WinScores[X + LO + FO] * os1];
				ChiVals[X + S1] = chitable[WinScores[X + F1] + WinScores[X + LO + F1] * os1];
			}
			else {

				*(ChiVals + X) = 0;
				*(ChiVals + X + SO) = 0;
				*(ChiVals + X + S1) = 0;
			}
		}

		for (X = 0; X < LenXoverSeq; X++) {
			if (ChiVals[X] > MChi)
				MChi = ChiVals[X];
		}
		for (X = 0; X < LenXoverSeq; X++) {
			if (ChiVals[X + SO] > MChi)
				MChi = ChiVals[X + SO];
		}
		for (X = 0; X < LenXoverSeq; X++) {
			if (ChiVals[X + S1] > MChi)
				MChi = ChiVals[X + S1];
		}
		return (MChi);
	}

	int MyMathFuncs::GrowMChiWin2P(int LO, int RO, int LenXoverSeq, int HWindowWidth, int TWin, int MaxY, int LS, int a, int c, int MaxFailCount, double *MPV, int *WinWin, double *MChi, int *TopL, int *TopR, int *TopLO, int *TopRO, unsigned char *Scores, unsigned char *mdmap) {
		//scores LS,2 
		int off1, FailCount;
		float A, C, E, B, D;
		double mchi, hMChi, tpv, lxdh, lx, h, t, htpv;
		off1 = MaxY*(LS + 1);
		FailCount = 0;

		lx = (double)(LenXoverSeq);
		h = (double)(HWindowWidth);
		lxdh = lx / h;
		A = (float)(a);
		C = (float)(c);
		mchi = *MChi;
		/*t = (double)(TWin);
		h = (double)(LenXoverSeq);
		h = h / t;*/
		//TW2 = TWin*2;
//#pragma loop(no_vector)
		htpv = mchi;// / lxdh;
		while (FailCount <= MaxFailCount) {
			A = A + (float)(Scores[LO + off1]);
			C = C + (float)(Scores[RO + off1]);
			B = (float)(TWin) - A;
			D = (float)(TWin) - C;
			if (B + D > 0) {
				E = A*D - B*C;
				E = E*E;
				E = E * 2;
				hMChi = (double)(E) / (double)((TWin*(A + C)*(B + D)));
				/*
				hMChi = (A * D - B * C);
				hMChi =	hMChi * hMChi;
				hMChi = hMChi*TWin*2;
				hMChi = hMChi / (A + B);
				hMChi = hMChi / (C + D);
				hMChi = hMChi / (A + C);
				hMChi = hMChi / (B + D);*/
				//tpv = hMChi;
				
				if (hMChi >= htpv) {
					htpv = hMChi;
					*TopLO = LO;
					*TopRO = RO;
					mchi = hMChi;
					*WinWin = TWin;
					//*MPV = tpv;
					FailCount = 0;
					*TopL = (int)(A);
					*TopR = (int)(C);
				}
				else {
					if (TWin == HWindowWidth && *TopL == 0 && *TopR == 0) {
						*TopL = (int)(A);
						*TopR = (int)(C);
					}
					FailCount++;
					if (FailCount > MaxFailCount) {
						a = (int)(A);
						c = (int)(C);
						*MChi = mchi;
						return(1);
					}
				}
			}
			else {
				FailCount++;
				if (FailCount > MaxFailCount) {
					a = (int)(A);
					c = (int)(C);
					*MChi = mchi;
					return(1);
				}
				hMChi = 0;
			}
			if (mdmap[RO] == 1 || mdmap[LO] == 1)
				break;
			RO++;
			LO--;

			if (mdmap[RO] == 1 || mdmap[LO] == 1)
				break;

			if (LO < 1)
				LO = LenXoverSeq;
			if (RO > LenXoverSeq * 2)
				RO = RO - LenXoverSeq * 2;
			if (RO > LenXoverSeq)
				RO = RO - LenXoverSeq;

			TWin++;
		}
		*MChi = mchi;
		a = (int)(A);
		c = (int)(C);
		return(1);
	}
	

	int MyMathFuncs::GrowMChiWin2P2(int MaxABWin, int LO, int RO, int LenXoverSeq, int HWindowWidth, int TWin, int MaxY, int LS, int a, int c, int MaxFailCount, double *MPV, int *WinWin, double *MChi, int *TopL, int *TopR, int *TopLO, int *TopRO, unsigned char *Scores, unsigned char *mdmap, float *chitable, int *chimap) {
		//scores LS,2 
		int off1, FailCount,aa,cc;
		float A, C, E, B, D;
		double mchi, hMChi, tpv, lxdh, lx, h, t, htpv;
		off1 = MaxY*(LS + 1);
		FailCount = 0;

		lx = (double)(LenXoverSeq);
		h = (double)(HWindowWidth);
		lxdh = lx / h;
		A = (float)(a);
		C = (float)(c);
		mchi = *MChi;
		/*t = (double)(TWin);
		h = (double)(LenXoverSeq);
		h = h / t;*/
		//TW2 = TWin*2;
		//#pragma loop(no_vector)
		htpv = mchi;
		while (FailCount <= MaxFailCount) {
			A = A + (float)(Scores[LO + off1]);
			C = C + (float)(Scores[RO + off1]);
			B = (float)(TWin)-A;
			D = (float)(TWin)-C;
			if (B + D > 0) {
				if (TWin <= MaxABWin) {
					aa = A;
					cc = C;
					//hx = chimap[TWin];
					hMChi = chitable[chimap[TWin] + aa + cc*(TWin + 1)];
					//tpv = ChiPVal(hMChi) * (h / t);
					if (hMChi >= htpv) {
						htpv = hMChi;
						*TopLO = LO;
						*TopRO = RO;
						mchi = hMChi;
						*WinWin = TWin;
						//*MPV = tpv;
						FailCount = 0;
						*TopL = A;
						*TopR = C;
					}
					else {
						if (TWin == HWindowWidth && *TopL == 0 && *TopR == 0) {
							*TopL = (int)(A);
							*TopR = (int)(C);
						}
						FailCount++;
						if (FailCount > MaxFailCount) {
							a = (int)(A);
							c = (int)(C);
							*MChi = mchi;
							return(1);
						}
					}
				}
				else {
					
					//if (A + C > 0 && B + D > 0) {
						E = A*D - B*C;
						E = E*E;
						E = E * 2;
						hMChi = (double)(E) / (double)((TWin*(A + C)*(B + D)));
						/*
						hMChi = (A * D - B * C);
						hMChi =	hMChi * hMChi;
						hMChi = hMChi*TWin*2;
						hMChi = hMChi / (A + B);
						hMChi = hMChi / (C + D);
						hMChi = hMChi / (A + C);
						hMChi = hMChi / (B + D);*/
						//tpv = hMChi;

						if (hMChi >= htpv) {
							htpv = hMChi;
							*TopLO = LO;
							*TopRO = RO;
							mchi = hMChi;
							*WinWin = TWin;
							//*MPV = tpv;
							FailCount = 0;
							*TopL = (int)(A);
							*TopR = (int)(C);
						}
						else {
							if (TWin == HWindowWidth && *TopL == 0 && *TopR == 0) {
								*TopL = (int)(A);
								*TopR = (int)(C);
							}
							FailCount++;
							if (FailCount > MaxFailCount) {
								a = (int)(A);
								c = (int)(C);
								*MChi = mchi;
								return(1);
							}
						}
					//}
					/*else {
						FailCount++;
						if (FailCount > MaxFailCount) {
							a = (int)(A);
							c = (int)(C);
							*MChi = mchi;
							return(1);
						}
						hMChi = 0;
					}*/
				}
			}
			else {
				FailCount++;
				if (FailCount > MaxFailCount) {
					a = (int)(A);
					c = (int)(C);
					*MChi = mchi;
					return(1);
				}
				hMChi = 0;
			}
			if (mdmap[RO] == 1 || mdmap[LO] == 1)
				break;
			RO++;
			LO--;

			if (mdmap[RO] == 1 || mdmap[LO] == 1)
				break;

			if (LO < 1)
				LO = LenXoverSeq;
			if (RO > LenXoverSeq * 2)
				RO = RO - LenXoverSeq * 2;
			if (RO > LenXoverSeq)
				RO = RO - LenXoverSeq;

			TWin++;
		}
		*MChi = mchi;
		a = (int)(A);
		c = (int)(C);
		return(1);
	}

	int  MyMathFuncs::CleanChiVals(int LenXoverSeq, int LenSeq, double *ChiVals) {
		int SO, SO2, X;
		SO = LenSeq + 1;
		SO2 = SO * 2;
		for (X = 0; X < LenXoverSeq; X++) 
				ChiVals[X] = 0;

		for (X = 0; X < LenXoverSeq; X++) 
			ChiVals[X + SO] = 0;

		for (X = 0; X < LenXoverSeq; X++)
			ChiVals[X + SO2] = 0;
			
		return(1);
	}
	int  MyMathFuncs::CleanChiVals2(int LenXoverSeq, int LenSeq, double *ChiVals) {
		int X;
		
		for (X = 0; X < LenXoverSeq; X++)
			ChiVals[X] = 0;

		return(1);
	}
	double  MyMathFuncs::CalcChiValsP(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals)
	{

		double ChiH, MChi;
		//double ChiH;
		float A, B, C, D, E;
		int X, LO, FO, SO, Y;
		FO = UBWS + 1;// LenSeq + (HWindowWidth)* 2 + 1;
		SO = LenSeq + 1;
		LO = HWindowWidth;
		//ChiH = 0;
		MChi = 0;
		
//#pragma omp parallel for private (Y, A, C, B, D, E, ChiH) 
		for (X = 0; X < LenXoverSeq; X++) {
			for (Y = 0; Y <= 2;Y++){
				A = (float)(*(WinScores + X + FO*Y));
				C = (float)(*(WinScores + X + LO + FO*Y));
				if (A - C > criticaldiff || A-C < -criticaldiff){
					B = HWindowWidth - A;
					D = HWindowWidth - C;

					if (A + C > 0 && B + D > 0){
						E = A*D - B*C;
						ChiH = (double)E*E*2/(double)(HWindowWidth*(A + C)*(B + D));
						*(ChiVals + X + SO*Y) = ChiH;
						if (ChiH > MChi)
							MChi = ChiH;
					}
					else
						*(ChiVals + X + SO*Y) = 0;
				}
				else
				*(ChiVals + X + SO*Y) = 0;//abs(A-C)/HWindowWidth;
			}

			/*
			//0

			A = *(WinScores + X);
			C = *(WinScores + X + LO);
			if (A - C >= criticaldiff || A - C <= -criticaldiff) {
				B = HWindowWidth - A;
				D = HWindowWidth - C;

				if (A + C > 0 && B + D > 0) {
					E = A*D - B*C;
					ChiH = (double)(E*E *2) / (double)(HWindowWidth*(A + C)*(B + D));
					*(ChiVals + X) = ChiH;
					if (ChiH > MChi)
						MChi = ChiH;
				}
				else
					*(ChiVals + X) = 0;
			}
			else
				*(ChiVals + X) = 0;//abs(A-C)/HWindowWidth;

								   //1
			A = *(WinScores + X + FO);
			C = *(WinScores + X + LO + FO);
			if (A - C >= criticaldiff || A - C <= -criticaldiff) {
				B = HWindowWidth - A;
				D = HWindowWidth - C;

				if (A + C > 0 && B + D > 0) {
					E = A*D - B*C;

					ChiH = (double)(E*E * 2) / (double)(HWindowWidth*(A + C)*(B + D));
					*(ChiVals + X + SO) = ChiH;
					if (ChiH > MChi)
						MChi = ChiH;
				}
				else
					*(ChiVals + X + SO) = 0;
			}
			else
				*(ChiVals + X + SO) = 0;//abs(A-C)/HWindowWidth;

										//2
			A = *(WinScores + X + FO * 2);
			C = *(WinScores + X + LO + FO * 2);
			if (A - C >= criticaldiff || A - C <= -criticaldiff) {
				B = HWindowWidth - A;
				D = HWindowWidth - C;

				if (A + C > 0 && B + D > 0) {
					E = A*D - B*C;
					ChiH = (double)(E*E * 2) / (double)(HWindowWidth*(A + C)*(B + D));
					*(ChiVals + X + SO * 2) = ChiH;
					if (ChiH > MChi)
						MChi = ChiH;
				}
				else
					*(ChiVals + X + SO * 2) = 0;
			}
			else
				*(ChiVals + X + SO * 2) = 0;//abs(A-C)/HWindowWidth;


				*/
		}
		
		return (MChi);
	}



	int MyMathFuncs::FindOverlapP(int lenseq, int BPos2, int EPos2, int *RSize, int *OLSeq) {
		int Z;
		int OLSize = 0;
		if (BPos2 < EPos2) {
			RSize[1] = EPos2 - BPos2 + 1;
			for (Z = BPos2; Z <= EPos2; Z++)
				OLSize = OLSize + OLSeq[Z];

		}
		else {
			RSize[1] = EPos2 + lenseq - BPos2 + 1;
			for (Z = 1; Z <= EPos2; Z++)
				OLSize = OLSize + OLSeq[Z];

			for (Z = BPos2; Z <= lenseq; Z++)
				OLSize = OLSize + OLSeq[Z];

		}
		return (OLSize);
	}

	float MyMathFuncs::FindOverlapP2(int lenseq, int BPos2, int EPos2, int *RSize, int *OLSeq) {
		int Z;
		float OLSize = 0;
		float Valid = 0;
		if (BPos2 < EPos2) {
			RSize[1] = EPos2 - BPos2 + 1;
			for (Z = BPos2; Z <= EPos2; Z++) {
				if (OLSeq[Z] == 1) {
					OLSize = OLSize + 1;
					Valid = Valid + 1;
				}
				else if (OLSeq[Z] == 0)
					Valid = Valid + 1;

			}
		}
		else {
			RSize[1] = EPos2 + lenseq - BPos2 + 1;
			for (Z = 1; Z <= EPos2; Z++) {
				if (OLSeq[Z] == 1) {
					OLSize = OLSize + 1;
					Valid = Valid + 1;
				}
				else if (OLSeq[Z] == 0)
					Valid = Valid + 1;
			}
			for (Z = BPos2; Z <= lenseq; Z++) {
				if (OLSeq[Z] == 1) {
					OLSize = OLSize + 1;
					Valid = Valid + 1;
				}
				else if (OLSeq[Z] == 0)
					Valid = Valid + 1;
			}

		}
		if (Valid > 0)
			OLSize = OLSize / Valid;
		else
			OLSize = 0;
		return (OLSize);
	}

	int MyMathFuncs::TransferDistP(int NSeqs, int cr, int Reps, float *tFMat, float *DstMat) {

		int Z, A, off0, off1, off2, off3, off4, off5;
		off0 = Reps + 1;
		off1 = NSeqs + 1;
		off2 = off0*off1;

		for (Z = 0; Z < NSeqs; Z++) {
			off4 = Z*off2;
			off3 = Z*off1;

			for (A = Z + 1; A <= NSeqs; A++) {
				off5 = A*off1;
				tFMat[Z + off5] = DstMat[cr + A*off0 + off4];
				tFMat[A + off3] = tFMat[Z + off5];
			}
		}
		return (1);
	}

	int MyMathFuncs::TreeGroupsXP(int NextNo, char *THolder, int TLen, int NLen, char *TMatch, float *DLen) {

		int X, Y, Z, Cnt, SeqID, NCount, Miss, Hit, off1;
		char *TArray;
		char *DoneNode;

		TArray = (char*)calloc((NextNo +1)*(NextNo + 1), sizeof(char));
		DoneNode = (char*)calloc((NextNo + 1), sizeof(char));
		



		off1 = NextNo + 1;
		NCount = -1;
		for (X = 1; X <= TLen; X++) {
			if (THolder[X] == 40) {// Then 'ie (
				Cnt = 1;
				NCount++;
				Y = X + 1;
				while (Cnt > 0) {
					if (THolder[Y] == 40)// Then 'ie (
						Cnt++;
					else if (THolder[Y] == 41)// Then 'ie )
						Cnt--;
					else if (THolder[Y] == 83) {//Then 'ie S
						SeqID = 0;
						for (Z = 1; Z <= NLen; Z++)
							SeqID = SeqID + (int)(0.1 + (THolder[Y + Z] - 48) * pow(10, NLen - Z));

						TArray[NCount + SeqID*off1] = 1;
					}
					Y++;
				}
			}

		}

		X = 0;
		for (X = 0; X <= NextNo; X++) {
			for (Y = 0; Y <= NextNo; Y++) {
				if (DoneNode[Y] == 0) {
					Miss = 0;
					Hit = 0;
					for (Z = 0; Z <= NextNo; Z++) {
						if (TArray[Y + Z*off1] == TMatch[X + Z*off1])
							Hit++;
						else
							Miss++;

						if (Miss > 0 && Hit > 0)
							break;
					}
					if (Miss == 0 || Hit == 0) {
						DLen[X] = DLen[X] + 1;
						DoneNode[Y] = 1;

					}
				}
			}
		}


		free(TArray);
		free(DoneNode);
		return(1);
	}

	int MyMathFuncs::TreeReps(int NSeqs, int Reps,  int BSRndNumSeed, int NameLen, float *DstMat, int LTI, int *LTree, char *tMatch, float *DL) {
		float *tFMat, *DLen;
		char *FHolder;
		int Y, Z, Dummy;

		tFMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));
		

         for( Y = 1; Y <= Reps; Y++){
			
			DLen = (float*)calloc((NSeqs + 1), sizeof(float));
			FHolder = (char*)calloc(NSeqs * 40 * 2 +1, sizeof(char));
			Dummy = TransferDistP(NSeqs, Y, Reps, tFMat, DstMat);
			LTree[LTI] = rdp_nj_clearcut_compat(0, NSeqs, 1, 100, BSRndNumSeed, 1, NSeqs, tFMat, FHolder);
			Dummy = TreeGroupsXP(NSeqs, FHolder, LTree[LTI], NameLen, tMatch, DLen);
			for (Z = 0; Z <= NSeqs; Z++)
				DL[Z] = DL[Z] + DLen[Z];
			free(DLen);
			free(FHolder);
        }
		
		
		
		free(tFMat);
		

		return(1);
	}

	int MyMathFuncs::TreeRepsP(int NSeqs, int Reps, int BSRndNumSeed, int NameLen, float *DstMat, int LTI, int *LTree, char *tMatch, float *DL) {
		
#pragma omp parallel
		{
			float *tFMat, *DLen;
			char *FHolder;
			int Y, Z, Dummy;

			tFMat = (float*)calloc((NSeqs + 1)*(NSeqs + 1), sizeof(float));
			/*DLen = (float*)calloc((NSeqs + 1), sizeof(float));
			FHolder = (char*)calloc(NSeqs * 40 * 2 + 1, sizeof(char));*/
#pragma omp for
			for (Y = 1; Y <= Reps; Y++) {
				DLen = (float*)calloc((NSeqs + 1), sizeof(float));
				FHolder = (char*)calloc(NSeqs * 40 * 2 + 1, sizeof(char));
				Dummy = TransferDistP(NSeqs, Y, Reps, tFMat, DstMat);
				LTree[LTI] = rdp_nj_clearcut_compat(0, NSeqs, 1, 100, BSRndNumSeed, 1, NSeqs, tFMat, FHolder);
				Dummy = TreeGroupsXP(NSeqs, FHolder, LTree[LTI], NameLen, tMatch, DLen);
#pragma omp critical
				{
					for (Z = 0; Z <= NSeqs; Z++)
						DL[Z] = DL[Z] + DLen[Z];
				}

				free(DLen);
				free(FHolder);
			}


			free(tFMat);
			
		}

		return(1);
	}


	int MyMathFuncs::MakeOLSeqP(int LS, int BPos, int EPos, int *RSize, int *OLSeq) {
		int X;

		if (BPos < EPos) {
			for (X = 0; X < BPos; X++)
				OLSeq[X] = 0;

			for (X = EPos + 1; X <= LS; X++)
				OLSeq[X] = 0;

			RSize[0] = EPos - BPos + 1;
			for (X = BPos; X <= EPos; X++)
				OLSeq[X] = 1;
		}
		else {
			for (X = EPos + 1; X < BPos; X++)
				OLSeq[X] = 0;

			RSize[0] = EPos + LS - BPos + 1;
			for (X = 1; X <= EPos; X++)
				OLSeq[X] = 1;
			for (X = BPos; X <= LS; X++)
				OLSeq[X] = 1;
		}
		return(1);
	}
	int MyMathFuncs::MakeOLSeqP2(int LS, int BPos, int EPos, int *RSize, int *OLSeq, unsigned char *MDSeq) {
		int X;

		if (BPos < EPos) {
			for (X = 0; X < BPos; X++) {
				if (MDSeq[X] == 0)
					OLSeq[X] = 0;
				else
					OLSeq[X] = 2;
			}

			for (X = EPos + 1; X <= LS; X++) {
				if (MDSeq[X] == 0)
					OLSeq[X] = 0;
				else
					OLSeq[X] = 2;
			}

			RSize[0] = EPos - BPos + 1;
			for (X = BPos; X <= EPos; X++) {
				if (MDSeq[X] == 0)
					OLSeq[X] = 1;
				else
					OLSeq[X] = 2;
			}
		}
		else {
			for (X = EPos + 1; X < BPos; X++) {
				if (MDSeq[X] == 0)
					OLSeq[X] = 0;
				else
					OLSeq[X] = 2;
			}
			RSize[0] = EPos + LS - BPos + 1;
			for (X = 1; X <= EPos; X++) {
				if (MDSeq[X] == 0)
					OLSeq[X] = 1;
				else
					OLSeq[X] = 2;
			}
			for (X = BPos; X <= LS; X++) {
				if (MDSeq[X] == 0)
					OLSeq[X] = 1;
				else
					OLSeq[X] = 2;
			}
		}
		return(1);
	}
int MyMathFuncs::MakeCollecteventsC(int NextNo, int lenstrainseq0, int WinPP, int *RSize, int *OLSeq, int UBCM, int *CompMat, int UBRL, int *RList, int *RNum, int Addnum, int UBSM, float *SMatSmall, int *ISeqs, int *Trace, short int *PCurrentXOver, int UBPXO, XOVERDEFINE *PXOList, int UBCE,XOVERDEFINE *collectevents) {

//	Dim Dummy As Long, OLSize As Long, GoOn As Long, CTest As Long, CPar As Long, Hits(1) As Single, TotS(1) As Single, NearPair As Long, tWinPP As Long, CSeq(1) As Long, A As Long, Z As Long, X As Long, Y As Long, SQ() As Long, oldY As Long, FoundOne() As Long
//Dim tMatch(1) As Double, tMatchX() As Long, BPos2 As Long, EPos2 As Long
//Dim ZZZX As Long, ZZX As Long, BestOL() As Double, BestPV() As Double

	double *BestOL, *BestPV, *tMatch;
	int BadMatch, GoOn, A,osbpv,X,Y, rlo, cm0, cm1, Z, BPos2, EPos2,pxo;
	int *SQ, *CSeq;
	unsigned char *TList;
	float OLSize,rsh;

	if (RNum[WinPP] == -1)
		return(0);
	osbpv = RNum[WinPP] + 1;
	rlo = UBRL + 1;
	BestOL = (double*)calloc(osbpv*(Addnum * 2 + 1), sizeof(double));
	BestPV = (double*)calloc(osbpv*(Addnum * 2 + 1), sizeof(double));
	tMatch = (double*)calloc(3, sizeof(double));
	SQ = (int*)calloc(3, sizeof(int));
	CSeq = (int*)calloc(3, sizeof(int));
	TList = (unsigned char*)calloc(3 * (NextNo + 1), sizeof(unsigned char));

	cm0 = CompMat[WinPP];
	cm1 = CompMat[WinPP + UBCM+1];
//ReDim BestOL(RNum(WinPP), AddNum * 2), BestPV(RNum(WinPP), AddNum * 2), , sizeof(double))

	for (X = 0; X <= RNum[WinPP]; X++) {
		for (Y = 0; Y <= Addnum; Y++) 
			BestPV[X + Y*osbpv] = 1;
	}

	for (Y = 0; Y <= RNum[WinPP]; Y++)
		TList[WinPP + 3 * RList[WinPP + Y*rlo]] = 1;
	
	for (Y = 0; Y <= RNum[cm0]; Y++)
		TList[cm0 + RList[cm0 + Y*rlo]*3] = 1;
	
	for (Y = 0; Y <= RNum[cm1]; Y++)
		TList[cm1 + RList[cm1 + Y*rlo]*3] = 1;

				
	for (Z = 0; Z <= NextNo; Z++){
		X = Z;
		for (Y = 1; Y <= PCurrentXOver[Z]; Y++){
			BadMatch = 0;
			pxo = X + Y*(UBPXO + 1);
			SQ[0] = PXOList[pxo].Daughter;
			if (TList[WinPP + SQ[0]*3] == 1 || TList[cm0 + SQ[0]*3] == 1 || TList[cm1 + SQ[0]*3] == 1){
                        
				SQ[1] = PXOList[pxo].MajorP;
				if (SQ[1] <= NextNo){ 
					if (TList[WinPP + SQ[1]*3] == 1 || TList[cm0 + SQ[1]*3] == 1 || TList[cm1 + SQ[1]*3] == 1){                                
						SQ[2] = PXOList[pxo].MinorP;
						if (SQ[2] <= NextNo){
							if (TList[WinPP + SQ[2]*3] == 1 || TList[cm0 + SQ[2] * 3] == 1 || TList[cm1 + SQ[2]*3] == 1) {
								tMatch[0] = 3;
								BadMatch = 0;
								if (TList[WinPP + SQ[0] * 3] == 0 && TList[WinPP + SQ[1] * 3] == 0 && TList[WinPP + SQ[2] * 3]  == 0) {
									BadMatch = 1;
								}
								else if (TList[cm0 + SQ[0] * 3] == 0 && TList[cm0 + SQ[1] * 3] == 0 && TList[cm0 + SQ[2] * 3] == 0) {
									BadMatch = 1;
								}
								else if (TList[cm1 + SQ[0] * 3] == 0 && TList[cm1 + SQ[1] * 3] == 0 && TList[cm1 + SQ[2] * 3] == 0) {
									BadMatch = 1;
								}
								//ie an event involving a potentially recombinant sequence is found.
								//	'check for region overlap
								BPos2 = PXOList[pxo].Beginning;
								EPos2 = PXOList[pxo].Ending;
                                            
								OLSize = (float)(FindOverlapP(lenstrainseq0, BPos2, EPos2, &RSize[0], &OLSeq[0]));

								if (OLSize > 0) {
									rsh = (float)(RSize[0] + RSize[1]);
									tMatch[1] = (OLSize * 2) / rsh;
									
								}
								else
									tMatch[1] = 0;
									
                                           
								tMatch[1] = round(tMatch[1] * 100000) / 100000;
                                            
								if ((tMatch[1] > 0.1 && BadMatch == 0) || (tMatch[1] > 0.5 && BadMatch == 1)){
                                                
                                                
									GoOn = 0;
									for (A = 0; A <= RNum[WinPP]; A++){
										if (RList[WinPP + A*rlo] == SQ[0] || RList[WinPP + A*rlo] == SQ[1] || RList[WinPP + A*rlo] == SQ[2]){
											CSeq[1] = A;
											GoOn = 1;
											break;
										}
									}
                                                    
									if (GoOn == 1 && (BestPV[CSeq[1] + PXOList[pxo].ProgramFlag*osbpv] > PXOList[pxo].Probability || BestPV[CSeq[1] + PXOList[pxo].ProgramFlag*osbpv] == 0)) {
										BestPV[CSeq[1] + PXOList[pxo].ProgramFlag*osbpv] = PXOList[pxo].Probability;
										BestOL[CSeq[1] + PXOList[pxo].ProgramFlag*osbpv] = tMatch[1];
										collectevents[CSeq[1] + PXOList[pxo].ProgramFlag*(UBCE+1)] = PXOList[pxo];


									}

                                            
								}
							}
						}
					}
				}
			}
		}
	}
	//need to make sure that the event pointed to by trace(0) and trace(1) is represented in collectevents
	for (A = 0; A <= RNum[WinPP]; A++){
		if (RList[WinPP + A*rlo] == ISeqs[0] || RList[WinPP + A*rlo] == ISeqs[1] || RList[WinPP + A*rlo] == ISeqs[2]) {
			CSeq[1] = A;
			break;
		}
	}

	collectevents[CSeq[1] + PXOList[Trace[0] + Trace[1] * (UBPXO + 1)].ProgramFlag*(UBCE + 1)] = PXOList[Trace[0] + Trace[1] * (UBPXO + 1)];
	free (BestOL);
	free(BestPV);
	free(SQ);
	free(TList);
	free(tMatch);
	free(CSeq);
	return(1);
}



	double MyMathFuncs::CalcChiValsP2(int UBWS, int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, float *chitable)
	{

		double ChiH, MChi;
		int A, B, C, D;
		int  X, LO, FO, SO, bp1, os1, os2, S1, F1;
		//float ;
		FO = UBWS + 1;// LenSeq + (HWindowWidth)* 2 + 1;
		SO = LenSeq + 1;
		S1 = SO * 2;
		F1 = FO * 2;
		LO = HWindowWidth;
		ChiH = 0;
		MChi = 0;

		os1 = HWindowWidth + 1;

		for (X = 0; X < LenXoverSeq; X++) {
			/*bp1 = X - HWindowWidth;
			if (bp1 < 1)
				bp1 = bp1 + LenXoverSeq;*/
			//if (BanWins[X] == 0 && BanWins[bp1] == 0) {


				ChiVals[X] = chitable[WinScores[X] + WinScores[X + LO] * os1];
				ChiVals[X + SO] = chitable[WinScores[X + FO] + WinScores[X + LO + FO] * os1];
				ChiVals[X + S1] = chitable[WinScores[X + F1] + WinScores[X + LO + F1] * os1];
			/*}
			else {

				*(ChiVals + X) = 0;
				*(ChiVals + X + SO) = 0;
				*(ChiVals + X + S1) = 0;
			}*/
		}

		for (X = 0; X < LenXoverSeq; X++) {
			if (ChiVals[X] > MChi)
				MChi = ChiVals[X];
		}
		for (X = 0; X < LenXoverSeq; X++) {
			if (ChiVals[X + SO] > MChi)
				MChi = ChiVals[X + SO];
		}
		for (X = 0; X < LenXoverSeq; X++) {
			if (ChiVals[X + S1] > MChi)
				MChi = ChiVals[X + S1];
		}
		return (MChi);
	}


	double MyMathFuncs::GCGetHiPValP(int lseq, int LenXoverSeq, int *FragCount, double *PVals, int *MaxY, int *MaxX, int *highenough) {

		double MPV;
		int X, Y, os, os2,mx,my;

		//os = LenXoverSeq+1 ;
		os = lseq + 1;
		MPV = 100;
		mx = -1;
		my = -1;
		//had trouble ensuring that the actual max value gets returned
//#pragma omp parallel for private (X, os2, Y)
		for (X = 0; X < 6; X++) {

			if (highenough[X] == 1) {
				os2 = X*os;

				for (Y = 0; Y <= FragCount[X]; Y++) {

					if (PVals[Y + os2] < MPV) {
						
//#pragma omp critical
//						{
							mx = X;
							my = Y;
//#pragma omp flush(MPV)
							MPV = PVals[Y + os2];
//						}
					}
				}
			}
		}
		*MaxX = mx;
		*MaxY = my;
		return(MPV);
	}

	int MyMathFuncs::DelPValsP(short int GCMaxOverlapFrags, int Y, int X, int LS, double *PVals, int *FragCount, int *FragSt, int *FragEn, int *MaxScorePos, int *DeleteArray) {
		//pvals, fragsst/fragen/maxscorepos - lenseq,6
		int Z, off1, off2, off3, GoOn;
		off1 = FragSt[Y + X*(LS + 1)];
		off2 = FragEn[MaxScorePos[Y + X*(LS + 1)] + X*(LS + 1)];
		off3 = Y + X*(LS + 1);
		GoOn = 1;
		if (off1 < off2) {
			for (Z = off1; Z <= off2; Z++) {
				if (DeleteArray[Z] >= GCMaxOverlapFrags) {
					GoOn = 0;
					PVals[off3] = 100;
					break;
				}
			}
		}
		else {
			for (Z = off1; Z <= FragCount[X]; Z++) {
				if (DeleteArray[Z] >= GCMaxOverlapFrags) {
					GoOn = 0;
					PVals[off3] = 100;
					break;
				}
			}
			if (Z >= FragCount[X]) {
				for (Z = 1; Z <= off2; Z++) {
					if (DeleteArray[Z] >= GCMaxOverlapFrags) {
						GoOn = 0;
						PVals[off3] = 100;
						break;
					}
				}
			}
		}


		return(GoOn);
	}

	int MyMathFuncs::MakeDeleteArrayP(int FragSt, int FragEn, int FragCount, int *DeleteArray) {

		int B;



		if (FragSt < FragEn) {
			for (B = FragSt; B <= FragEn; B++)
				DeleteArray[B] = DeleteArray[B] + 1;
		}
		else {
			for (B = FragSt; B <= FragCount; B++)
				DeleteArray[B] = DeleteArray[B] + 1;
			for (B = 0; B <= FragEn; B++)
				DeleteArray[B] = DeleteArray[B] + 1;

		}
		return(1);


	}

	int MyMathFuncs::FindMissingP(int LS, int Seq1, int Seq2, int Seq3, int Z, int En, unsigned char *MissingData) {
		int B, off1, off2, off3, GoOn;
		off1 = Seq1*(LS + 1);
		off2 = Seq2*(LS + 1);
		off3 = Seq3*(LS + 1);

		if (Z < En) {
			for (B = En; B >= Z; B--) {
				if (MissingData[B + off1] == 1 || MissingData[B + off2] == 1 || MissingData[B + off3] == 1)
					break;
			}
		}
		else {
			GoOn = 0;
			for (B = En; B >= 1; B--) {
				if (MissingData[B + off1] == 1 || MissingData[B + off2] == 1 || MissingData[B + off3] == 1) {
					GoOn = 1;
					break;
				}
			}
			if (GoOn == 0) {
				for (B = LS; B >= Z; B--) {
					if (MissingData[B + off1] == 1 || MissingData[B + off2] == 1 || MissingData[B + off3] == 1)
						break;
				}
			}
		}
		return(B);
	}

	int  MyMathFuncs::CheckSplitP(int step,int LS, int Be, int En, int Seq1, int Seq2, int Seq3, int *Split, unsigned char *MissingData) {
		int Z, off1, off2, off3;
		off1 = Seq1*(LS + 1);
		off2 = Seq2*(LS + 1);
		off3 = Seq3*(LS + 1);
		if (Be < En) {
			for (Z = Be; Z <= En; Z=Z+step) {
				if (MissingData[Z + off1] == 1 || MissingData[Z + off2] == 1 || MissingData[Z + off3] == 1) {
					*Split = 1;
					break;
				}
			}
		}
		else {
			for (Z = Be; Z <= LS; Z = Z + step) {
				if (MissingData[Z + off1] == 1 || MissingData[Z + off2] == 1 || MissingData[Z + off3] == 1) {
					*Split = 1;
					break;
				}
			}
			if (*Split == 0) {
				for (Z = 1; Z <= En; Z = Z + step) {
					if (MissingData[Z + off1] == 1 || MissingData[Z + off2] == 1 || MissingData[Z + off3] == 1) {
						*Split = 1;
						break;
					}
				}
			}

		}
		return(Z);
	}

	double MyMathFuncs::MakeSubProbP(int X, int LS, int LenXoverSeq, int BTarget, int ETarget, char *SubSeq, double *LL, double *KMax, double *MissPen, double *critval) {
		int C, t1, t2;
		t1 = 0;
		t2 = 0;
		double FragSc, Polys, LKLen, KAScore, PV, THld, warn;
		FragSc = 0.0;
		if (X<3) {
			for (C = BTarget; C <= ETarget; C++)
				FragSc = FragSc + SubSeq[C + X*(LS + 1)];
		}
		else {
			if (X == 3) {
				t1 = 0;
				t2 = 1;
			}
			else if (X == 4) {
				t1 = 0;
				t2 = 2;
			}
			else if (X == 5) {
				t1 = 1;
				t2 = 2;
			}

			for (C = BTarget; C <= ETarget; C++)
				FragSc = FragSc + SubSeq[C + t1*(LS + 1)] + SubSeq[C + t2*(LS + 1)];

		}

		if (BTarget < ETarget)
			Polys = ETarget - BTarget + 1;
		else
			Polys = ETarget + LenXoverSeq - BTarget + 1;

		FragSc = FragSc - (Polys - FragSc) * MissPen[X];

		if (FragSc > critval[X]) {
			LKLen = log(KMax[X] * LenXoverSeq);
			KAScore = (LL[X] * FragSc) - LKLen;
			if (KAScore > 0) {
				if (KAScore < 32) {

					THld = exp(-KAScore);
					PV = 1 - exp(-THld);

				}
				else {
					warn = 0;
					if (KAScore > 700) {

						warn = KAScore;
						KAScore = 701;

					}

					THld = exp(-KAScore);
					if (warn != 0) {
						KAScore = warn - 700;
						THld = THld / KAScore;
					}
					PV = THld;
				}

			}
			else
				PV = 1;
			return (PV);
		}
		else
			return (1);
	}

	double MyMathFuncs::ChiPVal2P(double X) {

		long double PValHolder;
		double ChiPValx;
		if (X == 0)
			ChiPValx = 1;
		else {

			PValHolder = (NormalZ(-sqrt(X)));

			if (PValHolder == 0) {//< 0.0000000001){
								  //if (X > 35)
								  //	Y=X;
								  //else
								  //	X=35.1;
				PValHolder = 0.0;
				PValHolder = 0.0000000001;
				PValHolder = PValHolder / (X - 34);
			}

			ChiPValx = (double)(PValHolder);
		}
		return(ChiPValx);

	}

	int MyMathFuncs::GetACP(int LenXoverSeq, int LS, int MaxY, int MaxX, int TWin, int *A, int *C,unsigned char *Scores) {
		int X, os;
		os = MaxY*(LS + 1);
		//scores ls,2
		//Calculate first A and first C for window size Twin
		for (X = ((MaxX + 1) - TWin); X <= MaxX; X++) {
			if (X > 0) {
				if (X <= LenXoverSeq)
					*A = *A + Scores[X + os];
				else
					*A = *A + Scores[X - LenXoverSeq + os];
			}
			else
				*A = *A + Scores[LenXoverSeq + X + os];

		}
		for (X = (MaxX + 1); X <= MaxX + TWin; X++) {
			if (X > 0) {
				if (X <= LenXoverSeq)
					*C = *C + Scores[X + os];
				else
					*C = *C + Scores[X - LenXoverSeq + os];
			}
			else
				*C = *C + Scores[LenXoverSeq + X + os];
		}
		return(1);
	}

	int MyMathFuncs::ClearDeleteArray(int ls, int *da)
	{
		int X;
		for (X = 0; X <= ls; X++)
			da[X] = 0;
		return(1);

	}

	int MyMathFuncs::ClearDeleteArrayB(int ls, unsigned char *da)
	{
		int X;
		for (X = 0; X <= ls; X++)
			da[X] = 0;
		return(1);

	}
	int MyMathFuncs::FindSubSeqGCAP(char gcindelflag, int LSeq, int seq1, int seq2, int seq3, short int *SeqNum, char *SubSeq, int *XPosDiff, int *XDiffPos, int *NDiff) {

		int Z, X, Y, SX;
		int s1, s2, s3, os1, os2, os3, mp, s1s, s2s, s3s, off1, off2, off3, GoOn;

		Y = 0;

		os1 = seq1 * (LSeq + 1);
		os2 = seq2 * (LSeq + 1);
		os3 = seq3 * (LSeq + 1);
		off1 = LSeq + 1;
		off2 = (LSeq + 1) * 2;
		off3 = (LSeq + 1) * 6;
		X = 0;
		//for (X = 0; X <= LSeq; X++)
		XPosDiff[X] = 0;
		XDiffPos[X] = 0;

		if (gcindelflag == 0) {

			for (X = 1; X <= LSeq; X++) {

				XPosDiff[X] = Y;
				s1 = SeqNum[X + os1];
				s2 = SeqNum[X + os2];
				s3 = SeqNum[X + os3];

				if (s1 != s2 || s1 != s3) {
					if (s1 != 46) {
						if (s2 != 46) {
							if (s3 != 46) {

								Y++;
								XPosDiff[X] = Y;
								XDiffPos[Y] = X;
								SubSeq[Y] = (s1 == s2);
								SubSeq[Y + off1] = (s1 == s3);
								SubSeq[Y + off2] = (s2 == s3);
							}
						}
					}

					SubSeq[Y + off3] = (SubSeq[Y] == 0 && SubSeq[Y + off1] == 0 && SubSeq[Y + off2] == 0);

				}

			}
		}
		else if (gcindelflag == 1) {
			for (X = 1; X <= LSeq; X++) {


				s1 = SeqNum[X + os1];
				s2 = SeqNum[X + os2];
				s3 = SeqNum[X + os3];
				GoOn = 0;
				if (s1 != s2) {

					if (s1 != 46) {
						if (s2 != 46) {
							if (s3 != 46) {
								Y++;
								XPosDiff[X] = Y;
								XDiffPos[Y] = X;
								SubSeq[Y] = 0;
								SubSeq[Y + off1] = (s1 == s3);
								SubSeq[Y + off2] = (s2 == s3);
							}
							else
								GoOn = 1;
						}
						else
							GoOn = 1;
					}
					else
						GoOn = 1;

					if (GoOn == 1) {
						XPosDiff[X] = Y;
						SX = X;


						if (s1 == s2) {
							s1s = 1;
							s2s = 0;
							s3s = 0;
						}
						else {
							s1s = 0;
							if (s1 == s3) {
								s2s = 1;
								s3s = 0;
							}
							else {
								s2s = 0;
								if (s2 == s3)
									s3s = 1;
								else
									s3s = 0;
							}
						}


						X++;

						if (X <= LSeq) {


							s1 = SeqNum[X + os1];
							s2 = SeqNum[X + os2];
							s3 = SeqNum[X + os3];

							//for (Z = X; Z <= LSeq; Z++){
							while (s1 == 46 || s2 == 46 || s3 == 46) {
								//	s1 = SeqNum[Z + os1];
								//	s2 = SeqNum[Z + os2];
								//	s3 = SeqNum[Z + os3];

								//if (s1 != 46 && s2 != 46 && s3 != 46) 
								//	break;	
								if (s1s != 0) {
									if (s1 != s2)
										s1s = 0;
								}

								if (s2s != 0) {
									if (s1 != s3)
										s2s = 0;
								}

								if (s3s != 0) {
									if (s2 != s3)
										s3s = 0;

								}
								X++;
								if (X > LSeq)
									break;

								s1 = SeqNum[X + os1];

								s2 = SeqNum[X + os2];
								s3 = SeqNum[X + os3];

							}
							//X=Z;
						}



						X--;
						Y++;

						mp = (int)(SX + (X - SX) / 2);

						for (Z = SX + 1; Z <= X; Z++)
							XPosDiff[Z] = 0;

						XPosDiff[mp] = Y;
						XDiffPos[Y] = mp;

						SubSeq[Y] = s1s;
						SubSeq[Y + off1] = s2s;
						SubSeq[Y + off2] = s3s;
					}

					SubSeq[Y + off3] = (SubSeq[Y] == 0 && SubSeq[Y + off1] == 0 && SubSeq[Y + off2] == 0);
				}
				else if (s1 != s3) {

					if (s1 != 46) {
						if (s2 != 46) {
							if (s3 != 46) {
								Y++;
								XPosDiff[X] = Y;
								XDiffPos[Y] = X;
								SubSeq[Y] = 1;
								SubSeq[Y + off1] = 0;
								SubSeq[Y + off2] = (s2 == s3);
							}
							else
								GoOn = 1;
						}
						else
							GoOn = 1;
					}
					else
						GoOn = 1;

					if (GoOn == 1) {
						XPosDiff[X] = Y;
						SX = X;


						if (s1 == s2) {
							s1s = 1;
							s2s = 0;
							s3s = 0;
						}
						else {
							s1s = 0;
							if (s1 == s3) {
								s2s = 1;
								s3s = 0;
							}
							else {
								s2s = 0;
								if (s2 == s3)
									s3s = 1;
								else
									s3s = 0;
							}
						}


						X++;

						if (X <= LSeq) {


							s1 = SeqNum[X + os1];
							s2 = SeqNum[X + os2];
							s3 = SeqNum[X + os3];

							//for (Z = X; Z <= LSeq; Z++){
							while (s1 == 46 || s2 == 46 || s3 == 46) {
								//	s1 = SeqNum[Z + os1];
								//	s2 = SeqNum[Z + os2];
								//	s3 = SeqNum[Z + os3];

								//if (s1 != 46 && s2 != 46 && s3 != 46) 
								//	break;	
								if (s1s != 0) {
									if (s1 != s2)
										s1s = 0;
								}

								if (s2s != 0) {
									if (s1 != s3)
										s2s = 0;
								}

								if (s3s != 0) {
									if (s2 != s3)
										s3s = 0;

								}
								X++;
								if (X > LSeq)
									break;

								s1 = SeqNum[X + os1];
								s2 = SeqNum[X + os2];
								s3 = SeqNum[X + os3];

							}
							//X=Z;
						}



						X--;
						Y++;

						mp = (int)(SX + (X - SX) / 2);

						for (Z = SX + 1; Z <= X; Z++)
							XPosDiff[Z] = 0;

						XPosDiff[mp] = Y;
						XDiffPos[Y] = mp;

						SubSeq[Y] = s1s;
						SubSeq[Y + off1] = s2s;
						SubSeq[Y + off2] = s3s;
					}

					SubSeq[Y + off3] = (SubSeq[Y] == 0 && SubSeq[Y + off1] == 0 && SubSeq[Y + off2] == 0);
				}
				else


					XPosDiff[X] = Y;

			}
		}
		else if (gcindelflag == 2) {
			for (X = 1; X <= LSeq; X++) {

				XPosDiff[X] = Y;
				s1 = SeqNum[X + os1];
				s2 = SeqNum[X + os2];
				s3 = SeqNum[X + os3];

				if (s1 != s2 || s1 != s3) {
					if (s1 != 46 && s2 != 46 && s3 != 46) {
						Y++;
						XPosDiff[X] = Y;
						XDiffPos[Y] = X;
						SubSeq[Y] = (s1 == s2);
						SubSeq[Y + LSeq + 1] = (s1 == s3);
						SubSeq[Y + (LSeq + 1) * 2] = (s2 == s3);
					}

					else {
						Y++;
						XPosDiff[X] = Y;
						XDiffPos[Y] = X;
						SubSeq[Y] = (s1 == s2);
						SubSeq[Y + (LSeq + 1)] = (s1 == s3);
						SubSeq[Y + (LSeq + 1) * 2] = (s2 == s3);
					}


					SubSeq[Y + (LSeq + 1) * 6] = (SubSeq[Y] == 0 && SubSeq[Y + (LSeq + 1)] == 0 && SubSeq[Y + (LSeq + 1) * 2] == 0);


				}

			}

		}
		for (X = 1; X <= Y; X++) {
			NDiff[0] = NDiff[0] + SubSeq[X];
			NDiff[1] = NDiff[1] + SubSeq[X + off1];
			NDiff[2] = NDiff[2] + SubSeq[X + off2];
		}

		XDiffPos[Y + 1] = LSeq;

		Z = XDiffPos[Y] + 1;

		for (X = Z; X <= LSeq; X++)
			XPosDiff[X] = Y;

		return(Y);
	}
	
	int MyMathFuncs::FindSubSeqGCAP6(int ubcs,unsigned char *cs, int ubfss, unsigned char *fssgc, char gcindelflag, int LSeq, int seq1, int seq2, int seq3, char *SubSeq, int *NDiff) {

		int h,Z, X, Y, SX, osf;
		int s1, s2, s3, os1, os2, os3, os4, os5,mp, s1s, s2s, s3s, off1, off2, off3, GoOn, ofs;

		Y = 0;

		os1 = seq1 * (ubcs + 1);
		os2 = seq2 * (ubcs + 1);
		os3 = seq3 * (ubcs + 1);
		os4 = (ubfss + 1) * 4;
		os5 = (ubfss + 1)*(ubfss + 1) * 4;
		off1 = LSeq + 1;
		off2 = (LSeq + 1) * 2;
		off3 = (LSeq + 1) * 6;
		X = 0;
		//for (X = 0; X <= LSeq; X++)
		
		//if (gcindelflag == 0) {

		for (X = 1; X <= ubcs; X++) {


			s1 = cs[X + os1];
			s2 = cs[X + os2];
			s3 = cs[X + os3];
			osf = s1 * 4 + os4 * s2 + os5 * s3;
			if (fssgc[3 + osf] != 0) {
				for (int z = 0; z <= 2; z++) {
					
					h = fssgc[z + osf];
					if (h > 0) {
						h--;
						Y++;
						SubSeq[Y + h*(LSeq + 1)] = 1;
					}
				}

			}
			/*if (s1 != s2 || s1 != s3) {
				if (s1 != 46) {
					if (s2 != 46) {
						if (s3 != 46) {

							Y++;

							SubSeq[Y] = (s1 == s2);
							SubSeq[Y + off1] = (s1 == s3);
							SubSeq[Y + off2] = (s2 == s3);
						}
					}
				}

				SubSeq[Y + off3] = (SubSeq[Y] == 0 && SubSeq[Y + off1] == 0 && SubSeq[Y + off2] == 0);

			}

		}*/
		}
		//}
		
		for (X = 1; X <= Y; X++) {
			NDiff[0] = NDiff[0] + SubSeq[X];
			NDiff[1] = NDiff[1] + SubSeq[X + off1];
			NDiff[2] = NDiff[2] + SubSeq[X + off2];
		}

		
		

		return(Y);
	}
	int MyMathFuncs::FindSubSeqGCAP7(int ubcs, unsigned char *cs, int ubfss, unsigned char *fssgc, char gcindelflag, int LSeq, int seq1, int seq2, int seq3, char *SubSeq, int *NDiff, int *XDiffPos, int *XPosDiff) {

		int h, Z, X, Y, SX, osf;
		int s1, s2, s3, os1, os2, os3, os4, os5, mp, s1s, s2s, s3s, off1, off2, off3, GoOn, ofs;

		Y = 0;

		XPosDiff[0] = 0;
		XDiffPos[0] = 0;
		os1 = seq1 * (ubcs + 1);
		os2 = seq2 * (ubcs + 1);
		os3 = seq3 * (ubcs + 1);
		os4 = (ubfss + 1) * 4;
		os5 = (ubfss + 1)*(ubfss + 1) * 4;
		off1 = LSeq + 1;
		off2 = (LSeq + 1) * 2;
		off3 = (LSeq + 1) * 6;
		X = 0;
		//for (X = 0; X <= LSeq; X++)

		//if (gcindelflag == 0) {

		for (X = 1; X <= ubcs; X++) {


			s1 = cs[X + os1];
			s2 = cs[X + os2];
			s3 = cs[X + os3];
			osf = s1 * 4 + os4 * s2 + os5 * s3;
			if (fssgc[3 + osf] != 0) {
				for (int z = 0; z <= 2; z++) {
					
					h = fssgc[z + osf];
					if (h > 0) {
						h--;
						Y++;
						XPosDiff[z+(X-1)*3+1] = Y;
						XDiffPos[Y] = z + (X-1) * 3+1;
						SubSeq[Y + h*(LSeq + 1)] = 1;
					}
					else {
						XPosDiff[z + ((X-1) * 3)+1] = Y;
					}

				}

			
			}
			else{
				XPosDiff[(X-1)*3+1] = Y;
				XPosDiff[(X-1)*3+2] = Y;
				XPosDiff[(X-1)*3+3] = Y;
			
			}
				
			/*if (s1 != s2 || s1 != s3) {
			if (s1 != 46) {
			if (s2 != 46) {
			if (s3 != 46) {

			Y++;

			SubSeq[Y] = (s1 == s2);
			SubSeq[Y + off1] = (s1 == s3);
			SubSeq[Y + off2] = (s2 == s3);
			}
			}
			}

			SubSeq[Y + off3] = (SubSeq[Y] == 0 && SubSeq[Y + off1] == 0 && SubSeq[Y + off2] == 0);

			}

			}*/
		}
		//}


		XDiffPos[Y + 1] = LSeq;

		Z = XDiffPos[Y] + 1;

		for (X = Z; X <= LSeq; X++)
			XPosDiff[X] = Y;

		

		for (X = 1; X <= Y; X++) {
			NDiff[0] = NDiff[0] + SubSeq[X];
			NDiff[1] = NDiff[1] + SubSeq[X + off1];
			NDiff[2] = NDiff[2] + SubSeq[X + off2];
		}




		return(Y);
	}
	int MyMathFuncs::FindSubSeqGCAP5(char gcindelflag, int LSeq, int seq1, int seq2, int seq3, short int *SeqNum, char *SubSeq, int *NDiff) {

		int Z, X, Y, SX;
		int s1, s2, s3, os1, os2, os3, mp, s1s, s2s, s3s, off1, off2, off3, GoOn;

		Y = 0;

		os1 = seq1 * (LSeq + 1);
		os2 = seq2 * (LSeq + 1);
		os3 = seq3 * (LSeq + 1);
		off1 = LSeq + 1;
		off2 = (LSeq + 1) * 2;
		off3 = (LSeq + 1) * 6;
		X = 0;
		//for (X = 0; X <= LSeq; X++)

		if (gcindelflag == 0) {

			for (X = 1; X <= LSeq; X++) {


				s1 = SeqNum[X + os1];
				s2 = SeqNum[X + os2];
				s3 = SeqNum[X + os3];

				if (s1 != s2 || s1 != s3) {
					if (s1 != 46) {
						if (s2 != 46) {
							if (s3 != 46) {

								Y++;

								SubSeq[Y] = (s1 == s2);
								SubSeq[Y + off1] = (s1 == s3);
								SubSeq[Y + off2] = (s2 == s3);
							}
						}
					}

					SubSeq[Y + off3] = (SubSeq[Y] == 0 && SubSeq[Y + off1] == 0 && SubSeq[Y + off2] == 0);

				}

			}
		}
		else if (gcindelflag == 1) {
			for (X = 1; X <= LSeq; X++) {


				s1 = SeqNum[X + os1];
				s2 = SeqNum[X + os2];
				s3 = SeqNum[X + os3];
				GoOn = 0;
				if (s1 != s2) {

					if (s1 != 46) {
						if (s2 != 46) {
							if (s3 != 46) {
								Y++;

								SubSeq[Y] = 0;
								SubSeq[Y + off1] = (s1 == s3);
								SubSeq[Y + off2] = (s2 == s3);
							}
							else
								GoOn = 1;
						}
						else
							GoOn = 1;
					}
					else
						GoOn = 1;

					if (GoOn == 1) {

						SX = X;


						if (s1 == s2) {
							s1s = 1;
							s2s = 0;
							s3s = 0;
						}
						else {
							s1s = 0;
							if (s1 == s3) {
								s2s = 1;
								s3s = 0;
							}
							else {
								s2s = 0;
								if (s2 == s3)
									s3s = 1;
								else
									s3s = 0;
							}
						}


						X++;

						if (X <= LSeq) {


							s1 = SeqNum[X + os1];
							s2 = SeqNum[X + os2];
							s3 = SeqNum[X + os3];

							//for (Z = X; Z <= LSeq; Z++){
							while (s1 == 46 || s2 == 46 || s3 == 46) {
								//	s1 = SeqNum[Z + os1];
								//	s2 = SeqNum[Z + os2];
								//	s3 = SeqNum[Z + os3];

								//if (s1 != 46 && s2 != 46 && s3 != 46) 
								//	break;	
								if (s1s != 0) {
									if (s1 != s2)
										s1s = 0;
								}

								if (s2s != 0) {
									if (s1 != s3)
										s2s = 0;
								}

								if (s3s != 0) {
									if (s2 != s3)
										s3s = 0;

								}
								X++;
								if (X > LSeq)
									break;

								s1 = SeqNum[X + os1];

								s2 = SeqNum[X + os2];
								s3 = SeqNum[X + os3];

							}
							//X=Z;
						}



						X--;
						Y++;

						mp = (int)(SX + (X - SX) / 2);


						SubSeq[Y] = s1s;
						SubSeq[Y + off1] = s2s;
						SubSeq[Y + off2] = s3s;
					}

					SubSeq[Y + off3] = (SubSeq[Y] == 0 && SubSeq[Y + off1] == 0 && SubSeq[Y + off2] == 0);
				}
				else if (s1 != s3) {

					if (s1 != 46) {
						if (s2 != 46) {
							if (s3 != 46) {
								Y++;

								SubSeq[Y] = 1;
								SubSeq[Y + off1] = 0;
								SubSeq[Y + off2] = (s2 == s3);
							}
							else
								GoOn = 1;
						}
						else
							GoOn = 1;
					}
					else
						GoOn = 1;

					if (GoOn == 1) {

						SX = X;


						if (s1 == s2) {
							s1s = 1;
							s2s = 0;
							s3s = 0;
						}
						else {
							s1s = 0;
							if (s1 == s3) {
								s2s = 1;
								s3s = 0;
							}
							else {
								s2s = 0;
								if (s2 == s3)
									s3s = 1;
								else
									s3s = 0;
							}
						}


						X++;

						if (X <= LSeq) {


							s1 = SeqNum[X + os1];
							s2 = SeqNum[X + os2];
							s3 = SeqNum[X + os3];

							//for (Z = X; Z <= LSeq; Z++){
							while (s1 == 46 || s2 == 46 || s3 == 46) {
								//	s1 = SeqNum[Z + os1];
								//	s2 = SeqNum[Z + os2];
								//	s3 = SeqNum[Z + os3];

								//if (s1 != 46 && s2 != 46 && s3 != 46) 
								//	break;	
								if (s1s != 0) {
									if (s1 != s2)
										s1s = 0;
								}

								if (s2s != 0) {
									if (s1 != s3)
										s2s = 0;
								}

								if (s3s != 0) {
									if (s2 != s3)
										s3s = 0;

								}
								X++;
								if (X > LSeq)
									break;

								s1 = SeqNum[X + os1];
								s2 = SeqNum[X + os2];
								s3 = SeqNum[X + os3];

							}
							//X=Z;
						}



						X--;
						Y++;

						mp = (int)(SX + (X - SX) / 2);


						SubSeq[Y] = s1s;
						SubSeq[Y + off1] = s2s;
						SubSeq[Y + off2] = s3s;
					}

					SubSeq[Y + off3] = (SubSeq[Y] == 0 && SubSeq[Y + off1] == 0 && SubSeq[Y + off2] == 0);
				}


			}
		}
		else if (gcindelflag == 2) {
			for (X = 1; X <= LSeq; X++) {


				s1 = SeqNum[X + os1];
				s2 = SeqNum[X + os2];
				s3 = SeqNum[X + os3];

				if (s1 != s2 || s1 != s3) {
					if (s1 != 46 && s2 != 46 && s3 != 46) {
						Y++;

						SubSeq[Y] = (s1 == s2);
						SubSeq[Y + LSeq + 1] = (s1 == s3);
						SubSeq[Y + (LSeq + 1) * 2] = (s2 == s3);
					}

					else {
						Y++;

						SubSeq[Y] = (s1 == s2);
						SubSeq[Y + (LSeq + 1)] = (s1 == s3);
						SubSeq[Y + (LSeq + 1) * 2] = (s2 == s3);
					}


					SubSeq[Y + (LSeq + 1) * 6] = (SubSeq[Y] == 0 && SubSeq[Y + (LSeq + 1)] == 0 && SubSeq[Y + (LSeq + 1) * 2] == 0);


				}

			}

		}
		for (X = 1; X <= Y; X++) {
			NDiff[0] = NDiff[0] + SubSeq[X];
			NDiff[1] = NDiff[1] + SubSeq[X + off1];
			NDiff[2] = NDiff[2] + SubSeq[X + off2];
		}




		return(Y);
	}

	int MyMathFuncs::FindSubSeqGCAP2(char gcindelflag, int LSeq, int seq1, int seq2, int seq3, short int *SeqNum, char *SubSeq, int *XPosDiff, int *XDiffPos, int *NDiff, unsigned char *binarray) {

		int Z, X, Y, SX;
		int s1, s2, s3, os1, os2, os3, mp, s1s, s2s, s3s, off1, off2, off3, GoOn, target, se3os, se2os, ah0,ah1,ah2, SY, os1x,os3x, os2x;
		unsigned char  ba2, ba3;
		ah0 = 0;
		ah1 = 0;
		ah2 = 0;

		Y = 0;

		os1 = seq1 * (LSeq + 1);
		os2 = seq2 * (LSeq + 1);
		os3 = seq3 * (LSeq + 1);
		off1 = LSeq + 1;
		off2 = (LSeq + 1) * 2;
		off3 = (LSeq + 1) * 6;
		
		
		//SubSeq(Len(StrainSeq(0)), 6)
		//BinArray(Len(StrainSeq(0)), Nextno)
		target = off1 * 7;
		for (X = 1; X < target; X++)
			SubSeq[X] = 0;
		X = 0;
		//for (X = 0; X <= LSeq; X++)
		XPosDiff[X] = 0;

		if (gcindelflag == 0) {
			
			for (X = 1; X <= LSeq; X++) {

				se2os = X + os2;
				
				if (binarray[se2os] == 1 ) {
					
					se3os = X + os3;
					if (binarray[se3os] == 1 ) {


							XDiffPos[++Y] = X;
							if (SeqNum[se2os] == SeqNum[se3os]) {
								SubSeq[Y + off2] = 1;
								ah2 ++;
							}
							else
								SubSeq[Y + off3] = 1;
					}
					else if (binarray[se3os] == 0) {
						XDiffPos[++Y] = X;
						SubSeq[Y + off1] = 1;//(s1 == s3);
						ah1++;
					}
					

				}
				else if (binarray[X + os3] == 1) {
					if (binarray[se2os] == 0 ) {

						
						XDiffPos[++Y] = X;
						SubSeq[Y] = 1;// (s1 == s2);
						ah0++;
					}

				}
				XPosDiff[X] = Y;
			}
			NDiff[0] = ah0; //NDiff[0] + SubSeq[X];
			NDiff[1] = ah1;// NDiff[1] + SubSeq[X + off1];
			NDiff[2] = ah2;// NDiff[2] + SubSeq[X + off2];

		}
		else if (gcindelflag == 1) {
			for (X = 1; X <= LSeq; X++) {

				
				ba2 =binarray[X + os2];
				ba3 = binarray[X + os3];
				GoOn = 0;
				if (ba2 != 0) {
					
					if (ba2 == 1) {
						if (ba3 == 0) {
							
								
								XDiffPos[++Y] = X;
								XPosDiff[X] = Y;
								SubSeq[Y + off1] = 1;
								ah1++;
							
						}
						else if (ba3 == 1) {
							
							XDiffPos[++Y] = X;
							XPosDiff[X] = Y;
							s2s = (SeqNum[X + os2] == SeqNum[X + os3]);
							SubSeq[Y + off2] = s2s;
							ah2 += s2s;
							SubSeq[Y + off3] = (s2s == 0);
						}
						else
							GoOn = 1;
					}
					else
						GoOn = 1;

					if (GoOn == 1) {
						XPosDiff[X] = Y;
						SX = X;
						if (ba3 == 0) {
							s2s = 1;
							
							X++;

							if (X <= LSeq) {


								s1 = SeqNum[X + os1];
								s3 = SeqNum[X + os3];

								//for (Z = X; Z <= LSeq; Z++){
								while (s1 == 46 || SeqNum[X + os2] == 46 || s3 == 46) {//(binarray[se2os] > 1 || binarray[se3os] > 1) {
																		  //	s1 = SeqNum[Z + os1];

									if (s2s != 0) {
										if (s1 != s3)
											s2s = 0;
									}

									
									X++;
									if (X > LSeq)
										break;

									s1 = SeqNum[X + os1];
									s3 = SeqNum[X + os3];

								}
								//X=Z;
							}
							ah1 += s2s;
							SubSeq[++Y + off1] = s2s;
							SubSeq[Y + off3] = (s2s == 0);
						}
						else {
							
							s3s = (SeqNum[X + os2] == SeqNum[X + os3]);
							if (s3s == 1) {
								X++;

								if (X <= LSeq) {


									
									s2 = SeqNum[X + os2];
									s3 = SeqNum[X + os3];

									//for (Z = X; Z <= LSeq; Z++){
									while (SeqNum[X + os1] == 46 || s2 == 46 || s3 == 46) {//(binarray[se2os] > 1 || binarray[se3os] > 1) {
										
										if (s3s != 0) {
											if (s2 != s3)
												s3s = 0;

										}
										X++;
										if (X > LSeq)
											break;

										
										s2 = SeqNum[X + os2];
										s3 = SeqNum[X + os3];

									}
									
								}
								
								ah2 += s3s;

								
								SubSeq[++Y + off2] = s3s;
								SubSeq[Y + off3] = (s3s == 0);
							}
							else{
								X++;

								if (X <= LSeq) {


									

									//for (Z = X; Z <= LSeq; Z++){
									while (SeqNum[X + os1] == 46 || SeqNum[X + os2] == 46 || SeqNum[X + os3] == 46) {//(binarray[se2os] > 1 || binarray[se3os] > 1) {
										X++;
										if (X > LSeq)
											break;

										
									}
									//X=Z;
								}
								
								SubSeq[++Y + off3] = 1;
								
							}
								
						}
						
						


						X--;
						

						mp = (int)(SX + (X - SX) / 2);

						for (Z = SX + 1; Z <= X; Z++)
							XPosDiff[Z] = 0;

						XPosDiff[mp] = Y;
						XDiffPos[Y] = mp;
						
						
					}

					
				}
				else if (ba3 != 0) {

					if (ba3 == 1) {
						
							
								
						XDiffPos[++Y] = X;
						XPosDiff[X] = Y;
						SubSeq[Y] = 1;
						ah0++;
								
					}
					else
						GoOn = 1;

					if (GoOn == 1) {
						XPosDiff[X] = Y;
						SX = X;
						s1s = 1;
						X++;

						if (X <= LSeq) {


							s1 = SeqNum[X + os1];
							s2 = SeqNum[X + os2];
							

							//for (Z = X; Z <= LSeq; Z++){
							while (s1 == 46 || s2 == 46 || SeqNum[X + os3] == 46) {
								//	s1 = SeqNum[Z + os1];
								//	s2 = SeqNum[Z + os2];
								//	s3 = SeqNum[Z + os3];

								//if (s1 != 46 && s2 != 46 && s3 != 46) 
								//	break;	
								if (s1s != 0) {
									if (s1 != s2)
										s1s = 0;
								}
								X++;
								if (X > LSeq)
									break;

								s1 = SeqNum[X + os1];
								s2 = SeqNum[X + os2];
								

							}
							//X=Z;
						}



						X--;
						mp = (int)(SX + (X - SX) / 2);

						for (Z = SX + 1; Z <= X; Z++)
							XPosDiff[Z] = 0;
						
						XDiffPos[++Y] = mp;
						XPosDiff[mp] = Y;
						ah0 += s1s;
						SubSeq[Y] = s1s;
						SubSeq[Y + off3] = (s1s == 0);
						
					}

					
				}
				else


					XPosDiff[X] = Y;

			}
			//for (X = 1; X <= Y; X++) {
			NDiff[0] = ah0; //NDiff[0] + SubSeq[X];
			NDiff[1] = ah1;// NDiff[1] + SubSeq[X + off1];
			NDiff[2] = ah2;// NDiff[2] + SubSeq[X + off2];
			//}
		}
				
		
		else if (gcindelflag == 2) {
			for (X = 1; X <= LSeq; X++) {

				XPosDiff[X] = Y;
				s1 = SeqNum[X + os1];
				s2 = SeqNum[X + os2];
				s3 = SeqNum[X + os3];

				if (s1 != s2 || s1 != s3) {
					if (s1 != 46 && s2 != 46 && s3 != 46) {
						Y++;
						XPosDiff[X] = Y;
						XDiffPos[Y] = X;
						SubSeq[Y] = (s1 == s2);
						SubSeq[Y + LSeq + 1] = (s1 == s3);
						SubSeq[Y + (LSeq + 1) * 2] = (s2 == s3);
					}

					else {
						Y++;
						XPosDiff[X] = Y;
						XDiffPos[Y] = X;
						SubSeq[Y] = (s1 == s2);
						SubSeq[Y + (LSeq + 1)] = (s1 == s3);
						SubSeq[Y + (LSeq + 1) * 2] = (s2 == s3);
					}


					SubSeq[Y + (LSeq + 1) * 6] = (SubSeq[Y] == 0 && SubSeq[Y + (LSeq + 1)] == 0 && SubSeq[Y + (LSeq + 1) * 2] == 0);


				}

			}
			for (X = 1; X <= Y; X++) {
				NDiff[0] = NDiff[0] + SubSeq[X];
				NDiff[1] = NDiff[1] + SubSeq[X + off1];
				NDiff[2] = NDiff[2] + SubSeq[X + off2];
			}

		}
		

		XDiffPos[Y + 1] = LSeq;

		Z = XDiffPos[Y] + 1;

		for (X = Z; X <= LSeq; X++)
			XPosDiff[X] = Y;

		return(Y);
	}



	int MyMathFuncs::FindSubSeqGCAP3(int UBND, int UBXPD1, int UBSS1, int UBSS2, int elementnum, int *lxos, char gcindelflag, int LSeq, int seq1, int seq2, int *elementseq, short int *SeqNum, char *SubSeq, int *XPosDiff, int *XDiffPos, int *NDiff, unsigned char *binarray) {

		int Z, X, Y, SX;
		int MinDiff, MaxDiff, seq3, e, s1, s2, s3, os1, os2, os3, mp, s1s, s2s, s3s, off1, off2, off3, GoOn, target, se3os, se2os, ah0, ah1, ah2, SY, os1x, os3x, os2x, xpdos,  xpdos2, ssos,ssos2, ndos, ndos2;
		unsigned char  ba2, ba3;
		

		
		xpdos = UBXPD1 + 1;
		ssos = (UBSS1 + 1)*(UBSS2 + 1);
		ndos = UBND + 1;
		os1 = seq1 * (LSeq + 1);
		os2 = seq2 * (LSeq + 1);
		
		off1 = LSeq + 1;
		off2 = (LSeq + 1) * 2;
		off3 = (LSeq + 1) * 6;


		//SubSeq(Len(StrainSeq(0)), 6)
		//BinArray(Len(StrainSeq(0)), Nextno)
		target =ssos*(elementnum+1);
		for (X = 1; X < target; X++)
				SubSeq[X] = 0;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private (e, Y, ah0, ah1, ah2, xpdos2, ssos2, ndos2, seq3, os3, X, se2os, se3os, ba2, ba3, s1s, s2s, s3s, mp, SX, Z, s1, s2, s3, MaxDiff, MinDiff)		
		for (e = 0; e <= elementnum; e++) {
			Y = 0;
			ah0 = 0;
			ah1 = 0;
			ah2 = 0;
			xpdos2 = e*xpdos;
			ssos2 = ssos*e;
			ndos2 = ndos*e;
			seq3 = elementseq[e];
			os3 = seq3 * (LSeq + 1);

			X = 0;
			//for (X = 0; X <= LSeq; X++)
			XPosDiff[X + xpdos2] = 0;

			if (gcindelflag == 0) {

				for (X = 1; X <= LSeq; X++) {

					se2os = X + os2;

					if (binarray[se2os] == 1) {

						se3os = X + os3;
						if (binarray[se3os] == 1) {


							XDiffPos[++Y + xpdos2] = X;
							if (SeqNum[se2os] == SeqNum[se3os]) {
								SubSeq[Y + off2 + ssos2] = 1;
								ah2++;
							}
							else
								SubSeq[Y + off3 + ssos2] = 1;
						}
						else if (binarray[se3os] == 0) {
							XDiffPos[++Y + xpdos2] = X;
							SubSeq[Y + off1 + ssos2] = 1;//(s1 == s3);
							ah1++;
						}


					}
					else if (binarray[X + os3] == 1) {
						if (binarray[se2os] == 0) {


							XDiffPos[++Y + xpdos2] = X;
							SubSeq[Y + ssos2] = 1;// (s1 == s2);
							ah0++;
						}

					}
					XPosDiff[X + xpdos2] = Y;
				}
				NDiff[ndos2] = ah0; //NDiff[0] + SubSeq[X];
				NDiff[1 + ndos2] = ah1;// NDiff[1] + SubSeq[X + off1];
				NDiff[2 + ndos2] = ah2;// NDiff[2] + SubSeq[X + off2];

			}
			else if (gcindelflag == 1) {
				for (X = 1; X <= LSeq; X++) {


					ba2 = binarray[X + os2];
					ba3 = binarray[X + os3];
					GoOn = 0;
					if (ba2 != 0) {

						if (ba2 == 1) {
							if (ba3 == 0) {


								XDiffPos[++Y + xpdos2] = X;
								XPosDiff[X + xpdos2] = Y;
								SubSeq[Y + off1 + ssos2] = 1;
								ah1++;

							}
							else if (ba3 == 1) {

								XDiffPos[++Y + xpdos2] = X;
								XPosDiff[X + xpdos2] = Y;
								s2s = (SeqNum[X + os2] == SeqNum[X + os3]);
								SubSeq[Y + off2 + ssos2] = s2s;
								ah2 += s2s;
								SubSeq[Y + off3 + ssos2] = (s2s == 0);
							}
							else
								GoOn = 1;
						}
						else
							GoOn = 1;

						if (GoOn == 1) {
							XPosDiff[X + xpdos2] = Y;
							SX = X;
							if (ba3 == 0) {
								s2s = 1;

								X++;

								if (X <= LSeq) {


									s1 = SeqNum[X + os1];
									s3 = SeqNum[X + os3];

									//for (Z = X; Z <= LSeq; Z++){
									while (s1 == 46 || SeqNum[X + os2] == 46 || s3 == 46) {//(binarray[se2os] > 1 || binarray[se3os] > 1) {
																						   //	s1 = SeqNum[Z + os1];

										if (s2s != 0) {
											if (s1 != s3)
												s2s = 0;
										}


										X++;
										if (X > LSeq)
											break;

										s1 = SeqNum[X + os1];
										s3 = SeqNum[X + os3];

									}
									//X=Z;
								}
								ah1 += s2s;
								SubSeq[++Y + off1 + ssos2] = s2s;
								SubSeq[Y + off3 + ssos2] = (s2s == 0);
							}
							else {

								s3s = (SeqNum[X + os2] == SeqNum[X + os3]);
								if (s3s == 1) {
									X++;

									if (X <= LSeq) {



										s2 = SeqNum[X + os2];
										s3 = SeqNum[X + os3];

										//for (Z = X; Z <= LSeq; Z++){
										while (SeqNum[X + os1] == 46 || s2 == 46 || s3 == 46) {//(binarray[se2os] > 1 || binarray[se3os] > 1) {

											if (s3s != 0) {
												if (s2 != s3)
													s3s = 0;

											}
											X++;
											if (X > LSeq)
												break;


											s2 = SeqNum[X + os2];
											s3 = SeqNum[X + os3];

										}

									}

									ah2 += s3s;


									SubSeq[++Y + off2 + ssos2] = s3s;
									SubSeq[Y + off3 + ssos2] = (s3s == 0);
								}
								else {
									X++;

									if (X <= LSeq) {




										//for (Z = X; Z <= LSeq; Z++){
										while (SeqNum[X + os1] == 46 || SeqNum[X + os2] == 46 || SeqNum[X + os3] == 46) {//(binarray[se2os] > 1 || binarray[se3os] > 1) {
											X++;
											if (X > LSeq)
												break;


										}
										//X=Z;
									}

									SubSeq[++Y + off3 + ssos2] = 1;

								}

							}




							X--;


							mp = (int)(SX + (X - SX) / 2);

							for (Z = SX + 1; Z <= X; Z++)
								XPosDiff[Z + xpdos2] = 0;

							XPosDiff[mp + xpdos2] = Y;
							XDiffPos[Y + xpdos2] = mp;


						}


					}
					else if (ba3 != 0) {

						if (ba3 == 1) {



							XDiffPos[++Y + xpdos2] = X;
							XPosDiff[X + xpdos2] = Y;
							SubSeq[Y + ssos2] = 1;
							ah0++;

						}
						else
							GoOn = 1;

						if (GoOn == 1) {
							XPosDiff[X + xpdos2] = Y;
							SX = X;
							s1s = 1;
							X++;

							if (X <= LSeq) {


								s1 = SeqNum[X + os1];
								s2 = SeqNum[X + os2];


								//for (Z = X; Z <= LSeq; Z++){
								while (s1 == 46 || s2 == 46 || SeqNum[X + os3] == 46) {
									//	s1 = SeqNum[Z + os1];
									//	s2 = SeqNum[Z + os2];
									//	s3 = SeqNum[Z + os3];

									//if (s1 != 46 && s2 != 46 && s3 != 46) 
									//	break;	
									if (s1s != 0) {
										if (s1 != s2)
											s1s = 0;
									}
									X++;
									if (X > LSeq)
										break;

									s1 = SeqNum[X + os1];
									s2 = SeqNum[X + os2];


								}
								//X=Z;
							}



							X--;
							mp = (int)(SX + (X - SX) / 2);

							for (Z = SX + 1; Z <= X; Z++)
								XPosDiff[Z + xpdos2] = 0;

							XDiffPos[++Y + xpdos2] = mp;
							XPosDiff[mp + xpdos2] = Y;
							ah0 += s1s;
							SubSeq[Y + ssos2] = s1s;
							SubSeq[Y + off3 + ssos2] = (s1s == 0);

						}


					}
					else


						XPosDiff[X + xpdos2] = Y;

				}
				//for (X = 1; X <= Y; X++) {
				NDiff[ndos2] = ah0; //NDiff[0] + SubSeq[X];
				NDiff[1 + ndos2] = ah1;// NDiff[1] + SubSeq[X + off1];
				NDiff[2 + ndos2] = ah2;// NDiff[2] + SubSeq[X + off2];
							   //}
			}


			else if (gcindelflag == 2) {
				for (X = 1; X <= LSeq; X++) {

					XPosDiff[X + xpdos2] = Y;
					s1 = SeqNum[X + os1];
					s2 = SeqNum[X + os2];
					s3 = SeqNum[X + os3];

					if (s1 != s2 || s1 != s3) {
						if (s1 != 46 && s2 != 46 && s3 != 46) {
							Y++;
							XPosDiff[X + xpdos2] = Y;
							XDiffPos[Y + xpdos2] = X;
							SubSeq[Y + ssos2] = (s1 == s2);
							SubSeq[Y + LSeq + 1 + ssos2] = (s1 == s3);
							SubSeq[Y + (LSeq + 1) * 2 + ssos2] = (s2 == s3);
						}

						else {
							Y++;
							XPosDiff[X + xpdos2] = Y;
							XDiffPos[Y + xpdos2] = X;
							SubSeq[Y + ssos2] = (s1 == s2);
							SubSeq[Y + (LSeq + 1) + ssos2] = (s1 == s3);
							SubSeq[Y + (LSeq + 1) * 2 + ssos2] = (s2 == s3);
						}


						SubSeq[Y + (LSeq + 1) * 6 + ssos2] = (SubSeq[Y + ssos2] == 0 && SubSeq[Y + (LSeq + 1) + ssos2] == 0 && SubSeq[Y + (LSeq + 1) * 2 + ssos2] == 0);


					}

				}
				for (X = 1; X <= Y; X++) {
					NDiff[ndos2] = NDiff[0 + ndos2] + SubSeq[X + ssos2];
					NDiff[1 + ndos2] = NDiff[1 + ndos2] + SubSeq[X + off1 + ssos2];
					NDiff[2 + ndos2] = NDiff[2 + ndos2] + SubSeq[X + off2 + ssos2];
				}

			}


			XDiffPos[Y + 1 + xpdos2] = LSeq;

			Z = XDiffPos[Y + xpdos2] + 1;

			for (X = Z; X <= LSeq; X++)
				XPosDiff[X + xpdos2] = Y;
			
			if (NDiff[ndos2] == Y || NDiff[1+ndos2] == Y || NDiff[2+ndos2] == Y)
				lxos[e] = 0;
				//for outer frags (ie matches instead of differences)
			else {
				NDiff[3 + ndos2] = NDiff[ndos2] + NDiff[1 + ndos2]; //seq1
				NDiff[4 + ndos2] = NDiff[ndos2] + NDiff[2 + ndos2];  //seq2
				NDiff[5 + ndos2] = NDiff[1 + ndos2] + NDiff[2 + ndos2];  //seq3


					//for inner fargs (ie genuine differences)
				NDiff[ndos2] = Y - NDiff[ndos2];
				NDiff[1 + ndos2] = Y - NDiff[1 + ndos2];
				NDiff[2 + ndos2] = Y - NDiff[2 + ndos2];



				MaxDiff = 0;
				MinDiff = LSeq;
				for (X = 0; X <= 5; X++) {
					if (MinDiff > NDiff[X + ndos2])
						MinDiff = NDiff[X + ndos2];
					if (MaxDiff < NDiff[X + ndos2])
						MaxDiff = NDiff[X + ndos2];

				}

				if (MinDiff < 3 && MaxDiff > MinDiff * 10) 
					lxos[e] = 0;
				else {
					if (NDiff[3 + ndos2] == 0)
						NDiff[3 + ndos2] = 1;
					if (NDiff[4 + ndos2] == 0)
						NDiff[4 + ndos2] = 1;
					if (NDiff[5 + ndos2] == 0)
						NDiff[5 + ndos2] = 1;

					SubSeq[Y + 1 + ssos2] = 0;
					SubSeq[Y + 1 + off1 + ssos2] = 0;
					SubSeq[Y + 1 + off2 + ssos2] = 0;
					SubSeq[Y + 1 + off3 + ssos2] = 0;

					lxos[e] = Y;
				}
					
			}
		}
		omp_set_num_threads(2);
		return(1);
	}


	int MyMathFuncs::GetBestMatch(int Nextno, int NumSeeds, int UBD,float *Dist, int *BestMatch) {
		int x, Y, off;
		float Lowest;
		off = UBD + 1;
		for (x = 0; x <= Nextno; x++) {
			Lowest = (float)(Nextno);
				
			for (Y = 0; Y <= NumSeeds; Y++) {
				if (Lowest > Dist[x + Y*off]) {
					Lowest = Dist[x + Y*off];
					BestMatch[x] = Y;

				}
			}
		}
		return(1);

	}

	int MyMathFuncs::GetBestMatch2(int Nextno, int NumSeeds, int UBD, float *Dist, int *BestMatch, int *NIY) {
		int x, Y, off;
		float Lowest;
		off = UBD + 1;
		for (x = 0; x <= Nextno; x++) {
			Lowest = (float)(Nextno);

			for (Y = 0; Y <= NumSeeds; Y++) {
				if (Lowest > Dist[x + Y*off]) {
					Lowest = Dist[x + Y*off];
					BestMatch[x] = Y;
				}
				else if (Lowest == Dist[x + Y*off]) {
					if (NIY[BestMatch[x]] > NIY[Y]) {
						//NIY[BestMatch[x]] = NIY[BestMatch[x]] - 1;
						BestMatch[x] = Y;
					}
				}
			}
			NIY[BestMatch[x]] = NIY[BestMatch[x]] + 1;
		}
		return(1);

	}

	int MyMathFuncs::GetClosestTo(int A, int Nextno, int UBD,int *Done,float *ClosestTo, float *Dist) {
		int Y, x,off;
		off = UBD + 1;
		for (Y = 0; Y <= Nextno; Y++) {
			if (Done[Y] == 0) {
				
				for (x = 0; x < A; x++) {
					if (Dist[Y + x*off] < ClosestTo[Y])
						ClosestTo[Y] = Dist[Y + x*off];
				}

			}
		}
		return(1);
	}







	int MyMathFuncs::FindSubSeqGCAP4(int UBND, int UBXPD1, int UBSS1, int UBSS2, int elementnum, int *lxos, char gcindelflag, int LSeq, int seq1, int *elementseq2, int *elementseq, short int *SeqNum, char *SubSeq, int *XPosDiff, int *XDiffPos, int *NDiff, unsigned char *binarray) {

		int Z, X, Y, SX;
		int MinDiff, MaxDiff, seq2, seq3, e, s1, s2, s3, os1, os2, os3, mp, s1s, s2s, s3s, off1, off2, off3, GoOn, target, se3os, se2os, ah0, ah1, ah2, SY, os1x, os3x, os2x, xpdos, xpdos2, ssos, ssos2, ndos, ndos2;
		unsigned char  ba2, ba3;



		xpdos = UBXPD1 + 1;
		ssos = (UBSS1 + 1)*(UBSS2 + 1);
		ndos = UBND + 1;
		os1 = seq1 * (LSeq + 1);
		
		off1 = LSeq + 1;
		off2 = (LSeq + 1) * 2;
		off3 = (LSeq + 1) * 6;


		//SubSeq(Len(StrainSeq(0)), 6)
		//BinArray(Len(StrainSeq(0)), Nextno)
		target = ssos*(elementnum + 1);
		for (X = 1; X < target; X++)
			SubSeq[X] = 0;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private (e, Y, ah0, ah1, ah2, xpdos2, ssos2, ndos2, seq2, seq3, os2, os3, X, se2os, se3os, ba2, ba3, s1s, s2s, s3s, mp, SX, Z, s1, s2, s3, MaxDiff, MinDiff)		
		for (e = 0; e <= elementnum; e++) {
			Y = 0;
			ah0 = 0;
			ah1 = 0;
			ah2 = 0;
			xpdos2 = e*xpdos;
			ssos2 = ssos*e;
			ndos2 = ndos*e;


			seq2 = elementseq2[e];
			os2 = seq2 * (LSeq + 1);
			seq3 = elementseq[e];
			os3 = seq3 * (LSeq + 1);

			X = 0;
			//for (X = 0; X <= LSeq; X++)
			XPosDiff[X + xpdos2] = 0;

			if (gcindelflag == 0) {

				for (X = 1; X <= LSeq; X++) {

					se2os = X + os2;

					if (binarray[se2os] == 1) {

						se3os = X + os3;
						if (binarray[se3os] == 1) {


							XDiffPos[++Y + xpdos2] = X;
							if (SeqNum[se2os] == SeqNum[se3os]) {
								SubSeq[Y + off2 + ssos2] = 1;
								ah2++;
							}
							else
								SubSeq[Y + off3 + ssos2] = 1;
						}
						else if (binarray[se3os] == 0) {
							XDiffPos[++Y + xpdos2] = X;
							SubSeq[Y + off1 + ssos2] = 1;//(s1 == s3);
							ah1++;
						}


					}
					else if (binarray[X + os3] == 1) {
						if (binarray[se2os] == 0) {


							XDiffPos[++Y + xpdos2] = X;
							SubSeq[Y + ssos2] = 1;// (s1 == s2);
							ah0++;
						}

					}
					XPosDiff[X + xpdos2] = Y;
				}
				NDiff[ndos2] = ah0; //NDiff[0] + SubSeq[X];
				NDiff[1 + ndos2] = ah1;// NDiff[1] + SubSeq[X + off1];
				NDiff[2 + ndos2] = ah2;// NDiff[2] + SubSeq[X + off2];

			}
			else if (gcindelflag == 1) {
				for (X = 1; X <= LSeq; X++) {


					ba2 = binarray[X + os2];
					ba3 = binarray[X + os3];
					GoOn = 0;
					if (ba2 != 0) {

						if (ba2 == 1) {
							if (ba3 == 0) {


								XDiffPos[++Y + xpdos2] = X;
								XPosDiff[X + xpdos2] = Y;
								SubSeq[Y + off1 + ssos2] = 1;
								ah1++;

							}
							else if (ba3 == 1) {

								XDiffPos[++Y + xpdos2] = X;
								XPosDiff[X + xpdos2] = Y;
								s2s = (SeqNum[X + os2] == SeqNum[X + os3]);
								SubSeq[Y + off2 + ssos2] = s2s;
								ah2 += s2s;
								SubSeq[Y + off3 + ssos2] = (s2s == 0);
							}
							else
								GoOn = 1;
						}
						else
							GoOn = 1;

						if (GoOn == 1) {
							XPosDiff[X + xpdos2] = Y;
							SX = X;
							if (ba3 == 0) {
								s2s = 1;

								X++;

								if (X <= LSeq) {


									s1 = SeqNum[X + os1];
									s3 = SeqNum[X + os3];

									//for (Z = X; Z <= LSeq; Z++){
									while (s1 == 46 || SeqNum[X + os2] == 46 || s3 == 46) {//(binarray[se2os] > 1 || binarray[se3os] > 1) {
																						   //	s1 = SeqNum[Z + os1];

										if (s2s != 0) {
											if (s1 != s3)
												s2s = 0;
										}


										X++;
										if (X > LSeq)
											break;

										s1 = SeqNum[X + os1];
										s3 = SeqNum[X + os3];

									}
									//X=Z;
								}
								ah1 += s2s;
								SubSeq[++Y + off1 + ssos2] = s2s;
								SubSeq[Y + off3 + ssos2] = (s2s == 0);
							}
							else {

								s3s = (SeqNum[X + os2] == SeqNum[X + os3]);
								if (s3s == 1) {
									X++;

									if (X <= LSeq) {



										s2 = SeqNum[X + os2];
										s3 = SeqNum[X + os3];

										//for (Z = X; Z <= LSeq; Z++){
										while (SeqNum[X + os1] == 46 || s2 == 46 || s3 == 46) {//(binarray[se2os] > 1 || binarray[se3os] > 1) {

											if (s3s != 0) {
												if (s2 != s3)
													s3s = 0;

											}
											X++;
											if (X > LSeq)
												break;


											s2 = SeqNum[X + os2];
											s3 = SeqNum[X + os3];

										}

									}

									ah2 += s3s;


									SubSeq[++Y + off2 + ssos2] = s3s;
									SubSeq[Y + off3 + ssos2] = (s3s == 0);
								}
								else {
									X++;

									if (X <= LSeq) {




										//for (Z = X; Z <= LSeq; Z++){
										while (SeqNum[X + os1] == 46 || SeqNum[X + os2] == 46 || SeqNum[X + os3] == 46) {//(binarray[se2os] > 1 || binarray[se3os] > 1) {
											X++;
											if (X > LSeq)
												break;


										}
										//X=Z;
									}

									SubSeq[++Y + off3 + ssos2] = 1;

								}

							}




							X--;


							mp = (int)(SX + (X - SX) / 2);

							for (Z = SX + 1; Z <= X; Z++)
								XPosDiff[Z + xpdos2] = 0;

							XPosDiff[mp + xpdos2] = Y;
							XDiffPos[Y + xpdos2] = mp;


						}


					}
					else if (ba3 != 0) {

						if (ba3 == 1) {



							XDiffPos[++Y + xpdos2] = X;
							XPosDiff[X + xpdos2] = Y;
							SubSeq[Y + ssos2] = 1;
							ah0++;

						}
						else
							GoOn = 1;

						if (GoOn == 1) {
							XPosDiff[X + xpdos2] = Y;
							SX = X;
							s1s = 1;
							X++;

							if (X <= LSeq) {


								s1 = SeqNum[X + os1];
								s2 = SeqNum[X + os2];


								//for (Z = X; Z <= LSeq; Z++){
								while (s1 == 46 || s2 == 46 || SeqNum[X + os3] == 46) {
									//	s1 = SeqNum[Z + os1];
									//	s2 = SeqNum[Z + os2];
									//	s3 = SeqNum[Z + os3];

									//if (s1 != 46 && s2 != 46 && s3 != 46) 
									//	break;	
									if (s1s != 0) {
										if (s1 != s2)
											s1s = 0;
									}
									X++;
									if (X > LSeq)
										break;

									s1 = SeqNum[X + os1];
									s2 = SeqNum[X + os2];


								}
								//X=Z;
							}



							X--;
							mp = (int)(SX + (X - SX) / 2);

							for (Z = SX + 1; Z <= X; Z++)
								XPosDiff[Z + xpdos2] = 0;

							XDiffPos[++Y + xpdos2] = mp;
							XPosDiff[mp + xpdos2] = Y;
							ah0 += s1s;
							SubSeq[Y + ssos2] = s1s;
							SubSeq[Y + off3 + ssos2] = (s1s == 0);

						}


					}
					else


						XPosDiff[X + xpdos2] = Y;

				}
				//for (X = 1; X <= Y; X++) {
				NDiff[ndos2] = ah0; //NDiff[0] + SubSeq[X];
				NDiff[1 + ndos2] = ah1;// NDiff[1] + SubSeq[X + off1];
				NDiff[2 + ndos2] = ah2;// NDiff[2] + SubSeq[X + off2];
									   //}
			}


			else if (gcindelflag == 2) {
				for (X = 1; X <= LSeq; X++) {

					XPosDiff[X + xpdos2] = Y;
					s1 = SeqNum[X + os1];
					s2 = SeqNum[X + os2];
					s3 = SeqNum[X + os3];

					if (s1 != s2 || s1 != s3) {
						if (s1 != 46 && s2 != 46 && s3 != 46) {
							Y++;
							XPosDiff[X + xpdos2] = Y;
							XDiffPos[Y + xpdos2] = X;
							SubSeq[Y + ssos2] = (s1 == s2);
							SubSeq[Y + LSeq + 1 + ssos2] = (s1 == s3);
							SubSeq[Y + (LSeq + 1) * 2 + ssos2] = (s2 == s3);
						}

						else {
							Y++;
							XPosDiff[X + xpdos2] = Y;
							XDiffPos[Y + xpdos2] = X;
							SubSeq[Y + ssos2] = (s1 == s2);
							SubSeq[Y + (LSeq + 1) + ssos2] = (s1 == s3);
							SubSeq[Y + (LSeq + 1) * 2 + ssos2] = (s2 == s3);
						}


						SubSeq[Y + (LSeq + 1) * 6 + ssos2] = (SubSeq[Y + ssos2] == 0 && SubSeq[Y + (LSeq + 1) + ssos2] == 0 && SubSeq[Y + (LSeq + 1) * 2 + ssos2] == 0);


					}

				}
				for (X = 1; X <= Y; X++) {
					NDiff[ndos2] = NDiff[0 + ndos2] + SubSeq[X + ssos2];
					NDiff[1 + ndos2] = NDiff[1 + ndos2] + SubSeq[X + off1 + ssos2];
					NDiff[2 + ndos2] = NDiff[2 + ndos2] + SubSeq[X + off2 + ssos2];
				}

			}


			XDiffPos[Y + 1 + xpdos2] = LSeq;

			Z = XDiffPos[Y + xpdos2] + 1;

			for (X = Z; X <= LSeq; X++)
				XPosDiff[X + xpdos2] = Y;

			if (NDiff[ndos2] == Y || NDiff[1 + ndos2] == Y || NDiff[2 + ndos2] == Y)
				lxos[e] = 0;
			//for outer frags (ie matches instead of differences)
			else {
				NDiff[3 + ndos2] = NDiff[ndos2] + NDiff[1 + ndos2]; //seq1
				NDiff[4 + ndos2] = NDiff[ndos2] + NDiff[2 + ndos2];  //seq2
				NDiff[5 + ndos2] = NDiff[1 + ndos2] + NDiff[2 + ndos2];  //seq3


																		 //for inner fargs (ie genuine differences)
				NDiff[ndos2] = Y - NDiff[ndos2];
				NDiff[1 + ndos2] = Y - NDiff[1 + ndos2];
				NDiff[2 + ndos2] = Y - NDiff[2 + ndos2];



				MaxDiff = 0;
				MinDiff = LSeq;
				for (X = 0; X <= 5; X++) {
					if (MinDiff > NDiff[X + ndos2])
						MinDiff = NDiff[X + ndos2];
					if (MaxDiff < NDiff[X + ndos2])
						MaxDiff = NDiff[X + ndos2];

				}

				if (MinDiff < 3 && MaxDiff > MinDiff * 10)
					lxos[e] = 0;
				else {
					if (NDiff[3 + ndos2] == 0)
						NDiff[3 + ndos2] = 1;
					if (NDiff[4 + ndos2] == 0)
						NDiff[4 + ndos2] = 1;
					if (NDiff[5 + ndos2] == 0)
						NDiff[5 + ndos2] = 1;

					SubSeq[Y + 1 + ssos2] = 0;
					SubSeq[Y + 1 + off1 + ssos2] = 0;
					SubSeq[Y + 1 + off2 + ssos2] = 0;
					SubSeq[Y + 1 + off3 + ssos2] = 0;

					lxos[e] = Y;
				}

			}
		}
		omp_set_num_threads(2);
		return(1);
	}

	int MyMathFuncs::CalcKMaxP(short int GCMissmatchPen, int XOLen, short int MCFlag, int MCCorrection, double LowestProb, double *pco, int *HiFragScore, double *critval, double *MissPen, double *LL, double *KMax, int *NDiff, int *highenough) {

		int GoOn, X;
		double PVM1, KMTL, P, Q, LenXoverSeq, mP, k, ll0, Z, zm, Yy, zdel, mx, d1, d2, d3, PCO;

		LenXoverSeq = (double)(XOLen);
		for (X = 0; X < 6; X++) {
			LL[X] = 0;
			if (highenough[X] == 1) {
				if (NDiff[X] > 0 && NDiff[X] < LenXoverSeq) {

					P = NDiff[X] / LenXoverSeq;
					Q = 1 - P;

					/* If the mismatch penalty is infinite: */
					if (GCMissmatchPen == 0) {
						LL[X] = -log(Q);
						KMax[X] = P;
					}
					else {
						mP = MissPen[X] * P;
						mx = MissPen[X];
						k = 0.0;
						ll0 = log(mP / Q);
						ll0 = ll0 / (mx + 1);
						Z = exp(2 * ll0);
						zdel = 1;
						Yy = 1;
						while (fabs(zdel) > 0.000001 || fabs(Yy) > 0.000001) {

							zm = pow(Z, -mx);
							Yy = Q * Z + P * zm - 1;
							zdel = Yy / (Q - mP * zm / Z);
							Z = Z - zdel;

						}

						LL[X] = log(Z);
						d1 = exp(LL[X]);
						d1 = d1 - 1;
						d2 = -(mx + 1) * LL[X];
						d3 = exp(d2);
						KMax[X] = d1 * (Q - (mP * d3));

					}
				}
			}
		}




			
		
		if (MCFlag == 0)
			PCO = LowestProb / MCCorrection;
		else
			PCO = LowestProb;
			

		if (PCO > 1)
			return(0);
			



		for (X = 0; X <= 5; X++) {
			if (KMax[X] > 0) {
				KMTL = KMax[X] * XOLen;
				KMTL = log(KMTL);
				PVM1 = (1 - PCO);

				if (PVM1 > 0) {
					PVM1 = -log(PVM1);
					if (PVM1 > 0) {
						PVM1 = -log(PVM1);
						critval[X] = (KMTL + PVM1) / LL[X];
						if (critval[X] < 4)
							critval[X] = 4;
					}
				}
				else 
					critval[X] = 4;
				
			}
			else
				critval[X] = 4;
					
			
		}

		GoOn = 0;
		for (X = 0; X <= 5; X++){
			highenough[X] = 0;
			if (HiFragScore[X] > 3 && HiFragScore[X] > critval[X]){
				highenough[X] = 1;
				GoOn = 1;
			}
		}
			
		*pco = PCO;

		return(GoOn);

	}

	int MyMathFuncs::GetMaxFragScoreP(int LenXoverSeq, int lseq, short int CircularFlag, short int GCMissmatchPen, double *MissPen, int *MaxScorePos, int *FragMaxScore, int *FragScore, int *FragCount, int *hiscore) {

		int fcx, fms, X, Y, Z, os, os2, Polys, Diffs, fs, ts, os3, os5, os6, os7, hsx, msp;
		float MPen;

		os = lseq + 1;
		os2 = lseq + 1;

		if (GCMissmatchPen > 0) {

			//Find the max score per frag by joing up frags, basically start at leftmost frag and
			//expand right, finding the max score including that frag -these are inner frags
			
			//this works terribly with 6 threads but ok with 3
			/*int procs;
			procs = omp_get_num_procs();
			procs = procs/2 - 1;
			if (procs < 3)
				procs = 3;
			omp_set_num_threads(procs);*/
			/*int procs;
			procs = omp_get_num_procs();
			if (procs > 1)
				omp_set_num_threads(2);
			else
				omp_set_num_threads(1);*/
			if (CircularFlag == 0) {
				
//

//#pragma omp parallel for private(X, hsx, os3, os6, MPen, fcx, Y, os5, Polys, fms, Diffs,Z, ts, msp, fs)
				for (X = 0; X < 6; X++) {
					hsx = 0;
					os3 = X*os2;
					os6 = X*os;
					MPen = (float)(MissPen[X]);
					fcx = FragCount[X];

					for (Y = 0; Y <= fcx; Y++) {
						os5 = Y + os3;

						if (FragScore[os5] > 0) {

							
							Polys = FragScore[os5];
							fms = Polys;
							Diffs = 0;

							
							msp = Y;

							for (Z = Y + 1; Z <= fcx; Z++) {


								fs = FragScore[Z + os3];
								if (fs <= 0) {
									fs = -fs;
									Diffs += fs;
								}

								Polys += fs;

								ts = (int)((Polys - Diffs) - (Diffs * MPen));
								if (ts < 0)
									break;
								else if (ts >= fms) {
									fms = ts;
									 msp = Z;
								}


							}
							MaxScorePos[os5] = msp;
							FragMaxScore[os5] = fms;
							if (fms > hsx) 
								hsx = fms;

							
							
						}
						else
							FragMaxScore[os5] = 0;

					}
					hiscore[X] = hsx;
				}
			}
			else {

//#pragma omp parallel for private(X, hsx, os3, MPen, fcx, Y, os5, Polys, fms, Diffs,Z, ts, msp, fs)		

				for (X = 0; X < 6; X++) {
					hsx = 0;
					os3 = X*os2;
					MPen = (float)(MissPen[X]);
					fcx = FragCount[X];

					for (Y = 0; Y <= fcx; Y++) {
						os5 = Y + os3;

						if (FragScore[os5] > 0) {
							Polys = FragScore[os5];
							Diffs = 0;
							fms = Polys;
							
							msp = Y;
							for (Z = Y + 1; Z <= fcx; Z++) {


								fs = FragScore[Z + os3];
								if (fs <= 0) {
									fs = -fs;
									Diffs += fs;
								}

								Polys += fs;

								ts = (int)((Polys - Diffs) - (Diffs * MPen));
								if (ts < 0)
									break;
								else if (ts >= fms) {
									fms = ts;
									msp = Z;
								}


							}
							MaxScorePos[os5] = msp;
							FragMaxScore[os5] = fms;
							if (fms > hsx) 
								hsx = fms;

								
						}
						else
							FragMaxScore[os5] = 0;


					}
					hiscore[X] = hsx;
				}


			}
		}
		else {
			for (X = 0; X <= 5; X++) {
				os3 = X*os2;
				os7 = X*os;

				for (Y = 0; Y <= FragCount[X]; Y++) {
					os = Y + os3;
					FragMaxScore[Y + os7] = FragScore[os];
					MaxScorePos[os] = Y;
				}
			}
		}
		//omp_set_num_threads(2);
		return (1);
	}


	int MyMathFuncs::GetMaxFragScoreP2(int elementnum, int *LenXoverSeq, int lseq, short int CircularFlag, short int GCMissmatchPen, double *MissPen, int *MaxScorePos, int *FragMaxScore, int *FragScore, int *FragCount, int *hiscore,  int *NDiffG) {

		int osfs, osfs2, e, osfc2,  osfms, osfms2, osmsp, osmsp2, osmp2, fcx, fms, X, Y, Z, os, os2, Polys, Diffs, fs, ts, os3, os5, os6, os7, hsx, msp;
		float MPen, LTG;
		//ReDim HiFragScoreG(5, GroupSize), FragMaxScoreG(GCDimSize, 5, GroupSize), MaxScorePosG(GCDimSize, 5, GroupSize), MissPenG(5, GroupSize), LTGG(GroupSize), NDiffG(6, GroupSize),FragScoreG(GCDimSize, 6, GroupSize), FragCountG(6, GroupSize),
		os = lseq + 1;
		os2 = lseq + 1;
		
		osfms = (lseq + 1) * 6;//
		osmsp = (lseq + 1) * 6;//
		
		osfs = (lseq + 1) * 7;//

		int procs;
		procs = omp_get_num_procs();
		if (procs > 1)
			omp_set_num_threads(2);
		else
			omp_set_num_threads(1);

#pragma omp parallel for private(e, osmp2, osfc2, osfms2, osmsp2, osfs2, LTG, X, hsx, os3, os6, MPen, fcx, Y, os5, Polys, fms, Diffs,Z, ts, msp, fs)
							  //#pragma omp parallel for private(X, hsx, os3,      MPen, fcx, Y, os5, Polys, fms, Diffs,Z, ts, msp, fs)		
		for (e = 0; e <= elementnum; e++) {
			osmp2 = e * 6;
			
			osfc2 = 7 * e;
			
			
			osfms2 = osfms*e;
			osmsp2 = osmsp*e;
			osfs2 = osfs*e;
			LTG = (float)(LenXoverSeq[e] * GCMissmatchPen);
			MissPen[osmp2] = (int)(LTG / NDiffG[osfc2]) + 1;
			MissPen[1 + osmp2] = (int)(LTG / NDiffG[1 + osfc2]) + 1;
			MissPen[2 + osmp2] = (int)(LTG / NDiffG[2 + osfc2]) + 1;
			MissPen[3 + osmp2] = (int)(LTG / NDiffG[3 + osfc2]) + 1;
			MissPen[4 + osmp2] = (int)(LTG / NDiffG[4 + osfc2]) + 1;
			MissPen[5 + osmp2] = (int)(LTG / NDiffG[5 + osfc2]) + 1;

			hiscore[osmp2] = 0;
			hiscore[1 + osmp2] = 0;
			hiscore[2 + osmp2] = 0;
			hiscore[3 + osmp2] = 0;
			hiscore[4 + osmp2] = 0;
			hiscore[5 + osmp2] = 0;



				if (GCMissmatchPen > 0) {

					//Find the max score per frag by joing up frags, basically start at leftmost frag and
					//expand right, finding the max score including that frag -these are inner frags

					//this works terribly with 6 threads but ok with 3
					int procs;
					procs = omp_get_num_procs();
					if (procs > 1)
						omp_set_num_threads(2);
					else
						omp_set_num_threads(1);

					if (CircularFlag == 0) {

						//


						for (X = 0; X < 6; X++) {
							hsx = 0;
							os3 = X*os2;
							os6 = X*os;
							MPen = (float)(MissPen[X + osmp2]);
							fcx = FragCount[X + osfc2];

							for (Y = 0; Y <= fcx; Y++) {
								os5 = Y + os3;

								if (FragScore[os5 + osfs2] > 0) {


									Polys = FragScore[os5 + osfs2];
									fms = Polys;
									Diffs = 0;


									msp = Y;

									for (Z = Y + 1; Z <= fcx; Z++) {


										fs = FragScore[Z + os3 + osfs2];
										if (fs <= 0) {
											fs = -fs;
											Diffs += fs;
										}

										Polys += fs;

										ts = (int)((Polys - Diffs) - (Diffs * MPen));
										if (ts < 0)
											break;
										else if (ts >= fms) {
											fms = ts;
											msp = Z;
										}


									}
									MaxScorePos[os5 + osmsp2] = msp;
									FragMaxScore[os5 + osfms2] = fms;
									if (fms > hsx)
										hsx = fms;



								}
								else
									FragMaxScore[os5 + osfms2] = 0;

							}
							hiscore[X+osmp2] = hsx;
						}
					}
					else {



						for (X = 0; X < 6; X++) {
							hsx = 0;
							os3 = X*os2;
							MPen = (float)(MissPen[X + osmp2]);
							fcx = FragCount[X + osfc2];

							for (Y = 0; Y <= fcx; Y++) {
								os5 = Y + os3;

								if (FragScore[os5 + osfs2] > 0) {
									Polys = FragScore[os5 + osfs2];
									Diffs = 0;
									fms = Polys;

									msp = Y;
									for (Z = Y + 1; Z <= fcx; Z++) {


										fs = FragScore[Z + os3 + osfs2];
										if (fs <= 0) {
											fs = -fs;
											Diffs += fs;
										}

										Polys += fs;

										ts = (int)((Polys - Diffs) - (Diffs * MPen));
										if (ts < 0)
											break;
										else if (ts >= fms) {
											fms = ts;
											msp = Z;
										}


									}
									MaxScorePos[os5 + osmsp2] = msp;
									FragMaxScore[os5 + osfms2] = fms;
									if (fms > hsx)
										hsx = fms;


								}
								else
									FragMaxScore[os5 + osfms2] = 0;


							}
							hiscore[X + osmp2] = hsx;
						}


					}
				}
				else {
					for (X = 0; X <= 5; X++) {
						os3 = X*os2;
						os7 = X*os;

						for (Y = 0; Y <= FragCount[X + osfc2]; Y++) {
							os = Y + os3;
							FragMaxScore[Y + os7 + osfms2] = FragScore[os + osfs2];
							MaxScorePos[os + osmsp2] = Y;
						}
					}
				}
		}
		//omp_set_num_threads(2);
		return (1);
	}


	int MyMathFuncs::XOHomologyP(short int inlyer, int lenstrainseq, int lenxoverseq, short int xoverwindow, char *xoverseqnumw, int *xoverhomologynum)
	{
		int limit = xoverwindow * 2 + 1;
		int limit2 = limit * 2;
		int lenstrainseq2, xoverwindow2, xoverwindow4, lenxoverseq2, ll1, ll2, off3, off4;
		int g = 0;
		int  z, x, off1, off2, h, of1, of2, of3, t1, t2, t3;



		//XOverHomologyNum(LenXOverSeq + XoverWindow * 2, 2)
		//empty some space in xoverhomologynum
		t1 = (1 + lenxoverseq + xoverwindow * 2)*3;
		for (x = 0; x <= t1; x++)
			xoverhomologynum[x] = 0;

		lenstrainseq2 = lenstrainseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		lenxoverseq2 = lenxoverseq * 2;

		off1 = lenstrainseq + xoverwindow2;
		off2 = lenstrainseq2 + xoverwindow4;

		t1 = 0;
		t2 = 0;
		t3 = 0;


		for (z = 1; z <= limit; z++) {
			t1 += *(xoverseqnumw + z);
			t2 += *(xoverseqnumw + z + off1);
			t3 += *(xoverseqnumw + z + off2);
			
		}



		ll1 = lenstrainseq + limit-1;
		ll2 = lenstrainseq2 + limit2-2;

		*(xoverhomologynum + 1) = t1;
		*(xoverhomologynum + 1 + ll1) = t2;
		*(xoverhomologynum + 1 + ll2) = t3;


		off3 = off1 + xoverwindow2;
		off4 = off2 + xoverwindow2;

		off1--;
		off2--;

//#pragma omp parallel num_threads(3)
//		{
//
//#pragma omp sections private(h, of1, of2, of3, x)
//		{
//#pragma omp section
//		{



			h = t1;
			of1 = -1;
			of2 = xoverwindow2;
			of3 = 0;


			for (x = 2; x <= lenxoverseq + 1; x++) {

				h = h - xoverseqnumw[x + of1] + xoverseqnumw[x + of2];
				xoverhomologynum[x + of3] = h;
			}

//		}
//#pragma omp section
//		{



			h = t2;
			of1 = off1;
			of2 = off3;
			of3 = ll1;


			for (x = 2; x <= lenxoverseq + 1; x++) {

				h = h - xoverseqnumw[x + of1] + xoverseqnumw[x + of2];
				xoverhomologynum[x + of3] = h;
			}

//		}
//#pragma omp section
//		{



			h = t3;
			of1 = off2;
			of2 = off4;
			of3 = ll2;


			for (x = 2; x <= lenxoverseq + 1; x++) {

				h = h - xoverseqnumw[x + of1] + xoverseqnumw[x + of2];
				xoverhomologynum[x + of3] = h;
			}

//		}
//		}

//		}




		return(1);
	}


	int MyMathFuncs::XOHomologyP2(short int inlyer, int lenstrainseq, int lenxoverseq, short int xoverwindow, char *xoverseqnumw, int *xoverhomologynum)
	{
		int limit = xoverwindow * 2 + 1;
		int limit2 = limit * 2;
		int lenstrainseq2, xoverwindow2, xoverwindow4, lenxoverseq2, ll1, ll2, off3, off4, thresh;
		int g = 0;
		int  z, x, off1, off2, h, of1, of2, of3, t1, t2, t3;

		thresh = 1;

		//XOverHomologyNum(LenXOverSeq + XoverWindow * 2, 2)
		//empty some space in xoverhomologynum
		t1 = (1 + lenxoverseq + xoverwindow * 2) * 3;
		for (x = 0; x <= t1; x++)
			xoverhomologynum[x] = 0;

		lenstrainseq2 = lenstrainseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		lenxoverseq2 = lenxoverseq * 2;

		off1 = lenstrainseq + xoverwindow2;
		off2 = lenstrainseq2 + xoverwindow4;

		t1 = 0;
		t2 = 0;
		t3 = 0;


		for (z = 1; z <= limit; z++) {
			t1 += *(xoverseqnumw + z);
			t2 += *(xoverseqnumw + z + off1);
			t3 += *(xoverseqnumw + z + off2);

		}



		ll1 = lenstrainseq + limit - 1;
		ll2 = lenstrainseq2 + limit2 - 2;

		*(xoverhomologynum + 1) = t1;
		*(xoverhomologynum + 1 + ll1) = t2;
		*(xoverhomologynum + 1 + ll2) = t3;


		off3 = off1 + xoverwindow2;
		off4 = off2 + xoverwindow2;

		off1--;
		off2--;

		//#pragma omp parallel num_threads(3)
		//		{
		//
		//#pragma omp sections private(h, of1, of2, of3, x)
		//		{
		//#pragma omp section
		//		{



		h = t1;
		of1 = -1;
		of2 = xoverwindow2;
		of3 = 0;


		for (x = 2; x <= lenxoverseq + 1; x++) {

			h = h - xoverseqnumw[x + of1] + xoverseqnumw[x + of2];
			xoverhomologynum[x + of3] = h;
		}

		//		}
		//#pragma omp section
		//		{



		h = t2;
		of1 = off1;
		of2 = off3;
		of3 = ll1;


		for (x = 2; x <= lenxoverseq + 1; x++) {

			h = h - xoverseqnumw[x + of1] + xoverseqnumw[x + of2];
			xoverhomologynum[x + of3] = h;
		}

		//		}
		//#pragma omp section
		//		{



		h = t3;
		of1 = off2;
		of2 = off4;
		of3 = ll2;


		for (x = 2; x <= lenxoverseq + 1; x++) {

			h = h - xoverseqnumw[x + of1] + xoverseqnumw[x + of2];
			xoverhomologynum[x + of3] = h;
		}

		//		}
		//		}

		//		}
		of1 = (off1+1)*(inlyer-1);
		thresh = 0;
		for (x = 1; x <= lenxoverseq + 1; x++) {
			//t1 = xoverhomologynum[x];
			//t2 = xoverhomologynum[x + off1+1];
			//t3 = xoverhomologynum[x + (off1+1)*2];
			if (xoverhomologynum[x + of1] <= xoverwindow) {
				thresh = x;
				break;
			}
		}


		return(thresh);
	}


	int MyMathFuncs::XOHomologyPB5(short int inlyer, int lenstrainseq, int lenxoverseq, short int xoverwindow, char *xoverseqnumw, int *xoverhomologynum)
	{
		int limit = xoverwindow * 2 + 1;
		int limit2 = limit * 2;
		int lenstrainseq2, xoverwindow2, xoverwindow4, lenxoverseq2, ll1, ll2, off3, off4;
		int g = 0;
		int  z, x, off1, off2, h, of1, of2, of3, t1, t2, t3;



		//XOverHomologyNum(LenXOverSeq + XoverWindow * 2, 2)
		//empty some space in xoverhomologynum
		t1 = (1 + lenxoverseq + xoverwindow * 2) * 3;
		for (x = 0; x <= t1; x++)
			xoverhomologynum[x] = 0;

		lenstrainseq2 = lenstrainseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		lenxoverseq2 = lenxoverseq * 2;

		off1 = lenstrainseq + xoverwindow2;
		off2 = lenstrainseq2 + xoverwindow4;

		t1 = 0;
		t2 = 0;
		t3 = 0;


		for (z = 1; z <= limit; z++) {
			t1 += *(xoverseqnumw + z);
			t2 += *(xoverseqnumw + z + off1);
			t3 += *(xoverseqnumw + z + off2);

		}



		ll1 = lenstrainseq + limit - 1;
		ll2 = lenstrainseq2 + limit2 - 2;

		*(xoverhomologynum + 1) = t1;
		*(xoverhomologynum + 1 + ll1) = t2;
		*(xoverhomologynum + 1 + ll2) = t3;


		off3 = off1 + xoverwindow2;
		off4 = off2 + xoverwindow2;

		off1--;
		off2--;

#pragma omp parallel num_threads(3)
		{

#pragma omp sections private(h, of1, of2, of3, x)
			{
#pragma omp section
				{



					h = t1;
					of1 = -1;
					of2 = xoverwindow2;
					of3 = 0;


					for (x = 2; x <= lenxoverseq + 1; x++) {

						h = h - xoverseqnumw[x + of1] + xoverseqnumw[x + of2];
						xoverhomologynum[x + of3] = h;
					}

				}
#pragma omp section
				{



					h = t2;
					of1 = off1;
					of2 = off3;
					of3 = ll1;


					for (x = 2; x <= lenxoverseq + 1; x++) {

						h = h - xoverseqnumw[x + of1] + xoverseqnumw[x + of2];
						xoverhomologynum[x + of3] = h;
					}

				}
#pragma omp section
				{



					h = t3;
					of1 = off2;
					of2 = off4;
					of3 = ll2;


					for (x = 2; x <= lenxoverseq + 1; x++) {

						h = h - xoverseqnumw[x + of1] + xoverseqnumw[x + of2];
						xoverhomologynum[x + of3] = h;
					}

				}
			}

		}




		return(1);
	}
	


	int MyMathFuncs::MakeImageDataP(int bkr, int bkg, int bkb, int SX, int SY, int PosE1, int PosE0, int PosS1, int PosS0, int StS, int StSX, int CurScale, float XAD, float Min, float MR, int UBID1, int UBID2, int UBID3, int UBRM1, int UBHM1, int *HeatMap, float *RegionMat, unsigned char *ImageData) {

		int Y, X, YP, XP, R, G, B, N, M, Band, os1, os2, os3, os4, os5, os6, os7, H3, ColPix, os8;
		float  RM, H2;
		unsigned char Rx, Gx, Bx;

		os1 = UBRM1 + 1;
		os3 = UBHM1 + 1;
		os4 = UBID1 + 1;
		os5 = UBID2 + 1;
		if (XAD>1)
			Band = (int)(XAD / 1 + 1 + 1);
		else
			Band = 0;
		//Band=5;
		
		
//#pragma omp parallel for private( YP, os2, os6, X, XP, RM, H2, H3, ColPix, R, G, B, Rx, Gx, Bx, M, os7, N, os8)

		for (Y = SY; Y <= PosE1; Y += StS) {
			
			YP = (int)((Y - PosS1) * XAD);
			//if (YP != (int)(((Y+1) - PosS1) * XAD)) {
				if (YP <= UBID3 && YP >= 0) {
					os2 = os1*Y;
					os6 = YP*os4*os5;
					for (X = SX; X <= PosE0; X += StSX) {
						XP = (int)((X - PosS0) * XAD);
						//if (XP != (int)((X - PosS0 + StSX) * XAD)) {
							if (XP <= UBID2 && XP >= 0) {
								//'Z = GetColPix(CurScale, MR, Y, Pict, PosS(0), PosS(1), XAD, PosE(0), StSX, X, UBound(RegionMat, 1), UBound(HeatMap, 1), Min, RegionMat(0, 0), ColPix, HeatMap(0, 0))
								RM = RegionMat[X + os2];
								if (RM >= Min) {
									H2 = (RM - Min) / MR;

									if (H2 <= 1)
										H3 = (int)(H2 * 1020);
									else
										H3 = 1020;

									ColPix = HeatMap[CurScale + H3*os3];
									R = (ColPix / 65536);
									G = ((ColPix - R * 65536) / 256);
									B = (ColPix - R * 65536 - G * 256);
									Rx = (unsigned char)(R);
									Gx = (unsigned char)(G);
									Bx = (unsigned char)(B);


								}
								else {
									Rx = bkr;
									Gx = bkg;
									Bx = bkb;
								}
								if (Band == 0) {
									os7 = XP*os4 + os6;
									ImageData[os7] = Rx;
									ImageData[1 + os7] = Gx;
									ImageData[2 + os7] = Bx;
								}
								else {
									for (M = YP - Band; M <= YP + Band; M++) {
										if (M <= UBID3 && M >= 0) {
											os8 = M*os4*os5;

											for (N = XP - Band; N <= XP + Band; N++) {
												if (N <= UBID2 && N >= 0) {
													os7 = N*os4 + os8;
													ImageData[os7] = Rx;
													ImageData[1 + os7] = Gx;
													ImageData[2 + os7] = Bx;
												}
											}
										}

									}


								}


							}
						//}
					}
				}
			//}
		}



		return(1);
	}


	int MyMathFuncs::ConvSimToDistP(int SLen, int Nextno, int UBDistance, int UBPermvalid, int UBFubvalid, int UBSubvalid, short int *RedoDist, float *Distance, float *FMat, float *SMat, float *PermValid, float *PermDiffs, float *Fubvalid, float *Fubdiffs, float *SubValid, float *SubDiffs)
	{

		double DistXY, ValidXY, TempVal1;
		int X, Y, osUBP, osUBF, osUBS, osUBD, off1;

		osUBP = (UBPermvalid + 1);

		osUBS = (UBSubvalid + 1);
		osUBF = (UBFubvalid + 1);
		osUBD = (UBDistance + 1);
		//this changes similarities to JC distances - it needs to be done in c

		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel for private(X, off1, Y, DistXY, ValidXY, TempVal1)
		for (X = 0; X <= Nextno; X++) {

			off1 = osUBP*X;
			for (Y = X + 1; Y <= Nextno; Y++) {
				//do distance
				if (RedoDist[X] + RedoDist[Y] > 0.000001) {
					DistXY = (double)(SubValid[Y + osUBS*X] + Fubvalid[Y + osUBF*X]);
					PermValid[Y + off1] = (float)(DistXY);
					PermValid[X + osUBP*Y] = (float)(DistXY);

					DistXY = (double)(SubDiffs[Y + osUBS*X] + Fubdiffs[Y + osUBF*X]);
					PermDiffs[Y + off1] = (float)(DistXY);
					PermDiffs[X + osUBP*Y] = (float)(DistXY);


					if (PermDiffs[Y + osUBP*X] > 0) {
						DistXY = (double)(PermDiffs[Y + osUBS*X] / PermDiffs[Y + osUBS*X]);
						Distance[Y + osUBD*X] = (float)(DistXY);
						Distance[X + osUBD*Y] = (float)(DistXY);
					}
					else {
						Distance[Y + osUBD*X] = (float)(0.0);
						Distance[X + osUBD*Y] = (float)(0.0);
					}

				}
				else {


					DistXY = (double)(PermValid[Y + osUBP*X] - SubValid[Y + osUBS*X]);
					Fubvalid[Y + osUBF*X] = (float)(DistXY);
					Fubvalid[X + osUBF*Y] = (float)(DistXY);

					DistXY = (double)(PermDiffs[Y + osUBP*X] - SubDiffs[Y + osUBS*X]);
					Fubdiffs[Y + osUBF*X] = (float)(DistXY);
					Fubdiffs[X + osUBF*Y] = (float)(DistXY);

				}
				//Do FMat
				ValidXY = (double)(Fubvalid[Y + osUBF*X]); //PermValidx(X, Y) - SubValid(X, Y)
				if (ValidXY > 0.0000001) {
					DistXY = (double)((ValidXY - Fubdiffs[Y + osUBF*X]) / ValidXY);
					if (DistXY > 0.25) {
						TempVal1 = ((4 * DistXY - 1) / 3);
						DistXY = (log(TempVal1));
						DistXY = (-0.75 * DistXY);
					}
					else
						DistXY = 10.0;
				}

				else
					DistXY = 10.0;

				FMat[Y + osUBF*X] = (float)(DistXY);
				FMat[X + osUBF*Y] = (float)(DistXY);

				//Do SMat
				ValidXY = (double)(SubValid[Y + osUBS*X]);
				if (ValidXY > 0.0000001) {
					DistXY = double((ValidXY - SubDiffs[Y + osUBS*X]) / ValidXY);
					if (DistXY > 0.25) {
						TempVal1 = ((4 * DistXY - 1) / 3);
						DistXY = (log(TempVal1));
						DistXY = (-0.75 * DistXY);
					}
					else
						DistXY = 10.0;

				}
				else
					DistXY = 10.0;

				SMat[Y + osUBS*X] = (float)(DistXY);
				SMat[X + osUBS*Y] = (float)(DistXY);

			}
		}



		for (X = 0; X <= Nextno; X++) {
			Distance[X + osUBD*X] = 1;
			FMat[X + osUBF*X] = 0;
			PermValid[X + osUBD*X] = (float)(SLen);
			if (X <= UBPermvalid)
				PermDiffs[X + osUBP*X] = 0;

			SubValid[X + osUBS*X] = (float)(SLen / 2);
			SubDiffs[X + osUBS*X] = 0;
		}
		omp_set_num_threads(2);
		return(1);
	}


	int MyMathFuncs::EraseEmptiesP(int Nextno, int UB, int UBFM, int SCO, int *ISeqs, float *FMat, float *FubValid, float *SMat, float *SubValid) {
		int X, Y, EraseF, off1, off2;

		off1 = UBFM + 1;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private(X, off2, Y, EraseF)
		for (X = 0; X <= Nextno; X++) {
			off2 = X*off1;
			if (FMat[X + off2] != 3) {

				EraseF = 0;

				for (Y = 0; Y <= 2; Y++) {
					if (ISeqs[Y] <= Nextno) {
						if (X != ISeqs[0] && X != ISeqs[1] && X != ISeqs[2]) {

							if (FubValid[ISeqs[Y] + off2] < SCO || SubValid[ISeqs[Y] + off2] < SCO) {
								EraseF = 1;
								break;
							}

						}
					}
				}

				if (EraseF == 1) {

					for (Y = 0; Y <= Nextno; Y++) {
						FMat[X + Y*off1] = 3.0;
						FMat[Y + off2] = 3.0;
						if (UB > 0) {
							SMat[X + Y*off1] = 3.0;
							SMat[Y + off2] = 3.0;
						}


					}

				}

			}
		}
		omp_set_num_threads(2);
		return(1);
		
	}

	int MyMathFuncs::MakeBigMap(int IStart, int XRes, int YRes, int MBN, int TType, int TNum, int TSH, float XSize, float TSingle, int UBMB1, int UBMB2, int UBMB3,int UBMB4, HDC pict, float *MapBlocks) {
		int X, off1, off2, off3, off4, off5, off6, off7, os2, os3, hold, os4, t1,t2;
		float AH1;
		//CPen newpen;

		off1 = UBMB1 + 1;
		off2 = UBMB2 + 1;
		off3 = UBMB3 + 1;
		off4 = off1*off2;
		off5 = off4*off3;
		off6 = TType + TNum*off1;
		os2 = 2 * off4;
		os3 = 3 * off4;
		t1 = (-IStart / TSingle) -13;
		t2 = (TSH-IStart) / TSingle;
		HGDIOBJ original = NULL;
		HGDIOBJ original2 = NULL;

		original = SelectObject(pict, GetStockObject(DC_PEN));
		original2 = SelectObject(pict, GetStockObject(DC_BRUSH));
		SelectObject(pict, GetStockObject(DC_PEN));
		SelectObject(pict, GetStockObject(DC_BRUSH));
//#pragma omp parallel for private(X, off7, AH1, hold)
		for (X = 1; X <= MBN; X++) {
			if (X <= UBMB4) {
				off7 = X*off5;
				os4 = off6 + off7;
				AH1 = MapBlocks[os2 + os4];

				if  (AH1 >= t1) {

					

					if (AH1 <= t2) {
						hold = MapBlocks[os3 + os4];
//#pragma omp critical
	//					{
						SetDCBrushColor(pict, hold);
						SetDCPenColor(pict, hold);

						Rectangle(pict, (5 + MapBlocks[os4] * XSize) * XRes, (IStart + (AH1 + 4) * TSingle) * YRes, (5 + MapBlocks[off4 + os4] * XSize) * XRes, (IStart + (AH1 + 13) * TSingle) * YRes);
						//MoveToEx(pict, (5 + MapBlocks[off6 + off7] * XSize) * XRes, (IStart + (AH1 + 4) * TSingle) * YRes, 0);

						//LineTo(pict, (5 + MapBlocks[off6 + off4 + off7] * XSize) * XRes, (IStart + (AH1 + 13) * TSingle) * YRes);
						//Form2.Picture5.Line(((5 + MapBlocks[off6 + off7] * XSize) * XRes), (IStart + (AH1 + 4) * TSingle) * YRes) - ((5 + MapBlocks[off6 + off4 + off7] * XSize) * XRes, (IStart + (AH1 + 13) * TSingle) * YRes), MapBlocks[off6 + os3 + off7], BF
//						}
					}
				}
			}

		}
		
		SelectObject(pict, original);
		SelectObject(pict, original2);

		return(1);
	}


	int MyMathFuncs::SeqColBlocks(HDC Pict, int UBSL, double tTYF, float TCA, double XConA, int NumSeqLines, int Targ, int VSV, int *SeqLines) {
		int x, hold, X1, X2, Y1, Y2, os, os2;
		//float TCA;
		HGDIOBJ original = NULL;
		HGDIOBJ original2 = NULL;
		os = UBSL + 1;
		original = SelectObject(Pict, GetStockObject(DC_PEN));
		original2 = SelectObject(Pict, GetStockObject(DC_BRUSH));

		SelectObject(Pict, GetStockObject(DC_PEN));
		SelectObject(Pict, GetStockObject(DC_BRUSH));

		//TCA = (float)(tTYF * XConA);



//#pragma omp parallel for private(os2, Y1, Y2, hold, X1, X2)
		for (x = 0; x < NumSeqLines; x++) {
			os2 = x*os;
			if (SeqLines[os2] == 0) { //rectangle

				Y1 = SeqLines[3 + os2] * tTYF - VSV;
				Y2 = SeqLines[5 + os2] * tTYF - VSV;
				if ((Y1 >= 0 && Y2 <= Targ) || (Y1 <= 0 && Y2 >= Targ) || (Y1 >= 0 && Y1 <= Targ) || (Y2 >= 0 && Y2 <= Targ)) {

					hold = SeqLines[1 + os2];
					SetDCBrushColor(Pict, hold);
					SetDCPenColor(Pict, hold);


					X1 = SeqLines[2 + os2] * TCA;
					X2 = SeqLines[4 + os2] * TCA;


					Rectangle (Pict, X1, Y1, X2, Y2);

				}
			}
			else
				break;

		}

		SelectObject(Pict, original);
		SelectObject(Pict, original2);

		return(x);

	}

	int MyMathFuncs::SeqColBlocksP(HDC Pict, int UBSL, double tTYF, float TCA, double XConA, int NumSeqLines, int Targ, int VSV, int *SeqLines) {
		int x, hold, X1, X2, Y1, Y2, os, os2;
		//float TCA;
		HGDIOBJ original = NULL;
		HGDIOBJ original2 = NULL;
		os = UBSL + 1;
		original = SelectObject(Pict, GetStockObject(DC_PEN));
		original2 = SelectObject(Pict, GetStockObject(DC_BRUSH));

		SelectObject(Pict, GetStockObject(DC_PEN));
		SelectObject(Pict, GetStockObject(DC_BRUSH));

		//TCA = (float)(tTYF * XConA);
		x = 0;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

// The source GUI drawing loop exits on its first non-rectangle row.  GCC's
// OpenMP front end correctly rejects `break` in a parallel loop, so keep this
// presentation-only routine serial when the native method passes use OpenMP.
#if !defined(NEXT_RDP_USE_REAL_OPENMP)
#pragma omp parallel for private(os2, Y1, Y2, hold, X1, X2)
#endif
		for (x = 0; x < NumSeqLines; x++) {
			os2 = x*os;
			if (SeqLines[os2] == 0) { //rectangle

				Y1 = SeqLines[3 + os2] * tTYF - VSV;
				Y2 = SeqLines[5 + os2] * tTYF - VSV;
				if ((Y1 >= 0 && Y2 <= Targ) || (Y1 <= 0 && Y2 >= Targ) || (Y1 >= 0 && Y1 <= Targ) || (Y2 >= 0 && Y2 <= Targ)) {

					hold = SeqLines[1 + os2];
					SetDCBrushColor(Pict, hold);
					SetDCPenColor(Pict, hold);


					X1 = SeqLines[2 + os2] * TCA;
					X2 = SeqLines[4 + os2] * TCA;


					Rectangle(Pict, X1, Y1, X2, Y2);

				}
			}
			else
				break;

		}

		SelectObject(Pict, original);
		SelectObject(Pict, original2);
		omp_set_num_threads(2);
		return(x);
		
	}


	


	int MyMathFuncs::ExtraRemovalsP(int Nextno, int UBF, int UBS, int *ISeqs, int *ExtraRemove, float *FMat, float *SMat) {

		int MaxRemove, Z, X, Y, osf, oss, WinX;

		osf = UBF + 1;
		oss = UBS + 1;

		for (Z = 0; Z <= Nextno; Z++) {
			if (FMat[Z + Z*osf] < 3) {
				for (X = 0; X <= Nextno; X++)
					ExtraRemove[X] = 0;


				for (X = 0; X <= Nextno; X++) {

					if (FMat[X + X*osf] < 3) {
						for (Y = X + 1; Y <= Nextno; Y++) {
							if (FMat[Y + Y*osf] < 3) {
								if (FMat[X + Y*osf] >= 3) {
									ExtraRemove[X] = ExtraRemove[X] + 1;
									ExtraRemove[Y] = ExtraRemove[Y] + 1;
								}

								if (SMat[X + Y*oss] >= 3) {
									ExtraRemove[X] = ExtraRemove[X] + 1;
									ExtraRemove[Y] = ExtraRemove[Y] + 1;
								}
							}
						}
					}

				}

				MaxRemove = 0;

				WinX = -1;
				for (X = 0; X <= Nextno; X++) {
					if (FMat[X + X*osf] < 3) {
						if (ExtraRemove[X] > MaxRemove) {
							if (X != ISeqs[0] && X != ISeqs[1] && X != ISeqs[2]) {

								MaxRemove = ExtraRemove[X];
								WinX = X;

							}

						}
					}
				}
				if (WinX == -1)
					break;

				if (MaxRemove > 0) {
					FMat[WinX + WinX*osf] = float(3.0);
					SMat[WinX + WinX*oss] = float(3.0);
				}
			}

		}

		return(1);
	}

	int MyMathFuncs::PrintSeqs(int X1, HDC Pict, int UBST, int LOS, int Targ, int UBSL, int NumSeqLines, int StartX, int VSV, int SLFS, int SeqSpaceIncrement, int FirstSeq, int *SeqLines, char *SeqText) {
		int X, os, os2, Y1, len, os3, Z;
		//LPCSTR test = "test";
		//long *enterorder;
		// enterorder = (long *)malloc(numsp*sizeof(long));
		//s = LOS * (UBST + 1);
		/*if (LOS < 2)
			LOS = 2;*/
		wchar_t *pwcs;
		pwcs = (wchar_t *)malloc(LOS*sizeof(wchar_t));
		
		//len = mbstowcs(pwcs, *SeqText, MB_CUR_MAX);
		

		os2 = UBSL + 1;
		for (X = StartX; X < NumSeqLines; X++) {
		
			
			os = X*os2;
			if (SeqLines[os] != 0) {

				
				Y1 = SLFS + (X - FirstSeq) * SeqSpaceIncrement;
				if (Y1 > -30 && Y1 < Targ) {

					//std::wstring stemp = s2ws(SeqText[SeqLines[1 + os]]);
					//LPCWSTR result = stemp.c_str();
					//return(Y1);
					os3 = LOS*SeqLines[1 + os];

					for (Z = 0; Z < LOS; Z++)
						pwcs[Z] = SeqText[os3+Z];
						//pwcs[Z] = SeqText[Z+SeqLines[1 + os]*os3];
					
					//TextOutA(Pict, 0, Y1, SeqText[SeqLines[1 + os]], LOS);
					TextOut(Pict, X1, Y1, pwcs,LOS);
					//return(Y1);
					//TextOut(Pict, 0, Y1, text, LOS);
				}
			}

		}
		free(pwcs);
		return(1);
	}

	int MyMathFuncs::PrintSeqsP(int X1, HDC Pict, int UBST, int LOS, int Targ, int UBSL, int NumSeqLines, int StartX, int VSV, int SLFS, int SeqSpaceIncrement, int FirstSeq, int *SeqLines, char *SeqText) {
		int X, os, os2, Y1, len, os3, Z, rv;
		//LPCSTR test = "test";
		//long *enterorder;
		// enterorder = (long *)malloc(numsp*sizeof(long));
		//s = LOS * (UBST + 1);
		wchar_t *pwcs;
		

		//len = mbstowcs(pwcs, *SeqText, MB_CUR_MAX);


		os2 = UBSL + 1;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private(X, os, Y1, os3, Z, pwcs, rv)
		for (X = StartX; X < NumSeqLines; X++) {


			os = X*os2;
			if (SeqLines[os] != 0) {


				Y1 = SLFS + (X - FirstSeq) * SeqSpaceIncrement;
				if (Y1 > -30 && Y1 < Targ) {

					//std::wstring stemp = s2ws(SeqText[SeqLines[1 + os]]);
					//LPCWSTR result = stemp.c_str();
					//return(Y1);
					os3 = LOS*SeqLines[1 + os];
					
					pwcs = (wchar_t *)malloc(LOS*sizeof(wchar_t));
//#pragma omp parallel for 
					for (Z = 0; Z < LOS; Z++)
						pwcs[Z] = SeqText[os3 + Z];
					//pwcs[Z] = SeqText[Z+SeqLines[1 + os]*os3];

					//TextOutA(Pict, 0, Y1, SeqText[SeqLines[1 + os]], LOS);
					rv = 0;
					while (rv == 0) {
						rv = TextOut(Pict, X1, Y1, pwcs, LOS);
					}
					

					free(pwcs);
					//return(Y1);
					//TextOut(Pict, 0, Y1, text, LOS);
				}
			}

		}
		omp_set_num_threads(2);
		return(1);
	}

	int MyMathFuncs::MakeBigMapB(int IStart, int XRes, int YRes, int MBN, int TType, int TNum, int TSH, float XSize, float TSingle, int UBMB1, int UBMB2, int UBMB3, int UBMB4, HDC pict, float *MapBlocks) {
		int X, off1, off2, off3, off4, off5, off6, off7, off8, os2, os3, hold;
		float AH1;
		//CPen newpen;

		off1 = UBMB1 + 1;
		off2 = UBMB2 + 1;
		off3 = UBMB3 + 1;
		off4 = off1*off2;
		off5 = off4*off3;
		//off6 = TType + TNum*off1;
		off8 = TType*off4 + TNum*off5;
		//os2 = 2 * off4;
		//os3 = 3 * off4;

		HGDIOBJ original = NULL;
		HGDIOBJ original2 = NULL;
		original = SelectObject(pict, GetStockObject(DC_PEN));
		original2 = SelectObject(pict, GetStockObject(DC_BRUSH));
		SelectObject(pict, GetStockObject(DC_PEN));
		SelectObject(pict, GetStockObject(DC_BRUSH));

		for (X = 1; X <= MBN; X++) {
			if (X <= UBMB2) {
				off7 = X*off1 + off8;
				AH1 = MapBlocks[2 + off7];

				if (IStart + (AH1 + 13) * TSingle >= 0) {



					if (IStart + AH1 * TSingle <= TSH) {
						hold = MapBlocks[3 + off7];
						SetDCBrushColor(pict, hold);
						SetDCPenColor(pict, hold);
						Rectangle(pict, (5 + MapBlocks[off7] * XSize) * XRes, (IStart + (AH1 + 4) * TSingle) * YRes, (5 + MapBlocks[1 + off7] * XSize) * XRes, (IStart + (AH1 + 13) * TSingle) * YRes);
						
						//Form2.Picture5.Line (((5 + MapBlocks(0,x,TType, TNum) * XSize) * XRes), (IStart + (AH1 + 4) * TSingle) * YRes)-((5 + MapBlocks(1,x,TType, TNum) * XSize) * XRes, (IStart + (AH1 + 13) * TSingle) * YRes), MapBlocks(3, X,TType, TNum ), BF
						//MoveToEx(pict, (5 + MapBlocks[off6 + off7] * XSize) * XRes, (IStart + (AH1 + 4) * TSingle) * YRes, 0);

						//LineTo(pict, (5 + MapBlocks[off6 + off4 + off7] * XSize) * XRes, (IStart + (AH1 + 13) * TSingle) * YRes);
						//Form2.Picture5.Line(((5 + MapBlocks[off6 + off7] * XSize) * XRes), (IStart + (AH1 + 4) * TSingle) * YRes) - ((5 + MapBlocks[off6 + off4 + off7] * XSize) * XRes, (IStart + (AH1 + 13) * TSingle) * YRes), MapBlocks[off6 + os3 + off7], BF

					}
				}
			}

		}

		SelectObject(pict, original);
		SelectObject(pict, original2);

		return(1);
	}


	float MyMathFuncs::GetMaxXPos(int CharLen, int TNum, int TType, int TDL0, int UBON, int UBTD1, int UBTD2, int UBTD3, int UBTD4, float PRat, short int *ONameLen, float *TreeDrawB) {
		int X, off1, off2, off3, off4, off5, AH1;
		float MaxXPos, CXP;
		MaxXPos = 0.0;
		off1 = UBTD1 + 1;
		off2 = (UBTD2 + 1)*off1;
		off3 = (UBTD3 + 1)*off2;
		off4 = (UBTD4 + 1)*off3;
		off5 = TNum*off2 + TType*off3;
		for (X = 0; X <= TDL0; X++) {

			CXP = TreeDrawB[X*off1 + off5] / PRat;

			AH1 = TreeDrawB[2 + X*off1 + off5];
			if (AH1 > -1) {
				if (AH1 <= UBON)
					CXP = CXP + CharLen * (ONameLen[AH1] + 2);

			}

			if (CXP > MaxXPos)
				MaxXPos = CXP;
			
		}
		return(MaxXPos);
	}


	int MyMathFuncs::MakeTreeDrawB2(int UBA, int UBB, int UBC, int UBD, int UBE, float *TreeDraw, float *TreeDrawB)
	{
		int A, B, C, D, E, off1, off2, off3, off4, off5, off6, off7, o1, o2, o3, o4, o5, o6, o7, o8;

		off1 = UBA + 1;
		off2 = (UBB + 1)*off1;
		off3 = (UBC + 1)*off2;
		off4 = (UBD + 1)*off3;

		o1 = UBD + 1;
		o2 = (UBE + 1)*o1;
		o3 = (UBA + 1)*o2;
		o4 = (UBB + 1)*o3;


		for (A = 0; A <= UBA; A++) {
			o5 = A*o2;
			for (B = 0; B <= UBB; B++) {
				o6 = B*o3 + o5;
				off5 = B*off1 + A;
				for (C = 0; C <= UBC; C++) {
					o7 = C*o4 + o6;
					off6 = C*off2 + off5;
					for (D = 0; D <= UBD; D++) {
						o8 = D + o7;
						off7 = D*off3 + off6;
						for (E = 0; E <= UBE; E++) {
							TreeDrawB[E*o1 + o8] = TreeDraw[E*off4 + off7];
							//TreeDrawB[D + E*o1 + A*o2 + B*o3 + C*o4] = TreeDraw[E*off4 + off7];
						}
					}
				}
			}
		}
		return(1);
	}

	float MyMathFuncs::GetMaxXPosB(int UBTTS1, int UBTTS2, int UBTT, int CharLen, int TNum, int TType, int TDL0, int UBON, int UBTD1, int UBTD2, int UBTD3, int UBTD4, float PRat, short int *ONameLen, float *TreeDrawB, int *TreeTraceSeqs, int *TreeTrace) {
		int X, off1, off2, off3, off4, off5, AH1, AH2, o1;
		float MaxXPos, CXP;
		MaxXPos = 0.0;
		off1 = UBTD1 + 1;
		off2 = (UBTD2 + 1)*off1;
		off3 = (UBTD3 + 1)*off2;
		off4 = (UBTD4 + 1)*off3;
		off5 = TNum*off2 + TType*off3;
		o1 = UBTTS1 + 1;
		for (X = 0; X <= TDL0; X++) {

			CXP = TreeDrawB[X*off1 + off5] / PRat;

			AH1 = TreeDrawB[2 + X*off1 + off5];
			if (AH1 > -1) {
				if (AH1 <= UBON) {

					if (UBTTS1 > 0) {
						if (UBTTS2 >= X) {
							AH2 = TreeTraceSeqs[1 + AH1*o1];
							if (UBTT >= AH2) {
								if (TreeTrace[AH2] <= UBON)
									CXP = CXP + CharLen * (ONameLen[TreeTrace[AH2]] + 2);
							}
						}
					}
				}

			}

			if (CXP > MaxXPos)
				MaxXPos = CXP;

		}
		return(MaxXPos);
	}

	int MyMathFuncs::DrawTreeLines(HDC Pict, int IStart, int TSHx, int TargetA, int TNum, int TType, int TDL1, int UBTD1, int UBTD2, int UBTD3, int UBTD4, float PRat, float TSingle, int *OS, float *TreeDrawB) {
		int X, GoOn, AH1, AH2, off1, off2, off3, off4, off5, os, off6, hold;
		HGDIOBJ original = NULL;
		HGDIOBJ original2 = NULL;

		off1 = UBTD1 + 1;
		off2 = (UBTD2 + 1)*off1;
		off3 = (UBTD3 + 1)*off2;
		off4 = (UBTD4 + 1)*off3;
		off5 = TNum*off2 + TType*off3 + off4;

		original = SelectObject(Pict, GetStockObject(DC_PEN));
		original2 = SelectObject(Pict, GetStockObject(DC_BRUSH));
		SelectObject(Pict, GetStockObject(DC_PEN));
		SelectObject(Pict, GetStockObject(DC_BRUSH));

		GoOn = 0;
		os = 0;
		AH1 = -1;
		for (X = 0; X <= TDL1; X++) {
			off6 = X*off1 + off5;
			AH1 = (int)(TreeDrawB[1 + off6]);
				
			AH2 = (int)(TreeDrawB[3 +off6]);

			if (AH1 >= TargetA || AH2 >= TargetA) {


				if (AH1 <= TSHx || AH2 <= TSHx) {
					hold = (int)(TreeDrawB[4 + off6]);
					SetDCPenColor(Pict, hold);
					
					
					MoveToEx(Pict, TreeDrawB[off6] / PRat, IStart + AH1 * TSingle, 0);
						
					LineTo(Pict, TreeDrawB[2 + off6] / PRat, IStart + AH2 * TSingle);
						
				}
				else{

					if (AH1 > os) {
						os = AH1;
						GoOn = 1;
					}
				}
			}
		}
		*OS = os;

		SelectObject(Pict, original);
		SelectObject(Pict, original2);

		return(GoOn);
	}


	int MyMathFuncs::DoAABlocksP(int xRes, int yRes, int UBPC, int UBIX, int UBID22, int UBID23, int UBID12, int UBID13, int *PCount, int *ImageX, unsigned char *ImageData, unsigned char *ImageData2) {
		int X, Y, Z, A, B, off1, off2, off3, off4, off6, off5, off7, off8, off9, off10, off11, off12;
		float XRes, YRes;
		XRes = (float)(xRes);
		YRes = (float)(yRes);
		off1 = UBPC + 1;
		off2 = 3 * (UBIX + 1);
		off3 = 4 * (UBID22 + 1);
		off4 = 4 * (UBID12 + 1);
//#pragma omp parallel for private(X, A, off5, off6, Y, B, off9, off7, off8, Z)
		
		for (Y = 0; Y <= UBID23; Y++) {
			B = (int)(Y / YRes);
			off9 = B*off1;
			for (X = 0; X <= UBID22; X++) {

				//x = (float)(X);
				A = (int)(X / XRes);
				off5 = A * 3;
				off6 = X * 4;
				PCount[A + off9] += 1;
				off7 = B*off2 + off5;
				off8 = off6 + Y*off3;
				for (Z = 0; Z <= 2; Z++)

//				{
//#pragma omp atomic
					ImageX[Z + off7] += (int)(ImageData2[Z + off8]);
//				}
			}
		}

//#pragma omp parallel for private(Y, off10, off11, off12, X, off6, off8, off5, off7, off9, Z)
		for (Y = 0; Y <= UBID13; Y++) {
			off10 = Y*off1;
			off11 = Y*off2;
			off12 = Y*off3;
			for (X = 0; X <= UBID12; X++) {
				off6 = X * 3;
				off8 = X * 4;
			
				off5 = X + off10;
				off7 = off6 + off11;
				off9 = off8 + off12;
//#pragma omp critical
//				{
					for (Z = 0; Z <= 2; Z++)
						ImageData[Z + off9] = (unsigned char)(ImageX[Z + off7] / PCount[off5]);

//				}

			}
		}
		return(1);
	}
	int MyMathFuncs::MakeSeqCatCount2P(int Nextno, int LSeq, int UBSN1, int UBSCC1, int StartPosInAlign, int EndPosInAlign, int *SeqCatCount, int *AA, short int *SeqNum, unsigned char *NucMat, unsigned char *SeqSpace, unsigned char  *NucMatB, unsigned char  *NucMat2, unsigned char *flp, unsigned char *ml, unsigned char *nl) {

		int X, Y, os1, os2, ANum, NumNucs, off3, off2, aapos;
		os1 = UBSN1 + 1;
		os2 = UBSCC1 + 1;
		ANum = 0;
		SeqCatCount[0] = 0;
		SeqCatCount[1] = 0;
		SeqCatCount[2] = 0;
		SeqCatCount[3] = 0;
		SeqCatCount[4] = 0;
		SeqCatCount[5] = 0;
		SeqCatCount[6] = 0;
		SeqCatCount[7] = 0;
		SeqCatCount[8] = 0;
		SeqCatCount[9] = 0;
		off3 = Nextno + 1;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
		if (StartPosInAlign <= EndPosInAlign) {


			for (X = StartPosInAlign; X <= EndPosInAlign; X++) {
				off2 = off3*X;
				AA[0] = 0;
				AA[1] = 0;
				AA[2] = 0;
				AA[3] = 0;
				AA[4] = 0;
				
#pragma omp parallel for private(Y, aapos)
				for (Y = 0; Y <= Nextno; Y++) {
					aapos = NucMat[SeqNum[X + Y*os1]];
					if (AA[aapos] == 0)
						AA[aapos] = 1;
				}
				ANum = (int)(AA[1] + AA[2] + AA[3] + AA[4]);

				SeqSpace[X] = AA[0];
				NumNucs = (int)(AA[0]);
				SeqCatCount[NumNucs + ANum*os2] = SeqCatCount[NumNucs + ANum*os2] + 1;

				for (Y = 1; Y <= 4; Y++) {
					if (AA[Y] == 1) {
						NucMat2[NucMatB[Y]] = NumNucs;
						NumNucs = NumNucs + 1;
					}
				}
				ml[X] = AA[0];
				nl[X] = ANum;

#pragma omp parallel for private(Y)
				for (Y = 0; Y <= Nextno; Y++)
					flp[Y + off2] = NucMat2[(int)(SeqNum[X + Y*os1])];
			}
		}
		else if (StartPosInAlign > EndPosInAlign) {
			for (X = StartPosInAlign; X <= LSeq; X++) {
				off2 = off3*X;
				AA[0] = 0;
				AA[1] = 0;
				AA[2] = 0;
				AA[3] = 0;
				AA[4] = 0;
#pragma omp parallel for private(Y, aapos)
				for (Y = 0; Y <= Nextno; Y++) {
					aapos = NucMat[SeqNum[X + Y*os1]];
					if (AA[aapos] == 0)
						AA[aapos] = 1;
				}

				ANum = (int)(AA[1] + AA[2] + AA[3] + AA[4]);

				SeqSpace[X] = AA[0];
				NumNucs = (int)(AA[0]);
				SeqCatCount[NumNucs + ANum*os2] = SeqCatCount[NumNucs + ANum*os2] + 1;
				NumNucs = AA[0];
				for (Y = 1; Y <= 4; Y++) {
					if (AA[Y] == 1) {
						NucMat2[NucMatB[Y]] = NumNucs;
						NumNucs = NumNucs + 1;
					}
				}
				ml[X] = AA[0];
				nl[X] = ANum;
#pragma omp parallel for private(Y)
				for (Y = 0; Y <= Nextno; Y++)
					flp[Y + off2] = NucMat2[(int)(SeqNum[X + Y*os1])];
			}
			for (X = 1; X <= EndPosInAlign; X++) {
				off2 = off3*X;
				AA[0] = 0;
				AA[1] = 0;
				AA[2] = 0;
				AA[3] = 0;
				AA[4] = 0;
#pragma omp parallel for private(Y, aapos)
				for (Y = 0; Y <= Nextno; Y++) {
					aapos = NucMat[SeqNum[X + Y*os1]];
					if (AA[aapos] == 0)
						AA[aapos] = 1;
				}

				ANum = (int)(AA[1] + AA[2] + AA[3] + AA[4]);

				SeqSpace[X] = AA[0];
				NumNucs = (int)(AA[0]);
				SeqCatCount[NumNucs + ANum*os2] = SeqCatCount[NumNucs + ANum*os2] + 1;
				NumNucs = AA[0];
				for (Y = 1; Y <= 4; Y++) {
					if (AA[Y] == 1) {
						NucMat2[NucMatB[Y]] = NumNucs;
						NumNucs = NumNucs + 1;
					}
				}
				ml[X] = AA[0];
				nl[X] = ANum;
#pragma omp parallel for private(Y)
				for (Y = 0; Y <= Nextno; Y++)
					flp[Y + off2] = NucMat2[(int)(SeqNum[X + Y*os1])];
			}

		}
		omp_set_num_threads(2);
		return(1);
	}


	int MyMathFuncs::FindSubSeqP(int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int seq3, short int spacerno, short int *seqnum, short int *xoverseqnum, char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer)

		

	{

		int hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4, xow, wmx, wpx, target;
		int s1, s2, s3, sz, os1, os2, se1, se2, se3, ah0, ah1, ah2;//,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int g = 0;

		se1 = seq1*lenseq;
		se2 = seq2*lenseq;
		se3 = seq3*lenseq;
		oc = 0;
		hoc = 0;
		xow = xoverwindow;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		os1 = xoverwindow + lenseq + xoverwindow2;
		os2 = xoverwindow + lenseq2 + xoverwindow4;
		ah0 = 0;
		ah1 = 0;
		ah2 = 0;



		if (spacerflag == 0) {
			//so1 = seq1*lenseq;
			//so2 = seq2*lenseq;
			//so3 = seq3*lenseq;
			target = (lenseq + xoverwindow2) * 3;
			for (x = 1; x < target; x++)
				xoverseqnumw[x] = 0;

			hoc = xoverwindow * 10;

			for (x = 1; x < lenseq; x++) {



				s1 = seqnum[x + se1];

				if (s1 != 46) {
					s2 = seqnum[x + se2];
					if (s2 != 46) {
						s3 = seqnum[x + se3];
						if (s3 != 46) {
							if (s1 != s2) {
								//if (s1 == s2 || s1 == s3 || s2 == s3){


								if (s1 == s3) {
									y++;
									//xoverseqnumw[y + xoverwindow] = 0;
									xoverseqnumw[y + os1] = 1;
									//xoverseqnumw[y + os2] = 0;
									ah1++;
									xdiffpos[y] = x;
									//xposdiff[x] = y;
								}
								else if (s2 == s3) {
									y++;
									//xoverseqnumw[y + xoverwindow] = 0;
									//xoverseqnumw[y + os1] = 0;
									xoverseqnumw[y + os2] = 1;
									ah2++;
									xdiffpos[y] = x;
									//xposdiff[x] = y;
								}


								//}
							}
							else if (s1 != s3) {
								//if (s1 == s2 || s1 == s3 || s2 == s3){


								if (s1 == s2) {
									y++;
									xoverseqnumw[y + xoverwindow] = 1;
									//xoverseqnumw[y + os1] = 0;
									//xoverseqnumw[y + os2] = 0;
									ah0++;
									xdiffpos[y] = x;
									//xposdiff[x] = y;
								}
								else if (s2 == s3) {
									y++;
									//xoverseqnumw[y + xoverwindow] = 0;
									//xoverseqnumw[y + os1] = 0;
									xoverseqnumw[y + os2] = 1;
									ah2++;
									xdiffpos[y] = x;
									//xposdiff[x] = y;
								}

								//}
							}
						}
					}
				}

				xposdiff[x] = y;

			}
			//return(ah1);	
		}

		else if (spacerflag == 1) {
			for (x = 1; x < lenseq; x++) {

				*(xposdiff + x) = y;
				s1 = *(seqnum + x + se1);
				s2 = *(seqnum + x + se2);
				s3 = *(seqnum + x + se3);

				if (s1 != s2 || s1 != s3) {

					if (s1 == s2 || s1 == s3 || s2 == s3) {

						if (s1 != 46) {
							if (s2 != 46) {
								if (s3 != 46) {
									if (s1 != s2 && s1 != s3) {
										//If seq1 is odd one


										if (outlyer == seq1) {
											if (oc > 0) {
												if (oc > hoc)
													hoc = oc;
												oc--;
											}
											if (s1 == s2) {
												y++;
												xoverseqnumw[y + xoverwindow] = 1;
												xoverseqnumw[y + os1] = 0;
												xoverseqnumw[y + os2] = 0;
												ah0++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}
											else if (s1 == s3) {
												y++;
												xoverseqnumw[y + xoverwindow] = 0;
												xoverseqnumw[y + os1] = 1;
												xoverseqnumw[y + os2] = 0;
												ah1++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}
											else if (s2 == s3) {
												y++;
												xoverseqnumw[y + xoverwindow] = 0;
												xoverseqnumw[y + os1] = 0;
												xoverseqnumw[y + os2] = 1;
												ah2++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}

										}
										else {
											for (z = 1; z <= spacerno; z++) {
												sz = *(seqnum + x + *(spacerseqs + z) * lenseq);
												if (sz == s1) {
													//If difference is legitimate
													if (*(xdiffpos + y) != x) {
														oc += 2;
														if (s1 == s2) {
															y++;
															xoverseqnumw[y + xoverwindow] = 1;
															xoverseqnumw[y + os1] = 0;
															xoverseqnumw[y + os2] = 0;
															ah0++;
															xdiffpos[y] = x;
															xposdiff[x] = y;
														}
														else if (s1 == s3) {
															y++;
															xoverseqnumw[y + xoverwindow] = 0;
															xoverseqnumw[y + os1] = 1;
															xoverseqnumw[y + os2] = 0;
															ah1++;
															xdiffpos[y] = x;
															xposdiff[x] = y;
														}
														else if (s2 == s3) {
															y++;
															xoverseqnumw[y + xoverwindow] = 0;
															xoverseqnumw[y + os1] = 0;
															xoverseqnumw[y + os2] = 1;
															ah2++;
															xdiffpos[y] = x;
															xposdiff[x] = y;
														}
														z = spacerno;

													}

												}

											}

										}
									}
									else if (s2 != s1  && s2 != s3) {
										//If seq2 is odd one

										if (outlyer == seq2) {
											if (oc >0) {
												if (oc > hoc)
													hoc = oc;
												oc--;
											}
											if (s1 == s2) {
												y++;
												xoverseqnumw[y + xoverwindow] = 1;
												xoverseqnumw[y + os1] = 0;
												xoverseqnumw[y + os2] = 0;
												ah0++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}
											else if (s1 == s3) {
												y++;
												xoverseqnumw[y + xoverwindow] = 0;
												xoverseqnumw[y + os1] = 1;
												xoverseqnumw[y + os2] = 0;
												ah1++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}
											else if (s2 == s3) {
												y++;
												xoverseqnumw[y + xoverwindow] = 0;
												xoverseqnumw[y + os1] = 0;
												xoverseqnumw[y + os2] = 1;
												ah2++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}

										}
										else {
											for (z = 1; z <= spacerno; z++) {
												sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
												if (s2 == sz) {
													//If difference is legitimate

													if (*(xdiffpos + y) != x) {
														oc += 2;
														if (s1 == s2) {
															y++;
															xoverseqnumw[y + xoverwindow] = 1;
															xoverseqnumw[y + os1] = 0;
															xoverseqnumw[y + os2] = 0;
															ah0++;
															xdiffpos[y] = x;
															xposdiff[x] = y;
														}
														else if (s1 == s3) {
															y++;
															xoverseqnumw[y + xoverwindow] = 0;
															xoverseqnumw[y + os1] = 1;
															xoverseqnumw[y + os2] = 0;
															ah1++;
															xdiffpos[y] = x;
															xposdiff[x] = y;
														}
														else if (s2 == s3) {
															y++;
															xoverseqnumw[y + xoverwindow] = 0;
															xoverseqnumw[y + os1] = 0;
															xoverseqnumw[y + os2] = 1;
															ah2++;
															xdiffpos[y] = x;
															xposdiff[x] = y;
														}
														z = spacerno;

													}
												}
											}
										}
									}
									else if (s3 != s1  && s3 != s2) {
										//If seq3 is odd one
										if (outlyer == seq3) {
											if (oc >0) {
												if (oc > hoc)
													hoc = oc;
												oc--;
											}
											if (s1 == s2) {
												y++;
												xoverseqnumw[y + xoverwindow] = 1;
												xoverseqnumw[y + os1] = 0;
												xoverseqnumw[y + os2] = 0;
												ah0++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}
											else if (s1 == s3) {
												y++;
												xoverseqnumw[y + xoverwindow] = 0;
												xoverseqnumw[y + os1] = 1;
												xoverseqnumw[y + os2] = 0;
												ah1++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}
											else if (s2 == s3) {
												y++;
												xoverseqnumw[y + xoverwindow] = 0;
												xoverseqnumw[y + os1] = 0;
												xoverseqnumw[y + os2] = 1;
												ah2++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}


										}
										else {
											for (z = 1; z <= spacerno; z++) {
												sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
												if (s3 == sz) {

													//If difference is legitimate
													if (*(xdiffpos + y) != x) {
														oc += 2;
														if (s1 == s2) {
															y++;
															xoverseqnumw[y + xoverwindow] = 1;
															xoverseqnumw[y + os1] = 0;
															xoverseqnumw[y + os2] = 0;
															ah0++;
															xdiffpos[y] = x;
															xposdiff[x] = y;
														}
														else if (s1 == s3) {
															y++;
															xoverseqnumw[y + xoverwindow] = 0;
															xoverseqnumw[y + os1] = 1;
															xoverseqnumw[y + os2] = 0;
															ah1++;
															xdiffpos[y] = x;
															xposdiff[x] = y;
														}
														else if (s2 == s3) {
															y++;
															xoverseqnumw[y + xoverwindow] = 0;
															xoverseqnumw[y + os1] = 0;
															xoverseqnumw[y + os2] = 1;
															ah2++;
															xdiffpos[y] = x;
															xposdiff[x] = y;
														}
														z = spacerno;
													}
												}
											}
										}

									}
								}
							}
						}
					}
				}

			}

		}
		else if (spacerflag>1) {


			for (x = 1; x < lenseq; x++) {
				s1 = seqnum[x + se1];
				s2 = seqnum[x + se2];
				s3 = seqnum[x + se3];
				xposdiff[x] = y;

				if (s1 != s2 || s1 != s3) {

					if (s1 != s2 && s1 != s3 && s2 != s3)
						g++;

					else {
						if (s1 != 46) {
							if (s2 != 46) {
								if (s3 != 46) {

									if (s1 != s2  && s1 != s3) {
										//If seq1 is odd one
										for (z = 1; z <= spacerno; z++) {
											sz = seqnum[x + spacerseqs[z] * lenseq];
											if (sz == s1) {
												//If difference is legitimate
												if (seq1 != outlyer)
													oc += 2;
												else {
													if (oc >0) {
														if (oc > hoc)
															hoc = oc;
														oc--;
													}
												}


												if (s1 == s2) {
													y++;
													xoverseqnumw[y + xoverwindow] = 1;
													xoverseqnumw[y + os1] = 0;
													xoverseqnumw[y + os2] = 0;
													ah0++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												else if (s1 == s3) {
													y++;
													xoverseqnumw[y + xoverwindow] = 0;
													xoverseqnumw[y + os1] = 1;
													xoverseqnumw[y + os2] = 0;
													ah1++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												else if (s2 == s3) {
													y++;
													xoverseqnumw[y + xoverwindow] = 0;
													xoverseqnumw[y + os1] = 0;
													xoverseqnumw[y + os2] = 1;
													ah2++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												break;
											}
										}
									}

									else if (s2 != s1  && s2 != s3) {
										//If seq2 is odd one
										for (z = 1; z <= spacerno; z++) {
											sz = seqnum[x + spacerseqs[z] * lenseq];

											if (s2 == sz) {
												//If difference is legitimate



												if (seq2 != outlyer)
													oc += 2;
												else {
													if (oc >0) {
														if (oc > hoc)
															hoc = oc;
														oc--;
													}
												}
												if (s1 == s2) {
													y++;
													xoverseqnumw[y + xoverwindow] = 1;
													xoverseqnumw[y + os1] = 0;
													xoverseqnumw[y + os2] = 0;
													ah0++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												else if (s1 == s3) {
													y++;
													xoverseqnumw[y + xoverwindow] = 0;
													xoverseqnumw[y + os1] = 1;
													xoverseqnumw[y + os2] = 0;
													ah1++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												else if (s2 == s3) {
													y++;
													xoverseqnumw[y + xoverwindow] = 0;
													xoverseqnumw[y + os1] = 0;
													xoverseqnumw[y + os2] = 1;
													ah2++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												break;
											}
										}
									}
									else if (s3 != s1  && s3 != s2) {
										//If seq3 is odd one
										for (z = 1; z <= spacerno; z++) {
											sz = seqnum[x + spacerseqs[z] * lenseq];

											if (s3 == sz) {

												//If difference is legitimate

												if (seq3 != outlyer)
													oc += 2;
												else {
													if (oc >0) {
														if (oc > hoc)
															hoc = oc;
														oc--;
													}
												}

												if (s1 == s2) {
													y++;
													xoverseqnumw[y + xoverwindow] = 1;
													xoverseqnumw[y + os1] = 0;
													xoverseqnumw[y + os2] = 0;
													ah0++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												else if (s1 == s3) {
													y++;
													xoverseqnumw[y + xoverwindow] = 0;
													xoverseqnumw[y + os1] = 1;
													xoverseqnumw[y + os2] = 0;
													ah1++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												else if (s2 == s3) {
													y++;
													xoverseqnumw[y + xoverwindow] = 0;
													xoverseqnumw[y + os1] = 0;
													xoverseqnumw[y + os2] = 1;
													ah2++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												break;
											}
										}
									}
								}
							}

						}
					}
				}
			}

		}

		wmx = y - xow;
		wpx = y + xow;
		for (b = 1; b <= xow; b++) {
			wmx++;
			wpx++;
			xoverseqnumw[b] = xoverseqnumw[wmx + xoverwindow];
			xoverseqnumw[b + lenseq + xoverwindow2] = xoverseqnumw[wmx + os1];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
			xoverseqnumw[b + lenseq2 + xoverwindow4] = xoverseqnumw[wmx + os2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

			xoverseqnumw[wpx] = xoverseqnumw[b + xoverwindow];//XOverSeqNum(X, 0)
			xoverseqnumw[wpx + lenseq + xoverwindow2] = xoverseqnumw[b + os1];
			xoverseqnumw[wpx + lenseq2 + xoverwindow4] = xoverseqnumw[b + os2];
		}
		//for (b = y+1; b < y + xoverwindow; b++)
		//	*(xdiffpos + b) = 0;


		//for (b = lenseq+1; b < lenseq + xoverwindow; b++)
		//	*(xposdiff + b) = 0;

		ah[0] = ah0;
		ah[1] = ah1;
		ah[2] = ah2;

		if (hoc < oc)
			hoc = oc;

		if (hoc < xoverwindow * 2)
			return (-y);
		else
			return(y);
	}

	
	int MyMathFuncs::FindSubSeqPB(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos,  char *xoverseqnumw, int *xdiffpos, int *xposdiff, unsigned char *FSSRDP)
	{

		int hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4, xow, wmx, wpx, target;
		int s1, s2, s3, sz, os1, os2, os3, se1, se2, se3, ah0, ah1, ah2,s1o,s2o,s3o, osf;//,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int g = 0;
		int lenseq;
		int h;
		int xh;

		lenseq = lenstrainseq0;
		oc = 0;
		hoc = 0;
		xow = xoverwindow;
		os3 = ubxos + 1;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		os1 = xoverwindow + lenseq + xoverwindow2;
		os2 = xoverwindow + lenseq2 + xoverwindow4;
		
		ah0 = 0;
		ah1 = 0;
		ah2 = 0;


		target = (lenseq + xoverwindow2) * 3;
		for (x = 1; x < target; x++)
			xoverseqnumw[x] = 0;

		hoc = xoverwindow * 10;
		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);
		xh = 0;
		for (x = 1; x <= ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];
			osf = se1*3 + (UBFSS + 1)*3*se2 + (UBFSS + 1)*(UBFSS + 1)*3*se3;
			for (z = 0; z <= 2; z++) {

				holder = FSSRDP[z + osf];
				xh++;
				if (holder > 0) {
					y = y + 1;
					xdiffpos[y] = xh;
					h = holder - 1;
					ah[h] = ah[h] + 1;
					//return(xow);
					xoverseqnumw[y + xow +h*os3] = 1;
					xposdiff[xh] = y;
				}
				else
					xposdiff[xh] = y;
			}
		}


		

		wmx = y - xow;
		wpx = y + xow;
		for (b = 1; b <= xow; b++) {
			wmx++;
			wpx++;
			xoverseqnumw[b] = xoverseqnumw[wmx + xow];
			xoverseqnumw[b + os3] = xoverseqnumw[wmx + + xow + os3];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
			xoverseqnumw[b + os3*2] = xoverseqnumw[wmx + xow + os3*2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

			xoverseqnumw[wpx] = xoverseqnumw[b + xow];//XOverSeqNum(X, 0)
			xoverseqnumw[wpx + os3] = xoverseqnumw[b + xow + os3];
			xoverseqnumw[wpx + os3*2] = xoverseqnumw[b + xow + os3*2];
		}

		if (hoc < oc)
			hoc = oc;

		if (hoc < xoverwindow * 2)
			return (-y);
		else
			return(y);
	}
	
	int MyMathFuncs::FindSubSeqPB4(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, int *xdiffpos, int *xposdiff, unsigned char *FSSRDP)
	{

		int hoc, oc, b, x, z, lenseq2;
		int s1, s2, s3, sz, os1, os2, os3, os4, os5, se1, se2, se3, ah0, ah1, ah2, s1o, s2o, s3o, osf;//,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int g = 0;
		int lenseq;
		int h;
		int xh;

		lenseq = lenstrainseq0;
		oc = 0;
		hoc = 0;
		os3 = ubxos + 1;
		lenseq2 = lenseq * 2;

		ah0 = 0;
		ah1 = 0;
		ah2 = 0;

		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);
		xh = 0;

		os4 = (UBFSS + 1) * 4;
		os5 = (UBFSS + 1)*(UBFSS + 1) * 4;

		for (x = 1; x <= ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];
			osf = se1 * 4 + os4 * se2 + os5 * se3;
			for (z = 0; z <= 2; z++) {

				holder = FSSRDP[z + osf];
				xh++;
				if (holder > 0) {
					y = y + 1;
					xdiffpos[y] = xh;
					xposdiff[xh] = y;
				}
				else
					xposdiff[xh] = y;
			}
		}




		return(y);
	}

	int MyMathFuncs::CleanXOSNW(int lenxoseq, int xoverwindow, int UBXO1, char *xoverseqnumw)
	{
		int target, x, y, os1, xoverwindow2;
		//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
		xoverwindow2 = xoverwindow * 2;
		target = (lenxoseq + xoverwindow2);
		
		for (x = 0; x < target; x++)
			xoverseqnumw[x] = 0;

		os1 = UBXO1+1;
		for (x = 0; x < target; x++)
			xoverseqnumw[x+os1] = 0;
		os1 = os1 * 2;
		for (x = 0; x < target; x++)
			xoverseqnumw[x+os1] = 0;


		return(1);
	}
	int MyMathFuncs::FindSubSeqPB2(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, unsigned char *FSSRDP)
	{

		int hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4, xow, wmx, wpx;
		int s1, s2, s3, sz, os1, os2, os3, os4, os5, se1, se2, se3, ah0, ah1, ah2, s1o, s2o, s3o, osf;//,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int g = 0;
		int lenseq;
		int h;
		

		lenseq = lenstrainseq0;
		oc = 0;
		hoc = 0;
		xow = xoverwindow;
		os3 = ubxos + 1;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		os1 = xoverwindow + lenseq + xoverwindow2;
		os2 = xoverwindow + lenseq2 + xoverwindow4;

		ah0 = 0;
		ah1 = 0;
		ah2 = 0;

		hoc = xoverwindow * 10;
		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);
		
		os4 = (UBFSS + 1) * 3;
		os5 = (UBFSS + 1)*(UBFSS + 1) * 3;
		for (x = 1; x <= ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];
			osf = se1 * 3 + os4 * se2 + os5 * se3;
			for (z = 0; z <= 2; z++) {

				holder = FSSRDP[z + osf];
				
				if (holder > 0) {
					y = y + 1;
					h = holder - 1;
					ah[h] = ah[h] + 1;
					//return(xow);
					xoverseqnumw[y + xow + h*os3] = 1;
				}
				
			}
		}




		wmx = y - xow;
		wpx = y + xow;
		for (b = 1; b <= xow; b++) {
			wmx++;
			wpx++;
			xoverseqnumw[b] = xoverseqnumw[wmx + xow];
			xoverseqnumw[b + os3] = xoverseqnumw[wmx + +xow + os3];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
			xoverseqnumw[b + os3 * 2] = xoverseqnumw[wmx + xow + os3 * 2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

			xoverseqnumw[wpx] = xoverseqnumw[b + xow];//XOverSeqNum(X, 0)
			xoverseqnumw[wpx + os3] = xoverseqnumw[b + xow + os3];
			xoverseqnumw[wpx + os3 * 2] = xoverseqnumw[b + xow + os3 * 2];
		}

		if (hoc < oc)
			hoc = oc;

		if (hoc < xoverwindow * 2)
			return (-y);
		else
			return(y);
	}

	int MyMathFuncs::FindSubSeqPB6(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, unsigned char *FSSRDP)
	{

		int b, x, z, lenseq2, xoverwindow2, xoverwindow4, xow, wmx, wpx;
		int  os1, os2, os3, os4, os5, se1, se2, se3, s1o, s2o, s3o, osf;//s1, s2, s3,hoc, oc,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int lenseq;
		int h;


		lenseq = lenstrainseq0;
		//oc = 0;
		//hoc = 0;
		xow = xoverwindow;
		os3 = ubxos + 1;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		os1 = xoverwindow + lenseq + xoverwindow2;
		os2 = xoverwindow + lenseq2 + xoverwindow4;

		//ah0 = 0;
		//ah1 = 0;
		//ah2 = 0;

		//hoc = xoverwindow * 10;
		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);

		os4 = (UBFSS + 1) * 4;
		os5 = (UBFSS + 1)*(UBFSS + 1) * 4;
		for (int x = 1; x <= ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];

			osf = se1 * 4 + os4 * se2 + os5 * se3;
			//if (FSSRDP[3 + osf] > 0) {
				holder = FSSRDP[3 + osf];
				//if (holder > 0) {
				for (int z = 0; z < holder; z++) {

					h = FSSRDP[z + osf]-1;

					//if (h > 0) {
						y++;
						//h--;
						ah[h]++;
						xoverseqnumw[y + xow + h*os3] = 1;
					//}
				}
			//}
		}




		wmx = y - xow;
		wpx = y + xow;
		for (b = 1; b <= xow; b++) {
			wmx++;
			wpx++;
			xoverseqnumw[b] = xoverseqnumw[wmx + xow];
			xoverseqnumw[b + os3] = xoverseqnumw[wmx + +xow + os3];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
			xoverseqnumw[b + os3 * 2] = xoverseqnumw[wmx + xow + os3 * 2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

			xoverseqnumw[wpx] = xoverseqnumw[b + xow];//XOverSeqNum(X, 0)
			xoverseqnumw[wpx + os3] = xoverseqnumw[b + xow + os3];
			xoverseqnumw[wpx + os3 * 2] = xoverseqnumw[b + xow + os3 * 2];
		}

		/*if (hoc < oc)
		hoc = oc;

		if (hoc < xoverwindow * 2)
		return (-y);
		else*/
		return(y);
	}
	//int MyMathFuncs::FindSubSeqPB6(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, unsigned char *FSSRDP)
	//{

	//	int hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4, xow, wmx, wpx;
	//	int s1, s2, s3, sz, os1, os2, os3, os4, os5, se1, se2, se3, ah0, ah1, ah2, s1o, s2o, s3o, osf;//,so1,so2,so3;
	//	int holder = 0;
	//	int y = 0;
	//	int g = 0;
	//	int lenseq;
	//	int h;


	//	lenseq = lenstrainseq0;
	//	oc = 0;
	//	hoc = 0;
	//	xow = xoverwindow;
	//	os3 = ubxos + 1;
	//	lenseq2 = lenseq * 2;
	//	xoverwindow2 = xoverwindow * 2;
	//	xoverwindow4 = xoverwindow * 4;
	//	os1 = xoverwindow + lenseq + xoverwindow2;
	//	os2 = xoverwindow + lenseq2 + xoverwindow4;

	//	ah0 = 0;
	//	ah1 = 0;
	//	ah2 = 0;

	//	hoc = xoverwindow * 10;
	//	s1o = seq1*(ubcs1 + 1);
	//	s2o = seq2*(ubcs1 + 1);
	//	s3o = seq3*(ubcs1 + 1);

	//	os4 = (UBFSS + 1) * 4;
	//	os5 = (UBFSS + 1)*(UBFSS + 1) * 4;
	//	for (x = 1; x <= ubcs1; x++) {
	//		se1 = CS[x + s1o];
	//		se2 = CS[x + s2o];
	//		se3 = CS[x + s3o];
	//		//return(s1o);
	//		osf = se1 * 4 + os4 * se2 + os5 * se3;
	//		if (FSSRDP[3 + osf] > 0) {
	//			for (z = 0; z <= 2; z++) {

	//				holder = FSSRDP[z + osf];

	//				if (holder > 0) {
	//					y = y + 1;
	//					h = holder - 1;
	//					ah[h] = ah[h] + 1;
	//					//return(xow);
	//					xoverseqnumw[y + xow + h*os3] = 1;
	//				}

	//			}
	//		}
	//	}




	//	wmx = y - xow;
	//	wpx = y + xow;
	//	for (b = 1; b <= xow; b++) {
	//		wmx++;
	//		wpx++;
	//		xoverseqnumw[b] = xoverseqnumw[wmx + xow];
	//		xoverseqnumw[b + os3] = xoverseqnumw[wmx + +xow + os3];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
	//		xoverseqnumw[b + os3 * 2] = xoverseqnumw[wmx + xow + os3 * 2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

	//		xoverseqnumw[wpx] = xoverseqnumw[b + xow];//XOverSeqNum(X, 0)
	//		xoverseqnumw[wpx + os3] = xoverseqnumw[b + xow + os3];
	//		xoverseqnumw[wpx + os3 * 2] = xoverseqnumw[b + xow + os3 * 2];
	//	}

	//	if (hoc < oc)
	//		hoc = oc;

	//	if (hoc < xoverwindow * 2)
	//		return (-y);
	//	else
	//		return(y);
	//}

	

	int MyMathFuncs::FindSubSeqPB3(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, unsigned char *FSSRDP)
	{

		int b, x, z, lenseq2, xoverwindow2, xoverwindow4, xow, wmx, wpx;
		int  os1, os2, os3, os4, os5, se1, se2, se3, s1o, s2o, s3o, osf;//s1, s2, s3,hoc, oc,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int lenseq;
		int h;


		lenseq = lenstrainseq0;
		//oc = 0;
		//hoc = 0;
		xow = xoverwindow;
		os3 = ubxos + 1;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		os1 = xoverwindow + lenseq + xoverwindow2;
		os2 = xoverwindow + lenseq2 + xoverwindow4;

		//ah0 = 0;
		//ah1 = 0;
		//ah2 = 0;

		//hoc = xoverwindow * 10;
		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);

		os4 = (UBFSS + 1) * 4;
		os5 = (UBFSS + 1)*(UBFSS + 1) * 4;
		for (int x = 1; x <= ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];
			
			osf = se1 * 4 + os4 * se2 + os5 * se3;
			if (FSSRDP[3 + osf] > 0) {
			//holder = FSSRDP[3 + osf];
			//if (holder > 0) {
				for (int z = 0; z <= 2; z++) {

					h = FSSRDP[z + osf];

					if (h > 0) {
						y++;
						h--;
						ah[h]++;
						xoverseqnumw[y + xow + h*os3] = 1;
					}
				}
			}
		}




		wmx = y - xow;
		wpx = y + xow;
		for (b = 1; b <= xow; b++) {
			wmx++;
			wpx++;
			xoverseqnumw[b] = xoverseqnumw[wmx + xow];
			xoverseqnumw[b + os3] = xoverseqnumw[wmx + +xow + os3];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
			xoverseqnumw[b + os3 * 2] = xoverseqnumw[wmx + xow + os3 * 2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

			xoverseqnumw[wpx] = xoverseqnumw[b + xow];//XOverSeqNum(X, 0)
			xoverseqnumw[wpx + os3] = xoverseqnumw[b + xow + os3];
			xoverseqnumw[wpx + os3 * 2] = xoverseqnumw[b + xow + os3 * 2];
		}

		/*if (hoc < oc)
			hoc = oc;

		if (hoc < xoverwindow * 2)
			return (-y);
		else*/
			return(y);
	}

	int MyMathFuncs::FindSubSeqPB5(int *ah, int UBFSS, int xoverwindow, int ubcs1, int lenstrainseq0, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int ubxos, char *xoverseqnumw, unsigned char *FSSRDP)
	{

		int hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4, xow, wmx, wpx;
		int s1, s2, s3, sz, os1, os2, os3, os4, os5, se1, se2, se3, ah0, ah1, ah2, s1o, s2o, s3o, osf;//,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int g = 0;
		int lenseq;
		int h;


		lenseq = lenstrainseq0;
		oc = 0;
		hoc = 0;
		xow = xoverwindow;
		os3 = ubxos + 1;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		os1 = xoverwindow + lenseq + xoverwindow2;
		os2 = xoverwindow + lenseq2 + xoverwindow4;

		ah0 = 0;
		ah1 = 0;
		ah2 = 0;

		hoc = xoverwindow * 10;
		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);

		os4 = (UBFSS + 1) * 4;
		os5 = (UBFSS + 1)*(UBFSS + 1) * 4;
		for (x = 1; x <= ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];
			osf = se1 * 4 + os4 * se2 + os5 * se3;
			if (FSSRDP[3 + osf] > 0) {
				for (z = 0; z <= 2; z++) {

					holder = FSSRDP[z + osf];

					if (holder > 0) {
						y = y + 1;
						//h = holder - 1;
						ah[holder - 1] = ah[holder - 1] + 1;
						//return(xow);
						xoverseqnumw[y + xow] = holder;
					}

				}
			}
		}
		wmx = y - xow;
		wpx = y + xow;
		for (b = 1; b <= xow; b++) {
			wmx++;
			wpx++;
			xoverseqnumw[b] = xoverseqnumw[wmx + xow];
			//xoverseqnumw[b + os3] = xoverseqnumw[wmx + +xow + os3];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
			//xoverseqnumw[b + os3 * 2] = xoverseqnumw[wmx + xow + os3 * 2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

			xoverseqnumw[wpx] = xoverseqnumw[b + xow];//XOverSeqNum(X, 0)
			//xoverseqnumw[wpx + os3] = xoverseqnumw[b + xow + os3];
			//xoverseqnumw[wpx + os3 * 2] = xoverseqnumw[b + xow + os3 * 2];
		}

		if (hoc < oc)
			hoc = oc;

		if (hoc < xoverwindow * 2)
			return (-y);
		else
			return(y);
	}


	


	int MyMathFuncs::FindSubSeqMCPB(int UBFSS, int ubcs1, int nextno, int seq1, int seq2, int seq3, unsigned char *CS, int *xdiffpos, int *xposdiff, unsigned char *FSSRDP)
	{

		int x, z;
		int se1, se2, se3, s1o, s2o, s3o, osf;//,so1,so2,so3;
		int y = 0;
		int xh;
		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);
		xh = 0;
		for (x = 1; x <= ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];
			osf = se1 * 3 + (UBFSS + 1) * 3 * se2 + (UBFSS + 1)*(UBFSS + 1) * 3 * se3;
			for (z = 0; z <= 2; z++) {
				xh++;
				if (FSSRDP[z + osf]) {
					y++;
					xdiffpos[y] = xh;
					xposdiff[xh] = y;
				}
				else
					xposdiff[xh] = y;
			}
		}




			return(y);
	}

	int MyMathFuncs::FindSubSeqMCPB2(int UBFSS, int ubcs1, int nextno, int seq1, int seq2, int seq3, unsigned char *CS,  unsigned char *FSSRDP, int *XDiffPos)
	{

		int x, z;
		int se1, se2, se3, s1o, s2o, s3o, osf;//,so1,so2,so3;
		int y = 0;
		int xh;
		s1o = seq1*(ubcs1 + 1);
		s2o = seq2*(ubcs1 + 1);
		s3o = seq3*(ubcs1 + 1);
		xh = 0;
		for (x = 1; x <= ubcs1; x++) {
			se1 = CS[x + s1o];
			se2 = CS[x + s2o];
			se3 = CS[x + s3o];
			osf = se1 * 3 + (UBFSS + 1) * 3 * se2 + (UBFSS + 1)*(UBFSS + 1) * 3 * se3;
			for (z = 0; z <= 2; z++) {
				xh++;
				if (FSSRDP[z + osf]) {
					y++;
					XDiffPos[y] = xh;
				
				}
				
			}
		}




		return(y);
	}

	int MyMathFuncs::FindSubSeqP4(int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int seq3, short int spacerno, short int *seqnum,  unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray)
	{

		int se2os, se3os, hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4,  wmx, wpx, target;
		int s1, s2, s3, sz,   ah0, ah1, ah2;//,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int g = 0;

		const int se1 = seq1*lenseq;
		const int se2 = seq2*lenseq;
		const int se3 = seq3*lenseq;
		oc = 0;
		hoc = 0;
		const int xow = xoverwindow;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		const int os1 = xoverwindow + lenseq + xoverwindow2;
		const int os2 = xoverwindow + lenseq2 + xoverwindow4;
		
		ah0 = 0;
		ah1 = 0;
		ah2 = 0;

		target = (lenseq + xoverwindow2 + 1) * 3;
		
		for (x = 1; x < target; x++)
			xoverseqnumw[x] = 0;

		if (spacerflag == 0) {
			//so1 = seq1*lenseq;
			//so2 = seq2*lenseq;
			//so3 = seq3*lenseq;
			
			y = 0;
			hoc = xoverwindow * 10;
			for (x = 1; x < lenseq; x++) {
				se2os = se2 + x;
				if (binarray[se2os] == 1) {//seq1 and seq2 are different
					se3os = se3 + x;	
					if (binarray[se3os] == 1) {//if seq1 is also different to seq3{
										//are seq2 and seq3 the same?
						if (seqnum[se2os] == seqnum[se3os]) {
							/*if (seqnum[x + se2] != 46) {
								if (seqnum[x + se3] != 46) {*/
									xdiffpos[++y] = x;
									xoverseqnumw[os2 + y] = 1;
									ah2++;
								/*}
							}*/

						}
					}

					else if (binarray[se3os] == 0) {//seq1 and seq3 are the same

						xdiffpos[++y] = x;
						xoverseqnumw[os1 + y] = 1;
						ah1++;

					}


				}
				else if (binarray[se3+x] == 1 ) {//seq1 and seq3 are different but seq1=seq2
					if (binarray[se2os] == 0) {
						xdiffpos[++y] = x;
						xoverseqnumw[xow + y] = 1;
						ah0++;
					}

				}
				xposdiff[x] = y;

			}
			//int hold;
			//for (x = 1; x < lenseq; x++) {
			//	if (binarray[x + se2]) {//seq1 and seq2 are different
			//		
			//		if (binarray[x + se3]) {//if seq1 is also different to seq3{
			//								are seq2 and seq3 the same?
			//			hold = (int)(seqnum[x + se2] == seqnum[x + se3]);
			//			y+= hold;
			//			xdiffpos[y] = x;							
			//			xoverseqnumw[y + os2]=hold;
			//			ah2+=hold;
			//				
			//			
			//		}

			//		else {//seq1 and seq3 are the same
			//			
			//			xdiffpos[y++] = x;
			//			xoverseqnumw[y + os1]=1;
			//			ah1++;
			//			
			//		}


			//	}
			//	else if (binarray[x + se3]) {//seq1 and seq3 are different but seq1=seq2
			//		xdiffpos[y++] = x;
			//		xoverseqnumw[y + xoverwindow]=1;
			//		ah0++;
			//		
			//	}
			//	xposdiff[x] = y;

			//}
			//
			
			
		}

		else if (spacerflag == 1) {
			for (x = 1; x < lenseq; x++) {
				*(xposdiff + x) = y;
				if (binarray[x + se2] || binarray[x + se3]) {

					s1 = *(seqnum + x + se1);
					s2 = *(seqnum + x + se2);
					s3 = *(seqnum + x + se3);

					//if (s1 != s2 || s1 != s3) {

						if (s1 == s2 || s1 == s3 || s2 == s3) {

							//if (s1 != 46) {
							//	if (s2 != 46) {
							//		if (s3 != 46) {
										if (binarray[x + se2] && binarray[x + se3]) {
											//If seq1 is odd one


											if (outlyer == seq1) {
												if (oc > 0) {
													if (oc > hoc)
														hoc = oc;
													oc--;
												}
												
												
												if (s2 == s3) {
													y++;
													
													xoverseqnumw[y + os2] = 1;
													ah2++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}

											}
											else {
												for (z = 1; z <= spacerno; z++) {
													sz = *(seqnum + x + *(spacerseqs + z) * lenseq);
													if (sz == s1) {
														//If difference is legitimate
														if (*(xdiffpos + y) != x) {
															oc += 2;
															
															
															 if (s2 == s3) {
																y++;
																
																xoverseqnumw[y + os2] = 1;
																ah2++;
																xdiffpos[y] = x;
																xposdiff[x] = y;
															}
															z = spacerno;

														}

													}

												}

											}
										}
										else if (s2 != s1  && s2 != s3) {
											//If seq2 is odd one

											if (outlyer == seq2) {
												if (oc > 0) {
													if (oc > hoc)
														hoc = oc;
													oc--;
												}
												
												if (s1 == s3) {
													y++;
													
													xoverseqnumw[y + os1] = 1;
													
													ah1++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												

											}
											else {
												for (z = 1; z <= spacerno; z++) {
													sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
													if (s2 == sz) {
														//If difference is legitimate

														if (*(xdiffpos + y) != x) {
															oc += 2;
															
															if (s1 == s3) {
																y++;
																
																xoverseqnumw[y + os1] = 1;
																
																ah1++;
																xdiffpos[y] = x;
																xposdiff[x] = y;
															}
															
															z = spacerno;

														}
													}
												}
											}
										}
										else if (s3 != s1  && s3 != s2) {
											//If seq3 is odd one
											if (outlyer == seq3) {
												if (oc > 0) {
													if (oc > hoc)
														hoc = oc;
													oc--;
												}
												if (s1 == s2) {
													y++;
													xoverseqnumw[y + xoverwindow] = 1;
													
													ah0++;
													xdiffpos[y] = x;
													xposdiff[x] = y;
												}
												


											}
											else {
												for (z = 1; z <= spacerno; z++) {
													sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
													if (s3 == sz) {

														//If difference is legitimate
														if (*(xdiffpos + y) != x) {
															oc += 2;
															if (s1 == s2) {
																y++;
																xoverseqnumw[y + xoverwindow] = 1;
																
																ah0++;
																xdiffpos[y] = x;
																xposdiff[x] = y;
															}
															
															z = spacerno;
														}
													}
												}
											}

										}
									//}
								//}
							//}
						}
					//}
				}
			}

		}
		else if (spacerflag>1) {


			for (x = 1; x < lenseq; x++) {
				xposdiff[x] = y;
				if (binarray[x + se2] || binarray[x + se3]) {
					s1 = seqnum[x + se1];
					s2 = seqnum[x + se2];
					s3 = seqnum[x + se3];


					//if (s1 != s2 || s1 != s3) {

						if (s1 != s2 && s1 != s3 && s2 != s3)
							g++;

						else {
							//if (s1 != 46) {
								//if (s2 != 46) {
									//if (s3 != 46) {

										if (s1 != s2  && s1 != s3) {
											//If seq1 is odd one
											for (z = 1; z <= spacerno; z++) {
												sz = seqnum[x + spacerseqs[z] * lenseq];
												if (sz == s1) {
													//If difference is legitimate
													if (seq1 != outlyer)
														oc += 2;
													else {
														if (oc > 0) {
															if (oc > hoc)
																hoc = oc;
															oc--;
														}
													}


													if (s2 == s3) {
														y++;
														
														xoverseqnumw[y + os2] = 1;
														ah2++;
														xdiffpos[y] = x;
														xposdiff[x] = y;
													}
													break;
												}
											}
										}

										else if (s2 != s1  && s2 != s3) {
											//If seq2 is odd one
											for (z = 1; z <= spacerno; z++) {
												sz = seqnum[x + spacerseqs[z] * lenseq];

												if (s2 == sz) {
													//If difference is legitimate



													if (seq2 != outlyer)
														oc += 2;
													else {
														if (oc > 0) {
															if (oc > hoc)
																hoc = oc;
															oc--;
														}
													}
													if (s1 == s3) {
														y++;
														
														xoverseqnumw[y + os1] = 1;
														
														ah1++;
														xdiffpos[y] = x;
														xposdiff[x] = y;
													}
													
													break;
												}
											}
										}
										else if (s3 != s1  && s3 != s2) {
											//If seq3 is odd one
											for (z = 1; z <= spacerno; z++) {
												sz = seqnum[x + spacerseqs[z] * lenseq];

												if (s3 == sz) {

													//If difference is legitimate

													if (seq3 != outlyer)
														oc += 2;
													else {
														if (oc > 0) {
															if (oc > hoc)
																hoc = oc;
															oc--;
														}
													}

													if (s1 == s2) {
														y++;
														xoverseqnumw[y + xoverwindow] = 1;
														
														ah0++;
														xdiffpos[y] = x;
														xposdiff[x] = y;
													}
													
													break;
												}
											}
										//}
									//}
								//}

							}
						}
					//}
				}
			}

		}

		wmx = y - xow;
		wpx = y + xow;
		for (b = 1; b <= xow; b++) {
			wmx++;
			wpx++;
			xoverseqnumw[b] = xoverseqnumw[wmx + xoverwindow];
			xoverseqnumw[b + lenseq + xoverwindow2] = xoverseqnumw[wmx + os1];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
			xoverseqnumw[b + lenseq2 + xoverwindow4] = xoverseqnumw[wmx + os2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

			xoverseqnumw[wpx] = xoverseqnumw[b + xoverwindow];//XOverSeqNum(X, 0)
			xoverseqnumw[wpx + lenseq + xoverwindow2] = xoverseqnumw[b + os1];
			xoverseqnumw[wpx + lenseq2 + xoverwindow4] = xoverseqnumw[b + os2];
		}

		
		//for (b = y+1; b < y + xoverwindow; b++)
		//	*(xdiffpos + b) = 0;


		//for (b = lenseq+1; b < lenseq + xoverwindow; b++)
		//	*(xposdiff + b) = 0;

		ah[0] = ah0;
		ah[1] = ah1;
		ah[2] = ah2;

		if (hoc < oc)
			hoc = oc;

		if (hoc < xoverwindow * 2)
			return (-y);
		else
			return(y);
	}


	int MyMathFuncs::FillSetsP(int SZ1, int lseq, int Nextno, int UBCX, int UBXO1, int UBRL1, int UBS, int *RNum, int *RList, short int *CurrentXOver, XOVERDEFINE *XOverlist, int *RI,  char *OLSeq, unsigned char *Sets) {
		int X, Y, UB, os, os2, os3, os4;
		int A, Z, TS, OL, SZ2;
		float sz1h, sz2h, olh;
		os3 = UBRL1 + 1;
		os4 = UBS + 1;
		os = UBXO1 + 1;

		if (Nextno <= UBCX)
			UB = Nextno;
		else
			UB = UBCX;
		for (X = 0; X <= UB; X++) {

			
				
			for (Y = 1; Y <= CurrentXOver[X]; Y++) {
				os2 = X + Y*os;
				RI[4] = XOverlist[os2].Beginning;
				RI[5] = XOverlist[os2].Ending;
					
				//GoOn = DoSetsB(SZ1, LSeq, RI(0), OLSeq(0));
				
				OL = 0;
				if (RI[4] < RI[5]) {
					SZ2 = RI[5] - RI[4] + 1;
					for (A = RI[4]; A <= RI[5]; A++)
						OL = OL + OLSeq[A];

				}
				else {
					SZ2 = RI[5] + lseq - RI[4] + 1;
					for (A = RI[4]; A <= lseq; A++)
						OL = OL + OLSeq[A];
					for (A = 1; A <= RI[5]; A++)
						OL = OL + OLSeq[A];

				}
				//return(OL);
				sz1h = (float)(SZ1);
				sz2h = (float)(SZ2);
				olh = (float)(OL);

				if (olh / ((sz1h + sz2h) / 2) > 0.3) {



					RI[0] = XOverlist[os2].Daughter;
					RI[1] = XOverlist[os2].MajorP;
					RI[2] = XOverlist[os2].MinorP;
						

						
					for (Z = 0; Z <= 2; Z++) {

						for (A = 0; A <= RNum[Z]; A++) {
							TS = RList[Z + A*os3];
							if (RI[0] == TS || RI[1] == TS || RI[2] == TS) {
								Sets[Z + RI[0]*os4] = 1;
								Sets[Z + RI[1] * os4] = 1;
								Sets[Z + RI[2] * os4] = 1;
							}
						}
					}
						
				}
			}
				
		}
		return(1);
	}

	int MyMathFuncs::FillSetsP2(int UBXO2, int SZ1, int lseq, int Nextno, int UBCX, int UBXO1, int UBRL1, int UBS, int *RNum, int *RList, short int *CurrentXOver, XOVERDEFINE *XOverlist, int *RI, char *OLSeq, unsigned char *Sets) {
		int X, Y, UB, os, os2, os3, os4, os5;
		int A, Z, TS, OL, SZ2;
		float sz1h, sz2h, olh;
		os3 = UBRL1 + 1;
		os4 = UBS + 1;
		os = UBXO1 + 1;

		if (Nextno <= UBCX)
			UB = Nextno;
		else
			UB = UBCX;
		for (Y = 1; Y <= UBXO2; Y++) {
			os5 = Y*os;
			for (X = 0; X <= UB; X++) {
				if (Y <= CurrentXOver[X]) {



					os2 = X + os5;
					RI[4] = XOverlist[os2].Beginning;
					RI[5] = XOverlist[os2].Ending;

					//GoOn = DoSetsB(SZ1, LSeq, RI(0), OLSeq(0));

					OL = 0;
					if (RI[4] < RI[5]) {
						SZ2 = RI[5] - RI[4] + 1;
						for (A = RI[4]; A <= RI[5]; A++)
							OL = OL + OLSeq[A];

					}
					else {
						SZ2 = RI[5] + lseq - RI[4] + 1;
						for (A = RI[4]; A <= lseq; A++)
							OL = OL + OLSeq[A];
						for (A = 1; A <= RI[5]; A++)
							OL = OL + OLSeq[A];

					}
					//return(OL);
					sz1h = (float)(SZ1);
					sz2h = (float)(SZ2);
					olh = (float)(OL);

					if (olh / ((sz1h + sz2h) / 2) > 0.3) {



						RI[0] = XOverlist[os2].Daughter;
						RI[1] = XOverlist[os2].MajorP;
						RI[2] = XOverlist[os2].MinorP;



						for (Z = 0; Z <= 2; Z++) {

							for (A = 0; A <= RNum[Z]; A++) {
								TS = RList[Z + A*os3];
								if (RI[0] == TS || RI[1] == TS || RI[2] == TS) {
									Sets[Z + RI[0] * os4] = 1;
									Sets[Z + RI[1] * os4] = 1;
									Sets[Z + RI[2] * os4] = 1;
								}
							}
						}

					}
				}
			}

		}
		return(1);
	}
	int MyMathFuncs::FillSetsP3(int BE, int EN,int UBXO2, int SZ1, int lseq, int Nextno, int UBCX, int UBXO1, int UBRL1, int UBS, int *RNum, int *RList, short int *CurrentXOver, XOVERDEFINE *XOverlist, int *RI, char *OLSeq, unsigned char *Sets) {
		int OL2, X, Y, UB, os, os2, os3, os4, os5;
		int A, Z, TS, OL, SZ2;
		float sz1h, sz2h, olh;
		os3 = UBRL1 + 1;
		os4 = UBS + 1;
		os = UBXO1 + 1;

		if (Nextno <= UBCX)
			UB = Nextno;
		else
			UB = UBCX;


		if (BE < EN) {
			for (Y = 1; Y <= UBXO2; Y++) {
				os5 = Y*os;
				for (X = 0; X <= UB; X++) {
					if (Y <= CurrentXOver[X]) {



						os2 = X + os5;
						RI[4] = XOverlist[os2].Beginning;
						RI[5] = XOverlist[os2].Ending;

						//GoOn = DoSetsB(SZ1, LSeq, RI(0), OLSeq(0));

						OL = 0;
						//OL2 = 0;
						if (RI[4] < RI[5]) {
							SZ2 = RI[5] - RI[4] + 1;
							if (RI[5] >= BE){
								if (RI[4] <= EN) {
									if (EN >= RI[5] && BE <= RI[4])
										OL = SZ2;
									else if (EN <= RI[5] && BE >= RI[4])
										OL = SZ1;
									else if (EN <= RI[5])
										OL = EN - RI[4] + 1;
									else if (BE >= RI[4])
										OL = RI[5] - BE + 1;
								}
							}

							//SZ2 = RI[5] - RI[4] + 1;
							/*for (A = RI[4]; A <= RI[5]; A++)
								OL2 = OL2 + OLSeq[A];

							if (OL != OL2)
								OL = OL2;*/

						}
						else {
							SZ2 = RI[5] + lseq - RI[4] + 1;
							if (RI[5] >= BE){
								
								if (EN >= RI[5])
									OL = RI[5] - BE + 1;
								else 
									OL = SZ1;
							}
							if (RI[4] <= EN) {

								if (BE <= RI[4])
									OL = EN - RI[4] + 1;
								else
									OL = SZ1;
							}
							
							/*for (A = RI[4]; A <= lseq; A++)
								OL2 = OL2 + OLSeq[A];
							for (A = 1; A <= RI[5]; A++)
								OL2 = OL2 + OLSeq[A];

							if (OL != OL2)
								OL = OL2; */
						}
						//return(OL);
						sz1h = (float)(SZ1);
						sz2h = (float)(SZ2);
						olh = (float)(OL);

						if (olh / ((sz1h + sz2h) / 2) > 0.3) {



							RI[0] = XOverlist[os2].Daughter;
							RI[1] = XOverlist[os2].MajorP;
							RI[2] = XOverlist[os2].MinorP;



							for (Z = 0; Z <= 2; Z++) {

								for (A = 0; A <= RNum[Z]; A++) {
									TS = RList[Z + A*os3];
									if (RI[0] == TS || RI[1] == TS || RI[2] == TS) {
										Sets[Z + RI[0] * os4] = 1;
										Sets[Z + RI[1] * os4] = 1;
										Sets[Z + RI[2] * os4] = 1;
									}
								}
							}

						}
					}
				}

			}
		}
		else {
			for (Y = 1; Y <= UBXO2; Y++) {
				os5 = Y*os;
				for (X = 0; X <= UB; X++) {
					if (Y <= CurrentXOver[X]) {



						os2 = X + os5;
						RI[4] = XOverlist[os2].Beginning;
						RI[5] = XOverlist[os2].Ending;

						//GoOn = DoSetsB(SZ1, LSeq, RI(0), OLSeq(0));

						OL = 0;
						//OL2 = 0;
						if (RI[4] < RI[5]) {
							SZ2 = RI[5] - RI[4] + 1;


							if (EN >= RI[4]){
								
								if (EN >= RI[5])
									OL = SZ2;
								else 
									OL = EN - RI[4] + 1;;
							}
							if (BE <= RI[5]) {

								if (BE <= RI[4])
									OL = OL + SZ2;
								else
									OL=OL+RI[5] - BE + 1;
							}



							/*for (A = RI[4]; A <= RI[5]; A++)
								OL2 = OL2 + OLSeq[A];

							if (OL != OL2)
								OL = OL2; */

						}
						else if (RI[4] > RI[5]) {
							
							SZ2 = RI[5] + lseq - RI[4] + 1;
							if (RI[5] <= EN)
								OL = RI[5];
							else
								OL = EN;

							if (RI[4] <= BE)
								OL = OL + lseq - BE + 1;
							else
								OL = OL + lseq - RI[4] + 1;

							/*for (A = RI[4]; A <= lseq; A++)
								OL2 = OL2 + OLSeq[A];
							for (A = 1; A <= RI[5]; A++)
								OL2 = OL2 + OLSeq[A];

							if (OL != OL2)
								OL = OL2;*/

						}
						//return(OL);
						sz1h = (float)(SZ1);
						sz2h = (float)(SZ2);
						olh = (float)(OL);

						if (olh / ((sz1h + sz2h) / 2) > 0.3) {



							RI[0] = XOverlist[os2].Daughter;
							RI[1] = XOverlist[os2].MajorP;
							RI[2] = XOverlist[os2].MinorP;



							for (Z = 0; Z <= 2; Z++) {

								for (A = 0; A <= RNum[Z]; A++) {
									TS = RList[Z + A*os3];
									if (RI[0] == TS || RI[1] == TS || RI[2] == TS) {
										Sets[Z + RI[0] * os4] = 1;
										Sets[Z + RI[1] * os4] = 1;
										Sets[Z + RI[2] * os4] = 1;
									}
								}
							}

						}
					}
				}

			}
		}

		return(1);
	}
	int MyMathFuncs::DoSetsAP(int Nextno, int UBCX, int UBXO1, int SZ1, int lseq, int *RI, char *OLSeq, char *Sets,  short int *CurrentXOver, XOVERDEFINE *XOverlist, unsigned char *DoIt, int *ISeqs) {
		//sets -2,nextno
		int X, Y, Z, A, OL, SZ2, UB, GoOn, os, os2;
		unsigned int ri0, ri1, ri2;
		float sz1h, sz2h, olh;
		os = UBXO1 + 1;
		
		if (Nextno <= UBCX)
			UB = Nextno;
		else
			UB = UBCX;

		for (X = 0; X <= UB; X++) {
			
			for (Y = 1; Y <= CurrentXOver[X]; Y++) {
				os2 = X + Y*os;
				RI[0] = XOverlist[os2].Daughter;
				RI[1] = XOverlist[os2].MajorP;
				RI[2] = XOverlist[os2].MinorP;
				GoOn = 0;
				DoIt[0] = 0; 
				DoIt[1] = 0; 
				DoIt[2] = 0;
				for (Z = 0; Z <= 2; Z++) {

					if (RI[0] == ISeqs[Z] || RI[1] == ISeqs[Z] || RI[2] == ISeqs[Z]) {
						GoOn = 1;
						DoIt[Z] = 1;
					}
				}
				if (GoOn == 1) {

					RI[4] = XOverlist[os2].Beginning;
					RI[5] = XOverlist[os2].Ending;
					ri0 = RI[0] * 3;
					ri1 = RI[1] * 3;
					ri2 = RI[2] * 3;
					
					OL = 0;
					if (RI[4] < RI[5]) {
						SZ2 = RI[5] - RI[4] + 1;
						for (A = RI[4]; A <= RI[5]; A++)
							OL = OL + OLSeq[A];

					}
					else {
						SZ2 = RI[5] + lseq - RI[4] + 1;
						for (A = RI[4]; A <= lseq; A++)
							OL = OL + OLSeq[A];
						for (A = 1; A <= RI[5]; A++)
							OL = OL + OLSeq[A];

					}
					sz1h = (float)(SZ1);
					sz2h = (float)(SZ2);
					olh = (float)(OL);
					if (olh / ((sz1h + sz2h) / 2) > 0.3) {
						for (Z = 0; Z < 3; Z++) {
							if (DoIt[Z] == 1) {
								Sets[Z + ri0] = 1;
								Sets[Z + ri1] = 1;
								Sets[Z + ri2] = 1;
							}
						}
					}
				}
			}
		}
		return(1);
	}


	//int MyMathFuncs::DoSetsAP2(int Nextno, int UBCX, int UBXO1, int SZ1, int lseq, int *RI, char *OLSeq, char *Sets, char *doit, short int *CurrentXOver, XOVERDEFINE *XOverlist, unsigned char *DoIt, int *ISeqs) {
	//	//sets -2,nextno
	//	int X, Y, Z, A, OL, SZ2, UB, GoOn, os, os2;
	//	unsigned int ri0, ri1, ri2;
	//	float sz1h, sz2h, olh;
	//	os = UBXO1 + 1;

	//	if (Nextno <= UBCX)
	//		UB = Nextno;
	//	else
	//		UB = UBCX;

	//	for (X = 0; X <= UB; X++) {

	//		for (Y = 1; Y <= CurrentXOver[X]; Y++) {
	//			os2 = X + Y*os;
	//			RI[0] = XOverlist[os2].Daughter;
	//			RI[1] = XOverlist[os2].MajorP;
	//			RI[2] = XOverlist[os2].MinorP;
	//			GoOn = 0;
	//			DoIt[0] = 0;
	//			DoIt[1] = 0;
	//			DoIt[2] = 0;
	//			for (Z = 0; Z <= 2; Z++) {

	//				if (RI[0] == ISeqs[Z] || RI[1] == ISeqs[Z] || RI[2] == ISeqs[Z]) {
	//					GoOn = 1;
	//					DoIt[Z] = 1;
	//				}
	//			}
	//			if (GoOn == 1) {

	//				RI[4] = XOverlist[os2].Beginning;
	//				RI[5] = XOverlist[os2].Ending;
	//				ri0 = RI[0] * 3;
	//				ri1 = RI[1] * 3;
	//				ri2 = RI[2] * 3;

	//				OL = 0;
	//				if (RI[4] < RI[5]) {
	//					SZ2 = RI[5] - RI[4] + 1;
	//					for (A = RI[4]; A <= RI[5]; A++)
	//						OL = OL + OLSeq[A];

	//				}
	//				else {
	//					SZ2 = RI[5] + lseq - RI[4] + 1;
	//					for (A = RI[4]; A <= lseq; A++)
	//						OL = OL + OLSeq[A];
	//					for (A = 1; A <= RI[5]; A++)
	//						OL = OL + OLSeq[A];

	//				}
	//				sz1h = (float)(SZ1);
	//				sz2h = (float)(SZ2);
	//				olh = (float)(OL);
	//				if (olh / ((sz1h + sz2h) / 2) > 0.3) {
	//					for (Z = 0; Z < 3; Z++) {
	//						if (doit[Z] == 1) {
	//							Sets[Z + ri0] = 1;
	//							Sets[Z + ri1] = 1;
	//							Sets[Z + ri2] = 1;
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//	return(1);
	//}
	int MyMathFuncs::CleanRedo(int RWNN, int x, int NextNo, int *RedoLS, int UBRL, int *RedoList) {
		int Z, A, GO, RedoListSize, os;
		Z = 0;
		os = UBRL + 1;
		RedoListSize = *RedoLS;
        while (Z <= RedoListSize){
			GO = 0;
			for (A = 1; A <= 3; A++) {
				if (RedoList[A + Z*os] == x) {
					GO = 1;
					break;
				}
			}
			if (GO == 1) {
				for (A = 0; A <= 3; A++)
					RedoList[A + Z*os] = RedoList[A + RedoListSize*os];
				RedoListSize--;
			}
            else
				Z++;
        }
		if (RWNN == 1) {
			for (Z = 0; Z <= RedoListSize; Z++) {
				for (A = 1; A <= 3; A++) {
					if (RedoList[A + Z*os] == NextNo) {
						RedoList[A + Z*os] = x;
						break;
					}
				}
			}
		}
		*RedoLS = RedoListSize;
		return(1);
	}


	int MyMathFuncs::MakeMoveDist(int NextNo, float *MoveDistF,float *MoveDistS, int UBFM, float *FMat, int UBSM, float *SMat) {
		int x, Y,osf, osf2, oss2, oss;
		osf = UBFM + 1;
		oss = UBSM + 1;

		for (x = 0; x <=NextNo; x++){
			oss2 = oss*x;
			osf2 = osf*x;
			for (Y = 0; Y <= NextNo; Y++){
				MoveDistF[x] = MoveDistF[x] + FMat[Y+osf2];
				MoveDistS[x] = MoveDistS[x] + SMat[Y + oss2];
			}
		}
		return(1);
	}


	int MyMathFuncs::MakeMatchMatX2P(int LSeq, int NextNo, int XX, char *ContainSite, float *SMat, float *MatchMat, float *BMatch, int *BPMatch, short int *SeqNum, int *iseqs) {
		//bpmatch 2,1,nexto
		//bmatch 2,nextno
		//containsite lseq,nextno
		//matchmat 2,nextno,nextno

		int X, AA, A, Y, Z, off1, ST, EN, s1, s2;
		float Diffs, Valid, th1, th2, th3;
		X = XX;
		off1 = LSeq + 1;
#pragma omp parallel for private(Y, ST, EN, Z)
		for (Y = 0; Y <= NextNo; Y++) {

			if (BMatch[X + Y * 3] > 0) {
				ST = BPMatch[X + Y * 6];
				if (ST >LSeq)
					ST = LSeq;
				EN = BPMatch[X + 3 + Y * 6];
				if (ST < EN) {
					for (Z = ST; Z <= EN; Z++)
						ContainSite[Z + Y*off1] = 1;
				}
				else {
					for (Z = ST; Z <= LSeq; Z++)
						ContainSite[Z + Y*off1] = 1;
					for (Z = 1; Z <= EN; Z++)
						ContainSite[Z + Y*off1] = 1;

				}
			}
		}
		for (AA = 0; AA <= 2; AA++) {
			Y = iseqs[AA];
			if (BMatch[X + Y * 3] > 0) {
#pragma omp parallel for private(Z, Diffs, Valid, ST, EN, A,s1, s2,th1, th2, th3)
				for (Z = 0; Z <= NextNo; Z++) {
					Diffs = 0.0;
					Valid = 0.0;
					ST = BPMatch[X + Y * 6];
					EN = BPMatch[X + 3 + Y * 6];
					if (BMatch[X + Z * 3] > 0) {

						if (ST < EN) {
							for (A = ST; A <= EN; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (ContainSite[A + Z*off1] == 1) {
											if (ContainSite[A + Y*off1] == 1) {
												if (s1 != s2)
													Diffs++;

												Valid++;
											}
										}
									}
								}
							}
						}
						else {
							for (A = ST; A <= LSeq; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (ContainSite[A + Z*off1] == 1) {
											if (ContainSite[A + Y*off1] == 1) {
												if (s1 != s2)
													Diffs++;

												Valid++;
											}
										}
									}
								}
							}
							for (A = 1; A <= EN; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (ContainSite[A + Z*off1] == 1) {
											if (ContainSite[A + Y*off1] == 1) {
												if (s1 != s2)
													Diffs++;

												Valid++;
											}
										}
									}
								}
							}

						}
					}

					else {
						if (ST < EN) {
							for (A = ST; A <= EN; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (s1 != s2)
											Diffs++;

										Valid++;

									}
								}
							}
						}
						else {
							for (A = ST; A <= LSeq; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (s1 != s2)
											Diffs++;

										Valid++;

									}
								}
							}
							for (A = 1; A <= EN; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (s1 != s2)
											Diffs++;

										Valid++;

									}
								}
							}

						}
					}
					if (Valid <= 30)
						MatchMat[X + AA * 3 + Z * 9] = 3;
					else if (Diffs / Valid >= 0.75)
						MatchMat[X + AA * 3 + Z * 9] = 3;
					else {
						th1 = 1 - (Diffs / Valid);
						th2 = (float)((4.0 * th1 - 1.0) / 3.0);
						th3 = (float)(log(th2));
						MatchMat[X + AA * 3 + Z * 9] = (float)(-0.75 * th3);
					}
					//MatchMat[X + Z*3 + Y*3*(NextNo+1)] = MatchMat[X + Y*3 + Z*3*(NextNo+1)];
				}
			}
			else {
				//this is smatsmall and not smat
				for (Z = 0; Z <= NextNo; Z++) {
					if (BMatch[X + Z * 3] == 0) {
						MatchMat[X + AA * 3 + Z * 9] = SMat[AA + Z * 3];
						//MatchMat[X + Z*3 + Y*9] = MatchMat[X + Y*3 + Z*9];
					}
				}
			}
		}

		return(1);

	}

	int MyMathFuncs::MakeMatchMatX2P2(int LSeq, int NextNo, int XX, char *ContainSite, float *SMat, float *MatchMat, float *BMatch, int *BPMatch, short int *SeqNum, int *iseqs) {
		//bpmatch 2,1,nexto
		//bmatch 2,nextno
		//containsite lseq,nextno
		//matchmat 2,nextno,nextno

		int X, AA, A, Y, Z, off1, ST, EN, s1, s2;
		float Diffs, Valid, th1, th2, th3;
		X = XX;
		off1 = LSeq + 1;
//#pragma omp parallel for private(Y, ST, EN, Z)
		for (Y = 0; Y <= NextNo; Y++) {

			if (BMatch[X + Y * 3] > 0) {
				ST = BPMatch[X + Y * 6];
				if (ST >LSeq)
					ST = LSeq;
				EN = BPMatch[X + 3 + Y * 6];
				if (ST < EN) {
					for (Z = ST; Z <= EN; Z++)
						ContainSite[Z + Y*off1] = 1;
				}
				else {
					for (Z = ST; Z <= LSeq; Z++)
						ContainSite[Z + Y*off1] = 1;
					for (Z = 1; Z <= EN; Z++)
						ContainSite[Z + Y*off1] = 1;

				}
			}
		}
		omp_set_num_threads(3);
#pragma omp parallel for private(AA, Y, Z, Diffs, Valid, ST, EN, A,s1, s2,th1, th2, th3)
		for (AA = 0; AA <= 2; AA++) {
			Y = iseqs[AA];
			if (BMatch[X + Y * 3] > 0) {

				for (Z = 0; Z <= NextNo; Z++) {
					Diffs = 0.0;
					Valid = 0.0;
					ST = BPMatch[X + Y * 6];
					EN = BPMatch[X + 3 + Y * 6];
					if (BMatch[X + Z * 3] > 0) {

						if (ST < EN) {
							for (A = ST; A <= EN; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (ContainSite[A + Z*off1] == 1) {
											if (ContainSite[A + Y*off1] == 1) {
												if (s1 != s2)
													Diffs++;

												Valid++;
											}
										}
									}
								}
							}
						}
						else {
							for (A = ST; A <= LSeq; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (ContainSite[A + Z*off1] == 1) {
											if (ContainSite[A + Y*off1] == 1) {
												if (s1 != s2)
													Diffs++;

												Valid++;
											}
										}
									}
								}
							}
							for (A = 1; A <= EN; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (ContainSite[A + Z*off1] == 1) {
											if (ContainSite[A + Y*off1] == 1) {
												if (s1 != s2)
													Diffs++;

												Valid++;
											}
										}
									}
								}
							}

						}
					}

					else {
						if (ST < EN) {
							for (A = ST; A <= EN; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (s1 != s2)
											Diffs++;

										Valid++;

									}
								}
							}
						}
						else {
							for (A = ST; A <= LSeq; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (s1 != s2)
											Diffs++;

										Valid++;

									}
								}
							}
							for (A = 1; A <= EN; A++) {
								s1 = SeqNum[A + Z*off1];
								if (s1 != 46) {
									s2 = SeqNum[A + Y*off1];
									if (s2 != 46) {

										if (s1 != s2)
											Diffs++;

										Valid++;

									}
								}
							}

						}
					}
					if (Valid <= 30)
						MatchMat[X + AA * 3 + Z * 9] = 3;
					else if (Diffs / Valid >= 0.75)
						MatchMat[X + AA * 3 + Z * 9] = 3;
					else {
						th1 = 1 - (Diffs / Valid);
						th2 = (float)((4.0 * th1 - 1.0) / 3.0);
						th3 = (float)(log(th2));
						MatchMat[X + AA * 3 + Z * 9] = (float)(-0.75 * th3);
					}
					//MatchMat[X + Z*3 + Y*3*(NextNo+1)] = MatchMat[X + Y*3 + Z*3*(NextNo+1)];
				}
			}
			else {
				//this is smatsmall and not smat
				for (Z = 0; Z <= NextNo; Z++) {
					if (BMatch[X + Z * 3] == 0) {
						MatchMat[X + AA * 3 + Z * 9] = SMat[AA + Z * 3];
						//MatchMat[X + Z*3 + Y*9] = MatchMat[X + Y*3 + Z*9];
					}
				}
			}
		}

		return(1);

	}

	int MyMathFuncs::MakeSDMP(int NextNo, int SLen, int *SP, int *EP, int *ISeqs, int *CompMat, unsigned char *MissingData, short int *SeqNum, double *SDM, double *DistMat) {
		int A, X, Y, Z, S0, S1, off0, off1, off2, offy;
		short int n0, n1, n2;
		double D0, D1, D2, V0, V1, V2, T;

//#pragma omp parallel for private(Seq2, se2, X, S1, S2)

		//compmat 2,1
		for (X = 0; X <3; X++) {
			S0 = ISeqs[CompMat[X]];
			S1 = ISeqs[CompMat[X + 3]];
			off2 = ISeqs[X] * (SLen + 1);
			off0 = S0*(SLen + 1);
			off1 = S1*(SLen + 1);

			for (Z = 0; Z < 5; Z++) {
				for (Y = 0; Y <= NextNo; Y++) {
					T = 0.0;
					D0 = 0.0;
					D1 = 0.0;
					D2 = 0.0;
					V0 = 0.0;
					V1 = 0.0;
					V2 = 0.0;
					A = SP[Z];
					offy = Y*(SLen + 1);
					while (A != EP[Z]) {
						if (MissingData[A + off0] == 0) {
							if (MissingData[A + off1] == 0) {
								if (MissingData[A + off2] == 0) {
									n0 = SeqNum[A + off0];
									if (n0 != 46) {
										n1 = SeqNum[A + off1];
										if (n1 != 46) {
											n2 = SeqNum[A + offy];
											if (n2 != 46) {

												if (n2 != n0 || n2 != n1) {
													if (n2 == n0) {
														D0 = D0 + 1;
														T = T + 1;
													}
													else if (n2 == n1) {
														D1 = D1 + 1;
														T = T + 1;
													}
													else if (n0 == n1) {
														D2 = D2 + 1;
														T = T + 1;
													}
												}
											}

										}

									}
									if (Z == 1 || Z == 2 || Z == 4) {
										if (SeqNum[A + off2] != 46) {
											if (SeqNum[A + off2] != SeqNum[A + offy])
												V1 = V1 + 1;

											V0 = V0 + 1;
										}
									}
								}
							}
						}
						A = A + 1;
						if (A > SLen)
							A = 1;

					}

					//sdm -2,2,nextno
					if (V0 > 0) {
						if (Z == 1)
							SDM[X + Y * 9] = V1 / V0;
						else if (Z == 2)
							SDM[X + 3 + Y * 9] = V1 / V0;
						else if (Z == 4)
							SDM[X + 6 + Y * 9] = V1 / V0;

					}
					else {
						if (Z == 1)
							SDM[X + Y * 9] = 10;
						else if (Z == 2)
							SDM[X + 3 + Y * 9] = 10;
						else if (Z == 4)
							SDM[X + 6 + Y * 9] = 10;

					}
					//distmat 2,4,nextno,2

					if (T > 0) {
						DistMat[X + Z * 3 + Y * 15] = D0 / T;

						DistMat[X + Z * 3 + Y * 15 + 15 * (NextNo + 1)] = D1 / T;

						DistMat[X + Z * 3 + Y * 15 + 30 * (NextNo + 1)] = D2 / T;
					}
					else {
						DistMat[X + Z * 3 + Y * 15] = 10;

						DistMat[X + Z * 3 + Y * 15 + 15 * (NextNo + 1)] = 10;

						DistMat[X + Z * 3 + Y * 15 + 30 * (NextNo + 1)] = 10;

					}

				}
			}
		}
		return (1);

	}
	int MyMathFuncs::MakeSDMP2(int NextNo, int SLen, int *SP, int *EP, int *ISeqs, int *CompMat, unsigned char *MissingData, short int *SeqNum, double *SDM, double *DistMat) {
		int A, X, Y, Z, S0, S1, off0, off1, off2, offy;
		short int n0, n1, n2;
		double D0, D1, D2, V0, V1, V2, T;
		omp_set_num_threads(3);
#pragma omp parallel for private(A, X, Y, Z, S0, S1, off0, off1, off2, offy, n0, n1, n2, D0, D1, D2, V0, V1, V2, T)
		for (X = 0; X <3; X++) {
			S0 = ISeqs[CompMat[X]];
			S1 = ISeqs[CompMat[X + 3]];
			off2 = ISeqs[X] * (SLen + 1);
			off0 = S0*(SLen + 1);
			off1 = S1*(SLen + 1);

			for (Z = 0; Z < 5; Z++) {
				for (Y = 0; Y <= NextNo; Y++) {
					T = 0.0;
					D0 = 0.0;
					D1 = 0.0;
					D2 = 0.0;
					V0 = 0.0;
					V1 = 0.0;
					V2 = 0.0;
					A = SP[Z];
					offy = Y*(SLen + 1);
					while (A != EP[Z]) {
						if (MissingData[A + off0] == 0) {
							if (MissingData[A + off1] == 0) {
								if (MissingData[A + off2] == 0) {
									n0 = SeqNum[A + off0];
									if (n0 != 46) {
										n1 = SeqNum[A + off1];
										if (n1 != 46) {
											n2 = SeqNum[A + offy];
											if (n2 != 46) {

												if (n2 != n0 || n2 != n1) {
													if (n2 == n0) {
														D0 = D0 + 1;
														T = T + 1;
													}
													else if (n2 == n1) {
														D1 = D1 + 1;
														T = T + 1;
													}
													else if (n0 == n1) {
														D2 = D2 + 1;
														T = T + 1;
													}
												}
											}

										}

									}
									if (Z == 1 || Z == 2 || Z == 4) {
										if (SeqNum[A + off2] != 46) {
											if (SeqNum[A + off2] != SeqNum[A + offy])
												V1 = V1 + 1;

											V0 = V0 + 1;
										}
									}
								}
							}
						}
						A = A + 1;
						if (A > SLen)
							A = 1;

					}

					//sdm -2,2,nextno
					if (V0 > 0) {
						if (Z == 1)
							SDM[X + Y * 9] = V1 / V0;
						else if (Z == 2)
							SDM[X + 3 + Y * 9] = V1 / V0;
						else if (Z == 4)
							SDM[X + 6 + Y * 9] = V1 / V0;

					}
					else {
						if (Z == 1)
							SDM[X + Y * 9] = 10;
						else if (Z == 2)
							SDM[X + 3 + Y * 9] = 10;
						else if (Z == 4)
							SDM[X + 6 + Y * 9] = 10;

					}
					//distmat 2,4,nextno,2

					if (T > 0) {
						DistMat[X + Z * 3 + Y * 15] = D0 / T;

						DistMat[X + Z * 3 + Y * 15 + 15 * (NextNo + 1)] = D1 / T;

						DistMat[X + Z * 3 + Y * 15 + 30 * (NextNo + 1)] = D2 / T;
					}
					else {
						DistMat[X + Z * 3 + Y * 15] = 10;

						DistMat[X + Z * 3 + Y * 15 + 15 * (NextNo + 1)] = 10;

						DistMat[X + Z * 3 + Y * 15 + 30 * (NextNo + 1)] = 10;

					}

				}
			}
		}
		omp_set_num_threads(2);
		return (1);

	}
	int MyMathFuncs::MakeLenFrag(int LenStrainSeq0, int NextNo,int ABPos, int AEPos, int *BCycle, int *BoundX, int UBSN21, short int *SeqNum2, int UBSN, short int *SeqNum) {
		int LenFrag, x, sn2, sn, Y;
		x = BoundX[0];
		LenFrag = 0;
		sn2 = UBSN21+1;
		sn = UBSN + 1;
		while (x >= 0){
			LenFrag++;

			for (Y = 0; Y <= NextNo; Y++)
				SeqNum2[LenFrag + Y*sn2] = SeqNum[x + Y*sn];
			
			x++;
			if (x == ABPos && BCycle[0] <= 0)
				break;
			if (x > LenStrainSeq0) {
				x = 1;
				BCycle[0] = BCycle[0] - 1;
				if (x == ABPos && BCycle[0] <= 0)
					break;
			}
		}
		while (x >= 0) {
			LenFrag++;
			for (Y = 0; Y <= NextNo; Y++)
				SeqNum2[LenFrag + Y*sn2] = SeqNum[x + Y*sn];
			x++;
			if (x == BoundX[1] && BCycle[1] <= 0)
				break;
			if (x > LenStrainSeq0) {
				x = 1;
				BCycle[1] = BCycle[1] - 1;
				if (x == BoundX[1] && BCycle[1] <= 0)
					break;
			}
		}
		x = BoundX[2];
		while (x >= 0) {
			LenFrag++;
			for (Y = 0; Y <= NextNo; Y++)
				SeqNum2[LenFrag + Y*sn2] = SeqNum[x + Y*sn];
			x++;
			if (x == AEPos && BCycle[2] == 0)
				break;
			if (x > LenStrainSeq0) {
				x = 1;
				BCycle[2] = BCycle[2] - 1;
				if (x == AEPos && BCycle[2] == 0)
					break;
			}
		}
		while (x >= 0) {
			LenFrag++;
			for (Y = 0; Y <= NextNo; Y++)
				SeqNum2[LenFrag + Y*sn2] = SeqNum[x + Y*sn];
			x++;
			if (x == BoundX[3] && BCycle[3] == 0)
				break;
			if (x > LenStrainSeq0) {
				x = 1;
				BCycle[3] = BCycle[3] - 1;
				if (x == BoundX[3] && BCycle[3] == 0)
					break;
			}
		}
		return (LenFrag);
	}

	int MyMathFuncs::MaketFSMat(int Nextno, int UBFM,int UBtFM, float *FMat, float *tFMat) {
		int X, Y, A, B, os, os2;
		os = UBFM + 1;
		os2 = UBtFM + 1;
		A = 0;
		B = 0;
		for (X = 0; X <= Nextno; X++) {

			if (FMat[X + X*os] != 3) {
				tFMat[A + A*os2] = 0;
				B = A + 1;

				for (Y = X + 1; Y <= Nextno; Y++) {
					if (FMat[Y + Y*os] != 3) {
						tFMat[A + B*os2] = FMat[X + Y*os];
						tFMat[B + A*os2] = FMat[X + Y*os];
						B++;
					}

				}
				A++;

			}
		}

		return(1);
	}
	int MyMathFuncs::MaketFSMatL(int Nextno, int UBFM, int UBtFM, float *FMat, float *tFMat, int *LR) {
		int X, Y, A, B, os, os2;
		os = UBFM + 1;
		os2 = UBtFM + 1;
		A = 0;
		B = 0;
		for (X = 0; X <= Nextno; X++) {

			if (LR[X] == 0) {
				tFMat[A + A*os2] = 0;
				B = A + 1;

				for (Y = X + 1; Y <= Nextno; Y++) {
					if (LR[Y]==0) {
						tFMat[A + B*os2] = FMat[X + Y*os];
						tFMat[B + A*os2] = FMat[X + Y*os];
						B++;
					}

				}
				A++;

			}
		}

		return(1);
	}
	int MyMathFuncs::MakeNodeDepthC(int Nextno, int PermNextno, int UBND1, int UBDM1, int UBDD, unsigned char *DoneDist, short int *NodeDepth, float *DMat, float *TraceBak) {
		int X, Y, A, B, os1, UB, os2, os3, hold, os4, os5, FirstA;
		
		os1 = UBND1 + 1;
		os2 = UBDD + 1;
		os3 = UBDM1 + 1;
		for (X = 0; X <= Nextno + 1; X++) {
			for (Y = 0; Y <= Nextno + 1; Y++) {
				NodeDepth[Y + X*os1] = -1;
			}
		}




		if (PermNextno > UBDM1)
			UB = UBDM1;
		else
			UB = PermNextno;

			
		for (X = 0; X <= UB; X++) {
			os4 = X*os2;
			os5 = X*os3;
			for (Y = X + 1; Y <= UB; Y++) {
				
				hold = (int)(DMat[Y + os5] * 1000);
				DoneDist[hold + os4] = 1;
				DoneDist[hold + Y*os2] = 1;

					
			}
		}
		A = -1;
		B = 0;

			
		for (X = 0; X <= PermNextno * 2; X++) {
			FirstA = 0;
			
			for (Y = 0; Y <= PermNextno; Y++) {
				if (DoneDist[X + Y*os2] == 1) {
					if (FirstA == 0) {
						FirstA = 1;
						A++;
						B = 0;
					}
					NodeDepth[A + B*os1] = Y;
					TraceBak[A] = ((float)(X)) / 1000;
					B++;
				}
			}
		}
		return(1);
	}
	int MyMathFuncs::FindSubSeqP8(int UBXO1, int UBXO2, int UBXO3, int *lenxoverseq, int en, unsigned char *goong, int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int *elementseq2, int *elementseq, short int spacerno, short int *seqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray)
	{

		int seq2, seq3, se2os, se3os, hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4, wmx, wpx, target, g2, xonos, xonos1;
		int s1, s2, s3, sz, ah0, ah1, ah2, xpdos, xpdos1;//,so1,so2,so3;
		int holder;
		int y;
		int g;

		const int se1 = seq1*lenseq;
		int se2;
		int se3;

		const int xow = xoverwindow;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		const int os1 = xoverwindow + lenseq + xoverwindow2;
		const int os2 = xoverwindow + lenseq2 + xoverwindow4;

		xpdos = lenseq + 200;
		xonos = (UBXO1 + 1)*(UBXO2 + 1);
		target = (UBXO1 + 1)*(UBXO2 + 1)*(en + 1);

		for (x = 0; x <= target; x++)
			xoverseqnumw[x] = 0;

		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel for private(g2, seq2, seq3, se2, se3,oc,hoc,ah0,ah1,ah2,y,g,xonos1,xpdos1,holder,x,se3os,se2os,b, wmx,wpx)
		for (g2 = 0; g2 <= en; g2++) {


			seq2 = elementseq2[g2];
			se2 = seq2*lenseq;
			seq3 = elementseq[g2];
			se3 = seq3*lenseq;

			if (goong[g2] == 1) {
				oc = 0;
				hoc = 0;
				ah0 = 0;
				ah1 = 0;
				ah2 = 0;
				y = 0;
				g = 0;
				xonos1 = g2*xonos;
				xpdos1 = g2*xpdos;
				holder = 0;
				if (spacerflag == 0) {
					//so1 = seq1*lenseq;
					//so2 = seq2*lenseq;
					//so3 = seq3*lenseq;


					y = 0;
					hoc = xoverwindow * 10;
					for (x = 1; x < lenseq; x++) {
						se2os = se2 + x;
						if (binarray[se2os] == 1) {//seq1 and seq2 are different
							se3os = se3 + x;
							if (binarray[se3os] == 1) {//if seq1 is also different to seq3{
													   //are seq2 and seq3 the same?
								if (seqnum[se2os] == seqnum[se3os]) {
									/*if (seqnum[x + se2] != 46) {
									if (seqnum[x + se3] != 46) {*/
									xdiffpos[++y + xpdos1] = x;
									xoverseqnumw[os2 + y + xonos1] = 1;
									ah2++;
									/*}
									}*/

								}
							}

							else if (binarray[se3os] == 0) {//seq1 and seq3 are the same

								xdiffpos[++y + xpdos1] = x;
								xoverseqnumw[os1 + y + xonos1] = 1;
								ah1++;

							}


						}
						else if (binarray[se3 + x] == 1) {//seq1 and seq3 are different but seq1=seq2
							if (binarray[se2os] == 0) {
								xdiffpos[++y + xpdos1] = x;
								xoverseqnumw[xow + y + xonos1] = 1;
								ah0++;
							}

						}
						xposdiff[x + xpdos1] = y;

					}
					//int hold;
					//for (x = 1; x < lenseq; x++) {
					//	if (binarray[x + se2]) {//seq1 and seq2 are different
					//		
					//		if (binarray[x + se3]) {//if seq1 is also different to seq3{
					//								are seq2 and seq3 the same?
					//			hold = (int)(seqnum[x + se2] == seqnum[x + se3]);
					//			y+= hold;
					//			xdiffpos[y] = x;							
					//			xoverseqnumw[y + os2]=hold;
					//			ah2+=hold;
					//				
					//			
					//		}

					//		else {//seq1 and seq3 are the same
					//			
					//			xdiffpos[y++] = x;
					//			xoverseqnumw[y + os1]=1;
					//			ah1++;
					//			
					//		}


					//	}
					//	else if (binarray[x + se3]) {//seq1 and seq3 are different but seq1=seq2
					//		xdiffpos[y++] = x;
					//		xoverseqnumw[y + xoverwindow]=1;
					//		ah0++;
					//		
					//	}
					//	xposdiff[x] = y;

					//}
					//


				}

				else if (spacerflag == 1) {
					for (x = 1; x < lenseq; x++) {
						*(xposdiff + x + xpdos1) = y;
						if (binarray[x + se2] || binarray[x + se3]) {

							s1 = *(seqnum + x + se1);
							s2 = *(seqnum + x + se2);
							s3 = *(seqnum + x + se3);

							//if (s1 != s2 || s1 != s3) {

							if (s1 == s2 || s1 == s3 || s2 == s3) {

								//if (s1 != 46) {
								//	if (s2 != 46) {
								//		if (s3 != 46) {
								if (binarray[x + se2] && binarray[x + se3]) {
									//If seq1 is odd one


									if (outlyer == seq1) {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}


										if (s2 == s3) {
											y++;

											xoverseqnumw[y + os2 + xonos1] = 1;
											ah2++;
											xdiffpos[y + xpdos1] = x;
											xposdiff[x + xpdos1] = y;
										}

									}
									else {
										for (z = 1; z <= spacerno; z++) {
											sz = *(seqnum + x + *(spacerseqs + z) * lenseq);
											if (sz == s1) {
												//If difference is legitimate
												if (*(xdiffpos + y + xpdos1) != x) {
													oc += 2;


													if (s2 == s3) {
														y++;

														xoverseqnumw[y + os2 + xonos1] = 1;
														ah2++;
														xdiffpos[y + xpdos1] = x;
														xposdiff[x + xpdos1] = y;
													}
													z = spacerno;

												}

											}

										}

									}
								}
								else if (s2 != s1  && s2 != s3) {
									//If seq2 is odd one

									if (outlyer == seq2) {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}

										if (s1 == s3) {
											y++;

											xoverseqnumw[y + os1 + xonos1] = 1;

											ah1++;
											xdiffpos[y + xpdos1] = x;
											xposdiff[x + xpdos1] = y;
										}


									}
									else {
										for (z = 1; z <= spacerno; z++) {
											sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
											if (s2 == sz) {
												//If difference is legitimate

												if (*(xdiffpos + y + xpdos1) != x) {
													oc += 2;

													if (s1 == s3) {
														y++;

														xoverseqnumw[y + os1 + xonos1] = 1;

														ah1++;
														xdiffpos[y + xpdos1] = x;
														xposdiff[x + xpdos1] = y;
													}

													z = spacerno;

												}
											}
										}
									}
								}
								else if (s3 != s1  && s3 != s2) {
									//If seq3 is odd one
									if (outlyer == seq3) {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}
										if (s1 == s2) {
											y++;
											xoverseqnumw[y + xoverwindow + xonos1] = 1;

											ah0++;
											xdiffpos[y + xpdos1] = x;
											xposdiff[x + xpdos1] = y;
										}



									}
									else {
										for (z = 1; z <= spacerno; z++) {
											sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
											if (s3 == sz) {

												//If difference is legitimate
												if (*(xdiffpos + y + xpdos1) != x) {
													oc += 2;
													if (s1 == s2) {
														y++;
														xoverseqnumw[y + xoverwindow + xonos1] = 1;

														ah0++;
														xdiffpos[y + xpdos1] = x;
														xposdiff[x + xpdos1] = y;
													}

													z = spacerno;
												}
											}
										}
									}

								}
								//}
								//}
								//}
							}
							//}
						}
					}

				}
				else if (spacerflag > 1) {


					for (x = 1; x < lenseq; x++) {
						xposdiff[x + xpdos1] = y;
						if (binarray[x + se2] || binarray[x + se3]) {
							s1 = seqnum[x + se1];
							s2 = seqnum[x + se2];
							s3 = seqnum[x + se3];


							//if (s1 != s2 || s1 != s3) {

							if (s1 != s2 && s1 != s3 && s2 != s3)
								g++;

							else {
								//if (s1 != 46) {
								//if (s2 != 46) {
								//if (s3 != 46) {

								if (s1 != s2  && s1 != s3) {
									//If seq1 is odd one
									for (z = 1; z <= spacerno; z++) {
										sz = seqnum[x + spacerseqs[z] * lenseq];
										if (sz == s1) {
											//If difference is legitimate
											if (seq1 != outlyer)
												oc += 2;
											else {
												if (oc > 0) {
													if (oc > hoc)
														hoc = oc;
													oc--;
												}
											}


											if (s2 == s3) {
												y++;

												xoverseqnumw[y + os2 + xonos1] = 1;
												ah2++;
												xdiffpos[y + xpdos1] = x;
												xposdiff[x + xpdos1] = y;
											}
											break;
										}
									}
								}

								else if (s2 != s1  && s2 != s3) {
									//If seq2 is odd one
									for (z = 1; z <= spacerno; z++) {
										sz = seqnum[x + spacerseqs[z] * lenseq];

										if (s2 == sz) {
											//If difference is legitimate



											if (seq2 != outlyer)
												oc += 2;
											else {
												if (oc > 0) {
													if (oc > hoc)
														hoc = oc;
													oc--;
												}
											}
											if (s1 == s3) {
												y++;

												xoverseqnumw[y + os1 + xonos1] = 1;

												ah1++;
												xdiffpos[y + xpdos1] = x;
												xposdiff[x + xpdos1] = y;
											}

											break;
										}
									}
								}
								else if (s3 != s1  && s3 != s2) {
									//If seq3 is odd one
									for (z = 1; z <= spacerno; z++) {
										sz = seqnum[x + spacerseqs[z] * lenseq];

										if (s3 == sz) {

											//If difference is legitimate

											if (seq3 != outlyer)
												oc += 2;
											else {
												if (oc > 0) {
													if (oc > hoc)
														hoc = oc;
													oc--;
												}
											}

											if (s1 == s2) {
												y++;
												xoverseqnumw[y + xoverwindow + xonos1] = 1;

												ah0++;
												xdiffpos[y + xpdos1] = x;
												xposdiff[x + xpdos1] = y;
											}

											break;
										}
									}
									//}
									//}
									//}

								}
							}
							//}
						}
					}

				}

				wmx = y - xow;
				wpx = y + xow;
				for (b = 1; b <= xow; b++) {
					wmx++;
					wpx++;
					xoverseqnumw[b + xonos1] = xoverseqnumw[wmx + xoverwindow + xonos1];
					xoverseqnumw[b + lenseq + xoverwindow2 + xonos1] = xoverseqnumw[wmx + os1 + xonos1];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
					xoverseqnumw[b + lenseq2 + xoverwindow4 + xonos1] = xoverseqnumw[wmx + os2 + xonos1];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

					xoverseqnumw[wpx + xonos1] = xoverseqnumw[b + xoverwindow + xonos1];//XOverSeqNum(X, 0)
					xoverseqnumw[wpx + lenseq + xoverwindow2 + xonos1] = xoverseqnumw[b + os1 + xonos1];
					xoverseqnumw[wpx + lenseq2 + xoverwindow4 + xonos1] = xoverseqnumw[b + os2 + xonos1];
				}


				//for (b = y+1; b < y + xoverwindow; b++)
				//	*(xdiffpos + b) = 0;


				//for (b = lenseq+1; b < lenseq + xoverwindow; b++)
				//	*(xposdiff + b) = 0;

				ah[0 + g2 * 3] = ah0;
				ah[1 + g2 * 3] = ah1;
				ah[2 + g2 * 3] = ah2;

				if (hoc < oc)
					hoc = oc;

				if (hoc < xoverwindow * 2)
					lenxoverseq[g2] = -y;
				else
					lenxoverseq[g2] = y;
			}
			else
				lenxoverseq[g2] = 0;

		}
		omp_set_num_threads(2);
		return(1);
	}

	int MyMathFuncs::CompressTE(int lseq, unsigned char *DecompressSeq, unsigned char *TEString, int *Decompress) {
		int A;
		for (A = 1; A <= lseq; A++)
			DecompressSeq[A-1] = TEString[Decompress[A]-1];
			
		return(1);
	}

	int MyMathFuncs::FakeMissing(int seq1,int lseq, int UBSN1, int UBTSN1, int *SeqnumBak, short int *Seqnum, short int *tSN) {
		int A, off, off2;
		off = seq1*(UBSN1 + 1);
		off2 = 3 * (UBTSN1 + 1);
		for (A = 0; A <= lseq; A++){
			if (SeqnumBak[A] < 50)
				Seqnum[A + off] = 46;
			tSN[A + off2] = Seqnum[A + off];
		}
		return(1);
	}


	int MyMathFuncs::MakeVarSites(int lseq, int BPos, int EPos, int sa, int SB, int SX, int SY, int UBSN1, short int *SeqNum, int *VXPos, short int *VarSiteMap, int *VSBE){
		int x, sao, sbo, sxo, syo, LenVarSeq;
		sao = sa*(UBSN1 + 1);
		sbo = SB*(UBSN1 + 1);
		sxo = SX*(UBSN1 + 1);
		syo = SY*(UBSN1 + 1);
		LenVarSeq = 0;
		for (x = 0; x<= lseq;x++){
            if (SeqNum[x + sao] != 46){
                if (SeqNum[x + sbo] != 46){
                    if (SeqNum[x + sxo] != 46){
                        if (SeqNum[x+ sao] != SeqNum[x + sbo] || SeqNum[x + sao] != SeqNum[x + sxo]){
                                
							LenVarSeq++;
                                
							VXPos[LenVarSeq] = x;
							if (SeqNum[x + sxo] == SeqNum[x + syo])
								VarSiteMap[LenVarSeq] = 2; 
                            else if (SeqNum[x + syo] != 46){
								if ((SeqNum[x + syo] == SeqNum[x + sbo] && SeqNum[x + sxo] == SeqNum[x + sao]) || (SeqNum[x + sxo] == SeqNum[x + sbo] && SeqNum[x + syo] == SeqNum[x + sao]))
									VarSiteMap[LenVarSeq] = -1;
								else if ((SeqNum[x + syo] == SeqNum[x+ sbo] && SeqNum[x+ syo] != SeqNum[x + sao]) || (SeqNum[x + syo] != SeqNum[x + sbo] && SeqNum[x + syo] == SeqNum[x+ sao]))
									VarSiteMap[LenVarSeq] = 0;
								else if (SeqNum[x+ syo] != SeqNum[x + sbo] && SeqNum[x + syo] != SeqNum[x + sao] && SeqNum[x + syo] != SeqNum[x + sao])
									VarSiteMap[LenVarSeq] = 1;
                                
                            }
                           
                        }
                    }
                }
            }
			if (x == BPos)
				VSBE[0] = LenVarSeq;
            
			if (x == EPos)
				VSBE[1] = LenVarSeq;
            
		}
			return(LenVarSeq);
	}
	double MyMathFuncs::MakeRCompatP(int *ISeqs, int *CompMat, int WinPP, int Nextno, int *RCompat, int *RCompatB, int *InPen, int *RCats, int *RNum, int *NRNum, int *GoodC, int *DoneX, int *Rlist, int *NRList, float *FAMat, double *LDist) {
		//goodc nextno,1
		int X, Y, Z, RL1, RL2, RL3, nCats, s0, s1, ds0, ds1;

		//get non-recombinant list
		s0 = ISeqs[CompMat[WinPP]];
		s1 = ISeqs[CompMat[WinPP + 3]];
		DoneX[s0] = 1;
		DoneX[s1] = 1;
		for (X = 0; X <= RNum[WinPP]; X++) {
			for (Y = 0; Y <= Nextno; Y++) {
				if (DoneX[Y] == 0) {
					if (GoodC[Y] == 1 || GoodC[Y + Nextno + 1] == 1) {

						if (FAMat[Rlist[WinPP + X * 3] + Y*(Nextno + 1)] < LDist[WinPP]) {
							//check and see if it is recombinant

							for (Z = 0; Z <= RNum[WinPP]; Z++) {

								if (Y == Rlist[WinPP + Z * 3])
									break;
							}

							if (Z == RNum[WinPP] + 1) {
								//ie it is non-recombinant
								DoneX[Y] = 1;
								NRList[WinPP + NRNum[WinPP] * 3] = Y;
								NRNum[WinPP] = NRNum[WinPP] + 1;
							}
						}
					}
				}
			}
		}

		NRNum[WinPP] = NRNum[WinPP] - 1;

		RCompat[WinPP] = 0;
		for (X = 0; X <= RNum[WinPP]; X++) {
			//ReDim RCats(Nextno * 3)
			for (Y = 0; Y <= Nextno * 3; Y++)
				RCats[Y] = 0;

			//get categories
			RL1 = Rlist[WinPP + X * 3];
			if (NRNum[WinPP] > -1) {
				for (Y = 0; Y <= NRNum[WinPP]; Y++) {
					RL2 = NRList[WinPP + Y * 3];
					RL3 = (long)((FAMat[RL1 + RL2*(Nextno + 1)] * 1000) + 0.0000001);
					//	RL4 = FAMat[RL1 + RL2*(Nextno+1)] * 1000;
					//	RL3 = long(RL4);
					RCats[RL3] = 1;
				}
			}

			//add the other iseqs
			if (FAMat[RL1 + s0*(Nextno + 1)] < LDist[WinPP]) {
				for (Y = 0; Y <= RNum[WinPP]; Y++) {
					RL2 = Rlist[WinPP + Y * 3];
					RL3 = (long)((FAMat[RL2 + s0*(Nextno + 1)] * 1000) + 0.0000001);
					RCats[RL3] = 1;
				}
			}
			if (FAMat[RL1 + s1*(Nextno + 1)] < LDist[WinPP]) {
				for (Y = 0; Y <= RNum[WinPP]; Y++) {
					RL2 = Rlist[WinPP + Y * 3];
					RL3 = (long)((FAMat[RL2 + s1*(Nextno + 1)] * 1000) + 0.0000001);
					RCats[RL3] = 1;
				}
			}
			//count the categories
			nCats = 0;
			for (Y = 0; Y <= Nextno; Y++)
				nCats = nCats + RCats[Y];

			if (nCats > RCompat[WinPP])
				RCompat[WinPP] = nCats;
		}
		ds0 = 0;
		ds1 = 0;
		for (X = 0; X <= RNum[WinPP]; X++) {
			//add iseqs to nrlist
			if (FAMat[Rlist[WinPP + X * 3] + s0*(Nextno + 1)] < LDist[WinPP] && ds0 == 0) {
				NRNum[WinPP] = NRNum[WinPP] + 1;
				NRList[WinPP + NRNum[WinPP] * 3] = s0;
				ds0 = 1;
			}
			if (FAMat[Rlist[WinPP + X * 3] + s1*(Nextno + 1)] < LDist[WinPP] && ds1 == 0) {
				NRNum[WinPP] = NRNum[WinPP] + 1;
				NRList[WinPP + NRNum[WinPP] * 3] = s1;
				ds1 = 1;
			}
		}
		RCompatB[WinPP] = 0;
		if (NRNum[WinPP] > -1) {
			for (X = 0; X <= NRNum[WinPP]; X++) {
				for (Y = 0; Y <= Nextno * 3; Y++)
					RCats[Y] = 0;
				//get categories
				RL1 = NRList[WinPP + X * 3];

				for (Y = 0; Y <= RNum[WinPP]; Y++) {
					RL2 = Rlist[WinPP + Y * 3];
					RL3 = (long)((FAMat[RL1 + RL2*(Nextno + 1)] * 1000) + 0.0000001);
					//RL4 = FAMat[RL1 + RL2*(Nextno+1)] * 1000;
					//RL3 = long(RL4);
					RCats[RL3] = 1;
				}

				/*//add the other iseqs
				if (FAMat[RL1 + s0*(Nextno+1)] < LDist[WinPP]){
				for (Y = 0; Y <= RNum[WinPP]; Y++){
				RL2 = Rlist[WinPP + Y*3];
				RL3 = (long)((FAMat[RL2 + s0*(Nextno+1)] * 1000)+0.0000001);
				RCats[RL3] = 1;
				}
				}
				if (FAMat[RL1 + s1*(Nextno+1)] < LDist[WinPP]){
				for (Y = 0; Y <= RNum[WinPP]; Y++){
				RL2 = Rlist[WinPP + Y*3];
				RL3 = (long)((FAMat[RL2 + s1*(Nextno+1)] * 1000)+0.0000001);
				RCats[RL3] = 1;
				}
				}*/
				//count them
				nCats = 0;
				for (Y = 0; Y <= Nextno; Y++)
					nCats = nCats + RCats[Y];

				nCats = nCats - 1;

				if (nCats > RCompatB[WinPP])
					RCompatB[WinPP] = nCats;
			}
		}
		if (NRNum[WinPP] > -1) {
			if (RCompatB[WinPP] < RCompat[WinPP])
				RCompat[WinPP] = RCompatB[WinPP];
		}

		if (RCompat[WinPP] > RNum[WinPP])
			RCompat[WinPP] = RNum[WinPP];

		if (RCompat[WinPP] > 0)
			RCompat[WinPP] = RCompat[WinPP] + InPen[WinPP];// 'penalise the inversions




		return(1);
	}


	int MyMathFuncs::FindSubSeqP7(int UBXO1, int UBXO2, int UBXO3, int *lenxoverseq, int en, unsigned char *goong,int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int *elementseq, short int spacerno, short int *seqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray)
	{

		int seq3, se2os, se3os, hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4, wmx, wpx, target,  g2, xonos, xonos1;
		int s1, s2, s3, sz, ah0, ah1, ah2, xpdos, xpdos1;//,so1,so2,so3;
		int holder;
		int y;
		int g;

		const int se1 = seq1*lenseq;
		const int se2 = seq2*lenseq;
		int se3; 
		
		const int xow = xoverwindow;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		const int os1 = xoverwindow + lenseq + xoverwindow2;
		const int os2 = xoverwindow + lenseq2 + xoverwindow4;

		xpdos = lenseq + 200;
		xonos = (UBXO1+1)*(UBXO2+1);
		target = (UBXO1 + 1)*(UBXO2 + 1)*(en+1);
		
		for (x = 0; x <= target; x++)
			xoverseqnumw[x] = 0;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private(g2, seq3,se3,oc,hoc,ah0,ah1,ah2,y,g,xonos1,xpdos1,holder,x,se3os,se2os,b, wmx,wpx, s1, s2, s3,sz, z)
		for (g2 = 0; g2 <= en; g2++) {
			
			seq3 = elementseq[g2];
			se3 = seq3*lenseq;
			
			if (goong[g2] == 1) {
				oc = 0;
				hoc = 0;
				ah0 = 0;
				ah1 = 0;
				ah2 = 0;
				y = 0;
				g = 0;
				xonos1 = g2*xonos;
				xpdos1 = g2*xpdos;
				holder = 0;
				if (spacerflag == 0) {
					//so1 = seq1*lenseq;
					//so2 = seq2*lenseq;
					//so3 = seq3*lenseq;


					y = 0;
					hoc = xoverwindow * 10;
					for (x = 1; x < lenseq; x++) {
						se2os = se2 + x;
						if (binarray[se2os] == 1) {//seq1 and seq2 are different
							se3os = se3 + x;
							if (binarray[se3os] == 1) {//if seq1 is also different to seq3{
													   //are seq2 and seq3 the same?
								if (seqnum[se2os] == seqnum[se3os]) {
									/*if (seqnum[x + se2] != 46) {
									if (seqnum[x + se3] != 46) {*/
									xdiffpos[++y + xpdos1] = x;
									xoverseqnumw[os2 + y + xonos1] = 1;
									ah2++;
									/*}
									}*/

								}
							}

							else if (binarray[se3os] == 0) {//seq1 and seq3 are the same

								xdiffpos[++y + xpdos1] = x;
								xoverseqnumw[os1 + y + xonos1] = 1;
								ah1++;

							}


						}
						else if (binarray[se3 + x] == 1) {//seq1 and seq3 are different but seq1=seq2
							if (binarray[se2os] == 0) {
								xdiffpos[++y + xpdos1] = x;
								xoverseqnumw[xow + y + xonos1] = 1;
								ah0++;
							}

						}
						xposdiff[x + xpdos1] = y;

					}
					//int hold;
					//for (x = 1; x < lenseq; x++) {
					//	if (binarray[x + se2]) {//seq1 and seq2 are different
					//		
					//		if (binarray[x + se3]) {//if seq1 is also different to seq3{
					//								are seq2 and seq3 the same?
					//			hold = (int)(seqnum[x + se2] == seqnum[x + se3]);
					//			y+= hold;
					//			xdiffpos[y] = x;							
					//			xoverseqnumw[y + os2]=hold;
					//			ah2+=hold;
					//				
					//			
					//		}

					//		else {//seq1 and seq3 are the same
					//			
					//			xdiffpos[y++] = x;
					//			xoverseqnumw[y + os1]=1;
					//			ah1++;
					//			
					//		}


					//	}
					//	else if (binarray[x + se3]) {//seq1 and seq3 are different but seq1=seq2
					//		xdiffpos[y++] = x;
					//		xoverseqnumw[y + xoverwindow]=1;
					//		ah0++;
					//		
					//	}
					//	xposdiff[x] = y;

					//}
					//


				}

				else if (spacerflag == 1) {
					for (x = 1; x < lenseq; x++) {
						*(xposdiff + x + xpdos1) = y;
						if (binarray[x + se2] || binarray[x + se3]) {

							s1 = *(seqnum + x + se1);
							s2 = *(seqnum + x + se2);
							s3 = *(seqnum + x + se3);

							//if (s1 != s2 || s1 != s3) {

							if (s1 == s2 || s1 == s3 || s2 == s3) {

								//if (s1 != 46) {
								//	if (s2 != 46) {
								//		if (s3 != 46) {
								if (binarray[x + se2] && binarray[x + se3]) {
									//If seq1 is odd one


									if (outlyer == seq1) {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}


										if (s2 == s3) {
											y++;

											xoverseqnumw[y + os2 + xonos1] = 1;
											ah2++;
											xdiffpos[y + xpdos1] = x;
											xposdiff[x + xpdos1] = y;
										}

									}
									else {
										for (z = 1; z <= spacerno; z++) {
											sz = *(seqnum + x + *(spacerseqs + z) * lenseq);
											if (sz == s1) {
												//If difference is legitimate
												if (*(xdiffpos + y + xpdos1) != x) {
													oc += 2;


													if (s2 == s3) {
														y++;

														xoverseqnumw[y + os2 + xonos1] = 1;
														ah2++;
														xdiffpos[y + xpdos1] = x;
														xposdiff[x + xpdos1] = y;
													}
													z = spacerno;

												}

											}

										}

									}
								}
								else if (s2 != s1  && s2 != s3) {
									//If seq2 is odd one

									if (outlyer == seq2) {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}

										if (s1 == s3) {
											y++;

											xoverseqnumw[y + os1 + xonos1] = 1;

											ah1++;
											xdiffpos[y + xpdos1] = x;
											xposdiff[x + xpdos1] = y;
										}


									}
									else {
										for (z = 1; z <= spacerno; z++) {
											sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
											if (s2 == sz) {
												//If difference is legitimate

												if (*(xdiffpos + y + xpdos1) != x) {
													oc += 2;

													if (s1 == s3) {
														y++;

														xoverseqnumw[y + os1 + xonos1] = 1;

														ah1++;
														xdiffpos[y + xpdos1] = x;
														xposdiff[x + xpdos1] = y;
													}

													z = spacerno;

												}
											}
										}
									}
								}
								else if (s3 != s1  && s3 != s2) {
									//If seq3 is odd one
									if (outlyer == seq3) {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}
										if (s1 == s2) {
											y++;
											xoverseqnumw[y + xoverwindow + xonos1] = 1;

											ah0++;
											xdiffpos[y + xpdos1] = x;
											xposdiff[x + xpdos1] = y;
										}



									}
									else {
										for (z = 1; z <= spacerno; z++) {
											sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
											if (s3 == sz) {

												//If difference is legitimate
												if (*(xdiffpos + y + xpdos1) != x) {
													oc += 2;
													if (s1 == s2) {
														y++;
														xoverseqnumw[y + xoverwindow + xonos1] = 1;

														ah0++;
														xdiffpos[y + xpdos1] = x;
														xposdiff[x + xpdos1] = y;
													}

													z = spacerno;
												}
											}
										}
									}

								}
								//}
								//}
								//}
							}
							//}
						}
					}

				}
				else if (spacerflag > 1) {


					for (x = 1; x < lenseq; x++) {
						xposdiff[x + xpdos1] = y;
						if (binarray[x + se2] || binarray[x + se3]) {
							s1 = seqnum[x + se1];
							s2 = seqnum[x + se2];
							s3 = seqnum[x + se3];


							//if (s1 != s2 || s1 != s3) {

							if (s1 != s2 && s1 != s3 && s2 != s3)
								g++;

							else {
								//if (s1 != 46) {
								//if (s2 != 46) {
								//if (s3 != 46) {

								if (s1 != s2  && s1 != s3) {
									//If seq1 is odd one
									for (z = 1; z <= spacerno; z++) {
										sz = seqnum[x + spacerseqs[z] * lenseq];
										if (sz == s1) {
											//If difference is legitimate
											if (seq1 != outlyer)
												oc += 2;
											else {
												if (oc > 0) {
													if (oc > hoc)
														hoc = oc;
													oc--;
												}
											}


											if (s2 == s3) {
												y++;

												xoverseqnumw[y + os2 + xonos1] = 1;
												ah2++;
												xdiffpos[y + xpdos1] = x;
												xposdiff[x + xpdos1] = y;
											}
											break;
										}
									}
								}

								else if (s2 != s1  && s2 != s3) {
									//If seq2 is odd one
									for (z = 1; z <= spacerno; z++) {
										sz = seqnum[x + spacerseqs[z] * lenseq];

										if (s2 == sz) {
											//If difference is legitimate



											if (seq2 != outlyer)
												oc += 2;
											else {
												if (oc > 0) {
													if (oc > hoc)
														hoc = oc;
													oc--;
												}
											}
											if (s1 == s3) {
												y++;

												xoverseqnumw[y + os1 + xonos1] = 1;

												ah1++;
												xdiffpos[y + xpdos1] = x;
												xposdiff[x + xpdos1] = y;
											}

											break;
										}
									}
								}
								else if (s3 != s1  && s3 != s2) {
									//If seq3 is odd one
									for (z = 1; z <= spacerno; z++) {
										sz = seqnum[x + spacerseqs[z] * lenseq];

										if (s3 == sz) {

											//If difference is legitimate

											if (seq3 != outlyer)
												oc += 2;
											else {
												if (oc > 0) {
													if (oc > hoc)
														hoc = oc;
													oc--;
												}
											}

											if (s1 == s2) {
												y++;
												xoverseqnumw[y + xoverwindow + xonos1] = 1;

												ah0++;
												xdiffpos[y + xpdos1] = x;
												xposdiff[x + xpdos1] = y;
											}

											break;
										}
									}
									//}
									//}
									//}

								}
							}
							//}
						}
					}

				}

				wmx = y - xow;
				wpx = y + xow;
				for (b = 1; b <= xow; b++) {
					wmx++;
					wpx++;
					xoverseqnumw[b + xonos1] = xoverseqnumw[wmx + xoverwindow + xonos1];
					xoverseqnumw[b + lenseq + xoverwindow2 + xonos1] = xoverseqnumw[wmx + os1 + xonos1];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
					xoverseqnumw[b + lenseq2 + xoverwindow4 + xonos1] = xoverseqnumw[wmx + os2 + xonos1];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

					xoverseqnumw[wpx + xonos1] = xoverseqnumw[b + xoverwindow + xonos1];//XOverSeqNum(X, 0)
					xoverseqnumw[wpx + lenseq + xoverwindow2 + xonos1] = xoverseqnumw[b + os1 + xonos1];
					xoverseqnumw[wpx + lenseq2 + xoverwindow4 + xonos1] = xoverseqnumw[b + os2 + xonos1];
				}


				//for (b = y+1; b < y + xoverwindow; b++)
				//	*(xdiffpos + b) = 0;


				//for (b = lenseq+1; b < lenseq + xoverwindow; b++)
				//	*(xposdiff + b) = 0;

				ah[0 + g2*3] = ah0;
				ah[1 + g2*3] = ah1;
				ah[2 + g2*3] = ah2;

				if (hoc < oc)
					hoc = oc;

				if (hoc < xoverwindow * 2)
					lenxoverseq[g2] = -y;
				else
					lenxoverseq[g2] = y;
			}
			else
				lenxoverseq[g2] = 0;

		}
		omp_set_num_threads(2);
		return(1);
	}
	

	


	int MyMathFuncs::MakeISeq4P(int Nextno, int UBNS, int UBIS4, short int *SeqCompressor4, short int *ISeq4, char *NumSeq) {
		//SeqCompressor4(4,4,4,4), ISeq4(UBIS4,Nextno)
		int X, Y, A, B, C, D, StepPos, os1, os2, os3, os4;
		os1 = UBNS + 1;
		os3 = UBIS4 + 1;
		for (X = 0; X <= Nextno; X++) {
			StepPos = 0;
			os2 = os1*X;
			os4 = os3*X;
			for (Y = 1; Y <= UBNS - 4; Y += 4) {
				StepPos++;
				A = (int)(NumSeq[Y + os2]);
				B = (int)(NumSeq[Y + 1 + os2]);
				C = (int)(NumSeq[Y + 2 + os2]);
				D = (int)(NumSeq[Y + 3 + os2]);
				ISeq4[StepPos + os4] = SeqCompressor4[A + B * 5 + C * 25 + D * 125];
			}
			StepPos++;
			A = 0;
			B = 0;
			C = 0;
			D = 0;
			if (Y <= UBNS) {
				A = (int)(NumSeq[Y + os2]);
				if (Y + 1 <= UBNS) {
					B = (int)(NumSeq[Y + 1 + os2]);
					if (Y + 2 <= UBNS) {
						C = (int)(NumSeq[Y + 2 + os2]);
						if (Y + 3 <= UBNS)
							D = (int)(NumSeq[Y + 3 + os2]);

					}
				}
			}
			ISeq4[StepPos + os4] = SeqCompressor4[A + B * 5 + C * 25 + D * 125];
		}
		return(1);
	}

	int MyMathFuncs::MakeISeq3P(int Nextno, int UBNS, int UBIS4, short int *SeqCompressor4, short int *ISeq4, char *NumSeq) {
		//SeqCompressor4(4,4,4), ISeq4(UBIS4,Nextno)
		int X, Y, A, B, C, StepPos, os1, os2, os3, os4;
		os1 = UBNS + 1;
		os3 = UBIS4 + 1;
		for (X = 0; X <= Nextno; X++) {
			StepPos = 0;
			os2 = os1*X;
			os4 = os3*X;
			for (Y = 1; Y <= UBNS - 3; Y += 3) {
				StepPos++;
				A = (int)(NumSeq[Y + os2]);
				B = (int)(NumSeq[Y + 1 + os2]);
				C = (int)(NumSeq[Y + 2 + os2]);
				
				ISeq4[StepPos + os4] = SeqCompressor4[A + B * 5 + C * 25];
			}
			StepPos++;
			A = 0;
			B = 0;
			C = 0;
			
			if (Y <= UBNS) {
				A = (int)(NumSeq[Y + os2]);
				if (Y + 1 <= UBNS) {
					B = (int)(NumSeq[Y + 1 + os2]);
					if (Y + 2 <= UBNS) {
						C = (int)(NumSeq[Y + 2 + os2]);
						
					}
				}
			}
			ISeq4[StepPos + os4] = SeqCompressor4[A + B * 5 + C * 25];
		}
		return(1);
	}

	int MyMathFuncs::MakeNumSeqP(int Nextno, int SLen, int UBNS, int StartPosInAlign, int EndPosInAlign, unsigned char *ConvNumSeq, short int *SeqNum, unsigned char *NumSeq) {
		int X, Y, Offset, os1, os2, os3, os4, os5, os6;

		os1 = SLen + 1;
		os5 = UBNS + 1;


		ConvNumSeq[66] = 1;
		ConvNumSeq[68] = 2;
		ConvNumSeq[72] = 3;
		ConvNumSeq[85] = 4;

		if (StartPosInAlign < EndPosInAlign) {
			for (X = 0; X <= Nextno; X++) {
				os4 = os1*X;
				os6 = os5*X;
				for (Y = StartPosInAlign; Y <= EndPosInAlign; Y++) {
					os2 = Y + os4;
					os3 = Y + os6 - StartPosInAlign;
					NumSeq[os3] = ConvNumSeq[SeqNum[os2]];//for some reason without this offset the wrong number is read from seqnum
				}
			}
		}
		else {
			Offset = SLen - StartPosInAlign + 1;
			for (X = 0; X <= Nextno; X++) {
				os4 = os1*X;
				os6 = os5*X;
				for (Y = StartPosInAlign; Y <= SLen; Y++) {
					os2 = Y + os4;
					os3 = Y + os6 - StartPosInAlign ;
					NumSeq[os3] = ConvNumSeq[SeqNum[os2]];

				}

				for (Y = 1; Y <= EndPosInAlign; Y++) {
					os2 = Y + os4;
					os3 = Y + os6 + Offset;
					NumSeq[os3] = ConvNumSeq[SeqNum[os2]];

				}

			}
		}

		return(1);
	}

	
	int MyMathFuncs::MakeBinArrayP4(int Seq1, int  Nextno, int UBIS4, int UBBC, int UBBA, short int *Maskseq, short int *ISeq4, unsigned char *BinArray, unsigned char *BinConverter4) {
		int X, Y, os1,s1os, os2, os3, os4, os5;
		os1 = UBIS4 + 1;
		os2 = UBBC + 1;
		os4 = UBBA + 1;
		s1os = os1*Seq1;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private(X, Y, os5, os3)
		for (X = Seq1 + 1; X <= Nextno; X++) {
			
			if (Maskseq[X] == 0) {
				os5 = os4*X;
				os3 = X*os1;
				for (Y = 1; Y <= UBIS4; Y++)
					BinArray[Y + os5] = BinConverter4[ISeq4[Y + s1os] + ISeq4[Y + os3]*os2];

			}
		}
		omp_set_num_threads(2);
		return(1);
	}


	int MyMathFuncs::MarkRemovalsP(int Nextno, int WinPP, int Redolistsize, int *RedoList, int *RNum, int *Rlist, unsigned char *DoPairs) {
		int X, Z, Y;
		for (X = 0; X <= Redolistsize; X++) {
			for (Y = 1; Y < 4; Y++) {
				for (Z = 0; Z <= RNum[WinPP]; Z++) {
					if (Rlist[WinPP + Z * 3] == RedoList[Y + X * 4]) {
						if (Y == 1) {
							if (DoPairs[RedoList[2 + X * 4] + RedoList[3 + X * 4] * (Nextno + 1)] == 1)
								break;
						}
						else if (Y == 2) {
							if (DoPairs[RedoList[1 + X * 4] + RedoList[3 + X * 4] * (Nextno + 1)] == 1)
								break;
						}
						else if (Y == 3) {
							if (DoPairs[RedoList[1 + X * 4] + RedoList[2 + X * 4] * (Nextno + 1)] == 1)
								break;
						}

					}
				}
				if (Z <= RNum[WinPP]) {
					RedoList[X * 4] = -1;
					break;
				}
			}
		}

		return(1);
	}

	int MyMathFuncs::MakePairsP(int Nextno, int Da, int Ma, int Mi, int WinPP, int *RNum, int *Rlist, unsigned char *DoPairs) {
		int WinPPY, off1;
		for (WinPPY = 0; WinPPY <= RNum[WinPP]; WinPPY++) {
			off1 = Rlist[WinPP + WinPPY * 3];
			if (off1 == Da || off1 == Ma || off1 == Mi) {

				DoPairs[Mi + Ma*(Nextno + 1)] = 1;
				DoPairs[Ma + Mi*(Nextno + 1)] = 1;

				DoPairs[Da + Ma*(Nextno + 1)] = 1;
				DoPairs[Ma + Da*(Nextno + 1)] = 1;

				DoPairs[Mi + Da*(Nextno + 1)] = 1;
				DoPairs[Da + Mi*(Nextno + 1)] = 1;

				break;
			}

		}
		return(WinPPY);
	}

	int MyMathFuncs::MakeBinArray2P(int UBPV1, float *permvalid, int UBDP1, unsigned char *dopairs, int Seq1, int LSeq, int Nextno, short int *Maskseq, short int *SeqNum, unsigned char *BinArray, int *slookupnum, int *slookup, unsigned char *isin, int *tracesub, int *actualsize, int MinSeqSize) {
		int Seq2, X, os, se1, se2, target, S1, S2, g, os2, os3;
		os = LSeq + 1;
		se1 = os*Seq1;
		target = LSeq + Nextno*os;
		os2 = UBPV1 + 1;
		os3 = UBDP1 + 1;
		for (Seq2 = 0; Seq2 <= target; Seq2++)
			BinArray[Seq2] = 0;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private(g, Seq2, se2, X, S1, S2)
		for (g = 1; g <= slookupnum[0]; g++) {
			//for (Seq2 = Seq1 + 1; Seq2 <= Nextno; Seq2++){
			Seq2 = slookup[g * 2];
			if (actualsize[Seq2] > MinSeqSize) {
				if (permvalid[Seq1 + Seq2*os2] > MinSeqSize) {
					if (isin[Seq2] == 0) {
						if (dopairs[Seq1 + Seq2*os3] == 1){
							if (tracesub[Seq1] != tracesub[Seq2]) {
								//if (Maskseq[Seq2] == 0){
								se2 = os*Seq2;
								for (X = 1; X <= LSeq; X++) {
									S1 = SeqNum[X + se1];
									S2 = SeqNum[X + se2];
									if (S1 != S2) {

										if (S1 != 46) {
											if (S2 != 46)
												BinArray[X + se2] = 1;
											else
												BinArray[X + se2] = 3;
										}
										else
											BinArray[X + se2] = 2;

									}
								}
							}
						}
					}
				}
			}
		}
		omp_set_num_threads(2);
		return(1);
	}

	int MyMathFuncs::MakeBinArray3P(int SNextNo, int UBDP1, unsigned char *dopairs, int Seq1, int LSeq, int Nextno, short int *Maskseq, short int *SeqNum, unsigned char *BinArray, int *slookupnum, int *slookup, int *tracesub, int *actualsize, int MinSeqSize) {
		int Seq2, X, os, se1, se2, target, S1, S2, g, os2, os3, tSeq2;
		os = LSeq + 1;
		se1 = os*Seq1;
		target = LSeq + Nextno*os;
		
		os3 = UBDP1 + 1;
		for (Seq2 = 0; Seq2 <= target; Seq2++)
			BinArray[Seq2] = 0;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private (g, Seq2, tSeq2, se2, S1, S2)
		for (g = 1; g <= slookupnum[0]; g++) {
			//for (Seq2 = Seq1 + 1; Seq2 <= Nextno; Seq2++){
			Seq2 = slookup[g * 2];
			if (actualsize[Seq2] > MinSeqSize) {
				//if (permvalid[Seq1 + Seq2*os2] > MinSeqSize) {
				if (Seq2 > SNextNo)
					tSeq2 = tracesub[Seq2];
				else
					tSeq2 = Seq2;
				
				//If DoPairs(TraceSub(Seq1), tSeq2) = 1 Then
				//if (isin[Seq2] == 0) {
						if (dopairs[tracesub[Seq1] + tSeq2*os3] == 1) {
							if (tracesub[Seq1] != tSeq2) {
								//if (Maskseq[Seq2] == 0){
								se2 = os*Seq2;
								for (X = 1; X <= LSeq; X++) {
									S1 = SeqNum[X + se1];
									S2 = SeqNum[X + se2];
									if (S1 != S2) {

										if (S1 != 46) {
											if (S2 != 46)
												BinArray[X + se2] = 1;
											else
												BinArray[X + se2] = 3;
										}
										else
											BinArray[X + se2] = 2;

									}
								}
							}
						}
					//}
				//}
			}
		}
		omp_set_num_threads(2);
		return(1);
	}
	
		
	int MyMathFuncs::FindSubSeqP5(int UBXSN, int UBBA,  int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int seq3, short int spacerno, short int *seqnum, short int *xoverseqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray3, unsigned char *binarray4, unsigned char *PairwiseTripArray, unsigned char *SS255RDP)
	{

		int hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4, wmx, wpx, target, D, A, B, ptaos1, ptaos2, osba, os5,os6;
		int s1, s2, s3, s2ba,s3ba,sz, ah0, ah1, ah2;//,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int g = 0;

		const int se1 = seq1*lenseq;
		const int se2 = seq2*lenseq;
		const int se3 = seq3*lenseq;
		oc = 0;
		hoc = 0;
		const int xow = xoverwindow;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		const int os1 = xoverwindow + lenseq + xoverwindow2;
		const int os2 = xoverwindow + lenseq2 + xoverwindow4;
		os5 = UBXSN + 1;
		ah0 = 0;
		ah1 = 0;
		ah2 = 0;
		ptaos1 = 256;
		ptaos2 = 256 * 256;
		target = (lenseq + xoverwindow2 + 1) * 3;
		osba = UBBA + 1;
		//s1ba = seq1*osba;
		s2ba = seq2*osba;
		s3ba = seq3*osba;
		for (x = 1; x < target; x++)
			xoverseqnumw[x] = 0;

		if (spacerflag == 0) {
			//so1 = seq1*lenseq;
			//so2 = seq2*lenseq;
			//so3 = seq3*lenseq;

			y = 0;
			x = -1;

			for (D = 1; D <= UBBA; D++){
				os6 = D + s3ba;
				A = PairwiseTripArray[binarray3[D + s2ba] + binarray3[os6] * ptaos1 + binarray4[os6] * ptaos2];
				if (A < 255) {
					B = SS255RDP[A];
					x++;
					if (B != 3) {
						xdiffpos[++y] = x;
						ah[B]++;
						xoverseqnumw[xow + y + B*os5] = 1;
					}
					xposdiff[x] = y;
					B = SS255RDP[A + 256];
					x++;
					if (B != 3) {

						xdiffpos[++y] = x;
						ah[B]++;
						xoverseqnumw[xow + y + B*os5] = 1;
					}
					xposdiff[x] = y;
					B = SS255RDP[A + 512];
					x++;
					if (B != 3) {
						xdiffpos[++y] = x;
						ah[B]++;
						xoverseqnumw[xow + y + B*os5] = 1;
					}
					xposdiff[x] = y;
					B = SS255RDP[A + 768];
					x++;
					if (B != 3){								
						xdiffpos[++y] = x;
						ah[B]++;
						xoverseqnumw[xow + y + B*os5] = 1;
					}
					xposdiff[x] = y;
				}
				else {

					xposdiff[x + 1] = y;
					xposdiff[x + 2] = y;
					xposdiff[x + 3] = y;
					xposdiff[x + 4] = y;
					x = x + 4;
				}

			}
				//LenXOverSeq = Y 
			ah0 = ah[0];
			ah1 = ah[1];
			ah2 = ah[2];

			hoc = xoverwindow * 10;
			


		}

		else if (spacerflag == 1) {
			for (x = 1; x < lenseq; x++) {
				*(xposdiff + x) = y;
				if (binarray3[x + se2] || binarray3[x + se3]) {

					s1 = *(seqnum + x + se1);
					s2 = *(seqnum + x + se2);
					s3 = *(seqnum + x + se3);

					//if (s1 != s2 || s1 != s3) {

					if (s1 == s2 || s1 == s3 || s2 == s3) {

						//if (s1 != 46) {
						//	if (s2 != 46) {
						//		if (s3 != 46) {
						if (binarray3[x + se2] && binarray3[x + se3]) {
							//If seq1 is odd one


							if (outlyer == seq1) {
								if (oc > 0) {
									if (oc > hoc)
										hoc = oc;
									oc--;
								}


								if (s2 == s3) {
									y++;

									xoverseqnumw[y + os2] = 1;
									ah2++;
									xdiffpos[y] = x;
									xposdiff[x] = y;
								}

							}
							else {
								for (z = 1; z <= spacerno; z++) {
									sz = *(seqnum + x + *(spacerseqs + z) * lenseq);
									if (sz == s1) {
										//If difference is legitimate
										if (*(xdiffpos + y) != x) {
											oc += 2;


											if (s2 == s3) {
												y++;

												xoverseqnumw[y + os2] = 1;
												ah2++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}
											z = spacerno;

										}

									}

								}

							}
						}
						else if (s2 != s1  && s2 != s3) {
							//If seq2 is odd one

							if (outlyer == seq2) {
								if (oc > 0) {
									if (oc > hoc)
										hoc = oc;
									oc--;
								}

								if (s1 == s3) {
									y++;

									xoverseqnumw[y + os1] = 1;

									ah1++;
									xdiffpos[y] = x;
									xposdiff[x] = y;
								}


							}
							else {
								for (z = 1; z <= spacerno; z++) {
									sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
									if (s2 == sz) {
										//If difference is legitimate

										if (*(xdiffpos + y) != x) {
											oc += 2;

											if (s1 == s3) {
												y++;

												xoverseqnumw[y + os1] = 1;

												ah1++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}

											z = spacerno;

										}
									}
								}
							}
						}
						else if (s3 != s1  && s3 != s2) {
							//If seq3 is odd one
							if (outlyer == seq3) {
								if (oc > 0) {
									if (oc > hoc)
										hoc = oc;
									oc--;
								}
								if (s1 == s2) {
									y++;
									xoverseqnumw[y + xoverwindow] = 1;

									ah0++;
									xdiffpos[y] = x;
									xposdiff[x] = y;
								}



							}
							else {
								for (z = 1; z <= spacerno; z++) {
									sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
									if (s3 == sz) {

										//If difference is legitimate
										if (*(xdiffpos + y) != x) {
											oc += 2;
											if (s1 == s2) {
												y++;
												xoverseqnumw[y + xoverwindow] = 1;

												ah0++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}

											z = spacerno;
										}
									}
								}
							}

						}
						//}
						//}
						//}
					}
					//}
				}
			}

		}
		else if (spacerflag>1) {


			for (x = 1; x < lenseq; x++) {
				xposdiff[x] = y;
				if (binarray3[x + se2] || binarray3[x + se3]) {
					s1 = seqnum[x + se1];
					s2 = seqnum[x + se2];
					s3 = seqnum[x + se3];


					//if (s1 != s2 || s1 != s3) {

					if (s1 != s2 && s1 != s3 && s2 != s3)
						g++;

					else {
						//if (s1 != 46) {
						//if (s2 != 46) {
						//if (s3 != 46) {

						if (s1 != s2  && s1 != s3) {
							//If seq1 is odd one
							for (z = 1; z <= spacerno; z++) {
								sz = seqnum[x + spacerseqs[z] * lenseq];
								if (sz == s1) {
									//If difference is legitimate
									if (seq1 != outlyer)
										oc += 2;
									else {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}
									}


									if (s2 == s3) {
										y++;

										xoverseqnumw[y + os2] = 1;
										ah2++;
										xdiffpos[y] = x;
										xposdiff[x] = y;
									}
									break;
								}
							}
						}

						else if (s2 != s1  && s2 != s3) {
							//If seq2 is odd one
							for (z = 1; z <= spacerno; z++) {
								sz = seqnum[x + spacerseqs[z] * lenseq];

								if (s2 == sz) {
									//If difference is legitimate



									if (seq2 != outlyer)
										oc += 2;
									else {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}
									}
									if (s1 == s3) {
										y++;

										xoverseqnumw[y + os1] = 1;

										ah1++;
										xdiffpos[y] = x;
										xposdiff[x] = y;
									}

									break;
								}
							}
						}
						else if (s3 != s1  && s3 != s2) {
							//If seq3 is odd one
							for (z = 1; z <= spacerno; z++) {
								sz = seqnum[x + spacerseqs[z] * lenseq];

								if (s3 == sz) {

									//If difference is legitimate

									if (seq3 != outlyer)
										oc += 2;
									else {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}
									}

									if (s1 == s2) {
										y++;
										xoverseqnumw[y + xoverwindow] = 1;

										ah0++;
										xdiffpos[y] = x;
										xposdiff[x] = y;
									}

									break;
								}
							}
							//}
							//}
							//}

						}
					}
					//}
				}
			}

		}

		wmx = y - xow;
		wpx = y + xow;
		for (b = 1; b <= xow; b++) {
			wmx++;
			wpx++;
			xoverseqnumw[b] = xoverseqnumw[wmx + xoverwindow];
			xoverseqnumw[b + lenseq + xoverwindow2] = xoverseqnumw[wmx + os1];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
			xoverseqnumw[b + lenseq2 + xoverwindow4] = xoverseqnumw[wmx + os2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

			xoverseqnumw[wpx] = xoverseqnumw[b + xoverwindow];//XOverSeqNum(X, 0)
			xoverseqnumw[wpx + lenseq + xoverwindow2] = xoverseqnumw[b + os1];
			xoverseqnumw[wpx + lenseq2 + xoverwindow4] = xoverseqnumw[b + os2];
		}


		//for (b = y+1; b < y + xoverwindow; b++)
		//	*(xdiffpos + b) = 0;


		//for (b = lenseq+1; b < lenseq + xoverwindow; b++)
		//	*(xposdiff + b) = 0;

		ah[0] = ah0;
		ah[1] = ah1;
		ah[2] = ah2;

		if (hoc < oc)
			hoc = oc;

		if (hoc < xoverwindow * 2)
			return (-y);
		else
			return(y);
	}

	int MyMathFuncs::FindSubSeqP6(int UBXSN, int UBBA, int *ah, short int spacerflag, short int outlyer, int xoverwindow, int lenseq, int nextno, int seq1, int seq2, int seq3, short int spacerno, short int *seqnum, short int *xoverseqnum, unsigned char *xoverseqnumw, short int *spacerseqs, int *xdiffpos, int *xposdiff, short int *validspacer, unsigned char *binarray3, unsigned char *binarray4, unsigned char *binarray5, unsigned char *PairwiseTripArray, unsigned char *SS255RDP)
	{

		int hoc, oc, b, x, z, lenseq2, xoverwindow2, xoverwindow4, wmx, wpx, target, D, A, B, ptaos1, ptaos2, osba, os5, os6;
		int s1, s2, s3, s2ba, s3ba, sz, ah0, ah1, ah2;//,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int g = 0;

		const int se1 = seq1*lenseq;
		const int se2 = seq2*lenseq;
		const int se3 = seq3*lenseq;
		oc = 0;
		hoc = 0;
		const int xow = xoverwindow;
		lenseq2 = lenseq * 2;
		xoverwindow2 = xoverwindow * 2;
		xoverwindow4 = xoverwindow * 4;
		const int os1 = xoverwindow + lenseq + xoverwindow2;
		const int os2 = xoverwindow + lenseq2 + xoverwindow4;
		os5 = UBXSN + 1;
		ah0 = 0;
		ah1 = 0;
		ah2 = 0;

		ptaos1 = 43;
		ptaos2 = 43 * 43;
		
		
		target = (lenseq + xoverwindow2 + 1) * 3;
		osba = UBBA + 1;
		//s1ba = seq1*osba;
		s2ba = seq2*osba;
		s3ba = seq3*osba;
		for (x = 1; x < target; x++)
			xoverseqnumw[x] = 0;

		if (spacerflag == 0) {
			//so1 = seq1*lenseq;
			//so2 = seq2*lenseq;
			//so3 = seq3*lenseq;

			y = 0;
			x = -1;

			for (D = 1; D <= UBBA; D++) {
				os6 = D + s3ba;
				A = PairwiseTripArray[binarray3[D + s2ba] + binarray3[os6] * ptaos1 + binarray4[os6] * ptaos2];
				
				//if (A < 63) {
					//*ah[0] += A;
					B = SS255RDP[A];
					x++;
					if (B != 3) {
						xdiffpos[++y] = x;
						ah[B]++;
						xoverseqnumw[xow + y + B*os5] = 1;
					}
					xposdiff[x] = y;
					B = SS255RDP[A + 64];
					x++;
					if (B != 3) {

						xdiffpos[++y] = x;
						ah[B]++;
						xoverseqnumw[xow + y + B*os5] = 1;
					}
					xposdiff[x] = y;
					B = SS255RDP[A + 128];
					x++;
					if (B != 3) {
						xdiffpos[++y] = x;
						ah[B]++;
						xoverseqnumw[xow + y + B*os5] = 1;
					}
					xposdiff[x] = y;
			//	}
				//else {
			//
			//		xposdiff[x + 1] = y;
			//		xposdiff[x + 2] = y;
			//		xposdiff[x + 3] = y;
					
			//		x = x + 3;
			//	}

			}
			//LenXOverSeq = Y 
			ah0 = ah[0];
			ah1 = ah[1];
			ah2 = ah[2];

			hoc = xoverwindow * 10;



		}

		else if (spacerflag == 1) {
			for (x = 1; x < lenseq; x++) {
				*(xposdiff + x) = y;
				if (binarray3[x + se2] || binarray3[x + se3]) {

					s1 = *(seqnum + x + se1);
					s2 = *(seqnum + x + se2);
					s3 = *(seqnum + x + se3);

					//if (s1 != s2 || s1 != s3) {

					if (s1 == s2 || s1 == s3 || s2 == s3) {

						//if (s1 != 46) {
						//	if (s2 != 46) {
						//		if (s3 != 46) {
						if (binarray3[x + se2] && binarray3[x + se3]) {
							//If seq1 is odd one


							if (outlyer == seq1) {
								if (oc > 0) {
									if (oc > hoc)
										hoc = oc;
									oc--;
								}


								if (s2 == s3) {
									y++;

									xoverseqnumw[y + os2] = 1;
									ah2++;
									xdiffpos[y] = x;
									xposdiff[x] = y;
								}

							}
							else {
								for (z = 1; z <= spacerno; z++) {
									sz = *(seqnum + x + *(spacerseqs + z) * lenseq);
									if (sz == s1) {
										//If difference is legitimate
										if (*(xdiffpos + y) != x) {
											oc += 2;


											if (s2 == s3) {
												y++;

												xoverseqnumw[y + os2] = 1;
												ah2++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}
											z = spacerno;

										}

									}

								}

							}
						}
						else if (s2 != s1  && s2 != s3) {
							//If seq2 is odd one

							if (outlyer == seq2) {
								if (oc > 0) {
									if (oc > hoc)
										hoc = oc;
									oc--;
								}

								if (s1 == s3) {
									y++;

									xoverseqnumw[y + os1] = 1;

									ah1++;
									xdiffpos[y] = x;
									xposdiff[x] = y;
								}


							}
							else {
								for (z = 1; z <= spacerno; z++) {
									sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
									if (s2 == sz) {
										//If difference is legitimate

										if (*(xdiffpos + y) != x) {
											oc += 2;

											if (s1 == s3) {
												y++;

												xoverseqnumw[y + os1] = 1;

												ah1++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}

											z = spacerno;

										}
									}
								}
							}
						}
						else if (s3 != s1  && s3 != s2) {
							//If seq3 is odd one
							if (outlyer == seq3) {
								if (oc > 0) {
									if (oc > hoc)
										hoc = oc;
									oc--;
								}
								if (s1 == s2) {
									y++;
									xoverseqnumw[y + xoverwindow] = 1;

									ah0++;
									xdiffpos[y] = x;
									xposdiff[x] = y;
								}



							}
							else {
								for (z = 1; z <= spacerno; z++) {
									sz = *(seqnum + x + *(spacerseqs + z)*lenseq);
									if (s3 == sz) {

										//If difference is legitimate
										if (*(xdiffpos + y) != x) {
											oc += 2;
											if (s1 == s2) {
												y++;
												xoverseqnumw[y + xoverwindow] = 1;

												ah0++;
												xdiffpos[y] = x;
												xposdiff[x] = y;
											}

											z = spacerno;
										}
									}
								}
							}

						}
						//}
						//}
						//}
					}
					//}
				}
			}

		}
		else if (spacerflag>1) {


			for (x = 1; x < lenseq; x++) {
				xposdiff[x] = y;
				if (binarray3[x + se2] || binarray3[x + se3]) {
					s1 = seqnum[x + se1];
					s2 = seqnum[x + se2];
					s3 = seqnum[x + se3];


					//if (s1 != s2 || s1 != s3) {

					if (s1 != s2 && s1 != s3 && s2 != s3)
						g++;

					else {
						//if (s1 != 46) {
						//if (s2 != 46) {
						//if (s3 != 46) {

						if (s1 != s2  && s1 != s3) {
							//If seq1 is odd one
							for (z = 1; z <= spacerno; z++) {
								sz = seqnum[x + spacerseqs[z] * lenseq];
								if (sz == s1) {
									//If difference is legitimate
									if (seq1 != outlyer)
										oc += 2;
									else {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}
									}


									if (s2 == s3) {
										y++;

										xoverseqnumw[y + os2] = 1;
										ah2++;
										xdiffpos[y] = x;
										xposdiff[x] = y;
									}
									break;
								}
							}
						}

						else if (s2 != s1  && s2 != s3) {
							//If seq2 is odd one
							for (z = 1; z <= spacerno; z++) {
								sz = seqnum[x + spacerseqs[z] * lenseq];

								if (s2 == sz) {
									//If difference is legitimate



									if (seq2 != outlyer)
										oc += 2;
									else {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}
									}
									if (s1 == s3) {
										y++;

										xoverseqnumw[y + os1] = 1;

										ah1++;
										xdiffpos[y] = x;
										xposdiff[x] = y;
									}

									break;
								}
							}
						}
						else if (s3 != s1  && s3 != s2) {
							//If seq3 is odd one
							for (z = 1; z <= spacerno; z++) {
								sz = seqnum[x + spacerseqs[z] * lenseq];

								if (s3 == sz) {

									//If difference is legitimate

									if (seq3 != outlyer)
										oc += 2;
									else {
										if (oc > 0) {
											if (oc > hoc)
												hoc = oc;
											oc--;
										}
									}

									if (s1 == s2) {
										y++;
										xoverseqnumw[y + xoverwindow] = 1;

										ah0++;
										xdiffpos[y] = x;
										xposdiff[x] = y;
									}

									break;
								}
							}
							//}
							//}
							//}

						}
					}
					//}
				}
			}

		}

		wmx = y - xow;
		wpx = y + xow;
		for (b = 1; b <= xow; b++) {
			wmx++;
			wpx++;
			xoverseqnumw[b] = xoverseqnumw[wmx + xoverwindow];
			xoverseqnumw[b + lenseq + xoverwindow2] = xoverseqnumw[wmx + os1];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
			xoverseqnumw[b + lenseq2 + xoverwindow4] = xoverseqnumw[wmx + os2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

			xoverseqnumw[wpx] = xoverseqnumw[b + xoverwindow];//XOverSeqNum(X, 0)
			xoverseqnumw[wpx + lenseq + xoverwindow2] = xoverseqnumw[b + os1];
			xoverseqnumw[wpx + lenseq2 + xoverwindow4] = xoverseqnumw[b + os2];
		}


		//for (b = y+1; b < y + xoverwindow; b++)
		//	*(xdiffpos + b) = 0;


		//for (b = lenseq+1; b < lenseq + xoverwindow; b++)
		//	*(xposdiff + b) = 0;

		ah[0] = ah0;
		ah[1] = ah1;
		ah[2] = ah2;

		if (hoc < oc)
			hoc = oc;

		if (hoc < xoverwindow * 2)
			return (-y);
		else
			return(y);
	}
	int MyMathFuncs::MakeXPD2(int lenxoseq,  int *xdiffpos, int *xposdiff)
	{
		int x, y, val;
			for (x = 0; x < lenxoseq; x++) {
				val = xdiffpos[x];
				for (y = xdiffpos[x]; y < xdiffpos[x + 1]; y++)
					xposdiff[y] = val;

			}
			
		return(1);
	}


	int MyMathFuncs::MakeBinArrayP(int Seq1, int LSeq, int Nextno, short int *Maskseq, short int *SeqNum, unsigned char *BinArray) {
		int Seq2, X, os, se1, se2, target, S1, S2;
		os = LSeq + 1;
		se1 = os*Seq1;
		target = LSeq + Nextno*os;
		for (Seq2 = 0; Seq2 <= target; Seq2++)
			BinArray[Seq2] = 0;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private(Seq2, se2, X, S1, S2)
		for (Seq2 = Seq1 + 1; Seq2 <= Nextno; Seq2++) {
			if (Maskseq[Seq2] == 0) {
				se2 = os*Seq2;
//#pragma omp parallel for private(X, S1, S2)
				for (X = 1; X <= LSeq; X++) {
					S1 = SeqNum[X + se1];
					S2 = SeqNum[X + se2];
					if (S1 != S2) {

						if (S1 != 46) {
							if (S2 != 46)
								BinArray[X + se2] = 1;
							else
								BinArray[X + se2] = 3;
						}
						else
							BinArray[X + se2] = 2;

					}
				}
			}
		}
		omp_set_num_threads(2);
		return(1);
	}

	

	
	//int  MyMathFuncs::MakeBinArray2P(int Seq1, int LSeq, int Nextno, short int *Maskseq, short int *SeqNum, unsigned char *BinArray, int *slookupnum, int *slookup, unsigned char *isin, int *tracesub, int *actualsize, int MinSeqSize) {
	//	int Seq2, X, os, se1, se2, target, S1, S2, g;
	//	os = LSeq + 1;
	//	se1 = os*Seq1;
	//	target = LSeq + Nextno*os;
	//	for (Seq2 = 0; Seq2 <= target; Seq2++)
	//		BinArray[Seq2] = 0;

	//	for (g = 1; g <= slookupnum[0]; g++) {
	//		//for (Seq2 = Seq1 + 1; Seq2 <= Nextno; Seq2++){
	//		Seq2 = slookup[g * 2];
	//		if (actualsize[Seq2] > MinSeqSize) {
	//			if (isin[Seq2] == 0) {
	//				if (tracesub[Seq1] != tracesub[Seq2]) {
	//					//if (Maskseq[Seq2] == 0){
	//					se2 = os*Seq2;
	//					for (X = 1; X <= LSeq; X++) {
	//						S1 = SeqNum[X + se1];
	//						S2 = SeqNum[X + se2];
	//						if (S1 != S2) {

	//							if (S1 != 46) {
	//								if (S2 != 46)
	//									BinArray[X + se2] = 1;

	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//	return(1);
	//}

	

	int MyMathFuncs::FindSubSeqP2(int UBXSN, int UBVO, int UBCS, int XoverWindow, int lenseq, int A, int B, int C, int *AH, unsigned char *CompressedSeqs3, unsigned char *XoverSeqNumW, int *XDP, int *XPD, unsigned char *SkipTrip, int *FindSS0)



	{

		int  x, lenseq2, xoverwindow2, xoverwindow4,target;
		int   os1, os2, ah0, ah1, ah2, os3, os5;//,so1,so2,so3;
		int holder = 0;
		int y = 0;
		int g = 0;




		lenseq2 = lenseq * 2;
		xoverwindow2 = XoverWindow * 2;
		xoverwindow4 = XoverWindow * 4;
		os1 = XoverWindow + lenseq + xoverwindow2;
		os2 = XoverWindow + lenseq2 + xoverwindow4;
		ah0 = 0;
		ah1 = 0;
		ah2 = 0;

		//if (spacerflag == 0) {
			//so1 = seq1*lenseq;
			//so2 = seq2*lenseq;
			//so3 = seq3*lenseq;
		target = (lenseq + xoverwindow2 + 1) * 3;
		for (x = 0; x < target; x++)
			XoverSeqNumW[x] = 0;




		int Z, Y, V1, V2, V3, aoff, boff, coff, voff, voff2, voff3, vo, FS1, xoff, SkipT, YV, voff6, newsize;
		Z = 0;
		xoff = UBXSN + 1;
		aoff = A*(UBCS + 1);
		boff = B*(UBCS + 1);
		coff = C*(UBCS + 1);
		voff = UBVO + 1;
		voff2 = voff*voff;
		voff3 = voff*voff2;
		voff6 = 2 * voff3;
		os3 = voff3 * 3;
		os5 = voff3 * 5;
		newsize = (int)(lenseq / 3) - 1;
		for (Y = 1; Y <= UBCS; Y++) {
			V1 = CompressedSeqs3[Y + aoff];
			V2 = CompressedSeqs3[Y + boff];
			V3 = CompressedSeqs3[Y + coff];
			vo = V1 + V2*voff + V3*voff2;
			SkipT = SkipTrip[vo];
			YV = (Y - 1) * 3 + 1;
			if (SkipT > 0) {
				

				if (FindSS0[vo] == 1) {
					Z++;
					FS1 = FindSS0[vo + voff3];
					AH[FS1] = AH[FS1] + 1;
					XoverSeqNumW[Z + XoverWindow + FS1*UBXSN] = 1;
					XDP[Z] = YV;
					SkipT = SkipT - 1;
				}
				XPD[YV] = Z;

				if (SkipT > 0) {
					if (FindSS0[vo + voff6] == 1) {
						Z++;
						FS1 = FS1 = FindSS0[vo + os3];
						AH[FS1] = AH[FS1] + 1;
						XoverSeqNumW[Z + XoverWindow + FS1*UBXSN] = 1;
						XDP[Z] = YV + 1;
						SkipT = SkipT - 1;

					}
					
					XPD[YV+1] = Z;
					if (SkipT > 0) {

						Z++;
						FS1 = FindSS0[vo + os5];
						AH[FS1] = AH[FS1] + 1;
						XoverSeqNumW[Z + XoverWindow + FS1*UBXSN] = 1;
						XDP[Z] = YV + 2;

					
					}
					XPD[YV + 2] = Z;
				}
				else {
					XPD[YV+1] = Z;
					XPD[YV + 2] = Z;

				}

			}
			else {
				XPD[YV] = Z;
				XPD[YV+1] = Z;
				XPD[YV + 2] = Z;

			}

		}





		int WMX, WPX;
		WMX = Z - XoverWindow;
		WPX = Z + XoverWindow;
		for (Y = 1; Y <= XoverWindow; Y++){
			WMX = WMX + 1;
			WPX = WPX + 1;
			XoverSeqNumW[Y] = XoverSeqNumW[WMX + XoverWindow];
			XoverSeqNumW[Y + xoff] = XoverSeqNumW[WMX + XoverWindow + xoff];
			XoverSeqNumW[Y + xoff * 2] = XoverSeqNumW[WMX + XoverWindow + xoff * 2];

			XoverSeqNumW[WPX] = XoverSeqNumW[Y + XoverWindow];
			XoverSeqNumW[WPX + xoff] = XoverSeqNumW[Y + XoverWindow + xoff];
			XoverSeqNumW[WPX + xoff*2] = XoverSeqNumW[Y + XoverWindow + xoff*2];
		}
			


		//wmx = Z - XoverWindow;
		//wpx = Z + XoverWindow;
		//for (b = 1; b <= XoverWindow; b++) {
		//	wmx++;
		//	wpx++;





		//	XoverSeqNumW[b] = XoverSeqNumW[wmx + XoverWindow];
		//	XoverSeqNumW[b + lenseq + xoverwindow2] = XoverSeqNumW[wmx + os1];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 1)
		//	XoverSeqNumW[b + lenseq2 + xoverwindow4] = XoverSeqNumW[wmx + os2];//XOverSeqNum(LenXOverSeq - XOverWindow + X, 2)

		//	XoverSeqNumW[wpx] = XoverSeqNumW[b + XoverWindow];//XOverSeqNum(X, 0)
		//	XoverSeqNumW[wpx + lenseq + xoverwindow2] = XoverSeqNumW[b + os1];
		//	XoverSeqNumW[wpx + lenseq2 + xoverwindow4] = XoverSeqNumW[b + os2];
		//}

		
		
		return(Z);
	}


	int  MyMathFuncs::FindFirstCOP(int x, int MedHomol, int  HighHomol, int LenXOverSeq, int UBXOHN, int *XOverHomologyNum) 
	{
		int os, X;
		os = UBXOHN + 1;
		X = x;
		while (XOverHomologyNum[X + (MedHomol - 1)*os] > XOverHomologyNum[X + (HighHomol - 1)*os]) {
			X = X + 1;
		}
		if (XOverHomologyNum[LenXOverSeq + (MedHomol - 1)*os] <= XOverHomologyNum[LenXOverSeq + (HighHomol - 1)*os]) {
			XOverHomologyNum[LenXOverSeq + (MedHomol - 1)*os] = XOverHomologyNum[1 +(MedHomol - 1)*os];
			XOverHomologyNum[LenXOverSeq + (HighHomol - 1)*os] = XOverHomologyNum[1 + (HighHomol - 1)*os];
		}
			return(X);
	}


	double MyMathFuncs::DefineEventP(int ShortOutFlag, int LongWindedFlag, int MedHomol, int HighHomol, int LowHomol, int TargetX, int CircularFlag, int XX, int  XOverWindow, int  lenseq, int  LenXoverSeq, int  SeqDaughter, int  SeqMinorP, int *EndFlag, int  *Be, int  *En, int  *NCommon, int  *XOverLength, char *XOverSeqNum, int *XDiffPos, int *XOverHomologyNum)
	{
		//xoverseqnum lenseq,2
		//xoverhomologynum LenXoverSeq + XOverWindow * 2, 2
		int X, off1, off2, off3, off4, off6, Storex, Store, TC, NC;
		TC = 0;

		off1 = LenXoverSeq + XOverWindow * 2 + 1;
		off2 = (MedHomol - 1)*off1;
		off3 = (HighHomol - 1)*off1;
		off6 = (LowHomol - 1)*off1;

		if ((SeqDaughter == 0 && SeqMinorP == 1) || (SeqDaughter == 1 && SeqMinorP == 0))
			TC = 0;
		else if ((SeqDaughter == 0 && SeqMinorP == 2) || (SeqDaughter == 2 && SeqMinorP == 0))
			TC = 1;
		else if ((SeqDaughter == 2 && SeqMinorP == 1) || (SeqDaughter == 1 && SeqMinorP == 2))
			TC = 2;

		off4 = XOverWindow + TC*(lenseq + 1 + XOverWindow * 2);




		X = XX;
		if (X > 1) {
			while (X != 1) {
				X--;

				if (XOverSeqNum[X + off4] == 0) {// || s1 == 64){
					X++;
					break;
				}
			}
		}

		NC = 0;
		if (XOverSeqNum[X + off4] == 0) {
			while (X <= LenXoverSeq + 1) {
				X++;

				if (XOverSeqNum[X + off4])// && s1 != 64)
					break;
				if (X == LenXoverSeq + 1) {
					if (NC == 1)
						break;
					X = 0;
					NC++;
					*EndFlag = 1;
				}
			}
		}



		if (CircularFlag == 0) {
			if (X == 1) {
				if (XDiffPos[X] < TargetX)
					*Be = 1;
				else
					*Be = XDiffPos[X];
			}
			else
				*Be = XDiffPos[X];
		}
		else
			*Be = XDiffPos[X];



		while (X <= LenXoverSeq) {


			if (XOverSeqNum[X + off4])
				*NCommon = *NCommon + 1;

			*XOverLength = *XOverLength + 1;

			if (*XOverLength >= LenXoverSeq)
				break;

			X++;


			//if (s1 == 64)
			//	break;

			if (X > LenXoverSeq) {

				if (CircularFlag == 1) {
					X = 0;
					*EndFlag = 1;
				}
				else {
					X = LenXoverSeq + 2 * XOverWindow + 1;
					break;
				}
			}

			if (XOverHomologyNum[X + off2] < XOverHomologyNum[X + off3] || XOverHomologyNum[X + off2] < XOverHomologyNum[X + off6]) {
				if (XOverSeqNum[X + off4] == 0)
					break;
			}

		}

		Storex = X;
		X--;
		Store = 0;

		while (X > 0) {
			if (XOverSeqNum[X + off4] == 0) {// || XOverSeqNum[X + off4] == 64){
				if (X < LenXoverSeq) {
					if (X > 0) {
						X--;
						Store = Store + 1;
						*XOverLength = *XOverLength - 1;
						if (*XOverLength == 1)
							break;
					}
					else {

						X = LenXoverSeq + 1;


						break;
					}
				}
				else
					break;
			}
			else
				break;
		}



		if (X == LenXoverSeq && CircularFlag == 0) {
			if (ShortOutFlag == 0 || ShortOutFlag == 6 || ShortOutFlag == 10)
				*En = lenseq;
			else
				*En = XDiffPos[LenXoverSeq];

		}
		else {
			if (X >= LenXoverSeq) {
				if (LongWindedFlag == 0)
					XDiffPos[X] = 0;
				else
					X = LenXoverSeq + 1;

			}
			else if (X < 1)
				X = LenXoverSeq + 1;

			*En = XDiffPos[X];
		}

		X = Storex + 1;
		if (*En == 0) {

			if (ShortOutFlag == 0 || ShortOutFlag == 6 || ShortOutFlag == 10)
				*En = lenseq - 1;
			else
				*En = XDiffPos[LenXoverSeq];

		}
		return((double)(X));

	}
	
	int MyMathFuncs::GoRightP(int Seq1, int Seq2, int Seq3, int CircularFlag, int startpos, int LS, int UBMD, unsigned char *MissingData) {
		int CycleX, Z, os1, os2, os3, X;
		os1 = Seq1*(UBMD + 1);
		os2 = Seq2*(UBMD + 1);
		os3 = Seq3*(UBMD + 1);
		Z = startpos;
		CycleX = 0;
		for (X = 1; X <= LS; X++) {
			Z = Z + 1;

			if (Z > LS) {
				if (CircularFlag == 0) {
					Z = 1;
					break;
				}
				else {
					Z = 1;
					CycleX = CycleX + 1;
					if (CycleX == 3)
						break;

				}
			}
			if (MissingData[Z + os1] == 0 && MissingData[Z + os2] == 0 && MissingData[Z + os3] == 0)
				break;


		}
		return (Z);
	}

	int MyMathFuncs::GoLeftP(int Seq1, int Seq2, int Seq3, int CircularFlag, int startpos, int LS, int UBMD, unsigned char *MissingData) {
		int X, CycleX, Z, os1, os2, os3;
		os1 = Seq1*(UBMD + 1);
		os2 = Seq2*(UBMD + 1);
		os3 = Seq3*(UBMD + 1);
		Z = startpos;
		CycleX = 0;
		for (X = 1; X <= LS; X++) {
			Z = Z - 1;

			if (Z < 1) {
				if (CircularFlag == 0) {
					Z = LS;
					break;
				}
				else {
					Z = LS;
					CycleX = CycleX + 1;
					if (CycleX == 3)
						break;

				}
			}
			if (MissingData[Z + os1] == 0 && MissingData[Z + os2] == 0 && MissingData[Z + os3] == 0)
				break;


		}
		return (Z);
	}
	
	int MyMathFuncs::DefineEventP2(int UBXOHN,int ShortOutFlag, int LongWindedFlag, int MedHomol, int HighHomol, int LowHomol, int TargetX, int CircularFlag, int XX, int  XOverWindow, int  lenseq, int  LenXoverSeq, int  SeqDaughter, int  SeqMinorP, int *EndFlag, int  *Be, int  *En, int  *NCommon, int  *XOverLength, char *XOverSeqNum, int *XOverHomologyNum) 
	{
		//xoverseqnum lenseq,2
		//xoverhomologynum LenXoverSeq + XOverWindow * 2, 2
		int X, off1, off2, off3, off4, off6, off1b, off2b, off3b, off4b, off6b, Storex,  TC, NC;
		TC = 0;

		off1 = LenXoverSeq + XOverWindow * 2 + 1;
		off2 = (MedHomol - 1)*off1;
		off3 = (HighHomol - 1)*off1;
		off6 = (LowHomol - 1)*off1;

		off1b = UBXOHN + 1;
		off2b = (MedHomol - 1)*off1b;
		off3b = (HighHomol - 1)*off1b;
		off6b = (LowHomol - 1)*off1b;

		if ((SeqDaughter == 0 && SeqMinorP == 1) || (SeqDaughter == 1 && SeqMinorP == 0))
			TC = 0;
		else if ((SeqDaughter == 0 && SeqMinorP == 2) || (SeqDaughter == 2 && SeqMinorP == 0))
			TC = 1;
		else if ((SeqDaughter == 2 && SeqMinorP == 1) || (SeqDaughter == 1 && SeqMinorP == 2))
			TC = 2;

		off4 = XOverWindow + TC*(lenseq + 1 + XOverWindow * 2);




		X = XX;
		if (X > 1) {
			while (X != 1) {
				X--;

				if (XOverSeqNum[X + off4] == 0) {// || s1 == 64){
					X++;
					break;
				}
			}
		}

		NC = 0;
		if (XOverSeqNum[X + off4] == 0) {
			while (X <= LenXoverSeq + 1) {
				X++;

				if (XOverSeqNum[X + off4])// && s1 != 64)
					break;
				if (X == LenXoverSeq + 1) {
					if (NC == 1)
						break;
					X = 0;
					NC++;
					*EndFlag = 1;
				}
			}
		}


		*Be = X;
		


		while (X <= LenXoverSeq) {


			if (XOverSeqNum[X + off4])
				*NCommon = *NCommon + 1;

			*XOverLength = *XOverLength + 1;

			if (*XOverLength >= LenXoverSeq)
				break;

			X++;


			//if (s1 == 64)
			//	break;

			if (X > LenXoverSeq) {

				if (CircularFlag == 1) {
					X = 0;
					*EndFlag = 1;
				}
				else {
					X = LenXoverSeq + 2 * XOverWindow + 1;
					break;
				}
			}

			if (XOverHomologyNum[X + off2b] < XOverHomologyNum[X + off3b] || XOverHomologyNum[X + off2b] < XOverHomologyNum[X + off6b]) {
				if (XOverSeqNum[X + off4] == 0)
					break;
			}

		}

		Storex = X;
		X--;
		
		while (X > 0) {
			if (XOverSeqNum[X + off4] == 0) {// || XOverSeqNum[X + off4] == 64){
				if (X < LenXoverSeq) {
					if (X > 0) {
						X--;
						
						*XOverLength = *XOverLength - 1;
						if (*XOverLength == 1)
							break;
					}
					else {

						X = LenXoverSeq + 1;


						break;
					}
				}
				else
					break;
			}
			else
				break;
		}

		
		*En = X;
		

		X = Storex + 1;
		
		
		return((int)(X));

	}

	int MyMathFuncs::CountVSites(int x, int lseq, int SA, int SB, int SX, int Epos, int UBSN, short int *Seqnum) {
		
		int VSiteNum, osa, osb,osx;
		VSiteNum = 0;
		osa = (UBSN + 1)*SA;
		osb = (UBSN + 1)*SB;
		osx = (UBSN + 1)*SX;
		while (x != Epos){

			if (Seqnum[x + osa] != 46) {
				if (Seqnum[x + osb] != 46) {
					if (Seqnum[x + osx] != 46) {
						if (Seqnum[x + osa] != Seqnum[x + osb] || Seqnum[x + osa] != Seqnum[x + osx])
							VSiteNum++;

					}
				}
			}

			x++;
			if (x == Epos)
				return(VSiteNum);
			if (x > lseq) {
				x = 1;
				if (x == Epos)
					return(VSiteNum);
			}

		}

		return(VSiteNum);
	}


	int MyMathFuncs::FindNextP(int UBXOHN, int StartPosX, int HighHomol, int MedHomol, int LowHomol, int LenXoverSeq, int xoverwindow, int *XOverHomologyNum)
	{
		int X, off1, off2, off3;
		int limit = xoverwindow * 2 + 1;

		//*(xoverhomologynum+1) = t1/limit;
		//*(xoverhomologynum+1+limit+lenxoverseq) = t2/limit;
		//*(xoverhomologynum+1+limit*2+lenxoverseq*2 ) = t3/limit;
		off1 = (MedHomol - 1)*(UBXOHN + 1);
		off2 = (HighHomol - 1)*(UBXOHN + 1);
		off3 = (LowHomol - 1)*(UBXOHN + 1);

		for (X = StartPosX; X <= LenXoverSeq; X++) {

			if (XOverHomologyNum[X + off1] > XOverHomologyNum[X + off2]) {
				if (XOverHomologyNum[X + off1] > XOverHomologyNum[X + off3])
					
						return X;
			}

		}
		return -1;
	}

	int MyMathFuncs::FindNextPB(int UBXOHN, int StartPosX, int HighHomol, int MedHomol, int LowHomol, int LenXoverSeq, int xoverwindow, int *XOverHomologyNum)
	{
		int X, off1, off2, off3;
		int limit = xoverwindow * 2 + 1;

		//*(xoverhomologynum+1) = t1/limit;
		//*(xoverhomologynum+1+limit+lenxoverseq) = t2/limit;
		//*(xoverhomologynum+1+limit*2+lenxoverseq*2 ) = t3/limit;
		off1 = (MedHomol - 1)*(UBXOHN + 1);
		off2 = (HighHomol - 1)*(UBXOHN + 1);
		off3 = (LowHomol - 1)*(UBXOHN + 1);

		for (X = StartPosX; X <= LenXoverSeq; X = X + 2) {

			if (XOverHomologyNum[X + off1] > XOverHomologyNum[X + off2]) {
				if (XOverHomologyNum[X + off1] > XOverHomologyNum[X + off3])
					if (X > 1) {
						if (XOverHomologyNum[X - 1 + off1] > XOverHomologyNum[X - 1 + off2]) {
							if (XOverHomologyNum[X - 1 + off1] > XOverHomologyNum[X - 1 + off3])
								return(X - 1);
							else
								return(X);
						}
						else
							return(X);
					}
					else
						return X;
			}

		}
		return -1;
	}

	int  MyMathFuncs::MarkDones(int Nextno, int lseq,int STA, int ENA, int A1, int A2, int A3, int UBDS1, int UBPXO,  unsigned char *DoneSeq, short int *PCurrentXOver, XOVERDEFINE *PXOList)
	{
		int x, Y, dso, dso2, pxo, pxo2, b1, b2, b3, MatchesX, STB, ENB, C;
		//double CPVal;
		dso = UBDS1 + 1;
		pxo = UBPXO + 1;
		C = 0;
		for (x = 0; x <= Nextno; x++){
                    
			for (Y = 1; Y <= PCurrentXOver[x]; Y++) {
                        
				dso2 = x + dso*Y;
                if (DoneSeq[dso2] == 0){
					pxo2 = x + pxo*Y;
					//CPVal = PXOList[pxo2].Probability;
					if (PXOList[pxo2].Probability > 0) {
						if (PXOList[pxo2].Beginning != PXOList[pxo2].Ending){
							b1 = PXOList[pxo2].Daughter;
							b2 = PXOList[pxo2].MinorP;
							b3 = PXOList[pxo2].MajorP;
							MatchesX = 0;

							if (A1 == b1 || A1 == b2 || A1 == b3) 
								MatchesX++;
									
							if (A2 == b1 || A2 == b2 || A2 == b3)
								MatchesX++;
									
							if (A3 == b1 || A3 == b2 || A3 == b3)
								MatchesX++;
									
							if (MatchesX > 1){
								STB = PXOList[pxo2].Beginning;
								ENB = PXOList[pxo2].Ending;
								if (STA < ENA) {
									if (STB < ENB) {
										if (abs(STB - STA) + abs(ENA - ENB) < 20) {
											DoneSeq[dso2] = 1;
											C++;
										}
									}
																		
									else{//STB>ENB meaning there is wrapping around the end 
										
										if (abs(STB - lseq - STA) + abs(ENA - ENB) < 20 || abs(STB - STA) + abs(lseq - ENA - ENB) < 20) {
											DoneSeq[dso2] = 1;
											C++;
										}
										
									}	
								}
								else {//STA>ENA
									if (STB < ENB) {
										if (abs(STA - lseq - STB) + abs(ENA - ENB) < 20 || abs(STB - STA) + abs(lseq - ENB - ENA) < 20) {
											DoneSeq[dso2] = 1;
											C++;
										}
									}
									else {//there is wrapping around the end for both
										if (abs(STB - STA) + abs(ENA - ENB) < 20) {
											DoneSeq[dso2] = 1;
											C++;
										}
										
									}
								}
							}
						}
					}
					else
						DoneSeq[dso2] = 1;
                            
				}
			}
        }
		return(C);	
	}

	int  MyMathFuncs::UpdatePlotsCP(int UBAD,float ff, HDC Pict, int LSeq, short int P1, short int P2, short int P3, short int P4, int StepSize, float XFactor, float oDMax, float oPMax, int MaxHits, int *Decompress, float *PDistPlt, float *ProbPlt, int *HitPlt) {
		int os;
		float a, b, c, hc, mh;
		

		a = 0;
		b = 0;
		c = 0;


//		//omp_set_num_threads(3);
//#pragma omp parallel 
//		{
//#pragma omp sections private (Pict)
//			{
//#pragma omp section 
//				{
					if (oDMax > 0) {
						os = P2 - P1 - 10;
						MoveToEx(Pict, 30 + Decompress[(int)(1/ff)] * XFactor + XFactor, P2 - 5 - (PDistPlt[1] / oDMax) * os, 0);
						for (a = 2; a <= UBAD; a = a + StepSize) {
							
							LineTo(Pict, 30 + Decompress[(int)(a / ff)] * XFactor + XFactor, P2 - 5 - (PDistPlt[(int)(a)] / oDMax)* os);
						}
					}
//				}
//#pragma omp section
//				{
					if (oPMax > 0) {
						os = P3 - P2 - 10;
						MoveToEx(Pict, 30 + Decompress[(int)(1 / ff)] * XFactor + XFactor, P3 - 5 - (ProbPlt[1] / oPMax) * os, 0);
						for (b = 2; b <= UBAD; b = b + StepSize) {
							
							LineTo(Pict, 30 + Decompress[(int)(b / ff)] * XFactor + XFactor, P3 - 5 - (ProbPlt[(int)(b)] / oPMax) * os);
						}
//					}
//				}
//#pragma omp section
//				{
					mh = (float)(MaxHits);
					if (MaxHits > 0) {
						os = P4 - P3 - 10;
						hc = (float)(HitPlt[1]);
						MoveToEx(Pict, 30 + Decompress[(int)(1 / ff)] * XFactor + XFactor, P4 - 5 - (hc / mh) * os, 0);
						for (c = 2; c <= UBAD; c = c + StepSize) {
							hc = (float)(HitPlt[(int)(c)]);
							LineTo(Pict, 30 + Decompress[(int)(c / ff)] * XFactor + XFactor, P4 - 5 - (hc / mh) * os);
						}
					}
				}
//			}
//		}
		return(1);
	}


	int MyMathFuncs::FindBestRecSignalP(char DoneTarget, int NextNo, int UB, int UB2, double *LowP, char *DoneSeq, int *Trace, short int *PCurrentXOver, XOVERDEFINE *PXOList) {
		int  X, Y, TotalNoRecombinants,os, os2;
		double LowPX, CPVal;
		TotalNoRecombinants = 0;
		LowPX = *LowP;
		for (X = 0; X <= NextNo; X++)
			TotalNoRecombinants = TotalNoRecombinants + PCurrentXOver[X];
		for (Y = 1; Y <= UB2; Y++) {
			os = Y*(UB + 1);
			os2 = Y*(NextNo + 1);
			for (X = 0; X <= NextNo; X++) {
				if (Y <= PCurrentXOver[X]) {

					if (DoneSeq[X + os2] == DoneTarget) {
						//CPVal = PXOList(X, Y).Probability

						CPVal = PXOList[X + os].Probability;
						//*LowP = (double)((PXOList[1+ 5*(NextNo+1)].PermPVal));
						//TotalNoRecombinants = (PXOList[1+ 5*(NextNo+1)].Eventnumber);
						//return(TotalNoRecombinants);
						if (CPVal > 0 && CPVal < LowPX) {
							if (PXOList[X + os].Beginning != PXOList[X + os].Ending) {
								LowPX = CPVal;
								Trace[0] = X;
								Trace[1] = Y;

							}


						}
						else if (CPVal == LowPX && Trace[0]>X) {
							if (PXOList[X + os].Beginning != PXOList[X + os].Ending) {
								LowPX = CPVal;
								Trace[0] = X;
								Trace[1] = Y;

							}

						}
					}
				}
			}
		}
		*LowP = LowPX;
		return(TotalNoRecombinants);
	}

	int MyMathFuncs::FindBestRecSignalP2(char DoneTarget, int NextNo, int UB, int UB2, double *LowP, char *DoneSeq, int *Trace, short int *PCurrentXOver, double *TestPVs) {
		int  X, Y, TotalNoRecombinants, os, os2, TV;
		double LowPX, CPVal;
		TotalNoRecombinants = 0;
		LowPX = *LowP;
		for (X = 0; X <= NextNo; X++)
			TotalNoRecombinants = TotalNoRecombinants + PCurrentXOver[X];
		
		if (UB2 < NextNo )
			TV = UB2;
		else
			TV = NextNo;

		for (X = 0; X <= TV; X++) {
			os = X*(UB + 1);
			
			for (Y = 1; Y <= PCurrentXOver[X]; Y++) {
				//if (Y <= PCurrentXOver[X]) {
					os2 = Y*(NextNo + 1);

					if (DoneSeq[X + os2] == DoneTarget) {
						//CPVal = PXOList(X, Y).Probability
						
						CPVal = TestPVs[Y + os];
						//*LowP = (double)((PXOList[1+ 5*(NextNo+1)].PermPVal));
						//TotalNoRecombinants = (PXOList[1+ 5*(NextNo+1)].Eventnumber);
						//return(TotalNoRecombinants);
						if (CPVal > 0 && CPVal < LowPX) {
							//if (PXOList[X + os].Beginning != PXOList[X + os].Ending) {
								LowPX = CPVal;
								Trace[0] = X;
								Trace[1] = Y;

							//}


						}
						else if (CPVal == LowPX && Trace[0]>X) {
							//if (PXOList[X + os].Beginning != PXOList[X + os].Ending) {
								LowPX = CPVal;
								Trace[0] = X;
								Trace[1] = Y;

							//}

						}
					}
				//}
			}
		}
		*LowP = LowPX;
		return(TotalNoRecombinants);
	}


	int MyMathFuncs::MakeTestPVs(int UBDS, unsigned char *DoneSeq, int NextNo, int UB, int UB2, short int *PCurrentXOver, XOVERDEFINE *PXOList, double *TestPVs) {
		int x, Y, os, os2, osTP, osPV, osDS;
		os = UB + 1;
		os2 = UB2 + 1;
		osDS = UBDS + 1;
		for (x = 0; x <= NextNo; x++) {
			osTP = x*os2;
			for (Y = 1; Y <= PCurrentXOver[x]; Y++) {
				osPV = Y*os;
				if (PXOList[x + osPV].Beginning != PXOList[x + osPV].Ending)
					TestPVs[Y + osTP] = PXOList[x + osPV].Probability;
				if (PXOList[x + osPV].OutsideFlag == 1)
					DoneSeq[x + Y*osDS] = 1;
			}
		}
		
				

			
			
		

		return(1);
	}



	int MyMathFuncs::AddjustCXO(int NextNo, int WinPP, double LowestProb,int UBDS1, int UBDS2, unsigned char *DoneSeq, int UBTD1, int UBTD2, unsigned char *TempDone, int *oRecombNo, int *RNum, int*RList, unsigned char *DoPairs, int UBTS, int *TraceSub, short int *tCurrentxover, int UBTXOL1, int UBTXOL2,  XOVERDEFINE *TempXOList,short int *PCurrentXOver, int UBPXO1, int UBPXO2, XOVERDEFINE *PXOList) {

		int x,Y,PCXO,DA,Mi,Ma, os, os2, WinPPY, os3, os4, os5, os6;
		os = UBPXO1 + 1;
		os3 = UBTXOL1 + 1;
		os4 = UBTD1 + 1;
		os5 = UBDS1 + 1;
		for (x = 0; x <= NextNo;x++){
			PCXO = PCurrentXOver[x];
			for (Y = 1; Y <= PCXO; Y++){
				if (Y <= UBPXO2){
					os2 = x + Y*os;
					DA = PXOList[os2].Daughter;
					Mi = PXOList[os2].MinorP;
					Ma = PXOList[os2].MajorP;
					if (DA > NextNo)
						DA = TraceSub[DA];
					
					if (Mi > UBTS)
						Mi = 0;
					
					if (Ma > UBTS)
						Ma = 0;
					
					if (Mi > NextNo)
						Mi = TraceSub[Mi];
					
					if (Ma > NextNo)
						Ma = TraceSub[Ma];
					
					WinPPY = MakePairsP(NextNo, DA, Ma, Mi, WinPP, RNum, RList, DoPairs);
					if (PXOList[os2].Probability <= LowestProb){
						if (WinPPY == RNum[WinPP]+ 1){
							if (tCurrentxover[DA] <= tCurrentxover[Mi] && tCurrentxover[DA] <= tCurrentxover[Ma]){
								tCurrentxover[x] = tCurrentxover[x] + 1;
								if (tCurrentxover[x] < UBTXOL2) {
									TempXOList[x + tCurrentxover[x] * os3] = PXOList[os2];
									if (tCurrentxover[x] <= UBTD2 && Y <= UBDS2)
										TempDone[x + tCurrentxover[x] * os4] = DoneSeq[x + Y*os5];
								}
								else
									tCurrentxover[x] = tCurrentxover[x] - 1;
							}
							else if (tCurrentxover[Mi] <= tCurrentxover[DA] && tCurrentxover[Mi] <= tCurrentxover[Ma]){
								tCurrentxover[Mi] = tCurrentxover[Mi] + 1;
								if (tCurrentxover[Mi] <= UBTXOL2){
									os6 = Mi + tCurrentxover[Mi] * os3;
									TempXOList[os6] = PXOList[os2];
									TempXOList[os6].Daughter = Mi;
									TempXOList[os6].MinorP = DA;
									if (tCurrentxover[Mi] <= UBTD2 && Y <= UBDS2)
										TempDone[Mi + tCurrentxover[Mi] * os4] = DoneSeq[x + Y*os5];
									
								}
								else
									tCurrentxover[Mi] = tCurrentxover[Mi] - 1;
								
							}
							else{
								tCurrentxover[Ma] = tCurrentxover[Ma] + 1;
								if (tCurrentxover[Ma] <= UBTXOL2) {
									os6 = Ma + tCurrentxover[Ma] * os3;
									TempXOList[os6] = PXOList[os2];
									TempXOList[os6].Daughter = Ma;
									TempXOList[os6].MajorP = DA;
									if (tCurrentxover[Ma] <= UBTD2 && Y <= UBDS2)
										TempDone[Ma + tCurrentxover[Ma] * os4] = DoneSeq[x + Y*os5];

								}
								else
									tCurrentxover[Ma] = tCurrentxover[Ma] - 1;
								
							}
							oRecombNo[100] = oRecombNo[100] + 1;
							oRecombNo[PXOList[os2].ProgramFlag] = oRecombNo[PXOList[os2].ProgramFlag] + 1;
						}
					}
					PCurrentXOver[x] = PCurrentXOver[x] - 1;
				}
			}
		}
		return(1);
	}

	int MyMathFuncs::CheckYannP(int SEN, int NextNo,int LenStrainSeq0, int BPos, int Epos, int *ISeqs, int UBSN, short int *SeqNum, unsigned char *IsPresent, int *TraceSub, int UBXH, int UBXHMi, int UBXHMa, unsigned char *ExtraHits, unsigned char *ExtraHitsMi, unsigned char *ExtraHitsMa, int *A, int *b) {
		int Y, Z, x, S1, s2, S3, os, n1, n2, N3, oeh1, oeh2, oeh3;
		S1 = ISeqs[0];
		s2 = ISeqs[1];
		S3 = ISeqs[2];
		os = UBSN + 1;
		oeh1 = UBXH + 1;
		oeh2 = UBXHMi + 1;
		oeh3 = UBXHMa + 1;


		if (BPos < Epos){
			for (x = 1; x <= BPos - 1; x++){
				n1 = SeqNum[x + S1*os];
				n2 = SeqNum[x + s2*os];
				N3 = SeqNum[x + S3*os];
        
				if (n1 > 46 && n2 > 46 && N3 > 46){
					if (n1 != n2 || n1 != N3){
						//check if this site is parsimony informative wrt recombination (i.e. recombinant and all its corecombinnats are considered as a single sequence as parents and their co-parenst are all considered single sequneces)
						IsPresent[n1] = 0;
						IsPresent[n2] = 0;
						IsPresent[N3] = 0;
						for (Y = 0; Y <= NextNo; Y++){
							Z = TraceSub[Y];
							if (Z <= UBXH && Z <= UBXHMi && Z <= UBXHMa){
								if (ExtraHits[Z + SEN*oeh1] == 0 && ExtraHitsMi[Z + SEN*oeh2] == 0 && ExtraHitsMa[Z + SEN*oeh3] == 0)
									IsPresent[SeqNum[x + Y*os]] = 1;
								
							}
						}
						if (IsPresent[n1] == 1 && IsPresent[n2] == 1 && IsPresent[N3] == 1){//this polymorphic site is informative about this ercombination event
							if (n1 == n2)
								A[0] = A[0] + 1;
							else if (n1 == N3)
								A[1] = A[1] + 1;
							else if (n2 == N3)
								A[2] = A[2] + 1;
							
						
						}
					}
				}
			}
			
			
			for (x = Epos+1; x <= LenStrainSeq0; x++) {
				n1 = SeqNum[x + S1*os];
				n2 = SeqNum[x + s2*os];
				N3 = SeqNum[x + S3*os];

				if (n1 > 46 && n2 > 46 && N3 > 46) {
					if (n1 != n2 || n1 != N3) {
						//check if this site is parsimony informative wrt recombination (i.e. recombinant and all its corecombinnats are considered as a single sequence as parents and their co-parenst are all considered single sequneces)
						IsPresent[n1] = 0;
						IsPresent[n2] = 0;
						IsPresent[N3] = 0;
						for (Y = 0; Y <= NextNo; Y++) {
							Z = TraceSub[Y];
							if (Z <= UBXH && Z <= UBXHMi && Z <= UBXHMa){
								if (ExtraHits[Z + SEN*oeh1] == 0 && ExtraHitsMi[Z + SEN*oeh2] == 0 && ExtraHitsMa[Z + SEN*oeh3] == 0) 
									IsPresent[SeqNum[x + Y*os]] = 1;

							}
						}
						if (IsPresent[n1] == 1 && IsPresent[n2] == 1 && IsPresent[N3] == 1) {//this polymorphic site is informative about this ercombination event
							if (n1 == n2)
								A[0] = A[0] + 1;
							else if (n1 == N3)
								A[1] = A[1] + 1;
							else if (n2 == N3)
								A[2] = A[2] + 1;


						}
					}
				}
			}
			
			for (x = BPos; x <= Epos; x++) {
				n1 = SeqNum[x + S1*os];
				n2 = SeqNum[x + s2*os];
				N3 = SeqNum[x + S3*os];

				if (n1 > 46 && n2 > 46 && N3 > 46) {
					if (n1 != n2 || n1 != N3) {
						//check if this site is parsimony informative wrt recombination (i.e. recombinant and all its corecombinnats are considered as a single sequence as parents and their co-parenst are all considered single sequneces)
						IsPresent[n1] = 0;
						IsPresent[n2] = 0;
						IsPresent[N3] = 0;
						for (Y = 0; Y <= NextNo; Y++) {
							Z = TraceSub[Y];
							if (Z <= UBXH && Z <= UBXHMi && Z <= UBXHMa) {
								if (ExtraHits[Z + SEN*oeh1] == 0 && ExtraHitsMi[Z + SEN*oeh2] == 0 && ExtraHitsMa[Z + SEN*oeh3] == 0)
									IsPresent[SeqNum[x + Y*os]] = 1;

							}
						}
						if (IsPresent[n1] == 1 && IsPresent[n2] == 1 && IsPresent[N3] == 1) {//this polymorphic site is informative about this ercombination event
							if (n1 == n2)
								b[0] = b[0] + 1;
							else if (n1 == N3)
								b[1] = b[1] + 1;
							else if (n2 == N3)
								b[2] = b[2] + 1;


						}
					}
				}
			}
		}
		else{
			for (x = Epos + 1; x <= BPos - 1;x++){
				n1 = SeqNum[x + S1*os];
				n2 = SeqNum[x + s2*os];
				N3 = SeqNum[x + S3*os];

				if (n1 > 46 && n2 > 46 && N3 > 46) {
					if (n1 != n2 || n1 != N3) {
						//check if this site is parsimony informative wrt recombination (i.e. recombinant and all its corecombinnats are considered as a single sequence as parents and their co-parenst are all considered single sequneces)
						IsPresent[n1] = 0;
						IsPresent[n2] = 0;
						IsPresent[N3] = 0;
						for (Y = 0; Y <= NextNo; Y++) {
							Z = TraceSub[Y];
							if (Z <= UBXH && Z <= UBXHMi && Z <= UBXHMa){
								if (ExtraHits[Z + SEN*oeh1] == 0 && ExtraHitsMi[Z + SEN*oeh2] == 0 && ExtraHitsMa[Z + SEN*oeh3] == 0) 
									IsPresent[SeqNum[x + Y*os]] = 1;

							}
						}
						if (IsPresent[n1] == 1 && IsPresent[n2] == 1 && IsPresent[N3] == 1) {//this polymorphic site is informative about this ercombination event
							if (n1 == n2)
								A[0] = A[0] + 1;
							else if (n1 == N3)
								A[1] = A[1] + 1;
							else if (n2 == N3)
								A[2] = A[2] + 1;


						}
					}
				}
			}
			
			for (x = 1; x <= Epos; x++) {
				n1 = SeqNum[x + S1*os];
				n2 = SeqNum[x + s2*os];
				N3 = SeqNum[x + S3*os];

				if (n1 > 46 && n2 > 46 && N3 > 46) {
					if (n1 != n2 || n1 != N3) {
						//check if this site is parsimony informative wrt recombination (i.e. recombinant and all its corecombinnats are considered as a single sequence as parents and their co-parenst are all considered single sequneces)
						IsPresent[n1] = 0;
						IsPresent[n2] = 0;
						IsPresent[N3] = 0;
						for (Y = 0; Y <= NextNo; Y++) {
							Z = TraceSub[Y];
							if (Z <= UBXH && Z <= UBXHMi && Z <= UBXHMa) {
								if (ExtraHits[Z + SEN*oeh1] == 0 && ExtraHitsMi[Z + SEN*oeh2] == 0 && ExtraHitsMa[Z + SEN*oeh3] == 0)
									IsPresent[SeqNum[x + Y*os]] = 1;

							}
						}
						if (IsPresent[n1] == 1 && IsPresent[n2] == 1 && IsPresent[N3] == 1) {//this polymorphic site is informative about this ercombination event
							if (n1 == n2)
								b[0] = b[0] + 1;
							else if (n1 == N3)
								b[1] = b[1] + 1;
							else if (n2 == N3)
								b[2] = b[2] + 1;


						}
					}
				}
			}
			
			for (x = BPos; x <= LenStrainSeq0; x++) {
				n1 = SeqNum[x + S1*os];
				n2 = SeqNum[x + s2*os];
				N3 = SeqNum[x + S3*os];

				if (n1 > 46 && n2 > 46 && N3 > 46) {
					if (n1 != n2 || n1 != N3) {
						//check if this site is parsimony informative wrt recombination (i.e. recombinant and all its corecombinnats are considered as a single sequence as parents and their co-parenst are all considered single sequneces)
						IsPresent[n1] = 0;
						IsPresent[n2] = 0;
						IsPresent[N3] = 0;
						for (Y = 0; Y <= NextNo; Y++) {
							Z = TraceSub[Y];
							if (Z <= UBXH && Z <= UBXHMi && Z <= UBXHMa) {
								if (ExtraHits[Z + SEN*oeh1] == 0 && ExtraHitsMi[Z + SEN*oeh2] == 0 && ExtraHitsMa[Z + SEN*oeh3] == 0)
									IsPresent[SeqNum[x + Y*os]] = 1;

							}
						}
						if (IsPresent[n1] == 1 && IsPresent[n2] == 1 && IsPresent[N3] == 1) {//this polymorphic site is informative about this ercombination event
							if (n1 == n2)
								b[0] = b[0] + 1;
							else if (n1 == N3)
								b[1] = b[1] + 1;
							else if (n2 == N3)
								b[2] = b[2] + 1;


						}
					}
				}
			}
		}

		return(1);
	}

	int MyMathFuncs::MarkOutsides(int UBDS, unsigned char *DoneSeq, int NextNo, int UB,  short int *PCurrentXOver, XOVERDEFINE *PXOList) {
		int x, Y, os, osDS;
		os = UB + 1;
		osDS = UBDS + 1;
		for (x = 0; x <= NextNo; x++) {
			for (Y = 1; Y <= PCurrentXOver[x]; Y++) {
				
				if (DoneSeq[x + Y*osDS] == 1)
					PXOList[x + Y*os].OutsideFlag = 1;
					
			}
		}
		return(1);
	}

	int  MyMathFuncs::UpdatePlotsCP2(int UBAD,float ff, HDC Pict, int LSeq, short int P1, short int P2, short int P3, short int P4, int StepSize, float XFactor, float oDMax, float oPMax, int MaxHits, int *Decompress, float *PDistPlt, float *ProbPlt, int *HitPlt, float *ll1, float *ll2, float *ll3) {
		int  os, p;
		float a, b, c,hc, mh;


		a = 0;
		b = 0;
		c = 0;
		p = 0;

				omp_set_num_threads(3);
		#pragma omp parallel 
				{
		#pragma omp sections private (p)
					{
		#pragma omp section 
						{
		if (oDMax > 0) {
			os = P2 - P1 - 10;
			p = 0;
			for (a = 1; a <= UBAD; a = a + StepSize) {
				p++;
				//LineTo(Pict, 30 + Decompress[a] * XFactor + XFactor, P2 - 5 - (PDistPlt[a] / oDMax)* os);
				ll1[p*2] = 30 + Decompress[(int)(a / ff)] * XFactor + XFactor;
				ll1[1+p*2] = P2 - 5 - (PDistPlt[(int)(a)] / oDMax)* os;
			}
			p++;
			ll1[p*2] = 30 + Decompress[(int)(UBAD / ff)] * XFactor + XFactor;
			ll1[1 + p * 2] = P2 - 5 - (PDistPlt[UBAD] / oDMax)* os;
			p++;
			ll1[p * 2] = 30 + Decompress[(int)(UBAD/ ff)] * XFactor + XFactor;
			ll1[1 + p * 2] = P2 - 5 ;
			p++;
			ll1[p * 2] = 30 + Decompress[(int)(1 / ff)] * XFactor + XFactor;
			ll1[1 + p * 2] = P2 - 5;
			p++;
			ll1[p * 2] = 30 + Decompress[(int)(1 / ff)] * XFactor + XFactor;
			ll1[1 + p * 2] = P2 - 5 - (PDistPlt[1] / oDMax)* os;
		}
						}
		#pragma omp section
						{
		if (oPMax > 0) {
			os = P3 - P2 - 10;
			p = 0;
			for (b = 1; b <= UBAD; b = b + StepSize) {
				p++;
				//LineTo(Pict, 30 + Decompress[b] * XFactor + XFactor, P3 - 5 - (ProbPlt[b] / oPMax) * os);
				
				if (ProbPlt[(int)(b)] > oPMax)
					oPMax = ProbPlt[(int)(b)];

				ll2[p*2] = 30 + Decompress[(int)(b / ff)] * XFactor + XFactor;
				ll2[1 + p * 2] = P3 - 5 - (ProbPlt[(int)(b)] / oPMax) * os;
			}
			p--;
			ll2[p * 2] = 30 + Decompress[(int)(UBAD / ff)] * XFactor + XFactor;
			ll2[1 + p * 2] = P3 - 5 - (ProbPlt[UBAD] / oPMax)* os;
			p++;
			ll2[p * 2] = 30 + Decompress[(int)(UBAD / ff)] * XFactor + XFactor;
			ll2[1 + p * 2] = P3 - 5 ;
			p++;
			ll2[p * 2] = 30 + Decompress[(int)(1 / ff)] * XFactor + XFactor;
			ll2[1 + p * 2] = P3 - 5 ;
			p++;
			ll2[p * 2] = 30 + Decompress[(int)(1 / ff)] * XFactor + XFactor;
			ll2[1 + p * 2] = P3 - 5 - (ProbPlt[1] / oPMax)* os;
								}
							}
			#pragma omp section
							{
			mh = (float)(MaxHits);
			if (MaxHits > 0) {
				os = P4 - P3 - 10;
				p = 0;
				for (c = 1; c <= UBAD; c = c + StepSize) {
					p++;
					hc = (float)(HitPlt[(int)(c)]);
					//LineTo(Pict, 30 + Decompress[c] * XFactor + XFactor, P4 - 5 - (hc / mh) * os);
					ll3[p*2] = 30 + Decompress[(int)(c / ff)] * XFactor + XFactor;
					ll3[1+p*2] = P4 - 5 - (hc / mh) * os;
				}
				p++;
				hc = (float)(HitPlt[UBAD]);
				ll3[p * 2] = 30 + Decompress[(int)(UBAD / ff)] * XFactor + XFactor;
				ll3[1 + p * 2] = P4 - 5 - (hc / mh)* os;
				p++;
				ll3[p * 2] = 30 + Decompress[(int)(UBAD / ff)] * XFactor + XFactor;
				ll3[1 + p * 2] = P4 - 5;
				p++;
				ll3[p * 2] = 30 + Decompress[(int)(1 / ff)] * XFactor + XFactor;
				ll3[1 + p * 2] = P4 - 5;
				p++;
				hc = (float)(HitPlt[1]);
				ll3[p * 2] = 30 + Decompress[(int)(1 / ff)] * XFactor + XFactor;
				ll3[1 + p * 2] = P4 - 5 - (hc / mh)* os;
			}
		}
					}
				}
		return(1);
	}
	double  MyMathFuncs::UpdateDonePVCO(double NPVal, double LPV,int Prg,int s1, int SIP, int UBXOL1, int UBDPV, short int *CurrentXOver, XOVERDEFINE *XoverList, double *DonePVCO) {
		int x,os, os2, os3;
		os = UBXOL1 + 1;
		os3 = Prg+(UBDPV + 1)*s1;
		for (x = 1; x <= CurrentXOver[s1]; x++){

			if (x != SIP) {
				os2 = x*os;
				if (XoverList[s1 + os2].ProgramFlag == Prg) {

					if (XoverList[s1 + os2].Probability > LPV) {
						LPV = XoverList[s1 + os2].Probability;
						DonePVCO[os3] = x;
					}
				}
			}
			else{
				if (NPVal > LPV) {
					LPV = NPVal;
					DonePVCO[os3] = x;
				}
					
			}
			
		}
		return(LPV);
	}

	int  MyMathFuncs::GetFragsP(short int CircularFlag, int LenXoverSeq, int LSeq, int maxcount, char *SubSeq, int *FragSt, int *FragEn, int *FragScore, int *FragCount) {

		int SX, A, X, Y, Z, GoOn, os, os2, os3, FE, FS, RR, FSc, FC6, os4, os5, target;

		GoOn = 0;
		os = maxcount + 1;
		RR = LenXoverSeq + 1;
		FC6 = 0;
		os4 = LSeq + 1;
		os5 = maxcount + 1;

		target = 6 * os;
		X = 0;
		Y = 0;
		Z = 0;



//
//		omp_set_num_threads(3);
//#pragma omp parallel 
//		{
//#pragma omp sections private (X)
//			{
//#pragma omp section
//				{
					for (X = 0; X <= target; X++)
						FragScore[X] = 0;

//				}
//#pragma omp section
//				{
					for (X = 0; X <= target; X++)

						FragSt[X] = 0;
//				}
//
//#pragma omp section
//				{
					for (X = 0; X <= target; X++)

						FragEn[X] = 0;
				/*}
			}
			
		}*/

		//target = 
		/*for (Y = 0; Y <= 6; Y++){
		osx = Y*os;
		for (X = 0; X < RR; X++)
		FragScore[X+osx] = 0;
		}
		for (Y = 0; Y <= 6; Y++){
		osx = Y*os;
		for (X = 0; X <= RR; X++){
		FragSt[X+osx] = 0;
		FragEn[X+osx] = 0;

		}
		}*/


		for (X = 1; X < RR; X++) {
			for (Y = 0; Y < 3; Y++) {
				if (SubSeq[X + Y*os4] == 1)
					break;
			}
			if (Y == 3)
				Y = 6;

			os2 = Y*os4;
			os = FragCount[Y] + Y*os5;
			FragSt[os] = X;
			FS = X;
			SX = X;
			X++;
			if (CircularFlag == 0) {
				while (SubSeq[os2 + X] == 1) {
					X++;
					if (X > LenXoverSeq)
						break;
				}
			}
			else {

				while (SubSeq[os2 + X] == 1) {
					X++;
					if (X > LenXoverSeq) {
						X = 1;
						GoOn = 1;
						break;

					}
				}
				if (X == 1) {
					while (SubSeq[os2 + X] == 1 && X != SX) {
						X++;

					}
				}

				if (X == SX)
					return(0);
			}



			X--;

			FragEn[os] = X;
			FE = X;
			if (FE > FS)
				FSc = FE - FS + 1;
			else if (FE < FS)
				FSc = FE - FS + RR;
			else
				FSc = 1;

			FragScore[os] = FSc;


			/*for (Z = 0; Z < 3; Z++){
			if (Y != Z){
			os3 = FragCount[Z] + Z*os5;
			FragSt[os3] = FS;
			FragEn[os3] = FE;
			FragScore[os3] = -FSc;
			FragCount[Z] = FragCount[Z] + 1;
			}
			}*/


			if (Y != 0) {
				os3 = FragCount[0];
				FragSt[os3] = FS;
				FragEn[os3] = FE;
				FragScore[os3] = -FSc;
				FragCount[0] = FragCount[0] + 1;
			}

			if (Y != 1) {
				os3 = FragCount[1] + os5;
				FragSt[os3] = FS;
				FragEn[os3] = FE;
				FragScore[os3] = -FSc;
				FragCount[1] = FragCount[1] + 1;
			}
			if (Y != 2) {
				os3 = FragCount[2] + 2 * os5;
				FragSt[os3] = FS;
				FragEn[os3] = FE;
				FragScore[os3] = -FSc;
				FragCount[2] = FragCount[2] + 1;
			}


			if (Y == 0)
				Z = 5;
			else if (Y == 1)
				Z = 4;
			else if (Y == 2)
				Z = 3;
			else if (Y == 6)
				Z = 6;

			if (Z < 6) {
				os3 = FragCount[Z] + Z*os5;
				if (FragCount[Z] > 0) {
					if (FragScore[os3 - 1] > 0) {
						FragEn[os3 - 1] = FE;
						FragScore[os3 - 1] = FragScore[os3 - 1] + FSc;
					}
					else {
						FragSt[os3] = FS;
						FragEn[os3] = FE;
						FragScore[os3] = FSc;
						FragCount[Z] = FragCount[Z] + 1;
					}
				}
				else {
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = FSc;
					FragCount[Z] = FragCount[Z] + 1;
				}

				for (A = 3; A < 6; A++) {
					if (A != Z) {
						os3 = FragCount[A] + A*os5;
						FragSt[os3] = FS;
						FragEn[os3] = FE;
						FragScore[os3] = -FSc;
						FragCount[A] = FragCount[A] + 1;
					}
				}

				/*if (Z != 3){
				os3 = FragCount[3] + 3*os5;
				FragSt[os3] = FS;
				FragEn[os3] = FE;
				FragScore[os3] = -FSc;
				FragCount[3] = FragCount[3] + 1;
				}
				if (Z != 4){
				os3 = FragCount[4] + 4*os5;
				FragSt[os3] = FS;
				FragEn[os3] = FE;
				FragScore[os3] = -FSc;
				FragCount[4] = FragCount[4] + 1;
				}
				if (Z != 5){
				os3 = FragCount[5] + 5*os5;
				FragSt[os3] = FS;
				FragEn[os3] = FE;
				FragScore[os3] = -FSc;
				FragCount[5] = FragCount[5] + 1;
				}*/


			}
			else {
				if (FragCount[3] > 0) {

					for (A = 3; A < 6; A++) {
						os3 = FragCount[A] + A*os5;
						if (FragScore[os3 - 1] > 0) {

							FragEn[os3 - 1] = FE;
							FragScore[os3 - 1] = FragScore[os3 - 1] + FSc;
						}
						else {
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = FSc;
							FragCount[A] = FragCount[A] + 1;
						}
					}//Next A
				}
				else {

					for (A = 3; A < 6; A++) {
						os3 = FragCount[A] + A*os5;
						FragSt[os3] = FS;
						FragEn[os3] = FE;
						FragScore[os3] = FSc;
						FragCount[A] = FragCount[A] + 1;
					}//Next A
				}
			}

			/*if (Y <= 2)
				FragCount[Y] = FragCount[Y] + 1;*/

			if (Y <= 2) 
					FragCount[Y] = FragCount[Y] + 1;
			

			for (A = 0; A < 6; A++) {
				if (FragCount[A] >= maxcount-1) {
					GoOn = 1;
					break;
				}


			}


			

			FC6++;
			if (GoOn == 1 || FC6 >= maxcount-1)
							break;

			

		}//next X
		FragCount[6] = FC6;


		for (X = 0; X <= 6; X++) {
			os = FragCount[X] + X*os5;
			FragSt[os] = LenXoverSeq;
			FragEn[os] = LenXoverSeq;
		}
		//omp_set_num_threads(2);
		return(1);
	}

	int  MyMathFuncs::MakeWindowSizeP(int BEP, int ENP, int *CriticalDiff, int LenXoverSeq, double CWinFract, int CWinSize, int *HWindowWidth, int lHWindowWidth, short int CProportionFlag) {

	//Public Sub MakeWindowSize(BEP, ENP, CriticalDiff As Long, LenXoverSeq As Long, CWinFract As Double, CWinSize As Long, HWindowWidth As Long, lHWindowWidth As Long, GoOn As Byte, FindallFlag As Byte, CProportionFlag As Integer, XPosDiff() As Long)
		float WindowWidth;
		int EN, BE;
		
		
		if (CProportionFlag == 0){
			WindowWidth = CWinSize;
			*HWindowWidth = (int)((WindowWidth / 2.0)+0.51);
		}
		else{
			if ((int)(CWinFract * LenXoverSeq) > 20 && (int)(CWinFract * LenXoverSeq) < (LenXoverSeq / 1.5))
				WindowWidth = ((int)((CWinFract * LenXoverSeq) / 2) * 2 - 2);
			else if ((int)(CWinFract * LenXoverSeq) <= 20){
    
				if (LenXoverSeq > 15)
					WindowWidth = 20;
				else
					return(0);
			}
			else if (int(CWinFract * LenXoverSeq) >= (LenXoverSeq / 1.5)){
    
				if (LenXoverSeq / 1.5 > 10)
					WindowWidth = ((int)((LenXoverSeq / 1.5) / 2) * 2 - 2);
				else
					return(0);
			}
			*HWindowWidth = (int)((WindowWidth / 2) + 0.51);
        
			*CriticalDiff = 2;
			if (CWinSize != *HWindowWidth * 2 && CProportionFlag == 0)
				CWinSize = *HWindowWidth * 2;
			
		}
		
		if(*HWindowWidth * 2 > LenXoverSeq)
			*HWindowWidth = (int)(((LenXoverSeq * 0.75) / 2)+0.51) - 1;
		if (*HWindowWidth <= *CriticalDiff)
			*HWindowWidth = (int)((LenXoverSeq / 2)+0.51) - 1;
		if (*HWindowWidth < 6) {
			*HWindowWidth = lHWindowWidth;
			return(0);
		}

		
		return(1);
	}
	int  MyMathFuncs::AddPVal(int Prog, double *mtP,int HWindowWidth,int LenXoverSeq,int MCCorrection,short int MCFlag, double MChi, double LowestProb) {

		//Public Sub AddPVal(MCFlag As Integer, MChi As Double, MCCorrection As Long, LenXoverSeq As Long, HWindowWidth As Long, mtP() As Double, Prog As Long)
		if (MCFlag == 0) {
			if (ChiPVal(MChi) * MCCorrection * ((float)(LenXoverSeq) / (float)(HWindowWidth)) <= mtP[Prog])
				mtP[Prog] = ChiPVal(MChi) * MCCorrection * ((float)(LenXoverSeq) / (float)(HWindowWidth));

		}
		else {
			if (ChiPVal(MChi) * (LenXoverSeq / (HWindowWidth)) <= mtP[Prog])
				mtP[Prog] = ChiPVal(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth));

		}
		return(0);
	}

	int  MyMathFuncs::MakeTWinP(unsigned char FindallFlag, int HWindowWidth, int *TWin, int LenXoverSeq) {

	//Public Sub MakeTWin(FindallFlag, HWindowWidth As Long, TWin As Long, LenXoverSeq As Long)
		if (FindallFlag == 0){
			*TWin = int(((float)(HWindowWidth) / 4.0) + 0.51);
			if (*TWin < 6)
				*TWin = 6;
			
			if (*TWin > HWindowWidth)
				*TWin = HWindowWidth;
			
			if (*TWin > int((float)(LenXoverSeq) / 2.0))
				*TWin = int((float)(LenXoverSeq) / 2.0);
			
		}
		else
			*TWin = HWindowWidth;
		
		return(1);
	}

	int  MyMathFuncs::DestroyPeakP(int MaxY, int LS, int RO, int LO, int LenXoverSeq, double *LOT, double *SmoothChi, double *ChiVals) {
		//smootchchi - ls,2
		int Circuit1, EraseAll, X;
		EraseAll = 0;
		if (RO < 1) {
			RO = LenXoverSeq - RO;
			if (RO < 2)
				return(0);
		}
		Circuit1 = 0;
		if (RO > LenXoverSeq * 2)
			RO = RO - LenXoverSeq * 2;
		if (RO > LenXoverSeq) {
			RO = RO - LenXoverSeq;
			if (RO >= LenXoverSeq)
				return(0);
			Circuit1 = 1;
		}

		else
			Circuit1 = 0;

		if (RO > LenXoverSeq - 2) {
			if (RO == LenXoverSeq) {
				LOT[0] = SmoothChi[1 + MaxY*(LS + 1)];
				LOT[1] = SmoothChi[2 + MaxY*(LS + 1)];
			}
			else if (RO == LenXoverSeq - 1) {
				LOT[0] = SmoothChi[LenXoverSeq + MaxY*(LS + 1)];
				LOT[1] = SmoothChi[1 + MaxY*(LS + 1)];
			}
		}
		else {
			LOT[0] = SmoothChi[RO + 1 + MaxY*(LS + 1)];
			LOT[1] = SmoothChi[RO + 2 + MaxY*(LS + 1)];
		}
		while (SmoothChi[RO + MaxY*(LS + 1)] > 0 && (SmoothChi[RO + MaxY*(LS + 1)] >= LOT[0] || SmoothChi[RO + MaxY*(LS + 1)] >= LOT[1])) {
			RO++;
			if (RO > LenXoverSeq) {
				RO = 1;
				if (Circuit1 == 1) {
					EraseAll = 1;
					break;
				}
				else
					Circuit1 = 1;
			}

			if (RO > LenXoverSeq - 2) {
				if (RO == LenXoverSeq) {
					LOT[0] = SmoothChi[1 + MaxY*(LS + 1)];
					LOT[1] = SmoothChi[2 + MaxY*(LS + 1)];
				}
				else if (RO == LenXoverSeq - 1) {
					LOT[0] = SmoothChi[LenXoverSeq + MaxY*(LS + 1)];
					LOT[1] = SmoothChi[1 + MaxY*(LS + 1)];
				}
			}
			else {
				LOT[0] = SmoothChi[RO + 1 + MaxY*(LS + 1)];
				LOT[1] = SmoothChi[RO + 2 + MaxY*(LS + 1)];
			}
		}
		if (EraseAll == 0) {
			if (LO < 1) {
				LO = LenXoverSeq - RO;
				if (LO < 2)
					return(0);

			}
			Circuit1 = 0;
			if (LO > LenXoverSeq) {
				LO = LO - LenXoverSeq;
				if (LO >= LenXoverSeq)
					return(0);
				Circuit1 = 1;
			}
			else
				Circuit1 = 0;

			if (LO < 3) {
				if (LO == 1) {
					LOT[0] = SmoothChi[LenXoverSeq + MaxY*(LS + 1)];
					LOT[1] = SmoothChi[LenXoverSeq - 1 + MaxY*(LS + 1)];
				}
				else if (LO == 2) {
					LOT[0] = SmoothChi[1 + MaxY*(LS + 1)];
					LOT[1] = SmoothChi[LenXoverSeq + MaxY*(LS + 1)];
				}
			}
			else {
				LOT[0] = SmoothChi[LO - 1 + MaxY*(LS + 1)];
				LOT[1] = SmoothChi[LO - 2 + MaxY*(LS + 1)];
			}
			while (SmoothChi[LO + MaxY*(LS + 1)] > 0 && (SmoothChi[LO + MaxY*(LS + 1)] >= LOT[0] || SmoothChi[LO + MaxY*(LS + 1)] >= LOT[1])) {
				LO--;
				if (LO < 1) {
					LO = LenXoverSeq;
					if (Circuit1 == 1) {
						EraseAll = 1;
						break;
					}
					else
						Circuit1 = 1;

				}
				if (LO < 3) {
					if (LO == 1) {
						LOT[0] = SmoothChi[LenXoverSeq + MaxY*(LS + 1)];
						LOT[1] = SmoothChi[LenXoverSeq - 1 + MaxY*(LS + 1)];
					}
					else if (LO == 2) {
						LOT[0] = SmoothChi[1 + MaxY*(LS + 1)];
						LOT[1] = SmoothChi[LenXoverSeq + MaxY*(LS + 1)];
					}
				}
				else {
					LOT[0] = SmoothChi[LO - 1 + MaxY*(LS + 1)];
					LOT[1] = SmoothChi[LO - 2 + MaxY*(LS + 1)];
				}
			}
		}
		//chivals ls,2
		if (EraseAll == 1) {
			for (X = 0; X <= LenXoverSeq; X++)
				ChiVals[X + MaxY*(LS + 1)] = 0;
		}
		else {
			if (LO < RO) {
				for (X = LO; X <= RO; X++)
					ChiVals[X + MaxY*(LS + 1)] = 0;
			}

			else {
				for (X = 0; X <= RO; X++)
					ChiVals[X + MaxY*(LS + 1)] = 0;

				for (X = LO; X <= LenXoverSeq; X++)
					ChiVals[X + MaxY*(LS + 1)] = 0;

			}
		}
		return(1);
	}


	int  MyMathFuncs::FastRecCheckMC(int SEN, int LongWindedFlag, double *BQPV, int EarlyBale, double UCTHresh, short int MCFlag, short int ShortOutFlag, int MCCorrection, double LowestProb, short int CircularFlag, int NextNo, int MaxABWin, int HWindowWidthX, int lHWindowWidth, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainSeq0, int CriticalDiffX, unsigned char FindallFlag,  int UBFSSMC, int UBCS, int Seq1, int Seq2, int Seq3, unsigned char *CS, unsigned char *FSSMC, short int *SeqNum, int UBWS, unsigned char *Scores, int *Winscores, int *XDiffPos, double *Chivals,  int *BanWin, unsigned char *MDMap, float *ChiTable2, int *Chimap, double *mtP, double *SmoothChi) {
		int tb,te,x, Redox, MaxX, maxz, TWin, MaxFailCount, WinWin, A, C, pMaxX, TopL, TopR, TopLO, TopRO, LO, RO;
		double mPrb,MChi, xMPV, MPV;
		double *LOT;
		short int MaxY;
		LOT = (double *)calloc(2, sizeof(double));
		tb = 0;
		te = 0;
		int	HWindowWidth, CriticalDiff, Dummy;
		int WasteOfTime, LenXoverSeq,  GoOn;

		GoOn = 0;
		HWindowWidth = HWindowWidthX;
		CriticalDiff = CriticalDiffX;

		LenXoverSeq = FindSubSeqMCPB2(UBFSSMC, UBCS, NextNo, Seq1, Seq2, Seq3, CS, FSSMC, XDiffPos);

		if (LenXoverSeq < CriticalDiff * 2 || LenXoverSeq < 7) {
			free(LOT);
			return(0);
		}
		

		GoOn = MakeWindowSizeP(tb, te, &CriticalDiff, LenXoverSeq, MCWinFract, MCWinSize, &HWindowWidth, lHWindowWidth,  MCProportionFlag);
		//                (int BEP, int ENP, int *CriticalDiff, int LenXoverSeq, double CWinFract, int CWinSize, int *HWindowWidth, int lHWindowWidth, unsigned char *FindallFlag, short int CProportionFlag, int *XPosDiff)
		if (GoOn == 0) {
			free(LOT);
			return(0);
		}



		Dummy = WinScoreCalcP(UBWS, CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0 + 1, Seq1, Seq2, Seq3, Scores, XDiffPos, SeqNum, Winscores);

		

		//return(0);
		if (FindallFlag == 0 && (CircularFlag == 0 )){
			/*if (SEN > 0) {

				if (SEN == 1)
					Dummy = ClearDeleteArray(LenStrainSeq0, BanWin);

				Dummy = MakeBanWinP(LenStrainSeq0, Seq1, Seq2, Seq3, HWindowWidth, LenStrainSeq0, LenXoverSeq, BanWin, MDMap, MissingData, XPosDiff, XDiffPos)
			}*/

        
			Dummy = ClearDeleteArrayB(LenStrainSeq0, MDMap);
			Dummy = ClearDeleteArray(LenStrainSeq0 + HWindowWidth * 2, BanWin);
        
    
			if (CircularFlag == 0){
				//Ban windows spanning the ends of the alignment
				MDMap[1] = 1;
				MDMap[LenXoverSeq] = 1;
				for (x = (LenXoverSeq - HWindowWidth + 2); x <=LenXoverSeq; x++)
					BanWin[x] = 1;
			}

			if (HWindowWidth > 4 && HWindowWidth <= MaxABWin)

				MChi = CalcChiVals4P3(UBWS, CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0, Winscores, Chivals, BanWin, &ChiTable2[Chimap[HWindowWidth]]);
			else
				MChi = CalcChiVals4P(UBWS, CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0, Winscores, Chivals, BanWin);

		}
		else {
			if (SEN == 1) {
				Dummy = ClearDeleteArrayB(LenStrainSeq0, MDMap);
				Dummy = ClearDeleteArray(LenStrainSeq0 + HWindowWidth * 2, BanWin);
			}
			//ReDim MDMap(LenXOverSeq), BanWin(LenXOverSeq + HWindowWidth * 2)
			MChi = CalcChiValsP(UBWS, CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0, Winscores, Chivals);

		}
		
		
		
		//End If
		//return(0);
		if (MCProportionFlag == 0){
			if ((ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 3) > LowestProb) {
				free(LOT);
				Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
				return(0);
			}
		
		}
		else {
			if (((ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 6) * MCCorrection) > LowestProb) {
				free(LOT);
				Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
				return(0);
			}
		}


		

		if (ShortOutFlag == 3)
			Dummy = AddPVal(3, mtP, HWindowWidth, LenXoverSeq, MCCorrection, MCFlag, MChi, LowestProb);
			//int  MyMathFuncs::AddPVal(int Prog, double *mtP,int HWindowWidth,int LenXoverSeq,int MCCorrection,short int MCFlag, double MChi, double LowestProb) {



			//Smooth to find peaks (uses a window of positions)


			//test if this p version is quicker
			//SmoothChiValsP LenXOverSeq, Len(StrainSeq(0)), Chivals(0, 0), SmoothChi(0, 0)

			//The p - version is slower

		Dummy = SmoothChiValsP(LenXoverSeq, LenStrainSeq0, Chivals, SmoothChi);
		//return(0);
		MChi = 0;

		
		WasteOfTime = 0; //Keeps track of how many cycles have been wasted looking for insignificant hit

	//Now find recombination events


		Redox = 0;

		


		while (Redox <=100){

			Redox ++;
    
			Dummy = FindMChiP(LenStrainSeq0, LenXoverSeq, &MaxX, &MaxY, &MChi, Chivals);

			

    
			if (MaxX == -1 || MaxY == -1 || Redox > 100) {
				free(LOT);
				Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
				Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
				return(0);
			}
    
			//The p-version of this seems slower than the old one
			mPrb = ChiPVal2P(MChi);
			if (MCProportionFlag == 0) {

				if ((mPrb * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 3.0) > LowestProb || mPrb == 1) {
					free(LOT);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(0);
				}
			
			}
			else{
				if (((mPrb * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 6.0) * MCCorrection) > LowestProb || mPrb == 1) {
					free(LOT);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(0);
				}
				
			}

			

			if (mPrb < LowestProb){
				maxz = MaxX;
        
				if (MaxX == 0)
					MaxX = 1;
		
				Dummy = MakeTWinP(FindallFlag, HWindowWidth, &TWin, LenXoverSeq);
				//Call MakeTWin(FindallFlag, HWindowWidth, TWin, LenXoverSeq)
        
        
				MaxFailCount = HWindowWidth * 2;
        
				if (MaxFailCount > int((LenXoverSeq - TWin * 2) / 2))
					MaxFailCount = int((LenXoverSeq - TWin * 2) / 2);
		
				if (MaxFailCount == 0)
					MaxFailCount = 1;

				WinWin = HWindowWidth; 
				A = 0; 
				C = 0;
				Dummy = GetACP(LenXoverSeq, LenStrainSeq0, MaxY, MaxX, TWin, &A, &C, Scores);
				//Now find the other breakpoint.
				pMaxX = MaxX;
				if (MaxX < 1){
					if (CircularFlag == 1)
						MaxX = LenXoverSeq + MaxX;
					else
						MaxX = 1;
				
				}
				else if (MaxX > LenXoverSeq){
					if (CircularFlag == 1)
						MaxX = MaxX - LenXoverSeq;
					else
						MaxX = LenXoverSeq - 1;
				
				}
        
				if (TWin >= HWindowWidth) {
					TopL = A;
					TopR = C;
				}
				else {
					TopL = 0;
					TopR = 0;
				}
        
				TWin++;
        
				//Store best values from the initial screen
        
				MPV = mPrb * ((float)(LenXoverSeq) / (float)(HWindowWidth));
        
				WinWin = HWindowWidth;
				TopLO = MaxX - HWindowWidth + 1;
				TopRO = MaxX + HWindowWidth;
        
				LO = (MaxX - TWin + 1);
				if (LO < 0)
					LO = LenXoverSeq + LO;
				RO = (MaxX + TWin);
				if (RO > LenXoverSeq * 2)
					RO = RO - LenXoverSeq * 2;
				if (RO >= LenXoverSeq){
					if (MDMap[LenXoverSeq] == 0)
						RO = RO - LenXoverSeq;
					else{
						TWin = TWin - 1;
						RO = LenXoverSeq;
						LO = (MaxX - TWin + 1);
						if (LO < 0)
							LO = LenXoverSeq + LO;
                        
					}
				}
				
				if (FindallFlag == 0){
            
            
					if (LongWindedFlag == 1 && (CircularFlag == 0 || SEN > 1)){
						//better use of floats might help here - also the p-version gives a slightly different result to the non-p version
						//These routines give slightly different results to the vc5 versions
                
						Dummy = GrowMChiWin2P2(MaxABWin, LO, RO, LenXoverSeq, HWindowWidth, TWin, MaxY, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, MDMap, ChiTable2, Chimap);
						/*if (WinWin > 1000)
							WinWin = 1000;*/

						if (WinWin < HWindowWidth)
							MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(WinWin));
						else
							MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth));
					
					}
					else{
                        
						Dummy = GrowMChiWinP2(MaxABWin, LO, RO, LenXoverSeq, HWindowWidth, TWin, MaxY, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, ChiTable2, Chimap);
						/*if (WinWin > 1000)
							WinWin = 1000;*/
						if (WinWin < HWindowWidth)
							MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(WinWin));
						else{
							if (MChi < 20000)
								MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth));
							else
								MPV = pow(10,-200);
						
						}
					}
           
				}
				else {
					Dummy = GrowMChiWinP2(MaxABWin, LO, RO, LenXoverSeq, HWindowWidth, TWin, MaxY, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, ChiTable2, Chimap);
					/*if (WinWin > 1000)
						WinWin = 1000;*/
					if (WinWin < HWindowWidth)
						MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(WinWin));
					else {
						if (MChi < 20000)
							MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth));
						else
							MPV = pow(10, -200);

					}


				}

				
				

				xMPV = MPV;
				if (MCProportionFlag == 0)
					xMPV = xMPV * 3;
				else
					xMPV = xMPV * 3;
			
				if (xMPV < *BQPV) {

					*BQPV = xMPV;
					

				}
        
				if (xMPV < UCTHresh && EarlyBale == 1) {
					free(LOT);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(1);
				}
				if (MCFlag == 0)
					xMPV = xMPV * MCCorrection;
        
			
				MPV = xMPV;
			
        
				if (ShortOutFlag == 3){
					if (MPV <= mtP[3])
						mtP[3] = MPV;
				
				}
        
				LO = MaxX - WinWin;
				RO = MaxX + WinWin - 1;
				if (RO > LenXoverSeq * 2)
					RO = RO - LenXoverSeq * 2;
				if (RO > LenXoverSeq){
					if (MDMap[RO - LenXoverSeq] == 0)
						RO = RO - LenXoverSeq;
					else
						RO = LenXoverSeq;
				}
				if (LO < 1)
					LO = LO + LenXoverSeq;
			
				if (MPV < LowestProb) {

					//*BQPV = MPV;
					free(LOT);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(1);

				}
				else{
					WasteOfTime = WasteOfTime + 1;
            
					if (WasteOfTime == 3) {
						free(LOT);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
						return(0);
					}
            
					RO = pMaxX;
					LO = pMaxX;

					//destroy the maxx peak
            
					Dummy = DestroyPeakP(MaxY, LenStrainSeq0, RO, LO, LenXoverSeq, &LOT[0], SmoothChi, Chivals);
            
					if (maxz == -1 || MaxY == -1) {
						free(LOT);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
						return(0);
					
					}
					Chivals[maxz + MaxY*(LenStrainSeq0 + 1)] = 0;
					if (Dummy == 0) {
						free(LOT);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
						return(0);
					}
            
				}
				//if (Redox > 3) {
					//free(LOT);
					//Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
					//Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
					//return(0);
				//}


			}
			else
				break;
		
				
		}
		free(LOT);
		Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
		Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
		return(0);

	}

	int  MyMathFuncs::FastRecCheckMC2(int SEN, int LongWindedFlag, double *BQPV, int EarlyBale, double UCTHresh, short int MCFlag, short int ShortOutFlag, int MCCorrection, double LowestProb, short int CircularFlag, int NextNo, int MaxABWin, int HWindowWidthX, int lHWindowWidth, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainSeq0, int CriticalDiffX, unsigned char FindallFlag, int UBFSSMC, int UBCS, int Seq1, int Seq2, int Seq3, unsigned char *CS, unsigned char *FSSMC, short int *SeqNum, unsigned char *MissingData, int UBWS, unsigned char *Scores, int *Winscores, int *XDiffPos,  int *XPosDiff, double *Chivals, int *BanWin, unsigned char *MDMap, float *ChiTable2, int *Chimap, double *mtP, double *SmoothChi) {
		int cxv1, tb, te, x, Redox, MaxX, maxz, TWin, MaxFailCount, WinWin, A, C, pMaxX, TopL, TopR, TopLO, TopRO, LO, RO;
		double mPrb, MChi, xMPV, MPV;
		double *LOT;
		short int MaxY;
		LOT = (double *)calloc(2, sizeof(double));
		tb = 0;
		te = 0;
		int	HWindowWidth, CriticalDiff, Dummy;
		int WasteOfTime, LenXoverSeq, GoOn;

		GoOn = 0;
		HWindowWidth = HWindowWidthX;
		CriticalDiff = CriticalDiffX;
		if (SEN == 0)
			LenXoverSeq = FindSubSeqMCPB2(UBFSSMC, UBCS, NextNo, Seq1, Seq2, Seq3, CS, FSSMC, XDiffPos);
		else
			LenXoverSeq = FindSubSeqCP(LenStrainSeq0 + 1, NextNo, Seq1, Seq2, Seq3, SeqNum, XDiffPos, XPosDiff);


		if (LenXoverSeq < CriticalDiff * 2 || LenXoverSeq < 7) {
			free(LOT);
			return(0);
		}


		GoOn = MakeWindowSizeP(tb, te, &CriticalDiff, LenXoverSeq, MCWinFract, MCWinSize, &HWindowWidth, lHWindowWidth, MCProportionFlag);
		//                (int BEP, int ENP, int *CriticalDiff, int LenXoverSeq, double CWinFract, int CWinSize, int *HWindowWidth, int lHWindowWidth, unsigned char *FindallFlag, short int CProportionFlag, int *XPosDiff)
		if (GoOn == 0) {
			free(LOT);
			return(0);
		}

		Dummy = WinScoreCalcP(UBWS, CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0 + 1, Seq1, Seq2, Seq3, Scores, XDiffPos, SeqNum, Winscores);

		//return(0);
		/*if (Seq1 == 0 && Seq2 == 105 && Seq3 == 149)
			cxv1 = 1;*/
		if (FindallFlag == 0 && ((LongWindedFlag == 1 && SEN > 0) || CircularFlag == 0)) {
			if (SEN > 0) {
				if (SEN == 1) {
					Dummy = ClearDeleteArray(LenStrainSeq0 + HWindowWidth * 2, BanWin);
					Dummy = ClearDeleteArrayB(LenStrainSeq0, MDMap);
				}

				Dummy = MakeBanWinP(LenStrainSeq0, Seq1, Seq2, Seq3, HWindowWidth, LenStrainSeq0, LenXoverSeq, BanWin, MDMap, MissingData, XPosDiff, XDiffPos);
			}
			else {

				Dummy = ClearDeleteArrayB(LenStrainSeq0, MDMap);
				Dummy = ClearDeleteArray(LenStrainSeq0 + HWindowWidth * 2, BanWin);
				//Dummy = MakeBanWinP(LenStrainSeq0, Seq1, Seq2, Seq3, HWindowWidth, LenStrainSeq0, LenXoverSeq, BanWin, MDMap, MissingData, XPosDiff, XDiffPos);*/
			}

			if (CircularFlag == 0) {
				//Ban windows spanning the ends of the alignment
				MDMap[1] = 1;
				MDMap[LenXoverSeq] = 1;
				for (x = (LenXoverSeq - HWindowWidth + 2); x <= LenXoverSeq; x++)
					BanWin[x] = 1;
			}

			if (HWindowWidth > 4 && HWindowWidth <= MaxABWin)

				MChi = CalcChiVals4P3(UBWS, CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0, Winscores, Chivals, BanWin, &ChiTable2[Chimap[HWindowWidth]]);
			else
				MChi = CalcChiVals4P(UBWS, CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0, Winscores, Chivals, BanWin);

		}
		else {
			if (SEN == 1) {
				Dummy = ClearDeleteArrayB(LenStrainSeq0, MDMap);
				Dummy = ClearDeleteArray(LenStrainSeq0 + HWindowWidth * 2, BanWin);
			}
			//ReDim MDMap(LenXOverSeq), BanWin(LenXOverSeq + HWindowWidth * 2)
			MChi = CalcChiValsP(UBWS, CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0, Winscores, Chivals);

		}
		//End If
		//return(0);
		if (MCProportionFlag == 0) {
			if ((ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 3) > LowestProb) {
				free(LOT);
				Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
				return(0);
			}

		}
		else {
			if (((ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 6) * MCCorrection) > LowestProb) {
				free(LOT);
				Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
				return(0);
			}
		}




		if (ShortOutFlag == 3)
			Dummy = AddPVal(3, mtP, HWindowWidth, LenXoverSeq, MCCorrection, MCFlag, MChi, LowestProb);
		//int  MyMathFuncs::AddPVal(int Prog, double *mtP,int HWindowWidth,int LenXoverSeq,int MCCorrection,short int MCFlag, double MChi, double LowestProb) {



		//Smooth to find peaks (uses a window of positions)


		//test if this p version is quicker
		//SmoothChiValsP LenXOverSeq, Len(StrainSeq(0)), Chivals(0, 0), SmoothChi(0, 0)

		//The p - version is slower

		Dummy = SmoothChiValsP(LenXoverSeq, LenStrainSeq0, Chivals, SmoothChi);
		//return(0);
		MChi = 0;


		WasteOfTime = 0; //Keeps track of how many cycles have been wasted looking for insignificant hit

						 //Now find recombination events


		Redox = 0;

		while (Redox <= 100) {

			Redox++;

			Dummy = FindMChiP(LenStrainSeq0, LenXoverSeq, &MaxX, &MaxY, &MChi, Chivals);

			if (MaxX == -1 || MaxY == -1 || Redox > 100) {
				free(LOT);
				Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
				return(0);
			}

			//The p-version of this seems slower than the old one
			mPrb = ChiPVal2P(MChi);
			if (MCProportionFlag == 0) {

				if ((mPrb * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 3.0) > LowestProb || mPrb == 1) {
					free(LOT);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(0);
				}

			}
			else {
				if (((mPrb * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 6.0) * MCCorrection) > LowestProb || mPrb == 1) {
					free(LOT);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(0);
				}

			}



			if (mPrb < LowestProb) {
				maxz = MaxX;

				if (MaxX == 0)
					MaxX = 1;

				Dummy = MakeTWinP(FindallFlag, HWindowWidth, &TWin, LenXoverSeq);
				//Call MakeTWin(FindallFlag, HWindowWidth, TWin, LenXoverSeq)


				MaxFailCount = HWindowWidth * 2;

				if (MaxFailCount > int((LenXoverSeq - TWin * 2) / 2))
					MaxFailCount = int((LenXoverSeq - TWin * 2) / 2);

				if (MaxFailCount == 0)
					MaxFailCount = 1;

				WinWin = HWindowWidth;
				A = 0;
				C = 0;
				Dummy = GetACP(LenXoverSeq, LenStrainSeq0, MaxY, MaxX, TWin, &A, &C, Scores);
				//Now find the other breakpoint.
				pMaxX = MaxX;
				if (MaxX < 1) {
					if (CircularFlag == 1)
						MaxX = LenXoverSeq + MaxX;
					else
						MaxX = 1;

				}
				else if (MaxX > LenXoverSeq) {
					if (CircularFlag == 1)
						MaxX = MaxX - LenXoverSeq;
					else
						MaxX = LenXoverSeq - 1;

				}

				if (TWin >= HWindowWidth) {
					TopL = A;
					TopR = C;
				}
				else {
					TopL = 0;
					TopR = 0;
				}

				TWin++;

				//Store best values from the initial screen

				MPV = mPrb * ((float)(LenXoverSeq) / (float)(HWindowWidth));

				WinWin = HWindowWidth;
				TopLO = MaxX - HWindowWidth + 1;
				TopRO = MaxX + HWindowWidth;

				LO = (MaxX - TWin + 1);
				if (LO < 0)
					LO = LenXoverSeq + LO;
				RO = (MaxX + TWin);
				if (RO > LenXoverSeq * 2)
					RO = RO - LenXoverSeq * 2;
				if (RO >= LenXoverSeq) {
					if (MDMap[LenXoverSeq] == 0)
						RO = RO - LenXoverSeq;
					else {
						TWin = TWin - 1;
						RO = LenXoverSeq;
						LO = (MaxX - TWin + 1);
						if (LO < 0)
							LO = LenXoverSeq + LO;

					}
				}
				/*if (Seq1 == 0 && Seq2 == 105 && Seq3 == 149)
					cxv1 = 1;*/

				if (FindallFlag == 0) {


					if (LongWindedFlag == 1 && (CircularFlag == 0 || SEN > 0)) {
						//better use of floats might help here - also the p-version gives a slightly different result to the non-p version
						//These routines give slightly different results to the vc5 versions

						Dummy = GrowMChiWin2P2(MaxABWin, LO, RO, LenXoverSeq, HWindowWidth, TWin, MaxY, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, MDMap, ChiTable2, Chimap);
						//Dummy = GrowMChiWin2P(LO, RO, LenXoverSeq, HWindowWidth, TWin, 0, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, MDMap);





					}
					else {

						Dummy = GrowMChiWinP2(MaxABWin, LO, RO, LenXoverSeq, HWindowWidth, TWin, MaxY, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, ChiTable2, Chimap);
						/*if (WinWin > 1000)
						WinWin = 1000;*/
						
					}

				}
				else {
					Dummy = GrowMChiWinP2(MaxABWin, LO, RO, LenXoverSeq, HWindowWidth, TWin, MaxY, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, ChiTable2, Chimap);
					/*if (WinWin > 1000)
					WinWin = 1000;*/


				}

				if (WinWin < HWindowWidth)
					MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(WinWin));
				else {
					if (MChi < 20000)
						MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth));
					else
						MPV = pow(10, -200);

				}


				xMPV = MPV;
				if (MCProportionFlag == 0)
					xMPV = xMPV * 3;
				else
					xMPV = xMPV * 3;

				if (xMPV < *BQPV) {

					*BQPV = xMPV;


				}

				if (xMPV < UCTHresh && EarlyBale == 1) {
					free(LOT);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(1);
				}
				if (MCFlag == 0)
					xMPV = xMPV * MCCorrection;


				MPV = xMPV;


				if (ShortOutFlag == 3) {
					if (MPV <= mtP[3])
						mtP[3] = MPV;

				}

				LO = MaxX - WinWin;
				RO = MaxX + WinWin - 1;
				if (RO > LenXoverSeq * 2)
					RO = RO - LenXoverSeq * 2;
				if (RO > LenXoverSeq) {
					if (MDMap[RO - LenXoverSeq] == 0)
						RO = RO - LenXoverSeq;
					else
						RO = LenXoverSeq;
				}
				if (LO < 1)
					LO = LO + LenXoverSeq;

				if (MPV < LowestProb) {

					//*BQPV = MPV;
					free(LOT);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(1);

				}
				else {
					WasteOfTime = WasteOfTime + 1;

					if (WasteOfTime == 3) {
						free(LOT);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
						return(0);
					}

					RO = pMaxX;
					LO = pMaxX;

					//destroy the maxx peak

					Dummy = DestroyPeakP(MaxY, LenStrainSeq0, RO, LO, LenXoverSeq, &LOT[0], SmoothChi, Chivals);

					if (maxz == -1 || MaxY == -1) {
						free(LOT);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
						return(0);

					}
					Chivals[maxz + MaxY*(LenStrainSeq0 + 1)] = 0;
					if (Dummy == 0) {
						free(LOT);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
						Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
						return(0);
					}

				}
			}
			else
				break;


		}
		free(LOT);
		Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, Chivals);
		Dummy = CleanChiVals(LenXoverSeq, LenStrainSeq0, SmoothChi);
		return(0);

	}

	int  MyMathFuncs::AddToMapP(int A, int S, double Win, int LS, int *APos, short int *Map) {
		int B, C, os3, SX, EX;
		os3 = (LS + 1);
		SX = (S - (int)((Win) / 2));
		EX = (S + (int)((Win) / 2));
		//return (EX);                            
		for (B = SX; B <= EX; B++) {
			if (B < 1)
				C = APos[LS] + B;
			else if (B > APos[LS])
				C = B - APos[LS];
			else
				C = B;


			Map[C + A*os3] = Map[C + A*os3] + 1;

		}
		return(1);
	}

	int  MyMathFuncs::WinScoreCalc4P(int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int Seq1, int Seq2, int Seq3, unsigned char *Scores, int *XDiffPos, short int *SeqNum, int *WinScores)
	{
		int goon, X, LO, RO, TO;
		TO = LenSeq;
		goon = 0;


		//Calculate scores per position
		for (X = 1; X <= LenXoverSeq; X++) {

			if (*(SeqNum + XDiffPos[X] + LenSeq*Seq1) == *(SeqNum + XDiffPos[X] + LenSeq*Seq2)) {
				*(Scores + X) = 1;
			}
			else {
				*(Scores + X) = 0;
			}

			*(WinScores + X) = 0;
		}



		//Calculate scores for first window

		*(WinScores) = 0;

		for (X = (LenXoverSeq - HWindowWidth + 1); X <= LenXoverSeq; X++) {
			*(WinScores) = *(WinScores)+*(Scores + X);
		}

		//Calculate scores for windows traversing the left end
		for (X = 1; X <= HWindowWidth; X++) {
			LO = ((LenXoverSeq - HWindowWidth) + X);
			*(WinScores + X) = *(WinScores + X - 1) - *(Scores + LO) + *(Scores + X);
		}

		//Calculate scores for internal windows
		for (X = HWindowWidth + 1; X <= LenXoverSeq; X++) {
			LO = X - HWindowWidth;
			*(WinScores + X) = *(WinScores + X - 1) - *(Scores + LO) + *(Scores + X);
		}

		//Calculate scores for windows traversing the right end
		for (X = LenXoverSeq + 1; X < LenXoverSeq + HWindowWidth; X++) {
			LO = (X - HWindowWidth);
			RO = X - LenXoverSeq;
			*(WinScores + X) = *(WinScores + X - 1) - *(Scores + LO) + *(Scores + RO);
		}

		for (X = 1; X < LenXoverSeq + HWindowWidth; X++) {
			if (*(WinScores + X) - *(WinScores + X + HWindowWidth - 1) > criticaldiff || *(WinScores + X) - *(WinScores + X + HWindowWidth - 1) < -criticaldiff) {
				goon = 1;
				break;
			}
		}

		return (goon);
	}

	int  MyMathFuncs::WinScoreCalc4P2(int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int Seq1, int Seq2, int Seq3, unsigned char *Scores, int *XDiffPos, short int *SeqNum, int *WinScores)
	{
		int goon, X, LO, RO, TO, os1, os2, h1;
		TO = LenSeq;
		goon = 0;
		os1 = LenSeq*Seq1;
		os2 = LenSeq*Seq2;
		//Calculate scores per position
		for (X = 1; X <= LenXoverSeq; X++) {
			h1 = XDiffPos[X];
			*(Scores + X) = (*(SeqNum + h1 + os1) == *(SeqNum + h1 + os2));
			*(WinScores + X) = 0;
		}



		//Calculate scores for first window

		*(WinScores) = 0;

		for (X = (LenXoverSeq - HWindowWidth + 1); X <= LenXoverSeq; X++) {
			*(WinScores) = *(WinScores)+*(Scores + X);
		}

		//Calculate scores for windows traversing the left end
		for (X = 1; X <= HWindowWidth; X++) {
			LO = ((LenXoverSeq - HWindowWidth) + X);
			*(WinScores + X) = *(WinScores + X - 1) - *(Scores + LO) + *(Scores + X);
		}

		//Calculate scores for internal windows
		for (X = HWindowWidth + 1; X <= LenXoverSeq; X++) {
			LO = X - HWindowWidth;
			*(WinScores + X) = *(WinScores + X - 1) - *(Scores + LO) + *(Scores + X);
		}

		//Calculate scores for windows traversing the right end
		for (X = LenXoverSeq + 1; X < LenXoverSeq + HWindowWidth; X++) {
			LO = (X - HWindowWidth);
			RO = X - LenXoverSeq;
			*(WinScores + X) = *(WinScores + X - 1) - *(Scores + LO) + *(Scores + RO);
		}

		for (X = 1; X < LenXoverSeq + HWindowWidth; X++) {
			if (*(WinScores + X) - *(WinScores + X + HWindowWidth - 1) > criticaldiff || *(WinScores + X) - *(WinScores + X + HWindowWidth - 1) < -criticaldiff) {
				goon = 1;
				break;
			}
		}

		return (goon);
	}
	double MyMathFuncs::CalcChiVals3P(int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals)
	{

		double A, B, C, D, E, ChiH, MChi;
		int X, Y, LO, FO, SO;
		FO = LenSeq + HWindowWidth * 2 + 1;
		SO = LenSeq + 1;
		LO = HWindowWidth;
		ChiH = 0;
		MChi = 0;
		Y = 0;
		for (X = 0; X < LenXoverSeq; X++) {

			A = *(WinScores + X);
			C = *(WinScores + X + LO);
			B = HWindowWidth - A;
			D = HWindowWidth - C;

			if (A + C > 0 && B + D >0) {
				E = A*D - B*C;
				ChiH = E*E * 2 / (HWindowWidth*(A + C)*(B + D));
				*(ChiVals + X) = ChiH;
				if (ChiH > MChi) {
					MChi = ChiH;
					//ma=A;
				}
			}
			else
				*(ChiVals + X) = 0;

		}
		return (MChi);
		//return (ma);
	}
	int  MyMathFuncs::SmoothChiVals3P(int LenXoverSeq, int LenSeq, double *ChiVals, double *SmoothChi)
	{
		short int Y;
		int  RO, X;
		double RunCount;
		int qWindowSize = 5;
		RO = qWindowSize * 2 + 1;
		Y = 0;
		RunCount = 0;
		for (X = 0 - qWindowSize; X <= 1 + qWindowSize; X++) {
			if (X < 1)
				RunCount += *(ChiVals + LenXoverSeq + X);
			else
				RunCount += *(ChiVals + X);
		}
		*(SmoothChi) = RunCount / RO;
		for (X = 1 - qWindowSize; X < LenXoverSeq - qWindowSize; X++) {
			if (X > 0 && X + RO <= LenXoverSeq)
				RunCount = RunCount - *(ChiVals + X) + *(ChiVals + X + RO);
			else if (X + RO > LenXoverSeq)
				RunCount = RunCount - *(ChiVals + X) + *(ChiVals + X + RO - LenXoverSeq);
			else if (X < 1)
				RunCount = RunCount - *(ChiVals + LenXoverSeq + X) + *(ChiVals + X + RO);

			*(SmoothChi + X + qWindowSize) = RunCount / RO;
		}



		return 1;
	}

	int  MyMathFuncs::FindMChi3P(int LenSeq, int LenXoverSeq, int *MaxX, short int *MaxY, double *MChi, double *ChiVals)
	{
		int X, SO, tMaxX;
		short int Y, tMaxY;
		double ChiV, tMChi;
		tMaxX = -1;
		tMaxY = -1;
		tMChi = 0;
		ChiV = 0;
		SO = LenSeq + 1;
		Y = 0;
		for (X = 0; X < LenXoverSeq; X++) {

			ChiV = *(ChiVals + X);
			if (ChiV > tMChi) {
				tMChi = ChiV;
				tMaxX = X;
				tMaxY = Y;
			}

		}
		*MChi = tMChi;
		*MaxX = tMaxX;
		*MaxY = tMaxY;
		return 1;
	}

	int  MyMathFuncs::FastRecCheckChim(unsigned char *MissingData, int *XPosDiff, int *LXOS, int YP, int SEN, int LongWindedFlag, double *BQPV, int EarlyBale, double UCTHresh, short int MCFlag, short int ShortOutFlag, int MCCorrection, double LowestProb, short int CircularFlag, int NextNo, int MaxABWin, int HWindowWidthX, int lHWindowWidth, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainSeq0, int CriticalDiffX, unsigned char FindallFlag, int UBFSSRDP, int UBCS, int Seq1, int Seq2, int Seq3, unsigned char *CS, unsigned char *FSSRDP, short int *SeqNum, int UBWS, unsigned char *Scores, int *Winscores, int UBXDP, int *XDP, double *Chivals, int *BanWin, unsigned char *MDMap, float *ChiTable2, int *Chimap, double *mtP, double *SmoothChi) {
		int dx1, tb, te, x, Redox, MaxX, maxz, TWin, MaxFailCount, WinWin, A, C, pMaxX, TopL, TopR, TopLO, TopRO, LO, RO;
		double mPrb, MChi, xMPV, MPV;
		double *LOT;
		short int MaxY;
		LOT = (double *)calloc(2, sizeof(double));
		tb = 0;
		te = 0;
		int	HWindowWidth, CriticalDiff, Dummy;
		int WasteOfTime, LenXoverSeq, GoOn;

		GoOn = 0;
		HWindowWidth = HWindowWidthX;
		CriticalDiff = CriticalDiffX;
		LenXoverSeq = 0;
		//LenXoverSeq = FindSubSeqDP(UBFSSMC, UBCS, NextNo, Seq1, Seq2, Seq3, CS, FSSMC, XDiffPos);
		//LenXoverSeq = FindSubSeqDP2(LenStrainSeq0 + 1, NextNo, Seq1, Seq2, Seq3, SeqNum, XDiffPos);




		if (SEN == 0) {
			if (YP == 0) {
				Dummy = FindSubSeqDP3(UBFSSRDP, UBCS, FSSRDP, CS, LenStrainSeq0 + 1, NextNo, Seq1, Seq2, Seq3, LXOS, UBXDP, XDP);
				
			}
			LenXoverSeq = LXOS[YP];
		}
		else {

			if (YP == 0) {
				if (SEN == 0)
					Dummy = FindSubSeqDP3(UBFSSRDP, UBCS, FSSRDP, CS, LenStrainSeq0 + 1, NextNo, Seq1, Seq2, Seq3, LXOS, UBXDP, XDP);
				else
					Dummy = FindSubSeqDP6(UBFSSRDP, UBCS, FSSRDP, CS, LenStrainSeq0 + 1, NextNo, Seq1, Seq2, Seq3, LXOS, UBXDP, XDP, XPosDiff);

			}
			LenXoverSeq = LXOS[YP];
			//LenXoverSeq = FindSubSeqDP(LenStrainSeq0 + 1, NextNo, Seq1, Seq2, Seq3, SeqNum, &XDP[YP*(UBXDP + 1)], &XPosDiff[YP*(UBXDP + 1)]);
			//LenXoverSeq = FindSubSeqDP(LenStrainSeq0 + 1, NextNo, Seq1, Seq2, Seq3, SeqNum, XDP, XPosDiff);

		}
		
		if (LenXoverSeq < CriticalDiff * 2 || LenXoverSeq < 7) {
			free(LOT);
			return(0);
		}


		GoOn = MakeWindowSizeP(tb, te, &CriticalDiff, LenXoverSeq, MCWinFract, MCWinSize, &HWindowWidth, lHWindowWidth, MCProportionFlag);
		//                (int BEP, int ENP, int *CriticalDiff, int LenXoverSeq, double CWinFract, int CWinSize, int *HWindowWidth, int lHWindowWidth, unsigned char *FindallFlag, short int CProportionFlag, int *XPosDiff)
		if (GoOn == 0) {
			free(LOT);
			return(0);
		}

		Dummy = WinScoreCalc4P2(CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0 + 1, Seq1, Seq2, Seq3, Scores, &XDP[YP*(UBXDP+1)], SeqNum, Winscores);

		//return(0);
		//If FindallFlag = 0 And ((LongWindedFlag = 1 And SEventNumber > 0) Or CircularFlag = 0) Then
		if (FindallFlag == 0 && ((LongWindedFlag == 1 && SEN > 0) || CircularFlag == 0)) {
			if (SEN > 0){
				if (SEN == 1){
					//Dummy = ClearDeleteArray(LenStrainSeq0, BanWin);
					Dummy = ClearDeleteArray(LenStrainSeq0 + HWindowWidth * 2, BanWin);
					Dummy = ClearDeleteArrayB(LenStrainSeq0, MDMap);
				}
        
				Dummy = MakeBanWinP(LenStrainSeq0 + HWindowWidth * 2, Seq1, Seq2, Seq3, HWindowWidth, LenStrainSeq0, LenXoverSeq, BanWin, MDMap, MissingData,  &XPosDiff[YP*(UBXDP + 1)], &XDP[YP*(UBXDP + 1)]);
			}
			else {
				//Dummy = ClearDeleteArray(LenStrainSeq0, BanWin);
				Dummy = ClearDeleteArray(LenStrainSeq0 + HWindowWidth * 2, BanWin);
				Dummy = ClearDeleteArrayB(LenStrainSeq0, MDMap);
			}



			if (CircularFlag == 0) {
				//Ban windows spanning the ends of the alignment
				MDMap[1] = 1;
				MDMap[LenXoverSeq] = 1;
				for (x = (LenXoverSeq - HWindowWidth + 2); x <= LenXoverSeq; x++)
					BanWin[x] = 1;
			}

			MChi = CalcChiVals5P(CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0, Winscores, Chivals, BanWin);
			//MChi = CalcChiVals3P(CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0, Winscores, Chivals);
		}
		else {
			if (SEN == 1) {
				Dummy = ClearDeleteArrayB(LenStrainSeq0, MDMap);
				Dummy = ClearDeleteArray(LenStrainSeq0 + HWindowWidth * 2, BanWin);
			}
			//ReDim MDMap(LenXOverSeq), BanWin(LenXOverSeq + HWindowWidth * 2)
			MChi = CalcChiVals3P(CriticalDiff, HWindowWidth, LenXoverSeq, LenStrainSeq0, Winscores, Chivals);

		}
		//End If
		//return(0);
		if (MCProportionFlag == 0) {
			if ((ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 3) > LowestProb) {
				free(LOT);
				Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
				Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
				return(0);
			}

		}
		else {
			if (((ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 6) * MCCorrection) > LowestProb) {
				free(LOT);
				Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
				Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
				return(0);
			}
		}

		if (ShortOutFlag == 3)
			Dummy = AddPVal(4, mtP, HWindowWidth, LenXoverSeq, MCCorrection, MCFlag, MChi, LowestProb);
		//int  MyMathFuncs::AddPVal(int Prog, double *mtP,int HWindowWidth,int LenXoverSeq,int MCCorrection,short int MCFlag, double MChi, double LowestProb) {



		//Smooth to find peaks (uses a window of positions)


		//test if this p version is quicker
		//SmoothChiValsP LenXOverSeq, Len(StrainSeq(0)), Chivals(0, 0), SmoothChi(0, 0)

		//The p - version is slower

		Dummy = SmoothChiVals3P(LenXoverSeq, LenStrainSeq0, Chivals, SmoothChi);
		//return(0);
		MChi = 0;

		double oxChi;
		WasteOfTime = 0; //Keeps track of how many cycles have been wasted looking for insignificant hit

						 //Now find recombination events


		Redox = 0;
		double lmprb;
		lmprb = 0;


		while (Redox <= 100) {

			Redox++;
			
			Dummy = FindMChi3P(LenStrainSeq0, LenXoverSeq, &MaxX, &MaxY, &MChi, Chivals);
			
			if (MaxX == -1 || MaxY == -1 || Redox > 100 || MChi < 0) {
				free(LOT);
				Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
				Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
				return(0);
			}

			//The p-version of this seems slower than the old one
			mPrb = ChiPVal2P(MChi);
			if (mPrb == lmprb) {
				free(LOT);
				Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
				Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
				return(0);

			}
			lmprb = mPrb;
			if (MCProportionFlag == 0) {

				if ((mPrb * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 3.0) > LowestProb || mPrb == 1) {
					free(LOT);
					Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(0);
				}

			}
			else {
				if (((mPrb * ((float)(LenXoverSeq) / (float)(HWindowWidth)) * 6.0) * MCCorrection) > LowestProb || mPrb == 1) {
					free(LOT);
					Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(0);
				}

			}
			if (mPrb < LowestProb && mPrb != 1) {
				
				maxz = MaxX;

				if (MaxX == 0)
					MaxX = 1;

				Dummy = MakeTWinP(FindallFlag, HWindowWidth, &TWin, LenXoverSeq);
				//Call MakeTWin(FindallFlag, HWindowWidth, TWin, LenXoverSeq)


				MaxFailCount = HWindowWidth * 2;

				if (MaxFailCount > int((LenXoverSeq - TWin * 2) / 2))
					MaxFailCount = int((LenXoverSeq - TWin * 2) / 2);

				if (MaxFailCount == 0)
					MaxFailCount = 1;

				WinWin = HWindowWidth;
				A = 0;
				C = 0;
				Dummy = GetACP(LenXoverSeq, LenStrainSeq0, 0, MaxX, TWin, &A, &C, Scores);
				//Now find the other breakpoint.
				pMaxX = MaxX;
				if (MaxX < 1) {
					if (CircularFlag == 1)
						MaxX = LenXoverSeq + MaxX;
					else
						MaxX = 1;

				}
				else if (MaxX > LenXoverSeq) {
					if (CircularFlag == 1)
						MaxX = MaxX - LenXoverSeq;
					else
						MaxX = LenXoverSeq - 1;

				}

				if (TWin >= HWindowWidth) {
					TopL = A;
					TopR = C;
				}
				else {
					TopL = 0;
					TopR = 0;
				}

				TWin++;

				//Store best values from the initial screen

				MPV = mPrb * ((float)(LenXoverSeq) / (float)(HWindowWidth));

				WinWin = HWindowWidth;
				TopLO = MaxX - HWindowWidth + 1;
				TopRO = MaxX + HWindowWidth;

				LO = (MaxX - TWin + 1);
				if (LO < 0)
					LO = LenXoverSeq + LO;
				RO = (MaxX + TWin);
				if (RO > LenXoverSeq * 2)
					RO = RO - LenXoverSeq * 2;
				if (RO >= LenXoverSeq) {
					if (MDMap[LenXoverSeq] == 0)
						RO = RO - LenXoverSeq;
					else {
						TWin = TWin - 1;
						RO = LenXoverSeq;
						LO = (MaxX - TWin + 1);
						if (LO < 0)
							LO = LenXoverSeq + LO;

					}
				}
				oxChi = MChi;
				if (FindallFlag == 0) {


					if (LongWindedFlag == 1 && (CircularFlag == 0 || SEN > 0)) {
						//better use of floats might help here - also the p-version gives a slightly different result to the non-p version
						//These routines give slightly different results to the vc5 versions
						Dummy = GrowMChiWin2P(LO, RO, LenXoverSeq, HWindowWidth, TWin, 0, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, MDMap);
                

						

					}
					else {
						Dummy = GrowMChiWinP(LO, RO, LenXoverSeq, HWindowWidth, TWin, 0, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores);

						//Dummy = GrowMChiWinP2(MaxABWin, LO, RO, LenXoverSeq, HWindowWidth, TWin, 0, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, ChiTable2, Chimap);
						/*if (WinWin > 1000)
							WinWin = 1000;*/
						
					}

				}
				else {
					Dummy = GrowMChiWinP(LO, RO, LenXoverSeq, HWindowWidth, TWin, 0, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores);
					//Dummy = GrowMChiWinP2(MaxABWin, LO, RO, LenXoverSeq, HWindowWidth, TWin, MaxY, LenStrainSeq0, A, C, MaxFailCount, &MPV, &WinWin, &MChi, &TopL, &TopR, &TopLO, &TopRO, Scores, ChiTable2, Chimap);
					/*if (WinWin > 1000)
						WinWin = 1000;*/
					


				}

				if (MChi < 0)
					MChi = oxChi;
				if (WinWin < HWindowWidth)
					MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(WinWin));
				else {
					if (MChi < 20000)
						MPV = ChiPVal2P(MChi) * ((float)(LenXoverSeq) / (float)(HWindowWidth));
					else
						MPV = pow(10, -200);

				}

				xMPV = MPV;
				if (MCProportionFlag == 0)
					xMPV = xMPV * 3;
				else
					xMPV = xMPV * 3;

				if (xMPV < *BQPV) {

					*BQPV = xMPV;


				}

				if (xMPV < UCTHresh && EarlyBale == 1) {
					free(LOT);
					Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(1);
				}
				if (MCFlag == 0)
					xMPV = xMPV * MCCorrection;


				MPV = xMPV;


				if (ShortOutFlag == 3) {
					if (MPV <= mtP[3])
						mtP[3] = MPV;

				}

				LO = MaxX - WinWin;
				RO = MaxX + WinWin - 1;
				if (RO > LenXoverSeq * 2)
					RO = RO - LenXoverSeq * 2;
				if (RO > LenXoverSeq) {
					if (MDMap[RO - LenXoverSeq] == 0)
						RO = RO - LenXoverSeq;
					else
						RO = LenXoverSeq;
				}
				if (LO < 1)
					LO = LO + LenXoverSeq;

				if (MPV < LowestProb) {

					//*BQPV = MPV;
					free(LOT);
					Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
					Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
					return(1);

				}
				else {
					WasteOfTime = WasteOfTime + 1;

					if (WasteOfTime == 3) {
						free(LOT);
						Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
						Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
						return(0);
					}

					RO = pMaxX;
					LO = pMaxX;

					//destroy the maxx peak

					Dummy = DestroyPeakP(0, LenStrainSeq0, RO, LO, LenXoverSeq, &LOT[0], SmoothChi, Chivals);

					if (maxz == -1 || MaxY == -1) {
						free(LOT);
						Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
						Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
						return(0);

					}
					Chivals[maxz + MaxY*(LenStrainSeq0 + 1)] = 0;
					if (Dummy == 0) {
						free(LOT);
						Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
						Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
						return(0);
					}

				}
			}
			else
				break;

		}
		free(LOT);
		Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, Chivals);
		Dummy = CleanChiVals2(LenXoverSeq, LenStrainSeq0, SmoothChi);
		return(0);

	}


	double  MyMathFuncs::CalcChiVals5P(int criticaldiff, int HWindowWidth, int LenXoverSeq, int LenSeq, int *WinScores, double *ChiVals, int *BanWins)
	{

		double A, B, C, D, E, ChiH, MChi;
		int X, LO, FO, SO, bp1;
		FO = LenSeq + HWindowWidth * 2 + 1;
		SO = LenSeq + 1;
		LO = HWindowWidth;
		ChiH = 0;
		MChi = 0;
		for (X = 0; X < LenXoverSeq; X++) {
			bp1 = X - HWindowWidth;
			if (bp1 < 1)
				bp1 = bp1 + LenXoverSeq;
			if (BanWins[X] == 0 && BanWins[bp1] == 0) {
				A = *(WinScores + X);
				C = *(WinScores + X + LO);
				if (A - C > criticaldiff || A - C < -criticaldiff) {
					B = HWindowWidth - A;
					D = HWindowWidth - C;
					if (A + C > 0 && B + D > 0) {
						E = A*D - B*C;
						ChiH = E*E * 2 / (HWindowWidth*(A + C)*(B + D));
						*(ChiVals + X) = ChiH;
						if (ChiH > MChi)
							MChi = ChiH;
					}
					else
						*(ChiVals + X) = 0;
				}
				else
					*(ChiVals + X) = 0;
			}
			else
				*(ChiVals + X) = 0;
		}
		return (MChi);
		//return (ma);
	}

	int  MyMathFuncs::GetFragsP3(short int CircularFlag, int LenXoverSeq, int LSeq, int maxcount, char *SubSeq, int *FragSt, int *FragEn, int *FragScore, int *FragCount) {

		int SX, A, X, Y, Z, GoOn, os, os2, os3, FE, FS, RR, FSc, FC6, os4, os5, target;

		GoOn = 0;
		os = maxcount + 1;
		RR = LenXoverSeq + 1;
		FC6 = 0;
		os4 = LSeq + 1;
		os5 = maxcount + 1;

		target = 6 * os;
		X = 0;
		Y = 0;
		Z = 0;



		//
		//		omp_set_num_threads(3);
		//#pragma omp parallel 
		//		{
		//#pragma omp sections private (X)
		//			{
		//#pragma omp section
		//				{
		for (X = 0; X <= target; X++)
			FragScore[X] = 0;

		//				}
		//#pragma omp section
		//				{
		for (X = 0; X <= target; X++)

			FragSt[X] = 0;
		//				}
		//
		//#pragma omp section
		//				{
		for (X = 0; X <= target; X++)

			FragEn[X] = 0;
		/*}
		}

		}*/

		//target = 
		/*for (Y = 0; Y <= 6; Y++){
		osx = Y*os;
		for (X = 0; X < RR; X++)
		FragScore[X+osx] = 0;
		}
		for (Y = 0; Y <= 6; Y++){
		osx = Y*os;
		for (X = 0; X <= RR; X++){
		FragSt[X+osx] = 0;
		FragEn[X+osx] = 0;

		}
		}*/

		if (CircularFlag == 0) {
			for (X = 1; X < RR; X++) {
				for (Y = 0; Y < 3; Y++) {
					if (SubSeq[X + Y*os4] == 1)
						break;
				}
				if (Y == 3)
					Y = 6;

				os2 = Y*os4;
				os = FragCount[Y] + Y*os5;
				FragSt[os] = X;
				FS = X;
				SX = X;
				X++;
				
				while (SubSeq[os2 + X] == 1) {
					X++;
					if (X > LenXoverSeq)
						break;
				}
				
				X--;

				FragEn[os] = X;
				FE = X;
				if (FE > FS)
					FSc = FE - FS + 1;
				else if (FE < FS)
					FSc = FE - FS + RR;
				else
					FSc = 1;

				FragScore[os] = FSc;

				if (Y != 0) {
					os3 = FragCount[0];
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = -FSc;
					FragCount[0] = FragCount[0] + 1;
				}

				if (Y != 1) {
					os3 = FragCount[1] + os5;
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = -FSc;
					FragCount[1] = FragCount[1] + 1;
				}
				if (Y != 2) {
					os3 = FragCount[2] + 2 * os5;
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = -FSc;
					FragCount[2] = FragCount[2] + 1;
				}


				if (Y == 0)
					Z = 5;
				else if (Y == 1)
					Z = 4;
				else if (Y == 2)
					Z = 3;
				else if (Y == 6)
					Z = 6;

				if (Z < 6) {
					os3 = FragCount[Z] + Z*os5;
					if (FragCount[Z] > 0) {
						if (FragScore[os3 - 1] > 0) {
							FragEn[os3 - 1] = FE;
							FragScore[os3 - 1] = FragScore[os3 - 1] + FSc;
						}
						else {
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = FSc;
							FragCount[Z] = FragCount[Z] + 1;
						}
					}
					else {
						FragSt[os3] = FS;
						FragEn[os3] = FE;
						FragScore[os3] = FSc;
						FragCount[Z] = FragCount[Z] + 1;
					}

					for (A = 3; A < 6; A++) {
						if (A != Z) {
							os3 = FragCount[A] + A*os5;
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = -FSc;
							FragCount[A] = FragCount[A] + 1;
						}
					}

					


				}
				else {
					if (FragCount[3] > 0) {

						for (A = 3; A < 6; A++) {
							os3 = FragCount[A] + A*os5;
							if (FragScore[os3 - 1] > 0) {

								FragEn[os3 - 1] = FE;
								FragScore[os3 - 1] = FragScore[os3 - 1] + FSc;
							}
							else {
								FragSt[os3] = FS;
								FragEn[os3] = FE;
								FragScore[os3] = FSc;
								FragCount[A] = FragCount[A] + 1;
							}
						}//Next A
					}
					else {

						for (A = 3; A < 6; A++) {
							os3 = FragCount[A] + A*os5;
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = FSc;
							FragCount[A] = FragCount[A] + 1;
						}//Next A
					}
				}

				/*if (Y <= 2)
				FragCount[Y] = FragCount[Y] + 1;*/

				if (Y <= 2)
					FragCount[Y] = FragCount[Y] + 1;


				for (A = 0; A < 6; A++) {
					if (FragCount[A] >= maxcount - 1) {
						GoOn = 1;
						break;
					}


				}




				FC6++;
				if (GoOn == 1 || FC6 >= maxcount - 1)
					break;



			}//next X
		}
		else

		{
			for (X = 1; X < RR; X++) {
				for (Y = 0; Y < 3; Y++) {
					if (SubSeq[X + Y*os4] == 1)
						break;
				}
				if (Y == 3)
					Y = 6;

				os2 = Y*os4;
				os = FragCount[Y] + Y*os5;
				FragSt[os] = X;
				FS = X;
				SX = X;
				X++;
				
				

				while (SubSeq[os2 + X] == 1) {
					X++;
					if (X > LenXoverSeq) {
						X = 1;
						GoOn = 1;
						break;

					}
				}
				if (X == 1) {
					while (SubSeq[os2 + X] == 1 && X != SX) {
						X++;

					}
				}

				if (X == SX)
					return(0);
				



				X--;

				FragEn[os] = X;
				FE = X;
				if (FE > FS)
					FSc = FE - FS + 1;
				else if (FE < FS)
					FSc = FE - FS + RR;
				else
					FSc = 1;

				FragScore[os] = FSc;


				/*for (Z = 0; Z < 3; Z++){
				if (Y != Z){
				os3 = FragCount[Z] + Z*os5;
				FragSt[os3] = FS;
				FragEn[os3] = FE;
				FragScore[os3] = -FSc;
				FragCount[Z] = FragCount[Z] + 1;
				}
				}*/


				if (Y != 0) {
					os3 = FragCount[0];
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = -FSc;
					FragCount[0] = FragCount[0] + 1;
				}

				if (Y != 1) {
					os3 = FragCount[1] + os5;
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = -FSc;
					FragCount[1] = FragCount[1] + 1;
				}
				if (Y != 2) {
					os3 = FragCount[2] + 2 * os5;
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = -FSc;
					FragCount[2] = FragCount[2] + 1;
				}


				if (Y == 0)
					Z = 5;
				else if (Y == 1)
					Z = 4;
				else if (Y == 2)
					Z = 3;
				else if (Y == 6)
					Z = 6;

				if (Z < 6) {
					os3 = FragCount[Z] + Z*os5;
					if (FragCount[Z] > 0) {
						if (FragScore[os3 - 1] > 0) {
							FragEn[os3 - 1] = FE;
							FragScore[os3 - 1] = FragScore[os3 - 1] + FSc;
						}
						else {
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = FSc;
							FragCount[Z] = FragCount[Z] + 1;
						}
					}
					else {
						FragSt[os3] = FS;
						FragEn[os3] = FE;
						FragScore[os3] = FSc;
						FragCount[Z] = FragCount[Z] + 1;
					}

					for (A = 3; A < 6; A++) {
						if (A != Z) {
							os3 = FragCount[A] + A*os5;
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = -FSc;
							FragCount[A] = FragCount[A] + 1;
						}
					}

					


				}
				else {
					if (FragCount[3] > 0) {

						for (A = 3; A < 6; A++) {
							os3 = FragCount[A] + A*os5;
							if (FragScore[os3 - 1] > 0) {

								FragEn[os3 - 1] = FE;
								FragScore[os3 - 1] = FragScore[os3 - 1] + FSc;
							}
							else {
								FragSt[os3] = FS;
								FragEn[os3] = FE;
								FragScore[os3] = FSc;
								FragCount[A] = FragCount[A] + 1;
							}
						}//Next A
					}
					else {

						for (A = 3; A < 6; A++) {
							os3 = FragCount[A] + A*os5;
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = FSc;
							FragCount[A] = FragCount[A] + 1;
						}//Next A
					}
				}

				/*if (Y <= 2)
				FragCount[Y] = FragCount[Y] + 1;*/

				if (Y <= 2)
					FragCount[Y] = FragCount[Y] + 1;


				for (A = 0; A < 6; A++) {
					if (FragCount[A] >= maxcount - 1) {
						GoOn = 1;
						break;
					}


				}




				FC6++;
				if (GoOn == 1 || FC6 >= maxcount - 1)
					break;



			}//next X
		}
		FragCount[6] = FC6;


		for (X = 0; X <= 6; X++) {
			os = FragCount[X] + X*os5;
			FragSt[os] = LenXoverSeq;
			FragEn[os] = LenXoverSeq;
		}
		//omp_set_num_threads(2);
		return(1);
	}


	int  MyMathFuncs::GetFragsP2(char *goon, int elementnum, int UBFC, int UBFS1, int UBFS2, int UBSS1, int UBSS2, short int CircularFlag, int *LenXoverSeq, int LSeq, int maxcount, char *SubSeq, int *FragSt, int *FragEn, int *FragScore, int *FragCount) {

		int e, SX, A, X, Y, Z, GoOn, os, os2, os3, FE, FS, RR, FSc, FC6, os4, os5, target, osfc, osss, osfs, osfc2, osss2, osfs2, os6, lxos;

		
		os = maxcount + 1;
		
		
		os4 = LSeq + 1;
		os5 = maxcount + 1;

		osss = (UBSS1 + 1)*(UBSS2 + 1);
		osfc = UBFC + 1;
		osfs = (UBFS1 + 1)*(UBFS2 + 1);

		target = (UBFS1+1)*(UBFS2+1)*(elementnum+1);
		



omp_set_num_threads(3);
#pragma omp parallel 
		{
#pragma omp sections private (X)
		{
#pragma omp section
		{
			for (X = 0; X < target; X++)
				FragScore[X] = 0;

		}
#pragma omp section
		{
			for (X = 0; X < target; X++)

				FragSt[X] = 0;
		}

#pragma omp section
		{
			for (X = 0; X < target; X++)

				FragEn[X] = 0;
		}
		}

		}

		

		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);


#pragma omp parallel for private(e, osss2, osfc2, osfs2, lxos, X, Y, Z, os6, FC6, RR, os2, os, FS, SX, FE, FSc, os3, A, GoOn)		
		for (e = 0; e <= elementnum; e++) {
			
			
			osss2 = osss*e;
			osfc2 = osfc*e;

			FragCount[osfc2] = 0;
			FragCount[1 + osfc2] = 0;
			FragCount[2 + osfc2] = 0;
			FragCount[3 + osfc2] = 0;
			FragCount[4 + osfc2] = 0;
			FragCount[5 + osfc2] = 0;
			FragCount[6 + osfc2] = 0;
			osfs2 = osfs*e;
			lxos = LenXoverSeq[e];

			goon[e] = 1;
			GoOn = 0;
			X = 0;
			Y = 0;
			Z = 0;
			FC6 = 0;

			RR = lxos + 1;
			for (X = 1; X < RR; X++) {
				os6 = X + osss2;
				for (Y = 0; Y < 3; Y++) {
					if (SubSeq[os6 + Y*os4] == 1)
						break;
				}
				if (Y == 3)
					Y = 6;

				os2 = Y*os4 + osss2;
				os = FragCount[Y + osfc2] + Y*os5 + osfs2;
				FragSt[os] = X;
				FS = X;
				SX = X;
				X++;
				if (CircularFlag == 0) {
					while (SubSeq[os2 + X] == 1) {
						X++;
						if (X > lxos)
							break;
					}
				}
				else {

					while (SubSeq[os2 + X] == 1) {
						X++;
						if (X > lxos) {
							X = 1;
							GoOn = 1;
							break;

						}
					}
					if (X == 1) {
						while (SubSeq[os2 + X] == 1 && X != SX) {
							X++;

						}
					}

					if (X == SX) {
						goon[e] = 0;
						break;
					}
				}



				X--;

				FragEn[os] = X;
				FE = X;
				if (FE > FS)
					FSc = FE - FS + 1;
				else if (FE < FS)
					FSc = FE - FS + RR;
				else
					FSc = 1;

				FragScore[os] = FSc;


				/*for (Z = 0; Z < 3; Z++){
				if (Y != Z){
				os3 = FragCount[Z] + Z*os5;
				FragSt[os3] = FS;
				FragEn[os3] = FE;
				FragScore[os3] = -FSc;
				FragCount[Z] = FragCount[Z] + 1;
				}
				}*/


				if (Y != 0) {
					os3 = FragCount[osfc2 ] + osfs2;
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = -FSc;
					//if (FragCount[osfc2] < maxcount)
						FragCount[osfc2] = FragCount[osfc2] + 1;
				}

				if (Y != 1) {
					os3 = FragCount[1 + osfc2] + os5 + osfs2;
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = -FSc;
					//if (FragCount[1 + osfc2] < maxcount)
						FragCount[1 + osfc2] = FragCount[1 + osfc2] + 1;
				}
				if (Y != 2) {
					os3 = FragCount[2 + osfc2] + 2 * os5 + osfs2;
					FragSt[os3] = FS;
					FragEn[os3] = FE;
					FragScore[os3] = -FSc;
					//if (FragCount[2 + osfc2] < maxcount)
						FragCount[2 + osfc2] = FragCount[2 + osfc2] + 1;
				}


				if (Y == 0)
					Z = 5;
				else if (Y == 1)
					Z = 4;
				else if (Y == 2)
					Z = 3;
				else if (Y == 6)
					Z = 6;

				if (Z < 6) {
					os3 = FragCount[Z + osfc2] + Z*os5 + osfs2;
					if (FragCount[Z + osfc2] > 0) {
						if (FragScore[os3 - 1] > 0) {
							FragEn[os3 - 1] = FE;
							FragScore[os3 - 1] = FragScore[os3 - 1] + FSc;
						}
						else {
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = FSc;
							//if (FragCount[Z + osfc2] < maxcount)
								FragCount[Z + osfc2] = FragCount[Z + osfc2] + 1;
						}
					}
					else {
						FragSt[os3] = FS;
						FragEn[os3] = FE;
						FragScore[os3] = FSc;
						//if (FragCount[Z + osfc2] < maxcount)
							FragCount[Z + osfc2] = FragCount[Z + osfc2] + 1;
					}

					for (A = 3; A < 6; A++) {
						if (A != Z) {
							os3 = FragCount[A + osfc2] + A*os5 + osfs2;
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = -FSc;
							//if (FragCount[A + osfc2] < maxcount)
								FragCount[A + osfc2] = FragCount[A + osfc2] + 1;
						}
					}



				}
				else {
					if (FragCount[3 + osfc2] > 0) {

						for (A = 3; A < 6; A++) {
							os3 = FragCount[A + osfc2] + A*os5 + osfs2;
							if (FragScore[os3 - 1] > 0) {

								FragEn[os3 - 1] = FE;
								FragScore[os3 - 1] = FragScore[os3 - 1] + FSc;
							}
							else {
								FragSt[os3] = FS;
								FragEn[os3] = FE;
								FragScore[os3] = FSc;
								//if (FragCount[A + osfc2] < maxcount)
									FragCount[A + osfc2] = FragCount[A + osfc2] + 1;
							}
						}//Next A
					}
					else {

						for (A = 3; A < 6; A++) {
							os3 = FragCount[A + osfc2] + A*os5 + osfs2;
							FragSt[os3] = FS;
							FragEn[os3] = FE;
							FragScore[os3] = FSc;
							//if (FragCount[A + osfc2] < maxcount)
								FragCount[A + osfc2] = FragCount[A + osfc2] + 1;
						}//Next A
					}
				}

				if (Y <= 2) {
					//
						//if (FragCount[Y + osfc2] < maxcount)
							FragCount[Y + osfc2] = FragCount[Y + osfc2] + 1;
				}
				//if (FC6 < maxcount)
					FC6++;
					for (A = 0; A < 6; A++) {
						if (FragCount[A + osfc2] >= maxcount-1) {
							GoOn = 1;
							break;
						}


					}
					

				if (GoOn == 1 || FC6 >= maxcount-1)
					break;

				

			}//next X
			if (goon[e] == 1) {
				FragCount[6 + osfc2] = FC6;
				for (X = 0; X <= 6; X++) {
					os = FragCount[X + osfc2] + X*os5 + osfs2;
					FragSt[os] = lxos;
					FragEn[os] = lxos;
				}
			}

			

		}

		//omp_set_num_threads(2);
		return(1);
	}

	double MyMathFuncs::FastGC(int lseq, int Nextno, double PCO, int UBSN, int UBTP, int UBDP, short int *SeqNum, int *VarSites, unsigned char *Mask, unsigned char *TestPairs, unsigned char *DP) {
		double PV;
		int VSPos, Z, X, Y, sno, snox, snoy, tpo, dpo;
		VSPos = -1;
		sno = UBSN + 1;
		tpo = UBTP + 1;
		dpo = UBDP + 1;
		for (Z = 1; Z <=lseq; Z++){
        
			for (X = 0; X< Nextno; X++){
				snox = Z+sno*X;
				if (SeqNum[snox] != 46){
					for (Y = X + 1; Y<= Nextno; Y++){
						snoy = Z + sno*Y;
						if (SeqNum[snoy] != 46){
							if (SeqNum[snoy] != SeqNum[snox]) {
								VSPos++;
								VarSites[VSPos] = Z;
								X = Nextno;
								break;
							}
						}
					}
				}
			}
        
		}
		

		int *xy;
		int c;
		xy = (int *)calloc((Nextno+1)*Nextno, sizeof(int));
		c = 0;
		for (X = 0; X < Nextno; X++) {
			if (Mask[X] == 0) {

				for (Y = X + 1; Y <= Nextno; Y++) {
					if (Mask[Y] == 0) {
						if (TestPairs[X + Y*tpo] == 1) {
							xy[c] = X;
							xy[c+1] = Y;
							c = c + 2;
						}
					}
				}
			}
		}

		c = c - 2;
		c = c / 2;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);



#pragma omp parallel 
		{
			int Z, X, Y, dpo,tpo, sno, snox, snoy,d;


			sno = UBSN + 1;
			tpo = UBTP + 1;
			dpo = UBDP + 1;

			int FC6, SP, NDiff, A, B, InC, OutC;


			unsigned char *SS;


			int *InSt, *OutSt, *InEn, *OutEn;
			float *InSc, *OutSc;
			int NSame;
			float LTG, MissPen, MissPen2;

			float  FS, HSX, FMS, Polys, Diffs, TS;
			int MSP;
			double InHS, OutHS;

			int *OutMSP, *InMSP;
			float *OutMS, *InMS;

			double OutKMax, InKMax, D1, D2, D3, ZM, ZDel, Yy, LZ, LL0, mP, mX, InLL, OutLL, P, Q;
			double KMTL, Critval, PVM1, tHLD;
			double InKAScore, InPV;
			
			FC6 = 0;
#pragma omp for// private(X,Y,Z,A,B,SP,tHLD, FS, TS,HSX,Polys,InHS, InPV, InKAScore,Diffs,FMS, MSP NDiff, FC6, snox,snoy,SS, NSame,InLL,PVM1,Critval,,P,Q,mP,mX, LL0, LZ,ZDel, Yy, ZM, D1,D2,D3,InKMax,KMTL,LTG,MissPen,MissPen2,InSt,OutSt,InEN,OutEn,InSc,OutSc,InC,OutC,InMS, InMSP)
			for (d = 0; d <= c; d++){
				X = xy[d * 2];
				Y = xy[1 + d * 2];
			//for (X = 0; X < Nextno; X++) {
			//	if (Mask[X] == 0) {

			//		for (Y = X + 1; Y <= Nextno; Y++) {
			//			if (Mask[Y] == 0) {
			//				if (TestPairs[X + Y*tpo] == 1) {
								SS = (unsigned char*)calloc(VSPos, sizeof(unsigned char));
								SP = -1;
								NDiff = 0;
								FC6 = 0;
								for (Z = 0; Z <= VSPos; Z++) {
									snox = VarSites[Z] + X*sno;
									if (SeqNum[snox] != 46) {
										snoy = VarSites[Z] + Y*sno;
										if (SeqNum[snoy] != 46) {
											SP++;
											//VSBak[SP] = VarSites(Z)
											if (SeqNum[snox] != SeqNum[snoy]) {
												SS[SP] = 1;
												NDiff++;
											}
											else
												SS[SP] = 0;

										}
									}
								}
								if (NDiff > 2 && SP > 1) {

									NSame = SP - NDiff;
									LTG = SP * 1;//1=gcmissmatchpen or a g-scale of 1
									MissPen = (float)((int)(LTG / NDiff) + 1);
									MissPen2 = (float)((int)(LTG / NSame) + 1);


									InSt = (int*)calloc(VSPos, sizeof(int));
									OutSt = (int*)calloc(VSPos, sizeof(int));
									InEn = (int*)calloc(VSPos, sizeof(int));
									OutEn = (int*)calloc(VSPos, sizeof(int));
									InSc = (float*)calloc(VSPos, sizeof(float));
									OutSc = (float*)calloc(VSPos, sizeof(float));

									InC = -1;
									OutC = -1;
									for (Z = 0; Z <= SP; Z++) {
										if (SS[Z] == 0) {
											InC++;
											InSt[InC] = Z;

											for (A = Z + 1; A < SP * 2; A++) {
												if (A > SP)
													B = A - SP - 1;
												else
													B = A;

												if (SS[B] == 1)
													break;
											}
											Z = A - 1;
											InEn[InC] = A - 1;
											InSc[InC] = InEn[InC] - InSt[InC] + 1;

											OutC++;
											OutSt[OutC] = InSt[InC];
											OutEn[OutC] = InEn[InC];
											OutSc[OutC] = -InSc[InC];
										}

										else {

											OutC++;
											OutSt[OutC] = Z;
											for (A = Z + 1; A < SP * 2; A++) {
												if (A > SP)
													B = A - SP - 1;
												else
													B = A;

												if (SS[B] == 0)
													break;
											}
											Z = A - 1;
											OutEn[OutC] = A - 1;
											OutSc[OutC] = OutEn[OutC] - OutSt[OutC] + 1;

											InC++;
											InSt[InC] = OutSt[OutC];
											InEn[InC] = OutEn[OutC];
											InSc[InC] = -OutSc[OutC];
										}
										FC6++;
										if (OutC >= SP - 1)
											break;

										if (InC >= SP - 1)
											break;



									}



									InMS = (float*)calloc(VSPos, sizeof(float));
									//OutMS = (float*)calloc(VSPos, sizeof(float));
									InMSP = (int*)calloc(VSPos, sizeof(int));
									//OutMSP = (int*)calloc(VSPos, sizeof(int));
									HSX = 0;
									for (Z = 0; Z <= InC; Z++) {

										if (InSc[Z] > 0) {
											Polys = InSc[Z];
											FMS = Polys;
											Diffs = 0.0;
											MSP = Z;
											for (B = Z + 1; B <= InC * 2; B++) {
												if (B > InC)
													A = B - InC - 1;
												else
													A = B;

												FS = (float)(InSc[A]);
												if (FS <= 0) {
													FS = FS*-1;
													Diffs = Diffs + FS;
												}
												Polys = Polys + FS;
												TS = ((Polys - Diffs) - (Diffs * MissPen));
												if (TS < 0)
													break;
												//TS = (int)(TS + 0.5);
												if (TS >= FMS) {
													FMS = TS;
													MSP = A;
												}
											}
											InMSP[Z] = MSP;
											InMS[Z] = FMS;

											if (FMS > HSX)
												HSX = FMS;

										}
										else
											InMS[Z] = 0;

									}
									InHS = (double)(HSX);


									InLL = 0;

									if (InHS > 3) {
										P = (double)((double)(NDiff) / (double)(SP));
										Q = (double)(1) - P;
										mP = (double)(MissPen)* P;
										mX = (double)(MissPen);
										LL0 = log(mP / Q);
										LL0 = LL0 / (mX + 1);
										LZ = exp(2 * LL0);
										ZDel = 1.0;
										Yy = 1.0;
										while (fabs(ZDel) > 0.000001 || fabs(Yy) > 0.000001) {
											ZM = pow(LZ, -mX);
											Yy = Q * LZ + P * ZM - 1;
											ZDel = Yy / (Q - mP * ZM / Z);
											LZ = LZ - ZDel;


										}
										InLL = log(LZ);
										D1 = exp(InLL);
										D1 = D1 - 1;
										D2 = -(mX + 1) * InLL;
										D3 = exp(D2);
										InKMax = D1 * (Q - (mP * D3));

										if (InKMax > 0) {
											KMTL = InKMax * SP;
											KMTL = log(KMTL);
											PVM1 = (1 - PCO);
											if (PVM1 > 0) {
												PVM1 = -log(PVM1);
												Critval = (KMTL + PVM1) / InLL;
												if (Critval < 4)
													Critval = 4;
											}
											else
												Critval = 4;

										}
										else
											Critval = 4;

										if (InHS > Critval) {
											

											InKAScore = InLL * InHS - KMTL;
											if (InKAScore >= 30)
												InKAScore = 30;
											if (InKAScore > 0 && InKAScore <= 30) {
												tHLD = exp(-InKAScore);
												InPV = 1 - exp(-tHLD);

											}
											if (InPV <= PCO) {
												DP[X + Y*dpo] = 1;
												DP[Y + X*dpo] = 1;
											}
										}

									}

									/*if (DP[X + Y*dpo] == 0){
										OutLL = 0;
										HSX = 0;
										for (Z = 0; Z<= OutC; Z++){
											if (OutSc[Z] > 0){
												Polys = OutSc[Z];
												FMS = Polys;
												Diffs = 0;
												MSP = Z;
												for (B = Z + 1; B <= OutC * 2; B++) {
													if (B > OutC)
														A = B - OutC - 1;
													else
														A = B;

													FS = OutSc[A];
													if (FS <= 0) {
														FS = -FS;
														Diffs = Diffs + FS;
													}
													Polys = Polys + FS;
													TS = ((Polys - Diffs) - (Diffs * MissPen2));
													if (TS < 0)
														break;
													TS = (int)(TS + 0.5);
													if (TS >= FMS){
														FMS = TS;
														MSP = A;
													}
												}
												OutMSP[Z] = MSP;
												OutMS[Z] = FMS;
												if (FMS > HSX)
													HSX = FMS;

											}
											else
												OutMS[Z] = 0;

										}
										OutHS = HSX;

										if (OutHS > 3){
											//P = NSame / SP;
											P = (double)((double)(NSame) / (double)(SP));
											Q = 1 - P;
											mP = MissPen2 * P;
											mX = MissPen2;
											LL0 = log(mP / Q);
											LL0 = LL0 / (mX + 1);
											LZ = exp(2 * LL0);
											ZDel = 1;
											Yy = 1;
											while (fabs(ZDel) > 0.000001 || fabs(Yy) > 0.000001) {
												ZM = pow(LZ, -mX);
												Yy = Q * LZ + P * ZM - 1;
												ZDel = Yy / (Q - mP * ZM / Z);
												LZ = LZ - ZDel;
											}
											OutLL = log(LZ);

											D1 = exp(OutLL);
											D1 = D1 - 1;
											D2 = -(mX + 1) * OutLL;
											D3 = exp(D2);
											OutKMax = D1 * (Q - (mP * D3));
											if (OutKMax > 0){
												KMTL = OutKMax * SP;
												KMTL = log(KMTL);
												PVM1 = (1 - PCO);
												if (PVM1 > 0) {
													PVM1 = -log(PVM1);
													Critval = (KMTL + PVM1) / OutLL;
													if (Critval < 4)
														Critval = 4;
												}
												else
													Critval = 4;

											}

											else
												Critval = 4;

											if (OutHS > Critval){
												double OutKAScore, OutPV;
												OutKAScore = OutLL * OutHS - KMTL;
												if (OutKAScore >= 30)
													OutKAScore = 30;
												if (OutKAScore > 0 && OutKAScore <= 30) {
													tHLD = exp(-OutKAScore);
													OutPV = 1 - exp(-tHLD);

												}
												if (OutPV <= PCO) {
													DP[X + Y*dpo] = 1;
													DP[Y + X*dpo] = 1;
												}
											}
										}
									}*/
									free(InSc);
									free(OutSc);
									free(InEn);
									free(OutEn);
									free(InSt);
									free(OutSt);
									free(InMS);
									//free(OutMS);
									free(InMSP);
									//free(OutMSP);
								}
								free(SS);
							//}
						//}
					//}
				//}
			}
			
		}
		free(xy);
		omp_set_num_threads(2);
		return (1);
	}


	int MyMathFuncs::cleanss(int y,int UBSS, char *SubSeq) {
		int x, z, off1, off2, off6;
		off1 = (UBSS + 1);
		off2 = (UBSS + 1)*2;
		off6 = (UBSS + 1) * 6;
		for (x = 0; x <= y; x++)
			SubSeq[x] = 0;
		for (x = 0; x <= y; x++)
			SubSeq[x+off1] = 0;
		for (x = 0; x <= y; x++)
			SubSeq[x+off2] = 0;
		for (x = 0; x <= y; x++)
			SubSeq[x+off6] = 0;
		
		return(1);
	}
	double MyMathFuncs::GCXoverDP2(double *BQPV, int ubcs, unsigned char *cs, int ubfss, unsigned char *fssgc, int MCFlag, int UBPV, double *PVals,double LowestProb, int MCCorrection,int ShortOutFlag, int CircularFlag, int GCDimSize, int lenstrainseq0, short int GCMissmatchPen, char GCIndelFlag, int Seq1, int Seq2, int Seq3, int UBFST, int *FragSt, int *FragEn, int UBFS, int *FragScore, int UBSS, char *SubSeq, int UBMSP, int *MaxScorePos, int UBFMS, int *FragMaxScore, int *HighEnough) {
		int x, LenXoverSeq, MaxDiff, MinDiff, GoOn, dummy;
		int *NDiff, *HiFragScore, *FragCount;
		double *LL, *KMax;
		float LTG;
		double *MissPen, MaxScore, PCO, *Critval;
		
		NDiff = (int*)calloc(8, sizeof(int));
		LL = (double*)calloc(8, sizeof(double));
		KMax = (double*)calloc(8, sizeof(double));
		Critval = (double*)calloc(8, sizeof(double));
		FragCount = (int*)calloc(8, sizeof(int));
		/*if (GCIndelFlag == 1)
			LenXoverSeq = FindSubSeqGCAP(GCIndelFlag, lenstrainseq0, Seq1, Seq2, Seq3, SeqNum(0, 0), SubSeq, XPosDiff(0), XDiffPos(0), NDiff);
		else*/
			LenXoverSeq = FindSubSeqGCAP6(ubcs, cs, ubfss, fssgc,GCIndelFlag, lenstrainseq0, Seq1, Seq2, Seq3,  SubSeq, NDiff);


		if (NDiff[0] == LenXoverSeq || NDiff[1] == LenXoverSeq || NDiff[2] == LenXoverSeq) {
			dummy = cleanss(LenXoverSeq, UBSS, SubSeq);
			free(KMax);
			free(LL);
			free(NDiff);
			free(FragCount);
			free(Critval);
			return(0);
		}


		HiFragScore = (int*)calloc(7, sizeof(int));

		SubSeq[LenXoverSeq + 1] = 0;
		SubSeq[LenXoverSeq + 1 + UBSS + 1] = 0;
		SubSeq[LenXoverSeq + 1 + 2 * (UBSS + 1)] = 0;
		SubSeq[LenXoverSeq + 1 + 6 * (UBSS + 1)] = 0;

		//for outer frags (ie matches instead of differences)
		NDiff[3] = NDiff[0] + NDiff[1]; //seq1
		NDiff[4] = NDiff[0] + NDiff[2]; //seq2
		NDiff[5] = NDiff[1] + NDiff[2]; //seq3

			//for inner fargs (ie genuine differences)
		NDiff[0] = LenXoverSeq - NDiff[0];
		NDiff[1] = LenXoverSeq - NDiff[1];
		NDiff[2] = LenXoverSeq - NDiff[2];

		LTG = LenXoverSeq * GCMissmatchPen;

		MissPen = (double*)calloc(7, sizeof(double));
			
		MissPen[0] = (int)(LTG / (float)(NDiff[0])) + 1;
		MissPen[1] = (int)(LTG / (float)(NDiff[1])) + 1; 
		MissPen[2] = (int)(LTG / (float)(NDiff[2])) + 1;


		MaxDiff = 0;
		MinDiff = lenstrainseq0;

		for (x = 0; x <= 5; x++) {
			if(MinDiff > NDiff[x])
				MinDiff = NDiff[x];
			if (MaxDiff < NDiff[x])
				MaxDiff = NDiff[x];
		}

	if (MinDiff < 3 && MaxDiff > MinDiff * 10) {
		free(KMax);
		free(LL);
		free(MissPen);
		free(HiFragScore);
		free(NDiff);
		free(FragCount);
		free(Critval);
		dummy = cleanss(LenXoverSeq, UBSS, SubSeq);
		return(0);
	}


	if (NDiff[3] == 0)
		NDiff[3] == 1;
	if (NDiff[4] == 0)
		NDiff[4] = 1;
	if (NDiff[5] == 0)
		NDiff[5] = 1;

	MissPen[3] = (int)(LTG / (int)(NDiff[3])) + 1;
	MissPen[4] = (int)(LTG / (int)(NDiff[4])) + 1;
	MissPen[5] = (int)(LTG / (int)(NDiff[5])) + 1;
	
	GoOn = GetFragsP(CircularFlag, LenXoverSeq, lenstrainseq0, GCDimSize, SubSeq, FragSt, FragEn, FragScore, FragCount);

		if (GoOn == 0 && ShortOutFlag != 3) {
			free(KMax);
			free(LL);
			free(MissPen);
			free(HiFragScore);
			free(NDiff);
			free(FragCount);
			free(Critval);
			dummy = cleanss(LenXoverSeq, UBSS, SubSeq);
			return(0);
		}

	dummy =GetMaxFragScoreP(LenXoverSeq, GCDimSize, CircularFlag, GCMissmatchPen, MissPen, MaxScorePos, FragMaxScore, FragScore, FragCount, HiFragScore);

	/*if (Seq1 == 1 && Seq2 == 17 && Seq3 == 34)
		x = 0;*/
	for (x = 0; x <= 5; x++) {
		if (HiFragScore[x] > 3)
			HighEnough[x] = 1;
		else
			HighEnough[x] = 0;
	}




			
	GoOn = CalcKMaxP(GCMissmatchPen, LenXoverSeq, MCFlag, MCCorrection, LowestProb, &PCO, HiFragScore, Critval, MissPen, LL, KMax, NDiff, HighEnough);

	if (GoOn == 0) {
		free(KMax);
		free(LL);
		free(MissPen);
		free(HiFragScore);
		free(NDiff);
		free(FragCount);
		free(Critval);
		dummy = cleanss(LenXoverSeq, UBSS, SubSeq);
		return(0);
	}




	MaxScore = GCCalcPValP2(GCDimSize, LenXoverSeq, FragMaxScore, PVals, FragCount, KMax, LL, HighEnough, Critval);

	
	free(KMax);
	free(LL);
	free(MissPen);
	free(HiFragScore);
	free(NDiff);
	free(FragCount);
	free(Critval);
	
	dummy = cleanss(LenXoverSeq, UBSS, SubSeq);
	*BQPV = MaxScore;
	if (MaxScore <= PCO)
		return(1);
	else
		return(0);
		
	}


	double MyMathFuncs::GCXoverDP(int MCFlag, int UBPV, double *PVals, double LowestProb, int MCCorrection, int ShortOutFlag, int CircularFlag, int GCDimSize, int lenstrainseq0, short int GCMissmatchPen, char GCIndelFlag, int Seq1, int Seq2, int Seq3, int UBFST, int *FragSt, int *FragEn, int UBFS, int *FragScore, short int *SeqNum, int UBSS, char *SubSeq, int UBMSP, int *MaxScorePos, int UBFMS, int *FragMaxScore, int *HighEnough) {
		int x, LenXoverSeq, MaxDiff, MinDiff, GoOn, dummy;
		int *NDiff, *HiFragScore, *FragCount;
		double *LL, *KMax;
		float LTG;
		double *MissPen, MaxScore, PCO, *Critval;

		NDiff = (int*)calloc(8, sizeof(int));
		LL = (double*)calloc(8, sizeof(double));
		KMax = (double*)calloc(8, sizeof(double));
		Critval = (double*)calloc(8, sizeof(double));
		FragCount = (int*)calloc(8, sizeof(int));

		LenXoverSeq = FindSubSeqGCAP5(GCIndelFlag, lenstrainseq0, Seq1, Seq2, Seq3, SeqNum, SubSeq, NDiff);




		if (NDiff[0] == LenXoverSeq || NDiff[1] == LenXoverSeq || NDiff[2] == LenXoverSeq) {
			free(NDiff);
			free(LL);
			free(KMax);
			free(FragCount);
			return(0);
		}


		HiFragScore = (int*)calloc(7, sizeof(int));

		SubSeq[LenXoverSeq + 1] = 0;
		SubSeq[LenXoverSeq + 1 + UBSS + 1] = 0;
		SubSeq[LenXoverSeq + 1 + 2 * (UBSS + 1)] = 0;
		SubSeq[LenXoverSeq + 1 + 6 * (UBSS + 1)] = 0;

		//for outer frags (ie matches instead of differences)
		NDiff[3] = NDiff[0] + NDiff[1]; //seq1
		NDiff[4] = NDiff[0] + NDiff[2]; //seq2
		NDiff[5] = NDiff[1] + NDiff[2]; //seq3

										//for inner fargs (ie genuine differences)
		NDiff[0] = LenXoverSeq - NDiff[0];
		NDiff[1] = LenXoverSeq - NDiff[1];
		NDiff[2] = LenXoverSeq - NDiff[2];

		LTG = LenXoverSeq * GCMissmatchPen;

		MissPen = (double*)calloc(7, sizeof(double));

		MissPen[0] = (int)(LTG / (float)(NDiff[0])) + 1;
		MissPen[1] = (int)(LTG / (float)(NDiff[1])) + 1;
		MissPen[2] = (int)(LTG / (float)(NDiff[2])) + 1;


		MaxDiff = 0;
		MinDiff = lenstrainseq0;

		for (x = 0; x <= 5; x++) {
			if (MinDiff > NDiff[x])
				MinDiff = NDiff[x];
			if (MaxDiff < NDiff[x])
				MaxDiff = NDiff[x];
		}

		if (MinDiff < 3 && MaxDiff > MinDiff * 10) {
			free(KMax);
			free(LL);
			free(MissPen);
			free(HiFragScore);
			free(NDiff);
			free(FragCount);
			return(0);
		}


		if (NDiff[3] == 0)
			NDiff[3] = 1;
		if (NDiff[4] == 0)
			NDiff[4] = 1;
		if (NDiff[5] == 0)
			NDiff[5] = 1;

		MissPen[3] = (int)(LTG / (int)(NDiff[3])) + 1;
		MissPen[4] = (int)(LTG / (int)(NDiff[4])) + 1;
		MissPen[5] = (int)(LTG / (int)(NDiff[5])) + 1;

		GoOn = GetFragsP(CircularFlag, LenXoverSeq, lenstrainseq0, GCDimSize, SubSeq, FragSt, FragEn, FragScore, FragCount);

		if (GoOn == 0 && ShortOutFlag != 3) {
			free(KMax);
			free(LL);
			free(MissPen);
			free(HiFragScore);
			free(NDiff);
			free(FragCount);
			return(0);
		}

		dummy = GetMaxFragScoreP(LenXoverSeq, GCDimSize, CircularFlag, GCMissmatchPen, MissPen, MaxScorePos, FragMaxScore, FragScore, FragCount, HiFragScore);


		for (x = 0; x <= 5; x++) {
			if (HiFragScore[x] > 3)
				HighEnough[x] = 1;
			else
				HighEnough[x] = 0;
		}





		GoOn = CalcKMaxP(GCMissmatchPen, LenXoverSeq, MCFlag, MCCorrection, LowestProb, &PCO, HiFragScore, Critval, MissPen, LL, KMax, NDiff, HighEnough);

		if (GoOn == 0) {
			free(KMax);
			free(LL);
			free(MissPen);
			free(HiFragScore);
			free(NDiff);
			free(FragCount);
			return(0);
		}




		MaxScore = GCCalcPValP2(GCDimSize, LenXoverSeq, FragMaxScore, PVals, FragCount, KMax, LL, HighEnough, Critval);


		free(KMax);
		free(LL);
		free(MissPen);
		free(HiFragScore);
		free(NDiff);
		free(FragCount);

		if (MaxScore <= PCO)
			return(1);
		else
			return(0);

	}

	double MyMathFuncs::GCCalcPValP2(int lseq, long LXover, int *FragMaxScore, double *PVals, int *FragCount, double *KMax, double *LL, int *highenough, double *critval) {
		int X, Y, os, os2;
		double MaxScore, LenXoverSeq, THld;
		float KAScore, LKLen, warn;
		os = lseq + 1;
		LenXoverSeq = (double)(LXover);
		MaxScore = 10;
		for (X = 0; X <= 5; X++) {
			if (highenough[X] == 1) {

				//10^2=100,log100 = 2; exp (10) = e^10 
				//only calculating scores over a critical maximum will massively speed this up
				if (KMax[X] > 0) {

					//1-exp(-exp(-(ll(x)*score-log(kmax*lenxoverseq))))
					// work out which score corresponds with a particular p val
					//pval = 1-exp(-exp(-(ll(x)*score-log(kmax*lenxoverseq))))
					//1-pval = exp(-exp(-(ll(x)*score-log(kmax*lenxoverseq))))
					//-log(1-pval) = exp(-(ll(x)*score-log(kmax*lenxoverseq)))
					//-log(-log(1-pval)) = ll(x)*score-log(kmax*lenxoverseq)
					//log(kmax*lenxoverseq)-log(-log(1-pval))=ll(x)*score
					//(log(kmax*lenxoverseq)-log(-log(1-pval)))/ll(x) = score
					LKLen = (float)(log(KMax[X] * LenXoverSeq));
					for (Y = 0; Y <= FragCount[X]; Y++) {
						os2 = Y + X*os;
						if (FragMaxScore[os2] > critval[X]) {
							KAScore = (float)((LL[X] * FragMaxScore[os2]) - LKLen);
							if (KAScore > 0) {
								if (KAScore < 32) {

									THld = exp((double)(-KAScore));
									PVals[os2] = 1 - exp(-THld);
								}
								else {
									warn = 0;
									if (KAScore > 700) {

										warn = KAScore;
										KAScore = 701;

									}

									THld = exp((double)(-KAScore));
									if (warn != 0) {
										KAScore = (float)(warn - 700);
										THld = THld /(double)(KAScore);
									}
									PVals[os2] = THld;
								}
							}
							else
								PVals[os2] = 1;//THld;


							if (PVals[os2] < MaxScore)
								MaxScore = PVals[os2];
						}
						else
							PVals[os2] = 1;
					}
				}
				else
					highenough[X] = 0;
			}
		}
		return (MaxScore);
	}

	int MyMathFuncs::AEFirstRDP(int Seq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DP, unsigned char *DP2, int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseq, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS,  unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact) {
	
		
	int redonum;	
		
//#pragma omp parallel 
		
//		{
			int Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo,  dp12, dp13, dp23, pv12, pv23, pv13;
			unsigned char  NewOneFound;
			double BQPV;
			int UBXOHN, UBXSNW, FRC;
			redonum = 0;
			rlo = UBRL + 1;
			dpo = UBDP + 1;
			//AH = (int*)calloc(4, sizeof(int));

			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XoverSeqNum(Len(StrainSeq(0)), 2)
			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XOverHomologyNum(Len(StrainSeq(0)) + XoverWindow * 2, 2)
			UBXOHN = LenStrainseq0 + XoverWindow * 2;
			UBXSNW = LenStrainseq0 + (int)(XOverWindowX / 2) * 2;
			short int *XoverSeqNum;
			char *XoverSeqNumW;
			int *XOverHomologyNum;

			XoverSeqNum = (short int*)calloc((LenStrainseq0 + 1) * 3, sizeof(short int));
			XoverSeqNumW = (char*)calloc((UBXSNW + 1) * 3, sizeof(char));
			XOverHomologyNum = (int*)calloc((UBXOHN + 1) * 3, sizeof(int));
//#pragma omp for
			for (Seq2 = oNextno + 1; Seq2 <= NextNo - 1; Seq2++) {
				if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
					s12o = (int)(Seq1 + (UBPV + 1)*Seq2);
					pv12 = (int)(PermValid[s12o]);
					if (PermDiffs[s12o] > MinDIffs && pv12 > MinSeqSize) {
						dp12 = (int)(DP[Seq1 + dpo*Seq2]);
						for (Seq3 = Seq2 + 1; Seq3 <= NextNo; Seq3++) {
							s13o = (int)(Seq1 + (UBPV + 1)*Seq3);
							s23o = (int)(Seq2 + (UBPV + 1)*Seq3);
							pv13 = (int)(PermValid[s13o]);
							pv23 = (int)(PermValid[s23o]);
							dpos3 = Seq3*dpo;
							dp13 = (int)(DP[Seq1 + dpos3]);
							dp23 = (int)(DP[Seq2 + dpos3]);
							if (Seq3 <= UBPV && tMaskseq[Seq3] == 0 && (dp23 + dp13 + dp12) > 0) {
								if (pv13 > MinSeqSize && pv23 > MinSeqSize) {
									if (PermDiffs[s23o] > (float)(MinDIffs)) {
										NewOneFound = 0;
										BQPV = 0;

										//int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0,int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, int UBXSNW, unsigned char *CS,short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag,int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact
										FRC = FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
										if (FRC == 1) {
											//Call XOver(Distance(), XPosDiff(), XDiffPos(), CurrentXOver(), XoverList(), SeqNum(), Seq1, Seq2, Seq3, 0)
//#pragma omp critical
//											{
												redonum++;
												RL[redonum*rlo] = Seq1;
												RL[1 + redonum*rlo] = Seq2;
												RL[2 + redonum*rlo] = Seq3;
											//}
											NewOneFound = 1;

										}

										if (BQPV > 0) {
											if (BQPV < SubThresh)
												NewOneFound = 1;

										}




										if (NewOneFound == 1) {
//#pragma omp critical
//											{

												DP2[Seq1 + dpo*Seq3] = 1;
												DP2[Seq3 + dpo*Seq1] = 1;
												DP2[Seq2 + dpo*Seq1] = 1;
												DP2[Seq1 + dpo*Seq2] = 1;
//											}
										}
									}
								}
							}
						}
					}
				}
			}

			free(XoverSeqNum);
			free(XoverSeqNumW);
			free(XOverHomologyNum);
//		}
	return(redonum);
}

	int MyMathFuncs::AEFirstRDP2(int Seq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DP, unsigned char *DP2, int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseq, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact) {


		int redonum;
		
		redonum = 0;
		int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, pv12, pv23, pv13, dsa,dsb,dsc,dsd,dse,dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int UBXOHN, UBXSNW, FRC;

			
			dsb = (UBCS + 1)*(NextNo + 1);
			dsa = (NextNo + 1)*(NextNo + 1);
			dsc = 4 * (UBFSSRDP + 1)*(UBFSSRDP + 1)*(UBFSSRDP + 1);
			dsd = (UBFact3x3 + 1)*(UBFact3x3 + 1)*(UBFact3x3 + 1);
			dse = 172;
			dsf = 172 * 172 * 51;
			rlo = UBRL + 1;
			dpo = UBDP + 1;

			UBXOHN = LenStrainseq0 + XoverWindow * 2;
			UBXSNW = LenStrainseq0 + (int)(XOverWindowX / 2) * 2;
			short int *XoverSeqNum;
			char *XoverSeqNumW;
			int *XOverHomologyNum;

			XoverSeqNum = (short int*)calloc((LenStrainseq0 + 1) * 3, sizeof(short int));
			XoverSeqNumW = (char*)calloc((UBXSNW + 1) * 3, sizeof(char));
			XOverHomologyNum = (int*)calloc((UBXOHN + 1) * 3, sizeof(int));
			
			//unsigned char *DP;
			//unsigned char *CS;
			////float *PermDiffs, *PermValid, *Distance, *TreeDistance;
			//unsigned char *tMaskseq,*FSSRDP; 
			
			//double *Fact3X3, *Fact, *ProbEstimate;*/

			/*PermDiffs = (float*)calloc(dsa, sizeof(float));
			PermValid = (float*)calloc(dsa, sizeof(float));
			Distance = (float*)calloc(dsa, sizeof(float));
			TreeDistance = (float*)calloc(dsa, sizeof(float));*/
			//DP = (unsigned char*)calloc(dsa, sizeof(unsigned char));
			//CS = (unsigned char*)calloc(dsb, sizeof(unsigned char));
			//FSSRDP = (unsigned char*)calloc(dsc, sizeof(unsigned char));
			//tMaskseq = (unsigned char*)calloc(NextNo+1, sizeof(unsigned char));
			/*Fact3X3 = (double*)calloc(dsd, sizeof(double));
			Fact = (double*)calloc(dse, sizeof(double));
			ProbEstimate = (double*)calloc(dsf, sizeof(double));*/

//#pragma omp critical
//{
				//for (x = 0; x <= dsf; x++)
				//	ProbEstimate[x] = ProbEstimateO[x];
				//for (x = 0; x <= dsd; x++)
				//	Fact3X3[x] = Fact3X3O[x];
				//for (x = 0; x <= dse; x++)
				//	Fact[x] = FactO[x];
				/*for (x = 0; x <= NextNo + 1; x++)
					tMaskseq[x] = tMaskseqO[x];
				for (x = 0; x <= dsb; x++)
					CS[x] = CSO[x];
				for (x = 0; x <= dsc; x++)
					FSSRDP[x] = FSSRDPO[x];*/
				//for (x = 0; x <= dsa; x++)
				//	DP[x] = DPO[x];
				//for (x = 0; x <= dsa; x++)
				//	PermDiffs[x] = PermDiffsO[x];
				//for (x = 0; x <= dsa; x++)
				//	PermValid[x] = PermValidO[x];
				//for (x = 0; x <= dsa; x++)
				//	Distance[x] = DistanceO[x];
				//for (x = 0; x <= dsa; x++)
				//	TreeDistance[x] = TreeDistanceO[x];


//}
			
			
			//AH = (int*)calloc(4, sizeof(int));

			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XoverSeqNum(Len(StrainSeq(0)), 2)
			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XOverHomologyNum(Len(StrainSeq(0)) + XoverWindow * 2, 2)
			
#pragma omp for
			for (Seq2 = oNextno + 1; Seq2 <= NextNo - 1; Seq2++) {
				if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
					s12o = (int)(Seq1 + (UBPV + 1)*Seq2);
					pv12 = (int)(PermValid[s12o]);
					if (PermDiffs[s12o] > MinDIffs && pv12 > MinSeqSize) {
						dp12 = (int)(DP[Seq1 + dpo*Seq2]);
						for (Seq3 = Seq2 + 1; Seq3 <= NextNo; Seq3++) {
							s13o = (int)(Seq1 + (UBPV + 1)*Seq3);
							s23o = (int)(Seq2 + (UBPV + 1)*Seq3);
							pv13 = (int)(PermValid[s13o]);
							pv23 = (int)(PermValid[s23o]);
							dpos3 = Seq3*dpo;
							dp13 = (int)(DP[Seq1 + dpos3]);
							dp23 = (int)(DP[Seq2 + dpos3]);
							if (Seq3 <= UBPV && tMaskseq[Seq3] == 0 && (dp23 + dp13 + dp12) > 0) {
								if (pv13 > MinSeqSize && pv23 > MinSeqSize) {
									if (PermDiffs[s23o] > (float)(MinDIffs)) {
										NewOneFound = 0;
										BQPV = 0;

										//int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0,int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, int UBXSNW, unsigned char *CS,short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag,int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact
										FRC = FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
										if (FRC == 1) {
											//Call XOver(Distance(), XPosDiff(), XDiffPos(), CurrentXOver(), XoverList(), SeqNum(), Seq1, Seq2, Seq3, 0)
#pragma omp critical
											{
												redonum++;
												RL[redonum*rlo] = Seq1;
												RL[1 + redonum*rlo] = Seq2;
												RL[2 + redonum*rlo] = Seq3;
											}
											NewOneFound = 1;

										}

										if (BQPV > 0) {
											if (BQPV < SubThresh)
												NewOneFound = 1;

										}




										if (NewOneFound == 1) {
#pragma omp critical
											{

												DP2[Seq1 + dpo*Seq3] = 1;
												DP2[Seq3 + dpo*Seq1] = 1;
												DP2[Seq2 + dpo*Seq1] = 1;
												DP2[Seq1 + dpo*Seq2] = 1;
											}
										}
									}
								}
							}
						}
					}
				}
			}

			free(XoverSeqNum);
			free(XoverSeqNumW);
			free(XOverHomologyNum);

			//free(PermDiffs);
			//free(PermValid);
			//free(Distance);
			//free(TreeDistance);
			//free(DP);
			//free(CS);
			//free(FSSRDP);
			//free(tMaskseq);
			//free(Fact3X3);
			//free(Fact);
			//free(ProbEstimate);

		}
		omp_set_num_threads(2);
		return(redonum);
	}

	int MyMathFuncs::AEFirstRDP3(int oSeq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DP, unsigned char *DP2, int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseq, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact) {


		int redonum;
		int procs;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int UBXOHN, UBXSNW, FRC;


			dsb = (UBCS + 1)*(NextNo + 1);
			dsa = (NextNo + 1)*(NextNo + 1);
			dsc = 4 * (UBFSSRDP + 1)*(UBFSSRDP + 1)*(UBFSSRDP + 1);
			dsd = (UBFact3x3 + 1)*(UBFact3x3 + 1)*(UBFact3x3 + 1);
			dse = 172;
			dsf = 172 * 172 * 51;
			rlo = UBRL + 1;
			dpo = UBDP + 1;

			UBXOHN = LenStrainseq0 + XoverWindow * 2;
			UBXSNW = LenStrainseq0 + (int)(XOverWindowX / 2) * 2;
			short int *XoverSeqNum;
			char *XoverSeqNumW;
			int *XOverHomologyNum;

			XoverSeqNum = (short int*)calloc((LenStrainseq0 + 1) * 3, sizeof(short int));
			XoverSeqNumW = (char*)calloc((UBXSNW + 1) * 3, sizeof(char));
			XOverHomologyNum = (int*)calloc((UBXOHN + 1) * 3, sizeof(int));

			//unsigned char *DP;
			//unsigned char *CS;
			////float *PermDiffs, *PermValid, *Distance, *TreeDistance;
			//unsigned char *tMaskseq,*FSSRDP; 

			//double *Fact3X3, *Fact, *ProbEstimate;*/

			/*PermDiffs = (float*)calloc(dsa, sizeof(float));
			PermValid = (float*)calloc(dsa, sizeof(float));
			Distance = (float*)calloc(dsa, sizeof(float));
			TreeDistance = (float*)calloc(dsa, sizeof(float));*/
			//DP = (unsigned char*)calloc(dsa, sizeof(unsigned char));
			//CS = (unsigned char*)calloc(dsb, sizeof(unsigned char));
			//FSSRDP = (unsigned char*)calloc(dsc, sizeof(unsigned char));
			//tMaskseq = (unsigned char*)calloc(NextNo+1, sizeof(unsigned char));
			/*Fact3X3 = (double*)calloc(dsd, sizeof(double));
			Fact = (double*)calloc(dse, sizeof(double));
			ProbEstimate = (double*)calloc(dsf, sizeof(double));*/

			//#pragma omp critical
			//{
			//for (x = 0; x <= dsf; x++)
			//	ProbEstimate[x] = ProbEstimateO[x];
			//for (x = 0; x <= dsd; x++)
			//	Fact3X3[x] = Fact3X3O[x];
			//for (x = 0; x <= dse; x++)
			//	Fact[x] = FactO[x];
			/*for (x = 0; x <= NextNo + 1; x++)
			tMaskseq[x] = tMaskseqO[x];
			for (x = 0; x <= dsb; x++)
			CS[x] = CSO[x];
			for (x = 0; x <= dsc; x++)
			FSSRDP[x] = FSSRDPO[x];*/
			//for (x = 0; x <= dsa; x++)
			//	DP[x] = DPO[x];
			//for (x = 0; x <= dsa; x++)
			//	PermDiffs[x] = PermDiffsO[x];
			//for (x = 0; x <= dsa; x++)
			//	PermValid[x] = PermValidO[x];
			//for (x = 0; x <= dsa; x++)
			//	Distance[x] = DistanceO[x];
			//for (x = 0; x <= dsa; x++)
			//	TreeDistance[x] = TreeDistanceO[x];


			//}


			//AH = (int*)calloc(4, sizeof(int));

			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XoverSeqNum(Len(StrainSeq(0)), 2)
			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XOverHomologyNum(Len(StrainSeq(0)) + XoverWindow * 2, 2)

#pragma omp for
			for (Seq1 = 0; Seq1 <= oNextno; Seq1++){
				for (Seq2 = oNextno + 1; Seq2 <= NextNo - 1; Seq2++) {
					if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
						s12o = (int)(Seq1 + (UBPV + 1)*Seq2);
						pv12 = (int)(PermValid[s12o]);
						if (PermDiffs[s12o] > MinDIffs && pv12 > MinSeqSize) {
							dp12 = (int)(DP[Seq1 + dpo*Seq2]);
							for (Seq3 = Seq2 + 1; Seq3 <= NextNo; Seq3++) {
								s13o = (int)(Seq1 + (UBPV + 1)*Seq3);
								s23o = (int)(Seq2 + (UBPV + 1)*Seq3);
								pv13 = (int)(PermValid[s13o]);
								pv23 = (int)(PermValid[s23o]);
								dpos3 = Seq3*dpo;
								dp13 = (int)(DP[Seq1 + dpos3]);
								dp23 = (int)(DP[Seq2 + dpos3]);
								if (Seq3 <= UBPV && tMaskseq[Seq3] == 0 && (dp23 + dp13 + dp12) > 0) {
									if (pv13 > MinSeqSize && pv23 > MinSeqSize) {
										if (PermDiffs[s23o] > (float)(MinDIffs)) {
											NewOneFound = 0;
											BQPV = 0;

											//int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0,int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, int UBXSNW, unsigned char *CS,short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag,int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact
											FRC = FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
											if (FRC == 1) {
												//Call XOver(Distance(), XPosDiff(), XDiffPos(), CurrentXOver(), XoverList(), SeqNum(), Seq1, Seq2, Seq3, 0)
	#pragma omp critical
												{
													redonum++;
													RL[redonum*rlo] = Seq1;
													RL[1 + redonum*rlo] = Seq2;
													RL[2 + redonum*rlo] = Seq3;
												}
												NewOneFound = 1;

											}

											if (BQPV > 0) {
												if (BQPV < SubThresh)
													NewOneFound = 1;

											}




											if (NewOneFound == 1) {
	#pragma omp critical
												{

													DP2[Seq1 + dpo*Seq3] = 1;
													DP2[Seq3 + dpo*Seq1] = 1;
													DP2[Seq2 + dpo*Seq1] = 1;
													DP2[Seq1 + dpo*Seq2] = 1;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			free(XoverSeqNum);
			free(XoverSeqNumW);
			free(XOverHomologyNum);

			//free(PermDiffs);
			//free(PermValid);
			//free(Distance);
			//free(TreeDistance);
			//free(DP);
			//free(CS);
			//free(FSSRDP);
			//free(tMaskseq);
			//free(Fact3X3);
			//free(Fact);
			//free(ProbEstimate);

		}
		omp_set_num_threads(2);
		return(redonum);
	}

	int MyMathFuncs::AESecondRDP3(int oSeq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DP, unsigned char *DP2, int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseq, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact) {


		int redonum;
		int procs;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int UBXOHN, UBXSNW, FRC;


			dsb = (UBCS + 1)*(NextNo + 1);
			dsa = (NextNo + 1)*(NextNo + 1);
			dsc = 4 * (UBFSSRDP + 1)*(UBFSSRDP + 1)*(UBFSSRDP + 1);
			dsd = (UBFact3x3 + 1)*(UBFact3x3 + 1)*(UBFact3x3 + 1);
			dse = 172;
			dsf = 172 * 172 * 51;
			rlo = UBRL + 1;
			dpo = UBDP + 1;

			UBXOHN = LenStrainseq0 + XoverWindow * 2;
			UBXSNW = LenStrainseq0 + (int)(XOverWindowX / 2) * 2;
			short int *XoverSeqNum;
			char *XoverSeqNumW;
			int *XOverHomologyNum;

			XoverSeqNum = (short int*)calloc((LenStrainseq0 + 1) * 3, sizeof(short int));
			XoverSeqNumW = (char*)calloc((UBXSNW + 1) * 3, sizeof(char));
			XOverHomologyNum = (int*)calloc((UBXOHN + 1) * 3, sizeof(int));

			//unsigned char *DP;
			//unsigned char *CS;
			////float *PermDiffs, *PermValid, *Distance, *TreeDistance;
			//unsigned char *tMaskseq,*FSSRDP; 

			//double *Fact3X3, *Fact, *ProbEstimate;*/

			/*PermDiffs = (float*)calloc(dsa, sizeof(float));
			PermValid = (float*)calloc(dsa, sizeof(float));
			Distance = (float*)calloc(dsa, sizeof(float));
			TreeDistance = (float*)calloc(dsa, sizeof(float));*/
			//DP = (unsigned char*)calloc(dsa, sizeof(unsigned char));
			//CS = (unsigned char*)calloc(dsb, sizeof(unsigned char));
			//FSSRDP = (unsigned char*)calloc(dsc, sizeof(unsigned char));
			//tMaskseq = (unsigned char*)calloc(NextNo+1, sizeof(unsigned char));
			/*Fact3X3 = (double*)calloc(dsd, sizeof(double));
			Fact = (double*)calloc(dse, sizeof(double));
			ProbEstimate = (double*)calloc(dsf, sizeof(double));*/

			//#pragma omp critical
			//{
			//for (x = 0; x <= dsf; x++)
			//	ProbEstimate[x] = ProbEstimateO[x];
			//for (x = 0; x <= dsd; x++)
			//	Fact3X3[x] = Fact3X3O[x];
			//for (x = 0; x <= dse; x++)
			//	Fact[x] = FactO[x];
			/*for (x = 0; x <= NextNo + 1; x++)
			tMaskseq[x] = tMaskseqO[x];
			for (x = 0; x <= dsb; x++)
			CS[x] = CSO[x];
			for (x = 0; x <= dsc; x++)
			FSSRDP[x] = FSSRDPO[x];*/
			//for (x = 0; x <= dsa; x++)
			//	DP[x] = DPO[x];
			//for (x = 0; x <= dsa; x++)
			//	PermDiffs[x] = PermDiffsO[x];
			//for (x = 0; x <= dsa; x++)
			//	PermValid[x] = PermValidO[x];
			//for (x = 0; x <= dsa; x++)
			//	Distance[x] = DistanceO[x];
			//for (x = 0; x <= dsa; x++)
			//	TreeDistance[x] = TreeDistanceO[x];


			//}


			//AH = (int*)calloc(4, sizeof(int));

			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XoverSeqNum(Len(StrainSeq(0)), 2)
			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XOverHomologyNum(Len(StrainSeq(0)) + XoverWindow * 2, 2)

#pragma omp for
			for (Seq1 = 0; Seq1 <= oNextno; Seq1++) {
				if (Seq1 <= UBPV && tMaskseq[Seq1] == 0) {
					for (Seq2 = Seq1 + 1; Seq2 <= oNextno; Seq2++) {
						if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
							s12o = (int)(Seq1 + (UBPV + 1)*Seq2);
							pv12 = (int)(PermValid[s12o]);
							if (PermDiffs[s12o] > MinDIffs && pv12 > MinSeqSize) {
								//dp12 = (int)(DP[Seq1 + dpo*Seq2]);
								for (Seq3 = oNextno + 1; Seq3 <= NextNo; Seq3++) {
									s13o = (int)(Seq1 + (UBPV + 1)*Seq3);
									s23o = (int)(Seq2 + (UBPV + 1)*Seq3);
									pv13 = (int)(PermValid[s13o]);
									pv23 = (int)(PermValid[s23o]);
									dpos3 = Seq3*dpo;
									dp13 = (int)(DP[Seq1 + dpos3]);
									dp23 = (int)(DP[Seq2 + dpos3]);
									dp213 = (int)(DP2[Seq1 + dpos3]);
									dp223 = (int)(DP2[Seq2 + dpos3]);
									if (Seq3 <= UBPV && tMaskseq[Seq3] == 0 ) {
										if (pv13 > MinSeqSize && pv23 > MinSeqSize && ((dp23 + dp13) > 0 || (dp223 + dp213)>0)) {
											if (PermDiffs[s23o] > (float)(MinDIffs) && PermDiffs[s13o] > (float)(MinDIffs)) {
												NewOneFound = 0;
												BQPV = 0;

												//int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0,int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, int UBXSNW, unsigned char *CS,short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag,int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact
												FRC = FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
												if (FRC == 1) {
													//Call XOver(Distance(), XPosDiff(), XDiffPos(), CurrentXOver(), XoverList(), SeqNum(), Seq1, Seq2, Seq3, 0)
#pragma omp critical
													{
														redonum++;
														RL[redonum*rlo] = Seq1;
														RL[1 + redonum*rlo] = Seq2;
														RL[2 + redonum*rlo] = Seq3;
													}
													

												}

												




												
											}
										}
									}
								}
							}
						}
					}
				}
			}
			free(XoverSeqNum);
			free(XoverSeqNumW);
			free(XOverHomologyNum);

			//free(PermDiffs);
			//free(PermValid);
			//free(Distance);
			//free(TreeDistance);
			//free(DP);
			//free(CS);
			//free(FSSRDP);
			//free(tMaskseq);
			//free(Fact3X3);
			//free(Fact);
			//free(ProbEstimate);

		}
		omp_set_num_threads(2);
		return(redonum);
	}

	int MyMathFuncs::PrimaryRDP3(int oSeq1, int MinDIffs, int MinSeqSize, int oNextno, int NextNo, double SubThresh, int UBRL, int *RL, int UBDP, unsigned char *DP, unsigned char *DP2, int UBPV, float *PermDiffs, float *PermValid, unsigned char *tMaskseq, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact) {


		int redonum;
		int procs;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int UBXOHN, UBXSNW, FRC;


			dsb = (UBCS + 1)*(NextNo + 1);
			dsa = (NextNo + 1)*(NextNo + 1);
			dsc = 4 * (UBFSSRDP + 1)*(UBFSSRDP + 1)*(UBFSSRDP + 1);
			dsd = (UBFact3x3 + 1)*(UBFact3x3 + 1)*(UBFact3x3 + 1);
			dse = 172;
			dsf = 172 * 172 * 51;
			rlo = UBRL + 1;
			dpo = UBDP + 1;

			UBXOHN = LenStrainseq0 + XoverWindow * 2;
			UBXSNW = LenStrainseq0 + (int)(XOverWindowX / 2) * 2;
			short int *XoverSeqNum;
			char *XoverSeqNumW;
			int *XOverHomologyNum;

			XoverSeqNum = (short int*)calloc((LenStrainseq0 + 1) * 3, sizeof(short int));
			XoverSeqNumW = (char*)calloc((UBXSNW + 1) * 3, sizeof(char));
			XOverHomologyNum = (int*)calloc((UBXOHN + 1) * 3, sizeof(int));

			//unsigned char *DP;
			//unsigned char *CS;
			////float *PermDiffs, *PermValid, *Distance, *TreeDistance;
			//unsigned char *tMaskseq,*FSSRDP; 

			//double *Fact3X3, *Fact, *ProbEstimate;*/

			/*PermDiffs = (float*)calloc(dsa, sizeof(float));
			PermValid = (float*)calloc(dsa, sizeof(float));
			Distance = (float*)calloc(dsa, sizeof(float));
			TreeDistance = (float*)calloc(dsa, sizeof(float));*/
			//DP = (unsigned char*)calloc(dsa, sizeof(unsigned char));
			//CS = (unsigned char*)calloc(dsb, sizeof(unsigned char));
			//FSSRDP = (unsigned char*)calloc(dsc, sizeof(unsigned char));
			//tMaskseq = (unsigned char*)calloc(NextNo+1, sizeof(unsigned char));
			/*Fact3X3 = (double*)calloc(dsd, sizeof(double));
			Fact = (double*)calloc(dse, sizeof(double));
			ProbEstimate = (double*)calloc(dsf, sizeof(double));*/

			//#pragma omp critical
			//{
			//for (x = 0; x <= dsf; x++)
			//	ProbEstimate[x] = ProbEstimateO[x];
			//for (x = 0; x <= dsd; x++)
			//	Fact3X3[x] = Fact3X3O[x];
			//for (x = 0; x <= dse; x++)
			//	Fact[x] = FactO[x];
			/*for (x = 0; x <= NextNo + 1; x++)
			tMaskseq[x] = tMaskseqO[x];
			for (x = 0; x <= dsb; x++)
			CS[x] = CSO[x];
			for (x = 0; x <= dsc; x++)
			FSSRDP[x] = FSSRDPO[x];*/
			//for (x = 0; x <= dsa; x++)
			//	DP[x] = DPO[x];
			//for (x = 0; x <= dsa; x++)
			//	PermDiffs[x] = PermDiffsO[x];
			//for (x = 0; x <= dsa; x++)
			//	PermValid[x] = PermValidO[x];
			//for (x = 0; x <= dsa; x++)
			//	Distance[x] = DistanceO[x];
			//for (x = 0; x <= dsa; x++)
			//	TreeDistance[x] = TreeDistanceO[x];


			//}


			//AH = (int*)calloc(4, sizeof(int));

			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XoverSeqNum(Len(StrainSeq(0)), 2)
			//ReDim XoverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
			//ReDim XOverHomologyNum(Len(StrainSeq(0)) + XoverWindow * 2, 2)

#pragma omp for
			for (Seq1 = oSeq1; Seq1 <= oNextno; Seq1++) {
				if (Seq1 <= UBPV && tMaskseq[Seq1] == 0) {
					for (Seq2 = Seq1 + 1; Seq2 <= NextNo-1; Seq2++) {
						if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
							s12o = (int)(Seq1 + (UBPV + 1)*Seq2);
							pv12 = (int)(PermValid[s12o]);
							if (PermDiffs[s12o] > MinDIffs && pv12 > MinSeqSize) {
								//dp12 = (int)(DP[Seq1 + dpo*Seq2]);
								for (Seq3 = Seq2 + 1; Seq3 <= NextNo; Seq3++) {
									s13o = (int)(Seq1 + (UBPV + 1)*Seq3);
									s23o = (int)(Seq2 + (UBPV + 1)*Seq3);
									pv13 = (int)(PermValid[s13o]);
									pv23 = (int)(PermValid[s23o]);
									
									if (Seq3 <= UBPV && tMaskseq[Seq3] == 0) {
										if (pv13 > MinSeqSize && pv23 > MinSeqSize){// && ((dp23 + dp13) > 0 || (dp223 + dp213)>0)) {
											if (PermDiffs[s23o] > (float)(MinDIffs) && PermDiffs[s13o] > (float)(MinDIffs)) {
												NewOneFound = 0;
												BQPV = 0;

												//int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0,int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, int UBXSNW, unsigned char *CS,short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag,int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact
												FRC = FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
												if (FRC == 1) {
													//Call XOver(Distance(), XPosDiff(), XDiffPos(), CurrentXOver(), XoverList(), SeqNum(), Seq1, Seq2, Seq3, 0)
#pragma omp critical
													{
														redonum++;
														RL[redonum*rlo] = Seq1;
														RL[1 + redonum*rlo] = Seq2;
														RL[2 + redonum*rlo] = Seq3;
													}


												}







											}
										}
									}
								}
							}
						}
					}
				}
			}
			free(XoverSeqNum);
			free(XoverSeqNumW);
			free(XOverHomologyNum);

			//free(PermDiffs);
			//free(PermValid);
			//free(Distance);
			//free(TreeDistance);
			//free(DP);
			//free(CS);
			//free(FSSRDP);
			//free(tMaskseq);
			//free(Fact3X3);
			//free(Fact);
			//free(ProbEstimate);

		}
		omp_set_num_threads(2);
		return(redonum);
	}

	int MyMathFuncs::AlistRDP3(short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact) {


		int redonum;
		int procs;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = procs/2;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int UBXOHN, UBXSNW, FRC;


			dsb = (UBCS + 1)*(NextNo + 1);
			dsa = (NextNo + 1)*(NextNo + 1);
			dsc = 4 * (UBFSSRDP + 1)*(UBFSSRDP + 1)*(UBFSSRDP + 1);
			dsd = (UBFact3x3 + 1)*(UBFact3x3 + 1)*(UBFact3x3 + 1);
			dse = 172;
			dsf = 172 * 172 * 51;
			

			UBXOHN = LenStrainseq0 + XoverWindow * 2;
			UBXSNW = LenStrainseq0 + (int)(XOverWindowX / 2) * 2;
			short int *XoverSeqNum;
			char *XoverSeqNumW;
			int *XOverHomologyNum;

			XoverSeqNum = (short int*)calloc((LenStrainseq0 + 1) * 3, sizeof(short int));
			XoverSeqNumW = (char*)calloc((UBXSNW + 1) * 3, sizeof(char));
			XOverHomologyNum = (int*)calloc((UBXOHN + 1) * 3, sizeof(int));


#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				Seq1 = AList[y * 3];
				Seq2 = AList[1 + y * 3];
				Seq3 = AList[2 + y * 3];

				NewOneFound = 0;
				BQPV = 0;

				FRC = FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, SubThresh, LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
				if (FRC == 1)
					RL[y] = 1;
				
							
			}
			free(XoverSeqNum);
			free(XoverSeqNumW);
			free(XOverHomologyNum);

			//free(PermDiffs);
			//free(PermValid);
			//free(Distance);
			//free(TreeDistance);
			//free(DP);
			//free(CS);
			//free(FSSRDP);
			//free(tMaskseq);
			//free(Fact3X3);
			//free(Fact);
			//free(ProbEstimate);

		}
		omp_set_num_threads(2);
		return(redonum);
	}

	int MyMathFuncs::AlistRDP4(int ubslpv, double *StoreLPV, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact) {


		int redonum, slpvo;
		int procs;
		slpvo = ubslpv + 1;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = procs / 2;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int UBXOHN, UBXSNW, FRC;


			dsb = (UBCS + 1)*(NextNo + 1);
			dsa = (NextNo + 1)*(NextNo + 1);
			dsc = 4 * (UBFSSRDP + 1)*(UBFSSRDP + 1)*(UBFSSRDP + 1);
			dsd = (UBFact3x3 + 1)*(UBFact3x3 + 1)*(UBFact3x3 + 1);
			dse = 172;
			dsf = 172 * 172 * 51;


			UBXOHN = LenStrainseq0 + XoverWindow * 2;
			UBXSNW = LenStrainseq0 + (int)(XOverWindowX / 2) * 2;
			short int *XoverSeqNum;
			char *XoverSeqNumW;
			int *XOverHomologyNum;

			XoverSeqNum = (short int*)calloc((LenStrainseq0 + 1) * 3, sizeof(short int));
			XoverSeqNumW = (char*)calloc((UBXSNW + 1) * 3, sizeof(char));
			XOverHomologyNum = (int*)calloc((UBXOHN + 1) * 3, sizeof(int));


#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				Seq1 = AList[y * 3];
				Seq2 = AList[1 + y * 3];
				Seq3 = AList[2 + y * 3];

				NewOneFound = 0;
				BQPV = 1;

				FRC = FastRecCheckPB(CircularFlag, 0, MCCorrection, MCFlag, 1, SubThresh, LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
				if (FRC == 1) {
					//BQPV = BQPV*(double)(MCCorrection);
					if (BQPV < StoreLPV[Seq1*slpvo] || BQPV < StoreLPV[Seq2*slpvo]|| BQPV < StoreLPV[Seq3*slpvo])
						RL[y] = 1;
					else
						RL[y] = 2;
				}

			}
			free(XoverSeqNum);
			free(XoverSeqNumW);
			free(XOverHomologyNum);

			//free(PermDiffs);
			//free(PermValid);
			//free(Distance);
			//free(TreeDistance);
			//free(DP);
			//free(CS);
			//free(FSSRDP);
			//free(tMaskseq);
			//free(Fact3X3);
			//free(Fact);
			//free(ProbEstimate);

		}
		omp_set_num_threads(2);
		return(redonum);
	}
	int MyMathFuncs::AlistGC(char GCIndelFlag, short int GCMissmatchPen, int GCDimSize,short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag,int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC) {
		//     
			
			//ReDim SubSeq(Len(StrainSeq(0)), 6)
			
		int redonum;
		int procs;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = (procs / 2)-1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int  FRC;
			int *HighEnough;
			HighEnough = (int*)calloc(10, sizeof(int));
			int UBSS, UBMSP, UBFMS;
			char *SubSeq;
			int *MaxScorePos, *FragMaxScore;
			UBSS = LenStrainseq0;
			UBMSP = GCDimSize;
			UBFMS = GCDimSize;
			SubSeq = (char*)calloc((LenStrainseq0 + 1) * 7, sizeof(char));
			//ReDim FragMaxScore(GCDimSize, 5)
			//ReDim MaxScorePos(GCDimSize, 5)
			MaxScorePos = (int*)calloc((GCDimSize + 1) * 6, sizeof(int));
			FragMaxScore = (int*)calloc((GCDimSize + 1) * 6, sizeof(int));

			//ReDim PVals(GCDimSize, 5)
			double *PVals;
			int UBPV;
			UBPV = GCDimSize;
			PVals = (double*)calloc((GCDimSize + 1) * 6, sizeof(double));
			
			int *FragST, *FragEN, *FragScore, UBFS, UBFST;
			UBFST = GCDimSize;
			UBFS = GCDimSize;
			FragST = (int*)calloc((GCDimSize + 1) * 7, sizeof(int));
			FragEN = (int*)calloc((GCDimSize + 1) * 7, sizeof(int));
			FragScore = (int*)calloc((GCDimSize + 1) * 7, sizeof(int));

			//ReDim DeleteArray(Len(StrainSeq(0)) + 1)


#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				Seq1 = AList[y * 3];
				Seq2 = AList[1 + y * 3];
				Seq3 = AList[2 + y * 3];

				NewOneFound = 0;
				BQPV = 0;

				FRC = GCXoverDP2(&BQPV, UBCS, CS, UBFSSGC, FSSGC, MCFlag, UBPV, PVals, LowestProb,  MCCorrection,  ShortOutFlag,  CircularFlag,  GCDimSize, LenStrainseq0,  GCMissmatchPen, GCIndelFlag,  Seq1,  Seq2,  Seq3,  UBFST, FragST, FragEN,  UBFS, FragScore,  UBSS, SubSeq, UBMSP, MaxScorePos, UBFMS, FragMaxScore, HighEnough);
					//FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
				if (FRC == 1)
					RL[y] = 1;


			}
			free(HighEnough);
			free(SubSeq);
			free(MaxScorePos);
			free(FragMaxScore);
			free(PVals);
			free(FragST);
			free(FragEN);
			free(FragScore);

		}


		omp_set_num_threads(2);
		return(redonum);
	}

	//int MyMathFuncs::AlistGC                             (char GCIndelFlag, short int GCMissmatchPen, int GCDimSize, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC) {

	int MyMathFuncs::AlistGC2(int ubslpv, double *StoreLPV, char GCIndelFlag, short int GCMissmatchPen, int GCDimSize, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC) {
		//     

		//ReDim SubSeq(Len(StrainSeq(0)), 6)

		int redonum, slpvo;
		int procs;
		redonum = 0;
		slpvo = ubslpv + 1;
		procs = omp_get_num_procs();
		procs = (procs / 2) - 1;
		
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int  FRC;
			int *HighEnough;
			HighEnough = (int*)calloc(10, sizeof(int));
			int UBSS, UBMSP, UBFMS;
			char *SubSeq;
			int *MaxScorePos, *FragMaxScore;
			UBSS = LenStrainseq0;
			UBMSP = GCDimSize;
			UBFMS = GCDimSize;
			SubSeq = (char*)calloc((LenStrainseq0 + 1) * 7, sizeof(char));
			//ReDim FragMaxScore(GCDimSize, 5)
			//ReDim MaxScorePos(GCDimSize, 5)
			MaxScorePos = (int*)calloc((GCDimSize + 1) * 6, sizeof(int));
			FragMaxScore = (int*)calloc((GCDimSize + 1) * 6, sizeof(int));

			//ReDim PVals(GCDimSize, 5)
			double *PVals;
			int UBPV;
			UBPV = GCDimSize;
			PVals = (double*)calloc((GCDimSize + 1) * 6, sizeof(double));

			int *FragST, *FragEN, *FragScore, UBFS, UBFST;
			UBFST = GCDimSize;
			UBFS = GCDimSize;
			FragST = (int*)calloc((GCDimSize + 1) * 7, sizeof(int));
			FragEN = (int*)calloc((GCDimSize + 1) * 7, sizeof(int));
			FragScore = (int*)calloc((GCDimSize + 1) * 7, sizeof(int));

			//ReDim DeleteArray(Len(StrainSeq(0)) + 1)


#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				Seq1 = AList[y * 3];
				Seq2 = AList[1 + y * 3];
				Seq3 = AList[2 + y * 3];
				/*if (y==635)
					BQPV = 0;*/
					
				NewOneFound = 0;
				BQPV = 0;

				FRC = GCXoverDP2(&BQPV, UBCS, CS, UBFSSGC, FSSGC, MCFlag, UBPV, PVals, LowestProb, MCCorrection, ShortOutFlag, CircularFlag, GCDimSize, LenStrainseq0, GCMissmatchPen, GCIndelFlag, Seq1, Seq2, Seq3, UBFST, FragST, FragEN, UBFS, FragScore, UBSS, SubSeq, UBMSP, MaxScorePos, UBFMS, FragMaxScore, HighEnough);
				//FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
				if (FRC == 1) {
					BQPV = BQPV*(double)(MCCorrection);
					if (BQPV < StoreLPV[1 + Seq1*slpvo] || BQPV < StoreLPV[1 + Seq2*slpvo] || BQPV < StoreLPV[1 + Seq3*slpvo])
						RL[y] = 1;
					else
						RL[y] = 2;

				}
			}
			free(HighEnough);
			free(SubSeq);
			free(MaxScorePos);
			free(FragMaxScore);
			free(PVals);
			free(FragST);
			free(FragEN);
			free(FragScore);

		}


		omp_set_num_threads(2);
		return(redonum);
	}
	
	int MyMathFuncs::MakeAListIS(int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs) {

		int A, b, x, ALC, Seq1, Seq2, Seq3, Spos, alo, anlo, WinPPY, malc;
		ALC = -1;

		int rp;

		malc = TripListLen + RNum[WinPP] * 3;
		rp = rs[0];
		rs[0] = -1;
		for (x = rp; x <= TripListLen; x++) {
			if (Worthwhilescan[x] > 0) {
				anlo = (UBAnL + 1)*x;
				Seq1 = Analysislist[anlo];
				Seq2 = Analysislist[anlo + 1];
				Seq3 = Analysislist[anlo + 2];

				for (WinPPY = 0; WinPPY <= RNum[WinPP]; WinPPY++) {
					if (ALC >= malc)
						return (ALC);
					A = RList[WinPP + WinPPY*(UBRL + 1)];
					if (A > PermNextno)
						b = TraceSub[A];
					else
						b = A;

					if (Seq1 == b || Seq2 == b || Seq3 == b) {
						if (b == Seq1)
							Seq1 = A;
						else if (b == Seq2)
							Seq2 = A;
						else if (b == Seq3)
							Seq3 = A;

						if (ActualSeqSize[Seq1] > MinSeqSize) {
							if (ActualSeqSize[Seq2] > MinSeqSize) {
								if (ActualSeqSize[Seq3] > MinSeqSize) {
									if (DoPairs[Seq1 + Seq2*(UBDP + 1)] == 1 && DoPairs[Seq1 + Seq3*(UBDP + 1)] == 1 && DoPairs[Seq2 + Seq3*(UBDP + 1)] == 1) {

										if (ProgBinRead[Worthwhilescan[x] * (UBPB + 1)] == 1) {

											ALC++;


											alo = (int)(ALC*(UBAL + 1));
											AList[alo] = Seq1;
											AList[1 + alo] = Seq2;
											AList[2 + alo] = Seq3;
										}

									}

								}
							}
						}
					}

				}
			}

		}
		return(ALC);
	}
	
	
	int MyMathFuncs::AlistMC2(unsigned char *wws, int StartP, int EndP, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth, int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int NextNo, int ubslpv, double *StoreLPV, short int *AList, int AListLen, unsigned char *RL, short int CircularFlag, int MCCorrection, short int MCFlag, double UCTHresh, double LowestProb, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainseq0, int UBCS, unsigned char *CS, int UBFSSMC, unsigned char *FSSMC, short int *SeqNum, int *Chimap, float *ChiTable2) {

		//     

		//ReDim SubSeq(Len(StrainSeq(0)), 6)

		int redonum, slpvo;
		int procs;
		redonum = 0;
		slpvo = ubslpv + 1;
		procs = omp_get_num_procs();
		procs = (procs / 2) - 1;

		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int  FRC;





			int UBWS;
			UBWS = LenStrainseq0 + HWindowWidth * 2;
			int *BanWin, *Winscores, *XDiffPos;
			BanWin = (int*)calloc(LenStrainseq0*2 + HWindowWidth * 2 + 1, sizeof(int));
			Winscores = (int*)calloc((LenStrainseq0 + HWindowWidth * 2 + 1) * 3, sizeof(int));
			XDiffPos = (int*)calloc(LenStrainseq0 + 201, sizeof(int));
			

			unsigned char *Scores, *MDMap;
			MDMap = (unsigned char*)calloc(LenStrainseq0 + 1, sizeof(unsigned char));
			Scores = (unsigned char*)calloc((LenStrainseq0 + 1) * 3, sizeof(unsigned char));

			double *Chivals, *SmoothChi, *mtP;
			Chivals = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			SmoothChi = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			mtP = (double*)calloc(101, sizeof(double));





#pragma omp for
			for (y = StartP; y <= EndP; y++) {

				/*if (Seq1 == 15 && Seq2 == 44 && Seq3 == 119) {
				BQPV = 2;
				Seq1 = 15;
				Seq2 = 44;
				Seq3 = 119;
				}*/
				//FRC = GCXoverDP2(&BQPV, UBCS, CS, UBFSSGC, FSSGC, MCFlag, UBPV, PVals, LowestProb, MCCorrection, ShortOutFlag, CircularFlag, GCDimSize, LenStrainseq0, GCMissmatchPen, GCIndelFlag, Seq1, Seq2, Seq3, UBFST, FragST, FragEN, UBFS, FragScore, UBSS, SubSeq, UBMSP, MaxScorePos, UBFMS, FragMaxScore, HighEnough);
				FRC = 0;
				//while (FRC == 0) {
					Seq1 = AList[y * 3];
					Seq2 = AList[1 + y * 3];
					Seq3 = AList[2 + y * 3];
					/*if (y == 635)
					BQPV = 0;*/

					NewOneFound = 0;
					BQPV = 1;
					//note that findallflag is forced to 1 here
					FRC = FastRecCheckMC(0, LongWindedFlag, &BQPV, 1, UCTHresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, MCWinFract, MCWinSize, MCProportionFlag, LenStrainseq0, CriticalDiff, 1, UBFSSMC, UBCS, Seq1, Seq2, Seq3, CS, FSSMC, SeqNum, UBWS, Scores, Winscores, XDiffPos, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);
					//FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
					/*if (BQPV <= 0) {
					BQPV = 1;
					}*/

					//FRC = 0;

					if (FRC == 1) {
						BQPV = BQPV*(double)(MCCorrection);
						if (BQPV < StoreLPV[3 + Seq1*slpvo] || BQPV < StoreLPV[3 + Seq2*slpvo] || BQPV < StoreLPV[3 + Seq3*slpvo])
							RL[y] = 1;
						else
							RL[y] = 2;

					}
					/*else {
						if (wws[y] == 0) {
							if (AList[2 + y * 3] < NextNo && y < AListLen) {
								if (AList[2 + y * 3] + 1 < AList[2 + (y + 1) * 3])
									AList[2 + y * 3] = AList[2 + y * 3] + 1;
								else
									break;
							}
							else
								break;
						}
						else
							break;
					}*/
				//}
			}
			free(BanWin);
			free(Winscores);
			free(XDiffPos);
			free(MDMap);
			free(Scores);
			free(Chivals);
			free(SmoothChi);
			free(mtP);

		}


		omp_set_num_threads(2);
		return(redonum);
	}






	int MyMathFuncs::AlistMC3(int SEN, unsigned char *wws, int StartP, int EndP, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth,  int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int NextNo,int ubslpv, double *StoreLPV, short int *AList, int AListLen,  unsigned char *RL, short int CircularFlag, int MCCorrection, short int MCFlag, double UCTHresh, double LowestProb, double MCWinFract, int MCWinSize, short int MCProportionFlag,  int LenStrainseq0,  int UBCS, unsigned char *CS, int UBFSSMC, unsigned char *FSSMC, short int *SeqNum, unsigned char *MissingData, int *Chimap, float *ChiTable2) {
		                  
		//     

		//ReDim SubSeq(Len(StrainSeq(0)), 6)

		int redonum, slpvo, dumdum;
		int procs;
		redonum = 0;
		slpvo = ubslpv + 1;
		procs = omp_get_num_procs();
		procs = (procs / 2) - 1;
		dumdum = 0;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int  FRC;

			int UBWS;
			UBWS = LenStrainseq0 + HWindowWidth * 2;
			int *BanWin, *Winscores, *XDiffPos, *XPosDiff;
			BanWin = (int*)calloc(LenStrainseq0*2 + HWindowWidth * 2+ 1, sizeof(int));
			Winscores = (int*)calloc((LenStrainseq0 + HWindowWidth * 2 + 1) * 3, sizeof(int));
			XDiffPos = (int*)calloc(LenStrainseq0 + 201, sizeof(int));
			XPosDiff = (int*)calloc(LenStrainseq0 + 201, sizeof(int));
				
			unsigned char *Scores, *MDMap;
			MDMap = (unsigned char*)calloc(LenStrainseq0 + 1, sizeof(unsigned char));
			Scores = (unsigned char*)calloc((LenStrainseq0 + 1)*3, sizeof(unsigned char));
				
			double *Chivals,  *SmoothChi, *mtP;
			Chivals = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			SmoothChi = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			mtP = (double*)calloc(101, sizeof(double));

#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				
				/*if (Seq1 == 15 && Seq2 == 44 && Seq3 == 119) {
					BQPV = 2;
					Seq1 = 15;
					Seq2 = 44;
					Seq3 = 119;
				}*/
				//FRC = GCXoverDP2(&BQPV, UBCS, CS, UBFSSGC, FSSGC, MCFlag, UBPV, PVals, LowestProb, MCCorrection, ShortOutFlag, CircularFlag, GCDimSize, LenStrainseq0, GCMissmatchPen, GCIndelFlag, Seq1, Seq2, Seq3, UBFST, FragST, FragEN, UBFS, FragScore, UBSS, SubSeq, UBMSP, MaxScorePos, UBFMS, FragMaxScore, HighEnough);
				FRC = 0;
				//while (FRC == 0) {
					Seq1 = AList[y * 3];
					Seq2 = AList[1 + y * 3];
					Seq3 = AList[2 + y * 3];
					/*if (y == 635)
					BQPV = 0;*/

					NewOneFound = 0;
					/*if (Seq1 == 177 && Seq2 == 207 && Seq3 == 261 && SEN==88)
						dumdum = 1;*/

					BQPV = 1;
					FRC = FastRecCheckMC2(SEN, LongWindedFlag, &BQPV, 0, UCTHresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, MCWinFract, MCWinSize, MCProportionFlag, LenStrainseq0, CriticalDiff,FindallFlag, UBFSSMC, UBCS, Seq1, Seq2, Seq3, CS, FSSMC, SeqNum, MissingData, UBWS, Scores, Winscores, XDiffPos, XPosDiff, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);
					//FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
					/*if (BQPV <= 0) {
						BQPV = 1;
					}*/
					if (FRC == 1) {
						BQPV = BQPV*(double)(MCCorrection);
						if (BQPV < StoreLPV[3 + Seq1*slpvo] || BQPV < StoreLPV[3 + Seq2*slpvo] || BQPV < StoreLPV[3 + Seq3*slpvo])
							RL[y] = 1;
						else
							RL[y] = 2;

					}
					/*else {
						if (wws[y] == 0) {
							if (AList[2 + y * 3] < NextNo && y < AListLen) {
								if (AList[2 + y * 3] + 1 < AList[2 + (y + 1) * 3])
									AList[2 + y * 3] = AList[2 + y * 3] + 1;
								else
									break;
							}
							else
								break;
						}
						else
							break;
					}*/
				//}
			}
			free(BanWin);
			free(Winscores);
			free(XDiffPos);
			free(XPosDiff);
			free(MDMap);
			free(Scores);
			free(Chivals);
			free(SmoothChi);
			free(mtP);

		}


		omp_set_num_threads(2);
		return(redonum);
	}
	
	int MyMathFuncs::MakeAListISP2(int prg,int *rs, int UBPB, unsigned char *ProgBinRead,int *TraceSub,int WinPP, int *RNum, int UBRL, int *RList, int UBAnL , short int *Analysislist, int TripListLen,unsigned char *Worthwhilescan,int *ActualSeqSize,int PermNextno, int NextNo, int MinSeqSize, int UBAL, int UBAL2, short int *AList, int UBDP, unsigned char *DoPairs) {

		int A, b, x, ALC, Seq1, Seq2, Seq3, Spos, alo, anlo, WinPPY,malc;
		ALC = -1;
		
		int rp;

		malc = TripListLen + RNum[WinPP] * 3;
		rp = rs[0];
		rs[0] = -1;
		if (rp == -1)
			rp = 0;
        for (x = rp; x<= TripListLen; x++){
            if (Worthwhilescan[x] > 0){
				anlo = (UBAnL + 1)*x;
				Seq1 = Analysislist[anlo];
				Seq2 = Analysislist[anlo + 1];
                Seq3 = Analysislist[anlo + 2];
                
                for (WinPPY = 0; WinPPY <= RNum[WinPP]; WinPPY++){
					if (ALC >= malc)
						return (ALC);
					A = RList[WinPP + WinPPY*(UBRL + 1)];
					if (A > PermNextno)
						b = TraceSub[A];
					else
						b = A;
                    
                    if (Seq1 == b || Seq2 == b || Seq3 == b){
						if (b == Seq1)
							Seq1 = A;
						else if (b == Seq2)
							Seq2 = A;
						else if (b == Seq3)
							Seq3 = A;
                        
                        if (ActualSeqSize[Seq1] > MinSeqSize){
							if (ActualSeqSize[Seq2] > MinSeqSize) {
								if (ActualSeqSize[Seq3] > MinSeqSize) {
                                    if (DoPairs[Seq1 + Seq2*(UBDP+1)] == 1 && DoPairs[Seq1 + Seq3*(UBDP + 1)] == 1 && DoPairs[Seq2 + Seq3*(UBDP + 1)] == 1){
                                    
                                         if (ProgBinRead[prg+Worthwhilescan[x]*(UBPB+1)] == 1){
											 
                                            ALC++;
											if (ALC > UBAL2) {
												ALC--;
												return(ALC);
											}

											alo = (int)(ALC*(UBAL + 1));
											AList[alo] = Seq1;
											AList[1 + alo] = Seq2;
											AList[2 + alo] = Seq3;
                                         }
                                         
                                    }
                                
                               }
                            }
                        }
                    }
                    
                }
            }
            
        }
	return(ALC);
}

	int MyMathFuncs::MakeAListISP3(float probdo, int prg, int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize, int UBAL, int UBAL2, short int *AList, int UBDP, unsigned char *DoPairs) {

		int A, b, x, ALC, Seq1, Seq2, Seq3, Spos, alo, anlo, WinPPY, malc;
		float curprop;
		ALC = -1;

		int rp;

		malc = TripListLen + RNum[WinPP] * 3;
		rp = rs[0];
		rs[0] = -1;
		if (rp == -1)
			rp = 0;
		curprop = 1;
		for (x = rp; x <= TripListLen; x++) {
			if (Worthwhilescan[x] > 0) {
				anlo = (UBAnL + 1)*x;
				Seq1 = Analysislist[anlo];
				Seq2 = Analysislist[anlo + 1];
				Seq3 = Analysislist[anlo + 2];

				for (WinPPY = 0; WinPPY <= RNum[WinPP]; WinPPY++) {
					if (ALC >= malc)
						return (ALC);
					A = RList[WinPP + WinPPY*(UBRL + 1)];
					if (A > PermNextno)
						b = TraceSub[A];
					else
						b = A;

					if (Seq1 == b || Seq2 == b || Seq3 == b) {
						if (b == Seq1)
							Seq1 = A;
						else if (b == Seq2)
							Seq2 = A;
						else if (b == Seq3)
							Seq3 = A;

						if (ActualSeqSize[Seq1] > MinSeqSize) {
							if (ActualSeqSize[Seq2] > MinSeqSize) {
								if (ActualSeqSize[Seq3] > MinSeqSize) {
									if (DoPairs[Seq1 + Seq2*(UBDP + 1)] == 1 && DoPairs[Seq1 + Seq3*(UBDP + 1)] == 1 && DoPairs[Seq2 + Seq3*(UBDP + 1)] == 1) {

										if (ProgBinRead[prg + Worthwhilescan[x] * (UBPB + 1)] == 1) {
											curprop = curprop + probdo;
											if (curprop > 1) {
												curprop = curprop - 1;
												ALC++;
												if (ALC > UBAL2) {
													ALC--;
													return(ALC);
												}

												alo = (int)(ALC*(UBAL + 1));
												AList[alo] = Seq1;
												AList[1 + alo] = Seq2;
												AList[2 + alo] = Seq3;
											}
										}
									}

								}
							}
						}
					}

				}
			}

		}
		return(ALC);
	}


	int MyMathFuncs::MakeBALISP3(int prg, int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize,  int UBDP, unsigned char *DoPairs) {

		int A, b, x, ALC, Seq1, Seq2, Seq3, Spos, alo, anlo, WinPPY, malc;
		ALC = -1;

		int rp;

		malc = TripListLen + RNum[WinPP] * 3;
		rp = rs[0];
		rs[0] = -1;
		if (rp == -1)
			rp = 0;
		for (x = rp; x <= TripListLen; x++) {
			if (Worthwhilescan[x] > 0) {
				anlo = (UBAnL + 1)*x;
				Seq1 = Analysislist[anlo];
				Seq2 = Analysislist[anlo + 1];
				Seq3 = Analysislist[anlo + 2];

				for (WinPPY = 0; WinPPY <= RNum[WinPP]; WinPPY++) {
					if (ALC >= malc)
						return (ALC);
					A = RList[WinPP + WinPPY*(UBRL + 1)];
					if (A > PermNextno)
						b = TraceSub[A];
					else
						b = A;

					if (Seq1 == b || Seq2 == b || Seq3 == b) {
						if (b == Seq1)
							Seq1 = A;
						else if (b == Seq2)
							Seq2 = A;
						else if (b == Seq3)
							Seq3 = A;

						if (ActualSeqSize[Seq1] > MinSeqSize) {
							if (ActualSeqSize[Seq2] > MinSeqSize) {
								if (ActualSeqSize[Seq3] > MinSeqSize) {
									if (DoPairs[Seq1 + Seq2*(UBDP + 1)] == 1 && DoPairs[Seq1 + Seq3*(UBDP + 1)] == 1 && DoPairs[Seq2 + Seq3*(UBDP + 1)] == 1) {

										if (ProgBinRead[prg + Worthwhilescan[x] * (UBPB + 1)] == 1) {

											ALC++;
											
										}

									}

								}
							}
						}
					}

				}
			}

		}
		return(ALC);
	}



	int MyMathFuncs::MakeAListISP(int prg, int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs) {

		int A, b, x, ALC, Seq1, Seq2, Seq3, Spos, alo, anlo, WinPPY, malc;
		ALC = -1;

		int rp;

		malc = TripListLen + RNum[WinPP] * 3;
		rp = rs[0];
		rs[0] = -1;
		if (rp == -1)
			rp = 0;
		for (x = rp; x <= TripListLen; x++) {
			if (Worthwhilescan[x] > 0) {
				anlo = (UBAnL + 1)*x;
				Seq1 = Analysislist[anlo];
				Seq2 = Analysislist[anlo + 1];
				Seq3 = Analysislist[anlo + 2];

				for (WinPPY = 0; WinPPY <= RNum[WinPP]; WinPPY++) {
					if (ALC >= malc)
						return (ALC);
					A = RList[WinPP + WinPPY*(UBRL + 1)];
					if (A > PermNextno)
						b = TraceSub[A];
					else
						b = A;

					if (Seq1 == b || Seq2 == b || Seq3 == b) {
						if (b == Seq1)
							Seq1 = A;
						else if (b == Seq2)
							Seq2 = A;
						else if (b == Seq3)
							Seq3 = A;

						if (ActualSeqSize[Seq1] > MinSeqSize) {
							if (ActualSeqSize[Seq2] > MinSeqSize) {
								if (ActualSeqSize[Seq3] > MinSeqSize) {
									if (DoPairs[Seq1 + Seq2*(UBDP + 1)] == 1 && DoPairs[Seq1 + Seq3*(UBDP + 1)] == 1 && DoPairs[Seq2 + Seq3*(UBDP + 1)] == 1) {

										if (ProgBinRead[prg + Worthwhilescan[x] * (UBPB + 1)] == 1) {

											ALC++;

											alo = (int)(ALC*(UBAL + 1));
											AList[alo] = Seq1;
											AList[1 + alo] = Seq2;
											AList[2 + alo] = Seq3;
										}

									}

								}
							}
						}
					}

				}
			}

		}
		return(ALC);
	}
	int MyMathFuncs::AlistChi(int SEN, unsigned char *MissingData, unsigned char *wws, int StartP, int EndP, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth, int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int NextNo, int ubslpv, double *StoreLPV, short int *AList, int AListLen, unsigned char *RL, short int CircularFlag, int MCCorrection, short int MCFlag, double UCTHresh, double LowestProb, double CWinFract, int CWinSize, short int CProportionFlag, int LenStrainseq0, int UBCS, unsigned char *CS, int UBFSSRDP, unsigned char *FSSRDP, short int *SeqNum, int *Chimap, float *ChiTable2) {

		int redonum, slpvo;
		int procs;
		redonum = 0;
		slpvo = ubslpv + 1;
		procs = omp_get_num_procs();
		procs = (procs / 2) - 1;

		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int dx1, x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int  FRC, tSeq1, tSeq2, tSeq3,tvx;

			int UBWS;
			UBWS = LenStrainseq0 + HWindowWidth * 2;
			int *BanWin, *Winscores, *LXOS, *XDP, *XPD;//, *XDiffPos

			//UBound(WinScoresX, 1), ScoresX(0), WinScoresX(0),, ChiValsX(0), BanWin(0), MDMap(0),  mtP(0), SmoothChiX(0)

			// UBound(XDP, 1), XDP(0, 0)
			//Dim LXOS() As Long
			//	ReDim LXOS(3)
			//	Dim XDP() As Long
			//	ReDim XDP(Len(StrainSeq(0)) + 200, 2)
			LXOS = (int*)calloc(3, sizeof(int));
			XDP = (int*)calloc((LenStrainseq0 + 201) * 3, sizeof(int));
			XPD = (int*)calloc((LenStrainseq0 + 201) * 3, sizeof(int));
			BanWin = (int*)calloc(LenStrainseq0*2 + HWindowWidth * 2 + 1, sizeof(int));
			Winscores = (int*)calloc((LenStrainseq0 + HWindowWidth * 2 + 1) * 3, sizeof(int));
			//XDiffPos = (int*)calloc(LenStrainseq0 + 201, sizeof(int));


			unsigned char *Scores, *MDMap;
			MDMap = (unsigned char*)calloc(LenStrainseq0 + 1, sizeof(unsigned char));
			Scores = (unsigned char*)calloc((LenStrainseq0 + 1) * 3, sizeof(unsigned char));

			double *Chivals, *SmoothChi, *mtP;
			Chivals = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			SmoothChi = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			mtP = (double*)calloc(101, sizeof(double));

#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				/*if (y == 4229)
					dx1 = 1;*/
				FRC = 0;
				//while (FRC == 0) {
					Seq1 = AList[y * 3];
					Seq2 = AList[1 + y * 3];
					Seq3 = AList[2 + y * 3];


					
					

					tSeq1 = Seq1;
					tSeq2 = Seq2;
					tSeq3 = Seq3;
					NewOneFound = 0;
					BQPV = 1;
					RL[y] = 0;
					//FRC = FastRecCheckMC(0, LongWindedFlag, &BQPV, 1, UCTHresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, MCWinFract, MCWinSize, MCProportionFlag, LenStrainseq0, CriticalDiff, 1, UBFSSMC, UBCS, Seq1, Seq2, Seq3, CS, FSSMC, SeqNum, UBWS, Scores, Winscores, XDiffPos, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);
					FRC = FastRecCheckChim(MissingData, XPD, LXOS, 0,SEN, LongWindedFlag, &BQPV, 1, UCTHresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, CWinFract,CWinSize, CProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSRDP, UBCS, tSeq1, tSeq2, tSeq3, CS, FSSRDP, SeqNum, UBWS, Scores, Winscores, LenStrainseq0 + 200, XDP, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);

					if (FRC == 1) {
						BQPV = BQPV*(double)(MCCorrection);
						if (BQPV < StoreLPV[3 + tSeq1*slpvo] || BQPV < StoreLPV[3 + tSeq2*slpvo] || BQPV < StoreLPV[3 + tSeq3*slpvo])
							RL[y] = 1;
						else
							RL[y] = 2;

					}

					tSeq1 = Seq2;
					tSeq2 = Seq3;
					tSeq3 = Seq1;
					NewOneFound = 0;
					BQPV = 1;
					//FRC = FastRecCheckMC(0, LongWindedFlag, &BQPV, 1, UCTHresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, MCWinFract, MCWinSize, MCProportionFlag, LenStrainseq0, CriticalDiff, 1, UBFSSMC, UBCS, Seq1, Seq2, Seq3, CS, FSSMC, SeqNum, UBWS, Scores, Winscores, XDiffPos, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);
					FRC = FastRecCheckChim(MissingData, XPD, LXOS, 1, SEN, LongWindedFlag, &BQPV, 1, UCTHresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, CWinFract, CWinSize, CProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSRDP, UBCS, tSeq1, tSeq2, tSeq3, CS, FSSRDP, SeqNum, UBWS, Scores, Winscores, LenStrainseq0 + 200, XDP, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);

					if (FRC == 1) {
						BQPV = BQPV*(double)(MCCorrection);
						if (BQPV < StoreLPV[3 + tSeq1*slpvo] || BQPV < StoreLPV[3 + tSeq2*slpvo] || BQPV < StoreLPV[3 + tSeq3*slpvo])
							RL[y] = RL[y]+4;
						else
							RL[y] = RL[y]+8;

					}
					
					tSeq1 = Seq3;
					tSeq2 = Seq1;
					tSeq3 = Seq2;
					NewOneFound = 0;
					BQPV = 1;
					//FRC = FastRecCheckMC(0, LongWindedFlag, &BQPV, 1, UCTHresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, MCWinFract, MCWinSize, MCProportionFlag, LenStrainseq0, CriticalDiff, 1, UBFSSMC, UBCS, Seq1, Seq2, Seq3, CS, FSSMC, SeqNum, UBWS, Scores, Winscores, XDiffPos, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);
					FRC = FastRecCheckChim(MissingData, XPD, LXOS, 2, SEN, LongWindedFlag, &BQPV, 1, UCTHresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, CWinFract, CWinSize, CProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSRDP, UBCS, tSeq1, tSeq2, tSeq3, CS, FSSRDP, SeqNum, UBWS, Scores, Winscores, LenStrainseq0 + 200, XDP, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);
					/*if (tSeq1 == 5 && tSeq2 == 0 && tSeq3 == 1)
						tvx=tSeq2;*/
					if (FRC == 1) {
						BQPV = BQPV*(double)(MCCorrection);
						if (BQPV < StoreLPV[3 + tSeq1*slpvo] || BQPV < StoreLPV[3 + tSeq2*slpvo] || BQPV < StoreLPV[3 + tSeq3*slpvo])
							RL[y] = RL[y]+16;
						else
							RL[y] = RL[y]+32;

					}
					//if (SEN > 0)
					//	break;
					
					/*if (wws[y] == 0) {
						if (AList[2 + y * 3] < NextNo && y < AListLen) {
							if (AList[2 + y * 3] + 1 < AList[2 + (y + 1) * 3])
								AList[2 + y * 3] = AList[2 + y * 3] + 1;
							else
								break;
						}
						else
							break;
					}
					else
						break;*/

					
					Seq1 = tSeq1;
					Seq2 = tSeq2;
					Seq3 = tSeq3;


					
				//}
			}
			free(LXOS);
			free(XDP);
			free(XPD);
			free(BanWin);
			free(Winscores);
			//free(XDiffPos);
			free(MDMap);
			free(Scores);
			free(Chivals);
			free(SmoothChi);
			free(mtP);

		}


		omp_set_num_threads(2);
		return(redonum);
	}


	int MyMathFuncs::MakeLowCI(int Y, int TargetNum, int oPermNum, int PNA, int UBMS, float *MapS, int UBPVM, float *PValMap){
		int x, Z, ZDP, os, os2, WinPos;
		float TopS;
		double ZZ;
		os = UBMS + 1;
		os2 = UBPVM + 1;
		for (x = 0; x<= TargetNum; x++){
			TopS = 0.0;
			for (Z = 1; Z <= oPermNum; Z++){
				ZZ = (double)(Z) / (double)(PNA);
				ZDP = round(ZZ);
				if (MapS[Y + ZDP*os] >= TopS) {
					TopS = MapS[Y + ZDP*os];
					WinPos = ZDP;
				}
			}
			PValMap[Y + x*os2] = MapS[Y + WinPos*os];
			MapS[Y + WinPos*os] = -MapS[Y + WinPos*os] - 1;
		}
		return (1);
	}

	int MyMathFuncs::MakeHighCI(int PermNum, int Y, int TargetNum, int oPermNum, int PNA, int UBMS, float *MapS, int UBPVM, float *PValMap) {
		int x, Z, ZDP, os, os2, WinPos;
		float TopS;
		double ZZ;
		os = UBMS + 1;
		os2 = UBPVM + 1;
		for (x = PermNum; x >= (PermNum - TargetNum); x--) {
			TopS = 1000000.0;
			for (Z = 1; Z <= oPermNum; Z++) {
				ZZ = (double)(Z) / (double)(PNA);
				ZDP = round(ZZ);
				if (MapS[Y + ZDP*os] <= TopS) {
					if (MapS[Y + ZDP*os] >= 0) {
						TopS = MapS[Y + ZDP*os];
						WinPos = ZDP;
					}
				}
			}
			PValMap[Y + x*os2] = MapS[Y + WinPos*os];
			MapS[Y + WinPos*os] = -MapS[Y + WinPos*os] - 1;
		}
		return (1);
	}

	int MyMathFuncs::MakeAListOSP(int UBAL2,int BusyWithExcludes,int UBSV, float *SubValid, int sNextno, int UBTS1, int prg, int *rs, int UBPB, unsigned char *ProgBinRead, int *TraceSub, int WinPP, int *RNum, int UBRL, int *RList, int UBAnL, short int *Analysislist, int TripListLen, unsigned char *Worthwhilescan, int *ActualSeqSize, int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs) {

		int  A, b, x, ALC, Seq1, Seq2, Seq3, Spos, alo, anlo, WinPPY, malc, GoOn, tSeq1, tSeq2, tSeq3, svo;
		ALC = -1;

		int rp;
		svo = UBSV + 1;
		malc = TripListLen + RNum[WinPP] * 3;
		rp = rs[0];
		rs[0] = -1;
		if (rp < 0)
			rp = 0;
		for (x = rp; x <= TripListLen; x++) {
			if (ALC == UBAL2)
				return (ALC);
			if (Worthwhilescan[x] > 0) {
				anlo = (UBAnL + 1)*x;
				Seq1 = Analysislist[anlo];
				
				
				//CurrentTripListNum = x
				GoOn = 1;
				
				if (Seq1 <= UBTS1) {
					if (Seq1 > sNextno)
						tSeq1 = TraceSub[Seq1];
					else
						tSeq1 = Seq1;

				}
				else
					GoOn = 0;
				
				Seq2 = Analysislist[anlo + 1];
				if (Seq2 <= UBTS1) {
					if (Seq2 > sNextno)
						tSeq2 = TraceSub[Seq2];
					else
						tSeq2 = Seq2;

				}
				else
					GoOn = 0;

				if (DoPairs[tSeq1 + tSeq2*(UBDP + 1)] == 1) {
					Seq3 = Analysislist[anlo + 2];
					if (Seq3 <= UBTS1) {
						if (Seq3 > sNextno)
							tSeq3 = TraceSub[Seq3];
						else
							tSeq3 = Seq3;

					}
					else
						GoOn = 0;

					if (DoPairs[tSeq1 + tSeq3*(UBDP + 1)] == 1 && DoPairs[tSeq2 + tSeq3*(UBDP + 1)] == 1) {
						if (GoOn == 1) {
							for (WinPPY = NextNo - RNum[WinPP]; WinPPY <= NextNo; WinPPY++) {
								A = WinPPY;
								if (A > PermNextno)
									b = TraceSub[A];
								else
									b = A;

								if (Seq1 == b || Seq2 == b || Seq3 == b) {
									if (b == Seq1)
										Seq1 = A;
									else if (b == Seq2)
										Seq2 = A;
									else if (b == Seq3)
										Seq3 = A;

									if (ActualSeqSize[Seq1] > MinSeqSize) {
										if (ActualSeqSize[Seq2] > MinSeqSize) {
											if (ActualSeqSize[Seq3] > MinSeqSize) {
												if (SubValid[tSeq1 + tSeq3*svo] > 20) {
													if (SubValid[tSeq1 + tSeq2*svo] > 20) {
														if (ProgBinRead[prg + Worthwhilescan[x] * (UBPB + 1)] == 1 || BusyWithExcludes == 1) {

															
															if (ALC >= UBAL2) {
																ALC = UBAL2;
																return(ALC);
															}
															ALC++;
															alo = (int)(ALC*(UBAL + 1));
															AList[alo] = Seq1;
															AList[1 + alo] = Seq2;
															AList[2 + alo] = Seq3;
														}

													}
												}
											}
										}
									}
								}
							}
						}


					}
				}
				
			}

		}
		return(ALC);
	}
	
int MyMathFuncs::MakeAListASEF(int BAL,int *rs, int oNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DP, unsigned char *tMaskseq,int UBPV, float *PermValid, float *PermDiffs) {
	int ALC, Seq1, Seq2, Seq3,alo;
	int rp;
	rp = rs[0];
	rs[0] = -1;
	ALC = -1;

   for (Seq1 = rp; Seq1 <= oNextno; Seq1++){
        if (Seq1 <= UBPV && tMaskseq[Seq1] == 0){           
			for (Seq2 = oNextno + 1; Seq2 <= NextNo; Seq2++) {
                    if (Seq2 <= UBPV && tMaskseq[Seq2] == 0){
                        if (PermDiffs[Seq1 + Seq2*(UBPV+1)] > 2 && PermValid[Seq1 + Seq2*(UBPV + 1)] > MinSeqSize){
                            for (Seq3 = Seq2 + 1; Seq3 <= NextNo; Seq3++){
                                if (Seq3 <= UBPV && tMaskseq[Seq3] == 0 && ((int)(DP[Seq2 +Seq3*(UBDP+1)]) + (int)(DP[Seq1 + Seq3*(UBDP + 1)]) + (int)(DP[Seq1 + Seq2*(UBDP + 1)]) > 0))
									if (PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize) {
										if (PermDiffs[Seq2 + Seq3*(UBPV + 1)] > 2) {
											ALC++;
											if (ALC > BAL) {
												rs[0] = Seq1;
												return(ALC - 1);
											}

											alo = (int)(ALC*(UBAL + 1));
											AList[alo] = Seq1;
											AList[1 + alo] = Seq2;
											AList[2 + alo] = Seq3;
										}
									}
                                }
                            }
                        }
                    }
				}
				if (ALC > 1000000) {
					rs[0] = Seq1+1;
					break;
				}
					
            }
        

	 return(ALC);
}

int MyMathFuncs::MakeAListASEF2(float propdo,int BAL, int *rs, int oNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DP, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs) {
	int ALC, Seq1, Seq2, Seq3, alo;
	int rp;
	float curprop;
	rp = rs[0];
	rs[0] = -1;
	ALC = -1;
	curprop = 0;
	for (Seq1 = rp; Seq1 <= oNextno; Seq1++) {
		if (Seq1 <= UBPV && tMaskseq[Seq1] == 0) {
			for (Seq2 = oNextno + 1; Seq2 <= NextNo; Seq2++) {
				if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
					if (PermDiffs[Seq1 + Seq2*(UBPV + 1)] > 2 && PermValid[Seq1 + Seq2*(UBPV + 1)] > MinSeqSize) {
						for (Seq3 = Seq2 + 1; Seq3 <= NextNo; Seq3++) {
							if (Seq3 <= UBPV && tMaskseq[Seq3] == 0 && ((int)(DP[Seq2 + Seq3*(UBDP + 1)]) + (int)(DP[Seq1 + Seq3*(UBDP + 1)]) + (int)(DP[Seq1 + Seq2*(UBDP + 1)]) > 0))
								if (PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize) {
									if (PermDiffs[Seq2 + Seq3*(UBPV + 1)] > 2) {
										curprop = curprop + propdo;
										if (curprop > 1) {
											curprop = curprop - 1;
											ALC++;


											if (ALC > BAL) {
												rs[0] = Seq1;
												return(ALC - 1);
											}

											alo = (int)(ALC*(UBAL + 1));
											AList[alo] = Seq1;
											AList[1 + alo] = Seq2;
											AList[2 + alo] = Seq3;
										}
									}
								}
						}
					}
				}
			}
		}
		if (ALC > 1000000) {
			rs[0] = Seq1 + 1;
			break;
		}

	}


	return(ALC);
}

int MyMathFuncs::GetBALF(int oNextno, int NextNo, int MinSeqSize, int UBDP, unsigned char *DP, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs) {
	int Seq1, Seq2, Seq3, alo, BAL;
	int rp;
	rp = 0;
	BAL = 0;
	for (Seq1 = rp; Seq1 <= oNextno; Seq1++) {
		if (Seq1 <= UBPV && tMaskseq[Seq1] == 0) {
			for (Seq2 = oNextno + 1; Seq2 <= NextNo; Seq2++) {
				if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
					if (PermDiffs[Seq1 + Seq2*(UBPV + 1)] > 2 && PermValid[Seq1 + Seq2*(UBPV + 1)] > MinSeqSize) {
						for (Seq3 = Seq2 + 1; Seq3 <= NextNo; Seq3++) {
							if (Seq3 <= UBPV && tMaskseq[Seq3] == 0 && ((int)(DP[Seq2 + Seq3*(UBDP + 1)]) + (int)(DP[Seq1 + Seq3*(UBDP + 1)]) + (int)(DP[Seq1 + Seq2*(UBDP + 1)]) > 0))
								if (PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize) {
									if (PermDiffs[Seq2 + Seq3*(UBPV + 1)] > 2) {
										BAL++;
									}
								}
						}
					}
				}
			}
		}
		

	}
	return(BAL);
}

int MyMathFuncs::MakeSubSeqPermP(int LS, int D, int P1, int P2, short int *PermSeqNum, int *XDiffPos, int *XPosDiff) {
	int LSSeq, B, osd, osp1, osp2;
	LSSeq = 0;
	osd = D*(LS + 1);
	osp1 = P1*(LS + 1);
	osp2 = P2*(LS + 1);

	for (B = 0; B <= LS + 200; B++)
		XDiffPos[B] = 0;

	for (B = 0; B <= LS + 200; B++)
		XPosDiff[B] = 0;

	for (B = 1; B <= LS; B++) {

		XPosDiff[B] = LSSeq;
		if (PermSeqNum[B + osd] != 46) {
			if (PermSeqNum[B + osp1] != 46) {
				if (PermSeqNum[B + osp2] != 46) {
					if (PermSeqNum[B + osd] != PermSeqNum[B + osp1]) {
						LSSeq++;
						XPosDiff[B] = LSSeq;
						XDiffPos[LSSeq] = B;
					}
					else if (PermSeqNum[B + osd] != PermSeqNum[B + osp2]) {
						LSSeq++;
						XPosDiff[B] = LSSeq;
						XDiffPos[LSSeq] = B;


					}
				}
			}

		}
	}

	return(LSSeq);
}

//int MyMathFuncs::MakeAListASES(int *rs, int oNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DP, unsigned char *DP2, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs) {
//	int ALC, Seq1, Seq2, Seq3, alo;
//	ALC = -1;
//	int rp;
//	rp = rs[0];
//	rs[0] = -1;
//	for (Seq1 = rp; Seq1 <= oNextno; Seq1++){
//		if (Seq1 <= UBPV && tMaskseq[Seq1] == 0) {
//			for (Seq2 = Seq1 + 1; Seq2 <= oNextno; Seq2++) {
//				if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
//					if (PermDiffs[Seq1 + Seq2*(UBPV + 1)] > 2 && PermValid[Seq1 + Seq2*(UBPV + 1)] > MinSeqSize) {
//						for (Seq3 = oNextno + 1; Seq3 <= NextNo; Seq3++) {
//							if (tMaskseq[Seq3] == 0) {
//								if (PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize) {
//									if ((DP[Seq1 + Seq3*(UBDP + 1)] + DP[Seq2 + Seq3*(UBDP + 1)] > 0) || (DP2[Seq1 + Seq3*(UBDP + 1)] + DP2[Seq2 + Seq3*(UBDP + 1)]) == 2) {
//										ALC++;
//										alo = (int)(ALC*(UBAL + 1));
//										AList[alo] = Seq1;
//										AList[1 + alo] = Seq2;
//										AList[2 + alo] = Seq3;
//									}
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		if (ALC > 1000000) {
//			rs[0] = Seq1 + 1;
//			break;
//		}
//	}
//	return(ALC);
//}
int MyMathFuncs::MakeAListASES(int BAL,int *rs, int oNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DP, unsigned char *DP2, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs) {
	int ALC, Seq1, Seq2, Seq3, alo, maxadd;
	ALC = -1;
	int rp;
	rp = rs[0];
	rs[0] = -1;
	for (Seq1 = rp; Seq1 <= oNextno; Seq1++){
		if (Seq1 <= UBPV && tMaskseq[Seq1] == 0) {

			for (Seq2 = Seq1 + 1; Seq2 <= oNextno; Seq2++) {
				if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
					if (PermDiffs[Seq1 + Seq2*(UBPV + 1)] > 2 && PermValid[Seq1 + Seq2*(UBPV + 1)] > MinSeqSize) {
						for (Seq3 = oNextno + 1; Seq3 <= NextNo; Seq3++) {
							if (tMaskseq[Seq3] == 0) {
								if (PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize) {
									if ((DP[Seq1 + Seq3*(UBDP + 1)] + DP[Seq2 + Seq3*(UBDP + 1)] > 0) || (DP2[Seq1 + Seq3*(UBDP + 1)] + DP2[Seq2 + Seq3*(UBDP + 1)]) == 2) {
										ALC++;
										if (ALC > BAL) {
											rs[0] = Seq1;
											return(ALC - 1);
										}
										alo = (int)(ALC*(UBAL + 1));
										AList[alo] = Seq1;
										AList[1 + alo] = Seq2;
										AList[2 + alo] = Seq3;
									}
								}
							}
						}
					}
				}
			}
			//needed to ensure that ALC doesnt overshoot the bounds of alist in the next cycle
			//without it the program will end up rescanning some triplets
			if (Seq1 < oNextno-1) {
				maxadd = (oNextno - (Seq1 + 1))*(NextNo - oNextno);
				if (ALC + maxadd > BAL) {
					rs[0] = Seq1+1;
					return(ALC);
				}
			}
		}
		if (ALC > 1000000) {
			rs[0] = Seq1 + 1;
			break;
		}
	}
	return(ALC);
}

int MyMathFuncs::MakeAListASES2(float propdo, int BAL, int *rs, int oNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DP, unsigned char *DP2, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs) {
	int ALC, Seq1, Seq2, Seq3, alo, maxadd;
	float curprop;
	ALC = -1;
	int rp;
	rp = rs[0];
	rs[0] = -1;

	curprop = 1;
	for (Seq1 = rp; Seq1 <= oNextno; Seq1++) {
		if (Seq1 <= UBPV && tMaskseq[Seq1] == 0) {

			for (Seq2 = Seq1 + 1; Seq2 <= oNextno; Seq2++) {
				if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
					if (PermDiffs[Seq1 + Seq2*(UBPV + 1)] > 2 && PermValid[Seq1 + Seq2*(UBPV + 1)] > MinSeqSize) {
						for (Seq3 = oNextno + 1; Seq3 <= NextNo; Seq3++) {
							if (tMaskseq[Seq3] == 0) {
								if (PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize) {
									if ((DP[Seq1 + Seq3*(UBDP + 1)] + DP[Seq2 + Seq3*(UBDP + 1)] > 0) || (DP2[Seq1 + Seq3*(UBDP + 1)] + DP2[Seq2 + Seq3*(UBDP + 1)]) == 2) {
										curprop = curprop + propdo;
										if (curprop > 1) {
											curprop = curprop - 1;
											ALC++;
											if (ALC > BAL) {
												rs[0] = Seq1;
												return(ALC - 1);
											}
											alo = (int)(ALC*(UBAL + 1));
											AList[alo] = Seq1;
											AList[1 + alo] = Seq2;
											AList[2 + alo] = Seq3;
										}
									}
								}
							}
						}
					}
				}
			}
			//needed to ensure that ALC doesnt overshoot the bounds of alist in the next cycle
			//without it the program will end up rescanning some triplets
			if (Seq1 < oNextno - 1) {
				maxadd = (oNextno - (Seq1 + 1))*(NextNo - oNextno);
				if (ALC + maxadd > BAL) {
					rs[0] = Seq1 + 1;
					return(ALC);
				}
			}
		}
		if (ALC > 1000000) {
			rs[0] = Seq1 + 1;
			break;
		}
	}
	return(ALC);
}


int MyMathFuncs::GetBALS(int oNextno, int NextNo, int MinSeqSize, int UBDP, unsigned char *DP, unsigned char *DP2, unsigned char *tMaskseq, int UBPV, float *PermValid, float *PermDiffs) {
	int Seq1, Seq2, Seq3, maxadd, BAL;
	
	int rp;
	rp = 0;
	BAL = 0;
	for (Seq1 = rp; Seq1 <= oNextno; Seq1++) {
		if (Seq1 <= UBPV && tMaskseq[Seq1] == 0) {

			for (Seq2 = Seq1 + 1; Seq2 <= oNextno; Seq2++) {
				if (Seq2 <= UBPV && tMaskseq[Seq2] == 0) {
					if (PermDiffs[Seq1 + Seq2*(UBPV + 1)] > 2 && PermValid[Seq1 + Seq2*(UBPV + 1)] > MinSeqSize) {
						for (Seq3 = oNextno + 1; Seq3 <= NextNo; Seq3++) {
							if (tMaskseq[Seq3] == 0) {
								if (PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize) {
									if ((DP[Seq1 + Seq3*(UBDP + 1)] + DP[Seq2 + Seq3*(UBDP + 1)] > 0) || (DP2[Seq1 + Seq3*(UBDP + 1)] + DP2[Seq2 + Seq3*(UBDP + 1)]) == 2) {
										BAL++;
										
									}
								}
							}
						}
					}
				}
			}
			
		}
		
	}
	return(BAL);
}

int MyMathFuncs::MakeAListP(int PropTrips, int NextNo, short int *MaskSeq, int UBAL1, short int *Analysislist) {

	int x, Y, CP, Z, CurPos, os1;
	os1 = UBAL1 + 1;
	CurPos = 0;
	CP = -1;
	for (x = 0; x <= NextNo - 2; x++) {
		if (MaskSeq[x] == 0) {
			for (Y = x + 1; Y <= NextNo - 1; Y++) {
				if (MaskSeq[Y] == 0) {
					//CurPos = PropTrips - 1;
					for (Z = Y + 1; Z <= NextNo; Z++) {
						if (MaskSeq[Z] == 0) {
							CurPos = CurPos + 1;
							if (CurPos == PropTrips) {// || Z == NextNo){
								CurPos = 0;
								CP++;
								Analysislist[CP*os1] = x;
								Analysislist[1 + CP*os1] = Y;
								Analysislist[2 + CP*os1] = Z;
							}
						}
					}

				}
			}
		}
	}
	return(CP);
}
int MyMathFuncs::MakeAListP2(float PropTrips, int NextNo, short int *MaskSeq, int UBAL1, short int *Analysislist) {

	int x, Y, CP, Z,  os1;
	float CurPos, PT;
	os1 = UBAL1 + 1;
	CurPos = 1;
	CP = -1;
	for (x = 0; x <= NextNo - 2; x++) {
		if (MaskSeq[x] == 0) {
			for (Y = x + 1; Y <= NextNo - 1; Y++) {
				if (MaskSeq[Y] == 0) {
					//CurPos = PropTrips - 1;
					for (Z = Y + 1; Z <= NextNo; Z++) {
						if (MaskSeq[Z] == 0) {
							CurPos = CurPos + PropTrips;
							if (CurPos > 1) {// || Z == NextNo){
								CurPos = CurPos - 1;
								CP++;
								Analysislist[CP*os1] = x;
								Analysislist[1 + CP*os1] = Y;
								Analysislist[2 + CP*os1] = Z;
							
							}
							
						}
					}

				}
			}
		}
	}
	return(CP);
}

int MyMathFuncs::MakeAListISE(int *rs,  int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs, int UBPV, float *PermValid) {
	int ALC, Seq1, Seq2, Seq3, Spos, alo;
	ALC = -1;
	int rp;
	rp = rs[0];
	rs[0] = -1;
    for (Seq1 = rp; Seq1 <= PermNextno; Seq1++){
		for (Seq2 = Seq1 + 1; Seq2 <= NextNo; Seq2++) {
            if (Seq1 != Seq2){
				if (DoPairs[Seq1 + Seq2*(UBDP+1)] == 1){
				
					if (Seq2 > PermNextno)
						Spos = Seq2 + 1;
					else
						Spos = PermNextno + 1;
					
					for (Seq3 = Spos; Seq3 <= NextNo; Seq3++){
						if (DoPairs[Seq1 + Seq3*(UBDP+1)] == 1 && DoPairs[Seq2 + Seq3*(UBDP+1)] == 1){
							if (Seq1 <= UBPV && Seq2 <= UBPV && Seq3 <= UBPV){
								if (PermValid[Seq1 + Seq2*(UBPV+1)] > MinSeqSize && PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize){
									if (Seq1 <= NextNo && Seq2 <= NextNo && Seq3 <= NextNo) {
										ALC++;
										alo = (int)(ALC*(UBAL + 1));
										AList[alo] = Seq1;
										AList[1 + alo] = Seq2;
										AList[2 + alo] = Seq3;
									}
								}
							}
                            
						}
					}
				}
            }
        }
		if (ALC > 1000000) {
			rs[0] = Seq1+1;
			break;
		}
    }
	return(ALC);

	}
int MyMathFuncs::MakeAListISE2(int BAL, float propdo,int *rs, int PermNextno, int NextNo, int MinSeqSize, int UBAL, short int *AList, int UBDP, unsigned char *DoPairs, int UBPV, float *PermValid) {
	int ALC, Seq1, Seq2, Seq3, Spos, alo;
	float propcur;
	ALC = -1;
	int rp;
	rp = rs[0];
	rs[0] = -1;
	propcur = 1.0;
	for (Seq1 = rp; Seq1 <= PermNextno; Seq1++) {
		for (Seq2 = Seq1 + 1; Seq2 <= NextNo; Seq2++) {
			if (Seq1 != Seq2) {
				if (DoPairs[Seq1 + Seq2*(UBDP + 1)] == 1) {

					if (Seq2 > PermNextno)
						Spos = Seq2 + 1;
					else
						Spos = PermNextno + 1;

					for (Seq3 = Spos; Seq3 <= NextNo; Seq3++) {
						if (DoPairs[Seq1 + Seq3*(UBDP + 1)] == 1 && DoPairs[Seq2 + Seq3*(UBDP + 1)] == 1) {
							if (Seq1 <= UBPV && Seq2 <= UBPV && Seq3 <= UBPV) {
								if (PermValid[Seq1 + Seq2*(UBPV + 1)] > MinSeqSize && PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize) {
									if (Seq1 <= NextNo && Seq2 <= NextNo && Seq3 <= NextNo) {
										propcur = propcur + propdo;
										if (propcur > 1) {
											propcur = propcur - 1;
											ALC++;
											
											alo = (int)(ALC*(UBAL + 1));
											AList[alo] = Seq1;
											AList[1 + alo] = Seq2;
											AList[2 + alo] = Seq3;
											if (ALC == BAL)
												break;
										}
									}
								}
							}

						}
					}
					if (ALC == BAL)
						break;
				}
			}
		}
		if (ALC >= BAL/2) {
			rs[0] = Seq1 + 1;
			break;
		}
	}
	return(ALC);

}




int MyMathFuncs::MakeBALISE(int PermNextno, int NextNo, int MinSeqSize, int UBDP, unsigned char *DoPairs, int UBPV, float *PermValid) {
	int BAL, Seq1, Seq2, Seq3, Spos, alo;
	
	int rp;
	rp = 0;
	BAL = 0;
	for (Seq1 = rp; Seq1 <= PermNextno; Seq1++) {
		for (Seq2 = Seq1 + 1; Seq2 <= NextNo; Seq2++) {
			if (Seq1 != Seq2) {
				if (DoPairs[Seq1 + Seq2*(UBDP + 1)] == 1) {

					if (Seq2 > PermNextno)
						Spos = Seq2 + 1;
					else
						Spos = PermNextno + 1;

					for (Seq3 = Spos; Seq3 <= NextNo; Seq3++) {
						if (DoPairs[Seq1 + Seq3*(UBDP + 1)] == 1 && DoPairs[Seq2 + Seq3*(UBDP + 1)] == 1) {
							if (Seq1 <= UBPV && Seq2 <= UBPV && Seq3 <= UBPV) {
								if (PermValid[Seq1 + Seq2*(UBPV + 1)] > MinSeqSize && PermValid[Seq1 + Seq3*(UBPV + 1)] > MinSeqSize && PermValid[Seq2 + Seq3*(UBPV + 1)] > MinSeqSize) {
									if (Seq1 <= NextNo && Seq2 <= NextNo && Seq3 <= NextNo) {
										BAL++;
									}
								}
							}

						}
					}
				}
			}
		}
		
	}
	return(BAL);

}



//int MyMathFuncs::AEFirstAlistRDP3(int UBDP,unsigned char *DP2, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact) {

int MyMathFuncs::MakeCompressSeqP(int NextNo, int UBR, unsigned char *Recoded, int UBCS, unsigned char *CompressSeq, int UBCR1, int UBCR2, unsigned char *CompressorRDP){
		//Compression As Long, Y As Long, x As Long, Z As Long
		int Y, x, Z, ro, uc3, uc2;
		uc3 = (UBCR1 + 1)*(UBCR2 + 1);
		uc2 = (UBCR1 + 1);
		//ReDim CompressSeq((UBound(Recoded, 1) / Compression) + 1, NextNo)
		for (Y = 0; Y<= NextNo; Y++){
			Z = 0;
			for (x = 1; x <= UBR - 3; x=x+3){
				Z++;
				ro = x + Y*(UBR + 1);
				CompressSeq[Z + Y*(UBCS+1)] = CompressorRDP[Recoded[ro] + Recoded[ro+1]*uc2 + Recoded[ro+2] * uc3];
			}
		}
		return(1);
}

	int MyMathFuncs::AEFirstAlistRDP3(int UBDP,unsigned char *DP2, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, unsigned char *CS, short int *Seqnum, int XoverWindow, short int XOverWindowX, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact) {


		int redonum;
		int procs;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = procs/2;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int UBXOHN, UBXSNW, FRC;

			dpo = UBDP + 1;
			dsb = (UBCS + 1)*(NextNo + 1);
			dsa = (NextNo + 1)*(NextNo + 1);
			dsc = 4 * (UBFSSRDP + 1)*(UBFSSRDP + 1)*(UBFSSRDP + 1);
			dsd = (UBFact3x3 + 1)*(UBFact3x3 + 1)*(UBFact3x3 + 1);
			dse = 172;
			dsf = 172 * 172 * 51;


			UBXOHN = LenStrainseq0 + XoverWindow * 2;
			UBXSNW = LenStrainseq0 + (int)(XOverWindowX / 2) * 2;
			short int *XoverSeqNum;
			char *XoverSeqNumW;
			int *XOverHomologyNum;

			XoverSeqNum = (short int*)calloc((LenStrainseq0 + 1) * 3, sizeof(short int));
			XoverSeqNumW = (char*)calloc((UBXSNW + 1) * 3, sizeof(char));
			XOverHomologyNum = (int*)calloc((UBXOHN + 1) * 3, sizeof(int));


#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				Seq1 = AList[y * 3];
				Seq2 = AList[1 + y * 3];
				Seq3 = AList[2 + y * 3];

				NewOneFound = 0;
				BQPV = 0;

				FRC = FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
				if (FRC == 1) {
					RL[y] = 1;
					NewOneFound = 1;
				}



				if (BQPV > 0) {
					if (BQPV < SubThresh)
						NewOneFound = 1;
				}




				if (NewOneFound == 1) {
#pragma omp critical
					{
						if (DP2[Seq1 + dpo*Seq3] == 0) {
							DP2[Seq1 + dpo*Seq3] = 1;
							DP2[Seq3 + dpo*Seq1] = 1;
						}
						if (DP2[Seq2 + dpo*Seq1] == 0) {
							DP2[Seq2 + dpo*Seq1] = 1;
							DP2[Seq1 + dpo*Seq2] = 1;
						}
					}
					
				}

			}
			free(XoverSeqNum);
			free(XoverSeqNumW);
			free(XOverHomologyNum);

			//free(PermDiffs);
			//free(PermValid);
			//free(Distance);
			//free(TreeDistance);
			//free(DP);
			//free(CS);
			//free(FSSRDP);
			//free(tMaskseq);
			//free(Fact3X3);
			//free(Fact);
			//free(ProbEstimate);

		}
		omp_set_num_threads(2);
		return(redonum);
	}




	//NumRedos = AEFirstAlistMC(LongWindedFlag, ShortOutFlag, MaxABWin, HWindowWidth, lHWindowWidth, CriticalDiff, FindallFlag, UBound(DP2), DP2(0, 0), AList(0, 0), ALC, Y, EPX, NextNo, SubThresh, RedoL3(0), CircularFlag,                                                                                                           MCCorrection, MCFlag, LowestProb,                MCWinFract, MCWinSize, MCProportionFlag, Len(StrainSeq(0)),                       ShortOutFlag, UBound(CompressSeq, 1), CompressSeq(0, 0), UBound(FSSMC, 2), FSSMC(0, 0, 0, 0), SeqNum(0, 0), Chimap(0), ChiTable2(0))

	//int MyMathFuncs::AlistMC2(unsigned char *wws, int StartP, int EndP, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth,  int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int NextNo,int ubslpv, double *StoreLPV, short int *AList, int AListLen,  unsigned char *RL, short int CircularFlag, int MCCorrection, short int MCFlag, double UCTHresh, double LowestProb, double MCWinFract, int MCWinSize, short int MCProportionFlag,  int LenStrainseq0,  int UBCS, unsigned char *CS, int UBFSSMC, unsigned char *FSSMC, short int *SeqNum, int *Chimap, float *ChiTable2) {
	
	int MyMathFuncs::AEFirstAlistMC(int SEN, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth, int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int UBDP, unsigned char *DP2, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, double MCWinFract, int MCWinSize, short int MCProportionFlag, int LenStrainseq0,  int UBCS, unsigned char *CS, int UBFSSMC, unsigned char *FSSMC, short int *SeqNum, unsigned char *MissingData, int *Chimap, float *ChiTable2) {
		//int MyMathFuncs::AlistGC(char GCIndelFlag, short int GCMissmatchPen, int GCDimSize,short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag,int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC) {


		int redonum;
		int procs;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = procs / 2;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int bsss,x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int  FRC;
			bsss = 0;
			dpo = UBDP + 1;

			int UBWS;
			UBWS = LenStrainseq0 + HWindowWidth * 2;
			int *BanWin, *Winscores, *XDiffPos, *XPosDiff;
			BanWin = (int*)calloc(LenStrainseq0*2 + HWindowWidth * 2 + 1, sizeof(int));
			Winscores = (int*)calloc((LenStrainseq0 + HWindowWidth * 2 + 1) * 3, sizeof(int));
			XDiffPos = (int*)calloc(LenStrainseq0 + 201, sizeof(int));
			XPosDiff = (int*)calloc(LenStrainseq0 + 201, sizeof(int));
			unsigned char *Scores, *MDMap;
			MDMap = (unsigned char*)calloc(LenStrainseq0 + 1, sizeof(unsigned char));
			Scores = (unsigned char*)calloc((LenStrainseq0 + 1) * 3, sizeof(unsigned char));

			double *Chivals, *SmoothChi, *mtP;
			Chivals = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			SmoothChi = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			mtP = (double*)calloc(101, sizeof(double));




#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				Seq1 = AList[y * 3];
				Seq2 = AList[1 + y * 3];
				Seq3 = AList[2 + y * 3];
				//48,229,234
				//if (Seq1 == 48 && Seq2 == 229 && Seq3 == 234)
				//	bsss = 1;

				NewOneFound = 0;
				BQPV = 1;
				

				//FRC = FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
				//FRC = GCXoverDP2(&BQPV, UBCS, CS, UBFSSGC, FSSGC, MCFlag, UBPV, PVals, LowestProb, MCCorrection, ShortOutFlag, CircularFlag, GCDimSize, LenStrainseq0, GCMissmatchPen, GCIndelFlag, Seq1, Seq2, Seq3, UBFST, FragST, FragEN, UBFS, FragScore, UBSS, SubSeq, UBMSP, MaxScorePos, UBFMS, FragMaxScore, HighEnough);

				FRC = FastRecCheckMC2(SEN, LongWindedFlag, &BQPV, 0, SubThresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, MCWinFract, MCWinSize, MCProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSMC, UBCS, Seq1, Seq2, Seq3, CS, FSSMC, SeqNum, MissingData, UBWS, Scores, Winscores, XDiffPos, XPosDiff, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);


				if (FRC == 1) {
					//BQPV = BQPV*(double)(MCCorrection);
					//if (BQPV < StoreLPV[3 + Seq1*slpvo] || BQPV < StoreLPV[3 + Seq2*slpvo] || BQPV < StoreLPV[3 + Seq3*slpvo]) {
						RL[y] = 1;
						NewOneFound = 1;
					//}
				}


				//BQPV = BQPV*(double)(MCCorrection);
				if (BQPV > 0) {
					if (BQPV < SubThresh)
						NewOneFound = 1;
				}




				if (NewOneFound == 1) {
#pragma omp critical
					{
						if (DP2[Seq1 + dpo*Seq3] == 0) {
							DP2[Seq1 + dpo*Seq3] = 1;
							DP2[Seq3 + dpo*Seq1] = 1;
						}
						if (DP2[Seq2 + dpo*Seq1] == 0) {
							DP2[Seq2 + dpo*Seq1] = 1;
							DP2[Seq1 + dpo*Seq2] = 1;
						}
					}

				}

			}
			free(BanWin);
			free(Winscores);
			free(XDiffPos);
			free(XPosDiff);
			free(MDMap);
			free(Scores);
			free(Chivals);
			free(SmoothChi);
			free(mtP);



		}
		omp_set_num_threads(2);
		return(redonum);
	}

	
	int MyMathFuncs::AEFirstAlistChi(int SEN, int LongWindedFlag, short int ShortOutFlag, int MaxABWin, int HWindowWidth, int lHWindowWidth, int CriticalDiff, unsigned char FindallFlag, int UBDP, unsigned char *DP2, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, double CWinFract, int CWinSize, short int CProportionFlag, int LenStrainseq0, int UBCS, unsigned char *CS, int UBFSSRDP, unsigned char *FSSRDP, short int *SeqNum, unsigned char *MissingData, int *Chimap, float *ChiTable2) {
		//int MyMathFuncs::AlistGC(char GCIndelFlag, short int GCMissmatchPen, int GCDimSize,short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag,int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC) {


		int redonum;
		int procs;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = procs / 2;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int bsss, x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int  FRC, tSeq1, tSeq2, tSeq3;
			bsss = 0;
			dpo = UBDP + 1;

			int UBWS;
			UBWS = LenStrainseq0 + HWindowWidth * 2;

			int *BanWin, *Winscores, *XDiffPos, *XPosDiff, *LXOS, *XDP, *XPD;
			LXOS = (int*)calloc(3, sizeof(int));
			XDP = (int*)calloc((LenStrainseq0 + 201) * 3, sizeof(int));
			XPD = (int*)calloc((LenStrainseq0 + 201) * 3, sizeof(int));
			BanWin = (int*)calloc(LenStrainseq0*2 + HWindowWidth * 2 + 1, sizeof(int));
			Winscores = (int*)calloc((LenStrainseq0 + HWindowWidth * 2 + 1) * 3, sizeof(int));
			XDiffPos = (int*)calloc(LenStrainseq0 + 201, sizeof(int));
			XPosDiff = (int*)calloc(LenStrainseq0 + 201, sizeof(int));
			unsigned char *Scores, *MDMap;
			MDMap = (unsigned char*)calloc(LenStrainseq0 + 1, sizeof(unsigned char));
			Scores = (unsigned char*)calloc((LenStrainseq0 + 1) * 3, sizeof(unsigned char));

			double *Chivals, *SmoothChi, *mtP;
			Chivals = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			SmoothChi = (double*)calloc((LenStrainseq0 + 1) * 3, sizeof(double));
			mtP = (double*)calloc(101, sizeof(double));




#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				Seq1 = AList[y * 3];
				Seq2 = AList[1 + y * 3];
				Seq3 = AList[2 + y * 3];
				//48,229,234
				//if (Seq1 == 48 && Seq2 == 229 && Seq3 == 234)
				//	bsss = 1;

				tSeq1 = Seq1;
				tSeq2 = Seq2;
				tSeq3 = Seq3;
				NewOneFound = 0;
				BQPV = 1;




				//                          FRC = FastRecCheckMC2(SEN, LongWindedFlag, &BQPV, 0, SubThresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, MCWinFract, MCWinSize, MCProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSMC, UBCS, Seq1, Seq2, Seq3, CS, FSSMC, SeqNum, MissingData, UBWS, Scores, Winscores, XDiffPos, XPosDiff, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);
				FRC = FastRecCheckChim(MissingData, XPD, LXOS, 0, SEN, LongWindedFlag, &BQPV, 0, SubThresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, CWinFract, CWinSize, CProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSRDP, UBCS, tSeq1, tSeq2, tSeq3, CS, FSSRDP, SeqNum, UBWS, Scores, Winscores, LenStrainseq0 + 200, XDP, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);


				if (FRC == 1) {
					
					RL[y] = 1;
					NewOneFound = 1;
					
				}


				
				if (BQPV > 0) {
					if (BQPV < SubThresh)
						NewOneFound = 1;
				}

				tSeq1 = Seq2;
				tSeq2 = Seq3;
				tSeq3 = Seq1;
				
				BQPV = 1;

				//                          FRC = FastRecCheckMC2(SEN, LongWindedFlag, &BQPV, 0, SubThresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, MCWinFract, MCWinSize, MCProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSMC, UBCS, Seq1, Seq2, Seq3, CS, FSSMC, SeqNum, MissingData, UBWS, Scores, Winscores, XDiffPos, XPosDiff, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);
				FRC = FastRecCheckChim(MissingData, XPD, LXOS, 1, SEN, LongWindedFlag, &BQPV, 0, SubThresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, CWinFract, CWinSize, CProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSRDP, UBCS, tSeq1, tSeq2, tSeq3, CS, FSSRDP, SeqNum, UBWS, Scores, Winscores, LenStrainseq0 + 200, XDP, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);


				if (FRC == 1) {

					RL[y] = RL[y] + 4;
					NewOneFound = 1;

				}

				if (BQPV > 0) {
					if (BQPV < SubThresh)
						NewOneFound = 1;
				}

				tSeq1 = Seq3;
				tSeq2 = Seq1;
				tSeq3 = Seq2;
				
				BQPV = 1;

				//                          FRC = FastRecCheckMC2(SEN, LongWindedFlag, &BQPV, 0, SubThresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, MCWinFract, MCWinSize, MCProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSMC, UBCS, Seq1, Seq2, Seq3, CS, FSSMC, SeqNum, MissingData, UBWS, Scores, Winscores, XDiffPos, XPosDiff, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);
				FRC = FastRecCheckChim(MissingData, XPD, LXOS, 2, SEN, LongWindedFlag, &BQPV, 0, SubThresh, MCFlag, ShortOutFlag, MCCorrection, LowestProb, CircularFlag, NextNo, MaxABWin, HWindowWidth, lHWindowWidth, CWinFract, CWinSize, CProportionFlag, LenStrainseq0, CriticalDiff, FindallFlag, UBFSSRDP, UBCS, tSeq1, tSeq2, tSeq3, CS, FSSRDP, SeqNum, UBWS, Scores, Winscores, LenStrainseq0 + 200, XDP, Chivals, BanWin, MDMap, ChiTable2, Chimap, mtP, SmoothChi);


				if (FRC == 1) {

					RL[y] = RL[y] + 16;
					NewOneFound = 1;

				}

				if (BQPV > 0) {
					if (BQPV < SubThresh)
						NewOneFound = 1;
				}


				if (NewOneFound == 1) {
#pragma omp critical
					{
						if (DP2[Seq1 + dpo*Seq3] == 0) {
							DP2[Seq1 + dpo*Seq3] = 1;
							DP2[Seq3 + dpo*Seq1] = 1;
						}
						if (DP2[Seq2 + dpo*Seq1] == 0) {
							DP2[Seq2 + dpo*Seq1] = 1;
							DP2[Seq1 + dpo*Seq2] = 1;
						}
					}

				}

			}
			free(BanWin);
			free(Winscores);
			free(XDiffPos);
			free(XPosDiff);
			free(MDMap);
			free(Scores);
			free(Chivals);
			free(SmoothChi);
			free(mtP);



		}
		omp_set_num_threads(2);
		return(redonum);
	}

	int MyMathFuncs::AEFirstAlistGC(char GCIndelFlag, short int GCMissmatchPen, int GCDimSize, int UBDP, unsigned char *DP2, short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC) {
	//int MyMathFuncs::AlistGC(char GCIndelFlag, short int GCMissmatchPen, int GCDimSize,short int *AList, int AListLen, int StartP, int EndP, int NextNo, double SubThresh, unsigned char *RL, int CircularFlag, int MCCorrection, int MCFlag, double LowestProb, int TargetX, int LenStrainseq0, int ShortOutFlag,int UBFSSGC, int UBCS, unsigned char *CS, unsigned char *FSSGC) {
	

		int redonum;
		int procs;
		redonum = 0;
		procs = omp_get_num_procs();
		procs = procs / 2;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);

#pragma omp parallel //private (Seq1, MinDIffs, MinSeqSize, oNextno, NextNo, SubThresh, UBRL,  UBDP,   UBPV,  CircularFlag,  MCCorrection,  MCFlag,  LowestProb,  TargetX,  LenStrainseq0,  ShortOutFlag,  UBD,  UBTD,  UBFSSRDP, UBCS, XoverWindow, XOverWindowX,  ProbEstimateInFileFlag, UBPE1, UBPE2,  UBFact3x3)

		{
			int x, y, Seq1, Seq2, Seq3, s12o, s13o, s23o, dpo, dpos3, rlo, dp12, dp13, dp23, dp213, dp223, pv12, pv23, pv13, dsa, dsb, dsc, dsd, dse, dsf;
			unsigned char  NewOneFound;
			double BQPV;
			int UBXOHN, UBXSNW, FRC;

			dpo = UBDP + 1;
			dsb = (UBCS + 1)*(NextNo + 1);
			dsa = (NextNo + 1)*(NextNo + 1);
			dsc = 4 * (UBFSSGC + 1)*(UBFSSGC + 1)*(UBFSSGC + 1);
			
			dse = 172;
			dsf = 172 * 172 * 51;

			int *HighEnough;
			HighEnough = (int*)calloc(10, sizeof(int));
			int UBSS, UBMSP, UBFMS;
			char *SubSeq;
			int *MaxScorePos, *FragMaxScore;
			UBSS = LenStrainseq0;
			UBMSP = GCDimSize;
			UBFMS = GCDimSize;
			SubSeq = (char*)calloc((LenStrainseq0 + 1) * 7, sizeof(char));
			//ReDim FragMaxScore(GCDimSize, 5)
			//ReDim MaxScorePos(GCDimSize, 5)
			MaxScorePos = (int*)calloc((GCDimSize + 1) * 6, sizeof(int));
			FragMaxScore = (int*)calloc((GCDimSize + 1) * 6, sizeof(int));

			//ReDim PVals(GCDimSize, 5)
			double *PVals;
			int UBPV;
			UBPV = GCDimSize;
			PVals = (double*)calloc((GCDimSize + 1) * 6, sizeof(double));

			int *FragST, *FragEN, *FragScore, UBFS, UBFST;
			UBFST = GCDimSize;
			UBFS = GCDimSize;
			FragST = (int*)calloc((GCDimSize + 1) * 7, sizeof(int));
			FragEN = (int*)calloc((GCDimSize + 1) * 7, sizeof(int));
			FragScore = (int*)calloc((GCDimSize + 1) * 7, sizeof(int));
			

			


#pragma omp for
			for (y = StartP; y <= EndP; y++) {
				Seq1 = AList[y * 3];
				Seq2 = AList[1 + y * 3];
				Seq3 = AList[2 + y * 3];

				NewOneFound = 0;
				BQPV = 0;
				

				//FRC = FastRecCheckP(CircularFlag, 0, MCCorrection, MCFlag, 1, (LowestProb / MCCorrection), LowestProb, NextNo, TargetX, LenStrainseq0, ShortOutFlag, UBD, Distance, UBTD, TreeDistance, UBFSSRDP, UBCS, UBXSNW, CS, Seqnum, Seq1, Seq2, Seq3, LenStrainseq0 + 1, XoverWindow, XOverWindowX, XoverSeqNum, XoverSeqNumW, UBXOHN, XOverHomologyNum, FSSRDP, ProbEstimateInFileFlag, UBPE1, UBPE2, ProbEstimate, UBFact3x3, Fact3X3, Fact, &BQPV);
				FRC = GCXoverDP2(&BQPV, UBCS, CS, UBFSSGC, FSSGC, MCFlag, UBPV, PVals, LowestProb, MCCorrection, ShortOutFlag, CircularFlag, GCDimSize, LenStrainseq0, GCMissmatchPen, GCIndelFlag, Seq1, Seq2, Seq3, UBFST, FragST, FragEN, UBFS, FragScore, UBSS, SubSeq, UBMSP, MaxScorePos, UBFMS, FragMaxScore, HighEnough);

				if (FRC == 1) {
					RL[y] = 1;
					NewOneFound = 1;
				}



				if (BQPV > 0) {
					if (BQPV < SubThresh)
						NewOneFound = 1;
				}




				if (NewOneFound == 1) {
#pragma omp critical
					{
						if (DP2[Seq1 + dpo*Seq3] == 0) {
							DP2[Seq1 + dpo*Seq3] = 1;
							DP2[Seq3 + dpo*Seq1] = 1;
						}
						if (DP2[Seq2 + dpo*Seq1] == 0) {
							DP2[Seq2 + dpo*Seq1] = 1;
							DP2[Seq1 + dpo*Seq2] = 1;
						}
					}

				}

			}
			free(HighEnough);
			free(SubSeq);
			free(MaxScorePos);
			free(FragMaxScore);
			free(PVals);
			free(FragST);
			free(FragEN);
			free(FragScore);

			

		}
		omp_set_num_threads(2);
		return(redonum);
	}

	float MyMathFuncs::DoHMMCycles(int nseed, int SLen, int HMMCycles, int LenStrainSeq0, int NumberXY, int NumberABC, unsigned char *RecodeB, double *TransitionBak, double *EmissionBak,  double *InitPBak, int *LaticePathBak) {
		double bestLike;

		int procs;
		procs = omp_get_num_procs();
		procs = procs / 2-1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel 
		{


			int NumCycles, x, Y, A, Dummy, Z;
			float iVal;


			double *TransitionM2, *EmissionM2, *InitP, *StateCount, *TransitionCount, *LaticeAB, *LaticeXY, *OptXY;
			int *DoneImba, *LaticePath;
			double Imballance, TotCount, Fudge, PathLike, MaxL, rd, PathMax;
			int Imballance2, di, Maxiterations;
			///ReDim TransitionM2(NumberXY - 1, NumberXY - 1), EmissionM2(NumberABC - 1, NumberXY - 1), InitP(NumberXY - 1);
			TransitionM2 = (double*)calloc(NumberXY*NumberXY, sizeof(double));
			EmissionM2 = (double*)calloc(NumberABC*NumberXY, sizeof(double));
			LaticeAB = (double*)calloc((SLen + 1)*NumberXY, sizeof(double));
			LaticeXY = (double*)calloc((SLen + 1)*NumberXY, sizeof(double));
			LaticePath = (int*)calloc((SLen + 4), sizeof(int));
			InitP = (double*)calloc(NumberXY, sizeof(double));
			PathMax = 0;
			srand(nseed + omp_get_thread_num()); // makes sure random numbers are not the same for all threads;
			//Dim DoneImba() As Long, Imballance As Double, Imballance2 As Long
			//Dim LaticePathBak() As Long, StateCount() As Double, TotCount As Double, TransitionCount() As Double, TotXY(), LaticeAB() As Double,
			//dim LaticeXY() As Double, LaticePath() As Long, Fudge As Double, Maxiterations As Long, PathLike As Double, MaxL As Double, PathMax As Double
			bestLike = -1000000;
#pragma omp for 
			for (NumCycles = 0; NumCycles <= HMMCycles; NumCycles++) {
				iVal = 5 / (double)(LenStrainSeq0);
				if (NumberXY > 1) {
					for (x = 0; x <= NumberXY - 1; x++) {
						for (Y = 0; Y <= NumberXY - 1; Y++) {
							if (x == Y)
								TransitionM2[x + Y*NumberXY] = log((double)(1 - iVal));
							else
								TransitionM2[x + Y*NumberXY] = log((double)(iVal) / (double)(NumberXY - 1));

						}
					}
				}
				else
					TransitionM2[0] = log(1);

				rd = rand();
				rd = rd / RAND_MAX;
				Imballance = ((int)(((double)(NumberABC)* rd) + 1)) / 10.0;



				DoneImba = (int*)calloc(NumberABC + 1, sizeof(int));
				di = 0;
				for (x = 0; x <= NumberXY - 1; x++) {

					while (di == 0) {
						rd = rand();
						rd = rd / RAND_MAX;
						Imballance2 = (int)(((NumberXY + 3) * rd)) - 2;
						if (Imballance2 >= 0 && Imballance2 < NumberXY) {
							if (DoneImba[Imballance2] == 0) {
								for (Y = 0; Y <= NumberABC - 1; Y++) {
									if (Y == Imballance2)
										EmissionM2[Y + x*NumberABC] = 1.0 / (double)(NumberABC)+Imballance * 2;
									else
										EmissionM2[Y + x*NumberABC] = 1.0 / (double)(NumberABC)-Imballance;

								}

								DoneImba[Imballance2] = 1;
								break;
							}
						}
					}
				}


				for (x = 0; x <= NumberXY - 1; x++) {

					for (Y = 0; Y <= NumberABC - 1; Y++)
						EmissionM2[Y + x*NumberABC] = log(EmissionM2[Y + x*NumberABC]);
				}


				//double testv2;
				for (x = 0; x <= NumberXY - 1; x++) {
					InitP[x] = log(1.0 / (double)(NumberXY));
					//testv2= InitP[x];

				}

				//LaticeXY = (double*)calloc((SLen+1)*NumberXY, sizeof(double));
				//LaticePath = (int*)calloc((SLen + 4), sizeof(int));
				//ReDim  LaticePath(SLen + 3)

				Fudge = 0.01;

				Maxiterations = 100;




				for (A = 1; A <= Maxiterations; A++) {

					//double *OptXY;
					//Dim Dummy As Long, x As Long, OptXY() As Double
					//ReDim OptXY(NumberXY - 1, NumberXY - 1)
					OptXY = (double*)calloc(NumberXY*NumberXY, sizeof(double));
					
					//ReDim LaticeAB(SLen, NumberXY - 1)
					//ReDim LaticeXY(SLen, NumberXY - 1)
					//double testv;

					for (x = 0; x <= NumberXY - 1; x++) {
						LaticeXY[x*(SLen + 1)] = EmissionM2[RecodeB[0] + x*NumberABC] + InitP[x];
						//testv = LaticeXY[x*(SLen + 1)];

					}



					Dummy = ViterbiCP(SLen, NumberABC, NumberXY, OptXY, RecodeB, LaticeXY, TransitionM2, EmissionM2, LaticeAB);

					MaxL = GetLaticePathP(SLen, NumberXY, LaticeXY, LaticeAB, LaticePath);

					if (PathMax == MaxL) {
						free(OptXY);
						
						break;

					}
					else
						PathMax = MaxL;




					if (NumberXY > 1) {

						TransitionCount = (double*)calloc(NumberXY*NumberXY, sizeof(double));
						StateCount = (double*)calloc(NumberABC*NumberXY, sizeof(double));
						//ReDim TransitionCount(NumberXY - 1, NumberXY - 1)

						//ReDim StateCount(NumberABC - 1, NumberXY - 1)


						Dummy = UpdateCountsP(SLen, NumberABC, NumberXY, LaticePath, RecodeB, TransitionCount, StateCount);

						for (x = 0; x <= NumberXY - 1; x++) {
							TotCount = 0;
							for (Y = 0; Y <= NumberXY - 1; Y++)
								TotCount = TotCount + TransitionCount[x + Y*NumberXY];

							TotCount = TotCount + (Fudge * NumberXY);
							for (Y = 0; Y <= NumberXY - 1; Y++) {
								TransitionCount[x + Y*NumberXY] = TransitionCount[x + Y*NumberXY] + Fudge;
								TransitionM2[x + Y*NumberXY] = log((TransitionCount[x + Y*NumberXY] / TotCount));
							}
						}

						for (Y = 0; Y <= NumberXY - 1; Y++) {
							TotCount = 0;
							for (Z = 0; Z <= NumberABC - 1; Z++)
								TotCount = TotCount + StateCount[Z + Y*NumberABC];

							TotCount = TotCount + (Fudge * NumberABC);
							for (Z = 0; Z <= NumberABC - 1; Z++)
								EmissionM2[Z + Y*NumberABC] = log(((StateCount[Z + Y*NumberABC] + Fudge) / TotCount));
						}
						free(TransitionCount);
						free(StateCount);
					}
					else
						TransitionM2[0] = log(1.0);

					free(OptXY);
					//free(LaticeAB);
					//free(LaticeXY);
				}

				if (MaxL > bestLike) {
#pragma omp critical
					{
						bestLike = MaxL;
						//ReDim TransitionBak(NumberXY - 1, NumberXY - 1), EmissionBak(NumberABC - 1, NumberXY - 1), InitPBak(NumberXY - 1), LaticePathBak(SLen)
						for (x = 0; x <= NumberXY - 1; x++) {
							InitPBak[x] = InitP[x];
							for (Y = 0; Y <= NumberXY - 1; Y++)
								TransitionBak[x + Y*NumberXY] = TransitionM2[x + Y*NumberXY];

							for (Y = 0; Y <= NumberABC - 1; Y++)
								EmissionBak[Y + x*NumberABC] = EmissionM2[Y + x*NumberABC];

						}
						/*for (x = 0; x <= SLen; x++)
							LaticePathBak[x] = LaticePath[x];*/

						
						memcpy(LaticePathBak, LaticePath, (SLen+1) * sizeof(int));

					}
				}
				
				free(DoneImba);

			}
			free(LaticePath);
			free(LaticeAB);
			free(LaticeXY);
			free(TransitionM2);
			free(EmissionM2);
			free(InitP);
		}
		omp_set_num_threads(2);
		return(bestLike);
	}

	float MyMathFuncs::DoHMMCyclesDetermin(int nseed, int SLen, int HMMCycles, int LenStrainSeq0, int NumberXY, int NumberABC, unsigned char *RecodeB, double *TransitionBak, double *EmissionBak, double *InitPBak, int *LaticePathBak) {
		double bestLike;

		int procs;
		procs = omp_get_num_procs();
		procs = procs / 2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel
		{


			int NumCycles, x, Y, A, Dummy, Z;
			float iVal;


			double *TransitionM2, *EmissionM2, *InitP, *StateCount, *TransitionCount, *LaticeAB, *LaticeXY, *OptXY;
			int *DoneImba, *LaticePath;
			double Imballance, TotCount, Fudge, PathLike, MaxL, rd, PathMax;
			int Imballance2, di, Maxiterations;
			///ReDim TransitionM2(NumberXY - 1, NumberXY - 1), EmissionM2(NumberABC - 1, NumberXY - 1), InitP(NumberXY - 1);



			//Dim DoneImba() As Long, Imballance As Double, Imballance2 As Long
			//Dim LaticePathBak() As Long, StateCount() As Double, TotCount As Double, TransitionCount() As Double, TotXY(), LaticeAB() As Double,
			//dim LaticeXY() As Double, LaticePath() As Long, Fudge As Double, Maxiterations As Long, PathLike As Double, MaxL As Double, PathMax As Double
			bestLike = -1000000;
			PathMax = 0;
#pragma omp for 
			for (NumCycles = 0; NumCycles <= HMMCycles; NumCycles++) {
				
				TransitionM2 = (double*)calloc(NumberXY*NumberXY, sizeof(double));
				EmissionM2 = (double*)calloc(NumberABC*NumberXY, sizeof(double));
				LaticeAB = (double*)calloc((SLen + 1)*NumberXY, sizeof(double));
				LaticeXY = (double*)calloc((SLen + 1)*NumberXY, sizeof(double));
				LaticePath = (int*)calloc((SLen + 4), sizeof(int));
				InitP = (double*)calloc(NumberXY, sizeof(double));
				srand((unsigned int )(nseed + NumCycles)); // makes sure random numbers are not the same for all threads;
				iVal = 5 / (double)(LenStrainSeq0);
				if (NumberXY > 1) {
					for (x = 0; x <= NumberXY - 1; x++) {
						for (Y = 0; Y <= NumberXY - 1; Y++) {
							if (x == Y)
								TransitionM2[x + Y*NumberXY] = log((double)(1 - iVal));
							else
								TransitionM2[x + Y*NumberXY] = log((double)(iVal) / (double)(NumberXY - 1));

						}
					}
				}
				else
					TransitionM2[0] = log(1);

				rd = rand();
				rd = rd / RAND_MAX;
				Imballance = ((int)(((double)(NumberABC)* rd) + 1)) / 10.0;



				DoneImba = (int*)calloc(NumberABC + 1, sizeof(int));
				di = 0;
				for (x = 0; x <= NumberXY - 1; x++) {

					while (di == 0) {
						rd = rand();
						rd = rd / RAND_MAX;
						Imballance2 = (int)(((NumberXY + 3) * rd)) - 2;
						if (Imballance2 >= 0 && Imballance2 < NumberXY) {
							if (DoneImba[Imballance2] == 0) {
								for (Y = 0; Y <= NumberABC - 1; Y++) {
									if (Y == Imballance2)
										EmissionM2[Y + x*NumberABC] = 1.0 / (double)(NumberABC)+Imballance * 2;
									else
										EmissionM2[Y + x*NumberABC] = 1.0 / (double)(NumberABC)-Imballance;

								}

								DoneImba[Imballance2] = 1;
								break;
							}
						}
					}
				}


				for (x = 0; x <= NumberXY - 1; x++) {

					for (Y = 0; Y <= NumberABC - 1; Y++)
						EmissionM2[Y + x*NumberABC] = log(EmissionM2[Y + x*NumberABC]);
				}


				//double testv2;
				for (x = 0; x <= NumberXY - 1; x++) {
					InitP[x] = log(1.0 / (double)(NumberXY));
					//testv2= InitP[x];

				}

				//LaticeXY = (double*)calloc((SLen+1)*NumberXY, sizeof(double));
				//LaticePath = (int*)calloc((SLen + 4), sizeof(int));
				//ReDim  LaticePath(SLen + 3)

				Fudge = 0.01;

				Maxiterations = 100;




				for (A = 1; A <= Maxiterations; A++) {

					//double *OptXY;
					//Dim Dummy As Long, x As Long, OptXY() As Double
					//ReDim OptXY(NumberXY - 1, NumberXY - 1)
					OptXY = (double*)calloc(NumberXY*NumberXY, sizeof(double));

					//ReDim LaticeAB(SLen, NumberXY - 1)
					//ReDim LaticeXY(SLen, NumberXY - 1)
					//double testv;

					for (x = 0; x <= NumberXY - 1; x++) {
						LaticeXY[x*(SLen + 1)] = EmissionM2[RecodeB[0] + x*NumberABC] + InitP[x];
						//testv = LaticeXY[x*(SLen + 1)];

					}



					Dummy = ViterbiCP(SLen, NumberABC, NumberXY, OptXY, RecodeB, LaticeXY, TransitionM2, EmissionM2, LaticeAB);

					MaxL = GetLaticePathP(SLen, NumberXY, LaticeXY, LaticeAB, LaticePath);

					if (PathMax == MaxL) {
						free(OptXY);

						break;

					}
					else
						PathMax = MaxL;




					if (NumberXY > 1) {

						TransitionCount = (double*)calloc(NumberXY*NumberXY, sizeof(double));
						StateCount = (double*)calloc(NumberABC*NumberXY, sizeof(double));
						//ReDim TransitionCount(NumberXY - 1, NumberXY - 1)

						//ReDim StateCount(NumberABC - 1, NumberXY - 1)


						Dummy = UpdateCountsP(SLen, NumberABC, NumberXY, LaticePath, RecodeB, TransitionCount, StateCount);

						for (x = 0; x <= NumberXY - 1; x++) {
							TotCount = 0;
							for (Y = 0; Y <= NumberXY - 1; Y++)
								TotCount = TotCount + TransitionCount[x + Y*NumberXY];

							TotCount = TotCount + (Fudge * NumberXY);
							for (Y = 0; Y <= NumberXY - 1; Y++) {
								TransitionCount[x + Y*NumberXY] = TransitionCount[x + Y*NumberXY] + Fudge;
								TransitionM2[x + Y*NumberXY] = log((TransitionCount[x + Y*NumberXY] / TotCount));
							}
						}

						for (Y = 0; Y <= NumberXY - 1; Y++) {
							TotCount = 0;
							for (Z = 0; Z <= NumberABC - 1; Z++)
								TotCount = TotCount + StateCount[Z + Y*NumberABC];

							TotCount = TotCount + (Fudge * NumberABC);
							for (Z = 0; Z <= NumberABC - 1; Z++)
								EmissionM2[Z + Y*NumberABC] = log(((StateCount[Z + Y*NumberABC] + Fudge) / TotCount));
						}
						free(TransitionCount);
						free(StateCount);
					}
					else
						TransitionM2[0] = log(1.0);

					free(OptXY);
					//free(LaticeAB);
					//free(LaticeXY);
				}

				if (MaxL > bestLike) {
#pragma omp critical
					{
						bestLike = MaxL;
						//ReDim TransitionBak(NumberXY - 1, NumberXY - 1), EmissionBak(NumberABC - 1, NumberXY - 1), InitPBak(NumberXY - 1), LaticePathBak(SLen)
						for (x = 0; x <= NumberXY - 1; x++) {
							InitPBak[x] = InitP[x];
							for (Y = 0; Y <= NumberXY - 1; Y++)
								TransitionBak[x + Y*NumberXY] = TransitionM2[x + Y*NumberXY];

							for (Y = 0; Y <= NumberABC - 1; Y++)
								EmissionBak[Y + x*NumberABC] = EmissionM2[Y + x*NumberABC];

						}
						/*for (x = 0; x <= SLen; x++)
						LaticePathBak[x] = LaticePath[x];*/


						memcpy(LaticePathBak, LaticePath, (SLen + 1) * sizeof(int));

					}
				}

				free(DoneImba);
				free(LaticePath);
				free(LaticeAB);
				free(LaticeXY);
				free(TransitionM2);
				free(EmissionM2);
				free(InitP);
			}

		}
		omp_set_num_threads(2);
		return(bestLike);
	}
	

	float MyMathFuncs::DoHMMCyclesSerial(int nseed, int SLen, int HMMCycles, int LenStrainSeq0, int NumberXY, int NumberABC, unsigned char *RecodeB, double *TransitionBak, double *EmissionBak, double *InitPBak, int *LaticePathBak) {
		double bestLike;

		/*int procs;
		procs = omp_get_num_procs();
		procs = procs / 2 - 1;
		if (procs < 3)
		procs = 3;
		omp_set_num_threads(procs);*/
		//#pragma omp parallel
		//		{


		int NumCycles, x, Y, A, Dummy, Z;
		float iVal;


		double *TransitionM2, *EmissionM2, *InitP, *StateCount, *TransitionCount, *LaticeAB, *LaticeXY, *OptXY;
		int *DoneImba, *LaticePath;
		double Imballance, TotCount, Fudge, PathLike, MaxL, rd, PathMax;
		int Imballance2, di, Maxiterations;
		///ReDim TransitionM2(NumberXY - 1, NumberXY - 1), EmissionM2(NumberABC - 1, NumberXY - 1), InitP(NumberXY - 1);



		//Dim DoneImba() As Long, Imballance As Double, Imballance2 As Long
		//Dim LaticePathBak() As Long, StateCount() As Double, TotCount As Double, TransitionCount() As Double, TotXY(), LaticeAB() As Double,
		//dim LaticeXY() As Double, LaticePath() As Long, Fudge As Double, Maxiterations As Long, PathLike As Double, MaxL As Double, PathMax As Double
		bestLike = -1000000;
		PathMax = 0;
		//#pragma omp for 
		TransitionM2 = (double*)calloc(NumberXY*NumberXY, sizeof(double));
		EmissionM2 = (double*)calloc(NumberABC*NumberXY, sizeof(double));
		LaticeAB = (double*)calloc((SLen + 1)*NumberXY, sizeof(double));
		LaticeXY = (double*)calloc((SLen + 1)*NumberXY, sizeof(double));
		LaticePath = (int*)calloc((SLen + 4), sizeof(int));
		InitP = (double*)calloc(NumberXY, sizeof(double));
		srand((unsigned int)(nseed));
		for (NumCycles = 0; NumCycles <= HMMCycles; NumCycles++) {


			//srand((unsigned int)(nseed + NumCycles)); // make sure random numbers are not the same for all threads;
			iVal = 5 / (double)(LenStrainSeq0);
			if (NumberXY > 1) {
				for (x = 0; x <= NumberXY - 1; x++) {
					for (Y = 0; Y <= NumberXY - 1; Y++) {
						if (x == Y)
							TransitionM2[x + Y*NumberXY] = log((double)(1 - iVal));
						else
							TransitionM2[x + Y*NumberXY] = log((double)(iVal) / (double)(NumberXY - 1));

					}
				}
			}
			else
				TransitionM2[0] = log(1);

			rd = rand();
			rd = rd / RAND_MAX;
			Imballance = ((int)(((double)(NumberABC)* rd) + 1)) / 10.0;



			DoneImba = (int*)calloc(NumberABC + 1, sizeof(int));
			di = 0;
			for (x = 0; x <= NumberXY - 1; x++) {

				while (di == 0) {
					rd = rand();
					rd = rd / RAND_MAX;
					Imballance2 = (int)(((NumberXY + 3) * rd)) - 2;
					if (Imballance2 >= 0 && Imballance2 < NumberXY) {
						if (DoneImba[Imballance2] == 0) {
							for (Y = 0; Y <= NumberABC - 1; Y++) {
								if (Y == Imballance2)
									EmissionM2[Y + x*NumberABC] = 1.0 / (double)(NumberABC)+Imballance * 2;
								else
									EmissionM2[Y + x*NumberABC] = 1.0 / (double)(NumberABC)-Imballance;

							}

							DoneImba[Imballance2] = 1;
							break;
						}
					}
				}
			}


			for (x = 0; x <= NumberXY - 1; x++) {

				for (Y = 0; Y <= NumberABC - 1; Y++)
					EmissionM2[Y + x*NumberABC] = log(EmissionM2[Y + x*NumberABC]);
			}


			//double testv2;
			for (x = 0; x <= NumberXY - 1; x++) {
				InitP[x] = log(1.0 / (double)(NumberXY));
				//testv2= InitP[x];

			}

			//LaticeXY = (double*)calloc((SLen+1)*NumberXY, sizeof(double));
			//LaticePath = (int*)calloc((SLen + 4), sizeof(int));
			//ReDim  LaticePath(SLen + 3)

			Fudge = 0.01;

			Maxiterations = 100;




			for (A = 1; A <= Maxiterations; A++) {

				//double *OptXY;
				//Dim Dummy As Long, x As Long, OptXY() As Double
				//ReDim OptXY(NumberXY - 1, NumberXY - 1)
				OptXY = (double*)calloc(NumberXY*NumberXY, sizeof(double));

				//ReDim LaticeAB(SLen, NumberXY - 1)
				//ReDim LaticeXY(SLen, NumberXY - 1)
				//double testv;

				for (x = 0; x <= NumberXY - 1; x++) {
					LaticeXY[x*(SLen + 1)] = EmissionM2[RecodeB[0] + x*NumberABC] + InitP[x];
					//testv = LaticeXY[x*(SLen + 1)];

				}



				Dummy = ViterbiCP(SLen, NumberABC, NumberXY, OptXY, RecodeB, LaticeXY, TransitionM2, EmissionM2, LaticeAB);

				MaxL = GetLaticePathP(SLen, NumberXY, LaticeXY, LaticeAB, LaticePath);

				if (PathMax == MaxL) {
					free(OptXY);

					break;

				}
				else
					PathMax = MaxL;




				if (NumberXY > 1) {

					TransitionCount = (double*)calloc(NumberXY*NumberXY, sizeof(double));
					StateCount = (double*)calloc(NumberABC*NumberXY, sizeof(double));
					//ReDim TransitionCount(NumberXY - 1, NumberXY - 1)

					//ReDim StateCount(NumberABC - 1, NumberXY - 1)


					Dummy = UpdateCountsP(SLen, NumberABC, NumberXY, LaticePath, RecodeB, TransitionCount, StateCount);

					for (x = 0; x <= NumberXY - 1; x++) {
						TotCount = 0;
						for (Y = 0; Y <= NumberXY - 1; Y++)
							TotCount = TotCount + TransitionCount[x + Y*NumberXY];

						TotCount = TotCount + (Fudge * NumberXY);
						for (Y = 0; Y <= NumberXY - 1; Y++) {
							TransitionCount[x + Y*NumberXY] = TransitionCount[x + Y*NumberXY] + Fudge;
							TransitionM2[x + Y*NumberXY] = log((TransitionCount[x + Y*NumberXY] / TotCount));
						}
					}

					for (Y = 0; Y <= NumberXY - 1; Y++) {
						TotCount = 0;
						for (Z = 0; Z <= NumberABC - 1; Z++)
							TotCount = TotCount + StateCount[Z + Y*NumberABC];

						TotCount = TotCount + (Fudge * NumberABC);
						for (Z = 0; Z <= NumberABC - 1; Z++)
							EmissionM2[Z + Y*NumberABC] = log(((StateCount[Z + Y*NumberABC] + Fudge) / TotCount));
					}
					free(TransitionCount);
					free(StateCount);
				}
				else
					TransitionM2[0] = log(1.0);

				free(OptXY);
				//free(LaticeAB);
				//free(LaticeXY);
			}

			if (MaxL > bestLike) {
				//#pragma omp critical
				//					{
				bestLike = MaxL;
				//ReDim TransitionBak(NumberXY - 1, NumberXY - 1), EmissionBak(NumberABC - 1, NumberXY - 1), InitPBak(NumberXY - 1), LaticePathBak(SLen)
				for (x = 0; x <= NumberXY - 1; x++) {
					InitPBak[x] = InitP[x];
					for (Y = 0; Y <= NumberXY - 1; Y++)
						TransitionBak[x + Y*NumberXY] = TransitionM2[x + Y*NumberXY];

					for (Y = 0; Y <= NumberABC - 1; Y++)
						EmissionBak[Y + x*NumberABC] = EmissionM2[Y + x*NumberABC];

				}
				/*for (x = 0; x <= SLen; x++)
				LaticePathBak[x] = LaticePath[x];*/


				memcpy(LaticePathBak, LaticePath, (SLen + 1) * sizeof(int));

			}
			//				}

			free(DoneImba);
		}

		free(LaticePath);
		free(LaticeAB);
		free(LaticeXY);
		free(TransitionM2);
		free(EmissionM2);
		free(InitP);

		//		}
		//		omp_set_num_threads(2);
		return(bestLike);
	}

	float MyMathFuncs::DoHMMCyclesPara(int nseed, int SLen, int HMMCycles, int LenStrainSeq0, int NumberXY, int NumberABC, unsigned char *RecodeB, double *TransitionBak, double *EmissionBak, double *InitPBak, int *LaticePathBak) {
		double bestLike;

		int procs;
		procs = omp_get_num_procs();
		procs = procs / 2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel
		{


			int NumCycles, x, Y, A, Dummy, Z;
			float iVal;


			double *TransitionM2, *EmissionM2, *InitP, *StateCount, *TransitionCount, *LaticeAB, *LaticeXY, *OptXY;
			int *DoneImba, *LaticePath;
			double Imballance, TotCount, Fudge, PathLike, MaxL, rd, PathMax;
			int Imballance2, di, Maxiterations;
			///ReDim TransitionM2(NumberXY - 1, NumberXY - 1), EmissionM2(NumberABC - 1, NumberXY - 1), InitP(NumberXY - 1);

			std::mt19937 rng;
			std::uniform_int_distribution<int> distribution(0, RAND_MAX);
			//Dim DoneImba() As Long, Imballance As Double, Imballance2 As Long
			//Dim LaticePathBak() As Long, StateCount() As Double, TotCount As Double, TransitionCount() As Double, TotXY(), LaticeAB() As Double,
			//dim LaticeXY() As Double, LaticePath() As Long, Fudge As Double, Maxiterations As Long, PathLike As Double, MaxL As Double, PathMax As Double
			bestLike = -1000000;
			PathMax = 0;
			
			TransitionM2 = (double*)calloc(NumberXY*NumberXY, sizeof(double));
			EmissionM2 = (double*)calloc(NumberABC*NumberXY, sizeof(double));
			LaticeAB = (double*)calloc((SLen + 1)*NumberXY, sizeof(double));
			LaticeXY = (double*)calloc((SLen + 1)*NumberXY, sizeof(double));
			LaticePath = (int*)calloc((SLen + 4), sizeof(int));
			InitP = (double*)calloc(NumberXY, sizeof(double));

#pragma omp for 
			

			for (NumCycles = 0; NumCycles <= HMMCycles; NumCycles++) {

				
				
				//srand((unsigned int)(nseed + NumCycles)); // make sure random numbers are not the same for all threads;
				rng.seed((unsigned int)(nseed + NumCycles));
				iVal = 5 / (double)(LenStrainSeq0);
				if (NumberXY > 1) {
					for (x = 0; x <= NumberXY - 1; x++) {
						for (Y = 0; Y <= NumberXY - 1; Y++) {
							if (x == Y)
								TransitionM2[x + Y*NumberXY] = log((double)(1 - iVal));
							else
								TransitionM2[x + Y*NumberXY] = log((double)(iVal) / (double)(NumberXY - 1));

						}
					}
				}
				else
					TransitionM2[0] = log(1);

				//rd = rand();
				//rd = rng();
				//
				rd = distribution(rng);
				rd = rd / RAND_MAX;
				Imballance = ((int)(((double)(NumberABC)* rd) + 1)) / 10.0;



				DoneImba = (int*)calloc(NumberABC + 1, sizeof(int));
				di = 0;
				for (x = 0; x <= NumberXY - 1; x++) {

					while (di == 0) {
						//rd = rand();
						//rd = rng();
						//rd = rd / RAND_MAX;
						
						rd = distribution(rng);
						rd = rd / RAND_MAX;
						Imballance2 = (int)(((NumberXY + 3) * rd)) - 2;
						if (Imballance2 >= 0 && Imballance2 < NumberXY) {
							if (DoneImba[Imballance2] == 0) {
								for (Y = 0; Y <= NumberABC - 1; Y++) {
									if (Y == Imballance2)
										EmissionM2[Y + x*NumberABC] = 1.0 / (double)(NumberABC)+Imballance * 2;
									else
										EmissionM2[Y + x*NumberABC] = 1.0 / (double)(NumberABC)-Imballance;

								}

								DoneImba[Imballance2] = 1;
								break;
							}
						}
					}
				}


				for (x = 0; x <= NumberXY - 1; x++) {

					for (Y = 0; Y <= NumberABC - 1; Y++)
						EmissionM2[Y + x*NumberABC] = log(EmissionM2[Y + x*NumberABC]);
				}


				//double testv2;
				for (x = 0; x <= NumberXY - 1; x++) {
					InitP[x] = log(1.0 / (double)(NumberXY));
					//testv2= InitP[x];

				}

				//LaticeXY = (double*)calloc((SLen+1)*NumberXY, sizeof(double));
				//LaticePath = (int*)calloc((SLen + 4), sizeof(int));
				//ReDim  LaticePath(SLen + 3)

				Fudge = 0.01;

				Maxiterations = 100;




				for (A = 1; A <= Maxiterations; A++) {

					//double *OptXY;
					//Dim Dummy As Long, x As Long, OptXY() As Double
					//ReDim OptXY(NumberXY - 1, NumberXY - 1)
					OptXY = (double*)calloc(NumberXY*NumberXY, sizeof(double));

					//ReDim LaticeAB(SLen, NumberXY - 1)
					//ReDim LaticeXY(SLen, NumberXY - 1)
					//double testv;

					for (x = 0; x <= NumberXY - 1; x++) {
						LaticeXY[x*(SLen + 1)] = EmissionM2[RecodeB[0] + x*NumberABC] + InitP[x];
						//testv = LaticeXY[x*(SLen + 1)];

					}



					Dummy = ViterbiCP(SLen, NumberABC, NumberXY, OptXY, RecodeB, LaticeXY, TransitionM2, EmissionM2, LaticeAB);

					MaxL = GetLaticePathP(SLen, NumberXY, LaticeXY, LaticeAB, LaticePath);

					if (PathMax == MaxL) {
						free(OptXY);

						break;

					}
					else
						PathMax = MaxL;




					if (NumberXY > 1) {

						TransitionCount = (double*)calloc(NumberXY*NumberXY, sizeof(double));
						StateCount = (double*)calloc(NumberABC*NumberXY, sizeof(double));
						//ReDim TransitionCount(NumberXY - 1, NumberXY - 1)

						//ReDim StateCount(NumberABC - 1, NumberXY - 1)


						Dummy = UpdateCountsP(SLen, NumberABC, NumberXY, LaticePath, RecodeB, TransitionCount, StateCount);

						for (x = 0; x <= NumberXY - 1; x++) {
							TotCount = 0;
							for (Y = 0; Y <= NumberXY - 1; Y++)
								TotCount = TotCount + TransitionCount[x + Y*NumberXY];

							TotCount = TotCount + (Fudge * NumberXY);
							for (Y = 0; Y <= NumberXY - 1; Y++) {
								TransitionCount[x + Y*NumberXY] = TransitionCount[x + Y*NumberXY] + Fudge;
								TransitionM2[x + Y*NumberXY] = log((TransitionCount[x + Y*NumberXY] / TotCount));
							}
						}

						for (Y = 0; Y <= NumberXY - 1; Y++) {
							TotCount = 0;
							for (Z = 0; Z <= NumberABC - 1; Z++)
								TotCount = TotCount + StateCount[Z + Y*NumberABC];

							TotCount = TotCount + (Fudge * NumberABC);
							for (Z = 0; Z <= NumberABC - 1; Z++)
								EmissionM2[Z + Y*NumberABC] = log(((StateCount[Z + Y*NumberABC] + Fudge) / TotCount));
						}
						free(TransitionCount);
						free(StateCount);
					}
					else
						TransitionM2[0] = log(1.0);

					free(OptXY);
					//free(LaticeAB);
					//free(LaticeXY);
				}

				if (MaxL > bestLike) {
#pragma omp critical
					{
						bestLike = MaxL;
						//ReDim TransitionBak(NumberXY - 1, NumberXY - 1), EmissionBak(NumberABC - 1, NumberXY - 1), InitPBak(NumberXY - 1), LaticePathBak(SLen)
						for (x = 0; x <= NumberXY - 1; x++) {
							InitPBak[x] = InitP[x];
							for (Y = 0; Y <= NumberXY - 1; Y++)
								TransitionBak[x + Y*NumberXY] = TransitionM2[x + Y*NumberXY];

							for (Y = 0; Y <= NumberABC - 1; Y++)
								EmissionBak[Y + x*NumberABC] = EmissionM2[Y + x*NumberABC];

						}
						/*for (x = 0; x <= SLen; x++)
						LaticePathBak[x] = LaticePath[x];*/


						memcpy(LaticePathBak, LaticePath, (SLen + 1) * sizeof(int));

					}
				}

				free(DoneImba);
				
			}
			free(LaticePath);
			free(LaticeAB);
			free(LaticeXY);
			free(TransitionM2);
			free(EmissionM2);
			free(InitP);
		}
		omp_set_num_threads(2);
		return(bestLike);
	}

	int MyMathFuncs::RecodeNucsLong(int Y, int LSeq, int UBRecoded, int UBReplace, short int *tSeqnum, unsigned char *NucMat,unsigned char *Replace, unsigned char *Recoded)

	{
		int x, off1, off2;
		unsigned char NN;
		off1 = UBRecoded + 1;
		off2 = UBReplace + 1;
		for (x = 1; x <= LSeq; x++) {

			NN = tSeqnum[x];
			NN = NucMat[NN];
			Recoded[x + Y*off1] = Replace[x + NN*off2];
		}
		return(1);
	}
	int MyMathFuncs::RecodeNucs(int NextNo,  int LS, int UBNC, int *NucCount, int UBR, unsigned char *Replace)
	{
		int x, H, W, Y, CycleX, os1, os2;
		os1 = UBNC + 1;
		os2 = UBR + 1;
		for (x = 1; x <= LS; x++) {
			CycleX = 1;
			do {
				H = 0;
				W = 0;
				for (Y = 1; Y <= 4; Y++) {
					if (H < NucCount[x + os1*Y]) {
						H = NucCount[x + os1*Y];
						W = Y;
					}
				}
				if (H == 0) {
					break;
				}
				NucCount[x + os1*W] = -1;
				Replace[x + os2*W] = CycleX;
				if (H == NextNo + 1) {
					break;
				}
				CycleX = CycleX + 1;
			} while (true);
		}
			return(1);

	}

int MyMathFuncs::FastRecCheckP(int CircularFlag, int DoAllFlag, int MCCorrection, int MCFlag,int EarlyBale,double UCThresh, double LowestProb, int NextNo, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, int UBXSNW, unsigned char *CS,short int *Seqnum,int Seq1, int Seq2, int Seq3, int LenStrainSeq, int XoverWindow, short int XOverWindowX, short int *XoverSeqNum, char *XoverSeqNumW, int UBXOHN, int *XOverHomologyNum, unsigned char *FSSRDP, int ProbEstimateInFileFlag,int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact, double *BQPV){

	//Dim GoOn As Long, NumSame As Long
	//Dim XOverLen As Long, A As Long, Dummy As Variant, OldX As Long

	//Dim AFact As Double, B As Long, BE As Long, EN As Long, NCommon As Long, StartPosX As Long
	//Dim NextPosX As Long
	//Dim oProbXOver As Double
	//Dim FindCycle As Integer, SeqDaughter As Long, SeqMinorP As Long, Temp As Integer, EndFlag As Long, HighHomol As Long, MedHomol As Long, LowHomol As Long
	//Dim x As Long, NumDifferent As Long

	//Dim SLen As Long
	//Dim AH(2) As Long

	int x, Dummy, tds12, tds23, tds13, Temp;
	int HighHomol, MedHomol, LowHomol, ActiveSeq, ActiveMajorP, ActiveMinorP, SeqDaughter, SeqMinorP;
	int StartPosX, FindCycle, OldX, NextPosX, NumDifferent, NumSame, GoOn;// , test1, test2, test3;
	double AFact, ProbabilityXOver,IndProb, ah1,ah2,ah3;
	int *AH, *NCommon, *XOverLength, *BE, *EN, *EndFlag;
	float tt;
	double *AvHomol;
	AH = (int*)calloc(4, sizeof(int));
	
	IndProb = 0.0;
	/*test1 = CS[1];
	test2 = CS[500];
	test3 = CS[UBCS + 2];*/



	int LenXoverSeq;
    
	LenXoverSeq = FindSubSeqPB3 (&AH[0], UBFSSRDP, XoverWindow, UBCS, LenStrainseq0, NextNo, Seq1, Seq2, Seq3, &CS[0], UBXSNW, &XoverSeqNumW[0], &FSSRDP[0]);
	//return(LenXoverSeq);
    
	if (LenXoverSeq < XoverWindow * 2) {
		free(AH);

		CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
			
		return(0);
	}
    
	if (AH[0] < XoverWindow / 3 || AH[1] < XoverWindow / 3 || AH[2] < XoverWindow / 3) {
		free(AH);
		
		CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
		return(0);
	}
    if (AH[0] == 0 || AH[1] == 0 || AH[2] == 0 ){
		CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
		free(AH);
		
		return(0);
    }
   
	
	AvHomol = (double*)calloc(4, sizeof(double));

	AvHomol[1] = ((double)(AH[0])) / ((double)(LenXoverSeq));
    AvHomol[2] = ((double)(AH[1])) / ((double)(LenXoverSeq));
    AvHomol[3] = ((double)(AH[2])) / ((double)(LenXoverSeq));
	ah1 = AvHomol[1];
	ah2 = AvHomol[2];
	ah3 = AvHomol[3];
    //AvHomol(3) = (CLng(AvHomol(3) * 10000)) / 10000 
    //AvHomol(1) = (CLng(AvHomol(1) * 10000)) / 10000 
    //AvHomol(2) = (CLng(AvHomol(2) * 10000)) / 10000 
    
    /*If MCFlag = 0 Then
        MCC = MCCorrection
    Else
        MCC = 1
    End If*/
    
	LenXoverSeq = abs(LenXoverSeq) - 1;
    

    
    //Work out identities (7.1/21)

    
    if (UBTD > 0){
		tds12 = Seq1 + Seq2*(UBTD + 1);
		tds13 = Seq1 + Seq3*(UBTD + 1);
		tds23 = Seq2 + Seq3*(UBTD + 1);
		if (TreeDistance[tds12] >= TreeDistance[tds13] && TreeDistance[tds12] >= TreeDistance[tds23])
			HighHomol = 1;
		else if (TreeDistance[tds13] >= TreeDistance[tds12] && TreeDistance[tds13] >= TreeDistance[tds23])
			HighHomol = 2;
		else if (TreeDistance[tds23] >= TreeDistance[tds12] && TreeDistance[tds23] >= TreeDistance[tds13])
			HighHomol = 3;
        
	}
    else{
		tds12 = Seq1 + Seq2*(UBD + 1);
		tds13 = Seq1 + Seq3*(UBD + 1);
		tds23 = Seq2 + Seq3*(UBD + 1);
		if (Distance[tds12] >= Distance[tds13] && Distance[tds12] >= Distance[tds23])
			HighHomol = 1;
		else if (Distance[tds13] >= Distance[tds12] && Distance[tds13] >= Distance[tds23])
			HighHomol = 2;
		else if (Distance[tds23] >= Distance[tds12] && Distance[tds23] >= Distance[tds13])
			HighHomol = 3;
   }

	StartPosX = XOHomologyP2(HighHomol, LenStrainSeq, LenXoverSeq, XoverWindow, XoverSeqNumW, XOverHomologyNum);


	if (StartPosX == 0 && ShortOutFlag != 3) {
		CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
		free(AH);
		free(AvHomol);
		return(0);
	}
    
    if (UBD > 0){
		if (Seq1 <= UBD && Seq2 <= UBD && Seq3 <= UBD){
			tds12 = Seq1 + Seq2*(UBD + 1);
			tds13 = Seq1 + Seq3*(UBD + 1);
			tds23 = Seq2 + Seq3*(UBD + 1);
			if (AvHomol[1] == AvHomol[2] && AvHomol[1] == AvHomol[3]){
				if (Distance[tds12] >= Distance[tds13] && Distance[tds12] >= Distance[tds23]){
					if (Distance[tds13] > Distance[tds23]) {
						AvHomol[2] = AvHomol[2] - 0.00001;
						AvHomol[3] = AvHomol[3] - 0.00002;
					}
					else {
						AvHomol[2] = AvHomol[2] - 0.00002;
						AvHomol[3] = AvHomol[3] - 0.00001;
					}
				}
				else if (Distance[tds13] >= Distance[tds12] && Distance[tds13] >= Distance[tds23]){
					if (Distance[tds12] > Distance[tds23]) {
						AvHomol[1] = AvHomol[1] - 0.00001;
						AvHomol[3] = AvHomol[3] - 0.00002;
					}
					else {
						AvHomol[1] = AvHomol[1] - 0.00002;
						AvHomol[3] = AvHomol[3] - 0.00001;
					}
				}
				else{
					if (Distance[tds12] > Distance[tds13]) {
						AvHomol[1] = AvHomol[1] - 0.00001;
						AvHomol[2] = AvHomol[2] - 0.00002;
					}
					else {
						AvHomol[1] = AvHomol[1] - 0.00002;
						AvHomol[2] = AvHomol[2] - 0.00001;
					}
				}
			}
			else if (AvHomol[1] == AvHomol[2]){
				if (Distance[tds12] > Distance[tds13])
					AvHomol[2] = AvHomol[2] - 0.00001;
				else
					AvHomol[1] = AvHomol[1] - 0.00001;
			}
			else if (AvHomol[1] == AvHomol[3]){
				if (Distance[tds12] > Distance[tds23])
					AvHomol[3] = AvHomol[3] - 0.00001;
				else
					AvHomol[1] = AvHomol[1] - 0.00001;
				
			}
			else if (AvHomol[2] == AvHomol[3]){
				if (Distance[tds13] > Distance[tds23])
					AvHomol[3] = AvHomol[3] - 0.00001;
				else
					AvHomol[2] = AvHomol[2] - 0.00001;
				
			}
    
		}
	}
    
    
    
    if (AvHomol[1] >= AvHomol[2] && AvHomol[1] >= AvHomol[3]){
		HighHomol = 1;

		if (AvHomol[2] >= AvHomol[3]) {
			MedHomol = 2;
			LowHomol = 3;
			ActiveSeq = Seq1;
			ActiveMajorP = Seq2;
			ActiveMinorP = Seq3;
			SeqDaughter = 0;
			SeqMinorP = 2;
		}
        
		else {
			MedHomol = 3; 
			LowHomol = 2;
			ActiveSeq = Seq2;
			ActiveMajorP = Seq1;
			ActiveMinorP = Seq3;
			SeqDaughter = 1;
			SeqMinorP = 2;

		}
	}
    else if (AvHomol[2] >= AvHomol[1] && AvHomol[2] >= AvHomol[3]){
		HighHomol = 2;

		if (AvHomol[1] >= AvHomol[3]) {
			MedHomol = 1;
			LowHomol = 3;
			ActiveSeq = Seq1;
			ActiveMajorP = Seq3;
			ActiveMinorP = Seq2;
			SeqDaughter = 0; 
			SeqMinorP = 1;
		}
		else {
			MedHomol = 3; 
			LowHomol = 1;
			ActiveSeq = Seq3;
			ActiveMajorP = Seq1;
			ActiveMinorP = Seq2;
			SeqDaughter = 2; 
			SeqMinorP = 1;

		}
	}
    else if ( AvHomol[3] >= AvHomol[1] && AvHomol[3] >= AvHomol[2] ){
		HighHomol = 3;
   
		if (AvHomol[1] >= AvHomol[2]) {
			MedHomol = 1; 
			LowHomol = 2;
			ActiveSeq = Seq2; 
			ActiveMajorP = Seq3;
			ActiveMinorP = Seq1;
			SeqDaughter = 1;
			SeqMinorP = 0;
		}
		else {
			MedHomol = 2;
			LowHomol = 1;
			ActiveSeq = Seq3;
			ActiveMajorP = Seq2;
			ActiveMinorP = Seq1;
			SeqDaughter = 2;
			SeqMinorP = 0;

		}

    }
   
	FindCycle = 0;
    
    
	OldX = -1;
    
	NCommon = (int*)calloc(1, sizeof(int));
	XOverLength = (int*)calloc(1, sizeof(int));
	BE = (int*)calloc(1, sizeof(int));
	EN = (int*)calloc(1, sizeof(int));
	EndFlag = (int*)calloc(1, sizeof(int));
	*EndFlag = 0;
	*BE = 0;
	*EN = 0;
	*XOverLength = 0;
	*NCommon = 0;

	x =-1;
    while (FindCycle < 4){
        
		NextPosX = 1;
        
        
		while (FindCycle < 4) {
            
			ProbabilityXOver = 0.0;
            
            
			x = FindNextP(UBXOHN, NextPosX, HighHomol, MedHomol, LowHomol, LenXoverSeq, XoverWindow, XOverHomologyNum);
                
                
           
            if (x > -1 && x != OldX){
				OldX = x;
                
				if (CircularFlag == 1 && XOverHomologyNum[x + (MedHomol - 1)*(UBXOHN + 1)] > XOverHomologyNum[x + (HighHomol - 1)*(UBXOHN + 1)] && x == 1)

					x = FindFirstCOP(x, MedHomol, HighHomol, LenXoverSeq, UBXOHN, XOverHomologyNum);
                    
                else{
                    
					NCommon[0] = 0;
					XOverLength[0] = 0;
                    //3246,0,1,2,3,1,322(targetx),1,x,15,3216,121,0,2,
					x = DefineEventP2(UBXOHN, ShortOutFlag, 1, MedHomol, HighHomol, LowHomol, TargetX, CircularFlag, x, XoverWindow, LenStrainseq0, LenXoverSeq, SeqMinorP, SeqDaughter, EndFlag, BE, EN, NCommon, XOverLength, XoverSeqNumW, XOverHomologyNum);
                    
                    if (XOverLength[0] > 2 && *EN != *BE && (*EN > *BE || CircularFlag == 1) ){

						NumDifferent = XOverLength[0] - NCommon[0];
                        
                        if (NCommon[0] > NumDifferent * 0.8){
                    
                        
							if (XOverLength[0] >= 170) {
								AFact =(double)(((double)(XOverLength[0])) / ((double)(169)));
								NumDifferent = round(NumDifferent * 169 / XOverLength[0]);
								XOverLength[0] = 169;
								NumSame = XOverLength[0] - NumDifferent;
							}
							else {

								AFact = (double)(1);
								NumSame = NCommon[0];
							}
                        
                        
							IndProb = AvHomol[MedHomol];
							GoOn = 0;
							if (ProbEstimateInFileFlag == 0) {
								if (ProbEstimate[XOverLength[0] + NumSame*(UBPE1+1) + ((int)(IndProb * 50))*(UBPE1+1)*(UBPE2+1)] < LowestProb)
									GoOn = 1;
							}
							else
								GoOn = 1;
                            
                            if (GoOn == 1){
                                
								if (XOverLength[0] <= UBFact3x3)

									ProbabilityXOver = ProbCalcP2(Fact3X3, UBFact3x3, XOverLength[0], NumSame, IndProb, LenXoverSeq);
								else

									ProbabilityXOver = ProbCalcP(Fact, XOverLength[0], XOverLength[0] - NumDifferent, IndProb, LenXoverSeq);
                                
                           
                        
                        
                            
                                if (AFact > 1.000001){
									if (ProbabilityXOver > 0)
										ProbabilityXOver = pow(ProbabilityXOver, AFact);
									else
										ProbabilityXOver = 0.05;
                                   
                                }
								
								if (ProbabilityXOver < pow(10, -300))
									ProbabilityXOver = pow(10, -300);


								*BQPV = ProbabilityXOver;

								if (ProbabilityXOver < UCThresh && ProbabilityXOver > 0) {

									if (EarlyBale == 1) {
										
										CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
										free(AH);
										free(NCommon);
										free(XOverLength);
										free(BE);
										free(EN);
										free(EndFlag);
										free(AvHomol);
										return(1);
									}

									if (MCFlag == 0)
										ProbabilityXOver = ProbabilityXOver * (double)(MCCorrection);



									if (ProbabilityXOver < LowestProb && ProbabilityXOver > 0) {
										*BQPV = ProbabilityXOver;
										CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
										free(AH);
										free(NCommon);
										free(XOverLength);
										free(BE);
										free(EN);
										free(EndFlag);
										free(AvHomol);
										return(1);

									}

								}
                           
                            }
                            
                            
                        }
                        
                    
                    }
                    
                }
                
               
				if (EndFlag[0] == 1) {
					EndFlag[0] = 0;
					x = LenXoverSeq;
				}
                
				if (x < LenXoverSeq + 1 && x > NextPosX)
					NextPosX = x + 1;

				else
					break;
                
                
                
            }    
            else
				break;
        }
        
        
        
        
		if (FindCycle == 0) {
			Temp = MedHomol;
			MedHomol = LowHomol;
			LowHomol = Temp;
		}
        else if (FindCycle == 1){
            
			tt = 0.7;
            if (AvHomol[HighHomol] < tt || DoAllFlag == 1){
				Temp = HighHomol;
				HighHomol = LowHomol;
				LowHomol = MedHomol;
				MedHomol = Temp;
			}
            else
                break;
            
		}
        else
            break;
        

		if (HighHomol == 1 && MedHomol == 2 && LowHomol == 3) {
			SeqDaughter = 0;
			SeqMinorP = 2;
		}
		else if (HighHomol == 1 && MedHomol == 3 && LowHomol == 2) {

			SeqDaughter = 1;
			SeqMinorP = 2;
		}
		else if(HighHomol == 2 && MedHomol == 1 && LowHomol == 3) {

			SeqDaughter = 0;
			SeqMinorP = 1;
		}
		else if (HighHomol == 2 && MedHomol == 3 && LowHomol == 1) {

			SeqDaughter = 2; 
			SeqMinorP = 1;
		}
		else if (HighHomol == 3 && MedHomol == 1 && LowHomol == 2) {

			SeqDaughter = 1; 
			SeqMinorP = 0;
		}
		else if (HighHomol == 3 && MedHomol == 2 && LowHomol == 1) {

			SeqDaughter = 2; 
			SeqMinorP = 0;
		}

        FindCycle++;
        
        
    }
	CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
    
	free(AH);
	free(NCommon);
	free(XOverLength);
	free(BE);
	free(EN);
	free(EndFlag);
	free(AvHomol);
	return(0);
}

int MyMathFuncs::FastRecCheckPB(int CircularFlag, int DoAllFlag, int MCCorrection, int MCFlag, int EarlyBale, double UCThresh, double LowestProb, int NextNo, int TargetX, int LenStrainseq0, int ShortOutFlag, int UBD, float *Distance, int UBTD, float *TreeDistance, int UBFSSRDP, int UBCS, int UBXSNW, unsigned char *CS, short int *Seqnum, int Seq1, int Seq2, int Seq3, int LenStrainSeq, int XoverWindow, short int XOverWindowX, short int *XoverSeqNum, char *XoverSeqNumW, int UBXOHN, int *XOverHomologyNum, unsigned char *FSSRDP, int ProbEstimateInFileFlag, int UBPE1, int UBPE2, double *ProbEstimate, int UBFact3x3, double *Fact3X3, double *Fact, double *BQPV) {

	//Dim GoOn As Long, NumSame As Long
	//Dim XOverLen As Long, A As Long, Dummy As Variant, OldX As Long

	//Dim AFact As Double, B As Long, BE As Long, EN As Long, NCommon As Long, StartPosX As Long
	//Dim NextPosX As Long
	//Dim oProbXOver As Double
	//Dim FindCycle As Integer, SeqDaughter As Long, SeqMinorP As Long, Temp As Integer, EndFlag As Long, HighHomol As Long, MedHomol As Long, LowHomol As Long
	//Dim x As Long, NumDifferent As Long

	//Dim SLen As Long
	//Dim AH(2) As Long

	int x, Dummy, tds12, tds23, tds13, Temp;
	int HighHomol, MedHomol, LowHomol, ActiveSeq, ActiveMajorP, ActiveMinorP, SeqDaughter, SeqMinorP;
	int StartPosX, FindCycle, OldX, NextPosX, NumDifferent, NumSame, GoOn;// , test1, test2, test3;
	double AFact, ProbabilityXOver, IndProb, ah1, ah2, ah3;
	int *AH, *NCommon, *XOverLength, *BE, *EN, *EndFlag;
	float tt;
	double *AvHomol;
	AH = (int*)calloc(4, sizeof(int));

	IndProb = 0.0;
	/*test1 = CS[1];
	test2 = CS[500];
	test3 = CS[UBCS + 2];*/



	int LenXoverSeq;

	LenXoverSeq = FindSubSeqPB3(&AH[0], UBFSSRDP, XoverWindow, UBCS, LenStrainseq0, NextNo, Seq1, Seq2, Seq3, &CS[0], UBXSNW, &XoverSeqNumW[0], &FSSRDP[0]);
	//return(LenXoverSeq);

	if (LenXoverSeq < XoverWindow * 2) {
		free(AH);

		CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);

		return(0);
	}

	if (AH[0] < XoverWindow / 3 || AH[1] < XoverWindow / 3 || AH[2] < XoverWindow / 3) {
		free(AH);

		CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
		return(0);
	}
	if (AH[0] == 0 || AH[1] == 0 || AH[2] == 0) {
		CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
		free(AH);

		return(0);
	}


	AvHomol = (double*)calloc(4, sizeof(double));

	AvHomol[1] = ((double)(AH[0])) / ((double)(LenXoverSeq));
	AvHomol[2] = ((double)(AH[1])) / ((double)(LenXoverSeq));
	AvHomol[3] = ((double)(AH[2])) / ((double)(LenXoverSeq));
	ah1 = AvHomol[1];
	ah2 = AvHomol[2];
	ah3 = AvHomol[3];
	//AvHomol(3) = (CLng(AvHomol(3) * 10000)) / 10000 
	//AvHomol(1) = (CLng(AvHomol(1) * 10000)) / 10000 
	//AvHomol(2) = (CLng(AvHomol(2) * 10000)) / 10000 

	/*If MCFlag = 0 Then
	MCC = MCCorrection
	Else
	MCC = 1
	End If*/

	LenXoverSeq = abs(LenXoverSeq) - 1;



	//Work out identities (7.1/21)


	if (UBTD > 0) {
		tds12 = Seq1 + Seq2*(UBTD + 1);
		tds13 = Seq1 + Seq3*(UBTD + 1);
		tds23 = Seq2 + Seq3*(UBTD + 1);
		if (TreeDistance[tds12] >= TreeDistance[tds13] && TreeDistance[tds12] >= TreeDistance[tds23])
			HighHomol = 1;
		else if (TreeDistance[tds13] >= TreeDistance[tds12] && TreeDistance[tds13] >= TreeDistance[tds23])
			HighHomol = 2;
		else if (TreeDistance[tds23] >= TreeDistance[tds12] && TreeDistance[tds23] >= TreeDistance[tds13])
			HighHomol = 3;

	}
	else {
		tds12 = Seq1 + Seq2*(UBD + 1);
		tds13 = Seq1 + Seq3*(UBD + 1);
		tds23 = Seq2 + Seq3*(UBD + 1);
		if (Distance[tds12] >= Distance[tds13] && Distance[tds12] >= Distance[tds23])
			HighHomol = 1;
		else if (Distance[tds13] >= Distance[tds12] && Distance[tds13] >= Distance[tds23])
			HighHomol = 2;
		else if (Distance[tds23] >= Distance[tds12] && Distance[tds23] >= Distance[tds13])
			HighHomol = 3;
	}

	StartPosX = XOHomologyP2(HighHomol, LenStrainSeq, LenXoverSeq, XoverWindow, XoverSeqNumW, XOverHomologyNum);


	if (StartPosX == 0 && ShortOutFlag != 3) {
		CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
		free(AH);
		free(AvHomol);
		return(0);
	}

	if (UBD > 0) {
		if (Seq1 <= UBD && Seq2 <= UBD && Seq3 <= UBD) {
			tds12 = Seq1 + Seq2*(UBD + 1);
			tds13 = Seq1 + Seq3*(UBD + 1);
			tds23 = Seq2 + Seq3*(UBD + 1);
			if (AvHomol[1] == AvHomol[2] && AvHomol[1] == AvHomol[3]) {
				if (Distance[tds12] >= Distance[tds13] && Distance[tds12] >= Distance[tds23]) {
					if (Distance[tds13] > Distance[tds23]) {
						AvHomol[2] = AvHomol[2] - 0.00001;
						AvHomol[3] = AvHomol[3] - 0.00002;
					}
					else {
						AvHomol[2] = AvHomol[2] - 0.00002;
						AvHomol[3] = AvHomol[3] - 0.00001;
					}
				}
				else if (Distance[tds13] >= Distance[tds12] && Distance[tds13] >= Distance[tds23]) {
					if (Distance[tds12] > Distance[tds23]) {
						AvHomol[1] = AvHomol[1] - 0.00001;
						AvHomol[3] = AvHomol[3] - 0.00002;
					}
					else {
						AvHomol[1] = AvHomol[1] - 0.00002;
						AvHomol[3] = AvHomol[3] - 0.00001;
					}
				}
				else {
					if (Distance[tds12] > Distance[tds13]) {
						AvHomol[1] = AvHomol[1] - 0.00001;
						AvHomol[2] = AvHomol[2] - 0.00002;
					}
					else {
						AvHomol[1] = AvHomol[1] - 0.00002;
						AvHomol[2] = AvHomol[2] - 0.00001;
					}
				}
			}
			else if (AvHomol[1] == AvHomol[2]) {
				if (Distance[tds12] > Distance[tds13])
					AvHomol[2] = AvHomol[2] - 0.00001;
				else
					AvHomol[1] = AvHomol[1] - 0.00001;
			}
			else if (AvHomol[1] == AvHomol[3]) {
				if (Distance[tds12] > Distance[tds23])
					AvHomol[3] = AvHomol[3] - 0.00001;
				else
					AvHomol[1] = AvHomol[1] - 0.00001;

			}
			else if (AvHomol[2] == AvHomol[3]) {
				if (Distance[tds13] > Distance[tds23])
					AvHomol[3] = AvHomol[3] - 0.00001;
				else
					AvHomol[2] = AvHomol[2] - 0.00001;

			}

		}
	}



	if (AvHomol[1] >= AvHomol[2] && AvHomol[1] >= AvHomol[3]) {
		HighHomol = 1;

		if (AvHomol[2] >= AvHomol[3]) {
			MedHomol = 2;
			LowHomol = 3;
			ActiveSeq = Seq1;
			ActiveMajorP = Seq2;
			ActiveMinorP = Seq3;
			SeqDaughter = 0;
			SeqMinorP = 2;
		}

		else {
			MedHomol = 3;
			LowHomol = 2;
			ActiveSeq = Seq2;
			ActiveMajorP = Seq1;
			ActiveMinorP = Seq3;
			SeqDaughter = 1;
			SeqMinorP = 2;

		}
	}
	else if (AvHomol[2] >= AvHomol[1] && AvHomol[2] >= AvHomol[3]) {
		HighHomol = 2;

		if (AvHomol[1] >= AvHomol[3]) {
			MedHomol = 1;
			LowHomol = 3;
			ActiveSeq = Seq1;
			ActiveMajorP = Seq3;
			ActiveMinorP = Seq2;
			SeqDaughter = 0;
			SeqMinorP = 1;
		}
		else {
			MedHomol = 3;
			LowHomol = 1;
			ActiveSeq = Seq3;
			ActiveMajorP = Seq1;
			ActiveMinorP = Seq2;
			SeqDaughter = 2;
			SeqMinorP = 1;

		}
	}
	else if (AvHomol[3] >= AvHomol[1] && AvHomol[3] >= AvHomol[2]) {
		HighHomol = 3;

		if (AvHomol[1] >= AvHomol[2]) {
			MedHomol = 1;
			LowHomol = 2;
			ActiveSeq = Seq2;
			ActiveMajorP = Seq3;
			ActiveMinorP = Seq1;
			SeqDaughter = 1;
			SeqMinorP = 0;
		}
		else {
			MedHomol = 2;
			LowHomol = 1;
			ActiveSeq = Seq3;
			ActiveMajorP = Seq2;
			ActiveMinorP = Seq1;
			SeqDaughter = 2;
			SeqMinorP = 0;

		}

	}

	FindCycle = 0;


	OldX = -1;

	NCommon = (int*)calloc(1, sizeof(int));
	XOverLength = (int*)calloc(1, sizeof(int));
	BE = (int*)calloc(1, sizeof(int));
	EN = (int*)calloc(1, sizeof(int));
	EndFlag = (int*)calloc(1, sizeof(int));
	*EndFlag = 0;
	*BE = 0;
	*EN = 0;
	*XOverLength = 0;
	*NCommon = 0;

	x = -1;
	while (FindCycle < 4) {

		NextPosX = 1;


		while (FindCycle < 4) {

			ProbabilityXOver = 0.0;


			x = FindNextP(UBXOHN, NextPosX, HighHomol, MedHomol, LowHomol, LenXoverSeq, XoverWindow, XOverHomologyNum);



			if (x > -1 && x != OldX) {
				OldX = x;

				if (CircularFlag == 1 && XOverHomologyNum[x + (MedHomol - 1)*(UBXOHN + 1)] > XOverHomologyNum[x + (HighHomol - 1)*(UBXOHN + 1)] && x == 1)

					x = FindFirstCOP(x, MedHomol, HighHomol, LenXoverSeq, UBXOHN, XOverHomologyNum);

				else {

					NCommon[0] = 0;
					XOverLength[0] = 0;
					//3246,0,1,2,3,1,322(targetx),1,x,15,3216,121,0,2,
					x = DefineEventP2(UBXOHN, ShortOutFlag, 1, MedHomol, HighHomol, LowHomol, TargetX, CircularFlag, x, XoverWindow, LenStrainseq0, LenXoverSeq, SeqMinorP, SeqDaughter, EndFlag, BE, EN, NCommon, XOverLength, XoverSeqNumW, XOverHomologyNum);

					if (XOverLength[0] > 2 && *EN != *BE && (*EN > *BE || CircularFlag == 1)) {

						NumDifferent = XOverLength[0] - NCommon[0];

						if (NCommon[0] > NumDifferent * 0.8) {


							if (XOverLength[0] >= 170) {
								AFact = (double)(((double)(XOverLength[0])) / ((double)(169)));
								NumDifferent = round(NumDifferent * 169 / XOverLength[0]);
								XOverLength[0] = 169;
								NumSame = XOverLength[0] - NumDifferent;
							}
							else {

								AFact = (double)(1);
								NumSame = NCommon[0];
							}


							IndProb = AvHomol[MedHomol];
							GoOn = 0;
							if (ProbEstimateInFileFlag == 0) {
								if (ProbEstimate[XOverLength[0] + NumSame*(UBPE1 + 1) + ((int)(IndProb * 50))*(UBPE1 + 1)*(UBPE2 + 1)] < LowestProb)
									GoOn = 1;
							}
							else
								GoOn = 1;

							if (GoOn == 1) {

								if (XOverLength[0] <= UBFact3x3)

									ProbabilityXOver = ProbCalcP2(Fact3X3, UBFact3x3, XOverLength[0], NumSame, IndProb, LenXoverSeq);
								else

									ProbabilityXOver = ProbCalcP(Fact, XOverLength[0], XOverLength[0] - NumDifferent, IndProb, LenXoverSeq);





								if (AFact > 1.000001) {
									if (ProbabilityXOver > 0)
										ProbabilityXOver = pow(ProbabilityXOver, AFact);
									else
										ProbabilityXOver = 0.05;

								}

								if (ProbabilityXOver < pow(10, -300))
									ProbabilityXOver = pow(10, -300);

								/*if (*BQPV > ProbabilityXOver && ProbabilityXOver>0)
									*BQPV = ProbabilityXOver;*/

								if (ProbabilityXOver < UCThresh && ProbabilityXOver > 0) {

									/*if (EarlyBale == 1) {

										CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
										free(AH);
										free(NCommon);
										free(XOverLength);
										free(BE);
										free(EN);
										free(EndFlag);
										free(AvHomol);
										return(1);
									}*/

									if (MCFlag == 0)
										ProbabilityXOver = ProbabilityXOver * (double)(MCCorrection);



									if (ProbabilityXOver < LowestProb && ProbabilityXOver > 0) {
										if (*BQPV > ProbabilityXOver)
											*BQPV = ProbabilityXOver;
										/*CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);
										free(AH);
										free(NCommon);
										free(XOverLength);
										free(BE);
										free(EN);
										free(EndFlag);
										free(AvHomol);
										return(1);*/

									}

								}

							}


						}


					}

				}


				if (EndFlag[0] == 1) {
					EndFlag[0] = 0;
					x = LenXoverSeq;
				}

				if (x < LenXoverSeq + 1 && x > NextPosX)
					NextPosX = x + 1;

				else
					break;



			}
			else
				break;
		}




		if (FindCycle == 0) {
			Temp = MedHomol;
			MedHomol = LowHomol;
			LowHomol = Temp;
		}
		else if (FindCycle == 1) {

			tt = 0.7;
			if (AvHomol[HighHomol] < tt || DoAllFlag == 1) {
				Temp = HighHomol;
				HighHomol = LowHomol;
				LowHomol = MedHomol;
				MedHomol = Temp;
			}
			else
				break;

		}
		else
			break;


		if (HighHomol == 1 && MedHomol == 2 && LowHomol == 3) {
			SeqDaughter = 0;
			SeqMinorP = 2;
		}
		else if (HighHomol == 1 && MedHomol == 3 && LowHomol == 2) {

			SeqDaughter = 1;
			SeqMinorP = 2;
		}
		else if (HighHomol == 2 && MedHomol == 1 && LowHomol == 3) {

			SeqDaughter = 0;
			SeqMinorP = 1;
		}
		else if (HighHomol == 2 && MedHomol == 3 && LowHomol == 1) {

			SeqDaughter = 2;
			SeqMinorP = 1;
		}
		else if (HighHomol == 3 && MedHomol == 1 && LowHomol == 2) {

			SeqDaughter = 1;
			SeqMinorP = 0;
		}
		else if (HighHomol == 3 && MedHomol == 2 && LowHomol == 1) {

			SeqDaughter = 2;
			SeqMinorP = 0;
		}

		FindCycle++;


	}
	CleanXOSNW(LenXoverSeq + XoverWindow * 2, XoverWindow, UBXSNW, XoverSeqNumW);

	free(AH);
	free(NCommon);
	free(XOverLength);
	free(BE);
	free(EN);
	free(EndFlag);
	free(AvHomol);

	if (*BQPV < LowestProb && *BQPV > 0) 
		return(1);
	else
		return(0);
}

	double MyMathFuncs::GCCalcPValP(int lseq, long LXover, long *FragMaxScore, double *PVals, long *FragCount, double *KMax, double *LL, int *highenough, double *critval) {
		int X, Y, os, os2;
		double MaxScore, LenXoverSeq, THld;
		float KAScore, LKLen, warn;
		os = lseq + 1;
		LenXoverSeq = (double)(LXover);
		MaxScore = 10;
		for (X = 0; X <= 5; X++) {
			if (highenough[X] == 1) {

				//10^2=100,log100 = 2; exp (10) = e^10 
				//only calculating scores over a critical maximum will massively speed this up
				if (KMax[X] > 0) {

					//1-exp(-exp(-(ll(x)*score-log(kmax*lenxoverseq))))
					// work out which score corresponds with a particular p val
					//pval = 1-exp(-exp(-(ll(x)*score-log(kmax*lenxoverseq))))
					//1-pval = exp(-exp(-(ll(x)*score-log(kmax*lenxoverseq))))
					//-log(1-pval) = exp(-(ll(x)*score-log(kmax*lenxoverseq)))
					//-log(-log(1-pval)) = ll(x)*score-log(kmax*lenxoverseq)
					//log(kmax*lenxoverseq)-log(-log(1-pval))=ll(x)*score
					//(log(kmax*lenxoverseq)-log(-log(1-pval)))/ll(x) = score
					LKLen = (float)(log(KMax[X] * LenXoverSeq));
					for (Y = 0; Y <= FragCount[X]; Y++) {
						os2 = Y + X*os;
						if (FragMaxScore[os2] > critval[X]) {
							KAScore = (float)((LL[X] * FragMaxScore[os2]) - LKLen);
							if (KAScore > 0) {
								if (KAScore < 32) {

									THld = exp((double)(-KAScore));
									PVals[os2] = 1 - exp(-THld);
								}
								else {
									warn = 0;
									if (KAScore > 700) {

										warn = KAScore;
										KAScore = 701;

									}

									THld = exp((double)(-KAScore));
									if (warn != 0) {
										KAScore = (float)(warn - 700);
										THld = THld / (double)(KAScore);
									}
									PVals[os2] = THld;
								}
							}
							else
								PVals[os2] = 1;//THld;


							if (PVals[os2] < MaxScore)
								MaxScore = PVals[os2];
						}
						else
							PVals[os2] = 1;
					}
				}
				else
					highenough[X] = 0;
			}
		}
		return (MaxScore);
	}


	int MyMathFuncs::SignalCountC(int Nextno, int UBXO1, int UBXO2, int AddNum, double LowestProb, short int *CurrentXOver, XOVERDEFINE *XOverlist, int *oRecombNo) {
		int X, Y, Pr,os;
		double Prob;
		os = UBXO1 + 1;

		for (X = 0; X <= Nextno; X++) {
			for (Y = 1; Y <= CurrentXOver[X]; Y++) {
				if (Y <= UBXO2) {
					Prob = XOverlist[X + Y*os].Probability;
					if (Prob < LowestProb && Prob > 0) 
						oRecombNo[XOverlist[X + Y*os].ProgramFlag]++;
					
				}
				else
					CurrentXOver[X]--;
				
			}
		}
		for (X = 0; X < AddNum; X++)
			oRecombNo[100] = oRecombNo[100] + oRecombNo[X];
			
		return(1);
	}
int MyMathFuncs::CMaxD2P(int incnum, int Seq1, int Seq2, int Seq3, int SBP, int EBP, int Nextno, int SLen, short int *SeqNum, short int *SeqnumX, int *IdenticalR, int *IdenticalF, unsigned char *NucMat, int *IncSeq3, unsigned char *IncSeq2, unsigned char *IncSeq, float *E, float *d, float *VScoreMat, float *AvDist, int *TotCount) {
	int  Se1, Se2, Se3, Se4, osx, X, Y, SBPM, EBPP;
	unsigned int v1, ofe1, ofe2, S1, S2, S3, SeX1, SeX2, SeX3, SeX4;
	//int A, B, C, D;
	char GoOn, go, go1;
	float e0, e1, e2, d0, d1, d2, FS, Dist1, Dist2, Dist3, Dist4, Dist, d3;


	if (SBP > 0)
		SBPM = IdenticalF[SBP - 1];
	else
		SBPM = IdenticalF[SBP];


	if (EBP < SLen)
		EBPP = IdenticalF[EBP + 1];
	else
		EBPP = IdenticalF[EBP];


	SBP = IdenticalF[SBP];
	EBP = IdenticalF[EBP];
	d3 = (float)(1 / 3);
	ofe2 = SLen + 1;
	for (X = 0; X <= Nextno; X++) {
		ofe1 = X*ofe2;
		for (Y = 0; Y <= SLen; Y++) {
			SeqnumX[Y + ofe1] = NucMat[SeqNum[Y + ofe1]];
		}
	}
	//VScoreMat(4, 4, 4, 4, 2)
	//	for (A = 0; A <= incnum-3; A++){

	int procs;
	procs = omp_get_num_procs();
	procs = procs/2 - 1;
	if (procs < 3)
		procs = 3;
	omp_set_num_threads(procs);

	for (Se1 = 0; Se1 <= Nextno - 3; Se1++) {
		//		Se1 = IncSeq3[A];
		if (IncSeq2[Se1] == 1) {
			//			for (B = A+1; B <= incnum-2; B++){
			SeX1 = Se1*ofe2;
			for (Se2 = Se1 + 1; Se2 <= Nextno - 2; Se2++) {
				//				Se2 = IncSeq3[B];
				if (IncSeq2[Se2] == 1) {
					//					for (C = B+1; C <= incnum-1; C++){
					go1 = IncSeq[Se1] + IncSeq[Se2];
					SeX2 = Se2*ofe2;
					for (Se3 = Se2 + 1; Se3 <= Nextno - 1; Se3++) {
						//						Se3 = IncSeq3[C];
						if (IncSeq2[Se3] == 1) {
							//							for (D = C+1; D <= incnum; D++){
							go = go1 + IncSeq[Se3];
							SeX3 = Se3*ofe2;

#pragma omp parallel for private (Se4, SeX4, X, S1, S2, S3, osx, v1, e0, e1, e2, d0, d1, d2, FS, GoOn, Dist, Dist1, Dist2, Dist3, Dist4)
							for (Se4 = Se3 + 1; Se4 <= Nextno; Se4++) {
								//								Se4 = IncSeq3[D];
								if (IncSeq2[Se4] == 1) {
									GoOn = go + IncSeq[Se4];

									if (GoOn > 0) {// Then 'Seq1 = Se1 Or Seq2 = Se1 Or Seq3 = Se1 Or Seq1 = Se2 Or Seq2 = Se2 Or Seq3 = Se2 Or Seq1 = Se3 Or Seq2 = Se3 Or Seq3 = Se3 Or Seq1 = Se4 Or Seq2 = Se4 Or Seq3 = Se4 Then

										e0 = 0;
										e1 = 0;
										e2 = 0;
										d0 = 0;
										d1 = 0;
										d2 = 0;


										SeX4 = Se4*ofe2;
										if (SBP < EBP) {

											
										//{

												//#pragma omp section
												//{
													
												for (X = 1; X <= SBPM; X++) {
													//if (Identical[X] == 0){
													osx = IdenticalR[X];
													S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
													S2 = SeqnumX[osx + SeX2];
													S3 = SeqnumX[osx + SeX3];


													v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

													e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													e2 = e2 + VScoreMat[v1 + 1250];
													/*
													if (S1 != S2 ){
													//S4 = SeqnumX[osx + SeX4];


													//maybe use a lookup table for these
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

													e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													e2 = e2 + VScoreMat[v1 + 1250];

													}
													else if (S1 != S3){

													//	S4 = SeqnumX[osx + SeX4];



													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

													e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													e2 = e2 + VScoreMat[v1 + 1250];

													}
													*/
													//}
												}
										
												//}
												//#pragma omp section
												//{
												
												
												for (X = EBPP; X <= IdenticalF[SLen]; X++) {
													//if (Identical[X] == 0){
													osx = IdenticalR[X];
													S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
													S2 = SeqnumX[osx + SeX2];
													S3 = SeqnumX[osx + SeX3];


													v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

													e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													e2 = e2 + VScoreMat[v1 + 1250];

													/*	if (S1 != S2 ){
													//S4 = SeqnumX[osx + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

													e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													e2 = e2 + VScoreMat[v1 + 1250];

													}
													else if ( S1 != S3){
													//S4 = SeqnumX[osx + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

													e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													e2 = e2 + VScoreMat[v1 + 1250];

													}
													*/
													//}
												}
												//}
											//}
											
											FS = e0 + e1 + e2;
											if (FS > 0) {
												e0 = e0 / FS;
												e1 = e1 / FS;
												e2 = e2 / FS;




												for (X = SBP; X <= EBP; X++) {
													//if (Identical[X] == 0){
													osx = IdenticalR[X];
													S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
													S2 = SeqnumX[osx + SeX2];
													S3 = SeqnumX[osx + SeX3];
													v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];
													/*	if (S1 != S2 ){
													//S4 = SeqnumX[osx + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}
													else if ( S1 != S3){
													//S4 = SeqnumX[osx + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}
													//}*/
												}

												FS = d0 + d1 + d2;
												if (FS > 0) {
													d0 = d0 / FS;
													d1 = d1 / FS;
													d2 = d2 / FS;
												}
												else {

													d0 = d3;
													d1 = d3;
													d2 = d3;
													e0 = d3;
													e1 = d3;
													e2 = d3;
												}
											}
											else {

												e0 = d3;
												e1 = d3;
												e2 = d3;
												d0 = d3;
												d1 = d3;
												d2 = d3;
											}
										}
										else {


											for (X = EBPP; X <= SBPM; X++) {
												//if (Identical[X] == 0){
												osx = IdenticalR[X];
												S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
												S2 = SeqnumX[osx + SeX2];
												S3 = SeqnumX[osx + SeX3];
												v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

												e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												e2 = e2 + VScoreMat[v1 + 1250];

												/*
												if (S1 != S2 || S1 != S3){
												//S4 = SeqnumX[osx + SeX4];
												v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

												e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												e2 = e2 + VScoreMat[v1 + 1250];

												}*/
												//}
											}
											FS = e0 + e1 + e2;
											if (FS > 0) {
												e0 = e0 / FS;
												e1 = e1 / FS;
												e2 = e2 / FS;



												for (X = 1; X <= EBP; X++) {
													//if (Identical[X] == 0){
													osx = IdenticalR[X];
													S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
													S2 = SeqnumX[osx + SeX2];
													S3 = SeqnumX[osx + SeX3];
													v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;
													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];
													/*if (S1 != S2 ){
													//	S4 = SeqnumX[osx  + SeX4];


													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}
													else if (S1 != S3){
													//S4 = SeqnumX[osx  + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx  + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}*/
													//}
												}
												for (X = SBP; X <= IdenticalF[SLen]; X++) {
													//if (Identical[X] == 0){
													osx = IdenticalR[X];
													S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
													S2 = SeqnumX[osx + SeX2];
													S3 = SeqnumX[osx + SeX3];
													v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];
													/*	if (S1 != S2 ){
													//S4 = SeqnumX[osx  + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx  + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}
													else if ( S1 != S3){
													//S4 = SeqnumX[osx  + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx  + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}*/
													//}
												}
												FS = d0 + d1 + d2;
												if (FS > 0) {
													d0 = d0 / FS;
													d1 = d1 / FS;
													d2 = d2 / FS;
												}
												else {

													d0 = d3;
													d1 = d3;
													d2 = d3;
													e0 = d3;
													e1 = d3;
													e2 = d3;
												}
											}
											else {

												e0 = d3;
												e1 = d3;
												e2 = d3;
												d0 = d3;
												d1 = d3;
												d2 = d3;
											}
										}
										//it doesn't matter what the actual distance is - the relative distance is what matters. or does it?
										if (d0 != d3 || d1 != d3) {
											Dist1 = (float)(fabs(d0 - e0));
											//Dist1 = Dist1*Dist1;// ^ 2
											Dist2 = (float)(fabs(d1 - e1));
											//Dist2 = Dist2*Dist2;
											Dist3 = (float)(fabs(d2 - e2));
											//Dist3 = Dist3*Dist3;
											Dist4 = Dist1 + Dist2 + Dist3;
											Dist = (float)Dist4;//(pow(Dist4,0.5));
							#pragma omp critical
											{
												if (Se1 == Seq1 || Se2 == Seq1 || Se3 == Seq1 || Se4 == Seq1) {
//#pragma omp atomic
													AvDist[0] += Dist;
//#pragma omp atomic
													TotCount[0] ++;

												}
												if (Se1 == Seq2 || Se2 == Seq2 || Se3 == Seq2 || Se4 == Seq2) {
//#pragma omp atomic
													AvDist[1] += Dist;
//#pragma omp atomic
													TotCount[1]++;

												}
												if (Se1 == Seq3 || Se2 == Seq3 || Se3 == Seq3 || Se4 == Seq3) {
//#pragma omp atomic
													AvDist[2] += Dist;
//#pragma omp atomic
													TotCount[2]++;

												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}


	omp_set_num_threads(2);
	return(1);
}


int MyMathFuncs::CMaxD2P3(int incnum, int Seq1, int Seq2, int Seq3, int SBP, int EBP, int Nextno, int SLen, short int *SeqNum, short int *SeqnumX, int *IdenticalR, int *IdenticalF, unsigned char *NucMat, int *IncSeq3, unsigned char *IncSeq2, unsigned char *IncSeq, float *E, float *d, float *VScoreMat, float *AvDist, int *TotCount) {
	int  Se1, Se2, Se3, Se4, osx, X, Y, SBPM, EBPP;
	unsigned int v1, ofe1, ofe2, S1, S2, S3, SeX1, SeX2, SeX3, SeX4;
	//int A, B, C, D;
	char GoOn, go, go1;
	float e0, e1, e2, d0, d1, d2, FS, Dist1, Dist2, Dist3, Dist4, Dist, d3;


	if (SBP > 0)
		SBPM = IdenticalF[SBP - 1];
	else
		SBPM = IdenticalF[SBP];


	if (EBP < SLen)
		EBPP = IdenticalF[EBP + 1];
	else
		EBPP = IdenticalF[EBP];


	SBP = IdenticalF[SBP];
	EBP = IdenticalF[EBP];
	d3 = (float)(1 / 3);
	ofe2 = SLen + 1;
	for (X = 0; X <= Nextno; X++) {
		ofe1 = X*ofe2;
		for (Y = 0; Y <= SLen; Y++) {
			SeqnumX[Y + ofe1] = NucMat[SeqNum[Y + ofe1]];
		}
	}
	//VScoreMat(4, 4, 4, 4, 2)
	//	for (A = 0; A <= incnum-3; A++){

	int procs;
	procs = omp_get_num_procs();
	procs = procs / 2 -1;
	if (procs < 3)
		procs = 3;
	omp_set_num_threads(procs);
	int b, c, dimsize;
	int *xy;
	dimsize = (incnum + 2)*(incnum+1)*(incnum)*(incnum - 1);
	dimsize = dimsize / 3;
	//dimsize = dimsize * 4;
	xy = (int *)calloc(dimsize, sizeof(int));
	c = 0;
	for (Se1 = 0; Se1 <= Nextno - 3; Se1++) {
		if (IncSeq2[Se1] == 1) {
			SeX1 = Se1*ofe2;
			for (Se2 = Se1 + 1; Se2 <= Nextno - 2; Se2++) {
				if (IncSeq2[Se2] == 1) {
					go1 = IncSeq[Se1] + IncSeq[Se2];
					SeX2 = Se2*ofe2;
					for (Se3 = Se2 + 1; Se3 <= Nextno - 1; Se3++) {
						if (IncSeq2[Se3] == 1) {
							go = go1 + IncSeq[Se3];
							SeX3 = Se3*ofe2;
							for (Se4 = Se3 + 1; Se4 <= Nextno; Se4++) {
								if (IncSeq2[Se4] == 1) {
									GoOn = go + IncSeq[Se4];
									if (GoOn > 0) {
										xy[c] = Se1;
										xy[c+1] = Se2;
										xy[c+2] = Se3;
										xy[c+3] = Se4;
										c = c + 4;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	c = c - 4;
	c = c / 4;
#pragma omp parallel for private (b, Se1, Se2, Se3, Se4, SeX1, SeX2, SeX3, SeX4, X, S1, S2, S3, osx, v1, e0, e1, e2, d0, d1, d2, FS, GoOn, Dist, Dist1, Dist2, Dist3, Dist4)
	for (b = 0; b <= c; b++) {
		Se1 = xy[b*4];
		Se2 = xy[1+ b*4];
		Se3 = xy[2 + b*4];
		Se4 = xy[3 +b*4];
		SeX1 = Se1*ofe2;
		SeX2 = Se2*ofe2;
		SeX3 = Se3*ofe2;
		SeX4 = Se4*ofe2;
		e0 = 0;
		e1 = 0;
		e2 = 0;
		d0 = 0;
		d1 = 0;
		d2 = 0;

		if (SBP < EBP) {


			for (X = 1; X <= SBPM; X++) {
				osx = IdenticalR[X];
				S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
				S2 = SeqnumX[osx + SeX2];
				S3 = SeqnumX[osx + SeX3];

				v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

				e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
				e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
				e2 = e2 + VScoreMat[v1 + 1250];
				
			}



			for (X = EBPP; X <= IdenticalF[SLen]; X++) {
				osx = IdenticalR[X];
				S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
				S2 = SeqnumX[osx + SeX2];
				S3 = SeqnumX[osx + SeX3];


				v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

				e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
				e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
				e2 = e2 + VScoreMat[v1 + 1250];

			}
			FS = e0 + e1 + e2;
			if (FS > 0) {
				e0 = e0 / FS;
				e1 = e1 / FS;
				e2 = e2 / FS;

				for (X = SBP; X <= EBP; X++) {
					//if (Identical[X] == 0){
					osx = IdenticalR[X];
					S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
					S2 = SeqnumX[osx + SeX2];
					S3 = SeqnumX[osx + SeX3];
					v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

					d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
					d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
					d2 = d2 + VScoreMat[v1 + 1250];
					
				}

				FS = d0 + d1 + d2;
				if (FS > 0) {
					d0 = d0 / FS;
					d1 = d1 / FS;
					d2 = d2 / FS;
				}
				else {

					d0 = d3;
					d1 = d3;
					d2 = d3;
					e0 = d3;
					e1 = d3;
					e2 = d3;
				}
			}
			else {

				e0 = d3;
				e1 = d3;
				e2 = d3;
				d0 = d3;
				d1 = d3;
				d2 = d3;
			}
		}
		else {


			for (X = EBPP; X <= SBPM; X++) {
				//if (Identical[X] == 0){
				osx = IdenticalR[X];
				S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
				S2 = SeqnumX[osx + SeX2];
				S3 = SeqnumX[osx + SeX3];
				v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

				e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
				e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
				e2 = e2 + VScoreMat[v1 + 1250];

				
			}
			FS = e0 + e1 + e2;
			if (FS > 0) {
				e0 = e0 / FS;
				e1 = e1 / FS;
				e2 = e2 / FS;



				for (X = 1; X <= EBP; X++) {
					//if (Identical[X] == 0){
					osx = IdenticalR[X];
					S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
					S2 = SeqnumX[osx + SeX2];
					S3 = SeqnumX[osx + SeX3];
					v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;
					d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
					d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
					d2 = d2 + VScoreMat[v1 + 1250];
					
				}
				for (X = SBP; X <= IdenticalF[SLen]; X++) {
					//if (Identical[X] == 0){
					osx = IdenticalR[X];
					S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
					S2 = SeqnumX[osx + SeX2];
					S3 = SeqnumX[osx + SeX3];
					v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

					d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
					d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
					d2 = d2 + VScoreMat[v1 + 1250];
					
				}
				FS = d0 + d1 + d2;
				if (FS > 0) {
					d0 = d0 / FS;
					d1 = d1 / FS;
					d2 = d2 / FS;
				}
				else {

					d0 = d3;
					d1 = d3;
					d2 = d3;
					e0 = d3;
					e1 = d3;
					e2 = d3;
				}
			}
			else {

				e0 = d3;
				e1 = d3;
				e2 = d3;
				d0 = d3;
				d1 = d3;
				d2 = d3;
			}
		}
		//it doesn't matter what the actual distance is - the relative distance is what matters. or does it?
		if (d0 != d3 || d1 != d3) {
			Dist1 = (float)(fabs(d0 - e0));
			//Dist1 = Dist1*Dist1;// ^ 2
			Dist2 = (float)(fabs(d1 - e1));
			//Dist2 = Dist2*Dist2;
			Dist3 = (float)(fabs(d2 - e2));
			//Dist3 = Dist3*Dist3;
			Dist4 = Dist1 + Dist2 + Dist3;
			Dist = (float)Dist4;//(pow(Dist4,0.5));
#pragma omp critical
			{
				if (Se1 == Seq1 || Se2 == Seq1 || Se3 == Seq1 || Se4 == Seq1) {
					//#pragma omp atomic
					AvDist[0] += Dist;
					//#pragma omp atomic
					TotCount[0] ++;

				}
				if (Se1 == Seq2 || Se2 == Seq2 || Se3 == Seq2 || Se4 == Seq2) {
					//#pragma omp atomic
					AvDist[1] += Dist;
					//#pragma omp atomic
					TotCount[1]++;

				}
				if (Se1 == Seq3 || Se2 == Seq3 || Se3 == Seq3 || Se4 == Seq3) {
					//#pragma omp atomic
					AvDist[2] += Dist;
					//#pragma omp atomic
					TotCount[2]++;

				}
			}
		}

						
	}


	omp_set_num_threads(2);
	free(xy);
	return(1);
}


int MyMathFuncs::CMaxD2P2(int incnum, int Seq1, int Seq2, int Seq3, int SBP, int EBP, int Nextno, int SLen, short int *SeqNum, short int *SeqnumX, const int *IdenticalR, int *IdenticalF, unsigned char *NucMat, int *IncSeq3, unsigned char *IncSeq2, unsigned char *IncSeq, float *E, float *d, float *VScoreMat, float *AvDist, int *TotCount) {
	int  Se1, Se2, Se3, Se4, osx, X, Y, sbpm, ebpp;
	unsigned int v1, ofe1, ofe2, S1, S2, S3, SeX1, SeX2, SeX3, SeX4;
	//int A, B, C, D;
	char GoOn, go, go1;
	float e0, e1, e2, d0, d1, d2, FS, Dist1, Dist2, Dist3, Dist4, Dist, d3;

	if (SBP > 0)
		sbpm = IdenticalF[SBP - 1];
	else
		sbpm = IdenticalF[SBP];


	if (EBP < SLen)
		ebpp = IdenticalF[EBP + 1];
	else
		ebpp = IdenticalF[EBP];

	const int SBPM = sbpm; 
	const int EBPP = ebpp;

	SBP = IdenticalF[SBP];
	EBP = IdenticalF[EBP];
	d3 = (float)(1 / 3);
	ofe2 = SLen + 1;
	for (X = 0; X <= Nextno; X++) {
		ofe1 = X*ofe2;
		for (Y = 0; Y <= SLen; Y++) {
			SeqnumX[Y + ofe1] = NucMat[SeqNum[Y + ofe1]];
		}
	}
	//VScoreMat(4, 4, 4, 4, 2)
	//	for (A = 0; A <= incnum-3; A++){
	for (Se1 = 0; Se1 <= Nextno - 3; Se1++) {
		//		Se1 = IncSeq3[A];
		if (IncSeq2[Se1] == 1) {
			//			for (B = A+1; B <= incnum-2; B++){
			SeX1 = Se1*ofe2;
			for (Se2 = Se1 + 1; Se2 <= Nextno - 2; Se2++) {
				//				Se2 = IncSeq3[B];
				if (IncSeq2[Se2] == 1) {
					//					for (C = B+1; C <= incnum-1; C++){
					go1 = IncSeq[Se1] + IncSeq[Se2];
					SeX2 = Se2*ofe2;
					for (Se3 = Se2 + 1; Se3 <= Nextno - 1; Se3++) {
						//						Se3 = IncSeq3[C];
						if (IncSeq2[Se3] == 1) {
							//							for (D = C+1; D <= incnum; D++){
							go = go1 + IncSeq[Se3];
							SeX3 = Se3*ofe2;

//#pragma omp parallel for private (Se4, SeX4, X, S1, S2, S3, osx, v1, e0, e1, e2, d0, d1, d2, FS, GoOn, Dist, Dist1, Dist2, Dist3, Dist4)
							for (Se4 = Se3 + 1; Se4 <= Nextno; Se4++) {
								//								Se4 = IncSeq3[D];
								if (IncSeq2[Se4] == 1) {
									GoOn = go + IncSeq[Se4];

									if (GoOn > 0) {// Then 'Seq1 = Se1 Or Seq2 = Se1 Or Seq3 = Se1 Or Seq1 = Se2 Or Seq2 = Se2 Or Seq3 = Se2 Or Seq1 = Se3 Or Seq2 = Se3 Or Seq3 = Se3 Or Seq1 = Se4 Or Seq2 = Se4 Or Seq3 = Se4 Then

										e0 = 0;
										e1 = 0;
										e2 = 0;
										d0 = 0;
										d1 = 0;
										d2 = 0;


										SeX4 = Se4*ofe2;
										if (SBP < EBP) {


											//{

											//#pragma omp section
											//{
//#pragma
//#pragma simd
											for (X = 1; X <= SBPM; X++) {
												//if (Identical[X] == 0){
												//osx = IdenticalR[X];
												S1 = SeqnumX[(int)(IdenticalR[X] + SeX1)]; //'85,72,66,68
												S2 = SeqnumX[(int)(IdenticalR[X] + SeX2)];
												S3 = SeqnumX[(int)(IdenticalR[X] + SeX3)];


												//v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

												//e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												//e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												//e2 = e2 + VScoreMat[v1 + 1250];
												/*
												if (S1 != S2 ){
												//S4 = SeqnumX[osx + SeX4];


												//maybe use a lookup table for these
												v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

												e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												e2 = e2 + VScoreMat[v1 + 1250];

												}
												else if (S1 != S3){

												//	S4 = SeqnumX[osx + SeX4];



												v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

												e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												e2 = e2 + VScoreMat[v1 + 1250];

												}
												*/
												//}
											}

											//}
											//#pragma omp section
											//{

//#pragma simd
											for (X = EBPP; X <= IdenticalF[SLen]; X++) {
												//if (Identical[X] == 0){
												osx = IdenticalR[X];
												S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
												S2 = SeqnumX[osx + SeX2];
												S3 = SeqnumX[osx + SeX3];


												v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

												e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												e2 = e2 + VScoreMat[v1 + 1250];

												/*	if (S1 != S2 ){
												//S4 = SeqnumX[osx + SeX4];
												v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

												e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												e2 = e2 + VScoreMat[v1 + 1250];

												}
												else if ( S1 != S3){
												//S4 = SeqnumX[osx + SeX4];
												v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

												e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												e2 = e2 + VScoreMat[v1 + 1250];

												}
												*/
												//}
											}
											//}
											//}

											FS = e0 + e1 + e2;
											if (FS > 0) {
												e0 = e0 / FS;
												e1 = e1 / FS;
												e2 = e2 / FS;



//#pragma simd
												for (X = SBP; X <= EBP; X++) {
													//if (Identical[X] == 0){
													osx = IdenticalR[X];
													S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
													S2 = SeqnumX[osx + SeX2];
													S3 = SeqnumX[osx + SeX3];
													v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];
													/*	if (S1 != S2 ){
													//S4 = SeqnumX[osx + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}
													else if ( S1 != S3){
													//S4 = SeqnumX[osx + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}
													//}*/
												}

												FS = d0 + d1 + d2;
												if (FS > 0) {
													d0 = d0 / FS;
													d1 = d1 / FS;
													d2 = d2 / FS;
												}
												else {

													d0 = d3;
													d1 = d3;
													d2 = d3;
													e0 = d3;
													e1 = d3;
													e2 = d3;
												}
											}
											else {

												e0 = d3;
												e1 = d3;
												e2 = d3;
												d0 = d3;
												d1 = d3;
												d2 = d3;
											}
										}
										else {

//#pragma simd
											for (X = EBPP; X <= SBPM; X++) {
												//if (Identical[X] == 0){
												osx = IdenticalR[X];
												S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
												S2 = SeqnumX[osx + SeX2];
												S3 = SeqnumX[osx + SeX3];
												v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

												e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												e2 = e2 + VScoreMat[v1 + 1250];

												/*
												if (S1 != S2 || S1 != S3){
												//S4 = SeqnumX[osx + SeX4];
												v1 = S1 + S2*5 + S3*25 + SeqnumX[osx + SeX4]*125;

												e0 = e0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
												e1 = e1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
												e2 = e2 + VScoreMat[v1 + 1250];

												}*/
												//}
											}
											FS = e0 + e1 + e2;
											if (FS > 0) {
												e0 = e0 / FS;
												e1 = e1 / FS;
												e2 = e2 / FS;


//#pragma simd
												for (X = 1; X <= EBP; X++) {
													//if (Identical[X] == 0){
													osx = IdenticalR[X];
													S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
													S2 = SeqnumX[osx + SeX2];
													S3 = SeqnumX[osx + SeX3];
													v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;
													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];
													/*if (S1 != S2 ){
													//	S4 = SeqnumX[osx  + SeX4];


													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}
													else if (S1 != S3){
													//S4 = SeqnumX[osx  + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx  + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}*/
													//}
												}
//#pragma simd
												for (X = SBP; X <= IdenticalF[SLen]; X++) {
													//if (Identical[X] == 0){
													osx = IdenticalR[X];
													S1 = SeqnumX[osx + SeX1]; //'85,72,66,68
													S2 = SeqnumX[osx + SeX2];
													S3 = SeqnumX[osx + SeX3];
													v1 = S1 + S2 * 5 + S3 * 25 + SeqnumX[osx + SeX4] * 125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];
													/*	if (S1 != S2 ){
													//S4 = SeqnumX[osx  + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx  + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}
													else if ( S1 != S3){
													//S4 = SeqnumX[osx  + SeX4];
													v1 = S1 + S2*5 + S3*25 + SeqnumX[osx  + SeX4]*125;

													d0 = d0 + VScoreMat[v1];  //1=2:3=4 = 1;1=2:3<>4 = 0.5;1<>2:3=4 = 0.5
													d1 = d1 + VScoreMat[v1 + 625]; //1=3:2=4 = 1; 1=3:2<>4 = 0.5
													d2 = d2 + VScoreMat[v1 + 1250];

													}*/
													//}
												}
												FS = d0 + d1 + d2;
												if (FS > 0) {
													d0 = d0 / FS;
													d1 = d1 / FS;
													d2 = d2 / FS;
												}
												else {

													d0 = d3;
													d1 = d3;
													d2 = d3;
													e0 = d3;
													e1 = d3;
													e2 = d3;
												}
											}
											else {

												e0 = d3;
												e1 = d3;
												e2 = d3;
												d0 = d3;
												d1 = d3;
												d2 = d3;
											}
										}
										//it doesn't matter what the actual distance is - the relative distance is what matters. or does it?
										if (d0 != d3 || d1 != d3) {
											Dist1 = (float)(fabs(d0 - e0));
											//Dist1 = Dist1*Dist1;// ^ 2
											Dist2 = (float)(fabs(d1 - e1));
											//Dist2 = Dist2*Dist2;
											Dist3 = (float)(fabs(d2 - e2));
											//Dist3 = Dist3*Dist3;
											Dist4 = Dist1 + Dist2 + Dist3;
											Dist = (float)Dist4;//(pow(Dist4,0.5));
#pragma omp critical
											{
												if (Se1 == Seq1 || Se2 == Seq1 || Se3 == Seq1 || Se4 == Seq1) {
													//#pragma omp atomic
													AvDist[0] += Dist;
													//#pragma omp atomic
													TotCount[0] ++;

												}
												if (Se1 == Seq2 || Se2 == Seq2 || Se3 == Seq2 || Se4 == Seq2) {
													//#pragma omp atomic
													AvDist[1] += Dist;
													//#pragma omp atomic
													TotCount[1]++;

												}
												if (Se1 == Seq3 || Se2 == Seq3 || Se3 == Seq3 || Se4 == Seq3) {
													//#pragma omp atomic
													AvDist[2] += Dist;
													//#pragma omp atomic
													TotCount[2]++;

												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}



	return(1);
}

double MyMathFuncs::ViterbiCP(int SLen, int NumberAB, int NumberXY, double *OptXY, unsigned char *RecodeB, double *LaticeXY, double *TransitionM2, double *EmissionM2, double *LaticeAB)
{


	//transition(numberxy-1, numberxy-1)
	//optxy(numberxy-1, numberxy-1)
	//emmision(numberab-1,numberxy-1)
	//laticeab(slen,numberxy-1)
	//laticexy(slen,numberxy-1)
	int A, B, X, off1, off2, off3, offa, offa2;

	off1 = NumberXY;
	off2 = NumberAB;
	off3 = SLen + 1;
	for (X = 1; X <= SLen; X++) {

		for (A = 0; A < NumberXY; A++) { //first state
			offa = X - 1 + A*off3;
			for (B = 0; B < NumberXY; B++) { //second state
				OptXY[A + B*off1] = LaticeXY[offa] + TransitionM2[A + B*off1] + EmissionM2[RecodeB[X] + B*off2];

			}
		}



		for (A = 0; A < NumberXY; A++) { //first state
			offa = X + A*off3;
			offa2 = A*off1;
			LaticeXY[offa] = -10000000000000000;
			for (B = 0; B < NumberXY; B++) { //second state
				if (LaticeXY[offa] < OptXY[B + offa2]) {
					LaticeXY[offa] = OptXY[B + offa2];
					LaticeAB[offa] = B;
				}
			}
		}

	}
	return(1);
}


float MyMathFuncs::ViterbiCPF(int SLen, int NumberAB, int NumberXY, float *OptXY, unsigned char *RecodeB, float *LaticeXY, float *TransitionM2, float *EmissionM2, float *LaticeAB)
{


	//transition(numberxy-1, numberxy-1)
	//optxy(numberxy-1, numberxy-1)
	//emmision(numberab-1,numberxy-1)
	//laticeab(slen,numberxy-1)
	//laticexy(slen,numberxy-1)
	int A, B, X, off1, off2, off3, offa, offa2;

	off1 = NumberXY;
	off2 = NumberAB;
	off3 = SLen + 1;
	for (X = 1; X <= SLen; X++) {

		for (A = 0; A < NumberXY; A++) { //first state
			offa = X - 1 + A*off3;
			for (B = 0; B < NumberXY; B++) { //second state
				OptXY[A + B*off1] = LaticeXY[offa] + TransitionM2[A + B*off1] + EmissionM2[RecodeB[X] + B*off2];

			}
		}



		for (A = 0; A < NumberXY; A++) { //first state
			offa = X + A*off3;
			offa2 = A*off1;
			LaticeXY[offa] = -10000000000000000;
			for (B = 0; B < NumberXY; B++) { //second state
				if (LaticeXY[offa] < OptXY[B + offa2]) {
					LaticeXY[offa] = OptXY[B + offa2];
					LaticeAB[offa] = B;
				}
			}
		}

	}
	return(1);
}


double  MyMathFuncs::ForwardCP(int SLen, int NumberAB, int NumberXY, double *ValXY, double *OptXY, unsigned char *RecodeB, double *LaticeXY, double *TransitionM2, double *EmissionM2)
{


	//transition(numberxy-1, numberxy-1)
	//optxy(numberxy-1, numberxy-1)
	//emmision(numberab-1,numberxy-1)
	//laticeab(slen,numberxy-1)
	//laticexy(slen,numberxy-1)
	int A, B, X, off1, off2, off3;
	double MinV;

	off1 = NumberXY;
	off2 = NumberAB;
	off3 = SLen + 1;

	for (X = 1; X <= SLen; X++) {

		for (A = 0; A < NumberXY; A++) { //first state
			for (B = 0; B < NumberXY; B++) { //second state
				OptXY[A + B*off1] = LaticeXY[X - 1 + A*off3] + TransitionM2[A + B*off1] + EmissionM2[RecodeB[X] + B*off2];

			}
		}


		//ReDim ValXY(NumberXY - 1)
		for (A = 0; A < NumberXY; A++)
			ValXY[A] = 0;

		for (A = 0; A < NumberXY; A++) {
			MinV = -10000000000000000;
			for (B = 0; B < NumberXY; B++) {
				if (MinV < OptXY[B + A*off1])
					MinV = OptXY[B + A*off1];

			}
			LaticeXY[X + A*off3] = 0;
			for (B = 0; B < NumberXY; B++) {
				ValXY[B] = OptXY[B + A*off1] - MinV;
				ValXY[B] = exp(ValXY[B]);
				LaticeXY[X + A*off3] = LaticeXY[X + A*off3] + ValXY[B];

			}
			LaticeXY[X + A*off3] = log(LaticeXY[X + A*off3]);
			LaticeXY[X + A*off3] = LaticeXY[X + A*off3] + MinV;

		}

	}
	return(1);
}

float  MyMathFuncs::ForwardCPF(int SLen, int NumberAB, int NumberXY, float *ValXY, float *OptXY, unsigned char *RecodeB, float *LaticeXY, float *TransitionM2, float *EmissionM2)
{


	//transition(numberxy-1, numberxy-1)
	//optxy(numberxy-1, numberxy-1)
	//emmision(numberab-1,numberxy-1)
	//laticeab(slen,numberxy-1)
	//laticexy(slen,numberxy-1)
	int A, B, X, off1, off2, off3;
	double MinV;

	off1 = NumberXY;
	off2 = NumberAB;
	off3 = SLen + 1;

	for (X = 1; X <= SLen; X++) {

		for (A = 0; A < NumberXY; A++) { //first state
			for (B = 0; B < NumberXY; B++) { //second state
				OptXY[A + B*off1] = LaticeXY[X - 1 + A*off3] + TransitionM2[A + B*off1] + EmissionM2[RecodeB[X] + B*off2];

			}
		}


		//ReDim ValXY(NumberXY - 1)
		for (A = 0; A < NumberXY; A++)
			ValXY[A] = 0;

		for (A = 0; A < NumberXY; A++) {
			MinV = -10000000000000000;
			for (B = 0; B < NumberXY; B++) {
				if (MinV < OptXY[B + A*off1])
					MinV = OptXY[B + A*off1];

			}
			LaticeXY[X + A*off3] = 0;
			for (B = 0; B < NumberXY; B++) {
				ValXY[B] = OptXY[B + A*off1] - MinV;
				ValXY[B] = exp(ValXY[B]);
				LaticeXY[X + A*off3] = LaticeXY[X + A*off3] + ValXY[B];

			}
			LaticeXY[X + A*off3] = log(LaticeXY[X + A*off3]);
			LaticeXY[X + A*off3] = LaticeXY[X + A*off3] + MinV;

		}

	}
	return(1);
}

double MyMathFuncs::ReverseCP(int SLen, int NumberAB, int NumberXY, double *ValXY, double *OptXY, unsigned char *RecodeB, double *LaticeXY, double *TransitionM2, double *EmissionM2)
{


	//transition(numberxy-1, numberxy-1)
	//optxy(numberxy-1, numberxy-1)
	//emmision(numberab-1,numberxy-1)
	//laticeab(slen,numberxy-1)
	//laticexy(slen,numberxy-1)
	int A, B, X, off1, off2, off3;
	double MinV;

	off1 = NumberXY;
	off2 = NumberAB;
	off3 = SLen + 1;

	for (X = SLen - 1; X >= 0; X--) {

		for (A = 0; A < NumberXY; A++) { //first state
			for (B = 0; B < NumberXY; B++) { //second state
				OptXY[A + B*off1] = LaticeXY[X + 1 + B*off3] + TransitionM2[A + B*off1] + EmissionM2[RecodeB[X + 1] + B*off2];

			}
		}


		//ReDim ValXY(NumberXY - 1)
		for (A = 0; A < NumberXY; A++)
			ValXY[A] = 0;

		for (A = 0; A < NumberXY; A++) {
			MinV = -10000000000000000;
			for (B = 0; B < NumberXY; B++) {
				if (MinV < OptXY[A + B*off1])
					MinV = OptXY[A + B*off1];

			}
			LaticeXY[X + A*off3] = 0;
			for (B = 0; B < NumberXY; B++) {
				ValXY[B] = OptXY[A + B*off1] - MinV;
				ValXY[B] = exp(ValXY[B]);
				LaticeXY[X + A*off3] = LaticeXY[X + A*off3] + ValXY[B];

			}
			LaticeXY[X + A*off3] = log(LaticeXY[X + A*off3]);
			LaticeXY[X + A*off3] = LaticeXY[X + A*off3] + MinV;

		}

	}
	return(1);
}

float MyMathFuncs::ReverseCPF(int SLen, int NumberAB, int NumberXY, float *ValXY, float *OptXY, unsigned char *RecodeB, float *LaticeXY, float *TransitionM2, float *EmissionM2)
{


	//transition(numberxy-1, numberxy-1)
	//optxy(numberxy-1, numberxy-1)
	//emmision(numberab-1,numberxy-1)
	//laticeab(slen,numberxy-1)
	//laticexy(slen,numberxy-1)
	int A, B, X, off1, off2, off3;
	double MinV;

	off1 = NumberXY;
	off2 = NumberAB;
	off3 = SLen + 1;

	for (X = SLen - 1; X >= 0; X--) {

		for (A = 0; A < NumberXY; A++) { //first state
			for (B = 0; B < NumberXY; B++) { //second state
				OptXY[A + B*off1] = LaticeXY[X + 1 + B*off3] + TransitionM2[A + B*off1] + EmissionM2[RecodeB[X + 1] + B*off2];

			}
		}


		//ReDim ValXY(NumberXY - 1)
		for (A = 0; A < NumberXY; A++)
			ValXY[A] = 0;

		for (A = 0; A < NumberXY; A++) {
			MinV = -10000000000000000;
			for (B = 0; B < NumberXY; B++) {
				if (MinV < OptXY[A + B*off1])
					MinV = OptXY[A + B*off1];

			}
			LaticeXY[X + A*off3] = 0;
			for (B = 0; B < NumberXY; B++) {
				ValXY[B] = OptXY[A + B*off1] - MinV;
				ValXY[B] = exp(ValXY[B]);
				LaticeXY[X + A*off3] = LaticeXY[X + A*off3] + ValXY[B];

			}
			LaticeXY[X + A*off3] = log(LaticeXY[X + A*off3]);
			LaticeXY[X + A*off3] = LaticeXY[X + A*off3] + MinV;

		}

	}
	return(1);
}




int MyMathFuncs::FillRmat(int Y, int Nextno, int UBRM1, int UBRM2, int UBDM1, int UBDM2, int UBDM3, double *RMat, double *DistMat, unsigned char *ZP) {
	int X, Z, osd1, osd2, osd3, osd12, osd123, osr1, osr2, osr12, os1, os2, os3, os4, os5, os6,os7;

	osd1 = UBDM1 + 1;
	osd2 = UBDM2 + 1;
	osd3 = UBDM3 + 1;
	osd12 = osd1*osd2;
	osd123 = osd12*osd3;
	osr1 = UBRM1 + 1;
	osr2 = UBRM2 + 1;
	osr12 = osr1*osr2;
	for (X = 0; X <= 1; X++) {
		if (X == 0) {

			ZP[0] = 0;
			ZP[1] = 1;
		}
		else {

			ZP[0] = 3;
			ZP[1] = 2;
		}

		for (Z = 0; Z <= Nextno; Z++) {
			os1 = Z*osr12;
			os2 = Z*osd12;
			os3 = Y + ZP[0] * osd1 + os2;
			os4 = Y + ZP[1] * osd1 + os2;
			os5 = Y + 4 * osd1 + os2;
			os6 = X + os1;
			os7 = 2 + os1;
			RMat[os6] = DistMat[os3];
			RMat[os6 + osr1 ] = DistMat[os3 + osd123];
			RMat[os6 + 2 * osr1] = DistMat[os3 + 2 * osd123];
			RMat[os6 + 3 * osr1] = DistMat[os4];
			RMat[os6 + 4 * osr1] = DistMat[os4 + osd123];
			RMat[os6 + 5 * osr1] = DistMat[os4 + 2 * osd123];


			RMat[os7] += DistMat[os3] / 2;
			RMat[os7 + osr1] += DistMat[os3 + osd123] / 2;
			RMat[os7 + 2 * osr1] += DistMat[os3 + 2 * osd123] / 2;

			RMat[os7 + 3 * osr1] += DistMat[os5] / 2;
			RMat[os7 + 4 * osr1 ] += DistMat[os5 + osd123] / 2;
			RMat[os7 + 5 * osr1] += DistMat[os5 + 2 * osd123] / 2;


		}
	}
	return(1);
}

int MyMathFuncs::FillIntTD(int UB, float *mindist, float *mintdist, float *adjustd, float *adjusttd, short int *IntTD, float *Distance, float *TreeDistance) {
	int X, Y, os1, os2;
	float MaxDist, MaxTDist, MinDist, MinTDist, IntervalD, IntervalTD, AdjustD, AdjustTD;
	MaxDist = -10000.0;
	MinDist = 10000.0;
	MaxTDist = -10000.0;
	MinTDist = 10000.0;
	os1 = UB + 1;
	for (Y = 0; Y <= UB; Y++) {
		os2 = os1*Y;
		for (X = Y+1; X <= UB; X++) {
			if (Distance[X + os2] > MaxDist) 
				MaxDist = Distance[X + os2];
			if (Distance[X + os2] < MinDist) 
				MinDist = Distance[X + os2];
			if (TreeDistance[X + os2] > MaxTDist)
				MaxTDist = TreeDistance[X + os2];
			if (TreeDistance[X + os2] < MinTDist)
				MinTDist = TreeDistance[X + os2];
		}
	}

	IntervalD = (float)(MaxDist - MinDist);

	AdjustD = (float)(64000.0 / IntervalD);

	IntervalTD = (float)(MaxTDist - MinTDist);

	AdjustTD = (float)(64000.0 / IntervalTD);
	int procs;
	procs = omp_get_num_procs();
	procs = procs/2 - 1;
	if (procs < 3)
		procs = 3;
	omp_set_num_threads(procs);
#pragma omp parallel for private (Y, X, os2)
	for (Y = 0; Y <= UB; Y++){
		os2 = os1*Y;
		for (X = Y+1; X <= UB; X++) 
			IntTD[X + os2] = (int)(((Distance[X + os2] - MinDist) * AdjustD) - 32000);
				
		
	}
#pragma omp parallel for private (Y, X, os2)
	for (Y = 1; Y <= UB; Y++){
		os2 = os1*Y;
		for (X = 0; X < Y; X++)
		
			IntTD[X + os2] = (int)(((TreeDistance[X + os2] - MinTDist) * AdjustTD) - 32000);
			
	}

	*mindist = MinDist;
	*mintdist = MinTDist;
	*adjustd = AdjustD;
	*adjusttd = AdjustTD;
	omp_set_num_threads(2);
	return(1);
}


int MyMathFuncs::ReadIntTD(int UB, float MinDist, float MinTDist, float AdjustD, float AdjustTD, short int *IntTD, float *Distance, float *TreeDistance) {
	int X, Y, os1, os2;
	
	os1 = UB + 1;
	int procs;
		procs = omp_get_num_procs();
		procs = procs/2 - 1;
		if (procs < 3)
			procs = 3;
		omp_set_num_threads(procs);
#pragma omp parallel for private (Y, X, os2)
	for (Y = 0; Y <= UB; Y++) {
		os2 = os1*Y;
		for (X = Y + 1; X <= UB; X++) {
			//Distance[X + os2] = (int)(((Distance[X + os2] - MinDist) * AdjustD) - 32000);
			Distance[X + os2] = (float)((IntTD[X + os2]) + 32000) / AdjustD + MinDist;
			Distance[Y + X*os1] = Distance[X + os2];
				//((CLng(IntTD(X, Y)) + 32000) / AdjustD) + MinDist
		}

	}
#pragma omp parallel for private (Y, X, os2)
	for (Y = 1; Y <= UB; Y++) {
		os2 = os1*Y;
		for (X = 0; X < Y; X++) {

			TreeDistance[X + os2] = (float)((IntTD[X + os2]) + 32000) / AdjustTD + MinTDist;
			TreeDistance[Y + X*os1] = TreeDistance[X + os2];
		}

	}

	omp_set_num_threads(2);
	return(1);
}

int MyMathFuncs::DoPermsXP(int LS, int SSWinLen, int SSNumPerms, char *PScores, char *VRandTemplate, char *VRandConv, int *PermPScores)
{
	int Z, B;
	int HN, HN1, HN2, HN3;
	int vo;
	int os2, os3;
	os2 = (SSNumPerms + 1);
	vo = LS + 1;
	for (Z = 1; Z <= SSWinLen; Z++) {
		HN = (int)(PScores[Z]);
		HN = HN*os2;
		PermPScores[HN] = PermPScores[HN] + 1;
	}

#pragma omp parallel for private (Z, os3, B, HN1, HN2, HN3)
	for (Z = 1; Z <= SSNumPerms; Z++) {
		os3 = Z*vo;
		for (B = 1; B <= SSWinLen; B++) {
			HN1 = (int)(PScores[B]);
			HN2 = (int)(VRandTemplate[B + os3]);
			HN2 = HN1 + HN2 * 16;
			HN3 = VRandConv[HN2];
			HN3 = Z + HN3*os2;
//#pragma  atomic
			PermPScores[HN3]++;
			
		}
	}
	return 1;
}

int MyMathFuncs::DoPerms3P(int LS, int SSWinLen, int SSNumPerms, int SSNumPerms2, int *PScores, char *VRandTemplate, char *VRandConv, int *PermPScores)
{
	int Z, B;
	int a, c, HN1, os, os2, os3, os4;
	int vo;
	vo = LS + 1;

	//PermPScores[0] = PScores[0];
	//PermPScores[1] = PScores[1]; 
	//PermPScores[15] = PScores[15];
	//PermPScores -numprems,15
	os = (SSNumPerms + 1);

	for (B = 2; B <= 14; B++) {
		os2 = B*os;
		for (Z = 0; Z <= SSNumPerms2; Z++)
			PermPScores[Z + os2] = 0;



	}

	for (Z = 0; Z <= 15; Z++)
		PermPScores[Z*os] = PScores[Z];


	os3 = 15 * os;

	c = 0;
	for (B = 2; B <= 14; B++) {
		if (PScores[B] > 0) {

			HN1 = PScores[B];
			for (a = 0; a < HN1; a++) {
				os4 = c - 1;
				for (Z = 1; Z <= SSNumPerms2; Z++)
					PermPScores[Z + os*VRandConv[B + 16 * VRandTemplate[Z + os4]]] += 1;


				c = c + SSNumPerms2;
			}
		}
	}
	return c;
}

double MyMathFuncs::UpdateCountsP(int SLen, int NumberABC, int NumberXY, int *LaticePath, unsigned char *RecodeB, double *TransitionCount, double *StateCount)
{

	//statecount(numberabc-1, numberxy)
	//transitioncount(numberxy-1,numberxy-1)
	int Z, off1, off2;

	off1 = NumberXY;
	off2 = NumberABC;
	for (Z = 0; Z < SLen; Z++)
		TransitionCount[LaticePath[Z] + LaticePath[Z + 1] * off1] = TransitionCount[LaticePath[Z] + LaticePath[Z + 1] * off1] + 1;



	for (Z = 0; Z <= SLen; Z++)
		StateCount[RecodeB[Z] + LaticePath[Z] * off2] = StateCount[RecodeB[Z] + LaticePath[Z] * off2] + 1;

	return(1);
}

float MyMathFuncs::UpdateCountsPF(int SLen, int NumberABC, int NumberXY, int *LaticePath, unsigned char *RecodeB, float *TransitionCount, float *StateCount)
{

	//statecount(numberabc-1, numberxy)
	//transitioncount(numberxy-1,numberxy-1)
	int Z, off1, off2;

	off1 = NumberXY;
	off2 = NumberABC;
	for (Z = 0; Z < SLen; Z++)
		TransitionCount[LaticePath[Z] + LaticePath[Z + 1] * off1] = TransitionCount[LaticePath[Z] + LaticePath[Z + 1] * off1] + 1;



	for (Z = 0; Z <= SLen; Z++)
		StateCount[RecodeB[Z] + LaticePath[Z] * off2] = StateCount[RecodeB[Z] + LaticePath[Z] * off2] + 1;

	return(1);
}

double MyMathFuncs::GetLaticePathP(int SLen, int NumberXY, double *LaticeXY, double *LaticeAB, int *LaticePath)
{

	int Y, X, off1;
	double MaxL;

	off1 = SLen + 1;

	MaxL = -10000000000000000;

	for (Y = 0; Y < NumberXY; Y++) {

		if (LaticeXY[SLen + Y*off1] > MaxL) {
			MaxL = LaticeXY[SLen + Y*off1];
			LaticePath[SLen] = (int)(LaticeAB[SLen + Y*off1]);
		}

	}

	for (X = SLen - 1; X >= 0; X--) {
		//for (Y = 0; Y< NumberXY; Y++)

		LaticePath[X] = (int)(LaticeAB[X + LaticePath[X + 1] * off1]);


	}
	return(MaxL);
}

float MyMathFuncs::GetLaticePathPF(int SLen, int NumberXY, float *LaticeXY, float *LaticeAB, int *LaticePath)
{

	int Y, X, off1;
	float MaxL;

	off1 = SLen + 1;

	MaxL = -100000000000;

	for (Y = 0; Y < NumberXY; Y++) {

		if (LaticeXY[SLen + Y*off1] > MaxL) {
			MaxL = LaticeXY[SLen + Y*off1];
			LaticePath[SLen] = (int)(LaticeAB[SLen + Y*off1]);
		}

	}

	for (X = SLen - 1; X >= 0; X--) {
		//for (Y = 0; Y< NumberXY; Y++)

		LaticePath[X] = (int)(LaticeAB[X + LaticePath[X + 1] * off1]);


	}
	return(MaxL);
}

int MyMathFuncs::CopyCharArray(int ub1, int ub2, unsigned char *From, unsigned char *To) {

	int x, y, target;
	target = ub1 + (ub1 + 1)*ub2;
	

	for (x = 0; x <= target; x++) 
			To[x] = From[x];
	

	return(1);

}
int MyMathFuncs::CollapseNodesXP(int NextNo, float CutOff, int *T, float *DLen, float *CMat, float *AMat, float *TraceBak) {

	int X, Y, Z, A, GoOn, off1, off2, off3;
	float LODist;

	off1 = NextNo + 1;
	for (X = 0; X <= NextNo; X++) {
		if (DLen[X] < CutOff) {

			//First find two sequences with this dist
			GoOn = 0;
			LODist = 100000;
			for (Y = 0; Y < NextNo; Y++) {
				off2 = Y*off1;
				for (Z = Y + 1; Z <= NextNo; Z++) {

					if (AMat[Z + off2] == TraceBak[X]) {
						GoOn = 1;
						T[0] = Y;
						T[1] = Z;
						off3 = Z*off1;
						for (A = 0; A <= NextNo; A++) {
							if (AMat[A + off2] == AMat[A + off3]) {
								if (AMat[A + off2] > TraceBak[X] && AMat[A + off2] < LODist)
									LODist = AMat[A + off2];

							}
						}
					}
				}
				//If GoOn = 1 Then Exit For
			}

			//Find next lowest dist
			if (GoOn == 1) {

				if (LODist < 100000) {
					for (Y = 0; Y < NextNo; Y++) {
						off2 = Y*off1;
						for (Z = Y + 1; Z <= NextNo; Z++) {

							if (AMat[Z + off2] == TraceBak[X]) {
								AMat[Y + Z*off1] = (float)(LODist);
								AMat[Z + off2] = (float)(LODist);

							}
						}
					}
				}
				else {
					for (Y = 0; Y < NextNo; Y++) {
						off2 = Y*off1;
						for (Z = Y + 1; Z <= NextNo; Z++) {

							if (AMat[Z + off2] == TraceBak[X]) {
								CMat[Y + Z*off1] = TraceBak[X];
								CMat[Z + off2] = TraceBak[X];

							}

						}
					}
				}

			}
		}
		else {
			for (Y = 0; Y < NextNo; Y++) {
				off2 = Y*off1;
				for (Z = Y + 1; Z <= NextNo; Z++) {
					if (AMat[Z + off2] == TraceBak[X]) {
						CMat[Y + Z*off1] = TraceBak[X];
						CMat[Z + off2] = TraceBak[X];

					}
				}
			}
		}

	}
	return(1);

}

int MyMathFuncs::CollapseNodesXP2(int NextNo, float CutOff, int *T, float *DLen, float *CMat, float *PAMat, float *TraceBak) {

	int X, Y, Z, A, GoOn, off1, off2, off3;
	float LODist;
	float *AMat;
	off1 = NextNo + 1;

	AMat =  (float *)calloc((NextNo+1)*(NextNo+1), sizeof(float));

	for (X = 0; X <= NextNo; X++){
		for (Y = 0; Y <= NextNo; Y++)
			AMat[X + Y*off1] = PAMat[X + Y*off1];
	}


	
	for (X = 0; X <= NextNo; X++) {
		if (DLen[X] < CutOff) {

			//First find two sequences with this dist
			GoOn = 0;
			LODist = 100000;
			for (Y = 0; Y < NextNo; Y++) {
				off2 = Y*off1;
				for (Z = Y + 1; Z <= NextNo; Z++) {

					if (AMat[Z + off2] == TraceBak[X]) {
						GoOn = 1;
						T[0] = Y;
						T[1] = Z;
						off3 = Z*off1;
						for (A = 0; A <= NextNo; A++) {
							if (AMat[A + off2] == AMat[A + off3]) {
								if (AMat[A + off2] > TraceBak[X] && AMat[A + off2] < LODist)
									LODist = AMat[A + off2];

							}
						}
					}
				}
				//If GoOn = 1 Then Exit For
			}

			//Find next lowest dist
			if (GoOn == 1) {

				if (LODist < 100000) {
					for (Y = 0; Y < NextNo; Y++) {
						off2 = Y*off1;
						for (Z = Y + 1; Z <= NextNo; Z++) {

							if (AMat[Z + off2] == TraceBak[X]) {
								AMat[Y + Z*off1] = (float)(LODist);
								AMat[Z + off2] = (float)(LODist);

							}
						}
					}
				}
				else {
					for (Y = 0; Y < NextNo; Y++) {
						off2 = Y*off1;
						for (Z = Y + 1; Z <= NextNo; Z++) {

							if (AMat[Z + off2] == TraceBak[X]) {
								CMat[Y + Z*off1] = TraceBak[X];
								CMat[Z + off2] = TraceBak[X];

							}

						}
					}
				}

			}
		}
		else {
			for (Y = 0; Y < NextNo; Y++) {
				off2 = Y*off1;
				for (Z = Y + 1; Z <= NextNo; Z++) {
					if (AMat[Z + off2] == TraceBak[X]) {
						CMat[Y + Z*off1] = TraceBak[X];
						CMat[Z + off2] = TraceBak[X];

					}
				}
			}
		}

	}
	free(AMat);
	return(1);

}


int MyMathFuncs::CollapseNodesXP3(int NextNo, float CutOff, int *T, float *DLen, float *CMat, float *PAMat, float *TraceBak) {

	int X, Y, Z, A, GoOn, off1, off2, off3;
	float LODist;
	float *AMat, *cm;
	//short int *sp;
	off1 = NextNo + 1;

	cm = (float *)calloc((NextNo + 1), sizeof(float));
	AMat = (float *)calloc((NextNo + 1)*(NextNo + 1), sizeof(float));
	//sp = (short int *)calloc(((NextNo + 1)*((NextNo + 1)/2), sizeof(short int));

	for (X = 0; X <= NextNo; X++)
		cm [X]= 10000;
	//LODist = 1000;
	for (X = 0; X <= NextNo; X++) {
		for (Y = 0; Y <= NextNo; Y++) {
			if (cm[X] > PAMat[X + Y*off1])
				cm[X] = PAMat[X + Y*off1];
			AMat[X + Y*off1] = PAMat[X + Y*off1];
		}
	}


	GoOn=0;
	for (X = 0; X <= NextNo; X++) {
		if (DLen[X] < CutOff) {//if branch x has less support than the bootstrap cutoff
								//traceback x indicates the branch number in amat that corresponds to branch x

			//First find two sequences with this dist
			GoOn = 0;
			LODist = 100000;
			for (Y = 0; Y < NextNo; Y++) {
				off2 = Y*off1;
				for (Z = Y + 1; Z <= NextNo; Z++) {

					if (AMat[Z + off2] == TraceBak[X]) {	//if branch x is immediately below the node (tree is midpoint rooted) 
															//of the MRCA of sequences Y and Z 
						GoOn = 1;
						/*T[0] = Y;
						T[1] = Z;*/
						off3 = Z*off1;
						for (A = 0; A <= NextNo; A++) {     
							if (AMat[A + off2] == AMat[A + off3]) { //if the MRCA node of sequence A and Z is the same as the MRCA node 
																	// of sequence A and Y (i.e. the A-Z/A-Y MRCA node is below the Y-Z MRCA node
																	//i.e sequence A is on the other side of branch X to seqs Z and Y 
								if (AMat[A + off2] > TraceBak[X] && AMat[A + off2] < LODist) { //this finds the closest outlyer on the other side of branch x.
																	// or rather it finds the number of the node (stored in lodist) of the node
																	// on the other side of branch x
									LODist = AMat[A + off2];
									//GoOn = 1;
									
								}

							}
						}
						//if (GoOn==1)
							break;
					}
					
				}
				if (GoOn == 1)
					break;
			}

			//Find next lowest dist
			if (GoOn == 1) {

				if (LODist < 100000) {//if there was a node below branch x (i.e. if branch x is not the "root branch")
					for (Y = 0; Y < NextNo; Y++) {
						off2 = Y*off1;
						for (Z = Y + 1; Z <= NextNo; Z++) {

							if (AMat[Z + off2] == TraceBak[X]) {//if the MRCA of sequence Z and Y is the node above branch x then
																//change the MRCA of sequnece Z and Y to the node below branch X
								AMat[Y + Z*off1] = (float)(LODist);
								AMat[Z + off2] = (float)(LODist);

							}
						}
					}
				}
				else {//if branch x is the root branch then dont collapse it and transfer the node labels for sequence pairs with this node as a
						//a mrca from AMat to CMat
					for (Y = 0; Y < NextNo; Y++) {
						off2 = Y*off1;
						for (Z = Y + 1; Z <= NextNo; Z++) {

							if (AMat[Z + off2] == TraceBak[X]) {
								CMat[Y + Z*off1] = TraceBak[X];
								CMat[Z + off2] = TraceBak[X];

							}

						}
					}
				}

			}
		}
		else {//if the branch is supported then transfer the node labels for sequence pairs with node as an MRCA from AMat to CMat
			for (Y = 0; Y < NextNo; Y++) {
				off2 = Y*off1;
				for (Z = Y + 1; Z <= NextNo; Z++) {
					if (AMat[Z + off2] == TraceBak[X]) {
						CMat[Y + Z*off1] = TraceBak[X];
						CMat[Z + off2] = TraceBak[X];

					}
				}
			}
		}

	}
	free(AMat);
	free(cm);
	return(1);

}


DMAT *
NJ_parse_distance_matrix(float *dists,int nextno) {

	DMAT *dmat = NULL;

	int state, dmat_type;
	int row;
	int fltcnt;
	int x, y, i, xx, j;
	int numvalread;
	int expectedvalues = -1;
	float val;
	int first_state = 0;


	/* allocate our distance matrix and token structure */
	dmat = (DMAT *)calloc(1, sizeof(DMAT));
	
	dmat->ntaxa = nextno + 1;

	/* set our initial working size according to the # of taxa */
	dmat->size = dmat->ntaxa;

	/* allocate space for the distance matrix values here */
	dmat->val =
		(float *)calloc(NJ_NCELLS(dmat->ntaxa), sizeof(float));
	

	/*  taxa names */
	dmat->taxaname = (char **)calloc(dmat->ntaxa, sizeof(char *));
	
	for (xx = 0; xx <= nextno; xx++) {
		dmat->taxaname[xx] = (char *)calloc(3, sizeof(char));
	}

	/* set the initial state of our state machine */
	dmat_type = NJ_PARSE_SYMMETRIC;
	row = -1;
	fltcnt = 0;
	numvalread = 0;


	//this is where clearcut parses the distance matrix file


	//dmat-> val = dists;

	for (i= 0; i < nextno; i++){
		for (j = i+1; j<=nextno;j++)
			dmat->val[NJ_MAP(i, j, dmat->size)] = dists[i+j*(nextno+1)];
	}




	/* now lets allocate space for the r and r2 columns */
	dmat->r = (float *)calloc(dmat->ntaxa, sizeof(float));
	dmat->r2 = (float *)calloc(dmat->ntaxa, sizeof(float));
	
	/* track some memory addresses */
	dmat->rhandle = dmat->r;
	dmat->r2handle = dmat->r2;
	dmat->valhandle = dmat->val;

	return(dmat);



	/* clean up our partial progress */

}


void
NJ_init_r(DMAT *dmat) {

	long int i, j, size;
	long int index;
	float *r, *r2, *val;
	long int size1;
	float size2;

	r = dmat->r;
	r2 = dmat->r2;
	val = dmat->val;
	size = dmat->size;
	size1 = size - 1;
	size2 = (float)(size - 2);

	index = 0;
	for (i = 0; i<size1; i++) {
		index++;
		for (j = i + 1; j<size; j++) {
			r[i] += val[index];
			r[j] += val[index];
			index++;
		}

		r2[i] = r[i] / size2;
	}

	return;
}

NJ_VERTEX *
NJ_init_vertex(DMAT *dmat) {

	long int i;
	NJ_VERTEX *vertex;

	/* allocate the vertex here */
	vertex = (NJ_VERTEX *)calloc(1, sizeof(NJ_VERTEX));

	/* allocate the nodes in the vertex */
	vertex->nodes = (NJ_TREE **)calloc(dmat->ntaxa, sizeof(NJ_TREE *));
	vertex->nodes_handle = vertex->nodes;

	/* initialize our size and active variables */
	vertex->nactive = dmat->ntaxa;
	vertex->size = dmat->ntaxa;

	/* initialize the nodes themselves */
	for (i = 0; i<dmat->ntaxa; i++) {

		vertex->nodes[i] = (NJ_TREE *)calloc(1, sizeof(NJ_TREE));

		vertex->nodes[i]->left = NULL;
		vertex->nodes[i]->right = NULL;

		vertex->nodes[i]->taxa_index = i;
	}

	return(vertex);
}

float
NJ_min_transform(DMAT *dmat,
	long int *ret_i,
	long int *ret_j) {

	long int i, j;   /* indices used for looping        */
	long int tmp_i = 0;/* to limit pointer dereferencing  */
	long int tmp_j = 0;/* to limit pointer dereferencing  */
	float smallest;  /* track the smallest trans. dist  */
	float curval;    /* the current trans. dist in loop */

	float *ptr;      /* pointer into distance matrix    */
	float *r2;       /* pointer to r2 matrix for computing transformed dists */

	smallest = (float)HUGE_VAL;

	/* track these here to limit pointer dereferencing in inner loop */
	ptr = dmat->val;
	r2 = dmat->r2;

	/* for every row */
	for (i = 0; i<dmat->size; i++) {
		ptr++;  /* skip diagonal */
		for (j = i + 1; j<dmat->size; j++) {   /* for every column */

											   /* find transformed distance in matrix at i, j */
			curval = *(ptr++) - (r2[i] + r2[j]);

			/* if the transformed distanance is less than the known minimum */
			if (curval < smallest) {

				smallest = curval;
				tmp_i = i;
				tmp_j = j;
			}
		}
	}

	/* pass back (by reference) the coords of the min. transformed distance */
	*ret_i = tmp_i;
	*ret_j = tmp_j;

	return(smallest);  /* return the min transformed distance */
}


NJ_TREE *
NJ_decompose(DMAT *dmat,
	NJ_VERTEX *vertex,
	long int x,
	long int y,
	int last_flag) {

	NJ_TREE *new_node;
	float x2clade, y2clade;

	/* compute the distance from the clade components to the new node */
	if (last_flag) {
		x2clade =
			(dmat->val[NJ_MAP(x, y, dmat->size)]);
	}
	else {
		x2clade =
			(dmat->val[NJ_MAP(x, y, dmat->size)]) / 2 +
			((dmat->r2[x] - dmat->r2[y]) / 2);
	}

	vertex->nodes[x]->dist = x2clade;

	if (last_flag) {
		y2clade =
			(dmat->val[NJ_MAP(x, y, dmat->size)]);
	}
	else {
		y2clade =
			(dmat->val[NJ_MAP(x, y, dmat->size)]) / 2 +
			((dmat->r2[y] - dmat->r2[x]) / 2);
	}

	vertex->nodes[y]->dist = y2clade;

	/* allocate new node to connect two sub-clades */
	new_node = (NJ_TREE *)calloc(1, sizeof(NJ_TREE));

	new_node->left = vertex->nodes[x];
	new_node->right = vertex->nodes[y];
	new_node->taxa_index = NJ_INTERNAL_NODE;  /* this is not a terminal node, no taxa index */

	if (last_flag) {
		return(new_node);
	}

	vertex->nodes[x] = new_node;
	vertex->nodes[y] = vertex->nodes[0];

	vertex->nodes = &(vertex->nodes[1]);

	vertex->nactive--;

	return(new_node);
}

static inline
void
NJ_compute_r(DMAT *dmat,
	long int a,
	long int b) {

	long int i;         /* a variable used in indexing */
	float *ptrx, *ptry; /* pointers into the distance matrix */

						/* some variables to limit pointer dereferencing in loop */
	long int size;
	float *r, *val;

	/* to limit pointer dereferencing */
	size = dmat->size;
	val = dmat->val;
	r = dmat->r + a + 1;

	/*
	* Loop through the rows and decrement the stored r values
	* by the distances stored in the rows and columns of the distance
	* matrix which are being removed post-join.
	*
	* We do the rows altogether in order to benefit from cache locality.
	*/
	ptrx = &(val[NJ_MAP(a, a + 1, size)]);
	ptry = &(val[NJ_MAP(b, b + 1, size)]);

	for (i = a + 1; i<size; i++) {
		*r -= *(ptrx++);

		if (i>b) {
			*r -= *(ptry++);
		}

		r++;
	}

	/* Similar to the above loop, we now do the columns */
	ptrx = &(val[NJ_MAP(0, a, size)]);
	ptry = &(val[NJ_MAP(0, b, size)]);
	r = dmat->r;
	for (i = 0; i<b; i++) {
		if (i<a) {
			*r -= *ptrx;
			ptrx += size - i - 1;
		}

		*r -= *ptry;
		ptry += size - i - 1;
		r++;
	}

	return;
}


static inline
void
NJ_collapse(DMAT *dmat,
	NJ_VERTEX *vertex,
	long int a,
	long int b) {


	long int i;     /* index used for looping */
	long int size;  /* size of dmat --> reduce pointer dereferencing */
	float a2clade;  /* distance from a to the new node that joins a and b */
	float b2clade;  /* distance from b to the new node that joins a and b */
	float cval;     /* stores distance information during loop */
	float *vptr;    /* pointer to elements in first row of dist matrix */
	float *ptra;    /* pointer to elements in row a of distance matrix */
	float *ptrb;    /* pointer to elements in row b of distance matrix */

	float *val, *r, *r2;  /* simply used to limit pointer dereferencing */


	

	/* some shortcuts to help limit dereferencing */
	val = dmat->val;
	r = dmat->r;
	r2 = dmat->r2;
	size = dmat->size;

	/* compute the distance from the clade components (a, b) to the new node */
	a2clade =
		((val[NJ_MAP(a, b, size)]) + (dmat->r2[a] - dmat->r2[b])) / 2.0;
	b2clade =
		((val[NJ_MAP(a, b, size)]) + (dmat->r2[b] - dmat->r2[a])) / 2.0;


	r[a] = 0.0;  /* we are removing row a, so clear dist. in r */

				 /*
				 * Fill the horizontal part of the "a" row and finish computing r and r2
				 * we handle the horizontal component first to maximize cache locality
				 */
	ptra = &(val[NJ_MAP(a, a + 1, size)]);   /* start ptra at the horiz. of a  */
	ptrb = &(val[NJ_MAP(a + 1, b, size)]);   /* start ptrb at comparable place */
	for (i = a + 1; i<size; i++) {

		/*
		* Compute distance from new internal node to others in
		* the distance matrix.
		*/
		cval =
			((*ptra - a2clade) +
			(*ptrb - b2clade)) / 2.0;

		/* incr.  row b pointer differently depending on where i is in loop */
		if (i<b) {
			ptrb += size - i - 1;  /* traverse vertically  by incrementing by row */
		}
		else {
			ptrb++;            /* traverse horiz. by incrementing by column   */
		}

		/* assign the newly computed distance and increment a ptr by a column */
		*(ptra++) = cval;

		/* accumulate the distance onto the r vector */
		r[a] += cval;
		r[i] += cval;

		/* scale r2 on the fly here */
		r2[i] = r[i] / (float)(size - 3);
	}

	/* fill the vertical part of the "a" column and finish computing r and r2 */
	ptra = val + a;  /* start at the top of the columb for "a" */
	ptrb = val + b;  /* start at the top of the columb for "b" */
	for (i = 0; i<a; i++) {

		/*
		* Compute distance from new internal node to others in
		* the distance matrix.
		*/
		cval =
			((*ptra - a2clade) +
			(*ptrb - b2clade)) / 2.0;

		/* assign the newly computed distance and increment a ptr by a column */
		*ptra = cval;

		/* accumulate the distance onto the r vector */
		r[a] += cval;
		r[i] += cval;

		/* scale r2 on the fly here */
		r2[i] = r[i] / (float)(size - 3);

		/* here, always increment by an entire row */
		ptra += size - i - 1;
		ptrb += size - i - 1;
	}


	/* scale r2 on the fly here */
	r2[a] = r[a] / (float)(size - 3);



	/*
	* Copy row 0 into row b.  Again, the code is structured into two
	* loops to maximize cache locality for writes along the horizontal
	* component of row b.
	*/
	vptr = val;
	ptrb = val + b;
	for (i = 0; i<b; i++) {
		*ptrb = *(vptr++);
		ptrb += size - i - 1;
	}
	vptr++;  /* skip over the diagonal */
	ptrb = &(val[NJ_MAP(b, b + 1, size)]);
	for (i = b + 1; i<size; i++) {
		*(ptrb++) = *(vptr++);
	}

	/*
	* Collapse r here by copying contents of r[0] into r[b] and
	* incrementing pointer to the beginning of r by one row
	*/
	r[b] = r[0];
	dmat->r = r + 1;


	/*
	* Collapse r2 here by copying contents of r2[0] into r2[b] and
	* incrementing pointer to the beginning of r2 by one row
	*/
	r2[b] = r2[0];
	dmat->r2 = r2 + 1;

	/* increment dmat pointer to next row */
	dmat->val += size;

	/* decrement the total size of the distance matrix by one row */
	dmat->size--;

	return;
}




NJ_TREE *NJ_neighbor_joining(DMAT *dmat, int outlyer) {


	NJ_TREE   *tree = NULL;
	NJ_VERTEX *vertex = NULL;

	long int a, b;
	float min;


	/* initialize the r and r2 vectors */
	NJ_init_r(dmat);

	/* allocate and initialize our vertex vector used for tree construction */
	vertex = NJ_init_vertex(dmat);
	

	/* we iterate until the working distance matrix has only 2 entries */
	while (vertex->nactive > 2) {

		/*
		* Find the global minimum transformed distance from the distance matrix
		*/
		min = NJ_min_transform(dmat, &a, &b);

		/*
		* Build the tree by removing nodes a and b from the vertex array
		* and inserting a new internal node which joins a and b.  Collapse
		* the vertex array similarly to how the distance matrix and r and r2
		* are compacted.
		*/
		NJ_decompose(dmat, vertex, a, b, 0);

		/* decrement the r and r2 vectors by the distances corresponding to a, b */
		NJ_compute_r(dmat, a, b);

		/* compact the distance matrix and the r and r2 vectors */
		NJ_collapse(dmat, vertex, a, b);
	}

	/* Properly join the last two nodes on the vertex list */
	tree = NJ_decompose(dmat, vertex, 0, 1, NJ_LAST);

	/* return the computed tree to the calling function */
	return(tree);
}


static inline
void
NJ_permute(long int *perm,
	long int size) {

	long int i;     /* index used for looping */
	long int swap;  /* we swap values to generate permutation */
	long int tmp;   /* used for swapping values */
	double K;

					/* check to see if vector of long ints is valid */
	if (!perm) {
		
		exit(-1);
	}

	/* init permutation as an ordered list of integers */
	for (i = 0; i<size; i++) {
		perm[i] = i;
	}

	/*
	* Iterate across the array from i = 0 to size -1, swapping ith element
	* with a randomly chosen element from a changing range of possible values
	*/
	for (i = 0; i<size; i++) {

		/* choose which element we will swap with */
		K = rand();
		swap = i + (int)((K / RAND_MAX)*(size - 1));
		//swap = i + NJ_genrand_int31_top(size - i);

		/* swap elements here */
		if (i != swap) {
			tmp = perm[swap];
			perm[swap] = perm[i];
			perm[i] = tmp;
		}
	}

	return;
}

static inline
float
NJ_find_hmin(DMAT *dmat,
	long int a,
	long int *min,
	long int *hmincount) {

	long int i;     /* index variable for looping                    */
	int size;       /* current size of distance matrix               */
	int mindex = 0; /* holds the current index to the chosen minimum */
	float curval;   /* used to hold current transformed values       */
	float hmin;     /* the value of the transformed minimum          */

	float *ptr, *r2, *val;  /* pointers used to reduce dereferencing in inner loop */

							/* values used for stochastic selection among multiple minima */
	float p, x;
	long int smallcnt;

	/* initialize the min to something large */
	hmin = (float)HUGE_VAL;

	/* setup some pointers to limit dereferencing later */
	r2 = dmat->r2;
	val = dmat->val;
	size = dmat->size;

	/* initialize values associated with minima tie breaking */
	p = 1.0;
	smallcnt = 0;


	ptr = &(val[NJ_MAP(a, a + 1, size)]);   /* at the start of the horiz. part */
	for (i = a + 1; i<size; i++) {

		curval = *(ptr++) - (r2[a] + r2[i]);  /* compute transformed distance */

		if (NJ_FLT_EQ(curval, hmin)) {  /* approx. equal */

			smallcnt++;

			p = 1.0 / (float)smallcnt;
			x = rand()/RAND_MAX;
			//x = genrand_real2();

			/* select this minimum in a way which is proportional to
			the number of minima found along the row so far */
			if (x < p) {
				mindex = i;
			}

		}
		else if (curval < hmin) {

			smallcnt = 1;
			hmin = curval;
			mindex = i;
		}
	}

	/* save off the the minimum index to be returned via reference */
	*min = mindex;

	/* save off the number of minima */
	*hmincount = smallcnt;

	/* return the value of the smallest tranformed distance */
	return(hmin);
}


static inline
float
NJ_find_vmin(DMAT *dmat,
	long int a,
	long int *min,
	long int *vmincount) {

	long int i;         /* index variable used for looping */
	long int size;      /* track the size of the matrix    */
	long int mindex = 0;/* track the index to the minimum  */
	float curval;       /* track value of current transformed distance  */
	float vmin;         /* the index to the smallest "vertical" minimum */

						/* pointers which are used to reduce pointer dereferencing in inner loop */
	float *ptr, *r2, *val;

	/* values used in stochastically breaking ties */
	float p, x;
	long int smallcnt;

	/* initialize the vertical min to something really big */
	vmin = (float)HUGE_VAL;

	/* save off some values to limit dereferencing later */
	r2 = dmat->r2;
	val = dmat->val;
	size = dmat->size;

	p = 1.0;
	smallcnt = 0;

	/* start on the first row and work down */
	ptr = &(val[NJ_MAP(0, a, size)]);
	for (i = 0; i<a; i++) {

		curval = *ptr - (r2[i] + r2[a]);  /* compute transformed distance */

		if (NJ_FLT_EQ(curval, vmin)) {  /* approx. equal */

			smallcnt++;

			p = 1.0 / (float)smallcnt;
			x = rand()/RAND_MAX;
			//x = genrand_real2();

			/* break ties stochastically to avoid systematic bias */
			if (x < p) {
				mindex = i;
			}

		}
		else if (curval < vmin) {

			smallcnt = 1;
			vmin = curval;
			mindex = i;
		}

		/* increment our working pointer to the next row down */
		ptr += size - i - 1;
	}

	/* pass back the index to the minimum found so far (by reference) */
	*min = mindex;

	/* pass back the number of minima along the vertical */
	*vmincount = smallcnt;

	/* return the value of the smallest transformed distance */
	return(vmin);
}

static inline
int
NJ_check_additivity(DMAT *dmat,
	long int a,
	long int b) {

	float a2clade, b2clade;
	float clade_dist;
	long int target;


	/* determine target taxon here */
	if (b == dmat->size - 1) {
		/* if we can't do a row here, lets do a column */
		if (a == 0) {
			if (b == 1) {
				target = 2;
			}
			else {
				target = 1;
			}
		}
		else {
			target = 0;
		}
	}
	else {
		target = b + 1;
	}


	/* distance between a and the root of clade (a,b) */
	a2clade =
		((dmat->val[NJ_MAP(a, b, dmat->size)]) +
		(dmat->r2[a] - dmat->r2[b])) / 2.0;

	/* distance between b and the root of clade (a,b) */
	b2clade =
		((dmat->val[NJ_MAP(a, b, dmat->size)]) +
		(dmat->r2[b] - dmat->r2[a])) / 2.0;

	/* distance between the clade (a,b) and the target taxon */
	if (b<target) {

		/* compute the distance from the clade root to the target */
		clade_dist =
			((dmat->val[NJ_MAP(a, target, dmat->size)] - a2clade) +
			(dmat->val[NJ_MAP(b, target, dmat->size)] - b2clade)) / 2.0;

		/*
		* Check to see that distance from clade root to target + distance from
		*  b to clade root are equal to the distance from b to the target
		*/
		if (NJ_FLT_EQ(dmat->val[NJ_MAP(b, target, dmat->size)],
			(clade_dist + b2clade))) {
			return(1);  /* join is legitimate   */
		}
		else {
			return(0);  /* join is illigitimate */
		}

	}
	else {

		/* compute the distance from the clade root to the target */
		clade_dist =
			((dmat->val[NJ_MAP(target, a, dmat->size)] - a2clade) +
			(dmat->val[NJ_MAP(target, b, dmat->size)] - b2clade)) / 2.0;

		/*
		* Check to see that distance from clade root to target + distance from
		*  b to clade root are equal to the distance from b to the target
		*/
		if (NJ_FLT_EQ(dmat->val[NJ_MAP(target, b, dmat->size)],
			(clade_dist + b2clade))) {
			return(1);  /* join is legitimate   */
		}
		else {
			return(0);  /* join is illegitimate */
		}
	}
}

static inline
int
NJ_check(int RJ, DMAT *dmat,
	long int a,
	long int b,
	float min,
	int additivity) {


	long int i, size;
	float *ptr, *val, *r2;


	/* some aliases for speed and readability reasons */
	val = dmat->val;
	r2 = dmat->r2;
	size = dmat->size;


	/* now determine if joining a, b will result in broken distances */
	if (additivity) {
		if (!NJ_check_additivity(dmat, a, b)) {
			return(0);
		}
	}

	/* scan the horizontal of row b, punt if anything < min */
	ptr = &(val[NJ_MAP(b, b + 1, size)]);
	for (i = b + 1; i<size; i++) {
		if (NJ_FLT_LT((*ptr - (r2[b] + r2[i])), min)) {
			return(0);
		}
		ptr++;
	}

	/* scan the vertical component of row a, punt if anything < min */
	if (RJ == 1) {  /* if we are doing random joins, we checked this */
		ptr = val + a;
		for (i = 0; i<a; i++) {
			if (NJ_FLT_LT((*ptr - (r2[i] + r2[a])), min)) {
				return(0);
			}
			ptr += size - i - 1;
		}
	}

	/* scan the vertical component of row b, punt if anything < min */
	ptr = val + b;
	for (i = 0; i<b; i++) {
		if (NJ_FLT_LT((*ptr - (r2[i] + r2[b])), min) && i != a) {
			return(0);
		}
		ptr += size - i - 1;
	}

	return(1);
}


void
NJ_free_vertex(NJ_VERTEX *vertex) {

	if (vertex) {
		if (vertex->nodes_handle) {
			free(vertex->nodes_handle);
		}
		free(vertex);
	}

	return;
}

NJ_TREE *
NJ_relaxed_nj(int RJ, DMAT *dmat) {


	NJ_TREE *tree;
	NJ_VERTEX *vertex;
	long int a, b, t, bh, bv, i;
	float hmin, vmin, hvmin;
	float p, q, x;
	int join_flag;
	int additivity_mode;
	long int hmincount, vmincount;
	long int *permutation = NULL;



	/* initialize the r and r2 vectors */
	NJ_init_r(dmat);

	additivity_mode = 1;

	
	/* allocate and initialize our vertex vector used for tree construction */
	vertex = NJ_init_vertex(dmat);

	/* loop until there are only 2 nodes left to join */
	while (vertex->nactive > 2) {

		switch (RJ) {

			/* RANDOMIZED JOINS */
		case 0:

			join_flag = 0;

			NJ_permute(permutation, dmat->size - 1);
			for (i = 0; i<dmat->size - 1 && (vertex->nactive>2); i++) {

				a = permutation[i];

				/* find min trans dist along horiz. of row a */
				hmin = NJ_find_hmin(dmat, a, &bh, &hmincount);
				if (a) {
					/* find min trans dist along vert. of row a */
					vmin = NJ_find_vmin(dmat, a, &bv, &vmincount);
				}
				else {
					vmin = hmin;
					bv = bh;
					vmincount = 0;
				}

				if (NJ_FLT_EQ(hmin, vmin)) {

					/*
					* The minima along the vertical and horizontal are
					* the same.  Compute the proportion of minima along
					* the horizonal (p) and the proportion of minima
					* along the vertical (q).
					*
					* If the same minima exist along the horizonal and
					* vertical, we break the tie in a way which is
					* non-biased.  That is, we break the tie based on the
					* proportion of horiz. minima versus vertical minima.
					*
					*/
					p = (float)hmincount / ((float)hmincount + (float)vmincount);
					q = 1.0 - p;
					x = rand()/RAND_MAX;// genrand_real2();

					if (x < p) {
						hvmin = hmin;
						b = bh;
					}
					else {
						hvmin = vmin;
						b = bv;
					}
				}
				else if (NJ_FLT_LT(hmin, vmin)) {
					hvmin = hmin;
					b = bh;
				}
				else {
					hvmin = vmin;
					b = bv;
				}

				if (NJ_check(RJ, dmat, a, b, hvmin, additivity_mode)) {

					/* swap a and b, if necessary, to make sure a < b */
					if (b < a) {
						t = a;
						a = b;
						b = t;
					}

					join_flag = 1;

					/* join taxa from rows a and b */
					NJ_decompose(dmat, vertex, a, b, 0);

					/* collapse matrix */
					NJ_compute_r(dmat, a, b);
					NJ_collapse(dmat, vertex, a, b);

					NJ_permute(permutation, dmat->size - 1);
				}
			}

			/* turn off additivity if go through an entire cycle without joining */
			if (!join_flag) {
				additivity_mode = 0;
			}

			break;



			/* DETERMINISTIC JOINS */
		case 1:

			join_flag = 0;

			for (a = 0; a<dmat->size - 1 && (vertex->nactive > 2);) {

				/* find the min along the horizontal of row a */
				hmin = NJ_find_hmin(dmat, a, &b, &hmincount);

				if (NJ_check(RJ, dmat, a, b, hmin, additivity_mode)) {

					join_flag = 1;

					/* join taxa from rows a and b */
					NJ_decompose(dmat, vertex, a, b, 0);

					/* collapse matrix */
					NJ_compute_r(dmat, a, b);
					NJ_collapse(dmat, vertex, a, b);

					if (a) {
						a--;
					}

				}
				else {
					a++;
				}
			}

			/* turn off additivity if go through an entire cycle without joining */
			if (!join_flag) {
				additivity_mode = 0;
			}

			break;
		}

	}  /* WHILE */

	   /* Join the last two nodes on the vertex list */
	tree = NJ_decompose(dmat, vertex, 0, 1, NJ_LAST);

	if (vertex) {
		NJ_free_vertex(vertex);
	}

	

	return(tree);
}

void
NJ_output_tree2(int nlen, int *tpos, char *otreex,
	NJ_TREE *tree,
	NJ_TREE *root,
	DMAT *dmat) {


	int   x2, x3, n2, n3, mod;
	float x;

	if (!tree) {
		return;
	}

	if (tree->taxa_index != NJ_INTERNAL_NODE) {

		//print the taxon name (actually its number starting with an S)
		*tpos = *tpos + 1;
		*(otreex + *tpos) = 83;

		mod = nlen;
		n2 = tree->taxa_index;
		n3 = (int)(n2 / mod);
		*tpos = *tpos + 1;
		*(otreex + *tpos) = (char)(48 + n3);
		while (mod > 1) {

			n2 -= n3*mod;
			n3 = (int)(n2 / (mod / 10));
			*tpos = *tpos + 1;
			*(otreex + *tpos) = (char)(48 + n3);
			mod /= 10;
		}
		//Print the branchlength
		*tpos = *tpos + 1;
		*(otreex + *tpos) = 58;//first print the colon
		
		x = tree->dist;
		if (x < 0.0) {
			x *= -1;
			//*tpos = *tpos + 1;
			//*(otreex + *tpos) = 45;//a minus sign
		}
		if (x < 1.0) {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = 48;// a zero
			*tpos = *tpos + 1;
			*(otreex + *tpos) = 46;//a dot
		}
		else {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = (int)(x)+48; //a number>0
			*tpos = *tpos + 1;
			*(otreex + *tpos) = 46;//a dot
		}

		x -= (int)(x);

		mod = 100000;
		x2 = (int)(x*mod);
		x3 = (int)(x2 / (mod / 10));
		while (mod > 10) {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = char(48 + x3);
			mod /= 10;
			x2 -= x3*mod;
			x3 = x2 / (mod / 10);
			//if (x2/10000000 < 1 && (double)(x2/10000000) > 0.1){
			//	tpos++;
			//	*(otree + tpos) = 46;
			//}
		}

//			fprintf(fp, "%s:%f",
//				dmat->taxaname[tree->taxa_index],
//				tree->dist);//print length of terminal branch in normal notation
		

	}
	else {


		if (tree->left && tree->right) {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = char(40);//open bracket
			//fprintf(fp, "(");
		}
		if (tree->left) {
			NJ_output_tree2(nlen, tpos, otreex, tree->left, root, dmat);
		}

		if (tree->left && tree->right) {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = char(44);
			//fprintf(fp, ",");
		}
		if (tree->right) {
			NJ_output_tree2(nlen, tpos, otreex, tree->right, root, dmat);
		}

		if (tree != root->left) {
			if (tree->left && tree->right) {
				if (tree != root) {
					//print close bracket and colon
					*tpos = *tpos + 1;
					*(otreex + *tpos) = char(41);
					*tpos = *tpos + 1;
					*(otreex + *tpos) = char(58);
					//now print the branch length
					x = tree->dist;
					if (x < 0.0) {
						x *= -1;
						//*tpos = *tpos + 1;
						//*(otreex + *tpos) = 45;//a minus sign
					}
					if (x < 1.0) {
						*tpos = *tpos + 1;
						*(otreex + *tpos) = 48;// a zero
						*tpos = *tpos + 1;
						*(otreex + *tpos) = 46;//a dot
					}
					else {
						*tpos = *tpos + 1;
						*(otreex + *tpos) = (int)(x)+48; //a number>0
						*tpos = *tpos + 1;
						*(otreex + *tpos) = 46;//a dot
					}

					x -= (int)(x);

					mod = 100000;
					x2 = (int)(x*mod);
					x3 = (int)(x2 / (mod / 10));
					while (mod > 10) {
						*tpos = *tpos + 1;
						*(otreex + *tpos) = char(48 + x3);
						mod /= 10;
						x2 -= x3*mod;
						x3 = x2 / (mod / 10);
						//if (x2/10000000 < 1 && (double)(x2/10000000) > 0.1){
						//	tpos++;
						//	*(otree + tpos) = 46;
						//}
					}

					//fprintf(fp, "):%f", tree->dist);//close bracket and print branch length in normal notation
					
				}
				else {
					*tpos = *tpos + 1;
					*(otreex + *tpos) = char(41);
					//fprintf(fp, ")");
				}
			}
		}
		else {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = char(41);
			//fprintf(fp, ")");
		}
	}

	return;
}

void
NJ_search_tree(int nlen, int *tpos, char *otreex,
	NJ_TREE *tree,
	NJ_TREE *root,
	DMAT *dmat, NJ_TREE *target, int outlyer) {


	int   x2, x3, n2, n3, mod;
	float x;

	if (!tree) {
		return;
	}

	if (tree->taxa_index == outlyer) {

		target = tree;


	}
	else {


		
		if (tree->left) {
			NJ_search_tree(nlen, tpos, otreex, tree->left, root, dmat,target,outlyer);
		}

		
		if (tree->right) {
			NJ_search_tree(nlen, tpos, otreex, tree->right, root, dmat, target,outlyer);
		}

		
	}

	return;
}

int
NJ_output_tree(int nlen, int *tpos, NJ_TREE *tree,
	DMAT *dmat,
	long int count, char *outtree,int outlyer) {
	//find the outlyer node
	NJ_TREE *target;

	target = (NJ_TREE *)calloc(1, sizeof(NJ_TREE));
	if (outlyer > 0) {
		NJ_search_tree(nlen, tpos, outtree, tree, tree, dmat, target, outlyer-1);
		NJ_output_tree2(nlen, tpos, outtree, target, tree, dmat);
	}
	else
		NJ_output_tree2(nlen, tpos, outtree, tree, tree, dmat);
	
	*tpos = *tpos + 1;
	*(outtree + *tpos) = char(59);
	// add the final ";"
	//fprintf(fp, ";\n");

	
	return(0);
}

void
NJ_free_tree(NJ_TREE *node) {

	if (!node) {
		return;
	}

	if (node->left) {
		NJ_free_tree(node->left);
	}

	if (node->right) {
		NJ_free_tree(node->right);
	}

	free(node);

	return;
}

void
NJ_free_dmat(DMAT *dmat) {

	long int i;

	if (dmat) {

		if (dmat->taxaname) {

			for (i = 0; i<dmat->ntaxa; i++) {
				if (dmat->taxaname[i]) {
					free(dmat->taxaname[i]);
				}
			}

			free(dmat->taxaname);
		}

		if (dmat->valhandle) {
			free(dmat->valhandle);
		}

		if (dmat->rhandle) {
			free(dmat->rhandle);
		}

		if (dmat->r2handle) {
			free(dmat->r2handle);
		}

		free(dmat);
	}

	return;
}



//All this clearcut code is written by Luke Sheneman
/*
* clearcut.c
*
* $Id: clearcut.c,v 1.2 2006/08/25 03:58:45 sheneman Exp $
*
*****************************************************************************
*
* Copyright (c) 2004,  Luke Sheneman
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
*  + Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*  + Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in
*    the documentation and/or other materials provided with the
*    distribution.
*  + The names of its contributors may not be used to endorse or promote
*    products derived  from this software without specific prior
*    written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************
*
* An implementation of the Relaxed Neighbor-Joining algorithm
*  of Evans, J., Sheneman, L., and Foster, J.
*
*
* AUTHOR:
*
*   Luke Sheneman
*   sheneman@cs.uidaho.edu
*
*/
float MyMathFuncs::Clearcut(int outlyer, int NextNo, int treetype, int nlen2, int nseed, int RJ, int UBD, float *dists, char *outtree) {
	//treetype = 1 for normal NJ, 0 for rapidNJ
	//RJ = 0 for randomized joins and 1 for deterministic joins during relaxed NJ construction
	int lt, nlen;
	int tpos;
	DMAT *dmat;         /* The working distance matrix */
	//DMAT *dmat_backup = NULL;/* A backup distance matrix    */
	NJ_TREE *tree;      /* The phylogenetic tree       */
	
	long int i;
	tpos = 0;
	/* some variables for tracking time */
	//struct timeval tv;
	unsigned long long startUs, endUs;


	/* check and parse supplied command-line arguments */
	//nj_args = NJ_handle_args(argc, argv); -these now get passed through the function call

	/* Initialize Mersenne Twister PRNG */
	//init_genrand(seed); already initialised
	srand(nseed);


	//switch (nj_args->input_mode) {

		/* If the input type is a distance matrix */
	///case NJ_INPUT_MODE_DISTANCE: input type always distances

		/* parse the distance matrix */
		dmat = NJ_parse_distance_matrix(dists, NextNo);
		

		/* If the input type is a multiple sequence alignment */
	
	//}

	/*
	* If we are going to generate multiple trees from
	* the same distance matrix, we need to make a backup
	* of the original distance matrix.
	*/
//	if (nj_args->ntrees > 1) {
//		dmat_backup = NJ_dup_dmat(dmat);
//	}

	/* process n trees */
		i = 0;
		

		/* RECORD THE PRECISE TIME OF THE START OF THE NEIGHBOR-JOINING */
		//gettimeofday(&tv, NULL);
		startUs = 10;//((unsigned long long) tv.tv_sec * 1000000ULL)
					 //  + ((unsigned long long) tv.tv_usec);


					 /*
					 * Invoke either the Relaxed Neighbor-Joining (treetype=0)
					 * or the "traditional" Neighbor-Joining algorithm(treetype=1)
					 */
		if (treetype == 1) {
			tree = NJ_neighbor_joining(dmat, outlyer);
		}
		else {
			tree = NJ_relaxed_nj(RJ, dmat);
		}

		

		/* RECORD THE PRECISE TIME OF THE END OF THE NEIGHBOR-JOINING */
		// gettimeofday(&tv, NULL);
		endUs = 20;// ((unsigned long long) tv.tv_sec * 1000000ULL)
				   //  + ((unsigned long long) tv.tv_usec);

		if (NextNo < 100)
			nlen = 10;
		else if (NextNo < 1000)
			nlen = 100;
		else
			nlen = 1000;

		/* Output the neighbor joining tree here */
		NJ_output_tree(nlen, &tpos, tree, dmat, i, outtree,outlyer);

		NJ_free_tree(tree);  /* Free the tree */
		NJ_free_dmat(dmat);  /* Free the working distance matrix */

		lt = tpos;

	

	return(lt);
}

}//end of MyMathFuncs

